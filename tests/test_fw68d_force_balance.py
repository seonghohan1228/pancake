"""Check native FW68D driver acceptance and cache isolation without running CFD."""

import csv
import json
import math
import os
from pathlib import Path
import runpy
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
driver = runpy.run_path(str(ROOT / "cases/02d_fw68d_steady_six/run_steady.py"))
Campaign = driver["Campaign"]
np = driver["np"]


def write_solution(run, case, clearance, force_error=0., marker="CONVERGED", pressure_change=0.):
    """Write a small real VTK/diagnostics fixture with a solved non-seed position."""
    ambient = float(case["tc0_boundary_temperature_k"])
    load = float(case["load_mean_n"])
    psi = math.radians(320.)
    position = [-.45 * clearance * math.cos(psi), -.45 * clearance * math.sin(psi), 0.]
    fields = {"p": [101325., 101325.], "film_content": [1., 1.],
              "h": [20e-6, 30e-6], "T": [ambient + 1] * 2,
              "T_bearing": [ambient + 1] * 2, "T_journal": [ambient + 1] * 2,
              "temperature_housing_outer": [ambient + 1] * 2,
              "mu": [.02, .02], "q_bearing": [0., 0.], "q_journal": [0., 0.],
              "friction_torque": [.01, .01],
              "F": [load * force_error, -load, 0.] * 2, "xB": position * 2}
    flat = run / "results/flat"
    flat.mkdir(parents=True, exist_ok=True)
    vtk = ET.Element("VTKFile")
    cell = ET.SubElement(ET.SubElement(vtk, "Piece"), "CellData")
    for name, values in fields.items():
        ET.SubElement(cell, "DataArray", Name=name).text = " ".join(map(str, values))
    ET.ElementTree(vtk).write(flat / "rank0.vts")
    (flat / "frame.pvts").write_text('<VTKFile><Piece Source="rank0.vts"/></VTKFile>')
    (flat / "results.pvd").write_text('<VTKFile><DataSet timestep="0" file="frame.pvts"/></VTKFile>')
    diagnostic = {"liquid_residual": 0., "cavitated_fraction": 0.,
                  "liquid_residual_kg_s": 0., "liquid_balance_scale_kg_s": 1.,
                  "liquid_residual_rel": 0., "energy_balance_residual_W": 0.,
                  "energy_balance_scale_W": 1., "energy_viscous_heat_W": 0.,
                  "energy_wall_heat_W": 0., "energy_axial_enthalpy_out_W": 0.,
                  "energy_groove_enthalpy_in_W": 0., "energy_mass_imbalance_correction_W": 0.,
                  "energy_feed_mass_kg_s": 0., "cross_film_profile_t_min_K": ambient + 1,
                  "cross_film_profile_t_max_K": ambient + 2, "cross_film_profile_residual": 0.,
                  "coupling_R_p": pressure_change, "coupling_R_T": 0., "coupling_converged": 1,
                  "coupling_restored_best": 0, "converged": 1}
    with (run / "results/diagnostics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostic))
        writer.writeheader()
        writer.writerow(diagnostic)
    log = (f"fbns iter=0 |G|=1e-13 |FB|=1e-13\n"
           f"Equilibrium: {marker} R_F=0.001 x={position[0]} y={position[1]}\n"
           f"Step 0 | F=(0, {-load}, 0) N | M=0.01\n")
    (run / "run.log").write_text(log)
    return log


def main():
    with TemporaryDirectory() as directory:
        work = Path(directory)
        args = SimpleNamespace(materials=driver["SETUP"], solver=None, work=directory,
                               ranks=1, coupling_max_iters=2000, coupling_relaxation=1.,
                               trial=None, rankfile=None, max_thermal_iterations=30)
        campaign = Campaign(args)
        case = campaign.cases["E"]
        # Separate physical supply and cooling, including cache/restart identities.
        def parsed(text):
            return dict(line.split(" = ", 1) for line in text.splitlines() if " = " in line)
        baseline = parsed(campaign.config_for(case, .15, 345., 60., 100., work))
        amb = case["tc0_boundary_temperature_k"]
        for keys in (("bc_z_south_temperature", "bc_z_north_temperature", "groove_supply_temperature"),
                     ("bush_outer_temperature", "shaft_ambient_temperature"),
                     ("bc_z_south_temperature",), ("bc_z_north_temperature",),
                     ("groove_supply_temperature",), ("bush_outer_temperature",),
                     ("shaft_ambient_temperature",)):
            changed_case = {**case, "thermal_boundaries_K": dict.fromkeys(keys, amb + 1)}
            config_text = campaign.config_for(changed_case, .15, 345., 60., 100., work)
            changed = parsed(config_text)
            expected = {"inlet_groove" if key == "groove_supply_temperature" else key for key in keys}
            assert {key for key in baseline if baseline[key] != changed[key]} == expected
            assert driver["config_identity"](config_text) != driver["config_identity"](
                campaign.config_for(case, .15, 345., 60., 100., work))
            assert driver["warm_start_identity"](config_text) != driver["warm_start_identity"](
                campaign.config_for(case, .15, 345., 60., 100., work))
        for invalid in ({"unknown": amb}, {"bush_outer_temperature": float("nan")},
                        {"bush_outer_temperature": -1}, {"bush_outer_temperature": True}):
            try:
                driver["boundary_temperatures"]({**case, "thermal_boundaries_K": invalid})
            except ValueError:
                pass
            else:
                raise AssertionError("invalid thermal boundary accepted")
        # Distinct baths exercise the actual coefficient closure, not just rendering.
        bath_case = {**case, "thermal_boundaries_K": {
            "bush_outer_temperature": amb + 1, "shaft_ambient_temperature": amb + 2}}
        bush_target = driver["housing_coefficient"](amb + 6, amb + 1, campaign.points, campaign.materials)[0]
        shaft_target = driver["shaft_coefficient"](amb + 7, amb + 2,
            case["angular_speed_magnitude_rad_s"], campaign.points, campaign.materials)[0]
        probe = Campaign.__new__(Campaign)
        probe.args = SimpleNamespace(max_thermal_iterations=2)
        probe.native_equilibrium = True
        probe.h_warm = {}
        probe.points, probe.materials = campaign.points, campaign.materials
        probe.find_completed = lambda *a, **kw: None
        def fake_single(c, eps, psi, hb, hs, label):
            assert math.isclose(hb, bush_target) and math.isclose(hs, shaft_target)
            folder = work / label
            folder.mkdir()
            return {"run_dir": str(folder), "eps": .3, "psi_deg": 270.,
                    "T_outer_mean_K": amb + 6, "T_shaft_mean_K": amb + 7}
        probe.single = fake_single
        closed = probe.evaluate(bath_case, .15, 345., "bath_check")
        assert closed["status"] == "load_and_thermal_closed" and closed["h_relative_mismatch"] == 0
        # The reporting split reproduces analytic uniform-film torque and axial flux.
        sys.path.insert(0, str(ROOT / "cases/02d_fw68d_steady_six"))
        from summarize_sensitivity import field_metrics, execution_rankfile
        rankfile = work / "parallel.rankfile"
        rankfile.write_text("rank 0=afml-r6 slot=0:4\n")
        job = {"scenario":"cooling_plus_1K", "case":"E", "status":"completed", "returncode":0,
               "rankfile":str(rankfile), "rankfile_sha256":driver["sha256"](rankfile)}
        assert execution_rankfile({"jobs":[job]}, "cooling_plus_1K", "E") == rankfile
        for jobs in ([], [job,job], [{**job,"status":"running"}], [{**job,"returncode":1}],
                     [{**job,"rankfile_sha256":"wrong"}]):
            try:
                execution_rankfile({"jobs":jobs}, "cooling_plus_1K", "E")
            except ValueError:
                pass
            else:
                raise AssertionError("missing, failed, ambiguous or changed parallel execution accepted")
        nt, nz, nr = 8, 2, 2
        h, mu, rho = 40e-6, .1, 930.
        radius, length, omega = .0215, .04, 20.
        shape = (nz, nt)
        uniform = lambda value: np.full(shape, value).ravel()
        analytical = 2 * np.pi * mu * omega * radius**3 * length / h
        raw = {"theta_rad": np.tile((np.arange(nt)+.5)*2*np.pi/nt, nz),
               "z_m": np.repeat((np.arange(nz)+.5)*length/nz, nt),
               "p": uniform(101335.), "film_content": uniform(1.),
               "rho_liquid_solution": uniform(rho), "cross_film_f0": uniform(h/mu),
               "cross_film_f1": uniform(h*h/(2*mu)), "cross_film_f2": uniform(h**3/(12*mu)),
               "cross_film_couette_h": uniform(h), "friction_torque": uniform(analytical),
               "T": uniform(amb+5), "mu": uniform(mu), "q_bearing": uniform(0.),
               "q_journal": uniform(0.), "T_journal": uniform(amb+5),
               "temperature_housing_outer": uniform(amb+5),
               **{f"cross_film_solid_{k}": uniform(amb+5) for k in range(nr)}}
        cfg = {**baseline, "n_theta_global": str(nt), "n_z_global": str(nz),
               "bush_radial_cells": str(nr), "bush_axial_cells": str(nz),
               "omega": str(omega), "inlet_groove": f"90 0 101325 {amb}"}
        split = field_metrics(cfg, raw, campaign.points, campaign.materials)
        expected_out = 2*radius*2*np.pi/(length/nz)*rho*h**3/(12*mu)*10
        assert math.isclose(split["couette_Nm"], analytical, rel_tol=1e-12)
        assert abs(split["pressure_shear_Nm"]) < 1e-14
        assert math.isclose(split["south_out_kg_s"], expected_out, rel_tol=1e-12)
        assert math.isclose(split["north_out_kg_s"], expected_out, rel_tol=1e-12)
        assert split["north_in_kg_s"] == split["groove_in_kg_s"] == 0
        load = float(case["load_mean_n"])
        clearance = float(campaign.experiment["geometry"]["radial_clearance_m"])
        assert np.array_equal(driver["applied_load"](load), [0., load])
        run = work / "E/native"
        log = write_solution(run, case, clearance)
        def check(log):
            return driver["check_solution"](run, log, 0, load, 0.,
                float(case["tc0_boundary_temperature_k"]), campaign.points, True, clearance)
        cold_case = {**case, "tc0_boundary_temperature_k": amb - 3}
        cold_log = write_solution(run, cold_case, clearance)
        cold_result = driver["check_solution"](run, cold_log, 0, load, 0., amb - 3,
                                                campaign.points, True, clearance)
        assert cold_result["T_min_K"] < amb
        try:
            check(cold_log)
        except ValueError:
            pass
        else:
            raise AssertionError("sub-reservoir temperature accepted")
        log = write_solution(run, case, clearance)
        result = check(log)
        assert math.isclose(result["eps"], .45)
        assert math.isclose(result["psi_deg"], 320.)
        assert result["force_error_rel"] == 0.
        # Inner pressure acceptance does not add a stricter pressure-change gate.
        check(write_solution(run, case, clearance, pressure_change=3e-4))
        for invalid_log in (log.replace("|G|=1e-13", "|G|=1.1e-6"),
                            log.replace("|FB|=1e-13", "|FB|=1.1e-6")):
            try:
                check(invalid_log)
            except ValueError:
                pass
            else:
                raise AssertionError("failed inner pressure equation accepted")
        for force_error, marker in ((.003, "CONVERGED"), (0., "NOT CONVERGED")):
            bad_log = write_solution(run, case, clearance, force_error, marker)
            try:
                check(bad_log)
            except ValueError:
                pass
            else:
                raise AssertionError("native failure or independent VTK force defect accepted")
        log = write_solution(run, case, clearance)
        config = campaign.config_for(case, .15, 345., 60., 100., run / "results")
        assert "motion_model = STEADY_EQUILIBRIUM" in config
        assert "bearing_x,bearing_y,bearing_z" in config
        (run / "config.txt").write_text(config)
        record = {"case_id": "E", "returncode": 0, "status": "valid_inner_solve",
                  "source_sha256": campaign.source_hash, "executable_sha256": None,
                  "driver_sha256": campaign.driver_hash, "rankfile_sha256": None,
                  "native_equilibrium": True, "mpi_ranks": 1, **campaign.input_hashes,
                  "config_sha256": driver["sha256"](run / "config.txt"),
                  "seed_eps": .15, "seed_psi_deg": 345., **result}
        assert campaign.matching_provenance(record)
        checked = campaign.recheck(run, record, case)
        assert checked is not None and math.isclose(checked["eps"], .45)
        for key in ("driver_sha256", "rankfile_sha256", "native_equilibrium"):
            old = dict(record)
            old.pop(key)
            assert not campaign.matching_provenance(old), key
        # Preserve the original wrong-load regression with current provenance.
        old_run = work / "E/old_downward_load"
        (old_run / "results").mkdir(parents=True)
        old_config = config.replace("external_load_y = " + str(load), "external_load_y = " + str(-load))
        (old_run / "config.txt").write_text(old_config)
        (old_run / "results/steady_state.dat").write_text("old orientation state")
        old = {**record, "config_sha256": driver["sha256"](old_run / "config.txt")}
        campaign.completed["E"] = [(old_run, old)]
        assert campaign.recheck(old_run, old, case) is None
        assert campaign.find_completed(case, .15, 345., 60., 100.) is None
        assert campaign.recover_checkpoint(case) is None
        # A preserved failed result can seed a fresh solve, but never satisfy its acceptance.
        (run / "results/steady_state.dat").write_text("initial guess only")
        (run / "metrics.json").write_text(json.dumps({**record, "status": "failed",
                                                     "source_sha256": "previous", "driver_sha256": "previous"}))
        restart_args = vars(args).copy()
        restart_args.update(work=str(work / "fresh"), initial_state_root=str(work))
        warm = Campaign(SimpleNamespace(**restart_args))
        assert warm.find_completed(case, .15, 345., 60., 100.) is None
        assert warm.recover_checkpoint(case) == run / "results/steady_state.dat"
        assert warm.restart_positions["E"] == (result["eps"], result["psi_deg"])
        assert warm.restart_provenance["E"]["initial_guess_only"]
        assert warm.restart_provenance["E"]["prior_status"] == "failed"
        assert warm.h_warm["E"] == (60., 100.)
        assert not warm.completed["E"]
        calls = []
        def evaluate(case, eps, psi, label):
            calls.append((eps, psi, label))
            return {**result, "status": "load_and_thermal_closed", "run_dir": str(run)}
        mock = SimpleNamespace(work=work, evaluate=evaluate)
        balanced = Campaign.balance(mock, case, .15, 345.)
        assert calls == [(.15, 345., "equilibrium")]
        assert balanced["status"] == "load_and_thermal_closed"
        saved = list((work / "E").glob("result_*.json"))
        assert len(saved) == 1 and math.isclose(json.loads(saved[0].read_text())["eps"], .45)
        stop_run = work / "stop_check"
        stop_run.mkdir()
        child = ("from pathlib import Path\nimport time\n"
                 "while not Path('results/STOP').exists(): time.sleep(.01)\n"
                 "Path('results/steady_state.dat').write_text('saved at stop')\nraise SystemExit(5)\n")
        returncode, deadline = driver["run_solver"]([sys.executable, "-c", child], stop_run,
                                                    dict(os.environ), timeout=.1, stop_grace=10)
        assert deadline and returncode == 5
        assert (stop_run / "results/steady_state.dat").read_text() == "saved at stop"
    print("FW68D native driver passed: force/inner-pressure acceptance, cache isolation, warm-start-only recovery and cooperative timeout")


if __name__ == "__main__":
    main()
