"""Check native FW68D driver acceptance and cache isolation without running CFD."""

import csv
import json
import math
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
driver = runpy.run_path(str(ROOT / "cases/02d_fw68d_steady_six/run_steady.py"))
Campaign = driver["Campaign"]
np = driver["np"]


def write_solution(run, case, clearance, force_error=0., marker="CONVERGED"):
    """Write a small real VTK/diagnostics fixture with a solved non-seed position."""
    ambient = float(case["tc0_boundary_temperature_k"])
    load = float(case["load_mean_n"])
    psi = math.radians(320.)
    position = [-.45 * clearance * math.cos(psi), -.45 * clearance * math.sin(psi), 0.]
    fields = {"p": [101325., 101325.], "film_content": [1., 1.],
              "h": [20e-6, 30e-6], "T": [ambient + 1] * 2,
              "T_bearing": [ambient + 1] * 2, "T_journal": [ambient + 1] * 2,
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
                  "coupling_R_p": 0., "coupling_R_T": 0., "coupling_converged": 1,
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
        load = float(case["load_mean_n"])
        clearance = float(campaign.experiment["geometry"]["radial_clearance_m"])
        assert np.array_equal(driver["applied_load"](load), [0., load])
        run = work / "E/native"
        log = write_solution(run, case, clearance)
        def check(log):
            return driver["check_solution"](run, log, 0, load, 0.,
                float(case["tc0_boundary_temperature_k"]), campaign.points, True, clearance)
        result = check(log)
        assert math.isclose(result["eps"], .45)
        assert math.isclose(result["psi_deg"], 320.)
        assert result["force_error_rel"] == 0.
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
    print("FW68D native driver passed: VTK position/force acceptance, single equilibrium call, provenance and wrong-load restart rejection")


if __name__ == "__main__":
    main()
