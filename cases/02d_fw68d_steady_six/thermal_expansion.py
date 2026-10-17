"""Free-expansion geometry, solid-volume temperatures and preserved A-F campaigns.

Run one case per disjoint OpenMPI rankfile. Every geometry starts with fresh
fields; only the native driver's same-geometry cooling loop reuses checkpoints.
"""

import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
import run_steady as driver


def expanded_geometry(shaft_c, bearing_c):
    """Expand both bodies from their 20 C dimensions, never from a previous solve."""
    if any(isinstance(t, bool) or not isinstance(t, (int, float)) or
           not math.isfinite(t) or t <= -273.15 for t in (shaft_c, bearing_c)):
        raise ValueError("solid temperatures must be finite Celsius above absolute zero")
    shaft_scale = 1 + 16e-6 * (shaft_c - 20)
    bearing_scale = 1 + 16e-6 * (bearing_c - 20)
    ds, db = .043 * shaft_scale, .043085 * bearing_scale
    geometry = {"journal_radius_m": ds / 2, "radial_clearance_m": (db - ds) / 2,
                "bearing_length_m": .04 * bearing_scale,
                "bush_outer_radius_m": .0425 * bearing_scale,
                "sensor_radius_m": .026 * bearing_scale}
    driver.validate_geometry(geometry)
    return geometry


def starting_geometry(row):
    """Use the supplied full-precision dimensions without rerounding their inputs."""
    return dict(zip(("journal_radius_m", "radial_clearance_m", "bearing_length_m",
                     "bush_outer_radius_m", "sensor_radius_m"),
                    (row["Rshaft_m"], row["clearance_m"], row["Lbearing_m"],
                     row["outerR_m"], row["sensor_radius_m"])))


def geometry_change(old, new):
    """Report proposed updates in metres and apply the specified numerical gates."""
    def dimensions(g):
        return {"Ds": 2*g["journal_radius_m"],
                "Db": 2*(g["journal_radius_m"]+g["radial_clearance_m"]),
                "L": g["bearing_length_m"], "Do": 2*g["bush_outer_radius_m"],
                "c": g["radial_clearance_m"]}
    a, b = dimensions(old), dimensions(new)
    changes = {key: b[key]-a[key] for key in a}
    return {"delta_m": changes, "converged": all(abs(v) <= (1e-9 if k == "c" else 1e-8)
                                                  for k, v in changes.items())}


def solid_observables(config, raw, sensor_radius):
    """Verify repeated solid slabs and sample moving sensors in log-radius."""
    nt, nz, nr, ns = (int(config[k]) for k in
                      ("n_theta_global", "n_z_global", "bush_radial_cells", "bush_axial_cells"))
    if nz % ns:
        raise ValueError("solid averaging requires equally repeated axial slabs")
    order = np.lexsort((np.asarray(raw["theta_rad"]).ravel(), np.asarray(raw["z_m"]).ravel()))
    solid = np.stack([np.asarray(raw[f"cross_film_solid_{k}"]).ravel()[order].reshape(nz, nt)
                      for k in range(nr)])
    repeated = solid.reshape(nr, ns, nz // ns, nt)
    if not np.isfinite(solid).all() or not np.array_equal(repeated, np.repeat(repeated[:, :, :1], nz//ns, axis=2)):
        raise ValueError("invalid solid temperatures or changed axial-slab mapping")
    shaft = np.asarray(raw["T_journal"])
    if not np.isfinite(shaft).all() or np.ptp(shaft) > 1e-10:
        raise ValueError("shaft is not isothermal")
    radius, outer = float(config["R"]), float(config["bush_outer_radius"])
    edges = np.linspace(radius, outer, nr+1)
    weights = np.diff(edges**2) / (outer**2-radius**2)
    centers = .5 * (edges[:-1]+edges[1:])
    if not centers[0] <= sensor_radius <= centers[-1]:
        raise ValueError("sensor interpolation would extrapolate solid temperatures")
    mid = .5*(solid[:, nz//2-1]+solid[:, nz//2])
    curve = [np.interp(np.log(sensor_radius), np.log(centers), mid[:, i]) for i in range(nt)]
    angles = np.degrees(np.asarray(raw["theta_rad"]).ravel()[order].reshape(nz, nt)[0])
    sensors = np.interp([240, 270, 300], angles, curve, period=360)-273.15
    return {"Ts_C": float(shaft.mean()-273.15),
            "Tb_C": float(weights @ repeated[:, :, 0].mean(axis=(1, 2))-273.15),
            "TC_C": sensors.tolist(), "sensor_radius_m": sensor_radius,
            "sensor_z_m": float(config["L"])/2,
            "solid_shell_edges_m": edges.tolist(), "solid_shell_weights": weights.tolist(),
            "solid_axial_repeats": nz//ns}


def parse_config(text):
    return dict(line.split(" = ", 1) for line in text.splitlines() if " = " in line)


def solve(args, label, geometry, tight, seed):
    """Run and independently recheck one complete load/thermal/cooling equilibrium."""
    destination = Path(args.work) / args.case / label
    destination.mkdir(parents=True, exist_ok=False)
    setup = json.loads(driver.SETUP.read_text())
    setup["additional_cases"] = []
    if geometry is not None:
        driver.validate_geometry(geometry)
        setup["geometry_by_case"] = {args.case: geometry}
    if tight:
        setup["solver_tolerances"] = {"steady_force_tolerance": 2e-5,
                                       "coupling_tolerance": 1e-8, "fbns_tol": 1e-8}
        setup["cooling_relative_tolerance"] = 1e-5
        setup["cooling_temperature_tolerance_K"] = 1e-4
    setup["force_balance"]["initial_positions"][args.case] = list(seed)
    setup_path = destination / "setup.json"
    driver.save_json(setup_path, setup)
    nominal_config = json.loads((Path(args.baseline) / f"{args.case}_config.json").read_text())
    options = SimpleNamespace(solver=args.solver, work=str(destination / "runs"),
        materials=str(setup_path), ranks=4, rankfile=args.rankfile, trial=None,
        timeout=7200, max_thermal_iterations=30, coupling_max_iters=int(nominal_config["coupling_max_iters"]),
        coupling_relaxation=1.0, initial_state_root=None)
    campaign = driver.Campaign(options)
    campaign.mpi_version = subprocess.check_output([campaign.command()[0], "--version"], text=True).splitlines()[:2]
    start = time.monotonic()
    result = campaign.balance(campaign.cases[args.case], *seed)
    run = Path(result["run_dir"])
    checked = campaign.recheck(run, result, campaign.cases[args.case])
    if checked is None:
        raise ValueError("completed result failed independent provenance/solution recheck")
    config = parse_config((run / "config.txt").read_text())
    raw = driver.read_fields(run)
    observables = solid_observables(config, raw, .026 if geometry is None else geometry["sensor_radius_m"])
    requested = setup.get("solver_tolerances", {})
    if (result["force_error_rel"] >= requested.get("steady_force_tolerance", .002) or
        float(result["diagnostics"]["coupling_R_T"]) >= requested.get("coupling_tolerance", 1e-6) or
        result["cross_film_profile_residual"] >= requested.get("coupling_tolerance", 1e-6) or
        result["energy_balance_residual_rel"] >= requested.get("coupling_tolerance", 1e-6)):
        raise ValueError("completed solution failed its requested tighter tolerance")
    summary = {"label": label, "tight": tight, "geometry": geometry,
               "elapsed_seconds": time.monotonic()-start, "config": config,
               "observables": observables, "result": result,
               "proposed_geometry": expanded_geometry(observables["Ts_C"], observables["Tb_C"])}
    if geometry is not None:
        summary["feedback"] = geometry_change(geometry, summary["proposed_geometry"])
    np.savez_compressed(destination / "fields.npz", **raw)
    driver.save_json(destination / "accepted.json", summary)
    print(f"{args.case} {label}: torque={result['torque_Nm']:.10g}, Ts={observables['Ts_C']:.8f}, Tb={observables['Tb_C']:.8f}", flush=True)
    return summary


def run_case(args):
    rows = json.loads(Path(args.geometry_input).read_text())["cases"]
    row = next(row for row in rows if row["case"] == args.case)
    geometry = starting_geometry(row)
    nominal = json.loads((Path(args.baseline) / f"{args.case}_config.json").read_text())
    seed = (float(nominal["e"])/float(nominal["c"]), float(nominal["attitude_angle_deg"]))
    record = {"case": args.case, "geometry_input_sha256": driver.sha256(args.geometry_input),
              "worker_sha256": driver.sha256(__file__), "status": "running", "solves": []}
    root = Path(args.work) / args.case
    root.mkdir(parents=True, exist_ok=False)
    def save():
        driver.save_json(root / "history.json", record)
    save()
    try:
        baseline = solve(args, "nominal_tight", None, True, seed)
        record["solves"].append(baseline)
        save()
        updates = 0
        tight = False
        while True:
            label = f"geometry_{updates:02d}" + ("_tight" if tight else "")
            current = solve(args, label, geometry, tight, seed)
            record["solves"].append(current)
            save()
            seed = current["result"]["eps"], current["result"]["psi_deg"]
            if not current["feedback"]["converged"] and updates < 5:
                geometry = current["proposed_geometry"]
                updates += 1
                continue
            if not tight:
                tight = True
                continue
            record["status"] = "complete" if current["feedback"]["converged"] else "feedback_unresolved_at_cap"
            record["geometry_updates"] = updates
            save()
            return
    except Exception as error:
        record["status"] = "failed"
        record["error"] = str(error)
        save()
        raise


def self_check():
    zero = expanded_geometry(20., 20.)
    assert math.isclose(zero["radial_clearance_m"], 42.5e-6, abs_tol=1e-17)
    for t in (20., 45., 80.):
        g = expanded_geometry(t, t)
        assert math.isclose(g["radial_clearance_m"], 42.5e-6*(1+16e-6*(t-20)), abs_tol=1e-17)
        assert math.isclose(2*(g["journal_radius_m"]+g["radial_clearance_m"]), .043085*(1+16e-6*(t-20)), abs_tol=1e-17)
    for bad in ((float("nan"), 30.), (True, 30.), (-274., 30.), (1000., 20.)):
        try:
            expanded_geometry(*bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted invalid temperatures {bad}")
    assert geometry_change(zero, zero)["converged"]
    setup = json.loads(driver.SETUP.read_text())
    exp = json.loads(driver.EXPERIMENT.read_text())
    case = copy.deepcopy(exp["cases"][0])
    g = expanded_geometry(65.7, 65.1)
    case["geometry"] = g
    text = driver.render_config(case, exp, setup["specified_boundaries"], setup["material_settings"],
                                .6, 320., 60., 100., Path("test-results"), driver.TABLE)
    cfg = parse_config(text)
    assert float(cfg["R"]) == g["journal_radius_m"] and float(cfg["c"]) == g["radial_clearance_m"]
    assert float(cfg["L"]) == g["bearing_length_m"] and float(cfg["e"]) == .6*g["radial_clearance_m"]
    assert float(cfg["bush_outer_radius"]) == g["bush_outer_radius_m"]
    assert float(cfg["shaft_ambient_area_ratio"]) == g["journal_radius_m"]/g["bearing_length_m"]
    assert float(cfg["inlet_groove"].split()[1]) == setup["specified_boundaries"]["groove"]["width_deg"]
    material = driver.case_materials(case, setup["material_settings"])
    assert material["shaft_exposed_end_radius_m"] == g["journal_radius_m"]
    assert material["housing_outer_diameter_m"] == 2*g["bush_outer_radius_m"]
    assert math.isclose(g["sensor_radius_m"]/.026, g["bearing_length_m"]/.04)
    # A nonuniform synthetic solid distinguishes volume weighting from a surface proxy.
    cfg.update(n_theta_global="4", n_z_global="4", bush_radial_cells="2", bush_axial_cells="2")
    raw = {"theta_rad": np.tile(np.arange(4)*np.pi/2, 4), "z_m": np.repeat(np.arange(4), 4),
           "T_journal": np.full(16, 310.), "cross_film_solid_0": np.full(16, 300.),
           "cross_film_solid_1": np.full(16, 320.)}
    sensor = g["journal_radius_m"]+.3*(g["bush_outer_radius_m"]-g["journal_radius_m"])
    obs = solid_observables(cfg, raw, sensor)
    edges = np.linspace(float(cfg["R"]), float(cfg["bush_outer_radius"]), 3)
    expected = np.diff(edges**2) @ np.array([300., 320.])/np.diff(edges[[0,-1]]**2)[0]-273.15
    assert math.isclose(obs["Tb_C"], expected, abs_tol=1e-12)
    raw["cross_film_solid_0"][4] += 1
    try:
        solid_observables(cfg, raw, sensor)
    except ValueError:
        pass
    else:
        raise AssertionError("changed axial mapping accepted")
    print("thermal expansion self-check passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--case", choices=list("ABCDEF"))
    parser.add_argument("--geometry-input")
    parser.add_argument("--baseline")
    parser.add_argument("--work")
    parser.add_argument("--solver")
    parser.add_argument("--rankfile")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        if not all((args.case, args.geometry_input, args.baseline, args.work, args.solver, args.rankfile)):
            parser.error("case, geometry-input, baseline, work, solver and rankfile are required")
        run_case(args)
