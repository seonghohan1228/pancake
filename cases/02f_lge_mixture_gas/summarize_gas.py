"""Recheck and compare matched native gas branches, preserving experimental coordinates."""

import argparse
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_gas as worker


def accepted(root, key, branch, phase="production"):
    """Revalidate numerical acceptance and provenance; never promote a partial run."""
    destination = root / phase / branch / key
    saved = json.loads((destination / "accepted.json").read_text())
    tight, layers = phase.startswith("refinement"), 32 if phase == "refinement_layers32" else 16
    if (saved["status"] != "accepted_conditional_native_coupled" or saved["gas_branch"] != branch or
            saved["tight"] != tight or saved["cross_film_layers"] != layers or
            saved["worker_sha256"] != worker.driver.sha256(worker.__file__)):
        raise ValueError("accepted scenario metadata differs")
    expected_setup = copy.deepcopy(json.loads((root / "inputs" / f"{key}_{branch}_setup.json").read_text()))
    expected_setup["material_settings"]["cross_film_layers"] = layers
    if tight:
        expected_setup["solver_tolerances"] = dict(steady_force_tolerance=2e-5, coupling_tolerance=1e-8, fbns_tol=1e-12)
        expected_setup["cooling_relative_tolerance"] = 1e-5
        expected_setup["cooling_temperature_tolerance_K"] = 1e-4
    if expected_setup != json.loads((destination / "setup.json").read_text()):
        raise ValueError("case setup differs from frozen matched inputs")
    options = worker.baseline.options(root / "inputs", destination, None, None)
    command = saved["result"]["command"]
    options.solver = command[command.index("--config")-1]
    if not Path(options.solver).is_file():
        options.solver = str(root / "source/build-linux/pancake")
    options.rankfile = None
    options.table = str(root / "inputs/fw68d.csv")
    options.ranks = saved["result"]["mpi_ranks"]
    campaign = worker.Campaign(options)
    # Rehash the synchronized rankfile when its original absolute path is unavailable.
    if "--rankfile" in command:
        rankfile = Path(command[command.index("--rankfile")+1])
        if not rankfile.is_file():
            rankfile = root / rankfile.name
        campaign.rankfile_hash = worker.driver.sha256(rankfile)
    run = Path(saved["result"]["run_dir"])
    if not run.is_dir():
        run = destination / "runs" / key / run.name
    case = campaign.cases[key]
    checked = campaign.recheck(run, saved["result"], case)
    if checked is None:
        raise ValueError("native result provenance or numerical recheck failed")
    metrics = json.loads((run / "metrics.json").read_text())
    if any(saved["result"].get(name) != value for name, value in metrics.items()):
        raise ValueError("accepted result differs from retained solver metrics")
    if worker.driver.sha256(destination / "fields.npz") != saved["fields_sha256"]:
        raise ValueError("accepted field snapshot hash changed")
    config = worker.parse_config((run / "config.txt").read_text())
    if saved["config"] != config or saved["case"] != case["measured"]:
        raise ValueError("accepted config or authoritative measured row differs")
    expected = worker.parse_config(campaign.config_for(case, checked["seed_eps"], checked["seed_psi_deg"],
        checked["h_bush_W_m2_K"], checked["h_shaft_W_m2_K"], run / "results"))
    for name in expected.keys() | config.keys():
        if name not in ("output_dir", "solution_table_file") and expected.get(name) != config.get(name):
            raise ValueError(f"accepted configuration control differs: {name}")
    if checked["status"] != "load_and_thermal_closed":
        raise ValueError("saved result is not load and cooling closed")
    ambient = float(case["tc0_boundary_temperature_k"])
    material = worker.driver.case_materials(case, campaign.materials)
    hb, _ = worker.driver.housing_coefficient(checked["T_outer_mean_K"], ambient, campaign.points, material)
    hs, _ = worker.driver.shaft_coefficient(checked["T_shaft_mean_K"], ambient,
        worker.driver.angular_speed(case), campaign.points, material)
    mismatch = max(abs(hb/float(config["bush_outer_heat_transfer"])-1),
                   abs(hs/float(config["shaft_ambient_heat_transfer"])-1))
    if (mismatch >= campaign.h_tolerance or checked["h_relative_mismatch"] >= campaign.h_tolerance or
            checked["T_change_K"] >= campaign.wall_tolerance):
        raise ValueError("saved cooling closure fails independent correlations")
    fields = worker.driver.read_fields(run)
    with np.load(destination / "fields.npz", allow_pickle=False) as snapshot:
        if set(snapshot.files) != set(fields):
            raise ValueError("field snapshot members changed")
        for name, values in fields.items():
            if not np.array_equal(snapshot[name], values):
                raise ValueError(f"field snapshot differs from raw VTK: {name}")
    observables = worker.baseline.solid_observables(config, fields, .026)
    if observables != saved["observables"]:
        raise ValueError("saved TC observables differ from raw solids")
    return saved, checked


def summarize(root, output, old_report=None, density_report=None):
    """Emit paired results only when both native branches independently pass."""
    output.mkdir(parents=True, exist_ok=False)
    experiment = json.loads((root / "inputs/experiment.json").read_text())
    rows, errors = [], []
    for case in experiment["cases"]:
        key, measured = case["case_id"], case["measured"]
        row = {**measured, "case_id": key, "S_density_reference_TC2": case["S_density_reference_TC2"],
               "paired_status": "unavailable", "verification_error": ""}
        for name in ("gas_off_torque_nm", "gas_on_torque_nm", "gas_off_f", "gas_on_f",
                     "gas_on_minus_off_torque_nm", "gas_on_minus_off_f", "gas_on_minus_off_torque_percent",
                     "gas_on_alpha_max", "gas_on_dissolved_min", "gas_off_max_temperature_K", "gas_on_max_temperature_K"):
            row[name] = None
        for branch in worker.BRANCHES:
            for name in ("torque_error_nm", "f_error", "torque_error_percent", "f_error_percent"):
                row[f"{branch}_{name}"] = None
            for index in range(1, 4):
                row[f"{branch}_tc{index}_cfd_c"] = None
                row[f"{branch}_tc{index}_error_k"] = None
        unavailable = dict(row)
        try:
            off_saved, off = accepted(root, key, "gas_off")
            on_saved, on = accepted(root, key, "gas_on")
            sign = math.copysign(1., float(off_saved["config"]["omega"]))
            off_torque, on_torque = sign*off["torque_Nm"], sign*on["torque_Nm"]
            if min(off_torque, on_torque) <= 0:
                raise ValueError("invalid projected torque")
            factor = float(measured["force_n_mean"])*.0215
            row.update(paired_status="accepted_conditional_matched_native_comparison",
                gas_off_torque_nm=off_torque, gas_on_torque_nm=on_torque,
                gas_off_f=off_torque/factor, gas_on_f=on_torque/factor,
                gas_on_minus_off_torque_nm=on_torque-off_torque,
                gas_on_minus_off_f=(on_torque-off_torque)/factor,
                gas_on_minus_off_torque_percent=100*(on_torque/off_torque-1),
                gas_on_alpha_max=on["alpha_gas_max"], gas_on_dissolved_min=on["dissolved_gas_min"],
                gas_off_max_temperature_K=off["cross_film_profile_t_max_K"],
                gas_on_max_temperature_K=on["cross_film_profile_t_max_K"])
            for branch, saved, torque in (("gas_off", off_saved, off_torque), ("gas_on", on_saved, on_torque)):
                measured_torque, measured_f = float(measured["bearing_residual_nm"]), float(measured["f"])
                row[f"{branch}_torque_error_nm"] = torque-measured_torque
                row[f"{branch}_f_error"] = torque/factor-measured_f
                row[f"{branch}_torque_error_percent"] = 100*(torque/measured_torque-1) if measured_torque else None
                row[f"{branch}_f_error_percent"] = 100*(torque/factor/measured_f-1) if measured_f else None
                for index, temperature in enumerate(saved["observables"]["TC_C"], 1):
                    row[f"{branch}_tc{index}_cfd_c"] = temperature
                    row[f"{branch}_tc{index}_error_k"] = temperature-float(measured[f"tc{index}_c_mean"])
        except (OSError, ValueError, KeyError, IndexError, worker.ET.ParseError) as error:
            row = unavailable
            row["verification_error"] = str(error)
            errors.append({"case_id": key, "error": str(error)})
        rows.append(row)
    worker.write_rows(output / "matched_comparison.csv", rows)
    focused = []
    check_cases = experiment.get("numerical_check_cases", {"refinement": ["M02", "M16"], "refinement_layers32": ["M02"]})
    for phase, keys, reference in (("refinement", check_cases["refinement"], "production"),
                                    ("refinement_layers32", check_cases["refinement_layers32"], "refinement")):
        for key in keys:
            for branch in worker.BRANCHES:
                item = dict(case_id=key, branch=branch, phase=phase, reference=reference,
                    status="unavailable", error="", torque_difference_nm=None, torque_difference_percent=None,
                    tc1_difference_k=None, tc2_difference_k=None, tc3_difference_k=None)
                try:
                    saved, checked = accepted(root, key, branch, phase)
                    base, original = accepted(root, key, branch, reference)
                    sign = math.copysign(1., float(saved["config"]["omega"]))
                    delta = sign*(checked["torque_Nm"]-original["torque_Nm"])
                    item.update(status="verified", torque_difference_nm=delta,
                        torque_difference_percent=100*(checked["torque_Nm"]/original["torque_Nm"]-1))
                    for index, (value, old) in enumerate(zip(saved["observables"]["TC_C"], base["observables"]["TC_C"]), 1):
                        item[f"tc{index}_difference_k"] = value-old
                except (OSError, ValueError, KeyError, IndexError, worker.ET.ParseError) as error:
                    item["error"] = str(error)
                focused.append(item)
    worker.write_rows(output / "focused_checks.csv", focused)
    upstream = {}
    for label, path in (("old_fixed30", old_report), ("density_only", density_report)):
        upstream[label] = ({r["selected_id"]: r for r in worker.load_rows(path)} if path else {})
    historical = []
    for case, paired in zip(experiment["cases"], rows):
        key = case["case_id"]
        item = dict(case_id=key, status="upstream_comparison_unavailable", error="",
            old_fixed30_torque_nm=None, density_only_torque_nm=None,
            density_minus_old_torque_nm=None, density_minus_old_torque_percent=None,
            matched_off_minus_density_torque_nm=None)
        for index in range(1, 4):
            item[f"density_minus_old_tc{index}_k"] = None
        try:
            old, density = upstream["old_fixed30"][key], upstream["density_only"][key]
            if any(r["status"] != "accepted_conditional_surrogate" for r in (old, density)):
                raise ValueError("upstream report row is not accepted")
            for row in (old, density):
                for name, value in case["measured"].items():
                    if name in row and row[name] != value:
                        raise ValueError(f"upstream measured coordinate differs: {name}")
            ot, dt = float(old["bearing_torque_cfd_nm"]), float(density["bearing_torque_cfd_nm"])
            if not all(math.isfinite(v) and v > 0 for v in (ot, dt)):
                raise ValueError("invalid upstream torque")
            item.update(status="upstream_reported_accepted", old_fixed30_torque_nm=ot, density_only_torque_nm=dt,
                density_minus_old_torque_nm=dt-ot, density_minus_old_torque_percent=100*(dt/ot-1),
                matched_off_minus_density_torque_nm=(paired["gas_off_torque_nm"]-dt
                    if paired["paired_status"].startswith("accepted") else None))
            for index in range(1, 4):
                item[f"density_minus_old_tc{index}_k"] = float(density[f"tc{index}_cfd_c"])-float(old[f"tc{index}_cfd_c"])
        except (ValueError, KeyError) as error:
            item["error"] = str(error)
        historical.append(item)
    worker.write_rows(output / "old_vs_density_only.csv", historical)
    for filename in ("pure_coordinates.csv", "existing_cfd_coordinates.csv"):
        (output / filename).write_bytes((root / "inputs" / filename).read_bytes())
    worker.driver.save_json(output / "report.json", {
        "status": ("partial_unaccepted_pairs_preserved" if errors else "complete"
            if all(item["status"] == "verified" for item in focused) else "production_complete_focused_checks_unavailable"),
        "accepted_pairs": len(rows)-len(errors), "required_pairs": len(rows), "errors": errors,
        "focused_checks_verified": sum(item["status"] == "verified" for item in focused),
        "focused_checks_required": len(focused),
        "upstream_reports": {label: {"path": str(path), "sha256": worker.driver.sha256(path)}
            for label, path in (("old_fixed30", old_report), ("density_only", density_report)) if path},
        "upstream_basis": "Historical comparisons use explicitly supplied accepted upstream report rows and preserve their CSV hashes; they do not revalidate historical solver outputs here.",
        "attribution": "Gas-on minus matched COUPLED gas-off only; comparison to the separate FBNS STRIATED fixed30 baseline cannot isolate gas release.",
        "coordinates": "Measured rows and saved S unchanged; rho(TC2) reference coordinate added;103pure-oil+6CFD coordinates copied unchanged.",
        "limits": "Conditional chart bases, instantaneous phase equilibrium, pure-oil density and declared caloric proxies; numerical acceptance is not experimental validation."})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--old-report", type=Path, help="verified old campaign case_results.csv")
    parser.add_argument("--density-report", type=Path, help="verified density-only case_results.csv")
    args = parser.parse_args()
    summarize(args.root, args.out, args.old_report, args.density_report)
