"""Recheck the preserved conditional mixture campaign before reporting or plotting it.

Requires --root, --execution-manifest, --plot-style-dir and a new --out directory.
Use --self-check alone for calculation tests. No solver is launched by this file.
"""

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import zipfile

import numpy as np


IDS = [f"M{i:02d}" for i in range(1, 18)]
RESULT_COLUMNS = ["bearing_torque_cfd_nm", "f_cfd", "bearing_error_nm", "bearing_abs_error_nm", "bearing_error_percent",
                  "total_motor_diagnostic_nm", "motor_error_nm", "motor_abs_error_nm", "motor_error_percent",
                  "signed_native_torque_nm", "omega_rad_s", "tc1_cfd_c", "tc2_cfd_c", "tc3_cfd_c",
                  "tc1_error_k", "tc2_error_k", "tc3_error_k", "force_error_rel", "h_min_m", "eps",
                  "psi_deg", "p_min_Pa", "p_max_Pa", "T_max_K", "cross_film_profile_t_max_K",
                  "film_content_min", "film_content_max", "cavitated_area_fraction",
                  "energy_balance_residual_rel", "liquid_residual_rel", "h_bush_W_m2_K",
                  "h_shaft_W_m2_K", "h_relative_mismatch", "T_change_K", "run_dir",
                  "pressure_boundary_pa_abs", "density_kg_m3", "cp_J_kg_K", "conductivity_W_m_K",
                  "beta_1_K", "p_cav_Pa_abs", "film_layers",
                  "refinement_status", "refinement_torque_difference_nm", "refinement_torque_difference_percent",
                  "refinement_tc1_difference_k", "refinement_tc2_difference_k", "refinement_tc3_difference_k",
                  "layers32_status", "layers32_vs_tight16_torque_difference_nm",
                  "layers32_vs_tight16_torque_difference_percent", "layers32_vs_tight16_tc1_difference_k",
                  "layers32_vs_tight16_tc2_difference_k", "layers32_vs_tight16_tc3_difference_k"]


def sha256(path):
    """Hash an existing evidence file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def load_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def compare_values(actual, expected):
    """Compare saved observables to independent field sampling, including their schema."""
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys()
        for key in expected:
            compare_values(actual[key], expected[key])
    elif isinstance(expected, (list, float, int)):
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    else:
        assert actual == expected


def residuals(row, signed_torque, omega):
    """Project native journal shear torque onto rotation, using saved S unchanged."""
    torque = float(signed_torque) * math.copysign(1., float(omega))
    load, radius = float(row["force_n_mean"]), .0215
    bearing, motor = float(row["bearing_residual_nm"]), float(row["directional_torque_nm_mean"])
    if not all(math.isfinite(v) for v in (torque, load, bearing, motor, omega)) or min(torque, load, bearing, motor) <= 0 or omega == 0:
        raise ValueError("Invalid torque projection or comparison denominator")
    total = torque + float(row["rpm_drag_reference_nm"])
    return {"bearing_torque_cfd_nm": torque, "f_cfd": torque/(load*radius),
            "bearing_error_nm": torque-bearing, "bearing_abs_error_nm": abs(torque-bearing),
            "bearing_error_percent": 100*(torque-bearing)/bearing,
            "total_motor_diagnostic_nm": total, "motor_error_nm": total-motor,
            "motor_abs_error_nm": abs(total-motor),
            "motor_error_percent": 100*(total-motor)/motor,
            "signed_native_torque_nm": signed_torque, "omega_rad_s": omega}


def empty_result(row, original_columns, status, error=""):
    """Export authoritative experimental columns; omit provisional preflight proposals."""
    return {**{key: row[key] for key in original_columns}, "status": status,
            "verification_error": error, **{key: None for key in RESULT_COLUMNS}}


def self_check():
    row = {"S": "0.12345678901234567", "force_n_mean": "1000", "bearing_residual_nm": ".2",
           "directional_torque_nm_mean": ".9", "rpm_drag_reference_nm": ".7"}
    positive, negative = residuals(row, .25, 10), residuals(row, -.25, -10)
    np.testing.assert_allclose(positive["bearing_error_percent"], 25)
    np.testing.assert_allclose(positive["motor_error_percent"], 100*.05/.9)
    np.testing.assert_allclose(positive["bearing_abs_error_nm"], abs(positive["bearing_error_nm"]))
    np.testing.assert_allclose(positive["motor_abs_error_nm"], abs(positive["motor_error_nm"]))
    below = residuals(row, .15, 10)
    assert below["bearing_error_nm"] < 0 and below["motor_error_nm"] < 0
    assert below["bearing_abs_error_nm"] == -below["bearing_error_nm"]
    assert below["motor_abs_error_nm"] == -below["motor_error_nm"]
    assert positive["bearing_torque_cfd_nm"] == negative["bearing_torque_cfd_nm"]
    assert positive["f_cfd"] == negative["f_cfd"] and row["S"] == "0.12345678901234567"
    preflight = {**row, "pressure_assignment_status": "unresolved", "solver_omega_rad_s": "null",
                 "mixture_cp_j_kg_k": "null", "proposed_radius_m": ".0215"}
    missing = empty_result(preflight, row, "missing_unrun")
    assert missing["S"] == row["S"] and missing["bearing_torque_cfd_nm"] is None
    assert not (set(preflight)-set(row)) & missing.keys()
    assert all(missing[key] is None for key in ("pressure_boundary_pa_abs", "density_kg_m3", "cp_J_kg_K",
                                               "conductivity_W_m_K", "beta_1_K", "p_cav_Pa_abs", "film_layers"))
    heterogeneous = [empty_result(preflight, row, "selfcheck_only") for _ in range(17)]
    heterogeneous[1].update(pressure_boundary_pa_abs=1234567., density_kg_m3=950., cp_J_kg_K=2010.,
                            conductivity_W_m_K=.13, beta_1_K=.7284/950., p_cav_Pa_abs=0., film_layers=16)
    heterogeneous[1].update(layers32_status="selfcheck_only",
                            layers32_vs_tight16_torque_difference_nm=.001,
                            layers32_vs_tight16_torque_difference_percent=.4)
    for index in range(1, 4):
        heterogeneous[1][f"layers32_vs_tight16_tc{index}_difference_k"] = .001
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(heterogeneous[0]))
    writer.writeheader()
    writer.writerows(heterogeneous)
    serialized = list(csv.DictReader(io.StringIO(stream.getvalue())))
    assert len(serialized) == 17 and all(saved["S"] == row["S"] for saved in serialized)
    assert all(not (set(preflight)-set(row)) & saved.keys() for saved in serialized)
    assert serialized[0]["cp_J_kg_K"] == "" and serialized[1]["cp_J_kg_K"] == "2010.0"
    assert serialized[0]["layers32_status"] == "" and serialized[1]["layers32_status"] == "selfcheck_only"
    try:
        residuals(row, -.25, 10)
    except ValueError:
        pass
    else:
        raise AssertionError("Opposite-sign torque accepted")


def verify_case(root, directory, row, execution, worker, tight=False, layers=16):
    """Re-read native outputs and require independent provenance, physics and cooling closure."""
    key = row["selected_id"]
    job_key = f"refinement_layers32_{key}" if layers == 32 else f"refinement_{key}" if tight else key
    job = execution["jobs"][job_key]
    assert job["state"] == "complete" and job["returncode"] == 0, "Execution job is not successful"
    accepted = load_json(directory / "accepted.json")
    assert accepted["status"] == "accepted_conditional_surrogate" and accepted["tight"] == tight
    assert accepted.get("cross_film_layers", 16) == layers
    assert accepted["case"] == row, "Saved case differs from authoritative matrix"
    assert accepted["worker_sha256"] == sha256(worker.__file__)
    assert accepted["fields_sha256"] == sha256(directory / "fields.npz")
    result = accepted["result"]
    run = Path(result["run_dir"]).resolve()
    assert run.is_relative_to((directory / "runs" / key).resolve()), "Run lies outside named case"
    config_path = run / "config.txt"
    config = worker.parse_config(config_path.read_text())
    assert config == accepted["config"]
    assert sha256(config_path) == result["config_sha256"]
    metrics = load_json(run / "metrics.json")
    for name, value in metrics.items():
        assert result.get(name) == value, f"Accepted record changed native metrics: {name}"
    checkpoint = run / "results/steady_state.dat"
    assert checkpoint.is_file() and sha256(checkpoint) == result["output_state_sha256"]
    if result.get("initial_state_file"):
        assert sha256(result["initial_state_file"]) == result["initial_state_sha256"]
    rankfile = None
    command = result["command"]
    if "--rankfile" in command:
        rankfile = Path(command[command.index("--rankfile")+1])
    assert rankfile and rankfile.is_file() and sha256(rankfile) == result["rankfile_sha256"]
    inputs = root / "inputs"
    expected_setup = load_json(inputs / f"{key}_setup.json")
    if tight:
        expected_setup["solver_tolerances"] = dict(steady_force_tolerance=2e-5, coupling_tolerance=1e-8, fbns_tol=1e-12)
        expected_setup["cooling_relative_tolerance"] = 1e-5
        expected_setup["cooling_temperature_tolerance_K"] = 1e-4
    expected_setup["material_settings"]["cross_film_layers"] = layers
    assert load_json(directory / "setup.json") == expected_setup, "Unexpected scenario settings"
    options = worker.options(inputs, directory, root / "source/build-linux/pancake", rankfile)
    campaign = worker.driver.Campaign(options)
    case = campaign.cases[key]
    assert case["measured"] == row
    assert campaign.matching_provenance(result), "Source, executable, input or MPI provenance mismatch"
    checked = campaign.recheck(run, result, case)
    assert checked is not None, "Native output recheck rejected result"
    closed = campaign.find_completed(case, result["seed_eps"], result["seed_psi_deg"], closed=True)
    assert closed is not None and Path(closed["run_dir"]).resolve() == run, "Cooling closure recheck rejected named result"
    coupling_tolerance = expected_setup.get("solver_tolerances", {}).get("coupling_tolerance", 1e-6)
    assert checked["force_error_rel"] <= expected_setup.get("solver_tolerances", {}).get("steady_force_tolerance", .002)
    assert max(float(checked["diagnostics"]["coupling_R_T"]), checked["cross_film_profile_residual"],
               checked["energy_balance_residual_rel"]) < coupling_tolerance
    fbns = re.findall(r"fbns iter=\d+ \|G\|=" + worker.driver.NUMBER + r" \|FB\|=" + worker.driver.NUMBER,
                      (run / "run.log").read_text(errors="replace"))
    assert fbns and max(map(float, fbns[-1])) < expected_setup.get("solver_tolerances", {}).get("fbns_tol", 1e-6)
    pressure, tc0 = float(row["pressure_bar_abs_mean"])*100000, float(row["tc0_c_mean"])+273.15
    numeric = {"R": .0215, "L": .04, "c": 42.5e-6, "bush_outer_radius": .0425,
               "n_theta_global": 720, "n_z_global": 192, "cross_film_layers": layers,
               "bush_radial_cells": 16, "bush_axial_cells": 16, "rho": 950,
               "external_load_x": 0, "external_load_y": float(row["force_n_mean"]),
               "bc_z_south_val": pressure, "bc_z_north_val": pressure,
               "liquid_specific_heat": expected_setup["material_settings"]["oil_cp_J_kg_K"],
               "thermal_conductivity": expected_setup["material_settings"]["oil_conductivity_W_m_K"],
               "density_temperature_intercept": 950,
               "density_temperature_slope": -950*expected_setup["material_settings"]["oil_thermal_expansion_1_K"],
               "p_cav": expected_setup["material_settings"]["cavitation_pressure_Pa"],
               "dissolved_gas_initial": 0, "dissolved_gas_henry_coeff": 0, "viscosity_pressure_coeff": 0}
    for name in ("bc_z_south_temperature", "bc_z_north_temperature", "shaft_ambient_temperature",
                 "bush_outer_temperature", "thermal_ambient_temperature"):
        numeric[name] = tc0
    for name, value in numeric.items():
        assert math.isclose(float(config[name]), value, rel_tol=1e-13, abs_tol=1e-14), name
    for name, value in {"density_model": "CONSTANT", "fluid_property_model": "SINGLE_PHASE",
                        "liquid_viscosity_model": "KINEMATIC_TABLE", "solution_table_clamp": "false",
                        "kinematic_viscosity_log_interpolation": "true", "motion_model": "STEADY_EQUILIBRIUM",
                        "groove_mixing_model": "FLUX_BALANCED"}.items():
        assert config[name] == value, name
    np.testing.assert_allclose(list(map(float, config["inlet_groove"].split())),
                               [90, 13.298312524795712, pressure, tc0], rtol=1e-13)
    assert math.isclose(float(config["omega"]), float(row["rpm"])*math.pi/30, rel_tol=1e-13)
    assert float(row["speed_rpm_mean"]) < 0 and float(row["signed_requested_rpm"]) < 0
    assert sha256(config["solution_table_file"]) == sha256(inputs / "solution_viscosity_surrogate.csv")
    raw = worker.driver.read_fields(run)
    with np.load(directory / "fields.npz", allow_pickle=False) as saved:
        assert set(saved.files) == set(raw)
        for name, array in raw.items():
            np.testing.assert_array_equal(saved[name], array, err_msg=name)
    observables = worker.solid_observables(config, raw, .026)
    compare_values(accepted["observables"], observables)
    assert observables["sensor_z_m"] == .02
    values = residuals(row, checked["torque_Nm"], float(config["omega"]))
    values.update(pressure_boundary_pa_abs=float(config["bc_z_south_val"]), density_kg_m3=float(config["rho"]),
                  cp_J_kg_K=float(config["liquid_specific_heat"]), conductivity_W_m_K=float(config["thermal_conductivity"]),
                  beta_1_K=float(expected_setup["material_settings"]["oil_thermal_expansion_1_K"]),
                  p_cav_Pa_abs=float(config["p_cav"]), film_layers=int(config["cross_film_layers"]))
    for index in range(1, 4):
        temperature = observables["TC_C"][index-1]
        values[f"tc{index}_cfd_c"] = temperature
        values[f"tc{index}_error_k"] = temperature-float(row[f"tc{index}_c_mean"])
    for name in ("force_error_rel", "h_min_m", "eps", "psi_deg", "p_min_Pa", "p_max_Pa", "T_max_K",
                 "cross_film_profile_t_max_K", "film_content_min", "cavitated_area_fraction",
                 "energy_balance_residual_rel", "h_bush_W_m2_K", "h_shaft_W_m2_K", "h_relative_mismatch", "T_change_K"):
        values[name] = checked[name]
    values.update(film_content_max=float(np.max(raw["film_content"])),
                  liquid_residual_rel=float(checked["diagnostics"]["liquid_residual_rel"]), run_dir=str(run))
    hashes = {str(path): sha256(path) for path in (directory / "accepted.json", directory / "fields.npz",
              directory / "setup.json", config_path, run / "metrics.json", run / "run.log", checkpoint)}
    return values, {"accepted": accepted, "independently_checked": checked, "observables": observables,
                    "artifact_sha256": hashes, "execution_job": job}


def plot_comparison(output, preflight, rows, style_directory):
    """Plot saved experimental coordinates and independently accepted CFD without fits."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(style_directory))
    from plot_style import apply_style, clean_log_ticks, save_figure, smart_legend
    pure = load_rows(preflight / "pure_coordinates.csv")
    historical = load_rows(preflight / "existing_cfd_coordinates.csv")
    assert len(pure) == 103 and len(historical) == 6 and len(rows) == 17
    apply_style()
    figure, axes = plt.subplots(figsize=(8.8, 6.8))
    for rpm, color in ((200, "#124D87"), (1500, "#A7CEE2")):
        for data, x_column, y_column, marker, label in (
                (pure, "S_TC2", "new_f", "o", f"Pure FW68D, {rpm} rpm"),
                (historical, "S_TC2", "unchanged_cfd_f", "D", f"Existing CFD, {rpm} rpm")):
            speed_column = "rpm_setting" if data is pure else "rpm"
            selected = [row for row in data if float(row[speed_column]) == rpm]
            x, y = [float(row[x_column]) for row in selected], [float(row[y_column]) for row in selected]
            line, = axes.plot(x, y, marker, ls="none", color=color, ms=4.7 if data is pure else 6.5, label=label)
            np.testing.assert_array_equal(line.get_xdata(), x)
            np.testing.assert_array_equal(line.get_ydata(), y)
    x, y = np.array([[float(row[key]) for row in rows] for key in ("S", "f")])
    xlow, xhigh, ylow, yhigh = [np.array([float(row[key]) for row in rows]) for key in
                              ("S_1min_min", "S_1min_max", "f_1min_min", "f_1min_max")]
    axes.errorbar(x, y, xerr=[x-np.minimum(x, xlow), np.maximum(x, xhigh)-x],
                  yerr=[y-np.minimum(y, ylow), np.maximum(y, yhigh)-y], fmt="^", ls="none",
                  color="#D55E00", ecolor="#D55E00", ms=6, elinewidth=.8, capsize=2,
                  label="Mixture experiment, chart 30%")
    accepted = [row for row in rows if row["status"] == "accepted_conditional_surrogate"]
    if accepted:
        line, = axes.plot([float(row["S"]) for row in accepted], [row["f_cfd"] for row in accepted],
                          "s", ls="none", color="#6A3D9A", ms=6, label="Mixture CFD, fixed chart 30%")
        np.testing.assert_array_equal(line.get_xdata(), [float(row["S"]) for row in accepted])
    axes.set(xscale="log", yscale="log", xlabel="$S$", ylabel="$f$")
    clean_log_ticks(axes)
    smart_legend(axes, loc="best", fontsize=11)
    assert not axes.get_title() and all(line.get_linestyle() == "None" for line in axes.lines)
    save_figure(figure, str(output / "mixture_stribeck_comparison"))
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--execution-manifest", type=Path)
    parser.add_argument("--plot-style-dir", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    self_check()
    if args.self_check:
        print("Mixture reporter calculation self-check passed")
        return
    if not all((args.root, args.execution_manifest, args.plot_style_dir, args.out)):
        parser.error("reporting requires root, execution-manifest, plot-style-dir, and new out directory")
    root, output = args.root.resolve(), args.out.resolve()
    preflight = root / "preflight_verified"
    execution = load_json(args.execution_manifest)
    manifest = load_json(preflight / "manifest.json")
    for name, expected in manifest["output_sha256"].items():
        assert sha256(preflight / name) == expected, f"Preflight hash mismatch: {name}"
    for name, expected in manifest["input_sha256"].items():
        assert sha256(root / "evidence" / name) == expected, f"Authoritative evidence hash mismatch: {name}"
    for name, expected in load_json(root / "inputs/manifest.json").items():
        assert sha256(root / "inputs" / name) == expected, f"Prepared input hash mismatch: {name}"
    matrix = load_rows(preflight / "case_matrix.csv")
    assert [row["selected_id"] for row in matrix] == IDS
    original = load_rows(root / "evidence/LGE-2026-10-18-mixture-stribeck-points.csv")
    assert len(original) == 17
    for row, saved in zip(matrix, original):
        assert all(row[key] == value for key, value in saved.items()), "Experimental columns changed"
    with zipfile.ZipFile(root / "evidence/LGE-2026-10-18-mixture-stribeck-evidence.zip") as archive:
        for output_name, member in (("pure_coordinates.csv", "unchanged-pureoil-plotted-points.csv"),
                                    ("existing_cfd_coordinates.csv", "unchanged-cfd-plotted-points.csv")):
            assert (preflight / output_name).read_bytes() == archive.read(member)
    sys.path.insert(0, str(root / "source/cases/02e_lge_mixture"))
    import run_mixture as worker
    rows, verified, rejected = [], {}, []
    for row in matrix:
        key, directory = row["selected_id"], root / "production" / row["selected_id"]
        job = execution.get("jobs", {}).get(key, {})
        failure = directory / "failure.json"
        status = ("failed_preserved" if failure.exists() or job.get("state") == "failed" or job.get("returncode", 0)
                  else "missing_accepted_artifact" if job.get("state") == "complete"
                  else "running_unfinished" if job.get("state") == "running" else "missing_unrun")
        result = empty_result(row, original[0], status, load_json(failure).get("error", "") if failure.exists() else "")
        result.update(execution_state=job.get("state"), execution_returncode=job.get("returncode"))
        if (directory / "accepted.json").exists():
            try:
                values, evidence = verify_case(root, directory, row, execution, worker)
                result.update(values, status="accepted_conditional_surrogate")
                verified[key] = evidence
            except Exception as error:
                result.update(status="rejected_verification", verification_error=f"{type(error).__name__}: {error}")
                rejected.append(key)
        refined = None
        refinement = root / "refinement" / key
        if refinement.exists():
            result["refinement_status"] = "failed_or_missing_preserved"
            if (refinement / "accepted.json").exists():
                try:
                    refined, evidence = verify_case(root, refinement, row, execution, worker, tight=True)
                    verified[f"refinement_{key}"] = evidence
                    result["refinement_status"] = "accepted_tight_tolerance_check"
                    if result["status"] == "accepted_conditional_surrogate":
                        difference = refined["bearing_torque_cfd_nm"]-result["bearing_torque_cfd_nm"]
                        result["refinement_torque_difference_nm"] = difference
                        result["refinement_torque_difference_percent"] = 100*difference/result["bearing_torque_cfd_nm"]
                        for index in range(1, 4):
                            result[f"refinement_tc{index}_difference_k"] = refined[f"tc{index}_cfd_c"]-result[f"tc{index}_cfd_c"]
                except Exception as error:
                    result["refinement_status"] = f"rejected_verification: {type(error).__name__}: {error}"
                    rejected.append(f"refinement_{key}")
        layers32 = root / "refinement_layers32" / key
        if layers32.exists():
            result["layers32_status"] = "failed_or_missing_preserved"
            if (layers32 / "accepted.json").exists():
                try:
                    fine, evidence = verify_case(root, layers32, row, execution, worker, tight=True, layers=32)
                    verified[f"refinement_layers32_{key}"] = evidence
                    result["layers32_status"] = "accepted_tight_32_layer_check"
                    if refined is not None:
                        difference = fine["bearing_torque_cfd_nm"]-refined["bearing_torque_cfd_nm"]
                        result["layers32_vs_tight16_torque_difference_nm"] = difference
                        result["layers32_vs_tight16_torque_difference_percent"] = 100*difference/refined["bearing_torque_cfd_nm"]
                        for index in range(1, 4):
                            result[f"layers32_vs_tight16_tc{index}_difference_k"] = fine[f"tc{index}_cfd_c"]-refined[f"tc{index}_cfd_c"]
                except Exception as error:
                    result["layers32_status"] = f"rejected_verification: {type(error).__name__}: {error}"
                    rejected.append(f"refinement_layers32_{key}")
        assert result["S"] == row["S"] and result["f"] == row["f"]
        rows.append(result)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "case_results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_comparison(output, preflight, rows, args.plot_style_dir)
    accepted = [row for row in rows if row["status"] == "accepted_conditional_surrogate"]
    methods = [f"Accepted conditional CFD: {len(accepted)}/17. All 17 authoritative experimental rows and their original columns are retained; missing, failed and verification-rejected cases have blank CFD values. Control columns describe independently verified accepted normal runs. Provisional preflight proposal columns are omitted from this comparison CSV and remain in the preserved preflight evidence.",
               "Scenario assumptions are frozen in inputs/decisions.json, including when no CFD cases pass verification. The fixed 30% chart curve assumes density rho = 950 kg/m³. Concentration basis and in-run composition remain unverified. No dissolved/free refrigerant or release is calculated. Nominal geometry has no thermal expansion.",
               "Recorded chamber pressure in bar absolute is converted using 1 bar = 100000 Pa and applied to the groove and both axial reservoirs. The selected 0 Pa oil vapor limit suppresses refrigerant release. It is neither a measured mixture bubble pressure nor R410A partial pressure.",
               "TC0 is the assumed feed and liquid-bath temperature proxy. Film, shaft and housing temperatures are solved. TC1, TC2 and TC3 are sampled at r = 0.026 m, z = 0.020 m, and theta = 240°, 270° and 300°, respectively. TC0 has no independent prediction error. The thermal values c_p = 2010 J/(kg K), k = 0.13 W/(m K), buoyancy beta = (0.7284/950) K⁻¹, and steel k = 15 W/(m K) are inherited analogues, not measured mixture properties.",
               "S is copied unchanged from each saved TC2 experimental coordinate. Friction is f_CFD = M/(W R), with R = 0.0215 m. Positive comparison torque is native cylindrical-journal shear torque projected onto configured rotation. It includes Couette and pressure-gradient shear, but not shaft-end, motor, seal or unresolved groove-recirculation drag. The inherited theta frame follows rotation. Negative logger RPM alone does not establish laboratory CW/CCW.",
               "Bearing error is CFD minus saved bearing_residual_nm; its percentage denominator is that saved bearing residual. The separate motor diagnostic adds the saved RPM parasitic reference to CFD torque; its percentage denominator is measured directional motor torque. bearing_error_nm and motor_error_nm are signed dimensional differences. bearing_abs_error_nm and motor_abs_error_nm are their absolute magnitudes. Transfer of the parasitic reference to pressurized mixture is unverified. Temperature errors are signed CFD minus measured values, reported in K.",
               "Every accepted result is rechecked against raw native metrics/log/config/VTK, saved NPZ/checkpoint, source/driver/worker/executable/input hashes, successful execution status and cooling closure. Normal and tight pressure solves use FBNS tolerances of 1e-10 and 1e-12, respectively. These tighter solves follow the preserved initial M02 pilot, whose relative mass residual was -4.553754639e-8 despite a final FBNS G residual of 6.023062e-9. The required mass closure remains 1e-8; acceptance and physics are unchanged. Capped/restored/unconverged trials are excluded. The six historical CFD markers are preserved references, not new verified mixture results.",
               "Plot bars show only the saved one-minute temporal envelopes. They are not confidence intervals or composition, calibration or numerical uncertainty. Focused comparisons use normal versus tight tolerances at 16 film layers, then tight tolerances at 16 versus 32 film layers. Available differences are reported separately and do not replace normal CFD points or establish mesh independence. Housing resolution remains 16 radial by 16 axial cells.",
               "Agreement alone does not validate the mixture model or prove steady composition/friction.",
               "| Case | Status | Bearing error (%) | TC1 error (K) | TC2 error (K) | TC3 error (K) |", "|---|---|---:|---:|---:|---:|"]
    methods.extend(f"| {row['selected_id']} | {row['status']} | {'' if row['bearing_error_percent'] is None else format(row['bearing_error_percent'], '.5g')} | " +
                   " | ".join("" if row[f"tc{index}_error_k"] is None else format(row[f"tc{index}_error_k"], '+.5g')
                              for index in range(1, 4)) + " |" for row in rows)
    (output / "methods_results.md").write_text("\n\n".join(methods[:9])+"\n\n"+"\n".join(methods[9:])+"\n", encoding="utf-8")
    report = {"status": "verified_conditional_comparison" if accepted else "no_accepted_cfd",
              "accepted_count": len(accepted), "rejected": rejected, "cases": rows, "verified_records": verified,
              "execution": execution, "scenario_decisions": load_json(root / "inputs/decisions.json"),
              "source_sha256": worker.driver.source_sha256(), "worker_sha256": sha256(worker.__file__),
              "reporter_sha256": sha256(__file__), "plot_style_sha256": sha256(args.plot_style_dir / "plot_style.py"),
              "analysis_source_sha256": {str(path): sha256(path) for path in
                  (Path(worker.driver.__file__), root / "source/cases/02d_fw68d_steady_six/thermal_expansion.py")},
              "evidence_sha256": {name: sha256(root / "evidence" / name) for name in manifest["input_sha256"]},
              "input_sha256": {str(path): sha256(path) for path in (args.execution_manifest,
                  preflight / "manifest.json", root / "inputs/manifest.json")},
              "output_sha256": {path.name: sha256(path) for path in sorted(output.iterdir())}}
    (output / "verified_comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"accepted": len(accepted), "rejected": rejected, "output": str(output)}))
    if rejected:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
