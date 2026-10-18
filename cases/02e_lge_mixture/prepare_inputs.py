"""Verify preserved M01-M17 evidence and write preparation metadata, never CFD inputs/runs.

Run with Python + NumPy and required --inputs/--out directory arguments;
an existing different output is refused so prior evidence cannot be overwritten.
"""

import argparse
import csv
from decimal import Decimal
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import zipfile

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FILES = {
    "prompt": "LGE-2026-10-18-mixture-pancake-run-prompt.txt",
    "points": "LGE-2026-10-18-mixture-stribeck-points.csv",
    "checks": "LGE-2026-10-18-mixture-stribeck-checks.json",
    "archive": "LGE-2026-10-18-mixture-stribeck-evidence.zip",
    "chart": "LGE-2026-10-13-fw68d-chart-source.json",
    "drag_reference": "DAQ-test-plots-torque_vs_rpm_stable_60s_summary.csv",
    "pure_points": "LGE-2026-10-02-stribeck-rpm-drag-points.csv",
    "cfd_cases": "LGE-2026-10-02-stribeck-rpm-drag-cases.csv",
}
BASELINE = ROOT / "cases/02d_fw68d_steady_six/steady_setup.json"


def sha256(data):
    """Hash bytes without changing the source."""
    return hashlib.sha256(data).hexdigest()


def csv_rows(data):
    """Keep authoritative values as strings; parse only explicit calculations."""
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def csv_bytes(rows):
    """Serialize rows in their existing column order, with explicit nulls."""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows({k: "null" if v is None else v for k, v in row.items()} for row in rows)
    return stream.getvalue().encode("utf-8")


def decode_axis(anchors, pixel, offset=0):
    """Decode between calibrated pixels in log(value + offset), without extrapolation."""
    anchors = sorted(anchors, key=lambda anchor: anchor["pixel"])
    if not math.isfinite(pixel):
        raise ValueError("Nonfinite chart pixel")
    for left, right in zip(anchors, anchors[1:]):
        if left["pixel"] <= pixel <= right["pixel"]:
            fraction = (pixel - left["pixel"]) / (right["pixel"] - left["pixel"])
            return math.exp(math.log(left["value"] + offset) + fraction *
                            math.log((right["value"] + offset) / (left["value"] + offset))) - offset
    raise ValueError("Outside calibrated chart axis")


def viscosity(temperature_c, nodes):
    """Interpolate ln(nu) linearly in Celsius; reject every out-of-table request."""
    temperatures = np.asarray(temperature_c, dtype=float)
    if not (np.isfinite(temperatures).all() and
            (temperatures >= nodes[0, 0]).all() and (temperatures <= nodes[-1, 0]).all()):
        raise ValueError("Temperature outside decoded 30% viscosity table")
    return np.exp(np.interp(temperatures, nodes[:, 0], np.log(nodes[:, 1])))


def validate_ids(points):
    """Reject omissions, duplicates, additions, or reordering of saved case IDs."""
    if [row["selected_id"] for row in points] != [f"M{i:02d}" for i in range(1, 18)]:
        raise ValueError("Expected exactly ordered unique M01-M17")


def expect_value_error(function):
    """Require rejection of an invalid preflight input."""
    try:
        function()
    except ValueError:
        return
    raise AssertionError("Invalid input was accepted")


def main():
    """Verify evidence before producing any prepared_unrun output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path, help="Preserved evidence directory")
    parser.add_argument("--out", required=True, type=Path, help="New preflight output directory")
    args = parser.parse_args()
    input_directory, output_directory = args.inputs.resolve(), args.out.resolve()
    inputs = {key: (input_directory / name).read_bytes() for key, name in FILES.items()}
    checks = json.loads(inputs["checks"])
    for key in ("chart", "drag_reference", "pure_points", "cfd_cases"):
        assert sha256(inputs[key]) == checks["input_sha256"][key], key
    points = csv_rows(inputs["points"])
    validate_ids(points)
    expect_value_error(lambda: validate_ids(points[:-1]))
    expect_value_error(lambda: validate_ids([points[0]] + points[:-1]))
    chart = json.loads(inputs["chart"])
    assert chart["calibration"]["reviewed"]
    assert chart["calibration"]["x"]["scale"] == "logK"
    assert chart["calibration"]["y"]["scale"] == "log"
    curves = [curve for curve in chart["curves"] if curve["kind"] == "concentration" and curve["value"] == 30]
    assert len(curves) == 1
    xanchors = chart["calibration"]["x"]["anchors"]
    yanchors = chart["calibration"]["y"]["anchors"]
    x_min, x_max = min(a["pixel"] for a in xanchors), max(a["pixel"] for a in xanchors)
    included = [point for point in curves[0]["points"] if x_min <= point["x"] <= x_max]
    excluded = [point for point in curves[0]["points"] if not x_min <= point["x"] <= x_max]
    nodes = np.array(sorted((decode_axis(xanchors, p["x"], 273.15),
                             decode_axis(yanchors, p["y"]), p["p_MPa"]) for p in included))
    assert np.isfinite(nodes).all() and (np.diff(nodes[:, 0]) > 0).all() and (nodes[:, 1] > 0).all()
    np.testing.assert_allclose(nodes, checks["chart"]["decoded_nodes_T_C_nu_mm2_s_pressure_MPa"], rtol=1e-13)
    np.testing.assert_allclose(viscosity([30, 40], nodes), [8.2435, 6.6992], rtol=1e-5)
    np.testing.assert_allclose(viscosity(nodes[:, 0], nodes), nodes[:, 1], rtol=1e-14)
    for value in (nodes[0, 0] - 1e-9, nodes[-1, 0] + 1e-9, float("nan"), float("inf")):
        expect_value_error(lambda value=value: viscosity(value, nodes))
    expect_value_error(lambda: decode_axis(xanchors, x_min - 1))
    expect_value_error(lambda: decode_axis(yanchors, float("nan")))
    np.testing.assert_allclose(decode_axis([{"pixel": 0, "value": 1}, {"pixel": 2, "value": 100}], 1), 10)
    archive = zipfile.ZipFile(io.BytesIO(inputs["archive"]))
    assert archive.testzip() is None
    archived_script = archive.read("plot_mixture.py")
    assert sha256(archived_script) == checks["script_sha256"]
    assert csv_rows(archive.read("selected-mixture-points.csv")) == points
    archive_checks = json.loads(archive.read("checks.json"))
    assert archive_checks == checks
    pure, cfd = csv_rows(inputs["pure_points"]), csv_rows(inputs["cfd_cases"])
    pure_coordinates = archive.read("unchanged-pureoil-plotted-points.csv")
    cfd_coordinates = archive.read("unchanged-cfd-plotted-points.csv")
    archived_pure, archived_cfd = csv_rows(pure_coordinates), csv_rows(cfd_coordinates)
    assert len(pure) == len(archived_pure) == 103 and len(cfd) == len(archived_cfd) == 6
    for originals, archived, columns in ((pure, archived_pure, ("rpm_setting", "S_TC2", "new_f")),
                                          (cfd, archived_cfd, ("case_id", "rpm", "S_TC2", "unchanged_cfd_f"))):
        for original, saved in zip(originals, archived):
            for column in columns:
                assert original[column] == saved[column], (column, original[column], saved[column])
    references = csv_rows(inputs["drag_reference"])
    assert len({row["target_rpm"] for row in references}) == len(references)
    drag = {int(float(row["target_rpm"])): float(row["mean_torque_magnitude_Nm"]) for row in references}
    baseline = json.loads(BASELINE.read_bytes())
    baseline_head_bytes = subprocess.check_output(
        ["git", "show", "HEAD:cases/02d_fw68d_steady_six/steady_setup.json"], cwd=ROOT)
    assert json.loads(baseline_head_bytes) == baseline, "Baseline metadata differs from source HEAD"
    sensors = baseline["comparison_sensor_locations"]
    geometry = baseline["common_solver_settings"]
    assert (geometry["R"], geometry["L"], geometry["c"]) == (.0215, .04, 42.5e-6)
    matrix, errors = [], []
    for row in points:
        rpm, load = int(row["rpm"]), float(row["force_n_mean"])
        measured_rpm, requested_rpm = float(row["speed_rpm_mean"]), float(row["signed_requested_rpm"])
        assert load > 0 and rpm > 0 and abs(requested_rpm) == rpm
        assert measured_rpm * requested_rpm > 0
        assert row["all_criteria_pass"] == "True"
        assert int(row["window_minutes"]) == int(row["minute_summaries"]) == 15
        for field in ("tc0_c", "tc2_c", "force_n", "pressure_bar_abs", "speed_rpm"):
            assert float(row[field + "_mean_criterion_fraction"]) <= 1
        assert float(row["mean_speed_command_error_fraction"]) <= .01
        nu = float(viscosity(float(row["tc2_c_mean"]), nodes))
        mu = 950 * nu * 1e-6
        residual = float(row["directional_torque_nm_mean"]) - drag[rpm]
        s_value = mu * (rpm / 60) * (.04 * .043) / load * (.0215 / 42.5e-6)**2
        friction = residual / (load * .0215)
        recomputed = {"nu_tc2_mm2_s": nu, "mu_tc2_pa_s": mu, "bearing_residual_nm": residual,
                      "rpm_drag_reference_nm": drag[rpm], "S": s_value, "f": friction}
        for key, value in recomputed.items():
            np.testing.assert_allclose(value, float(row[key]), rtol=1e-13, atol=1e-15)
        assert float(row["rho_kg_m3_assumed"]) == 950
        errors.append({"selected_id": row["selected_id"], **{key + "_difference": value-float(row[key]) for key, value in recomputed.items()}})
        pressure_pa = str(Decimal(row["pressure_bar_abs_mean"]) * Decimal(100000))
        prepared = {**row, "status": "prepared_unrun", "measured_pressure_pa_abs": pressure_pa,
                    "motor_rotation_sign": 1 if measured_rpm > 0 else -1,
                    "nominal_signed_motor_rpm": rpm if measured_rpm > 0 else -rpm,
                    "rotation_mapping_status": "signed motor RPM verified; CFD theta mapping unresolved",
                    "solver_omega_rad_s": None, "proposed_radius_m": .0215,
                    "proposed_diameter_m": .043, "proposed_length_m": .04,
                    "proposed_nominal_radial_clearance_m": 42.5e-6,
                    "proposed_housing_outer_radius_m": .0425,
                    "proposed_groove_center_deg": 90, "proposed_groove_arc_width_m": .005,
                    "proposed_groove_width_deg": 13.298312524795712,
                    "proposed_groove_axial_length_m": .04,
                    "proposed_applied_load_x_n": 0, "proposed_applied_load_y_n": row["force_n_mean"],
                    "proposed_tc0_supply_ambient_proxy_k": float(row["tc0_c_mean"]) + 273.15,
                    "tc0_proxy_status": "proposed feed/surrounding temperature assumption; walls solved",
                    "tc2_role": "saved plotting coordinate and housing comparison; not imposed temperature",
                    "pressure_sensor_location": None, "pressure_assignment_status": "unresolved",
                    "groove_pressure_pa_abs": None, "south_reservoir_pressure_pa_abs": None,
                    "north_reservoir_pressure_pa_abs": None, "cavity_pressure_pa_abs": None,
                    "mixture_cp_j_kg_k": None, "mixture_conductivity_w_m_k": None,
                    "mixture_thermal_expansion_k_inverse": None,
                    "surrounding_fluid_heat_transfer_model": None,
                    "concentration_basis": chart["concentrationBasis"],
                    "chart_pressure_basis": chart["pressureBasis"],
                    "fixed_composition_status": "conditional 30% chart scenario; in-run composition unverified",
                    "density_status": "conditional constant 950 kg/m3; not measured",
                    "proposed_n_theta": 720, "proposed_n_z": 192, "proposed_film_layers": 16,
                    "proposed_housing_radial_cells": 16, "proposed_housing_axial_cells": 16,
                    "geometry_sensor_reference": "cases/02d_fw68d_steady_six/steady_setup.json; inherited proposal"}
        for sensor in ("TC1", "TC2", "TC3"):
            for coordinate in ("radius_m", "z_m", "theta_deg"):
                prepared[f"proposed_{sensor.lower()}_{coordinate}"] = sensors[sensor][coordinate]
        matrix.append(prepared)
    ordered = sorted(points, key=lambda row: int(row["start_bin_id"]))
    assert all(int(a["end_bin_id"]) < int(b["start_bin_id"]) for a, b in zip(ordered, ordered[1:]))
    ordered_s = sorted(points, key=lambda row: float(row["S"]))
    pilots = [{"role": role, "selected_id": ordered_s[index]["selected_id"], "S": ordered_s[index]["S"],
               "rpm": ordered_s[index]["rpm"], "status": "proposed_unrun"}
              for role, index in (("low_S", 0), ("middle_S", len(ordered_s)//2), ("high_S", -1))]
    report = {"status": "prepared_unrun", "counts": {"mixture": 17, "pure": 103, "existing_cfd": 6},
              "checks_passed": ["input/archive hash checks", "independent chart axis decoding and nodes",
                                "30/40C references", "exact-node interpolation", "strict domain and nonfinite rejection",
                                "IDs/count/order rejection checks", "17 saved S/f and drag/viscosity recalculations",
                                "selected condition fractions and nonoverlap", "signed measured/requested rotation agreement",
                                "103+6 source coordinates equal archive; archived bytes preserved",
                                "exact decimal absolute-pressure conversion"],
              "not_recomputed": "Raw minute summaries are not supplied here; selection and temporal bars are retained evidence, not independently regenerated.",
              "nu_30c_mm2_s": float(viscosity(30, nodes)), "nu_40c_mm2_s": float(viscosity(40, nodes)),
              "valid_temperature_range_c": nodes[:, 0][[0, -1]].tolist(),
              "excluded_uncalibrated_curve_points": excluded, "recalculation_differences": errors,
              "rotation_note": "All motor signs checked from measured speed; CW/CCW and native omega mapping require the measurement-axis convention.",
              "method_reference": FILES["archive"] + "!plot_mixture.py::decode,nu,coords (inspected, never executed)",
              "units_note": "Original archive variable names D=.04,L=.043 were interchanged; their product is correct. This matrix uses physical D=.043,L=.040.",
              "unresolved": ["concentration basis and in-run composition", "CFD rotation mapping",
                             "pressure sensor location and groove/axial/cavity pressures",
                             "mixture cp/conductivity and exterior-fluid heat-transfer assumptions",
                             "supported inputs for any native mixture phase-release model"],
              "simulation": {"run": False, "executable_sha256": None, "mpi_ranks": None,
                             "return_code": None, "convergence": None}}
    solver_viscosity_rows = [{"c_mass_fraction": 0, "T_C": t, "P_MPa": "", "viscosity_mm2s": v}
                            for t, v, _ in nodes]
    solver_viscosity = csv_bytes(solver_viscosity_rows)
    solver_rows = csv_rows(solver_viscosity)
    assert all(row["c_mass_fraction"] == "0" and row["P_MPa"] == "" for row in solver_rows)
    np.testing.assert_array_equal([[float(row["T_C"]), float(row["viscosity_mm2s"])]
                                  for row in solver_rows], nodes[:, :2])
    report["solver_viscosity_note"] = (
        "solution_viscosity_surrogate.csv contains only decoded viscosity nodes. c_mass_fraction=0 "
        "means no tracked species in the proposed surrogate, not measured mixture composition. "
        "Blank P_MPa supplies no pressure law or boundary; chart pressure labels are not imposed.")
    report["checks_passed"].append("solver viscosity-only CSV row values match decoded nodes exactly")
    outputs = {"case_matrix.csv": csv_bytes(matrix), "pure_coordinates.csv": pure_coordinates,
               "existing_cfd_coordinates.csv": cfd_coordinates,
               "solution_viscosity_surrogate.csv": solver_viscosity,
               "viscosity_nodes.csv": csv_bytes([{"temperature_c": t, "nu_mm2_s": v,
                                                    "conditional_mu_pa_s": 950*v*1e-6,
                                                    "chart_label_pressure_mpa_not_boundary": p}
                                                   for t, v, p in nodes]),
               "pilot_plan.json": json.dumps({"status": "proposed_unrun", "selection": "min, median, max saved S; no torque selection", "pilots": pilots}, indent=2).encode(),
               "selfcheck_report.json": json.dumps(report, indent=2).encode()}
    manifest = {"status": "prepared_unrun", "input_sha256": {FILES[k]: sha256(v) for k, v in inputs.items()},
                "input_directory": str(input_directory), "source_root": str(ROOT),
                "source_sha256": {str(Path(__file__).resolve()): sha256(Path(__file__).read_bytes()),
                                  str(BASELINE): sha256(BASELINE.read_bytes())},
                "baseline_head_sha256": sha256(baseline_head_bytes),
                "baseline_matches_head_json": True,
                "archive_member_sha256": {name: sha256(archive.read(name)) for name in archive.namelist()},
                "output_sha256": {name: sha256(data) for name, data in outputs.items()},
                "runtime": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__, "platform": platform.platform()},
                "checkout_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "source_note": "Source checkout is inferred from this script, not the caller's working directory. Baseline JSON matches this checkout HEAD; HEAD is not a CFD executable identity.",
                "solver_viscosity_note": report["solver_viscosity_note"],
                "output_null_encoding": "CSV literal null; JSON null", "original_columns_preserved": list(points[0])}
    outputs["manifest.json"] = json.dumps(manifest, indent=2).encode()
    output_directory.mkdir(parents=True, exist_ok=True)
    for name, data in outputs.items():
        target = output_directory / name
        if target.exists() and target.read_bytes() != data:
            raise FileExistsError(f"Refusing to replace different prior preflight evidence: {target}")
    for name, data in outputs.items():
        (output_directory / name).write_bytes(data)
    assert csv_rows(outputs["case_matrix.csv"]) == [{key: "null" if value is None else str(value) for key, value in row.items()} for row in matrix]
    print(json.dumps({"status": "prepared_unrun", "output": str(output_directory), "cases": len(matrix), "pilots": pilots}))


if __name__ == "__main__":
    main()
