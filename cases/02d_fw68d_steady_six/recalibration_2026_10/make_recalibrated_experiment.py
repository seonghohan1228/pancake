"""Build experiment_recalibrated_20261023.json from experiment.json.

Applies two measurement corrections and keeps every original value:

1. Radial load: the Oct 20 cylinder recalibration replaces the provisional
   F = 2752.246 (V - 0.0130312) N with F = K V, K = 2068.137 N/V (bounds
   1844.74-2291.53 N/V), V the CTT output in volts. Each saved force is first
   converted back to its CTT voltage with the provisional law.
2. Parasitic torque: the Oct 21 sweep with the test bearing removed (dry,
   unloaded, atmospheric). The up branch is the primary reference; the down
   branch and the Sep 11 references are kept for comparison.

Only load_mean_n changes the CFD input. The torque fields are comparison
targets. Run from this folder: python make_recalibrated_experiment.py
"""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE_DIR = HERE.parent
SOURCE = CASE_DIR / "experiment.json"
TARGET = CASE_DIR / "experiment_recalibrated_20261023.json"

K_OLD, V0_OLD = 2752.246065784855, 0.013031204699778745
calibration = json.loads((HERE / "load_cell_calibration_20261020.json").read_text())
K_MID, K_MIN, K_MAX = calibration["K_mid"], calibration["K_min"], calibration["K_max"]

reference = {}
with open(HERE / "parasitic_torque_20261021_by_rpm.csv", newline="") as f:
    for row in csv.DictReader(f):
        reference[round(float(row["rpm"]))] = row

LOAD_FIELDS = ("load_mean_n", "load_raw_min_n", "load_raw_max_n",
               "load_1min_mean_min_n", "load_1min_mean_max_n")


def recalibrate(force_old, k=K_MID):
    """Saved force with the provisional law -> force with the Oct 20 law [N]."""
    return k * (force_old / K_OLD + V0_OLD)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    for case in data["cases"]:
        old_load = case["load_mean_n"]
        for key in LOAD_FIELDS:
            case[key + "_provisional_2752"] = case[key]
            case[key] = recalibrate(case[key])
        case["ctt_v_mean_derived"] = old_load / K_OLD + V0_OLD
        case["load_mean_n_bounds_k_min_max"] = [recalibrate(old_load, K_MIN), recalibrate(old_load, K_MAX)]
        case["projected_pressure_pa_provisional_2752"] = case["projected_pressure_pa"]
        case["projected_pressure_pa"] = case["load_mean_n"] / (case["journal_diameter_m"] * case["bearing_length_m"])
        case["sommerfeld_number_provisional_2752"] = case["sommerfeld_number"]
        n_s = case["rpm_nominal"] / 60.0
        ratio = (case["journal_diameter_m"] / 2.0 / case["radial_clearance_m"]) ** 2
        case["sommerfeld_number"] = ratio * case["mu_pa_s"] * n_s / case["projected_pressure_pa"]

        row = reference[round(case["rpm_nominal"])]
        motor = case["motor_torque_mean_nm"]
        case["parasitic_torque_sep27_common_nm"] = case.pop("parasitic_torque_nm")
        case["apparent_bearing_torque_sep27_common_nm"] = case.pop("apparent_bearing_torque_nm")
        case["friction_coefficient_sep27_common_provisional_2752"] = case.pop("friction_coefficient")
        case["parasitic_torque_sep11_rpm_nm"] = float(row["old_sep11_nm"])
        case["parasitic_torque_20261021_up_nm"] = float(row["new_up_nm"])
        case["parasitic_torque_20261021_down_nm"] = float(row["new_down_nm"])
        radius = case["journal_diameter_m"] / 2.0
        for tag in ("sep11_rpm", "20261021_up", "20261021_down"):
            torque = motor - case[f"parasitic_torque_{tag}_nm"]
            case[f"apparent_bearing_torque_{tag}_nm"] = torque
            case[f"friction_coefficient_{tag}"] = torque / (case["load_mean_n"] * radius)

    data["description"] = ("Six A-F cases with the Oct 20 load-cell recalibration applied to all load "
                           "fields and the Oct 21 parasitic-torque sweep added as comparison targets. "
                           "Original values are kept with suffixes _provisional_2752 and _sep27_common.")
    data["recalibration_2026_10"] = {
        "applied": "2026-10-23 at Seongho's request for the CFD rerun",
        "load_law_new": f"F = {K_MID} * V_CTT  [N, V]",
        "load_law_new_bounds_n_per_v": [K_MIN, K_MAX],
        "load_law_provisional": f"F = {K_OLD} * (V_CTT - {V0_OLD})",
        "load_law_source": "recalibration_2026_10/load_cell_calibration_20261020.json; RECALIBRATION_2026-10.md",
        "parasitic_reference_primary": "parasitic_torque_20261021_up_nm",
        "parasitic_source": "recalibration_2026_10/parasitic_torque_20261021_by_rpm.csv",
        "source_experiment_sha256": sha256(SOURCE),
        "calibration_sha256": sha256(HERE / "load_cell_calibration_20261020.json"),
        "parasitic_sha256": sha256(HERE / "parasitic_torque_20261021_by_rpm.csv"),
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    for case in data["cases"]:
        print(f"{case['case_id']}: load {case['load_mean_n_provisional_2752']:.1f} -> {case['load_mean_n']:.1f} N"
              f" (bounds {case['load_mean_n_bounds_k_min_max'][0]:.1f}-{case['load_mean_n_bounds_k_min_max'][1]:.1f}),"
              f" S {case['sommerfeld_number']:.5f}, M_b(up) {case['apparent_bearing_torque_20261021_up_nm']:.4f} N m")


if __name__ == "__main__":
    main()
