"""Write a recalibrated copy of a verified 02e preflight directory.

Load: the 2026-10-20 load-cell recalibration replaces the provisional overlay
F_old = K_old (V - V0_old) with F = K V. Each recorded force maps affinely,
F_new = K (F_old / K_old + V0_old), so means, sample standard deviations, sums
and sums of squares transform exactly. K defaults to the mid value K_mid; the
K_min/K_max bounds are recorded as force_n_mean_lo/hi.

Torque: the 2026-10-21 bearing-removed parasitic sweep (up branch) replaces the
2026-09-11 reference. bearing_residual_nm becomes M_motor - M_p,up from the
supplied table; the old values are kept as *_sep11 columns.

S and f are rescaled to the new load and torque so 02e input checks hold. All
other columns, including the original window-selection criteria, are unchanged.
"""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

K_OLD = 2752.246065784855
V0_OLD = 0.013031204699778745
R_JOURNAL = 0.0215


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recalibrate_rows(rows, calibration, parasitic):
    """Return recalibrated rows; parasitic maps case id -> (bearing torque up, motor torque)."""
    k_mid, k_lo, k_hi = calibration["K_mid"], calibration["K_min"], calibration["K_max"]
    a, b = k_mid / K_OLD, k_mid * V0_OLD
    out = []
    for row in rows:
        key = row["selected_id"]
        if key not in parasitic:
            raise ValueError(f"no Oct 21 parasitic torque for {key}")
        torque_up, motor = parasitic[key]
        if abs(motor - float(row["directional_torque_nm_mean"])) > 1e-9 * max(1.0, motor):
            raise ValueError(f"motor torque mismatch for {key}")
        new = dict(row)
        n = float(row["force_n_count"])
        mean_old, total, squares = float(row["force_n_mean"]), float(row["force_n_sum"]), float(row["force_n_sq_sum"])
        mean_new = a * mean_old + b
        new.update(force_n_mean=repr(mean_new), force_n_std=repr(a * float(row["force_n_std"])),
                   force_n_sum=repr(a * total + b * n),
                   force_n_sq_sum=repr(a * a * squares + 2 * a * b * total + b * b * n),
                   force_n_mean_range_1min=repr(a * float(row["force_n_mean_range_1min"])))
        drag_old, torque_old = float(row["rpm_drag_reference_nm"]), float(row["bearing_residual_nm"])
        drag_new = motor - torque_up
        load_ratio = mean_old / mean_new
        new.update(rpm_drag_reference_nm=repr(drag_new), bearing_residual_nm=repr(torque_up),
                   S=repr(float(row["S"]) * load_ratio),
                   f=repr(torque_up / (mean_new * R_JOURNAL)),
                   f_1min_min=repr((float(row["f_1min_min"]) * mean_old * R_JOURNAL + drag_old - drag_new) / (mean_new * R_JOURNAL)),
                   f_1min_max=repr((float(row["f_1min_max"]) * mean_old * R_JOURNAL + drag_old - drag_new) / (mean_new * R_JOURNAL)),
                   force_n_mean_sep02_overlay=row["force_n_mean"], S_sep02_overlay=row["S"], f_sep11_parasitic=row["f"],
                   rpm_drag_reference_nm_sep11=row["rpm_drag_reference_nm"], bearing_residual_nm_sep11=row["bearing_residual_nm"],
                   force_n_mean_lo=repr(k_lo * (mean_old / K_OLD + V0_OLD)),
                   force_n_mean_hi=repr(k_hi * (mean_old / K_OLD + V0_OLD)),
                   load_recalibration_K_N_per_V=repr(k_mid))
        out.append(new)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True, help="LGE-2026-10-20-recal-calibration.json")
    parser.add_argument("--parasitic", type=Path, required=True, help="LGE-2026-10-21-mixture17_cfd_new_parasitic.csv")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    calibration = json.loads(args.calibration.read_text())
    parasitic = {r["selected_id"]: (float(r["Mexp_up"]), float(r["directional_torque_nm_mean"]))
                 for r in csv.DictReader(args.parasitic.open(encoding="utf-8"))}
    shutil.copytree(args.preflight, args.out)
    matrix = args.out / "case_matrix.csv"
    rows = list(csv.DictReader(matrix.open(encoding="utf-8")))
    new_rows = recalibrate_rows(rows, calibration, parasitic)
    fields = list(rows[0]) + [k for k in new_rows[0] if k not in rows[0]]
    with matrix.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(new_rows)
    record = dict(source_preflight=str(args.preflight), source_case_matrix_sha256=sha256(args.preflight / "case_matrix.csv"),
                  calibration=str(args.calibration), calibration_sha256=sha256(args.calibration),
                  parasitic=str(args.parasitic), parasitic_sha256=sha256(args.parasitic),
                  K_old=K_OLD, V0_old=V0_OLD, K_mid=calibration["K_mid"], K_min=calibration["K_min"], K_max=calibration["K_max"],
                  torque_reference="2026-10-21 bearing-removed parasitic sweep, up branch",
                  load_ratio_new_over_old=[float(n["force_n_mean"]) / float(o["force_n_mean"]) for o, n in zip(rows, new_rows)],
                  case_matrix_sha256=sha256(matrix))
    (args.out / "recalibration.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(dict(cases=len(new_rows), ratio_min=min(record["load_ratio_new_over_old"]),
                          ratio_max=max(record["load_ratio_new_over_old"]))))


if __name__ == "__main__":
    main()
