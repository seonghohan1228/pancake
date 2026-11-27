"""Check the 02e load/torque recalibration of a preflight case matrix."""

import math
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]
module = runpy.run_path(str(ROOT / "cases/02e_lge_mixture/recalibrate_preflight.py"))
recalibrate_rows, K_OLD, V0_OLD = module["recalibrate_rows"], module["K_OLD"], module["V0_OLD"]

samples = [1000.0, 1010.0, 990.0, 1005.0]
n = len(samples)
mean = sum(samples) / n
std = math.sqrt(sum((x - mean) ** 2 for x in samples) / (n - 1))
row = dict(selected_id="M01", rpm="200", force_n_count=str(n), force_n_sum=repr(sum(samples)),
           force_n_sq_sum=repr(sum(x * x for x in samples)), force_n_mean=repr(mean), force_n_std=repr(std),
           force_n_mean_range_1min="10.0", directional_torque_nm_mean="0.9", rpm_drag_reference_nm="0.6",
           bearing_residual_nm="0.3", S="0.01", f=repr(0.3 / (mean * 0.0215)),
           f_1min_min=repr(0.29 / (mean * 0.0215)), f_1min_max=repr(0.31 / (mean * 0.0215)))
calibration = dict(K_mid=2068.0, K_min=1845.0, K_max=2292.0)
new = recalibrate_rows([row], calibration, {"M01": (0.5, 0.9)})[0]

mapped = [2068.0 * (x / K_OLD + V0_OLD) for x in samples]
new_mean = sum(mapped) / n
assert math.isclose(float(new["force_n_mean"]), new_mean, rel_tol=1e-12)
assert math.isclose(float(new["force_n_sum"]), sum(mapped), rel_tol=1e-12)
assert math.isclose(float(new["force_n_sq_sum"]), sum(x * x for x in mapped), rel_tol=1e-12)
assert math.isclose(float(new["force_n_std"]),
                    math.sqrt(sum((x - new_mean) ** 2 for x in mapped) / (n - 1)), rel_tol=1e-12)
assert math.isclose(float(new["S"]), 0.01 * mean / new_mean, rel_tol=1e-14)
assert math.isclose(float(new["bearing_residual_nm"]), 0.5) and math.isclose(float(new["rpm_drag_reference_nm"]), 0.4)
assert math.isclose(float(new["f"]), 0.5 / (new_mean * 0.0215), rel_tol=1e-14)
assert math.isclose(float(new["f_1min_min"]), 0.49 / (new_mean * 0.0215), rel_tol=1e-12)
assert float(new["force_n_mean_lo"]) < new_mean < float(new["force_n_mean_hi"])
assert new["force_n_mean_sep02_overlay"] == row["force_n_mean"] and new["bearing_residual_nm_sep11"] == "0.3"
for bad in ({}, {"M01": (0.5, 0.8)}):
    try:
        recalibrate_rows([row], calibration, bad)
    except ValueError:
        pass
    else:
        raise AssertionError("missing or inconsistent parasitic torque must be rejected")
print("PASS: load affine map, sample statistics, S/f rescaling and parasitic replacement")
