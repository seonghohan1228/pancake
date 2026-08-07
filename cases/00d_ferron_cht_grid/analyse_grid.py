"""Score the G5/G6/G7 Ferron CHT equilibria and apply pre-registered checks."""

import csv
import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
SOURCE_CASE = CASE_DIR.parent / "01_ferron_wallbc_fields"
REF_PATH = CASE_DIR.parent / "ref_data" / "ferron1983_fig4_2000rpm_4kN_temperature.csv"
PRESSURE_REF_PATH = CASE_DIR.parent / "ref_data" / "ferron1983_fig4_2000rpm_4kN_pressure.csv"
P_AMBIENT = 101325.0
sys.path.insert(0, str(SOURCE_CASE))

import solver_io as S


LIMITS = {
    "eps_rel_pct": 0.5,
    "psi_deg": 0.25,
    "p_max_rel_pct": 1.0,
    "mean_error_K": 0.10,
    "rms_error_K": 0.10,
    "peak_temperature_K": 0.25,
    "peak_location_deg": 2.0,
}


def read_reference(path=REF_PATH):
    points = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row or row[0].lstrip().startswith("#"):
                continue
            points.append((float(row[0]), float(row[1]) + 273.15))
    return points


def circular_difference(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def periodic_interpolate(theta, values, target):
    pairs = sorted((x % 360.0, y) for x, y in zip(theta, values))
    xs = [x for x, _ in pairs]
    ys = [y for _, y in pairs]
    xs = [xs[-1] - 360.0, *xs, xs[0] + 360.0]
    ys = [ys[-1], *ys, ys[0]]
    target %= 360.0
    for i in range(len(xs) - 1):
        if xs[i] <= target <= xs[i + 1]:
            weight = (target - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + weight * (ys[i + 1] - ys[i])
    raise RuntimeError("periodic interpolation bracket not found")


def score_profile(run_dir, reference):
    profile = S.midplane_profile(CASE_DIR / run_dir, "T_bearing")
    if profile is None:
        raise RuntimeError(f"missing T_bearing profile in {run_dir}")
    theta, temperature = profile
    errors = []
    for x_ref, t_ref in reference:
        t_model = periodic_interpolate(theta, temperature, 90.0 + x_ref)
        errors.append(t_model - t_ref)
    peak_index = max(range(len(temperature)), key=temperature.__getitem__)
    return {
        "mean_error_K": sum(errors) / len(errors),
        "rms_error_K": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "max_abs_error_K": max(abs(value) for value in errors),
        "peak_temperature_K": temperature[peak_index],
        "peak_location_deg": (theta[peak_index] - 90.0) % 360.0,
        "n_reference": len(errors),
    }


def score_pressure(run_dir, reference):
    """Peak midplane gauge pressure and errors against the digitised peak."""
    theta, pressure = S.midplane_profile(CASE_DIR / run_dir, "p")
    peak_index = max(range(len(pressure)), key=pressure.__getitem__)
    model_gauge = pressure[peak_index] - P_AMBIENT
    model_location = (theta[peak_index] - 90.0) % 360.0
    exp_location, exp_gauge = max(reference, key=lambda item: item[1])
    return {
        "p_peak_gauge_MPa": model_gauge / 1.0e6,
        "p_peak_location_deg": model_location,
        "p_peak_error_pct": 100.0 * (model_gauge - exp_gauge) / exp_gauge,
        "p_peak_location_error_deg": circular_difference(model_location, exp_location),
    }


def relative_change_percent(coarse, fine):
    return abs(fine - coarse) / max(abs(fine), 1.0e-300) * 100.0


def analyse():
    equilibria = json.loads((CASE_DIR / "grid_equilibrium.json").read_text(encoding="utf-8"))
    reference = read_reference()
    pressure_reference = read_reference(PRESSURE_REF_PATH)
    out = {"grids": {}, "g6_to_g7": {}, "limits": LIMITS, "pass": False}

    for tag in ("G5", "G6", "G7"):
        result = equilibria[tag]["flow"]
        metrics = result["metrics"]
        row = {
            "eps": result["eps"],
            "psi_deg": result["psi"],
            "run_dir": result["dir"],
            "W_N": metrics["W"],
            "force_angle_deg": metrics["force_angle_deg"],
            "p_max_Pa": metrics["p_max"],
            "q_supply_m3s": equilibria[tag]["q_supply"],
            "q_supply_Lmin": equilibria[tag]["q_supply"] * 6.0e4,
            "coupling": metrics["coupling"],
        }
        row.update(score_profile(result["dir"], reference))
        row.update(score_pressure(result["dir"], pressure_reference))
        out["grids"][tag] = row

    g6, g7 = out["grids"]["G6"], out["grids"]["G7"]
    checks = {
        "eps_rel_pct": relative_change_percent(g6["eps"], g7["eps"]),
        "psi_deg": circular_difference(g6["psi_deg"], g7["psi_deg"]),
        "p_max_rel_pct": relative_change_percent(g6["p_max_Pa"], g7["p_max_Pa"]),
        "mean_error_K": abs(g7["mean_error_K"] - g6["mean_error_K"]),
        "rms_error_K": abs(g7["rms_error_K"] - g6["rms_error_K"]),
        "peak_temperature_K": abs(g7["peak_temperature_K"] - g6["peak_temperature_K"]),
        "peak_location_deg": circular_difference(g7["peak_location_deg"],
                                                   g6["peak_location_deg"]),
    }
    out["g6_to_g7"] = {
        key: {"value": value, "limit": LIMITS[key], "pass": value <= LIMITS[key]}
        for key, value in checks.items()
    }
    out["pass"] = all(item["pass"] for item in out["g6_to_g7"].values())
    (CASE_DIR / "grid_analysis.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("grid eps psi_deg p_max_MPa mean_K rms_K peak_C peak_x_deg")
    for tag, row in out["grids"].items():
        print(f"{tag} {row['eps']:.8f} {row['psi_deg']:.6f} "
              f"{row['p_max_Pa']/1e6:.6f} {row['mean_error_K']:+.4f} "
              f"{row['rms_error_K']:.4f} {row['peak_temperature_K']-273.15:.4f} "
              f"{row['peak_location_deg']:.4f}")
    print("G6->G7 checks:")
    for key, item in out["g6_to_g7"].items():
        print(f"  {key}: {item['value']:.6g} <= {item['limit']:.6g} "
              f"{'PASS' if item['pass'] else 'FAIL'}")
    print(f"VERDICT: {'PASS' if out['pass'] else 'FAIL'}")


if __name__ == "__main__":
    analyse()
