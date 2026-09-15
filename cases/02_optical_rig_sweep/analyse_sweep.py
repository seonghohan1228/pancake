"""Reduce case-02 sweep runs to sweep_results.csv and answer Q1-Q5.

Angle convention (case-02 house convention, see README):
    phi = solver theta in degrees, 0 <= phi < 360.
    phi = 0 is the line of centres (minimum film) because attitude_angle_deg = 0
    and h = c - e cos(theta); phi increases in the direction of journal rotation.
    The cavity is therefore expected in 0 < phi < 180 (diverging film) and the
    pressure peak in 180 < phi < 360 (converging film).
"""

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S

FULL_FILM = 0.999            # film_content below this counts as cavitated
VAPOROUS_FLOOR = 500.0       # [Pa abs] vaporous threshold used in the exclusion argument
BAR = 1.0e5


def plane(run_dir, field):
    """Final-frame field on the (z, theta) plane plus the axis vectors."""
    a = S._read_pvts(S.final_pvts(run_dir))
    phi = np.degrees(np.asarray(a["theta_rad"][0])) % 360.0
    z = np.asarray(a["z_m"][0])
    value = np.asarray(a[field][0])
    phis, zs = np.unique(phi), np.unique(z)
    out = np.full((len(zs), len(phis)), np.nan)
    out[np.searchsorted(zs, z), np.searchsorted(phis, phi)] = value
    if np.isnan(out).any():
        raise RuntimeError(f"incomplete {field} plane in {run_dir}")
    return phis, zs, out


def largest_wrapped_run(mask):
    """(start_index, length) of the longest circular True run; (None, 0) if empty."""
    n = len(mask)
    if not mask.any():
        return None, 0
    if mask.all():
        return 0, n
    starts = [i for i in range(n) if mask[i] and not mask[i - 1]]
    best, best_len = None, 0
    for start in starts:
        length = 0
        while mask[(start + length) % n]:
            length += 1
        if length > best_len:
            best, best_len = start, length
    return best, best_len


def cavity_geometry(run_dir):
    """Circumferential arc, its end angles, axial reach and peak void fraction."""
    phis, zs, film = plane(run_dir, "film_content")
    d_phi = 360.0 / len(phis)
    cavitated = film < FULL_FILM
    columns = cavitated.any(axis=0)
    start, length = largest_wrapped_run(columns)
    axial = cavitated.sum(axis=0).max() / len(zs) * 100.0 if columns.any() else 0.0
    geometry = {
        "cav_arc_deg": round(columns.sum() * d_phi, 3),
        "cav_arc_contiguous_deg": round(length * d_phi, 3),
        "cav_axial_coverage_pct": round(float(axial), 3),
        "void_max": round(float((1.0 - film).max()), 6),
    }
    if start is None:
        geometry["phi_start_deg"] = ""
        geometry["phi_end_deg"] = ""
    else:
        geometry["phi_start_deg"] = round(float(phis[start]), 3)
        geometry["phi_end_deg"] = round(float(phis[(start + length - 1) % len(phis)]), 3)
    return geometry


def row_of(name, entry):
    m = entry.get("metrics") or {}
    omega = entry["omega"]
    torque = m.get("friction_torque", float("nan"))
    row = {
        "case": name,
        "run": entry["run"],
        "eps": entry["eps"],
        "rpm": entry["rpm"],
        "p_sat_setting": entry["p_cav_key"],
        "p_cav_bar_abs": round(entry["p_cav"] / BAR, 4),
        "c_um": round(entry["c"] * 1e6, 1),
        "L_mm": round(entry["L"] * 1e3, 2),
        "p_end_bar_abs": round(entry["p_end"] / BAR, 3),
        "mu_Pa_s": entry["mu"],
        "W_N": round(m.get("W", float("nan")), 3),
        "Fx_N": round(m.get("Fx", float("nan")), 3),
        "Fy_N": round(m.get("Fy", float("nan")), 3),
        "attitude_deg": round((m.get("force_angle_deg", float("nan")) + 180.0) % 360.0, 3),
        "p_min_bar_abs": round(m.get("p_min", float("nan")) / BAR, 5),
        "p_max_bar_abs": round(m.get("p_max", float("nan")) / BAR, 5),
        "friction_torque_Nm": round(torque, 5),
        "dissipation_W": round(torque * omega, 3),
        "h_min_um": round(m.get("h_min", float("nan")) * 1e6, 3),
        "T_max_C": round(m.get("T_max", float("nan")) - 273.15, 3),
        "mu_min_Pa_s": m.get("mu_min", ""),
        "wall_s": entry["wall_s"],
    }
    row.update(cavity_geometry(CASE_DIR / entry["dir"]))
    return row


COLUMNS = ["case", "run", "eps", "rpm", "p_sat_setting", "p_cav_bar_abs", "c_um", "L_mm",
           "p_end_bar_abs", "mu_Pa_s", "W_N", "Fx_N", "Fy_N", "attitude_deg",
           "p_min_bar_abs", "p_max_bar_abs", "cav_arc_deg", "cav_arc_contiguous_deg",
           "phi_start_deg", "phi_end_deg", "cav_axial_coverage_pct", "void_max",
           "friction_torque_Nm", "dissipation_W", "h_min_um", "T_max_C",
           "mu_min_Pa_s", "wall_s"]


def answers(rows):
    """Q1-Q5 from the RUN1 subset."""
    r1 = [r for r in rows if r["run"] == "RUN1"]
    out = {}

    # Q1: lowest rpm with a cavity, per threshold setting.
    q1 = {}
    for setting in sorted({r["p_sat_setting"] for r in r1}):
        per_eps = {}
        for eps in sorted({r["eps"] for r in r1}):
            hits = [r["rpm"] for r in r1
                    if r["p_sat_setting"] == setting and r["eps"] == eps and r["cav_arc_deg"] > 0]
            per_eps[eps] = min(hits) if hits else None
        q1[setting] = per_eps
    out["Q1_onset_rpm"] = q1

    # Q2: P1 candidates -- degassed clean at the same (eps, rpm), arc >= 90 deg at 0.8 bar.
    clean = {(r["eps"], r["rpm"]) for r in r1
             if r["p_sat_setting"] == "degassed" and r["cav_arc_deg"] == 0}
    candidates = [{"eps": r["eps"], "rpm": r["rpm"], "arc_deg": r["cav_arc_deg"],
                   "phi_start_deg": r["phi_start_deg"], "phi_end_deg": r["phi_end_deg"],
                   "W_N": r["W_N"], "p_max_bar_abs": r["p_max_bar_abs"]}
                  for r in r1
                  if r["p_sat_setting"] == "p080" and r["cav_arc_deg"] >= 90.0
                  and (r["eps"], r["rpm"]) in clean]
    candidates.sort(key=lambda c: -c["arc_deg"])
    old = [r for r in r1 if r["eps"] == 0.5 and r["rpm"] == 1000]
    out["Q2_old_P1_still_valid"] = bool(
        (0.5, 1000) in clean
        and any(r["p_sat_setting"] == "p080" and r["cav_arc_deg"] >= 90.0 for r in old))
    out["Q2_P1_candidates"] = candidates

    # Q3: degassed margin over the vaporous floor.
    margins = [{"eps": r["eps"], "rpm": r["rpm"],
                "p_min_bar_abs": r["p_min_bar_abs"],
                "margin_x": round(r["p_min_bar_abs"] * BAR / VAPOROUS_FLOOR, 2)}
               for r in r1 if r["p_sat_setting"] == "degassed"]
    out["Q3_degassed_margin"] = sorted(margins, key=lambda m: m["margin_x"])
    out["Q3_points_below_20x"] = [m for m in margins if m["margin_x"] < 20.0]

    # Q4: titration sensitivity, arc per bar of threshold, 0.2 -> 0.8 bar.
    sensitivity = []
    for eps in sorted({r["eps"] for r in r1}):
        for rpm in sorted({r["rpm"] for r in r1}):
            def arc(setting):
                hit = [r["cav_arc_deg"] for r in r1
                       if r["eps"] == eps and r["rpm"] == rpm and r["p_sat_setting"] == setting]
                return hit[0] if hit else None
            low, high = arc("p020"), arc("p080")
            if low is None or high is None:
                continue
            sensitivity.append({"eps": eps, "rpm": rpm, "arc_020_deg": low, "arc_080_deg": high,
                                "d_arc_per_bar": round((high - low) / 0.6, 2)})
    out["Q4_titration_sensitivity"] = sensitivity

    # Q5: peak film load over the whole RUN1 matrix.
    peak = max(r1, key=lambda r: r["W_N"])
    out["Q5_W_max"] = {"W_N": peak["W_N"], "eps": peak["eps"], "rpm": peak["rpm"],
                       "p_sat_setting": peak["p_sat_setting"],
                       "deflection_pct_of_c_config1": round(0.00899 * peak["W_N"], 2),
                       "deflection_pct_of_c_config2": round(0.00683 * peak["W_N"], 2)}

    # Q6: the isopentane arm -- lowest eps and rpm that produce a cavity at
    # p_cav = 0.12 bar, which decides whether the vaporous arm exists at all.
    iso = [r for r in r1 if r["p_sat_setting"] == "iso"]
    iso_hits = [r for r in iso if r["cav_arc_deg"] > 0]
    out["Q6_isopentane"] = {
        "cavitating_cases": len(iso_hits),
        "total_cases": len(iso),
        "lowest_rpm": min((r["rpm"] for r in iso_hits), default=None),
        "lowest_eps": min((r["eps"] for r in iso_hits), default=None),
        "onset_by_eps": {eps: min((r["rpm"] for r in iso_hits if r["eps"] == eps),
                                  default=None)
                         for eps in sorted({r["eps"] for r in iso})},
        "vaporous_arm_exists": bool(iso_hits),
        "cases": sorted([{"eps": r["eps"], "rpm": r["rpm"], "arc_deg": r["cav_arc_deg"],
                          "p_min_bar_abs": r["p_min_bar_abs"]} for r in iso_hits],
                        key=lambda r: (r["eps"], r["rpm"])),
    }

    # Q7: transducer range and the rigid-bore check, over every completed run.
    p_max = max(r["p_max_bar_abs"] for r in rows)
    p_min = min(r["p_min_bar_abs"] for r in rows)
    out["Q7_transducer"] = {
        "p_min_bar_abs": p_min, "p_max_bar_abs": p_max,
        "covered_by_0_2_bar": p_max <= 2.0,
        "covered_by_0_7_bar": p_max <= 7.0,
        "covered_by_0_10_bar": p_max <= 10.0,
        "recommended_full_scale_bar": math.ceil(p_max * 1.25),
        "resolution_needed_pct_fs": round(p_min / math.ceil(p_max * 1.25) * 100.0, 4),
    }
    bore = [{"case": r["case"], "expansion_um": round(23.0 * (r["p_max_bar_abs"] - 1.0) / 10.0, 2),
             "pct_of_c": round(23.0 * (r["p_max_bar_abs"] - 1.0) / 10.0 / r["c_um"] * 100.0, 2)}
            for r in rows]
    out["rigid_bore"] = max(bore, key=lambda b: b["pct_of_c"])
    return out


def main():
    ledger = json.loads((CASE_DIR / "sweep_runs.json").read_text())
    rows = [row_of(name, entry) for name, entry in sorted(ledger.items())
            if entry.get("rc") == 0 and entry.get("metrics")]
    with (CASE_DIR / "sweep_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    (CASE_DIR / "sweep_analysis.json").write_text(json.dumps(answers(rows), indent=2))
    print(f"{len(rows)} runs -> sweep_results.csv, sweep_analysis.json")


if __name__ == "__main__":
    main()
