"""Case 01-CONJUGATE: pull the reported numbers out of the archived equilibria.

Writes conjugate_analysis.json and prints the tables that go into the summary.

The quantity under test is the BUSH-SIDE FILM SURFACE temperature, written by the
solver as `T_bearing` (`temperature_bearing`), which is the closest model quantity
to Ferron's thermocouples 0.5 mm below the bush inner surface. Against it sit two
references:

  * the 16 digitised Fig. 4 points (cases/ref_data/...temperature.csv), compared
    point by point and by arc segment, with the stated mapping
    solver_theta = (90 + x) mod 360 and no fitted shift;
  * variant A `walls_adiabatic` of cases/01_ferron_thd_crossfilm, which is the
    SAME model with the LOCAL 1-D bush closure. That comparison is the point of
    the case: whether 2-D bush conduction moves the predicted peak off the groove.
"""

import csv
import json
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S

REF = CASE_DIR.parent / "ref_data"
LOCAL_CASE = CASE_DIR.parent / "01_ferron_thd_crossfilm"
GROOVE_DEG = 90.0
P_AMB = 101325.0

VARIANTS = ("outer_fixed313", "outer_robin80")

# Same convention as the LOCAL case: the 18 deg groove is centred on x = 0, and
# its two edge cells carry a one-cell spike in the wall temperatures because the
# groove pins the film temperature over a single cell. Peaks are quoted with that
# band removed so a numerical edge artefact is never reported as a thermal crest.
GROOVE_EXCLUDE_DEG = 12.0

# Arc segments used for the error breakdown, copied from the LOCAL case summary so
# the two are directly comparable.
SEGMENTS = ((0.0, 140.0), (140.0, 270.0), (270.0, 360.0))


def read_csv(path):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]))
    return xs, ys


def peak(xs, ys):
    i = max(range(len(ys)), key=lambda k: ys[k])
    return xs[i], ys[i]


def peak_away_from_groove(xs, ys, band=GROOVE_EXCLUDE_DEG):
    keep = [(x, y) for x, y in zip(xs, ys) if band <= x <= 360.0 - band]
    if not keep:
        return peak(xs, ys)
    i = max(range(len(keep)), key=lambda k: keep[k][1])
    return keep[i]


def profile(run_dir, field, base=CASE_DIR):
    got = S.midplane_profile(Path(base) / run_dir, field)
    if got is None:
        return None
    theta, v = got
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def interp(xs, ys, x):
    """Linear interpolation on the periodic theta axis."""
    if x <= xs[0]:
        return ys[0]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            f = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + f * (ys[i] - ys[i - 1])
    return ys[-1]


def errors_vs_experiment(xs, t_celsius, xr, yr):
    """Point-by-point model - experiment, plus mean/rms and per-segment means."""
    d = [interp(xs, t_celsius, x) - y for x, y in zip(xr, yr)]
    out = {
        "n_points": len(d),
        "mean_K": sum(d) / len(d),
        "rms_K": (sum(v * v for v in d) / len(d)) ** 0.5,
        "max_abs_K": max(abs(v) for v in d),
        "per_point": [{"x_deg": x, "exp_C": y, "model_C": interp(xs, t_celsius, x),
                       "d_K": v} for x, y, v in zip(xr, yr, d)],
        "segments": {},
    }
    for lo, hi in SEGMENTS:
        sel = [v for x, v in zip(xr, d) if lo <= x < hi]
        out["segments"][f"{lo:.0f}-{hi:.0f}"] = {
            "n": len(sel),
            "mean_K": sum(sel) / len(sel) if sel else None,
        }
    return out


def local_reference():
    """Variant A of the LOCAL-bush case: same model, 1-D radial bush closure."""
    eq_path = LOCAL_CASE / "equilibrium.json"
    an_path = LOCAL_CASE / "thd_analysis.json"
    if not eq_path.exists():
        return None
    eq = json.loads(eq_path.read_text())
    if "walls_adiabatic" not in eq:
        return None
    rec = eq["walls_adiabatic"]["flow"]
    prof = profile(rec["dir"], "cross_film_t_bush", base=LOCAL_CASE)
    entry = {"dir": str(LOCAL_CASE.name) + "/" + rec["dir"],
             "eps": rec["eps"], "psi_deg": rec["psi"]}
    if an_path.exists():
        a = json.loads(an_path.read_text())
        v = a.get("variants", {}).get("walls_adiabatic", {}).get("flow", {})
        entry.update({k: v.get(k) for k in
                      ("T_bush_peak_nogroove_C", "T_bush_peak_nogroove_theta_deg",
                       "T_mid_peak_nogroove_C", "p_peak_gauge_MPa", "h_min_um",
                       "friction_torque_Nm", "T_shaft_C", "W_N")})
    return entry, prof


def main():
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    out = {"variants": {}, "experiment": {}, "local_bush_reference": {}, "runtime": {}}

    xp, yp = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv")
    xt, yt = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    out["experiment"] = {
        "p_peak_MPa": peak(xp, yp)[1] / 1e6,
        "p_peak_theta_deg": peak(xp, yp)[0],
        "T_peak_C": peak(xt, yt)[1],
        "T_peak_theta_deg": peak(xt, yt)[0],
        "T_inlet_C": yt[0],
        "n_temperature_points": len(xt),
    }

    ref = local_reference()
    if ref is not None:
        entry, prof = ref
        if prof is not None:
            xs, tb = prof
            entry["errors_vs_experiment"] = errors_vs_experiment(
                xs, [v - 273.15 for v in tb], xt, yt)
        out["local_bush_reference"] = entry

    for group in VARIANTS:
        if group not in eq:
            continue
        entry = {}
        for pas in ("pinned", "flow"):
            rec = eq[group][pas]
            m = rec["metrics"]
            xs, ps = profile(rec["dir"], "p")
            _, ts = profile(rec["dir"], "T")
            xb, tb = profile(rec["dir"], "T_bearing")
            _, tj = profile(rec["dir"], "T_journal")
            _, qb = profile(rec["dir"], "q_bearing")
            _, tsh = profile(rec["dir"], "cross_film_t_shaft")
            _, cr = profile(rec["dir"], "cross_film_core_rise")
            th_p, pk_p = peak(xs, ps)
            th_bg, pk_bg = peak_away_from_groove(xb, tb)
            th_tg, pk_tg = peak_away_from_groove(xs, ts)
            entry[pas] = {
                "dir": rec["dir"], "eps": rec["eps"], "psi_deg": rec["psi"],
                "attitude_classical_deg": rec["psi"] - m["force_angle_deg"],
                "W_N": m["W"], "force_angle_deg": m["force_angle_deg"],
                "W_closure_pct": 100.0 * (m["W"] - rec["W_target"]) / rec["W_target"],
                "p_peak_gauge_MPa": (pk_p - P_AMB) / 1e6,
                "p_peak_theta_deg": th_p,
                "T_bearing_peak_C": pk_bg - 273.15,
                "T_bearing_peak_theta_deg": th_bg,
                "T_bearing_min_C": min(tb) - 273.15,
                "T_bearing_max_C": max(tb) - 273.15,
                "T_mid_peak_nogroove_C": pk_tg - 273.15,
                "T_mid_peak_nogroove_theta_deg": th_tg,
                "T_journal_mean_C": sum(tj) / len(tj) - 273.15,
                "q_bearing_min_Wm2": min(qb), "q_bearing_max_Wm2": max(qb),
                "T_shaft_C": tsh[0] - 273.15,
                "core_rise_max_K": max(cr), "core_rise_min_K": min(cr),
                "T_field_max_C": m["T_max"] - 273.15,
                "T_min_C": m["T_min"] - 273.15,
                "mu_min": m.get("mu_min"), "mu_max": m.get("mu_max"),
                "h_min_um": m["h_min"] * 1e6, "cav_frac": m["cav_frac"],
                "friction_torque_Nm": m["friction_torque"],
                "q_supply_m3s": m.get("q_supply"),
                "q_supply_Lmin": (m.get("q_supply") or 0.0) * 6e4,
                "liquid_residual": m.get("liquid_residual"),
                "picard_iters": (m.get("coupling") or {}).get("iters"),
                "picard_R_p": (m.get("coupling") or {}).get("R_p"),
                "picard_R_T": (m.get("coupling") or {}).get("R_T"),
                "newton_steps": rec["newton_steps"],
                "bush_cfg": {k: rec["cfg"][k] for k in
                             ("bush_outer_boundary", "bush_htc")},
                "errors_vs_experiment": errors_vs_experiment(
                    xb, [v - 273.15 for v in tb], xt, yt),
            }
        entry["q_supply_m3s"] = eq[group]["q_supply"]
        out["variants"][group] = entry

    totals = {}
    for log_file in sorted(CASE_DIR.glob("*/fb_log.json")):
        log = json.loads(log_file.read_text())
        runs = log["iters"]
        ok = [r for r in runs if r.get("ok")]
        totals[log_file.parent.name] = {
            "solver_calls": len(runs),
            "wall_s": round(sum(r.get("wall_s", 0.0) for r in runs), 1),
            "picard_max": max(((r.get("coupling") or {}).get("iters", 0) for r in ok),
                              default=0),
            "worst_liquid_residual": max((abs(r.get("liquid_residual", 0.0)) for r in ok),
                                         default=0.0),
        }
    out["runtime"] = totals

    (CASE_DIR / "conjugate_analysis.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
