"""Case 01-THD-crossfilm: pull the reported numbers out of the archived equilibria.

Writes thd_analysis.json and prints the tables that go into the summary.

The new quantity relative to case 01-THD is the BUSH-SIDE FILM SURFACE
temperature `cross_film_t_bush`, which the cross-film model makes available for
the first time. That is the closest model quantity to Ferron's thermocouples,
which sit 0.5 mm below the bush inner surface; the depth-averaged film
temperature is also reported so the two can be compared directly.
"""

import csv
import json
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S

REF = CASE_DIR.parent / "ref_data"
GROOVE_DEG = 90.0
P_AMB = 101325.0

VARIANTS = ("walls_adiabatic", "bush40_Jadiabatic", "bush40_J40")

# The depth-averaged J-adiabatic equilibrium of case 01-THD, for the comparison
# row. Copied from cases/01_ferron_thd/thd_analysis.json (flow pass).
DEPTH_AVERAGED_REF = {
    "case": "01_ferron_thd fig4_Jadiabatic (flow), cross_film = NONE",
    "eps": 0.5738030720305524,
    "psi_deg": -36.32370169904824,
    "W_N": 3999.70,
    "T_field_max_C": 43.92,
    "T_mid_peak_C": 43.78,
    "T_mid_peak_theta_deg": 281.8,
    "h_min_um": 64.78,
    "p_peak_gauge_MPa": 1.2702,
    "friction_torque_Nm": 2.808,
    "mu_min": 0.0238061,
    "mu_max": 0.0274285,
}


def read_csv(path):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]))
    return xs, ys


# The 18 deg groove is centred on x = 0 (= 360). Its two edge cells carry a
# one-cell spike/dip in the wall temperatures -- the groove pins the film
# temperature, and the conjugate closure reacts to that step over a single cell.
# Reporting the raw field maximum would report that spike rather than the thermal
# crest, so every peak is also quoted with a band around the groove removed.
GROOVE_EXCLUDE_DEG = 12.0


def peak(xs, ys):
    i = max(range(len(ys)), key=lambda k: ys[k])
    return xs[i], ys[i]


def peak_away_from_groove(xs, ys, band=GROOVE_EXCLUDE_DEG):
    keep = [(x, y) for x, y in zip(xs, ys) if band <= x <= 360.0 - band]
    if not keep:
        return peak(xs, ys)
    i = max(range(len(keep)), key=lambda k: keep[k][1])
    return keep[i]


def range_away_from_groove(xs, ys, band=GROOVE_EXCLUDE_DEG):
    keep = [y for x, y in zip(xs, ys) if band <= x <= 360.0 - band]
    return (min(keep), max(keep)) if keep else (min(ys), max(ys))


def profile(run_dir, field):
    got = S.midplane_profile(CASE_DIR / run_dir, field)
    if got is None:
        return None
    theta, v = got
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def main():
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    out = {"variants": {}, "experiment": {}, "depth_averaged_reference": DEPTH_AVERAGED_REF,
           "runtime": {}}

    xp, yp = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv")
    xt, yt = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    out["experiment"]["p_peak_MPa"] = peak(xp, yp)[1] / 1e6
    out["experiment"]["p_peak_theta_deg"] = peak(xp, yp)[0]
    out["experiment"]["T_peak_C"] = peak(xt, yt)[1]
    out["experiment"]["T_peak_theta_deg"] = peak(xt, yt)[0]
    out["experiment"]["T_inlet_C"] = yt[0]

    for group in VARIANTS:
        if group not in eq:
            continue
        entry = {}
        for pas in ("pinned", "flow"):
            rec = eq[group][pas]
            m = rec["metrics"]
            xs, ps = profile(rec["dir"], "p")
            _, ts = profile(rec["dir"], "T")
            xb, tb = profile(rec["dir"], "cross_film_t_bush")
            _, tsh = profile(rec["dir"], "cross_film_t_shaft")
            _, cr = profile(rec["dir"], "cross_film_core_rise")
            th_p, pk_p = peak(xs, ps)
            th_t, pk_t = peak(xs, ts)
            th_b, pk_b = peak(xb, tb)
            th_tg, pk_tg = peak_away_from_groove(xs, ts)
            th_bg, pk_bg = peak_away_from_groove(xb, tb)
            tb_lo_g, tb_hi_g = range_away_from_groove(xb, tb)
            entry[pas] = {
                "dir": rec["dir"], "eps": rec["eps"], "psi_deg": rec["psi"],
                "attitude_classical_deg": rec["psi"] - m["force_angle_deg"],
                "W_N": m["W"], "force_angle_deg": m["force_angle_deg"],
                "W_closure_pct": 100.0 * (m["W"] - rec["W_target"]) / rec["W_target"],
                "p_peak_gauge_MPa": (pk_p - P_AMB) / 1e6,
                "p_peak_theta_deg": th_p,
                "p_max_gauge_MPa": (m["p_max"] - P_AMB) / 1e6,
                "T_mid_peak_C": pk_t - 273.15, "T_mid_peak_theta_deg": th_t,
                "T_bush_peak_C": pk_b - 273.15, "T_bush_peak_theta_deg": th_b,
                "T_bush_min_C": min(tb) - 273.15,
                # Groove edge cells removed; these are the reported numbers.
                "T_mid_peak_nogroove_C": pk_tg - 273.15,
                "T_mid_peak_nogroove_theta_deg": th_tg,
                "T_bush_peak_nogroove_C": pk_bg - 273.15,
                "T_bush_peak_nogroove_theta_deg": th_bg,
                "T_bush_min_nogroove_C": tb_lo_g - 273.15,
                "T_bush_max_nogroove_C": tb_hi_g - 273.15,
                "T_shaft_C": tsh[0] - 273.15,
                "core_rise_max_K": max(cr), "core_rise_min_K": min(cr),
                "T_field_max_C": m["T_max"] - 273.15,
                "T_field_max_theta_deg": m.get("T_max_theta_deg"),
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
                "wall_cfg": {k: rec["cfg"][k] for k in
                             ("shaft_htc", "bush_htc", "bush_k")},
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

    (CASE_DIR / "thd_analysis.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
