"""Case 01-THD: pull the reported numbers out of the archived equilibria.

Writes thd_analysis.json and prints the tables that go into SUMMARY.md.
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
ECC_RESCALE = 145.0 / 152.0


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


def profile(run_dir, field):
    theta, v = S.midplane_profile(CASE_DIR / run_dir, field)
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def main():
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    out = {"fig4": {}, "fig5": {}, "experiment": {}}

    xp, yp = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv")
    xt, yt = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    out["experiment"]["p_peak_MPa"], out["experiment"]["p_peak_theta_deg"] = \
        peak(xp, yp)[1] / 1e6, peak(xp, yp)[0]
    out["experiment"]["T_peak_C"], out["experiment"]["T_peak_theta_deg"] = \
        peak(xt, yt)[1], peak(xt, yt)[0]
    out["experiment"]["T_inlet_C"] = yt[0]

    for group in ("fig4_Jadiabatic", "fig4_Jiso40"):
        if group not in eq:
            continue
        entry = {}
        for pas in ("pinned", "flow"):
            rec = eq[group][pas]
            m = rec["metrics"]
            xs, ps = profile(rec["dir"], "p")
            xs_t, ts = profile(rec["dir"], "T")
            th_p, pk_p = peak(xs, ps)
            th_t, pk_t = peak(xs_t, ts)
            entry[pas] = {
                "dir": rec["dir"], "eps": rec["eps"], "psi_deg": rec["psi"],
                "attitude_classical_deg": rec["psi"] - m["force_angle_deg"],
                "W_N": m["W"], "force_angle_deg": m["force_angle_deg"],
                "p_peak_gauge_MPa": (pk_p - P_AMB) / 1e6,
                "p_peak_theta_deg": th_p,
                "p_max_gauge_MPa": (m["p_max"] - P_AMB) / 1e6,
                "T_mid_peak_C": pk_t - 273.15, "T_mid_peak_theta_deg": th_t,
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
            }
        entry["q_supply_m3s"] = eq[group]["q_supply"]
        out["fig4"][group] = entry

    xe, ye = read_csv(REF / "ferron1983_fig5_2000rpm_eccentricity.csv")
    out["experiment"]["fig5_rescaled"] = [
        {"W_N": x, "eps_145": y, "eps_rescaled": y * ECC_RESCALE}
        for x, y in zip(xe, ye)]

    for group, rec in eq.items():
        if group == "fig4_Jadiabatic":
            w = 4000.0
        elif group.startswith("fig5_W"):
            w = float(group[len("fig5_W"):])
        else:
            continue
        m = rec["flow"]["metrics"]
        out["fig5"][str(int(w))] = {
            "group": group, "W_N": w, "eps": rec["flow"]["eps"],
            "psi_deg": rec["flow"]["psi"],
            "attitude_classical_deg": rec["flow"]["psi"] - m["force_angle_deg"],
            "W_solved_N": m["W"], "force_angle_deg": m["force_angle_deg"],
            "h_min_um": m["h_min"] * 1e6,
            "p_max_gauge_MPa": (m["p_max"] - P_AMB) / 1e6,
            "T_field_max_C": m["T_max"] - 273.15,
            "cav_frac": m["cav_frac"],
            "friction_torque_Nm": m["friction_torque"],
            "q_supply_Lmin": rec["q_supply"] * 6e4,
            "picard_iters": (m.get("coupling") or {}).get("iters"),
            "newton_steps": rec["flow"]["newton_steps"],
        }

    # Wall clock and iteration health across every archived solver call.
    totals = {}
    for group_dir in sorted(CASE_DIR.glob("fig*/fb_log.json")):
        log = json.loads(group_dir.read_text())
        runs = [r for r in log["iters"]]
        totals[group_dir.parent.name] = {
            "solver_calls": len(runs),
            "wall_s": round(sum(r.get("wall_s", 0.0) for r in runs), 1),
            "picard_max": max((r.get("coupling") or {}).get("iters", 0)
                              for r in runs if r.get("ok")),
        }
    out["runtime"] = totals

    (CASE_DIR / "thd_analysis.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
