"""Case 01-THD: force-balance search for the Ferron 1983 Fig. 4 / Fig. 5 points.

Unknowns  : (eps = e/c, attitude_angle_deg) of a STATIC bearing.
Target    : Fp = (0, -W) N.  Same empirical sign convention as 01-prelim: the
            moving-bearing reference run carries external_load_y = +W and its
            pressure force sits in the -y half plane, so equilibrium needs
            Fp = -F_ext, i.e. magnitude W at -90 deg.
Method    : damped Newton with a forward-difference Jacobian; every solver call
            (including Jacobian probes) is archived as <group>/fb_iter_NN/.
Tolerance : |W| to 0.1 % and force direction to 0.05 deg.

Two-pass groove closure.  `groove_mixing_model = FLOW_WEIGHTED` needs the
inlet's `supply_flow_rate` Q_s, which the signed-off parameter set does not
supply and the Ferron paper (as digitised here) does not give.  Rather than
invent one, each group is solved twice:

  pass 1  groove_mixing_model = PINNED   (lambda = 1, needs no Q_s)
  pass 2  groove_mixing_model = FLOW_WEIGHTED with Q_s taken from pass 1's own
          steady mass balance, Q_s = inlet_mass_source / rho.  In steady state
          the net liquid injected by the groove penalty cells IS the fresh-oil
          charge (everything else that crosses the groove is recirculated film),
          so this is the solver's self-consistent value, not a fitted one.

Pass 2 restarts from the pass-1 equilibrium, so it normally costs 1-2 Newton
steps.  Both passes are reported; the difference is the hot-oil carry-over
sensitivity.

Usage: python force_balance.py [fig4 | fig5 | all]
"""

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solver_io as S

CASE_DIR = Path(__file__).resolve().parent
PROGRESS = CASE_DIR / "progress.log"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text()

CLEARANCE = 1.52e-4
N_THETA, N_Z = 720, 192
ANGLE_TARGET = -90.0          # deg; Fp direction that balances +W in +y

W_TOL_REL = 1.0e-3            # 0.1 %
ANGLE_TOL_DEG = 0.05

D_EPS = 5.0e-3                # Jacobian probe steps
D_PSI = 0.5

# Journal-side wall bracket (Tier 1, user-signed). The bush side is common:
# bearing_heat_transfer = 7200 W/(m^2 K) at 313.15 K.
JOURNAL = {"Jadiabatic": 0.0, "Jiso40": 10000.0}

Q_SUPPLY_PLACEHOLDER = 1.0e-4  # unused under PINNED; kept positive so the token parses


def log_path(group):
    return CASE_DIR / group / "fb_log.json"


def load_log(group):
    p = log_path(group)
    return json.loads(p.read_text()) if p.exists() else {"iters": []}


def save_log(group, log):
    p = log_path(group)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(log, indent=2))


def evaluate(group, log, cfgkeys, eps, psi, tag):
    """Run one STATIC steady case at (eps, psi); return its metrics record.

    Re-uses an archived iterate when the same (config keys, eps, psi) already ran.
    """
    key = json.dumps(cfgkeys, sort_keys=True)
    for rec in log["iters"]:
        if (rec.get("cfgkey") == key and abs(rec["eps"] - eps) < 1e-12
                and abs(rec["psi"] - psi) < 1e-12 and rec["ok"]):
            return rec
    n = len(log["iters"])
    run_dir = CASE_DIR / group / f"fb_iter_{n:02d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.txt").write_text(TEMPLATE.format(
        e=repr(eps * CLEARANCE), attitude=repr(psi),
        n_theta=N_THETA, n_z=N_Z,
        journal_htc=repr(cfgkeys["journal_htc"]),
        q_supply=repr(cfgkeys["q_supply"]),
        groove_mixing_model=cfgkeys["groove_mixing_model"],
        coupling_max_iters=200, coupling_relaxation=repr(1.0)))
    rc, wall = S.run_case(run_dir)
    m = S.frame_metrics(run_dir) if rc == 0 else None
    rec = {"dir": f"{group}/{run_dir.name}", "cfgkey": key, "eps": eps, "psi": psi,
           "tag": tag, "rc": rc, "wall_s": round(wall, 1),
           "ok": bool(m and "W" in m)}
    if rec["ok"]:
        rec.update({k: m[k] for k in
                    ("Fx", "Fy", "W", "force_angle_deg", "p_max", "h_min",
                     "cav_frac", "friction_torque", "T_max", "T_min",
                     "T_max_theta_deg", "mu_min", "mu_max", "q_supply",
                     "liquid_residual", "coupling") if k in m})
    log["iters"].append(rec)
    save_log(group, log)
    if rec["ok"]:
        cpl = rec.get("coupling") or {}
        S.progress(PROGRESS,
                   f"{rec['dir']} [{tag}] eps={eps:.6f} psi={psi:.4f} -> "
                   f"W={rec['W']:.2f} N angle={rec['force_angle_deg']:.4f} deg "
                   f"Tmax={rec.get('T_max', float('nan'))-273.15:.2f} C "
                   f"picard={cpl.get('iters')} ({wall:.0f} s)")
    else:
        S.progress(PROGRESS, f"{rec['dir']} [{tag}] FAILED rc={rc}")
    return rec


def converged(rec, w_target):
    d_ang = (rec["force_angle_deg"] - ANGLE_TARGET + 180.0) % 360.0 - 180.0
    return (abs(rec["W"] - w_target) / w_target < W_TOL_REL
            and abs(d_ang) < ANGLE_TOL_DEG), d_ang


def newton(group, label, cfgkeys, w_target, eps0, psi0, max_iter=12):
    log = load_log(group)
    target = (0.0, -w_target)
    eps, psi = eps0, psi0
    rec = evaluate(group, log, cfgkeys, eps, psi, f"{label}/x0")
    for it in range(max_iter):
        ok, d_ang = converged(rec, w_target)
        S.progress(PROGRESS, f"{group}/{label} iter {it}: eps={eps:.6f} psi={psi:.4f} "
                             f"W={rec['W']:.2f} (dW={rec['W']-w_target:+.2f} N) "
                             f"dAngle={d_ang:+.4f} deg")
        if ok:
            S.progress(PROGRESS, f"{group}/{label} CONVERGED eps={eps:.6f} psi={psi:.4f}")
            return {"group": group, "label": label, "cfg": cfgkeys, "W_target": w_target,
                    "eps": eps, "psi": psi, "dir": rec["dir"], "metrics": rec,
                    "newton_steps": it}
        r = (rec["Fx"] - target[0], rec["Fy"] - target[1])
        re_ = evaluate(group, log, cfgkeys, eps + D_EPS, psi, f"{label}/J_eps")
        rp_ = evaluate(group, log, cfgkeys, eps, psi + D_PSI, f"{label}/J_psi")
        if not (re_["ok"] and rp_["ok"]):
            raise RuntimeError(f"{group}/{label}: Jacobian probe failed")
        j = [[(re_["Fx"] - rec["Fx"]) / D_EPS, (rp_["Fx"] - rec["Fx"]) / D_PSI],
             [(re_["Fy"] - rec["Fy"]) / D_EPS, (rp_["Fy"] - rec["Fy"]) / D_PSI]]
        det = j[0][0] * j[1][1] - j[0][1] * j[1][0]
        if abs(det) < 1e-12:
            raise RuntimeError(f"{group}/{label}: singular Jacobian")
        d_e = (-r[0] * j[1][1] + r[1] * j[0][1]) / det
        d_p = (-j[0][0] * r[1] + j[1][0] * r[0]) / det
        # Damp so one step cannot leave the physical eps range or spin psi.
        scale = min(1.0, 0.15 / max(abs(d_e), 1e-12), 20.0 / max(abs(d_p), 1e-12))
        eps = min(0.95, max(0.02, eps + scale * d_e))
        psi = psi + scale * d_p
        rec = evaluate(group, log, cfgkeys, eps, psi, f"{label}/x{it+1}")
    raise RuntimeError(f"{group}/{label}: not converged in {max_iter} Newton steps")


def solve_point(group, journal_htc, w_target, eps0, psi0):
    """Two-pass solve: PINNED groove first, then FLOW_WEIGHTED with its own Q_s."""
    pinned = {"journal_htc": journal_htc, "groove_mixing_model": "PINNED",
              "q_supply": Q_SUPPLY_PLACEHOLDER}
    a = newton(group, "pinned", pinned, w_target, eps0, psi0)
    q_s = a["metrics"]["q_supply"]
    S.progress(PROGRESS, f"{group}: derived Q_supply = {q_s:.6e} m^3/s "
                         f"({q_s*6e4:.2f} L/min) from the PINNED steady mass balance")
    flow = {"journal_htc": journal_htc, "groove_mixing_model": "FLOW_WEIGHTED",
            "q_supply": q_s}
    b = newton(group, "flow", flow, w_target, a["eps"], a["psi"])
    return {"pinned": a, "flow": b, "q_supply": q_s}


def store(results):
    (CASE_DIR / "equilibrium.json").write_text(json.dumps(results, indent=2))


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out_path = CASE_DIR / "equilibrium.json"
    results = json.loads(out_path.read_text()) if out_path.exists() else {}

    if which in ("fig4", "all"):
        for variant, htc in JOURNAL.items():
            group = f"fig4_{variant}"
            if group in results:
                continue
            results[group] = solve_point(group, htc, 4000.0, 0.60, -40.0)
            store(results)

    if which in ("fig5", "all"):
        # J-adiabatic only. 4 kN is already solved as fig4_Jadiabatic and is
        # reused rather than re-run.
        guesses = {2000.0: (0.45, -47.0), 3000.0: (0.55, -43.0),
                   6000.0: (0.72, -33.0), 8000.0: (0.78, -29.0),
                   10000.0: (0.82, -26.0)}
        for w, (e0, p0) in sorted(guesses.items()):
            group = f"fig5_W{int(w)}"
            if group in results:
                continue
            results[group] = solve_point(group, JOURNAL["Jadiabatic"], w, e0, p0)
            store(results)

    print(json.dumps({k: {"eps": v["flow"]["eps"], "psi": v["flow"]["psi"]}
                      for k, v in results.items()}, indent=2))


if __name__ == "__main__":
    main()
