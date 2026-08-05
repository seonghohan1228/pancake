"""Case 01: 2D force-balance search for the Ferron 1983 Fig. 6 operating point.

Unknowns  : (eps = e/c, attitude_angle_deg) of a STATIC bearing.
Target    : the steady pressure-force resultant that balances the +6000 N
            y-load used by the moving-bearing run dt_T1_fbns, i.e.
            Fp = (0, -6000) N.  The sign is empirical: dt_T1_fbns runs under
            external_load_y = +6000 and its pressure force sits in the -y half
            plane (Fp = (-1001, -3046) N at t = 0.0135 s, angle -108 deg),
            so equilibrium requires Fp = -F_ext = (0, -6000) N.
Method    : damped Newton with a forward-difference Jacobian; every solver call
            (including Jacobian probes) is archived as fb_iter_NN/.
Tolerance : |W| to 0.1 % (6 N) and force direction to 0.05 deg.

Usage: python force_balance.py [mu40|mu55|all]
"""

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solver_io as S

CASE_DIR = Path(__file__).resolve().parent
PROGRESS = CASE_DIR / "progress.log"
LOG_JSON = CASE_DIR / "fb_log.json"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text()

CLEARANCE = 1.52e-4
W_TARGET = 6000.0
ANGLE_TARGET = -90.0          # deg; Fp direction that balances +6000 N in +y
TARGET = (0.0, -W_TARGET)     # Cartesian target for the Newton residual

W_TOL_REL = 1.0e-3            # 0.1 %
ANGLE_TOL_DEG = 0.05

D_EPS = 5.0e-3                # Jacobian probe steps
D_PSI = 0.5

VISCOSITIES = {"mu40": 0.0277, "mu55": 0.0166}


def load_log():
    return json.loads(LOG_JSON.read_text()) if LOG_JSON.exists() else {"iters": []}


def save_log(log):
    LOG_JSON.write_text(json.dumps(log, indent=2))


def evaluate(log, mu, eps, psi, tag):
    """Run one STATIC steady case at (eps, psi); return its metrics dict.

    Re-uses an archived iterate when the same (mu, eps, psi) was already run.
    """
    for rec in log["iters"]:
        if (rec["mu"] == mu and abs(rec["eps"] - eps) < 1e-12
                and abs(rec["psi"] - psi) < 1e-12 and rec["ok"]):
            return rec
    n = len(log["iters"])
    run_dir = CASE_DIR / f"fb_iter_{n:02d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.txt").write_text(
        TEMPLATE.format(e=repr(eps * CLEARANCE), attitude=repr(psi), mu=repr(mu)))
    rc, wall = S.run_case(run_dir)
    m = S.frame_metrics(run_dir) if rc == 0 else None
    rec = {"dir": run_dir.name, "mu": mu, "eps": eps, "psi": psi, "tag": tag,
           "rc": rc, "wall_s": round(wall, 1), "ok": bool(m and "W" in m)}
    if rec["ok"]:
        rec.update({k: m[k] for k in
                    ("Fx", "Fy", "W", "force_angle_deg", "p_max", "h_min",
                     "cav_frac", "friction_torque") if k in m})
    log["iters"].append(rec)
    save_log(log)
    if rec["ok"]:
        S.progress(PROGRESS, f"{run_dir.name} [{tag}] mu={mu} eps={eps:.6f} "
                             f"psi={psi:.4f} -> W={rec['W']:.2f} N "
                             f"angle={rec['force_angle_deg']:.4f} deg ({wall:.0f} s)")
    else:
        S.progress(PROGRESS, f"{run_dir.name} [{tag}] FAILED rc={rc}")
    return rec


def converged(rec):
    d_ang = (rec["force_angle_deg"] - ANGLE_TARGET + 180.0) % 360.0 - 180.0
    return (abs(rec["W"] - W_TARGET) / W_TARGET < W_TOL_REL
            and abs(d_ang) < ANGLE_TOL_DEG), d_ang


def solve(label, mu, eps0, psi0, max_iter=12):
    log = load_log()
    eps, psi = eps0, psi0
    rec = evaluate(log, mu, eps, psi, f"{label}/x0")
    for it in range(max_iter):
        ok, d_ang = converged(rec)
        S.progress(PROGRESS, f"{label} iter {it}: eps={eps:.6f} psi={psi:.4f} "
                             f"W={rec['W']:.2f} (dW={rec['W']-W_TARGET:+.2f} N) "
                             f"dAngle={d_ang:+.4f} deg")
        if ok:
            S.progress(PROGRESS, f"{label} CONVERGED at eps={eps:.6f} psi={psi:.4f}")
            return {"label": label, "mu": mu, "eps": eps, "psi": psi,
                    "dir": rec["dir"], "metrics": rec, "iters": it}
        r = (rec["Fx"] - TARGET[0], rec["Fy"] - TARGET[1])
        re_ = evaluate(log, mu, eps + D_EPS, psi, f"{label}/J_eps")
        rp_ = evaluate(log, mu, eps, psi + D_PSI, f"{label}/J_psi")
        if not (re_["ok"] and rp_["ok"]):
            raise RuntimeError(f"{label}: Jacobian probe failed")
        j = [[(re_["Fx"] - rec["Fx"]) / D_EPS, (rp_["Fx"] - rec["Fx"]) / D_PSI],
             [(re_["Fy"] - rec["Fy"]) / D_EPS, (rp_["Fy"] - rec["Fy"]) / D_PSI]]
        det = j[0][0] * j[1][1] - j[0][1] * j[1][0]
        if abs(det) < 1e-12:
            raise RuntimeError(f"{label}: singular Jacobian")
        d_e = (-r[0] * j[1][1] + r[1] * j[0][1]) / det
        d_p = (-j[0][0] * r[1] + j[1][0] * r[0]) / det
        # Damp so one step cannot leave the physical eps range or spin psi.
        scale = min(1.0, 0.15 / max(abs(d_e), 1e-12), 20.0 / max(abs(d_p), 1e-12))
        eps_new = min(0.98, max(0.02, eps + scale * d_e))
        psi_new = psi + scale * d_p
        eps, psi = eps_new, psi_new
        rec = evaluate(log, mu, eps, psi, f"{label}/x{it+1}")
    raise RuntimeError(f"{label}: not converged in {max_iter} Newton steps")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    targets = ["mu40", "mu55"] if which == "all" else [which]
    guesses = {"mu40": (0.50, -45.0), "mu55": (0.62, -40.0)}
    out_path = CASE_DIR / "equilibrium.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    for label in targets:
        e0, p0 = guesses[label]
        out[label] = solve(label, VISCOSITIES[label], e0, p0)
        out_path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
