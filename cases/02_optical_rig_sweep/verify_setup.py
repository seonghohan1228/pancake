"""Pre-flight gates V1 and V2. Both must pass before RUN1 is dispatched.

V1  grooveless ends-fed anchor. Every verified case in this repository feeds
    through an `inlet_groove`, so removing it is unvalidated. Compare a low-eps
    point against the Ocvirk short-bearing solution, which is the only closed
    form available for this configuration. At L/D ~ 1 short-bearing is a lower
    bound on pressure (it discards circumferential flow), so the gate is order
    of magnitude and sign, not equality.

V2  grid independence. G6 was frozen for Ferron geometry (R = 50 mm, c = 152 um,
    L = 80 mm), not this one. Re-run the worst corner at G6 and G7 and require
    the integral outputs to agree within 2 %.

Usage (from repo root):
    python3 cases/02_optical_rig_sweep/verify_setup.py
"""

import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S
import run_sweep as R
import analyse_sweep as A

VERIFY_ROOT = CASE_DIR / "runs" / "verify"
BAR = 1.0e5
V2_TOLERANCE_PCT = 2.0


def ocvirk(eps, rpm, mu, c=R.C_NOMINAL, L=R.L_NOMINAL, radius=R.R_JOURNAL):
    """Short-bearing peak gauge pressure [Pa] and load [N]."""
    omega = 2.0 * math.pi * rpm / 60.0
    speed = omega * radius

    def shape(angle):
        return eps * math.sin(angle) / (1.0 + eps * math.cos(angle)) ** 3

    peak_shape = max(shape(math.radians(deg)) for deg in range(1, 180))
    p_peak = 3.0 * mu * speed * L * L / (4.0 * radius * c * c) * peak_shape
    load = (mu * omega * radius * L ** 3 / c ** 2
            * eps / (2.0 * (1.0 - eps ** 2) ** 2)
            * math.sqrt(16.0 * eps ** 2 + math.pi ** 2 * (1.0 - eps ** 2)))
    return p_peak, load


def solve(name, eps, rpm, key, grid, mu=None):
    spec = R.point(name, eps, rpm, key, mu=mu if mu else R.MU_20, run="VERIFY")
    run_dir = VERIFY_ROOT / name
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.txt").write_text(R.TEMPLATE.format(
        c=spec["c"], e=spec["eps"] * spec["c"], L=spec["L"],
        n_theta=grid[0], n_z=grid[1], omega=spec["omega"], mu=spec["mu"],
        p_cav=spec["p_cav"], p_end=spec["p_end"], thermal_block=spec["thermal"]))
    rc, wall = S.run_case(run_dir)
    if rc != 0:
        raise SystemExit(f"{name}: solver returned {rc}")
    metrics = S.frame_metrics(run_dir)
    metrics.update(A.cavity_geometry(run_dir))
    metrics["wall_s"] = round(wall, 1)
    R.progress(f"verify {name} ({grid[0]}x{grid[1]}): W={metrics['W']:.1f} N "
               f"p_max={metrics['p_max'] / BAR:.3f} bar")
    return metrics


def relative_pct(a, b):
    return abs(a - b) / max(abs(b), 1e-30) * 100.0


def main():
    report = {}

    # V1 -- grooveless ends-fed anchor at eps = 0.5, 1000 rpm, degassed.
    numeric = solve("v1_g6", 0.5, 1000, "degassed", R.GRID)
    p_peak, load = ocvirk(0.5, 1000, R.MU_20)
    v1 = {
        "solver_W_N": round(numeric["W"], 2),
        "ocvirk_W_N": round(load, 2),
        "W_ratio": round(numeric["W"] / load, 3),
        "solver_p_max_gauge_bar": round((numeric["p_max"] - R.P_AMBIENT) / BAR, 4),
        "ocvirk_p_max_gauge_bar": round(p_peak / BAR, 4),
        "p_ratio": round((numeric["p_max"] - R.P_AMBIENT) / p_peak, 3),
    }
    # Short-bearing overestimates at L/D ~ 1, so the solver should sit below it
    # but within the same order: ratios in [0.2, 1.2].
    v1["pass"] = all(0.2 <= v1[k] <= 1.2 for k in ("W_ratio", "p_ratio"))
    report["V1_grooveless_anchor"] = v1

    # V2 -- grid independence at the worst corner.
    coarse = solve("v2_g6", 0.7, 3000, "degassed", (720, 192))
    fine = solve("v2_g7", 0.7, 3000, "degassed", (1080, 288))
    checks = {
        "W_rel_pct": relative_pct(coarse["W"], fine["W"]),
        "p_max_rel_pct": relative_pct(coarse["p_max"], fine["p_max"]),
        "cav_arc_rel_pct": relative_pct(coarse["cav_arc_deg"], fine["cav_arc_deg"]),
    }
    v2 = {
        "g6": {"W_N": round(coarse["W"], 2),
               "p_max_bar_abs": round(coarse["p_max"] / BAR, 4),
               "cav_arc_deg": coarse["cav_arc_deg"], "wall_s": coarse["wall_s"]},
        "g7": {"W_N": round(fine["W"], 2),
               "p_max_bar_abs": round(fine["p_max"] / BAR, 4),
               "cav_arc_deg": fine["cav_arc_deg"], "wall_s": fine["wall_s"]},
        "checks": {k: round(v, 4) for k, v in checks.items()},
        "limit_pct": V2_TOLERANCE_PCT,
        "pass": all(v <= V2_TOLERANCE_PCT for v in checks.values()),
    }
    report["V2_grid_independence"] = v2

    report["pass"] = v1["pass"] and v2["pass"]
    (CASE_DIR / "verify_setup.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not report["pass"]:
        raise SystemExit("verification gate FAILED -- do not dispatch RUN1")
    R.progress("verify: V1 and V2 both PASS, RUN1 cleared")


if __name__ == "__main__":
    main()
