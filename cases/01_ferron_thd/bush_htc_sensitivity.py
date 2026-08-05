"""Diagnostic ONLY - not a deliverable and not used in any reported figure.

The Tier-1 bush wall closes the film's heat balance with a single lumped
coefficient h_b = K_b / (R1 ln(R2/R1)).  The signed-off value uses
K_b = 250 W/(m K), giving h_b = 7213 W/(m^2 K).  The Part A result under-predicts
the measured film temperature rise, so this script re-runs the CONVERGED
J-adiabatic geometry (eps and attitude frozen at their force-balanced values) at
two lower bush conductivities to show how much of the gap the wall closure owns.

Nothing here is fitted: eps/psi are not re-solved, the load is therefore no
longer exactly 4 kN, and the numbers are reported as a sensitivity, not a result.

Usage: python bush_htc_sensitivity.py
"""

import json
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S

TEMPLATE = (CASE_DIR / "config_template.txt").read_text()
CLEARANCE = 1.52e-4
GEOM_FACTOR = 0.05 * 0.6931471805599453  # R1 ln(R2/R1) with R1 = 0.05, R2 = 0.10

# K_b [W/(m K)] -> label.  250 is the signed-off value; 120 and 50 bracket
# leaded bronze and steel, the usual bush materials.
K_BUSH = {250.0: "signed_off", 120.0: "bronze_like", 50.0: "steel_like"}


def main():
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    base = eq["fig4_Jadiabatic"]["flow"]
    eps, psi = base["eps"], base["psi"]
    q_s = eq["fig4_Jadiabatic"]["q_supply"]
    out = {}
    for k_b, label in K_BUSH.items():
        h_b = k_b / GEOM_FACTOR
        run_dir = CASE_DIR / "sensitivity_bush" / f"kb_{int(k_b)}"
        run_dir.mkdir(parents=True, exist_ok=True)
        cfg = TEMPLATE.format(
            e=repr(eps * CLEARANCE), attitude=repr(psi), n_theta=720, n_z=192,
            journal_htc=repr(0.0), q_supply=repr(q_s),
            groove_mixing_model="FLOW_WEIGHTED",
            coupling_max_iters=200, coupling_relaxation=repr(1.0))
        cfg = cfg.replace("bearing_heat_transfer = 7200",
                          f"bearing_heat_transfer = {h_b!r}")
        (run_dir / "config.txt").write_text(cfg)
        rc, wall = S.run_case(run_dir)
        m = S.frame_metrics(run_dir) if rc == 0 else None
        rec = {"K_b_W_mK": k_b, "h_b_W_m2K": h_b, "label": label, "rc": rc,
               "wall_s": round(wall, 1)}
        if m:
            rec.update({"W_N": m["W"], "T_max_C": m["T_max"] - 273.15,
                        "T_min_C": m["T_min"] - 273.15,
                        "T_rise_K": m["T_max"] - m["T_min"],
                        "p_max_gauge_MPa": (m["p_max"] - 101325.0) / 1e6,
                        "mu_min": m.get("mu_min"),
                        "picard_iters": (m.get("coupling") or {}).get("iters")})
        out[str(int(k_b))] = rec
        S.progress(CASE_DIR / "progress.log",
                   f"sensitivity_bush/kb_{int(k_b)}: h_b={h_b:.0f} W/(m^2 K) -> "
                   f"Tmax={rec.get('T_max_C', float('nan')):.2f} C "
                   f"W={rec.get('W_N', float('nan')):.0f} N ({wall:.0f} s)")
    (CASE_DIR / "sensitivity_bush.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
