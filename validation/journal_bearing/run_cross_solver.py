"""Task 1.4 -- pancake against an independent solver at the THD benchmark geometry.

Raimondi & Boyd validates pancake only at `L/D = 1/2` and `1` with `c/R = 1e-3`.
Ferron runs at `L/D = 0.80`, `c/R = 2.9e-3` and Brito at `L/D = 0.80`,
`c/R = 1.71e-3`, so until now nothing independent had ever been run at the
geometry where the load is over-predicted by 52-68%.

This compares pancake's isothermal Elrod-Adams solve against `christopherson.py`
-- a separate Python finite-volume solver with Christopherson clamping, audited
against the R&B charts to 1.3% in S and 0.80 deg in attitude angle.  The two
share no code.

Both are run ISOTHERMAL at the inlet viscosity so the comparison isolates the
Reynolds discretisation, cavitation treatment and force integration.  A
disagreement here would mean a solver defect; agreement moves the entire load
discrepancy onto the reference data or the thermal physics.

Usage:  $PY run_cross_solver.py <path-to-pancake.exe> [workdir]
"""
import json, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import christopherson as CH
import ferron as F
import brito as B
from run_ferron import run_case as ferron_run_case, P_AMB

# Christopherson resolution. The R&B audit runs at 1440x60; these geometries are
# far less anisotropic than the L/D = 8 audit case, so this converges quickly.
CH_NTH, CH_NZ = 1440, 60


def pancake_isothermal(exe, work, tag, *, eps, rpm, groove_deg, nth=720, nz=40):
    """One isothermal pancake solve at Ferron's geometry, via the Ferron harness
    so the config is byte-identical to the validation runs apart from thermal."""
    res = ferron_run_case(exe, work, tag, eps=eps, rpm=rpm, groove_deg=groove_deg,
                          thermal="ISOTHERMAL", h_wall=0.0, nth=nth, nz=nz)
    W = math.hypot(res["fx"], res["fy"])
    mid = res["nz"] // 2
    row = res["p"][mid]
    ip = max(range(len(row)), key=lambda k: row[k])
    return dict(W_N=W,
                attitude_deg=(180.0 - math.degrees(math.atan2(res["fy"], res["fx"]))) % 360.0,
                p_max_MPa=(row[ip] - P_AMB) / 1e6,
                p_max_angle_deg=(ip + 0.5) * 360.0 / res["nth"],
                theta_min=min(min(r) for r in res["theta"]),
                converged=res["converged"], R_p=res["R_p"])


def christopherson_case(*, eps, R, c, L, mu, rpm, groove_deg):
    res = CH.solve(eps=eps, R=R, c=c, L=L, mu=mu,
                   omega=2.0 * math.pi * rpm / 60.0, p_amb=P_AMB,
                   n_theta=CH_NTH, n_z=CH_NZ,
                   groove_deg=groove_deg, groove_width_deg=18.0)
    return dict(W_N=res["W"], attitude_deg=res["attitude_deg"],
                p_max_MPa=res["p_max_gauge"] / 1e6,
                p_max_angle_deg=res["p_max_angle_deg"],
                rupture_deg=res["rupture_deg"],
                converged=res["converged"], sweeps=res["sweeps"])


# The groove sits diametrically opposite the load line in Ferron's frame; the
# pancake harness solves a fixed point for it. Seeded from the committed
# ferron_results.json isothermal case so this script does not have to re-iterate.
FERRON_GROOVE_DEG = 307.32

CASES = [
    # (label, eps, rpm) -- the eccentricities the load comparison is quoted at.
    ("eps=0.47 (2 kN measured)", 0.47, 2000.0),
    ("eps=0.67 (4 kN measured)", 0.67, 2000.0),
    ("eps=0.78 (8 kN measured)", 0.78, 2000.0),
]


def main(exe, work):
    os.makedirs(work, exist_ok=True)
    print("=== pancake vs independent Christopherson solver, Ferron geometry ===")
    print(f"L/D = {F.L / (2 * F.R):.3f}, c/R = {F.C / F.R:.3e}, "
          f"mu = {F.MU_INLET:g} Pa.s, ISOTHERMAL, groove at {FERRON_GROOVE_DEG:.1f} deg")
    print(f"{'case':>26} {'W pancake':>10} {'W indep':>10} {'diff':>7} "
          f"{'p_max pk':>9} {'p_max in':>9} {'diff':>7}", flush=True)

    rows = []
    for label, eps, rpm in CASES:
        pk = pancake_isothermal(exe, work, f"x_eps{eps:.2f}", eps=eps, rpm=rpm,
                                groove_deg=FERRON_GROOVE_DEG)
        ch = christopherson_case(eps=eps, R=F.R, c=F.C, L=F.L, mu=F.MU_INLET,
                                 rpm=rpm, groove_deg=FERRON_GROOVE_DEG)
        dW = (pk["W_N"] - ch["W_N"]) / ch["W_N"]
        dp = (pk["p_max_MPa"] - ch["p_max_MPa"]) / ch["p_max_MPa"]
        rows.append(dict(label=label, eps=eps, rpm=rpm, pancake=pk,
                         independent=ch, dW=dW, dp_max=dp))
        print(f"{label:>26} {pk['W_N']:>10.1f} {ch['W_N']:>10.1f} {dW:>6.1%} "
              f"{pk['p_max_MPa']:>9.4f} {ch['p_max_MPa']:>9.4f} {dp:>6.1%}", flush=True)

    worst = max(abs(r["dW"]) for r in rows)
    print(f"\nworst load disagreement between the two solvers: {worst:.1%}")
    print("Interpretation: a few percent means the discretisation, cavitation and")
    print("force integration are sound at this geometry, and the 52-68% load gap")
    print("against Ferron's MEASURED eccentricity is not a solver defect.")

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "cross_solver_results.json")
    json.dump(rows, open(out, "w"), indent=1)
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "_work/cross")
