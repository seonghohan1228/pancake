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


def diagnose(exe, work, eps=0.78, rpm=2000.0):
    """Explain the residual pancake-minus-independent difference.

    MEASURED RESULT (eps = 0.78, 2000 rpm, recorded 2026-08-02).

    A. Grid is the whole story.  At fixed beta = 1e8 pancake converges
       monotonically *down* onto the reference under refinement:

           720x40   15010 N   +6.2%
          1440x60   14717 N   +4.2%
          2880x60   14540 N   +2.9%

       while the reference is already grid-converged by 1440x60 (14129 N,
       unchanged at 2880).  pancake converges more slowly in mesh because its
       Couette term is an upwinded flux on the film content (first order),
       whereas this reference solves for pressure directly with the wedge as a
       source term and has no convective discretisation at all.  Monotone
       convergence toward an independent solution under refinement is the
       signature of discretisation error, not of a defect.

    B. Raising beta makes it WORSE, which was not the expectation.  The guess was
       that pancake would approach the incompressible limit the reference solves;
       instead, at a fixed 2880x60 mesh:

          beta = 1e8   +2.9%   p_max/beta = 0.058
          beta = 3e8   +3.5%   p_max/beta = 0.019
          beta = 1e9   +9.4%   p_max/beta = 0.0059

       Raising beta shrinks `p_max/beta` toward the 4e-3 free-boundary resolution
       floor established by `run_raimondi_boyd.convergence_study()`, so it trades
       a small compressibility difference for a larger under-resolution error.
       The two effects pull in opposite directions and beta = 1e8 sits near the
       optimum for this case.  Do not "improve" a comparison by raising beta.
    """
    import run_ferron as RF

    print(f"\n=== Diagnosing the residual difference (eps = {eps}, {rpm:g} rpm) ===",
          flush=True)
    ch_ref = None
    print("\nA. Grid convergence of each solver at beta = 1e8", flush=True)
    print(f"{'grid':>12} {'W pancake':>10} {'W indep':>10} {'diff':>7} {'p/beta':>8}",
          flush=True)
    rows_grid = []
    for nth, nz in ((720, 40), (1440, 60), (2880, 60)):
        RF.BULK = 1.0e8
        pk = pancake_isothermal(exe, work, f"d_g{nth}", eps=eps, rpm=rpm,
                                groove_deg=FERRON_GROOVE_DEG, nth=nth, nz=nz)
        ch = CH.solve(eps=eps, R=F.R, c=F.C, L=F.L, mu=F.MU_INLET,
                      omega=2.0 * math.pi * rpm / 60.0, p_amb=P_AMB,
                      n_theta=nth, n_z=nz, groove_deg=FERRON_GROOVE_DEG,
                      groove_width_deg=18.0)
        ch_ref = ch
        d = (pk["W_N"] - ch["W"]) / ch["W"]
        rows_grid.append(dict(nth=nth, nz=nz, W_pancake=pk["W_N"], W_indep=ch["W"],
                              diff=d, p_over_beta=pk["p_max_MPa"] * 1e6 / 1.0e8,
                              pancake_converged=pk["converged"],
                              indep_converged=ch["converged"]))
        print(f"{f'{nth}x{nz}':>12} {pk['W_N']:>10.1f} {ch['W']:>10.1f} {d:>6.1%} "
              f"{pk['p_max_MPa'] * 1e6 / 1.0e8:>8.4f}", flush=True)

    print("\nB. Compressibility: raise beta at a fixed 2880x60 mesh", flush=True)
    print(f"{'beta':>10} {'W pancake':>10} {'W indep':>10} {'diff':>7} {'p/beta':>8} "
          f"{'resolved':>9}", flush=True)
    rows_beta = []
    for beta in (1.0e8, 3.0e8, 1.0e9):
        RF.BULK = beta
        pk = pancake_isothermal(exe, work, f"d_b{beta:g}", eps=eps, rpm=rpm,
                                groove_deg=FERRON_GROOVE_DEG, nth=2880, nz=60)
        pob = pk["p_max_MPa"] * 1e6 / beta
        d = (pk["W_N"] - ch_ref["W"]) / ch_ref["W"]
        rows_beta.append(dict(beta=beta, W_pancake=pk["W_N"], W_indep=ch_ref["W"],
                              diff=d, p_over_beta=pob, resolved=bool(pob >= 4e-3),
                              pancake_converged=pk["converged"]))
        print(f"{beta:>10.0e} {pk['W_N']:>10.1f} {ch_ref['W']:>10.1f} {d:>6.1%} "
              f"{pob:>8.4f} {'yes' if pob >= 4e-3 else 'NO':>9}", flush=True)
    RF.BULK = 1.0e8

    print("\nMeasured 2026-08-02: A shows monotone convergence onto the reference")
    print("under mesh refinement, which is discretisation error, not a defect.")
    print("B rises rather than falls -- raising beta pushes p_max/beta toward the")
    print("4e-3 free-boundary floor, trading compressibility for under-resolution.")
    print("beta = 1e8 is near the optimum here; do not raise it to flatter a result.")
    return dict(grid=rows_grid, beta=rows_beta, eps=eps, rpm=rpm)


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

    diag = diagnose(exe, work)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "cross_solver_results.json")
    json.dump(dict(cases=rows, diagnosis=diag), open(out, "w"), indent=1)
    print(f"wrote {out}")
    return rows, diag


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "_work/cross")
