"""Independent finite-volume Reynolds solver with Christopherson cavitation.

Why this exists
---------------
Every verdict in this directory compares pancake against *measured* data, so a
disagreement has at least three possible owners: the solver, the digitised
reference, and the modelling assumptions.  The Raimondi & Boyd charts separate
the first from the rest, but only at `L/D = 1/2` and `1` with `c/R = 1e-3`.  The
THD benchmarks run at `L/D = 0.8` with `c/R = 2.9e-3` (Ferron) and `1.71e-3`
(Brito), which is outside the regime where the solver has ever been checked
against anything independent.

This module closes that gap.  It solves the same PDE by an entirely separate
route -- Python, finite volume, red-black SOR, Christopherson clamping -- so
agreement with pancake at the benchmark geometry rules out the discretisation
and localises any remaining error in the reference data or the physics.

It is deliberately *not* a copy of pancake's algorithm.  pancake solves the
Elrod-Adams film content with a compressibility (`beta`) and a switch function;
this solves for pressure directly and clamps it, which is the classical
Christopherson (1941) route to the Swift-Stieber boundary condition.  The two
share no code and no failure mode: in particular this solver has no `p_max/beta`
resolution requirement, so it is a clean check on that known trap.

What it does and does not model
-------------------------------
  * Steady, isoviscous, incompressible, rigid, laminar.
  * Christopherson clamping gives the Swift-Stieber (Reynolds) rupture condition
    -- correct rupture boundary, zero pressure gradient there.
  * It is NOT mass-conserving in the cavitated zone: there is no reformation
    model, exactly like the charts it reproduces.  Do not use it to judge
    reformation, void fraction, or anything JFO-specific.  Use it for load,
    peak pressure, attitude angle and rupture position.

Self-validation
---------------
`audit()` checks this solver against two references it cannot have been tuned
to: the exact infinite-length Swift-Stieber solution (`swift_stieber_infinite`
in `raimondi_boyd.py`, derived analytically) and the transcribed Raimondi & Boyd
chart values.  Run it before trusting any comparison made with this module.

Usage:  $PY christopherson.py            # runs audit()
"""
import math

import numpy as np

from analytic import sommerfeld_number
from raimondi_boyd import rb_lookup, swift_stieber_infinite


def solve(*, eps, R, c, L, mu, omega, p_amb=1.0e5,
          n_theta=720, n_z=40, groove_deg=None, groove_width_deg=4.0,
          p_supply=None, omega_relax=1.85, tol=1e-9, max_sweeps=60000,
          check_every=1):
    """Solve the steady isoviscous Reynolds equation with Christopherson cavitation.

    Film thickness follows pancake's convention with `attitude_angle_deg = 180`:
    `h = c (1 + eps cos theta)`, so the minimum gap sits at theta = 180 deg and
    theta increases with the entrainment direction.

    `groove_deg` places a full-width Dirichlet supply groove at that angle (the
    fed-film condition the Raimondi & Boyd charts assume); leave it None for an
    ungrooved bearing.  `p_supply` defaults to ambient.

    Returns a dict with the pressure field [Pa absolute], load, attitude angle,
    peak pressure and its location, rupture angle, and the convergence evidence.
    """
    if p_supply is None:
        p_supply = p_amb

    d_theta = 2.0 * math.pi / n_theta
    d_z = L / n_z
    ds = R * d_theta                       # arc length step
    U = omega * R

    theta = (np.arange(n_theta) + 0.5) * d_theta
    h = c * (1.0 + eps * np.cos(theta))            # cell-centred film thickness
    h_e = 0.5 * (h + np.roll(h, -1))               # east face of cell i

    # Face conductances of the Poiseuille term, G = h^3/(12 mu) * (face area)/(distance).
    Ge = h_e ** 3 / (12.0 * mu) * d_z / ds                  # (n_theta,)
    Gz = h ** 3 / (12.0 * mu) * ds / d_z                    # interior z faces
    Gz_bnd = h ** 3 / (12.0 * mu) * ds / (0.5 * d_z)        # to the z boundary

    aE = np.repeat(Ge[:, None], n_z, axis=1)
    aW = np.repeat(np.roll(Ge, 1)[:, None], n_z, axis=1)
    aN = np.repeat(Gz[:, None], n_z, axis=1)
    aS = np.repeat(Gz[:, None], n_z, axis=1)
    aN[:, -1] = 0.0
    aS[:, 0] = 0.0

    # Couette wedge: net outflow (U/2)(h_e - h_w) dz moved to the right-hand side.
    couette = -(0.5 * U) * (h_e - np.roll(h_e, 1)) * d_z
    b = np.repeat(couette[:, None], n_z, axis=1)

    aP = aE + aW + aN + aS
    # Dirichlet ambient at both axial ends.
    aP[:, 0] += Gz_bnd
    aP[:, -1] += Gz_bnd
    b[:, 0] += Gz_bnd * p_amb
    b[:, -1] += Gz_bnd * p_amb

    groove_mask = np.zeros((n_theta, n_z), dtype=bool)
    if groove_deg is not None:
        dth = np.abs((np.degrees(theta) - groove_deg + 180.0) % 360.0 - 180.0)
        groove_mask[dth <= 0.5 * groove_width_deg, :] = True

    p = np.full((n_theta, n_z), p_amb, dtype=float)
    p[groove_mask] = p_supply

    # Red-black ordering so the Gauss-Seidel sweep vectorises.
    ii, jj = np.meshgrid(np.arange(n_theta), np.arange(n_z), indexing="ij")
    colors = [(ii + jj) % 2 == 0, (ii + jj) % 2 == 1]

    # Point SOR converges in O(n_theta) sweeps, and much worse when the mesh is
    # strongly anisotropic: at L/D = 8 the theta:z conductance ratio is ~4e3 and a
    # 1440x60 grid needs more than 4e4 sweeps. Long-bearing comparisons should use
    # a coarse z grid (the solution is nearly z-uniform there anyway) rather than a
    # looser tolerance, which would silently bias the rupture angle.
    sweeps, residual, p_prev = 0, float("inf"), p.copy()
    for sweeps in range(1, max_sweeps + 1):
        checking = (sweeps % check_every == 0) or sweeps == max_sweeps
        if checking:
            np.copyto(p_prev, p)
        for mask in colors:
            pE = np.roll(p, -1, axis=0)
            pW = np.roll(p, 1, axis=0)
            pN = np.empty_like(p); pN[:, :-1] = p[:, 1:]; pN[:, -1] = 0.0
            pS = np.empty_like(p); pS[:, 1:] = p[:, :-1]; pS[:, 0] = 0.0
            p_gs = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP
            upd = mask & ~groove_mask
            p[upd] = (1.0 - omega_relax) * p[upd] + omega_relax * p_gs[upd]
            # Christopherson: the film cannot sustain pressure below the cavity
            # value, and clamping inside the sweep (not after it) is what makes
            # the converged rupture boundary satisfy Swift-Stieber.
            np.maximum(p, p_amb, out=p)
            p[groove_mask] = p_supply

        if checking:
            residual = float(np.max(np.abs(p - p_prev))) / max(1.0, float(np.max(p)) - p_amb)
            if residual < tol:
                break

    p_gauge = p - p_amb
    fx = float(np.sum(p_gauge * np.cos(theta)[:, None]) * R * d_theta * d_z)
    fy = float(np.sum(p_gauge * np.sin(theta)[:, None]) * R * d_theta * d_z)
    W = math.hypot(fx, fy)

    mid = n_z // 2
    i_peak = int(np.argmax(p_gauge[:, mid]))

    # First angle past the minimum gap at which the mid-plane film cavitates.
    rupture = float("nan")
    for i in range(n_theta // 2, n_theta):
        if p_gauge[i, mid] <= 1e-9:
            rupture = (i + 0.5) * 360.0 / n_theta
            break

    return dict(p=p, theta=theta, h=h, fx=fx, fy=fy, W=W,
                attitude_deg=(180.0 - math.degrees(math.atan2(fy, fx))) % 360.0,
                p_max_gauge=float(np.max(p_gauge[:, mid])),
                p_max_angle_deg=(i_peak + 0.5) * 360.0 / n_theta,
                rupture_deg=rupture, sweeps=sweeps, residual=float(residual),
                converged=bool(residual < tol), n_theta=n_theta, n_z=n_z, L=L)


# --- self-validation --------------------------------------------------------
# Geometry of the Raimondi & Boyd harness, so the two agree by construction on
# everything except the algorithm (see run_raimondi_boyd.py).
RB_R, RB_C, RB_MU, RB_N_REV_S = 0.025, 25e-6, 0.01, 5.0


def audit(n_theta=1440, n_z=60, n_z_long=20, verbose=True):
    """Check this solver against references it cannot have been tuned to.

    Returns (rows, worst_S_error, worst_phi_error_deg).
    """
    omega = 2.0 * math.pi * RB_N_REV_S
    rows = []

    if verbose:
        print("=== Christopherson solver audit ===", flush=True)
        print("A. Exact infinite-length Swift-Stieber (analytic, raimondi_boyd.py)")
        print(f"{'eps':>5} {'S exact':>9} {'S here':>9} {'dS':>7} "
              f"{'phi exact':>10} {'phi here':>9} {'dphi':>7}", flush=True)

    # A. Long-bearing limit: L/D = 8 approaches the analytic infinite solution.
    # Coarse z on purpose -- see the anisotropy note in solve().
    for eps in (0.4, 0.6, 0.8):
        exact = swift_stieber_infinite(eps)
        res = solve(eps=eps, R=RB_R, c=RB_C, L=8.0 * 2.0 * RB_R, mu=RB_MU,
                    omega=omega, n_theta=n_theta, n_z=n_z_long, groove_deg=0.0)
        S = sommerfeld_number(mu=RB_MU, N_rev_s=RB_N_REV_S, W=res["W"],
                              L=res["L"], D=2.0 * RB_R, R=RB_R, c=RB_C)
        dS = (S - exact["S"]) / exact["S"]
        dphi = res["attitude_deg"] - exact["phi_deg"]
        rows.append(dict(case="infinite", eps=eps, S=S, S_ref=exact["S"],
                         dS=dS, phi=res["attitude_deg"], phi_ref=exact["phi_deg"],
                         dphi=dphi, converged=res["converged"]))
        if verbose:
            print(f"{eps:>5.1f} {exact['S']:>9.5f} {S:>9.5f} {dS:>6.1%} "
                  f"{exact['phi_deg']:>10.2f} {res['attitude_deg']:>9.2f} {dphi:>6.2f}", flush=True)

    # B. Finite-width charts.
    if verbose:
        print("\nB. Raimondi & Boyd chart values (transcribed, audited to 0.6%)")
        print(f"{'L/D':>5} {'eps':>5} {'S chart':>9} {'S here':>9} {'dS':>7} "
              f"{'phi chart':>10} {'phi here':>9} {'dphi':>7}")
    for ld_key, ld in (("0.5", 0.5), ("1", 1.0)):
        for eps in (0.4, 0.6, 0.8):
            ref = rb_lookup(ld_key, eps)
            if ref is None:
                continue
            S_ref, phi_ref = ref
            res = solve(eps=eps, R=RB_R, c=RB_C, L=ld * 2.0 * RB_R, mu=RB_MU,
                        omega=omega, n_theta=n_theta, n_z=n_z, groove_deg=0.0)
            S = sommerfeld_number(mu=RB_MU, N_rev_s=RB_N_REV_S, W=res["W"],
                                  L=res["L"], D=2.0 * RB_R, R=RB_R, c=RB_C)
            dS = (S - S_ref) / S_ref
            dphi = res["attitude_deg"] - phi_ref
            rows.append(dict(case=f"L/D={ld_key}", eps=eps, S=S, S_ref=S_ref,
                             dS=dS, phi=res["attitude_deg"], phi_ref=phi_ref,
                             dphi=dphi, converged=res["converged"]))
            if verbose:
                print(f"{ld_key:>5} {eps:>5.1f} {S_ref:>9.5f} {S:>9.5f} {dS:>6.1%} "
                      f"{phi_ref:>10.2f} {res['attitude_deg']:>9.2f} {dphi:>6.2f}", flush=True)

    finite = [r for r in rows if r["case"] != "infinite"]
    worst_S = max(abs(r["dS"]) for r in finite)
    worst_phi = max(abs(r["dphi"]) for r in finite)
    if verbose:
        print(f"\nworst finite-width error: {worst_S:.1%} in S, {worst_phi:.2f} deg "
              f"in attitude angle")
        print("(criteria used for pancake: 10% and 5 deg)")
        if not all(r["converged"] for r in rows):
            print("WARNING: not every audit case reached the SOR tolerance")
    return rows, worst_S, worst_phi


if __name__ == "__main__":
    audit()
