# Case 00c: JFO solver verification (correctness, not grid independence)

Triggered by 00a findings: rupture-front position drifts with grid at beta = 1e9
while Sommerfeld/Gumbel ladders converge cleanly. Question: is the front placement
(and BC/flux handling) correct, or only self-consistent? Case 00 grid work is
PAUSED until this closes.

## Tests

### 1. Exact 1D benchmark (this directory, no external data)

Infinitely long journal bearing emulated with NEUMANN z-boundaries and a tiny n_z:
theta-only JFO problem with a fixed-pressure groove. Ground truth from
`exact_1d_jfo.py`: incompressible-JFO semi-analytic solution (quadrature to
machine precision) for the same h(theta), groove, p_cav.
Compared: rupture angle, reformation angle, p(theta) profile, load.
Sweep: n_theta = 120..1920, beta = 1e9 / 5e9.
Runs in `oneD_n<ntheta>_b<beta>/`.

### 2. Published 2D contour comparison (waiting on user data)

Ausas et al. 2007, ASME J. Tribol. 129:868 grooved journal bearing:
2D pressure field + cavitation boundary contours, mass-conserving algorithm.
Alternative if easier to obtain: Woloszynski, Podsiadlo & Stachowiak 2015,
Tribol. Lett. 58:18 (FBNS), same benchmark family.
Needed: digitized contour(s) and/or centerline profiles per cases/REFERENCE_DATA.md.

### 3. Per-cell mass-residual output field (code addition, DONE - passes)

Optional output field `mass_residual` (add it to `output_fields`) exposing the
converged per-cell continuity residual in kg/s, built from the same flux
expressions `solve_elrod` assembles. Definition in PHYSICS.md §15.1, test in
`tests/test_mass_residual.cpp`. Uncommitted; no commit without user approval.

Result on `oneD_n120_b1e9` (converged steady state), scaled by the peak Couette
face flux 4.4e-2 kg/s:

| region | max abs(r) / Couette flux |
|---|---|
| groove cells (theta = 81-99 deg) | 0 (penalty rows, source reported separately) |
| cell upstream of the groove edge (theta jumps 0.259 -> 1.00007) | 6e-17 |
| rupture front / reformation front | 8e-11 |
| peak-pressure zone (theta = 175-185 deg) | 1.1e-7 |
| whole domain | 1.1e-7 |

The peak-pressure blob is 3.2e-11 of the local row conductance
`rho_0 beta h^3/(12 mu) A/d`, i.e. the Krylov residual profile at the default
`linear_rtol = 1e-8`, and it is smooth rather than localized on any interface.

**Verdict: PASSES the #3 criterion.** The residual is at machine/solver precision
everywhere except the deliberately modeled groove source. In particular the
groove edges and the free-boundary cells conserve mass to round-off. The
grid-divergent rupture front is therefore not caused by a non-conservative flux
or boundary treatment; the remaining candidates are scheme accuracy and the
free-boundary (switch-function) treatment, which test #1 addresses.

Not yet exercised by a run: the z-boundary rows, because the 00c 1D cases use
NEUMANN axial ends. `tests/test_mass_residual.cpp` test 2 covers Dirichlet ends
on a 2D case and finds them as tight as the interior (2e-14 of the row scale).

## Verdict criteria (user rules)

- #1: pancake rupture/reformation angles -> exact values as grid refines, at rates
  consistent with the scheme order; p(theta) L2 error -> 0. Any O(1) offset that
  survives refinement = solver defect.
- #2: contour-level agreement within digitization accuracy.
- #3: residual localized at machine precision everywhere except deliberately
  modeled sources (groove cells, z-boundary rows).
