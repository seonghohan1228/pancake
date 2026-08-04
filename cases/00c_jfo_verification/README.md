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
Runs in `oneD_n<ntheta>_b<beta>/`. The same configuration with
`cavitation_model = FBNS` runs in `oneD_n<ntheta>_fbns/` (bulk_modulus is
irrelevant there). `extract_front.py <case_dir> ...` tabulates both.

**Result (2026-08-04): the criterion FAILS for `JFO` and PASSES for `FBNS`.**
Front reported at the face between the last full-film and the first cavitated
cell of the mid-z row, same estimator for both solvers:

| n_theta | cell | JFO rupture err | cells | FBNS rupture err | cells | FBNS reform err | cells |
|---|---|---|---|---|---|---|---|
| 120 | 3.000 deg | -23.77 deg | -7.9 | -2.767 deg | -0.92 | +1.99 deg | +0.66 |
| 240 | 1.500 deg | -17.77 deg | -11.8 | -1.267 deg | -0.84 | +0.49 deg | +0.33 |
| 480 | 0.750 deg | -10.27 deg | -13.7 | -0.517 deg | -0.69 | +0.49 deg | +0.66 |
| 960 | 0.375 deg | -5.02 deg | -13.4 | -0.142 deg | -0.38 | +0.12 deg | +0.32 |

The JFO error shrinks in degrees but **grows** in cell widths: it is a
first-order error with a prefactor of roughly 7 cells, i.e. the front is not
converging with the mesh. Integrated quantities follow: JFO p_max is +50.4% /
+39.9% / +24.0% / +12.2% above the exact 7.597141 MPa and W' is +42.0% / +41.0% /
+31.1% / +17.8% above 495463.6 N/m, while FBNS is +0.54% / +0.28% / +0.14% /
+0.07% and +0.70% / +0.27% / +0.07% / -0.04%.

Diagnosis: the Elrod switch reads the free boundary through `p - p_cav`, which is
`O(beta^-1)` small there because `p = p_cav + beta ln(theta)`. Test #3 had
already ruled out a non-conservative flux or boundary treatment, so the free
boundary treatment itself was the remaining candidate — and replacing it with the
complementarity form fixes it without touching any other part of the operator.
FBNS converges in 4 Newton iterations at every resolution.

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

### 4. FBNS per-cell mass residual (code addition, DONE - passes)

`mass_residual` under `cavitation_model = FBNS` mirrors the FBNS operator
(`FBNS::cell_outflow`) rather than `solve_elrod`'s, because the face fluxes
differ: pressure-driven Poiseuille with no bulk modulus and no switch function,
always-upwind Couette, and a plain half-cell axial Dirichlet link. On the
converged `oneD_n480_fbns` mid row (peak Couette face flux 4.38e-2 kg/s):

| region | max abs(r) | / Couette flux |
|---|---|---|
| groove cells | 0 (source reported separately) | 0 |
| rupture / reformation front cells | 1.5e-15 kg/s | 3.3e-14 |
| every other cell | 1.5e-11 kg/s | 3.5e-10 |

The whole-domain worst case is the Krylov tolerance (`linear_rtol = 1e-8`), and
the free-boundary cells - the ones the JFO front error lives on - are at machine
precision. No clamp accounting is needed: the FBNS truncations act on Newton
iterates, not on the converged state.

## Verdict criteria (user rules)

- #1: pancake rupture/reformation angles -> exact values as grid refines, at rates
  consistent with the scheme order; p(theta) L2 error -> 0. Any O(1) offset that
  survives refinement = solver defect.
- #2: contour-level agreement within digitization accuracy.
- #3: residual localized at machine precision everywhere except deliberately
  modeled sources (groove cells, z-boundary rows).
