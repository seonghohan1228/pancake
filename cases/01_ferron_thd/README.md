# Case 01 THD: Ferron 4000 rpm / 6 kN with energy equation (STAGED - NOT RUN)

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Blocked on one parameter. Everything else is ready: FBNS phase 2 (commit 94edb45)
supports ENERGY_EQUATION; the groove MUST be thermally fed
(`inlet_groove = 90 18 1.713e5 313.15`) - a thermally transparent groove causes
unbounded thermal runaway under any cavitation model (see PHYSICS.md note).

## Parameter table (sign-off required before dispatch)

| parameter | value | status |
|---|---|---|
| geometry, speed, load, groove, supply p/T | as Case 00/01-prelim frozen set | VERIFIED (user, 08-04) |
| mu at 40 C | 0.0277 Pa s | VERIFIED |
| viscosity_temperature_coeff (Andrade b, in K) | ??? | **UNVERIFIED - BLOCKS RUN** |
| rho_cp = 1.72e6, k = 0.13 | Ferron oil | verified 08-04 with density/cp/k |
| journal/bearing wall treatment (HTC values or measured wall T) | ??? | UNVERIFIED - affects T comparison directly |

Needed from the user (from the Ferron paper):
1. The oil's viscosity at a second temperature (or the paper's stated
   thermoviscosity law) so b is fitted, not guessed. The linear coefficient
   0.034 1/K verified earlier is NOT the Andrade b (solver expects Kelvin units;
   b ~ mu drops by factor exp(b*(1/T1 - 1/T2))). From 0.034 1/K at 313 K,
   b ~ 0.034*313^2 ~ 3300 K is a consistent conversion - confirm or supply data.
2. Wall thermal boundary info: Ferron measured bush temperatures; the solver
   models walls via journal/bearing HTC or fixed wall temperatures. Choice +
   values need sign-off (candidates: measured shaft T from the paper, or
   conjugate treatment via cross_film model - Elrod-era module, check FBNS compat).

## Protocol once unblocked

Force-balance search (same driver as 01-prelim) at G6, FBNS + ENERGY_EQUATION +
SINGLE_PHASE + EMPIRICAL viscosity, fed groove at 313.15 K. Deliverables:
- p(theta) overlay vs ref_data/ferron1983_fig6_4000rpm_6kN_pressure.csv
  (test: does THD fix the early falling-flank rupture from 01-prelim?)
- T(theta) at midplane vs ref_data/ferron1983_fig6_4000rpm_6kN_temperature.csv
  (note: experiment is 0.5 mm sub-surface bush T - model equivalent must be
  stated, not silently equated to film T)
- equilibrium eps and attitude vs paper values.
