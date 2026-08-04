# Validation campaign plan (fresh start, 2026-08-04)

Goal: validate pancake for pressure distribution, load capacity, eccentricity under
external load (transient till convergence), THD, and gaseous cavitation against
high-quality reference data. The user judges pass/fail per case; agents run batches
(Opus 5, max 4 MPI ranks, sequential cases). This supersedes the root
`VALIDATION_PLAN.md` (previous attempt, kept for history only).

## Rules

- No batch runs on unverified parameters: every case gets a parameter table with
  verified/unverified flags and user sign-off before dispatch.
- All case inputs saved under `cases/<case>/`; results and `run.log` stay in place.
- User-digitized reference curves go in `cases/ref_data/` per `cases/REFERENCE_DATA.md`.
- Progress visible in `cases/<case>/progress.log`.

## Case ladder

| Case | Validates | Reference | Status |
|---|---|---|---|
| 00 grid/dt independence | discretization only (no physics comparison) | GCI, Celik et al. 2008 | RUNNING (params frozen: 4000 rpm, 6 kN cond., c = 152 um hot, 10 kg test mass for dt) |
| 01 Ferron benchmark | pressure p(theta) at midplane, eccentricity, attitude, THD temperature | Ferron, Frene & Boncompain 1983 (Fig. 6 = 4000 rpm, 6 kN, user-digitized) | waiting on 00 verdict + ref data |
| 02 static locus | load capacity vs eccentricity (Sommerfeld sweep) | Someya Databook; Raimondi & Boyd 1958 (L/D = 1, 0.5) | waiting on 00 |
| 03 transient to equilibrium | trajectory + converged endpoint under load step | endpoint vs 01/02; no classic transient experiment exists | after 01/02 |
| 04 THD secondary (optional) | temperature field on a second dataset | Boncompain et al. 1986 / Dowson 1966-67 | on request |
| 05 gaseous cavitation | cavity extent + pressure in cavitated zone | Etsion & Ludwig 1982; Braun & Hendricks 1984 | after 01 |

## Case 00 detail

`cases/00_grid_dt_independence/README.md`. 00a: 5 grids (60x16 to 480x128) steady,
fixed eps = 0.6, isothermal JFO + 2 THD confirmation grids. 00b: 5 timesteps
(4e-5 to 2.5e-6 s) transient 6 kN load step. Guidance: <1% change between
consecutive levels; user decides.

Beta finding (2026-08-04, runs grid_G{3,6}_beta5e9 and grid_G3_typediff): the grid
drift in p_max/attitude is beta-entangled (softer liquid -> front resolved
progressively), convection scheme is not a factor. USER DECISION: 00 stays at
beta = 1e9 (current value); every validation case must use the bulk modulus of the
paper's fluid (pin it during that case's parameter sign-off; if the paper does not
state one, propose a literature value for that oil and get sign-off).

## Frozen Ferron parameters (user-verified 2026-08-04)

R = 50 mm, L = 80 mm, c_hot = 152 um (145 cold + thermal expansion), groove 18 deg
wide at +90 deg (y-axis), load downward, supply 70 kPa gauge at 40 C,
mu(40 C) = 0.0277 Pa s, rho = 860, cp = 2000, k = 0.13, thermoviscosity 0.034 1/K.
Fig. 6 condition: 4000 rpm, 6 kN.

## Decision gates

1. 00 pass -> freeze grid + dt for all cases.
2. 01 pass -> pressure/eccentricity/THD credible; proceed 02/03 (cheap) and 05.
3. Each case: agent writes SUMMARY.md + comparison plots; user rules pass/fail;
   failures get a diagnosis case before any code fix is attempted.
