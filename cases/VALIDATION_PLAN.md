# Validation campaign plan (fresh start, 2026-08-04)

> **Where everything is: `cases/INDEX.md`.** That file is the front door - one
> table of every case directory, its stage, its status, its figures and whether
> the bulk data is local, on afml-r6 or both. This file is the ladder and the
> decision gates; the index is the map. Current scientific claims and their
> limits are controlled by `cases/VALIDATION_STATUS.md`.

Goal: validate pancake for pressure distribution, load capacity, eccentricity under
external load (transient till convergence), THD, and gaseous cavitation against
high-quality reference data. The user judges pass/fail per case; agents run batches
(Opus 5, max 4 MPI ranks, sequential cases). This supersedes the root
`VALIDATION_PLAN.md` (previous attempt, kept for history only).

## Rules

- No batch runs on unverified parameters: every case gets a parameter table with
  verified/unverified flags and user sign-off before dispatch.
- All case inputs saved under `cases/<case>/`; results and `run.log` stay in place.
  One case = one directory with `README.md` (what and why, written before the run)
  and `SUMMARY.md` (what came out). Every solver call archived as its own
  subdirectory with the `config.txt` that produced it. Nothing is ever deleted or
  overwritten - superseded output is renamed, never removed.
- User-digitized reference curves go in `cases/ref_data/` per `cases/REFERENCE_DATA.md`.
- Progress visible in `cases/<case>/progress.log`.
- DEFAULT EXECUTION HOST (user, 08-05): afml-r6 workstation (user@143.248.99.6),
  24 pinned ranks (OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1). Local Windows only
  for smoke tests, post-processing, and batches already in flight.
  Sync verification mandatory both ways: remote HEAD hash must equal the intended
  local commit (bundle transfer; no shared remote), untracked configs checksummed;
  results tar'd, sha256-verified across the transfer, extracted to results_remote/
  with one metric re-derived locally as an end-to-end check. Cross-platform gate
  on first use of each branch: one steady case must match Windows to 0.01%.
- Canonical workflow: develop local -> sync code to remote (bundle, hash-verified)
  + compile -> run remote -> copy results back (checksummed). cases/ additionally
  mirrored to `afml-r6:~/pancake/cases` via `cases/sync_to_r6.sh` after each batch
  (manifest-based, append-only, idempotent; it also reports remote-only files,
  which is how remote run output gets noticed).
- Solver branches: codex/fbns-cavitation (FBNS + THD) -> codex/full-film-predicate
  -> codex/fbns-gas (gas/outgassing). Stacked, tested, unmerged pending user.

## Directory layout (consolidated 2026-08-05)

Every case, probe and diagnostic run for this campaign lives under `cases/`.
Nothing is kept in a home directory, in `/tmp` or in an agent scratch tree; the
strays that were (`afml-r6:~/cases_run`, `afml-r6:/tmp/ksp_repro`, three loose
remote run directories) were folded in, and `cases/INDEX.md` records where each
one came from. Directory prefixes are the campaign stage and match the case
numbers in the ladder below - `00*` discretisation and solver verification,
`01*` Ferron family, `02*` static locus, `90*` diagnostics. No directory was
renamed during the consolidation; the reasoning is stated in the index.

Stage `90_diagnostics/` is new. It holds work that explains or defends solver
behaviour rather than physics, so a validation case never has to carry a solver
postmortem: `ksp_rank_robustness/` (the rank-dependent `DIVERGED_DTOL` that cost
the conjugate campaign a variant, closed by `f0e772d`) and `remote_harness/`
(reference copies of the loose afml-r6 driver scripts). Nothing in stage 90 gets
a pass/fail verdict.

## Case ladder

| Case | Validates | Reference | Status |
|---|---|---|---|
| 00 grid/dt independence | discretization only (no physics comparison) | GCI, Celik et al. 2008 | **PASSED by the user as an operational campaign resolution (08-05): G6 720x192, dt = 2e-5 s.** Scope limit: load/torque and the tested transient endpoint are converged; the original sweep did not show grid independence of cell-centred `p_max`, and its THD pair did not exercise the later Ferron thermoviscous/conjugate model. G6/dt is therefore a frozen campaign protocol, not a universal discretization-error certificate. **Reinforced 08-11 by case 01e's independent ladders, which do exercise the Ferron thermoviscous/conjugate model:** the fresh D-grid ladder runs G1--G8 and reports cell-centred `p_max` converging to 0.007 % across G7--G8, closing the original `p_max` gap for that model; the D-dt ladder extends the transient endpoint from 0.06 s to 0.10 s with 20/10/5 microseconds agreeing to 0.0012 % in final eccentricity. Both confirm the frozen G6 / 2e-5 s selection rather than replacing it. |
| 00c JFO verification | solver correctness | exact 1D JFO benchmark; mass_residual field | DONE: Elrod front error first-order/~7-cell prefactor -> replaced by FBNS (fd2b0ea): front <=1 cell, no beta. Phase 2 (94edb45) unblocks THD on FBNS (validated: 0.44 K vs Elrod on shared arc, W -12% thermal, sane). Gas on FBNS blocked on is_full_film refactor (docs/FBNS_GAS_COUPLING.md); refactor dispatched 08-05. Finding: thermally-transparent grooves cause thermal runaway under ANY model - THD grooves must be fed (t_supply) |
| 01 Ferron benchmark | pressure p(theta) at midplane, eccentricity, attitude, THD temperature | Ferron 1983, Figs. 4/5/6 (user-digitized in ref_data/) | PRELIM + THD(2000rpm/4kN) done: pressure peak +0.8%/1.25 deg (THD), eps(W) parallel within clearance ambiguity. THD temperature: CHT (2-D conjugate bush, outer h=10 to 40 C, shaft adiabatic) is verdict-ready - peak location 265 vs 258 deg, mean +0.02 K, rms 1.17 K vs Fig. 4 (01_ferron_wallbc_fields/SUMMARY.md); peak value -2.3 K remains. User verdict pending. ReynSim comparison in 01_ferron_vs_reynsim: This year inside band, Last year -17.6% peak |
| 02 static locus | load capacity vs eccentricity (Sommerfeld sweep) | R&B (textbook transcription found in validation/journal_bearing, scored 08-05: S within 1.8-2.5%, phi within 0.54 deg at L/D=1 and 0.5 - user verdict pending); Ferron Fig. 5 + Brito via S-eps collapse (`02_static_locus/COLLAPSE_ANALYSIS.md`): Ferron's discrepancy has no load-dependent structure over 5x load, but its absolute level remains exactly degenerate between a uniform clearance bias and a uniform load-model bias. Brito retains a +9.3% speed-dependent effective-clearance trend that no constant rescale removes; thermal clearance growth and under-predicted temperature remain hypotheses, not a unique diagnosis. One absolute film-thickness measurement is needed to break the level degeneracy. | R&B Part III / Someya still welcome for the L/D=1/4 block |
| 03 transient to equilibrium | trajectory + converged endpoint under load step | endpoint vs 01/02; no classic transient experiment exists | after 01/02 |
| 04 THD secondary (optional) | temperature field on a second dataset | Boncompain et al. 1986 / Dowson 1966-67 | on request |
| 05 gaseous cavitation | cavity extent + pressure in cavitated zone | Etsion & Ludwig 1982; Braun & Hendricks 1984 | **05a experimental extent comparison RUN 2026-08-05; user verdict pending. It does not validate the gaseous-cavitation model.** On all 16 scored runs `alpha_gas = 0`, so the comparison exercises single-phase FBNS/JFO under an air-release pressure datum. FBNS rupture mean \|error\| is 5.9 deg and the 1840-rpm ambient-pressure trend is captured, but reformation is 25-35 deg late at 3000 rpm. JFO fails. Spatial cavity pressure and finite-rate gas release remain unvalidated; 05b is blocked on a new pressure/desorption formulation. |

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
