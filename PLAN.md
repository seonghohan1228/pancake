# Development Plan

## 2026-09-15 - Solver audit and measured repairs

- [x] Trace the property, inlet, energy, linear-solver and MPI paths; consolidate
  inlet geometry and prescribed-value row handling without removing models.
- [x] Reuse final coupled states and eliminate duplicate EOS bracket evaluations.
- [x] Honor runtime stopping controls, preserve failed-iteration counts, reject
  undersized ghost exchanges and stop after PETSc initialization failure.
- [x] Correct incoming liquid pressure work under temperature-dependent density.
- [x] Complete and record 65 native and 65 workstation checks, matched output
  comparisons and measured FW68D runtime reductions at 1/4 MPI ranks.
- [ ] Discuss retention/labeling of approximate legacy gas closures with the user.
- [ ] Decide the validity policy for linear density with coupled caloric properties.
- [ ] Select a production workload and accuracy target before replacing the
  numerical Jacobian or implementing higher-order conservative transport.

See [the audit](docs/SOLVER_AUDIT_2026_09_15.md). Completed earlier work below
remains part of the development record.

## 2026-09-11 - Oil operating-window screening

- [x] Add an isolated, reproducible steady force-balanced load/rpm scan using
  the requested oil viscosity and density, 1 bar end/groove pressure, and the
  existing remote FBNS executable. Preserve unrelated local and remote work.
- [x] Cover 1--5 kN and 250--2000 rpm, record a lower-speed probe, export pressure
  and liquid-filling contours, and audit selected mesh refinements in case 02b.
- [x] Correct the geometry on user instruction to D = 43 mm, L = 40 mm and
  radial clearance = 42.5 micrometres; rerun in case 02c and mark 02b superseded.
- [x] Add requested temperature maps using the energy equation and two explicit
  provisional cooling scenarios, with constant viscosity and checked force closure.
- [ ] Replace provisional inlet/wall temperatures, cooling, heat capacity and
  conductivity with rig/oil data before interpreting the thermal maps as predictions.
- [ ] Confirm groove location/width and oil vapor
  pressure before treating this screening as an experiment
  prediction. Assess deformation and viscosity change at high-load points.

## 2026-09-08 model-selection discussion

- [x] Explain physical versus numerical load fluctuations, coupled energy
  changes, and model assumptions against current source.
- [x] Map pressure, phase, property, thermal, boundary and motion options to
  appropriate validation observables in [MODEL_SELECTION.md](docs/MODEL_SELECTION.md).
- [x] Re-read Ferron's full paper; correct reduced-solid-model equivalence
  claims and distinguish archived branch results from the new coupled path.
- [x] Record the required local-to-remote simulation workflow in AGENTS.md.
- [ ] Assess Ferron fresh mixing flow against pressure-derived inlet mass
  supply with independently specified hydraulic and thermal conditions.
- [ ] Quantify sensor-depth and reduced shaft/bush thermal-model errors before
  claiming Ferron temperature validation. Do not equate film-average and wall
  temperatures or apply a synthetic fixed-eccentricity load change to Ferron.
- [ ] Assess resolved transient statistics and numerical damping if physical
  unsteadiness is a validation target; smooth load alone is not evidence.

## 2026-09-07 coupled cavitation and heating implementation

- [x] Solve pressure, conserved oil/dissolved/free inventories and gas volume
  simultaneously; retain the gas-free JFO limit and exact accepted face fluxes.
- [x] Add native fixed-blend R410A density and caloric properties; check against
  the published Helmholtz EOS through independent CoolProp reference calls.
- [x] Predict film heating with conservative energy, compression work,
  viscous dissipation and dissolved/free energy differences. Verify release,
  compression, wall cooling and reservoir enthalpy independently.
- [x] Support separate axial and internal reservoir compositions and directional
  heat/species transport for groove, hole and absent internal supply.
- [x] Add concentration/temperature viscosity tables with strict bounds and
  a labeled FVC68B reference fixture. Do not relabel it as FW68D.
- [x] Pass 51 CTest checks on native Windows and local WSL Linux; build the GUI.
- [x] Complete paired oil-air, thermal and gas-cavity runs on the remote
  workstation after the face-gap transport and outflow-property repairs.
  Preserve failed R410A runs and distinguish numerical/domain failures.
- [x] Repair the outgoing R410A Newton-trial domain gate exposed at 20°C.
  Verify the formerly failing dissolved-solution case and complete the
  four-temperature 20-bar equilibrium matrix without weakening accepted-state
  phase checks. Retain the finite-rate condensation-domain failures.
- [x] Establish spatial, temporal and MPI refinement evidence for the
  synthetic acceptance cases. The workstation pilot matches 1/4/8-rank
  loads with a 3.70× speedup on eight ranks. The final thermal mesh doubling
  changes load by 0.234%; the finest independent JFO pressure error is 0.0432%.
- [ ] Extend accuracy/runtime evidence to experimentally characterized FW68D
  geometries and properties; these case-specific results are not universal.
- [ ] Resolve the Etsion gas-cavity mismatch; retain failed comparisons as
  evidence. A converged homogeneous model is not automatically validated.
- [ ] Obtain FW68D mixture solubility, actual-concentration viscosity, partial
  volume and caloric data, charge temperature and thermal boundary data.
- [ ] Assess separate refrigerant-rich liquid, blend fractionation, oil-vapor
  latent energy and cross-film/conjugate heating against the intended tests.

The older dated plan entries below are preserved as development history;
this section and the linked current physics documents control present scope.

## 2026-09-07 clarification: FW68D, pressure reservoirs and phase identity

- [x] Record FW68D, 20–80 °C, reservoirs to 20 bar absolute, R410A charge near
  10 bar, and the independent oil-air requirement in the
  [detailed review](docs/FW68D_CAVITATION_REVIEW.md).
- [x] Verify that pressure ends admit reverse flow and that groove, circular
  hole and no-internal-inlet geometries have distinct meanings. Add per-end
  inflow state and conservative bidirectional species/thermal tests.
- [x] Reject silent boundary-pressure replacement under `LOCAL_PSAT`.
- [x] Separate absorption rate from release and verify compression without
  fictitious gas loss. Do not transfer paraffin kinetics to FW68D.
- [x] Compare published mechanisms and recorded numerical results; correct
  the Etsion uncertainty-count and gas-validation overstatements. Quantify the
  pure-R410A ideal-gas error using a traceable EOS diagnostic.
- [ ] Obtain charge temperature and measured/sourced FW68D-air and
  FW68D-R410A material data. Do not substitute R290/PZ68S or R410A/FVS50.
- [x] Implement and verify coupled pressure/species/EOS/volume closure; use the
  1-ms air-bearing timestep failure below as the comparison case.
- [x] Implement R410A fixed-blend thermodynamics with explicit phase-domain
  limits and caloric properties; mixture experimental assessment remains open.
- [x] Generalize internal reservoir species treatment for reverse flow and
  independently characterized groove/hole reservoirs.
- [ ] Reproduce gaseous and vaporous bearing benchmarks with the completed
  formulation, then validate FW68D experiments and measure runtime at fixed
  accuracy. A property diagnostic is not a bearing validation.

## 2026-09-06 audit: oil-air and oil-R410A target

The [physics audit](docs/PHYSICS_AUDIT_2026-09.md) is the current assessment of
this target. The earlier roadmap and completed items are retained below.

- [x] Review accessible branches, reachable history and recorded validation;
  isolate work from the dirty Ferron campaign branches.
- [x] Repair old/new liquid storage, conservative gas transport, phase-exchange
  mass fractions, reservoir sources and actual gas boundary-flux accounting.
- [x] Remove transient gas relaxation and kinetic overshoot; correct steady
  flash mass fractions without claiming a steady transport solution.
- [x] Reject unsupported species and unphysical dissolved partial volume;
  unify mean/cross-film table viscosity and expose gas convergence/volume errors.
- [x] Fail the CLI on nonconvergence; preserve unrelated output-directory files.
- [x] Verify 40 tests on native Windows and Linux, exact 1-D FBNS fronts and
  integrated quantities, and a synthetic oil-air transient across MPI counts.
- [x] Compare pre-audit and repaired solvers on frozen steady and transient
  inputs, with direct inventory/time checks and publication-style figures.
  [Evidence](cases/90_diagnostics/physics_comparison/SUMMARY.md): steady curves
  unchanged; old closed reaction loses 11.548% of its gas. A 1 ms air transient
  fails the 10/5/2.5 us timestep check for pressure/load despite converged inner
  solves. Gas fraction alone is not a sufficient acceptance metric.
- [x] Close pressure, phase inventories, EOS and gap volume simultaneously.
  Require both conservation and volume residuals; no post-step mass-destroying
  theta/alpha clipping. Keep FBNS as the verified zero-gas limit. Use the
  archived air timestep failure as an acceptance case; require removal of the
  timestep-linked load alternation without averaging or lost species mass.
- [x] Fix transient timestamps, final partial intervals and moving-geometry time
  levels. Verify initial output, integrated heating and rigid-body acceleration
  against exact CLI reference cases, including zero-duration initial output.
- [ ] Rerun transient campaigns with corrected time levels and establish
  fluid-force/thermal/species splitting accuracy before fitting time constants.
- [ ] Obtain the oil grade, operating ranges and initial dissolved/entrained gas
  state. Validate oil-air pressure, flow, load and void history independently.
- [ ] Add oil-specific R32/R125 solubility/partitioning, liquid properties and
  vapor EOS for R410A. Justify any reduced blend approximation with data.
- [ ] Separate vapor and noncondensable inventories where both are present;
  include latent/dissolution enthalpy when it affects the thermal balance.
- [x] Complete discrete-flux consistency for energy, free-gas velocity and JFO
  transport; verify steady species transport rather than a local feed flash.
- [ ] Reconcile later Ferron thermal branches on their own task; preserve their
  measured-observable, boundary-condition and validation-eligibility caveats.
- [ ] Profile representative verified cases before optimizing KSP, assembly,
  allocation or MPI layout. The present small-case timings are not scaling data.

## Roadmap — refined 2026-08-05 (post-campaign consolidation)

Inputs: `cases/INDEX.md` (campaign map), `docs/UNFINISHED_BUSINESS.md`
(validation gaps, prioritized), `docs/CODE_AUDIT_2026-08.md` (code/GUI audit).
Current claims are controlled by `cases/VALIDATION_STATUS.md`; the two audit files
are dated inventories and may describe work that has since closed.
Old-plan reconciliation: the 2026-08-02 campaign's Phase 1 (reference-data
verification) closed via user-verified parameters (c_hot = 152 um et al.);
Phase 3 (grid/dt independence) PASSED but on FBNS, which made its beta-scaling
caveat obsolete (Elrod's beta-entangled front was root-caused and replaced,
00c); Phase 4 (validation matrix) is partially hit: Ferron pressure,
eccentricity (parallel within clearance ambiguity), THD temperature
(the formerly favorable CHT score was invalidated by inlet-source accounting;
corrected G6 is +3.82 K mean / 4.02 K RMS). The WP-1b
VOID_COUPLED convergence items and the five deferred modelling gaps
(2026-07-31) remain open and are now tracked in UNFINISHED_BUSINESS.

Priority queue (user gates in brackets):

0a. **Ferron Figure 6 measured-flow campaign** - **IN PROGRESS 2026-08-09**.
    At 4000 rpm / 6000 N, start the force-balance search from `eps=0.58` and
    hold the user-digitised Figure 10 supply fixed at `Qi=7.837 L/min`. Run four
    G6 thermal closures: depth-averaged, adiabatic cross-film, both-wall 40 C
    cross-film, and fixed-outer-40 C CHT. Do not rerun `h=80`; it represents
    deferred axial-face cooling. Return the Ubuntu archive, compact profiles,
    Ferron overlay, numerical verdicts, and exact remote provenance.

0. **Ferron CHT grid sensitivity** - **DONE, FAIL 2026-08-07**. The first pass
   exposed and fixed a groove-source diagnostic that grew as `1/dz`. The clean
   hash-verified rerun passes five of seven G6-to-G7 checks but misses the 0.10 K
   mean/RMS limits at 0.111/0.107 K. Corrected G6 is +3.82 K mean / 4.02 K RMS;
   the previous +0.02/1.17 K validation score is superseded. Evidence: case 00d.

1. ~~**Picard cap correctness**~~ - **DONE 2026-08-05**, see "In Progress" below.
   Return-best + non-convergence flag landed; under-relaxation measured
   (omega = 1 never converges, 0.95 converges in 48 sweeps, ADAPTIVE in 52);
   within that historical closure, the CHT reruns move no reported number beyond
   printed digits. Case 00d later superseded its physical validation score after
   correcting the inferred supply flow.
2. **Audit quick wins** - unknown-config-key hard error (the campaign shipped
   a live example of the silent-drop failure); VANLEER reclassified live and
   fixed/fenced; GUI anisotropic-stiffness load bug; config.txt default ->
   FBNS; energy/gas coefficient unification onto ReynoldsCoeff + A/B rerun.
3. **Load-capacity adjudication** - the S-eps collapse removes load-dependent
   structure from Ferron but cannot separate a uniform clearance bias from a
   uniform load-model bias. Do not label the load calculation exonerated without
   an absolute film-thickness or independently calibrated-clearance dataset.
4. **Gaseous-cavitation formulation** - Etsion extent case 05a has run, but all
   scored cases have `alpha_gas = 0`; it therefore does not validate gas release.
   Cavity pressure and finite-rate desorption need a new formulation before 05b.
   Supplements:
   Ochiai 2019 (meshable geometry), Reinke 2022 (gas-release rate constants),
   Engels 2024 (void fraction). [user picks scope]
5. **Case 03 transient-to-equilibrium** - nearly free (retarget dt_T3 config,
   2/4/6 kN), cross-checks moving-bearing vs static force balance.
6. **Ferron temperature residual** - resolve the fresh-oil groove-mixing flow
   against Ferron's measured inlet-flow evidence before STRIATED, shaft-Robin,
   or other thermal-model sensitivity work. The 2026-08-07 boundary audit also
   found that the current z-collapsed bush omits Ferron's axial-face \(h_b\)
   condition and the scalar shaft omits its r-z/43 C end condition. Independent
   outer-bush and shaft sink temperatures now permit bounded sensitivity runs;
   a full Ferron reproduction still requires r-theta-z bush and r-z shaft solids.
   Case `01e_ferron_terra_repro` is the isolated remote rerun in progress; its
   reduced-model numerical/grid verdicts must remain separate from any Ferron
   temperature-validation verdict. Its D add-on separately reruns the present
   CHT grid ladder and a 4 kN, 2000 rpm physical timestep ladder; neither is a
   substitute for the missing full-solid and independent-flow evidence. A
   repaired-legacy mode reproduces Claude's original GCT only for comparison;
   it is report-only because the nominal THD viscosity input is effectively inert.
   The fresh steady D-grid ladder covers G1--G8 and reports both G6--G7 and
   G7--G8 checks; its figures are generated locally from
   transferred field data. The physical timestep ladder also passes its
   registered 20/10/5 microsecond comparisons, selecting 20 microseconds; the
   80/160/320 microsecond stability probes and force-imbalanced 40 microsecond
   endpoint stay visible in the report rather than being counted as successful
   cases. That ladder now runs to 0.10 s of physical time rather than 0.06 s,
   and the 20/10/5 microsecond rungs agree there to 0.0012 percent in final
   eccentricity and to machine precision in load, so the 20 microsecond
   selection holds over the longer settling window.
7. **Integration decision** - local `main` contains the 45-commit Claude snapshot
   while `origin/main` is behind it. Reconcile and verify on a task branch before
   any merge or push; do not continue development directly on `main`.
8. **Structural** (audit S-items): GUI exposure of FBNS/cross-film/bush/KSP
   keys; test gaps (transient THD, MOVING_BEARING+FBNS, gas+cross-film,
   CONDUCTION_2D+FBNS); FBNS Schur conditioning near complementarity switch;
   VOID_COUPLED-under-FBNS p_void update site.

## In Progress

### Audit quick-win batch (`codex/audit-quickwins`) — priority-queue item 2

Complete 2026-08-05. Evidence and the deliberate exclusions are in
`docs/reports/QUICKWINS_2026-08.md`; per-item rationale in `CHANGES.md`.

| Item | Audit ref | Status |
|---|---|---|
| Unknown-config-key hard error, with typo suggestion | Q1 | **DONE 2026-08-05** — `tests/test_config_keys.cpp`; 105 case configs corrected |
| VANLEER far-upwind cell for negative face flux | C / F4 | **DONE 2026-08-05** (theta faces; z-face TVD still open, F4 stays partial) |
| GUI anisotropic-stiffness data loss on load | 3.5 / S9 | **DONE 2026-08-05** - normalize on new-case and on support edits only; load/save preserve |
| `config.txt` cavitation default -> FBNS | Q5 | **DONE 2026-08-05** - struct default stays JFO until the 11 default-reliant tests pin their model |
| Coefficient unification: energy axial BCs + gas transport | A / S1 | **DONE 2026-08-05** - A/B on the archived Ferron cross-film point: load/pressure/temperature move at the 6-digit output floor, `cross_film_core_rise` and `q_journal` by ~5e-5 relative |
### Picard cap correctness — DONE 2026-08-05

The STEADY_STATE outer loop returned whatever iterate `coupling_max_iters` landed
on and never said it had failed to meet `coupling_tolerance`. On the conjugate
Ferron cases that iterate was 6.3x worse than one the loop had already passed
through (`max(R_p,R_T)` bottoms out at 4.3e-5 near outer 56, then oscillates
between 5e-5 and 3.6e-4 to the cap). Landed:

- new `Coupling` module: best-iterate tracking, rollback **on cap exhaustion
  only** (a converging run is bit-identical to before), a whole-field snapshot
  taken before under-relaxation so the delivered state is a real end-of-sweep
  solve, loud `NOT CONVERGED` logging independent of `log_outer_iters`, and eight
  `coupling_*` columns appended to `diagnostics.csv`;
- `coupling_relaxation_mode = ADAPTIVE` (default off) plus
  `log_coupling_history` for the per-sweep residual trace.

Measured on afml-r6 at 24 pinned ranks: the loop gain is only marginally above 1
(omega = 0.95 already converges, in 48 sweeps against 200-and-never at 1.0), and
ADAPTIVE gets there in 52 with no tuned value. The archived answers survive: the
bush-side wall temperature profile moves <= 1e-3 K, W by 0.11 N in 4002 N,
p_max by 40 Pa in 1.405 MPa, T_max by 1 mK across the iterate-200 state, the
rolled-back best iterate and every converged relaxation setting. **The -2.3 K CHT
peak residual, mean +0.02 K and rms 1.17 K are unchanged within that closure**.
Case 00d later invalidated the inferred supply flow and therefore supersedes the
physical score, not the Picard comparison. Report:
`cases/01_ferron_wallbc_fields/PICARD_FIX_REPORT.md`.

Free by-product: `01_ferron_conjugate/outer_robin80` (h = 80), abandoned since
its `fb_iter_06` aborted on the KSP divergence, was re-run and **closes** -
eps = 0.602770, mean -0.20 K, rms 1.21 K over the 16 Fig. 4 points, within
0.26 K of the h = 10 production group. The h = 80 hole in the conjugate case is
shut, and the outer-coefficient sensitivity turns out to be small next to the
force-balance closure.

Remaining, deliberately not done here: HANDOFF **D8** is untouched - the
relaxation still acts on p and T while the Elrod primary state is theta; and the
best iterate on the conjugate case is still 43x above `coupling_tolerance`, so
the rollback is an honesty measure and under-relaxation is the fix. Whether the
frozen case protocol should adopt omega < 1 or ADAPTIVE is a user decision,
because it changes a protocol every archived conjugate number was produced under
(the numbers themselves do not move).

### Rank-dependent linear-solve divergence — DONE 2026-08-05

Solves diverged as a function of MPI rank count: the Ferron THD force-balance
iterate aborted at 4 ranks and converged at 1, 2 and 6; `outer_robin80` lost
`fb_iter_06` on 4 local ranks and `fb_iter_03` on 24 remote ranks; `test_fbns`
aborted at `-n 4`, which ctest never ran. Cause was block-Jacobi preconditioning
of convection-dominated operators on a theta-only decomposition — every rank
boundary cuts the flow direction, and the steady energy and Elrod equations have
no diffusion to spare. Default preconditioner is now ASM with one cell of overlap
(`linear_pc_type`, `linear_pc_overlap`), block Jacobi is kept on a single rank
where ASM degenerates but takes a worse ILU path, and a failed solve is retried on
`linear_solver_fallbacks` escalation rungs before aborting. See CHANGES.md
2026-08-05 for the iteration-count, wall-time and solution-agreement tables.
ctest now covers `-n 1` and `-n 4`, not only `-n 2` and `-n 6`. Verified at 1, 2,
4 and 6 ranks locally and at 24 ranks on afml-r6; solutions match the pre-fix
binary to all printed digits everywhere they converged before, byte-identically
in serial.

Remaining: the FBNS Schur complement is fragile whenever a complementarity pivot
sits on its threshold — a last-bit difference flips the pivot and leaves a
near-singular row. The fallback ladder covers it, but a conditioning fix (row
scaling, or a regularised pivot near the switch) would remove the need.

### Validation campaign restart (2026-08-04) — see `cases/INDEX.md` first

**Front door: `cases/INDEX.md`** — one table of every case directory, its campaign
stage, its status, its figures, and whether the bulk data is local, on afml-r6 or
both. `cases/VALIDATION_PLAN.md` holds the ladder and the decision gates;
`cases/REFERENCE_DATA.md` holds the digitisation format and what is still owed by
the user.

Consolidated 2026-08-05: every case, probe and diagnostic run now lives under
`cases/`, numbered by stage (`00*` discretisation and solver verification, `01*`
Ferron family, `02*` static locus, `90*` diagnostics). Strays that used to live
on the remote host (`~/cases_run`, `/tmp/ksp_repro`, three loose run directories)
were folded in; the repo-root `validation/` tree is a separate, pre-campaign
analytical benchmark harness (WP-7) and stays where it is.

Fresh validation campaign run from scratch at the user's direction; the earlier
`VALIDATION_PLAN.md` (root) is retained for history but no longer tracks status.
Case 00 (grid/dt independence on the Ferron geometry) is running.

Case 00c (JFO solver verification), test #1 — exact 1D benchmark: **CLOSED by a
new solver.** The Elrod-Adams rupture front is 7–13 *cell widths* upstream of the
exact solution at n_theta = 120..960 and does not converge relative to the mesh,
which is a free-boundary-treatment defect, not a conservation defect (test #3
already ruled that out). `cavitation_model = FBNS` (`src/fbns.{hpp,cpp}`,
PHYSICS.md §7.3, `tests/test_fbns.cpp`) solves the same JFO conditions as a
complementarity problem — no bulk modulus, no switch function — and lands the
rupture front within 1 cell and the reformation front within 0.66 cells at every
resolution, with p_max/W' errors of 0.5%→0.07%, in 4 grid-independent Newton
iterations. Runs in `cases/00c_jfo_verification/oneD_n*_fbns/`. Test #2 (published
2D contours) still waits on data.

FBNS **phase 2 — DONE 2026-08-04** (energy equation + `SINGLE_PHASE` properties;
PHYSICS.md §7.3.6/§7.3.7, `tests/test_fbns.cpp` tests 8–10, 27/27 suite).
`temperature_model = ENERGY_EQUATION` and `fluid_property_model = SINGLE_PHASE`
are now supported, so Cases 01-THD and 04 can run on FBNS. Coupling is
frozen-coefficient Picard around Newton, the same split the Elrod path uses: the
face mobilities read the per-cell `mu` field, held fixed inside one Newton solve,
with the `mu(T,p)` feedback closed by the existing outer loop. The `SINGLE_PHASE`
validation quirk (zero Henry slope and zero log-mixing exponent rejected even with
no gas species present) is fixed. Elrod/Gumbel byte-identical; isothermal FBNS
byte-identical to phase 1. Acceptance run on the `grid_G3_fbns` geometry: 9 Picard
iterations, W 9.96 → 8.76 kN thermal.

Two things came out of it that are **not** FBNS-specific and are still open:

- **Groove thermal runaway (both cavitation models).** A groove that pins pressure
  but is thermally transparent (`feeds_fresh_oil = false`) diverges to `T ~ 1e6 K`
  under a strong Andrade slope, because the Poiseuille dissipation scales as
  `1/mu` while `grad p` is held fixed by the pressure pin. Reproduces identically
  under JFO. Workaround: give the groove a supply temperature. A fix belongs in
  `energy.cpp` on its own branch (it changes the Elrod path).
- **Gaseous cavitation on FBNS — memo, `docs/FBNS_GAS_COUPLING.md`.** The bound
  `theta <= 1 - alpha_gas` composes cleanly with Fischer-Burmeister (~20 lines in
  `fbns.cpp`); the blocker was that "full film" is spelled `theta >= 1.0` at 15
  sites across four files, all of which would switch off at once. Prerequisite: a
  shared `is_full_film` predicate landed as a byte-identical refactor first.

  Prerequisite **DONE 2026-08-05** (`codex/full-film-predicate`): `FullFilm` in
  `src/full_film.hpp`. One ceiling, `FullFilm::bound`, returning exactly `1.0`
  today; all 15 sites routed through `is_full_film` / `is_full_film_face` /
  `gate` / `face_gate`. Byte-identical across five regression cases (41 output
  files, 41 sha256 matches) and `ctest` 27/27.

  FBNS **phase 3 — DONE 2026-08-05** (`codex/fbns-gas`; PHYSICS.md §7.3.7,
  `docs/FBNS_GAS_COUPLING.md` §8, `tests/test_fbns.cpp` tests 11-13, `ctest`
  27/27). `cavitation_model = FBNS` + `fluid_property_model = TWO_PHASE` runs:
  the bound is `(1 - alpha_gas) - theta >= 0` with `alpha_gas` lagged, so both
  Jacobian blocks and the Schur structure are unchanged. Complementarity against
  the moved bound is satisfied to 0 exactly and the per-cell mass residual is
  `5.1e-15` of the peak Couette flux. Elrod/Gumbel byte-identical to 28ac88a
  (54/54 output files over six configs, including JFO+TWO_PHASE and
  JFO+VOID_COUPLED+LOCAL_PSAT). The `alpha_gas <-> theta` feedback reuses the
  existing `void_coupling_relaxation` lever. **Case 05 can now be run on FBNS as
  well as JFO** — worth doing as a cross-model check, since the two treat the gas
  volume differently by design (inert on Elrod, displacing on FBNS).

  Still open from the memo:

  - `gas_pressure_coupling = VOID_COUPLED` under FBNS stays rejected. The
    per-cell plateau `p_void` is advanced only inside `Reynolds::solve_elrod`'s
    outer loop; FBNS would run it frozen. Needs its own update site, its own
    relaxation schedule and its own validation (memo §8.3).
  - Merging `Reynolds`'s `VoidParams::theta_full` onto `FullFilm::bound`. Not
    done on purpose: it changes the validated JFO gaseous path (case 05) and
    would have broken the byte-identity that made phase 3 checkable (memo §8.4).
  - `config_r290_pz68_quick.txt` under FBNS is stiff — a timestep-independent
    load swing from the saturated feed flashing ~23% of the gap per step. Every
    Newton solve converges and the gas mass stays bounded, but the first tens of
    steps are transient noise. The Elrod `VOID_COUPLED` path does not resolve it
    either (memo §8.5).
  - Cross-film Dowson coefficients in the FBNS operator — **DONE 2026-08-05**
    (`codex/fbns-crossfilm`); see the entry below.

Cross-film generalized coefficients under FBNS — **DONE 2026-08-05**
(`codex/fbns-crossfilm`; PHYSICS.md §7.3.8/§7.3.9, §13.1.3, §15.1, CHANGES.md,
`ctest` 27/27). `cavitation_model = FBNS` now accepts
`cross_film_thermal_model = POLYNOMIAL`: the face mobilities are `rho_l F2` and
`2 F1/F0`, i.e. `mu(T(y))` acts inside the film integrals rather than through a
depth-averaged `mu`. `CrossFilm::update` rebuilds the F-integrals once per outer
Picard iteration, before `FBNS::solve`, so within one Newton solve they are
frozen coefficients exactly as `mu` already was — the Jacobian blocks, the
linearity of the mass balance and the quadratic Newton rate are unchanged.

Same branch, second half: the **diagnostics operator now mirrors the assembled
one** for both cavitation paths. `src/reynolds_coefficients.hpp` holds the single
definition of the Poiseuille mobility, the Couette thickness/face flux, the FBNS
mobility, the Elrod `Gamma_base` and the `VoidParams` closure; `reynolds.cpp`,
`fbns.cpp` and `diagnostics.cpp` all assemble from it. Previously `Diagnostics`
used `h^3/(12 mu)` with the scalar `bulk_modulus` and plain `h` where the
assembly had moved to `F2`, `beta_bar` and `2 F1/F0`, so `mass_residual` reported
a fabricated imbalance on the cross-film and `VOID_COUPLED` runs it exists to
police. `EOS::switch_function` in the diagnostics also now reads
`VoidParams::theta_full` instead of a hard 1.

Verified: `tests/test_fbns.cpp` test 14 (F2 differs from `h^3/12mu` by 9.7 % and
`2F1/F0` from `h` by 8.8 % on the test bearing, so the branch genuinely engages;
FBNS and Elrod cross-film temperatures agree to 0.48 K over the shared full-film
arc; per-cell mass residual `2.3e-10` of the peak Couette flux over every
non-groove cell, fronts included). With `cross_film_thermal_model = NONE` the
FBNS path is byte-identical to `6dacec0` on the Ferron THD Fig. 4 configuration
(720 x 192, `mpiexec -n 4`, all four `.vts` pieces plus `diagnostics.csv`,
`mass_residual` included in the output field set).

Still open on this branch's topic:
  - `VoidParams::theta_full` and `FullFilm::bound` remain two ceilings (unchanged
    from the phase-3 note above); the diagnostics now follow the assembly on
    both, but merging them is still separate work.
  - §13.1.2's stated approximations are unchanged and now apply to FBNS runs too:
    the midplane collocation is Couette-only, the convective term is saturated
    against the dissipative one, and where that saturation binds the quartic
    degenerates towards a fixed flattened shape rather than resolving the wall
    layers.

Film-surface output fields — **DONE 2026-08-05** (`codex/fbns-crossfilm`;
PHYSICS.md §13.1.2 and the new §15.2, README.md, CHANGES.md, `ctest` 27/27).
`temperature_journal`, `temperature_bearing`, `heat_flux_journal` and
`heat_flux_bearing` expose the two ends of the cross-film profile next to the
depth-averaged `temperature`, which is what a wall thermocouple or a bush
heat-flux measurement can actually be compared against. Opt-in through
`output_fields`, allocated only when selected (unselected runs stay
byte-identical), and a selection under `cross_film_thermal_model = NONE` is a
validation error rather than a silent field of zeros. The fluxes are positive
INTO the fluid, the negative of the module's internal `q_shaft` / `q_bush`.

Conjugate bush, 2-D (r, theta) conduction — **DONE 2026-08-05**
(`codex/conjugate-bush`; `src/bush_conduction.{hpp,cpp}`, PHYSICS.md §13.1.2.1,
README.md, CHANGES.md, `tests/test_bush_conduction.cpp`, `ctest` 30/30).
`bush_thermal_model = CONDUCTION_2D` (default `LOCAL`) replaces the per-cell 1-D
radial bush closure with steady conduction in the bush annulus, which is Ferron's
own bush model and supplies the mechanism the wall-BC study of case
01-THD-crossfilm proved missing: circumferential redistribution inside the bush.
The film-side flux is affine in the bore temperature, so it enters as a Robin
condition and the annulus stays a linear system; it is solved directly (block
Thomas in theta plus one dense Schur solve for the periodic wrap, ~10^7 flops at
720 x 16) for both right-hand sides at once, so `T_b = alpha + beta T_shaft`
remains exact and the shaft balance still closes in one step. The bush is
z-collapsed as the paper's is, and the axial end faces are not modelled (~1 %).
`bush_outer_boundary` now expresses a fixed outer surface directly instead of as
a large-HTC limit, for both bush models. `LOCAL` verified byte-identical to
`b74db94` on two cross-film configs, whole `results/` tree.

Ferron polynomial viscosity — **DONE 2026-08-05** (`codex/viscosity-polynomial`;
PHYSICS.md §13.2.1 and §17.5, README.md, CHANGES.md, `ctest` 27/27).
`liquid_viscosity_model = POLYNOMIAL` selects
`mu(T) = mu (k0 - k1 Tbar + k2 Tbar^2)` with `Tbar = T[degC]/T_ref[degC]`, the
law Ferron, Frene & Boncompain (1983) actually used, so the validation cases can
run that paper's viscosity instead of an Andrade fit to it. New keys
`viscosity_poly_k0/k1/k2` default to `1, 0, 0` (identity), so the model is inert
until configured. **The nondimensional temperature is a CELSIUS ratio**, verified
against the paper text: its nomenclature reads "T = temperature degC", and its
Table 1 set (3.287, 3.064, 0.777) gives `mubar = 1.000` at 40 degC and `0.267` at
80 degC only under that reading — a Kelvin ratio gives 0.82, which is not the
oil. Verified live under both FBNS+`ENERGY_EQUATION` and Elrod+`ENERGY_EQUATION`
(`tests/test_fbns.cpp::test_polynomial_viscosity_couples`).

Still open, deliberately not done here:
  - No validation case is re-run on `POLYNOMIAL` yet. Cases 01-THD and 04 sit on
    the Ferron geometry and should be repeated with the paper's own viscosity
    law before their temperature/load gaps are attributed to anything else.
  - The law is temperature-only. Ferron's rig ran to ~70 bar, where a Barus term
    is not negligible; combining `POLYNOMIAL` with a pressure factor would need
    a new model, not a new coefficient.

Case 00c (JFO solver verification), test #3 — per-cell mass-residual output
field: **DONE**. `mass_residual` output field added (`src/diagnostics.{hpp,cpp}`,
PHYSICS.md §15.1, `tests/test_mass_residual.cpp`). The 00c converged state is
per-cell conservative at round-off away from the groove penalty cells, so the
grid-divergent rupture front is not a conservation defect. Uncommitted pending
user review. Test #1 (exact 1D benchmark) still open; test #2 waits on published
contour data.

### Validation campaign (2026-08-02) — see `VALIDATION_PLAN.md`

The full diagnosis of the Ferron/Brito temperature and load mismatches, and the
phased task list for other agents, live in `VALIDATION_PLAN.md`. Summary of
status; update both files together.

- **Phase 0 (guardrails) — DONE 2026-08-02.** KSP convergence guard (defect D6)
  with the `test_ksp_guard` regression test; config warnings for the silent
  THD-decoupling trap and for adiabatic walls; `t_supply` validation; deferred
  correction unstaled in the energy and gas paths (B2); iteration counts recorded
  by `run_brito.py`. 23/23 tests pass. No physics changed.
- **Phase 1 (reference-data verification) — TODO, highest value, no code
  changes.** Ferron L / radial clearance / the speed of the eccentricity-vs-load
  figure (the source panel does not label one; 1000 rpm fits far better than the
  assumed 2000). Brito Fig. 6 polar re-digitisation, hot running clearance, and
  wiring up the dead `solve_eps_for_load`. An R&B regression at Ferron's
  L/D = 0.8.
- **Phase 2 (energy-equation fixes) — 2.1-2.4 IMPLEMENTED 2026-08-02; 2.5 TODO.**
  2.1/2.2 theta-weighted the cavitated thermal inventory (D14) and measured **no
  steady effect** — the source and the mass flux scale with theta identically and
  cancel, so D14 is ruled out as a contributor. 2.3 added hot-oil carry-over at
  grooves (`groove_mixing_model`). 2.4 added the cross-film `T(y)` profile,
  conjugate shaft/bush walls and the generalized (Dowson) Reynolds equation as one
  switch, `cross_film_thermal_model` — new module `src/cross_film.{hpp,cpp}`, new
  `tests/test_cross_film.cpp`, 24/24 tests pass. Everything is off by default and
  bit-identical to the previous solver when off (verified byte-for-byte against a
  HEAD build on JFO/GUMBEL/transient/4-rank). Two approximations inside the
  order-4 collocation are forced by measured non-convergence, not by convenience,
  and are written up in PHYSICS.md §13.1.2 and CHANGES.md.
- **Phase 3 (grid and time-step independence) — TODO**, blocked on Phase 2.
  Must scale `beta` per grid level to hold `p_max/beta >= 4e-3`, or the rupture
  angle error masquerades as grid error.
- **Phase 4 (validation matrix) — TODO.** JFO (adds Engels 2024 measured void
  fraction, the direct Elrod observable), energy equation (Ferron + full Brito
  matrix), gaseous cavitation (blocked on reformulating the VOID_COUPLED
  coupling below).

### WP-1b — VOID_COUPLED made usable (2026-07-31)

Driver: the optical journal-bearing rig study (50 mm bore, c = 250 um,
L/D = 0.75, white mineral oil, air vs CO2 at ~1 atm) needs the released-gas
volume to influence the cavity, because that is the only mechanism by which the
dissolved species can matter. Under `NONE` the gas is volumetrically inert and
air and CO2 give bit-identical fields; `VOID_COUPLED` was the intended path but
diverged (10 831 bar) or hung.

Done:
- Bounded, closed-form plateau `p_void = L theta p_sat / ((1-theta) + L theta)`
  replacing the divergent explicit `m_g R_s T / ((1-theta) h)` (PHYSICS.md §16.3).
- `void_coupling_relaxation` damping the `alpha_gas <-> theta_full <-> theta`
  limit cycle in the inner Elrod loop. 22/22 tests pass, WP-1 gates unchanged.
- Analytical benchmarks in `validation/journal_bearing`: infinitely-long
  (Sommerfeld) and short-bearing (Ocvirk) limits, both converging to the closed
  forms at the expected rates.

Open — **VOID_COUPLED still does not converge and is not usable quantitatively**:
- The outer `equilibrate_gas` <-> `solve_elrod` Picard coupling exhausts
  `coupling_max_iters` at `R_p` between 0.57 and 7.2 (tolerance 1e-9). Answers
  are relaxation-dependent (w = 0.10/0.30/0.60 give p_max 12.1/6.29/6.42 bar at a
  fixed grid), which proves they are not fixed points. Fixing this is a
  numerical-methods change: replace the segregated Picard iteration with an
  implicit/Newton treatment, or fold the flash into the Elrod switch function so
  the void fraction and film content are solved together rather than alternately.
- The solver reports success on iteration-budget exhaustion. `ElrodStats` already
  carries `converged`; the run should surface a non-zero exit or an explicit
  warning when the coupling fails, so that "ran to completion" cannot be mistaken
  for "converged".
- `equilibrate_gas` still assigns free gas independently of film content, so
  `alpha_gas` is a flash diagnostic rather than a transported void fraction.
- Remaining validation targets (status as of 2026-07-31, see
  `validation/journal_bearing/README.md`):
  - **Raimondi & Boyd finite-width charts — DONE, passes.** With a fed film,
    S within 5.8% and attitude angle within 4.08 deg over eps = 0.2-0.9 at
    L/D = 1/2 and 1. Ungrooved, the same solver has no load capacity at all.
  - **Ferron 1983 measured pressure and temperature — DONE, partial.** See the
    README verdict; the geometry in circulation as L/D = 0.55, C/R = 1.5e-3 is
    wrong (it is L/D = 0.80, c = 145 um).
  - **Cole & Hughes cavity photographs — NOT SOURCED.** Paywalled with no numeric
    reproduction reachable; Jakobsson & Floberg 1957 was attempted as a
    substitute and is also unavailable. Independently, the reformation angle is
    not a prediction of this solver in either available inlet configuration, so
    that half of the target cannot be scored even with the data.
- **Next validation targets, located and screened 2026-07-31.** Verdict 2 (energy
  equation) currently rests on a single experiment, and the cavitation model has
  no void-fraction test at all. Three reachable datasets close both gaps:
  - **Brito et al., ASME J. Tribology 129(1), 2006, DOI 10.1115/1.2401206** —
    open access at repositorium.sdum.uminho.pt. 100 mm plain bearing, bronze
    bush, two axial grooves at ±90 deg to the load, diametral clearance
    171 ± 5 um at 20 C. Measures mid-plane pressure, oil–bush interface
    temperature, **oil–shaft interface temperature via a rotary joint**, flow
    rate and bush torque, swept over speed to 4000 rpm, 2 and 10 kN, three
    supply temperatures and three supply pressures. The shaft-side temperature
    is the quantity pancake's fixed-wall energy model cannot produce and which
    is the suspected cause of the Ferron temperature-peak error; the supply
    sweep addresses the inlet-mixing difference against Alakhramsing; and the
    clearance is quoted diametrally with tolerance and reference temperature,
    removing the ambiguity Ferron leaves. Full text verified, not just screened.
  - **Brito, PhD thesis, U. Minho 2009, hdl:1822/9884** — open access, likely
    the tabulated dataset behind the above.
  - **Engels et al., Lubricants 13(1), 2024, DOI 10.3390/lubricants13010007** —
    open access, **experimental void fraction** in a hydrodynamic bearing. This
    is the quantity the Elrod model actually predicts, and the gap Cole & Hughes
    was meant to fill. Screened by abstract only.
  Lund & Tonnesen 1984 (ASME 106(2):237) remains the other V6 target and is
  cited by Brito as the key prior two-groove work; not open access.
- **Elrod free-boundary accuracy depends on `p_max/bulk_modulus`.** Because
  pressure is carried inside the film content, the rupture angle can be 20 deg
  wrong at `beta = 1e9` on a 720-cell mesh while every convergence indicator
  reports success. Any new cavitated case needs a mesh/beta check, not just a
  residual check. `run_raimondi_boyd.convergence_study()` is the harness.

### Deferred modelling gaps (optical-rig study, 2026-07-31)

Recorded here so they are not silently inherited as validated physics.

1. **No vaporous-cavitation physics.** The "vaporous" case is produced only by
   setting `dissolved_gas_initial = 0` so the plateau falls back to the scalar
   `p_cav`. There is no latent heat, no phase-change kinetics and no vapour
   transport - it is a prescribed low-pressure plateau, not a vaporisation model.
   Braun & Hendricks (1984) argued oil vaporisation needs energy not available in
   normal bearing operation, so a real degassed rig should sit between the
   gaseous and "vaporous" predictions rather than at the 0.005 bar extreme. Do
   not quote the vaporous case as a validated prediction of vapour behaviour.
2. **Reabsorption kinetics are symmetric.** `update_gas_state` uses the single
   `gas_mass_transfer_rate` for both release and resorption, and the steady
   `EQUILIBRIUM` path makes both instantaneous. Physically resorption is
   diffusion-limited and much slower than release, which is why gaseous cavities
   persist downstream of reformation while vapour cavities collapse at once.
   Etsion & Ludwig (1982) proposed exactly this asymmetry as the mechanism behind
   the ~50 kPa pressure variation they measured *inside* the cavitation zone. A
   separate, slower resorption time constant would be needed to reproduce it.
3. **Cavity pressure is assumed uniform.** `LOCAL_PSAT` pins every cavitated cell
   at `p_sat`. Etsion & Ludwig measured ~50 kPa of variation across the cavity,
   which is ~10% of a 0.5 bar plateau. Predicted cavity *extent* is trustworthy;
   predicted cavity *pressure* is nominal.
4. **No axial re-flooding of a cavitated film.** Elrod-Adams sets the switch
   function to zero in cavitated cells, so there is no pressure-driven path for
   liquid to re-enter. An ungrooved bearing therefore falls irreversibly into the
   analytic fully-starved state `theta*h = const` (verified: rms 0.0065 against
   the closed form, `theta_min` 0.1766 vs 0.1765). Etsion & Ludwig photographed
   strong reverse flow re-flooding the cavity in a real submerged bearing, so
   this is a model limitation, not a physical prediction. Ungrooved geometries
   need an explicit inlet as a numerical stand-in.
5. **`p_sat` ceiling was geometry-specific.** The `p_sat <= 0.7 bar` limit found
   for the ungrooved case does not apply once a groove is present: at eps 0.5 /
   1500 rpm the grooved bearing runs cleanly at `p_sat = 1.0 bar` (p_max 3.12
   bar, arc 173 deg). Above `p_sat ~ 0.9 bar` the reformation front reaches the
   groove edge, so reformation angle stops being an independent observable.

## Completed Phases

### Phase 1 — Infrastructure
Mesh, Field, Communicator, IO, FVM operators, PETSc linear system. All working with MPI domain decomposition in $\theta$.

### Phase 2 — Physics Foundations
Film thickness (static eccentric), barotropic EOS ($p = p_{cav} + \beta\ln\theta$), VTK output.

### Phase 3 — Reynolds Solver (Gumbel)
Steady-state compressible Reynolds equation with lagged density. Gumbel cavitation (post-solve pressure clamp). BiCGStab + block Jacobi. All tests pass.

### Windows GUI Infrastructure
Native Win32 GUI target with tabbed parameter controls, output-field selection,
solver launch/stop, process log capture, and lightweight VTK heatmap preview.
This phase is intentionally an external orchestrator around `config.txt`,
`pancake.exe`, and generated VTK files so the numerical core remains shared.
The active Windows workflow is now MSYS2 MinGW64: one CMake preset builds
`pancake.exe`, `pancake_gui.exe`, tests, and copied runtime DLLs. The preset
also injects the MinGW bin directory into PATH for PowerShell builds. The main
`Workspace` tab uses a left/right layout: result preview and process log on the
left, and a scrollable one-row-per-parameter inspector on the right. Solver
backend controls are always visible in the top-right header across all tabs.
The preview now supports timestep selection, axes, a unit-labeled colorbar,
literal numeric edit preservation, raw Elrod/JFO `film_content` output, and a
`load_angle_deg` visualization reference measured from the positive y
axis. The GUI assembles multi-rank preview fields, groups axial boundary
controls by south/north side, uses structured inlet controls for groove/circular
inlet definitions, and provides per-value unit selectors while saving
solver-native config units. Both Windows executables embed the shared Pancake
logo icon so the GUI and backend present the same application identity. Elrod
axial boundaries are now pressure inputs in the GUI and config writer; the
solver derives film content internally from pressure, `p_cav`, and
`bulk_modulus`.

Follow-up correction (2026-06-08): deprecated manual `bc_z_*_theta` boundary
inputs were removed from generated configs and the GUI. The parser still
accepts old keys, but Reynolds, velocity/force post-processing, energy boundary
fluxes, and initialization use pressure-derived Elrod boundary values. Tests
now guard both the pressure-derived theta and the ignored legacy theta path.

Follow-up correction (2026-06-08): output-directory setup now fails with a
clear solver message when Windows cannot clear a locked `results/...` folder,
instead of terminating on an uncaught filesystem exception. Users should close
ParaView/preview handles or select a different `output_dir` before rerunning.

### GUI Modernization (Dear ImGui docking + ImPlot, 2026-06-11)
Rewrote the GUI per `GUI_MODERNIZATION_DIRECTIVE.md` as `src/gui/` on Win32 +
Direct3D 11: worker-thread solver subprocess with a live log-scale convergence
plot, ImPlot field heatmap with cavitation-zone overlay, run history with
overlay comparison, CSV/PNG export everywhere, validation that merges the
solver's own `SimulationConfig::validate()`, and a single statically linked
exe (embedded fonts, `%LOCALAPPDATA%\PancakeGui` settings, per-monitor-v2
DPI). See `AGENT_REPORT.md`, `BUILD.md`, `docs/GUI_ARCHITECTURE.md`.
Remaining follow-ups: human smoke test of the interactive workflow; delete
legacy `src/gui_win32.cpp` once this branch's pending edits to it are
settled; solver-side proposals (structured `PROGRESS` stdout line,
cooperative cancellation) remain open in `AGENT_REPORT.md` section 6.

### Windows GUI Redesign (2026-05-29)
Rebuilt the GUI UI layer around a single tab-less workspace (action/run bar,
collapsible input rail with a pinned live-summary card, result viewport, console)
to support fast tweak→run→inspect loops. Replaced group-box frames with painted
section headers/labels to remove the long-standing z-order/paint glitches. Added
live derived journal-bearing quantities and inline validation that gates `Run`,
a run-state chip with a stdout-driven progress bar, hybrid units (global SI ↔
engineering selector plus targeted per-field overrides), scenario presets, and
custom-config Open/Save via Windows file dialogs. `pancake_gui.exe` now links
statically into a single self-contained binary (no runtime DLLs deployed beside
it). The solver/preview/log/unit backends were preserved; `config.txt` format
and the solver contract are unchanged.

### Windows GUI Motion/Preview Update (2026-06-04)
Added structured Motion / Loads controls for moving-bearing selection,
pressure/motion/temperature time-method selection, external-load
magnitude/direction/z input, bearing initial state, support
stiffness/damping, and film-thickness safety settings. The GUI now keeps the
current contour stable while the solver is running and refreshes results only
on completion or explicit user refresh. Preview timestep switching now uses a
cached per-step list of `processor*` VTS files instead of recursively scanning
the whole results tree for each field/step change.

### Live Summary Numerical Diagnostics (2026-06-07)
Added a nominal circumferential Couette Courant number to the live summary:
$Co_\theta = |\omega|\Delta t/(2\Delta\theta)$ with
$\Delta\theta = 2\pi/N_\theta$. This is an input-summary diagnostic for the
surface-driven Reynolds advection term; pressure-driven Poiseuille velocities
remain post-solve quantities. Added `docs/PROJECT_CONTRIBUTIONS.md` to summarize
substantive solver, model, output, testing, and documentation contributions
without listing one-off repairs or GUI construction work.

### Elrod MPI Active-Set Synchronization (2026-06-07)
The Elrod-Adams outer active-set loop now exchanges solved `theta` ghost cells
after every outer solve, clamp, and near-full-film snap. This keeps
$g(\theta)$ current at MPI partition interfaces when running with multiple
ranks. Added a post-solve ghost-consistency regression and a 6-rank
`test_elrod_a2_mpi6` CTest entry because the visible contour artifact appeared
with six processor pieces.

### Solver Initialization Update (2026-06-04)
The timestep flooded guess now initializes `pressure` from the first configured
inlet `p_supply`, not `p_cav`. The Elrod `theta` field is initialized from the
same pressure through the EOS, with `p_cav` used only as the no-inlet fallback.

### Bearing Force Convention Update (2026-06-04)
Corrected macroscopic force integration to use one bearing-side convention:
`pressure_force_*` is the pressure force applied to the moving bearing surface,
including the bearing surface area vector
`((R+h)e_r - h_theta e_theta - (R+h)h_z e_z) dtheta dz`; `viscous_force_*` is
the bearing-side shear force; `fluid_force_*` is their sum. Legacy `load_*`
output requests map to `pressure_force_*`, but duplicate runtime `load_*`
fields are no longer allocated.

### Energy Equation Initial Implementation (2026-06-04)
Added `solution_mode = STEADY_STATE` for one-step operating-point calculations
and `temperature_model = ENERGY_EQUATION` for a first-pass film-averaged thermal
solve. The solver now carries `temperature` and `heat_generation` fields,
computes Couette plus pressure-flow viscous heat generation, applies wall heat
transfer to the journal and bearing, and exposes the thermal controls in the
Win32 GUI. The energy solve uses non-conservative upwind thermal advection to
avoid spurious `T div(F)` source terms near pressure-supply regions, treats
pressure supply inlets as pressure-only constraints, and uses
`temperature_reference` only for axial pressure-boundary thermal inflow.
Transient runs now assemble the `t=0` solve as steady and store that solved
state before applying transient pressure/temperature terms at `t=dt`.
Thermoviscous viscosity feedback, deformation, and coupled ETHD iterations
remain in Phase D.

Follow-up correction (2026-06-07): `heat_generation` now evaluates the
pressure-flow dissipation from active full-film faces instead of one-sided
cell-centered gradients at a cavitation front. This keeps the hard
Elrod-Adams/JFO pressure plateau intact while preventing inactive cavitated
faces from creating mesh-dependent corner heat spikes in the derived thermal
source.

### Oil/Propane Mixture Property Layer (2026-06-05)

Added the first JFO-compatible dissolved-gas implementation slice. The solver
now owns `src/fluid_properties.hpp/cpp`, with config-driven enums for
`fluid_property_model`, `dissolved_gas_species`, `oil_gas_solution_model`,
`density_model`, and `viscosity_model`. Defaults remain
`CONSTANT`/`PURE_OIL`, so existing pure-oil cases keep their old numerical
behavior.

Implemented fields:
- `rho_liquid_solution`, `mu_liquid_solution`, `cp_liquid_solution`,
  `rho_cp_liquid_solution`, `k_liquid_solution`, `mu`, and `rho_gas`.
- `dissolved_gas`, `free_gas_mass`, `alpha_gas`, and `gas_mass_transfer`.

Current algorithm:
1. Refresh liquid oil plus dissolved-gas properties from pressure,
   temperature, and dissolved gas.
2. Feed field-valued `rho_liquid_solution` and `mu` into Gumbel/Elrod,
   velocity reconstruction, force post-processing, heat generation, and thermal
   convection.
3. Solve the Reynolds/JFO step with `theta` still acting as the liquid-content
   variable.
4. For `GAS_CAVITATION_MIXTURE`, transport free gas with the film velocity and
   apply a finite-rate release/resorption update between dissolved gas and free
   gas.
5. Cap the reported `alpha_gas` by `gas_alpha_max`; do not delete free gas just
   because a cell is locally full-film under JFO.

This is deliberately a segregated JFO nonlinear corrector plus property/gas
coupling loop. It is not labelled PISO because the codebase does not yet have a
Navier-Stokes pressure-velocity backend.

Tests added/updated:
- `test_fluid_properties`: config aliases, solubility monotonicity, mixture
  density bounds, viscosity selection/table interpolation, the calibrated
  R290/PZ68 40 C / 10 bar point, saturated zero source, release, resorption,
  and bounded `alpha_gas`.
- `test_gas_transport`: dissolved-gas transport constant preservation,
  disabled-model and zero-timestep no-ops, bounded advection, and
  mass-weighted conservation.
- `test_io`: legacy `load_x` output selection writes `pressure_force_x` without
  a duplicate `load_x` array.
- `test_forces`: force diagnostics no longer require legacy `load_x`.

Deferred work:
- Boundary-condition-rich dissolved/free gas transport beyond the first axial
  open-boundary and pressure-inlet reservoir-composition contract; future work
  should add selectable gas-specific BCs where a case genuinely needs them.
- Additional empirical propane/oil table calibration for temperatures and
  pressures beyond the supplied chart point.
- Bubble-radius dynamics, slip, and a true pressure-velocity PISO backend.

Follow-up correction (2026-06-08): added shared CLI/GUI validation for
absolute-pressure operation and gas-model dependencies. Solver runs now fail
before launch when `p_cav` is left at gauge zero, when gas cavitation is paired
with a non-JFO cavitation model, or when required Henry/density/viscosity/gas
transfer coefficients are missing. Added `config_r290_pz68_quick.txt` as a
small smoke case using the same R290/PZ68 constants.

Follow-up correction (2026-06-08): corrected the first gas/free-surface coupling
semantics. Free gas is no longer clipped to the JFO void fraction before
resorption, which silently deleted gas in pressure-inlet/full-film cells without
a negative `gas_mass_transfer` source. `free_gas_mass` is now conservatively
transported with the film-averaged velocity, and `alpha_gas` is capped only by
`gas_alpha_max` while resorption is controlled by the finite-rate transfer law.

Follow-up correction (2026-06-08): added the first open gas-composition
boundary contract. Pressure supply/inlet cells impose the configured incoming
mixture (`dissolved_gas_initial`, zero `free_gas_mass`, zero `alpha_gas`, zero
local `gas_mass_transfer`). Axial pressure-boundary inflow carries the same
initial dissolved-gas composition with zero free gas, and axial outflow convects
released gas out of the domain. Focused tests now guard inlet composition,
axial dissolved-gas inflow, and axial free-gas outflow.

---

## Phase A — Mass-Conserving Cavitation (Elrod-Adams)

Replace the Gumbel pressure clamp with the Elrod-Adams formulation. The universal variable $\theta$ (fractional film content) becomes the primary unknown. Mass is conserved across the cavitation boundary by construction.

### A.0 — Codebase Audit and Preparation

**Status: complete** (this document).

Findings:
- `EOS::dp_dtheta()` already exists but is unused — will become the core of $\gamma_{eff}$.
- `FVM::divergence` with upwind/TVD is implemented and tested — will be used for the Couette convection term.
- `FVM::ddt` is implemented — needed for transient $\partial(\theta h)/\partial t$.
- `LinearSystem` supports 5-point stencil with periodic $\theta$ — no changes needed.
- Current `Reynolds::solve` operates on pressure as primary variable — must be refactored.

---

### A.1 — Universal Variable Reformulation

**Goal:** Reformulate the Reynolds equation with $\theta$ (film content) as the primary variable. No cavitation modelling yet — full-film everywhere ($\theta \ge 1$). Validate that this reproduces the current Gumbel-based pressure field.

**Mathematical formulation:**

Starting from the compressible Reynolds equation:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta} + \frac{\partial(\rho h)}{\partial t}$$

Substitute $\rho = \rho_0 \theta$ and apply the chain rule $\partial p / \partial x_i = (dp/d\theta)\,\partial\theta/\partial x_i$:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\Gamma \frac{\partial\theta}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\Gamma \frac{\partial\theta}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho_0 \theta h)}{\partial\theta} + \frac{\partial(\rho_0 \theta h)}{\partial t}$$

where the **effective diffusion coefficient** is:

$$\Gamma(\theta) = \frac{\rho_0 \theta \, h^3}{12\mu} \cdot \frac{dp}{d\theta} = \frac{\rho_0 \theta \, h^3}{12\mu} \cdot \frac{\beta}{\theta} = \frac{\rho_0 \beta \, h^3}{12\mu}$$

Note: for full-film ($\theta \ge 1$), $\Gamma$ is independent of $\theta$. This simplifies A.1 — the equation is linear in $\theta$.

The **Couette term** $\frac{\omega}{2}\frac{\partial(\rho_0 \theta h)}{\partial\theta}$ is a one-dimensional divergence of the flux $\frac{\omega R h}{2} \cdot \rho_0\theta$. In the FVM framework, this is a **convection** of $\theta$ with known face fluxes:

$$F_e^{(\theta)} = \frac{\omega R \, h_e}{2} \cdot \Delta z, \qquad F^{(z)} = 0$$

where $h_e$ is the film thickness interpolated to the east face. The `FVM::divergence` operator with upwind (or TVD) scheme handles this directly.

The **transient term** $\partial(\rho_0 \theta h)/\partial t = \rho_0 h \,\partial\theta/\partial t$ (for static $h$) is discretised as a weighted `ddt` with $h$ as a cell-volume multiplier.

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.hpp/cpp` | New function `Reynolds::solve_elrod(fields, sys, mesh, cfg)`. Keep `Reynolds::solve` as a fallback. |
| `config.hpp` | Add `enum CavitationModel { GUMBEL, ELROD_ADAMS }` and a config field to select it. |
| `main.cpp` | Add `"theta"` field. Initialise $\theta = 1$ (full film). Dispatch to the appropriate solver. |
| `fvm.hpp/cpp` | Add `FVM::ddt_weighted(sys, phi, weight, dt, mesh)` — like `ddt` but multiplies the cell volume by a per-cell weight field (for $h \cdot \partial\theta/\partial t$). |
| `equation_of_state.hpp` | No change — `dp_dtheta` already exists. |

**Solve sequence per timestep (A.1, full-film only):**

1. Compute $\Gamma = \rho_0 \beta h^3 / (12\mu)$ at physical + ghost cells.
2. Compute Couette face fluxes: $F_e(i,j) = \omega R \cdot h_e / 2 \cdot \Delta z$.
3. Reset system. Assemble:
   - `FVM::laplacian(sys, Gamma, mesh)` — Poiseuille diffusion
   - `FVM::divergence(sys, couette_flux_theta, zero_flux_z, theta, UPWIND, mesh)` — Couette convection
   - Dirichlet z-BC: $\theta = \theta(p_{cav}) = 1$ at $z = 0, L$
   - (Optional) `FVM::ddt_weighted(sys, theta, h_field, dt, mesh)` — transient
4. Solve $A\,\theta^{n+1} = b$.
5. Recover pressure: $p = p_{cav} + \beta \ln(\max(\theta, 1))$.
6. Update density: $\rho = \rho_0 \theta$.

**FVM sign convention check:**

The assembled equation is:

$$\underbrace{-\text{div}(\Gamma\nabla\theta)\cdot V}_{\text{laplacian}} + \underbrace{\text{div}(\mathbf{F}\,\theta)\cdot V}_{\text{divergence}} = \text{source}$$

The PDE $\text{div}(\Gamma\nabla\theta) = \text{div}(\mathbf{F}\,\theta) + \partial(\rho_0\theta h)/\partial t$ rearranges to:

$$-\text{div}(\Gamma\nabla\theta)\cdot V + \text{div}(\mathbf{F}\,\theta)\cdot V = \rho_0 h \frac{\partial\theta}{\partial t}\cdot V$$

The LHS matches. The transient goes to the RHS and is handled by `ddt_weighted` (which adds $h V/\Delta t$ to $a_P$ and $h V \theta^n / \Delta t$ to `source`).

**Validation:**

- **Regression test**: Run the static journal bearing from Phase 3. The $\theta$-formulation must produce the same pressure field as the Gumbel solver (since $\theta \ge 1$ everywhere in full-film regions; the Gumbel clamp has no effect when pressure is already above $p_{cav}$). Compare $\max|p_{Elrod} - p_{Gumbel}|$ < tolerance.
- **Unit test**: 1D channel with uniform $h$, $\Gamma = \text{const}$. Exact solution for $\theta(x)$ is linear between Dirichlet BCs. Verify 2nd-order convergence.

---

### A.2 — Elrod-Adams Switch Function and Cavitation

**Status: complete** (2026-04-12).

**Goal:** Enable the cavitation switch $g(\theta)$. In cavitated regions ($\theta < 1$), the Poiseuille (diffusion) term vanishes and only the Couette (convection) term transports liquid. This is the core of mass-conserving cavitation.

**Modified diffusion coefficient:**

$$\Gamma_g(\theta) = \frac{\rho_0 \beta \, h^3}{12\mu} \cdot g(\theta)$$

where the switch function $g$ satisfies the **Elrod-Adams complementarity conditions**:

$$g \in \{0, 1\}, \qquad \theta \ge 1, \qquad g(\theta - 1) = 0, \qquad p \ge p_{cav}$$

These imply:
- $g = 1 \Longleftrightarrow \theta \ge 1$ (full film, pressure active)
- $g = 0 \Longleftrightarrow \theta < 1$ (cavitated, pressure clamped to $p_{cav}$)

**Switch-function update strategy:**

Implement an **iterative flag update** with safeguards (Vijayaraghavan & Keith, 1989):

```
for outer_iter = 1 to max_outer:
    1. Compute g(i,j) from current theta:
       g(i,j) = 1  if theta(i,j) >= 1
       g(i,j) = 0  if theta(i,j) < 1
    
    2. Compute Gamma_g = rho0 * beta * h^3 / (12*mu) * g
    
    3. Assemble and solve for theta_new (same as A.1 but with Gamma_g)
    
    4. Post-process:
       - In full-film cells (g=1): if theta_new < 1, mark for switch
       - In cavitated cells (g=0): if theta_new > 1, mark for switch
    
    5. Check convergence:
       - No cells changed flag, AND
       - ||theta_new - theta_old|| / ||theta_new|| < tol
    
    6. Update theta = theta_new
```

This is simpler than Fischer-Burmeister or Murty pivoting. It works well for steady and mildly transient problems. If convergence is poor, Phase A.3 adds a more robust complementarity solver.

**Regularisation:**

To avoid division-by-zero in the cavitated region and improve convergence:
- $g_{reg}(\theta) = \tfrac{1}{2}(1 + \tanh((\theta - 1)/\epsilon))$ with $\epsilon \sim 10^{-3}$ can be used as a smooth approximation during initial iterations, then tightened to a hard switch.

**Boundary conditions in $\theta$-formulation:**

| Boundary | Current (pressure) | $\theta$-formulation |
|----------|-------------------|---------------------|
| $z = 0, L$ (Dirichlet) | $p = p_{cav}$ | $\theta = 1$ (since $\theta(p_{cav}) = 1$) |
| $\theta$ (periodic) | Ghost exchange for $p$ | Ghost exchange for $\theta$ |
| Supply groove (future) | $p = p_{supply}$ | $\theta = \exp((p_{supply} - p_{cav})/\beta)$ |

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Add outer iteration loop around the linear solve. Compute `g` field from $\theta$. Multiply $\Gamma$ by `g`. |
| `config.hpp` | Add `max_outer_iters`, `outer_tol`, `g_regularisation_eps`. |
| `equation_of_state.hpp` | Add `switch_function(theta)` and `switch_function_regularised(theta, eps)`. |
| Tests | New `test_elrod.cpp`. |

**Validation:**

- **1D inclined slider bearing** (Elrod & Adams, 1975): A linearly converging channel $h(x) = h_1 + (h_2 - h_1) x / L$ with $h_1 > h_2$. The analytical JFO solution gives:
  - Pressure rise in the converging zone ($\theta = 1$, $p > p_{cav}$)
  - Cavitation onset at $\theta = \theta^*$ where $\partial p/\partial x = 0$
  - Cavitated zone ($\theta < 1$, $p = p_{cav}$) with linear $\theta$ decrease
  - Reformation where $\theta$ returns to 1

  Validate: compare computed cavitation boundary position, peak pressure, and $\theta$ profile against the analytical solution.

- **Journal bearing**: Compare cavitation extent (angular range where $\theta < 1$) against published data (e.g., Vijayaraghavan & Keith 1989, Fesanghary & Khonsari 2011). Verify that the sum $\int \rho_0 \theta h \, dA$ is conserved to machine precision between timesteps (mass conservation check).

- **Mass conservation diagnostic**: At the end of each solve, compute and print:
  $$M_{error} = \left|\sum_f F_f \theta_f - \text{(transient contribution)}\right| \bigg/ \sum_f |F_f|$$
  This should be $O(10^{-12})$ (solver tolerance), not $O(1)$ as with Gumbel.

---

### A.3 — Robust Complementarity Solver (if needed)

**Goal:** If the iterative flag update from A.2 shows poor convergence (oscillating flags, slow outer iteration), implement a more robust complementarity strategy.

**Status: skipped** (A.2 converges well for current cases).

---

### A.4 — Text-based Configuration

**Goal:** Move simulation parameters from a hardcoded header to a text-based config file to avoid constant recompilation.

**Proposed format:** Simple key-value pairs (`key = value`).

**Code changes:**
| File | Change |
|------|--------|
| `config.hpp` | Add `load_from_file(std::string path)` method to `SimulationConfig`. |
| `main.cpp` | Parse command line for config path (default: `config.txt`). |

---

### A.5 — Fix Numerical Starvation (Cascade Collapse)

**Status: complete** (Flooded initialization implemented in solve_elrod).

---

### A.6 — Flexible Inlet Boundary Conditions

**Goal:** Allow users to define internal oil supply features (holes and grooves) via config.

**Status: complete** (Implemented via penalty method in linear assembly).

---

### A.7 — Staggered Grid Velocity Fields

**Status: complete** (Implemented in Reynolds namespace, interpolated for VTK output).

---

### A.8 — Configurable Outlet Boundary Conditions

**Goal:** Allow users to choose between Fixed Value (Dirichlet) and Zero Gradient (Neumann) at the axial ends.

**Status: complete** (Implemented in Reynolds solvers and velocity calculation).

---

## Phase B — Dissolved and Free Gas Transport

Model gaseous cavitation: dissolved gas comes out of solution when pressure drops, forming free gas bubbles that alter the mixture properties.

### B.1 — Fluid Property Model (Mixture EOS)

**Goal:** Replace the single-component barotropic EOS with a two-phase mixture model that accounts for liquid, dissolved gas, and free gas.

**Mixture state variables (per cell):**

| Variable | Symbol | Description |
|----------|--------|-------------|
| Liquid volume fraction | $\alpha_l$ | $= \theta$ from Phase A (for cavitated cells, $\alpha_l < 1$) |
| Free gas volume fraction | $\alpha_g$ | Bubbles in the gap |
| Dissolved gas mass fraction | $c_d$ | Gas dissolved in the liquid (by mass) |

Constraint: $\alpha_l + \alpha_g = 1$ (two-phase mixture; vapour cavitation is already handled by $\theta$).

**Mixture density:**

$$\rho_{mix} = \alpha_l \rho_l(p) + \alpha_g \rho_g(p)$$

where:
- $\rho_l(p) = \rho_{l,0} \exp((p - p_{ref})/\beta_l)$ — liquid (existing barotropic EOS)
- $\rho_g(p) = p / (R_g T)$ — ideal gas (isothermal)

**Mixture viscosity** (Einstein or Krieger-Dougherty model for bubbly flow):

$$\mu_{mix} = \mu_l (1 - \alpha_g)^{-[\mu] \alpha_{g,max}}$$

where $[\mu] \approx 2.5$ (Einstein) and $\alpha_{g,max} \approx 0.64$. For dilute gas ($\alpha_g \ll 1$), this simplifies to $\mu_{mix} \approx \mu_l (1 + 2.5\alpha_g)$.

**Henry's law** (dissolved gas equilibrium):

$$c_{d,eq}(p) = H \cdot p$$

where $H$ is the Henry's law constant (units: $\text{kg}_{gas}/(\text{kg}_{liquid}\cdot\text{Pa})$). At equilibrium, the dissolved gas concentration equals $c_{d,eq}$. When $p$ drops below saturation, dissolved gas is released as free gas.

**Code changes:**

| File | Change |
|------|--------|
| `src/equation_of_state.hpp` | Refactor into `src/fluid_properties.hpp`. Add `MixtureEOS` class with methods: `rho_liquid(p)`, `rho_gas(p, T)`, `rho_mix(alpha_g, p)`, `mu_mix(alpha_g)`, `henry_equilibrium(p)`. Keep the old `EOS` namespace as a thin wrapper for backward compatibility. |
| `config.hpp` | Add gas properties: $R_g$, $T$, $H$, $\mu_g$, $\alpha_{g,max}$. |
| Tests | `test_mixture_eos.cpp`: round-trip density, viscosity limits ($\alpha_g = 0 \to$ pure liquid, $\alpha_g \to \alpha_{g,max} \to$ high viscosity), Henry's law. |

**Validation:**
- Pure liquid limit ($\alpha_g = 0$): must recover Phase A results exactly.
- Pure gas limit ($\alpha_l = 0$): ideal gas EOS.
- Mixture density monotonicity in $p$.

---

### B.2 — Dissolved Gas Transport Equation

**Goal:** Add a transport equation for the dissolved gas mass fraction $c_d$ in the liquid phase.

**Governing equation:**

$$\frac{\partial}{\partial t}(\alpha_l \rho_l c_d h) + \text{div}(\alpha_l \rho_l c_d \mathbf{u}_l h) = -\dot{m}_{dg} h + D_{eff} \nabla^2(\alpha_l \rho_l c_d h)$$

where:
- $\mathbf{u}_l$ is the mean liquid velocity (Couette + Poiseuille contributions)
- $\dot{m}_{dg}$ is the mass transfer rate from dissolved to free gas (degassing source)
- $D_{eff}$ is an effective diffusion coefficient (molecular + turbulent mixing)

**Degassing/resorption source term:**

$$\dot{m}_{dg} = k_{dg} \, \alpha_l \rho_l \, (c_d - c_{d,eq}(p))$$

where $k_{dg}$ is a rate constant ($s^{-1}$). When $c_d > c_{d,eq}$ (supersaturated), gas is released ($\dot{m}_{dg} > 0$). When $c_d < c_{d,eq}$ (undersaturated), gas is resorbed ($\dot{m}_{dg} < 0$).

**Simplification for thin films:** In the lubrication limit, the mean liquid velocity in the $\theta$-direction is $\bar{u}_\theta = \omega R / 2 - h^2/(12\mu) \cdot \partial p / \partial s$. The convective flux uses this velocity.

**Code changes:**

| File | Change |
|------|--------|
| `src/gas_transport.hpp/cpp` | New module: `GasTransport::solve_dissolved(fields, sys, mesh, cfg)`. Assembles and solves the $c_d$ transport equation using existing FVM operators. |
| `main.cpp` | Add `"c_dissolved"` field. Call `GasTransport::solve_dissolved` after Reynolds solve. |
| `config.hpp` | Add `k_degassing`, `D_eff`, `c_d_initial`. |
| `fvm.hpp/cpp` | May need a general convection-diffusion helper, or assemble manually from existing operators. |

**Validation:**
- 1D advection test: uniform flow with known inlet $c_d$, verify advection speed.
- Equilibrium test: set $p = \text{const}$ everywhere, verify $c_d \to H \cdot p$ at steady state.
- Conservation: total dissolved mass $\int \alpha_l \rho_l c_d h \, dA$ plus transferred mass $\int \dot{m}_{dg} h \, dA \cdot \Delta t$ must balance.

---

### B.3 — Free Gas Volume Fraction Equation

**Goal:** Evolve the free gas volume fraction $\alpha_g$ based on degassing from the liquid and gas compressibility.

**Governing equation:**

$$\frac{\partial}{\partial t}(\alpha_g \rho_g h) + \text{div}(\alpha_g \rho_g \mathbf{u}_g h) = +\dot{m}_{dg} h$$

The free gas is transported by a mixture of Couette flow and pressure-driven expansion. For simplicity (no slip between gas and liquid), assume $\mathbf{u}_g = \mathbf{u}_l$ (homogeneous mixture).

**Code changes:**

| File | Change |
|------|--------|
| `src/gas_transport.hpp/cpp` | Add `GasTransport::solve_free_gas(fields, sys, mesh, cfg)`. Coupled to dissolved gas via $\dot{m}_{dg}$. |
| `main.cpp` | Add `"alpha_gas"` field. Initialise $\alpha_g = 0$ (no free gas initially). |

**Validation:**
- Mass conservation: total gas mass ($\int [\alpha_l \rho_l c_d + \alpha_g \rho_g] h \, dA$) must be constant in a closed domain (no flux BCs).
- Degassing test: uniform pressure drop at $t = 0$; verify exponential approach to equilibrium with rate $k_{dg}$.

---

### B.4 — Coupling to Reynolds Equation

**Goal:** Close the loop: the Reynolds equation uses mixture properties ($\rho_{mix}$, $\mu_{mix}$), and the gas transport equations use the pressure field from Reynolds.

**Coupled solve sequence per timestep:**

```
1. Sync ghost cells for all fields (theta, c_d, alpha_g, h)
2. Compute mixture properties:
   rho_mix(i,j) = (1 - alpha_g) * rho_l(p) + alpha_g * rho_g(p)
   mu_mix(i,j)  = mu_l * (1 - alpha_g)^(-2.5 * alpha_g_max)
3. Reynolds solve (Phase A solver, using rho_mix and mu_mix):
   Gamma_g = beta * h^3 / (12 * mu_mix) * g(theta)
   → solve for theta → recover p
4. Dissolved gas transport (B.2):
   → solve for c_d
5. Free gas transport (B.3):
   → solve for alpha_g
6. Check outer convergence:
   max(|theta_new - theta_old|, |c_d_new - c_d_old|, |alpha_g_new - alpha_g_old|) < tol
7. Repeat 2-6 if not converged (Picard outer iteration)
8. Store old-time values, advance timestep
```

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Accept `rho_mix` and `mu_mix` fields instead of scalar `cfg.rho` and `cfg.mu`. |
| `main.cpp` | Implement the coupled outer iteration loop. |
| `config.hpp` | Add `max_coupling_iters`, `coupling_tol`. |

**Validation:**
- **Decoupled limit**: set $\alpha_g = 0$ and $c_d = H \cdot p$ everywhere, disable gas transport. Must recover Phase A results exactly.
- **Journal bearing with gas release**: compare cavitation extent and peak pressure against published two-phase results (e.g., Ausas et al. 2009, Bayada & Chupin 2013).
- **Convergence study**: verify that outer iteration converges in $\le 10$ iterations for typical bearing conditions.

---

## Phase C - Moving Bearing Dynamics

**Status: complete** (2026-06-04, moving outer bearing with fixed shaft).

Allow the outer bearing surface to move in response to pressure and viscous
film forces, independent applied load, support stiffness, and support damping.
This turns the current static-eccentricity solve into a lagged fluid-structure
time integration while keeping the shaft fixed in space.

### C.1 - Time-dependent Film Geometry

**Goal:** Replace static eccentricity inputs with bearing-center state variables
when dynamic motion is enabled.

For bearing-center displacement $(x_b(t), y_b(t))$ measured from the fixed shaft
center, the nominal film thickness becomes:

$$h(\theta,t) = c + x_b(t)\cos\theta + y_b(t)\sin\theta$$

with optional existing misalignment terms retained as additive axial tilt:

$$h(\theta,z,t) = c + x_b(t)\cos\theta + y_b(t)\sin\theta
                 - (z-L/2)\alpha_x\cos\theta
                 - (z-L/2)\alpha_y\sin\theta$$

When `bearing_initial_from_attitude = true`, the initial moving-bearing
position is $x_b=-e\cos\psi,\ y_b=-e\sin\psi$, so the moving model starts from
the same gap as the static shaft-eccentricity model. The equivalent attitude
angle is updated from the bearing state as $\psi=\operatorname{atan2}(-y_b,-x_b)$.

The time derivative needed by the Reynolds equation is:

$$\frac{\partial h}{\partial t} =
    \dot{x}_b\cos\theta + \dot{y}_b\sin\theta$$

**Code changes:**

| File | Change |
|------|--------|
| `config.hpp` | Added `motion_model`, bearing initial state, support stiffness/damping, bearing mass, independent external load, safety limits, and time-method keys. |
| `journal_motion.hpp/cpp` | Added bearing state, force-balance integrators, attitude conversion, and uniform output field helpers. |
| `film_thickness.hpp/cpp` | Added moving-bearing film-thickness update and `dh_dt`. |
| `main.cpp` | Stores bearing state over time; updates `h` and `dh_dt` before each Reynolds solve; advances bearing motion from fluid plus external force. |
| `io.cpp` | Uses existing field-selection path for bearing position, external load, and result-force scalar fields. |
| `gui_win32.cpp` | Added clean Motion / Loads controls for the moving-bearing and time-method config keys, plus magnitude/direction input for the independent external load. |

**Validation:**
- Static limit: zero velocity and fixed journal state must reproduce the
  current static-eccentricity pressure, load, and torque.
- Initial attitude: moving-bearing initialization from `e` and
  `attitude_angle_deg` must exactly reproduce the static film thickness.
- Safety: runs must stop or reject input before $h_{min} \le 0$.

---

### C.2 - Transient Reynolds Coupling

**Goal:** Activate the transient squeeze-film contribution created by moving
geometry.

For the Elrod variable, the transient term expands as:

$$\frac{\partial(\rho_0\theta h)}{\partial t}
  = \rho_0 h\frac{\partial\theta}{\partial t}
  + \rho_0\theta\frac{\partial h}{\partial t}$$

`FVM::ddt_weighted(sys, theta, rho*h, dt, mesh)` handles the first term for
all transient moving-bearing pressure time methods. This capacity term is
mandatory for liquid-content storage; even `EULER_EXPLICIT` keeps it in the
linear solve. The second term is assembled as an explicit source using the old
or Picard-lagged value of $\theta$:

$$S_h = \rho_0\theta^n\frac{\partial h}{\partial t}$$

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Added optional transient assembly for dynamic `h`, including the explicit $\theta\,\partial h/\partial t$ source for Elrod and an explicit squeeze source for Gumbel. |
| `field.hpp/cpp` | `dh_dt` is stored as a regular scalar field with ghost updates through existing infrastructure. |
| Tests | Added `test_journal_motion` for dynamic film-thickness sign, attitude recovery, config aliases, and integrators. |

Follow-up correction (2026-06-07): `pressure_time_method = EULER_EXPLICIT` no
longer skips the Elrod weighted capacity term. It only leaves the
$\theta\,\partial h/\partial t$ squeeze source lagged. This prevents
post-startup artificial cavitation caused by applying squeeze without the
matching old liquid-content storage.

Startup correction (2026-06-07): the earlier forced steady `t=0` transient
startup was removed. Transient runs now begin from the initialized uniform
outlet-pressure/film-content fields and use transient terms from the first
solve. `omega_ramp_time` provides the intended no-rotation startup path by
ramping `omega` linearly from zero to the target speed over the configured
timespan. The Windows GUI exposes this ramp plus advanced fluid-property and
gaseous-cavitation controls.

**Validation:**
- Squeeze-film check: with $\omega = 0$ and prescribed radial approach,
  pressure must be generated only by $\partial h/\partial t$.
- Time-step refinement: transient load should converge as `dt` is reduced.

---

### C.3 - Bearing Equation of Motion

**Goal:** Couple pressure and viscous fluid force back to bearing motion.

The first implementation uses a three-degree-of-freedom rigid bearing support:

$$m_b\ddot{x}_b + c_x\dot{x}_b + k_x x_b =
  F_x^{fluid} + F_x^{external}$$

$$m_b\ddot{y}_b + c_y\dot{y}_b + k_y y_b =
  F_y^{fluid} + F_y^{external}$$

$$m_b\ddot{z}_b + c_z\dot{z}_b + k_z z_b =
  F_z^{fluid} + F_z^{external}$$

where $\mathbf{F}^{fluid}$ is the equal-and-opposite force to the integrated
pressure plus viscous force on the shaft surface. It is saved as
`fluid_force_x`, `fluid_force_y`, and `fluid_force_z`. The independent applied
load is saved separately as `external_load_x`, `external_load_y`, and
`external_load_z`.

**Implementation strategy:**
1. Use lagged fluid force from the latest Reynolds solve.
2. Support `EULER_EXPLICIT`, `EULER_IMPLICIT`, `CRANK_NICOLSON`, `RK2`, and
   `RK4` for motion. `CRANK_NICOLSON` is the formal name for the mixed
   implicit/explicit trapezoidal method; `SEMI_IMPLICIT` and `IMEX` are accepted
   aliases in config parsing.
3. Add a Picard loop around bearing motion and Reynolds solve when lagging is
   too weak for stiff supports or small clearances.

**Code changes:**

| File | Change |
|------|--------|
| `src/journal_motion.hpp/cpp` | New module for bearing state, force balance, and time integration. |
| `main.cpp` | Added motion/Reynolds coupling sequence. |
| `config.hpp` | Added motion-model selection: `STATIC`, `MOVING_BEARING`. |
| Tests | Added static-equivalent geometry, config alias, and constant-load integrator checks. |

**Validation:**
- Static equilibrium: with large support stiffness, the dynamic solve should
  converge to the configured static eccentricity.
- Free response: with hydrodynamic force disabled, the journal ODE must match
  the analytic damped oscillator.
- Linearized bearing coefficients: small perturbations about equilibrium should
  produce repeatable stiffness/damping estimates.

---

## Phase D - Thermal and ETHD Implementation

Add temperature, viscosity variation, heat generation, heat transfer, and elastic
film-thickness deformation. This phase should be introduced incrementally:
thermal hydrodynamics (THD) first, then elastic deformation, then fully coupled
elasto-thermo-hydrodynamics (ETHD).

### D.1 - Energy Equation and Temperature Field

**Goal:** Add a lubricant temperature field and solve a thin-film energy balance.

**Status:** First constant-property implementation is complete. Remaining D.1
work is validation breadth and optional higher-order transient thermal schemes.

A first-order averaged film model is:

$$\rho c_p h\left(\frac{\partial T}{\partial t}
  + \bar{u}_\theta\frac{1}{R}\frac{\partial T}{\partial\theta}
  + \bar{u}_z\frac{\partial T}{\partial z}\right)
  =
  \frac{1}{R^2}\frac{\partial}{\partial\theta}
  \left(k h\frac{\partial T}{\partial\theta}\right)
  + \frac{\partial}{\partial z}
  \left(k h\frac{\partial T}{\partial z}\right)
  + \Phi - Q_w$$

where $\Phi$ is viscous dissipation and $Q_w$ is heat loss to journal and
bearing walls. A practical initial dissipation model is:

$$\Phi \approx \frac{\mu U^2}{h}
       + \frac{h^3}{12\mu}\left[
         \left(\frac{1}{R}\frac{\partial p}{\partial\theta}\right)^2
         + \left(\frac{\partial p}{\partial z}\right)^2\right]$$

**Code changes:**

| File | Change |
|------|--------|
| `src/energy.hpp/cpp` | Added constant-property energy-equation module using existing FVM diffusion/convection operators. |
| `main.cpp` | Added `temperature` and `heat_generation` fields and calls the energy solve after velocity/force updates. |
| `config.hpp` | Added `solution_mode`, `temperature_model`, `temperature_initial`, `temperature_reference`, `rho_cp`, `thermal_conductivity`, wall heat-transfer coefficients, and wall temperatures. |
| `gui_win32.cpp` | Added Energy controls, steady-state solution selection, and preview units for thermal fields. |
| Tests | Added `test_energy` for config aliases, isothermal heat refresh, and steady Couette wall balance. |

**Validation:**
- Isothermal limit: disabling energy coupling must reproduce Phase A/B results.
- Pure diffusion: no velocity and no heat source must match a linear or
  manufactured temperature solution.
- Energy balance: generated heat minus wall losses must match the change in
  integrated thermal energy.

---

### D.2 - Thermoviscous Reynolds Coupling

**Goal:** Feed temperature-dependent viscosity and density back into Reynolds.

Start with a configurable viscosity law:

$$\mu(T,p) = \mu_{ref}
  \exp[-a_T(T - T_{ref}) + a_p(p - p_{ref})]$$

The Reynolds diffusion coefficient then becomes cell-wise:

$$\gamma = \frac{\rho(T,p)h^3}{12\mu(T,p)}$$

and the Elrod coefficient becomes:

$$\Gamma_g = \frac{\rho\beta h^3}{12\mu(T,p)}g(\theta_{film})$$

**Code changes:**

| File | Change |
|------|--------|
| `src/fluid_properties.hpp/cpp` | Extend the planned mixture-property module with viscosity-temperature and optional Barus pressure-viscosity laws. |
| `reynolds.cpp` | Accept field-valued `rho` and `mu` for both Gumbel and Elrod paths. |
| `main.cpp` | Add THD Picard loop: update properties, solve Reynolds, solve energy, check convergence. |

**Validation:**
- Constant-property limit recovers the current solver exactly.
- Hotter lubricant should reduce viscosity and lower generated pressure for
  the same geometry and speed.

---

### D.3 - Elastic Film Deformation

**Goal:** Add elastic and thermal deformation to the effective film thickness.

The effective film thickness is:

$$h_{eff} = h_{kinematic} + \delta_p + \delta_T$$

where $\delta_p$ is pressure-induced elastic deflection and $\delta_T$ is
thermal expansion. The first implementation should support a local linear
compliance model:

$$\delta_p(\theta,z) = C_p(\theta,z)\,[p(\theta,z)-p_{ref}]$$

and a local thermal expansion model:

$$\delta_T(\theta,z) = \alpha_s\,t_s\,[T_s(\theta,z)-T_{ref}]$$

Later versions can replace this with a convolution influence matrix or an
external structural solve without changing the Reynolds interface.

**Code changes:**

| File | Change |
|------|--------|
| `src/deformation.hpp/cpp` | New module for elastic and thermal deflection fields. |
| `film_thickness.cpp` | Combine kinematic, elastic, and thermal film-thickness contributions. |
| `config.hpp` | Add elastic compliance, thermal expansion, solid thickness, and deformation toggles. |
| Tests | Local-compliance response, thermal expansion response, and zero-compliance regression. |

**Validation:**
- Zero compliance and zero thermal expansion recover THD results.
- Positive pressure compliance must increase local film thickness where the
  configured sign convention says the bearing surface moves away from the
  journal.
- Deformation should converge under mesh refinement for smooth pressure fields.

---

### D.4 - Fully Coupled ETHD Loop

**Goal:** Couple Reynolds, energy, fluid properties, and deformation.

Recommended timestep sequence:

```
1. Update journal state and kinematic film thickness
2. Initialize Picard iteration from old p, theta, T, mu, h_eff
3. Compute deformation from pressure and temperature
4. Update effective film thickness
5. Update fluid properties from T, p, and gas state
6. Solve Reynolds / Elrod equation
7. Compute velocities and viscous heat generation
8. Solve energy equation
9. Check coupled residuals for h_eff, p/theta, T, and mu
10. Repeat 3-9 until converged or max iterations reached
11. Compute load, torque, diagnostics, and advance time
```

**Code changes:**

| File | Change |
|------|--------|
| `main.cpp` | Add coupled ETHD Picard driver with configurable convergence tolerances. |
| `config.hpp` | Add `enable_thermal`, `enable_deformation`, `max_ethd_iters`, and `ethd_tol`. |
| Tests | Decoupled-limit regression and coupled manufactured/sanity cases. |

**Validation:**
- Decoupled limit: disabling thermal and deformation coupling must recover the
  current Elrod/Phase B result.
- THD-only and EHD-only toggles must be independently usable.
- Coupled residuals should decrease monotonically for moderate load/speed test
  cases or report a clear non-convergence diagnostic.

---

## Phase Summary

| Phase | Description | Primary Variable | Key Output |
|-------|-------------|------------------|------------|
| 1-3 | Infrastructure + Gumbel Reynolds | $p$ | Pressure field, basic cavitation |
| A.1 | $\theta$-reformulation (full-film only) | $\theta$ | Equivalent to Phase 3 |
| A.2 | Elrod-Adams switch + cavitation | $\theta$ with $g(\theta)$ | Mass-conserving cavitation |
| A.3 | Robust complementarity (if needed) | $\theta$ | Improved convergence |
| B.1 | Mixture EOS | — | $\rho_{mix}$, $\mu_{mix}$ |
| B.2 | Dissolved gas transport | $c_d$ | Gas release/resorption |
| B.3 | Free gas volume fraction | $\alpha_g$ | Bubble evolution |
| B.4 | Coupled solver | $\theta, c_d, \alpha_g$ | Full two-phase bearing |
| C.1-C.3 | Moving bearing dynamics | $x_b, y_b, z_b, \dot{x}_b, \dot{y}_b, \dot{z}_b$ | Fluid force, external load, squeeze-film motion |
| D.1 | Energy equation | $T$ | Temperature and heat-generation fields |
| D.2-D.4 | Thermal coupling and ETHD solver | $\mu(T,p), h_{eff}$ | Temperature-dependent viscosity, deformation, coupled ETHD response |

---

## Phase E — Audit Remediation (2026-06-11; docs consolidated 2026-06-14)

A two-agent physics/plan review, a consolidated verdict, and a
literature-grounded audit (20+ primary sources with quoted claims) identified
model-level and numerical defects. **Their full content — the agent-executable
work packages with physics, citations, code touch-points, and test acceptance
gates — is inlined under "Phase E Work Package Specifications" below.** These
were previously split across `AUDIT.md` / `AUDIT_PLAN.md` / `REVIEW.md` /
`REVIEW_VERDICT.md`, now removed; the standalone originals are preserved in git
history at commit `ae567d6`. The `§A`/`§B`/`§C`/`§D` labels in the "Closes"
lines are the original audit's finding IDs.

WP-1/WP-2 complete what Phase B.4 originally specified (mixture density into
Reynolds + Picard coupling loop); the B.4 implementation that landed left free
gas volumetrically inert, which is the audit's headline defect. **WP-12 (new,
2026-06-14)** adds the 2-D-from-file property-table backend and the
`p_sat`-driven cavitation onset; its physics basis is
`docs/CAVITATION_OIL_REFRIGERANT.md`.

| WP | Title | Closes (AUDIT) | Status |
|----|-------|----------------|--------|
| WP-4 | Conservation/regime/convergence diagnostics + doc reconciliation | Part V.7, §B1 | **DONE** (2026-06-11) |
| WP-3 | Numerics honesty (type-differencing/TVD, dead time selectors, fixed-bearing transient capacity term) | §C1, §C5 | **DONE** (2026-06-11) |
| WP-5 | Mixture viscosity with gas-limit asymptote; solubility tables | §A3, §A4 | **DONE** (2026-06-11; measured isotherm table data still pending from user) |
| WP-2 | Outer Picard coupling; converged steady state; steady+gas guard | §C4 | **DONE** (2026-06-14); `test_coupling`; transient ≡ master, steady iterates to convergence |
| WP-1 | Couple released gas into film continuity/pressure; onset at p_sat(T,c_d) | §A1, §A2 | **DONE** (2026-06-14, Stage 1) `gas_pressure_coupling=VOID_COUPLED`; θ_full/p_void/β̄; `test_eos` + `test_gas_coupling` acceptance gates (0-D conservation, resorption, vaporous limit, coupled pressure-cap 679 kPa→p_sat). NONE ≡ master; acceptance tests caught/fixed the β̄ volume-fraction weighting and gas-free bit-exactness |
| WP-6 | θ-weighted energy capacity/advection and cavitated-zone shear | §A5 | **DONE** (2026-06-14) `cavitated_film_model=STRIATED`; Couette shear/heat-gen θ-weighted in cavitated cells; `test_coupling`; FULL_FILM ≡ master, STRIATED lowers torque |
| WP-7 | Quantitative validation campaign (Ausas/V&K/Giacopini/Grando/Ferron) | Part V | TODO (grid-convergence harness V7 seeded in `validation/grid_convergence/`) |
| WP-9 | Starved / mass-flux inlet | §B2 | TODO |
| WP-8 | Linearized K/C coefficients + whirl margin | §D1 | TODO |
| WP-10 | Strategic: FBNS/complementarity cavitation kernel | §C2, §C3 | TODO |
| WP-11 | Axial-boundary flux consistency (shared flux function; cavitated-cell reformation rate; INLET_OUTLET unification) | §B4 | **DONE** (2026-06-15) `consistent_boundary_flux` (default off); cavitated boundary cells re-flood at `cs·(p_bc−p_cav)/β`. The shared `SimulationConfig::elrod_boundary_outflow` is used by both `solve_elrod` and the diagnostics balance, so r_l stays small with the flag on (7e-8 vs 1.4e-3 before). Follow-up (2026-06-15): energy and dissolved-gas convection now carry the reformation reflood via the shared `reformation_boundary_velocity` — reflood liquid convects `temperature_reference` (energy, RESERVOIR) and `dissolved_gas_initial` at the full liquid density (gas); `gas_boundary_inflow` mirrors it; `test_reformation_boundary_reservoir_state` (test #4 of the WP-11 spec). |
| WP-12 | 2-D property tables from an editable file (P,T grid) + `p_sat`-driven cavitation onset threshold | §A3/§A4 (data), §A1 (onset) | **DONE** (2026-06-14); `test_property_tables`; R290 table smoke run verifies p_sat=1.0 MPa onset |

### Phase E implementation notes (2026-06-11, WP-4/WP-3/WP-5)

- New `src/diagnostics.hpp/cpp`: per-step global liquid/gas mass balances that
  mirror the assembled operators, with clamp/snap mass accounting
  (`cavitation_clamp_mass`, `gas_clamp_mass` fields), inlet penalty cells
  re-evaluated as physical residuals, and `results/<run>/diagnostics.csv`
  output (interval: `diagnostics_interval`). Closed/open-domain residuals sit
  at `linear_rtol` (new key); the PLAN A.2 O(1e-12) target is now actually
  checked in `tests/test_diagnostics.cpp`. `validate()` adds Taylor-number and
  p_sat saturation-context warnings (`SimulationConfig::saturation_pressure`).
- WP-3: `theta/thermal/gas_convection_scheme` keys activate the previously
  dead TVD machinery; `TYPE_DIFFERENCING` is implemented for theta; the
  Elrod capacity term is assembled for **all** transient runs (fixed-bearing
  transients were silently quasi-steady); dead
  `pressure_time_method`/`temperature_time_method` members removed (parser
  warns); the Gumbel Couette term is a face-flux divergence matched to the
  Elrod discretization (matched full-film agreement ~0.01%, was ~1.2%).
  `tests/test_schemes.cpp` + `test_schemes_mpi6` cover boundedness, order,
  matching, relaxation, MPI ghosts, and back-compat.
- WP-5: `gas_mixture_viscosity_model` (EINSTEIN_DILUTE default, DUKLER_VOID,
  MCADAMS_QUALITY, KRIEGER_DOUGHERTY, LINEAR_QUALITY) + `mu_gas`; Einstein is
  rejected above the 0.6 packing bound; R290 configs ship MCADAMS_QUALITY.
  Isotherm `*_table` data remains a flagged user-supplied input.
- Test-suite adjustments forced by the capacity-gate fix: `test_elrod_a1` and
  `test_elrod_a2` now declare `solution_mode = STEADY_STATE` (their original
  quasi-steady semantics); transient film-content storage has dedicated
  coverage in `test_schemes`.

---


## Phase E — Work Package Specifications (consolidated)

## WP-12 — 2-D property tables from an editable file + `p_sat`-driven cavitation onset

**Closes:** the "isotherm table data pending" flag on WP-5, and §A1 onset
decoupling for the gas-laden case. **Severity:** medium (data fidelity) + the
user-facing "the higher-release species ruptures first" physics. **Depends on:**
none for steps 1–5; step 6 (threshold) is a default-off precursor to WP-1.
**Added** 2026-06-14 from `docs/CAVITATION_OIL_REFRIGERANT.md`.

### Physics
The supplier PTSV chart (`data/R290_PZ68S/r290_pz68s.csv` — a P×T grid of
saturated solubility and kinematic viscosity) *is* the liquid-solution property
surface. Solubility doubles as the saturation locus: inverting
$c_{sat}(p,T)=c_d$ at fixed $T$ gives the local bubble-point $p_{sat}(T,c_d)$,
the gaseous-cavitation onset. Whichever species releases first (higher pressure)
governs onset, so the effective threshold is
$p_{cav,\mathrm{eff}} = \max(p_{cav,\mathrm{oil}},\, p_{sat}(T,c_d))$ —
general, not R290-specific. For the shipped config $p_{cav}=3.0\times10^4$ Pa
but $p_{sat}\approx 1.0$ MPa, so propane outgasses first and the current Elrod
solve ignores it.

### Code (steps land in order; each keeps `ctest` green with defaults)
1. `PropertyTable2D` + bilinear `interpolate_2d` (non-uniform axes, edge-clamped) in `config.hpp`/`fluid_properties`.
2. `src/property_table_io.{hpp,cpp}`: header-mapped CSV loader; SI unit convert (MPa→Pa, °C→K, %→fraction, mm²/s→m²/s); axes reconstructed from the file (**resolution editable post-compile, no recompile**); on-load validation (headers, monotone axes, rectangular grid, finite).
3. Config keys `solubility_table_file` / `viscosity_table_file` / `density_table_file` + `config_dir` resolution; load in `validate()`; a 2-D file wins over the inline 1-D table (warn); a missing file is fatal.
4. Kinematic→dynamic viscosity per cell `μ = ν·ρ_solution` in `update_solution_fields` (ρ from the density model; PURE_OIL is a documented degraded mode; optional density file/column override).
5. Generalize `SimulationConfig::saturation_pressure` to invert the 2-D isotherm by **first-crossing scan** (the CSV is real VLE data: monotone-then-flat per isotherm with rare 1-LSB blips, some isotherms pinned at 100% — a binary search is wrong); shared `invert_csat()` for the 1-D and 2-D branches.
6. `cavitation_threshold = SCALAR_PCAV | LOCAL_PSAT` (default `SCALAR_PCAV` = current behaviour). A `p_cav_eff` field is populated from `max(p_cav, p_sat(c_d,T))` and threaded through `reynolds.cpp` via a `cav_pressure()` helper at the easy `p_cav` sites (`:114,:189-190,:233,:239,:485`); the force gauge datum (`:601`) stays the **scalar** `p_cav` (onset ≠ force reference); `p_cav_eff` is ghost-exchanged. Reconciles with WP-1: `p_void = max(p_cav_eff, m_g R_g T/((1-θ)h))`, so `p_cav_eff` is the gas-free floor of WP-1's plateau, not a competing threshold.

### Tests
New `tests/test_property_tables.cpp` (interp; CSV load incl. the README
spot-check (1.0 MPa, 40 °C) → 21.746 %; malformed-file errors). Extend
`test_fluid_properties` (2-D solubility/viscosity, ν→μ, 2-D `p_sat` round-trip +
flat-tail/pinned-isotherm/noise edges), `test_io` (file keys/precedence),
`test_forces` (`SCALAR_PCAV`≡master to 1e-10; `LOCAL_PSAT` onset≈`p_sat`; `-n 6`
rank parity). **Decoupled-limit:** defaults reproduce master bit-for-bit.

## Ground rules for every agent (read first)

1. The audit *why* (literature, defect IDs §A/§B/§C) is summarized in `PHYSICS.md` and preserved in git history at `ae567d6`; this section is the *how*.
2. Honor `CLAUDE.md` conventions: snake_case, minimal comments (headers
   self-documenting), modular/swappable models, **new physics must be
   disablable at the highest level** (config enum/flag, default = old
   behavior), `std::filesystem::path`, GCC/Clang/MSVC-buildable.
3. **Decoupled-limit regressions are sacred.** Every WP keeps the full existing
   `ctest` suite green. When you add physics, the "off" setting must reproduce
   the previous results bit-for-bit (or to solver tolerance — state which).
4. Test harness: add `tests/test_<name>.cpp`, register with
   `add_pancake_test(...)` in `tests/CMakeLists.txt` (runs under `mpiexec -n 2`
   by default; add an `-n 6` variant if your change touches ghost exchange, cf.
   `test_elrod_a2_mpi6`).
5. On completion: append a dated entry to `CHANGES.md`, update the relevant
   `PHYSICS.md` section with the model equations actually implemented, mark the
   WP status in `PLAN.md` Phase E, and update `README.md` if user-facing
   options changed.
6. Notation matches `PHYSICS.md`: θ = fractional film content, β = bulk
   modulus, `p_cav` = cavitation/plateau pressure, c_d = dissolved gas mass
   fraction, α_g = free-gas volume fraction, m_g = free gas mass per unit area
   (`free_gas_mass`), ρ_g = gas density, h = gap, ω = shaft speed, R = radius.

**Dependency graph (do in this order unless parallelizing):**

```
WP-4 (diagnostics)  ──────────────┐  (no deps; gives observability for all)
WP-3 (numerics honesty) ──────────┤  (no deps)
WP-5 (property models) ──┐        │
WP-2 (outer coupling) ───┼──> WP-1 (two-phase closure) ──> WP-7 (validation, gas cases)
WP-6 (θ-weighted energy/shear)────┘        WP-2 ──> WP-8 (K/C extraction)
WP-9 (starved inlet)  (independent)
WP-10 (FBNS kernel)   (strategic; after WP-1 stabilizes)
WP-7 (validation, pure-oil cases V1–V4) can start immediately.
```

---

## WP-1 — Couple released gas into the film continuity / pressure problem

**Closes:** AUDIT §A1, §A2 (the two "contradicts-the-literature" defects).
**Severity:** MAJOR. **Depends on:** WP-2 (outer loop), WP-5 (mixture μ); WP-4
strongly recommended first (conservation logs make debugging this possible).

### Physics

**The defect.** Today `ρ = ρ_l·θ` (`fluid_properties.cpp:216`) and
`reynolds.cpp` never reads `alpha_gas`/`free_gas_mass`; the outgassing source
`Δṁ_g` updates diagnostics only. Released gas can neither open void nor relieve
pressure — the defining mechanism of gaseous cavitation
(Braun & Hannon 2010, DOI 10.1243/13506501JET772; Etsion & Ludwig 1982, ASME
J. Lubr. Technol. 104(2):157–163). Published magnitude of what this suppresses:
load −3% at ε = 0.3, **+21% at ε = 0.9** (Int. J. Refrigeration 2023,
S0140700723004577).

**Target closure (Σα = 1).** Volume fractions liquid/gas must satisfy

$$\alpha_l + \alpha_g = 1,\qquad \alpha_l \equiv \theta,\qquad
\alpha_g = \frac{m_g}{\rho_g(p,T)\,h},\qquad
\rho_g = \frac{p}{R_g T}.$$

Mixture density and the film mass per unit area become

$$\bar\rho = \rho_l\,\theta + \rho_g\,\alpha_g, \qquad
M = \bar\rho\,h = \rho_l\theta h + m_g .$$

(Homogeneous-mixture precedent: Grando, Priest & Prata 2006, Tribology Letters,
DOI 10.1007/s11249-006-9027-6 — "ρ̄ = φρ_g + (1−φ)ρ_l … enter the Reynolds
equation directly"; same closure in Int. J. Refrig. 2023 and in the
solubility-based compressible Reynolds equation of ASME J. Tribology
147(2):024101, DOI 10.1115/1.4066414.)

**Saturation pressure (gaseous onset).** Already implicit in the code's Henry
law — make it explicit and *use* it:

$$p_{sat}(T, c_d) = \frac{c_d}{H(T)},\qquad
H(T) = H_{ref}\exp\!\left[E_H\!\left(\tfrac1T - \tfrac1{T_{ref}}\right)\right]$$

(inversion of `equilibrium_dissolved_gas`, `fluid_properties.cpp:76–85`). For
the shipped R290/PZ68 config: p_sat = 0.2175/2.175e-7 = 1.0 MPa. Void onset for
gas-laden oil must track p_sat, not the vaporous `p_cav` (Grando 2006: gas
released "when saturation pressure is reached … two-phase flow thereafter";
experimental: Braun & Hendricks 1984, ASLE Trans. 27(1):1–14 — cavitation
extent grows with gas solubility; cavity contents are released gas, not vapor).

**Stage 1 — constrained-θ void coupling (incremental, keeps Elrod machinery).**
Implement behind `gas_pressure_coupling = NONE | VOID_COUPLED`
(default `NONE` = current behavior):

1. **Gas displaces liquid (void ceiling).** The full-film ceiling for θ is no
   longer 1 but

   $$\theta_{full}(i,j) = 1 - \alpha_g(i,j),$$

   used in the switch function `g(θ) = [θ ≥ θ_full]`, in the post-solve snap
   (`reynolds.cpp:333–337` snaps to θ_full, not 1.0), and in the supply/BC film
   content. A full-film cell that releases gas at p < p_sat acquires α_g > 0,
   its ceiling drops, it cavitates — **gaseous rupture at p_sat emerges
   mechanistically** instead of being imposed.
2. **Local void pressure.** The cavitated-zone plateau becomes cell-local:

   $$p_{void}(i,j) = \max\!\left(p_{cav},\; \frac{m_g R_g T}{(1-\theta)\,h}\right),$$

   i.e. the gas filling the void sets the plateau; at release equilibrium it
   approaches p_sat. Pressure recovery uses
   `p = p_void + β̄ ln(θ/θ_full)` in pressurized cells, `p = p_void` in
   cavitated cells. (Experimental support for non-fixed, near-saturation cavity
   pressures: Etsion & Ludwig 1982; Braun & Hendricks 1984 — cavity pressure
   0.137 → 0.034 MPa as speed rises.)
3. **Mixture compressibility.** In pressurized cells containing gas, the
   effective bulk modulus is the harmonic mixture (gas isothermal
   compressibility = 1/p):

   $$\frac{1}{\bar\beta} = \frac{\theta}{\beta} + \frac{\alpha_g}{p},$$

   so Γ in those cells uses β̄ — gas presence makes the film soft and caps
   pressure near p_void. This *is* the capping mechanism; do not add an ad-hoc
   clamp.
4. **Mass bookkeeping.** Liquid-solution continuity keeps θ as unknown with
   release sink, gas continuity is the existing `free_gas_mass` transport with
   the release source — both per unit area:

   $$\partial_t(\rho_l\theta h) + \nabla\!\cdot\!(\text{Couette/Poiseuille fluxes}) = -\Delta\dot m_g,
   \qquad
   \partial_t m_g + \nabla\!\cdot\!(m_g \mathbf{u}) = +\Delta\dot m_g,$$

   with the **existing** finite-rate kinetics (keep them — peer-reviewed
   precedent: Hao & Gu 2014, Tribology Int. 78:14–26, DOI
   10.1016/j.triboint.2014.04.028; experimental motivation: Etsion & Ludwig's
   p/p_s > 1 release/resorption loop):

   $$\Delta\dot m_g = k_{dg}\,(c_d - c_{eq}(p,T))\,\rho_l\,h\,\max(\theta,\theta_{min}).$$

5. **Forces/outputs** use ρ̄ and the WP-5 mixture viscosity.

**Stage 2 (strategic, pairs with WP-10).** Full homogeneous compressible
mixture EOS ρ̄(p, T, w) with equilibrium flash + kinetic relaxation, replacing
the hard switch entirely (Grando 2006; Bayada & Chupin 2013, ASME J. Tribology
135(4):041702). Do not start Stage 2 until Stage 1 is validated.

> **Design caution for the implementing agent.** Stage 1's assembly
> (θ_full ceiling + local p_void + β̄) is a synthesis consistent with the cited
> closures, not a transcription of a single published scheme. The 0-D and 1-D
> tests below are therefore *mandatory acceptance gates*, and the equilibrium
> limit must reproduce Grando's behavior (pressure plateau at p_sat; classical
> Reynolds-BC solution recovered when gas is disabled).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Add `enum class GasPressureCoupling { NONE, VOID_COUPLED }`, key `gas_pressure_coupling` (default NONE). Validation: `VOID_COUPLED` requires `GAS_CAVITATION_MIXTURE` and TRANSIENT (until WP-2's steady loop lands). Add `gas_constant_specific` if not present (R_g). |
| `fluid_properties.cpp/hpp` | Export `saturation_pressure(c_d, T, cfg)` (Henry inversion incl. TABLE path); set `rho = ρ_l θ + ρ_g α_g` when coupling on; keep `rho = ρ_l θ` when NONE. |
| `reynolds.cpp` (`solve_elrod`) | Per-cell `theta_full` field from α_g (lagged one outer iteration); switch, snap, Γ (with β̄), and Couette flux base density use ρ̄/θ_full; pressure recovery via local `p_void`. Release sink in the θ equation source. |
| `equation_of_state.hpp` | Generalize `switch_function(theta, theta_full)`, `pressure_from_theta(theta, theta_full, p_void, beta_bar)`; keep old signatures as θ_full=1, p_void=p_cav wrappers. |
| `gas_transport.cpp` | No structural change (transport already exists); ensure the release source is consistent with the θ-equation sink (one function, used by both). |
| `main.cpp` | Order inside WP-2's outer loop: properties → Reynolds(θ, p) → gas transport → exchange → properties → … (already close; the difference is iteration, WP-2). |
| `io.cpp` | Output `p_void`, `theta_full`, `p_sat` as optional fields. |

### Test cases (`tests/test_gas_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Decoupled regression**: `gas_pressure_coupling = NONE` on `config_r290_pz68_quick.txt` | Fields identical to current master to solver tolerance (≤ 1e-10 relative) |
| 2 | **Σα invariant**: random release/transport steps | `θ + α_g ≤ 1 + 1e-12` every cell, every step |
| 3 | **0-D closed-cell equilibrium**: single cell, no flux, p < p_sat initially, kinetics on | p → p_sat asymptotically; total propane mass (c_d ρ_l θ h + m_g) constant to 1e-12; equilibrium c_d = H p_sat |
| 4 | **0-D resorption**: start gassy at p > p_sat | m_g → 0, c_d → c_eq(p); no negative masses |
| 5 | **1-D pressure cap**: converging-diverging channel, supersaturated feed | In the gassy zone, p ∈ [p_void, p_sat·(1+5%)]; no cell sustains p ≪ p_sat with θ = θ_full while supersaturated |
| 6 | **Mass conservation (periodic θ, no-flux z)** | Global liquid and total-propane residual (WP-4 logger) ≤ 1e-10 per step |
| 7 | **Load-shift direction** (journal, ε = 0.9, R290 config) | Load with coupling ON > load with NONE (direction per IJR 2023 +21%; log the ratio, no hard % gate) |
| 8 | **Vaporous limit**: c_d = 0, no gas | Identical to plain Elrod (p_void ≡ p_cav, θ_full ≡ 1) |

---

## WP-2 — Outer coupling (Picard) iteration; real steady state; steady+gas guard

**Closes:** AUDIT §C4 (REVIEW N1, N7). **Severity:** MAJOR. **Depends on:** none
(but unblocks WP-1, WP-6 fidelity, WP-8).

### Physics

A segregated one-sweep "solution" is not a fixed point of the coupled system;
published transient cavitation codes iterate the coupling inside the step —
Ausas, Jai & Buscaglia 2009 (ASME J. Tribology 131(3):031702, DOI
10.1115/1.3142903) update even the journal motion "within the same relaxation
process." Fixed-point (Picard) iteration with under-relaxation:

$$\phi^{(k+1)} \leftarrow \omega_\phi\,\phi^{new} + (1-\omega_\phi)\,\phi^{(k)},
\qquad \phi \in \{\mu, \rho, T, c_d, m_g\},\ \omega_\phi \in (0,1].$$

Convergence on normalized residuals each outer iteration k:

$$R_p = \frac{\max_{ij}|p^{(k)}-p^{(k-1)}|}{p_{ref}},\quad
R_T = \frac{\max_{ij}|T^{(k)}-T^{(k-1)}|}{T_{ref}},\quad
R_c = \max_{ij}|c_d^{(k)}-c_d^{(k-1)}|,$$

stop when all < `coupling_tolerance` (default 1e-6) or `coupling_max_iters`
(default 30; warn on hitting the cap). **Report R_p in pressure space**, not
θ-space — β amplifies θ-noise by ~1e9 (AUDIT §C2/N4), so a θ-residual hides
pressure-level error.

`STEADY_STATE` = run this loop to convergence once (this is the operating
point). Steady + gas: currently `gas_dt = 0` silently freezes the gas model
(`main.cpp:243`, `fluid_properties.cpp:222`) — replace with either a hard
config error or `steady_gas_model = FROZEN | EQUILIBRIUM` where EQUILIBRIUM
sets c_d = c_eq(p,T) and α_g from the equilibrium flash each outer iteration.

### Code implementation

| File | Change |
|------|--------|
| `main.cpp` | Extract the per-step block (`main.cpp:220–260`) into `solve_coupled_step(fields, sys, mesh, cfg, bearing_state, dt)`; wrap in the outer loop with relaxation + residual log (rank-0, `Utils::log`). Steady mode: iterate to tolerance instead of `step == 0` single pass. |
| `config.hpp` | Keys: `coupling_max_iters`, `coupling_tolerance`, `coupling_relaxation` (one factor; per-field override later if needed), `steady_gas_model`. `validate()`: error on `GAS_CAVITATION_MIXTURE && STEADY_STATE && steady_gas_model == FROZEN` unless user sets an explicit `allow_frozen_gas_steady = true`. |
| `journal_motion.cpp` | (With WP-1/WP-8 in mind) optionally update bearing state inside the outer loop (Ausas-style) behind `motion_coupling = LAGGED | IN_LOOP`, default LAGGED. |
| `io.cpp`/logging | Per-step line: outer iters, final R_p/R_T/R_c, Elrod inner iters, flag flips (WP-4 shares this). |

### Test cases (`tests/test_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Constant-property limit**: `fluid_property_model = CONSTANT`, ISOTHERMAL | Converges in exactly 1 outer iteration; results identical to master |
| 2 | **Monotone convergence**: EMPIRICAL_CORRELATION μ(T) with strong dissipation | R_p, R_T decrease monotonically (allowing first-iter transient); converged < max iters |
| 3 | **Ordering independence**: permute sweep order (Reynolds↔energy) | Converged steady fields agree to `coupling_tolerance` (this is the definition of a fixed point — fails on master) |
| 4 | **Steady ≡ long transient**: steady solve vs transient run to t→∞, same config | max|p_steady − p_transient|/p_ref < 5·tolerance |
| 5 | **Guard**: GAS_CAVITATION_MIXTURE + STEADY_STATE + FROZEN | `validate()` throws/errors with actionable message |
| 6 | **Relaxation robustness**: ω = 0.5 and 1.0 on test 2 | Same converged answer |

---

## WP-3 — Numerics honesty: convection schemes, dead selectors, capacity-term gap

**Closes:** AUDIT §C1, §C5 (REVIEW N3, N6 + audit's new finding). **Severity:**
medium, cheap. **Depends on:** none.

### Physics

- First-order upwind *inside the cavitated zone* is standard — even FBNS uses
  it (Woloszynski et al. 2015, Tribology Letters, DOI 10.1007/s11249-015-0487-4:
  "in the cavitation region … a first-order upwind formula was employed").
  Upwinding *everywhere* is not: the canonical scheme is **type differencing** —
  central in the full film, upwind across/inside cavitation (Vijayaraghavan &
  Keith 1989, Tribology Trans. 32(2):225–233, DOI 10.1080/10402008908981882;
  San Andrés Notes 06). False diffusion of full upwind ≈ |u|Δx/2 (Patankar
  1980, ch. 5) lands exactly on the rupture/reformation/thermal/gas fronts.
- TVD face reconstruction (already coded, dead in `fvm.cpp:13–27,118–162`):
  $$\phi_f = \phi_U + \tfrac12\,\psi(r)\,(\phi_D-\phi_U),\quad
  r = \frac{\phi_U-\phi_{UU}}{\phi_D-\phi_U},\quad
  \psi_{vanLeer}(r)=\frac{r+|r|}{1+|r|}.$$
- Transient capacity: $\partial(\rho_0\theta h)/\partial t$ must be assembled
  for **all** transient runs; currently gated on `MOVING_BEARING`
  (`reynolds.cpp:261–262`), so a fixed-bearing transient (e.g. an
  `omega_ramp_time` startup) is quasi-steady in θ — film-history physics
  silently absent.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Key `theta_convection_scheme = UPWIND \| TVD_VANLEER \| TVD_MINMOD \| TYPE_DIFFERENCING` (default UPWIND for back-compat), similarly `thermal_convection_scheme`, `gas_convection_scheme`. **Delete** `pressure_time_method`/`temperature_time_method` members + GUI combos; parser accepts the old keys with a one-line warning ("ignored: field solves are backward-Euler"). |
| `reynolds.cpp` | Pass configured scheme to `FVM::divergence` (`:259`); TYPE_DIFFERENCING = per-face: central if `g=1` on both sides, upwind otherwise. Change the transient gate at `:261` to `cfg.solution_mode != STEADY_STATE` for the `ddt_weighted` capacity; squeeze source remains conditional on `dh_dt` presence. Align the Gumbel Couette discretization (`:146–147`) with the Elrod path's scheme so the A.1↔Gumbel limit compares like-for-like. |
| `energy.cpp`, `gas_transport.cpp` | Pass configured schemes (`energy.cpp:344`, `gas_transport.cpp:241,260`). |
| `gui_win32.cpp` | Remove dead time-method combos; add scheme dropdowns. |
| `fvm.cpp` | None expected — limiters/deferred correction exist; verify ghost-width (n_ghost = 2) suffices for the UU stencil at partition boundaries (it does; add the MPI test below anyway). |

### Test cases (extend `tests/test_fvm.cpp`, new `tests/test_schemes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **TVD boundedness**: 1-D step advection of θ | No new extrema (θ ∈ [min, max] of initial data); sharper front than upwind (L1 error strictly smaller) |
| 2 | **Order of accuracy**: smooth 1-D profile | Upwind ~1st order, TVD ≥ ~1.5 observed order on refinement |
| 3 | **A.1↔Gumbel match**: matched schemes both paths | Residual strictly below the current 1.2% (record value in test) |
| 4 | **Capacity-term fix**: fixed bearing, ω step at t=0 | θ relaxes over finite time (exponential-like), not instantaneously; steady final state matches steady solve |
| 5 | **Regression**: all schemes default UPWIND | Existing ctest suite untouched |
| 6 | **MPI**: TVD run, 2 vs 6 ranks | Rank-count independence of fields to 1e-12 (cf. `test_elrod_a2_mpi6`) |
| 7 | **Config back-compat**: old config with `pressure_time_method = CRANK_NICOLSON` | Parses, warns, runs |

Also deliverable: `validation/grid_convergence/` script running 60×20 / 120×40 /
240×80 on the journal case, reporting cavitation extent, h_min, load (archived
table — AUDIT §C1 requirement; feeds WP-7).

---

## WP-4 — Conservation, regime, and convergence diagnostics

**Closes:** AUDIT Part V item 7, §B1 logging (REVIEW N9, M3, N8-logging).
**Severity:** low effort / high value. **Depends on:** none. **Do first.**

### Physics

Per-step global balances on the physical domain (Ω = bearing surface,
dA = R dθ dz):

$$M_l = \int_\Omega \rho_l\,\theta\,h\,dA,\qquad
M_{gas,tot} = \int_\Omega \big(c_d\,\rho_l\,\theta\,h + m_g\big)\,dA,$$

$$r_l^{n+1} = \frac{M_l^{n+1}-M_l^{n} - \Delta t\,(\Phi_{in}-\Phi_{out}+S_{exch})}{\max(M_l^{n},\epsilon)},$$

likewise r_gas (exchange terms cancel in M_gas,tot — only boundary fluxes
remain). Boundary fluxes from the assembled face fluxes (Couette + Poiseuille at
z-boundaries and inlet penalties), not re-derived. Clamp/snap operations
(`reynolds.cpp:324–337`, gas clamps) must be **accounted as explicit source
terms** in the balance so their mass effect is visible (REVIEW N9: clamps are
not mass-neutral). Target per PLAN A.2: O(1e-12) in closed domains.

Regime guards at startup (config validate + rank-0 log):

- Couette Reynolds number and Taylor criterion (laminar validity):
  $$Re_c = \frac{\rho\,\omega R\,c}{\mu},\qquad
  Ta = Re_c\sqrt{c/R} \lesssim 41,$$
  warn above (Taylor-vortex onset; standard journal-bearing criterion, cf.
  Hamrock, *Fundamentals of Fluid Film Lubrication*; San Andrés Notes).
- Saturation context (AUDIT §B1): log `p_sat(T_init, c_d_init)`, ratio
  `bc_z_*_val / p_sat`, and `p_cav / p_sat`; warn when boundary oil enters
  saturated (ratio ≥ 1) so the current R290 configuration is visible.
- Per-step: Elrod outer iterations, flag flips, (after WP-2) coupling residuals.

### Code implementation

| File | Change |
|------|--------|
| `src/diagnostics.hpp/cpp` (new) | `Diagnostics::mass_balance(fields, mesh, cfg, dt)` returning {M_l, M_gas, r_l, r_gas, clamp_mass}; MPI_Allreduce inside; rank-0 CSV append (`results/<run>/diagnostics.csv`) + log line. |
| `reynolds.cpp` | Accumulate clamp/snap mass deltas into a fields-resident scalar (or returned struct) for the balance. |
| `config.hpp::validate()` | Re/Ta computation + warnings; p_sat ratio warnings; key `diagnostics_interval` (default 1 = every step). |
| `main.cpp` | Call after each step; write CSV. |

### Test cases (`tests/test_diagnostics.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Closed-domain conservation**: periodic θ, zero-gradient z, no inlets, transient | \|r_l\| ≤ 1e-12 per step (the PLAN goal, now actually checked) |
| 2 | **Open-domain accounting**: Dirichlet z ends | Balance closes: residual after flux accounting ≤ 1e-10 |
| 3 | **Clamp visibility**: initialize θ < θ_min somewhere | clamp_mass ≠ 0 reported, balance still closes |
| 4 | **Gas total**: WP-1 test 3 config | r_gas ≤ 1e-12 (exchange cancels) |
| 5 | **Re guard**: synthetic config with μ tiny | validate() warns (Ta > 41) |
| 6 | **Saturation log**: R290 config | startup log contains p_sat = 1.0e6 and ratio = 1.0 warning |

---

## WP-5 — Property models: mixture viscosity with gas-limit asymptote; solubility tables

**Closes:** AUDIT §A3, §A4 (REVIEW W3, W4). **Severity:** medium. **Depends
on:** none (WP-1 consumes the result).

### Physics

**Free-gas (bubbly/foam) viscosity.** Einstein μ = μ_l(1 + 2.5α_g) is a dilute
result (α ≲ 0.1) and is directionally wrong as α_g → 1 (gas must *thin* the
film). Replace with selectable models, all satisfying μ̄(0) = μ_l:

- `DUKLER_VOID` (linear void-weighted; the Grando 2006 density-side analogue):
  $$\bar\mu = \alpha_g\,\mu_g + (1-\alpha_g)\,\mu_l \;\to\; \mu_g \text{ as } \alpha_g\to1.$$
- `MCADAMS_QUALITY` (homogeneous two-phase standard; quality
  x_g = m_g/(m_g + ρ_lθh)):
  $$\frac{1}{\bar\mu} = \frac{x_g}{\mu_g} + \frac{1-x_g}{\mu_l}.$$
- `EINSTEIN_DILUTE` (keep, for back-compat) — but `validate()` errors when
  selected with `gas_alpha_max > 0.6` (Krieger–Dougherty packing bound):
  $$\bar\mu = \mu_l\,(1-\alpha_g/\alpha_{max})^{-2.5\,\alpha_{max}}\ \text{(KD, optional fourth model)}.$$

(Grando 2006 quote: "μ̄ = χμ_g + (1−χ)μ_l" — quality-weighted linear; offer it
as `LINEAR_QUALITY` if exact Grando replication is wanted for WP-7 V5.)

**Dissolved-gas (solution) viscosity.** The validated form for R290/oil pairs
is the Eyring-MTSM activity model (Int. J. Refrigeration 2023,
S0140700723002451: "viscosities of R-290/oil mixtures could be well represented
by the Eyring-MTSM model"; measured 303–348 K on PAG 68 / POE 75 / PVE 68).
Do **not** implement the activity model — feed measured isotherms through the
existing `TABLE` path (`viscosity_table`, `solubility_table`,
`density_table`), and document linear Henry + `exp(a_c c_d)` as a local tangent
valid near the (10 bar, 40 °C) calibration point only. Grando's reference case
(2 bar, w_sat = 7.13 %) shows how far the shipped 21.75 % case extrapolates.

> **Data task (flag to user, do not fabricate):** populating the R290/PZ68
> tables requires the actual isotherm data (e.g. the IJR 2023 dataset or
> supplier data for PZ68). Ship the *plumbing* + a clearly-labeled placeholder
> table; mark the case "table data pending" in config comments and README.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `enum class GasMixtureViscosityModel { EINSTEIN_DILUTE, DUKLER_VOID, MCADAMS_QUALITY, KRIEGER_DOUGHERTY, LINEAR_QUALITY }`, key `gas_mixture_viscosity_model` (default EINSTEIN_DILUTE = current behavior); `mu_gas` config key; validation pairing with `gas_alpha_max` as above. |
| `fluid_properties.cpp` | Replace hardcoded Einstein at `:202–207` with the dispatch; α_g clamp from `gas_alpha_max`, not the hardcoded 0.99. |
| `config_r290_pz68*.txt` | Set `gas_alpha_max = 0.6`, select `MCADAMS_QUALITY`; comment the validity statement. |

### Test cases (extend `tests/test_fluid_properties.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Asymptotes**: each model | μ̄(α=0) = μ_l exactly; DUKLER/MCADAMS/KD: μ̄ → μ_g (or finite KD limit) as α→α_max; monotone decreasing in α (except EINSTEIN which is documented-increasing) |
| 2 | **Einstein guard**: EINSTEIN + α_max = 1.0 | validate() error |
| 3 | **Regression**: default EINSTEIN, α_max as shipped | Existing test_fluid_properties values unchanged |
| 4 | **Table round-trip**: synthetic 3-point viscosity/solubility tables | Interpolation exact at nodes, linear between, clamped outside |
| 5 | **Calibration point**: R290 case (existing test) | Still passes with new defaults documented |

---

## WP-6 — θ-weighted energy capacity, enthalpy flux, and cavitated-zone shear

**Closes:** AUDIT §A5 inconsistency (REVIEW N5, part of W6). **Severity:**
medium. **Depends on:** none (interacts with WP-2 for converged THD).

### Physics

In a striated cavitated zone only the liquid fraction θ carries thermal mass,
enthalpy, and (predominantly) shear. Replace, behind
`cavitated_film_model = FULL_FILM | STRIATED` (default FULL_FILM = current):

$$\text{capacity: } \rho c_p\,\theta\,h\,\frac{\partial T}{\partial t},\qquad
\text{flux: } \rho c_p\,\theta_f\,h_f\,\bar u\,T,\qquad
\text{Couette shear: } \tau_\theta = \theta\,\frac{\mu\,\omega R}{h}\ (\text{cavitated cells}),$$

with full-film cells (θ = 1) unchanged. This restores internal consistency: the
gas transport already θ-weights its capacity (`gas_transport.cpp:31,193`).
Friction torque/power-loss are first-order THD deliverables (groove-design
practice reports power loss + h_min: ASME J. Tribology 144(12):121801, 2022) —
the FULL_FILM model overstates them by up to ~1/θ where θ → θ_min = 1e-6.
Cross-film profile and conjugate conduction remain Phase D scope; record that
the validated 1983–84 baseline includes them (Ferron 1983, ASME 105(3):422;
Lund & Tonnesen 1984, ASME 106(2):237–244 — 4th-order cross-film polynomial +
sleeve conduction) so temperature claims stay gated on WP-7 V6.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `cavitated_film_model` enum + key; validation: STRIATED requires ELROD_ADAMS. |
| `energy.cpp` | Capacity `:301` ×θ; flux loops `:308–340` ×θ_face (face value: upwind θ consistent with flux direction); wall-loss term: decide and document whether Q_w scales with θ (liquid contact) — recommended ×θ for consistency. |
| `reynolds.cpp` | `:519` and `:543`: Couette term ×(g_f + (1−g_f)·θ) i.e. full shear in full film, θ-scaled in cavitated cells. |

### Test cases (extend `tests/test_energy.cpp`, `tests/test_forces.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Full-film limit**: θ ≡ 1, STRIATED | Bit-identical to FULL_FILM |
| 2 | **Scaling**: uniform cavitated strip θ = 0.5 | Friction torque contribution of that strip = 0.5× FULL_FILM value; energy capacity likewise |
| 3 | **Energy balance**: WP-4 style closed thermal balance with θ-weighting | Heat generated − wall loss = dE/dt to 1e-10 |
| 4 | **Direction**: journal case with cavitation, STRIATED vs FULL_FILM | Power loss strictly lower with STRIATED; T_max not higher |
| 5 | **Regression**: default FULL_FILM | Suite unchanged |

---

## WP-7 — Validation campaign (archived, quantitative)

**Closes:** AUDIT Part V (REVIEW W11). **Depends on:** V1–V4 none; V5 needs
WP-1; V6 needs WP-2+WP-6.

### Structure

`validation/<case_id>/` containing `config*.txt`, `reference/` (digitized CSV +
provenance note: paper, figure/table number, digitization method), `run.ps1` +
`run.sh`, `compare.py` (tolerance-gated, exit code), `README.md`. Each case also
gets a coarse "smoke" variant wired into ctest (`add_test(NAME validation_<id>_smoke ...)`).

### Cases, references, acceptance

| ID | Case | Reference | Acceptance gate |
|----|------|-----------|-----------------|
| V1 | 1-D JFO slider (rupture + reformation), analytic | Elrod & Adams 1975; PLAN A.2 derivation | Rupture/reformation locations within 1 cell at 200 cells; p_peak ≤ 1 %; θ-profile L2 ≤ 1 % |
| V2 | Oscillatory squeeze film (exact solution) + dynamically loaded JB | Ausas, Jai & Buscaglia 2009, ASME JT 131(3):031702 (paper ships reference source code — use it for code-to-code) | p(t) L2 ≤ 2 % vs exact; JB orbit qualitative match + mass residual ≤ 1e-10 |
| V3 | Grooved journal bearing, flooded + (post-WP-9) starved | Vijayaraghavan & Keith 1989: Trib. Trans. 32(2):225–233 (algorithm cases, code-to-code) and Wear 134(2):377–397 (grooved, vs experiment) | Cavitation extent ±1 cell @ published grid; load ±5 % |
| V4 | LCP analytic suite (pocketed slider, squeeze damper) | Giacopini et al. 2010, ASME JT 132(4):041702 | p-profile L2 ≤ 2 % on the 1-D cases |
| V5 | Refrigerant equilibrium limit | Grando, Priest & Prata 2006 (Trib. Lett.) | With gas coupling ON and equilibrium kinetics (k_dg → large): recover classical Reynolds-BC solution for the moderately loaded case (their validation); plateau at p_sat |
| V6 | THD test bearing | Ferron, Frêne & Boncompain 1983, ASME 105(3):422; Lund & Tonnesen 1984, ASME 106(2):237 | Phase 1: qualitative T(θ) shape + documented gaps (no conjugate conduction). Hard gates deferred to Phase D; record deltas. |
| V7 | Grid convergence + conservation history | WP-3/WP-4 outputs | Monotone convergence of extent/load/h_min; archived plots |

**Citation hygiene for the implementing agent:** the V&K 1989 Trib. Trans.
paper is code-to-code only (its own abstract: comparisons "with the results
obtained using Elrod's algorithm"); the experimental comparison is the Wear
companion. Do not cite the former as experimental validation.

---

## WP-8 — Linearized dynamic coefficients (K_ij, C_ij) and whirl margin

**Closes:** AUDIT §D1 (REVIEW W9b). **Depends on:** WP-2 (needs a converged
operating point).

### Physics

Small-perturbation method (Lund 1987, ASME J. Tribology 109(1):37–41, DOI
10.1115/1.3261324 — the canonical reference; applied to gas-oil refrigerant
bearings via the perturbed compressible Reynolds equation in ASME JT
147(2):024101, DOI 10.1115/1.4066414):

$$K_{ij} = -\left.\frac{\partial F_i}{\partial x_j}\right|_{eq}
\approx -\frac{F_i(x_j^{eq}+\Delta x_j)-F_i(x_j^{eq}-\Delta x_j)}{2\,\Delta x_j},
\qquad
C_{ij} = -\left.\frac{\partial F_i}{\partial \dot x_j}\right|_{eq},$$

i, j ∈ {x, y}; force = fluid force on the journal/bearing per the existing
convention. Perturbation sizes: Δx = ε_K·c with ε_K ~ 1e-2 (linearity check:
halving Δx changes K by < 1 %); velocity perturbations Δẋ = ε_C·c·ω applied as
prescribed `dh/dt` with frozen position. Each force evaluation is a **converged
WP-2 operating point** (transient effects off; squeeze term from the prescribed
ẋ only).

Stability deliverable: for a rigid rotor of mass m per bearing, the
characteristic equation `det(m s² I + C s + K) = 0`; report eigenvalues, onset
speed (Re(s) = 0), and whirl frequency ratio Im(s)/ω (≈ 0.5 for plain bearings
— sanity anchor; San Andrés Notes 05). Lund's validity caveat (small motion
about equilibrium) goes in the output header.

### Code implementation

| File | Change |
|------|--------|
| `src/dynamic_coefficients.hpp/cpp` (new) | `compute_dynamic_coefficients(fields, sys, mesh, cfg, state)` → struct {K[2][2], C[2][2], whirl_ratio, onset_mass_or_speed}; loops ±perturbations calling the WP-2 steady driver. |
| `config.hpp` | `compute_dynamic_coefficients` (bool), `perturbation_position_fraction`, `perturbation_velocity_fraction`. |
| `main.cpp` | Post-run hook when enabled; CSV + log output. |
| `io.cpp` | Append K/C block to a summary file. |

### Test cases (`tests/test_dynamic_coefficients.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Linearity**: Δx and Δx/2 | K_ij relative change < 1 % |
| 2 | **Sign structure** (plain journal, moderate ε) | K_xx, K_yy > 0; K_xy·K_yx < 0 (cross-coupling signs per Lund 1987 / San Andrés Notes 05) |
| 3 | **Symmetry of C** (aligned, no gas) | C_xy ≈ C_yx within 5 % |
| 4 | **Consistency vs transient**: small free decay orbit vs linear model prediction | Frequency/decay within 10 % |
| 5 | **Whirl anchor**: light-load plain bearing | whirl ratio ∈ [0.45, 0.55] |

---

## WP-9 — Starved / mass-flux inlet option

**Closes:** AUDIT §B2 (REVIEW W8). **Depends on:** none.

### Physics

Starved feed = the supply cannot maintain a full film at the inlet; prescribe
inlet film fraction θ_in < 1 (or feed flow). JFO reformation downstream is
already handled by Elrod machinery. Precedent: starved cases inside the Elrod
framework since Vijayaraghavan & Keith 1989 (Wear 134(2):377–397 — "both
flooded and starved inlet conditions … various degrees of starvation"); groove
geometry as a design input per ASME JT 144(12):121801 (2022).

Two modes on `InletConfig`:

- `STARVED_THETA`: penalty-pin θ = θ_in (θ_in ∈ (θ_min, 1)); cavitated inlet
  cell (g = 0), Couette flux carries ρ_l θ_in h u into the domain.
- `MASS_FLUX`: prescribe feed ṁ [kg/s]; convert to equivalent θ_in over the
  inlet footprint from the local Couette flux:
  θ_in = ṁ / (ρ_l (ωR/2) h L_z,inlet), clamped to (θ_min, θ_supply(p)).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `InletConfig::Type` extends with `STARVED_THETA`, `MASS_FLUX` (parse `inlet_N_mode`, `inlet_N_theta`, `inlet_N_mdot`); validation: starved modes require ELROD_ADAMS. |
| `reynolds.cpp` (`apply_inlet_conditions`, `:60–91`) | Branch: pressure mode unchanged; starved modes pin θ_in via the same penalty; MASS_FLUX computes θ_in per cell from local h. |
| `energy.cpp`/`gas_transport.cpp` | Inlet thermal/composition contracts treat starved inlets like supply cells (existing reservoir-composition path). |
| `gui_win32.cpp` | Inlet mode selector + value field. |

### Test cases (`tests/test_inlets.cpp` extension)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Flooded equivalence**: θ_in = θ_supply(p_supply) | Matches pressure-inlet result to 1e-10 |
| 2 | **Monotonicity**: θ_in ∈ {1.0, 0.8, 0.5} | Load strictly decreasing; cavitation extent strictly increasing |
| 3 | **Flux fidelity** (MASS_FLUX) | Integrated inlet mass flux = prescribed ṁ within 1 % |
| 4 | **Conservation**: starved case under WP-4 logger | Balance closes ≤ 1e-10 |

---

## WP-10 — Strategic: complementarity (FBNS) cavitation kernel

**Closes:** AUDIT §C2, §C3 (β-stiffness, switch-iteration fragility; REVIEW
N4 + S3 caveat). **Depends on:** stabilize WP-1 first (formulation choice
interacts). **This retires `snap_tol` and the β-conditioning class of bugs.**

### Physics

JFO mass conservation as a complementarity problem (Giacopini et al. 2010, ASME
JT 132(4):041702, DOI 10.1115/1.4002215):

$$\bar p \ge 0,\qquad \vartheta \ge 0,\qquad \bar p\,\vartheta = 0,
\qquad \bar p = p - p_{cav},\ \ \vartheta = 1-\theta,$$

reformulated as a smooth unconstrained system via the Fischer–Burmeister
function (Woloszynski, Podsiadlo & Stachowiak 2015, Tribology Letters, DOI
10.1007/s11249-015-0487-4):

$$\phi_{FB}(\bar p_{ij}, \vartheta_{ij}) = \bar p_{ij} + \vartheta_{ij}
- \sqrt{\bar p_{ij}^2 + \vartheta_{ij}^2} = 0,$$

solved with Newton + Schur-complement reduction (FBNS) coupled to the
discretized Reynolds residual R(p, θ) = 0. Literature-reported performance: ~2
orders of magnitude faster than Elrod-Adams/LCP-pivoting baselines,
mesh-independent Newton counts (900 → 360 000 DOF), upwind θ-flux in the
cavitated region retained. Motivation for leaving the β-switch: solution
"strongly depend[s]" on β, β dominates the discrete system's stability
(eigenvalue analysis), softening tolerable only ~2 orders (3–40 % error beyond)
— ASME JT 139(3):031703 (2017); practitioner β-reduction practice per San
Andrés Notes 06. No β appears in the FBNS formulation; physical β returns only
through an optional compressible full-film EOS.

### Code implementation

| File | Change |
|------|--------|
| `src/reynolds_fbns.cpp/hpp` (new) | `Reynolds::solve_fbns(...)` parallel to `solve_elrod`; Newton loop with analytic Jacobian blocks (Reynolds w.r.t. p and θ; φ_FB diagonal blocks); reuse `LinearSystem`/PETSc (KSP inside Newton, or SNES if available). |
| `config.hpp` | `cavitation_model = GUMBEL \| ELROD_ADAMS \| FBNS`; `fbns_max_newton`, `fbns_tolerance`. |
| `main.cpp` | Dispatch. WP-1's mixture terms enter through ρ̄, p_void analogues — coordinate with whoever owns WP-1 Stage 2. |

### Test cases (`tests/test_fbns.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Equivalence**: V1 slider + journal case vs ELROD_ADAMS | p, θ agree within discretization tolerance (≤ 0.5 %) |
| 2 | **Complementarity**: converged fields | max(p̄_ij·ϑ_ij) ≤ 1e-10·p_ref |
| 3 | **Mesh-independent iterations** | Newton counts within ±2 across 60×20 → 240×80 |
| 4 | **No snap needed** | Flag-chatter test (the scenario behind `snap_tol`, reynolds.cpp:328–337 comment) converges without snapping |
| 5 | **Timing** | ≥ 5× faster than ELROD_ADAMS outer loop on the 240×80 case (log; literature says ~100×, don't gate on that) |

---

## WP-11 — Axial-boundary flux consistency (single source of truth)

**Closes:** AUDIT §B4 (line-level BC review, 2026-06-11). **Severity:** medium
(BC-3 live when cavitation reaches the ends; BC-1/2/4 latent — no shipped
config uses `INLET_OUTLET`). **Depends on:** none; coordinate with WP-1 (the
shared flux function is where `p_void` will later replace `p_cav`).

### Physics

The three solvers currently assemble the *same* axial boundary three different
ways. Required invariant: **one liquid mass flux per boundary face, used by
Elrod assembly, energy convection, gas convection, and the diagnostics
balance**, with direction-dependent temperature/composition layered on top.

Per face (south shown; north mirrored), with `p_bc = elrod_boundary_pressure(bc_val)`:

$$\dot m''_{bc} =
\begin{cases}
\dfrac{\rho_0\,\theta_P\,h^3}{12\mu}\,\dfrac{p_{bc}-p_P}{\Delta z/2} & \theta_P \ge 1 \ \text{(full film — EOS chain rule } \beta\,d\theta=\theta\,dp\text{)}\\[2ex]
\max\!\left(0,\ \dfrac{\rho_0\,h^3}{12\mu}\,\dfrac{p_{bc}-p_{cav}}{\Delta z/2}\right) & \theta_P < 1 \ \text{(cavitated — physical reformation rate)}
\end{cases}$$

The cavitated branch replaces the current ungated `Γ_base` Dirichlet link,
which overestimates re-flooding by ≈ `β(1−θ_P)/(p_bc−p_cav)` (~10²–10³×)
because the EOS `p = p_cav + β ln θ` is invalid where p ≡ p_cav. After WP-1,
`p_cav` in the cavitated branch becomes `p_void(i,j)`.

`INLET_OUTLET` semantics unified: for the **mass** equation it behaves as
DIRICHLET (bidirectional, like the Gumbel path already does); direction
dependence applies only to *what the flux carries* — inflow brings
`temperature_reference` (RESERVOIR mode) and `dissolved_gas_initial` / zero
free gas; outflow convects cell values (energy/gas already do exactly this).
The lagged `fields["pressure"] < bc_val` gate (stale by one timestep, frozen
across outer iterations) is deleted; direction comes from the sign of the
shared flux, evaluated with the current iterate.

Cavitated-cell reformation inflow must also carry reservoir enthalpy and
composition (today the g-gate zeroes both, so re-flooding liquid silently
adopts the cell's T and c_d).

> **Implemented (2026-06-15).** The shared flux lives on `SimulationConfig`, not a
> separate `boundary_flux.cpp`: `elrod_boundary_outflow` (liquid mass balance, used
> by `solve_elrod` + diagnostics) and `reformation_boundary_velocity` /
> `is_reformation_boundary` (the reflood Poiseuille velocity, used by `energy.cpp`,
> `gas_transport.cpp`, and `diagnostics.cpp`). The cavitated reflood carries
> `temperature_reference` (energy, RESERVOIR mode) and `dissolved_gas_initial` at the
> full liquid density $\rho_\ell$ (gas), resolving the g-gate-zeroes-both defect.
> Acceptance: `test_reformation_boundary_reservoir_state` (the spec's test #4) plus
> the existing `r_l ~7e-8` liquid-balance gate; the planned `test_boundary_fluxes.cpp`
> table below was folded into `test_diagnostics.cpp`/`test_coupling.cpp`.

### Code implementation

| File | Change |
|------|--------|
| `src/boundary_flux.hpp/cpp` (new) | `axial_boundary_mass_flux(fields, mesh, cfg, side, i)` implementing the two-branch formula; plus `axial_boundary_volume_flux` (free gas) and a struct return {mass_flux, is_inflow}. |
| `reynolds.cpp` (`solve_elrod`) | Full-film boundary cells: keep the Dirichlet link (equivalent). Cavitated boundary cells: drop the `Γ_base` link; add the reformation flux as an explicit θ-equation source `+ṁ''_{bc}·RΔθ/(ρ_0 h ...)` consistent with the assembled units. Delete the `INLET_OUTLET` lagged gate (`:388–414`); INLET_OUTLET assembles like DIRICHLET. |
| `energy.cpp` | `south/north_boundary_flux` call the shared function (×`c_p` face value) instead of the local g-gated expression; RESERVOIR/OPEN logic unchanged. |
| `gas_transport.cpp` | `south/north_boundary_mass_flux`/`_volume_flux` call the shared function; composition logic unchanged. |
| `diagnostics.cpp` | Balance uses the same shared function (removes any residual reimplementation). |
| `config.hpp` | Document NEUMANN = sealed/symmetry plane in the key comment + PHYSICS.md §6. |

### Test cases (`tests/test_boundary_fluxes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Triple consistency**: random full-film states | Elrod-assembled boundary flux, energy flux ÷ ρc_p-face, gas mass flux agree cell-by-cell to 1e-12 (by construction once shared) |
| 2 | **Cavitated reformation rate**: 1-D z-strip, cavitated interior, submerged end | Re-flood mass rate = `ρ₀h³(p_bc−p_cav)/(12μ·Δz/2)` analytic within 1 %; assert the old `Γ_base` link would give > 10× (regression demonstrating the fix) |
| 3 | **Outflow continuity**: interior p > p_bc with INLET_OUTLET | Liquid mass leaves; WP-4 balance closes ≤ 1e-10; energy/gas convect the same direction in the same step |
| 4 | **Reformation carries reservoir state**: cavitated end cell re-floods (RESERVOIR mode) | Cell c_d moves toward `dissolved_gas_initial`, T toward `temperature_reference` — not frozen |
| 5 | **Direction-flip robustness**: p_bc oscillating about interior pressure | Flux continuous through zero; no chatter; conservation holds |
| 6 | **Full-film regression**: shipped configs (ends full-film) | Results identical to current master to solver tolerance |
| 7 | **Gumbel parity**: INLET_OUTLET vs DIRICHLET in Gumbel | Identical (documents BC-4 resolution) |

## Phase E priority (consolidated from the review verdict)

- **Tier 1 (MAJOR).** WP-1 (two-phase closure — released gas is volumetrically inert to pressure) and WP-2 (outer Picard coupling — `STEADY_STATE` is a single segregated sweep, not a fixed point). **WP-12** supplies the `p_sat` data backend both depend on.
- **Tier 2 (medium).** WP-3 numerics (done), WP-5 property laws (done), WP-6 θ-weighted energy/shear, WP-11 boundary-flux consistency, WP-8 FSI / K-C extraction.
- **Tier 3 (low-effort honesty).** dead time selectors (done, WP-3), frozen-gas-in-steady (WP-2), config/doc drift (this consolidation), starvation option (WP-9), mass diagnostics (done, WP-4), regime guards (done, WP-4).
- **Tier 4.** WP-7 quantitative validation campaign.


## Phase F — Maintainer Handoff Remediation (2026-07-13)

`HANDOFF.md` is the retiring senior maintainer''s verified risk register. This
phase turns its D1-D16 defects and B1-B10 backlog into an execution order. It
does not replace the handoff: keep the defect rationale, invariants, and cited
code locations there, and record implementation evidence here.

### Operating rules

- Correctness and fail-loud behavior precede performance, refactoring, and new
  physics. Do not start WP-8/WP-9/WP-10 while F1-F6 are open.
- Use one `codex/<short-scope>` branch per row below. A row is DONE only after
  its focused tests, the full available platform suite, and required docs pass.
- Preserve the model-off equivalence gates from HANDOFF I6 and the conservation
  accounting from I4/I5. Numerical changes also require MPI and grid checks
  when their stencil, state, or boundary behavior is affected.
- Record commands, rank counts, residual tolerances, and any accepted numerical
  deltas in `CHANGES.md`; do not use "tests pass" without the test evidence.

### Ordered work

| ID | Work | Handoff source | Status |
|----|------|----------------|--------|
| F0 | Consolidate the authoritative worktree, establish the branch rule, and capture a clean Windows test baseline | I9, §7 | **DONE (2026-07-13)**; baseline evidence recorded with the consolidation commit |
| F1 | Make PETSc solves fail loudly: convergence reason, iteration count, finite-field guard | D6 | TODO |
| F2 | Initialize `theta` with the same effective cavitation pressure used by the EOS and solver | D3 | TODO |
| F3 | Replace every hard-coded `theta >= 1` decision with one shared per-cell full-film predicate | D2, I3 | TODO |
| F4 | Correct reversed-flow van Leer reconstruction and add north/south TVD correction | D1 | **PARTIAL (2026-08-05)** - reversed-flow reconstruction fixed and pinned by the directional-symmetry gate; the north/south TVD correction is still absent |
| F5 | Synchronize all properties consumed across partitions and add coupled-path 1-rank/2-rank parity | D5, B9, I1 | TODO |
| F6 | Add a fully populated config parse/serialize round-trip gate before extending the config surface | I7, B1a | TODO |
| F7 | Harden linear/nonlinear iteration: warm starts, `theta`-space relaxation, tolerance-scaled snapping, exact constraint rows | D7, D8, D10, D11 | TODO |
| F8 | Reconcile model and documentation consistency: force datum, Gumbel path, striated thermal transport, transient density storage, dead code | D9, D13-D16 | TODO |
| F9 | Implement ghost-exchanged per-cell `p_cav_eff(T,c_d)` and retain a deliberate scalar force datum | D4, D9 | TODO |
| F10 | Establish structural seams: `Simulation` step API, typed field names, split config implementation, generated config descriptors | B1, B3-B5 | TODO |
| F11 | Add restart and timestep controls after the step API exists: checkpoint/restart, adaptive `dt`, optional transient Picard | B2, B6, B7 | TODO |
| F12 | Decide and document the linear-algebra strategy: bless PETSc or hide it behind a backend boundary before any Trilinos work | B10 | TODO — owner decision required |
| F13 | Run the quantitative validation and R290 transfer-rate sensitivity campaigns | WP-7, B8 | TODO |
| F14 | Resume new capability work in dependency order: mass-flux inlet, K/C and whirl margin, then complementarity cavitation | WP-9, WP-8, WP-10 | TODO |

### Acceptance gates for F1-F6

1. **F1 — solver status:** force a non-converged KSP case and assert a nonzero,
   reason-bearing failure; normal solves report iterations and remain unchanged
   within tolerance. NaN/Inf injection must stop before output is written.
2. **F2 — EOS triangle:** for `SCALAR_PCAV` and `LOCAL_PSAT`, step-zero
   `pressure -> theta -> pressure` round-trips to solver tolerance. The shipped
   R290 case must not carry the previous approximately 0.97 MPa inconsistency.
3. **F3 — one film classifier:** no consumer-owned `theta >= 1.0` decisions
   remain in Reynolds post-processing, energy, gas transport, or boundary
   logic. A nonzero-void `VOID_COUPLED` case changes continuously as
   `alpha_gas -> 0`, while `NONE` remains equivalent to the baseline.
4. **F4 — directional symmetry:** mirrored one-dimensional advection with
   positive and negative flux produces mirrored results. An axial step becomes
   sharper under VANLEER than UPWIND without violating boundedness.
5. **F5 — MPI parity:** the quick TWO_PHASE transient agrees field-by-field for
   one and two ranks at approximately `linear_rtol`; retain the six-rank Elrod
   interface test.
6. **F6 — config closure:** a config populated with every non-deprecated key
   survives parse/serialize/parse with field-wise equality, including table
   paths and GUI-only edit cycles.

### Structural and science gates

- F7 must land before transient Picard or adaptive time stepping: both depend
  on trustworthy convergence detection and relaxation of the actual Elrod
  state.
- F10''s `Simulation` API precedes checkpoint/restart. The restart acceptance
  test is `continuous run == split run + restart` within solver tolerance.
- F12 must produce an explicit architecture decision record. Do not begin a
  backend migration while PETSc-specific behavior remains exposed to solver
  assembly code.
- F13 archives grid, timestep, scheme, rank, and `k_dg` sensitivity tables.
  Only after these are credible may F14 add new physical capability.

## References

- Elrod, H.G. & Adams, M.L. (1975). A computer program for cavitation and starvation problems. *Cavitation and related phenomena in lubrication*, 37-41.
- Vijayaraghavan, D. & Keith, T.G. (1989). Development and evaluation of a cavitation algorithm. *Tribology Transactions*, 32(2), 225-233.
- Fesanghary, M. & Khonsari, M.M. (2011). A modification of the switch function in the Elrod cavitation algorithm. *ASME J. Tribology*, 133(2), 024501.
- Ausas, R.F., Jai, M. & Buscaglia, G.C. (2009). A mass-conserving algorithm for dynamical lubrication problems with cavitation. *ASME J. Tribology*, 131(3), 031702.
- Bayada, G. & Chupin, L. (2013). Compressible fluid model for hydrodynamic lubrication cavitation. *ASME J. Tribology*, 135(4), 041702.

### Audit / two-phase / refrigerant-bearing references (DOIs)

Key DOIs (from the consolidated audit): Elrod-Adams numerics — 10.1080/10402008908981882,
10.1115/1.3142903, 10.1115/1.4002215, 10.1007/s11249-015-0487-4; refrigerant
two-phase — 10.1007/s11249-006-9027-6, S0140700723004577, 10.1115/1.4066414,
10.1016/j.triboint.2014.04.028; experiments — Etsion & Ludwig (ASME JLT
104(2):157), Braun & Hendricks (10.1080/05698198408981539); THD — Ferron 1983
(ASME 105(3):422), Lund & Tonnesen (10.1115/1.3260891); dynamics — Lund 1987
(10.1115/1.3261324); review — Braun & Hannon (10.1243/13506501JET772).


## FW68D / R410A properties (2026-09-15)

- [x] Archive the user's chart calibration and all 125 decoded property points.
- [x] Add solubility interpolation in absolute pressure and temperature, with the same isotherm inverse for bubble pressure.
- [x] Add kinematic viscosity interpolation in actual refrigerant mass fraction and temperature, including c=0 and c=1.
- [x] Expose linear/log viscosity interpolation and clamp/reject range handling in config and GUI.
- [x] Add configurable constant density (940 in the FW68D example) and linear density (954.9 - 0.7284 T_C).
- [x] Convert kinematic to dynamic viscosity using the selected density and document source limits.
- [x] Add reproducible importer checks and C++ source-node, interpolation, units, density, bounds and config regression checks.

## Branch consolidation (2026-09-15)

- [x] Combine the Ferron, optical-rig and FW68D property work.
- [x] Preserve the unfinished bubble solver as a separate optional target.
- [x] Test the merged state on native Windows and workstation afml-r6; prepare
  fast-forward publication. See `docs/INTEGRATION_2026_09_15.md` for evidence.

## 2026-09-27 - FW68D comparison plan and workstation MPI timing

- [x] Define a 24-point, two-speed FW68D-only cavitation comparison matrix from
  the supplied experimental Sommerfeld ranges, with four corner pilots and
  eight initial spatial/timestep verification repeats.
- [x] Review moving-bearing, pure-oil viscosity, film energy and reduced CHT
  compatibility. Identify cross-film wall-shear and speed-ramp consistency
  work required before quantitative experimental comparison.
- [x] Prepare a synthetic geometry-matched timing fixture without changing
  numerical source; build locally and pass native 1/2-rank short checks plus
  a full 360-by-96 smoke check before the isolated workstation build.
- [x] Complete and archive initial 360-by-96 workstation scaling and throughput:
  1 rank minimizes core cost, 14 two-rank jobs lead tested queue throughput,
  and 28 ranks minimize isolated latency for this short startup fixture.
- [ ] Freeze experimental load/reference-viscosity, thermal and support inputs
  before launching the physical 24-point comparison campaign.
- [ ] Correct cross-film wall forces/torque and any speed-ramp usage, verify
  timestep sensitivity, then execute the accepted physical campaign.

See [the campaign and benchmark plan](cases/90_diagnostics/fw68d_mpi_benchmark/README.md).
The short timing fixture is not an experimental validation or a mesh/time
independence study. Earlier completed tasks remain part of this record.

- [x] Reuse the prior accepted G6 720-by-192 mesh, 20-us timestep and 16-cell
  radial bush discretization in the FW68D plan and final MPI performance check.
- [ ] Resolve the G6 synthetic startup's minimum-film guard before step 3
  (reproduced at 24/28 ranks), then obtain an admissible G6 scaling comparison.

## 2026-09-27 - Replace the transient sweep with steady FW68D cases A–F

- [x] Replace the active 24-point transient matrix with the six user-selected
  measured cases, retaining G6 and adding two Case E temperature sensitivities.
- [x] Preserve exact experiment inputs and align the pure-oil table with the
  supplied viscosity nodes; independently verify mu, S and friction targets.
- [x] Specify STEADY_STATE with static inner geometry and an outer vector
  load-balance search; dt=20 us is inactive in the physical steady equations.
- [ ] Establish feed/load-frame, cavity-pressure and fluid/solid thermal inputs.
- [ ] Correct and verify cross-film wall shear before quantitative friction.
- [ ] Build/check locally, run the force-balanced steady cases remotely, return
  provenance/results and compare paired-speed temperature/friction/cavitation.

See [the current six-case plan](cases/02d_fw68d_steady_six/README.md).
The prior startup-motion investigation remains archived and is not a gate
on static equilibrium solves; their thermal and force convergence are checked
independently.


## 2026-09-27 - Record FW68D groove and housing boundaries

- [x] Set nominal 5 mm full-length axial groove opposite applied load and
  tie its pressure/temperature to the common bearing-end reservoir.
- [x] Record OD 85 mm, L 40 mm stainless housing and natural convection to
  still oil at TC0; estimate conditional h for all six reference temperatures.
- [x] Define applied force (0,-W) and fluid-on-bearing target (0,+W) explicitly.
- [ ] Resolve common reservoir pressure and gauge/absolute meaning of p_cav=0.
- [ ] Establish stainless conductivity, shaft cooling and FW68D thermal properties.
- [ ] Couple the natural-convection estimate to thermal equilibrium, assess
  omitted housing end-face cooling, and verify the selected cooling sensitivity.

These additions resolve groove placement and housing dimensions from the prior
entry. They do not mark physical runs or thermal convergence complete.


## 2026-09-27 - Enable and validate steady FW68D execution

- [x] Apply atmospheric end/groove pressure and the requested zero cavity plateau.
- [x] Select sourced ordinary/analogous stainless and PVE thermal properties.
- [x] Represent both exposed shaft-end areas and separate rotating-end cooling
  from natural convection at the stationary housing.
- [x] Correct cross-film wall forces/torque and verify native one/four-rank limits.
- [x] Implement the force/thermal campaign driver and pass its no-CFD self-check.
- [ ] Finish G6 native/remote preflight and measure a suitable MPI rank layout.
- [ ] Obtain the six accepted steady equilibria and two Case E sensitivities,
  return evidence, and compare temperature, friction and cavity regions.

The original discretization is retained. Non-converged trials are evidence for
numerical diagnosis only and do not complete an experimental comparison.

### FW68D steady execution follow-up (2026-09-27)

- [x] Rebuild f754147 in an isolated remote snapshot; focused remote tests pass.
- [x] Identify and preserve two invalid native timing inputs with inherited
  housing ambient override; correct generated reservoir references.
- [x] Retain 1e-6 coupling tolerances, use fixed relaxation and a sufficient cap.
- [x] Complete corrected native check and physical-core steady timing test.
- [x] Accept a converged fixed-position solve with tight FBNS tolerance before the load/thermal sweep.

- [x] Detect an inner FBNS stopping floor masking the outer pressure tolerance.
- [x] Repeat the fixed-position preflight with fbns_tol=1e-12; pressure crosses tolerance smoothly at iteration 484.
- [x] Assess localized pinned-groove wall-temperature undershoot explicitly; document baseline limitation without clamping.

- [x] Launch A-F on six disjoint four-rank slots, with two E temperature variants queued.
- [ ] Return completed production equilibria, inspect model limitations and compare with experimental targets.


## 2026-09-28 - Recover physically admissible FW68D equilibria

- [x] Diagnose the failed pressure/thermal fixed points and preserve baseline runs.
- [x] Repair groove enthalpy exchange and FBNS mass/energy transport consistency.
- [x] Add conservation and internal-profile admissibility checks without clipping.
- [x] Add validated reusable steady initial fields and bounded force backtracking.
- [ ] Rebuild the isolated workstation snapshot and validate representative CHT trials.
- [ ] Close thermal coefficients and bearing position for A-F and E temperature variants.
- [ ] Return results and provenance, inspect all acceptance checks, and merge verified code.

### September 28 validation progress

- [x] Pass the corrected workstation build, 14 focused tests and explicit
  four-rank I/O/diagnostics checks.
- [x] Verify the position-search algorithm from two distinct starts with a
  rejected full step and recovery by backtracking.
- [x] Pass Case E inner conservation and physical checks at its trial position.
- [ ] Finish external cooling and force closure; an inner pass is not an
  accepted experimental case.

### September 28 accepted-case progress

- [x] Complete Case E force and thermal closure and return the full evidence.
- [ ] Complete A/C position searches; latest accepted-position errors remain
  above the 0.2% target.
- [x] Start D on E's freed eight-core slot and queue F/E_low/E_high.
- [x] Resolve B's execution blocker through approved use of the already-deployed
  solver and inputs; start B on cores 24--27 without a new launcher transfer.
- [ ] Finish all remaining cases and experimental temperature/friction/cavity
  comparisons before declaring the campaign complete or merging the task.


## 2026-09-28 - Correct the confirmed upward bearing push

The earlier completed E milestone above is superseded for experimental use:
all downward-load cases solve the wrong confirmed physical configuration.

- [x] Identify the user-confirmed upward push and invalidate prior comparisons.
- [x] Stop affected workstation drivers and preserve all wrong-direction outputs.
- [x] Correct explicit load, force residuals, initial attitudes and restart checks.
- [x] Pass an independent force-sign regression and a short native CFD check.
- [x] Deploy a fresh isolated snapshot, rebuild and verify on the workstation.
- [ ] Recompute A-F and E sensitivities with both load and thermal closure.


### Final convergence correction before the upward-load campaign

- [x] Prevent skipped steady pressure updates from satisfying outer convergence.
- [x] Pass exact/perturbed warm-state and unchanged transient checks at 1/2/4/6 ranks.
- [x] Validate the final native and Linux snapshot, then launch independent cases
  on seven disjoint four-core slots; queue the eighth temperature variant.
- [ ] Accept only the same-position load, pressure, thermal and conservation closure.


### September28 thermal recovery

- [x] Reproduce and fix the absolute-Kelvin energy accuracy floor; pass28 native checks.
- [ ] Deploy the shifted solve and recover D/F without weakening acceptance.
- [ ] Repair the groove-entry cross-film cooling inconsistency and revalidate A.
- [ ] Recompute final corrected equilibria before reporting comparisons.


### September 28 - Integrated steady force equilibrium

- [x] Identify that the native steady CFD loop held position fixed while the external driver performed expensive fully converged thermal trials.
- [x] Stop the superseded workstation jobs and preserve every run and the unsuccessful groove experiment.
- [x] Add bearing-position Newton corrections within the steady pressure/cavitation/energy iteration, without physical transient terms.
- [x] Verify force sign, geometry/probe consistency, timestep independence, frozen temperature probes and rejected-step restoration on native Windows at 1/2/4 MPI ranks; pass 22 affected CTests plus the CLI cap regression.
- [x] Complete the equivalent Linux regression checks: 22 focused CTests, four-rank IO/diagnostics and driver acceptance checks passed.
- [x] Deploy isolated revision 23fcbe7 and restart at 16:47 KST; verify seven cases on 28 disjoint physical cores with the original resolution and strict force/thermal checks.
- [ ] Resolve any remaining thermal-profile failure and return accepted temperature, friction and cavitation comparisons.


### September 28 - User-selected inner pressure stop

- [x] Confirm that the requested 1e-6 criterion refers to the inner pressure-equation residual.
- [x] Terminate all seven active cases that already reached that threshold and verify zero remaining campaign processes.
- [x] Preserve and synchronize existing output; mark the cases stopped by the user rather than thermally converged.
- [x] Set the local campaign FBNS tolerance to 1e-6 for any later authorized run.
- [ ] Add a portable stop-and-export path before a campaign that needs intermediate field recovery.
- [ ] Resolve the cross-film physicality defect before accepting final temperature comparisons.
