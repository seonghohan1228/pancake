# Journal-bearing physics and implementation audit

Later implementation: [COUPLED conservation](COUPLED_CAVITATION.md),
[thermal energy](COUPLED_ENERGY.md) and [R410A](R410A_EOS.md) supersede the
capability gaps identified here. This dated audit preserves the original
branch/history findings and failures; it is not a current capability list.

Remote follow-up on 2026-09-07: `/home/user/pancake` on `afml-r6` is detached
at `f0e772d` and contains unfinished campaign output. Its 69 reachable commits
are all present in the already audited local object database. The remote
branch tips include `case05_79ddbee`, `codex/conjugate-bush`, `codex/fbns-gas`,
`codex/thermal-ksp-robustness` and `ksp-fix-f0e772d`. No remote Git history
was rewritten or fetched. New runs use a separate source snapshot; see
[the workstation results](../cases/90_diagnostics/coupled_validation/SUMMARY.md).

Date: 2026-09-06. Baseline: local `main` at `b9dcaf3`. Scope: accessible
branches, reachable history, test evidence, governing equations and the target
oil-air / oil-R410A experiments. Work was isolated from the dirty Ferron checkout.

## Decision

Keep FBNS and the reduced Reynolds formulation as the inexpensive hydrodynamic
baseline. The exact 1-D equilibrium-cavitation problem is well resolved. The
code is **not yet a validated model of simultaneous vaporous and gaseous
cavitation**, and R410A is not implemented. Passing regression tests or a
converged Newton solve does not change that assessment.

This audit repairs reproducible conservation, property, convergence-reporting
and output-safety defects. It does not invent oil-specific property data, fit
unmeasured transfer rates, or certify the existing gas-pressure closure.

## Branch and history reconciliation

All 105 baseline reachable commits and the local/cached remote branch graph were
inventoried. Remote heads were checked read-only with `git ls-remote --heads`;
no fetch, push, remote merge or remote deletion was performed.
The [complete baseline commit inventory](physics_audit_history_2026-09.txt)
records full hashes and subjects. An existing detached worktree at `b74db94`
was left intact as well.

| Ref at inspection | Head | Assessment |
|---|---|---|
| Local `main` | `b9dcaf3` | Audit base; 51 commits beyond remote main, including journal-bearing feature history |
| `origin/main` | `87df724` | Older published baseline |
| `origin/feature/journal_bearing` | `e6742b5` | Already contained in local main |
| `origin/feature/bubble` | `72736c6` | Separate spherical/stereographic soap-film rewrite; replaces bearing components and is not a journal-cavitation implementation |
| `codex/ferron-rerun-terra` | `da81b06` | Later Ferron campaign work |
| `codex/ferron-fig6-flow` | `446327c` | Later Ferron flow/thermal comparison, including separate wall-temperature controls |

The Ferron branches contain 19 post-main commits between main and the later
head, mostly campaign analysis/docs plus wall-temperature implementation and
tests. Their dirty and untracked campaign files in `C:/Dev/pancake` were inspected
read-only and preserved. They were not silently folded into this physics repair.

Important history: `fd2b0ea` introduced FBNS; `69a9011` moved the full-film bound
for gas; `e487c25` and `13cab1d` shared coefficients; `f0e772d` strengthened MPI
robustness; `d726408` exposed outer-coupling failure; `df0802b` corrected the
groove-source accounting that had invalidated a favorable thermal comparison.
This repair builds on those changes rather than reverting to an older solver.

## Defects repaired

| Defect | Consequence | Repair and evidence |
|---|---|---|
| Current density/gap multiplied old concentration or film content | Moving film and changing properties created or lost species inventory | Store old `rho_l`, `h`, `theta`; difference complete `rho_l h theta` and `c rho_l h theta`; variable-capacity and closed-bearing tests |
| Gas transport subtracted `c div(F)` from conservative convection | Dissolved gas was not conserved when carrier flux diverged | Conservative species equation; use the exact FBNS face mass flux in dissolved transport |
| Phase exchange changed `c` at fixed liquid mass | Wrong solution-mass denominator; oil and total gas not conserved | Advance dissolved, free and liquid inventories together; reconstruct `c`, density and film content |
| Explicit exchange and relaxation of only free gas | Overshoot past equilibrium and artificial gas loss | Analytic frozen-state exchange, availability bounds, no transient numerical relaxation |
| Reservoir composition imposed only after the gas solve | Downstream cells saw the wrong inlet composition; unmeasured gas source | Impose reservoir rows inside both solves; retain displaced equations for source accounting |
| Diagnostics reconstructed gas flux after reaction/property changes | Reported boundary flux differed from the flux actually transported | Save transport boundary fluxes; include equal/opposite liquid phase-transfer source |
| Unknown gas species silently selected propane | An `R410A` input could run with propane properties | Strict AIR / PROPANE / R290 parsing; reject R410A and typos |
| Dissolved partial volume fell back to ideal-gas volume | Unphysical liquid-solution density | Require a positive dissolved-species partial liquid density for MIXTURE |
| Cross-film viscosity helper ignored the 2-D table | Mean-film and cross-film mobility used different viscosity models | Apply kinematic-table conversion inside the shared viscosity helper |
| Coupling checked only pressure and temperature | Evolving gas state could be called steady | Include concentration and gas-volume change in `R_g`; reject nonfinite metrics |
| Failed nonlinear solves ended with success status | Campaign scripts could score unconverged fields | Exit code 4 and explicit failure message for failed inner/outer solves or minimum-film stop |
| Physics advanced at label t=0 and through an extra final interval; geometry lagged the output time | Startup histories and integrated heating used the wrong elapsed time | Save the initial state, advance geometry to each target time, use the actual final partial interval in all solves, and test against exact uniform heating and rigid-body acceleration |
| Output preparation recursively deleted the result directory | Reruns could destroy reports and unrelated inputs | Remove only solver-named products; preserve unrelated files and skip symlinks |

The old 35-test Windows suite passed before these repairs. The initial added
physics regressions then produced nine failures. That demonstrates why an
unchanged green test suite was insufficient. Existing tests with incorrect
frozen-inventory assumptions were replaced with conservation assertions; their
physical expectations were corrected rather than their tolerances loosened.

Mechanisms and equations are in [PHYSICS: conservative gas update](../PHYSICS.md#18-conservative-gas-update-and-experimental-scope-2026-09-06).

## Evidence reviewed and rerun

| Evidence | Result | Interpretation |
|---|---|---|
| Full native Windows CTest, MinGW GCC 16.1, MS-MPI, PETSc, Release | 40/40 pass | Includes new physics regressions on 1/2/4 ranks and forced CLI failure |
| Full Ubuntu 22.04 WSL CTest, GCC 11.4, OpenMPI, PETSc 3.15, Release | 40/40 pass | Native Linux source/build compatibility; 6-rank tests need oversubscription on this host |
| Exact 1-D FBNS benchmark, 120/240/480/960 theta cells | Both fronts within one cell at all four grids | Equilibrium free-boundary verification only |
| Finest 1-D benchmark | Peak pressure +0.0701%; integrated load -0.0361% | Reference discretization accuracy, not experiment |
| Synthetic 48x8 oil-air transient, 1/2/4 ranks | Gas balance residual below 2e-11; final fields agree within 6e-10 normalized at saved VTK precision | Conservative, reproducible mild transient; coefficients are synthetic |
| Gas-volume admissibility in that transient | Positive gap-volume excess, about 2.7e-4 of the gap | Pressure/gas splitting remains physically inconsistent; reported explicitly |

The 1-D native wall times were of order 0.1-0.5 s across these four grids,
including process startup and output. They establish a small-case execution
scale, not a speedup, scaling result or performance forecast for conjugate THD.
See [reproduction and machine-readable evidence](../cases/90_diagnostics/physics_audit/README.md).

Historical validation was reviewed, not all rerun. The authoritative August
ledger remains available in [Validation status](../cases/VALIDATION_STATUS.md).

- The corrected Ferron peak-pressure comparison is +7.31% and 1.25 degrees at
  its scored condition. The earlier very close temperature agreement used a
  groove source that grew as `1/dz` and cannot be used.
- The original G6/G7 conjugate thermal grid check failed its 0.1 K mean/RMS
  limits (0.111/0.107 K). Later local D-series campaign files report passing
  G6/G7 and G7/G8 checks, and passing 100 ms timestep checks. Those are distinct
  configurations, so the later results do not retroactively validate the old run.
- Specifically, the later D-series mean/RMS changes are 0.01108/0.00750 K and
  0.01321/0.00916 K. Its T2/T3 and T3/T4 eccentricity differences are 0.000323%
  and 0.001240%, with zero quoted maximum-temperature difference.
- Those later files still set `temperature_validation_eligible: false`:
  unresolved fresh supply and missing bush axial / shaft end cooling remain.
  A thermocouple 0.5 mm below the bore is also not the same observable as an
  interface temperature or a film average. Relevant untracked evidence was
  `cases/01e_ferron_terra_repro/{campaign_analysis,d_grid_campaign_analysis,
  d_discretization_100ms_analysis}.json` in the original checkout.
- All scored Etsion cavity-extent runs used `alpha_gas = 0`. They test the
  equilibrium plateau model, not gas release, gaseous void fraction or kinetics.
- Absolute clearance and uniform load-model bias remain degenerate in the
  Ferron collapse. A fitted clearance is not independent load validation.
- The older Elrod front error decreases first-order in degrees but spans many
  cells. A roughly fixed error measured in cells is compatible with first-order
  convergence; it is not proof of mathematical nonconvergence. FBNS is still
  substantially more accurate in this benchmark. This is separate from the
  demonstrated nonlinear nonconvergence of the VOID_COUPLED operating points.

## What remains physically unresolved

1. **Gas pressure and occupied volume must close together.** The current FBNS
   solve uses lagged `alpha_gas` in `theta <= 1-alpha_gas`, then transports and
   reacts gas. The post-step fields can violate this inequality. A gas-volume
   cap can also make `alpha_gas` disagree with the gas EOS inventory. The new
   `gas_volume_excess` and `gas_eos_defect` diagnostics expose both. Reducing a
   residual or clipping theta does not solve this; clipping loses liquid mass.
2. **Steady feed flash is a local approximation.** It preserves the local feed
   fraction with the correct lever rule but does not solve conservative steady
   species transport. Its convergence does not establish the correct residence
   history or free-gas distribution.
3. **Vapor and noncondensable gas are not separate inventories.** Classical JFO
   gives an equilibrium pressure floor and liquid content, not explicit vapor
   mass, nucleation, latent heat, or a gas/vapor mixture EOS. Air expansion,
   air desorption and refrigerant evaporation cannot be distinguished by the
   existing single `free_gas_mass` field.
4. **Closure data are not generic.** Saturated R290/PZ68 tables cannot be reused
   for R410A or for arbitrary nonequilibrium dissolved concentration. A saturated
   two-dimensional viscosity table is only a surface through composition space.
5. **Thermal and transport splitting remain approximations.** The transient
   driver performs one flow/gas/energy sweep. The energy equation has no
   refrigerant latent/dissolution enthalpy closure. JFO gas fluxes and the
   energy/velocity reconstruction still need a full discrete-flux consistency
   audit beyond the exact FBNS dissolved-carrier path repaired here.
6. **Transient fluid-force coupling is explicit.** Timing is now corrected:
   the initial VTK state is unadvanced, every solved state and geometry share a
   target time, and the last interval ends at `end_t`. Geometry is predicted
   with the previous fluid force, so stiff fluid-structure feedback still
   requires a timestep study or coupled iteration. Old transient campaigns must
   be rerun before comparing their startup phases or fitted time constants.

These are open model/driver work items, not certified behavior. Existing R290
examples may now fail earlier because failure exit codes and conservative
kinetics no longer conceal their coupling problems.

## Route to oil-air and oil-R410A experiments

First establish the oil grade, temperature/absolute-pressure range, speed,
clearance, supply geometry and saturation/degassing procedure. Measure or bound
initial dissolved concentration and entrained bubble content separately.

For oil-air, use measured solubility and oil viscosity/density, then identify
release/resorption time scales independently of bearing pressure. Air released
from solution and pre-existing air bubbles need separate initial conditions;
they can produce similar cavity images by different mechanisms.

R410A is a 50/50 **mass** blend of R32 and R125, not propane and not a gas constant
substitution ([NIST REFPROP mixture table](https://www.nist.gov/system/files/documents/srd/REFPROP8_manua3.htm)).
Track component partitioning in oil and vapor, or justify a reduced blend
approximation against the actual oil's VLE data. Use suitable vapor EOS data
and composition-dependent liquid properties. If air is also present, its
noncondensable inventory and partial pressure are additional state variables.
Refrigerant release produces condensable vapor; calling all free material
"gas" must not erase that thermodynamic distinction. The historical R290
terminology and max-threshold construction are not a general multicomponent law.

For a fast solver, retain the surface mesh while testing the lubrication
conditions: `h/R << 1`, `h/L << 1`, and a small reduced inertia scale
`rho U h^2 / (mu ell)` along each relevant pressure-variation length `ell`.
Grooves, jets and local acceleration can violate these conditions even when
the land film is thin. Use resolved 3-D CFD selectively for those regions and
mechanisms when measurements show the reduced model is insufficient.

Next implementation priorities: conservative pressure/species/volume coupling
with simultaneous residual checks; transient coupling accuracy studies; calibrated
oil-air verification; oil-R410A component properties and thermal terms; then
independent pressure, load, flow, temperature and void-fraction validation.
Only profile and optimize this verified formulation. A PETSc-to-Trilinos
migration does not remedy the missing physics and is not part of this repair.
