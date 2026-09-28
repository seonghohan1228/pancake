# Change Log

## 2026-09-11 - Correct the oil screening geometry

- Rerun the steady oil matrix as case 02c using the user's 43 mm journal
  diameter, 40 mm bearing length and 42.5 micrometre radial clearance. Preserve
  case 02b and explicitly mark its results and recommendations superseded.
- Use a new local case and remote output directory so no old equilibrium or
  field is reused after the geometry change. Keep the oil, end pressure,
  prescribed groove and vapor-floor assumptions otherwise identical.
- Reassess the load and speed cuts using cavity filling deficit as well as
  extent. The much smaller clearance changes eccentricity and cavity fraction;
  the original operating-point recommendation is not carried forward.
- Add the subsequently requested energy-equation temperature plots for nine
  equilibria under adiabatic and prescribed-cooling assumptions. Record
  provisional heat capacity, conductivity and temperatures, keep viscosity
  fixed, and explicitly separate this sensitivity calculation from coupled THD.

## 2026-09-11 - Screen steady oil cavitation across load and speed

- Add case 02b with a standalone remote sweep driver, vector force-balance
  continuation, coordinate-based VTK extraction, common-scale contour plots,
  and an audit of residuals, physical fields, force closure and mass balance.
- Use the requested 0.210 Pa s, 937 kg/m3 oil and 100000 Pa end/groove pressure.
  Record inherited geometry, assumed groove, 500 Pa cavity threshold, executable
  commit and binary hash so the scan can be reproduced without altering prior
  campaigns. Retain solver source and runtime unchanged.
- Add selected 480 by 128 repeats of the 240 by 64 matrix, a 100 rpm probe,
  physical nondimensional axes and explicit limits on vapor/thermal/deformation
  interpretation. Keep generated fields and images out of version control.

## Unreleased - Windows MSYS2 GUI Workflow

### 2026-09-08 - Clarify model selection and Ferron comparison scope

- Add a source-checked model guide covering all cavitation, material, thermal,
  boundary and motion choices and their compatibility. Explain why timestep
  consistency, not steady load, is the numerical criterion, and identify
  bubble dynamics absent from the homogeneous model.
- Derive the distinction between the older liquid-temperature equation and
  coupled mixture energy. Separate fixed-position load changes from the
  position changes needed to maintain an experimental applied load.
- Re-read Ferron et al. (1983) and correct claims that the reduced bush and
  shaft exactly reproduce its spatial heat-conduction models. Clarify end-face
  versus outer-cylinder cooling, sensor temperature, and hot clearance.
- Label the separate Ferron branch's archived thermal/mass-supply disagreement
  explicitly; do not present those values as new coupled-model results.
- Record the user's local edit/build/check, remote source/build/run, local
  result-return workflow. Documentation-only review; no solver behavior or
  campaign results changed and no CFD rerun claimed.

### 2026-09-07 - Couple gas volume, conserved species and film heating

- Add COUPLED pressure/oil/dissolved/free-species conservation and optional
  temperature in one sparse nonlinear solve. The old split update could
  converge pressure while leaving gas and liquid occupying incompatible
  volumes. Keep previous physical storage fixed during nonlinear iteration.
- Recover the pure-vapor JFO limit, permit finite release and separate
  absorption, and offer explicit local phase equilibrium. Preserve all
  inventories through the pressure/volume closure instead of capping gas
  volume after transport. Rationalize complementarity residuals to avoid
  cancellation near a disappearing phase.
- Use the actual conservative component fluxes for energy and published
  velocities. Evaluate known gap geometry at the face; upwind only the
  transported state. Keep pressure and Couette donors independent.
- Add internal energy with viscous dissipation, compression work, wall and
  in-plane heat transfer, reservoir enthalpy and desorption energy. Do not
  call the legacy liquid-only thermal update after the coupled solve.
- Make both pressure ends and internal groove/hole reservoirs bidirectional
  for composition and heat. Outflow extrapolates local properties as well
  as transported state; a downstream temperature must not impose a false
  cold discharge resistance. A draining groove does not pin film temperature.
- Implement the fixed-composition R410A Helmholtz EOS with density, internal
  energy and heat capacity. Reject a resolved subcritical free gas above its
  dew line. Allow film pressure above the 20-bar reservoir range; do not
  confuse reservoir operating limits with peak hydrodynamic pressure.
- Correct an outlet-domain check exposed by the remote 20°C run: numerical
  perturbations of an outgoing gas-free state need the same EOS continuation
  as interior trials. Rejecting them as real bubbles broke the Jacobian before
  iteration. Reproduce the failure with a high-pressure dissolved-solution
  test, repair the gate, and retain strict checks on accepted free phases.
- Add actual-concentration/temperature log-viscosity interpolation. Reject
  extrapolation and saturated pressure/temperature property substitutions in
  the coupled path. Label the DOE FVC68B table separately from FW68D.
- Expose the solver, material table, independent kinetics and thermal inputs
  in the GUI. Remove unsubstantiated mixture validation claims and inherited
  propane material values from the R410A preset.
- Correct the cavity-extent analysis to count total liquid deficit, including
  gas cavities. Increase VTK numeric precision so output rounding does not
  hide timestep differences. Update numeric output-time assertions.
- Pass all 51 regression tests on native Windows, local WSL Linux and the
  isolated remote Linux build. Run the paired and thermal campaigns on the
  workstation, recording machine/MPI provenance. An eight-rank pilot gives
  3.70× speedup with the same load as one rank. Retain the Etsion mismatch
  and failed R410A cases in the comparison report.
- Document the equations and limitations in COUPLED_CAVITATION.md,
  COUPLED_ENERGY.md and R410A_EOS.md. Full FW68D validation, separate condensed
  refrigerant, blend fractionation and general thermal evaporation remain open.

### 2026-09-07 - Correct reservoir states and separate gas absorption

- Give each axial pressure reservoir an independent inflow temperature,
  dissolved mass fraction and entrained free-gas mass fraction. Outflow uses
  the interior state. Record actual incoming free-species flux as the liquid
  carrier flux times `q/(1-q)`. Defaults reproduce the previous boundary state.
- Reject explicit end temperatures paired with thermal extrapolation, invalid
  reservoir fractions, and reservoir pressures below the effective cavity
  floor. Previously `LOCAL_PSAT` could silently raise a requested outlet to
  the initial charge's bubble point. Remove the advice equating atmospheric
  reservoir pressure with a universal cavitation pressure.
- Add independent `gas_absorption_rate`, including absorption-only and zero
  absorption limits, while retaining conservative analytical exchange.
  The historical equal-rate choice remains the default, not a calibrated law.
- Verify pressure-gradient reversal against discrete analytical species and
  energy profiles, gas inventory against boundary flux, independent absorption
  against the exact reaction, and compression against the ideal-gas volume law.
  All 40 tests execute successfully on Windows and Linux; retain the Linux
  six-rank launch failure and successful resource-setting retry in the evidence.
- Write the detailed FW68D/air/R410A physics review, distinguishing corrected
  code from unresolved pressure/species/volume coupling and missing material
  data. Add reproducible kinetic and pure-R410A EOS diagnostic figures.
  Correct two overstatements in the archived Etsion report and replot its
  measured/archived reformation angles with the stated uncertainty.
  No validated R410A bearing capability is claimed.

### 2026-09-06 - Measure the effect of the physics repairs

- Rebuild pre-audit `b9dcaf3` and compare with repaired `2bc4d3c` on identical
  frozen Ferron, optical-rig, synthetic oil-air and analytic timing inputs.
  Add a common direct-inventory reaction probe without modifying either core.
- Archive hashes, configurations, full diagnostics, portable profiles and
  histories, and five publication-style comparison figures in
  `cases/90_diagnostics/physics_comparison`. Both steady cases retain their
  pressure/temperature fields at saved precision. The old closed reaction
  loses 11.548% of total gas in 20 ms; the repair conserves it to roundoff.
- Align static transient histories by actual elapsed time and demonstrate
  the old excess heating and missing final motion state with exact solutions.
- Extend the air test to 1 ms, then halve the timestep twice. Inner solves
  converge and gas inventory is conserved, but load alternation persists and
  grows. Record this failed timestep check explicitly in PHYSICS, PLAN and
  validation status; do not turn a terminal-sample difference into a physical
  load claim. No pressure/volume-coupling repair is claimed by this comparison.
- Keep production output and the dirty campaign checkout untouched. Plot
  dependencies were installed under ignored build output, not globally.

### 2026-09-06 - Conserve cavitating-mixture inventories and audit readiness

- Audited all accessible branches, reachable history and validation evidence
  from isolated local main; preserved the dirty Ferron checkout. See
  `docs/PHYSICS_AUDIT_2026-09.md` for findings, provenance and unresolved work.
- Difference full old/new liquid and dissolved inventories so changing gap,
  density and solution mass do not create species. Replace nonconservative gas
  convection with conservative transport using exact FBNS carrier fluxes.
- Integrate frozen-state release/resorption analytically, update liquid mass
  and concentration together, and remove mass-losing transient gas relaxation.
  Correct the steady feed-flash lever rule and retain capacity-limited excess
  in solution. The steady flash remains an approximation to species transport.
- Impose gas reservoir composition inside both linear solves and account for
  displaced rows. Preserve the boundary flux actually used by transport and
  include phase transfer in liquid continuity diagnostics.
- Reject unknown species (including unsupported R410A), impossible mass
  fractions and missing dissolved partial-volume density. Share 2-D viscosity
  table conversion with cross-film property calls. A zero exchange rate is valid.
- Include gas change and nonfinite metrics in steady convergence; preserve
  history initialization in restored snapshots. Return exit 4 for failed
  nonlinear solves and minimum-film termination instead of reporting success.
- Add gas-volume and EOS-defect diagnostics without concealing volume error
  through nonconservative clipping. The pressure/gas/volume closure is still open.
- Correct the transient driver: save the unadvanced initial state, predict
  geometry at the target time from the previous fluid force, use a shortened
  final interval consistently in flow/gas/energy/motion, and always save the
  final state. Permit end_t = 0 for initial output. A CLI regression checks exact
  insulated Couette heating and rigid-body acceleration over a partial step.
- Replace recursive output-directory deletion with targeted solver-product
  cleanup, preserving unrelated reports, inputs and symlink targets.
- Add MPI-1/2/4 conservation regressions and a CLI failure contract, plus exact
  1-D and synthetic air reproduction scripts. Full native Windows and Ubuntu
  suites pass 40/40. Update README, PHYSICS, PLAN and validation claims; correct
  the reversed bubble-point wording and the front-convergence interpretation.

### 2026-08-11 - Ferron campaign branch integration

- Ran the case-02 optical-rig sweep at the confirmed lab-oil viscosity
  (ISO VG 46 white mineral oil, mu(20 C) = 0.068 Pa s by grade conversion),
  with the isopentane threshold `iso` = 0.12 bar added from Raoult and the
  eccentricity band tightened to the measured +/-0.02 stack-up. Added a `MU`
  sensitivity suite because cavity geometry is nonlinear in viscosity even
  though load and pressure are linear. Gated the dispatch behind `verify_setup.py`:
  V1 anchors the unvalidated grooveless ends-fed configuration against the
  Ocvirk short-bearing solution, V2 confirms G6 transfers to this geometry
  (0.004 % in load across G6 to G7). 156 of 159 points completed; RUN1 is
  complete at 75/75. The three failures are physical rather than infrastructural
  and are documented in `sweep_summary.md`.

- Completed the 0.10 s physical timestep ladder (40/20/10/5 us) that had been
  interrupted after its first rung. The 20/10/5 us rungs agree to 0.0012 % in
  final eccentricity, 0.0004 degrees in attitude and machine precision in load,
  so the frozen `dt = 2e-5 s` is now verified over a 0.10 s settling window
  rather than 0.06 s. The 40 us rung carries visible transient noise and stays
  a report-only stability probe. Report figures regenerated from the extended
  history.

- Folded the D-grid G8 rung (1620 x 432) into the local grid-convergence ledgers
  and regenerated `terra_d_grid_convergence.png`. G7 -> G8 passes every gate with
  the largest margin at 20 % of limit (eccentricity 0.050 % against 0.5 %, peak
  bearing temperature 0.012 K against 0.25 K), so G6 remains the frozen campaign
  resolution and the ladder is closed at both ends.

- Added case 02, the optical-rig operating-window sweep, for the Rev B geometry
  (D=49.3 mm, c=350 um, L/D=1, grooveless ends-fed, paraffin oil at 20 C). It
  runs pinned-eccentricity steady FBNS solves over eccentricity, speed and
  cavitation threshold, and reduces them to cavity arc, start/end angles, axial
  coverage, peak void, film load and dissipation. Two substitutions against the
  original specification are recorded in the case README: the prescribed
  saturation pressures are set through `SCALAR_PCAV` rather than `LOCAL_PSAT`,
  because `LOCAL_PSAT` derives `p_sat` from the Henry solubility model and
  degenerates to `p_cav` under constant fluid properties; and the L/D variation
  suite is dropped because the rig is now at L/D=1, above the 0.75 it proposed.

- Added a separate Ferron Figure 6 campaign at 4000 rpm and 6000 N using the
  user-digitised Figure 10 fresh supply, `Qi=7.837 L/min` at `eps=0.578`.
  It runs depth-averaged, adiabatic cross-film, fixed-40 C cross-film, and
  fixed-outer-40 C CHT models on G6. The axial-face `h_b=80` case is excluded
  because axial cooling remains outside the present reduced CHT model.
- Added compact Figure 6 profile export and a 300 dpi temperature comparison
  against Ferron's subsurface bearing thermocouple data.

- Reordered the Ferron thermal-model comparison and standardized its line map:
  depth-averaged grey dotted, adiabatic cross-film blue, 40 C cross-film
  orange, fixed-40 C CHT yellow, and outer-Robin CHT green.

- Expanded the 2000 rpm / 4 kN temperature figure to compare cross-film THD
  with isothermal walls, depth-averaged THD, outer-Robin heat transfer, CHT with
  fixed outer bush temperature, and adiabatic cross-film THD. The compact visual
  exporter now retains the depth-averaged temperature profile for each case.

- Added the 2000 rpm / 4 kN Ferron temperature comparison using the current
  R1-G7 bearing-interface temperature. No last-year circumferential temperature
  curve is drawn because the LGEP exports contain only unverified scalar time
  histories, not a comparable temperature field.

## Unreleased - Windows MSYS2 GUI Workflow

### 2026-08-07 - Ferron thermal-boundary audit corrections

- Added optional `shaft_ambient_temperature` and `bush_outer_temperature`.
  Positive values independently set the scalar shaft sink and bush outer-boundary
  temperatures; zero or negative values preserve the historical fallback to
  `thermal_ambient_temperature`. Parsing, serialization, validation, and the
  cross-film config round-trip test cover both keys.
- Corrected the Ferron documentation: the present bush model is z-collapsed and
  omits the paper's axial-face \(h_b\) condition; the scalar shaft closure omits
  Ferron's r-z shaft conduction and 43 C end temperatures. These changes expose
  boundary-condition sensitivity without claiming a full-solid reproduction.
- Added isolated case `01e_ferron_terra_repro` for the new Ubuntu campaign. It
  keeps smoke and production archives separate, distinguishes model-derived
  groove source from experimental flow, and pre-registers R0/R1/R2 plus the
  corrected reproduction of Claude's Fig. 4/Fig. 5 and wall-model matrix.
- Added a distinct discretisation add-on to case 01e: a fresh G5/G6/G7 steady
  force-balanced CHT grid ladder and a physical moving-bearing timestep ladder
  from 80 to 5 microseconds. This does not reuse Claude's invalid nominal THD
  grid pair or conflate steady grid convergence with physical-time refinement.
  Failed coarse stability probes are archived without aborting finer timestep
  legs.
- Added reproducible 300 dpi D-grid convergence, surface-temperature/pressure
  profile, and film-field contour figures from the transferred Ubuntu data.
- Completed the 24-rank Ubuntu timestep ladder. The 80 microsecond leg became
  unstable and stopped early, 40 microseconds reached the endpoint but missed
  force balance, and the registered 20/10/5 microsecond comparisons passed.
  Added full eccentricity histories and report-style grid, timestep, Ferron
  pressure, and CHT-temperature figures. Convergence panels reserve 10% of the
  data range above the maximum and 40% below the minimum.
- Extended the fresh force-balanced R1 grid sequence from G5--G7 to G1--G7 so
  the total-cell convergence plot includes the inexpensive coarse grids rather
  than borrowing points from Claude's fixed-eccentricity study.
- Extended the transient stability sweep above 80 microseconds with 160 and
  320 microsecond probes. Incomplete coarse histories remain visible as dashed
  curves and in a separate time-reached plot, but are excluded from the
  common-60-ms endpoint comparison.
- Added an explicit `--rerun-failed` recovery path for timestep launch errors;
  it cannot overwrite an archived solver run with `rc=0`.
- Added G8 (1620x432) to the fresh R1 grid study and a distinct 0.10 s
  40/20/10/5 microsecond ladder. The report averages eccentricity over the last
  50 ms, labels it only as mean eccentricity, omits diverged histories, fixes
  the mean panel to a 0.01 span with its maximum at 60% height, and applies the
  10%/40% padding rule to the time-reached panel.
- Added an unattended queue that waits for the independent G8 result before
  launching the four 0.10 s cases, preserving the one-24-rank-job policy.
- Added a separately labelled repaired-legacy reproduction of Claude's original
  G1--G5 isothermal grid sweep and nominal THD G3/G4 pair. Minimal validator
  repairs are machine-readable; non-converged coupling remains a failed verdict.

### 2026-08-07 - Ferron conjugate-THD grid protocol

- Added case 00d, a pre-registered G5/G6/G7 film-plane grid check of the final
  Ferron conjugate equilibrium. Every grid independently re-closes the two-pass
  groove model and force balance with converged adaptive pressure-temperature
  coupling.
- Acceptance limits target the quantities actually quoted in validation:
  eccentricity, attitude, cell-centred peak pressure, and the mean, RMS, peak,
  and peak location of the bush-temperature comparison.
- Fixed inlet-flow diagnostics at a full-length pressure groove. Inlet row
  replacement already discarded overlapping axial boundary links, but the
  diagnostic restored them and made the inferred supply flow grow as `1/dz`.
  Boundary and inlet-source accounting now use the same row precedence as the
  assembled FBNS and Elrod operators; a two-grid regression test pins it.
- The clean Ubuntu rerun removes the false 16.26 to 30.65 L/min divergence;
  corrected supply flow is 4.94 to 5.18 L/min. The registered grid verdict still
  fails narrowly on mean/RMS temperature change, and the corrected G6 validation
  score is +3.82 K mean / 4.02 K RMS rather than +0.02 K / 1.17 K.
- Added reproducible Ferron G5--G7 profile comparisons and G6 pressure and
  temperature contours. Remote Ubuntu VTK data are reduced to a portable NumPy
  archive, while the publication-style PNG figures are rendered on Windows.
- Added a direct profile and pointwise-difference figure comparing Claude's
  archived G6 conjugate result with the corrected G6 result.

### 2026-08-07 - Cross-platform baseline and validation-claim reconciliation

- Rebuilt the exact solver commit from a fresh Windows preset and verified all
  35 tests, then transferred the hash-verified commit to an isolated Ubuntu
  checkout, built in Release mode, and verified the same 35 tests there.
- Added `cases/VALIDATION_STATUS.md` as the authoritative claim ledger and marked
  the older root plan, first-campaign report, and unfinished-business audit as
  historical snapshots.
- Corrected three overclaims: the frozen campaign grid is not a universal
  `p_max`/coupled-THD certificate; the clearance collapse does not distinguish a
  uniform load bias; and Etsion case 05a has `alpha_gas = 0`, so it does not
  validate gaseous-cavitation physics.
- Made the Linux build command generator-neutral and Release by default, updated
  both platform test baselines, and removed the CMake FetchContent timestamp
  warning for the pinned GUI dependencies.

### 2026-08-05 - Energy axial boundaries and gas transport routed through ReynoldsCoeff (audit A / S1)

`src/reynolds_coefficients.hpp` exists so the three assemblies of the Reynolds
equation cannot drift apart. Two files never honoured it:

- **`energy.cpp`** branched correctly on cross-film for its interior faces and its
  dissipation term, but `south_boundary_flux` and `north_boundary_flux` hand-rolled
  `u_z = -g h^2 dp_dz / (12 mu)` with no cross-film branch at all. So with
  `cross_film_thermal_model = POLYNOMIAL` the energy equation transported heat
  through its interior faces with Dowson's `F2` and out through its axial ends
  with the classical law. The Ferron cross-film case has Dirichlet axial ends, so
  those are its oil outflow faces.
- **`gas_transport.cpp`** had zero references to cross-film anywhere: neither the
  Couette thickness `2 F1/F0` nor the Poiseuille mobility `F2` reached its face
  velocities, so gas was advected on a velocity field inconsistent with the
  pressure solve that produced it.

Both now call two new helpers:

- `ReynoldsCoeff::couette_face_velocity` - `U/2`, or `U (F1/F0)/h` generalized.
- `ReynoldsCoeff::poiseuille_face_velocity` - `-g h^2 (dp/ds)/(12 mu L)`, or
  `-g F2 (dp/ds)/(h L)` generalized, so that `h` times the velocity is exactly the
  `-g F2 grad p` the Reynolds solve conserved.

Both return the classical expression with its **original floating-point
association** when the model is off. That is not a claim: the full 35-case suite
was captured verbatim before and after and the only line that differs is
`test_ksp_guard`'s nondeterministic MPI abort text.

**A/B on a real cross-film point.** `cases/01_ferron_wallbc_fields/rerun_xf_adiabatic/fb_iter_13`
(Ferron 1983 Fig. 4, 2000 rpm / 4 kN, 720x192, FBNS + ENERGY_EQUATION +
`cross_film_thermal_model = POLYNOMIAL`, `bush_thermal_model = LOCAL`, Dirichlet
axial ends, the converged force-balance iterate) was rerun locally on 4 ranks with
the pre-change and post-change binaries. The old binary reproduces the archived
24-rank numbers (`p_max` 1.41105e6, `T_max` 321.768 K, `friction_torque` 2.59609,
`W` 4000.06 N), so the comparison is against the real archived point.

VTK output carries 6 significant digits, which is the resolution floor of this
table. Over 138 240 cells:

| field | scale | max abs diff | max rel diff | reading |
|---|---|---|---|---|
| `p` | 1.41105e6 Pa | 10 Pa | 7.1e-6 | one printed digit |
| `Fp` (load) | 4000.06 N | 0.01 N | 2.5e-6 | one printed digit |
| `T` | 321.768 K | 0.001 K | 3.1e-6 | one printed digit |
| `T_bearing` / `cross_film_t_bush` | 323.953 K | 0.001 K | 3.1e-6 | one printed digit |
| `friction_torque` | 2.59609 N m | 1e-5 | 3.9e-6 | one printed digit |
| `mu` | 0.0272947 Pa s | 2e-7 | 7.3e-6 | one printed digit |
| **`cross_film_core_rise`** | 1.06394 K | 5.6e-5 K | **5.3e-5** | above the floor |
| **`q_journal`** | 14467.1 W/m^2 | 0.6 | **4.1e-5** | above the floor |
| `film_content` | 1 | 3e-5 | 3.0e-5 | one front cell |
| `h`, `T_journal`, `cross_film_t_shaft` | - | 0 | 0 | unchanged |
| `mass_residual`, `q_bearing` | ~1e-10, ~3e-9 | - | >1 | noise about zero (adiabatic case) |

Printed summary: `T_shaft` 318.279 K and `T_bush_mean` 316.974 K unchanged;
`core_rise_max` 1.06394 -> 1.06395 K; the steady coupling loop converged in **15**
outer iterations instead of 16, both inside `coupling_tolerance = 1e-6`. Global
liquid mass moved in the 9th significant digit
(2.96795788e-3 -> 2.96795780e-3 kg).

**Reading.** Exactly what the audit predicted: small and boundary-driven. Nothing
in the load-capacity or peak-temperature comparison moves at the precision the
validation figures are drawn at, so no archived verdict changes. The two
quantities that do move beyond the print floor - the cross-film core temperature
rise and the journal-surface heat flux - are the two that are computed *from* the
wall/boundary energy balance, which is precisely the balance that was being closed
with two inconsistent mobilities.

**Not done here.** `energy.cpp`'s interior faces keep their own (already
cross-film-aware) expressions rather than being re-routed, because their
cross-film branch computes the volumetric flux `q` directly instead of `h*u` and
re-routing would perturb cross-film results for no correctness gain. And
`cfg.reformation_boundary_velocity`, the Elrod cavitated-end reflood closure,
still takes `(h_face, mu_face)` with no cross-film branch; it is only reachable on
the Elrod path with `consistent_boundary_flux`, and is recorded in
`docs/reports/QUICKWINS_2026-08.md` rather than changed unreviewed.

### 2026-08-05 - Shipped `config.txt` now defaults to FBNS (audit Q5)

`config.txt` is the file a new user starts from and it pointed at the superseded
Elrod-Adams path. Of the 338 case configs under `cases/`, 289 run `FBNS` and 41
run `JFO`; every one of the 338 sets the key explicitly, so this changes no case.

`cavitation_model = JFO` -> `FBNS`, and the surrounding comments were wrong in two
ways:

- They described the FBNS scope as *"Phase 1: ISOTHERMAL + fluid_property_model =
  CONSTANT"*. `validate()` has supported `ENERGY_EQUATION` with `CONSTANT`,
  `SINGLE_PHASE` and `TWO_PHASE`, plus the cross-film thermal model, for two
  phases now. The comment now states the real scope *and* what FBNS refuses
  (`gas_pressure_coupling = VOID_COUPLED`, `INLET_OUTLET` axial boundaries).
- `max_outer_iters` and `outer_tol` sat directly under the cavitation model with
  nothing saying they are Elrod-only. They are now labelled as not used on the
  FBNS path, which is the single most confusing thing about switching models.

**The struct default stays `JFO`.** Eleven test files never set
`cavitation_model` and would silently change solver path if it flipped:
`test_bush_conduction`, `test_cross_film`, `test_eos`, `test_fvm`,
`test_gas_transport`, `test_io`, `test_journal_motion`, `test_ksp_guard`,
`test_linear_system`, `test_property_tables`, `test_velocity`. Pinning their model
explicitly is a separate change with its own review; splitting it this way keeps
the test-path change auditable instead of implicit.

`test_config_keys` now also requires `config.txt` to pass `validate()`, not just
to parse - that is what holds the new default against the FBNS scope rules.

Left alone: the three `config_r290_*.txt` still run `JFO`. Those are gaseous
cavitation cases and the R290 work is where `VOID_COUPLED` lives, which FBNS
refuses; flipping them is a modelling decision, not a default cleanup.

### 2026-08-05 - The GUI no longer rewrites motion settings on load (audit 3.5 / S9)

`GuiApp::normalize_motion_config()` ran on **load** and again on **save**. Opening
an externally authored config overwrote `bearing_stiffness_y/z` and
`bearing_damping_y/z` with the x value, forced `bearing_initial_from_attitude =
true`, and zeroed all six `bearing_initial_*`. Saving then wrote the mutilated
version back. A config with `kx = 1e6, ky = 2e6` loaded and saved in the GUI came
back isotropic with no warning - silent data loss on a file the user did not
author in the GUI.

The isotropy contract itself is a legitimate design choice: the editor exposes one
support-stiffness and one support-damping row for all three directions. The defect
was applying it destructively to data that came from elsewhere. It is now applied
where it belongs:

| when | before | now |
|---|---|---|
| New case | normalize | normalize (unchanged - the GUI authors this config) |
| Load | normalize | preserve, and log what the editor cannot represent |
| Edit support row | - | mirror x onto y/z (`GuiApp::mirror_isotropic_support`) |
| Save | normalize | write what is in the config |

So a GUI-built case is still isotropic by construction, and an anisotropic or
explicitly initialised config now round-trips through load/save unchanged. On load
the user gets a log note naming what the editor will not show them, rather than
finding out from a diff. The same loop now also surfaces `parse_errors` from the
Q1 change in the GUI log.

Documentation updated to match: the stale "the GUI mirrors it onto the x/y/z
config members at save" comment in `schema.cpp`, and the motion sentence in
`docs/GUI_ARCHITECTURE.md`.

Not covered by ctest: `GuiApp` links Win32/ImGui and has no test target, so this
was verified by inspection and a clean GUI build, not by an automated test. That
gap is real and is recorded in the batch report.

### 2026-08-05 - van Leer far-upwind cell now follows the flux sign (D1 / F4)

`FVM::divergence` hard-coded the far-upwind cell regardless of the face flux
sign. On the east face it always differenced against `i-1` and, for `F_e < 0`,
inverted the positive-flux gradient ratio. For a negative face flux the upwind
cell is `i+1` and the far-upwind cell is `i+2`, so the scheme was reconstructing
from the wrong stencil with the wrong `r` - neither TVD nor second order. Same on
the west face.

The reconstruction is now written in terms of the face's own upwind triple
(far-upwind, upwind, downwind), which flips with the flux sign. The positive-flux
branch is algebraically and floating-point identical to the old code.

**Reclassified from latent to live.** `CHANGES.md` recorded D1 as "latent for
omega > 0" and `PLAN.md` had F4 as TODO. That caveat held only for the Elrod
theta equation, whose face flux is pure Couette and therefore carries the sign of
`omega`. It does not hold for the two paths the campaign added or extended:

- `energy.cpp` face flux is `0.5*omega*R*couette - g*F2*dp_dtheta/R`; the
  Poiseuille term can reverse the face with `omega > 0`, and
  `thermal_convection_scheme = VANLEER` reaches `FVM::divergence` unmodified.
- `gas_transport.cpp` has the same structure via `gas_convection_scheme`.

`omega < 0` is also legal (`validate()` never requires positivity), which puts
every Elrod face on the reversed branch.

**Test first, then fix.** `tests/test_schemes.cpp` gained a directional-symmetry
case (the F4 acceptance gate): the mirrored step profile advected with mirrored
flux is the same discrete problem up to a permutation, so its L1 error must match
the forward run's. Against the old code it fails outright -

| | forward | backward (old code) | backward (fixed) |
|---|---|---|---|
| L1 error | 0.0266008 | 0.0626316 | 0.0266008 |
| range | - | [0.140741, **1.00469**] | [0.2, 1.0] |

- 2.4x the error and a genuine overshoot past the initial maximum, i.e. the
reversed branch was not bounded either. After the fix the two agree to every
printed digit and the bounds are exact.

**Nothing else moved.** The full 35-case suite was run on the pre-fix and
post-fix binaries with verbose output captured and diffed: the only differences
are the three new directional-test lines (the test is registered at 1, 4 and 6
ranks) and nondeterministic MPI abort formatting in `test_ksp_guard`. No existing
test's numbers changed, which is the expected result - every currently tested
VANLEER path has non-negative face fluxes, so the fix is a no-op on them.

**Not fixed here, and now documented:** the deferred correction still covers the
circumferential faces only, so axial convection degrades to first-order upwind
under `VANLEER`. That is the second half of F4 and is stated in `PHYSICS.md 5.5`
rather than left for a reader to discover.

### 2026-08-05 - An unrecognized config key is now a hard error (audit Q1)

`SimulationConfig::load_from_stream` had no terminal `else`: any key the parser
did not recognize was silently discarded and the run continued on the default the
user never asked for. It now records the key in a new `parse_errors` list, which
`validate()` promotes into its error list, so every caller that already refuses to
run on a validation error - the solver's `main`, the GUI's config check - refuses
to run on a misspelled key too.

**This was not hypothetical.** `cases/01_ferron_thd/config_template.txt` set
`cross_film_model = NONE` while the parser has only ever known
`cross_film_thermal_model`. That line was dropped on every run of the Ferron THD
case and was replicated into 105 generated `fb_iter_*/config.txt` files across
`fig4_*`, `fig5_*` and the rest. It happened to be harmless - `NONE` is the
default - but nothing in the toolchain could have told anyone that. All 105 files
and the four mentions in `cases/01_ferron_thd/SUMMARY.md` now use the real key.

The error names the offending key and suggests the closest known one, using two
signals rather than one:

- **Edit distance**, length-scaled (`max(2, len/3)`), for ordinary typos:
  `bush_outer_radus` -> `bush_outer_radius`.
- **Shared underscore-separated words**, for the wrong *name* of the right knob.
  `cross_film_model` is eight edits from `cross_film_thermal_model` - far outside
  any sane distance threshold - but shares all three of its words. This is the
  case that actually shipped, and distance alone would have missed it.

An unrelated token (`totally_made_up_key`) matches neither and is reported with no
suggestion rather than a misleading one.

The deliberate legacy branches are untouched and pinned by test: the key aliases
(`viscosity_model`, `oil_gas_solution_model`, `gas_mixture_viscosity_model`), the
value aliases (`ELROD_ADAMS` -> `JFO`), and the warn-and-drop keys for removed
models (`pressure_time_method`, `temperature_time_method`,
`dissolved_gas_bunsen_coeff`) still parse as warnings, not errors.

New `tests/test_config_keys.cpp` (35 ctest cases now, was 34):

1. The shipped `cross_film_model = NONE` line, verbatim, must produce exactly one
   parse error naming the key and suggesting `cross_film_thermal_model`.
2. Every legacy branch still parses with no error and the expected warning count.
3. Every key in the new `known_config_keys()` manifest is accepted by the parser -
   this is what keeps the suggestion list from drifting away from the chain.
4. `to_config_text()` output re-parses with no unknown key, so the writer and the
   reader cannot diverge.
5. Every `config*.txt` in the repo root and under `cases/` parses with no unknown
   key. `*_template.txt` files are included, with their `{placeholder}` values
   substituted, because the templates are exactly where the dropped line lived.

Known limitation, unchanged by this work: a *recognized* key with a malformed
value still throws `std::invalid_argument` out of `std::stod` rather than being
reported. That is a separate defect and is not addressed here.
### 2026-08-05 - The Picard loop returns its best iterate and admits when it did not converge

**Defect.** The STEADY_STATE outer (Picard) loop in `main.cpp` returned whatever
iterate `coupling_max_iters` happened to land on, and said nothing about having
failed to meet `coupling_tolerance`. On the conjugate Ferron cases
(`bush_thermal_model = CONDUCTION_2D` with a weak outer Robin surface) that is a
limit cycle: `max(R_p, R_T)` decays geometrically to 4.3e-5 at outer 56, then
oscillates between 5e-5 and 3.6e-4 for the remaining 144 sweeps. The delivered
iterate 200 was **6.3x worse** than iterate 56, on every one of the eight
`rerun_cht_robin10` force-balance calls, and on four of the six `outer_robin80`
calls that ran to completion (the other two scraped in at 187 and 190 sweeps).

**Fix, in a new `Coupling` module** (`src/coupling.hpp`, `src/coupling.cpp`) so
the loop control is testable without a solver:

- `Coupling::Controller` ranks iterates by `max(R_p, R_T)` and, **only when the
  loop exhausts the cap without meeting the tolerance**, asks for a rollback to
  the lowest-residual sweep. A loop that meets the tolerance always returns its
  last iterate, so every converging case is bit-identical to before.
- `Coupling::Snapshot` deep-copies every registered field, ghost layers and
  `old_data` included, so the restored state is one self-consistent iterate
  (p, theta, T, wall state, properties, forces) and needs no ghost exchange. The
  snapshot is taken *before* under-relaxation, so what is delivered is a genuine
  end-of-sweep solve rather than a blend.
- A non-converged steady run now prints `NOT CONVERGED`, the cap it hit and the
  best-vs-final comparison **in red regardless of `log_outer_iters`**, and
  `results/diagnostics.csv` gains `coupling_iters, coupling_R_p, coupling_R_T,
  coupling_converged, coupling_best_iter, coupling_best_residual,
  coupling_restored_best, coupling_relaxation`. New columns are appended, so
  header-keyed readers are unaffected.

**Under-relaxation is the actual cure, and it was measured, not assumed.** Sweeps
to `coupling_tolerance = 1e-6` on the `rerun_cht_robin10` equilibrium, 24 pinned
ranks on afml-r6:

| `coupling_relaxation` | 1.0 | 0.95 | 0.9 | 0.8 | 0.7 | 0.5 | 0.3 | ADAPTIVE |
|---|---|---|---|---|---|---|---|---|
| outer sweeps | **200, never** | 48 | 53 | 64 | 76 | 113 | 197 | **52** |

Even 0.95 converges, so the loop gain is only marginally above 1 - which is why
`coupling_relaxation_mode = ADAPTIVE` (new, default off) works: it starts at the
ceiling, halves the factor the first time the residual rises and lets it recover
at 1.1x per falling sweep. It reached the tolerance in 52 sweeps with no tuned
value, ending at omega = 0.974 having dipped to 0.5. Backoff 0.5 was chosen over
0.7 (70 sweeps) and 0.8 (69 sweeps) by measurement.

**The limit cycle did not corrupt any archived number, and that is the answer to
the question the finding raised.** Comparing the iterate-200 state against the
rolled-back best iterate and against every converged relaxation setting, on the
same equilibrium:

| quantity | spread across all 12 runs | rollback vs iterate 200 |
|---|---|---|
| bush-side wall temperature `T_bearing`, full midplane profile | <= 1e-3 K | 0 K (identical) |
| `T_bearing` peak | 46.711-46.712 C | identical |
| W | 4002.37-4002.48 N | 0.01 N |
| p_max | 1405340-1405380 Pa | 0 Pa |
| T_max | 320.154-320.155 K | 0 K |
| h_min, cavitated fraction | identical to all printed digits | identical |

The **-2.3 K CHT peak residual stands unchanged**; so do mean +0.02 K and rms
1.17 K. The one mobile quantity is the *argmax cell* of a nearly flat maximum,
which hops +-1 cell (0.5 deg). No figure needed regenerating.

New config keys, all inert by default: `coupling_keep_best_iterate` (true;
rollback only ever fires on cap exhaustion), `coupling_relaxation_mode`
(FIXED/ADAPTIVE), `coupling_relaxation_min`, `coupling_relaxation_backoff`,
`coupling_relaxation_growth`, `log_coupling_history`. `coupling_max_iters`,
`coupling_tolerance` and `coupling_relaxation` are now validated.

`cases/01_ferron_conjugate/force_balance.py` and
`cases/01_ferron_wallbc_fields/force_balance.py` take the coupling schedule from
`PANCAKE_COUPLING_MAX_ITERS`, `PANCAKE_COUPLING_RELAXATION` and
`PANCAKE_COUPLING_EXTRA`, all defaulting to the frozen protocol so an unset
environment writes a byte-identical config.

Tests: `test_coupling` gains a monotone-convergence case (best == last, no
rollback), a synthetic oscillating history proving the best iterate is returned
and the run is called NON-CONVERGED, a slow-but-monotone case needing no
rollback, an adaptive-schedule check and an exact `Snapshot` round-trip including
ghosts and `old_data`. 34/34 green on Windows (MinGW) and on afml-r6 (GCC).

**Reruns (24 pinned ranks on afml-r6, one job at a time).** `rerun_cht_robin10`
under the frozen protocol reproduces its archive to 9 significant figures
(eps 0.603313784 -> 0.603313785, psi to 5.3e-6 deg) with the rollback firing on
all eight calls; scored against the 16 Fig. 4 points the archive, the rolled-back
best iterate and a fully converged ADAPTIVE re-closure are **identical to three
decimals** (mean +0.017 K, rms 1.173 K, peak 46.712 C, residual -2.33 K).
`outer_fixed313` converges anyway and reproduces its archived -2.971 / 3.581
exactly. Fully converged, the same batch costs 388 outer sweeps instead of 1600.

**New result out of the same batch:** `01_ferron_conjugate/outer_robin80` - the
h = 80 variant that had been abandoned when `fb_iter_06` aborted on the
rank-dependent KSP divergence - now closes: eps = 0.602770, psi = -37.7046, mean
error **-0.20 K**, rms **1.21 K**, peak 46.450 C. That is within 0.26 K of the
h = 10 production group and 2.8 K better than the -1.31 K / 1.76 K on the books
for it, which came from a never-converged near-equilibrium state. Closing the
force balance, not the outer coefficient, is what supplies the heat. Re-closed
under ADAPTIVE it lands on the same equilibrium to six digits from 501 sweeps
instead of 2178.

Full report, residual traces and the run-by-run comparison:
`cases/01_ferron_wallbc_fields/PICARD_FIX_REPORT.md`;
figure `cases/01_ferron_wallbc_fields/ferron_picard_convergence.png`.
Physics: `PHYSICS.md` §7.4 (new) and §13.1.2 (rewritten - it previously said the
limit cycle was "reported rather than tuned away", which is now superseded).
Case docs: `01_ferron_wallbc_fields/{README,SUMMARY}.md`,
`01_ferron_conjugate/{README,SUMMARY}.md` (the abandoned-`outer_robin80` notes
are now closed), `cases/INDEX.md`.

**Operational note recorded because it cost an hour.** afml-r6 was shared with
another session's 24-rank batch. Three 24-rank jobs on 56 cores is not a 1.3x
slowdown but a **45x** one under MPI busy-wait (692 s against ~99 s for the same
call). Lanes are serialized here; one 24-rank job at a time is the policy.

### 2026-08-05 - Validation cases consolidated into one place, with an index

Documentation and directory structure only; no solver code, no runs, no results
changed. Every case, probe and diagnostic run for the validation campaign now
lives under `cases/`, and `cases/INDEX.md` is the front door: one table of every
case directory, its campaign stage, its status, its figures, and whether the bulk
data is local, on afml-r6, or both.

**Strays folded in.** Four run trees lived outside `cases/` on afml-r6, three of
them in `/tmp` or the home directory - i.e. in the two places most likely to be
cleaned - while being the evidence for claims made in this file and in case
summaries. All were moved (never copied-and-deleted, never overwritten):

| was | now |
|---|---|
| `~/cases_run/01_ferron_conjugate` (24 ranks, 421 files) | `~/pancake/cases/01_ferron_conjugate/r6_24rank_rerun/` |
| `/tmp/ksp_repro` (153 files) | `~/pancake/cases/90_diagnostics/ksp_rank_robustness/` |
| `~/diag_h10` | `~/pancake/cases/01_ferron_conjugate/diag_h10_frozen_geometry/` |
| `~/gate_fixed313` | `~/pancake/cases/01_ferron_conjugate/gate_fixed313/` |
| `~/smoke` | `~/pancake/cases/01_ferron_wallbc_fields/smoke_coarse_180x24/` |

The conjugate rerun is **not** a duplicate of the synced 4-rank archive: its
`outer_fixed313` leg reproduces it to every printed digit (eps = 0.593397,
psi = -37.5036), which is the rank-independence evidence, and its `outer_robin80`
leg aborted three force-balance iterates earlier - the failure that started the
KSP investigation. Both archives are kept whole.

**New stage 90.** `cases/90_diagnostics/` holds work that explains or defends
solver behaviour rather than physics, so a validation case never has to carry a
solver postmortem. `ksp_rank_robustness/` is the minimal reproducer behind the
ASM + fallback-ladder fix below, with the six run configs and full solver logs
mirrored locally and the 137 MB of VTK left on afml-r6. `remote_harness/` holds
reference copies of the loose afml-r6 driver and analysis scripts, each mapped to
the case it drove.

**No directory was renamed.** The existing `00/00c/01/02` prefixes already encode
the campaign stage and are the vocabulary `PLAN.md`, `VALIDATION_PLAN.md`, this
file and every case summary already use; about sixty archived `config.txt` and
`run.log` files name their own case in a header comment, so renaming would leave
the archives self-inconsistent for no functional gain. The mapping to the
stage scheme is stated in `cases/INDEX.md` instead. Because nothing moved, the
inter-case Python references (`02_static_locus` -> `01_ferron_prelim`,
`01_ferron_vs_reynsim` -> both, `01_ferron_conjugate` ->
`01_ferron_thd_crossfilm`) are untouched; all 20 analysis and plot scripts under
`cases/` were run standalone as a check and all 20 exit 0.

Two provenance findings came out of that check, both recorded in the index rather
than papered over: `00_grid_dt_independence/dt_fbns_analysis.json` records
`end_t = 0.15` for T2/T3/T4 while the configs now on disk say `0.06` (the JSON is
the record of the run as executed and was left alone), and
`01_ferron_conjugate/conjugate_analysis.json` was missing and has been
regenerated from its own archived run directories.

New documentation: `cases/INDEX.md`; `README.md` for `01_ferron_prelim`,
`01_ferron_vs_reynsim`, `01_ferron_wallbc_fields`, `02_static_locus` and both
stage-90 directories; `SUMMARY.md` for `01_ferron_conjugate` (its README ended in
a dangling pointer to a summary that was never written).

### 2026-08-05 - Rank-robust linear solves: ASM preconditioning and a fallback ladder

Linear solves diverged as a function of the MPI rank count. The same
configuration converged at 1, 2 and 6 ranks and aborted at 4 with
`DIVERGED_DTOL`; the conjugate-bush campaign lost `outer_robin80` at force-balance
iterate 6 on 4 ranks and at iterate 3 on 24, and `tests/test_fbns` aborted under
`mpiexec -n 4` while passing the `-n 2` and `-n 6` variants ctest actually ran.

**Root cause.** The mesh is decomposed in theta only, so every rank boundary cuts
the flow direction. The two operators that failed are the ones with essentially
no diffusion to fall back on: in `STEADY_STATE` the energy equation carries no
time-derivative term, leaving convection against a conduction coefficient
weighted by the liquid film depth (`k * theta * h ~ 1e-5` against a convective
flux of order `rho_cp * u * h * dz`), and the Elrod film-content equation is
hyperbolic by construction. Both therefore have cell Peclet numbers of order
1e4-1e6, and their inverse is global along theta. `PCBJACOBI` discards the
coupling across every rank boundary, so the preconditioner truncates exactly the
direction the solution propagates along, and its quality degrades as the slabs
get thinner. BiCGStab's non-monotone residual then turns a poor preconditioner
into a divergence rather than a slowdown. Nothing about this is stochastic: each
failure is deterministic for a given decomposition, which is why it looked like a
rank lottery from the outside.

**Fix, part 1 - additive Schwarz.** The default preconditioner is now
`PCASM` with one cell of overlap and the PETSc default ILU(0) subdomain solves
(`linear_pc_type`, `linear_pc_overlap`). The overlap restores the theta coupling
across rank boundaries. Measured on the Ferron THD force-balance iterate
(`cases/01_ferron_vs_reynsim/diag`, 720 x 192, `ENERGY_EQUATION` + FBNS), total
Krylov iterations over the whole run:

| ranks | block Jacobi | ASM(1) |
|-------|--------------|--------|
| 1     | 3723         | 3723   |
| 2     | 4010         | 3755   |
| 4     | **diverged** | 3761   |
| 6     | 3961         | 3735   |

The iteration count is now flat in rank count instead of drifting upward and
falling over. On the isothermal steady G3-FBNS case the same measurement goes
from 416/454/439/452 to 416/409/412/420.

**On one rank ASM is not used.** With a single rank additive Schwarz degenerates
to one block covering the whole matrix, which *is* block Jacobi — but PETSc
reaches the subdomain ILU by a different route (`MatCreateSubMatrices` rather
than the diagonal block) and the resulting ILU(0) factors are not identical. On
this operator the ASM route is measurably the weaker one: the serial THD iterate
took 2.8x longer and needed 25 fallback solves. `LinearSystem` therefore selects
`PCBJACOBI` when `mesh.size == 1`, which also keeps serial runs bit-for-bit on
the pre-fix Krylov path.

**Fix, part 2 - a fallback ladder.** Removing the rank dependence is necessary
but not sufficient, and the 24-rank runs show it directly. Both the FBNS Schur
complement and the Elrod film-content operator go near-singular whenever a
complementarity pivot sits on the fence: a last-bit difference flips the pivot
and leaves one near-singular row behind, and which Krylov path survives that is
not predictable from the rank count. A failed solve is now retried on
`linear_solver_fallbacks` (default 2) escalation rungs - BiCGStab(ell) with
ASM(2)/ILU(1), then BiCGStab(ell) with block Jacobi/ILU(1) - before the run
aborts. The rungs are deliberately *different* paths rather than tighter versions
of one path, because that is what the evidence supports: the system that defeated
ASM in serial was solved in 15 iterations by block Jacobi, and vice versa at 4
ranks. Retries cost nothing on a solve that converges, and every retry and
recovery is logged on rank 0 so a production log shows when the ladder carried a
run. `linear_solver_fallbacks = 0` restores abort-on-first-failure.

**BiCGStab(ell) is not made the default, and the measurement is why.** On
`outer_robin80` fb_iter_03 at 24 ranks the one temperature solve that diverges
does so after exactly 25 iterations under *both* block Jacobi and ASM(1) - the
preconditioner does not decide that solve's fate, the Krylov method does, and
`linear_ksp_type = BCGSL` removes the divergence outright. But paying for
BiCGStab(ell) on all ~450 solves to avoid one costs more than the wasted 25
iterations: against the pre-fix baseline, `BCGSL` runs the Ferron THD iterate at
+9.4% / +16.1% / +26.7% wall time at 1 / 2 / 6 ranks and G3-FBNS at +7.7% to
+14.4%, so at 6 ranks it leaves the 20% budget. The ladder buys the same
robustness at the cost of the one solve that needs it. `linear_ksp_type = BCGSL`
remains available for a run that would rather not lose a solve at all.

**Verification.**

| case | ranks | before | after |
|------|-------|--------|-------|
| `01_ferron_vs_reynsim/diag` | 1 / 2 / 4 / 6 | OK / OK / **diverged** / OK | OK / OK / OK / OK |
| `01_ferron_conjugate/outer_robin80` fb_iter_03 | 4 | OK | OK |
| `01_ferron_conjugate/outer_robin80` fb_iter_06 | 4 | **diverged** | OK |
| `01_ferron_conjugate/outer_robin80` fb_iter_03 | 24 (afml-r6) | **diverged** | OK, one recovery |
| `01_ferron_conjugate/outer_robin80` fb_iter_06 | 24 (afml-r6) | OK | OK |
| `tests/test_fbns` | 1 / 2 / 4 / 6 | OK / OK / **diverged** / OK | OK / OK / OK / OK |

No fallback rung fired in any of the local runs: there the ASM change is what
carries them. The 24-rank fb_iter_03 case is the exception and the one that
justifies the ladder - the temperature solve still diverged after 25 iterations
under ASM, and rung 1 (BiCGStab(ell) + ASM(2)/ILU(1)) finished it in 24. Krylov
iteration totals at 24 ranks: 43909 with ASM against 47636 with block Jacobi on
fb_iter_03, 43689 against 47374 on fb_iter_06 - ASM is 8% cheaper in iterations
at that width. Wall time 51 s against 49 s and 52 s against 49 s.

The 24-rank runs were made on afml-r6 with `OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1`, both binaries built from the same tree at c3f7b09 and at
this commit. The pre-fix binary aborts fb_iter_03 with `linear solve
'temperature' failed after 25 iterations: KSP diverged (DIVERGED_DTOL, reason
-4)`, which is the failure the conjugate-bush campaign reported. Every integral
metric of fb_iter_06 at 24 ranks - including `liquid_residual` - is identical
between the two binaries to all printed digits.

**Solutions are unchanged.** Serial output is byte-identical to the previous
binary. At 2 and 6 ranks every integral metric of the THD iterate - `W`, `Fx`,
`Fy`, `force_angle_deg`, `p_max`, `T_max`, `T_max_theta_deg`, `cav_frac`,
`film_content_min`, `h_min`, `mu_min/max`, `friction_torque`, `q_supply` - agrees
to all printed digits, and the fields agree to one unit in the last digit of the
six-significant-figure VTK output (1e-3 K on temperature, 10 Pa on pressure,
1e-6 on film content). Global liquid mass agrees to 1.4e-11 relative. The same
holds for the conjugate `fb_iter_03` case at 4 ranks. Byte-identity is not
expected in parallel: the Krylov path changed, so the two solutions differ within
the solver tolerance, and that is what the numbers show.

**Cost.** Wall time on healthy cases, worst case +12.8%, inside the 20% budget:

| case | ranks | before | after | delta |
|------|-------|--------|-------|-------|
| G3-FBNS isothermal | 1 / 2 / 4 / 6 | 1059 / 598 / 465 / 428 ms | 1027 / 606 / 456 / 444 ms | -3.0% / +1.3% / -1.9% / +3.7% |
| Ferron THD iterate | 1 / 2 / 6 | 54.9 / 32.2 / 17.5 s | 54.1 / 33.8 / 19.6 s | -1.4% / +5.0% / +11.8% |
| conjugate fb_iter_03 | 4 | 244.8 s | 276.1 s | +12.8% |

ASM needs *fewer* Krylov iterations (40151 against 42329 on the conjugate case)
but each one costs more - the overlap scatter and the larger subdomain factor.
That trade is the price of rank-independence.

Both platforms are green: 34/34 on Windows (MSYS2 MINGW64, MS-MPI) and 34/34 on
Linux (afml-r6, OpenMPI).

**Test coverage.** ctest ran the MPI tests only at `-n 2` and `-n 6`, which is
why a defect that fires at `-n 4` survived several releases. Added `-n 1` and
`-n 4` variants of `test_fbns`, and `-n 4` variants of `test_schemes`,
`test_energy` and `test_elrod_a2`. Suite is 34 tests, green.

`test_energy`'s two analytic wall-balance identities asserted 1e-6 K on a 305 K
field - 3e-9 relative - while asking the solver for `linear_rtol = 1e-8`. How far
inside its tolerance a solve lands is a property of the preconditioner, so the
assertion was measuring the Krylov path rather than the balance it names. Both
now set `linear_rtol = 1e-13`, matching what `test_diagnostics`, `test_fbns` and
`test_mass_residual` already do for tight identities.

### 2026-08-05 - Conjugate bush: optional 2-D (r, theta) conduction in the annulus

The conjugate bush closed **locally**: per cell, 1-D radial conduction through
the wall in series with the outer surface, `q_bush = G (T_b - T_amb)`. Every watt
entering the bore therefore left through the wall directly underneath it, so the
bush surface temperature could only peak where the film was hottest. On the
Ferron rig that put the predicted bush peak at 347.7 deg from the groove against
a measured 258.2 deg, and case 01-THD-crossfilm identified the missing mechanism
by name: circumferential redistribution inside the bush.

`bush_thermal_model = CONDUCTION_2D` (default `LOCAL`) now solves Ferron's own
bush equation instead — steady `div(k_b grad T) = 0` in the annulus
`[R, bush_outer_radius] x [0, 2pi)`, z-collapsed exactly as the paper's is.

**The coupling stays linear and non-iterative.** The cross-film profile makes the
film-side flux affine in the bore temperature,
`q_bush = A0 + A1 T_shaft - B T_bush` with `B = (k/h)(1 + b_slope/2) > 0`, which
is the same relation the local closure inverts against `G (T_b - T_amb)`. Here it
is imposed as a Robin condition on the bore, so the annulus is a linear system;
it is solved for both right-hand sides on one factorisation, giving
`T_b(theta) = alpha + beta T_shaft` exactly, and the global shaft balance closes
`T_shaft` in one step exactly as before.

**New module `src/bush_conduction.{hpp,cpp}`.** Deliberately free of MPI and of
`Fields`: it takes the z-averaged flux law on the global theta grid and returns
the bore temperature there, which makes it unit-testable on its own.

- Finite volume, uniform in theta on the solver's own grid, `bush_radial_cells`
  (default 16) cells in r. Every radial face uses the **exact logarithmic shell
  conductance**, so a purely radial problem is reproduced to round-off at *any*
  radial resolution — measured 3e-13 relative at n_r = 1, 4, 16 and 64.
- The bore condition is put in series with the bore-to-first-node half cell,
  which removes `T_b` from the linear system without losing accuracy.
- **Direct solve, no tolerance to tune.** Block tridiagonal in theta with
  n_r x n_r blocks and one periodic wrap: block Thomas on blocks 0..n-2 carrying
  n_r extra RHS columns, then one dense n_r x n_r Schur solve closes the wrap.
  About 10^7 flops at 720 x 16. An alternating line-relaxation solver was written
  first and **rejected**: the bush grid is near-isotropic (radial and
  circumferential face spacings within a factor of a few), and it stalled at
  1.5e-5 relative residual after 500 sweeps on the very first coupled test.
- Serial and deterministic. The three coupling arrays are all-reduced onto the
  global theta grid and every rank then runs the same arithmetic on the same
  input, so there is no broadcast and no rank-0 asymmetry to get wrong. Covered
  at `-n 2` and `-n 4`.

**Outer boundary, for both bush models.** `bush_outer_boundary = CONVECTION`
(default) is the Robin condition with `bush_outer_heat_transfer`, with `0` an
adiabatic outer surface; `FIXED_TEMPERATURE` holds that surface at
`thermal_ambient_temperature`. The second is the `h -> infinity` limit of the
same series resistance and is written as such in both closures, so a fixed-wall
study no longer has to be expressed as a large-HTC limit.

**The Newton linearisation of the film energy equation takes the bush as
frozen.** `dT_b/dT_bar` is nonlocal under `CONDUCTION_2D` — one film cell barely
moves a bush its neighbours also feed — so `dT_b/dT_bar = 0`, the fully-sunk-wall
coefficient. Two alternatives were written and measured before settling on it:
the bore admittance with neighbours frozen reproduces the same number to 5 %
(`dT_b/dT_bar = 0.053` on the Ferron rig) but **grows without bound as
`bush_radial_cells` is refined**, which would make the film's damping depend on
the bush mesh; and the series-to-ambient conductance is mesh-free but so small on
a nearly adiabatic bush that the clamp zeroes the wall coefficient entirely, and
with the wall term fully explicit **the temperature KSP diverged** on the Ferron
rig at 720 x 192. Only the path to the fixed point is affected; the converged
answer always uses the solved `T_b`.

**`LOCAL` is untouched.** Verified, not asserted: `b74db94` and this commit were
both built and run on the same two cross-film configs (adiabatic walls, and the
fixed-wall limit), and the complete `results/` trees — flat VTK and
`diagnostics.csv` — are byte-identical.

Tests: new `tests/test_bush_conduction.cpp` (analytic annulus at four radial
resolutions and both outer conditions; exact reduction to the LOCAL closure when
nothing varies with theta; mirror symmetry; a localised hot patch spread
circumferentially — peak cut to 28 %, effective support widened 2.5x, bore
temperature smooth where the flux is a step; zero net heat with nonzero local
flux at an adiabatic outer surface; exactness of the `alpha + beta T_shaft`
split). `tests/test_cross_film.cpp` gains the end-to-end coupled checks
(z-uniformity of the bore temperature, the direct solve's residual, peak cut and
step smoothing against `LOCAL` on identical input, zero net bush heat when the
outer surface is adiabatic, and byte-identity of every field under `LOCAL`) and
now also runs at `-n 4`. 30/30 green.

Docs: PHYSICS.md §13.1.2.1 (equations, boundary conditions, discretisation,
solver, MPI, the z-collapse and the neglected end faces), README.md.

### 2026-08-05 - Film-surface temperatures and wall heat fluxes as optional output fields

The written `temperature` field is the depth average `T_bar`. Under
`cross_film_thermal_model = POLYNOMIAL` that is the mean of a profile whose two
ends differ, so it is not the quantity a wall thermocouple measures, and nothing
in the output exposed the wall heat flux at all. Four optional fields now do,
selected in `output_fields`:

| key | VTK name | definition | unit |
| --- | --- | --- | --- |
| `temperature_journal` | `T_journal` | `T(eta=0)`, shaft-side film surface | K |
| `temperature_bearing` | `T_bearing` | `T(eta=1)`, bush-side film surface | K |
| `heat_flux_journal` | `q_journal` | wall heat flux at `eta=0` | W/m^2 |
| `heat_flux_bearing` | `q_bearing` | wall heat flux at `eta=1` | W/m^2 |

`temperature` is unchanged — still `T_bar`, still absolute kelvin, as are the two
new temperatures.

**The fluxes are positive INTO the fluid.** That is the negative of the module's
internal `q_shaft = (k/h)(a+b)` and `q_bush = (k/h)(b-a)`, which are positive for
heat *leaving* the film (PHYSICS.md §13.1.2). A field named after a surface
should be positive into the fluid it belongs to, so the flip happens once, at the
write, and is stated in the field docs, PHYSICS.md §15.2, the config comment and
the GUI tooltip.

**Zero cost when unused.** The fields are allocated only when the model is active
*and* the key appears in `output_fields`; `CrossFilm::update` skips the stores
otherwise. A run that does not ask for them is byte-identical, which
`tests/test_cross_film.cpp` checks on the raw field storage, ghosts included.

**Asking for them without the model is an error, not zeros.** Under
`cross_film_thermal_model = NONE` there is no profile to evaluate, so
`SimulationConfig::validate` rejects the selection by name instead of writing a
plausible-looking field of zeros.

The temperatures are the profile evaluated at the walls rather than a copy of the
conjugate wall temperatures, so a future change to the shape functions cannot
leave the two disagreeing; the test asserts both the closed-form equality and the
agreement with `cross_film_t_shaft` / `cross_film_t_bush`. The flux check
compares against `-k dT/dy` taken by a 5-point stencil on the profile (exact for
a quartic), so the sign convention is verified independently of the expression
that produces it. Two limit cases pin the physics: an adiabatic bush gives an
identically zero `heat_flux_bearing` cell by cell, while an adiabatic shaft gives
a locally nonzero `heat_flux_journal` whose area integral vanishes — the upstream
heat path the conjugate shaft exists for.

`ctest` 27/27.

### 2026-08-05 - Cross-film generalized coefficients under FBNS, and one definition of the Reynolds coefficients

Two changes on one branch, because the second is what makes the first checkable.

**1. `cavitation_model = FBNS` now accepts `cross_film_thermal_model = POLYNOMIAL`.**

The validation guard that refused the combination is gone. The FBNS face
mobilities are now the generalized (Dowson) ones the Elrod path already used:
Poiseuille `rho_l F2` instead of `rho_l h^3/(12 mu)`, and Couette thickness
`2 F1/F0` instead of `h`. That is the point of the model — `mu(T(y))` acts
*inside* the film integrals rather than through a depth-averaged `mu`, and
because `mu` is convex in `T` while the core is the hottest part of the film, the
effective viscosity implied by `F2` is the lower one (PHYSICS.md §13.1.3). On the
test bearing `F2` differs from `h^3/(12 mu)` by 9.7 % and `2F1/F0` from `h` by
8.8 %.

The storage term keeps the geometric `h`: it measures a volume of liquid, not a
flux.

**The Newton structure is unchanged, deliberately.** `CrossFilm::update` rebuilds
`F0/F1/F2` once per outer (Picard) iteration, *before* `FBNS::solve` is entered,
so inside one Newton solve they are frozen coefficients in exactly the sense `mu`
already was. `A` and `B` are still built once at the top of the solve, the mass
balance is still linear in `(p, alpha)`, every nonlinearity still lives in the
Fischer-Burmeister equation, and the residual and the Jacobian read the *same*
frozen mobilities — which is the condition for the quadratic local rate. A
Jacobian carrying `d(F2)/dT` against a residual that did not would be the limping
case and is not what was done.

**2. `src/reynolds_coefficients.hpp` — one definition of the operator coefficients.**

`Reynolds::solve_elrod`, `Reynolds::solve`, `FBNS::solve` and `Diagnostics` had
four local copies of the same coefficients. `Diagnostics` *measures* the other
three, so a drifted copy there is not a shortcut but a wrong number, and it had
drifted:

| quantity | assembly used | diagnostics used |
|---|---|---|
| Elrod `Gamma_base` | `rho_0 beta_bar F2` | `rho_0 * cfg.bulk_modulus * h^3/(12 mu)` |
| Couette face flux | `2F1/F0` face average | plain `h` face average |
| Elrod switch function | `theta >= VoidParams::theta_full` | `theta >= 1` |
| FBNS boundary mobility | `rho_l F2` | `rho_l h^3/(12 mu)` |

So on exactly the cross-film and `gas_pressure_coupling = VOID_COUPLED` runs that
`mass_residual` exists to police, it reported a fabricated imbalance of the size
of those differences. All four call sites now assemble from the shared header.

**Changed**

- `src/reynolds_coefficients.hpp` (new): `poiseuille` / `poiseuille_coefficient`,
  `couette_thickness`, `couette_face_flux`, `fbns_mobility`, `elrod_gamma_base`,
  and the `VoidParams` closure (`void_params`, `gas_inventory_plateau`) moved out
  of `reynolds.cpp`'s anonymous namespace.
- `src/fbns.cpp`: `poiseuille_gamma` and `couette_face_flux` delegate to the
  shared header, so the generalized branch reaches the Jacobian blocks, the
  residual `cell_outflow`, the cold-start solve and the axial Dirichlet links
  together.
- `src/diagnostics.cpp`: local `gamma_base`, `couette_face_flux` and the two
  boundary-mobility lambdas replaced by the shared definitions; the Elrod switch
  function now reads `VoidParams::theta_full`.
- `src/reynolds.cpp`: local copies removed; `Gamma_base` and the Couette face
  flux built from the shared header. Bit-identical (same expressions, same
  floating-point association).
- `src/config.hpp`: the `FBNS requires cross_film_thermal_model = NONE` error is
  removed, with a comment recording why the combination is now sound.
- `src/fbns.hpp`, `PHYSICS.md` §7.3.8 / new §7.3.9 / §13.1.3 / §15.1, `PLAN.md`.
- `tests/test_fbns.cpp`: the scope test now asserts the combination is *accepted*;
  `ThermalCase` registers the cross-film fields and `solve_coupled` calls
  `CrossFilm::update` with the same ghost exchange `main.cpp` performs; new test 14.

**Verified**

- `ctest` 27/27 (`mpiexec -n 2` and `-n 6`).
- Test 14: the generalized coefficients genuinely engage (9.7 % / 8.8 %, 4.5 K
  core-to-mean rise); FBNS and Elrod driven through the identical Picard loop with
  the identical conjugate-wall closure agree to 0.48 K over 454 shared full-film
  cells, and the conjugate bush surface temperature to 0.50 K — within the
  1.05x of the depth-averaged difference it is affine in; per-cell mass residual
  `2.3e-10` of the peak Couette flux over every non-groove cell, `8.8e-12` on the
  16 rupture/reformation front cells.
- **Cross-film off is byte-identical to `6dacec0`** for FBNS: the Ferron THD
  Fig. 4 configuration (720 x 192, `mpiexec -n 4`, `ENERGY_EQUATION` +
  `POLYNOMIAL` viscosity + fed groove) reproduces all four `.vts` pieces and
  `diagnostics.csv` byte for byte, with `mass_residual` added to the output field
  set so the diagnostics change is inside the comparison.

**Not changed, and still true**

The §13.1.2 approximations of the cross-film profile now apply to FBNS runs as
well: the midplane collocation is Couette-only (the pressure-gradient term was
removed because it limit-cycled against the rupture position), and the convective
term is saturated against the dissipative one, so where that saturation binds the
quartic degenerates towards a fixed flattened shape rather than resolving the
wall layers. Quote any cross-film temperature together with that caveat.

### 2026-08-05 - Ferron polynomial viscosity (`liquid_viscosity_model = POLYNOMIAL`)

A third selectable liquid-viscosity law, the quadratic used by Ferron, Frene &
Boncompain, *ASME J. Lubr. Technol.* **105** (1983) 422:

```
mu(T) = mu * (k0 - k1*Tbar + k2*Tbar^2),   Tbar = T[degC] / T_ref[degC]
```

Motivation: the shipped benchmark bearing in `tests/test_fbns.cpp` already *is*
Ferron's (R = 50 mm, L = 80 mm, c = 152 um, mu = 0.0277 Pa s, rho = 860), and
comparing against that paper's THD results needs that paper's viscosity law, not
an Andrade fit to it.

**The units question, settled from the paper text.** `Tbar` is a ratio of
**CELSIUS** temperatures. The paper's nomenclature reads "T = temperature degC"
and "Tbar = T/T_0 nondimensional temperature", and its Table 1 coefficients only
work that way:

| reading | mubar(40 C) | mubar(80 C) | mu(80 C) |
|---|---|---|---|
| Celsius ratio (correct) | `3.287 - 3.064 + 0.777 = 1.000` | `0.267` | 7.4 mPa s (measured VG32) |
| Kelvin ratio | 1.000 (trivially) | `0.820` | 22.7 mPa s (not an oil) |

`property_reference_temperature` stays in Kelvin like every other temperature in
the config; the conversion is explicit inside the law
(`SimulationConfig::polynomial_reference_temperature_celsius`), so there is
exactly one place where 273.15 appears.

**Changed**

- `src/config.hpp`: `LiquidViscosityModel::POLYNOMIAL` (tokens `POLYNOMIAL`,
  `POLY`, `FERRON`); new keys `viscosity_poly_k0` / `k1` / `k2`, defaults
  `1, 0, 0` (the identity `mu(T) = mu`, so the model is inert unless configured);
  `polynomial_reference_temperature_celsius()`; config writer emits the three
  keys with the Celsius caveat in a comment.
- `src/fluid_properties.cpp`: the `POLYNOMIAL` branch of
  `liquid_solution_viscosity`. A **pure temperature law** — it deliberately
  applies neither the Barus pressure factor nor the dissolved-gas factor, which
  belong to `EMPIRICAL`. `mu` is floored at `1e-3 * mu` with a single
  rank-0 warning per run, because a polynomial extrapolated past its fit crosses
  zero (a `k2 = 0` set always does) and a negative `mu` makes the Reynolds
  operator indefinite.
- `src/gui/schema.cpp`, `src/gui/panels_setup.cpp`: the three coefficients and
  the new model entry, with the Celsius convention in the tooltip.

**Validation** (`validate()`): errors if
`property_reference_temperature <= 273.15 K` (the Celsius ratio has no usable
normalisation at or below 0 degC) or if any coefficient is non-finite. Warns —
does not error — if `k0 - k1 + k2` deviates from 1 by more than 5 %, since `mu`
is defined *at* the reference temperature and the polynomial should return 1
there; and warns if `k1 = k2 = 0`, which is `CONSTANT` by another name.

**Tests**

- `tests/test_fluid_properties.cpp`: `test_ferron_polynomial_viscosity` pins
  `mu(40 C) = 0.0277` and `mu(80 C) = 0.0074` on Ferron's coefficients, asserts
  the Kelvin reading is *not* what the code produces, checks monotone thinning
  40 -> 80 C, checks that pressure and dissolved gas do not move the result, and
  checks the negative-viscosity floor. `test_ferron_polynomial_validation`
  covers the four validation paths plus the config round-trip.
- `tests/test_fbns.cpp`: `test_polynomial_viscosity_couples` runs the coupled
  FBNS+`ENERGY_EQUATION` and Elrod+`ENERGY_EQUATION` Picard loops on Ferron's own
  bearing and coefficients. Both converge; the polynomial moves `mu` by 44 %
  (FBNS) and 48 % (Elrod) of `mu_ref`, and the pressure moves away from the same
  run with `mu` frozen by 2.2e-1 and 1.3e0 of `p_ref` respectively — i.e. the
  `mu` field reaches the Reynolds operator on both paths. `POLYNOMIAL` and
  `EMPIRICAL` differ by 3.0e-2 `mu_ref` and 1.2e-2 `p_ref`, so they are not the
  same model in disguise.

Full suite: 27/27 green (MSYS2 MinGW Release, MS-MPI, `-n 2` and `-n 6`).

### 2026-08-05 - FBNS phase 3: gaseous cavitation (`FBNS` + `TWO_PHASE`)

`cavitation_model = FBNS` now runs with `fluid_property_model = TWO_PHASE`. This
is step 2 of `docs/FBNS_GAS_COUPLING.md`, on the shared `FullFilm` predicate that
landed in 28ac88a.

**The model.** A resolved free-gas phase occupies a fraction `alpha_gas` of the
gap, so the liquid can only fill what is left and the complementarity bound moves
from `1 - theta >= 0` to `(1 - alpha_gas) - theta >= 0`. `alpha_gas` is lagged —
the gas path is operator-split and runs after the flow solve — so within one
Newton solve it is a coefficient, the bound stays affine in `theta` with the same
derivative `-1`, and **both Jacobian blocks, the column swap, the
`cond(A_F) <= 4 + 2 sqrt(2)` bound and the Schur reduction are unchanged**. See
PHYSICS.md §7.3.7.

**Changed**

- `src/full_film.hpp`: `bound()` returns `max(1 - alpha_gas, 0)`, via two new
  helpers — `gas_displaces_liquid(cfg)` (the model gate) and
  `gas_volume_fraction(fields, cfg, i, j)` (the clamped free-gas fraction). Every
  full-film test in the solver reads it, so the ceiling moved in all seventeen
  places at once. No call-site edits.
- `src/fbns.cpp`: the cavity fraction becomes `alpha = theta_full - theta` in the
  residual, the update and the clamp ceiling; the cold-start fill writes
  `theta_full` per cell.
- `src/config.hpp`: `elrod_boundary_outflow` and `is_reformation_boundary` now take
  the ceiling as an explicit parameter (no default) instead of comparing against a
  literal `1.0`. These are the two sites the step-1 refactor could not route,
  because they take film-content values and live in the header `full_film.hpp`
  includes. Seven call sites in `reynolds.cpp`, `diagnostics.cpp`, `energy.cpp` and
  `gas_transport.cpp` pass `FullFilm::bound(...)`. A no-op today; correct wiring
  for whenever the Elrod ceiling merges.
- `src/fluid_properties.cpp`: `gas_state_relaxation` under-relaxes the free-gas
  state whenever the gas displaces liquid, not only under `VOID_COUPLED`. Same
  lever (`void_coupling_relaxation`, default 0.3), same reason — the
  `alpha_gas <-> theta` feedback limit-cycles without it.
- `src/main.cpp`: `alpha_gas` ghost cells are synced after the gas update under
  `TWO_PHASE`, so the energy, velocity and diagnostics passes in the same outer
  iteration evaluate the ceiling on a current field rather than a stale one.

**Guards lifted**: `fluid_property_model = TWO_PHASE` under FBNS, and the
converse guard that pinned TWO_PHASE to `cavitation_model = JFO`. Gas transport
and the gas state machine needed no changes — they were already model-agnostic.

**Guard kept**: `gas_pressure_coupling = VOID_COUPLED` under FBNS. It makes the
cavitation datum a per-cell field (`p_void`) that is advanced and under-relaxed
**only** inside `Reynolds::solve_elrod`'s outer flag-update loop. FBNS has no such
loop, so it would run the plateau frozen at `p_cav` while reporting
`VOID_COUPLED`. The error message says exactly that. `cross_film_thermal_model`
and `INLET_OUTLET` stay rejected for their existing, unrelated reasons.

**`p_cav` is not a field on this path, and does not need to be.** Under
`cavitation_threshold = LOCAL_PSAT` the datum is evaluated at the *feed* state
(`dissolved_gas_initial`, `temperature_initial` — both config scalars), so it is
uniform over the domain by construction, deliberately, so that it injects no net
force. The Fischer-Burmeister argument `a = (p - p_cav)/p_ref` already reads the
correct local value. `p_void` under `VOID_COUPLED` is the only genuinely per-cell
cavitation pressure in the codebase, and it is guarded off (above).

**Bug found and fixed while landing this.** The memo's §2 says the cold-start test
generalizes to `theta >= theta_full`. It must not. That test is a *sentinel*
comparison against the value `main.cpp` initialises `theta` to, not a full-film
test: an uncavitated gassy cell sits at `theta = 1 - alpha_gas` by definition, so
generalizing it declares a cold start on every gassy step, discards the warm start
and refills `theta` down to the new ceiling in one go. Measured effect before the
fix: pressure to 2.9e8 Pa on `config_r290_pz68_quick.txt`. The literal `1` is kept
and the reason is written at the site.

**Tests** (`tests/test_fbns.cpp` 11-13, 27/27 suite green including the `-n 6`
variants):

| Check | Result |
|---|---|
| complementarity vs the moved bound, 5 transient steps on the 00c bearing | `max(-(p-p_cav)/p_ref)`, `max(theta-(1-alpha_gas))`, `max` product all **0 exactly** |
| per-cell mass residual / peak Couette flux, away from the groove | `5.1e-15` |
| free gas released, FBNS vs JFO | `alpha_gas` max `2.50e-2` vs `8.44e-2` (same order; asserted within 10x, fronts logged only) |
| `alpha_gas = 0` vs phase-2 `SINGLE_PHASE`, 3 steps | bit-identical `p` and `theta` |
| `FullFilm::bound` over 24 model combinations | exactly `1.0` except FBNS + `TWO_PHASE` |
| Elrod / Gumbel vs 28ac88a, 6 configs | **54/54 output files byte-identical** |

The regression set covers JFO+TWO_PHASE, JFO+VOID_COUPLED+LOCAL_PSAT, JFO
isothermal, JFO+SINGLE_PHASE+STRIATED, GUMBEL and FBNS+SINGLE_PHASE.

**Known stiffness, reported not hidden.** `config_r290_pz68_quick.txt` with
`cavitation_model = FBNS` converges every Newton solve (61/61 steps, 2-11
iterations, mean 6.3, 3.9 s wall) and keeps the gas mass bounded
(3.99e-5 -> 2.76e-5 kg, monotone out through the axial ends), but its load history
is transient noise for tens of steps. The feed enters at its bubble point and
flashes enough propane to fill ~23% of the gap per step; the ceiling drops by that
much at once and the load swings between 5.8e2 N and 2.4e5 N on alternate steps,
decaying ~30% per step pair. The swing is **timestep-independent** (identical at
dt = 5e-6, 5e-8, 5e-9 s) because the released fraction is itself proportional to
dt, so it is the model's displacement rate, not a resolution artefact. The Elrod
path with `VOID_COUPLED` does not resolve this case either — it drifts
monotonically instead of oscillating. Lower `void_coupling_relaxation` to damp it.

### 2026-08-05 - One shared full-film predicate (`FullFilm`, `src/full_film.hpp`)

Byte-identical refactor. The question "does this cell still hold an unbroken
liquid film?" was spelled as the literal `theta >= 1.0` (and its `< 1.0` and
face-pair variants) at fifteen sites across `reynolds.cpp`, `energy.cpp`,
`gas_transport.cpp` and `diagnostics.cpp`. All fifteen now call one predicate.

**Why now.** This is the stated prerequisite of `docs/FBNS_GAS_COUPLING.md` §4,
step 1. With a resolved free-gas phase the ceiling on the film content drops from
1 to `1 - alpha_gas`, and a ceiling that has to move in fifteen places moves
wrong: the instant any gas came out of solution every cell in the domain would
fail `theta >= 1` and every Poiseuille switch — momentum, enthalpy, dissolved-gas
transport — would turn off at once, silently, with a converged Newton residual.
Routing every user through `FullFilm::bound` makes that a one-line change.
The gas-coupled bound is deliberately *not* part of this change.

**New header `src/full_film.hpp`** (namespace `FullFilm`):

- `bound(fields, cfg, i, j)` — the ceiling. Returns exactly `1.0` today for every
  model, with the future one-line body written out in the docstring.
- `holds(theta_value, bound_value, tol = 0)` — value-level core, for callers that
  hold a film content rather than a cell index, or that need a tolerance band.
- `is_full_film(fields, cfg, i, j)` and `is_full_film_face(fields, cfg, ...)` —
  cell and face predicates; a face is full-film only if both cells it joins are.
- `gate` / `face_gate` — the same tests as the 1/0 switch `g` the flux
  expressions multiply through.

**Sites routed** (17, the 15 the memo lists plus two adjacent ones in the same
files): `reynolds.cpp` 7 — the four `calculate_velocities` face switches, the
`calculate_macroscopic_properties` shear switch, the `consistent_boundary_flux`
branch in `add_boundary_link`, and the Elrod flag-flip counter (which compared
against a bare `1.0` while the snap beside it already used
`VoidParams::theta_full`); `energy.cpp` 4, deleting the file-local
`full_film_switch` / `full_film_face_switch`; `gas_transport.cpp` 3, deleting its
own copy of `full_film_switch`; `diagnostics.cpp` 3 — the TYPE_DIFFERENCING face
test, both axial boundary flux gates, and the `cavitated_fraction` count.
`calculate_velocities` no longer binds a `theta` reference at all.

**Byte identity.** Five cases, each run on `mpiexec -n 2` from a clean worktree
build of `94edb45` and from this branch, sha256 over every `.vts`, `.pvts`,
`.pvd` and `diagnostics.csv`: 41 files, 41 identical. JFO steady
(`grid_G3`, 240x64, 45% of cells cavitated), FBNS steady (`grid_G3_fbns`), JFO +
`ENERGY_EQUATION`, JFO + `ENERGY_EQUATION` + `cavitated_film_model = STRIATED`
(the liquid-fraction weighting), and a transient `TWO_PHASE` JFO + energy run
(`config_r290_pz68_quick.txt`) — the only path that executes
`GasTransport::solve`, i.e. the three `gas_transport.cpp` sites. Full `ctest`
green, 27/27, unchanged from baseline.

**Deliberately left as literals**, none of them a per-cell full-film test:
`SimulationConfig::elrod_boundary_outflow` and `is_reformation_boundary`
(`config.hpp`) take a film-content *value*, and `config.hpp` cannot include
`full_film.hpp` because `full_film.hpp` includes `config.hpp`;
`FBNS::solve`'s cold-start detection reads a globally reduced minimum, and the
memo's §2 sketch changes it as part of the FBNS gas work itself; and
`EOS::dp_dtheta` is a value-level EOS derivative with no live caller in the
solver.

### 2026-08-04 - FBNS phase 2: energy equation and single-phase properties

Lifts the phase-1 scope guards that kept `cavitation_model = FBNS` isothermal and
constant-property, so the THD validation cases (01-THD, 04) can run on the solver
that actually resolves the free boundary. `FULL_SOMMERFELD`, `GUMBEL` and `JFO`
are untouched and verified byte-identical (sha256 of the flat VTK product and
`diagnostics.csv` on the `grid_G3_fbns` geometry at 120x16, before and after);
isothermal FBNS is likewise byte-identical to phase 1.

**What was already right.** No new coupling loop was needed: the STEADY_STATE
outer (Picard) loop in `main.cpp` already ran properties -> cavitation solve ->
properties -> gas -> velocities -> energy for every model, and `FBNS::solve`
already built its face mobilities from the per-cell `mu` and
`rho_liquid_solution` *fields*, not from `cfg.mu`. The guards were the only thing
stopping a thermal FBNS run. `theta` also carries identical semantics under both
mass-conserving models, so the D14 liquid-fraction weighting of heat generation,
stored energy, conduction and convected enthalpy under `cavitated_film_model =
STRIATED` applies to FBNS through the same code path, unchanged. Full-film FBNS
cells hold `theta = 1.0` exactly (the Newton update clamps the cavity fraction at
zero), so the `theta >= 1` full-film test needs no tolerance - checked on the
converged 120x16 field: 1100 of 1920 cells exactly 1.0, none in `(0.999, 1)`.

**Jacobian / viscosity consistency.** Frozen-coefficient Picard around Newton,
stated explicitly in `fbns.hpp` and PHYSICS.md §7.3.6. `FBNS::solve` reads `mu`
and `rho_liquid_solution` once per call, builds `A`/`B` from them and holds them
fixed for that Newton iteration; no `d(mu)/dT` or `d(mu)/dp` term enters the
Jacobian, which is exactly what keeps it consistent with the residual
`cell_outflow` evaluates (both read the same frozen mobility) and preserves the
quadratic rate and finite termination. The `mu(T,p)` feedback closes in the outer
loop. This matches what the Elrod path does. Tested exactly rather than by
inequality: replaying the converged coupled `mu` field through an isothermal FBNS
solve as a prescribed coefficient reproduces the coupled pressure to
`1.1e-9 * p_ref` and `theta` to `5.8e-9`, against a `mu` field that the Andrade
law moved by 41%. A silent fallback to the scalar `cfg.mu` would show up as a
~40% pressure difference.

**Energy equation, FBNS-specific (`src/energy.cpp`).** Two changes, both gated on
FBNS so the Elrod path is bit-for-bit unchanged:

- `solved_boundary_pressure` now uses `mass_conserving_cavitation()` instead of
  testing for JFO. Both models clamp the axial Dirichlet value at `p_cav`; FBNS
  was previously handed the raw configured value, so the boundary enthalpy flux
  was computed against a gradient the flow solve never saw.
- New `poiseuille_face_activity`: under FBNS the Poiseuille enthalpy flux is
  *not* gated by the full-film switch. Elrod carries a barotropic
  `p = p_cav + beta ln(theta)` through the ruptured zone, which is not the
  physical pressure, so it must gate. FBNS pins a ruptured cell at `p = p_cav` by
  complementarity, so its gradients are physical: between two ruptured cells the
  flux vanishes on its own, and at a rupture front or an open axial end it is
  precisely the reformation flux the FBNS mass balance conserved. Gating it left
  the energy equation convecting enthalpy on a different velocity field from the
  one carrying mass - the defect `consistent_boundary_flux` patches on the Elrod
  path, which FBNS gets for free. The switch is deliberately *kept* on the
  Poiseuille dissipation term (a striated-film modelling claim, not a
  conservation requirement, reaching only the single front cells), so
  `heat_generation` output is unchanged for isothermal FBNS runs.

**Config validation (`src/config.hpp`).**

- `TWO_PHASE` is the only `fluid_property_model` still refused under FBNS, with
  the reason named (the bound moves to `theta <= 1 - alpha_gas`).
  `gas_pressure_coupling`, `cross_film_thermal_model` and `INLET_OUTLET` remain
  refused for their original reasons.
- **Defect fixed:** `SINGLE_PHASE` no longer demands invented gas numbers. A
  pure-oil THD run has no dissolved species, so `dissolved_gas_henry_coeff = 0`
  and `solution_viscosity_gas_coeff = 0` are the *correct* inputs, and both were
  hard errors. They are now errors only when a gas species actually participates
  (`TWO_PHASE`, `cavitation_threshold = LOCAL_PSAT`, or
  `dissolved_gas_initial > 0`). `EMPIRICAL` viscosity with all three exponents
  zero now warns that it is identical to `CONSTANT`, which is the real mistake
  that check was reaching for.

**Verification** (`tests/test_fbns.cpp`, new tests 8-10; full suite 27/27 on 2 and
6 ranks). Effectively-1D grid, 720x4, Neumann ends, Andrade `E_mu = 3700` K:
`SINGLE_PHASE` with zero coefficients reproduces `CONSTANT` bit for bit; FBNS and
Elrod temperatures agree to 0.44 K over 465 shared full-film cells at
eccentricity 0.3; Picard counts 10 (FBNS) and 9 (Elrod); load capacity
1.653e4 N thermal vs 1.807e4 N isothermal.

The **fronts are not compared against Elrod**, and the brief's "fronts within 1
cell" acceptance line cannot be met by any implementation: Elrod's rupture front
sits 12-18 deg from the exact 1D solution on this bearing and that gap is its
barotropic smearing, not a grid error - it stays near 13 deg from n_theta = 240 to
720, and raising `bulk_modulus` to sharpen the switch makes it *worse* and stops
the Elrod outer loop converging (200+ iterations at beta = 1e10 against FBNS's
13). The exact benchmark of PHYSICS.md §7.3.5 is the only front reference worth
asserting against, and it is untouched. The residual `dT` is that misplaced front
leaking into Elrod's pressure field and hence into the Poiseuille share of the
enthalpy transport (~25% here); it scales with the thermal rise (0.24 K at
eccentricity 0.2, 0.44 at 0.3, 0.70 at 0.4, 1.42 at 0.6), which is why the
comparison is made at 0.3.

**Acceptance run** (`grid_G3_fbns` geometry, 240x64, `EMPIRICAL` viscosity,
`viscosity_temperature_coeff = 3700`, `temperature_reference =
property_reference_temperature = 313.15`, 2 ranks): converges in **9 Picard
iterations** to `R_p = 0`, `R_T = 2.9e-9`. `T` in [313.15, 325.9] K,
`W = 8.76 kN` against `9.96 kN` isothermal (-12.1%), friction torque 4.99 Nm
against 6.34 Nm - both in the sanity direction for thinner oil.

**Instability found, pre-existing, affects both cavitation models.** That run
requires the groove to be thermally fed. A groove that pins pressure but is left
thermally transparent (`feeds_fresh_oil = false`, which is what
`cases/00_grid_dt_independence/grid_G3_fbns/config.txt` has) diverges to
`T ~ 1e6 K` under a strong Andrade slope, under FBNS *and* under JFO, at every
`coupling_relaxation` from 1.0 down to 0.15. Mechanism: the Poiseuille
dissipation `h^3 |grad p|^2 / (12 mu)` scales as `1/mu`, and at the groove edge
`grad p` is held fixed by the pressure pin, so hotter -> thinner -> hotter is
unbounded. Giving the groove its supply temperature removes the artefact (a fed
groove is a reservoir, not a sheared film). Not fixed here: any fix changes the
Elrod path too and belongs in its own scoped change.

**Gaseous cavitation: memo, not code** (`docs/FBNS_GAS_COUPLING.md`). The design
question was whether `theta <= 1 - alpha_gas` composes with the Fischer-Burmeister
function. It does, cleanly: with `alpha_gas` lagged (which the existing
operator-split gas path already guarantees) the new argument is affine in `theta`
with the same derivative, so `B`, the column swap, the `cond(A_F)` bound and the
Schur complement are all unchanged and the whole change is ~20 lines of
`fbns.cpp`. The blocker is elsewhere: the codebase identifies "full film" by the
literal test `theta >= 1.0` at 15 sites across `reynolds.cpp`, `energy.cpp`,
`gas_transport.cpp` and `diagnostics.cpp`. Move the bound and every one of those
switches turns off at once the moment any gas comes out of solution - a silently
fully-cavitated bearing with a converged Newton residual. Lifting the guard needs
a shared `is_full_film` predicate landed first as a byte-identical refactor on its
own branch; the memo sets that out step by step.

### 2026-08-04 - FBNS complementarity cavitation (`cavitation_model = FBNS`)

New module `src/fbns.{hpp,cpp}` implementing the Fischer-Burmeister-Newton-Schur
algorithm of Woloszynski, Podsiadlo & Stachowiak, *Efficient Solution to the
Cavitation Problem in Hydrodynamic Lubrication*, Tribol. Lett. **58**:18 (2015).
Selected with `cavitation_model = FBNS`; the `FULL_SOMMERFELD`, `GUMBEL` and
`JFO` paths are untouched and verified bit-identical (byte-compare of
`diagnostics.csv` and every VTK/PVD product of `cases/00c_jfo_verification/
oneD_n120_b1e9` on 4 ranks, before and after the whole change).

**Why.** `cases/00c_jfo_verification` established that the Elrod-Adams rupture
front carries a first-order error with a ~7-cell prefactor: at n_theta = 120..960
the front sits 7-13 *cell widths* upstream of the exact 1D solution, and the
error does not shrink relative to the mesh. The 2026-08-04 `mass_residual` work
already ruled out a non-conservative flux or boundary treatment - the converged
state is per-cell conservative at round-off. What is left is the free-boundary
treatment itself: Elrod reads the switch through the degenerate `p - p_cav`
signal, which is `O(beta^-1)` small near the front because `p = p_cav + beta
ln(theta)`. FBNS removes the mechanism rather than tuning it.

**Formulation** (PHYSICS.md §7.3). Incompressible liquid; unknowns `p` and
`theta` per cell. The JFO conditions `p - p_cav >= 0`, `1 - theta >= 0`,
`(p - p_cav)(1 - theta) = 0` are imposed directly through the smooth,
unconstrained Fischer-Burmeister equation `FB(a,b) = a + b - sqrt(a^2 + b^2)`
with `a = (p - p_cav)/p_ref` and `b = 1 - theta`, paired with the discrete
Reynolds mass balance (pressure-driven Poiseuille, first-order upwind Couette
transport of `theta`, backward-Euler `d(rho h theta)/dt`). Newton's step uses the
paper's column swap - the unknown with the larger FB derivative goes into the
first block, bounding `cond(A_F) <= 4 + 2 sqrt(2)` - and the Schur reduction, so
one 5-point sparse solve per iteration through the existing `LinearSystem`/PETSc
KSP. **No artificial bulk modulus and no switch function** anywhere: a cavitated
cell has `p = p_cav`, so its pressure gradient - and hence its Poiseuille flux -
vanishes on its own, and the reformation jump is carried by the downstream cell's
pressure gradient.

**Scaling is not optional.** FB weights its two arguments equally; raw pascals
against a film content in [0,1] make the pressure branch dominate by six orders
of magnitude, drive `J_F,p -> 0` everywhere, and degenerate the swap. Pressure is
nondimensionalized by `p_ref` (`fbns_pressure_scale`, 0 = the Couette pressure
scale `6 mu |omega| R^2 / c^2` raised to cover any supplied overpressure); the
mass balance is scaled by the largest assembled row sum, which only sets what
`fbns_tol` means.

**Verified** against `cases/00c_jfo_verification/exact_1d_jfo.py` (exact rupture
302.766608 deg, reformation 79.006725 deg, p_max 7.597141 MPa, W' 495463.6 N/m)
in new sibling cases `oneD_n{120,240,480,960}_fbns`, 4 ranks, front reported at
the face between the last full-film and the first cavitated cell:

| n_theta | cell | FBNS rupture err | cells | JFO rupture err | cells | FBNS reform err | FBNS p_max | FBNS W' |
|---|---|---|---|---|---|---|---|---|
| 120 | 3.000 deg | -2.767 deg | -0.92 | -23.77 deg | -7.9 | +1.99 deg (+0.66) | +0.54% | +0.70% |
| 240 | 1.500 deg | -1.267 deg | -0.84 | -17.77 deg | -11.8 | +0.49 deg (+0.33) | +0.28% | +0.27% |
| 480 | 0.750 deg | -0.517 deg | -0.69 | -10.27 deg | -13.7 | +0.49 deg (+0.66) | +0.14% | +0.07% |
| 960 | 0.375 deg | -0.142 deg | -0.38 | -5.02 deg | -13.4 | +0.12 deg (+0.32) | +0.07% | -0.04% |

JFO `p_max` is +50.4% / +39.9% / +24.0% / +12.2% and `W'` +42.0% / +41.0% /
+31.1% / +17.8% on the same grids. FBNS converges in **4 Newton iterations at
every resolution**; the n = 480 residual history is
`6.40e-4 -> 6.02e-4 -> 3.93e-5 -> 1.23e-5 -> 3.49e-12`, the last step dropping
seven orders because once the active set stops changing the remaining system is
exactly linear and the semismooth Newton step terminates finitely.

**mass_residual.** `Diagnostics::mass_residual` gained an FBNS branch: the FBNS
face fluxes differ from `solve_elrod`'s (pressure-driven Poiseuille without beta
or `g(theta)`, always-upwind Couette, plain half-cell axial Dirichlet link), so
the field mirrors `FBNS::cell_outflow` - the same function the Newton residual
uses, so there is one operator, not two. No clamp term: the Newton truncations
act on iterates, not on the converged state. On the converged n = 480 run the
worst cell is 1.5e-11 kg/s = 3.5e-10 of the peak Couette face flux (the Krylov
tolerance), and 3e-14 of it on the rupture and reformation cells themselves;
groove cells are exactly zero, their imbalance being reported as
`inlet_mass_source`. `liquid_boundary_inflow` gained the matching FBNS branch.

**Boundaries.** Inlet cells replace *both* rows - mass balance becomes
`p* = p*_supply`, complementarity becomes `alpha = 0` with
`(J_F,p, J_F,alpha) = (0,1)`. Forcing that pair also forces the swap, which is
what keeps the Schur complement non-singular on those rows; no penalty
coefficient is used. Axial `DIRICHLET` is a half-cell Poiseuille link, so a
cavitated boundary cell re-floods at the physical rate that the Elrod path needs
`consistent_boundary_flux` to approximate. `NEUMANN` carries no term and needs no
`theta` boundary condition, because Couette transport is purely circumferential.

**Scope, enforced not assumed.** Phase 1 is `ISOTHERMAL` +
`fluid_property_model = CONSTANT`, steady or transient, static or moving bearing,
`DIRICHLET`/`NEUMANN` axial ends, MPI-parallel. The energy equation, the
dissolved-gas / two-phase layer, `gas_pressure_coupling`, the cross-film Dowson
coefficients and `INLET_OUTLET` each couple back into the Reynolds operator and
are **not** in the FBNS Jacobian, so `validate()` rejects them with a message
pointing at `cavitation_model = JFO`. `theta_convection_scheme` is warned as
ignored. `cavitated_film_model = STRIATED` now accepts `JFO` or `FBNS` (both
carry a mass-conserving `theta`).

**New config keys**: `fbns_max_iters` (100), `fbns_tol` (1e-9),
`fbns_relaxation` (1.0), `fbns_pressure_scale` (0 = auto). `bulk_modulus` is
unused on this path.

**Tests**: `tests/test_fbns.cpp` (2 and 6 ranks) covers pointwise
complementarity, per-cell mass conservation, the discrete cavitated-zone flux
invariant `h_face * theta = const` across rank boundaries (the rank-independence
check), free-boundary accuracy against the exact benchmark *including* the
requirement that Elrod be at least 4x worse on the same grid (so the benchmark
cannot silently stop discriminating), transient backward-Euler storage, Dirichlet
axial ends with the telescoping identity, and the phase-1 validation gates.
27/27 tests pass. `tests/test_coupling.cpp` gained the STRIATED + FBNS case.

### 2026-08-04 - Per-cell mass-residual output field (case 00c test #3)

New optional output field `mass_residual`: the per-cell liquid-continuity
residual of the converged JFO state, in kg/s (cell-integrated, not
area-normalized). Written only when listed in `output_fields`.

**Why.** `diagnostics.cpp` reports one number for the whole domain, and it closes
at solver tolerance. That is compatible with a genuinely conservative scheme and
equally compatible with large local errors that cancel on summation — which is
exactly the failure mode that could produce the grid-divergent rupture front case
00a found. The field makes the distinction observable in 2D.

**Definition** (PHYSICS.md §15.1). Net liquid mass flux out of the four faces
(Poiseuille + Couette) plus the transient storage and squeeze terms, minus the
mass the theta clamp/snap injected this step. Every term is taken from the
expressions `solve_elrod` assembles rather than re-derived: harmonically averaged
`Gamma_base * g(theta)` on interior faces, `cfg.elrod_boundary_outflow` on active
axial Dirichlet / `INLET_OUTLET` faces, and the configured `theta_convection_scheme`
face value. Inlet penalty cells carry only the clamp term: their flux imbalance is
the modelled supply source that the global balance already reports as
`inlet_mass_source`. Interior faces telescope, so the field sums exactly to the
global closure: `sum(r) == liquid_residual * liquid_mass / step`.

**Factoring.** `Diagnostics::elrod_inlet_cell_outflow` was already the per-cell
operator; it is now `elrod_cell_outflow` and is shared by the global balance and
the residual field. Two consistency fixes fell out of that sharing, both of which
made the assembled operator and the diagnostic agree where they previously did
not:

- the per-cell outflow only opened an axial boundary face under `DIRICHLET`,
  while `liquid_boundary_inflow` also opened it under `INLET_OUTLET` on inflow.
  An inlet cell sitting on an `INLET_OUTLET` axial row therefore had its boundary
  flux counted in the domain flux but not in its own inlet source. Both now share
  one `axial_link_active` predicate, matching `solve_elrod`.
- `couette_face_theta` mirrored only `UPWIND` and `TYPE_DIFFERENCING`; under
  `LINEAR` or `VANLEER` the diagnostic used upwind where the assembly used the
  deferred central / limited value. Both branches are now mirrored. (No effect on
  any shipped case: all of them run `UPWIND` or `TYPE_DIFFERENCING`.)

**Cost when unselected.** `main.cpp` registers the field only if
`cfg.output_field_enabled("mass_residual")`, and `Diagnostics::mass_residual`
returns immediately when the field is absent. Verified end to end: running
`cases/00c_jfo_verification/oneD_n120_b1e9` with and without the field gives
byte-identical `diagnostics.csv` and byte-identical VTK arrays for every other
field.

**Test** (`tests/test_mass_residual.cpp`, 2 ranks, joins ctest):

1. Steady grooved bearing (00c geometry): 464 clean full-film cells at
   7.1e-15 of the row conductance; field sums to the reported `r_l` to 2.7e-17.
2. Open Dirichlet axial ends, transient: boundary rows close as tightly as the
   interior (2.0e-14); sum identity holds to 4.3e-20 against a nonzero boundary
   flux and nonzero storage.
3. Perturbing one converged cell by 2% raises the field to 2.1e3 times the
   Couette flux on that cell and its four neighbours, leaves every other cell at
   its tolerance, and moves the global `r_l` by nothing at all. This is the
   feature's justification stated as an assertion.
4. Not registering the field is a no-op that leaves the global balance untouched.

**First result on case 00c** (`oneD_n120_b1e9`, converged steady state). Groove
cells are zero by construction; the cell immediately upstream of the groove edge
carries `theta = 0.259` next to a penalty cell pinned at `theta = 1.00007` and
still balances to 6e-17 of the Couette flux. The rupture and reformation fronts
sit at 8e-11. The largest value in the domain, 1.1e-7 of the Couette flux
(3.2e-11 of the row conductance), is a smooth blob over the peak-pressure zone —
the Krylov residual profile at the default `linear_rtol = 1e-8`, not a flux
defect. **The JFO assembly is per-cell conservative at this resolution, so the
grid-divergent front is not caused by a non-conservative flux or boundary
treatment.** Test #3 of case 00c is closed; the front's grid sensitivity has to
be explained by the scheme's accuracy or the free-boundary treatment, not by
conservation.

Files: `src/diagnostics.{hpp,cpp}`, `src/main.cpp`, `src/config.hpp`,
`tests/test_mass_residual.cpp`, `tests/CMakeLists.txt`, `PHYSICS.md`,
`README.md`, `cases/00c_jfo_verification/README.md`. Full suite 25/25 green.

### 2026-08-02 - Phase 2.3/2.4: groove carry-over, cross-film T(y), conjugate
walls, generalized Reynolds

Three coupled additions, all off by default and all bit-identical to the previous
solver when off. New module `src/cross_film.{hpp,cpp}`, new test
`tests/test_cross_film.cpp` (24/24 tests pass), physics written up in PHYSICS.md
sections 13.1.1-13.1.3.

**Why they are one change.** The wall heat flux, the effective viscosity in the
Reynolds coefficient, and the shape of the temperature profile across the film
are three readings of the same unknown T(y). Fixing any one of them alone means
inventing a closure for the other two, which is what the constant-h_wall Newton
cooling was.

**2.3 - hot-oil carry-over at grooves** (`groove_mixing_model`).
`T_groove = lambda*T_supply + (1-lambda)*T_upstream`, with lambda either given
per inlet (`FIXED_LAMBDA`) or computed from the supply flow against the solved
upstream Couette carry-over (`FLOW_WEIGHTED`). The upstream sample is taken per
axial row at the cell before the groove's leading edge, in the sense the shaft
entrains; that cell can live on any rank, so its temperature and flux are gathered
with a sum reduction over the (undecomposed) axial index. `PINNED` is the default
and reproduces the old hard pin exactly.

**2.4a - cross-film profile.** T(eta) = T_bar + a psi1 + b psi2 + c psi4 on three
zero-mean shape functions, so the depth average is the depth-averaged solve's own
answer by construction. Closed by the two wall values plus one midplane
collocation of the cross-film energy balance. Recovers the classical
conduction/dissipation parabola T = T_w + S eta(1-eta) exactly.

**2.4b - conjugate walls.** The shaft is circumferentially isothermal at a
temperature solved from the global film-to-shaft balance; with
`shaft_ambient_heat_transfer = 0` that balance is "net zero", so heat entering the
shaft in the loaded arc leaves it again upstream. This is the transport path the
prescribed-T_wall closure has no representation of at all. The bush closes per
cell through 1-D radial conduction to an outer convection boundary. Both reduce
to closed form; only four scalars are reduced across ranks. When active, the two
Newton-cooling coefficients and the two prescribed wall temperatures are ignored,
and the validator says so rather than letting them look effective.

**2.4c - generalized Reynolds.** Dowson's F0/F1/F2 replace h^3/(12 mu) and h/2:
q = -F2 grad p + U F1/F0. Not an effective-viscosity correction of h^3/12mu - the
integrals are the coefficient. The dissipation source follows exactly from the
same integrals, Phi = U^2/F0 + F2 |grad p|^2 (the Couette-Poiseuille cross term
integrates to zero identically), and the energy equation's convective flux uses
the same q the pressure solve conserved.

**Two approximations, both forced by measurement rather than convenience, both
documented in the source with the evidence:**

1. *The collocation's convective term is saturated against its dissipative term.*
   One midplane collocation cannot represent the thin wall layers that appear at
   reduced Peclet 10-30. Unbounded it produced a 494 K core rise on the Ferron
   rig, against a physical scale of about 3 K. Saturated, the core rise is 1.1 K
   and the conduction-dominated limit is still exact.
2. *The collocation uses the Couette part of the flow only.* A dp/ds dependence
   makes the profile shape jump whenever an Elrod cell flips between full film
   and cavitated, closing a loop profile -> viscosity -> Couette flux -> rupture
   position -> profile. Measured on the Ferron rig: the coupling residual R_p sat
   at 9.4 bar, *identical at 40 and at 120 outer iterations*, and under-relaxing
   the profile amplitude did not move it (9.4 / 2.3 / 7.6 / 2.2 bar at
   w = 1.0 / 0.5 / 0.3 / 0.15). Removing the pressure terms from the collocation
   alone took the same case to R_p = 8e-9 and changed the load by 0.7%. What is
   given up is about 6% of the dissipation and 14% of the midplane velocity, both
   feeding a shape correction that the saturation already caps; everything else
   still uses the full pressure field.

**Measured on Ferron (2000 rpm, eps = 0.534, 720x40, both legs converged to
R_p < 1e-3): the load prediction improves, the temperature-peak location does
not move at all.**

    leg                                   W/N   T_max/C  T_max err  T-peak X/pi  err/deg
    depth-averaged, adiabatic (before)   3971    47.48    -1.62 K       0.949     +71.8
    cross-film + conjugate (after)       3710    46.74    -2.36 K       0.949     +71.8

Load falls 6.6%, in the predicted direction and for the predicted reason: mu is
convex in T and the core is the hottest part of the film, so the effective
viscosity implied by F2 is below mu(T_mean). T_max falls with it, which moves
T_max slightly *away* from Ferron's measurement and outside the 2 K band.

**The 71 deg peak-location error does not move, to three decimals in X/pi.** The
diagnosis this task was built on — that a conjugate shaft carrying heat upstream
is what relocates the peak — is not supported by this implementation. Sized from
the run itself: the shaft exchange is (k/h) dT ~ 4 kW/m^2 against an advective
enthalpy divergence rho c h u dT/ds ~ 60 kW/m^2, so the wall path is ~6% of
advection, and a *circumferentially isothermal* shaft spreads what it does return
uniformly over 360 deg. Full detail, and the two candidate explanations (the
measurement is a bush sub-surface temperature; the shaft model is too weak),
are in VALIDATION_PLAN.md under "Phase 2 results".

Also recorded, two robustness findings that are **not** caused by this change but
that bound what can be validated:

- With a strongly cooled bush (`bush_outer_heat_transfer = 500`) the Ferron case
  aborts in the Elrod solve with `DIVERGED_BREAKDOWN`. The Phase 0.1 KSP guard
  catches it loudly rather than returning a plausible-looking field, which is the
  guard doing its job, but the wall sensitivity could not be swept as far as
  intended.
- The Ferron THD case at eps = 0.67, 720x40 **diverges at 4 MPI ranks and
  converges at 1 and 2**. Verified against a binary built from the pre-change
  HEAD: identical trace, identical abort. It is a pre-existing rank-dependence,
  not something introduced here, and it contradicts the rank-independence note in
  `mpi_launch.py`, which was established on the *isothermal* Ferron case. All
  Phase 2 validation runs were therefore done at 2 ranks.

**Bit-identity when off, verified rather than assumed.** Outputs compared
byte-for-byte against a binary built from the pre-change HEAD, on four
configurations: JFO steady, GUMBEL steady, transient, and 4-rank MPI. All four
identical. The classical branches keep their original floating-point association
deliberately (see `reynolds_diffusion`'s `prefactor` argument), because
`rho*(h*h*h/(12*mu))` and `rho*h*h*h/(12*mu)` are not the same double.

### 2026-08-04 - Brito Fig. 6 independently re-digitised: verdict re-scored

Seongho re-read Fig. 6 with digiplotter, blind to the values in `brito.py`. The
six points agree with the committed polar re-digitisation to |d eps| <= 0.013
(mean 0.005) and |d phi| <= 1.0 deg -- inside the stated uncertainties at every
point -- and refute the Cartesian alternative that sat ~0.15 lower. Recorded as
`ECC_ATTITUDE_INDEPENDENT`; the spread of the two readings is now the measured
reading uncertainty.

Consequently the load comparison is scoreable again, and it goes from
NOT SCOREABLE to **FAIL / OPEN**: the solver reaches the applied 10 kN at
eps = 0.648 / 0.426 / 0.367 (1000 / 3000 / 4000 rpm, c at 20 C) against measured
0.770 / 0.577 / 0.545 -- a deficit of 0.12-0.18 growing with speed, equivalently
+65-106% load at the measured eccentricities. Same direction and comparable
size as the Ferron discrepancy that Ferron's own theory shares. Brito's
eccentricity is likewise probe displacement over an assumed clearance, so the
clearance-systematics suspect carries over, but the digitisation no longer can.

### 2026-08-02 - Phase 2.1/2.2: consistent theta-weighting (measured: no steady effect)

Defect D14 fixed. Under `cavitated_film_model = STRIATED` every thermal term that
scales with liquid content is now weighted by the liquid fraction through one
shared helper: heat generation, stored energy, conduction, and convected enthalpy
in both directions. Previously only the Couette shear was weighted while
capacity, conduction and the convective flux carried the full gap, so a ruptured
cell stored and conducted as though it were solid oil. `FULL_FILM` remains the
default and is unchanged, so the cavitated zone is still not thermally modelled
unless asked for; 23/23 tests pass either way.

**Measured effect on Ferron (2000 rpm, adiabatic, 720x40): none.** T_max moves
0.003 K, load 0.2 N, peak pressure and rupture unchanged. The reason is an exact
cancellation, and it corrects an estimate recorded earlier in this campaign:

    dT/ds ~ Phi / (rho_cp theta h u) = (mu U^2 theta / h) / (rho_cp theta h U/2)
          = 2 mu U / (rho_cp h^2)          -- theta cancels

Weighting the heat source and the mass flux by the same fraction leaves the
steady temperature rise unchanged. The earlier claim that carrying the full film
through the cavitated arc "overstates the mass flux by ~2.5x and damps the
temperature rise" was wrong: it considered the flux in isolation and missed that
the source scales identically. The 0.003 K residual is conduction, which does not
cancel but is negligible against advection at this Peclet number.

So the fix is correct physics and makes the model self-consistent, but D14 is
**not** a contributor to the temperature discrepancy and is now ruled out as a
candidate. It should still matter where the cancellation breaks: transients
(the storage term), and any regime where conduction competes with advection.

### 2026-08-02 - Brito load comparison is not scoreable as it stands

Applied the Ferron lesson before repeating the mistake: instead of prescribing a
digitised eccentricity and calling the load a prediction, solve for the
eccentricity that carries the load Brito APPLIED and therefore knew exactly, so
the uncertain quantity becomes a check rather than an input. `solve_eps_for_load`
already existed in the harness and had never been called.

Eccentricity carrying 10 kN at Pf = 70 kPa, Tf = 40 C:

    N rpm    c at 20 C    c hot (+8.8%)   digitised (polar)   d eps
     1000        0.648          0.678               0.761    -0.113
     3000        0.426          0.461               0.576    -0.150
     4000        0.367          0.401               0.558    -0.191

- **The digitisation ambiguity is larger than the discrepancy.** `brito.py`
  records that the earlier Cartesian reading of Fig. 6 -- a polar chart -- gave
  eccentricities about 0.15 lower, i.e. near 0.61/0.43/0.41. The solver answers
  0.648/0.426/0.367, tracking the Cartesian reading to within ~0.04 at every
  speed. That does not prove the polar re-digitisation wrong -- using the solver
  to choose the digitisation would be circular -- but it does mean the Brito load
  comparison cannot discriminate until Fig. 6 is re-read from the source. It is
  now reported as NOT SCOREABLE rather than as a 65-106% error.
- The hot-clearance bracket closes 0.030-0.035, about a quarter of the gap at
  1000 rpm. It rests on assumed metal temperatures Brito did not measure, so it
  is a bracket and never a correction.
- The deficit grows monotonically with speed, the qualitative signature of a
  thermal mechanism, but not as U^2, so it does not by itself confirm the
  cross-film hypothesis. Eccentricity compresses a load error; re-test the
  scaling on load at fixed eps after the Phase 2 energy fixes.
- The two 1000 rpm rows did not converge and are marked in the JSON.

### 2026-08-02 - Validation harnesses can run on multiple ranks

`PANCAKE_VALIDATION_RANKS` selects the rank count (default 1). Verified exact:
1 and 4 ranks agree bit-identically at `linear_rtol = 1e-12`, and to 4.9e-5 at
the 1e-8 default, so the residual is the linear solve rather than the
partitioning.

Enabling it exposed a silent truncation in `read_field`: under N ranks the solver
writes one theta slice per rank to `processorN/`, and the reader only ever opened
`processor0`, returning 1/N of the circumference. It failed in the worst way --
the rank-independence check on LOAD passed throughout, because load comes from
`Fp`, a globally reduced uniform field that is correct in any single rank's file.
Only a downstream IndexError exposed it. `read_field` now stitches the per-rank
pieces (rank N writes `solution_<step>_N.vts`, so the rank-0 basename does not
exist elsewhere) and raises on a stitched extent that does not match `n_theta`.

### 2026-08-02 - Case data archived behind each verdict

`case_archive/` holds the verbatim `config.txt` for all 197 cases plus the 2-D
fields for the 38 that back a figure or a verdict, 11 MB total. Previously both
existed only in the gitignored 1.4 GB `_work/` tree, so a case was reproducible
only through a harness that may since have changed, and no contour could be drawn
without re-running the solver. `archive_cases.py --all` regenerates the full 2-D
set; it is derived data, not a source of truth.

### 2026-08-02 - Cross-solver residual explained; consolidated validation report

- **The few-percent pancake-minus-reference difference is mesh convergence, not a
  defect.** At fixed `beta = 1e8` pancake moves monotonically onto the
  independent solver under refinement: 15010 N (720x40), 14717 N (1440x60),
  14540 N (2880x60) against a reference already grid-converged at 14129 N. The
  slower mesh convergence is expected: pancake's Couette term is a first-order
  upwinded flux on the film content, while the reference solves for pressure
  directly with the wedge as a source term and has no convective discretisation
  at all. **No defect was found.**
- **Raising `bulk_modulus` toward incompressibility makes agreement worse, not
  better** -- 2.9%, 3.5%, 9.4% at beta = 1e8, 3e8, 1e9 on a fixed 2880x60 mesh.
  The expectation was the opposite and is recorded as refuted. Raising beta
  drives `p_max/beta` toward the 4e-3 free-boundary resolution floor, trading a
  small compressibility difference for a larger under-resolution error; the two
  oppose each other and beta = 1e8 is near optimal for this case. Do not raise
  beta to flatter a comparison.
- **Practical consequence:** the production 1440x60 Ferron mesh carries ~4% of
  discretisation error in load -- inside the Raimondi & Boyd tolerance and an
  order of magnitude below the discrepancy under investigation, but it must be
  stated whenever load is quoted to better than ~5%.
- **New `validation/VALIDATION_REPORT.md`**: one document holding every
  verification and validation verdict, with the reproducing command against each
  number, verification and validation kept separate, and the reference-data
  provenance caveats attached to the numbers they qualify.
- **New `make_report_figures.py`** producing `validation_cross_solver.png` and
  `validation_load_in_context.png` in the repository plot style. The literature
  comparison points are plotted as isolated markers rather than joined into
  curves, because they are secondhand single digitisations.

### 2026-08-02 - Phase 1: the load over-prediction is not a solver defect

Two hypotheses from `VALIDATION_PLAN.md` were tested and both are wrong. Recorded
because the wrong version is easy to re-derive.

- **Ferron's geometry in `ferron.py` was already correct**, confirmed against
  three independent open sources: L/D = 0.80, c = 145 um radial. The competing
  `L/D = 0.55, C/R = 1.5e-3` belongs to a *different experiment*, Tonnesen &
  Hansen (1981) -- also a 100 mm bearing loaded to ~9 kN, so the two are easy to
  conflate. 1.5e-3 is also what you get by misreading 145 um as diametral, so the
  wrong number has two independent routes into a config file.
- **The eccentricity-vs-load data is 2000 rpm.** Alakhramsing's Fig. 10a carries
  no speed label and the 1000 rpm hypothesis fitted the loads far better, but
  Lubricants 11(7) 308 (2023) Fig. 8 reproduces the identical dataset captioned
  "At 2000 rpm", agreeing within 0.03 in eps at every load. The hypothesis was
  wrong.
- **Ferron's own theory over-predicts his own measurements by ~75-80%**, worse
  than pancake's 52-68%: his theory reaches eps = 0.67 at ~7.2 kN where he
  measured it at 4 kN, and Alakhramsing's FE model is displaced the same way.
- **New `christopherson.py`**: an independent Python finite-volume Reynolds
  solver with Christopherson clamping, sharing no code with pancake and having no
  `p_max/beta` failure mode. Self-audited against the Raimondi & Boyd charts to
  1.3% in S and 0.80 deg in attitude angle -- tighter than pancake's own
  5.8%/4.08%. `audit()` regenerates the check.
- **New `run_cross_solver.py`** runs both solvers isothermally at Ferron's exact
  geometry, the first independent check ever made at L/D = 0.80, c/R = 2.9e-3:
  3248/7933/15010 N (pancake) against 3206/7591/14129 N (independent) at
  eps = 0.47/0.67/0.78, agreeing to 1.3-6.2%, peak pressure to 3-4%. Two
  unrelated discretisations agree while both exceed the measured load by 60-90%.
  **The discrepancy is in the reference data or the shared physics, not the
  discretisation.** pancake reads a few percent high and growing with eps, which
  is worth a separate look but is not the 50% effect.
- **The measured temperature is a bush temperature 0.5 mm below the surface**,
  not a film temperature, so it is not the quantity pancake reports. The two open
  reproductions also disagree by ~30 deg on where the peak sits inside a 250-285
  deg plateau, bounding how precisely the 71 deg peak-location error can be scored
  against this data. Recorded in `ferron.py`.

Consequence: stop chasing W(eps) against Ferron's measurements, which nobody
matches including Ferron. Validate on peak pressure, its position, and the
temperature profile; quote W(eps) as an open discrepancy shared with the
literature with the clearance-derived eccentricity as the leading suspect.

### 2026-08-02 - Validation plan and Phase 0 guardrails

Reason: the temperature and load mismatches against Ferron (1983) and Brito
(2006) were being investigated on a solver that could not report a failed linear
solve, and on configurations that could silently disable the very coupling under
test. Three audits (defect register, energy-equation implementation, an
independent re-solve of both load cases) were combined into `VALIDATION_PLAN.md`,
and its Phase 0 landed here. No physics changed in this entry.

Diagnosis recorded in `VALIDATION_PLAN.md`:

- **The load over-prediction is not a discretisation bug.** Load integration,
  film thickness, clearance handling, groove pressure and viscosity wiring all
  audit clean, and an independent Christopherson solver reproduces pancake's
  loads at the same inputs. It is inherited from reference data and from shared
  modelling assumptions.
- **The 71 deg temperature-peak error is structural**: with wall temperature a
  prescribed constant there is no path for heat to travel circumferentially
  against the flow, so no value of `h_wall` can move the peak.
- **The uniform-across-the-film temperature assumption is invalid on these
  rigs**: Br ~ 2.6, reduced Peclet ~ 10-30, unequal walls. Dissipation-driven
  core superheat is ~3 K (Ferron 2000 rpm) to ~6 K (Brito 4000 rpm) and scales
  as U^2, worth ~10-20% of load at high speed. Wall heat flux and effective
  viscosity both follow from the same missing `T(z)` profile, so they cannot be
  fixed independently. This gives a discriminator: after the reference data is
  settled, a residual error growing as U^2 is cross-film physics; a
  speed-independent one is geometry or property data.

Phase 0 changes:

- **Defect D6 fixed.** `LinearSystem::solve` now takes an equation label, checks
  `KSPGetConvergedReason`, scans the solution for non-finite values under an
  `MPI_Allreduce` so all ranks agree, and aborts with a named diagnosis. All five
  call sites labelled. Previously a diverged solve or a NaN was copied into the
  field and the run reported success.
- **New test `test_ksp_guard`** (23 tests total, all passing). It asserts on the
  guard's own message via `PASS_REGULAR_EXPRESSION` rather than on a nonzero exit,
  so a PETSc-internal failure that never reached the guard cannot pass it. The
  test aborts by design.
- **Silent THD decoupling now warns.** `temperature_model = ENERGY_EQUATION`
  combined with `fluid_property_model = CONSTANT` or
  `liquid_viscosity_model = CONSTANT` freezes viscosity, so the energy equation
  runs and writes a plausible temperature field with no effect on pressure or
  load. The shipped `config.txt` was in exactly this state; it now carries a
  comment saying so and the solver warns at startup. The example's physics was
  left unchanged deliberately, since enabling coupling requires the
  `EMPIRICAL` + `solution_viscosity_gas_coeff` workaround.
- **Adiabatic walls now warn.** Both wall heat-transfer coefficients at zero is a
  perfectly insulating film; temperature results must be quoted with the wall
  condition attached.
- **Fed inlets with non-positive `t_supply` are now rejected.** `t_supply`
  defaults to 0 K and `feeds_fresh_oil` is set by the mere presence of the
  optional token, so a mistyped value silently pinned the groove to absolute zero.
- **Deferred-correction staleness (B2) fixed in the energy and gas paths.**
  `FVM::divergence` reads `phi.old`, which on a live field is the previous
  timestep: stale across outer iterations and, in `STEADY_STATE` (one timestep,
  many outer iterations), never refreshed at all, silently degrading LINEAR and
  VANLEER to first-order upwind. `energy.cpp` and `gas_transport.cpp` now pass a
  lagged copy carrying the current iterate, the idiom `solve_elrod` already used
  for theta -- the Elrod path was correct and was left alone.
- **`run_brito.py` now records `elrod_iters`, `coupling_iters` and
  `coupling_capped`**, so a run that exhausted the coupling budget is no longer
  indistinguishable in the JSON from one that converged. `coupling_max_iters` is
  now sourced from a single constant instead of being hardcoded in the template.

### 2026-07-31 - Cavitation and energy-equation validation

Reason: the committed analytical benchmarks validated the Reynolds
discretisation only -- they ran `FULL_SOMMERFELD` and `GUMBEL`, so the
mass-conserving Elrod-Adams free boundary, the `LOCAL_PSAT` threshold and every
predicted rupture angle and void fraction were unvalidated. Three targets were
attempted; two produced usable verdicts.

- **New validation modules** in `validation/journal_bearing/`: `raimondi_boyd.py`
  (chart data plus an exact infinite-length Swift-Stieber solution used to audit
  it), `ferron.py` (Ferron's rig and digitised measurements with per-number
  provenance), `run_raimondi_boyd.py`, `run_ferron.py`, `make_tables.py` and
  `plot_validation.py`.
- **A free-boundary resolution requirement was found, and it is not
  self-announcing.** Elrod carries pressure inside the film content
  (`p = p_cav + beta ln theta`), so the accuracy of the rupture boundary is
  controlled by `p_max/beta` and degrades as the film becomes more
  incompressible. At `beta = 1e9` on a 720-cell mesh the rupture angle is
  **20.8 deg** from the exact Swift-Stieber value while the solver reports
  `mc = 0`, `ns = 0`, `R_p = 0` and exits zero. At `beta = 1e8` with 2880 cells
  it is within 0.23 deg. Recorded as `CONVERGENCE_MEASURED` and regenerable via
  `convergence_study()`.
- **Target A (Raimondi & Boyd) passes with a fed film.** With a full-film supply
  groove at the maximum gap -- the inlet condition the charts assume -- the
  Sommerfeld number is within 5.8% and the attitude angle within 4.08 deg over
  eps = 0.2-0.9 at L/D = 1/2 and 1, against 10% and 5 deg criteria. Every case
  converged on the Elrod tolerance.
- **Target A also shows the ungrooved bearing has no load capacity.** Without a
  supply, load falls three to six orders of magnitude below the chart at every
  eccentricity (0.007 N against kN at L/D = 1, eps = 0.2). Below eps ~ 0.6 the
  film is driven onto the `theta_min` clamp; above it the solution settles on the
  analytic starved state `theta h = const`. This follows from the Elrod switch
  function suppressing Poiseuille flow in cavitated cells, so nothing re-floods
  the film.
- **Target C (Cole & Hughes 1956) could not be sourced** and is recorded as such.
  Its reformation criterion is additionally untestable: ungrooved the film never
  reforms, grooved the reformation point is pinned at the groove by construction.

Reference data is cited per number and audited, not taken from memory: the
transcribed Raimondi & Boyd table reproduces the exact Swift-Stieber solution to
0.6% in S and 0.4 deg in attitude angle. Ferron's geometry is corrected against
the commonly circulated secondary values -- the rig is L/D = 0.80 with
c = 145 um, not L/D = 0.55 with C/R = 1.5e-3.

### 2026-07-31 - VOID_COUPLED: bounded plateau closure and coupling relaxation

Reason: `gas_pressure_coupling = VOID_COUPLED` was unusable on a sub-atmospheric
optical-rig case (50 mm bore, c = 250 um, L/D = 0.75, white mineral oil + air/CO2
at 1 atm). It produced a 10 831 bar peak pressure with air and never converged
with CO2, so the gas-species dependence of gaseous cavitation could not be
simulated at all. Two independent defects, both fixed here.

- **Unbounded plateau (divergence).** `void_params` evaluated the cavitated-zone
  plateau explicitly as `p_void = m_g R_s T / ((1 - theta) h)` with the void
  floored at `1e-6`. The flash in `equilibrate_gas` assigns free gas wherever
  `p < p_sat` regardless of film content, so a *full-film* cell could carry gas
  and divide it by a ~zero void, giving O(1e4) bar plateaus. Replaced with the
  closed-form fixed point of the same balance: requiring the gas released between
  `p_sat` and `p` to fill the void `(1 - theta)` at ideal-gas density gives
  `p = L theta p_sat / [(1 - theta) + L theta]` with `L = H rho_l R_s T` the
  Ostwald coefficient. This tends to `p_sat` as `theta -> 1` instead of
  diverging, never exceeds `p_sat`, and is monotone in `theta`. New helpers:
  `Reynolds::gas_inventory_plateau`, `SimulationConfig::ostwald_coefficient`,
  `SimulationConfig::gas_specific_constant_for_species`.
- **Limit cycle (non-convergence).** The two-way `alpha_gas <-> theta_full <->
  theta` feedback oscillated: the rupture front moved every iteration with zero
  flag flips, so the Elrod loop stalled at a finite mass residual (`mc ~ 0.13`)
  and ran to `max_outer_iters` forever. Added under-relaxation, weight
  `void_coupling_relaxation` (new key, default 0.3, active only under
  VOID_COUPLED): the lagged plateau is held in a new `p_void` field advanced once
  per Elrod outer iteration, and `alpha_gas` / `free_gas_mass` are relaxed toward
  the flash result. Relaxation changes only the iteration path, not the fixed
  point, so converged solutions are unchanged.
- A rejected intermediate fix is recorded because it looks correct and is not:
  capping free gas by the available void `(1 - theta)` also removes the
  divergence, but it suppresses outgassing in a supersaturated *full* film, which
  is precisely the physics WP-1 gate 4 (`test_gas_coupling`) asserts. Gas creates
  its own void by displacing liquid, so the cap must not be applied.
- All 22 tests pass, including the WP-1 acceptance gates (vaporous limit and
  pressure cap) unchanged.

**VOID_COUPLED is still not usable quantitatively.** The two defects above are
real and fixed, but they were not the whole story, and this must not be read as
the mode now working:

- Fixed: pressures are bounded (no more O(1e4) bar), and the *inner* Elrod loop
  converges cleanly (`mc -> 0`, zero flag flips) for the first solves.
- Not fixed: the *outer* coupling loop (`equilibrate_gas` <-> `solve_elrod`) does
  not converge on the rig case. It exhausts `coupling_max_iters` with
  `R_p = 0.57` (CO2, 180x30, w = 0.10) to `R_p = 7.2` (air, 360x60, w = 0.3)
  against a 1e-9 tolerance.
- The decisive evidence is relaxation dependence: at a fixed 180x30 grid, w =
  0.10 / 0.30 / 0.60 give p_max = 12.1 / 6.29 / 6.42 bar and cavity arc 64 / 100
  / 108 deg. Relaxation cannot move a fixed point, so these runs are not at one.
  A grid sequence 180x30 -> 360x60 likewise fails to converge (p_max 12.1, 8.9,
  8.6, 17.7 bar).
- Note the solver exits 0 and logs "Simulation finished successfully" when the
  iteration budget is exhausted, so process exit status is not a convergence
  test; `log_outer_iters = true` and the reported `R_p` are.
- Remaining work is a numerical-methods change, not a bug fix: the segregated
  Picard coupling between the flash and the Elrod solve needs replacing with an
  implicit/Newton treatment (or the flash needs folding into the Elrod switch
  itself). Until then use `gas_pressure_coupling = NONE` for quantitative work.

The closed-form plateau is independently useful even so, since it is exact:
at void fraction alpha the sustainable gaseous plateau is
`p = L p_sat / (alpha/(1-alpha) + L)`, which for p_sat = 1 bar and alpha = 0.4
gives 0.13 bar for air (L = 0.098, effectively a vaporous cavity) against
0.62 bar for CO2 (L = 1.08, genuinely gaseous).


### 2026-07-13 - Repository consolidation, agent Git rule, and Phase F plan

Reason: the latest GUI/config work and the maintainer handoff were left as a
large uncommitted change set on `main`, while stale feature branches no longer
represented the product state. Future agents also needed one explicit workflow
and an ordered response to the verified defects in `HANDOFF.md`.

- Consolidated the current GUI/config/property-table work and `HANDOFF.md` as
  the authoritative repository state; obsolete local feature branches are
  removed as part of the repository cleanup after validation and commit.
- Added the same compact Git rule to `AGENTS.md`, `CLAUDE.md`, and
  `CONTRIBUTING.md`: one task, one `codex/<scope>` branch, one verified merge,
  then delete the branch. It also prohibits mixing local agent state or
  generated artifacts with product work and requires explicit authorization
  for remote operations.
- Added `PLAN.md` Phase F, which converts `HANDOFF.md` D1-D16/B1-B10 into a
  dependency-ordered remediation program. Silent solver failure, EOS
  consistency, full-film classification, TVD directionality, MPI ghost parity,
  and config round-tripping precede refactoring or new physics.
- Reconciled `README.md` with the current GUI: expanded units, property-table
  file selection, cavitation controls, and contour color ranges.
- No numerical or physical model was changed by the cleanup itself.
- Validation evidence: `cmake --preset windows-native-mingw` and the release
  build completed successfully; CTest passed 22/22 tests in 24.28 s, including
  `test_elrod_a2_mpi6`, `test_schemes_mpi6`, and `test_gas_coupling`.

### 2026-07-07 - Maintainer handoff document (HANDOFF.md)

Reason: capture the current solver architecture, its invariants, and a freshly
verified defect list in one place for future contributors, before the original
context is lost.

- **Added `HANDOFF.md`** (repo root): one-page architecture/data-flow
  orientation, the nine invariants that keep the solver correct (ghost
  discipline, `old()` semantics, EOS triangle, single boundary-flux source,
  clamp accounting, bit-exactness gates, triplicated config surface,
  cross-platform contract, docs-move-with-code), a 16-item verified defect
  list (D1-D16) with file:line citations, the improvement backlog (B1-B10),
  doc-vs-code mismatches, and the verification gates (ctest, diagnostics
  residuals, grid-convergence harness, MPI, model-off equivalence).
- Headline new findings from this audit pass (all verified against current
  code, none fixed yet): the `FVM::divergence` van Leer TVD correction uses
  the wrong far-upwind cell for negative face fluxes and has no z-face branch
  (D1, latent for `omega > 0`); post-processing/energy/gas/boundary code still
  tests `theta >= 1.0` while the WP-1 solver uses `theta_full = 1 - alpha_g`
  (D2); initial `theta` at `main.cpp:179` uses raw `p_cav` instead of
  `effective_p_cav()`, so LOCAL_PSAT runs start ~`p_sat - p_cav` inconsistent
  (D3); thermal-property ghosts go stale between the gas and energy solves
  (D5, rank-dependent); `KSPSolve` convergence is never checked (D6); the
  Picard relaxation acts on `p`/`T` while the Elrod state is `theta` (D8); the
  force gauge datum drifted from PLAN.md WP-12's stated convention (D9);
  STRIATED weights shear only, not thermal capacity/advection (D14).
- No solver code changed in this entry.

### 2026-06-15 - Cavitation threshold/coupling in GUI, property-table loading in GUI, propane config fix

Reason: the propane case was cavitating at the oil value `p_cav = 30 kPa` instead of
the dissolved-propane bubble point `p_sat ~ 1 MPa` (PHYSICS.md S16.1), and the two
flags that fix it were not reachable from the GUI; separately, choosing a TABLE
property model errored with no way to supply the data file.

- **Cavitation physics exposed in the GUI.** The Case Setup > Cavitation panel now
  has **Cavitation threshold** (`SCALAR_PCAV` = oil `p_cav`, `LOCAL_PSAT` =
  `max(p_cav, p_sat)` = propane bubble point) and **Released-gas pressure coupling**
  (`NONE` / `VOID_COUPLED`), with hover text explaining each.
- **`gas_pressure_coupling` now serialized.** `to_config_text()` was parsing but not
  writing it, so a GUI save silently dropped it; added the line (config.hpp). The
  config round-trips and all 22 tests still pass.
- **Property tables usable from the GUI.** The GUI now links the solver's
  `property_table_io.cpp` and ingests the configured `*_table_file` CSVs into the
  in-memory 2-D tables (`GuiApp::reload_property_tables`, run on load / new / preset /
  table-file edit). TABLE models therefore validate exactly as the solver does, and
  the FluidProperties panel shows a **table-file path + Browse** for each TABLE model
  with a live "loaded N x M grid" / "not loaded" check. Verified: the supplier file
  `data/R290_PZ68S/r290_pz68s.csv` loads as 79 x 45 for both solubility and viscosity,
  0 errors, p_sat = 1.0 MPa.
- **`*_table_file` paths now serialized.** `to_config_text()` wrote the inline tables
  but not the `solubility_table_file` / `viscosity_table_file` / `density_table_file`
  paths, so a GUI Browse-and-save silently dropped the path and the solver got TABLE
  with no file (the recurring "viscosity_table needs ... entries" error, even after the
  GUI showed the table loaded). Added the three lines (config.hpp); round-trip verified,
  22/22 tests pass.
- **`data/` now ships next to the binaries.** The case configs reference the table CSVs
  by relative path (`data/R290_PZ68S/...`), which the solver resolves against the config's
  directory - so a TABLE run from a build dir failed validation ("viscosity_table needs ...
  entries") because `data/` was never copied there. CMake now copies `data/` next to both
  `pancake` and `pancake_gui` (new `pancake_copy_data_dir`), and `config_r290_pz68_table.txt`
  joins the copied case configs. Verified: the table config runs from inside the build dir.
  For configs saved elsewhere, Browse to an absolute CSV path (the GUI hint says so).
- **Shipped propane config corrected.** `config_r290_pz68.txt` now sets
  `cavitation_threshold = LOCAL_PSAT` (+ an explicit `gas_pressure_coupling = NONE`
  with both options documented). Verified by re-run: the cavitated-cell plateau moves
  30 kPa -> 1.0 MPa and the cavitated extent grows 68 -> 506 cells at the same step.
- PHYSICS.md S16.1 updated (the shipped config no longer carries the defect; the
  default `SCALAR_PCAV` still does, by design, and is now a labeled GUI choice).

### 2026-06-15 - GUI: unit options, motion simplification, contour colormap, species coefficient fill

Reason: usability pass on the Case Setup / Results panels. All edits are GUI-only
(`src/gui/`); `SimulationConfig` and the solver are untouched, so existing
`config_*.txt` files keep loading and the form still stores solver-native SI.

- **More display units.** `UnitOption` gained an affine `offset` (so a unit can be
  `scale + offset`, not just a multiplier) with shared `unit_to_display` /
  `unit_to_stored` helpers. New/extended families: pressure adds **bar**; bulk
  modulus gets its own Pa/kPa/MPa/**GPa** family; tilt slope adds **mrad**
  (small-angle: 1 m/m ~ 1 rad); temperature adds **degC** (offset -273.15);
  dissolved-gas mass fractions add **%**; bearing support gets **N/m, kN/m, MN/m**
  (stiffness) and **N.s/m, kN.s/m** (damping). Inlet supply temperature also uses
  the K/degC family. Conversions stay display-only; the config file is unchanged.
- **Cavitation model label.** "Full Sommerfeld (reference)" -> "Full Sommerfeld"
  (the "(reference)" tag was dropped; the hover description is unchanged).
- **Bearing motion simplified.** The form no longer exposes initial position
  (it comes from the geometry: eccentricity + attitude angle) or initial velocity
  (zero); per-axis stiffness/damping are collapsed to a single **isotropic**
  stiffness and damping. `GuiApp::normalize_motion_config()` enforces this contract
  (attitude-derived position, zero velocity, x-value mirrored to y/z) on new case,
  load, and save, so the written config stays consistent.
- **Field-map color scale.** New colormap selector — Viridis (sequential) or a
  **blue->white->red diverging** map (`diverging_colormap()`, registered once) — and
  an optional **manual min/max range** that overrides the auto data range.
- **Checkmark grouping.** In the heatmap/monitor toolbars and the Advanced toggle,
  each checkbox now sits tight to its own label with a wider gap before the next
  control, so a checkbox reads as belonging to its label.
- **Species coefficient fill.** Selecting a dissolved-gas species fills its
  coefficient block (`apply_dissolved_gas_species_preset`): **Propane** loads the
  calibrated R-290-in-PZ68S numbers from `config_r290_pz68.txt` (Henry slope
  2.175e-7 1/Pa, van't Hoff 1800 K, mu_gas 8.7e-6, etc., calibrated at the
  10 bar / 40 C reference state); **Air** loads a generic air-in-mineral-oil set.
  Only the coefficients are written; the property/solubility/density/viscosity
  model selectors stay as the user set them.

### 2026-06-15 - Config model-surface cleanup (rename / reorder / prune + LINEAR scheme)

Reason: the model selectors had accumulated redundant and mislabeled options. This
pass renames them for clarity, orders every enum simplest -> most accurate, prunes
options with no distinct justification, and documents each model in PHYSICS.md §17.
Old config tokens and keys remain accepted as aliases, so existing `config_*.txt`
files keep working; full suite stays green (22/22).

- **Renames** (enum types, values, config keys; old names still parse):
  - `cavitation_model`: `ELROD_ADAMS` -> **`JFO`** (the model is the Elrod-Adams
    algorithm for JFO cavitation; the internal `solve_elrod`/`elrod_*` code keeps its
    name). Order: `FULL_SOMMERFELD` -> `GUMBEL` -> `JFO`.
  - `fluid_property_model`: `OIL_DISSOLVED_GAS` -> **`SINGLE_PHASE`**,
    `GAS_CAVITATION_MIXTURE` -> **`TWO_PHASE`**.
  - `oil_gas_solution_model` -> key **`solubility_model`**.
  - `density_model`: `PURE_OIL` -> **`CONSTANT`**, `MASS_VOLUME_MIXING` -> **`MIXTURE`**
    (now documents both phase densities as explicit inputs: `rho` and
    `dissolved_gas_liquid_density`).
  - `viscosity_model` -> key **`liquid_viscosity_model`**; `PURE_OIL` -> **`CONSTANT`**,
    `EMPIRICAL_CORRELATION` -> **`EMPIRICAL`**.
  - `gas_mixture_viscosity_model` -> key **`mixture_viscosity_model`**; `*_QUALITY`/
    `*_VOID` suffixes dropped (`MCADAMS`, `LINEAR`, `DUKLER`).
  - `bc_z_*_thermal`: `OPEN` -> **`ZERO_GRADIENT`**, `RESERVOIR` -> **`CONSTANT`**.
- **Pruned** (with graceful alias-to-survivor on parse):
  - Liquid viscosity `LOG_MIXING` (it is `EMPIRICAL` with the T and p coefficients 0).
  - Solubility `BUNSEN` (the same linear-in-p law as Henry); the unused
    `dissolved_gas_bunsen_coeff` key now warns and is ignored.
  - Mixture viscosity `EINSTEIN_DILUTE` (thickens with gas, wrong as alpha->1) and
    `KRIEGER_DOUGHERTY` (suspension packing law). Survivors: `DUKLER` (legacy,
    kept for comparison), `LINEAR` (Grando), `MCADAMS`. **Default is now `MCADAMS`**
    (was `EINSTEIN_DILUTE`).
  - Convection `TVD_MINMOD` (van Leer dominates it); `TVD_VANLEER` -> `VANLEER`.
- **Added**: convection scheme **`LINEAR`** (2nd-order central, unbounded) as a
  reference scheme alongside `UPWIND`/`VANLEER`/`TYPE_DIFFERENCING`, in
  `couette_face_value` (reynolds.cpp) and `FVM::divergence` (fvm.cpp).
- **Kept as labeled references** (per the "keep the simplest model for comparison"
  principle): `FULL_SOMMERFELD` (already implemented as the no-clamp `Reynolds::solve`
  path) and all five `motion_time_method` integrators (the non-default ones tagged
  `[reference]` in the GUI).
- GUI (`panels_setup.cpp`) combos reordered/relabeled with the new names and hover
  text; `schema.cpp` drops the dead Bunsen entry. PHYSICS.md gains **§17 "Configurable
  Model Options"** (the canonical reference the GUI hover text mirrors). README and
  all example configs updated to the new tokens.

### 2026-06-15 - WP-11 follow-up: reformation reflood carries reservoir state

Reason: complete the WP-11 boundary-flux unification. The shared
`elrod_boundary_outflow` already re-floods a cavitated axial-boundary cell at the
physical reformation rate in the *liquid* (Elrod) balance, but the energy and
dissolved-gas convection still gated their boundary face by the full-film switch,
so the reflood inflow carried no temperature or composition — the reformation
liquid entered the liquid balance but was thermally/compositionally invisible.

- **Shared reflood velocity** (`SimulationConfig::is_reformation_boundary`,
  `reformation_boundary_velocity`): one Poiseuille reflood velocity driven from the
  cavitation plateau `(p_bc - p_cav)/(½Δz)`, matching the inflow
  `elrod_boundary_outflow` imposes on the liquid balance (the bulk modulus cancels,
  leaving `ρ_l h³/(12μ)·dp/dz`). The liquid, energy, and gas equations now see one
  reformation rate.
- **Energy** ([energy.cpp](src/energy.cpp)): a cavitated boundary cell under
  `consistent_boundary_flux` re-floods at the reformation rate; the inflowing
  reservoir liquid carries `temperature_reference` (RESERVOIR mode) through the
  existing inflow/outflow handling. Full-film and flag-off paths unchanged.
- **Dissolved gas** ([gas_transport.cpp](src/gas_transport.cpp)): the reflood
  convects the reservoir composition `dissolved_gas_initial` into the cell and uses
  the full liquid density `rho_liquid_solution` (not the cavitation-reduced Elrod
  density `ρ = ρ_l·θ`), so the gas reflood mass equals the Elrod liquid reflood.
- **Diagnostics** ([diagnostics.cpp](src/diagnostics.cpp)): `gas_boundary_inflow`
  mirrors the same reflood, so the reported `gas_boundary_flux` tracks the solver
  at a cavitated end instead of reading zero.
- **Test** (new `test_reformation_boundary_reservoir_state` in
  [test_diagnostics.cpp](tests/test_diagnostics.cpp)): with the flag off a cavitated
  boundary cell stays depleted (0.05) and the diagnostic gas flux is zero; with it
  on the cell re-saturates to the reservoir composition (0.2175) and the diagnostic
  mirrors the reflood as a positive inflow.
- Default off (`consistent_boundary_flux = false`) is bit-identical to before.
  Full suite: 22/22.

### 2026-06-15 - WP-1 acceptance gates + WP-11 shared boundary flux

Reason: harden the two newest model extensions per the AUDIT_PLAN acceptance
criteria.

- **WP-1 acceptance tests** (new tests/test_gas_coupling.cpp): 0-D closed-cell
  release conserves total propane and reaches a steady state; under-saturated
  resorption; the vaporous limit (no gas) reproduces the single-phase Elrod
  solve bit-for-bit; and a coupled pressure-cap gate where VOID_COUPLED relieves
  a supersaturated low-pressure arc up to the bubble point (679 kPa -> 1.0 MPa)
  while NONE leaves it depressed with volumetrically-inert gas. These gates
  caught and fixed two real bugs in the
  VOID_COUPLED kernel: the mixture bulk modulus must weight the LIQUID volume
  fraction (1-alpha_g), not the Elrod theta, so beta_bar -> beta exactly at
  alpha_g = 0; and a gas-free cell now short-circuits void_params to the exact
  single-phase parameters (avoiding a 1/(1/beta) round-off that broke the
  no-gas==single-phase identity).
- **WP-11 shared flux** (new SimulationConfig::elrod_boundary_outflow): the
  axial-boundary liquid mass flux is now defined once and used by both
  solve_elrod's assembly and the diagnostics mass balance. With
  consistent_boundary_flux on, the cavitated reformation rate is accounted in
  the balance, so the liquid residual drops from ~1.4e-3 to ~7e-8. Flag off is
  bit-identical to before. Full suite: 22/22.

### 2026-06-14 - WP-2/WP-1/WP-6/WP-11: coupling, two-phase pressure, striated zone, boundary flux

Reason: complete the oil-refrigerant cavitation model per PLAN.md Phase E. Every
work package is gated and defaults to the prior behaviour (full ctest stays
green bit-for-bit; the R290 table case runs end-to-end with all of them on).

- **WP-2 — outer Picard coupling.** STEADY_STATE now iterates the
  property/Reynolds/energy sweep to a pressure-space residual instead of a single
  segregated pass (`coupling_max_iters`/`coupling_tolerance`/`coupling_relaxation`).
  TRANSIENT advances one sweep per step exactly as before. The silent `gas_dt = 0`
  steady freeze is replaced by `steady_gas_model = FROZEN | EQUILIBRIUM` with a
  validation guard (`allow_frozen_gas_steady` to opt out); EQUILIBRIUM uses a new
  `FluidProperties::equilibrate_gas` local flash. `tests/test_coupling.cpp`.
- **WP-1 — released gas coupled into pressure** (`gas_pressure_coupling =
  NONE | VOID_COUPLED`, default NONE). The EOS is generalized to
  `p = p_void + beta_bar ln(theta/theta_full)` with single-phase back-compat
  overloads; the Elrod solve uses per-cell `theta_full = 1 - alpha_g` (void
  ceiling), `beta_bar` (mixture compressibility `1/beta_bar = theta/beta +
  alpha_g/p`), and a local `p_void` plateau. Gaseous rupture at p_sat then emerges
  mechanistically. NONE reproduces master bit-for-bit; VOID_COUPLED is stable and
  shifts the load as the released gas softens the film. `tests/test_eos.cpp`.
- **WP-6 — striated cavitated zone** (`cavitated_film_model = FULL_FILM |
  STRIATED`, default FULL_FILM). Under STRIATED only the liquid fraction theta
  carries Couette shear and heat generation in cavitated cells (friction torque
  and viscous dissipation are no longer overstated by ~1/theta). Requires
  ELROD_ADAMS. `tests/test_coupling.cpp`.
- **WP-11 (partial) — cavitated boundary reformation** (`consistent_boundary_flux`,
  default off). A cavitated axial-boundary cell re-floods at the physical
  Poiseuille rate `cs*(p_bc - p_cav)/beta` (an explicit inflow) rather than the
  barotropic Gamma_base Dirichlet link, which overstates re-flooding by
  ~beta/(p_bc-p_cav). Default off = identical assembly. REMAINING: a single shared
  boundary-flux function across the Reynolds/energy/gas solves and the diagnostics
  mass balance (the balance still uses the old link, so r_l rises with the flag on).
- New smoke config `config_r290_pz68_table.txt` exercises the table + LOCAL_PSAT +
  VOID_COUPLED + STRIATED stack end-to-end.

### 2026-06-14 - WP-12: 2-D PTSV property tables from file + p_sat cavitation onset

Reason: ingest the measured PZ68S/R290 PTSV data (a P-T grid) and make cavitation
onset track the dissolved-refrigerant bubble point, so the higher-release-pressure
species governs rupture (general physics, not hard-coded to R290). Implements the
new WP-12 in PLAN.md Phase E; default-off so existing results are unchanged.

- **2-D property tables** (`PropertyTable2D` in `config.hpp`, bilinear with edge
  clamp). New `src/property_table_io.{hpp,cpp}` loads a `P_MPa,T_C,<value>` CSV
  (header-mapped, BOM/CRLF/comment tolerant, SI unit conversion, grid resolution
  taken from the file so it is editable post-compile, on-load validation).
- Config keys `solubility_table_file` / `viscosity_table_file` /
  `density_table_file` (resolved against the config dir then CWD); loaded in
  `main.cpp` before validation; a 2-D file takes precedence over an inline 1-D
  table (warned). `data/R290_PZ68S/r290_pz68s.csv` is the shipped surface.
- The property layer consumes the 2-D tables: solubility `c_sat(p,T)`, density
  `rho(p,T)`, and viscosity as KINEMATIC `nu(p,T)` converted per cell to dynamic
  `mu = nu * rho_solution`. `SimulationConfig::saturation_pressure` inverts the
  2-D surface along the local isotherm via a first-crossing scan (robust to the
  saturated flat tail and measured-data noise); shared `invert_csat` for 1-D/2-D.
- **Cavitation threshold** `cavitation_threshold = SCALAR_PCAV | LOCAL_PSAT`
  (default SCALAR_PCAV). `SimulationConfig::effective_p_cav()` returns
  `max(p_cav, p_sat(c_d, T))` under LOCAL_PSAT (a uniform datum -> no spurious
  net force), threaded through the Elrod/Gumbel onset, recovery, boundary, inlet,
  and force-gauge sites. For the R290 case p_cav=1e5 but p_sat=1.0 MPa, so propane
  outgasses first.
- Tests: new `tests/test_property_tables.cpp` (interpolation, CSV load incl. the
  real data spot-check (1 MPa,40 C)->21.746%/6.843 mm2/s, malformed-file errors,
  2-D solubility/p_sat integration, effective_p_cav). Full suite 20/20 green
  (SCALAR_PCAV reproduces master). Smoke config `config_r290_pz68_table.txt`
  runs end-to-end with the table + LOCAL_PSAT (p_sat onset = 1.0 MPa).

### 2026-06-14 - Documentation consolidation; oil-refrigerant cavitation plan

Reason: the audit/review cycle left several standalone working docs at the repo
root, duplicating physics and plan content. Consolidated per the user's request.

- **Docs merged + removed.** `AUDIT.md`, `AUDIT_PLAN.md`, `REVIEW.md`,
  `REVIEW_VERDICT.md`, and the 2025 model-changes note were folded into the
  canonical files and deleted. Their verbatim originals are preserved in git
  history at commit `ae567d6` (one snapshot commit before removal).
  - `PLAN.md` Phase E now inlines the full WP-1..WP-11 specifications (physics,
    code touch-points, test acceptance gates), the review priority tiers, and
    the audit DOI reference list — previously only referenced from
    `AUDIT_PLAN.md`. Added **WP-12** (2-D property tables from an editable file
    + `p_sat`-driven cavitation onset).
  - `PHYSICS.md` gained **§16 Oil–Refrigerant Gaseous Cavitation** (two-threshold
    taxonomy, generalized JFO complementarity with `p_sat`, PTSV-data closures,
    model provenance and honest validity). The cited research basis stays in
    `docs/CAVITATION_OIL_REFRIGERANT.md`.
- **New research doc** `docs/CAVITATION_OIL_REFRIGERANT.md`: literature-verified
  survey of oil+R290 cavitation (gaseous/outgassing governed by the solution
  bubble-point), with status tags and a verification ledger.
- No code changes in this entry; WP-12/WP-1/WP-2/WP-11/WP-6 implementation
  follows per `PLAN.md` Phase E.

### 2026-06-12 - Light theme, regrouped form, physical Monitors panel

Reason: user review round 2.

- Light visual theme (neutral grays, blue accent); solver-log and status
  colors retuned for the light background.
- Parameter groups split (Mesh & time -> Time stepping + Mesh) and ordered
  by edit frequency; frequently-tweaked groups open at the top, the rest
  collapsed by default; the filter forces matching groups open.
- Convection scheme dropdown reduced to Upwind / Linear (type differencing)
  per user decision; TVD limiter values loaded from existing config files
  are shown and preserved; thermal/gas scheme rows removed from the form
  (keys still parsed/round-tripped). Raw config.txt tab removed; presets
  kept.
- Convergence panel replaced by Monitors: selectable physical series
  (h_min [um], T min/max [K], |F| [N], friction torque [N.m], cavitated
  area [%], liquid residual on a log axis) in stacked subplots with real
  engineering units and a linked time axis - different ranges live on
  separate axes instead of one squashed log plot. Data parsed live from
  DETAILED step lines and diagnostics.csv; merged CSV export.

### 2026-06-12 - GUI refinements: units, equal-aspect contour, option help

Reason: user feedback on the modernized GUI.

- Contour view clamped to the contour region exactly: the plot widget is
  sized to the domain aspect (axis decorations measured per frame), so no
  area outside the field is ever shown; the Fit button moved into the
  toolbar row beside Field / Cavitation zone.
- Inlet zone overlay toggle added next to the cavitation toggle: cells with
  inlet_indicator = 1 are tinted light blue, the hover readout reports
  "inlet region", and a hint appears when the case has inlets but the
  inlet_indicator field is not in output_fields.

- MPI process count moved next to Run/Cancel (visibility: more processes =
  faster solve); removed from the menu.
- Field map now renders in cell space with an equal-aspect constraint: the
  cell ratio shows as-is and window resizes zoom instead of stretching.
  Pan/zoom is constrained near the domain; a `Fit` button in the plot's
  upper-right corner restores the full-domain view. Tick labels remain
  physical (deg / m).
- Run History panel removed at user request (the key-results summary still
  uses the last run's record internally). Profile overlays went with it.
- Control labels moved to the left of every selector/checkbox.
- Per-field display units (Pa/kPa/MPa, deg/rad, rad/s-rpm, m/mm/um, s/ms)
  with in-place value conversion; the config file stays solver-native SI and
  choices persist in `%LOCALAPPDATA%\PancakeGui\settings.txt`.
- Hover descriptions on all dropdown options. Convection schemes expose
  upwind / linear (type differencing, central in the full film, V&K 1989) /
  TVD with accuracy and when-to-use notes. A standalone pure-central
  "linear" scheme for thermal/gas transport does not exist in the solver and
  would be a solver-side addition.

### 2026-06-11 - GUI modernization: Dear ImGui docking + ImPlot rewrite

Reason: executing `GUI_MODERNIZATION_DIRECTIVE.md` (Phase 0 scan, redesign,
standalone packaging). Full report: `AGENT_REPORT.md`; build instructions:
`BUILD.md`; design: `docs/GUI_ARCHITECTURE.md`.

- Replaced the single-file owner-drawn Win32 GUI (`src/gui_win32.cpp`, no
  longer built) with a docking ImGui + ImPlot application in `src/gui/`
  (Win32 + Direct3D 11 backends, fetched as pinned releases via CMake
  FetchContent).
- Solver subprocess management moved off the UI thread: a pipe-reader thread
  (ANSI parse, progress), a `diagnostics.csv` tail thread (live residuals),
  and an async results parser; cancellation terminates the subprocess with a
  partial-output warning.
- New capabilities required by the directive: live log-scale convergence
  plot, pressure heatmap with cavitation-zone overlay and hover readout,
  profile plots with multi-run overlay, in-session sortable run history,
  CSV export for all tables/curves, PNG export for all plots, key-results
  summary with copy-to-clipboard.
- Validation is now exactly the solver's: `validate_config()` merges
  `SimulationConfig::validate()` so Run is disabled iff the solver would
  reject the case; fields highlight in red with constraint tooltips.
- Standalone packaging: static CRT (`-static` / `/MT`), embedded Roboto +
  JetBrains Mono fonts (`cmake/embed_binary.cmake`), per-monitor-v2 DPI
  manifest + version resource, all runtime writes in
  `%LOCALAPPDATA%\PancakeGui\`. Import table verified to contain only
  Windows-shipped DLLs.
- New presets `release` (MSVC, untested here - no MSVC on the build machine)
  and `release-mingw` (verified); `windows-native-mingw` unchanged and still
  builds solver + GUI + tests.
- Solver source untouched (the GUI reuses header-only `src/config.hpp`).

### 2026-06-11 - Phase E remediation: WP-4 diagnostics, WP-3 numerics honesty, WP-5 mixture viscosity

Reason: first implementation pass of the audit remediation plan
(`AUDIT_PLAN.md`), executing the three no-dependency work packages in the
prescribed order.

**WP-4 — Conservation, regime, and convergence diagnostics (`src/diagnostics.hpp/cpp`, new)**
- Per-step global liquid and total-gas mass balances that mirror the assembled
  discrete operators: storage from `theta` vs `theta.old` with the assembled
  capacity weight and squeeze source, boundary fluxes from the same
  Dirichlet/INLET_OUTLET formulas, and inlet penalty cells re-evaluated as
  physical operator residuals (the injected mass). Results append to
  `<output_dir>/diagnostics.csv` every `diagnostics_interval` steps together
  with Elrod outer-iteration/flag-flip counts and the cavitated-cell fraction.
- Clamp/snap operations are no longer silently non-mass-neutral: the theta_min
  clamp and full-film snap accumulate their per-cell mass into
  `cavitation_clamp_mass`, and gas-transport clamps/supply-composition
  overwrites accumulate into `gas_clamp_mass`; both enter the balance.
- `Reynolds::solve_elrod` now returns an `ElrodStats` struct (outer iterations,
  cumulative flag flips, final max theta change, converged flag) used by the
  per-step log and CSV.
- `SimulationConfig::validate()` regime guards: Couette Reynolds number and
  Taylor-number warning above the laminar Taylor-vortex onset (Ta > 41), and
  gaseous-saturation context — `SimulationConfig::saturation_pressure(c_d, T)`
  (Henry inversion including van't Hoff slope; Bunsen and table inversions)
  is logged with `p_cav/p_sat`, and boundary or inlet feeds at/above `p_sat`
  warn that the oil enters saturated. For the shipped R290/PZ68 constants
  p_sat = 1.0 MPa, exactly the 10 bar feed.
- New `linear_rtol` config key (default 1e-8, back-compatible) threads through
  all Krylov solves; the balance residual sits at this tolerance (verified at
  1e-13 → residuals O(1e-12), the PLAN A.2 target).
- Tests: `tests/test_diagnostics.cpp` — closed-domain conservation,
  open-domain boundary-flux accounting, inlet-source accounting, clamp
  visibility, closed-cell gas-exchange closure (total propane conserved to
  1e-12), and the Taylor/saturation guard warnings.
- Documentation reconciliation (AUDIT §B1/B3): PHYSICS.md §6 now states the
  axial Dirichlet value is the configured absolute `bc_z_*_val` (shipped cases
  are submerged, not `p = p_cav` vented ends) and carries the
  constant-cavity-pressure caveat (Etsion & Ludwig 1982; Braun & Hendricks
  1984); the §13.2.1 R290 table now matches the shipped `p_cav = 3e4` /
  `bc = 1e6` values. New PHYSICS.md §15 documents the diagnostics.

**WP-3 — Numerics honesty (`fvm`, `reynolds`, `energy`, `gas_transport`, `config`, GUI)**
- Convection schemes are now selectable and live: `theta_convection_scheme`
  (UPWIND default, TVD_VANLEER, TVD_MINMOD, TYPE_DIFFERENCING),
  `thermal_convection_scheme`, `gas_convection_scheme` (UPWIND/TVD). The
  previously dead TVD machinery in `FVM::divergence` is wired up;
  TYPE_DIFFERENCING (Vijayaraghavan & Keith 1989) applies a deferred central
  correction on full-film faces and pure upwind across/inside cavitation,
  via face masks. Deferred corrections in the Elrod outer loop lag one outer
  iteration through a scratch carrier so `theta.old` stays the previous
  timestep for the capacity term.
- Dead time selectors removed: `pressure_time_method` and
  `temperature_time_method` never altered the backward-Euler field solves and
  are deleted from the config struct, writer, shipped configs, and GUI; the
  parser accepts the old keys with an "ignored" warning collected in
  `SimulationConfig::parse_warnings` and logged at startup.
- Capacity-term gate fixed (audit's silent-physics finding): the
  `rho0*h*dtheta/dt` liquid-content storage is assembled for every transient
  Elrod solve, fixed or moving bearing; previously a fixed-bearing transient
  (e.g. an `omega_ramp_time` startup) was quasi-steady in theta. The Gumbel
  squeeze source gate likewise no longer requires MOVING_BEARING.
- Gumbel/Elrod wedge alignment: the Gumbel Couette source is now a face-flux
  divergence with rho*h factorized into a central rho0*h face value times the
  scheme-interpolated content theta = rho/rho0, exactly mirroring the Elrod
  Couette flux. Matched full-film comparison agrees to ~0.01% peak pressure
  (historical central-vs-upwind mismatch was ~1.2%).
- Legacy steady-style tests (`test_elrod_a1`, `test_elrod_a2`) now declare
  `solution_mode = STEADY_STATE`, preserving their original quasi-steady
  semantics under the corrected transient gate; `test_elrod_bc` and
  `test_journal_motion` updated for the removed keys.
- Tests: `tests/test_schemes.cpp` (+ 6-rank `test_schemes_mpi6`) — TVD
  boundedness and sharpness vs upwind on step advection, observed refinement
  orders (~1.09 upwind, ~1.99 TVD van Leer), the matched Gumbel/Elrod
  full-film comparison, finite-relaxation/steady-limit checks for the capacity
  term, TVD ghost-layer consistency across ranks, and config back-compat.
- Deliverable: `validation/grid_convergence/` (config + run.ps1/run.sh +
  README) running 60x20/120x40/240x80 on a fed journal case; archived run
  shows monotone convergence of cavitation extent (5.17% -> 4.83% -> 4.80%)
  with balance residuals at ~1e-12.

**WP-5 — Mixture viscosity with gas-limit asymptote (`fluid_properties`, `config`, GUI)**
- New `gas_mixture_viscosity_model` dispatch replaces the hardcoded Einstein
  factor: EINSTEIN_DILUTE (back-compat default), DUKLER_VOID (void-weighted
  linear), MCADAMS_QUALITY (homogeneous two-phase standard, mass quality
  x = m_g/(m_g + rho_l*theta*h)), KRIEGER_DOUGHERTY (packing-limited), and
  LINEAR_QUALITY (Grando 2006 replication). All satisfy mu(alpha=0) = mu_l;
  the void/quality models thin toward the new `mu_gas` key. The free-gas
  volume-fraction clamp uses `gas_alpha_max` instead of a hardcoded 0.99.
- Validation: EINSTEIN_DILUTE is rejected with `gas_alpha_max > 0.6`
  (Krieger-Dougherty packing bound); DUKLER/MCADAMS/LINEAR require
  `mu_gas > 0`.
- Shipped R290/PZ68 configs select MCADAMS_QUALITY with `gas_alpha_max = 0.6`
  and saturated-propane-vapor `mu_gas = 8.7e-6 Pa s`, and document that the
  Henry/log-viscosity tangents are valid near (10 bar, 40 C) only — measured
  isotherm tables are pending (flagged, not fabricated).
- Tests: `tests/test_fluid_properties.cpp` extended with per-model liquid-limit
  and gas-limit asymptotes, monotonicity, and both validation guards.

GUI: removed the dead pressure/thermal time-step combos; added theta/thermal/
gas scheme dropdowns, the free-gas viscosity model combo, and the `mu_gas`
field with inline validation mirroring the CLI guards.

Full ctest suite green: 19/19 (16 prior + test_diagnostics + test_schemes +
test_schemes_mpi6).

### 2026-06-11 - Audit and remediation planning documents

Reason: a multi-agent physics/plan review identified model-level defects in the
gaseous-cavitation path (released gas volumetrically inert to the pressure
solve; fixed vaporous `p_cav` for a gas-saturated film), plus numerics and
documentation drift. The findings were verified against the literature and the
remediation was decomposed into agent-executable work packages.

**`REVIEW.md`, `REVIEW_VERDICT.md`, `AUDIT.md`, `AUDIT_PLAN.md`, `PLAN.md`**
- `REVIEW.md`: two-agent review thread (Reviewer A findings W1-W11/M1-M4;
  Reviewer B code-grounded verdicts plus new findings N1-N9).
- `REVIEW_VERDICT.md`: consolidated, prioritized issue list from both rounds.
- `AUDIT.md`: literature-grounded audit. Each issue carries verified file:line
  code facts, quoted primary sources (ASME J. Tribology, STLE Tribology
  Transactions, Wear, Tribology International, Int. J. Refrigeration, Proc.
  IMechE Part J), a verdict, and a prescribed fix. Headline: couple released
  gas into film continuity (published effect: load -3% at eccentricity 0.3,
  +21% at 0.9) and track gaseous void onset at p_sat(T, c_d) instead of a
  fixed sub-bar vaporous threshold. Also documents a fact missed by the
  reviews: the Elrod capacity term only assembles for MOVING_BEARING runs, so
  fixed-bearing transients are quasi-steady in theta.
- `AUDIT_PLAN.md`: ten work packages (WP-1..WP-10) with physics equations,
  citations, code touch-points, and test cases with acceptance gates.
- `PLAN.md`: added Phase E status table referencing `AUDIT_PLAN.md`; notes
  that WP-1/WP-2 complete the coupling Phase B.4 originally specified.

No solver code was changed in this pass.

### 2026-06-08 - Pressure-derived Elrod boundaries and R290/PZ68 rescue pass

Reason: the interrupted gaseous-cavitation implementation still exposed manual
axial film-content boundary inputs (`bc_z_*_theta`) and the default config
enabled a gas-density model without the physical constants needed to keep the
liquid density meaningful. That contradicted the pressure/bulk-modulus boundary
contract and made the default case too easy to corrupt.

**`src/config.hpp`, `src/reynolds.cpp`, `src/energy.cpp`**
- `ELROD_ADAMS` axial boundary film content is now derived internally from
  `bc_z_*_val`, `p_cav`, and `bulk_modulus`. Legacy `bc_z_*_theta` keys remain
  parser-compatible but no longer drive Reynolds, velocity reconstruction,
  force post-processing, energy boundary fluxes, or initial pressure.
- `GasTransport::solve` now returns immediately when `dt <= 0`, avoiding a
  singular boundaryless steady advection solve in steady operating-point runs.

**`src/gui_win32.cpp`, `config.txt`, `config_r290_pz68.txt`**
- Removed axial film-content boundary controls from the GUI and generated
  configs. The GUI now round-trips the full gas-property coefficient set:
  Henry temperature factor, dissolved-gas liquid density, dissolved-gas
  diffusivity, Andrade temperature viscosity coefficient, and Barus pressure
  viscosity coefficient.
- Reset `config.txt` to a conservative pure-oil default with absolute
  `p_cav = 1e5` instead of gauge zero. The calibrated R290/PZ68 starting case
  now uses absolute pressure consistently (`p_cav = 1e5`, 10 bar supply/reference)
  and writes under the ignored `results/` tree.
- Added `config_r290_pz68_quick.txt`, copied beside the solver and GUI at build
  time, for a short gaseous-cavitation smoke run.

**Tests and docs**
- Added shared `SimulationConfig::validate()` coverage and wired the validator
  into the CLI and GUI save/run path. Under-specified gas models now fail before
  launch instead of silently using blank constants.
- Output-directory setup now reports locked result folders cleanly instead of
  terminating on an uncaught `std::filesystem` exception. This is the expected
  failure mode when ParaView or another viewer still has `results/...` open.
- Corrected the gas bookkeeping so full-film pressure/inlet cells no longer
  delete `free_gas_mass` simply because JFO `theta >= 1`. Free gas is now
  transported conservatively as `free_gas_mass`, remains visible through
  `alpha_gas`, and disappears through negative `gas_mass_transfer` only when the
  liquid is undersaturated at the local pressure/temperature.
- Added gas-composition boundary handling for pressure inlets and axial open
  ends. Pressure supply/inlet cells now impose `dissolved_gas_initial`, zero
  `free_gas_mass`, zero `alpha_gas`, and zero local `gas_mass_transfer`.
  Axial pressure-boundary inflow carries the same initial dissolved-gas
  composition with zero free gas; axial outflow convects released gas out.
- Updated Elrod and energy boundary tests to assert pressure-derived theta and
  ignored legacy theta inputs.
- Added R290/PZ68 reference-point property coverage and a zero-timestep
  gas-transport no-op regression. Added regressions for inlet-composition
  enforcement, axial dissolved-gas inflow composition, and axial free-gas
  outflow removal.
- Updated README, PHYSICS, and PLAN to document the pressure-only boundary
  contract, the calibrated R290/PZ68 constants, and vectorized VTK output names.

### 2026-06-07 - Outlet initialization, omega ramp startup, and gas-property GUI controls

Reason: the transient startup should no longer force a steady first pressure
solve. Starting from the outlet-compatible initialized field and ramping shaft
speed is a clearer numerical experiment and avoids hiding startup dynamics.

**`src/config.hpp`, `src/main.cpp`**
- `SimulationConfig::initial_pressure()` now seeds uniform pressure from the
  axial outlet condition instead of the first inlet supply pressure. For
  `ELROD_ADAMS`, it recovers pressure from the outlet `bc_z_*_theta` film
  content through the EOS.
- Added `omega_ramp_time`; transient runs use
  `omega(t) = omega_target * clamp(t / omega_ramp_time, 0, 1)`. Steady-state
  runs ignore the ramp and use the target speed directly.
- Removed the forced steady `t=0` transient solve and the first-step flooded
  refill. Transient terms are active from the first transient solve.

**`src/gui_win32.cpp`, `config.txt`**
- Added GUI control for `omega_ramp_time`.
- Added a Fluid Properties section with selectable
  `fluid_property_model`, dissolved-gas species, solubility/density/viscosity
  models, gas release/resorption coefficients, bounds, reference state, and
  optional property tables.
- The live summary reports ramp status and the selected property/gas model
  controls, and warns about missing tabulated data or zero gas-transfer rate
  when a gas model is selected.

**Tests and docs**
- Updated inlet/startup regression coverage for outlet-based initialization and
  linear omega ramping.
- Updated `README.md`, `PHYSICS.md`, and `PLAN.md` to describe the reverted
  first-step startup behavior and the new ramp/property controls.

### 2026-06-07 - Preserve Elrod transient capacity for explicit pressure mode

Reason: moving-bearing transient runs with `pressure_time_method =
EULER_EXPLICIT` could develop a second, nonphysical cavitated region after the
groove even though the steady `t=0` solution was correct. The Elrod transient
assembly skipped the `rho*h*dtheta/dt` capacity term for explicit pressure mode
but still applied the lagged `theta*dh_dt` squeeze source. That is not a
conservative explicit discretization of `d(theta*h)/dt`; it lets the squeeze
source act without the old liquid-content storage term that should resist
instantaneous depletion.

**`src/reynolds.cpp`**
- Always assembles the Elrod weighted transient capacity term for non-steady
  moving-bearing solves. The squeeze term remains lagged with `theta.old`, so
  `EULER_EXPLICIT` still means the geometry source is explicit, not that
  pressure/film-content storage is removed.

**`tests/test_elrod_bc.cpp`**
- Added a regression for `pressure_time_method = EULER_EXPLICIT` with zero
  `dh_dt`, verifying the transient solve retains the old `theta` state through
  the capacity term instead of collapsing immediately to the steady boundary
  solution.

**`README.md`, `PHYSICS.md`, `PLAN.md`**
- Documented that transient Elrod capacity is mandatory for mass conservation
  and is assembled for all pressure time-method labels.

### 2026-06-07 - Elrod MPI active-set ghost synchronization

Reason: 6-rank runs showed rank-aligned pressure/film-content blocks. The
Elrod-Adams outer active-set loop solved new physical `theta` values, but the
next outer iteration recomputed the cavitation switch using stale cross-rank
`theta` ghost cells. That made the full-film/cavitated switch inconsistent at
MPI partition interfaces.

**`src/reynolds.cpp`**
- Synchronizes `theta` ghost cells after every Elrod outer solve, clamp, and
  near-full-film snap before recomputing the active-set switch or exiting the
  solver.

**`tests/test_elrod_a2.cpp`, `tests/CMakeLists.txt`**
- Added a post-solve ghost-consistency regression that compares `theta` ghost
  cells against neighboring ranks' current physical boundary cells.
- Added a 6-rank `test_elrod_a2_mpi6` CTest entry to cover the rank count that
  exposed the contour blocking.

**`README.md`, `PHYSICS.md`, `PLAN.md`**
- Documented that Elrod active-set ghost exchange happens inside the nonlinear
  outer loop, not only between timesteps.

### 2026-06-07 - Live-summary Courant number and project contribution summary

**`src/gui_win32.cpp`**
- The live summary now shows `Co_theta`, a nominal circumferential Couette
  Courant number computed from the current `omega`, `dt`, and `n_theta` inputs:
  `Co_theta = |omega| dt / (2 dtheta)`, with `dtheta = 2*pi/n_theta`.
- Increased the summary card height so the added derived quantity remains
  visible with the existing default summary rows.

**`docs/PROJECT_CONTRIBUTIONS.md` (new), `README.md`, `PHYSICS.md`,
`docs/WINDOWS_PORT.md`**
- Added a project-contribution summary that lists substantive solver, model,
  output, testing, and documentation work while excluding one-off fixes and GUI
  construction work.
- Documented that the displayed Courant number is a Couette-advection input
  diagnostic and does not include pressure-driven Poiseuille velocities.

### 2026-06-07 - Face-consistent cavitation-front heat generation

Reason: the Elrod-Adams pressure field intentionally has a hard cavitation
plateau (`p = p_cav` where `film_content < 1`), but the derived Poiseuille heat
source must not treat inactive cavitation-front faces as a valid cell-centered
pressure gradient. That produced mesh-dependent hot spots near the corners of
the full-film/cavitated interface. This is not phase-change heat; the default
case has gas transport disabled.

**`src/energy.cpp`**
- Changed `heat_generation` pressure-flow dissipation to average squared
  pressure gradients over active full-film faces. Cavitated faces contribute
  zero Poiseuille heat, while physical axial boundary cells use the adjacent
  interior face rather than the boundary pressure. This matches the face-gated
  pressure-driven velocity without reintroducing boundary heat spikes.

**`tests/test_energy.cpp`**
- Added a cavitation-front regression with two connected full-film cells
  surrounded by cavitated cells. The expected heat is the half-face active
  contribution on each full-film cell; the old one-sided gradient would be
  twice as large.

**`README.md`, `PHYSICS.md`, `PLAN.md`**
- Documented that the sharp Elrod pressure/film-content contour is the JFO
  active set, while `heat_generation` is now face-consistent post-processing
  rather than a gas phase-change source.

### 2026-06-07 - Dissolved-gas transport (gaseous cavitation, JFO Phase B)

**`src/gas_transport.{hpp,cpp}` (new)**
- Advection-diffusion transport of the dissolved-gas mass fraction `dissolved_gas`
  (c_d) using the SAME film mass flux (rho*h*u, gated Couette+Poiseuille) that the
  Elrod-Adams solve conserves, in non-conservative (constant-preserving) form.
  Operator splitting: transport, then `FluidProperties::update_gas_state` does the
  local Henry's-law release/resorption to/from free gas. A later 2026-06-08
  correction transports `free_gas_mass` and caps `alpha_gas` by `gas_alpha_max`
  instead of deleting free gas in full-film JFO cells. Wired into `main.cpp`
  after the property refresh and ghost sync.
- New config `dissolved_gas_diffusivity` [m^2/s] (effective Fickian diffusivity).

**`src/fluid_properties.cpp` -- (P,T)-dependent property models**
- Free gas now thickens the mixture viscosity (Einstein dilute suspension,
  mu = mu_l(1 + 2.5 alpha_g)); dissolved gas still thins it via the log-mixing law.
- Henry solubility gained a van't Hoff temperature factor
  (`dissolved_gas_henry_temp_coeff` E_H): c_sat = H_ref exp[E_H(1/T-1/T_ref)] p.
- EMPIRICAL_CORRELATION viscosity is now mu_oil*exp[E_mu(1/T-1/T_ref)] (Andrade T,
  `viscosity_temperature_coeff`) * exp[a_c c_d] (log-mixing) * exp[alpha_p(p-p_ref)]
  (Barus, `viscosity_pressure_coeff`).
- MASS_VOLUME_MIXING density uses a configurable dissolved-gas liquid-phase density
  (`dissolved_gas_liquid_density`) instead of the ideal-gas density (which is wrong
  for dissolved gas); 0 keeps the legacy behaviour.
- `config_r290_pz68.txt` pre-fills every constant; the (P,T) models reproduce the
  chart point (nu ~= 6.84 cSt vs 6.843 cSt at 40 C / 1 MPa, c_d=0.2175).

**`config_r290_pz68.txt` (new)**, **`tests/test_gas_transport.cpp` (new)**
- Calibrated R290/PZ68 example (Henry H=2.175e-7 /Pa, viscosity a_c=-11.40, 40 C /
  1 MPa) reproducing the chart point (solubility 21.75 %, nu=6.843 cSt). Tests cover
  constant-preservation, no-op for non-mixture models, and bounded mass-conserving
  advection.

### 2026-06-07 - Per-step logging, verbosity option, and matching GUI contour names

**`src/main.cpp`, `src/config.hpp`**
- The time loop now logs **every** timestep (decoupled from the VTK save cadence,
  which still follows `write_interval`); written steps are tagged `[saved]`.
- Added `output_verbosity = SIMPLIFIED | DETAILED`. SIMPLIFIED is a compact
  `Step N t=...` line; DETAILED adds `h_min`, temperature range, fluid force, and
  friction torque per step. Honored by both the CLI and (via the captured stdout)
  the GUI log; the GUI progress bar still tracks the `t=` token.

**`src/gui_win32.cpp`**
- Added an Output "log detail" combo (SIMPLIFIED/DETAILED) with config round-trip.
- The contour/preview field selector now lists the **output names** (`p`, `T`,
  `U_X/U_Y/U_Z`, `Fp_X`, `film_content`, `h`, ...) so they match the VTK arrays and
  ParaView; the reader resolves `<vector>_X/_Y/_Z` to the array + component and the
  preview field is no longer lower-cased (mixed-case `U`/`Fp`/`T` resolve correctly).

### 2026-06-07 - OpenFOAM-style VTK output names and Cartesian vector fields

Reason: the VTK output mixed long names (`pressure`, `temperature`) with redundant
per-component scalar arrays (`velocity_x/y/z/theta`, `*_force_x/y/z`), which is hard
to read in ParaView.

**`src/output_naming.hpp` (new)**
- Single source of truth mapping internal field names to output names
  (`pressure`->`p`, `temperature`->`T`, velocity->`U`) and grouping component triples
  into Cartesian vectors: `Fp` (pressure force), `Fv` (viscous force), `F` (fluid
  force), `Fext` (external load), `xB` (bearing position).

**`src/io.cpp`**
- Writes scalars under their OpenFOAM names and each enabled vector group as a single
  3-component Cartesian `DataArray`; per-component scalar arrays are no longer emitted.
  Flat output uses `U = (u_theta, u_z, 0)` (the unwrapped-plane Cartesian frame), curved
  output uses `U = (vx, vy, vz)`. `CellData`/`PCellData` mark `Vectors="U"`.

**`src/gui_win32.cpp`**
- Preview keeps its descriptive combo labels but resolves each selection to the new
  array name + component (`velocity_theta` -> `U[0]`, `pressure_force_x` -> `Fp[0]`, ...)
  so the heatmap preview still works against the renamed/vectorized output.

**`tests/test_io.cpp`**
- Asserts the unified `U`/`Fp` vectors and the absence of per-component scalar arrays.

### 2026-06-06 - Fix axial-boundary pressure consistency and add inlet thermal modes

Reason: the energy solve, velocity reconstruction, and z-force used the raw
`bc_z_*_val` as the axial boundary pressure, but the Elrod-Adams solve pins the
axial boundary on film content (`bc_z_*_theta`). With a submerged config
(`bc_z_*_val = 1.5e6`, `bc_z_*_theta = 1.0`) this fabricated a 1.5 MPa axial
gradient at the ends: `velocity_z` showed a non-physical inflow at both ends and
the near-boundary cells were dragged to `temperature_reference` (a cold drop).
Removing the phantom value then exposed a corner heat spike where the pinned
1.5 MPa supply groove met a vented end.

**`src/energy.cpp`, `src/reynolds.cpp`**
- Added `solved_boundary_pressure()`: for Elrod-Adams the axial boundary pressure
  is `pressure_from_theta(bc_z_*_theta)`; Gumbel keeps `bc_z_*_val`. Used in the
  energy boundary fluxes / `dp/dz` diagnostics, in `calculate_velocities`, and in
  the `calculate_macroscopic_properties` z-shear, so post-processing matches the
  solved field instead of the configured pressure value.

**`src/energy.cpp`, `src/config.hpp`**
- Added `ThermalInflowMode` (`OPEN` default / `RESERVOIR`) per axial side
  (`bc_z_*_thermal`). `OPEN` models a submerged / same-oil boundary: inflow is
  zero-gradient (no externally fixed temperature). `RESERVOIR` carries
  `temperature_reference` in on inflow (an actual fed inlet/outlet).
- Inlets gained an optional supply temperature
  (`inlet_* ... t_supply` -> `feeds_fresh_oil`). Open large-clearance regions
  (no `t_supply`) stay pressure-only and thermally transparent; fed inlets pin
  their cells to `t_supply` in the energy solve.

**`config.txt`**
- The submerged sample now encodes the 1.5 MPa bath at the axial ends through
  `bc_z_*_theta = 1.0015011` (= `exp(1.5e6/1e9)`) and `bc_z_*_thermal = OPEN`, so
  the ends and the open groove share the same bath pressure and recirculating
  oil. Result: temperature stays in a physical ~356-367 K band (previously
  301 K cold cells and a 4073 K spike), ends allow both leakage outflow and bath
  inflow without spurious heating/cooling.

**`tests/test_energy.cpp`, `tests/test_forces.cpp`**
- Tests that prescribed a pressure field consistent with `bc_z_*_val` now make
  the boundary pressure self-consistent with the solver (set `bc_z_*_theta` via
  the EOS, or use Gumbel) and select `RESERVOIR` where reservoir inflow is under
  test.
- Added `test_elrod_axial_boundary_uses_solved_pressure` regression: a full-film
  interior below `bc_z_*_val` must be axial outflow (no `temperature_reference`
  dragged in) under Elrod-Adams.

**`src/gui_win32.cpp`**
- Added Axial-Boundaries controls for `bc_z_south_thermal` / `bc_z_north_thermal`
  (OPEN / RESERVOIR) and Inlet controls for the open-region vs fed-inlet choice
  plus a supply-temperature field, with full config round-trip.

**`tests/test_inlets.cpp`, `tests/test_elrod_bc.cpp`**
- Recalibrated penalty-pin and theta-Dirichlet assertions to the linear-solver
  residual floor (rtol = 1e-8): inlet pressures are checked relative to
  `p_supply`, and `theta`/pressure boundary checks use a `beta`-scaled tolerance.
  The pins are accurate to ~5e-11 relative; the previous sub-floor absolute
  tolerances (`1e-6`, `1e-3`, `1e-8`) failed on the MinGW PETSc build although
  the solutions are correct. A real BC leak still trips them by many orders.

### 2026-06-05 - Thermal boundary regressions for energy solve

Reason: previous temperature calculations were known to be fragile near axial
pressure boundaries, and the new field-valued property path needs boundary
coverage too.

**`src/energy.cpp`**
- Boundary heat generation now treats `INLET_OUTLET` as a pressure boundary in
  `dp/dz` diagnostics, matching the thermal boundary-flux logic.

**`tests/test_energy.cpp`**
- Added north-boundary inflow coverage to mirror the existing south-boundary
  inflow test.
- Added a pure-outflow boundary case to verify `temperature_reference` is not
  imposed on outflow cells.
- Added an `INLET_OUTLET` boundary heat-generation test that also verifies
  field-valued `mu` is used at the boundary.

**`PHYSICS.md`**
- Documented the `INLET_OUTLET` boundary heat-generation treatment.

### 2026-06-05 - JFO-compatible propane/oil mixture-property layer

Reason: the solver needs oil plus dissolved-gas mixture properties before
gaseous cavitation transport, while preserving the Elrod-Adams/JFO
liquid-content variable as the mass-conserving cavitation backbone.

**`src/config.hpp` / `config.txt`**
- Added fluid-property enums and config keys:
  `fluid_property_model`, `dissolved_gas_species`,
  `oil_gas_solution_model`, `density_model`, and `viscosity_model`.
- Added dissolved-gas, finite-rate transfer, gas alpha limit, reference
  pressure/temperature, and `x:value` table inputs for solubility, density, and
  viscosity.
- Default remains `CONSTANT` pure oil. The default config and GUI no longer
  advertise RK pressure/temperature method semantics, though the parser still
  accepts old aliases.

**`src/fluid_properties.hpp` / `src/fluid_properties.cpp`**
- Added liquid-solution property evaluation for pure oil, Henry/Bunsen/table
  solubility, mass-volume or table density, and log/empirical/table viscosity.
- Added ideal-gas density for propane/R290 and air.
- Added finite-rate release/resorption between `dissolved_gas` and
  `free_gas_mass`, with `alpha_gas` capped by the configured free-gas limit.

**`src/main.cpp` / `src/reynolds.cpp` / `src/energy.cpp`**
- Added property and gas diagnostic fields:
  `rho_liquid_solution`, `mu_liquid_solution`, `cp_liquid_solution`,
  `rho_cp_liquid_solution`, `k_liquid_solution`, `rho_gas`, `mu`,
  `dissolved_gas`, `free_gas_mass`, `alpha_gas`, and `gas_mass_transfer`.
- Fed field-valued density/viscosity into Gumbel, Elrod-Adams, velocity,
  force/torque, heat generation, and thermal convection calculations.
- Kept the update as a segregated JFO/property/gas corrector, not PISO.

**`src/main.cpp` / `src/reynolds.cpp` / `src/config.hpp` / `src/gui_win32.cpp`**
- Removed runtime allocation/filling of legacy `load_x/load_y/load_z`.
- Legacy output requests map to `pressure_force_x/y/z`; the GUI exposes the new
  property/gas diagnostic output fields and preserves raw property config.

**`tests/test_fluid_properties.cpp` / `tests/test_io.cpp` / `tests/test_forces.cpp`**
- Added property-model tests for pure-oil fallback, solubility monotonicity,
  density bounds, viscosity model selection, table interpolation, saturated
  zero source, release, resorption, and bounded `alpha_gas`.
- Added an IO regression for legacy `load_x` mapping and updated force
  diagnostics to avoid legacy load fields.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the mixture-property API, equations, JFO compatibility, deferred
  transport work, and legacy output-field mapping.

### 2026-06-04 - Steady startup step for transient solves

Reason: transient pressure and temperature terms should not compare the first
physical solve against arbitrary initial guesses. The solved `t=0` state should
seed `old_time` before transient terms are applied.

**`src/main.cpp` / `src/reynolds.cpp`**
- The timestep loop now passes a step-local config with
  `solution_mode = STEADY_STATE` for `step == 0`, even when the user selected a
  transient run.
- Reynolds squeeze/transient assembly now respects `solution_mode`, so step 0
  omits pressure transient terms and the `t=dt` solve uses the solved step-0
  fields as previous time.

**`tests/test_energy.cpp`**
- Added a regression that sets an intentionally wrong previous temperature and
  verifies the steady initial energy solve ignores it.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the startup sequence: steady `t=0`, store solved old fields, then
  transient solve from `t=dt`.

### 2026-06-04 - Corrected energy boundary treatment

Reason: the first energy-equation implementation only tested uniform cases and
used conservative thermal advection. Near axial pressure boundaries and supply
grooves, the Reynolds solve has local mass source/sink behavior, so the
conservative `div(F*T)` form introduced a spurious `T div(F)` term and produced
large nonphysical hot/cold cells.

**`src/energy.cpp`**
- Added axial thermal boundary convection. Pressure-boundary inflow now uses
  `temperature_reference`; pressure-boundary outflow uses an upwind
  zero-gradient temperature treatment.
- Pressure supply inlet cells remain pressure-only constraints and no longer
  impose a temperature Dirichlet condition.
- Converted thermal advection to the documented non-conservative form by
  subtracting the flux-divergence term from the upwind conservative assembly.

**`tests/test_energy.cpp`**
- Added boundary-focused checks that pressure supply inlet cells are not fixed
  to a thermal Dirichlet value and that axial pressure-boundary inflow remains
  bounded by the reference temperature.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the corrected inlet/axial thermal boundary model and the
  non-conservative advection correction.

### 2026-06-04 - Active vector velocity output

Reason: ParaView should read `velocity` immediately as a vector field, while
users still need scalar components for direct plotting.

**`src/io.cpp`**
- Marked `velocity` as the active VTK cell vector with
  `CellData/PCellData Vectors="velocity"` in both curved and flat output.
- Kept the 3-component `velocity` array and added scalar component arrays.
  Curved output writes `velocity_x`, `velocity_y`, `velocity_z`, and
  `velocity_theta`; flat output writes `velocity_theta` and `velocity_z`.

**`tests/test_io.cpp` / `tests/CMakeLists.txt`**
- Added an output metadata regression that checks `.vts` and `.pvts` files for
  active vector metadata, a 3-component `velocity`, and component scalar arrays.

### 2026-06-04 - Energy equation and steady-state solution mode

Reason: the next development step is thermal hydrodynamics, and fast case
checks need a single steady operating-point solve instead of a full transient
loop.

**`config.txt` / `src/config.hpp` / `src/main.cpp`**
- Synced the root config with the build-directory operating conditions:
  `dt = 1e-6`, `write_interval = 1e-5`, RK2 time-method selections,
  moving-bearing mode, independent external load, and the smaller minimum film
  thickness limit.
- Added `solution_mode = TRANSIENT|STEADY_STATE`. `STEADY_STATE` runs one
  pressure/thermal step at `t=0` and skips motion advancement.
- Added `temperature_model = ISOTHERMAL|ENERGY_EQUATION`, thermal properties,
  wall temperatures, wall heat-transfer coefficients, and default output for
  `temperature` and `heat_generation`.
- Added `temperature` and `heat_generation` fields to the solver loop and call
  the energy solve after pressure, velocity, and force updates.

**`src/energy.hpp` / `src/energy.cpp`**
- Added a constant-property film-averaged energy equation using the existing
  FVM diffusion, upwind convection, weighted transient, and source assembly.
- Computes viscous heat generation from Couette shear plus pressure-flow shear
  and stores it as `heat_generation` in W/m^2.
- Applies journal and bearing wall heat-transfer sinks/sources for steady and
  transient thermal solves. `ISOTHERMAL` leaves `temperature` fixed while still
  refreshing `heat_generation`.

**`src/gui_win32.cpp`**
- Added a clean Energy section for thermal model, temperature, material, and
  wall heat-transfer inputs.
- Added `solution_mode` selection in Mesh & Time and thermal preview units for
  `temperature` and `heat_generation`.

**`tests/test_energy.cpp` / `tests/CMakeLists.txt`**
- Added regression coverage for thermal config aliases, isothermal heat refresh,
  and the steady Couette wall-balance scale.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the energy equation, steady-state mode, GUI thermal controls, and
  Phase D.1 implementation status.

### 2026-06-04 - Bearing-side force convention and geometry

Reason: the force outputs mixed shaft-side and bearing-side conventions. This
made `load_*`, `viscous_force_*`, and `fluid_force_*` hard to interpret and
could give surprising signs/scales when plotting or using them for motion.

**`src/reynolds.cpp` / `src/reynolds.hpp`**
- Reworked `calculate_macroscopic_properties` so pressure, viscous, and total
  fluid forces are all forces applied by the fluid to the moving bearing
  surface.
- Pressure integration now uses the bearing surface area vector
  `((R+h)e_r - h_theta e_theta - (R+h)h_z e_z) dtheta dz`, so circumferential
  and axial film-thickness slopes are represented.
- Viscous shear now uses bearing-side signs:
  `tau_theta = mu omega R / h - h/(2R) dp/dtheta` and
  `tau_z = -h/2 dp/dz` in full-film regions.
- `fluid_force_*` is now `pressure_force_* + viscous_force_*`, not a negated
  mixed shaft/bearing quantity.

**`src/config.hpp` / `src/main.cpp` / `src/gui_win32.cpp` / `config.txt`**
- Added `pressure_force_x`, `pressure_force_y`, and `pressure_force_z` output
  fields. Legacy `load_x`, `load_y`, and `load_z` remain aliases for the
  pressure-only resultant but are no longer part of the default output list.

**`tests/test_forces.cpp` / `tests/CMakeLists.txt`**
- Added regression coverage for pressure force surface geometry, axial viscous
  shear sign/area, and shaft friction-torque scaling.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the force names and bearing-side convention.

### 2026-06-04 - Inlet-pressure initialization

Reason: the pressure field's flooded initial guess should start from the inlet
supply pressure, not the cavitation pressure.

**`src/config.hpp`**
- Added `SimulationConfig::initial_pressure()`, returning the first configured
  inlet `p_supply` and falling back to `p_cav` only when no inlet exists.

**`src/main.cpp`**
- Initializes `pressure` from `cfg.initial_pressure()` at field creation and
  before each timestep solve.
- Initializes the Elrod `theta` guess from the same pressure through the EOS, so
  the pressure and film-content fields are consistent.

**`tests/test_inlets.cpp`**
- Added a regression assertion for the first-inlet initialization rule.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the inlet-pressure flooded guess and the no-inlet fallback.

### 2026-06-04 - GUI load controls and preview refresh behavior

Reason: the moving-bearing keys were available in `config.txt`, but not all of
them had clean GUI inputs. The preview also refreshed contours during solver
polling and recursively scanned result folders while switching steps, which made
interactive use noisy and slow.

**`src/gui_win32.cpp`**
- Added a Motion / Loads rail section for `motion_model`, pressure/motion/
  temperature time-method selection, external load magnitude/direction/z,
  bearing initial state, support stiffness/damping, mass, minimum film
  thickness, and stop-on-nonpositive-film settings.
- Converts the GUI external-load magnitude and in-plane direction into
  `external_load_x` and `external_load_y` on save, while preserving
  `external_load_z` as an independent component.
- Stopped refreshing the contour preview during solver timer polling. The
  current contour remains visible while a run is active and refreshes when the
  solver exits or when the user clicks Refresh.
- Replaced per-load recursive VTS discovery with a cached per-step list of all
  `processor*` files under `flat/` or the curved output root. Field and timestep
  switching now reuses that list.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented the Motion / Loads GUI inputs, load-vector conversion, stable
  preview behavior during runs, and cached preview file indexing.

### 2026-06-04 - Moving bearing dynamics and force-separated output

Reason: dynamic bearing studies need the film forces to move the bearing
housing while the shaft remains fixed. The previous roadmap used "journal
motion" wording, but the requested model is the opposite sign convention:
bearing-center displacement changes the gap as
$h = c + x_b\cos\theta + y_b\sin\theta$ and the equivalent attitude angle is
computed from the minimum-gap direction.

**`src/config.hpp` / `config.txt`**
- Added `MotionModel` with `STATIC` and `MOVING_BEARING`.
- Added `TimeSteppingMethod` with `EULER_EXPLICIT`, `EULER_IMPLICIT`,
  `CRANK_NICOLSON`, `RK2`, and `RK4`. The parser also accepts common aliases:
  `EXPLICIT`, `IMPLICIT`, `SEMI_IMPLICIT`, `IMEX`, `CRANK_NICHOLSON`,
  `MIDPOINT`, and `RUNGE_KUTTA_*`.
- Added moving-bearing state, mass, support stiffness/damping, independent
  external load, minimum film-thickness safety, and pressure/motion/temperature
  time-method config keys.
- Extended default output fields with `viscous_force_*`, `fluid_force_*`,
  `external_load_*`, and `bearing_*` fields.

**`src/journal_motion.hpp` / `src/journal_motion.cpp`**
- Added a moving-bearing state module. Initial state can be derived from
  `e`/`attitude_angle_deg` so `MOVING_BEARING` starts from exactly the same
  film thickness as the static eccentricity model.
- Implemented motion updates for explicit Euler, implicit Euler,
  Crank-Nicolson, RK2, and RK4 under lagged fluid force plus independent
  external load.
- Added equivalent eccentricity/attitude helpers and uniform output-field
  writers for bearing state and external load.

**`src/film_thickness.hpp` / `src/film_thickness.cpp`**
- Added `compute_moving_bearing`, using the fixed-shaft/moving-bearing sign
  convention and computing `dh_dt` for squeeze-film coupling.

**`src/reynolds.hpp` / `src/reynolds.cpp`**
- `calculate_macroscopic_properties` now returns a `ForceComponents` result.
- Kept legacy `load_x/load_y/load_z` as pressure-only force on the shaft.
- Added integrated viscous force from circumferential and axial shear.
- Added `fluid_force_x/y/z` as the equal-and-opposite total force applied to
  the moving bearing, separate from configured `external_load_x/y/z`.
- Added moving-gap squeeze-film source terms for Gumbel and Elrod; Elrod uses
  implicit `rho*h*dtheta/dt` for non-explicit pressure methods and lagged
  `theta*dh_dt` as an explicit source.

**`src/main.cpp`**
- Added moving-bearing state to the timestep loop. Each step updates `h` and
  `dh_dt`, solves Reynolds, computes forces, writes force/state outputs, and
  advances the bearing state for the next step.
- Added a minimum film-thickness stop condition to prevent running into
  nonphysical negative gaps.

**`src/gui_win32.cpp`**
- Preserves and emits the new time-method and moving-bearing config keys.
- Added output checkboxes and preview units for the new force and bearing-state
  fields.

**`tests/test_journal_motion.cpp` / `tests/CMakeLists.txt`**
- Added regression coverage for moving-bearing sign convention, attitude
  recovery, config aliases, and constant-load motion integrators.

**`README.md` / `PHYSICS.md` / `PLAN.md`**
- Documented moving-bearing dynamics, force separation, time-method names,
  squeeze-film coupling, and the updated Phase C status.

### 2026-05-31 - Added future journal-motion and ETHD roadmap

Reason: previous git history showed the energy equation only in the discarded
bubble-solver plan, while the journal-bearing roadmap did not include dynamic
journal motion or ETHD. The current physics documentation also only described
the static-eccentricity, isothermal Reynolds/Elrod model.

**`PLAN.md`**
- Added Phase C, Dynamic Journal Motion, covering time-dependent film geometry,
  squeeze-film transient Reynolds coupling, and a two-degree-of-freedom journal
  force-balance model.
- Added Phase D, Thermal and ETHD Implementation, covering the lubricant energy
  equation, thermoviscous Reynolds coupling, elastic/thermal film-thickness
  deformation, and the fully coupled ETHD Picard loop.
- Extended the phase summary table with the new dynamic-motion and ETHD phases.

**`PHYSICS.md`**
- Added planned dynamic journal-motion equations for $h(\theta,z,t)$,
  $\partial h/\partial t$, the expanded Elrod transient term, and rigid-journal
  force balance.
- Added planned ETHD equations for depth-averaged energy transport, viscous
  dissipation, wall heat loss, temperature-dependent viscosity, elastic
  compliance, thermal expansion, and coupled iteration order.

### 2026-05-29 - GUI redesign to a single workspace + single-file binary

Reason: the previous tab-based GUI felt dated and finicky. Live combo boxes were
wrapped in `BS_GROUPBOX` frames, which forced constant `bring_control_to_front`
/ z-order hacks (combos sometimes only painted after a click), the workflow was
split across `Workspace` / `Output` / `Raw Config` tabs, and there was no live
feedback — bad input was only reported by a `MessageBox` at save time. The goal
was a layout that lets engineers change parameters and run quickly while
preventing misunderstanding.

**`src/gui_win32.cpp`** (UI layer rewritten from scratch; solver-process,
VTS-preview, ANSI-log, and unit-conversion backends preserved):
- Replaced the tabbed shell with one tab-less workspace: action/run bar on top,
  a resizable input rail on the left with a pinned live-summary card, a result
  viewport in the center, and a resizable console at the bottom (two drag
  splitters).
- Removed all `BS_GROUPBOX` frames around live controls. Section headers and
  field labels are now painted in the rail's `WM_PAINT`; only interactive
  controls are child windows. This eliminates the z-order/paint bugs and the
  `raise_*` / `bring_control_to_front` workarounds entirely.
- Input rail: collapsible sections, a search filter, and a Basic/`Advanced`
  toggle for progressive disclosure.
- Live summary card recomputes derived journal-bearing quantities on every edit
  (`ε = e/c`, `h_min = c(1−ε)`, `h_max`, `c/R`, `L/D`, `U = ωR`, previewed-field
  mean/max).
- Inline validation: impossible/malformed inputs (`e ≥ c`, non-positive
  grid/clearance, `dt`/`write_interval`/`end_t` ordering, empty output names)
  paint the field red and disable `Run` until resolved, replacing the
  save-time-only `MessageBox`.
- Run state: owner-drawn `Run`/`Stop`, a state chip with elapsed time, and a
  progress bar driven by parsing `t=<value>` from solver stdout against `end_t`.
- Hybrid units: a global SI ↔ engineering (`mm`/`MPa`/`rpm`) selector sets
  defaults, with per-field override dropdowns kept only where multiple units are
  realistic; fixed-unit fields use a plain painted suffix (down from ~30 unit
  combos). Saved config stays in solver-native units.
- Scenario presets: Default bearing, High eccentricity, Misaligned, Circular
  oil hole.
- Custom config files: `config.txt` beside the exe is the default; `Open…` /
  `Save As…` use the standard Windows file dialogs (`comdlg32`). The loaded file
  name and a dirty `*` marker appear in the title bar.

**`CMakeLists.txt`**
- `pancake_gui` now links statically (`-static` under MinGW, static CRT under
  MSVC) and links `comdlg32` + `ole32`, so it ships as a single self-contained
  `pancake_gui.exe` with no runtime DLLs beside it. Dropped the
  `pancake_deploy_mingw_runtime(pancake_gui)` call (the solver still deploys its
  PETSc/MPI/MinGW DLLs).

### New files

**`src/gui_win32.cpp`**
- Added a native Win32 GUI for Windows builds.
- The GUI loads `config.txt` from the same directory as `pancake_gui.exe`, falls
  back to hardcoded `SimulationConfig` defaults when the file is missing, saves
  edits back to that file, and can launch `pancake.exe` from the same directory.
- Replaced the raw-editor-only prototype with workspace, output, and raw-config
  views.
- Added process launch/stop, captured solver logs, output folder/PVD open
  actions, and a GDI heatmap preview for the latest ASCII VTK scalar field.
- Reworked the main workspace into a split-pane engineering page:
  the result preview expands with the window, the process log stays below it,
  and a draggable horizontal splitter lets users tune the pane sizes.
- Lifted solver backend controls into the always-visible top-right header so
  the backend can be run/stopped from any tab.
- Moved the parameter editor into the right side of the main `Workspace` tab as
  a scrollable one-row-per-parameter inspector.
- Streamlined solver launch around the native `pancake.exe` backend built by
  the MSYS2 MinGW64 workflow.
- Preserved literal numeric edit text in GUI saves so typed or file-provided
  values such as `0.0010` stay visible as `0.0010`.
- Added preview timestep selection with previous/next controls.
- Added a heatmap colorbar with min/max labels for the current preview field,
  then expanded it with a field title, clearer number formatting, and plot axes.
- Fixed multi-rank preview loading by assembling the selected timestep from all
  matching `processor*` VTS files instead of treating rank 0 as a complete
  field.
- Added `load_angle_deg`, a visualization reference angle measured from the
  positive y axis, and rotated the preview display relative to that reference.
- Reworked axial boundary controls into South and North groups, and renamed the
  theta boundary inputs to film content value.
- Added per-value GUI unit selectors for length, time, angle, tilt, inlet
  pressure, and omega while still saving solver-native config units.
- Extended unit selectors across the full parameter inspector, including
  pressure fields, grid counts, dimensionless controls, and inlet supply
  pressure.
- Preserved raw inlet token text in the structured editor so values such as
  `1.5e6` are not rewritten as `1500000` after saving.
- Switched the process log to a RichEdit control that renders ANSI color codes
  instead of showing raw `[32m` / `[0m` fragments, with plain text fallback if
  RichEdit is unavailable.
- Normalized raw-config edit text to Windows CRLF newlines so multiline config
  text displays correctly.
- Replaced the raw multiline inlet editor with a structured groove/circular
  selector that comments out the unused inlet form when saving.
- Fixed the all-grey/low-contrast appearance: added a light modern palette and
  `WM_CTLCOLORSTATIC` / `WM_CTLCOLOREDIT` / `WM_CTLCOLORLISTBOX` handlers so
  labels, group boxes, and check boxes share one clean surface and editable
  fields stay white. Previously these controls fell back to the system 3D-face
  grey on the window background.
- Subclassed the tab control's `WM_ERASEBKGND` so its display area is painted
  with the page colour instead of the default grey slab behind the child
  controls.
- Forced every combo box above the interleaved group-box frames after creation,
  fixing the preview field/step selectors that only painted after a click.
- Reworked the preview annotations: wider right margin so colorbar min/max
  values are no longer truncated, the field name as a centred colorbar title,
  a stats subtitle at the top, and a decluttered axis footer.
- Moved preview title/stat annotations farther above the contour so min/max and
  grid metadata do not overlap the plotted field.
- Restored `film_content` to write and preview the raw Elrod/JFO universal
  variable instead of clamping it to `[0, 1]`, because pressure recovery depends
  on the same theta values and compressed full-film cells must remain visible.
- Added a regression check that the Elrod A.2 journal-bearing case contains
  both cavitated theta values below one and compressed full-film theta values
  above one.
- Centered the preview z-axis title, labeled the z-axis from `0` to `L`, and
  expanded the colorbar to include units plus min/max-inclusive tick labels.
- Fixed Elrod axial DIRICHLET boundaries so the solver uses the dimensionless
  `bc_z_south_theta` / `bc_z_north_theta` film-content values instead of
  converting the pressure-valued `bc_z_*_val` entries into theta. Added sanity
  coverage for this separation.
- Changed the default bulk modulus to `1e9` Pa so a 1.5 MPa supply pressure
  produces film-content values near 1.0015 rather than millions in the default
  Windows GUI workflow.
- Added blue/magenta to the ANSI log colour table so every code the solver
  emits (`\033[31..36m`, `\033[0m`) is colourised rather than only stripped.

**`src/pancake.ico` / `src/pancake_logo.svg` / `src/resource.h` / `src/pancake.rc`**
- Added the Pancake logo as a Windows application icon resource. The editable
  SVG mirrors the supplied blue film-profile logo, while the generated ICO is
  embedded into both `pancake.exe` and `pancake_gui.exe` for File Explorer,
  taskbar, and title-bar branding.

**`src/pancake_gui.manifest` / `src/pancake_gui.rc`**
- Added a `CREATEPROCESS_MANIFEST_RESOURCE_ID` manifest that declares the
  ComCtl32 v6 dependency. The previous `#pragma comment(linker, ...)` was guarded
  by `_MSC_VER`, so MinGW builds shipped without it and fell back to the classic
  (grey, Windows-95 era) common controls. Embedding the manifest via a compiled
  `.rc` enables modern visual styles under MinGW (and MSVC).

**`CMakePresets.json`**
- Kept a single Windows preset, `windows-native-mingw`, for building native
  `pancake.exe`, `pancake_gui.exe`, tests, and runtime DLL deployment from the
  MSYS2 MinGW64 PETSc/MS-MPI stack.

### Modified files

**`CMakeLists.txt` / `tests/CMakeLists.txt`**
- Split solver sources into `pancake_core`, kept `pancake` as the CLI target,
  and added the Windows-only `pancake_gui` target.
- Embedded the shared `pancake.ico` resource into both Windows executables, so
  the GUI and console backend present the same application logo.
- Enabled the `RC` language for the GUI target and added `src/pancake_gui.rc`
  (with its `INCLUDE_DIRECTORIES` source property) so the ComCtl32 v6 manifest
  is compiled in by windres/rc and visual styles work on every Windows compiler.
- Added PETSc discovery through `pkg-config`, including MSYS2 PETSc package
  names, while preserving the existing `PETSC_DIR` fallback.
- Restricted native Windows PETSc auto-discovery to MPI-enabled package names
  so CMake does not select the serial `petsc-sso` variant that collides with
  MS-MPI headers.
- Added MinGW runtime deployment for native Windows builds so required DLLs
  such as `libpetsc-dmo.dll`, `libstdc++-6.dll`, and `libopenblas.dll` are
  copied beside `pancake.exe`, `pancake_gui.exe`, and test executables.
- Updated the Windows CMake preset to prepend `C:\msys64\mingw64\bin` to PATH,
  preventing native compiler helper DLL lookup failures in plain PowerShell.

**`cmake/deploy_mingw_runtime.cmake`**
- Added a small post-build deployment helper that runs MSYS2 `ldd`, filters
  dependencies from the MinGW runtime directory, and copies them next to the
  target executable for File Explorer launches.

**`src/config.hpp` / `src/main.cpp`**
- Added config serialization for GUI defaults and same-directory config loading
  for the CLI default path.
- Added output persistence settings: `output_write_3d`, `output_write_flat`,
  and `output_fields`.
- Added `load_angle_deg` as a visualization-only load-reference angle.

**`src/io.cpp`**
- Honours the new output persistence settings by skipping disabled VTK/PVD
  products and disabled field arrays.

**`config.txt`**
- Added the new output persistence keys to the checked-in default config and
  included the axial inlet/outlet theta selectors emitted by config
  serialization.

**Markdown documentation**
- Rewrote the README, BUILDING guide, and Windows-port notes around the current
  MSYS2 MINGW64 workflow. Visual Studio is documented as unnecessary for the
  current native solver + GUI path.

**`.gitignore`**
- Removed stale Visual Studio GUI build-directory entries and kept the active
  `build-windows-mingw/` workflow ignored.

---

## 2026-05-26 — Cross-Platform Windows Port Assessment and Documentation

### New files

**`ASSESSMENT.md`**
- Added a Windows-port readiness assessment with file/line-specific blockers, warnings, dependency findings, and effort estimates.
- Documented the current PETSc implementation versus the stated Trilinos strategy as the primary dependency blocker.

**`BUILDING.md` / `ARCHITECTURE.md` / `CONTRIBUTING.md`**
- Added cross-platform build instructions, intended solver/GUI layering, and contribution conventions for Linux and Windows.
- Made `std::filesystem::path`, non-POSIX core code, explicit integer widths, and GCC/MSVC compatibility part of the contributor contract.

**`docs/WINDOWS_PORT.md`**
- Added an ordered checklist for the first native Windows build, starting with WSL-side cleanup and ending with Windows validation.

**`.gitattributes`**
- Added text/binary and line-ending rules for cross-platform development.

### Modified files

**`README.md` / `CLAUDE.md`**
- Updated build and platform guidance so future contributors and agents see the Linux/Windows dual-build strategy immediately.

**`.gitignore`**
- Added `build-linux/` and `build-windows/` as dedicated out-of-tree build directories.

---

## 2026-05-21 — Macroscopic Properties: Load Capacity and Friction Torque

### Modified files

**`src/reynolds.hpp` / `src/reynolds.cpp` — Post-processing calculations**
- Added `Reynolds::calculate_macroscopic_properties` to calculate the integrated load capacity vector ($F_x, F_y, F_z$) and the scalar friction torque $T$ over the full 2D mesh domain.
- Computes $F_z$ properly taking into account the surface normal variation caused by shaft tilt.
- Uses `MPI_Allreduce` to correctly sum the integrated values across all decomposed sub-domains.

**`src/main.cpp` — Global field allocation**
- Added four new uniform scalar fields (`load_x`, `load_y`, `load_z`, `friction_torque`) to the visualization output.
- Hooked up `calculate_macroscopic_properties` directly in the solver's main time loop.

**`PHYSICS.md` — Documentation**
- Appended a new section detailing the mathematical integration of the load capacity from the pressure field and the friction torque from the Couette and Poiseuille shear stress formulas.

---



## 2026-04-12 — Phase A.8: Configurable Outlet Boundary Conditions

### New files

**`tests/test_bc.cpp`**

Unit test for axial boundary conditions. Validates:
1. **Neumann (Zero Gradient)**: One end DIRICHLET and the other NEUMANN results in a uniform pressure field.
2. **Dirichlet (Fixed Value)**: Both ends DIRICHLET with different values results in a linear pressure profile.

### Modified files

**`src/config.hpp` — Boundary condition definitions**

- Added `BCType` enum (`DIRICHLET`, `NEUMANN`).
- Added axial boundary parameters to `SimulationConfig`: `bc_z_south_type`, `bc_z_south_val`, `bc_z_north_type`, `bc_z_north_val`.
- Updated `load_from_file` to support these parameters.

**`src/reynolds.cpp` — Boundary condition logic**

- Updated `Reynolds::solve` and `Reynolds::solve_elrod` to conditionally apply boundary terms based on the configured `BCType`.
- In `calculate_velocities`, the pressure gradient at boundaries is now correctly determined by the BC type (zero gradient enforced for `NEUMANN`).

---

## 2026-04-12 — Phase A.7: Staggered Grid Velocity Fields

### New files

**`tests/test_velocity.cpp`**

Unit test for the staggered grid velocity calculation. Verifies:
1. **Couette Flow**: In a concentric bearing with no pressure gradient, $u_\theta$ correctly matches $\omega R / 2$ and $u_z$ is zero.
2. **Poiseuille Flow**: A prescribed axial pressure gradient produces the analytical $u_z = -(h^2/12\mu) \partial p/\partial z$ response.

### Modified files

**`src/field.hpp` / `src/field.cpp` — Staggered grid support**

- `Field` constructor now correctly sizes physical dimensions for `FACE_THETA` and `FACE_Z`. 
- `FACE_Z` fields are allocated with $N_z + 1$ physical rows to account for the non-periodic axial faces.

**`src/io.cpp` — Staggered visualization**

- Completed `get_centered_data` for `FACE_THETA` and `FACE_Z`.
- Face-based velocities are now interpolated to cell centers during VTK output, ensuring they can be visualized alongside pressure and theta.

**`src/reynolds.hpp` / `src/reynolds.cpp` — Velocity calculation**

- Implemented `Reynolds::calculate_velocities`.
- $u_\theta$ is calculated at $\theta$-faces using a centered pressure gradient and averaged film thickness.
- $u_z$ is calculated at $z$-faces, including the bearing ends (Dirichlet $p_{cav}$ boundary).
- The Poiseuille (pressure-driven) component is only active in full-film regions ($\theta \ge 1$), consistent with the Elrod-Adams model.

**`src/main.cpp`**

- Added `velocity_theta` and `velocity_z` fields.
- Call `calculate_velocities` at each timestep after the Reynolds solve.

---

## 2026-04-12 — Phase A.6: Flexible Inlet Boundary Conditions

### New files

**`tests/test_inlets.cpp`**

Unit test for the internal Dirichlet (penalty) boundary conditions. Validates that both axial grooves and circular supply holes correctly force the prescribed supply pressure in the target cells for both GUMBEL and ELROD_ADAMS models.

### Modified files

**`src/config.hpp` — Inlet definitions**

- Added `InletConfig` struct.
- Added `std::vector<InletConfig> inlets` to `SimulationConfig`.
- Updated `load_from_file` to support multiple `inlet_circular` and `inlet_groove` definitions from the text config.

**`src/reynolds.cpp` — Penalty method for inlets**

- Implemented `apply_inlet_conditions` helper using a large-diagonal penalty method ($10^{20}$). This effectively locks inlet cells to their supply values without requiring invasive changes to the sparse matrix topology.
- Wired the helper into `Reynolds::solve` and the `Reynolds::solve_elrod` outer loop.

**`config.txt`**

Added example inlet definitions (one axial groove at 180° and one circular hole at 90°) to demonstrate usage.

---

## 2026-04-12 — Phase A.4: Text-based configuration and Phase A.5: Starvation fix

### New files

**`config.txt`**

Created a text-based configuration file containing all simulation parameters. This allows users to modify the simulation without recompiling. The format is a simple `key = value` structure with support for `#` comments.

### Modified files

**`src/config.hpp` — Added `load_from_file`**

Updated `SimulationConfig` with a `load_from_file(const std::string& path)` method. It parses the text file and updates the struct fields accordingly. This satisfies the user request for an editable config file.

**`src/main.cpp` — Load config and re-initialise fields**

- Added a call to `cfg.load_from_file("config.txt")` at the start of `main`.
- **Numerical Starvation Fix**: In the time loop, `theta` and `pressure` are now reset to full-film values ($\theta = 1$, $p = p_{cav}$) before every `solve_elrod` call. This provides a "flooded guess" that allows oil from the $z$-boundaries to reach the interior of the domain, preventing the "theta=0 everywhere" trap.

**`src/reynolds.cpp` — Flooded initialization in outer loop**

Modified the Elrod-Adams outer loop to force $g = 1$ for the very first iteration (`outer == 0`). This ensures that the solver starts with a fully active diffusion field, establishing a healthy pressure zone that "pulls" lubricant into the bearing from the boundaries. Subsequent iterations correctly use the `switch_function(theta)` to develop the cavitation zone.

**`PLAN.md` — Updated with Phase A.4 and A.5**

Updated the development plan to reflect the completion of the text-based configuration and the numerical starvation fix.

---

## 2026-04-12 — Phase A.2: Elrod-Adams switch function and mass-conserving cavitation

### New files

**`tests/test_elrod_a1.cpp` — updated regression test**

Added `cfg.max_outer_iters = 1` for the ELROD_ADAMS branch of `test_regression_vs_gumbel`. The intent of this test is Phase A.1: compare the full-film θ-formulation (g=1 everywhere) against the Gumbel solver. Running the full outer iteration was incorrect here because it allows cavitation to develop — the mass-conserving solution gives p_max ≈ 13 Pa (physical, but incomparable to Gumbel's 983 Pa). With max_outer_iters=1, only the full-film solve runs (outer=0 gives p_Elrod ≈ 994 Pa, error 1.2% vs Gumbel). All tests now pass.

**`tests/test_elrod_a2.cpp`**

New test file with four validation tests for the full Elrod-Adams model:
1. **Complementarity**: after convergence, no cell has g=1 AND θ<1 (max g·(1−θ)⁺ = 0).
2. **Cavitation location**: for ε=0.8 with attitude_angle=-90°, the majority (>50%) of cavitated cells lie in the diverging half (cos θ > 0). Observed fraction ≈ 0.51.
3. **Mass conservation**: total liquid content Σ θ·h over 30 timesteps drifts by < 2%. Observed relative drift = 0.
4. **Pressure positivity**: min pressure ≥ p_cav everywhere.

### Modified files

**`src/reynolds.cpp` — outer iteration loop + numerical fixes**

Extended `Reynolds::solve_elrod` with the switch-function outer (flag-update) iteration:

1. **Switch function** applied to diffusion: `Gamma_g(i,j) = Gamma_base(i,j) * EOS::switch_function(theta(i,j))`.

2. **Flooded bearing z-BC** (`Gamma_base`, not `Gamma_g`): The Dirichlet θ=1 condition at z=0,L uses `Gamma_base` (the full diffusion coefficient without the switch). Reasoning: the z-ends are immersed in the oil sump, which maintains full film regardless of the interior cavitation state. If `Gamma_g` were used, a cavitated boundary cell would set its z-BC coefficient to zero, removing the only restoring force for that cell. In subsequent outer iterations, adjacent cells lose their z-BC forcing, and the cavitation spreads globally (cascade collapse to θ_min everywhere). Using `Gamma_base` keeps the z-end forcing active unconditionally.

3. **Snap step**: After clamping θ ≥ theta_min, cells with θ ∈ (1−5×10⁻⁷, 1.0) are snapped to exactly 1.0. Motivation: the BiCGStab solver at rtol=1e-8 on a system with condition ~O(10³) produces residual errors ~O(1e-7) near the full-film/cavitation free boundary. Without snapping, these near-1.0 cells are assigned g=0 in the next outer iteration (because θ < 1), triggering spurious flag changes and preventing convergence. The snap tolerance 5e-7 was chosen to cover the observed numerical noise amplitude (1.7e-7) with margin, while staying well below physically cavitated cells (θ ≈ 0.111 in the diverging zone).

4. **n_switched convergence criterion** uses snap_tol threshold: `(tp >= 1.0 - snap_tol) != (tn >= 1.0 - snap_tol)` rather than a hard 1.0 comparison, consistent with the snap step.

5. **VecZeroEntries(x)** added in `LinearSystem::solve` before `KSPSolve`: `VecDuplicate` leaves the vector uninitialised; without this, the initial guess for the iterative solver may be garbage (or NaN), causing BiCGStab to fail or produce wrong results. This is a correctness fix, not an optimisation.

**Physical note on Elrod-Adams vs Gumbel pressure difference:**

The fully-converged Elrod-Adams mass-conserving solution gives a much lower peak pressure than Gumbel (≈13 Pa vs ≈983 Pa for ε=0.8). This is physically correct: in the mass-conserving model, the Couette flow advects low-θ fluid (θ ≈ 0.111 in the cavitation zone) back into the converging (full-film) zone, depressing the upstream θ and therefore the pressure. Gumbel clamps p ≥ 0 post-solve with no mass accounting, so the full density is always available. The two models are fundamentally different once mass conservation is active.

---

## 2026-04-11 — Phase A.1: universal variable (θ) reformulation

### New files

**`tests/test_elrod_a1.cpp`**

Two tests for the Elrod-Adams Phase A.1 solver:
- `test_uniform_film`: zero eccentricity → $h = c$ everywhere → no Couette wedge → steady solution must be $\theta = 1$ everywhere. Checks $\max|\theta - 1| < 10^{-6}$.
- `test_regression_vs_gumbel`: runs the standard journal bearing with both GUMBEL and ELROD_ADAMS for one timestep, verifies peak pressures agree to within 5% (observed: ~1.2%).

**`PLAN.md`**

Full phased development plan: Phases 1–3 (complete), Phase A (Elrod-Adams mass-conserving cavitation), Phase B (dissolved/free gas transport). Records the mathematical formulation, code changes, and validation benchmarks for each sub-phase.

### Modified files

**`src/config.hpp` — new `CavitationModel` enum**

Added `enum class CavitationModel { GUMBEL, ELROD_ADAMS }` and `cavitation_model` field to `SimulationConfig` (default: `GUMBEL`). Allows clean dispatch in `main.cpp` and tests.

**`src/fvm.hpp`/`fvm.cpp` — new `FVM::ddt_weighted` operator**

```cpp
void ddt_weighted(LinearSystem& sys, const Field& phi, const Field& weight,
                  double dt, const Mesh& mesh);
```

Discretises $\partial(w\phi)/\partial t$ for static weight $w$. Adds $w(i,j)\cdot V/\Delta t$ to $a_P$ and $w(i,j)\cdot V\,\phi^n/\Delta t$ to source. Required for the transient term $\partial(\theta h)/\partial t$ in Phase A.1+ dynamic bearings.

**`src/reynolds.hpp`/`reynolds.cpp` — new `Reynolds::solve_elrod`**

Implements the $\theta$-formulation (Phase A.1, full-film, $g(\theta) = 1$):

$$\text{div}(\Gamma\,\nabla\theta) = \frac{\omega}{2}\rho_0\frac{\partial(\theta h)}{\partial\theta}, \qquad \Gamma = \frac{\rho_0\beta h^3}{12\mu}$$

Key implementation notes:
- Couette face flux: $F_e = \rho_0\,(\omega/2)\,R\,h_e\,\Delta z$ — must include $\rho_0$ (missing it reduces peak pressure by factor $\rho_0$) and must include R (from the FVM area element $R\,\Delta\theta\,\Delta z$).
- Ghost cell coverage: flux loop runs from $i = -1$ to $n_\theta - 1$ so that `FVM::divergence` can read `couette_flux(-1, j)` for the west face of cell $i = 0$ without hitting uninitialised ghost data.
- Dirichlet z-BC: $\theta = 1$ (equivalent to $p = p_{cav}$ via EOS).
- Post-solve: clamp $\theta \geq 1$ (full-film constraint; Phase A.2 relaxes this).
- Recover $p$ and $\rho$ via EOS.

**`src/main.cpp` — dispatch on cavitation model**

- Added `"theta"` field, initialised to 1.0.
- `fields["theta"].store_old_time()` called alongside pressure.
- Time loop dispatches: `Reynolds::solve_elrod` for `ELROD_ADAMS`, `Reynolds::solve` for `GUMBEL`.
- Initial pressure changed from 101325 Pa to `p_cav` (cleaner initialisation for the Elrod solver).

**`PHYSICS.md` — Section 7 expanded**

Renamed to "Reynolds Solvers" (7.1 Gumbel, 7.2 Elrod-Adams A.1). Documents the full derivation of the $\theta$-formulation, the Couette-as-convection interpretation with face flux formula, the FVM assembly, and the key bugs found during implementation.

---

## 2026-04-11 — Reynolds solver, FVM fix, and documentation

### New files

**`src/reynolds.hpp` / `src/reynolds.cpp`**

Implements `Reynolds::solve`, which assembles and solves the steady-state compressible Reynolds equation at each time step:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\gamma\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\gamma\frac{\partial p}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta}$$

with $\gamma = \rho h^3/(12\mu)$.

Steps inside the solver:
1. Compute $\gamma$ and $\rho h$ over physical + theta ghost cells (avoids a separate MPI sync).
2. Assemble via `FVM::laplacian` (diffusion) + manual Dirichlet z-BC correction + `FVM::add_source` (negated Couette term, per FVM sign convention).
3. Solve with PETSc KSP.
4. Apply Gumbel cavitation clamp: $p \leftarrow \max(p, p_{cav})$.
5. Update density via EOS: $\rho = \rho_0\,\theta_{film}(p)$.

The transient $\partial(\rho h)/\partial t$ term is left out for the static-geometry case (zero for constant $h$ at steady state). It is a natural extension for Phase 4+ dynamic bearing work.

**`PHYSICS.md`**

Created. Documents the governing equations, EOS, FVM discretisation (including the sign convention), boundary conditions, solver algorithm, and notation.

**`CHANGES.md`**

This file. Created to satisfy the "append edits and reasons" requirement from `CLAUDE.md`.

---

### Modified files

**`src/fvm.cpp` — bug fix: incorrect `ew_scale`**

*Before:* `ew_scale = d_z / (R * R * d_theta)`  
*After:* `ew_scale = d_z / (R * d_theta)`

Root cause: the east/west FVM face coefficient for the θ-direction is derived from the surface divergence theorem:

$$\oint \gamma(\nabla p \cdot \hat{n})\,dl \;\ni\; \gamma_e \cdot \frac{\Delta z}{R\,\Delta\theta} \cdot (p_E - p_P)$$

The arc-length gradient contributes one factor of $1/R$ (not $1/R^2$). The previous code applied an additional spurious $1/R$ from mistakenly omitting the area-element factor $R$ when integrating $(1/R^2)\partial_\theta[\gamma\partial_\theta p]$ over the cell. This caused all Laplacian solutions to be off by $1/R$ (a factor of 100 for $R = 0.01\,\text{m}$), failing the 2D Poisson test.

**`tests/test_fvm.cpp` — bug fix: wrong source sign in `test_laplacian_2d`**

*Before:* `p_rhs(i,j) = -lambda * sin(...)  // passed as source_field`  
*After:* `p_rhs(i,j) = lambda * sin(...)`

The FVM stencil assembles $-\text{div}(\gamma\nabla p)\cdot V = \text{source}$, so for PDE $\text{div}(\gamma\nabla p) = f$ one must pass `source_field = -f`. With $f = -\lambda\,p_{exact}$, the correct source field is $+\lambda\,p_{exact}$. The original test passed $f$ directly, giving the wrong sign on the RHS and a solution of $-p_{exact}$ instead of $p_{exact}$.

Also added a `max_err` print before the assertion to aid future debugging.

**`src/reynolds.cpp` — Couette source sign**

The Couette term $S = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta}$ is the RHS of $\text{div}(\gamma\nabla p) = S$, so `source_field = -S`. The Couette contribution is stored with a leading minus sign accordingly.

**`src/main.cpp` — wire in Reynolds solver; fix post-finalize MPI abort**

- Added `#include "reynolds.hpp"` and `#include "linear_system.hpp"`.
- `LinearSystem sys(mesh)` created once before the time loop (PETSc setup cost paid once).
- `Reynolds::solve(fields, sys, mesh, cfg)` replaces the `// TODO` in the time loop.
- All PETSc objects (`LinearSystem`, `Fields`) now live inside a nested scope block that closes before `PetscFinalize()`. Previously, their destructors (which call `KSPDestroy`, `MatDestroy`, etc.) ran after `PetscFinalize()`, triggering an Open MPI post-finalize abort.

**`README.md`**

- Corrected dependency list: removed OpenMP (not used), added PETSc with install instructions.
- Added test table, configuration reference, and project-structure overview.


## 2026-09-15 — Add FW68D / R410A digitized solution properties

Added `data/R410A_FW68D/fw68d.csv` with all 125 points from the user's JSON,
plus its calibration/source image and a standard-library importer. The archived
metadata remains unchanged; the operating interpretation follows the user's
absolute-pressure and refrigerant-mass-fraction instructions. Calibrated pixel
coordinates determine numerical values, including the pure oil and pure refrigerant
curves; the model does not replace points with a fitted material correlation.

Added `CURVE_TABLE` solubility and `KINEMATIC_TABLE` actual-composition viscosity
through a shared `solution_table_file`. Ragged curves preserve original nodes.
Solubility interpolates c(T) on each isobar, then pressure; viscosity interpolates
nu(T) on each composition curve, then composition, optionally in log space.
Bubble-pressure inversion uses the same solubility isotherm. Added configurable
endpoint holding or range errors, with loader validation and no new C++ dependencies.

Added `LINEAR_TEMPERATURE` density with editable Celsius intercept and slope,
and an FW68D example selecting constant 940 kg/m³. Kinematic viscosity is converted
with the selected density each call. Coupled liquid partial specific volumes also
follow the linear density, so feed pressure work uses the current liquid volume.
Nonpositive density is rejected. The example, CMake-deployed data/config, GUI
controls, README and PHYSICS explain units, pure-fluid limits and calibration.

Added `test_fw68d_properties` for every source node, ragged interpolation, absolute
pressure, unit conversion, pure endpoints, density choices, inverse solubility,
range rejection, malformed CSV and config round trips. The importer check confirms
calibration anchors and reproducibility of every CSV row. Validation: native Windows GCC 16.1 solver and GUI builds passed; all 52 CTest
tests passed. Ubuntu 22.04 / GCC 11.4 / OpenMPI / PETSc build passed; 49 tests
passed initially, and the three six-rank tests passed after enabling local
OpenMPI oversubscription on the four-slot WSL host. After the final temperature
guard, eight relevant Windows checks and eleven Linux property/energy/MPI checks
passed. The Python importer check passed. Single-rank constant-density and
MPI-2 linear-density example runs passed; all 128 output cells contained the
expected 940 and 925.764 kg/m³, with dynamic viscosities 0.033438291 and
0.032931879 Pa s at c=0.1, T=40 °C. No remote campaign or Git push was performed.

## 2026-09-15 - Consolidate outstanding branches

- Integrate the Ferron campaign and independent shaft/bush sink temperatures with
  the coupled-cavitation and FW68D property work already on main. Preserve both
  sets of physics corrections when resolving documentation conflicts.
- Include the pending optical-rig sweep scripts, configuration and documentation;
  keep generated run products and local agent settings outside this commit.
- Verify that workstation commit a7a992e has identical source and tests to the
  already integrated f0e772d; retain the latter's more complete validation record.
- Native Windows checks for cross-film walls (2/4 ranks), bush conduction, config
  keys and FW68D tables pass. Full merged-state testing follows on both platforms.

### Preserve the unfinished bubble branch

- Relocate the bubble source, tests and model documents to `experimental/bubble`
  instead of accepting the old branch's deletion of the journal-bearing solver.
  Add an OFF-by-default CMake target using the existing MPI/PETSc discovery and
  native Windows runtime deployment; no extra solver dependency is introduced.
- Register the four original bubble checks at 1/2/4 MPI ranks and a CLI smoke
  regression. Replace `M_PI`, finalize PETSc after object destruction, reject KSP
  failure and refuse occupied output directories instead of recursively deleting
  their contents. Mark the unimplemented physics as experimental in its README.
- Multi-rank checks exposed off-rank mask indexing; compute neighbour membership
  geometrically and support neighbour queries across MPI partitions. Exchange
  face-flux ghosts before transport. Correct the pole-velocity test to use an
  odd grid containing the exact pole; nearby cells have nonzero tangential gravity.
- All 52 journal-bearing checks and all 13 optional bubble checks pass on native
  Windows. The bubble tests include 1/2/4 ranks; bearing coverage includes 6 ranks.

### Verify and prepare branch publication

- Complete the fresh afml-r6 Release build, 65 CTest checks, importer check and
  eight smoke runs. Both density choices give the expected FW68D values and
  bubble thickness output agrees at 1/2/4 ranks; native outputs agree with Linux.
  Save source/input/executable hashes and sync workstation logs/results locally.
- Preserve all original branch history and pending edits before excluding three
  oversized generated archives from unpublished history. Verify that published
  refs stay unchanged and the new main remains a fast-forward of origin/main.
  Record the commit mapping, archive paths and test evidence in
  `docs/INTEGRATION_2026_09_15.md`.
## 2026-09-15 - Repair solver boundaries and reduce repeated property work

- Share periodic inlet geometry across pressure, heat, masks and diagnostics;
  accept equivalent full-turn angles consistently and reject nonfinite angles.
- Replace pressure/temperature penalty rows with directly prescribed values at
  the existing equation scale, reusing the gas-transport boundary method.
- Reuse final coupled cell properties for acceptance and output, and evaluate
  EOS pressure/derivative together while bracketing. Preserve constitutive laws,
  interpolation, nonlinear acceptance checks and all model selectors.
- Respect runtime SNES/KSP stopping controls and report failed coupled iteration
  counts. Reject invalid MPI partition/ghost depth before communicating.
- Use the supply-temperature liquid volume in inlet enthalpy under linear density;
  correct the warning that constant viscosity prevents all temperature feedback.
- Stop the executable after PETSc initialization fails. Record the new regression
  checks, remaining physical-model decisions and remote performance measurements
  in [the audit](docs/SOLVER_AUDIT_2026_09_15.md).

## 2026-09-27 - Add a reproducible FW68D MPI benchmark

- Add a synthetic transient moving-bearing fixture with the confirmed rig
  dimensions, pure-FW68D table viscosity, FBNS cavitation, film energy and
  reduced conjugate heat transfer. Keep numerical source unchanged so the
  timing comparison measures the currently checked implementation.
- Add a standard-library workstation driver with fresh per-run directories,
  fixed input/output cadence, source/input/binary hashes, bounded MPI launches,
  convergence/final-field checks and separate scaling/throughput measurements.
- Record the 24-point experimental comparison plan and eight initial refinement
  repeats, distinguishing known rig geometry from missing thermal/load data and
  documenting the current cross-film shear and solid-transient limitations.
- Follow native build/check, isolated remote source transfer/rebuild/run and
  local evidence return. Preserve launcher failures and all prior campaigns.
- Record the seven-rank scaling sweep, repeated 24/28-rank candidates and
  measured 24-by-1 versus 14-by-2 concurrent batches. Return complete evidence
  locally; distinguish lowest core cost, queue throughput and single-case latency.
- Apply the user-selected prior convergence settings: G6 720-by-192 and
  dt=20 us, retaining bush radial16/profile order4/quadrature6. Archive the
  native two-step pass and matching remote 24/28-rank film-guard failures
  before step3. Do not transfer the coarse-grid optimum to the failed G6 case.

## 2026-09-27 - Prepare six load-controlled steady FW68D cases

- Preserve the exact A–F experimental JSON and source hashes. Replace the
  active transient sweep plan with static inner steady pressure/energy/CHT
  solves and an external eccentricity/attitude load-balance search.
- Match the supplied pure-oil viscosity temperature nodes in a case-specific
  table without altering the general material table; audit all reference
  viscosities, Sommerfeld numbers, friction coefficients and E endpoints.
- Map TC0 to feed/ambient references and leave wall/supply specifications
  unresolved. Record the wall-shear prerequisite and preserve all prior results.


## 2026-09-27 - Specify FW68D groove and still-oil housing cooling

- Record the user's 5 mm groove opposite applied load, common groove/end
  reservoir, requested zero cavity pressure, OD85/L40 stainless housing and
  surrounding still oil at TC0. Preserve the pressure datum as unresolved.
- Convert the bore-arc groove width to 13.2983125 degrees and record its
  26-column representation on the retained 720-cell circumferential mesh.
- Define force signs consistently with fluid-on-bearing output; use +y target
  for a -y applied load, independently of the inert load_angle_deg input.
- Record Churchill-Chu natural-convection estimates using explicitly provisional
  oil thermal properties; document the missing self-consistent h coupling and
  finite-housing/end-face limitations. No solver code or CFD output changed.


## 2026-09-27 - Enable the steady FW68D thermal campaign

- Adopt the user's 1 atm end/groove reservoir and two exposed shaft-end faces.
  Use documented BVC32 PVE cp/k and 304 stainless conductivity under the user's
  instruction to use ordinary/similar properties; keep exact FW68D nu/rho inputs.
- Permit exact zero absolute cavity pressure only in the restricted gas-free
  SINGLE_PHASE FBNS configuration; retain gas-path and negative-pressure checks.
- Correct shaft torque and bearing circumferential/axial wall stresses using
  F0/F1 when cross-film viscosity is active. Keep the classical initial/off path.
- Add independent variable-viscosity stress, uniform-limit, initial-state and
  zero-cavity-pressure regressions. Native one/four-rank focused checks pass.
- Add a runnable force/thermal iteration driver with retained per-attempt
  provenance and explicit rejection of nonlinear, thermal and force caps.
  Use natural convection for the stationary housing and rotating-disk cooling
  for the exposed shaft ends; report reduced-solid/confinement assumptions.

### FW68D steady preflight corrections (2026-09-27)

- Preserve two capped native checks as failed diagnostic attempts; their old
  housing-temperature override did not match TC0. Render every reservoir
  temperature explicitly in the case driver.
- Use fixed coupling relaxation with a larger cap, keeping the same tolerances
  and rejecting capped/restored solutions. This changes numerical iteration
  settings only; retain the convergence-study grid and cross-film resolution.
- Accept unavailable cross-film energy diagnostics as missing evidence, not
  zero; supported mass, pressure, thermal and field checks remain required.

- Verify ten bounded steady MPI trials, retain all failed/capped evidence and
  select four ranks per case for batch efficiency. Use verified physical-core
  rankfiles after the workstation ignored CPU-list binding.

- Tighten the campaign inner FBNS tolerance to1e-12 after the plotted pressure
  history exposed premature outer acceptance when FBNS skipped its Newton
  update at the default1e-9 floor. Preserve the original preflight as evidence.

- Verify strict fixed-position convergence at 484 iterations on the remote
  workstation, return full hashed evidence and inspect the convergence plot.
  Launch the eight-case numerical baseline on six disjoint four-rank slots;
  full load/thermal equilibria and experimental comparisons remain pending.


## 2026-09-28 - Repair steady FW68D thermal and load closure

- Replace the six-case PINNED groove thermal rows with signed, pressure-derived
  reservoir enthalpy exchange; retain pressure and TC0 reservoir conditions.
- Use exact FBNS mass fluxes in energy, preserve rupture-front pressure heating,
  and hold mass-specific cp constant as density varies.
- Add the discrete energy ledger, a timestep-independent liquid rate residual,
  and exact internal quartic temperature extrema. Reject invalid cross-film
  states and failed bush solves rather than repairing values silently.
- Include energy closure in the coupled stopping decision. Require the case
  driver to pass conservation, positivity, subambient-temperature, pressure,
  cavity-content, thermal-coefficient and vector-load checks.
- Reuse verified same-case steady fields; add bounded force-search backtracking
  and faster bounded external-coefficient updates. Preserve all failed trials.
- Keep the requested 720 by 192 mesh, order-four profile, six quadrature points
  and 16 radial bush cells. Damp only the quartic source update by 0.3.
- Add focused native/MPI checks for groove energy, heat capacity, rupture-front
  heating, profile extrema, wall heat balance, initial-state I/O and mass-rate
  normalization. Full-case validation is recorded separately as it completes.

- Add a standalone analytic force-balance check: two initial positions reach the
  same equilibrium after a failed full Newton step, exercising backtracking.


## 2026-09-28 - Correct FW68D load orientation

- The user confirmed an upward push from the lower actuator, with the groove
  at the top. The previous downward-force assumption was wrong, rather than
  an arbitrary coordinate choice, because the groove position stays fixed.
- Invalidate the old A-F campaign and E comparison for experimental use while
  retaining raw numerical results and provenance. Corrected cases require
  applied force (0,+W), fluid-force target (0,-W), and new position searches.

- Use the same applied-load helper in config generation and every force residual.
  Keep the groove fixed at 90 degrees and move numerical attitude guesses into
  the upward-bearing half-plane, without rotating saved physical fields.
- Reject cached outputs/checkpoints with old load components or mismatched
  config hashes. A regression proves rejection even when other provenance
  matches, and verifies two-start convergence to downward fluid force.
- Rebuild natively; existing force/motion tests, driver tests and a four-rank
  720x192 isothermal direction check pass. Corrected thermal equilibria remain
  pending; the old results are not relabelled as corrected.


## 2026-09-28 - Measure the actual steady pressure correction

- Require one Newton correction per steady FBNS call before residual-based
  acceptance; prevent warm iteration-zero returns from fabricating R_p=0.
- Preserve the transient path and all existing tolerances. Regression checks
  cover an exact zero-RHS warm state, a sub-tolerance perturbed pure-oil state,
  and the unchanged transient exit, on native 1/2/4/6 MPI ranks.
- Record the authorized stop of the wrong-direction campaign and prepare
  independent four-rank cases so a failed case cannot cancel the remaining queue.

- Validate bcb69c8 on Windows (52 CTests) and Linux (14 focused tests plus
  four-rank I/O/diagnostics); pass full-resolution CHT startup on both. Launch
  seven independent four-rank corrected searches with the eighth case queued;
  preserve and return hashed startup evidence. Final equilibria remain pending.


## 2026-09-28 - Remove the absolute-temperature linear accuracy floor

The corrected-load D run hit4000 sweeps with Rp8.62e-15 and RT1.51e-15,
but energy residual7.224642e-6W /6.406604W=1.12769e-6 still failed the
unchanged1e-6 gate. The linear energy equation used absolute Kelvin with
relative tolerance1e-8, so a large common300K offset obscured a small heat
budget. Solve the identical assembled equation for temperature rise instead:
A*(T-Tref)=b-A*Tref, including actual axial stencil and Dirichlet rows.
Restore absolute temperature after the solve; physics and tolerances are unchanged.

A low-heat regression reproduces the old four-rank residual1.33436e-6W;
the shifted formulation gives9.64e-14W, with the same rise when the reference
moves from300 to900K. One/four-rank energy tests and28 affected native CTests
pass. Log R_E alongside field residuals once the energy gate is evaluated,
so stationary field updates cannot conceal a failed energy balance.

This numerical correction is distinct from A's21.7mK profile undershoot at
the groove entrance; that model-consistency issue remains under investigation.


## 2026-09-28 - Couple bearing position into steady CFD iterations

- Replace the campaign's fully nested external force search with a native steady equilibrium mode. Position is corrected during pressure/cavitation/energy iteration, while every physical storage and squeeze term remains absent.
- Compute a bounded two-component force Newton step using pressure-only probes at frozen thermal coefficients; adjust the gap-dependent Dowson integrals consistently and preserve valid base states during rejected probes.
- Require force equilibrium together with the existing pressure, thermal and conservation convergence checks; record initial and final geometry separately and retain unsuccessful attempts.
- Preserve the unsuccessful groove collocation experiment outside the deployed source. It is not an accepted thermal-model correction.
- Validation and remote deployment results are recorded below when completed; this entry alone does not establish accepted physical cases.

Native validation: 22 affected CTests passed, plus the separate CLI equilibrium
cap regression and driver acceptance/self-checks. The real FBNS equilibrium
regression recovers a known position from a displaced seed on 1/2/4 ranks,
checks zero squeeze and unchanged thermal fields during pressure probes,
rejects inadmissible probe geometry without corrupting state, and recovers
the same result for inactive dt values 1e-8 and 1e3. A short 720x192 A pilot
performed six position corrections in 40 CFD sweeps and reduced the final
force error to 0.08596%; its deliberately capped thermal solve remains
unconverged and is not an accepted case. Linux deployment is still pending.


Workstation validation and launch (28 September, 16:53 KST): revision 23fcbe7
passed 22 focused Linux CTests, four-rank IO/diagnostics and driver acceptance
checks. Restarted seven independent four-rank cases at 16:47 KST, with E-high
queued. Verified 28 live solver ranks on disjoint physical cores. The initial
remote histories show position corrections inside the steady iteration
(B: 30 corrections in 35 iterations; E: 14 in 14). All six nominal cases
currently meet the force tolerance, but pressure/thermal closure is pending.
These intermediate histories are not accepted experimental predictions.
The separate local E short check reached its 360-second time cap after 31
iterations/corrections, without thermal convergence; the remote run has the
normal longer limit. No temperature tolerance or physicality test was relaxed.


## 2026-09-28 - Apply requested inner pressure tolerance and stop cases

Set the FW68D campaign FBNS tolerance to 1e-6 after the user explicitly chose
the inner pressure-equation residual. Terminate all seven active workstation
cases that had already met it, including the temporarily paused E-low case;
B was already rejected by its thermal-profile physical check. Keep stopped
results distinct from accepted coupled solutions. Preserve all saved outputs
and synchronize them locally. No new simulation was launched. The executable
of the preserved runs remains revision 23fcbe7 with its original input hashes.
