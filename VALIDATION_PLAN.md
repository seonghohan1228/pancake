# Validation & Fix Plan (2026-08-02)

Agent-facing plan to resolve the CFD-vs-experiment mismatches and complete the
validation matrix. Produced from three audits: the defect register
(HANDOFF.md D1-D16, all still open), a line-level audit of `src/energy.cpp`,
and an independent re-solve of the Ferron and Brito load cases with an external
Christopherson solver.

## 0. Diagnosis — read before doing anything

**The solver core is not the problem for load.** Load integration
(`reynolds.cpp:668-797`), film thickness (`film_thickness.cpp:27`), clearance
handling, groove pressure (absolute), viscosity wiring and omega were all
audited clean, and an independent isothermal Reynolds solve reproduces
pancake's loads to within a few percent at the same inputs. Raimondi & Boyd
passes (S within 5.8%). Therefore the +52-68% (Ferron) / +65-106% (Brito) load
over-prediction is inherited from **reference data and shared physics
assumptions**, not from the discretization.

**The energy equation is dimensionally correct but structurally incomplete.**
Depth-averaged form, convection velocity, Couette (`mu U^2/h`) and Poiseuille
(`h^3|grad p|^2/12mu`) dissipation are all correct per-area forms. The failures
are:

1. **No conjugate wall model** (`energy.cpp:128-144`): Newton cooling to a
   *prescribed constant* T_wall. No shaft or bush conduction, no
   circumferential heat spreading by the rotating shaft. This is the diagnosed
   cause of the 71 deg temperature-peak location error (predicted X/pi = +0.95
   vs measured +0.55) — with no upstream heat path, T can only climb
   monotonically groove-to-groove. No value of h_wall fixes this.
2. **Groove hard-pins T to t_supply** (`energy.cpp:224-266`, penalty 1e20):
   zero hot-oil carry-over. Real grooves mix 60-80% recirculated hot oil.
   Sized in-repo: +8 K inlet = -22% load. For Brito this reset happens twice
   per revolution.
3. **Cavitated-zone thermal inventory is not theta-weighted** (D14/B3):
   capacity `rho_cp*h`, conduction `k*h`, and convected flux use full `h`, not
   `theta*h` (`energy.cpp:316,318,334,352`); Couette dissipation charges the
   full film under the default `cavitated_film_model = FULL_FILM`
   (`energy.cpp:73-75`). With theta_min ~ 0.39 the liquid mass flux through
   the ~100 deg cavitated arc is over-stated ~2.5x, damping the temperature
   rise that feeds the next inlet.
4. **The uniform-T-across-the-film assumption is invalid on these rigs.**
   The energy equation is depth-averaged with mu(T_mean) and no cross-film
   profile; validity requires Br << 1, reduced Peclet << 1, and equal wall
   temperatures, and all three fail: Br = mu*U^2/(k*dT) ~ 2.6 (Ferron),
   Pe*eps^2 = rho*c*U*h^2/(k*R) ~ 10-30, walls unequal. Dissipation-driven
   core superheat mu*U^2/(8k) is ~3 K (Ferron 2000 rpm) to ~6 K (Brito
   4000 rpm), scaling as U^2, worth ~10-20% of load at high speed via the
   ~3.4%/K viscosity sensitivity. This is also why the constant-h_wall
   closure cannot be fixed independently: wall heat flux and effective
   viscosity both come from the same missing T(z) profile (see task 2.4).
5. Net effect of 1-4: the film in the loaded arc is **4-6 K colder than
   Ferron's measured bush wall** — physically impossible under adiabatic
   walls — inflating viscosity and contributing an estimated 10-25% to the
   load over-prediction (more at 4000 rpm where cross-film effects peak).
   The remaining ~30-50% points at reference data (Phase 1) and is partly
   shared with Alakhramsing's own FE model (+20-30% on the same rig).
   Discriminator: after Phase 1 settles the data, any residual error that
   grows as U^2 is cross-film physics; a speed-independent residual is
   geometry/property data.

**Config traps found (fix first, they poison every experiment):**

- **B1 — silent THD decoupling.** `fluid_property_model = CONSTANT` or
  `liquid_viscosity_model = CONSTANT` freezes mu at `cfg.mu`
  (`fluid_properties.cpp:148`) while `temperature_model = ENERGY_EQUATION`
  still runs and writes plausible-looking T fields. The shipped `config.txt`
  (lines 45/57/65) is in exactly this state. No warning is emitted.
- **B2 — steady-state scheme downgrade.** LINEAR/VANLEER deferred correction
  is built from `phi.old` (`fvm.cpp:147-193`), refreshed only at end of
  timestep (`main.cpp:405`). STEADY_STATE runs one timestep, so every outer
  iteration sees the uniform initial field: all schemes silently degrade to
  first-order upwind. Any steady-state scheme-sensitivity study is currently
  meaningless.
- **D6 — KSP convergence never checked** (`linear_system.cpp:119-122`);
  divergence/NaN propagates and the run exits 0.
- **Free-boundary resolution is governed by `p_max/beta`, not residuals**:
  keep `p_max/beta >= 4e-3` or the rupture angle is wrong while mc, ns, R_p
  all read zero (established by `run_raimondi_boyd.convergence_study`).
- `R_p` is normalized by `property_reference_pressure = 1e5`, so `R_p = 1`
  means 1 bar of per-iteration pressure movement. Brito 3000 rpm / 70 kPa has
  R_p = 3.56 — genuinely unconverged, do not quote it.

---

## Phase 0 — Guardrails — **DONE 2026-08-02** (23/23 tests pass, no physics changed)

| # | Task | Outcome |
|---|------|---------|
| 0.1 | KSP convergence + finiteness guard (defect D6) | **Done.** `LinearSystem::solve` takes an equation label, checks `KSPGetConvergedReason`, scans for non-finite values under `MPI_Allreduce`, aborts with a named message. All five call sites labelled. New `tests/test_ksp_guard.cpp` asserts on the guard's *own message* (`PASS_REGULAR_EXPRESSION`), not on a nonzero exit, so a PETSc-internal failure cannot pass it; verified firing with `DIVERGED_NANORINF`, exit 1. No existing test relied on a diverged solve |
| 0.2 | Warn on the silent THD-decoupling trap (B1) | **Done.** Warning fires on the shipped `config.txt`, confirming the trap was live. The example's physics was deliberately **not** changed (enabling coupling needs the `EMPIRICAL` + `solution_viscosity_gas_coeff` workaround, an unrequested behaviour change); it now carries an in-place comment. Also added: a warning when both wall heat-transfer coefficients are zero (adiabatic film) |
| 0.3 | Validate `t_supply` on fed inlets | **Done.** Non-positive `t_supply` with `feeds_fresh_oil` is now a hard error |
| 0.4 | Record iteration counts in `run_brito.py` | **Done.** `elrod_iters`, `coupling_iters`, `coupling_capped` now in the JSON; `coupling_max_iters` sourced from one constant instead of hardcoded in the config template |
| 0.5 | Fix deferred-correction staleness (B2) | **Done, narrower than planned.** `solve_elrod` was already correct (it copies the current iterate into `theta_lag`); the bug was only in `energy.cpp:365` and `gas_transport.cpp:281`, which passed the live field whose `.old` is the previous timestep. Both now use the same lagged-copy idiom. Ghost currency verified: both fields are exchanged immediately before their solves (`main.cpp:329`, `:314`) |

Do **not** "fix" `ddt_weighted`'s new-time capacity (B7 in the energy audit) —
it is consistent with the non-conservative `F·grad T` convection form used
deliberately; changing one without the other breaks the pair.

## Phase 1 — Reference-data verification (no code changes; can run parallel to Phase 0)

The single highest-value item: **one scalar of ~0.55 on `mu*omega*R^2/c^2`
explains the entire Ferron load error**, and the only unverified inputs in
that group are L, c, and the speed of the eps-vs-load figure.

| # | Task | Acceptance |
|---|------|-----------|
| 1.1 | Obtain primary Ferron/Frene/Boncompain 1983 values for **L, radial clearance c, and the speed of the eccentricity-vs-load data** (currently digitised from Alakhramsing Fig. 10a, whose panel does not label a speed — `ferron.py:82-87`). Test the 1000 rpm hypothesis: an independent solve at 1000 rpm gives 1604/3798/7069 N at eps = 0.47/0.67/0.78 vs measured 2000/4000/8000 — a far better fit than 2000 rpm | Speed and geometry pinned to the primary source; load table re-run; error either collapses or is confirmed real |
| 1.2 | Brito Fig. 6: the polar re-digitisation of eccentricity (`brito.py:105-137`) moved eps up ~0.15 and drives the whole +65-106%. Re-digitise independently; also run the **Cartesian reading as a documented bracket** (`solve_eps_for_load` at `run_brito.py:261-289` is dead code — wire it up to report eps-for-10kN as a third line) | Both readings + eps-for-load reported side by side with uncertainties |
| 1.3 | Brito hot running clearance: Cd = 171 um is quoted at 20 C; bulk rise ~20 K with bronze bush on steel shaft grows c by ~15-20% (W ~ 1/c^2, ~30% of load). Add a clearance-at-temperature sensitivity line to run_brito | Load bracket [c_20C, c_hot] reported per case |
| 1.4 | Bridge the validated regime to the benchmark regime: add an R&B-style fed-film regression at **L/D = 0.8, c/R ~ 2-3e-3** (Ferron/Brito ratios) against the independent Christopherson solver used in the audit | Agreement < 5% in W at fixed eps across eps = 0.4-0.8 |
| 1.5 | Re-run or drop Brito 3000 rpm / 70 kPa (R_p = 3.56) and mark 4000/70 (R_p = 0.42) marginal; never quote unconverged rows | brito_results.json has no rows with R_p > 1e-2 quoted in tables |

Note for 1.2/1.3: the 4000 rpm Brito cases cannot be rescued by eccentricity
at all (isothermal load at eps = 0.40 is still 14.8 kN vs 10 kN applied); they
require mu_eff ~ 0.0105 Pa.s (~71 C mean film). Those cases are the cleanest
evidence that the film temperature is under-predicted — keep them as the
sentinel metric for Phase 2.

## Phase 2 — Energy-equation fixes (sequential; size each with a before/after Ferron + Brito rerun)

| # | Task | Files | Sizing / acceptance |
|---|------|-------|---------------------|
| 2.1 | Theta-weight the cavitated thermal inventory: capacity `rho_cp*theta*h`, conduction `k*theta*h`, convected flux `rho_cp*theta*h*u` (D14/B3). Mirror the Reynolds-side weighting (`fluid_properties.cpp:258`) | `src/energy.cpp:316,318,334,352` | Rerun Ferron adiabatic: report delta T_max, delta T-profile through the cavitated arc |
| 2.2 | Scale Couette dissipation by theta in cavitated cells; either make STRIATED the default under JFO or add the theta factor to FULL_FILM (B4) | `src/energy.cpp:73-75`, `src/config.hpp:620` | Ferron STRIATED-vs-FULL_FILM delta documented (do 2.1 and 2.2 together; their errors partially cancel — see energy audit section 1) |
| 2.3 | Hot-oil carry-over at grooves: replace the hard pin with `T_groove = lambda*T_supply + (1-lambda)*T_upstream_film`, lambda from a mixing model (flow-weighted: supply flow vs upstream Couette flux is the standard closure) with a config override | `src/energy.cpp:224-266` | Ferron groove-exit T ~ 41.9 C (Alakhramsing) reproduced within ~0.5 K; load moves toward measurement; peak location expected NOT to move (in-repo sweep showed location insensitive to inlet T) |
| 2.4 | **Cross-film temperature profile + wall coupling** — the fix for the 71 deg peak-location error and the high-speed load residual. The film's uniform-T assumption fails all three validity criteria on these rigs (Br ~ 2.6, reduced Peclet ~ 10-30, unequal walls; core superheat ~3 K Ferron / ~6 K Brito 4000 rpm, scaling as U^2). Do it Lund-Tonnesen style: assume a low-order polynomial T(z) per cell (4th order), collocated with (a) the depth-averaged energy balance, (b) wall BCs from a minimal conjugate model — shaft circumferentially averaged (solved from the global film-shaft heat balance), 1D radial bush conduction to an outer convection BC. The same T(z) then supplies BOTH the wall heat flux (replacing the constant h_wall closure) and the viscosity z-integrals for a **generalized Reynolds equation** (Dowson F0/F1/F2) — do not fix one without the other; they come from the same missing profile. Full 2D bush conduction stays Phase D scope | new module + `src/energy.cpp:128-144`, `src/reynolds.cpp` (gamma assembly) | Ferron: T_max within 2 K **and** peak location within ~15 deg with one wall setting (no more bracketing); Brito: shaft interface T vs supply-T sweep qualitatively reproduced; Brito 4000 rpm sentinel load closes while 1000 rpm is untouched (cross-film effect must scale as U^2 — if the 1000 rpm gap moves equally, something else changed) |
| 2.5 | Sensitivity: enable `consistent_boundary_flux = true` (default false forces theta -> 1 on axial-boundary rows of the cavitated zone — artificial end re-flooding, `src/reynolds.cpp:490-494`) | config only | Cavity shape/side-flow delta documented |

Order matters: 2.1+2.2 first (cheap, local), then 2.3, then 2.4. Rerun the
Brito 4000 rpm sentinel after each step; success = predicted mean film
temperature rising toward the ~71 C the load data implies.

## Phase 3 — Grid and time-step independence (after Phase 2, before final validation)

Prereq: 0.1 (trustworthy convergence), 0.5 (schemes actually differ).

1. **Grid:** extend `validation/grid_convergence/` from 3 to 4 levels
   (60x20 → 480x160), energy equation **on**, JFO on, at Ferron geometry.
   Track W, p_max, rupture angle, T_max, T-peak location, cavitated fraction.
   Compute observed order + GCI (Richardson) per quantity. Constraint: scale
   beta per level to keep `p_max/beta >= 4e-3`, otherwise the rupture angle
   error masquerades as grid error. Acceptance: monotone convergence, observed
   order ~1-2 (upwind-dominated), GCI < 2% on W and T_max at the working grid.
2. **Time step:** transient case (start-up ramp `omega_ramp_time > 0`, or a
   step load under `motion_model` dynamics). dt sweep x8 range at fixed grid.
   BDF1 everywhere (D12): acceptance = observed temporal order ~1 and the
   working dt within 1% of the Richardson limit on eps(t) and T_max(t).
3. **Scheme sensitivity:** UPWIND vs LINEAR vs VANLEER on theta and thermal
   convection, only after 0.5. Note D1: van Leer far-upwind cell is wrong for
   F<0 and there is no z-face TVD branch — fix or restrict the study to
   positive-flux configurations.

## Phase 4 — Validation matrix (final; each row = one agent task with a pass/fail verdict written into validation/README.md)

### 4a. JFO / mass-conserving cavitation

| Case | Observable | Criterion |
|------|-----------|----------|
| V1: 1D JFO slider/pocket (analytic) | rupture + reformation position | exact to O(dx) |
| Raimondi & Boyd re-run (regression) | S, attitude | < 10%, < 5 deg (currently 5.8%, 4.1) |
| Beta/free-boundary study (exists) | rupture angle vs Swift-Stieber | < 1 deg at p_max/beta >= 4e-3 |
| Engels et al. 2024 (Lubricants 13(1), open access) | **measured void fraction** — the direct Elrod observable; nothing currently tests it | digitise + compare theta field; criterion TBD from paper's uncertainty |
| Vijayaraghavan & Keith grooved case (V3) | cavity boundaries with grooves | qualitative match |

### 4b. Energy equation / THD

| Case | Observable | Criterion |
|------|-----------|----------|
| Ferron 1983 (rerun after Phase 2) | p_max, peak location, T_max, **T-peak location**, load-at-eps, Q | p: keep current pass; T_max < 2 K; T-peak < 15 deg; load: re-judge after Phase 1 |
| Brito 2006 full matrix | p profile, p_max vs (N, P_supply), eccentricity locus, flow rate, **shaft interface T** vs T_supply/P_supply sweeps | p_max within digitisation + clearance bracket; interface T trends correct; sentinel: 4000 rpm cases |
| Ferron/Brito cross-check | same rig lineage — inconsistency between the two verdicts indicates a data problem, not a solver problem | verdicts agree in direction |

### 4c. Gaseous cavitation / dissolved gas

Blocked on a numerics change, not a bug fix: the VOID_COUPLED outer Picard
loop (`equilibrate_gas` <-> `solve_elrod`) does not converge (R_p 0.57-7.2;
answers relaxation- and grid-dependent). Order:

1. Reformulate the coupling: fold the flash/plateau `p_void = L*theta*p_sat /
   ((1-theta) + L*theta)` into the Elrod switch (single nonlinear solve) or
   Newton on the coupled residual. Fix D3 (initial theta from
   `effective_p_cav()`, `main.cpp:179`) and D2 (hardcoded `theta >= 1.0` vs
   per-cell ceiling) on the way; D5 (stale rho_cp/k ghosts between gas and
   energy solves) before any MPI-parallel validation run.
2. Convergence gate: R_p < 1e-6, grid- and relaxation-independent to < 1%.
3. Validate: Grando et al. refrigerant-oil equilibrium (V5); then the
   staged targets in docs/CAVITATION_OIL_REFRIGERANT.md (PZ68S + R290,
   p_cav -> local bubble point = D4/F9); Engels 2024 void fraction again if
   the rig has dissolved-gas character.

### Standing rules for all validation tasks

- Never quote a run with R_p > 1e-2 (in units of 1e5 Pa), mc > 1e-10, ns > 0,
  or the coupling-iteration cap hit.
- Always quote temperature results with the wall condition attached.
- Keep `property_reference_temperature` at the temperature where mu is quoted;
  moving it with supply temperature silently cancels the physics (documented
  trap, validation/journal_bearing/README.md).
- Reference data is audited per number: record digitisation uncertainty next
  to every digitised value (brito.py is the template).
- After each phase: update PLAN.md status rows, CHANGES.md, and this file.
