# HANDOFF — maintainer's guide for future contributors

Date: 2026-07-07. Written after a full audit pass over the solver core
(`reynolds.cpp`, `fvm.cpp`, `energy.cpp`, `gas_transport.cpp`, `main.cpp`,
`config.hpp`, `equation_of_state.hpp`). Every defect listed in §4 was verified
by reading the current code at the cited line, not inherited from older reviews.

**Who this is for:** anyone (human or model) making edits to this codebase
without the context of how it grew. Read §1–§3 before touching solver code.
§4 is the verified punch list — pick items from the top. §5 is the longer-term
backlog. §7 tells you how to prove you didn't break anything.

Companion documents:

- `PHYSICS.md` — governing equations and model derivations. Keep it in sync.
- `PLAN.md` — work-package history (WP-1…WP-12) and open items (WP-7/8/9/10).
- `AUDIT.md`, `REVIEW.md`, `REVIEW_VERDICT.md` — the 2026-06 literature-grounded
  audit that drove WP-1…WP-12. Historical; most findings are now closed.
- `docs/CAVITATION_OIL_REFRIGERANT.md` — why the R290/PZ68 case is *gaseous*
  (outgassing) cavitation and why the onset must track `p_sat(T, c_d)`.
- `BUILDING.md`, `docs/WINDOWS_PORT.md` — build contract (Linux + MSYS2 MinGW64).

---

## 1. What this code is, in one page

2-D compressible Reynolds (thin-film) solver on a cylindrical `theta × z` mesh,
C++20 + MPI (1-D decomposition in theta, periodic wrap), PETSc BiCGStab +
block-Jacobi for the linear solves, VTK/PVD output. A Win32/ImGui GUI
(`src/gui/`) orchestrates the solver as a subprocess through `config.txt`.

One transient timestep (`main.cpp:248-411`) does, in order:

1. `update_film_geometry` — `h`, `dh_dt` from static eccentricity/tilt or the
   moving-bearing state.
2. `FluidProperties::update_solution_fields` — liquid-solution `rho`, `mu`,
   `cp`, `k` from `(p, T, c_d)`; two-phase `mu` closure if free gas exists.
3. `Reynolds::solve_elrod` (JFO) or `Reynolds::solve` (Gumbel/Sommerfeld) —
   pressure/film content. Elrod runs its own inner active-set loop over the
   switch function `g(theta)`.
4. `GasTransport::solve` — advect `dissolved_gas` (material form) and
   `free_gas_mass` (conservative form) with the film-averaged velocity; then
   `FluidProperties::update_gas_state` — finite-rate release/resorption between
   dissolved and free gas (operator splitting).
5. `Reynolds::calculate_velocities` / `calculate_macroscopic_properties` —
   staggered face velocities, bearing forces, friction torque.
6. `Energy::solve` — film-averaged temperature, non-conservative (advective)
   form, with wall heat loss and axial boundary convection.
7. Diagnostics (mass balances to `diagnostics.csv`), VTK output, `store_old_time`
   on the six state fields, bearing ODE advance (lagged fluid force).

`STEADY_STATE` mode wraps steps 2–6 in a Picard loop
(`coupling_max_iters`, `coupling_tolerance`, `coupling_relaxation`) and runs
exactly one "timestep". TRANSIENT runs the sweep once per step — there is
currently **no** within-step outer iteration in transient mode (deliberate,
see WP-2 notes in `PLAN.md`).

State variables and their meaning:

| Field | Meaning |
|---|---|
| `theta` | Elrod universal variable: fractional film content. Full film when `theta >= theta_full` (see below), cavitated below. |
| `pressure` | Recovered from `theta` through the (generalized) barotropic EOS. Never independent of `theta` under JFO. |
| `rho` | Always `rho_liquid_solution * theta`. This is the *liquid* mass density per gap volume, not the two-phase mixture density. |
| `dissolved_gas` | Mass fraction `c_d` of gas dissolved in the liquid. |
| `free_gas_mass` | Released gas mass per unit bearing area (`kg/m^2`). |
| `alpha_gas` | Void fraction of free gas, derived in `fluid_properties.cpp`. |

The WP-1 generalized EOS (`equation_of_state.hpp:19-27`): full-film ceiling
`theta_full = 1 - alpha_g`, cavitated plateau `p_void = max(p_cav_eff,
m_g R_g T / ((1-theta) h))`, mixture bulk modulus `beta_bar`. Under the default
`gas_pressure_coupling = NONE` these collapse to `(1, p_cav_eff, beta)` and the
solve is bit-identical to the single-phase model. **This ceiling is the source
of the most important open bug family — see D2 in §4.**

---

## 2. Invariants — do not break these

These are the contracts that keep the solver correct. They are easy to violate
silently.

**I1 — Ghost discipline.** Any field read at `i-1`, `i+1`, `i-2`, `i+2` inside
an assembly loop crosses MPI partition boundaries. `n_ghost = 2` exists because
the TVD/deferred-correction stencils need `i±2`. If you add a field that a
solver reads at neighbor offsets, you must `comm.update_ghosts(...)` it after
every write and *before* the read. The existing sync points are hand-maintained
lists in `main.cpp:291-323` — fragile by design; see D5. The regression that
proves this matters is `test_elrod_a2_mpi6` (a visible contour artifact
appeared only at 6 ranks).

**I2 — `old()` semantics.** `Field::old_data` is *previous timestep* for
`ddt`/`ddt_weighted`. The deferred-correction schemes in `FVM::divergence`
also read `phi.old(...)` as the *lagged iterate*. `solve_elrod` resolves this
collision with a separate `theta_lag` field whose `old` is refreshed per outer
iteration (`reynolds.cpp:373-374`) while `theta.old` stays at the previous
timestep for the capacity term. If you add a deferred-corrected equation, copy
that pattern; never call `store_old_time()` on a state field mid-step.

**I3 — The EOS triangle.** `theta`, `pressure`, `rho` must stay mutually
consistent through `EOS::pressure_from_theta` / `theta_from_pressure` with the
*same* `(theta_full, p_void, beta_bar)` triple everywhere. Anything that writes
one of the three must reconcile the others. Violations of this triangle are
exactly bugs D3 and D8 below.

**I4 — One boundary-flux formula.** `SimulationConfig::elrod_boundary_outflow`
and `reformation_boundary_velocity` (`config.hpp:1248-1281`) are the single
source of truth for what crosses an axial boundary. `solve_elrod`, the energy
boundary convection, the gas boundary convection, and the diagnostics mass
balance all call them. If you change boundary physics, change it *there* and
nowhere else, or the mass-balance residual `r_l` in `diagnostics.csv` stops
closing (that residual sitting at ~`linear_rtol` is your proof of consistency;
WP-11 measured 7e-8 with `consistent_boundary_flux` on vs 1.4e-3 off).

**I5 — Clamps are accounted, not hidden.** Every clamp/snap/overwrite of a
conserved quantity must accumulate its implied mass into
`cavitation_clamp_mass` / `gas_clamp_mass` so `Diagnostics::mass_balance` can
close. See `reynolds.cpp:468-483` and `gas_transport.cpp:290-301` for the
pattern.

**I6 — Bit-exactness gates.** The acceptance style of this repo: a new model
must reduce *exactly* to the old one when disabled (`NONE ≡ master`,
`FULL_FILM ≡ master`, `alpha_g = 0 ⇒ beta_bar = beta`). Preserve this — write
the equivalence test before writing the feature. `test_gas_coupling` and
`test_coupling` show the pattern.

**I7 — Config surface is triplicated** (parser if-else chain
`config.hpp:696-882`, `to_config_text()` serializer, GUI `schema.cpp` specs).
Until this is unified (B1 in §5), any new key must be added in all three
places *plus* validation *plus* `PHYSICS.md`/GUI hover text. Two shipped bugs
(`gas_pressure_coupling` and `*_table_file` silently dropped on GUI save) came
from missing the serializer. A `parse(to_config_text(cfg)) == cfg` round-trip
test would have caught both — write it (B1a) before adding more keys.

**I8 — Cross-platform contract.** `std::filesystem::path` for paths, no
POSIX-only headers in core code, GCC/Clang/MSVC-buildable. See `ASSESSMENT.md`.

**I9 — Docs move with code.** Repo convention (CLAUDE.md): CHANGES.md gets an
entry with the *reason*, PHYSICS.md gets the mechanism/equations, PLAN.md
tracks status (never delete completed items).

---

## 3. Numerical scheme summary (what is actually implemented)

- **Diffusion** (`FVM::laplacian`): 5-point FVM, harmonic-mean face
  coefficients — correct choice for `h^3` contrast and for the Elrod `Gamma_g`
  which is exactly zero in cavitated cells (harmonic mean blocks the face).
- **Convection** (`FVM::divergence`): implicit first-order upwind everywhere,
  plus optional *deferred* corrections: `LINEAR` (central), `VANLEER` (TVD),
  `TYPE_DIFFERENCING` (central in full film, upwind across cavitation —
  Vijayaraghavan & Keith 1989). Deferred parts lag one outer
  iteration/timestep. **The TVD correction exists only for theta-direction
  faces and its reversed-flow branch is wrong — see D1.**
- **Time**: backward Euler (BDF1) for *all* field equations, via
  `ddt`/`ddt_weighted`. The Elrod capacity term `∂(ρ₀θh)/∂t` assembles for
  every transient run regardless of motion model (WP-3 fix, `reynolds.cpp:393`).
  The squeeze source `ρ₀ θⁿ ∂h/∂t` is explicit/lagged. Bearing rigid-body ODE
  supports up to RK4 — pointless accuracy next to the BDF1 film, see D12.
- **Elrod inner loop**: flag-update active set with flooded first iteration
  (`g=1` everywhere at `outer=0`), `theta_min` clamp, snap-to-full-film within
  `snap_tol = 5e-7` (hardcoded, derived from `rtol=1e-8` — see D10),
  convergence = `max|Δθ| < outer_tol` *and* zero flag flips.
- **Linear solver**: PETSc BiCGStab + block Jacobi, `rtol = linear_rtol`
  (config), `max_iters = 1000`, fresh zero initial guess every solve, **no
  convergence check** — see D6/D7.
- **Inlets**: penalty rows `a_p = 1e20` (`reynolds.cpp:152`,
  `energy.cpp:259`) — see D11.
- **Energy**: non-conservative (advective) form — conservative `div(F T)`
  assembled, then `T·div(F)` subtracted (`energy.cpp:268-291`). This kills
  spurious sources at supply cells but forfeits discrete energy conservation;
  the same trick is used for `dissolved_gas`. `free_gas_mass` is conservative.
- **Gas exchange**: operator-split, first-order kinetics
  `Δc = k_dg (c_d − c_eq) dt` (`fluid_properties.cpp:252-315`); instantaneous
  flash (`equilibrate_gas`) in steady EQUILIBRIUM mode. Hao & Gu (2014)
  precedent — the *kinetics form* is fine; the rate constant for R290 is an
  uncalibrated tuning parameter (see B8).

---

## 4. Verified defects and inconsistencies (this audit, 2026-07-07)

Ordered by (severity × cheapness). Each was confirmed against the code at the
cited lines. "Latent" = wrong code that current shipped configs don't hit.

### D1 — HIGH/latent: TVD (`VANLEER`) is wrong for reversed flow, and z-faces have no TVD at all

`fvm.cpp:171-193`. For a face with `F < 0` the far-upwind cell is `i+2`
(east face) — the code inverts the ratio `r` but keeps the increment
`(phi_up − phi.old(i-1,j))`, which for reversed flow spans the *downwind* side.
The identity `ψ(r) = r·ψ(1/r)` justifies the r-inversion only if the increment
is switched to `(phi_up − phi.old(i+2,j))` as well. Compare
`couette_face_value` (`reynolds.cpp:94-111`), where the caller orders
`phi_uu/phi_u/phi_d` by flow direction and gets it right — the two paths
disagree, so a `VANLEER` Gumbel run and a `VANLEER` Elrod run at `omega < 0`
diverge from each other. Currently latent because the Couette flux sign follows
`omega`, and every shipped config has `omega > 0`; but z-direction fluxes in
energy/gas transport change sign cell-to-cell — which is presumably why nobody
noticed that **the TVD branch simply has no north/south face correction**
(the loop body only handles east/west): axial thermal/gas advection silently
runs first-order upwind even with `VANLEER` selected.

Fix: mirror the direction-aware stencil of `couette_face_value` inside
`FVM::divergence` (helper lambda per face), add the z-face branch. Test: 1-D
step advection with `F < 0` must match the `F > 0` result mirrored; axial
advection of a temperature step must sharpen when switching UPWIND → VANLEER.

### D2 — HIGH/active-with-WP-1: hardcoded `theta >= 1.0` full-film tests contradict the VOID_COUPLED ceiling `theta_full = 1 − alpha_g`

The solver's switch function and snap use the per-cell `theta_full`
(`reynolds.cpp:361, 473-474`), but everything downstream still tests `1.0`:

- `calculate_velocities` Poiseuille gates: `reynolds.cpp:557, 582, 593, 601`
- `calculate_macroscopic_properties` shear gate: `reynolds.cpp:670`
- Elrod convergence flag counting: `reynolds.cpp:496` (`1.0 - snap_tol`)
- `Energy::full_film_switch` / face switch: `energy.cpp:11-17` (heat generation
  gates, thermal flux gates, boundary fluxes)
- `GasTransport::full_film_switch`: `gas_transport.cpp:23-25` and every face
  gate at `gas_transport.cpp:244, 263`
- `SimulationConfig::is_reformation_boundary`: `config.hpp:1266` (`theta_cell < 1.0`)

Consequence: under `gas_pressure_coupling = VOID_COUPLED`, a pressurized cell
with `theta ∈ [theta_full, 1)` is full-film to the Reynolds solve but
*cavitated* to the velocity reconstruction, force integration, dissipation,
thermal/gas face gating, and boundary classification — Poiseuille transport
and shear vanish exactly where the new model says pressure gradients exist.
`n_switched` also miscounts, so the Elrod loop can report non-convergence in
void-coupled runs.

Fix: one shared helper (e.g. `Reynolds::full_film(fields, cfg, i, j)` or a
precomputed `g` field exported from `solve_elrod`) used by all consumers.
Test: VOID_COUPLED case with nonzero `alpha_gas`; assert velocities/forces
change continuously as `alpha_gas → 0` (they currently jump).

### D3 — HIGH/active: initial `theta` uses raw `p_cav`, not `effective_p_cav()`

`main.cpp:179`: `theta_initial = theta_from_pressure(p_init, cfg.p_cav, beta)`
while every solve and recovery uses `cfg.effective_p_cav()`. Under
`cavitation_threshold = LOCAL_PSAT` (which the shipped
`config_r290_pz68.txt` now sets) the initial state is inconsistent: recovering
pressure from that `theta` through the solver's EOS gives
`p_sat + (p_init − p_cav)` — an overshoot of `p_sat − p_cav ≈ 0.97 MPa` in the
propane case. The Dirichlet boundaries pull it back after the first solve, but
the first transient step carries spurious stored mass through
`theta.old` in the capacity term, and step-0 diagnostics/output are wrong.

Fix: use `cfg.effective_p_cav()` at `main.cpp:179` (and audit
`test_elrod_bc.cpp` initial-pressure tests for the same substitution).
One-line fix; add a step-0 EOS-consistency assertion to `test_coupling`.

### D4 — MEDIUM: `LOCAL_PSAT` is a uniform feed-state constant; the plan and the physics call for a per-cell field

`config.hpp:1215-1219` computes `max(p_cav, p_sat(T_initial, c_d_initial))`
once. PLAN.md's own WP-12 spec (step 6) describes a ghost-exchanged
`p_cav_eff(i,j)` field from `max(p_cav, p_sat(c_d_ij, T_ij))` threaded through
`reynolds.cpp` — that part was not implemented. Physically the bubble point
varies strongly over the film: `T` rises tens of K around the circumference
and `c_d` is depleted downstream of rupture, so onset is currently wrong
wherever the local state departs from the feed state. All the ingredients
exist (`saturation_pressure(c_d, T)`, per-cell `VoidParams`, ghost exchange).
Note the WP-1 reconciliation rule when you do this:
`p_void = max(p_cav_eff_local, m_g R_g T/((1−θ)h))` — the local field is the
gas-free floor of the plateau, not a competing threshold. Keep the *force
gauge datum* scalar (see D9).

### D5 — MEDIUM/MPI: stale thermal-property ghosts between the gas solve and the energy solve

Within one coupling sweep: gas transport changes `c_d`
(`main.cpp:310-311`), `update_solution_fields` (`main.cpp:313`) then updates
`rho_cp_liquid_solution` / `k_liquid_solution` on *physical* cells, but only
`rho` and `mu` ghosts are re-synced (`main.cpp:314-315`) before
`Energy::solve` reads `rho_cp`/`k` at `i±1` face averages
(`energy.cpp:333-335, 351-353`). At theta-partition interfaces the energy
equation therefore uses pre-gas-transport capacities — a one-sweep lag that
makes results weakly rank-count-dependent. Same audit point: after
`update_gas_state`, `alpha_gas`/`free_gas_mass`/`dissolved_gas` ghosts are only
refreshed at the *next* sweep's blanket `comm.update_ghosts(fields)`
(`main.cpp:281`), which happens to be adequate for the lagged `void_params`
reads but is an accident of ordering, not a designed contract.

Fix: after `main.cpp:313`, sync the full property set (or bite the bullet and
give `Fields` a dirty-flag so `update_ghosts(fields)` is cheap and always
correct). Test: 1-rank vs 2-rank transient TWO_PHASE run, assert field-wise
agreement to ~`linear_rtol` (the existing MPI test covers theta only).

### D6 — HIGH (cheap): `KSPSolve` result is never checked

`linear_system.cpp:119-122`. If BiCGStab stalls/diverges, garbage enters
`theta`/`pressure`/`temperature` and the run continues, writes VTK, and
"finishes successfully". Add `KSPGetConvergedReason` + iteration count; on
divergence, abort with the reason (or, once adaptive dt exists, cut the step).
While there: the time loop has no finite-ness guard on solved fields — one
`std::isfinite` sweep per step (cheap vs a solve) turns silent NaN plumes into
an immediate, located failure.

### D7 — MEDIUM (cheap, big win): zero initial guess every solve

`linear_system.cpp:121` zeroes `x` before every `KSPSolve`. Inside the Elrod
active-set loop the previous outer iterate is an excellent guess; across
timesteps likewise. Seed `x` from the current field +
`KSPSetInitialGuessNonzero(ksp, PETSC_TRUE)`. Expect a severalfold cut in
Krylov iterations for transient runs. Also expose `KSPSetFromOptions` so PC
experiments (GAMG/hypre) need no recompile.

### D8 — MEDIUM: Picard relaxation acts on `p` and `T`, but the Elrod state is `theta`

`main.cpp:330-331` relaxes `pressure`/`temperature` after the sweep. Under JFO
the next Reynolds solve reassembles from `theta` (unrelaxed); relaxed `p`
re-enters only via `void_params` reads, INLET_OUTLET inflow decisions, and
type-differencing masks — and it breaks the EOS triangle (I3) because `p` no
longer equals `EOS(theta)`. So `coupling_relaxation < 1` under-damps exactly
the variable that matters and de-synchronizes the rest. Fix: relax `theta`
(and optionally `c_d`), then *recover* `p` through the EOS; or relax the
property fields instead. Test: WP-2's own acceptance case ("ω = 0.5 and 1.0
give the same converged answer") — extend it to a JFO+gas case where it
currently doesn't hold tightly.

### D9 — LOW/subtle: force gauge datum moved to `effective_p_cav()`, and the "zero net force" claim is only true for x/y

`reynolds.cpp:640` uses `p − effective_p_cav()`; PLAN.md WP-12 step 6
explicitly said the force datum should remain the *scalar* `p_cav`
("onset ≠ force reference"). For the closed circumferential surface the
uniform-datum x/y force vanishes analytically (the integrand is
`d/dθ[(R+h) sinθ]`-exact), but the *axial* component
`∫ −(R+h) h_z p_datum dθ dz` does **not** vanish for a tilted bearing
(`h_z ≠ 0`), so switching a tilted case to LOCAL_PSAT shifts `pressure_force_z`
by ~`(p_sat − p_cav) · ∫(R+h)h_z dA`. Decide the convention (recommend:
datum = scalar `p_cav`, per PLAN; or better, ambient/back pressure) — then fix
either the code or the comment at `reynolds.cpp:636-639`, and PLAN/PHYSICS
accordingly.

### D10 — MEDIUM: `snap_tol` is welded to an assumed solver tolerance

`reynolds.cpp:467`: `snap_tol = 5e-7`, justified by `rtol = 1e-8` and
condition ~1e3. `linear_rtol` is user-configurable; loosen it and the snap
stops absorbing solver noise near the free boundary → flag oscillation → the
outer loop hits `max_outer_iters`. Make it `snap_tol = C · linear_rtol · 50`
(documented constant) or a validated config key coupled to `linear_rtol`.

### D11 — MEDIUM: penalty rows `a_p = 1e20`

`reynolds.cpp:152`, `energy.cpp:259-261`. Scale-free 1e20 against stencil
coefficients of wildly different magnitude puts ~17 decades into the matrix;
block-Jacobi/BiCGStab tolerates it today, but it couples badly with D10 and
will bite any future preconditioner work. The off-diagonal links of penalty
rows are also left in place. Replace with exact row replacement
(`a_p = 1`, neighbor links zeroed, `source = value`) — `LinearSystem` already
owns the stencil arrays, so this is mechanical.

### D12 — LOW: temporal-order mismatch and time accumulation

All field equations are BDF1 while the bearing ODE offers RK4
(`journal_motion.cpp:75-94`) — the coupled system is first-order regardless;
either add BDF2 to `ddt_weighted` (three-level, small change) or document that
`dt`-refinement is the only accuracy control. Separately `main.cpp:408`
accumulates `t += dt` in floating point; use `t = step * dt`.

### D13 — LOW: Gumbel-path inconsistencies (fine if Gumbel stays a reference model)

- Transient Gumbel assembles only the squeeze source `−ρ ∂h/∂t`
  (`reynolds.cpp:249-258`); the compressibility storage `h ∂ρ/∂t` has no
  capacity term (Elrod has it). Small because `beta` is large — but it means
  Gumbel transients are not the `beta → beta` limit of Elrod transients.
- INLET_OUTLET degenerates to plain DIRICHLET in the Gumbel path
  (`reynolds.cpp:194-203` adds the link unconditionally), while the Elrod path
  checks flow direction (`reynolds.cpp:442-451`). Also note the Elrod
  INLET_OUTLET decision samples the *previous* iterate's pressure with no
  hysteresis — a boundary cell sitting at the threshold can flip-flop between
  outer iterations.
- Full Sommerfeld negative absolute pressures are knowingly unphysical
  (comment at `reynolds.cpp:275-276`) — fine for a reference model, but
  `validate()` should warn when FULL_SOMMERFELD is combined with gas models.

### D14 — MEDIUM/physics: STRIATED weights shear only; thermal capacity and advection still carry full `ρ·c_p·h` through cavitated cells

`cavitated_film_model = STRIATED` scales Couette shear and heat generation by
`theta` (`energy.cpp:73-75`, `reynolds.cpp:673-674`) — but the thermal capacity
`ρc_p h` (`energy.cpp:318`) and the convective flux `ρc_p h u` faces
(`energy.cpp:338, 355`) are not `theta`-weighted, so a cell that is 30% liquid
transports and stores heat as if it were full. WP-6's title promised
"θ-weighted energy capacity/advection *and* cavitated-zone shear"; only the
shear half landed. Same story for the liquid-solution `rho_cp`/`k` not having
any two-phase (void-weighted) closure. Fix behind the existing STRIATED flag
to preserve the `FULL_FILM ≡ master` gate.

### D15 — LOW/physics: transient capacity treats `ρ₀` as constant in time

`ddt_weighted(theta, rho0·h, dt)` with `rho0 = rho_liquid_solution(p,T,c_d)`
re-evaluated each step: the discrete storage term omits `θh · ∂ρ₀/∂t`. Density
drifts of a few % over a thermal transient leak that much mass-rate into the
balance. The diagnostics residual `r_l` partially hides it because the balance
uses the same convention. Document it in PHYSICS.md at minimum; the clean fix
is assembling `∂(ρ₀θh)/∂t` with `(ρ₀h)ⁿ⁺¹θⁿ⁺¹ − (ρ₀h)ⁿθⁿ` (store `ρ₀ⁿ`).

### D16 — LOW/cleanup: dead code and stale claims

- `EOS::dp_dtheta` (`equation_of_state.hpp:43-45`) and
  `switch_function_smooth` (`:60-62`) are unused — delete or mark reference-only.
- `validate()` clears the error vector, so table errors must be merged *after*
  it (`main.cpp:147-150` comment documents the trap) — make `validate()`
  append-only instead.
- Deprecated `bc_z_*_theta` keys are parsed and silently ignored — emit the
  one-line warning the parser already emits for other retired keys.
- The free-gas advection uses the liquid film-averaged velocity with a pure
  volume-flux face (`gas_transport.cpp:253, 271`) — homogeneous no-slip
  assumption; fine, but PHYSICS.md should say so explicitly next to the ideal-gas
  `rho_g` closure, because the pair is mass-inconsistent across sharp `p`
  fronts (watch `r_gas` in reformation-heavy cases).

---

## 5. Improvement backlog (beyond bug fixes)

Carried over from the 2026-07-07 codebase review, deduplicated against
PLAN.md's open WP-7/8/9/10.

- **B1 — single source of truth for the config surface.** One descriptor table
  (key, member, group, units, help, validation) generating parser + serializer
  + GUI schema. **B1a (do first, one afternoon):** round-trip test
  `parse(to_config_text(cfg))` == fully-populated `cfg`, field by field.
- **B2 — checkpoint/restart.** Binary dump of the six state fields + bearing
  state + `t`/`step` + config hash, every N steps; "restart ≡ continuous" test.
  Prerequisite: **B3**.
- **B3 — extract a `Simulation` class** from `main()` (init / advance_step /
  finalize). Gives restart, in-process GUI embedding, and testable stepping.
- **B4 — typed field names.** `constexpr std::string_view` constants
  (`field_names::pressure` …) used at registration and access; kills the
  runtime-typo class. Full registry redesign not needed.
- **B5 — split `config.hpp`** (1.6 kloc, all inline) into header + TU.
- **B6 — adaptive `dt`** driven by Elrod outer-iteration count and coupling
  residuals; pairs with D6's divergence detection.
- **B7 — transient Picard option** (`coupling_max_iters > 1` in TRANSIENT,
  default 1) for stiff void-coupled cases — with D8 fixed first, else it
  converges on the wrong variable.
- **B8 — `k_dg` sensitivity study** for R290 in the WP-7 harness (bracket
  2–3 decades, record cavitation extent/load response) before anyone invests
  in calibrating it.
- **B9 — MPI coverage for the coupled path**: 1-rank vs 2-rank bitwise(-ish)
  comparison on `config_r290_pz68_quick.txt`; catches D5-class regressions.
- **B10 — PETSc abstraction decision.** Either pimpl `LinearSystem` for the
  stated Trilinos target or amend `ASSESSMENT.md`/README to bless PETSc — the
  current doc/code contradiction misleads newcomers.
- **Longer-horizon physics** (post WP-7/8/9/10): per-cell `p_cav_eff` (D4),
  elastic shell compliance (even 1-D Winkler changes `h` at these loads),
  Patir–Cheng roughness flow factors, x–y cross-coupled bearing DOFs for whirl.

---

## 6. Doc-vs-code mismatches to reconcile (cheap, do alongside fixes)

| Where | Claim | Reality |
|---|---|---|
| PLAN.md WP-12 step 6 | `p_cav_eff` is a ghost-exchanged *field* | uniform scalar from feed state (`config.hpp:1215`) — D4 |
| PLAN.md WP-12 step 6 | force gauge datum stays scalar `p_cav` | code uses `effective_p_cav()` (`reynolds.cpp:640`) — D9 |
| PLAN.md WP-6 | "θ-weighted energy capacity/advection **and** shear" DONE | only shear is θ-weighted — D14 |
| `config.hpp:1213` comment | uniform datum "without injecting a net force" | true for x/y only; axial force shifts under tilt — D9 |
| PHYSICS.md (Elrod normalization) | `ρ = ρ_cav·θ` (canonical) | code uses `ρ = ρ₀(p,T,c_d)·θ`; O(p/β) difference — document, don't fix |
| ASSESSMENT.md | Trilinos is the fixed backend target | PETSc everywhere — B10 |

---

## 7. How to verify your changes

- **Build (Windows):** CMake preset `windows-native-mingw` (MSYS2 MinGW64)
  builds `pancake.exe`, `pancake_gui.exe`, tests, and copies runtime DLLs +
  `data/`. `ctest` needs the mingw64 `bin` on PATH. Linux: `cmake -S . -B
  build-linux && cmake --build build-linux -j`.
- **Unit gates:** the 22 ctest executables; the repo's habit is one focused
  test file per module — extend the matching `tests/test_*.cpp`, don't create
  grab-bag files.
- **Conservation gate:** run with `diagnostics_interval = 1`; `r_l` (and
  `r_gas` for TWO_PHASE) in `results/.../diagnostics.csv` must sit near
  `linear_rtol` in closed/consistent-flux configurations. If your change moves
  it, you broke I4 or I5.
- **Grid gate:** `validation/grid_convergence/run.ps1|sh` (60×20 / 120×40 /
  240×80 on the ε = 0.8 journal case) — cavitation extent, `h_min`, load must
  converge monotonically; archive the summary CSV.
- **MPI gate:** `test_elrod_a2_mpi6`; for anything touching ghosts, also do an
  ad-hoc 1-vs-2-rank field comparison until B9 exists.
- **Model-off gates (I6):** every optional model disabled must reproduce master
  bit-exactly; the acceptance tests in `test_gas_coupling`/`test_coupling`
  show how.
- **Shipped configs:** `config_r290_pz68_quick.txt` is the fast smoke for the
  gas path; the full `config_r290_pz68.txt` is the reference propane case
  (LOCAL_PSAT, JFO, TWO_PHASE — the config most sensitive to D2/D3/D4).

---

## 8. Parting advice

1. Fix in the order D6 → D3 → D2 → D1 → D5. Four of the five are small; all
   five convert silent wrongness into either correctness or loud failure.
   Do B1a the same week.
2. This codebase's superpower is its acceptance-gate culture (I6) and the
   mass-balance diagnostics (I4/I5). Every regression I found lives exactly
   where those gates don't reach: post-processing, initialization, ghost
   ordering, and reversed-flow branches. Extend the gates toward those
   corners rather than adding features first.
3. Distrust any `theta >= 1.0` literal you see — after WP-1, "full film" is a
   per-cell question. Grep for it after every merge (`rg "theta.*>=\s*1\.0"`).
4. The GUI is deliberately an external orchestrator around `config.txt` and
   the solver binary. Keep it that way; every time schema knowledge leaked
   into two places it produced a shipped bug (I7).
5. When physics and PLAN.md disagree with the code, assume the code drifted —
   but write the failing test *before* "fixing" either. Two of the mismatches
   in §6 (force datum, p_cav_eff field) are places where the plan's intent was
   right and the implementation quietly narrowed it.
