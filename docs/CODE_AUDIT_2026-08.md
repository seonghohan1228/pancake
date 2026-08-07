# Code Audit — August 2026

> **Dated audit inventory.** Several quick wins identified here were implemented
> after the audit. Use `PLAN.md` for implementation status and
> `cases/VALIDATION_STATUS.md` for current scientific claims; retain this file as
> the evidence and rationale behind remaining structural work.

Branch `codex/thermal-ksp-robustness` @ `f0e772d`. Read-only audit; no source
was modified.

**Scope.** What the recent campaign (FBNS complementarity cavitation, cross-film
thermal, conjugate bush conduction, KSP robustness) left behind: legacy and
clutter, maintainability, GUI coverage, and the cross-platform contract.

**Method.** `CLAUDE.md` conventions and `ARCHITECTURE.md` first, then the source
under `src/`, `tests/`, `src/gui/`, the shipped `config.txt`, and the 338 case
configs under `cases/`. Line references are to the audited commit.

---

## Verdict

The physics core is in better shape than the campaign's pace would suggest. The
`ReynoldsCoeff` unification is real, every header is guarded, there is not a
single `TODO`/`FIXME`/`HACK` marker in `src/` or `tests/`, and the removed models
(`TVD_MINMOD`, `BUNSEN`, `LOG_MIXING`, `pressure_time_method`) were retired with
warn-and-map parser branches and pinned by tests rather than left to rot.

The damage is concentrated in one place: **the newest and most valuable physics
is undiscoverable.** Cross-film thermal, conjugate bush conduction, the FBNS
numerics block, the coupling controls and the KSP robustness settings are absent
from the shipped `config.txt`, absent from the GUI, and — because the config
parser has no unknown-key branch — a user who misspells one of them gets silence.
The campaign shipped the models and skipped the surface.

Three findings are load-bearing:

| # | Finding | Severity |
|---|---|---|
| A | Coefficient unification is incomplete: `energy.cpp` axial boundaries and all of `gas_transport.cpp` still hand-roll `h^3/(12 mu)` with no cross-film branch | **High** |
| B | 43 of 146 config keys are in neither `config.txt` nor the GUI, and unknown keys parse silently | **High** |
| C | The superseded Elrod-Adams path is still the shipped default while 86% of the project's own cases run FBNS | **Medium-High** |

---

## 1. Legacy and clutter

### 1.1 The Elrod-Adams default contradicts the project's own practice

`cavitation_model` defaults to `JFO` in both places that matter:

- `src/config.hpp:887` — `CavitationModel cavitation_model = CavitationModel::JFO;`
- `config.txt:143` — `cavitation_model = JFO`

Against that, of the 338 case configs under `cases/`:

| Value | Count |
|---|---|
| `FBNS` | 289 |
| `JFO` | 41 |
| `GUMBEL` | 4 |
| `FULL_SOMMERFELD` | 4 |

Every one of the 338 sets the key explicitly, so **flipping the default breaks no
case.** The blast radius is entirely in the test suite: 11 test files never set
`cavitation_model` and would silently change solver path
(`test_bush_conduction`, `test_cross_film`, `test_eos`, `test_fvm`,
`test_gas_transport`, `test_io`, `test_journal_motion`, `test_ksp_guard`,
`test_linear_system`, `test_property_tables`, `test_velocity`).

`JFO` and `ELROD_ADAMS`/`ELROD` are aliases of the same enumerator
(`config.hpp:1123-1131`), not alternatives.

**Recommendation.** Change `config.txt:143` to `FBNS` now — that file is the new
user's starting point and it currently points at the superseded scheme. Leave the
struct default at `JFO` until the 11 tests above pin their model explicitly, then
flip it too. Splitting the change this way keeps the test-path change auditable
instead of implicit.

Also stale: `config.txt:142` still says *"Phase 1: ISOTHERMAL +
fluid_property_model = CONSTANT"* as the FBNS scope, but `config.hpp:1465-1470`
documents FBNS supporting `ENERGY_EQUATION` with `CONSTANT`/`SINGLE_PHASE`/
`TWO_PHASE`. The comment under-sells the solver by two phases.

### 1.2 `TimeSteppingMethod` — the enum is fine, its reachability is the issue

`config.hpp:48` declares five enumerators and `journal_motion.cpp:81-93`
dispatches **all five** to real implementations, with a defensive fallthrough.
Nothing is dead code.

What *is* dead is the whole subsystem's reachability: `JournalMotion::advance`
returns immediately unless the bearing moves (`journal_motion.cpp:133`), and
`motion_model` defaults to `STATIC`. So `motion_time_method` has zero effect on
a default run, and the GUI's "Reference only" label on RK2
(`src/gui/panels_setup.cpp:778`) describes the accuracy, not the status.

**This is not clutter.** The honest comment at `config.hpp:45-47` already says
field solves are backward-Euler and that the integrator applies to the bearing
EOM only. Leave it. The one gap worth closing is test coverage:
`TimeSteppingMethod::EULER_IMPLICIT`, the documented stable default, is never
selected by any test — only `RK4` is (`tests/test_journal_motion.cpp`).

### 1.3 The van Leer far-upwind bug is still present and is **no longer latent**

`src/fvm.cpp:174-193`. The far-upwind cell is hard-coded regardless of flux sign:

```cpp
// east face
double phi_e_ho = phi_e_up + 0.5 * limiter(r, scheme) * (phi_e_up - phi.old(i - 1, j));
// west face
double phi_w_ho = phi_w_up + 0.5 * limiter(r, scheme) * (phi_w_up - phi.old(i - 2, j));
```

For `F_e < 0` the upwind cell is `i+1`, so the far-upwind must be `i+2`, not
`i-1`. Separately, the `is_tvd` block covers east/west only — the north/south
faces get no TVD correction at all, so z-convection silently degrades to
first-order upwind under `VANLEER`.

`CHANGES.md:1412-1414` records this as D1 and calls it *"latent for omega > 0"*;
`PLAN.md:2098` still has it as F4 = TODO. **That caveat is now out of date.** It
holds only for the Elrod theta equation, where the face flux is pure Couette
(`reynolds_coefficients.hpp:104`) and therefore carries the sign of `omega`. On
the two paths the campaign added or extended:

- **Energy** (`energy.cpp:542`): `flux_theta = 0.5*omega*R*couette_face -
  g_flow*f2_face*dp_dtheta/R`. The Poiseuille term can flip the face flux
  negative **with `omega > 0`**, and `thermal_convection_scheme = VANLEER`
  passes straight through (`energy.cpp:580-582` downgrades only
  `TYPE_DIFFERENCING`).
- **Gas transport** (`gas_transport.cpp:283, :309`): same structure via
  `gas_convection_scheme`.

Also note `omega < 0` is supported (`config.hpp:1516` uses `std::abs(omega)`,
and `validate()` never requires positivity), which puts every Elrod face on the
buggy branch.

**This should be reclassified from latent to live before the next validation
run.**

### 1.4 `beta` and `bulk_modulus`

`beta` is **not a config key** — no `key == "beta"` branch exists anywhere in the
parser. It appears only as a mathematical symbol in comments and as unrelated
local names. A user writing `beta = 5e9` gets it silently dropped (see §2.2).

`bulk_modulus` is live for Elrod and genuinely ignored by FBNS: `fbns.cpp`
contains zero references to it, and `fbns.hpp:25` says so. Every consumer
(`theta_from_pressure`, `elrod_boundary_outflow`, `void_params::beta_bar`) is
reached only from `solve_elrod` or from `Diagnostics` under an explicit
`FBNS::active(cfg)` / `JFO` gate.

Two residues:

1. `main.cpp:183-185` special-cases it — `FBNS::active(cfg) ? 1.0 :
   EOS::theta_from_pressure(..., cfg.bulk_modulus)`.
2. `validate()` still runs `require_positive(bulk_modulus, ...)`
   (`config.hpp:1288`) **unconditionally**, so an FBNS run with
   `bulk_modulus = 0` errors out on a parameter it never uses. Gate that check on
   the cavitation model — small fix, removes a real confusion for FBNS users.

### 1.5 `BUNSEN` — a warning every single run emits

`config.hpp:1048-1052` warns that `dissolved_gas_bunsen_coeff` is ignored because
the BUNSEN solubility model was removed. The parser handling is correct. The
problem is that the key is still written into **every shipped config** —
`config.txt:76`, `config_r290_pz68.txt:65`, `config_r290_pz68_quick.txt:59` — and
~40 files under `cases/`, all with the value `0`. So the warning fires on every
run and trains users to ignore parse warnings.

**Quick win:** delete the key from the shipped configs. Keep the parser branch.

Stale doc to match: `PHYSICS.md:1432` still prints the BUNSEN formula as one of
three solubility branches when the code offers two (`SolubilityModel { HENRY,
TABLE }`, `config.hpp:63`). Same section shows `LOG_MIXING` and
`EMPIRICAL_CORRELATION` as distinct rows with distinct formulas; both now alias
to `LiquidViscosityModel::EMPIRICAL` (`config.hpp:642-645`), so selecting
`LOG_MIXING` gives the full Andrade x Barus law, not the reduced one the doc
shows.

### 1.6 Genuinely dead members

| Key | Status |
|---|---|
| `bc_z_south_theta` / `bc_z_north_theta` | Declared `config.hpp:876, :880`, assigned `:1113, :1121`, **never read anywhere**. Self-documented as deprecated and pinned by `tests/test_elrod_bc.cpp:114` and `tests/test_energy.cpp:558`. Deliberate, tested dead code — safe to delete the members and parser branches whenever the back-compat window closes. |
| `load_angle_deg` | The real leak. Declared `config.hpp:678`, assigned `:988`, echoed at `:2013`, rendered as an editable GUI spinbox (`src/gui/schema.cpp:92-93`) — and consumed by nothing. `PHYSICS.md:21` and `README.md:540` both describe it as *"a visualization-only reference angle"*; `ARCHITECTURE.md:27` claims the GUI "rotates the preview around a load-reference angle". No such rotation exists in the code. Either wire it into the preview or delete it and correct all three docs. |

---

## 2. Maintainability

### 2.1 The coefficient unification is incomplete — and the gap is the exact defect the header exists to prevent

`src/reynolds_coefficients.hpp:16-29` states its own purpose plainly: three call
sites must agree or the solver lies about itself, and when they drifted,
`mass_residual` "reported a fabricated imbalance on exactly the runs it exists to
police."

Three files honour it — `reynolds.cpp` (8 calls), `diagnostics.cpp` (8),
`fbns.cpp` (3). **Two do not.**

**`src/energy.cpp`.** Interior faces branch correctly on cross-film
(`:533-540`, `:566-568`), and so does the dissipation term (`:196-197`). But the
axial boundary flux helpers `south_boundary_flux` (`:245`) and
`north_boundary_flux` (`:270`) compute the boundary velocity as

```cpp
u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
```

with **no cross-film branch at all** — `grep cross_film_f src/energy.cpp`
returns `132, 196, 533, 539, 566, 567` and nothing in the 245-295 range. So with
`cross_film_thermal_model` active, the energy equation transports heat through
its interior faces with Dowson's `F2` and through its axial boundary faces with
the classical `h^3/(12 mu)`. The same inconsistency applies to
`cfg.reformation_boundary_velocity`, which also takes `(h_face, mu_face)`.

This matters for exactly the case the campaign was built for: the Ferron
cross-film run has Dirichlet axial ends, so the mismatched faces are the oil
outflow.

**`src/gas_transport.cpp`.** Zero occurrences of `CrossFilm` or `cross_film` in
the file. `:247-248` hand-rolls both terms:

```cpp
const double u_theta = 0.5 * cfg.omega * R
    - g_face * h_face * h_face * dp_dtheta / (12.0 * mu_face * R);
```

Neither the Couette thickness (`2 F1/F0`) nor the Poiseuille mobility (`F2`) is
cross-film aware, so gas is advected with a velocity field inconsistent with the
pressure solve that produced it. Currently unexercised only because gas transport
plus cross-film has no test (§2.4) — it is a trap, not yet a failure.

**Recommendation.** Route both files through `ReynoldsCoeff::couette_face_flux`
and `ReynoldsCoeff::poiseuille` / `poiseuille_coefficient`. The header's own
contract says the classical branch is bit-identical when the model is off, so
existing isothermal and non-cross-film results should not move. Verify that claim
with a bit-comparison run rather than assuming it.

### 2.2 `config.hpp` — 2322 lines, and it is a header

The single largest file in the tree, and it is `#include`d by 43 translation
units, pulling `<iostream>`, `<fstream>`, `<filesystem>` and `<iomanip>` into all
of them.

| Region | Lines | Note |
|---|---|---|
| Enums + `to_config_value` overloads | ~35-500 | 24 enums, one stringify overload each |
| `load_from_stream` | 956-1208 | **146** `key == "..."` branches, one flat `else if` chain |
| `validate()` | **1219-1747** | **529 lines in one function**, 67 `warn`/`error`/`require_*` calls |
| `to_config_text` | ~1960-2300 | key-by-key writer |

`validate()` is the sprawl. It is not badly written — the lambdas at `:1224-1244`
are the right idea — but 529 lines of straight-line cross-model consistency
checks in one function is the ceiling. It is also where the campaign's knowledge
lives: the warnings at `:1308`, `:1321`, `:1327`, `:1440`, `:1492` are the only
place a user learns that a THD-looking run with `CONSTANT` properties has no
thermal feedback, or that `theta_convection_scheme` is ignored under FBNS. That
knowledge is worth keeping; it just needs to be split by subsystem
(`validate_geometry`, `validate_cavitation`, `validate_thermal`,
`validate_gas`, `validate_linear_solver`) so each is reviewable.

**The parser has no terminal `else`.** `config.hpp:1206-1208` closes the chain
without a fallback, so any unrecognized key is silently discarded. This is the
root cause that lets the dead-key class regrow, and it already bites: three
warn-and-drop keys are handled explicitly, but everything else — a typo, a
removed key, a key from a newer branch — vanishes without a trace. Live example:
`cases/01_ferron_thd/config_template.txt:51` sets `cross_film_model = NONE`, but
the parser only knows `cross_film_thermal_model`. That line has been dropped on
every run, and is replicated across ~24 generated `fb_iter_*/config.txt` files.
Harmless today only because `NONE` is the default.

**Highest-value single fix in this audit:** add the terminal `else` that pushes
an "unrecognized key" entry into `parse_warnings`. Perhaps 10 lines. It makes
every future dead key self-reporting.

### 2.3 43 of 146 config keys are absent from the shipped `config.txt`

Parser: 146 distinct keys. Shipped `config.txt`: 103. The 43-key gap is not
random — it is precisely the campaign's output:

- **Cross-film (6):** `cross_film_thermal_model`, `cross_film_profile_order`,
  `cross_film_quadrature_points`, `cross_film_relaxation`,
  `shaft_ambient_heat_transfer`, `shaft_ambient_area_ratio`
- **Bush conduction (5):** `bush_thermal_model`, `bush_thermal_conductivity`,
  `bush_outer_radius`, `bush_outer_heat_transfer`, `bush_radial_cells`,
  `bush_outer_boundary`
- **Coupling / steady loop (3):** `coupling_max_iters`, `coupling_tolerance`,
  `coupling_relaxation`
- **KSP robustness (5):** `linear_ksp_type`, `linear_pc_type`,
  `linear_pc_overlap`, `linear_solver_fallbacks`, `linear_rtol`
- **Cavitation / FBNS scope (3):** `cavitated_film_model`,
  `cavitation_threshold`, `consistent_boundary_flux`
- Plus `groove_mixing_model`, `steady_gas_model`, `gas_pressure_coupling`,
  `thermal_ambient_temperature`, the `*_table_file` trio, and the
  `pressure_force_*` output names.

Combined with §2.2, a user who wants the conjugate bush model has no example to
copy, no GUI control, and no error if they guess the key name wrong.

### 2.4 Test coverage — 34 ctest cases, well built, with four structural holes

The suite is in good shape: 25 test sources, 25 registrations, no orphans; 33 of
34 cases run under MPI; FBNS is policed at 1, 2, 4 and 6 ranks. The
rank-independence work has real teeth.

Gaps, ranked by how much they would hurt:

| Gap | Evidence |
|---|---|
| **Transient THD — uncovered** | The `SolutionMode != STEADY_STATE` branch at `energy.cpp:592-594` calls `FVM::ddt_weighted`, which appears **nowhere** in `tests/`. All 13 `Energy::solve` call sites across `test_energy.cpp`, `test_cross_film.cpp` and `test_fbns.cpp` run `STEADY_STATE`. The thermal-capacity field built at `energy.cpp:491-505` is assembled but its time derivative is never exercised. |
| **`MOVING_BEARING` + any pressure solve — uncovered** | `test_journal_motion.cpp:49-108` unit-tests the geometry kernel with no solve. `test_elrod_bc.cpp:198-225` sets the flag but zeroes `dh_dt`, sets `omega = 0`, and never calls `compute_moving_bearing`. No FBNS test sets `motion_model` at all — `test_fbns.cpp:369-410` fakes the squeeze by hand under `STATIC`. The real loop (`main.cpp:109-114` → `fbns.cpp:186-187` → `main.cpp:475-479`) never runs. The combination is legal: FBNS scope validation (`config.hpp:1465-1499`) says nothing about `motion_model`. |
| **Gas transport + cross-film — uncovered, and structurally blocked** | `GasTransport::solve` requires `TWO_PHASE` **and** `dt > 0` (`gas_transport.cpp:201-203`), i.e. transient-only; every cross-film test is `STEADY_STATE`. Same root as the transient-THD gap. This is the untested combination that §2.1 identifies as already inconsistent. |
| **Bush `CONDUCTION_2D` + FBNS — uncovered** | FBNS + cross-film + **LOCAL** bush is well covered (`test_fbns.cpp:884-1006`, ranks 1/2/4/6). But `test_fbns.cpp` never sets `bush_thermal_model`, so the 2-D annulus is only tested standalone or against a bare `Energy::solve`. The all-reduce/redistribute at `cross_film.cpp:84, :313-331` never runs through a complementarity solve. |

Smaller items: `ConvectionScheme::LINEAR`, `DensityModel::TABLE`,
`OutputVerbosity` and `TimeSteppingMethod::EULER_IMPLICIT` are never selected;
`INLET_OUTLET` appears only as a rejected-scope case; `FULL_SOMMERFELD` is never
solved. `linear_solver_options_from_config` (`linear_system.hpp:40-57`) is called
by no test, and the non-default preconditioners are never rank-tested — despite
`config.hpp:105-106` recording that `BJACOBI` "is known to diverge at some rank
counts on the steady energy and Elrod equations," which is the exact regression
class `tests/CMakeLists.txt:42-45` was written to catch.

`src/main.cpp` has **no end-to-end test**: the outer Picard loop (`:293-393`),
the convergence criterion (`:386-391`), `omega_at_time` ramping (`:268`) and the
`stop_on_nonpositive_film` guard (`:266-276`) are never driven by ctest.
`test_fbns.cpp:611-676` hand-reimplements a mirror of the sweep, so the two can
drift silently.

There is no restart/checkpoint feature — not a coverage gap, an absent feature.
Worth noting given how long the validation runs are.

### 2.5 Error messages and header hygiene

**Messages are a strength.** The `validate()` warnings are specific, name the key,
and say what will happen — `"viscosity is frozen at mu, so the temperature field
has NO effect on pressure or load"` (`:1308`) is exactly right for a user who
thinks they are running THD. Keep this standard.

**Header hygiene is clean.** All 31 headers under `src/` and `src/gui/` carry
`#pragma once`; include counts are modest (1-6, except `config.hpp`'s 13). The
only issue is `config.hpp`'s size and reach (§2.2).

---

## 3. GUI (`src/gui/`)

### 3.1 Coverage: 88 keys exposed of 146 parsed

| Source | Keys |
|---|---|
| `schema.cpp:78-264` `double_specs()` | 57 |
| `schema.cpp:266-281` `int_specs()` | 4 |
| `schema.cpp:283-299` `bool_specs()` | 4 |
| `panels_setup.cpp:461-814` enum combos | 17 |
| `panels_setup.cpp:713-761, :816-921` paths + inlet editor | ~6 |

About 50 keys are unreachable; ~34 are real solver knobs (the rest is the
deliberate motion-vector isotropy collapse).

**Missing, by subsystem:**

- **Cross-film — 100% missing.** No widget for any of the six keys.
- **Bush conduction — 100% missing.** No widget for any of the six keys.
- **FBNS numerics — 100% missing.** `fbns_max_iters`, `fbns_tol`,
  `fbns_relaxation`, `fbns_pressure_scale`, `cavitated_film_model`,
  `consistent_boundary_flux`. This one is actively misleading: the GUI offers
  `FBNS` in the cavitation dropdown (`panels_setup.cpp:499`), but the only two
  convergence controls it exposes — `max_outer_iters` and `outer_tol` — are
  explicitly *not used* on the FBNS path (`config.hpp:901-903`). **A GUI user who
  selects FBNS is editing two dead controls and cannot reach the four live ones.**
- **Coupling — missing.** `coupling_max_iters`/`_tolerance`/`_relaxation` drive
  the steady THD Picard loop; none reachable.
- **KSP — only `linear_rtol` survives.** The ASM-vs-BJACOBI choice documented at
  `config.hpp:99-109` as the fix for rank-dependent divergence is invisible —
  precisely the knob a GUI user hitting a failed run would need.
- Also missing: `thermal_convection_scheme`, `gas_convection_scheme`
  (`panels_setup.cpp:783-786` acknowledges this, but that comment predates the
  THD and gas models), `output_fields`, `thermal_ambient_temperature`,
  `steady_gas_model`, `void_coupling_relaxation`.

### 3.2 Enum coverage: no stale lists, but 8 enums have no widget

Good news first — of the 17 enums that *do* have a combobox, **every one lists
all its enumerators.** `CavitationModel` includes FBNS, `TimeSteppingMethod` has
all five, `ConvectionScheme` all four. Nothing has gone stale.

Eight enums have no combobox at all: `CrossFilmModel`, `BushThermalModel`,
`BushOuterBoundary`, `GrooveMixingModel`, `KrylovMethod`, `LinearPreconditioner`,
`SteadyGasModel`, `CavitatedFilmModel`.

### 3.3 Output fields: viewing works, selecting does not

There is no `output_fields` editor anywhere. The string appears in
`panels_results.cpp` only inside warnings (`:464`, `:500`, `:560`, `:719`) that
tell the user a field is missing with no way to add it. Users get whatever
`default_output_fields()` (`config.hpp:1968-2001`) ships, which contains none of
the cross-film or bush fields.

Viewing is fine — the field picker is built dynamically from the VTK arrays
(`panels_results.cpp:232-261`). But the two hardcoded metadata chains,
`field_unit()` (`:47-65`) and `field_desc()` (`:68-102`), have already drifted:
they know `T_journal`/`T_bearing`/`q_journal`/`q_bearing` but not
`cross_film_t_shaft`, `cross_film_t_bush` or `cross_film_core_rise` — the names
`cross_film.cpp:108-111` actually writes and the crossfilm validation case
actually uses. Those three render with no unit and no tooltip, along with ~10
others (`p_void`, `rho_gas`, `free_gas_mass`, `dh_dt`, `cross_film_f0/f1/f2`, …).

### 3.4 Can a new user set up the Ferron THD case from the GUI alone? **No.**

`cases/01_ferron_thd/config_template.txt` — blocked by 5 keys:

| Key | Line | Why it blocks |
|---|---|---|
| `coupling_max_iters = 200` | `:29` | `force_balance.py:94` passes 200; the GUI is stuck at the default 30 (`config.hpp:918`). The steady THD Picard loop will not converge in 30. |
| `coupling_tolerance`, `coupling_relaxation` | `:30-31` | Not exposed |
| `groove_mixing_model` | `:129` | Case runs a two-pass `PINNED` → `FLOW_WEIGHTED` closure (`force_balance.py:12-18, 163-172`). No selector. |
| `output_fields` | `:135` | Needs `mu` and `friction_torque`; not editable either way |

Additionally the groove line at `:126` carries `mixing_lambda` and
`supply_flow_rate` as 5th and 6th tokens (`config.hpp:1204-1205`), and the GUI
inlet editor (`panels_setup.cpp:881-902`) exposes neither — so `FLOW_WEIGHTED`
would be unusable even if the model selector existed.

`cases/01_ferron_thd_crossfilm/` is blocked by 13 keys: everything above plus the
entire conjugate wall block, none of which has any GUI representation.

A GUI-only user can build the geometry, the Ferron polynomial viscosity
(`schema.cpp:208-218` even names the k0/k1/k2 values), FBNS cavitation and the
energy equation — then stalls at the coupling controls.

### 3.5 One data-loss bug

`normalize_motion_config()` (`gui_app.cpp:525-531`) is called on **load**
(`:562`). Opening an externally authored config overwrites
`bearing_stiffness_y/z` and `bearing_damping_y/z` with the x value, forces
`bearing_initial_from_attitude = true`, and zeroes all six `bearing_initial_*`.
Save (`:579`) then persists the mutilated version.

A config with anisotropic support (`kx = 1e6`, `ky = 2e6`) loaded and saved in
the GUI comes back isotropic, with no warning. The isotropy contract is a
documented design choice (`schema.cpp:244-247`), but applying it destructively on
load rather than only on GUI-originated edits is data loss.

### 3.6 Structural smell: three hand-synced parallel lists

Nothing detects a new config key. `double_specs`/`int_specs`/`bool_specs`
(`schema.cpp:79-299`) vs the parser (`config.hpp:983-1207`); `field_unit()` vs
`field_desc()` (`panels_results.cpp:47-102`); and the enum option lists
(`panels_setup.cpp:474-806`), which duplicate every enumerator by hand with no
`static_assert` or count check. This is exactly how the ~34 gaps in §3.1
accumulated, and it will happen again.

`validate_config()` (`schema.cpp:351-425`) inherits the problem:
`field_for_message()` (`:306-326`) can only highlight a field that appears in one
of the three spec lists, so a solver error about e.g. `bush_outer_radius` shows
as list-only text with no control to fix it.

`docs/GUI_ARCHITECTURE.md` (65 lines) has zero mentions of cross-film, bush
conduction, FBNS numerics or KSP settings, and no "known coverage gaps" section
recording that these are unexposed.

---

## 4. Cross-platform contract

**No drift introduced by the campaign.** This held up well.

- Zero POSIX-only headers anywhere in core `src/` (`unistd.h`, `dirent.h`,
  `sys/stat.h`, `pwd.h` — no hits outside `src/gui/`).
- `windows.h` appears only in `src/main.cpp` and `src/gui/`, correctly guarded
  with `NOMINMAX` (`main.cpp:29-33`).
- No hardcoded backslash separators.
- The new modules (`cross_film.cpp`, `bush_conduction.cpp`, `fbns.cpp`,
  `gas_transport.cpp`) do no file I/O and introduce no platform surface.

One **pre-existing, not campaign-introduced** inconsistency: `src/io.cpp` mixes
path styles. Newer code uses `std::filesystem` correctly (`:422`, `:440` —
`root / "results.pvd"`), but `update_pvd` (`:118`, `:125`, `:136`),
`update_flat_pvd` (`:147`) and the `write_*_vts` functions (`:170`, `:277`,
`:348`) build paths by `std::string` concatenation with `"/"`, because
`cfg.output_dir` is a `std::string` (`config.hpp:925`). This works on Windows
(Win32 accepts forward slashes) but contradicts `ARCHITECTURE.md:46`. Low
priority; worth folding into any future `io.cpp` work.

---

## 5. Prioritized recommendations

### Quick wins — hours, low risk

| # | Action | Effort | Risk | Payoff |
|---|---|---|---|---|
| Q1 | Add the terminal `else` to the config parser that pushes an "unrecognized key" warning (`config.hpp:1206-1208`) | ~1 h | **Very low** | Highest ratio in this audit. Makes every future dead key self-reporting; immediately surfaces the `cross_film_model` typo in ~25 case files |
| Q2 | Delete `dissolved_gas_bunsen_coeff` from `config.txt` and the three `config_r290_*.txt` | 15 min | Very low | Stops a spurious warning on every run |
| Q3 | Add the 43 missing keys to `config.txt` as commented defaults with one-line explanations | 2-3 h | Very low | Makes cross-film, bush and KSP discoverable at all |
| Q4 | Gate `require_positive(bulk_modulus)` on `!FBNS::active(cfg)` (`config.hpp:1288`) | 15 min | Very low | Removes a confusing hard error on an unused parameter |
| Q5 | Switch `config.txt:143` to `FBNS`; fix the stale Phase 1 comment at `:142` | 15 min | Very low — all 338 cases set the key explicitly | New users land on the validated path |
| Q6 | Fix `PHYSICS.md:1432` (BUNSEN) and the `LOG_MIXING`/`EMPIRICAL_CORRELATION` rows | 1 h | None | Removes a documented model that does not exist |
| Q7 | Add `cross_film_t_shaft`, `cross_film_t_bush`, `cross_film_core_rise` (+ the ~10 others) to `field_unit()`/`field_desc()` | 1 h | Very low | Cross-film results become readable in the GUI |
| Q8 | Resolve `load_angle_deg`: wire it into the preview or delete it and correct `PHYSICS.md:21`, `README.md:540`, `ARCHITECTURE.md:27` | 1-2 h | Low | Removes a documented feature that does not exist |

### Structural — days, needs care

| # | Action | Effort | Risk | Notes |
|---|---|---|---|---|
| S1 | **Route `energy.cpp` axial boundaries and `gas_transport.cpp` through `ReynoldsCoeff`** (§2.1) | 1-2 d | **Medium** — touches assembled fluxes | Do this first among the structural items. The header claims bit-identical reduction when cross-film is off; **verify that with a bit-comparison run, do not assume it**. Any change to non-cross-film results is a bug in the change, not the baseline |
| S2 | **Fix the van Leer far-upwind indexing and add the z-face TVD branch** (`fvm.cpp:174-193`, D1/F4) | 1-2 d | **Medium** | Reclassify from latent to live first (§1.3) — it is reachable on the energy and gas paths with `omega > 0`. Needs a reversed-flow test in `test_schemes.cpp` before the fix, so the fix is provably a fix |
| S3 | Add a transient-THD test | 0.5-1 d | Low | Closes the single largest coverage hole and unblocks S4 |
| S4 | Add gas-transport + cross-film and `CONDUCTION_2D` + FBNS integration tests | 1-2 d | Low | S4 is what would have caught S1 |
| S5 | Add a `MOVING_BEARING` + FBNS coupled test driving the real loop | 1 d | Low | Currently a legal, completely untested combination |
| S6 | Split `validate()` (529 lines) into per-subsystem functions | 1 d | Low — pure refactor | Preserve every warning verbatim; they are the campaign's institutional memory |
| S7 | Generate the GUI spec lists from a single key manifest shared with the parser, or add a build-time completeness check | 2-4 d | Medium | Fixes the root cause of §3.6 rather than the symptom. If too large, a test asserting "every parser key appears in a spec list or a documented exclusion list" gets 80% of the value in ~4 h |
| S8 | Add GUI panels for cross-film, bush conduction, FBNS numerics and coupling controls | 2-3 d | Low | After S7, so the new panels cannot drift. Enough to make the Ferron THD case GUI-reachable |
| S9 | Fix `normalize_motion_config()` to apply on edit, not on load (`gui_app.cpp:525-562`) | 0.5 d | Low | Data-loss bug |
| S10 | Add an end-to-end `main.cpp` test (outer loop, convergence, ramping, film guard) | 1-2 d | Low | Stops `test_fbns.cpp:611-676`'s hand-written mirror from drifting |
| S11 | Add odd-rank variants for `test_mass_residual`, `test_diagnostics`, `test_gas_transport` | 0.5 d | Low | These do global reductions and ghost-dependent flux sums — exactly the class the rank work targeted |

### Suggested order

`Q1` → `Q2, Q4, Q5, Q6` → `S2` (reclassify, test, fix) → `S3` → `S1` (with the
bit-comparison guard) → `S4` → `Q3, Q7` → `S7` → `S8`.

`Q1` first because it costs an hour and changes what every subsequent
investigation can see. `S2` before `S1` because the van Leer bug affects the same
flux paths `S1` touches, and fixing them in the other order makes the
bit-comparison ambiguous.
