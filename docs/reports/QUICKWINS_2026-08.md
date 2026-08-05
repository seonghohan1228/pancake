# Audit quick-win batch — August 2026

Branch `codex/audit-quickwins`, from `main` @ `79ddbee`. Five changes, five
commits, full `ctest` green (35/35) after each one and at the end.

Inputs: `docs/CODE_AUDIT_2026-08.md` (items Q1, Q5, findings A and C, §3.5) and
`PLAN.md` priority-queue item 2.

---

## Summary

| # | Change | Audit ref | Commit | Result |
|---|---|---|---|---|
| 1 | Unknown config key is a hard error, with a typo suggestion | Q1 | `9beb61d` | New test file; 105 case configs corrected |
| 2 | van Leer far-upwind cell follows the flux sign | C / F4 | `1c4e5fa` | Test fails on old code, passes on new; no existing test's numbers move |
| 3 | GUI stops rewriting motion settings on load and save | 3.5 / S9 | `3b0a3e3` | Data-loss bug closed; no test target exists for it |
| 4 | Shipped `config.txt` defaults to FBNS | Q5 | `e953c18` | Struct default deliberately left at `JFO` |
| 5 | Energy axial boundaries + gas transport onto `ReynoldsCoeff` | A / S1 | `13cab1d` | A/B run: moves at the output precision floor except two boundary-derived fields |

Full per-item rationale is in `CHANGES.md`. This report records the evidence and
the things deliberately left undone.

---

## 1. Unknown-config-key hard error (Q1)

The parser chain had no terminal `else`. Unknown keys now populate
`SimulationConfig::parse_errors`, which `validate()` promotes into its error
list, so both the solver's `main` and the GUI's config check refuse to run.

**The campaign had shipped a live example.**
`cases/01_ferron_thd/config_template.txt` set `cross_film_model = NONE` while the
parser has only ever known `cross_film_thermal_model`. The line was dropped on
every run and replicated into **105** generated `fb_iter_*/config.txt` files. The
value was `NONE`, the default, so no archived result is affected — but nothing in
the toolchain could have told anyone that. All 105 files and the four mentions in
`cases/01_ferron_thd/SUMMARY.md` now use the real key.

**Why the suggestion needs two signals.** `cross_film_model` is eight edits from
`cross_film_thermal_model`, far outside any sane distance threshold, but shares
all three of its underscore-separated words. Edit distance alone would have missed
the one case that actually shipped. The heuristic accepts a candidate when
`distance <= max(2, len/3)` **or** it shares at least half the unknown key's
words. `totally_made_up_key` matches neither and gets no suggestion.

**Test** (`tests/test_config_keys.cpp`, 5 cases): the shipped line verbatim; the
legacy branches; manifest-vs-parser agreement; `to_config_text()` round-trip; and
every `config*.txt` in the repo root and under `cases/` — templates included, with
`{placeholder}` values substituted, because the templates are where the dropped
line lived. `config.txt` is additionally required to pass `validate()`.

**Found while writing the manifest, not fixed:** a *recognized* key with a
malformed value still throws `std::invalid_argument` out of `std::stod` rather
than being reported. `cases/*/config_template.txt` files hit this, which is why
the test substitutes their placeholders. Worth a follow-up.

---

## 2. van Leer far-upwind cell (finding C / PLAN F4)

`FVM::divergence` hard-coded the far-upwind cell at `i-1` (east) and `i-2` (west)
regardless of the flux sign, and for a negative flux inverted the positive-flux
gradient ratio. For `F < 0` the upwind cell is the `+theta` neighbour and the
far-upwind cell is one further in that direction.

**Reclassified from latent to live.** The recorded "latent for `omega > 0`"
caveat held only for the Elrod theta equation, whose faces are pure Couette. The
energy (`energy.cpp`) and gas-transport (`gas_transport.cpp`) face fluxes carry a
Poiseuille term that can reverse a face with `omega > 0`, and both accept
`VANLEER` through their own scheme keys. `omega < 0` is also legal, which puts
every Elrod face on the reversed branch.

**Test written first, and it fails on the old code** — the PLAN F4 acceptance
gate ("mirrored one-dimensional advection with positive and negative flux
produces mirrored results"):

| | forward | backward, old code | backward, fixed |
|---|---|---|---|
| L1 error | 0.0266008 | 0.0626316 | 0.0266008 |
| range | — | [0.140741, **1.00469**] | [0.2, 1.0] |

2.4x the error and a genuine overshoot past the initial maximum: the reversed
branch was not bounded either. Registered at 1, 4 and 6 ranks.

**Did the fix change any green test's results? No.** The full 35-case suite was
run verbosely on the pre-fix and post-fix binaries and the outputs diffed. The
only differences are the three new directional-test lines and `test_ksp_guard`'s
nondeterministic MPI abort text. This is the expected answer: every currently
tested `VANLEER` path has non-negative face fluxes, so the fix is a no-op on them.

**Not fixed:** the deferred correction still covers circumferential faces only, so
axial convection degrades to first-order upwind under `VANLEER`. That is the
second half of F4. It is now stated in `PHYSICS.md` §5.5 rather than left
implicit, and F4 is marked PARTIAL in `PLAN.md`.

---

## 3. GUI motion-config data loss (§3.5 / S9)

`GuiApp::normalize_motion_config()` ran on **load** and again on **save**. Opening
an externally authored config overwrote `bearing_stiffness_y/z` and
`bearing_damping_y/z` with the x value, forced `bearing_initial_from_attitude =
true` and zeroed all six `bearing_initial_*`; saving persisted that. A config with
`kx = 1e6, ky = 2e6` came back isotropic with no warning.

The isotropy contract itself is legitimate — the editor exposes one support row
for all three directions. The defect was applying it destructively to data from
elsewhere.

| when | before | now |
|---|---|---|
| New case | normalize | normalize (the GUI authors this config) |
| Load | normalize | preserve; log what the editor cannot represent |
| Edit support row | — | mirror x onto y/z (`mirror_isotropic_support`) |
| Save | normalize | write what is in the config |

The load path also now surfaces `parse_errors` from item 1 in the GUI log.

**Verification gap, stated plainly:** `GuiApp` links Win32/ImGui and has no test
target, so this was verified by inspection and a clean GUI build, not by an
automated test. Nothing in `ctest` covers the GUI config round trip. Closing that
would mean extracting the config-mutation policy out of `GuiApp` into a testable
free function — worth doing, out of scope here.

---

## 4. Shipped default flipped to FBNS (Q5)

`config.txt` is the new user's starting point and pointed at the superseded
Elrod-Adams path while 289 of the project's 338 case configs run FBNS. All 338 set
the key explicitly, so no case changes.

Two comment defects fixed with it: the FBNS scope was described as *"Phase 1:
ISOTHERMAL + fluid_property_model = CONSTANT"*, two phases out of date; and
`max_outer_iters` / `outer_tol` sat under the model with nothing saying they are
Elrod-only.

**The struct default stays `JFO`, as instructed.** These eleven test files never
set `cavitation_model` and would silently change solver path if it flipped:

`test_bush_conduction`, `test_cross_film`, `test_eos`, `test_fvm`,
`test_gas_transport`, `test_io`, `test_journal_motion`, `test_ksp_guard`,
`test_linear_system`, `test_property_tables`, `test_velocity`.

(The new `test_config_keys` also relies on the default, but only through
`to_config_text()` — it exercises no solver path.)

Pinning those eleven is the prerequisite for flipping the struct default, and it
is a separate reviewable change.

**Left alone:** the three `config_r290_*.txt` still run `JFO`. Those are gaseous
cavitation cases and the R290 work is where `VOID_COUPLED` lives, which FBNS
refuses. Flipping them is a modelling decision, not a default cleanup.

---

## 5. Coefficient unification (finding A / S1)

`energy.cpp`'s `south_boundary_flux` / `north_boundary_flux` and every face of
`gas_transport.cpp` hand-rolled `h^3/(12 mu)` with no cross-film branch. Both now
call `ReynoldsCoeff::couette_face_velocity` and
`ReynoldsCoeff::poiseuille_face_velocity`.

**Bit-identity when the model is off is verified, not assumed.** The helpers keep
the classical expression's original floating-point association. Full 35-case suite
captured verbatim before and after: the only differing line is `test_ksp_guard`'s
nondeterministic abort text.

### A/B on the archived cross-film point

`cases/01_ferron_wallbc_fields/rerun_xf_adiabatic/fb_iter_13` — Ferron 1983
Fig. 4 (2000 rpm, 4 kN), 720x192, FBNS + `ENERGY_EQUATION` +
`cross_film_thermal_model = POLYNOMIAL`, `bush_thermal_model = LOCAL`, Dirichlet
axial ends, `groove_mixing_model = FLOW_WEIGHTED`; the converged force-balance
iterate (`Fx` = 0.028 N, `W` = 4000.07 N). Rerun locally on 4 ranks with the
pre-change and post-change binaries.

The pre-change binary reproduces the archived 24-rank numbers (`p_max` 1.41105e6,
`T_max` 321.768 K, `friction_torque` 2.59609, `W` 4000.06 N), so the A/B is
against the real archived point and not a locally drifted one.

VTK output carries 6 significant digits — that is the resolution floor of this
table. Over 138 240 cells:

| field | scale | max abs diff | max rel diff | reading |
|---|---|---|---|---|
| `p` | 1.41105e6 Pa | 10 Pa | 7.1e-6 | one printed digit |
| `Fp` (load) | 4000.06 N | 0.01 N | 2.5e-6 | one printed digit |
| `F` | 3985.55 N | 0.01 N | 2.5e-6 | one printed digit |
| `T` | 321.768 K | 0.001 K | 3.1e-6 | one printed digit |
| `T_bearing` / `cross_film_t_bush` | 323.953 K | 0.001 K | 3.1e-6 | one printed digit |
| `friction_torque` | 2.59609 N m | 1e-5 | 3.9e-6 | one printed digit |
| `mu` | 0.0272947 Pa s | 2e-7 | 7.3e-6 | one printed digit |
| **`cross_film_core_rise`** | 1.06394 K | 5.6e-5 K | **5.3e-5** | above the floor |
| **`q_journal`** | 14467.1 W/m^2 | 0.6 | **4.1e-5** | above the floor |
| `film_content` | 1 | 3e-5 | 3.0e-5 | one front cell |
| `h`, `T_journal`, `cross_film_t_shaft` | — | 0 | 0 | unchanged |
| `mass_residual` | ~1e-10 | — | >1 | noise about zero |
| `q_bearing` | ~3e-9 W | — | >1 | noise about zero (adiabatic bush) |

Printed summary: `T_shaft` 318.279 K and `T_bush_mean` 316.974 K unchanged;
`core_rise_max` 1.06394 -> 1.06395 K; the steady coupling loop converged in
**15** outer iterations instead of 16, both inside `coupling_tolerance = 1e-6`.
Global liquid mass moved in the 9th significant digit
(2.96795788e-3 -> 2.96795780e-3 kg).

**Reading.** Exactly what the audit predicted: small and boundary-driven. Nothing
in the load-capacity or peak-temperature comparison moves at the precision the
validation figures are drawn at, so **no archived verdict changes.** The two
fields that do move beyond the print floor — the cross-film core temperature rise
and the journal-surface heat flux — are the two computed *from* the wall and
boundary energy balance, which is precisely the balance that was previously being
closed with two inconsistent mobilities. That the discrepancy shows up there and
nowhere else is the strongest available evidence that the diagnosis was right.

### Deliberately not done

1. **`energy.cpp` interior faces** keep their own (already cross-film-aware)
   expressions. Their cross-film branch computes the volumetric flux `q` directly
   rather than `h*u`; re-routing would perturb cross-film results with no
   correctness gain and would have muddied this A/B.
2. **`cfg.reformation_boundary_velocity`** — the Elrod cavitated-end reflood
   closure — still takes `(h_face, mu_face)` with no cross-film branch. It is
   reachable only on the Elrod path under `consistent_boundary_flux`, so it is a
   real but narrower instance of the same defect. Left for a reviewed change.

---

## Test suite

35 ctest cases (was 34; `test_config_keys` is new). Green after every commit and
at the end of the batch. Built with the MSYS2 MinGW64 toolchain, MS-MPI, PETSc;
`ctest` runs at 1, 2, 4 and 6 ranks depending on the case.

## Working-environment note

Another session was operating in `C:\Dev\pancake` concurrently: it held
`pancake.exe` open (blocking the shared build directory) and switched the shared
worktree's branch mid-task, which briefly landed the first commit on `main`. That
was corrected — `main` is back at `79ddbee` with the other session's uncommitted
work intact, and all five commits are on `codex/audit-quickwins`. The rest of the
batch was built and run from a separate `git worktree` so the two sessions could
not collide again. No process was killed and nothing on `afml-r6` was touched.

## Follow-ups this batch surfaced

| item | where |
|---|---|
| A malformed *value* for a known key still throws out of `std::stod` | `config.hpp` parser |
| No test target covers the GUI config load/save policy | `src/gui/` |
| Axial-face TVD correction is still absent (`VANLEER` is theta-only) | `fvm.cpp`, PLAN F4 |
| `reformation_boundary_velocity` is not cross-film aware | `config.hpp` |
| The eleven tests that rely on the `cavitation_model` struct default | listed in §4 |
