# GUI for LGE bearing engineers: final report (2026-11-26)

Branch `lge-gui` (worktree `C:\Dev\pancake\worktrees\lge-gui`), not pushed. Plan: `docs/GUI_LGE_PLAN.md` (rev. 2 + rev. 3, approved). Glossary: `docs/GUI_GLOSSARY.md`. Screens: `docs/gui_screens/` (1600 × 1000 and 1280 × 800).

## What is where

| Class | Contents |
|---|---|
| **Basic** | D, L, c; feed type, angle, size; analysis type; speed; load and direction (on the shaft) or load table; ambient p and T; fluid pair; content or "saturated"; roughness; dissolved gas on/off; thermal on/off; comparison mode; mesh level; CPU cores (default 1) |
| **Expanded** | Tilt; hole axial position; fixed shaft position; gas model (plateau / gas transport); tolerances with defaults and meaning; time stepping; exact cells; cross-film layers; zone values (feed overpressure and temperature, housing, shaft) |
| **Auto** (read-only) | Solver form (FBNS / COUPLED, load balance, thermal form), cavity shear (FBNS: no shear from the bubble, solver key STRIATED), properties at the ambient state, end pressures and temperatures, output settings |
| **Removed** (keys kept on save) | Elrod and legacy cavitation forms and constants, contact / mixed lubrication, raw support stiffness and damping, PETSc options, gas EOS internals, debug outputs, convection schemes |

Screen map: Workflow tree (lights, estimates, onset chart) | Case setup panel + Live view (geometry, mesh) | Results, Run monitor, Sweep/batch, Solver log; dialogs: fluid library, boundary zones, thermal-transient confirmation, settings.

## Acceptance tests

| # | Test | Result |
|---|---|---|
| 1 | Round trip of 7 deck configs (O4, X4, optical rig, compressor GC+VC, VC, transient, P5g) | **Pass**: 201-203 solver keys each equal; original lines byte-identical (CRLF, afml-r6 paths) |
| 2 | Reproduction from GUI-written configs | **Pass**. Compressor GC+VC 3450 N: h_min 1.6876 µm (deck 1.6877), cavitated area 57.09 % (equal), torque −0.004 %. VC only: identical. Transient (3 rev, afml-r6, 24 ranks): third-revolution h_min 1.6417-1.7624 µm (deck 1.6406-1.7619). Tester 2610 N in the plateau form: h_min 3.323 µm vs P5g 3.242 (+2.5 %, model form R1), cavitated area 11.90 vs 11.83 %, p_max 89.8 vs 90.6 bar |
| 2' | STRIATED default (R3) | Compressor: h_min +0.02 %, area equal, **torque −18.2 %** vs the deck (FULL_FILM) |
| 3 | GUI run = terminal run (same config, 8 ranks) | **Pass**: h_min, area, torque, p_max equal to ~1e-10 |
| 4 | Fluid: pair from CSV, out-of-range warning, shipped examples | **Pass**: pair created and saved through the library API validates; content above the table turns Fluid yellow; example steady and transient run (10 s, 207 s on 4 cores). The dialog's click path was not driven (no input automation) |
| 5 | Estimates, tester | **Pass**: S = 0.0159 (deck 0.016), p_sat 9.65 bar (deck 9.7) |
| 6 | Angle references, compressor | **Pass**: feed zero puts the feed at 0 and h_min at θ_in = −160.25° (deck −160.25°); h_max zero puts h_min at 180° |
| 7 | New case 1 core; tolerances change iterations | **Pass**: 1 core; load-balance tolerance 1e-2 → 1e-4: 7 → 8 iterations (Newton converges quadratically) |
| 8 | Clean build; release folder only | **Pass**: 70/70 tests; GUI started from the release folder with a clean PATH and ran the example. Folder: exes, 12 runtime DLLs, `fluids/`, `examples/`, README (85 MB); MS-MPI installed separately |
| 9 | Screenshots | **Done**: 18 PNG; overlaps checked and fixed at both sizes |
| 10 | Run-time estimate within ×2 (this machine: 4 cores, 8 threads) | Isothermal steady **pass** (130 s est., 147 s meas.). Thermal steady **miss**: 486 s est., 202 s meas. (×2.4; the case converged in 36 outer iterations, near the low end of 0.8-8 core-h). Easy cases run faster than the typical value (example steady: ×4) |
| 11 | FBNS never shows a nonzero gas field | **Pass**: unit test; the real G_eps0.9300 output has alpha_gas = 0, dissolved_gas constant, no density, none offered |
| 12 | Thermal transient small test | **Pass**: solver runs energy + transient in the simplified forms (0.18-0.26 s/step on 180 × 48, 4 processes); the deck thermal form is steady-only (solver's reason) |

Solver changes (approved D2, D3): `STEADY_EQUILIBRIUM` accepts FBNS with dissolved gas (matches driver balances within 0.02 %); the config writer keeps the steady coupling keys.

## Open items

1. Thermal-steady run-time typical value: one local case is ×2.4 low; calibrate with more cases or the O-run timings.
2. Proposed, not built: `gas_molar_mass` key for a free ideal gas; unknown enum values as parse errors (X4's `HOMOGENEOUS` ran as FULL_FILM).
3. Not visually verified: Polar, Profiles, Transient and Comparison tabs with data; the fluid dialog with the gas block on; a batch end to end.
4. The deck torque and power-loss numbers were FULL_FILM; GUI numbers with STRIATED are lower (−18 % at 57 % cavitated area). Decide how to present this in the user guide.
5. afml-r6 clone left at `2e6b303` (was `55808eb`), `build-linux` rebuilt there.
6. Next: the user guide (KPI) from `docs/gui_screens/` and the glossary.
