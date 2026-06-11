# AGENT_REPORT — GUI Modernization (Phase 0 scan + plan)

Produced per `GUI_MODERNIZATION_DIRECTIVE.md` before any code change.
Updated at the end of the work with acceptance-checklist results (§7 at bottom).

---

## 1. File inventory (GUI-relevant)

| Path | Role |
|---|---|
| `src/gui_win32.cpp` | Entire existing GUI (3727 lines): raw Win32, owner-drawn controls, config editor, solver launcher, VTS preview canvas, console |
| `src/config.hpp` | Header-only `SimulationConfig`: every solver parameter, text `key = value` parser/serializer. Shared with solver; the GUI reuses it read-only |
| `src/pancake_gui.rc` | Resource script: icon + ComCtl32 v6 manifest |
| `src/pancake_gui.manifest` | ComCtl32 v6 + asInvoker manifest (no DPI declaration) |
| `src/resource.h` | Single icon id |
| `src/pancake.ico`, `src/pancake_logo.svg` | App icon assets |
| `src/pancake.rc` | Solver exe icon resource (solver-side; untouched) |
| `CMakeLists.txt` | Single tree: solver (`pancake`), GUI (`pancake_gui`), tests. GUI target already statically linked under MinGW |
| `CMakePresets.json` | One preset: `windows-native-mingw` (Ninja, MSYS2 mingw64) |
| `cmake/deploy_mingw_runtime.cmake` | Copies MinGW runtime DLLs next to the **solver** (GUI no longer uses it) |

## 2. Dependency audit

### Build-time artifacts found in the tree
- `build-windows-mingw/` — gitignored (covered by `.gitignore: build-windows-mingw/`). Not required at runtime.
- `build-gui-check/` — stray scratch build dir at repo root containing `CMakeCache.txt`, `CMakeFiles/`, `tmp/` with downloaded text files. **Not gitignored** → added to `.gitignore` as part of this work. Not required at runtime.
- No stray `.ninja_deps`/`.ninja_log`/`*.obj` outside build dirs.

### Runtime DLL imports of current exes (`objdump -p`, MinGW binutils; dumpbin unavailable — no MSVC on this machine)

`pancake_gui.exe` (current, already `-static`):
`COMCTL32.dll, comdlg32.dll, GDI32.dll, KERNEL32.dll, msvcrt.dll, SHELL32.dll, USER32.dll` — **all ship with Windows**. No vcruntime/msvcp, no MinGW DLLs.

`pancake.exe` (solver — out of scope, ships separately):
`KERNEL32.dll, msvcrt.dll` (Windows) + **external**: `libgcc_s_seh-1.dll, libstdc++-6.dll, msmpi.dll, libpetsc-dmo.dll` (which pulls OpenBLAS/metis/gfortran/winpthread, all deployed beside it by `deploy_mingw_runtime.cmake`). The two-exe folder therefore also contains the solver's DLLs; this is the solver's existing distribution model and is unchanged.

### Runtime file reads of the current GUI
- `config.txt` next to the GUI exe (loaded at startup if present; also via Open/Save dialogs) — config is a *document*, not an asset; this stays.
- Result files under the configured output dir: `*.vts` (XML structured grid), `results.pvd` (handed to ParaView via ShellExecute).
- **No fonts/images/shaders read from disk.** No registry, no settings persistence at all (window layout, solver path, MPI ranks are lost on exit).

## 3. Feature map (existing GUI) — KEEP / REWORK / DROP

| Feature | Verdict | Why |
|---|---|---|
| Parameter editor, ~80 fields grouped by physics domain | REWORK | Keep grouping + common/advanced split; rebuild in ImGui with per-field units and inline validation |
| Live validation (red fields, run disabled, error list) | KEEP (re-implement) | Already matches directive §4.3 |
| SI/Engineering unit toggle | REWORK | Directive §4.4 wants one consistent system → SI display with units on every field; engineering equivalents shown in tooltips |
| Presets (Default / High ecc / Misaligned / Oil hole) | KEEP | Cheap, useful as "Reset to defaults" + examples |
| Raw config.txt text editor with Sync/Apply | KEEP | Engineers hand-edit configs; low cost in ImGui |
| Open/Save/Save As/Reload config | KEEP | Core workflow |
| Run/Stop + MPI ranks + mpiexec discovery | KEEP (re-implement on worker threads) | Required; old version polled on a 400 ms UI timer → moves to reader thread |
| Progress from `t=` stdout parsing | REWORK | Keep defensive parse; add `Step N` + `diagnostics.csv` tail for residuals; structured-format proposal in §6 below |
| Console with ANSI color parsing | KEEP | Maps directly to log viewer requirement |
| VTS heat-map preview canvas (GDI) | REWORK | Becomes ImPlot heatmap with pan/zoom/hover, cavitation overlay, axes+units, PNG export |
| Step selector (Latest / step N, prev/next) | KEEP | Needed for transient runs |
| Open results folder / Open in ParaView | KEEP | One-line ShellExecute, high value |
| Summary card with derived quantities | REWORK | Becomes key-results panel (min h, peak p, attitude angle, friction power…) with copy-value |
| Owner-drawn Win32 theme/splitters | DROP | Replaced by ImGui docking |
| Run history / comparison | (absent) | NEW — required by §4.9/§5.1 |
| CSV/PNG export | (absent) | NEW — required by §4.8 |
| Settings persistence | (absent) | NEW — `%LOCALAPPDATA%\PancakeGui\` |
| DPI awareness | (absent) | NEW — per-monitor v2 manifest + font rebuild on `WM_DPICHANGED` |

## 4. Solver interface (discovered contract)

- **Invocation:** `pancake.exe [-c|--config] <config_file>`; with MPI: `mpiexec.exe -n <N> pancake.exe -c <config>`. Runs fine without mpiexec for 1 rank (MS-MPI singleton init; this is what the current GUI does).
- **Input:** plain-text `key = value` config with `#` comments; full schema lives in `src/config.hpp` (~90 keys: geometry, grid, time, physics, thermal, fluid-property models, bearing motion, BCs, cavitation, inlets, output). The GUI reuses `config.hpp` for parse/serialize — solver source untouched.
- **stdout:** ANSI-colored lines via rank-0 logger. Per step: `Step <n> t=<t>` (+ in DETAILED verbosity: `| h_min=… | T=[…] K | F=(…) N | M=… Nm | outer=… flips=… | r_l=… [r_gas=…]`), `[saved]` marker on write steps; `Simulation finished successfully.` on success; warnings/errors in yellow/red.
- **Exit codes:** 0 success, 1 config validation error, 3 output-dir setup failure, other = PETSc/MPI runtime error.
- **Result files:** `<output_dir>/results.pvd`, `solution_<step>.pvts` + `processor<r>/solution_<step>_<r>.vts` (XML structured grid, cell-centered fields: `p, h, film_content, T, rho, mu, F/Fp/Fv vectors, …`); optional `flat/` unwrapped copies; `diagnostics.csv` (step, time, outer_iters, residuals, masses, cavitated_fraction — written every `diagnostics_interval` steps).
- **Cancellation:** the solver has **no cooperative shutdown channel**. The GUI will TerminateProcess with a visible warning (see proposal in §6).

## 5. Proposed plan (ordered)

1. **Phase 0 report** (this file) — commit before code changes. ✅
2. **Scaffolding:** `third_party/` (ImGui docking + ImPlot, pinned via FetchContent; vendored fonts + stb_image_write), new `src/gui/` target with Win32+D3D11 ImGui backends, CMake presets (`release` = MSVC static CRT per §6.2 of the directive; `release-mingw` = fully static MinGW, used for verification on this machine — see Deviations). Minimal window builds + launches.
3. **App shell:** `AppState`, dark/neutral theme with single accent, embedded fonts (sans + mono), per-monitor-v2 DPI with font rebuild, settings + `imgui.ini` in `%LOCALAPPDATA%\PancakeGui\`, status bar (solver state / elapsed / last residual).
4. **Case setup panel:** parameter schema table (label, key, unit, default, bounds, common/advanced), generated widgets, live validation with red highlight + tooltip, Run disabled with reason, presets + reset-to-defaults, raw-config tab, Open/Save/Reload.
5. **Solver runner:** worker thread + CreateProcess pipes, log ring buffer with ANSI parsing, progress (`t=`/`Step`/finished), `diagnostics.csv` tailing for residuals, cancellation, solver discovery (GUI-exe dir → settings override; clear "not found" state), MPI ranks.
6. **Results:** worker-thread VTS + diagnostics parsing (malformed ⇒ readable error), ImPlot pressure heatmap with cavitation region (film_content < 1) overlay, h/p mid-plane profiles, convergence plot (residual vs step, log scale), key-scalar summary with copy buttons.
7. **History + export:** in-session run table (inputs + key results, sortable), overlay of ≥2 runs on profile plots, CSV export for tables/curves, PNG export of plots (backbuffer crop via stb_image_write).
8. **Packaging + docs:** version/icon `.rc`, `BUILD.md`, `docs/GUI_ARCHITECTURE.md`, remove `src/gui_win32.cpp` from the build, gitignore fixes, acceptance checklist run (§7 below).

## 6. Proposals for the solver maintainer (NOT implemented — solver code untouched)

- **Structured progress line:** emit `PROGRESS step=<n> t=<t> end_t=<T> res=<r>` once per step (machine-first, in addition to the human line). The GUI parses the existing format defensively but it is whitespace/locale-fragile.
- **Cooperative cancellation:** poll for a sentinel file `<config>.stop` (or handle `CTRL_BREAK_EVENT`) and exit with a distinct code after writing the current step. Until then the GUI must TerminateProcess, which can leave a partially written VTS file.

## 7. Decisions taken under latitude (§5.2) and deviations

- **Toolchain deviation (environment, not choice):** this machine has **no MSVC installation** (`vswhere` finds nothing; `C:\Program Files\Microsoft Visual Studio\18` is empty; no `cl.exe`). Installing VS Build Tools requires admin + multi-GB download — the directive's own user profile ("no admin rights") argues against depending on it. Resolution: the CMake project fully supports MSVC (preset `release`, `CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>"`, `/W4 /WX` on GUI sources) **and** MinGW (`release-mingw`, `-static`, `-Wall -Wextra -Werror`). Verification in this session uses the MinGW build; its import table meets the same acceptance bar (Windows-shipped DLLs only, no vcruntime/msvcp, no MinGW DLLs — already proven by the current GUI). When Build Tools are present, `cmake --preset release` produces the MSVC binary with zero source changes.
- **Solver exe name:** the solver is `pancake.exe`, not `solver.exe`. Discovery order: `pancake.exe` next to the GUI exe → user override in Settings (persisted). `solver.exe` is also accepted as a fallback name to honor the directive's wording.
- **App identity:** target/exe stays `pancake_gui.exe` (product name "Pancake GUI") — same two-file folder users already know.
- **Units:** SI display (m, Pa, rad/s, K, s) matching the config file exactly — no hidden conversions to corrupt round-trips; every field labeled; tooltips show common engineering equivalents (µm, MPa, rpm).
- **Layout:** single window, docking branch: left Case Setup, center tabbed Results/Convergence/Log, right Run History, persistent bottom status bar. Default layout built programmatically; user rearrangements persist via `imgui.ini` in `%LOCALAPPDATA%`.
- **Fonts:** vendored in `third_party/fonts/` (OFL/Apache licensed sans + mono), converted to C arrays at build time by a CMake script — no network needed at build for fonts, zero runtime file reads.
- **History persistence:** in-session only (directive leaves it optional); CSV export of the history table covers cross-session needs.
- **Code signing:** out of scope per directive; recommended as a future step to reduce SmartScreen friction.

## 8. Acceptance checklist (§7 of the directive) — results

Filled in at the end of the work; see bottom-of-file UPDATE section.
