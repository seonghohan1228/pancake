# Pancake GUI architecture

The GUI (`src/gui/`) is a Dear ImGui (docking) + ImPlot application on the
Win32 + Direct3D 11 backends. It drives the solver `pancake.exe` strictly as a
subprocess and never links solver-algorithm code; the shared sources are the
header-only `src/config.hpp` (`SimulationConfig` parse/serialize/validate) and
`src/property_table_io.cpp` (the `*_table_file` CSV loader), which guarantee the
GUI and solver always agree on the config format, on what constitutes a valid
case, and on how the property tables are read.

## File map

| File | Responsibility |
|---|---|
| `main_gui.cpp` | Platform layer: window, D3D11 device/swapchain, ImGui/ImPlot lifecycle, per-monitor-v2 DPI (font/style rebuild on `WM_DPICHANGED`), frame loop, backbuffer-region PNG capture |
| `gui_app.hpp/.cpp` | `GuiApp` - the single application-state object (no globals): lifecycle, dock layout, menu, status bar, run control, config file handling, CSV/PNG export glue |
| `panels_setup.cpp` | Case Setup window: run bar (Run/Cancel + MPI process count), schema-generated parameter form ordered by edit frequency with per-field display-unit selectors, enum/model selectors with per-option hover descriptions, inlet editor. Cavitation panel exposes the `cavitation_threshold` (SCALAR_PCAV / LOCAL_PSAT) and `gas_pressure_coupling` (NONE / VOID_COUPLED) selectors. Selecting a dissolved-gas species fills its solubility/gas/viscosity coefficient block (`apply_dissolved_gas_species_preset` — calibrated R-290-in-PZ68S, or generic air-in-oil); model selectors stay user-controlled. Each TABLE property model shows a `*_table_file` path + Browse field whose CSV is loaded via `GuiApp::reload_property_tables` (the linked `property_table_io.cpp`) for a live load check |
| `panels_results.cpp` | Results (equal-aspect field heatmap in cell space with cavitation/inlet overlays, Fit button, a colormap selector (Viridis / blue-red diverging) and an optional manual color range, profiles, key-results summary), Monitors (selectable physical traces: h_min, T range, force, torque, cavitated %, residual - stacked subplots in real units, parsed live from DETAILED step lines + diagnostics.csv) and Solver Log windows |
| `schema.hpp/.cpp` | Data-driven parameter descriptors (label, unit family, help, bounds, group, common/advanced) over `SimulationConfig` members; affine display-unit families (`UnitOption` = scale + offset; config stays solver-native): Pa/kPa/MPa/bar, Pa/MPa/GPa (bulk modulus), deg/rad, rad/s-rpm, m/mm/um, s/ms, m/m-mrad (tilt), K/degC, kg/kg-% (dissolved gas), N/m-kN/m-MN/m, N.s/m-kN.s/m; `validate_config()` merges the solver's own `SimulationConfig::validate()` so Run is disabled exactly when the solver would reject the case. Bearing motion is exposed as isotropic support stiffness/damping with the initial position taken from the geometry (attitude angle) and zero initial velocity (`GuiApp::normalize_motion_config`) |
| `solver_runner.hpp/.cpp` | Subprocess management: `CreateProcessW` with redirected stdout/stderr, a reader thread (ANSI-color parse, progress from step lines), a diagnostics.csv tail thread (residuals), cancellation via `TerminateProcess` |
| `results.hpp/.cpp` | Async, crash-proof parser for the solver's VTK XML output (`results.pvd` -> `.pvts` -> per-rank `.vts` pieces, assembled into a global cell grid); malformed files become readable error strings |
| `exports.hpp/.cpp` | CSV writers (column/row oriented), PNG encode (vendored stb_image_write), Win32 open/save dialogs |
| `settings.hpp/.cpp` | `%LOCALAPPDATA%\PancakeGui\` paths (settings.txt, imgui.ini) and persisted user settings |
| `theme.hpp/.cpp` | Neutral-gray theme with a single blue accent; embedded Roboto / JetBrains Mono fonts rebuilt per DPI scale |
| `win_util.hpp` | UTF-8 <-> UTF-16 and exe-path helpers |

## Threading model

The UI thread only draws; it never blocks on the solver or on files.

- **Solver reader thread** (`SolverRunner::reader_thread_main`): blocking
  `ReadFile` on the stdout/stderr pipe, line assembly, ANSI color parsing,
  progress extraction (`Step N t=...`, `finished successfully`), exit-code
  collection. Output lands in a mutex-guarded queue.
- **Diagnostics tail thread**: polls `diagnostics.csv` every ~300 ms while the
  solver runs and parses complete rows (residuals, cavitated fraction).
- **Results worker** (`ResultsLoader`): one worker at a time scans the output
  directory / parses a step grid; the UI polls a revision counter and swaps a
  `shared_ptr<const ResultGrid>` snapshot.
- The UI drains the queues once per frame (`drain_log` / `drain_diagnostics`)
  and copies a small `SolverProgress` struct under the mutex.

Cancellation: `TerminateProcess` (the solver has no cooperative shutdown
channel yet - see the proposal in `AGENT_REPORT.md` §6); the log warns that
the last output file may be incomplete.

## Solver contract used

- Launch: `pancake.exe -c <config>` (working dir = config dir), or
  `mpiexec -n <N> pancake.exe -c <config>` for MPI ranks > 1.
- Progress: indented `Step <n> t=<t>` stdout lines vs. `end_t`; residuals from
  `diagnostics.csv`.
- Results: `<output_dir>/results.pvd` + `solution_<step>.pvts` +
  `processor<r>/*.vts` (ASCII XML, cell-centered fields incl. `p`, `h`,
  `film_content`, `T`, vector groups `F/Fp/Fv/Fext/xB`, and the constant
  resultant fields `friction_torque`, `bearing_attitude_angle`).
- Exit codes: 0 ok, 1 config validation, 3 output-dir setup failure.

## Packaging

Static CRT (MSVC `/MT` via `CMAKE_MSVC_RUNTIME_LIBRARY`, MinGW `-static`),
fonts embedded as byte arrays generated by `cmake/embed_binary.cmake`,
icon + version info + per-monitor-v2 manifest in `src/gui/pancake_gui.rc`.
The exe reads zero asset files at startup and writes only under
`%LOCALAPPDATA%\PancakeGui\`. See `BUILD.md`.
