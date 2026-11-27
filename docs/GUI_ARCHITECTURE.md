# Pancake GUI architecture

The GUI (`src/gui/`) is a Dear ImGui (docking) + ImPlot application on the
Win32 + Direct3D 11 backends, built for compressor-bearing engineers
(`docs/GUI_LGE_PLAN.md`, `docs/GUI_GLOSSARY.md`). It drives the solver
`pancake.exe` strictly as a subprocess. The shared sources with the solver are
the header-only `src/config.hpp` (`SimulationConfig` parse, write and
`validate()`) and `src/property_table_io.cpp` (the table loader), so the GUI and
the solver agree on the config format, on what a valid case is, and on the
fluid tables.

## Data flow

```
panels edit CaseModel (engineer inputs, SI)
   -> to_config(case, pair, loaded config, case as loaded)   writes only changed key groups
   -> SimulationConfig::validate() + GUI checks              -> status lights (green/yellow/red)
   -> estimates (S, p_sat, margin, onset verdict, run time)
Save: patch_config_text(original text, ...)                  untouched lines stay byte-identical
Run:  save, then pancake.exe -c <case.txt> (mpiexec -n N ... for N > 1 cores)
```

The config file is the only case format. GUI-only settings (fluid pair id,
"saturated", mesh level, cores) are `#@gui key = value` comment lines, which the
solver ignores.

## File map

| File | Responsibility |
|---|---|
| `case_model.hpp/.cpp` | Engineer-level `CaseModel`; `to_config` / `from_config`; model selection (FBNS plateau, COUPLED gas transport, thermal forms); settings the GUI does not show. Toolkit-free, unit-tested |
| `config_text_patch.hpp/.cpp` | Text-preserving save of a loaded config; new-case text with `#@gui` lines. Toolkit-free |
| `fluid_library.hpp/.cpp` | Oil-refrigerant pairs in `fluids/` (`pair.txt` + solver solution table); optional gas block; property queries with range flags. Toolkit-free |
| `estimates.hpp/.cpp` | Sommerfeld number, saturation pressure, pressure margin, onset chart and verdict, run-time model (plan 0.2). Toolkit-free |
| `key_results.hpp/.cpp` | Key results of a solved step, design checks (gas release, Lambda), result fields by solver form (FBNS never shows gas fields). Toolkit-free |
| `case_status.hpp/.cpp` | Workflow-node lights and reasons, solver validation mapped to nodes, solver-form text. Toolkit-free |
| `angle_reference.hpp` | Angle zero (fixed, maximum film, minimum film, feed, opposite the feed) and display transforms |
| `units.hpp/.cpp` | Display units per field (SI stored) |
| `ui_widgets.hpp/.cpp` | Shared form rows, unit inputs, status dots, risk colours, tooltips |
| `gui_app.hpp/.cpp` | `GuiApp` state and shell: case files, fluid library, menus, dock layout, status bar, log, modals, screenshot series |
| `run_control.cpp` | Run queue (single, comparison, batch, two-step thermal route), solver launch and stage text, key results parsed off the UI thread, results loading, thermal-transient confirmation, scripted `--run` |
| `workflow_tree.cpp` | Run bar, workflow tree with lights, estimates box, panel window |
| `panel_*.cpp` | Geometry, fluid, operating conditions, boundary zones, mesh, solution panels |
| `view_geometry.cpp`, `view_mesh.cpp` | Live views: bearing section to scale (clearance exaggerated), unwrapped film with the mesh |
| `dialog_fluid_library.cpp` | Fluid library: list, copy, create, import, table editing with paste, plots |
| `onset_plot.cpp` | Onset chart widget |
| `results_window.cpp`, `results_tabs.cpp`, `result_data.hpp` | Results: key-results card, field maps, polar, profiles, transient histories, comparison, exports |
| `run_monitor.cpp` | Live monitors during a run, stage, progress, command line |
| `batch.cpp` | Sweeps and case lists, scripts, live table and curves |
| `case_report.cpp` | One-page HTML case report (inline SVG figures) |
| `solver_runner.hpp/.cpp` | Subprocess management, stdout and diagnostics tails |
| `results.hpp/.cpp` | Async, crash-proof VTK reader |
| `exports.hpp/.cpp`, `settings.hpp/.cpp`, `theme.hpp/.cpp`, `win_util.hpp`, `main_gui.cpp` | CSV/PNG export and dialogs, persisted settings, theme and fonts, string helpers, platform layer |

## Threading model

The UI thread only draws; it never blocks on the solver or on large files.
The solver reader thread and the diagnostics tail feed mutex-guarded queues;
`ResultsLoader` parses VTK on a worker; key results of a finished run are
parsed with `std::async`; the transient monitor tails `motion.csv` on a worker.

## Solver contract used

- Launch: `pancake.exe -c <case.txt>` (working directory = case folder), or
  `mpiexec -n <N> pancake.exe -c <case.txt>`.
- Output: `<case>_results/` with `flat/results.pvd` (unwrapped film, read by the
  GUI), `diagnostics.csv`, `equilibrium_history.csv` (steady load balance),
  `motion.csv` (transient), `iterations/NNNNNN/` snapshots every 10 outer
  iterations (steady live views).

## Command line

```
pancake_gui.exe [--open <case.txt>] [--run] [--fluids <dir>] [--screenshots <dir>] [--size WxH]
```

`--run` runs the opened case through the GUI's own Run path and writes
`gui_key_results.json` into its results folder. `--screenshots` saves every
screen as PNG (user-guide figures) and exits. `--fluids` sets the fluid
directory for this session only.

## Packaging

`tools/make_release.ps1` collects `pancake.exe`, `pancake_gui.exe`, the runtime
DLLs, `fluids/` and `examples/` into one folder; Microsoft MPI is a separate
install. `data/` and `tests/data/lge/` (LGE test data) are not shipped. See
`BUILD.md`.
