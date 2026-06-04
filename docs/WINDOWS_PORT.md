# Windows Development Notes

The current Windows workflow is native MSYS2 MINGW64. It builds both the solver
backend (`pancake.exe`) and the GUI (`pancake_gui.exe`) from the same CMake
preset.

Visual Studio is not required for the current workflow. Older Visual Studio
GUI-only presets were removed to keep the build path clear and reproducible.

## Required Installations

1. Microsoft MPI Runtime and SDK
   - Install both `msmpisetup.exe` and `msmpisdk.msi`.
   - Download: <https://www.microsoft.com/en-us/download/details.aspx?id=105289>

2. MSYS2
   - Download: <https://www.msys2.org/>
   - Use the **MSYS2 MINGW64** shell.

3. MSYS2 packages

   ```bash
   pacman -Syu
   pacman -S --needed git base-devel mingw-w64-x86_64-toolchain \
     mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja \
     mingw-w64-x86_64-pkgconf mingw-w64-x86_64-msmpi \
     mingw-w64-x86_64-petsc
   ```

## Build

From **MSYS2 MINGW64**:

```bash
cd /c/Dev/pancake
cmake --fresh --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

The same commands also work from PowerShell in `C:\Dev\pancake`. The preset
prepends `C:\msys64\mingw64\bin` to `PATH`, which keeps MinGW compiler helpers
and DLLs discoverable.

Expected CMake output:

```text
PETSc pkg-config module: petsc-dmo
```

Do not use the serial PETSc package `petsc-sso`; it collides with MPI headers
and causes errors around `MPI_Wtick`, `MPI_Pcontrol`, and
`PETSC_FUNCTION_NAME_CXX`.

## Outputs

```text
C:\Dev\pancake\build-windows-mingw\pancake.exe
C:\Dev\pancake\build-windows-mingw\pancake_gui.exe
```

The post-build runtime deploy step copies MinGW/PETSc DLLs beside `pancake.exe`
and the tests. This includes `libpetsc-dmo.dll`, `libstdc++-6.dll`,
`libopenblas.dll`, and related dependencies. `pancake_gui.exe` is statically
linked (`-static` under MinGW, static CRT under MSVC), so it ships as a single
self-contained executable that needs no runtime DLLs and is not part of the
deploy step.

## Run

```bash
cd /c/Dev/pancake/build-windows-mingw
./pancake.exe -c config.txt
mpiexec -n 2 ./pancake.exe -c config.txt
./pancake_gui.exe
```

The GUI expects `pancake.exe` in the same directory as `pancake_gui.exe` or in
`build-windows-mingw/`.

## GUI Layout

The GUI is a single tab-less workspace: an action bar on top, a resizable input
rail on the left with a pinned live-summary card beneath it, a result viewport
in the center, and a resizable console at the bottom.

- Action bar: owner-drawn `Run` / `Stop` buttons, MPI `ranks`, a run-state chip
  (`Idle` / `Running NN% · Ns` with a progress bar / `Done` / `Failed`), and the
  file actions `Open…`, `Save`, `Save As…`, `Reload`, `Reset`, `Solver…`.
- Input rail: a search filter, an `Advanced` toggle, a preset picker, and a
  unit-system selector, above collapsible sections (Geometry, Operating/Fluid,
  Mesh & Time, Cavitation, Axial Boundaries, Inlet, Output, Raw config.txt).
  Section headers and field labels are painted by the rail's `WM_PAINT`; only
  interactive controls are child windows. There are **no `BS_GROUPBOX` frames**
  around live combos, which removes the z-order/paint fights of the old layout.
- Live summary card: derived journal-bearing quantities (`ε = e/c`,
  `h_min = c(1−ε)`, `h_max`, `c/R`, `L/D`, `U = ωR`, previewed field mean/max)
  recomputed on every edit, plus a validation status line.
- Validation: malformed or physically impossible inputs (`e ≥ c`, non-positive
  grid/clearance, time-step ordering, empty output names) paint the field red
  via `WM_CTLCOLOREDIT` and disable `Run` until resolved.
- Configuration files: defaults to `config.txt` beside the executable; `Open…`
  and `Save As…` use `GetOpenFileNameW` / `GetSaveFileNameW` (comdlg32) for
  arbitrary `.txt` cases. The loaded file name and a dirty `*` marker show in
  the title bar.
- Units: a global SI ↔ engineering selector sets defaults; per-field override
  dropdowns remain only where multiple units are realistic. Saved config values
  remain in solver-native units.
- Result viewport: field/timestep selectors, previous/next/refresh, a heatmap
  with axes and a unit-labeled colorbar (min/max-inclusive ticks). Each preview
  timestep is assembled from all matching `processor*` VTS files.
- Progress: the run-state chip and progress bar are driven by parsing
  `t=<value>` from solver stdout against `end_t`; `finished successfully` marks
  completion.
- `load_angle_deg`: visualization reference angle from the positive y axis; the
  preview rotates the circumferential display relative to this load reference.
- Axial boundaries: for `ELROD_ADAMS`, `DIRICHLET` / `INLET_OUTLET` inflow use
  the dimensionless film-content rows rather than the pressure rows.
- `film_content`: raw Elrod/JFO universal variable; cavitated cells are below
  one and compressed cells can exceed one, unclamped in VTK and the preview.
- Process log: ANSI color codes from the backend render as color.
- Inlet editor: `GROOVE` / `CIRCULAR` modes show only the relevant fields; the
  unused config form is written as a comment.
- Visual styles: `src/pancake_gui.rc` embeds a ComCtl32 v6 manifest so the GUI
  uses modern themed controls under MinGW; `enable_language(RC)` is required in
  the CMake configure step for the GUI target.
- Application logo: `src/pancake.ico` is embedded in both `pancake.exe` and
  `pancake_gui.exe`; `src/pancake_logo.svg` is the editable source-style asset.

## Validation

```bash
cd /c/Dev/pancake
ctest --test-dir build-windows-mingw --output-on-failure
```

Known current status: `test_elrod_a1` fails a numerical pressure-regression
threshold. The native build, runtime DLL deployment, smoke run, and other tests
are working.

## Portability Notes

Linux remains a supported target from the same source tree. Keep platform build
directories separate (`build-linux/`, `build-windows-mingw/`) and avoid
platform-specific code in solver/core files unless it is guarded and documented.
