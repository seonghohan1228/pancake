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

The post-build runtime deploy step copies MinGW/PETSc DLLs beside the
executables and tests. This includes `libpetsc-dmo.dll`, `libstdc++-6.dll`,
`libopenblas.dll`, and related dependencies.

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

The GUI uses a persistent backend header plus one main workspace:

- Top-right header: solver executable, MPI ranks, Run, and Stop stay visible
  on every tab.
- `Workspace` left side: result preview above the process log. This side
  receives most extra width/height when the window is maximized.
- `Workspace` right side: scrollable parameter inspector with one item per row
  and units beside each value.
- Horizontal splitter: drag it to resize the result preview and log pane.
- Timestep selector: choose `Latest` or a specific VTS step, or use previous
  and next buttons to flip through written result files.
- Colorbar and axes: show the current field, units, min/max-inclusive ticks,
  and preview coordinates beside the heatmap.
- Multi-rank assembly: each preview timestep is assembled from all matching
  `processor*` VTS files.
- `load_angle_deg`: visualization reference angle measured from the positive y
  axis; the preview rotates the circumferential display relative to this load
  reference.
- Parameter units: per-value unit selectors affect GUI input while saved config
  values remain in solver-native units.
- Axial boundaries: for `ELROD_ADAMS`, the pressure rows are not used as
  theta boundary values. `DIRICHLET` and `INLET_OUTLET` inflow use the separate
  dimensionless film-content rows.
- `film_content`: visualization output is the raw Elrod/JFO universal variable.
  Cavitated cells are below one, while compressed full-film cells can exceed
  one and are not clamped in the VTK files or GUI preview.
- Process log: ANSI color codes from the backend are rendered as color instead
  of raw escape fragments.
- Inlet editor: `GROOVE` and `CIRCULAR` modes show only the relevant fields;
  the unused config form is written as a comment.
- Visual styles: `src/pancake_gui.rc` embeds a ComCtl32 v6 manifest so the GUI
  uses modern themed controls under MinGW. A light palette plus `WM_CTLCOLOR*`
  handlers and a tab-control background subclass keep every surface consistent
  instead of the default system grey. `enable_language(RC)` is required in the
  CMake configure step for the GUI target.
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
