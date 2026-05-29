# pancake

`pancake` is a 2D finite-volume solver for compressible thin-film lubrication
in a journal bearing. It solves the Reynolds equation on a cylindrical
`(theta, z)` mesh with MPI parallelism and writes VTK/PVD results for ParaView.

The current Windows development workflow uses **MSYS2 MINGW64**. Visual Studio
is not required to build or run the native solver or GUI in the current code.

## Current Windows Workflow

This workflow builds:

```text
C:\Dev\pancake\build-windows-mingw\pancake.exe
C:\Dev\pancake\build-windows-mingw\pancake_gui.exe
```

The build also copies the required MinGW/PETSc runtime DLLs beside the
executables, so File Explorer launches work after Microsoft MPI is installed.
Both Windows executables embed the shared logo from `src/pancake.ico`.

### 1. Install Prerequisites

Install Microsoft MPI Runtime and SDK:

- Download: <https://www.microsoft.com/en-us/download/details.aspx?id=105289>
- Install both `msmpisetup.exe` and `msmpisdk.msi`.

Install MSYS2:

- Download: <https://www.msys2.org/>
- Open **MSYS2 MINGW64** from the Start menu.

Update MSYS2:

```bash
pacman -Syu
```

If MSYS2 asks you to close the terminal, close it, reopen **MSYS2 MINGW64**,
then run:

```bash
pacman -Syu
```

Install the build packages:

```bash
pacman -S --needed git base-devel mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja \
  mingw-w64-x86_64-pkgconf mingw-w64-x86_64-msmpi \
  mingw-w64-x86_64-petsc
```

### 2. Clone

From **MSYS2 MINGW64**:

```bash
mkdir -p /c/Dev
cd /c/Dev
git clone https://github.com/seonghohan1228/pancake.git
cd /c/Dev/pancake
```

If the repository already exists:

```bash
cd /c/Dev/pancake
```

### 3. Configure And Build

From **MSYS2 MINGW64**:

```bash
cd /c/Dev/pancake
cmake --fresh --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

From PowerShell, use the same preset from `C:\Dev\pancake`:

```powershell
cmake --fresh --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

The preset prepends `C:\msys64\mingw64\bin` to `PATH`, which prevents the
silent `cc1plus.exe` missing-DLL failure that can happen in a plain PowerShell
session.

Expected configure output includes:

```text
PETSc pkg-config module: petsc-dmo
```

If it mentions `petsc-sso`, CMake found the serial PETSc package. Re-run the
fresh configure command above from **MSYS2 MINGW64**.

### 4. Run

Single process:

```bash
cd /c/Dev/pancake/build-windows-mingw
./pancake.exe -c config.txt
```

MPI run:

```bash
cd /c/Dev/pancake/build-windows-mingw
mpiexec -n 2 ./pancake.exe -c config.txt
```

GUI:

```bash
cd /c/Dev/pancake/build-windows-mingw
./pancake_gui.exe
```

You can also double-click `pancake_gui.exe` or `pancake.exe` in File Explorer
from `C:\Dev\pancake\build-windows-mingw`.

### 5. Test

```bash
cd /c/Dev/pancake
ctest --test-dir build-windows-mingw --output-on-failure
```

Current known status: the native build succeeds, the GUI smoke launch succeeds,
and all tests pass except `test_elrod_a1`, which currently fails a numerical
regression threshold rather than a build/runtime dependency.

## GUI Workflow

`pancake_gui.exe` reads and writes `config.txt` beside the executable. Solver
backend controls are always visible in the top-right header, so you can edit
parameters, output selections, or raw config without leaving the main workspace.
The `Workspace` tab places the result preview and process log on the left and a
scrollable parameter inspector on the right. Drag the horizontal splitter to
rebalance result view and log height. Maximizing the window gives the extra
space to the result preview while preserving the parameter column.

Use the tabs to choose which result fields to save, run or stop the backend
solver, inspect colorized process logs, edit raw config, and preview VTK
results. The preview has field and timestep selectors, previous/next step
buttons, axes, and a unit-labeled colorbar with min/max-inclusive ticks.
Multi-rank result previews are assembled from all `processor*` VTS files for
the selected timestep.

Numeric edit boxes preserve the text you typed. For example, `0.0010` remains
`0.0010` after saving.

The Windows application logo is embedded from `src/pancake.ico` in both
`pancake.exe` and `pancake_gui.exe`. The editable source-style logo is
`src/pancake_logo.svg`.

The parameter inspector puts one item per row with unit selectors directly
beside values: lengths (`m`, `cm`, `mm`, `um`), counts (`cells`), times (`s`,
`ms`), angles (`deg`, `rad`), tilt (`m/m`, `deg`, `rad`), pressure (`Pa`,
`kPa`, `MPa`), and rotation speed (`rad/s`, `rpm`). Config files are still
saved in solver-native units.

For `ELROD_ADAMS`, axial `DIRICHLET` and `INLET_OUTLET` inflow boundaries use
the dimensionless `film content` rows as theta boundary values. The pressure
rows remain pressure values for pressure solvers and inlet/outlet switching.
Pressure supply inlets are converted through the EOS as
`theta = exp((p_supply - p_cav) / bulk_modulus)`, so realistic bulk modulus
values are important for sane film-content magnitudes.

The inlet editor switches between `GROOVE` and `CIRCULAR`. Saving writes the
selected inlet as the active config line and keeps the unused inlet form as a
comment for quick switching later.

The `film_content` output is the Elrod/JFO universal variable used by the
solver. Cavitated cells are below `1`, and compressed full-film cells can be
above `1`; the GUI and VTK output do not clamp those values.

`load_angle_deg` is a visualization reference angle measured from the positive
y axis. The solver still stores cylindrical `theta` from the positive x axis;
the GUI preview rotates the displayed circumferential coordinate so the load
reference is easier to inspect. Positive `load_angle_deg` follows increasing
solver `theta`, so the preview's zero column corresponds to
`theta = 90 deg + load_angle_deg`.

## Output

Results are written under the configured `output_dir`, normally:

```text
C:\Dev\pancake\build-windows-mingw\results
```

Open `results.pvd` in ParaView for full visualization.

## Project Layout

```text
src/
  config.hpp            SimulationConfig and text config parsing
  mesh.hpp/cpp          Structured cylindrical grid, MPI decomposition
  field.hpp/cpp         Scalar fields with ghost layers
  communicator.hpp/cpp  MPI ghost-cell exchange
  fvm.hpp/cpp           FVM operators
  linear_system.hpp/cpp PETSc sparse system wrapper
  io.hpp/cpp            VTK/PVD writers
  reynolds.hpp/cpp      Reynolds equation assembly and post-processing
  gui_win32.cpp         Native Windows GUI shell around pancake.exe

cmake/
  deploy_mingw_runtime.cmake  Copies required MinGW/PETSc DLLs beside .exe files

tests/
  test_*.cpp            Solver and physics regression tests
```

## Linux

Linux remains supported from the same CMake source tree. Use `BUILDING.md` for
Linux and portability notes; use this README as the current Windows quickstart.
