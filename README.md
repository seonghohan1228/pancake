# pancake

`pancake` is a 2D finite-volume solver for compressible thin-film lubrication
in a journal bearing. It solves the Reynolds equation on a cylindrical
`(theta, z)` mesh with MPI parallelism and writes VTK/PVD results for ParaView.
The solver supports static eccentricity and an optional moving outer-bearing
model where the shaft remains fixed while pressure and viscous film forces move
the bearing housing. It also includes an optional first-pass film-averaged
energy equation for lubricant temperature and viscous heat generation.

The current Windows development workflow uses **MSYS2 MINGW64**. Visual Studio
is not required to build or run the native solver or GUI in the current code.

## Current Windows Workflow

This workflow builds:

```text
C:\Dev\pancake\build-windows-mingw\pancake.exe
C:\Dev\pancake\build-windows-mingw\pancake_gui.exe
```

`pancake_gui.exe` is a single self-contained executable: the C/C++ runtime is
statically linked and it depends only on Windows system DLLs, so it runs
standalone. `pancake.exe` (the solver) still needs Microsoft MPI plus its
PETSc/MinGW runtime DLLs, which the build copies beside it so File Explorer
launches work after Microsoft MPI is installed. Both executables embed the
shared logo from `src/pancake.ico`.

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

Current known status: the native build succeeds, the solver and GUI link, and
all checked-in tests pass.

## GUI Workflow

`pancake_gui.exe` opens on a single workspace built for fast parameter-tweak →
run → inspect loops. There are no tabs; everything is visible at once:

- **Action bar** (top): the primary `Run` / `Stop` buttons, an MPI `ranks`
  field, a live run-state chip (`Idle` / `Running NN% · Ns` with a progress bar
  / `Done` / `Failed`), and the file actions `Open…`, `Save`, `Save As…`,
  `Reload`, `Reset`, and `Solver…`.
- **Input rail** (left, resizable): a search box, an `Advanced` toggle, a
  preset picker, and a unit-system selector, above collapsible parameter
  sections (Geometry, Operating/Fluid, Energy, Mesh & Time, Motion/Loads,
  Cavitation, Axial Boundaries, Inlet, Output, and a collapsible Raw config.txt
  editor).
  Click a section header to collapse/expand it; type in the search box to
  filter rows across all sections.
- **Live summary** (pinned under the rail): derived journal-bearing quantities
  recomputed as you type — eccentricity ratio `ε = e/c`, `h_min = c(1−ε)`,
  `h_max`, `c/R`, `L/D`, surface speed `U = ωR`, plus the previewed field's
  mean/max after a run. A status line reports validation problems.
- **Results viewport** (center): field and timestep selectors, prev/next/refresh
  buttons, and the contour preview with axes and a unit-labeled colorbar.
  Multi-rank previews are assembled from all `processor*` VTS files for the
  selected timestep. Preview contours are not updated while the solver is
  running; the current contour stays stable until the run finishes or `Refresh`
  is clicked. `Open folder` and `Open in ParaView` are here.
- **Console** (bottom, resizable): the colorized solver process log.

**Configuration files.** The GUI defaults to `config.txt` beside the executable.
`Save` writes the current file (creating `config.txt` if none is loaded);
`Save As…` and `Open…` use the standard Windows file dialogs so you can keep a
library of named cases anywhere on disk. The loaded file name appears in the
title bar, with `*` while there are unsaved edits.

**Validation and gating.** Physically impossible or malformed inputs (`e ≥ c`,
non-positive grid/clearance, `dt`/`write_interval` ordering, empty
output names) highlight the offending field in red and are listed in the
summary card; `Run` stays disabled until the inputs are valid. Numeric edit
boxes preserve the exact text you typed.

**Units.** A global unit-system selector (`SI` ↔ engineering `mm/MPa/rpm`) sets
sensible defaults; individual fields keep a small override dropdown where more
than one unit is realistic — lengths (`m`, `mm`, `cm`, `um`), times (`s`,
`ms`), angles (`deg`, `rad`), tilt (`m/m`, `deg`, `rad`), pressure (`Pa`,
`kPa`, `MPa`), and rotation speed (`rad/s`, `rpm`). Fixed-unit fields show a
plain suffix. Config files are still saved in solver-native units.

**Presets.** The preset picker populates the form with a starting point
(Default bearing, High eccentricity, Misaligned, Circular oil hole) in one
click.

The Windows application logo is embedded from `src/pancake.ico` in both
`pancake.exe` and `pancake_gui.exe`. The editable source-style logo is
`src/pancake_logo.svg`.

For `ELROD_ADAMS`, axial `DIRICHLET` and `INLET_OUTLET` inflow boundaries use
the dimensionless `film content` rows as theta boundary values. The pressure
rows remain pressure values for pressure solvers and inlet/outlet switching.
Pressure supply inlets are converted through the EOS as
`theta = exp((p_supply - p_cav) / bulk_modulus)`, so realistic bulk modulus
values are important for sane film-content magnitudes.
The solver initializes the pressure field and flooded Elrod film-content guess
from the first configured inlet `p_supply`; it falls back to `p_cav` only when
no inlet is configured.

The inlet editor switches between `GROOVE` and `CIRCULAR`. Saving writes the
selected inlet as the active config line and keeps the unused inlet form as a
comment for quick switching later.

The `film_content` output is the Elrod/JFO universal variable used by the
solver. Cavitated cells are below `1`, and compressed full-film cells can be
above `1`; the GUI and VTK output do not clamp those values.

`motion_model = MOVING_BEARING` enables outer-bearing motion with the shaft
fixed in space. The initial moving-bearing position can be derived from
`e`/`attitude_angle_deg`, so the moving model starts from the same film
thickness as the static case. The hydrodynamic result used for motion is written
as `fluid_force_x`, `fluid_force_y`, and `fluid_force_z`. This is the total
fluid force applied to the moving bearing surface. Its components are
`pressure_force_x/y/z` for the pressure-only bearing-surface force and
`viscous_force_x/y/z` for the shear-only bearing-surface force. Legacy
`load_x/load_y/load_z` remain pressure-force aliases for older configs. The
independent applied load is written separately as `external_load_x`,
`external_load_y`, and `external_load_z`.
In the GUI, the applied load is entered as an in-plane magnitude/direction plus
an independent z component; saving converts those controls to
`external_load_x/y/z`.

The supported time-method names are `EULER_EXPLICIT`, `EULER_IMPLICIT`,
`CRANK_NICOLSON`, `RK2`, and `RK4`. `CRANK_NICOLSON` is the correct formal name
for the mixed implicit/explicit trapezoidal method; the config parser also
accepts aliases such as `SEMI_IMPLICIT`, `IMEX`, and `CRANK_NICHOLSON`.

`solution_mode = STEADY_STATE` runs one pressure/thermal solve at `t=0` and
skips bearing-motion advancement. This is intended for fast operating-point
checks. `TRANSIENT` keeps the normal time loop.

`temperature_model = ENERGY_EQUATION` solves the film-averaged thermal balance
after the Reynolds solve and force post-processing. The solver writes
`temperature` and `heat_generation` when those fields are enabled in
`output_fields`. `heat_generation` is viscous dissipation per bearing surface
area in `W/m^2`; the wall heat-transfer coefficients close the steady thermal
balance. `ISOTHERMAL` leaves the temperature field fixed but still refreshes
`heat_generation` for inspection.

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
  energy.hpp/cpp        Film-averaged energy equation and heat generation
  journal_motion.hpp/cpp Moving-bearing state and time integration
  linear_system.hpp/cpp PETSc sparse system wrapper
  io.hpp/cpp            VTK/PVD writers
  reynolds.hpp/cpp      Reynolds equation assembly and post-processing
  gui_win32.cpp         Native Windows GUI shell around pancake.exe (single-file, static)

cmake/
  deploy_mingw_runtime.cmake  Copies required MinGW/PETSc DLLs beside pancake.exe
                              (the GUI is statically linked and needs none)

tests/
  test_*.cpp            Solver and physics regression tests
```

## Linux

Linux remains supported from the same CMake source tree. Use `BUILDING.md` for
Linux and portability notes; use this README as the current Windows quickstart.
