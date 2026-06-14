# pancake

`pancake` is a 2D finite-volume solver for compressible thin-film lubrication
in a journal bearing. It solves the Reynolds equation on a cylindrical
`(theta, z)` mesh with MPI parallelism and writes VTK/PVD results for ParaView.
The solver supports static eccentricity and an optional moving outer-bearing
model where the shaft remains fixed while pressure and viscous film forces move
the bearing housing. It also includes an optional first-pass film-averaged
energy equation for lubricant temperature and viscous heat generation.

For a feature-level record of substantive project contributions, see
`docs/PROJECT_CONTRIBUTIONS.md`.

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

`pancake_gui.exe` is a Dear ImGui (docking) + ImPlot desktop app on Win32 +
Direct3D 11 (sources in `src/gui/`; see `BUILD.md` and
`docs/GUI_ARCHITECTURE.md`). It drives `pancake.exe` as a subprocess; the two
exes travel together in one folder. Dockable panels (layout is remembered in
`%LOCALAPPDATA%\PancakeGui\`):

- **Case Setup** (left): Run/Cancel bar with the MPI process count beside it
  (more processes = faster solve; needs MS-MPI above 1) and a stated reason
  while inputs are invalid, parameter filter, `Advanced` toggle, presets, and
  the grouped parameter form (labels left, tooltip help, live red-highlight
  validation, inlet editor). Groups are ordered by edit frequency — Operating
  conditions, Geometry, Lubricant, Cavitation, Time stepping, Mesh open at
  the top; Boundaries, Thermal, Motion, Fluid properties, Output, Numerics
  collapsed below. Numeric fields have display-unit selectors (`Pa/kPa/MPa`,
  `deg/rad`, `rad/s`/`rpm`, `m/mm/um`, `s/ms`) — values convert in place
  while the config file stays solver-native SI. Every dropdown option shows a
  hover description; the convection scheme offers Upwind and Linear (type
  differencing), with TVD values from existing config files preserved.
  Validation merges the solver's own `SimulationConfig::validate()`, so Run
  is enabled exactly when the solver would accept the case. `Enter`/`F5`
  runs.
- **Results** (center): field heatmap over the unwrapped bearing surface drawn
  with the true cell aspect ratio (equal-aspect; window resizes zoom rather
  than stretch, a `Fit` button in the plot corner restores the full domain),
  cavitation zone and inlet zone overlays, pan/zoom/hover readout, Viridis
  colorbar;
  mid-plane/arbitrary-slice pressure and film-thickness profiles; and a key
  results summary (min film thickness, peak pressure, attitude angle, friction
  torque/power, cavitated fraction) with one-click Copy per value.
- **Monitors**: live physical traces while the solver runs — minimum film
  thickness [µm], film temperature range [K], fluid force [N], friction
  torque [N·m], cavitated area [%], and the liquid mass residual (log) —
  selectable, each in its own stacked subplot with real engineering units
  and a shared time axis (sources: DETAILED step lines + `diagnostics.csv`).
- **Solver Log**: the colorized solver output (ANSI colors preserved),
  clipped/virtualized for long runs, with Copy-all.
- **Status bar**: solver state (idle / running step N / finished / failed /
  cancelled), elapsed time, last residual, progress vs `end_t`, and a
  persistent warning chip when `pancake.exe` is not found.

Every plot exports PNG and every table/curve exports CSV. `Results > Open in
ParaView` hands `results.pvd` to ParaView for full 3D inspection.

**Configuration files.** The GUI reads/writes the same `key = value` config
the solver consumes (shared `src/config.hpp`). It opens the last-used config
(or `config.txt` beside the exe) at startup; Run auto-saves the current case
first. Relative `output_dir` resolves next to the config file.

**Solver discovery.** `pancake.exe` next to the GUI exe is used by default;
a custom path can be set in *Tools > Settings* (persisted). MPI ranks > 1
launch through `mpiexec` (Microsoft MPI).

**Standalone binary.** The GUI is a single statically linked exe: fonts
(Roboto, JetBrains Mono) are embedded byte arrays, the import table contains
only Windows-shipped DLLs, all runtime writes go to `%LOCALAPPDATA%\PancakeGui\`,
and it runs from read-only locations. Per-monitor-v2 DPI aware.

The Windows application logo is embedded from `src/pancake.ico` in both
`pancake.exe` and `pancake_gui.exe`. The editable source-style logo is
`src/pancake_logo.svg`.

For `ELROD_ADAMS`, axial `DIRICHLET` and `INLET_OUTLET` inflow boundaries are
configured by pressure. The solver derives the boundary film content internally
from `bc_z_*_val`, `p_cav`, and `bulk_modulus`:
`theta = exp((max(p_bc, p_cav) - p_cav) / bulk_modulus)`. Legacy
`bc_z_*_theta` keys are still parsed for old files but are ignored by the solve
and are no longer emitted by the GUI. Pressure supply inlets use the same EOS
conversion, so realistic bulk modulus values are important for sane
film-content magnitudes. The solver initializes the pressure field uniformly
from the first non-Neumann outlet pressure, not from supply-groove pressure.
Supply inlets are still imposed during the solve.

The inlet editor switches between `GROOVE` and `CIRCULAR`. Saving writes the
selected inlet as the active config line and keeps the unused inlet form as a
comment for quick switching later.

The `film_content` output is the Elrod/JFO universal variable used by the
solver. Cavitated cells are below `1`, and compressed full-film cells can be
above `1`; the GUI and VTK output do not clamp those values.
In MPI runs, the Elrod active-set loop synchronizes solved `theta` ghost cells
between outer iterations so rank boundaries use current cavitation flags, not
stale timestep data.

`motion_model = MOVING_BEARING` enables outer-bearing motion with the shaft
fixed in space. The initial moving-bearing position can be derived from
`e`/`attitude_angle_deg`, so the moving model starts from the same film
thickness as the static case. The hydrodynamic result used for motion is written
as `fluid_force_x`, `fluid_force_y`, and `fluid_force_z`. This is the total
fluid force applied to the moving bearing surface. Its components are
`pressure_force_x/y/z` for the pressure-only bearing-surface force and
`viscous_force_x/y/z` for the shear-only bearing-surface force. Legacy
`load_x/load_y/load_z` output requests are mapped to `pressure_force_x/y/z`
for older configs, but the solver no longer allocates duplicate `load_*`
fields. The independent applied load is written separately as `external_load_x`,
`external_load_y`, and `external_load_z`.
In the GUI, the applied load is entered as an in-plane magnitude/direction plus
an independent z component; saving converts those controls to
`external_load_x/y/z`.

Field solves (pressure/theta, temperature, gas transport) are backward-Euler
in time; only the bearing equation of motion has a selectable integrator
(`motion_time_method`: `EULER_EXPLICIT`, `EULER_IMPLICIT`, `CRANK_NICOLSON`,
`RK2`, `RK4`). The old `pressure_time_method` / `temperature_time_method` keys
never changed the field solves and have been removed; old configs still parse
with a one-line "ignored" warning.

Convection face interpolation is selectable per transported variable with
`theta_convection_scheme`, `thermal_convection_scheme`, and
`gas_convection_scheme` (`UPWIND` default; `TVD_VANLEER`, `TVD_MINMOD`; theta
additionally accepts `TYPE_DIFFERENCING` — central across full-film faces,
upwind in cavitation, per Vijayaraghavan & Keith 1989). The Gumbel wedge term
uses the same scheme, so the two Reynolds paths discretize like-for-like.
`linear_rtol` controls the inner Krylov tolerance. `diagnostics_interval` sets
the cadence of the per-step global mass-balance / convergence CSV written to
`<output_dir>/diagnostics.csv` (masses, boundary fluxes, inlet sources,
clamp-mass accounting, normalized residuals, Elrod outer iterations and flag
flips; see PHYSICS.md §15). Startup validation also reports laminar-regime
(Taylor number) and gas-saturation (`p_sat`) context warnings. A
grid-convergence study harness lives in `validation/grid_convergence/`.

`solution_mode = STEADY_STATE` runs one pressure/thermal solve at `t=0` and
skips bearing-motion advancement. This is intended for fast operating-point
checks. `TRANSIENT` now starts from the initialized outlet-pressure fields and
uses transient terms from the first solve. Set `omega_ramp_time > 0` to ramp
rotation linearly from zero to the configured `omega` over that physical time;
`omega_ramp_time = 0` applies full speed immediately.
For every transient `ELROD_ADAMS` run — fixed or moving bearing — the
`rho*h*dtheta/dt` film-content capacity term is assembled in the theta solve;
the lagged `theta*dh_dt` squeeze source is added when geometry moves. Dropping
the capacity term is non-conservative, makes fixed-bearing transients
quasi-steady in `theta`, and can create artificial cavitation after startup.

`temperature_model = ENERGY_EQUATION` solves the film-averaged thermal balance
after the Reynolds solve and force post-processing. The solver writes
`temperature` and `heat_generation` when those fields are enabled in
`output_fields`. `heat_generation` is viscous dissipation per bearing surface
area in `W/m^2`. In `ELROD_ADAMS`, the pressure-flow part of this heat source is
computed only from active full-film faces; cavitated cells and inactive
cavitation-front faces do not add Poiseuille heat. Large values at the
full-film side of a cavitation boundary are pressure-gradient dissipation, not
gas phase-change heat unless the opt-in gas model and mass-transfer rate are
enabled. Supply inlets fix pressure; under Elrod-Adams that pressure is
converted to film content internally. They do not impose temperature. Axial pressure-boundary inflow uses
`temperature_reference` as the external reservoir temperature, while axial
outflow uses an upwind zero-gradient temperature treatment. The wall
heat-transfer coefficients close the steady thermal balance. `ISOTHERMAL`
leaves the temperature field fixed but still refreshes `heat_generation` for
inspection.

`fluid_property_model = CONSTANT` keeps the previous pure-oil behavior. The
new opt-in `OIL_DISSOLVED_GAS` and `GAS_CAVITATION_MIXTURE` modes add a
liquid oil plus dissolved gas property layer before the JFO cavitation solve.
For propane/R290, the configurable `oil_gas_solution_model`, `density_model`,
and `viscosity_model` keys can use Henry/Bunsen-style correlations or
`x:value` tables. All pressure inputs are absolute Pa; use `p_cav = 1e5` for an
atmospheric cavitation plateau, not `0` gauge. `config_r290_pz68.txt` is the
calibrated PZ68/R290 starting case for 40 C and 10 bar, and
`config_r290_pz68_quick.txt` is the same physics on a small smoke-test mesh.
Both use 10 bar supply/reference saturation with a 1 bar absolute sump/outlet
so the quick case actually releases gas in cavitated cells. They set Henry
solubility, dissolved-propane liquid-density mixing,
Andrade/Barus/log-mixing viscosity coefficients, finite-rate gas transfer, and
diffusivity explicitly. The solver validates required gas constants before a
run and rejects under-specified custom gas models.
The effective film viscosity with free gas present is selected by
`gas_mixture_viscosity_model`: `EINSTEIN_DILUTE` (back-compat default, valid
only for small void fractions and rejected when `gas_alpha_max > 0.6`),
`DUKLER_VOID`, `MCADAMS_QUALITY` (shipped R290 default), `KRIEGER_DOUGHERTY`,
or `LINEAR_QUALITY`; the void/quality-weighted models thin toward the
configured `mu_gas`. The Henry tangent and `exp(a_c c_d)` viscosity correction
are calibrated at (10 bar, 40 C) only. For the full operating map, measured PTSV
data loads from a runtime-editable CSV via `solubility_table_file` /
`viscosity_table_file` (a `P_MPa,T_C,solubility_pct,viscosity_mm2s` grid;
`data/R290_PZ68S/r290_pz68s.csv` ships one) with `oil_gas_solution_model = TABLE`
and `viscosity_model = TABLE`. The grid resolution is whatever the file holds (no
recompile); the kinematic table viscosity is converted to dynamic per cell with
the solution density, and the saturation surface is inverted to the local
bubble point. The solver writes
`rho_liquid_solution`, `mu_liquid_solution`, `mu`, `dissolved_gas`,
`free_gas_mass`, `alpha_gas`, and `gas_mass_transfer` diagnostics when
requested. Released free gas is transported as `free_gas_mass` with the
film-averaged velocity and resorbs only through the finite-rate transfer law;
it is no longer deleted just because a pressure inlet locally restores
full-film JFO `theta`. `theta` remains the liquid-content variable used by
Elrod-Adams. Pressure supply/inlet cells impose the configured incoming mixture
composition (`dissolved_gas_initial`, zero `free_gas_mass`, zero `alpha_gas`,
and zero local `gas_mass_transfer`). Axial pressure-boundary inflow uses the
same incoming mixture, while axial outflow convects dissolved gas and released
free gas out of the domain. The GUI exposes these selections and gas-property
constants in the Fluid Properties section with the Advanced toggle.

Several gated cavitation extensions complete the oil-refrigerant model (each
defaults to the prior behaviour; see PHYSICS.md §16 and PLAN.md Phase E).
`cavitation_threshold = LOCAL_PSAT` sets the onset to the local bubble point
`max(p_cav, p_sat(c_d,T))`, so the higher-release-pressure species (propane)
governs rupture rather than the oil `p_cav`. `gas_pressure_coupling =
VOID_COUPLED` lets released gas open void and relieve pressure (a void ceiling,
local plateau, and mixture compressibility) instead of being volumetrically
inert. `cavitated_film_model = STRIATED` weights Couette shear and heat
generation by the liquid fraction in cavitated cells. `consistent_boundary_flux`
re-floods cavitated axial-boundary cells at the physical Poiseuille rate. In
STEADY_STATE the coupled sweep iterates to a residual (`coupling_max_iters`,
`coupling_tolerance`) and `steady_gas_model` selects the steady gas treatment.
`config_r290_pz68_table.txt` exercises the table + LOCAL_PSAT + VOID_COUPLED +
STRIATED stack end-to-end.

Fast gaseous-cavitation smoke run:

```powershell
cd C:\Dev\pancake\build-windows-mingw
.\pancake.exe -c config_r290_pz68_quick.txt
```

When `velocity` is enabled in `output_fields`, VTK output writes `U` as the
active 3-component cell vector so ParaView reads it directly. Curved output
stores Cartesian `(U_X, U_Y, U_Z)` on the cylinder; flat output stores
`(u_theta, u_z, 0)` on the unwrapped plane. Resultant force/load/bearing
component fields are grouped into `Fp`, `Fv`, `F`, `Fext`, and `xB` vector
arrays rather than duplicate scalar component arrays.

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
If a run exits with `Output directory setup failed` on Windows, close ParaView
or any preview using that result folder, or change `output_dir` for the next
run. The solver clears the selected output directory before writing new VTK/PVD
files so stale timesteps cannot be mixed into a new case.

## Project Layout

```text
src/
  config.hpp            SimulationConfig and text config parsing
  mesh.hpp/cpp          Structured cylindrical grid, MPI decomposition
  field.hpp/cpp         Scalar fields with ghost layers
  communicator.hpp/cpp  MPI ghost-cell exchange
  fluid_properties.hpp/cpp
                         Oil/dissolved-gas properties and free-gas source terms
  fvm.hpp/cpp           FVM operators
  energy.hpp/cpp        Film-averaged energy equation and heat generation
  journal_motion.hpp/cpp Moving-bearing state and time integration
  linear_system.hpp/cpp PETSc sparse system wrapper
  io.hpp/cpp            VTK/PVD writers
  reynolds.hpp/cpp      Reynolds equation assembly and post-processing
  gui/                  ImGui + ImPlot desktop GUI around pancake.exe
                        (see docs/GUI_ARCHITECTURE.md; built per BUILD.md)
  gui_win32.cpp         Legacy Win32 GUI (no longer built; pending removal)

third_party/
  fonts/, stb/          Assets embedded into pancake_gui.exe at build time

cmake/
  deploy_mingw_runtime.cmake  Copies required MinGW/PETSc DLLs beside pancake.exe
                              (the GUI is statically linked and needs none)
  embed_binary.cmake          Converts vendored fonts into C array headers

tests/
  test_*.cpp            Solver and physics regression tests
```

## Linux

Linux remains supported from the same CMake source tree. Use `BUILDING.md` for
Linux and portability notes; use this README as the current Windows quickstart.
