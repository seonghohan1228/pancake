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
  collapsed below. Numeric fields have display-unit selectors (`Pa/kPa/MPa/bar`,
  `GPa` for bulk modulus, `K/degC`, `deg/rad/mrad`, `rad/s`/`rpm`, `m/mm/um`,
  `s/ms`, mass fraction/percent, and engineering stiffness/damping units) — values convert in place
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
  or blue-white-red color scale, and optional manual color limits;
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
first. Relative `output_dir` resolves next to the config file. Cavitation
threshold and released-gas pressure coupling are editable in the Cavitation
panel. TABLE property models expose CSV file pickers and a live loaded-grid
status; the selected paths are preserved when the GUI saves the case.

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

For `JFO`, axial `DIRICHLET` and `INLET_OUTLET` inflow boundaries are
configured by pressure. The solver derives the boundary film content internally
from `bc_z_*_val`, `p_cav`, and `bulk_modulus`:
`theta = exp((max(p_bc, p_cav) - p_cav) / bulk_modulus)`. Legacy
`bc_z_*_theta` keys are still parsed for old files but are ignored by the solve
and are no longer emitted by the GUI. Pressure supply inlets use the same EOS
conversion, so realistic bulk modulus values are important for sane
film-content magnitudes. The solver initializes the pressure field uniformly
from the first non-Neumann outlet pressure, not from supply-groove pressure.
Supply inlets are still imposed during the solve.

`cavitation_model = FBNS` solves the same JFO conditions as a complementarity
problem instead of through the Elrod-Adams bulk-modulus formulation. It imposes
`p - p_cav >= 0`, `theta_full - theta >= 0` and
`(p - p_cav)(theta_full - theta) = 0` directly via
the Fischer-Burmeister function and a Newton iteration with a Schur reduction
(Woloszynski, Podsiadlo & Stachowiak, Tribol. Lett. 58:18, 2015). There is no
artificial bulk modulus and no switch function, so `bulk_modulus` is unused on
this path. On the exact 1D benchmark of `cases/00c_jfo_verification` it places
the rupture front within one cell at every resolution from 120 to 960
circumferential cells, where `JFO` is 7-13 cells off, and it converges in four
Newton iterations independent of grid size. Scope: `temperature_model`
`ISOTHERMAL` or `ENERGY_EQUATION`, `fluid_property_model` `CONSTANT`,
`SINGLE_PHASE` (so THD runs with temperature- and pressure-dependent oil
properties are supported) or `TWO_PHASE`, `cavitated_film_model` `FULL_FILM` or
`STRIATED`, `DIRICHLET`/`NEUMANN` axial boundaries, `cross_film_thermal_model`
`POLYNOMIAL`; steady or transient, static or moving bearing, MPI-parallel.
Temperature reaches the solver through its per-cell Poiseuille mobility, which is
frozen inside each Newton solve and refreshed by the outer Picard loop, exactly
as on the Elrod path.

With `cross_film_thermal_model = POLYNOMIAL` that mobility is the generalized
(Dowson) one, `rho_l F2`, and the Couette flux is carried by `2 F1/F0` instead of
`h` — so `mu(T(y))` acts inside the film integrals rather than through a
depth-averaged viscosity. `CrossFilm::update` rebuilds the integrals once per
outer iteration, before the Newton solve, so they are frozen coefficients like
`mu` and both Jacobian blocks, the linearity of the mass balance and the
quadratic Newton rate are unchanged. See PHYSICS.md sections 7.3.9 and 13.1.3.

Under `TWO_PHASE` the full-film ceiling `theta_full` drops from 1 to
`1 - alpha_gas`, because the released gas occupies part of the gap and the liquid
can only fill what is left. `alpha_gas` is lagged (the gas path runs after the
flow solve), so it is a coefficient of the Newton solve rather than an unknown and
both Jacobian blocks keep their gas-free forms. The ceiling has one definition,
`FullFilm::bound` in `src/full_film.hpp`, which every full-film test in the solver
reads; it returns exactly 1 for every model other than FBNS + `TWO_PHASE`, so the
Elrod and Gumbel paths are unchanged bit for bit. The `alpha_gas <-> theta`
feedback this opens takes the same under-relaxation lever the Elrod void-coupled
path uses, `void_coupling_relaxation` (default 0.3); on a strongly supersaturated
feed it is worth lowering. See `docs/FBNS_GAS_COUPLING.md` for the limits.

Still rejected at config validation rather than silently ignored:
`gas_pressure_coupling = VOID_COUPLED` (it makes the cavitation plateau a per-cell
field that only the Elrod outer loop advances), the cross-film thermal model, and
`INLET_OUTLET` axial boundaries - use `cavitation_model = JFO` for
those. Tuning keys: `fbns_max_iters`,
`fbns_tol`, `fbns_relaxation`, and `fbns_pressure_scale` (0 = derived from the
operating point; it conditions the Fischer-Burmeister function and rarely needs
setting). See PHYSICS.md §7.3.

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
`gas_convection_scheme` (`UPWIND` default; `LINEAR` — 2nd-order central, unbounded;
`VANLEER` — 2nd-order TVD, bounded; theta additionally accepts `TYPE_DIFFERENCING`
— central across full-film faces, upwind in cavitation, per Vijayaraghavan & Keith 1989). The Gumbel wedge term
uses the same scheme, so the two Reynolds paths discretize like-for-like.
`linear_rtol` controls the inner Krylov tolerance. `linear_ksp_type` (`BCGS`
default, plus `BCGSL`, `GMRES`, `FGMRES`), `linear_pc_type` (`ASM` default, plus
`BJACOBI`, `JACOBI`, `NONE`) and `linear_pc_overlap` (default 1) select the
Krylov method and preconditioner for every inner solve. The defaults are chosen
for rank-independence, not for the lowest cost per solve: the mesh is decomposed
in theta only, so every rank boundary cuts the flow direction of the
convection-dominated operators, and block Jacobi — which discards the coupling
across exactly those cuts — gets weaker as the rank count rises. Before
2026-08-05 that made a steady THD or Elrod solve diverge at one rank count while
the identical configuration converged at another. Additive Schwarz with one cell
of overlap keeps that coupling; on the Ferron THD force-balance iterate the total
Krylov iteration count is now flat in rank count (3723/3755/3761/3735 at 1/2/4/6
ranks, against 3723/4010/diverged/3961 before). `linear_solver_fallbacks`
(default 2, 0 disables) is the number of escalation rungs a failed solve is
retried on — a different Krylov path, not merely a tighter one — before the run
aborts; retries cost nothing on a solve that converges and each one that fires is
logged. `diagnostics_interval` sets
the cadence of the per-step global mass-balance / convergence CSV written to
`<output_dir>/diagnostics.csv` (masses, boundary fluxes, inlet sources,
clamp-mass accounting, normalized residuals, Elrod outer iterations and flag
flips; see PHYSICS.md §15). Startup validation also reports laminar-regime
(Taylor number) and gas-saturation (`p_sat`) context warnings. A
grid-convergence study harness lives in `validation/grid_convergence/`.

Adding `mass_residual` to `output_fields` writes the per-cell liquid-continuity
residual of the converged mass-conserving state in `kg/s` (cell-integrated), which localizes
in 2D what the global balance can only report as one number: a defect in the
boundary or groove flux handling shows up as a hot cell, while the global
residual stays at solver tolerance because neighbouring errors cancel. Away from
groove cells, clamped cells and the rupture front the field must sit at round-off
of the barotropic conductance. The field is registered only when selected, so
runs that omit it are bit-identical. It is defined for `cavitation_model = JFO`
and `FBNS` (each mirroring its own solver's face fluxes); validation warns and
writes zeros otherwise. Both the solvers and this diagnostic assemble their
coefficients from one header, `src/reynolds_coefficients.hpp`, so the residual
reported is the residual of the operator actually assembled — including under the
cross-film model and `gas_pressure_coupling = VOID_COUPLED`, where local copies
had previously drifted. See PHYSICS.md §15.1.

`solution_mode = STEADY_STATE` runs one pressure/thermal solve at `t=0` and
skips bearing-motion advancement. This is intended for fast operating-point
checks. `TRANSIENT` now starts from the initialized outlet-pressure fields and
uses transient terms from the first solve. Set `omega_ramp_time > 0` to ramp
rotation linearly from zero to the configured `omega` over that physical time;
`omega_ramp_time = 0` applies full speed immediately.
For every transient `JFO` run — fixed or moving bearing — the
`rho*h*dtheta/dt` film-content capacity term is assembled in the theta solve;
the lagged `theta*dh_dt` squeeze source is added when geometry moves. Dropping
the capacity term is non-conservative, makes fixed-bearing transients
quasi-steady in `theta`, and can create artificial cavitation after startup.

`temperature_model = ENERGY_EQUATION` solves the film-averaged thermal balance
after the Reynolds solve and force post-processing. The solver writes
`temperature` and `heat_generation` when those fields are enabled in
`output_fields`. `heat_generation` is viscous dissipation per bearing surface
area in `W/m^2`. In `JFO`, the pressure-flow part of this heat source is
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

`cross_film_thermal_model = POLYNOMIAL` (default `NONE`) replaces the
uniform-across-the-gap assumption of the energy equation with a reconstructed
`T(y)` profile. It is a single switch for three coupled things: the cross-film
profile, a conjugate wall model (circumferentially isothermal shaft whose
temperature is solved from the global film-to-shaft heat balance, plus 1-D radial
conduction through the bush to an outer convection boundary), and a generalized
Reynolds equation using Dowson's `F0/F1/F2` viscosity integrals in place of
`h^3/(12 mu)` (under both `cavitation_model = JFO` and `FBNS`). When it is active, `journal_heat_transfer`,
`bearing_heat_transfer`, `journal_wall_temperature` and
`bearing_wall_temperature` are ignored — the conjugate model replaces them, and
the config validator warns rather than letting them look effective. Its own keys
are `cross_film_profile_order` (2 or 4), `cross_film_quadrature_points`,
`cross_film_relaxation`, `shaft_ambient_heat_transfer`,
`shaft_ambient_area_ratio`, `thermal_ambient_temperature`,
`shaft_ambient_temperature`, `bush_thermal_conductivity`, `bush_outer_radius`,
`bush_outer_heat_transfer` and `bush_outer_temperature`. The two optional
surface temperatures use `thermal_ambient_temperature` when set to `0` or a
negative value. With `NONE` nothing in the module runs and results
are bit-identical to the model-free solver. The equations as implemented, and
the two stated approximations inside the order-4 collocation, are in PHYSICS.md
sections 13.1.2 and 13.1.3.

`bush_thermal_model` selects how the bush closes that flux, and is read only when
the cross-film model is active:

| value | bush model |
| --- | --- |
| `LOCAL` (default) | per-cell 1-D radial conduction in series with the outer surface. Every watt entering the bore leaves through the wall underneath it, so the bush surface can only peak where the film is hottest. |
| `CONDUCTION_2D` | reduced steady 2-D `(r, theta)` conduction in the bush annulus `[R, bush_outer_radius]`. Heat picked up in the loaded arc can leave at a different angular position. |

`CONDUCTION_2D` adds `bush_radial_cells` (default 16). Both models share
`bush_outer_boundary`, which is `CONVECTION` (Robin with
`bush_outer_heat_transfer` to `bush_outer_temperature`, falling back to
`thermal_ambient_temperature`; `0` is an adiabatic outer surface) or
`FIXED_TEMPERATURE` (that surface held at the same selected temperature — the
large-HTC limit, expressed exactly). The 2-D bush is z-collapsed: it consumes the
z-averaged film flux and returns one bore temperature per theta. Ferron's solid
model retains the axial term and applies \(h_b\) on the axial faces; those faces
are absent here. The scalar shaft closure likewise is not Ferron's r-z shaft
conduction model with 43 C axial-end temperatures. Its linear system is solved
directly once per outer iteration and costs about 10^7 flops at 720 x 16. `LOCAL`
with `CONVECTION` is bit-identical to
the solver before the option existed. See PHYSICS.md §13.1.2.1.

While the model is active, four optional output fields expose the film surfaces
the depth-averaged `temperature` cannot show. Select them in `output_fields`:

| field | meaning | unit |
| --- | --- | --- |
| `temperature_journal` | `T(eta=0)`, shaft-side film surface temperature | K |
| `temperature_bearing` | `T(eta=1)`, bush-side film surface temperature | K |
| `heat_flux_journal` | wall heat flux at the shaft-side surface | W/m^2 |
| `heat_flux_bearing` | wall heat flux at the bush-side surface | W/m^2 |

Both fluxes are **positive into the fluid**, the opposite sign to the internal
`q_shaft`/`q_bush` of PHYSICS.md §13.1.2, which are positive for heat leaving the
film: a cooled bush therefore gives a negative `heat_flux_bearing`, and an
adiabatic wall gives an identically zero field. `temperature` is unchanged by
this — it stays the depth-averaged `T_bar`, in absolute kelvin like every
temperature in the solver. The fields are evaluations of the same profile the
`F`-integrals integrate, are allocated only when selected (so an unselected run
is byte-identical), and selecting one with `cross_film_thermal_model = NONE` is a
config error rather than a field of zeros. VTK writes them as `T_journal`,
`T_bearing`, `q_journal`, `q_bearing`.

`groove_mixing_model` (default `PINNED`) controls hot-oil carry-over at fed
inlets. `PINNED` holds groove cells at `t_supply`, the historical behaviour.
`FIXED_LAMBDA` and `FLOW_WEIGHTED` instead impose
`lambda*t_supply + (1-lambda)*T_upstream`, taking `lambda` either from the
inlet's own token or from its supply flow rate against the solved upstream
Couette carry-over. Both read two optional trailing tokens on the inlet line:
`inlet_groove = theta width p_supply t_supply [lambda] [q_supply]`.

`fluid_property_model = CONSTANT` keeps the previous pure-oil behavior. The
new opt-in `SINGLE_PHASE` (oil + dissolved gas, single-phase) and `TWO_PHASE`
(adds a free-gas phase) modes add a liquid oil plus dissolved gas property layer
before the JFO cavitation solve. For propane/R290, the configurable
`solubility_model`, `density_model`, and `liquid_viscosity_model` keys use a
Henry's-law correlation or measured tables. All pressure inputs are absolute Pa; use `p_cav = 1e5` for an
atmospheric cavitation plateau, not `0` gauge. `config_r290_pz68.txt` is the
calibrated PZ68/R290 starting case for 40 C and 10 bar, and
`config_r290_pz68_quick.txt` is the same physics on a small smoke-test mesh.
Both use 10 bar supply/reference saturation with a 1 bar absolute sump/outlet
so the quick case actually releases gas in cavitated cells. They set Henry
solubility, dissolved-propane liquid-density mixing,
Andrade/Barus viscosity coefficients, finite-rate gas transfer, and
diffusivity explicitly. The solver validates required gas constants before a
run and rejects under-specified custom gas models.
For a THD run on a single-grade oil, `liquid_viscosity_model = POLYNOMIAL`
selects the Ferron, Frene & Boncompain (1983) quadratic
`mu(T) = mu * (k0 - k1*Tbar + k2*Tbar^2)` with coefficients `viscosity_poly_k0`,
`viscosity_poly_k1`, `viscosity_poly_k2` (defaults `1, 0, 0`, the identity).
**`Tbar` is the ratio of temperatures in degrees CELSIUS,
`T[C] / property_reference_temperature[C]`, not Kelvin** — that is the paper's
own convention and the only reading under which its coefficients reproduce the
oil. The config key itself stays in Kelvin like every other temperature; the
conversion happens inside the viscosity law. Ferron's ISO VG32 set is
`viscosity_poly_k0 = 3.287`, `viscosity_poly_k1 = 3.064`,
`viscosity_poly_k2 = 0.777` with `mu = 0.0277` and
`property_reference_temperature = 313.15`, which gives exactly `mu` at 40 C and
`0.0074 Pa.s` at 80 C. The law carries no pressure or dissolved-gas term (use
`EMPIRICAL` for those); `mu` is floored at `1e-3 * mu` with one warning per run
if the polynomial is extrapolated far enough to go negative. See PHYSICS.md
§13.2.1.

The effective film viscosity with free gas present is selected by
`mixture_viscosity_model`: **`MCADAMS`** (default; homogeneous reciprocal blend by
mass quality, correct asymptotes), `LINEAR` (mass-quality linear blend, Grando
2006, refrigerant-oil validated), or `DUKLER` (legacy void-weighted linear blend,
kept for comparison). All thin toward the configured `mu_gas`. The Henry tangent and `exp(a_c c_d)` viscosity correction
are calibrated at (10 bar, 40 C) only. For the full operating map, measured PTSV
data loads from a runtime-editable CSV via `solubility_table_file` /
`viscosity_table_file` (a `P_MPa,T_C,solubility_pct,viscosity_mm2s` grid;
`data/R290_PZ68S/r290_pz68s.csv` ships one) with `solubility_model = TABLE`
and `liquid_viscosity_model = TABLE`. The grid resolution is whatever the file holds (no
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
  reynolds_coefficients.hpp
                         Single definition of the assembled operator coefficients
                         (Poiseuille/Dowson mobility, Couette thickness and face
                         flux, void/EOS parameters), shared by reynolds.cpp,
                         fbns.cpp and diagnostics.cpp
  fbns.hpp/cpp          FBNS complementarity cavitation solver (cavitation_model = FBNS)
  gui/                  ImGui + ImPlot desktop GUI around pancake.exe
                        (see docs/GUI_ARCHITECTURE.md; built per BUILD.md)

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
