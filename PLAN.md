# Development Plan

## Completed Phases

### Phase 1 — Infrastructure
Mesh, Field, Communicator, IO, FVM operators, PETSc linear system. All working with MPI domain decomposition in $\theta$.

### Phase 2 — Physics Foundations
Film thickness (static eccentric), barotropic EOS ($p = p_{cav} + \beta\ln\theta$), VTK output.

### Phase 3 — Reynolds Solver (Gumbel)
Steady-state compressible Reynolds equation with lagged density. Gumbel cavitation (post-solve pressure clamp). BiCGStab + block Jacobi. All tests pass.

### Windows GUI Infrastructure
Native Win32 GUI target with tabbed parameter controls, output-field selection,
solver launch/stop, process log capture, and lightweight VTK heatmap preview.
This phase is intentionally an external orchestrator around `config.txt`,
`pancake.exe`, and generated VTK files so the numerical core remains shared.
The active Windows workflow is now MSYS2 MinGW64: one CMake preset builds
`pancake.exe`, `pancake_gui.exe`, tests, and copied runtime DLLs. The preset
also injects the MinGW bin directory into PATH for PowerShell builds. The main
`Workspace` tab uses a left/right layout: result preview and process log on the
left, and a scrollable one-row-per-parameter inspector on the right. Solver
backend controls are always visible in the top-right header across all tabs.
The preview now supports timestep selection, axes, a unit-labeled colorbar,
literal numeric edit preservation, raw Elrod/JFO `film_content` output, and a
`load_angle_deg` visualization reference measured from the positive y
axis. The GUI assembles multi-rank preview fields, groups axial boundary
controls by south/north side, uses structured inlet controls for groove/circular
inlet definitions, and provides per-value unit selectors while saving
solver-native config units. Both Windows executables embed the shared Pancake
logo icon so the GUI and backend present the same application identity. Elrod
axial boundaries are now pressure inputs in the GUI and config writer; the
solver derives film content internally from pressure, `p_cav`, and
`bulk_modulus`.

Follow-up correction (2026-06-08): deprecated manual `bc_z_*_theta` boundary
inputs were removed from generated configs and the GUI. The parser still
accepts old keys, but Reynolds, velocity/force post-processing, energy boundary
fluxes, and initialization use pressure-derived Elrod boundary values. Tests
now guard both the pressure-derived theta and the ignored legacy theta path.

Follow-up correction (2026-06-08): output-directory setup now fails with a
clear solver message when Windows cannot clear a locked `results/...` folder,
instead of terminating on an uncaught filesystem exception. Users should close
ParaView/preview handles or select a different `output_dir` before rerunning.

### GUI Modernization (Dear ImGui docking + ImPlot, 2026-06-11)
Rewrote the GUI per `GUI_MODERNIZATION_DIRECTIVE.md` as `src/gui/` on Win32 +
Direct3D 11: worker-thread solver subprocess with a live log-scale convergence
plot, ImPlot field heatmap with cavitation-zone overlay, run history with
overlay comparison, CSV/PNG export everywhere, validation that merges the
solver's own `SimulationConfig::validate()`, and a single statically linked
exe (embedded fonts, `%LOCALAPPDATA%\PancakeGui` settings, per-monitor-v2
DPI). See `AGENT_REPORT.md`, `BUILD.md`, `docs/GUI_ARCHITECTURE.md`.
Remaining follow-ups: human smoke test of the interactive workflow; delete
legacy `src/gui_win32.cpp` once this branch's pending edits to it are
settled; solver-side proposals (structured `PROGRESS` stdout line,
cooperative cancellation) remain open in `AGENT_REPORT.md` section 6.

### Windows GUI Redesign (2026-05-29)
Rebuilt the GUI UI layer around a single tab-less workspace (action/run bar,
collapsible input rail with a pinned live-summary card, result viewport, console)
to support fast tweak→run→inspect loops. Replaced group-box frames with painted
section headers/labels to remove the long-standing z-order/paint glitches. Added
live derived journal-bearing quantities and inline validation that gates `Run`,
a run-state chip with a stdout-driven progress bar, hybrid units (global SI ↔
engineering selector plus targeted per-field overrides), scenario presets, and
custom-config Open/Save via Windows file dialogs. `pancake_gui.exe` now links
statically into a single self-contained binary (no runtime DLLs deployed beside
it). The solver/preview/log/unit backends were preserved; `config.txt` format
and the solver contract are unchanged.

### Windows GUI Motion/Preview Update (2026-06-04)
Added structured Motion / Loads controls for moving-bearing selection,
pressure/motion/temperature time-method selection, external-load
magnitude/direction/z input, bearing initial state, support
stiffness/damping, and film-thickness safety settings. The GUI now keeps the
current contour stable while the solver is running and refreshes results only
on completion or explicit user refresh. Preview timestep switching now uses a
cached per-step list of `processor*` VTS files instead of recursively scanning
the whole results tree for each field/step change.

### Live Summary Numerical Diagnostics (2026-06-07)
Added a nominal circumferential Couette Courant number to the live summary:
$Co_\theta = |\omega|\Delta t/(2\Delta\theta)$ with
$\Delta\theta = 2\pi/N_\theta$. This is an input-summary diagnostic for the
surface-driven Reynolds advection term; pressure-driven Poiseuille velocities
remain post-solve quantities. Added `docs/PROJECT_CONTRIBUTIONS.md` to summarize
substantive solver, model, output, testing, and documentation contributions
without listing one-off repairs or GUI construction work.

### Elrod MPI Active-Set Synchronization (2026-06-07)
The Elrod-Adams outer active-set loop now exchanges solved `theta` ghost cells
after every outer solve, clamp, and near-full-film snap. This keeps
$g(\theta)$ current at MPI partition interfaces when running with multiple
ranks. Added a post-solve ghost-consistency regression and a 6-rank
`test_elrod_a2_mpi6` CTest entry because the visible contour artifact appeared
with six processor pieces.

### Solver Initialization Update (2026-06-04)
The timestep flooded guess now initializes `pressure` from the first configured
inlet `p_supply`, not `p_cav`. The Elrod `theta` field is initialized from the
same pressure through the EOS, with `p_cav` used only as the no-inlet fallback.

### Bearing Force Convention Update (2026-06-04)
Corrected macroscopic force integration to use one bearing-side convention:
`pressure_force_*` is the pressure force applied to the moving bearing surface,
including the bearing surface area vector
`((R+h)e_r - h_theta e_theta - (R+h)h_z e_z) dtheta dz`; `viscous_force_*` is
the bearing-side shear force; `fluid_force_*` is their sum. Legacy `load_*`
output requests map to `pressure_force_*`, but duplicate runtime `load_*`
fields are no longer allocated.

### Energy Equation Initial Implementation (2026-06-04)
Added `solution_mode = STEADY_STATE` for one-step operating-point calculations
and `temperature_model = ENERGY_EQUATION` for a first-pass film-averaged thermal
solve. The solver now carries `temperature` and `heat_generation` fields,
computes Couette plus pressure-flow viscous heat generation, applies wall heat
transfer to the journal and bearing, and exposes the thermal controls in the
Win32 GUI. The energy solve uses non-conservative upwind thermal advection to
avoid spurious `T div(F)` source terms near pressure-supply regions, treats
pressure supply inlets as pressure-only constraints, and uses
`temperature_reference` only for axial pressure-boundary thermal inflow.
Transient runs now assemble the `t=0` solve as steady and store that solved
state before applying transient pressure/temperature terms at `t=dt`.
Thermoviscous viscosity feedback, deformation, and coupled ETHD iterations
remain in Phase D.

Follow-up correction (2026-06-07): `heat_generation` now evaluates the
pressure-flow dissipation from active full-film faces instead of one-sided
cell-centered gradients at a cavitation front. This keeps the hard
Elrod-Adams/JFO pressure plateau intact while preventing inactive cavitated
faces from creating mesh-dependent corner heat spikes in the derived thermal
source.

### Oil/Propane Mixture Property Layer (2026-06-05)

Added the first JFO-compatible dissolved-gas implementation slice. The solver
now owns `src/fluid_properties.hpp/cpp`, with config-driven enums for
`fluid_property_model`, `dissolved_gas_species`, `oil_gas_solution_model`,
`density_model`, and `viscosity_model`. Defaults remain
`CONSTANT`/`PURE_OIL`, so existing pure-oil cases keep their old numerical
behavior.

Implemented fields:
- `rho_liquid_solution`, `mu_liquid_solution`, `cp_liquid_solution`,
  `rho_cp_liquid_solution`, `k_liquid_solution`, `mu`, and `rho_gas`.
- `dissolved_gas`, `free_gas_mass`, `alpha_gas`, and `gas_mass_transfer`.

Current algorithm:
1. Refresh liquid oil plus dissolved-gas properties from pressure,
   temperature, and dissolved gas.
2. Feed field-valued `rho_liquid_solution` and `mu` into Gumbel/Elrod,
   velocity reconstruction, force post-processing, heat generation, and thermal
   convection.
3. Solve the Reynolds/JFO step with `theta` still acting as the liquid-content
   variable.
4. For `GAS_CAVITATION_MIXTURE`, transport free gas with the film velocity and
   apply a finite-rate release/resorption update between dissolved gas and free
   gas.
5. Cap the reported `alpha_gas` by `gas_alpha_max`; do not delete free gas just
   because a cell is locally full-film under JFO.

This is deliberately a segregated JFO nonlinear corrector plus property/gas
coupling loop. It is not labelled PISO because the codebase does not yet have a
Navier-Stokes pressure-velocity backend.

Tests added/updated:
- `test_fluid_properties`: config aliases, solubility monotonicity, mixture
  density bounds, viscosity selection/table interpolation, the calibrated
  R290/PZ68 40 C / 10 bar point, saturated zero source, release, resorption,
  and bounded `alpha_gas`.
- `test_gas_transport`: dissolved-gas transport constant preservation,
  disabled-model and zero-timestep no-ops, bounded advection, and
  mass-weighted conservation.
- `test_io`: legacy `load_x` output selection writes `pressure_force_x` without
  a duplicate `load_x` array.
- `test_forces`: force diagnostics no longer require legacy `load_x`.

Deferred work:
- Boundary-condition-rich dissolved/free gas transport beyond the first axial
  open-boundary and pressure-inlet reservoir-composition contract; future work
  should add selectable gas-specific BCs where a case genuinely needs them.
- Additional empirical propane/oil table calibration for temperatures and
  pressures beyond the supplied chart point.
- Bubble-radius dynamics, slip, and a true pressure-velocity PISO backend.

Follow-up correction (2026-06-08): added shared CLI/GUI validation for
absolute-pressure operation and gas-model dependencies. Solver runs now fail
before launch when `p_cav` is left at gauge zero, when gas cavitation is paired
with a non-JFO cavitation model, or when required Henry/density/viscosity/gas
transfer coefficients are missing. Added `config_r290_pz68_quick.txt` as a
small smoke case using the same R290/PZ68 constants.

Follow-up correction (2026-06-08): corrected the first gas/free-surface coupling
semantics. Free gas is no longer clipped to the JFO void fraction before
resorption, which silently deleted gas in pressure-inlet/full-film cells without
a negative `gas_mass_transfer` source. `free_gas_mass` is now conservatively
transported with the film-averaged velocity, and `alpha_gas` is capped only by
`gas_alpha_max` while resorption is controlled by the finite-rate transfer law.

Follow-up correction (2026-06-08): added the first open gas-composition
boundary contract. Pressure supply/inlet cells impose the configured incoming
mixture (`dissolved_gas_initial`, zero `free_gas_mass`, zero `alpha_gas`, zero
local `gas_mass_transfer`). Axial pressure-boundary inflow carries the same
initial dissolved-gas composition with zero free gas, and axial outflow convects
released gas out of the domain. Focused tests now guard inlet composition,
axial dissolved-gas inflow, and axial free-gas outflow.

---

## Phase A — Mass-Conserving Cavitation (Elrod-Adams)

Replace the Gumbel pressure clamp with the Elrod-Adams formulation. The universal variable $\theta$ (fractional film content) becomes the primary unknown. Mass is conserved across the cavitation boundary by construction.

### A.0 — Codebase Audit and Preparation

**Status: complete** (this document).

Findings:
- `EOS::dp_dtheta()` already exists but is unused — will become the core of $\gamma_{eff}$.
- `FVM::divergence` with upwind/TVD is implemented and tested — will be used for the Couette convection term.
- `FVM::ddt` is implemented — needed for transient $\partial(\theta h)/\partial t$.
- `LinearSystem` supports 5-point stencil with periodic $\theta$ — no changes needed.
- Current `Reynolds::solve` operates on pressure as primary variable — must be refactored.

---

### A.1 — Universal Variable Reformulation

**Goal:** Reformulate the Reynolds equation with $\theta$ (film content) as the primary variable. No cavitation modelling yet — full-film everywhere ($\theta \ge 1$). Validate that this reproduces the current Gumbel-based pressure field.

**Mathematical formulation:**

Starting from the compressible Reynolds equation:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta} + \frac{\partial(\rho h)}{\partial t}$$

Substitute $\rho = \rho_0 \theta$ and apply the chain rule $\partial p / \partial x_i = (dp/d\theta)\,\partial\theta/\partial x_i$:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\Gamma \frac{\partial\theta}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\Gamma \frac{\partial\theta}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho_0 \theta h)}{\partial\theta} + \frac{\partial(\rho_0 \theta h)}{\partial t}$$

where the **effective diffusion coefficient** is:

$$\Gamma(\theta) = \frac{\rho_0 \theta \, h^3}{12\mu} \cdot \frac{dp}{d\theta} = \frac{\rho_0 \theta \, h^3}{12\mu} \cdot \frac{\beta}{\theta} = \frac{\rho_0 \beta \, h^3}{12\mu}$$

Note: for full-film ($\theta \ge 1$), $\Gamma$ is independent of $\theta$. This simplifies A.1 — the equation is linear in $\theta$.

The **Couette term** $\frac{\omega}{2}\frac{\partial(\rho_0 \theta h)}{\partial\theta}$ is a one-dimensional divergence of the flux $\frac{\omega R h}{2} \cdot \rho_0\theta$. In the FVM framework, this is a **convection** of $\theta$ with known face fluxes:

$$F_e^{(\theta)} = \frac{\omega R \, h_e}{2} \cdot \Delta z, \qquad F^{(z)} = 0$$

where $h_e$ is the film thickness interpolated to the east face. The `FVM::divergence` operator with upwind (or TVD) scheme handles this directly.

The **transient term** $\partial(\rho_0 \theta h)/\partial t = \rho_0 h \,\partial\theta/\partial t$ (for static $h$) is discretised as a weighted `ddt` with $h$ as a cell-volume multiplier.

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.hpp/cpp` | New function `Reynolds::solve_elrod(fields, sys, mesh, cfg)`. Keep `Reynolds::solve` as a fallback. |
| `config.hpp` | Add `enum CavitationModel { GUMBEL, ELROD_ADAMS }` and a config field to select it. |
| `main.cpp` | Add `"theta"` field. Initialise $\theta = 1$ (full film). Dispatch to the appropriate solver. |
| `fvm.hpp/cpp` | Add `FVM::ddt_weighted(sys, phi, weight, dt, mesh)` — like `ddt` but multiplies the cell volume by a per-cell weight field (for $h \cdot \partial\theta/\partial t$). |
| `equation_of_state.hpp` | No change — `dp_dtheta` already exists. |

**Solve sequence per timestep (A.1, full-film only):**

1. Compute $\Gamma = \rho_0 \beta h^3 / (12\mu)$ at physical + ghost cells.
2. Compute Couette face fluxes: $F_e(i,j) = \omega R \cdot h_e / 2 \cdot \Delta z$.
3. Reset system. Assemble:
   - `FVM::laplacian(sys, Gamma, mesh)` — Poiseuille diffusion
   - `FVM::divergence(sys, couette_flux_theta, zero_flux_z, theta, UPWIND, mesh)` — Couette convection
   - Dirichlet z-BC: $\theta = \theta(p_{cav}) = 1$ at $z = 0, L$
   - (Optional) `FVM::ddt_weighted(sys, theta, h_field, dt, mesh)` — transient
4. Solve $A\,\theta^{n+1} = b$.
5. Recover pressure: $p = p_{cav} + \beta \ln(\max(\theta, 1))$.
6. Update density: $\rho = \rho_0 \theta$.

**FVM sign convention check:**

The assembled equation is:

$$\underbrace{-\text{div}(\Gamma\nabla\theta)\cdot V}_{\text{laplacian}} + \underbrace{\text{div}(\mathbf{F}\,\theta)\cdot V}_{\text{divergence}} = \text{source}$$

The PDE $\text{div}(\Gamma\nabla\theta) = \text{div}(\mathbf{F}\,\theta) + \partial(\rho_0\theta h)/\partial t$ rearranges to:

$$-\text{div}(\Gamma\nabla\theta)\cdot V + \text{div}(\mathbf{F}\,\theta)\cdot V = \rho_0 h \frac{\partial\theta}{\partial t}\cdot V$$

The LHS matches. The transient goes to the RHS and is handled by `ddt_weighted` (which adds $h V/\Delta t$ to $a_P$ and $h V \theta^n / \Delta t$ to `source`).

**Validation:**

- **Regression test**: Run the static journal bearing from Phase 3. The $\theta$-formulation must produce the same pressure field as the Gumbel solver (since $\theta \ge 1$ everywhere in full-film regions; the Gumbel clamp has no effect when pressure is already above $p_{cav}$). Compare $\max|p_{Elrod} - p_{Gumbel}|$ < tolerance.
- **Unit test**: 1D channel with uniform $h$, $\Gamma = \text{const}$. Exact solution for $\theta(x)$ is linear between Dirichlet BCs. Verify 2nd-order convergence.

---

### A.2 — Elrod-Adams Switch Function and Cavitation

**Status: complete** (2026-04-12).

**Goal:** Enable the cavitation switch $g(\theta)$. In cavitated regions ($\theta < 1$), the Poiseuille (diffusion) term vanishes and only the Couette (convection) term transports liquid. This is the core of mass-conserving cavitation.

**Modified diffusion coefficient:**

$$\Gamma_g(\theta) = \frac{\rho_0 \beta \, h^3}{12\mu} \cdot g(\theta)$$

where the switch function $g$ satisfies the **Elrod-Adams complementarity conditions**:

$$g \in \{0, 1\}, \qquad \theta \ge 1, \qquad g(\theta - 1) = 0, \qquad p \ge p_{cav}$$

These imply:
- $g = 1 \Longleftrightarrow \theta \ge 1$ (full film, pressure active)
- $g = 0 \Longleftrightarrow \theta < 1$ (cavitated, pressure clamped to $p_{cav}$)

**Switch-function update strategy:**

Implement an **iterative flag update** with safeguards (Vijayaraghavan & Keith, 1989):

```
for outer_iter = 1 to max_outer:
    1. Compute g(i,j) from current theta:
       g(i,j) = 1  if theta(i,j) >= 1
       g(i,j) = 0  if theta(i,j) < 1
    
    2. Compute Gamma_g = rho0 * beta * h^3 / (12*mu) * g
    
    3. Assemble and solve for theta_new (same as A.1 but with Gamma_g)
    
    4. Post-process:
       - In full-film cells (g=1): if theta_new < 1, mark for switch
       - In cavitated cells (g=0): if theta_new > 1, mark for switch
    
    5. Check convergence:
       - No cells changed flag, AND
       - ||theta_new - theta_old|| / ||theta_new|| < tol
    
    6. Update theta = theta_new
```

This is simpler than Fischer-Burmeister or Murty pivoting. It works well for steady and mildly transient problems. If convergence is poor, Phase A.3 adds a more robust complementarity solver.

**Regularisation:**

To avoid division-by-zero in the cavitated region and improve convergence:
- $g_{reg}(\theta) = \tfrac{1}{2}(1 + \tanh((\theta - 1)/\epsilon))$ with $\epsilon \sim 10^{-3}$ can be used as a smooth approximation during initial iterations, then tightened to a hard switch.

**Boundary conditions in $\theta$-formulation:**

| Boundary | Current (pressure) | $\theta$-formulation |
|----------|-------------------|---------------------|
| $z = 0, L$ (Dirichlet) | $p = p_{cav}$ | $\theta = 1$ (since $\theta(p_{cav}) = 1$) |
| $\theta$ (periodic) | Ghost exchange for $p$ | Ghost exchange for $\theta$ |
| Supply groove (future) | $p = p_{supply}$ | $\theta = \exp((p_{supply} - p_{cav})/\beta)$ |

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Add outer iteration loop around the linear solve. Compute `g` field from $\theta$. Multiply $\Gamma$ by `g`. |
| `config.hpp` | Add `max_outer_iters`, `outer_tol`, `g_regularisation_eps`. |
| `equation_of_state.hpp` | Add `switch_function(theta)` and `switch_function_regularised(theta, eps)`. |
| Tests | New `test_elrod.cpp`. |

**Validation:**

- **1D inclined slider bearing** (Elrod & Adams, 1975): A linearly converging channel $h(x) = h_1 + (h_2 - h_1) x / L$ with $h_1 > h_2$. The analytical JFO solution gives:
  - Pressure rise in the converging zone ($\theta = 1$, $p > p_{cav}$)
  - Cavitation onset at $\theta = \theta^*$ where $\partial p/\partial x = 0$
  - Cavitated zone ($\theta < 1$, $p = p_{cav}$) with linear $\theta$ decrease
  - Reformation where $\theta$ returns to 1

  Validate: compare computed cavitation boundary position, peak pressure, and $\theta$ profile against the analytical solution.

- **Journal bearing**: Compare cavitation extent (angular range where $\theta < 1$) against published data (e.g., Vijayaraghavan & Keith 1989, Fesanghary & Khonsari 2011). Verify that the sum $\int \rho_0 \theta h \, dA$ is conserved to machine precision between timesteps (mass conservation check).

- **Mass conservation diagnostic**: At the end of each solve, compute and print:
  $$M_{error} = \left|\sum_f F_f \theta_f - \text{(transient contribution)}\right| \bigg/ \sum_f |F_f|$$
  This should be $O(10^{-12})$ (solver tolerance), not $O(1)$ as with Gumbel.

---

### A.3 — Robust Complementarity Solver (if needed)

**Goal:** If the iterative flag update from A.2 shows poor convergence (oscillating flags, slow outer iteration), implement a more robust complementarity strategy.

**Status: skipped** (A.2 converges well for current cases).

---

### A.4 — Text-based Configuration

**Goal:** Move simulation parameters from a hardcoded header to a text-based config file to avoid constant recompilation.

**Proposed format:** Simple key-value pairs (`key = value`).

**Code changes:**
| File | Change |
|------|--------|
| `config.hpp` | Add `load_from_file(std::string path)` method to `SimulationConfig`. |
| `main.cpp` | Parse command line for config path (default: `config.txt`). |

---

### A.5 — Fix Numerical Starvation (Cascade Collapse)

**Status: complete** (Flooded initialization implemented in solve_elrod).

---

### A.6 — Flexible Inlet Boundary Conditions

**Goal:** Allow users to define internal oil supply features (holes and grooves) via config.

**Status: complete** (Implemented via penalty method in linear assembly).

---

### A.7 — Staggered Grid Velocity Fields

**Status: complete** (Implemented in Reynolds namespace, interpolated for VTK output).

---

### A.8 — Configurable Outlet Boundary Conditions

**Goal:** Allow users to choose between Fixed Value (Dirichlet) and Zero Gradient (Neumann) at the axial ends.

**Status: complete** (Implemented in Reynolds solvers and velocity calculation).

---

## Phase B — Dissolved and Free Gas Transport

Model gaseous cavitation: dissolved gas comes out of solution when pressure drops, forming free gas bubbles that alter the mixture properties.

### B.1 — Fluid Property Model (Mixture EOS)

**Goal:** Replace the single-component barotropic EOS with a two-phase mixture model that accounts for liquid, dissolved gas, and free gas.

**Mixture state variables (per cell):**

| Variable | Symbol | Description |
|----------|--------|-------------|
| Liquid volume fraction | $\alpha_l$ | $= \theta$ from Phase A (for cavitated cells, $\alpha_l < 1$) |
| Free gas volume fraction | $\alpha_g$ | Bubbles in the gap |
| Dissolved gas mass fraction | $c_d$ | Gas dissolved in the liquid (by mass) |

Constraint: $\alpha_l + \alpha_g = 1$ (two-phase mixture; vapour cavitation is already handled by $\theta$).

**Mixture density:**

$$\rho_{mix} = \alpha_l \rho_l(p) + \alpha_g \rho_g(p)$$

where:
- $\rho_l(p) = \rho_{l,0} \exp((p - p_{ref})/\beta_l)$ — liquid (existing barotropic EOS)
- $\rho_g(p) = p / (R_g T)$ — ideal gas (isothermal)

**Mixture viscosity** (Einstein or Krieger-Dougherty model for bubbly flow):

$$\mu_{mix} = \mu_l (1 - \alpha_g)^{-[\mu] \alpha_{g,max}}$$

where $[\mu] \approx 2.5$ (Einstein) and $\alpha_{g,max} \approx 0.64$. For dilute gas ($\alpha_g \ll 1$), this simplifies to $\mu_{mix} \approx \mu_l (1 + 2.5\alpha_g)$.

**Henry's law** (dissolved gas equilibrium):

$$c_{d,eq}(p) = H \cdot p$$

where $H$ is the Henry's law constant (units: $\text{kg}_{gas}/(\text{kg}_{liquid}\cdot\text{Pa})$). At equilibrium, the dissolved gas concentration equals $c_{d,eq}$. When $p$ drops below saturation, dissolved gas is released as free gas.

**Code changes:**

| File | Change |
|------|--------|
| `src/equation_of_state.hpp` | Refactor into `src/fluid_properties.hpp`. Add `MixtureEOS` class with methods: `rho_liquid(p)`, `rho_gas(p, T)`, `rho_mix(alpha_g, p)`, `mu_mix(alpha_g)`, `henry_equilibrium(p)`. Keep the old `EOS` namespace as a thin wrapper for backward compatibility. |
| `config.hpp` | Add gas properties: $R_g$, $T$, $H$, $\mu_g$, $\alpha_{g,max}$. |
| Tests | `test_mixture_eos.cpp`: round-trip density, viscosity limits ($\alpha_g = 0 \to$ pure liquid, $\alpha_g \to \alpha_{g,max} \to$ high viscosity), Henry's law. |

**Validation:**
- Pure liquid limit ($\alpha_g = 0$): must recover Phase A results exactly.
- Pure gas limit ($\alpha_l = 0$): ideal gas EOS.
- Mixture density monotonicity in $p$.

---

### B.2 — Dissolved Gas Transport Equation

**Goal:** Add a transport equation for the dissolved gas mass fraction $c_d$ in the liquid phase.

**Governing equation:**

$$\frac{\partial}{\partial t}(\alpha_l \rho_l c_d h) + \text{div}(\alpha_l \rho_l c_d \mathbf{u}_l h) = -\dot{m}_{dg} h + D_{eff} \nabla^2(\alpha_l \rho_l c_d h)$$

where:
- $\mathbf{u}_l$ is the mean liquid velocity (Couette + Poiseuille contributions)
- $\dot{m}_{dg}$ is the mass transfer rate from dissolved to free gas (degassing source)
- $D_{eff}$ is an effective diffusion coefficient (molecular + turbulent mixing)

**Degassing/resorption source term:**

$$\dot{m}_{dg} = k_{dg} \, \alpha_l \rho_l \, (c_d - c_{d,eq}(p))$$

where $k_{dg}$ is a rate constant ($s^{-1}$). When $c_d > c_{d,eq}$ (supersaturated), gas is released ($\dot{m}_{dg} > 0$). When $c_d < c_{d,eq}$ (undersaturated), gas is resorbed ($\dot{m}_{dg} < 0$).

**Simplification for thin films:** In the lubrication limit, the mean liquid velocity in the $\theta$-direction is $\bar{u}_\theta = \omega R / 2 - h^2/(12\mu) \cdot \partial p / \partial s$. The convective flux uses this velocity.

**Code changes:**

| File | Change |
|------|--------|
| `src/gas_transport.hpp/cpp` | New module: `GasTransport::solve_dissolved(fields, sys, mesh, cfg)`. Assembles and solves the $c_d$ transport equation using existing FVM operators. |
| `main.cpp` | Add `"c_dissolved"` field. Call `GasTransport::solve_dissolved` after Reynolds solve. |
| `config.hpp` | Add `k_degassing`, `D_eff`, `c_d_initial`. |
| `fvm.hpp/cpp` | May need a general convection-diffusion helper, or assemble manually from existing operators. |

**Validation:**
- 1D advection test: uniform flow with known inlet $c_d$, verify advection speed.
- Equilibrium test: set $p = \text{const}$ everywhere, verify $c_d \to H \cdot p$ at steady state.
- Conservation: total dissolved mass $\int \alpha_l \rho_l c_d h \, dA$ plus transferred mass $\int \dot{m}_{dg} h \, dA \cdot \Delta t$ must balance.

---

### B.3 — Free Gas Volume Fraction Equation

**Goal:** Evolve the free gas volume fraction $\alpha_g$ based on degassing from the liquid and gas compressibility.

**Governing equation:**

$$\frac{\partial}{\partial t}(\alpha_g \rho_g h) + \text{div}(\alpha_g \rho_g \mathbf{u}_g h) = +\dot{m}_{dg} h$$

The free gas is transported by a mixture of Couette flow and pressure-driven expansion. For simplicity (no slip between gas and liquid), assume $\mathbf{u}_g = \mathbf{u}_l$ (homogeneous mixture).

**Code changes:**

| File | Change |
|------|--------|
| `src/gas_transport.hpp/cpp` | Add `GasTransport::solve_free_gas(fields, sys, mesh, cfg)`. Coupled to dissolved gas via $\dot{m}_{dg}$. |
| `main.cpp` | Add `"alpha_gas"` field. Initialise $\alpha_g = 0$ (no free gas initially). |

**Validation:**
- Mass conservation: total gas mass ($\int [\alpha_l \rho_l c_d + \alpha_g \rho_g] h \, dA$) must be constant in a closed domain (no flux BCs).
- Degassing test: uniform pressure drop at $t = 0$; verify exponential approach to equilibrium with rate $k_{dg}$.

---

### B.4 — Coupling to Reynolds Equation

**Goal:** Close the loop: the Reynolds equation uses mixture properties ($\rho_{mix}$, $\mu_{mix}$), and the gas transport equations use the pressure field from Reynolds.

**Coupled solve sequence per timestep:**

```
1. Sync ghost cells for all fields (theta, c_d, alpha_g, h)
2. Compute mixture properties:
   rho_mix(i,j) = (1 - alpha_g) * rho_l(p) + alpha_g * rho_g(p)
   mu_mix(i,j)  = mu_l * (1 - alpha_g)^(-2.5 * alpha_g_max)
3. Reynolds solve (Phase A solver, using rho_mix and mu_mix):
   Gamma_g = beta * h^3 / (12 * mu_mix) * g(theta)
   → solve for theta → recover p
4. Dissolved gas transport (B.2):
   → solve for c_d
5. Free gas transport (B.3):
   → solve for alpha_g
6. Check outer convergence:
   max(|theta_new - theta_old|, |c_d_new - c_d_old|, |alpha_g_new - alpha_g_old|) < tol
7. Repeat 2-6 if not converged (Picard outer iteration)
8. Store old-time values, advance timestep
```

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Accept `rho_mix` and `mu_mix` fields instead of scalar `cfg.rho` and `cfg.mu`. |
| `main.cpp` | Implement the coupled outer iteration loop. |
| `config.hpp` | Add `max_coupling_iters`, `coupling_tol`. |

**Validation:**
- **Decoupled limit**: set $\alpha_g = 0$ and $c_d = H \cdot p$ everywhere, disable gas transport. Must recover Phase A results exactly.
- **Journal bearing with gas release**: compare cavitation extent and peak pressure against published two-phase results (e.g., Ausas et al. 2009, Bayada & Chupin 2013).
- **Convergence study**: verify that outer iteration converges in $\le 10$ iterations for typical bearing conditions.

---

## Phase C - Moving Bearing Dynamics

**Status: complete** (2026-06-04, moving outer bearing with fixed shaft).

Allow the outer bearing surface to move in response to pressure and viscous
film forces, independent applied load, support stiffness, and support damping.
This turns the current static-eccentricity solve into a lagged fluid-structure
time integration while keeping the shaft fixed in space.

### C.1 - Time-dependent Film Geometry

**Goal:** Replace static eccentricity inputs with bearing-center state variables
when dynamic motion is enabled.

For bearing-center displacement $(x_b(t), y_b(t))$ measured from the fixed shaft
center, the nominal film thickness becomes:

$$h(\theta,t) = c + x_b(t)\cos\theta + y_b(t)\sin\theta$$

with optional existing misalignment terms retained as additive axial tilt:

$$h(\theta,z,t) = c + x_b(t)\cos\theta + y_b(t)\sin\theta
                 - (z-L/2)\alpha_x\cos\theta
                 - (z-L/2)\alpha_y\sin\theta$$

When `bearing_initial_from_attitude = true`, the initial moving-bearing
position is $x_b=-e\cos\psi,\ y_b=-e\sin\psi$, so the moving model starts from
the same gap as the static shaft-eccentricity model. The equivalent attitude
angle is updated from the bearing state as $\psi=\operatorname{atan2}(-y_b,-x_b)$.

The time derivative needed by the Reynolds equation is:

$$\frac{\partial h}{\partial t} =
    \dot{x}_b\cos\theta + \dot{y}_b\sin\theta$$

**Code changes:**

| File | Change |
|------|--------|
| `config.hpp` | Added `motion_model`, bearing initial state, support stiffness/damping, bearing mass, independent external load, safety limits, and time-method keys. |
| `journal_motion.hpp/cpp` | Added bearing state, force-balance integrators, attitude conversion, and uniform output field helpers. |
| `film_thickness.hpp/cpp` | Added moving-bearing film-thickness update and `dh_dt`. |
| `main.cpp` | Stores bearing state over time; updates `h` and `dh_dt` before each Reynolds solve; advances bearing motion from fluid plus external force. |
| `io.cpp` | Uses existing field-selection path for bearing position, external load, and result-force scalar fields. |
| `gui_win32.cpp` | Added clean Motion / Loads controls for the moving-bearing and time-method config keys, plus magnitude/direction input for the independent external load. |

**Validation:**
- Static limit: zero velocity and fixed journal state must reproduce the
  current static-eccentricity pressure, load, and torque.
- Initial attitude: moving-bearing initialization from `e` and
  `attitude_angle_deg` must exactly reproduce the static film thickness.
- Safety: runs must stop or reject input before $h_{min} \le 0$.

---

### C.2 - Transient Reynolds Coupling

**Goal:** Activate the transient squeeze-film contribution created by moving
geometry.

For the Elrod variable, the transient term expands as:

$$\frac{\partial(\rho_0\theta h)}{\partial t}
  = \rho_0 h\frac{\partial\theta}{\partial t}
  + \rho_0\theta\frac{\partial h}{\partial t}$$

`FVM::ddt_weighted(sys, theta, rho*h, dt, mesh)` handles the first term for
all transient moving-bearing pressure time methods. This capacity term is
mandatory for liquid-content storage; even `EULER_EXPLICIT` keeps it in the
linear solve. The second term is assembled as an explicit source using the old
or Picard-lagged value of $\theta$:

$$S_h = \rho_0\theta^n\frac{\partial h}{\partial t}$$

**Code changes:**

| File | Change |
|------|--------|
| `reynolds.cpp` | Added optional transient assembly for dynamic `h`, including the explicit $\theta\,\partial h/\partial t$ source for Elrod and an explicit squeeze source for Gumbel. |
| `field.hpp/cpp` | `dh_dt` is stored as a regular scalar field with ghost updates through existing infrastructure. |
| Tests | Added `test_journal_motion` for dynamic film-thickness sign, attitude recovery, config aliases, and integrators. |

Follow-up correction (2026-06-07): `pressure_time_method = EULER_EXPLICIT` no
longer skips the Elrod weighted capacity term. It only leaves the
$\theta\,\partial h/\partial t$ squeeze source lagged. This prevents
post-startup artificial cavitation caused by applying squeeze without the
matching old liquid-content storage.

Startup correction (2026-06-07): the earlier forced steady `t=0` transient
startup was removed. Transient runs now begin from the initialized uniform
outlet-pressure/film-content fields and use transient terms from the first
solve. `omega_ramp_time` provides the intended no-rotation startup path by
ramping `omega` linearly from zero to the target speed over the configured
timespan. The Windows GUI exposes this ramp plus advanced fluid-property and
gaseous-cavitation controls.

**Validation:**
- Squeeze-film check: with $\omega = 0$ and prescribed radial approach,
  pressure must be generated only by $\partial h/\partial t$.
- Time-step refinement: transient load should converge as `dt` is reduced.

---

### C.3 - Bearing Equation of Motion

**Goal:** Couple pressure and viscous fluid force back to bearing motion.

The first implementation uses a three-degree-of-freedom rigid bearing support:

$$m_b\ddot{x}_b + c_x\dot{x}_b + k_x x_b =
  F_x^{fluid} + F_x^{external}$$

$$m_b\ddot{y}_b + c_y\dot{y}_b + k_y y_b =
  F_y^{fluid} + F_y^{external}$$

$$m_b\ddot{z}_b + c_z\dot{z}_b + k_z z_b =
  F_z^{fluid} + F_z^{external}$$

where $\mathbf{F}^{fluid}$ is the equal-and-opposite force to the integrated
pressure plus viscous force on the shaft surface. It is saved as
`fluid_force_x`, `fluid_force_y`, and `fluid_force_z`. The independent applied
load is saved separately as `external_load_x`, `external_load_y`, and
`external_load_z`.

**Implementation strategy:**
1. Use lagged fluid force from the latest Reynolds solve.
2. Support `EULER_EXPLICIT`, `EULER_IMPLICIT`, `CRANK_NICOLSON`, `RK2`, and
   `RK4` for motion. `CRANK_NICOLSON` is the formal name for the mixed
   implicit/explicit trapezoidal method; `SEMI_IMPLICIT` and `IMEX` are accepted
   aliases in config parsing.
3. Add a Picard loop around bearing motion and Reynolds solve when lagging is
   too weak for stiff supports or small clearances.

**Code changes:**

| File | Change |
|------|--------|
| `src/journal_motion.hpp/cpp` | New module for bearing state, force balance, and time integration. |
| `main.cpp` | Added motion/Reynolds coupling sequence. |
| `config.hpp` | Added motion-model selection: `STATIC`, `MOVING_BEARING`. |
| Tests | Added static-equivalent geometry, config alias, and constant-load integrator checks. |

**Validation:**
- Static equilibrium: with large support stiffness, the dynamic solve should
  converge to the configured static eccentricity.
- Free response: with hydrodynamic force disabled, the journal ODE must match
  the analytic damped oscillator.
- Linearized bearing coefficients: small perturbations about equilibrium should
  produce repeatable stiffness/damping estimates.

---

## Phase D - Thermal and ETHD Implementation

Add temperature, viscosity variation, heat generation, heat transfer, and elastic
film-thickness deformation. This phase should be introduced incrementally:
thermal hydrodynamics (THD) first, then elastic deformation, then fully coupled
elasto-thermo-hydrodynamics (ETHD).

### D.1 - Energy Equation and Temperature Field

**Goal:** Add a lubricant temperature field and solve a thin-film energy balance.

**Status:** First constant-property implementation is complete. Remaining D.1
work is validation breadth and optional higher-order transient thermal schemes.

A first-order averaged film model is:

$$\rho c_p h\left(\frac{\partial T}{\partial t}
  + \bar{u}_\theta\frac{1}{R}\frac{\partial T}{\partial\theta}
  + \bar{u}_z\frac{\partial T}{\partial z}\right)
  =
  \frac{1}{R^2}\frac{\partial}{\partial\theta}
  \left(k h\frac{\partial T}{\partial\theta}\right)
  + \frac{\partial}{\partial z}
  \left(k h\frac{\partial T}{\partial z}\right)
  + \Phi - Q_w$$

where $\Phi$ is viscous dissipation and $Q_w$ is heat loss to journal and
bearing walls. A practical initial dissipation model is:

$$\Phi \approx \frac{\mu U^2}{h}
       + \frac{h^3}{12\mu}\left[
         \left(\frac{1}{R}\frac{\partial p}{\partial\theta}\right)^2
         + \left(\frac{\partial p}{\partial z}\right)^2\right]$$

**Code changes:**

| File | Change |
|------|--------|
| `src/energy.hpp/cpp` | Added constant-property energy-equation module using existing FVM diffusion/convection operators. |
| `main.cpp` | Added `temperature` and `heat_generation` fields and calls the energy solve after velocity/force updates. |
| `config.hpp` | Added `solution_mode`, `temperature_model`, `temperature_initial`, `temperature_reference`, `rho_cp`, `thermal_conductivity`, wall heat-transfer coefficients, and wall temperatures. |
| `gui_win32.cpp` | Added Energy controls, steady-state solution selection, and preview units for thermal fields. |
| Tests | Added `test_energy` for config aliases, isothermal heat refresh, and steady Couette wall balance. |

**Validation:**
- Isothermal limit: disabling energy coupling must reproduce Phase A/B results.
- Pure diffusion: no velocity and no heat source must match a linear or
  manufactured temperature solution.
- Energy balance: generated heat minus wall losses must match the change in
  integrated thermal energy.

---

### D.2 - Thermoviscous Reynolds Coupling

**Goal:** Feed temperature-dependent viscosity and density back into Reynolds.

Start with a configurable viscosity law:

$$\mu(T,p) = \mu_{ref}
  \exp[-a_T(T - T_{ref}) + a_p(p - p_{ref})]$$

The Reynolds diffusion coefficient then becomes cell-wise:

$$\gamma = \frac{\rho(T,p)h^3}{12\mu(T,p)}$$

and the Elrod coefficient becomes:

$$\Gamma_g = \frac{\rho\beta h^3}{12\mu(T,p)}g(\theta_{film})$$

**Code changes:**

| File | Change |
|------|--------|
| `src/fluid_properties.hpp/cpp` | Extend the planned mixture-property module with viscosity-temperature and optional Barus pressure-viscosity laws. |
| `reynolds.cpp` | Accept field-valued `rho` and `mu` for both Gumbel and Elrod paths. |
| `main.cpp` | Add THD Picard loop: update properties, solve Reynolds, solve energy, check convergence. |

**Validation:**
- Constant-property limit recovers the current solver exactly.
- Hotter lubricant should reduce viscosity and lower generated pressure for
  the same geometry and speed.

---

### D.3 - Elastic Film Deformation

**Goal:** Add elastic and thermal deformation to the effective film thickness.

The effective film thickness is:

$$h_{eff} = h_{kinematic} + \delta_p + \delta_T$$

where $\delta_p$ is pressure-induced elastic deflection and $\delta_T$ is
thermal expansion. The first implementation should support a local linear
compliance model:

$$\delta_p(\theta,z) = C_p(\theta,z)\,[p(\theta,z)-p_{ref}]$$

and a local thermal expansion model:

$$\delta_T(\theta,z) = \alpha_s\,t_s\,[T_s(\theta,z)-T_{ref}]$$

Later versions can replace this with a convolution influence matrix or an
external structural solve without changing the Reynolds interface.

**Code changes:**

| File | Change |
|------|--------|
| `src/deformation.hpp/cpp` | New module for elastic and thermal deflection fields. |
| `film_thickness.cpp` | Combine kinematic, elastic, and thermal film-thickness contributions. |
| `config.hpp` | Add elastic compliance, thermal expansion, solid thickness, and deformation toggles. |
| Tests | Local-compliance response, thermal expansion response, and zero-compliance regression. |

**Validation:**
- Zero compliance and zero thermal expansion recover THD results.
- Positive pressure compliance must increase local film thickness where the
  configured sign convention says the bearing surface moves away from the
  journal.
- Deformation should converge under mesh refinement for smooth pressure fields.

---

### D.4 - Fully Coupled ETHD Loop

**Goal:** Couple Reynolds, energy, fluid properties, and deformation.

Recommended timestep sequence:

```
1. Update journal state and kinematic film thickness
2. Initialize Picard iteration from old p, theta, T, mu, h_eff
3. Compute deformation from pressure and temperature
4. Update effective film thickness
5. Update fluid properties from T, p, and gas state
6. Solve Reynolds / Elrod equation
7. Compute velocities and viscous heat generation
8. Solve energy equation
9. Check coupled residuals for h_eff, p/theta, T, and mu
10. Repeat 3-9 until converged or max iterations reached
11. Compute load, torque, diagnostics, and advance time
```

**Code changes:**

| File | Change |
|------|--------|
| `main.cpp` | Add coupled ETHD Picard driver with configurable convergence tolerances. |
| `config.hpp` | Add `enable_thermal`, `enable_deformation`, `max_ethd_iters`, and `ethd_tol`. |
| Tests | Decoupled-limit regression and coupled manufactured/sanity cases. |

**Validation:**
- Decoupled limit: disabling thermal and deformation coupling must recover the
  current Elrod/Phase B result.
- THD-only and EHD-only toggles must be independently usable.
- Coupled residuals should decrease monotonically for moderate load/speed test
  cases or report a clear non-convergence diagnostic.

---

## Phase Summary

| Phase | Description | Primary Variable | Key Output |
|-------|-------------|------------------|------------|
| 1-3 | Infrastructure + Gumbel Reynolds | $p$ | Pressure field, basic cavitation |
| A.1 | $\theta$-reformulation (full-film only) | $\theta$ | Equivalent to Phase 3 |
| A.2 | Elrod-Adams switch + cavitation | $\theta$ with $g(\theta)$ | Mass-conserving cavitation |
| A.3 | Robust complementarity (if needed) | $\theta$ | Improved convergence |
| B.1 | Mixture EOS | — | $\rho_{mix}$, $\mu_{mix}$ |
| B.2 | Dissolved gas transport | $c_d$ | Gas release/resorption |
| B.3 | Free gas volume fraction | $\alpha_g$ | Bubble evolution |
| B.4 | Coupled solver | $\theta, c_d, \alpha_g$ | Full two-phase bearing |
| C.1-C.3 | Moving bearing dynamics | $x_b, y_b, z_b, \dot{x}_b, \dot{y}_b, \dot{z}_b$ | Fluid force, external load, squeeze-film motion |
| D.1 | Energy equation | $T$ | Temperature and heat-generation fields |
| D.2-D.4 | Thermal coupling and ETHD solver | $\mu(T,p), h_{eff}$ | Temperature-dependent viscosity, deformation, coupled ETHD response |

---

## Phase E — Audit Remediation (2026-06-11; docs consolidated 2026-06-14)

A two-agent physics/plan review, a consolidated verdict, and a
literature-grounded audit (20+ primary sources with quoted claims) identified
model-level and numerical defects. **Their full content — the agent-executable
work packages with physics, citations, code touch-points, and test acceptance
gates — is inlined under "Phase E Work Package Specifications" below.** These
were previously split across `AUDIT.md` / `AUDIT_PLAN.md` / `REVIEW.md` /
`REVIEW_VERDICT.md`, now removed; the standalone originals are preserved in git
history at commit `ae567d6`. The `§A`/`§B`/`§C`/`§D` labels in the "Closes"
lines are the original audit's finding IDs.

WP-1/WP-2 complete what Phase B.4 originally specified (mixture density into
Reynolds + Picard coupling loop); the B.4 implementation that landed left free
gas volumetrically inert, which is the audit's headline defect. **WP-12 (new,
2026-06-14)** adds the 2-D-from-file property-table backend and the
`p_sat`-driven cavitation onset; its physics basis is
`docs/CAVITATION_OIL_REFRIGERANT.md`.

| WP | Title | Closes (AUDIT) | Status |
|----|-------|----------------|--------|
| WP-4 | Conservation/regime/convergence diagnostics + doc reconciliation | Part V.7, §B1 | **DONE** (2026-06-11) |
| WP-3 | Numerics honesty (type-differencing/TVD, dead time selectors, fixed-bearing transient capacity term) | §C1, §C5 | **DONE** (2026-06-11) |
| WP-5 | Mixture viscosity with gas-limit asymptote; solubility tables | §A3, §A4 | **DONE** (2026-06-11; measured isotherm table data still pending from user) |
| WP-2 | Outer Picard coupling; converged steady state; steady+gas guard | §C4 | **DONE** (2026-06-14); `test_coupling`; transient ≡ master, steady iterates to convergence |
| WP-1 | Couple released gas into film continuity/pressure; onset at p_sat(T,c_d) | §A1, §A2 | **DONE** (2026-06-14, Stage 1) `gas_pressure_coupling=VOID_COUPLED`; θ_full/p_void/β̄; `test_eos`; NONE ≡ master, VOID_COUPLED stable (load shift verified) |
| WP-6 | θ-weighted energy capacity/advection and cavitated-zone shear | §A5 | **DONE** (2026-06-14) `cavitated_film_model=STRIATED`; Couette shear/heat-gen θ-weighted in cavitated cells; `test_coupling`; FULL_FILM ≡ master, STRIATED lowers torque |
| WP-7 | Quantitative validation campaign (Ausas/V&K/Giacopini/Grando/Ferron) | Part V | TODO (grid-convergence harness V7 seeded in `validation/grid_convergence/`) |
| WP-9 | Starved / mass-flux inlet | §B2 | TODO |
| WP-8 | Linearized K/C coefficients + whirl margin | §D1 | TODO |
| WP-10 | Strategic: FBNS/complementarity cavitation kernel | §C2, §C3 | TODO |
| WP-11 | Axial-boundary flux consistency (shared flux function; cavitated-cell reformation rate; INLET_OUTLET unification) | §B4 | **PARTIAL** (2026-06-14) `consistent_boundary_flux` (default off) gives cavitated boundary cells the physical reformation rate `cs·(p_bc−p_cav)/β`; default ≡ master. REMAINING: share one flux fn across energy/gas/**diagnostics** (the mass balance still uses the old link, so r_l rises when the flag is on) |
| WP-12 | 2-D property tables from an editable file (P,T grid) + `p_sat`-driven cavitation onset threshold | §A3/§A4 (data), §A1 (onset) | **DONE** (2026-06-14); `test_property_tables`; R290 table smoke run verifies p_sat=1.0 MPa onset |

### Phase E implementation notes (2026-06-11, WP-4/WP-3/WP-5)

- New `src/diagnostics.hpp/cpp`: per-step global liquid/gas mass balances that
  mirror the assembled operators, with clamp/snap mass accounting
  (`cavitation_clamp_mass`, `gas_clamp_mass` fields), inlet penalty cells
  re-evaluated as physical residuals, and `results/<run>/diagnostics.csv`
  output (interval: `diagnostics_interval`). Closed/open-domain residuals sit
  at `linear_rtol` (new key); the PLAN A.2 O(1e-12) target is now actually
  checked in `tests/test_diagnostics.cpp`. `validate()` adds Taylor-number and
  p_sat saturation-context warnings (`SimulationConfig::saturation_pressure`).
- WP-3: `theta/thermal/gas_convection_scheme` keys activate the previously
  dead TVD machinery; `TYPE_DIFFERENCING` is implemented for theta; the
  Elrod capacity term is assembled for **all** transient runs (fixed-bearing
  transients were silently quasi-steady); dead
  `pressure_time_method`/`temperature_time_method` members removed (parser
  warns); the Gumbel Couette term is a face-flux divergence matched to the
  Elrod discretization (matched full-film agreement ~0.01%, was ~1.2%).
  `tests/test_schemes.cpp` + `test_schemes_mpi6` cover boundedness, order,
  matching, relaxation, MPI ghosts, and back-compat.
- WP-5: `gas_mixture_viscosity_model` (EINSTEIN_DILUTE default, DUKLER_VOID,
  MCADAMS_QUALITY, KRIEGER_DOUGHERTY, LINEAR_QUALITY) + `mu_gas`; Einstein is
  rejected above the 0.6 packing bound; R290 configs ship MCADAMS_QUALITY.
  Isotherm `*_table` data remains a flagged user-supplied input.
- Test-suite adjustments forced by the capacity-gate fix: `test_elrod_a1` and
  `test_elrod_a2` now declare `solution_mode = STEADY_STATE` (their original
  quasi-steady semantics); transient film-content storage has dedicated
  coverage in `test_schemes`.

---


## Phase E — Work Package Specifications (consolidated)

## WP-12 — 2-D property tables from an editable file + `p_sat`-driven cavitation onset

**Closes:** the "isotherm table data pending" flag on WP-5, and §A1 onset
decoupling for the gas-laden case. **Severity:** medium (data fidelity) + the
user-facing "the higher-release species ruptures first" physics. **Depends on:**
none for steps 1–5; step 6 (threshold) is a default-off precursor to WP-1.
**Added** 2026-06-14 from `docs/CAVITATION_OIL_REFRIGERANT.md`.

### Physics
The supplier PTSV chart (`data/R290_PZ68S/r290_pz68s.csv` — a P×T grid of
saturated solubility and kinematic viscosity) *is* the liquid-solution property
surface. Solubility doubles as the saturation locus: inverting
$c_{sat}(p,T)=c_d$ at fixed $T$ gives the local bubble-point $p_{sat}(T,c_d)$,
the gaseous-cavitation onset. Whichever species releases first (higher pressure)
governs onset, so the effective threshold is
$p_{cav,\mathrm{eff}} = \max(p_{cav,\mathrm{oil}},\, p_{sat}(T,c_d))$ —
general, not R290-specific. For the shipped config $p_{cav}=3.0\times10^4$ Pa
but $p_{sat}\approx 1.0$ MPa, so propane outgasses first and the current Elrod
solve ignores it.

### Code (steps land in order; each keeps `ctest` green with defaults)
1. `PropertyTable2D` + bilinear `interpolate_2d` (non-uniform axes, edge-clamped) in `config.hpp`/`fluid_properties`.
2. `src/property_table_io.{hpp,cpp}`: header-mapped CSV loader; SI unit convert (MPa→Pa, °C→K, %→fraction, mm²/s→m²/s); axes reconstructed from the file (**resolution editable post-compile, no recompile**); on-load validation (headers, monotone axes, rectangular grid, finite).
3. Config keys `solubility_table_file` / `viscosity_table_file` / `density_table_file` + `config_dir` resolution; load in `validate()`; a 2-D file wins over the inline 1-D table (warn); a missing file is fatal.
4. Kinematic→dynamic viscosity per cell `μ = ν·ρ_solution` in `update_solution_fields` (ρ from the density model; PURE_OIL is a documented degraded mode; optional density file/column override).
5. Generalize `SimulationConfig::saturation_pressure` to invert the 2-D isotherm by **first-crossing scan** (the CSV is real VLE data: monotone-then-flat per isotherm with rare 1-LSB blips, some isotherms pinned at 100% — a binary search is wrong); shared `invert_csat()` for the 1-D and 2-D branches.
6. `cavitation_threshold = SCALAR_PCAV | LOCAL_PSAT` (default `SCALAR_PCAV` = current behaviour). A `p_cav_eff` field is populated from `max(p_cav, p_sat(c_d,T))` and threaded through `reynolds.cpp` via a `cav_pressure()` helper at the easy `p_cav` sites (`:114,:189-190,:233,:239,:485`); the force gauge datum (`:601`) stays the **scalar** `p_cav` (onset ≠ force reference); `p_cav_eff` is ghost-exchanged. Reconciles with WP-1: `p_void = max(p_cav_eff, m_g R_g T/((1-θ)h))`, so `p_cav_eff` is the gas-free floor of WP-1's plateau, not a competing threshold.

### Tests
New `tests/test_property_tables.cpp` (interp; CSV load incl. the README
spot-check (1.0 MPa, 40 °C) → 21.746 %; malformed-file errors). Extend
`test_fluid_properties` (2-D solubility/viscosity, ν→μ, 2-D `p_sat` round-trip +
flat-tail/pinned-isotherm/noise edges), `test_io` (file keys/precedence),
`test_forces` (`SCALAR_PCAV`≡master to 1e-10; `LOCAL_PSAT` onset≈`p_sat`; `-n 6`
rank parity). **Decoupled-limit:** defaults reproduce master bit-for-bit.

## Ground rules for every agent (read first)

1. The audit *why* (literature, defect IDs §A/§B/§C) is summarized in `PHYSICS.md` and preserved in git history at `ae567d6`; this section is the *how*.
2. Honor `CLAUDE.md` conventions: snake_case, minimal comments (headers
   self-documenting), modular/swappable models, **new physics must be
   disablable at the highest level** (config enum/flag, default = old
   behavior), `std::filesystem::path`, GCC/Clang/MSVC-buildable.
3. **Decoupled-limit regressions are sacred.** Every WP keeps the full existing
   `ctest` suite green. When you add physics, the "off" setting must reproduce
   the previous results bit-for-bit (or to solver tolerance — state which).
4. Test harness: add `tests/test_<name>.cpp`, register with
   `add_pancake_test(...)` in `tests/CMakeLists.txt` (runs under `mpiexec -n 2`
   by default; add an `-n 6` variant if your change touches ghost exchange, cf.
   `test_elrod_a2_mpi6`).
5. On completion: append a dated entry to `CHANGES.md`, update the relevant
   `PHYSICS.md` section with the model equations actually implemented, mark the
   WP status in `PLAN.md` Phase E, and update `README.md` if user-facing
   options changed.
6. Notation matches `PHYSICS.md`: θ = fractional film content, β = bulk
   modulus, `p_cav` = cavitation/plateau pressure, c_d = dissolved gas mass
   fraction, α_g = free-gas volume fraction, m_g = free gas mass per unit area
   (`free_gas_mass`), ρ_g = gas density, h = gap, ω = shaft speed, R = radius.

**Dependency graph (do in this order unless parallelizing):**

```
WP-4 (diagnostics)  ──────────────┐  (no deps; gives observability for all)
WP-3 (numerics honesty) ──────────┤  (no deps)
WP-5 (property models) ──┐        │
WP-2 (outer coupling) ───┼──> WP-1 (two-phase closure) ──> WP-7 (validation, gas cases)
WP-6 (θ-weighted energy/shear)────┘        WP-2 ──> WP-8 (K/C extraction)
WP-9 (starved inlet)  (independent)
WP-10 (FBNS kernel)   (strategic; after WP-1 stabilizes)
WP-7 (validation, pure-oil cases V1–V4) can start immediately.
```

---

## WP-1 — Couple released gas into the film continuity / pressure problem

**Closes:** AUDIT §A1, §A2 (the two "contradicts-the-literature" defects).
**Severity:** MAJOR. **Depends on:** WP-2 (outer loop), WP-5 (mixture μ); WP-4
strongly recommended first (conservation logs make debugging this possible).

### Physics

**The defect.** Today `ρ = ρ_l·θ` (`fluid_properties.cpp:216`) and
`reynolds.cpp` never reads `alpha_gas`/`free_gas_mass`; the outgassing source
`Δṁ_g` updates diagnostics only. Released gas can neither open void nor relieve
pressure — the defining mechanism of gaseous cavitation
(Braun & Hannon 2010, DOI 10.1243/13506501JET772; Etsion & Ludwig 1982, ASME
J. Lubr. Technol. 104(2):157–163). Published magnitude of what this suppresses:
load −3% at ε = 0.3, **+21% at ε = 0.9** (Int. J. Refrigeration 2023,
S0140700723004577).

**Target closure (Σα = 1).** Volume fractions liquid/gas must satisfy

$$\alpha_l + \alpha_g = 1,\qquad \alpha_l \equiv \theta,\qquad
\alpha_g = \frac{m_g}{\rho_g(p,T)\,h},\qquad
\rho_g = \frac{p}{R_g T}.$$

Mixture density and the film mass per unit area become

$$\bar\rho = \rho_l\,\theta + \rho_g\,\alpha_g, \qquad
M = \bar\rho\,h = \rho_l\theta h + m_g .$$

(Homogeneous-mixture precedent: Grando, Priest & Prata 2006, Tribology Letters,
DOI 10.1007/s11249-006-9027-6 — "ρ̄ = φρ_g + (1−φ)ρ_l … enter the Reynolds
equation directly"; same closure in Int. J. Refrig. 2023 and in the
solubility-based compressible Reynolds equation of ASME J. Tribology
147(2):024101, DOI 10.1115/1.4066414.)

**Saturation pressure (gaseous onset).** Already implicit in the code's Henry
law — make it explicit and *use* it:

$$p_{sat}(T, c_d) = \frac{c_d}{H(T)},\qquad
H(T) = H_{ref}\exp\!\left[E_H\!\left(\tfrac1T - \tfrac1{T_{ref}}\right)\right]$$

(inversion of `equilibrium_dissolved_gas`, `fluid_properties.cpp:76–85`). For
the shipped R290/PZ68 config: p_sat = 0.2175/2.175e-7 = 1.0 MPa. Void onset for
gas-laden oil must track p_sat, not the vaporous `p_cav` (Grando 2006: gas
released "when saturation pressure is reached … two-phase flow thereafter";
experimental: Braun & Hendricks 1984, ASLE Trans. 27(1):1–14 — cavitation
extent grows with gas solubility; cavity contents are released gas, not vapor).

**Stage 1 — constrained-θ void coupling (incremental, keeps Elrod machinery).**
Implement behind `gas_pressure_coupling = NONE | VOID_COUPLED`
(default `NONE` = current behavior):

1. **Gas displaces liquid (void ceiling).** The full-film ceiling for θ is no
   longer 1 but

   $$\theta_{full}(i,j) = 1 - \alpha_g(i,j),$$

   used in the switch function `g(θ) = [θ ≥ θ_full]`, in the post-solve snap
   (`reynolds.cpp:333–337` snaps to θ_full, not 1.0), and in the supply/BC film
   content. A full-film cell that releases gas at p < p_sat acquires α_g > 0,
   its ceiling drops, it cavitates — **gaseous rupture at p_sat emerges
   mechanistically** instead of being imposed.
2. **Local void pressure.** The cavitated-zone plateau becomes cell-local:

   $$p_{void}(i,j) = \max\!\left(p_{cav},\; \frac{m_g R_g T}{(1-\theta)\,h}\right),$$

   i.e. the gas filling the void sets the plateau; at release equilibrium it
   approaches p_sat. Pressure recovery uses
   `p = p_void + β̄ ln(θ/θ_full)` in pressurized cells, `p = p_void` in
   cavitated cells. (Experimental support for non-fixed, near-saturation cavity
   pressures: Etsion & Ludwig 1982; Braun & Hendricks 1984 — cavity pressure
   0.137 → 0.034 MPa as speed rises.)
3. **Mixture compressibility.** In pressurized cells containing gas, the
   effective bulk modulus is the harmonic mixture (gas isothermal
   compressibility = 1/p):

   $$\frac{1}{\bar\beta} = \frac{\theta}{\beta} + \frac{\alpha_g}{p},$$

   so Γ in those cells uses β̄ — gas presence makes the film soft and caps
   pressure near p_void. This *is* the capping mechanism; do not add an ad-hoc
   clamp.
4. **Mass bookkeeping.** Liquid-solution continuity keeps θ as unknown with
   release sink, gas continuity is the existing `free_gas_mass` transport with
   the release source — both per unit area:

   $$\partial_t(\rho_l\theta h) + \nabla\!\cdot\!(\text{Couette/Poiseuille fluxes}) = -\Delta\dot m_g,
   \qquad
   \partial_t m_g + \nabla\!\cdot\!(m_g \mathbf{u}) = +\Delta\dot m_g,$$

   with the **existing** finite-rate kinetics (keep them — peer-reviewed
   precedent: Hao & Gu 2014, Tribology Int. 78:14–26, DOI
   10.1016/j.triboint.2014.04.028; experimental motivation: Etsion & Ludwig's
   p/p_s > 1 release/resorption loop):

   $$\Delta\dot m_g = k_{dg}\,(c_d - c_{eq}(p,T))\,\rho_l\,h\,\max(\theta,\theta_{min}).$$

5. **Forces/outputs** use ρ̄ and the WP-5 mixture viscosity.

**Stage 2 (strategic, pairs with WP-10).** Full homogeneous compressible
mixture EOS ρ̄(p, T, w) with equilibrium flash + kinetic relaxation, replacing
the hard switch entirely (Grando 2006; Bayada & Chupin 2013, ASME J. Tribology
135(4):041702). Do not start Stage 2 until Stage 1 is validated.

> **Design caution for the implementing agent.** Stage 1's assembly
> (θ_full ceiling + local p_void + β̄) is a synthesis consistent with the cited
> closures, not a transcription of a single published scheme. The 0-D and 1-D
> tests below are therefore *mandatory acceptance gates*, and the equilibrium
> limit must reproduce Grando's behavior (pressure plateau at p_sat; classical
> Reynolds-BC solution recovered when gas is disabled).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Add `enum class GasPressureCoupling { NONE, VOID_COUPLED }`, key `gas_pressure_coupling` (default NONE). Validation: `VOID_COUPLED` requires `GAS_CAVITATION_MIXTURE` and TRANSIENT (until WP-2's steady loop lands). Add `gas_constant_specific` if not present (R_g). |
| `fluid_properties.cpp/hpp` | Export `saturation_pressure(c_d, T, cfg)` (Henry inversion incl. TABLE path); set `rho = ρ_l θ + ρ_g α_g` when coupling on; keep `rho = ρ_l θ` when NONE. |
| `reynolds.cpp` (`solve_elrod`) | Per-cell `theta_full` field from α_g (lagged one outer iteration); switch, snap, Γ (with β̄), and Couette flux base density use ρ̄/θ_full; pressure recovery via local `p_void`. Release sink in the θ equation source. |
| `equation_of_state.hpp` | Generalize `switch_function(theta, theta_full)`, `pressure_from_theta(theta, theta_full, p_void, beta_bar)`; keep old signatures as θ_full=1, p_void=p_cav wrappers. |
| `gas_transport.cpp` | No structural change (transport already exists); ensure the release source is consistent with the θ-equation sink (one function, used by both). |
| `main.cpp` | Order inside WP-2's outer loop: properties → Reynolds(θ, p) → gas transport → exchange → properties → … (already close; the difference is iteration, WP-2). |
| `io.cpp` | Output `p_void`, `theta_full`, `p_sat` as optional fields. |

### Test cases (`tests/test_gas_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Decoupled regression**: `gas_pressure_coupling = NONE` on `config_r290_pz68_quick.txt` | Fields identical to current master to solver tolerance (≤ 1e-10 relative) |
| 2 | **Σα invariant**: random release/transport steps | `θ + α_g ≤ 1 + 1e-12` every cell, every step |
| 3 | **0-D closed-cell equilibrium**: single cell, no flux, p < p_sat initially, kinetics on | p → p_sat asymptotically; total propane mass (c_d ρ_l θ h + m_g) constant to 1e-12; equilibrium c_d = H p_sat |
| 4 | **0-D resorption**: start gassy at p > p_sat | m_g → 0, c_d → c_eq(p); no negative masses |
| 5 | **1-D pressure cap**: converging-diverging channel, supersaturated feed | In the gassy zone, p ∈ [p_void, p_sat·(1+5%)]; no cell sustains p ≪ p_sat with θ = θ_full while supersaturated |
| 6 | **Mass conservation (periodic θ, no-flux z)** | Global liquid and total-propane residual (WP-4 logger) ≤ 1e-10 per step |
| 7 | **Load-shift direction** (journal, ε = 0.9, R290 config) | Load with coupling ON > load with NONE (direction per IJR 2023 +21%; log the ratio, no hard % gate) |
| 8 | **Vaporous limit**: c_d = 0, no gas | Identical to plain Elrod (p_void ≡ p_cav, θ_full ≡ 1) |

---

## WP-2 — Outer coupling (Picard) iteration; real steady state; steady+gas guard

**Closes:** AUDIT §C4 (REVIEW N1, N7). **Severity:** MAJOR. **Depends on:** none
(but unblocks WP-1, WP-6 fidelity, WP-8).

### Physics

A segregated one-sweep "solution" is not a fixed point of the coupled system;
published transient cavitation codes iterate the coupling inside the step —
Ausas, Jai & Buscaglia 2009 (ASME J. Tribology 131(3):031702, DOI
10.1115/1.3142903) update even the journal motion "within the same relaxation
process." Fixed-point (Picard) iteration with under-relaxation:

$$\phi^{(k+1)} \leftarrow \omega_\phi\,\phi^{new} + (1-\omega_\phi)\,\phi^{(k)},
\qquad \phi \in \{\mu, \rho, T, c_d, m_g\},\ \omega_\phi \in (0,1].$$

Convergence on normalized residuals each outer iteration k:

$$R_p = \frac{\max_{ij}|p^{(k)}-p^{(k-1)}|}{p_{ref}},\quad
R_T = \frac{\max_{ij}|T^{(k)}-T^{(k-1)}|}{T_{ref}},\quad
R_c = \max_{ij}|c_d^{(k)}-c_d^{(k-1)}|,$$

stop when all < `coupling_tolerance` (default 1e-6) or `coupling_max_iters`
(default 30; warn on hitting the cap). **Report R_p in pressure space**, not
θ-space — β amplifies θ-noise by ~1e9 (AUDIT §C2/N4), so a θ-residual hides
pressure-level error.

`STEADY_STATE` = run this loop to convergence once (this is the operating
point). Steady + gas: currently `gas_dt = 0` silently freezes the gas model
(`main.cpp:243`, `fluid_properties.cpp:222`) — replace with either a hard
config error or `steady_gas_model = FROZEN | EQUILIBRIUM` where EQUILIBRIUM
sets c_d = c_eq(p,T) and α_g from the equilibrium flash each outer iteration.

### Code implementation

| File | Change |
|------|--------|
| `main.cpp` | Extract the per-step block (`main.cpp:220–260`) into `solve_coupled_step(fields, sys, mesh, cfg, bearing_state, dt)`; wrap in the outer loop with relaxation + residual log (rank-0, `Utils::log`). Steady mode: iterate to tolerance instead of `step == 0` single pass. |
| `config.hpp` | Keys: `coupling_max_iters`, `coupling_tolerance`, `coupling_relaxation` (one factor; per-field override later if needed), `steady_gas_model`. `validate()`: error on `GAS_CAVITATION_MIXTURE && STEADY_STATE && steady_gas_model == FROZEN` unless user sets an explicit `allow_frozen_gas_steady = true`. |
| `journal_motion.cpp` | (With WP-1/WP-8 in mind) optionally update bearing state inside the outer loop (Ausas-style) behind `motion_coupling = LAGGED | IN_LOOP`, default LAGGED. |
| `io.cpp`/logging | Per-step line: outer iters, final R_p/R_T/R_c, Elrod inner iters, flag flips (WP-4 shares this). |

### Test cases (`tests/test_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Constant-property limit**: `fluid_property_model = CONSTANT`, ISOTHERMAL | Converges in exactly 1 outer iteration; results identical to master |
| 2 | **Monotone convergence**: EMPIRICAL_CORRELATION μ(T) with strong dissipation | R_p, R_T decrease monotonically (allowing first-iter transient); converged < max iters |
| 3 | **Ordering independence**: permute sweep order (Reynolds↔energy) | Converged steady fields agree to `coupling_tolerance` (this is the definition of a fixed point — fails on master) |
| 4 | **Steady ≡ long transient**: steady solve vs transient run to t→∞, same config | max|p_steady − p_transient|/p_ref < 5·tolerance |
| 5 | **Guard**: GAS_CAVITATION_MIXTURE + STEADY_STATE + FROZEN | `validate()` throws/errors with actionable message |
| 6 | **Relaxation robustness**: ω = 0.5 and 1.0 on test 2 | Same converged answer |

---

## WP-3 — Numerics honesty: convection schemes, dead selectors, capacity-term gap

**Closes:** AUDIT §C1, §C5 (REVIEW N3, N6 + audit's new finding). **Severity:**
medium, cheap. **Depends on:** none.

### Physics

- First-order upwind *inside the cavitated zone* is standard — even FBNS uses
  it (Woloszynski et al. 2015, Tribology Letters, DOI 10.1007/s11249-015-0487-4:
  "in the cavitation region … a first-order upwind formula was employed").
  Upwinding *everywhere* is not: the canonical scheme is **type differencing** —
  central in the full film, upwind across/inside cavitation (Vijayaraghavan &
  Keith 1989, Tribology Trans. 32(2):225–233, DOI 10.1080/10402008908981882;
  San Andrés Notes 06). False diffusion of full upwind ≈ |u|Δx/2 (Patankar
  1980, ch. 5) lands exactly on the rupture/reformation/thermal/gas fronts.
- TVD face reconstruction (already coded, dead in `fvm.cpp:13–27,118–162`):
  $$\phi_f = \phi_U + \tfrac12\,\psi(r)\,(\phi_D-\phi_U),\quad
  r = \frac{\phi_U-\phi_{UU}}{\phi_D-\phi_U},\quad
  \psi_{vanLeer}(r)=\frac{r+|r|}{1+|r|}.$$
- Transient capacity: $\partial(\rho_0\theta h)/\partial t$ must be assembled
  for **all** transient runs; currently gated on `MOVING_BEARING`
  (`reynolds.cpp:261–262`), so a fixed-bearing transient (e.g. an
  `omega_ramp_time` startup) is quasi-steady in θ — film-history physics
  silently absent.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Key `theta_convection_scheme = UPWIND \| TVD_VANLEER \| TVD_MINMOD \| TYPE_DIFFERENCING` (default UPWIND for back-compat), similarly `thermal_convection_scheme`, `gas_convection_scheme`. **Delete** `pressure_time_method`/`temperature_time_method` members + GUI combos; parser accepts the old keys with a one-line warning ("ignored: field solves are backward-Euler"). |
| `reynolds.cpp` | Pass configured scheme to `FVM::divergence` (`:259`); TYPE_DIFFERENCING = per-face: central if `g=1` on both sides, upwind otherwise. Change the transient gate at `:261` to `cfg.solution_mode != STEADY_STATE` for the `ddt_weighted` capacity; squeeze source remains conditional on `dh_dt` presence. Align the Gumbel Couette discretization (`:146–147`) with the Elrod path's scheme so the A.1↔Gumbel limit compares like-for-like. |
| `energy.cpp`, `gas_transport.cpp` | Pass configured schemes (`energy.cpp:344`, `gas_transport.cpp:241,260`). |
| `gui_win32.cpp` | Remove dead time-method combos; add scheme dropdowns. |
| `fvm.cpp` | None expected — limiters/deferred correction exist; verify ghost-width (n_ghost = 2) suffices for the UU stencil at partition boundaries (it does; add the MPI test below anyway). |

### Test cases (extend `tests/test_fvm.cpp`, new `tests/test_schemes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **TVD boundedness**: 1-D step advection of θ | No new extrema (θ ∈ [min, max] of initial data); sharper front than upwind (L1 error strictly smaller) |
| 2 | **Order of accuracy**: smooth 1-D profile | Upwind ~1st order, TVD ≥ ~1.5 observed order on refinement |
| 3 | **A.1↔Gumbel match**: matched schemes both paths | Residual strictly below the current 1.2% (record value in test) |
| 4 | **Capacity-term fix**: fixed bearing, ω step at t=0 | θ relaxes over finite time (exponential-like), not instantaneously; steady final state matches steady solve |
| 5 | **Regression**: all schemes default UPWIND | Existing ctest suite untouched |
| 6 | **MPI**: TVD run, 2 vs 6 ranks | Rank-count independence of fields to 1e-12 (cf. `test_elrod_a2_mpi6`) |
| 7 | **Config back-compat**: old config with `pressure_time_method = CRANK_NICOLSON` | Parses, warns, runs |

Also deliverable: `validation/grid_convergence/` script running 60×20 / 120×40 /
240×80 on the journal case, reporting cavitation extent, h_min, load (archived
table — AUDIT §C1 requirement; feeds WP-7).

---

## WP-4 — Conservation, regime, and convergence diagnostics

**Closes:** AUDIT Part V item 7, §B1 logging (REVIEW N9, M3, N8-logging).
**Severity:** low effort / high value. **Depends on:** none. **Do first.**

### Physics

Per-step global balances on the physical domain (Ω = bearing surface,
dA = R dθ dz):

$$M_l = \int_\Omega \rho_l\,\theta\,h\,dA,\qquad
M_{gas,tot} = \int_\Omega \big(c_d\,\rho_l\,\theta\,h + m_g\big)\,dA,$$

$$r_l^{n+1} = \frac{M_l^{n+1}-M_l^{n} - \Delta t\,(\Phi_{in}-\Phi_{out}+S_{exch})}{\max(M_l^{n},\epsilon)},$$

likewise r_gas (exchange terms cancel in M_gas,tot — only boundary fluxes
remain). Boundary fluxes from the assembled face fluxes (Couette + Poiseuille at
z-boundaries and inlet penalties), not re-derived. Clamp/snap operations
(`reynolds.cpp:324–337`, gas clamps) must be **accounted as explicit source
terms** in the balance so their mass effect is visible (REVIEW N9: clamps are
not mass-neutral). Target per PLAN A.2: O(1e-12) in closed domains.

Regime guards at startup (config validate + rank-0 log):

- Couette Reynolds number and Taylor criterion (laminar validity):
  $$Re_c = \frac{\rho\,\omega R\,c}{\mu},\qquad
  Ta = Re_c\sqrt{c/R} \lesssim 41,$$
  warn above (Taylor-vortex onset; standard journal-bearing criterion, cf.
  Hamrock, *Fundamentals of Fluid Film Lubrication*; San Andrés Notes).
- Saturation context (AUDIT §B1): log `p_sat(T_init, c_d_init)`, ratio
  `bc_z_*_val / p_sat`, and `p_cav / p_sat`; warn when boundary oil enters
  saturated (ratio ≥ 1) so the current R290 configuration is visible.
- Per-step: Elrod outer iterations, flag flips, (after WP-2) coupling residuals.

### Code implementation

| File | Change |
|------|--------|
| `src/diagnostics.hpp/cpp` (new) | `Diagnostics::mass_balance(fields, mesh, cfg, dt)` returning {M_l, M_gas, r_l, r_gas, clamp_mass}; MPI_Allreduce inside; rank-0 CSV append (`results/<run>/diagnostics.csv`) + log line. |
| `reynolds.cpp` | Accumulate clamp/snap mass deltas into a fields-resident scalar (or returned struct) for the balance. |
| `config.hpp::validate()` | Re/Ta computation + warnings; p_sat ratio warnings; key `diagnostics_interval` (default 1 = every step). |
| `main.cpp` | Call after each step; write CSV. |

### Test cases (`tests/test_diagnostics.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Closed-domain conservation**: periodic θ, zero-gradient z, no inlets, transient | \|r_l\| ≤ 1e-12 per step (the PLAN goal, now actually checked) |
| 2 | **Open-domain accounting**: Dirichlet z ends | Balance closes: residual after flux accounting ≤ 1e-10 |
| 3 | **Clamp visibility**: initialize θ < θ_min somewhere | clamp_mass ≠ 0 reported, balance still closes |
| 4 | **Gas total**: WP-1 test 3 config | r_gas ≤ 1e-12 (exchange cancels) |
| 5 | **Re guard**: synthetic config with μ tiny | validate() warns (Ta > 41) |
| 6 | **Saturation log**: R290 config | startup log contains p_sat = 1.0e6 and ratio = 1.0 warning |

---

## WP-5 — Property models: mixture viscosity with gas-limit asymptote; solubility tables

**Closes:** AUDIT §A3, §A4 (REVIEW W3, W4). **Severity:** medium. **Depends
on:** none (WP-1 consumes the result).

### Physics

**Free-gas (bubbly/foam) viscosity.** Einstein μ = μ_l(1 + 2.5α_g) is a dilute
result (α ≲ 0.1) and is directionally wrong as α_g → 1 (gas must *thin* the
film). Replace with selectable models, all satisfying μ̄(0) = μ_l:

- `DUKLER_VOID` (linear void-weighted; the Grando 2006 density-side analogue):
  $$\bar\mu = \alpha_g\,\mu_g + (1-\alpha_g)\,\mu_l \;\to\; \mu_g \text{ as } \alpha_g\to1.$$
- `MCADAMS_QUALITY` (homogeneous two-phase standard; quality
  x_g = m_g/(m_g + ρ_lθh)):
  $$\frac{1}{\bar\mu} = \frac{x_g}{\mu_g} + \frac{1-x_g}{\mu_l}.$$
- `EINSTEIN_DILUTE` (keep, for back-compat) — but `validate()` errors when
  selected with `gas_alpha_max > 0.6` (Krieger–Dougherty packing bound):
  $$\bar\mu = \mu_l\,(1-\alpha_g/\alpha_{max})^{-2.5\,\alpha_{max}}\ \text{(KD, optional fourth model)}.$$

(Grando 2006 quote: "μ̄ = χμ_g + (1−χ)μ_l" — quality-weighted linear; offer it
as `LINEAR_QUALITY` if exact Grando replication is wanted for WP-7 V5.)

**Dissolved-gas (solution) viscosity.** The validated form for R290/oil pairs
is the Eyring-MTSM activity model (Int. J. Refrigeration 2023,
S0140700723002451: "viscosities of R-290/oil mixtures could be well represented
by the Eyring-MTSM model"; measured 303–348 K on PAG 68 / POE 75 / PVE 68).
Do **not** implement the activity model — feed measured isotherms through the
existing `TABLE` path (`viscosity_table`, `solubility_table`,
`density_table`), and document linear Henry + `exp(a_c c_d)` as a local tangent
valid near the (10 bar, 40 °C) calibration point only. Grando's reference case
(2 bar, w_sat = 7.13 %) shows how far the shipped 21.75 % case extrapolates.

> **Data task (flag to user, do not fabricate):** populating the R290/PZ68
> tables requires the actual isotherm data (e.g. the IJR 2023 dataset or
> supplier data for PZ68). Ship the *plumbing* + a clearly-labeled placeholder
> table; mark the case "table data pending" in config comments and README.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `enum class GasMixtureViscosityModel { EINSTEIN_DILUTE, DUKLER_VOID, MCADAMS_QUALITY, KRIEGER_DOUGHERTY, LINEAR_QUALITY }`, key `gas_mixture_viscosity_model` (default EINSTEIN_DILUTE = current behavior); `mu_gas` config key; validation pairing with `gas_alpha_max` as above. |
| `fluid_properties.cpp` | Replace hardcoded Einstein at `:202–207` with the dispatch; α_g clamp from `gas_alpha_max`, not the hardcoded 0.99. |
| `config_r290_pz68*.txt` | Set `gas_alpha_max = 0.6`, select `MCADAMS_QUALITY`; comment the validity statement. |

### Test cases (extend `tests/test_fluid_properties.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Asymptotes**: each model | μ̄(α=0) = μ_l exactly; DUKLER/MCADAMS/KD: μ̄ → μ_g (or finite KD limit) as α→α_max; monotone decreasing in α (except EINSTEIN which is documented-increasing) |
| 2 | **Einstein guard**: EINSTEIN + α_max = 1.0 | validate() error |
| 3 | **Regression**: default EINSTEIN, α_max as shipped | Existing test_fluid_properties values unchanged |
| 4 | **Table round-trip**: synthetic 3-point viscosity/solubility tables | Interpolation exact at nodes, linear between, clamped outside |
| 5 | **Calibration point**: R290 case (existing test) | Still passes with new defaults documented |

---

## WP-6 — θ-weighted energy capacity, enthalpy flux, and cavitated-zone shear

**Closes:** AUDIT §A5 inconsistency (REVIEW N5, part of W6). **Severity:**
medium. **Depends on:** none (interacts with WP-2 for converged THD).

### Physics

In a striated cavitated zone only the liquid fraction θ carries thermal mass,
enthalpy, and (predominantly) shear. Replace, behind
`cavitated_film_model = FULL_FILM | STRIATED` (default FULL_FILM = current):

$$\text{capacity: } \rho c_p\,\theta\,h\,\frac{\partial T}{\partial t},\qquad
\text{flux: } \rho c_p\,\theta_f\,h_f\,\bar u\,T,\qquad
\text{Couette shear: } \tau_\theta = \theta\,\frac{\mu\,\omega R}{h}\ (\text{cavitated cells}),$$

with full-film cells (θ = 1) unchanged. This restores internal consistency: the
gas transport already θ-weights its capacity (`gas_transport.cpp:31,193`).
Friction torque/power-loss are first-order THD deliverables (groove-design
practice reports power loss + h_min: ASME J. Tribology 144(12):121801, 2022) —
the FULL_FILM model overstates them by up to ~1/θ where θ → θ_min = 1e-6.
Cross-film profile and conjugate conduction remain Phase D scope; record that
the validated 1983–84 baseline includes them (Ferron 1983, ASME 105(3):422;
Lund & Tonnesen 1984, ASME 106(2):237–244 — 4th-order cross-film polynomial +
sleeve conduction) so temperature claims stay gated on WP-7 V6.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `cavitated_film_model` enum + key; validation: STRIATED requires ELROD_ADAMS. |
| `energy.cpp` | Capacity `:301` ×θ; flux loops `:308–340` ×θ_face (face value: upwind θ consistent with flux direction); wall-loss term: decide and document whether Q_w scales with θ (liquid contact) — recommended ×θ for consistency. |
| `reynolds.cpp` | `:519` and `:543`: Couette term ×(g_f + (1−g_f)·θ) i.e. full shear in full film, θ-scaled in cavitated cells. |

### Test cases (extend `tests/test_energy.cpp`, `tests/test_forces.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Full-film limit**: θ ≡ 1, STRIATED | Bit-identical to FULL_FILM |
| 2 | **Scaling**: uniform cavitated strip θ = 0.5 | Friction torque contribution of that strip = 0.5× FULL_FILM value; energy capacity likewise |
| 3 | **Energy balance**: WP-4 style closed thermal balance with θ-weighting | Heat generated − wall loss = dE/dt to 1e-10 |
| 4 | **Direction**: journal case with cavitation, STRIATED vs FULL_FILM | Power loss strictly lower with STRIATED; T_max not higher |
| 5 | **Regression**: default FULL_FILM | Suite unchanged |

---

## WP-7 — Validation campaign (archived, quantitative)

**Closes:** AUDIT Part V (REVIEW W11). **Depends on:** V1–V4 none; V5 needs
WP-1; V6 needs WP-2+WP-6.

### Structure

`validation/<case_id>/` containing `config*.txt`, `reference/` (digitized CSV +
provenance note: paper, figure/table number, digitization method), `run.ps1` +
`run.sh`, `compare.py` (tolerance-gated, exit code), `README.md`. Each case also
gets a coarse "smoke" variant wired into ctest (`add_test(NAME validation_<id>_smoke ...)`).

### Cases, references, acceptance

| ID | Case | Reference | Acceptance gate |
|----|------|-----------|-----------------|
| V1 | 1-D JFO slider (rupture + reformation), analytic | Elrod & Adams 1975; PLAN A.2 derivation | Rupture/reformation locations within 1 cell at 200 cells; p_peak ≤ 1 %; θ-profile L2 ≤ 1 % |
| V2 | Oscillatory squeeze film (exact solution) + dynamically loaded JB | Ausas, Jai & Buscaglia 2009, ASME JT 131(3):031702 (paper ships reference source code — use it for code-to-code) | p(t) L2 ≤ 2 % vs exact; JB orbit qualitative match + mass residual ≤ 1e-10 |
| V3 | Grooved journal bearing, flooded + (post-WP-9) starved | Vijayaraghavan & Keith 1989: Trib. Trans. 32(2):225–233 (algorithm cases, code-to-code) and Wear 134(2):377–397 (grooved, vs experiment) | Cavitation extent ±1 cell @ published grid; load ±5 % |
| V4 | LCP analytic suite (pocketed slider, squeeze damper) | Giacopini et al. 2010, ASME JT 132(4):041702 | p-profile L2 ≤ 2 % on the 1-D cases |
| V5 | Refrigerant equilibrium limit | Grando, Priest & Prata 2006 (Trib. Lett.) | With gas coupling ON and equilibrium kinetics (k_dg → large): recover classical Reynolds-BC solution for the moderately loaded case (their validation); plateau at p_sat |
| V6 | THD test bearing | Ferron, Frêne & Boncompain 1983, ASME 105(3):422; Lund & Tonnesen 1984, ASME 106(2):237 | Phase 1: qualitative T(θ) shape + documented gaps (no conjugate conduction). Hard gates deferred to Phase D; record deltas. |
| V7 | Grid convergence + conservation history | WP-3/WP-4 outputs | Monotone convergence of extent/load/h_min; archived plots |

**Citation hygiene for the implementing agent:** the V&K 1989 Trib. Trans.
paper is code-to-code only (its own abstract: comparisons "with the results
obtained using Elrod's algorithm"); the experimental comparison is the Wear
companion. Do not cite the former as experimental validation.

---

## WP-8 — Linearized dynamic coefficients (K_ij, C_ij) and whirl margin

**Closes:** AUDIT §D1 (REVIEW W9b). **Depends on:** WP-2 (needs a converged
operating point).

### Physics

Small-perturbation method (Lund 1987, ASME J. Tribology 109(1):37–41, DOI
10.1115/1.3261324 — the canonical reference; applied to gas-oil refrigerant
bearings via the perturbed compressible Reynolds equation in ASME JT
147(2):024101, DOI 10.1115/1.4066414):

$$K_{ij} = -\left.\frac{\partial F_i}{\partial x_j}\right|_{eq}
\approx -\frac{F_i(x_j^{eq}+\Delta x_j)-F_i(x_j^{eq}-\Delta x_j)}{2\,\Delta x_j},
\qquad
C_{ij} = -\left.\frac{\partial F_i}{\partial \dot x_j}\right|_{eq},$$

i, j ∈ {x, y}; force = fluid force on the journal/bearing per the existing
convention. Perturbation sizes: Δx = ε_K·c with ε_K ~ 1e-2 (linearity check:
halving Δx changes K by < 1 %); velocity perturbations Δẋ = ε_C·c·ω applied as
prescribed `dh/dt` with frozen position. Each force evaluation is a **converged
WP-2 operating point** (transient effects off; squeeze term from the prescribed
ẋ only).

Stability deliverable: for a rigid rotor of mass m per bearing, the
characteristic equation `det(m s² I + C s + K) = 0`; report eigenvalues, onset
speed (Re(s) = 0), and whirl frequency ratio Im(s)/ω (≈ 0.5 for plain bearings
— sanity anchor; San Andrés Notes 05). Lund's validity caveat (small motion
about equilibrium) goes in the output header.

### Code implementation

| File | Change |
|------|--------|
| `src/dynamic_coefficients.hpp/cpp` (new) | `compute_dynamic_coefficients(fields, sys, mesh, cfg, state)` → struct {K[2][2], C[2][2], whirl_ratio, onset_mass_or_speed}; loops ±perturbations calling the WP-2 steady driver. |
| `config.hpp` | `compute_dynamic_coefficients` (bool), `perturbation_position_fraction`, `perturbation_velocity_fraction`. |
| `main.cpp` | Post-run hook when enabled; CSV + log output. |
| `io.cpp` | Append K/C block to a summary file. |

### Test cases (`tests/test_dynamic_coefficients.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Linearity**: Δx and Δx/2 | K_ij relative change < 1 % |
| 2 | **Sign structure** (plain journal, moderate ε) | K_xx, K_yy > 0; K_xy·K_yx < 0 (cross-coupling signs per Lund 1987 / San Andrés Notes 05) |
| 3 | **Symmetry of C** (aligned, no gas) | C_xy ≈ C_yx within 5 % |
| 4 | **Consistency vs transient**: small free decay orbit vs linear model prediction | Frequency/decay within 10 % |
| 5 | **Whirl anchor**: light-load plain bearing | whirl ratio ∈ [0.45, 0.55] |

---

## WP-9 — Starved / mass-flux inlet option

**Closes:** AUDIT §B2 (REVIEW W8). **Depends on:** none.

### Physics

Starved feed = the supply cannot maintain a full film at the inlet; prescribe
inlet film fraction θ_in < 1 (or feed flow). JFO reformation downstream is
already handled by Elrod machinery. Precedent: starved cases inside the Elrod
framework since Vijayaraghavan & Keith 1989 (Wear 134(2):377–397 — "both
flooded and starved inlet conditions … various degrees of starvation"); groove
geometry as a design input per ASME JT 144(12):121801 (2022).

Two modes on `InletConfig`:

- `STARVED_THETA`: penalty-pin θ = θ_in (θ_in ∈ (θ_min, 1)); cavitated inlet
  cell (g = 0), Couette flux carries ρ_l θ_in h u into the domain.
- `MASS_FLUX`: prescribe feed ṁ [kg/s]; convert to equivalent θ_in over the
  inlet footprint from the local Couette flux:
  θ_in = ṁ / (ρ_l (ωR/2) h L_z,inlet), clamped to (θ_min, θ_supply(p)).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `InletConfig::Type` extends with `STARVED_THETA`, `MASS_FLUX` (parse `inlet_N_mode`, `inlet_N_theta`, `inlet_N_mdot`); validation: starved modes require ELROD_ADAMS. |
| `reynolds.cpp` (`apply_inlet_conditions`, `:60–91`) | Branch: pressure mode unchanged; starved modes pin θ_in via the same penalty; MASS_FLUX computes θ_in per cell from local h. |
| `energy.cpp`/`gas_transport.cpp` | Inlet thermal/composition contracts treat starved inlets like supply cells (existing reservoir-composition path). |
| `gui_win32.cpp` | Inlet mode selector + value field. |

### Test cases (`tests/test_inlets.cpp` extension)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Flooded equivalence**: θ_in = θ_supply(p_supply) | Matches pressure-inlet result to 1e-10 |
| 2 | **Monotonicity**: θ_in ∈ {1.0, 0.8, 0.5} | Load strictly decreasing; cavitation extent strictly increasing |
| 3 | **Flux fidelity** (MASS_FLUX) | Integrated inlet mass flux = prescribed ṁ within 1 % |
| 4 | **Conservation**: starved case under WP-4 logger | Balance closes ≤ 1e-10 |

---

## WP-10 — Strategic: complementarity (FBNS) cavitation kernel

**Closes:** AUDIT §C2, §C3 (β-stiffness, switch-iteration fragility; REVIEW
N4 + S3 caveat). **Depends on:** stabilize WP-1 first (formulation choice
interacts). **This retires `snap_tol` and the β-conditioning class of bugs.**

### Physics

JFO mass conservation as a complementarity problem (Giacopini et al. 2010, ASME
JT 132(4):041702, DOI 10.1115/1.4002215):

$$\bar p \ge 0,\qquad \vartheta \ge 0,\qquad \bar p\,\vartheta = 0,
\qquad \bar p = p - p_{cav},\ \ \vartheta = 1-\theta,$$

reformulated as a smooth unconstrained system via the Fischer–Burmeister
function (Woloszynski, Podsiadlo & Stachowiak 2015, Tribology Letters, DOI
10.1007/s11249-015-0487-4):

$$\phi_{FB}(\bar p_{ij}, \vartheta_{ij}) = \bar p_{ij} + \vartheta_{ij}
- \sqrt{\bar p_{ij}^2 + \vartheta_{ij}^2} = 0,$$

solved with Newton + Schur-complement reduction (FBNS) coupled to the
discretized Reynolds residual R(p, θ) = 0. Literature-reported performance: ~2
orders of magnitude faster than Elrod-Adams/LCP-pivoting baselines,
mesh-independent Newton counts (900 → 360 000 DOF), upwind θ-flux in the
cavitated region retained. Motivation for leaving the β-switch: solution
"strongly depend[s]" on β, β dominates the discrete system's stability
(eigenvalue analysis), softening tolerable only ~2 orders (3–40 % error beyond)
— ASME JT 139(3):031703 (2017); practitioner β-reduction practice per San
Andrés Notes 06. No β appears in the FBNS formulation; physical β returns only
through an optional compressible full-film EOS.

### Code implementation

| File | Change |
|------|--------|
| `src/reynolds_fbns.cpp/hpp` (new) | `Reynolds::solve_fbns(...)` parallel to `solve_elrod`; Newton loop with analytic Jacobian blocks (Reynolds w.r.t. p and θ; φ_FB diagonal blocks); reuse `LinearSystem`/PETSc (KSP inside Newton, or SNES if available). |
| `config.hpp` | `cavitation_model = GUMBEL \| ELROD_ADAMS \| FBNS`; `fbns_max_newton`, `fbns_tolerance`. |
| `main.cpp` | Dispatch. WP-1's mixture terms enter through ρ̄, p_void analogues — coordinate with whoever owns WP-1 Stage 2. |

### Test cases (`tests/test_fbns.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Equivalence**: V1 slider + journal case vs ELROD_ADAMS | p, θ agree within discretization tolerance (≤ 0.5 %) |
| 2 | **Complementarity**: converged fields | max(p̄_ij·ϑ_ij) ≤ 1e-10·p_ref |
| 3 | **Mesh-independent iterations** | Newton counts within ±2 across 60×20 → 240×80 |
| 4 | **No snap needed** | Flag-chatter test (the scenario behind `snap_tol`, reynolds.cpp:328–337 comment) converges without snapping |
| 5 | **Timing** | ≥ 5× faster than ELROD_ADAMS outer loop on the 240×80 case (log; literature says ~100×, don't gate on that) |

---

## WP-11 — Axial-boundary flux consistency (single source of truth)

**Closes:** AUDIT §B4 (line-level BC review, 2026-06-11). **Severity:** medium
(BC-3 live when cavitation reaches the ends; BC-1/2/4 latent — no shipped
config uses `INLET_OUTLET`). **Depends on:** none; coordinate with WP-1 (the
shared flux function is where `p_void` will later replace `p_cav`).

### Physics

The three solvers currently assemble the *same* axial boundary three different
ways. Required invariant: **one liquid mass flux per boundary face, used by
Elrod assembly, energy convection, gas convection, and the diagnostics
balance**, with direction-dependent temperature/composition layered on top.

Per face (south shown; north mirrored), with `p_bc = elrod_boundary_pressure(bc_val)`:

$$\dot m''_{bc} =
\begin{cases}
\dfrac{\rho_0\,\theta_P\,h^3}{12\mu}\,\dfrac{p_{bc}-p_P}{\Delta z/2} & \theta_P \ge 1 \ \text{(full film — EOS chain rule } \beta\,d\theta=\theta\,dp\text{)}\\[2ex]
\max\!\left(0,\ \dfrac{\rho_0\,h^3}{12\mu}\,\dfrac{p_{bc}-p_{cav}}{\Delta z/2}\right) & \theta_P < 1 \ \text{(cavitated — physical reformation rate)}
\end{cases}$$

The cavitated branch replaces the current ungated `Γ_base` Dirichlet link,
which overestimates re-flooding by ≈ `β(1−θ_P)/(p_bc−p_cav)` (~10²–10³×)
because the EOS `p = p_cav + β ln θ` is invalid where p ≡ p_cav. After WP-1,
`p_cav` in the cavitated branch becomes `p_void(i,j)`.

`INLET_OUTLET` semantics unified: for the **mass** equation it behaves as
DIRICHLET (bidirectional, like the Gumbel path already does); direction
dependence applies only to *what the flux carries* — inflow brings
`temperature_reference` (RESERVOIR mode) and `dissolved_gas_initial` / zero
free gas; outflow convects cell values (energy/gas already do exactly this).
The lagged `fields["pressure"] < bc_val` gate (stale by one timestep, frozen
across outer iterations) is deleted; direction comes from the sign of the
shared flux, evaluated with the current iterate.

Cavitated-cell reformation inflow must also carry reservoir enthalpy and
composition (today the g-gate zeroes both, so re-flooding liquid silently
adopts the cell's T and c_d).

### Code implementation

| File | Change |
|------|--------|
| `src/boundary_flux.hpp/cpp` (new) | `axial_boundary_mass_flux(fields, mesh, cfg, side, i)` implementing the two-branch formula; plus `axial_boundary_volume_flux` (free gas) and a struct return {mass_flux, is_inflow}. |
| `reynolds.cpp` (`solve_elrod`) | Full-film boundary cells: keep the Dirichlet link (equivalent). Cavitated boundary cells: drop the `Γ_base` link; add the reformation flux as an explicit θ-equation source `+ṁ''_{bc}·RΔθ/(ρ_0 h ...)` consistent with the assembled units. Delete the `INLET_OUTLET` lagged gate (`:388–414`); INLET_OUTLET assembles like DIRICHLET. |
| `energy.cpp` | `south/north_boundary_flux` call the shared function (×`c_p` face value) instead of the local g-gated expression; RESERVOIR/OPEN logic unchanged. |
| `gas_transport.cpp` | `south/north_boundary_mass_flux`/`_volume_flux` call the shared function; composition logic unchanged. |
| `diagnostics.cpp` | Balance uses the same shared function (removes any residual reimplementation). |
| `config.hpp` | Document NEUMANN = sealed/symmetry plane in the key comment + PHYSICS.md §6. |

### Test cases (`tests/test_boundary_fluxes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Triple consistency**: random full-film states | Elrod-assembled boundary flux, energy flux ÷ ρc_p-face, gas mass flux agree cell-by-cell to 1e-12 (by construction once shared) |
| 2 | **Cavitated reformation rate**: 1-D z-strip, cavitated interior, submerged end | Re-flood mass rate = `ρ₀h³(p_bc−p_cav)/(12μ·Δz/2)` analytic within 1 %; assert the old `Γ_base` link would give > 10× (regression demonstrating the fix) |
| 3 | **Outflow continuity**: interior p > p_bc with INLET_OUTLET | Liquid mass leaves; WP-4 balance closes ≤ 1e-10; energy/gas convect the same direction in the same step |
| 4 | **Reformation carries reservoir state**: cavitated end cell re-floods (RESERVOIR mode) | Cell c_d moves toward `dissolved_gas_initial`, T toward `temperature_reference` — not frozen |
| 5 | **Direction-flip robustness**: p_bc oscillating about interior pressure | Flux continuous through zero; no chatter; conservation holds |
| 6 | **Full-film regression**: shipped configs (ends full-film) | Results identical to current master to solver tolerance |
| 7 | **Gumbel parity**: INLET_OUTLET vs DIRICHLET in Gumbel | Identical (documents BC-4 resolution) |

## Phase E priority (consolidated from the review verdict)

- **Tier 1 (MAJOR).** WP-1 (two-phase closure — released gas is volumetrically inert to pressure) and WP-2 (outer Picard coupling — `STEADY_STATE` is a single segregated sweep, not a fixed point). **WP-12** supplies the `p_sat` data backend both depend on.
- **Tier 2 (medium).** WP-3 numerics (done), WP-5 property laws (done), WP-6 θ-weighted energy/shear, WP-11 boundary-flux consistency, WP-8 FSI / K-C extraction.
- **Tier 3 (low-effort honesty).** dead time selectors (done, WP-3), frozen-gas-in-steady (WP-2), config/doc drift (this consolidation), starvation option (WP-9), mass diagnostics (done, WP-4), regime guards (done, WP-4).
- **Tier 4.** WP-7 quantitative validation campaign.


## References

- Elrod, H.G. & Adams, M.L. (1975). A computer program for cavitation and starvation problems. *Cavitation and related phenomena in lubrication*, 37-41.
- Vijayaraghavan, D. & Keith, T.G. (1989). Development and evaluation of a cavitation algorithm. *Tribology Transactions*, 32(2), 225-233.
- Fesanghary, M. & Khonsari, M.M. (2011). A modification of the switch function in the Elrod cavitation algorithm. *ASME J. Tribology*, 133(2), 024501.
- Ausas, R.F., Jai, M. & Buscaglia, G.C. (2009). A mass-conserving algorithm for dynamical lubrication problems with cavitation. *ASME J. Tribology*, 131(3), 031702.
- Bayada, G. & Chupin, L. (2013). Compressible fluid model for hydrodynamic lubrication cavitation. *ASME J. Tribology*, 135(4), 041702.

### Audit / two-phase / refrigerant-bearing references (DOIs)

Key DOIs (from the consolidated audit): Elrod-Adams numerics — 10.1080/10402008908981882,
10.1115/1.3142903, 10.1115/1.4002215, 10.1007/s11249-015-0487-4; refrigerant
two-phase — 10.1007/s11249-006-9027-6, S0140700723004577, 10.1115/1.4066414,
10.1016/j.triboint.2014.04.028; experiments — Etsion & Ludwig (ASME JLT
104(2):157), Braun & Hendricks (10.1080/05698198408981539); THD — Ferron 1983
(ASME 105(3):422), Lund & Tonnesen (10.1115/1.3260891); dynamics — Lund 1987
(10.1115/1.3261324); review — Braun & Hannon (10.1243/13506501JET772).
