# Development Plan

## Completed Phases

### Phase 1 — Infrastructure
Mesh, Field, Communicator, IO, FVM operators, PETSc linear system. All working with MPI domain decomposition in $\theta$.

### Phase 2 — Physics Foundations
Film thickness (static eccentric), barotropic EOS ($p = p_{cav} + \beta\ln\theta$), VTK output.

### Phase 3 — Reynolds Solver (Gumbel)
Steady-state compressible Reynolds equation with lagged density. Gumbel cavitation (post-solve pressure clamp). BiCGStab + block Jacobi. All tests pass.

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

---

## References

- Elrod, H.G. & Adams, M.L. (1975). A computer program for cavitation and starvation problems. *Cavitation and related phenomena in lubrication*, 37-41.
- Vijayaraghavan, D. & Keith, T.G. (1989). Development and evaluation of a cavitation algorithm. *Tribology Transactions*, 32(2), 225-233.
- Fesanghary, M. & Khonsari, M.M. (2011). A modification of the switch function in the Elrod cavitation algorithm. *ASME J. Tribology*, 133(2), 024501.
- Ausas, R.F., Jai, M. & Buscaglia, G.C. (2009). A mass-conserving algorithm for dynamical lubrication problems with cavitation. *ASME J. Tribology*, 131(3), 031702.
- Bayada, G. & Chupin, L. (2013). Compressible fluid model for hydrodynamic lubrication cavitation. *ASME J. Tribology*, 135(4), 041702.
