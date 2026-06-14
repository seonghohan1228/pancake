# Physics and Numerical Methods

## 1. Problem Overview

`pancake` solves the **Reynolds equation** for thin-film lubrication in a **journal bearing** — a shaft of radius $R$ rotating at angular velocity $\omega$ inside a slightly larger bore, separated by a thin lubricant film. The geometry is described in cylindrical coordinates $(\theta, z)$ on the bearing surface, where $\theta \in [0, 2\pi)$ is the circumferential angle and $z \in [0, L]$ is the axial coordinate.

---

## 2. Film Geometry

The local film thickness for a statically eccentric shaft:

$$h(\theta) = c - e\cos(\theta - \psi)$$

where:
- $c$ = radial clearance (gap at concentric position)
- $e$ = eccentricity (offset of shaft centre from bore centre)
- $\psi$ = attitude angle (angular position of minimum gap)

The solver coordinate $\theta$ is measured from the positive x axis. The GUI
also exposes `load_angle_deg`, a visualization-only reference angle measured
from the positive y axis. It does not change the Reynolds solve; it rotates the
preview coordinate so bearing results can be inspected relative to the usual
load direction. Positive `load_angle_deg` follows increasing solver $\theta$,
so preview angle zero maps to $\theta = 90^\circ + \text{load_angle_deg}$.

The eccentricity ratio $\varepsilon = e/c \in [0,1)$ controls the bearing load: $\varepsilon = 0$ gives a concentric (uniform) film; $\varepsilon \to 1$ gives metal-to-metal contact.

---

## 3. Reynolds Equation

Starting from the Navier-Stokes equations integrated through the film thickness $h$ (lubrication approximation: $h \ll R$, $Re \ll 1$), the **compressible thin-film Reynolds equation** in arc-length coordinates $(s, z)$, with $s = R\theta$, is:

$$\frac{\partial}{\partial s}\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial s}\right] + \frac{\partial}{\partial z}\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial z}\right] = \frac{U}{2}\frac{\partial(\rho h)}{\partial s} + \frac{\partial(\rho h)}{\partial t}$$

where $U = \omega R$ is the shaft surface speed. Substituting $s = R\theta$:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\left[\frac{\rho h^3}{12\mu}\frac{\partial p}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta} + \frac{\partial(\rho h)}{\partial t}$$

Defining the **diffusion coefficient** $\gamma = \rho h^3 / (12\mu)$, this is:

$$\underbrace{\frac{1}{R^2}\frac{\partial}{\partial\theta}\left[\gamma\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\left[\gamma\frac{\partial p}{\partial z}\right]}_{\text{pressure-driven flow (Poiseuille)}} = \underbrace{\frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta}}_{\text{Couette source}} + \underbrace{\frac{\partial(\rho h)}{\partial t}}_{\text{transient / compressibility}}$$

The LHS is recognised as the surface divergence $\text{div}(\gamma\,\text{grad}\,p)$ on the cylindrical surface.

**Steady-state incompressible** (standard form): set $\rho = \text{const}$ and $\partial/\partial t = 0$, giving the classical half-Sommerfeld problem.

---

## 4. Barotropic Equation of State and Cavitation

A **barotropic EOS** links density to pressure and handles sub-ambient (cavitated) regions via a **universal variable** $\theta_{film}$ (fractional film content):

| Region | Condition | Pressure | Density |
|--------|-----------|----------|---------|
| Full film | $\theta_{film} \ge 1$ | $p = p_{cav} + \beta\ln\theta_{film}$ | $\rho = \rho_0\,\theta_{film}$ |
| Cavitated | $\theta_{film} < 1$ | $p = p_{cav}$ | $\rho = \rho_0\,\theta_{film}$ |

Inverse relations:

$$\theta_{film}(p) = \begin{cases} \exp\!\left(\dfrac{p - p_{cav}}{\beta}\right) & p > p_{cav} \\ 1 & p \le p_{cav} \end{cases}$$

where $\beta$ is the bulk modulus. The current solver uses the **Gumbel (half-Sommerfeld) cavitation condition**: after solving for pressure, all cells with $p < p_{cav}$ are clamped to $p_{cav}$. Density is then updated from the EOS.

The VTK/GUI field named `film_content` writes this Elrod/JFO universal variable
directly. It is not clamped for visualization: values below one indicate
cavitated cells, while values above one indicate compressed full-film cells and
are the same values used to recover pressure through the EOS.

For pressure-valued supply features, the imposed Elrod variable is
$\theta_{supply} = \exp((p_{supply}-p_{cav})/\beta)$. A pressure of
$1.5\,\text{MPa}$ with $\beta=10^5\,\text{Pa}$ implies
$\theta_{supply}\approx 3.27\times 10^6$, which is numerically valid for this
EOS but not physically representative of oil. A realistic oil bulk modulus is
typically closer to $10^9\,\text{Pa}$, where the same pressure gives
$\theta_{supply}\approx 1.0015$.

The pressure derivative $dp/d\theta_{film}$ (used for effective diffusion in a full $\theta_{film}$-formulation):

$$\frac{dp}{d\theta_{film}} = \begin{cases} \beta / \theta_{film} & \theta_{film} \ge 1 \\ 0 & \theta_{film} < 1 \end{cases}$$

---

## 5. Finite Volume Discretisation

### 5.1 Grid

The domain is a structured $N_\theta \times N_z$ grid of cell-centred control volumes. Each cell $(i, j)$ has:
- centre at $\theta_i = (i + \tfrac{1}{2})\Delta\theta$, $z_j = (j + \tfrac{1}{2})\Delta z$
- area $A_{cell} = R\,\Delta\theta\,\Delta z$

### 5.2 FVM Integral Form

Integrating $\text{div}(\gamma\,\text{grad}\,p) = f$ over cell area and applying the divergence theorem on the cylindrical surface:

$$\oint_{\partial\Omega} \gamma\,(\nabla p \cdot \hat{n})\,dl = \int_\Omega f\,dA$$

Face contributions (outward flux positive):

| Face | Coefficient | Expression |
|------|-------------|------------|
| East ($+\theta$) | $a_E$ | $\gamma_e\,\Delta z\,/\,(R\,\Delta\theta)$ |
| West ($-\theta$) | $a_W$ | $\gamma_w\,\Delta z\,/\,(R\,\Delta\theta)$ |
| North ($+z$) | $a_N$ | $\gamma_n\,R\,\Delta\theta\,/\,\Delta z$ |
| South ($-z$) | $a_S$ | $\gamma_s\,R\,\Delta\theta\,/\,\Delta z$ |

Face values of $\gamma$ use **harmonic averaging**: $\gamma_e = 2\gamma_i\gamma_{i+1}/(\gamma_i + \gamma_{i+1})$.

### 5.3 Assembled Stencil

$$a_P\,p_P - a_E\,p_E - a_W\,p_W - a_N\,p_N - a_S\,p_S = \text{source}$$

with $a_P = a_E + a_W + a_N + a_S$ (plus any boundary contributions).

The LHS sums face contributions $a_f(p_P - p_{\text{nbr}})$, which is the *negative* of the outward flux sum from the divergence theorem. The stencil therefore represents:

$$-\,\text{div}(\gamma\,\nabla p)\cdot A_\text{cell} = \text{source}$$

**Sign convention — `source_field` passed to `FVM::add_source` must be $-f$, not $+f$.**

For a PDE of the form $\text{div}(\gamma\,\nabla p) = f$, rearranging gives:

$$-\,\text{div}(\gamma\,\nabla p)\cdot A_\text{cell} = -f\cdot A_\text{cell}$$

so `FVM::add_source` must receive `source_field` $= -f$. Applied to the Reynolds equation with Couette term $S = \tfrac{\omega}{2}\partial_\theta(\rho h)$:

$$\text{div}(\gamma\,\nabla p) = S \quad\Longrightarrow\quad \text{source\_field} = -S$$

### 5.4 Time Derivative

The `FVM::ddt` operator discretises $\partial\phi/\partial t$ with a first-order backward Euler scheme, adding:

$$a_P \leftarrow a_P + \frac{A_\text{cell}}{\Delta t}, \qquad \text{source} \leftarrow \text{source} + \frac{A_\text{cell}}{\Delta t}\,\phi^n$$

### 5.5 Convection Schemes

`FVM::divergence` assembles $\text{div}(F\,\phi)$ for a known face flux field $F$. The upwind (first-order) scheme adds implicit contributions; higher-order schemes add an explicit deferred correction using $\phi^{\text{old}}$ (the previous timestep or, in the Elrod outer loop, the lagged iterate).

The face interpolation is selectable per transported variable
(`theta_convection_scheme`, `thermal_convection_scheme`,
`gas_convection_scheme`; default `UPWIND` for backward compatibility):

- `UPWIND` — first-order, unconditionally bounded, false diffusion $\approx |u|\Delta x/2$ (Patankar 1980) concentrated exactly on rupture/reformation/thermal/gas fronts.
- `TVD_VANLEER`, `TVD_MINMOD` — limited high-order deferred correction
  $\phi_f = \phi_U + \tfrac12\psi(r)(\phi_D-\phi_U)$ with
  $r = (\phi_U-\phi_{UU})/(\phi_D-\phi_U)$; sharper fronts, no new extrema.
- `TYPE_DIFFERENCING` (theta only) — Vijayaraghavan & Keith (1989): central
  interpolation across faces where both adjacent cells are full-film
  ($g = 1$), first-order upwind across and inside the cavitated region. This
  is the canonical Elrod-Adams discretization; full upwind everywhere is a
  conservative-but-diffusive default.

Validation (`tests/test_schemes.cpp`): step advection is bounded and sharper
under TVD than upwind; observed refinement orders are $\approx 1.09$ (upwind)
and $\approx 1.99$ (TVD van Leer) on a smooth profile with $\Delta t \propto \Delta x^2$.

---

## 6. Boundary Conditions

### Circumferential ($\theta$)
**Periodic**: ghost cells are exchanged via MPI `Sendrecv` between adjacent ranks (rank 0 ↔ rank $N_{ranks}-1$ for the wrap-around).

### Axial ($z$)
**Dirichlet** $p = p_{bc}$ at $z = 0$ and $z = L$, where $p_{bc}$ is the configured
absolute boundary pressure `bc_z_*_val`. Implemented by adding the boundary face
flux to the stencil using a half-cell distance:

$$a_P \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}, \qquad \text{source} \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}\,p_{bc}$$

A truly open (vented) end corresponds to $p_{bc} = p_{cav}$; the shipped
configurations are **submerged** ends with $p_{bc} > p_{cav}$ (e.g. a flooded
sump or a pressurized refrigerant housing), which the Elrod path converts to a
boundary film content $\theta_{bc} = \exp[(p_{bc}-p_{cav})/\beta] > 1$.
`NEUMANN` selects zero-gradient ends; `INLET_OUTLET` applies the Dirichlet
value only on inflow.

**Cavity-pressure caveat.** The model holds the cavitated-zone pressure at the
single constant `p_cav`. Measurements show the cavity pressure of gas-laden
oils is neither constant nor equal to the vapor pressure — it tracks the gas
release state and drops with speed (Etsion & Ludwig 1982; Braun & Hendricks
1984). A constant `p_cav` is therefore a model input to be chosen per case, and
the WP-1 void coupling (`PLAN.md` Phase E) is the planned mechanistic
replacement; see §16 for the oil–refrigerant gaseous-cavitation theory and
`docs/CAVITATION_OIL_REFRIGERANT.md` for the cited basis.

---

## 7. Reynolds Solvers

Two solvers are available, selected by `cfg.cavitation_model`.

Before the time loop, the pressure field is initialized from the first
non-Neumann axial outlet pressure, not from a supply groove. For
`ELROD_ADAMS`, the initialized and boundary film content are derived from that
pressure and the bulk modulus:

$$p^0 = \max(p_{bc},p_{cav}),\qquad
\theta^0 =
\max\!\left[\exp\!\left(\frac{p^0-p_{cav}}{\beta}\right),\theta_{min}\right]$$

Legacy `bc_z_*_theta` config keys are parsed for old files but are not used by
the solver. This keeps boundary and initial conditions pressure-based:
changing pressure or bulk modulus changes the corresponding Elrod film content
through the EOS. Supply inlets are still imposed during the linear solve.

### 7.1 Gumbel Solver (`Reynolds::solve`)

The solver follows a **lagged-coefficient** approach at each time step:

1. Compute $\gamma^n = \rho^n h^3 / (12\mu)$ at physical + theta ghost cells (using already-synced ghost values of $\rho$ and $h$).
2. Reset and assemble the linear system:
   - `FVM::laplacian` with $\gamma^n$
   - Dirichlet z-boundary corrections
   - Couette source: $-\frac{\omega}{2}\frac{\partial(\rho h)^n}{\partial\theta}$,
     evaluated as a face-flux divergence matched to the Elrod path: the product
     $\rho h$ is factorized into a central face value of $\rho_0 h$ times the
     transported content $\theta = \rho/\rho_0$ interpolated with the configured
     `theta_convection_scheme`. This keeps the Gumbel and Elrod wedge
     discretizations like-for-like (validated to $\approx 0.01\%$ peak-pressure
     agreement under `TYPE_DIFFERENCING` in `tests/test_schemes.cpp`).
3. Solve $A\,p^{n+1} = b$ via PETSc KSP (BiCGStab + block Jacobi).
4. Apply Gumbel cavitation: $p^{n+1} \leftarrow \max(p^{n+1},\, p_{cav})$.
5. Update $\rho^{n+1} = \rho_0\,\theta_{film}(p^{n+1})$ via EOS.

For transient runs an explicit squeeze source $-\rho\,\partial h/\partial t$ is
added whenever a `dh_dt` field is present; the term is omitted only for
`STEADY_STATE` operating-point solves. Gumbel remains a pressure clamp and is
not mass-conserving (a validation warning says so).

### 7.2 Elrod-Adams Solver — Phase A.1 (`Reynolds::solve_elrod`, full-film)

Uses the **universal variable** $\theta$ (fractional film content) as the primary unknown. This is the foundation for mass-conserving cavitation.

**Derivation of the $\theta$-formulation.** Substituting $\rho = \rho_0\theta$ and $\partial p/\partial x = (\beta/\theta)\,\partial\theta/\partial x$ (for full film, $\theta \ge 1$) into the Reynolds equation:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\frac{\rho_0\theta h^3}{12\mu}\cdot\frac{\beta}{\theta}\frac{\partial\theta}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\frac{\rho_0\theta h^3}{12\mu}\cdot\frac{\beta}{\theta}\frac{\partial\theta}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho_0\theta h)}{\partial\theta}$$

The $\theta$ cancels in the diffusion coefficient:

$$\underbrace{\text{div}(\Gamma\,\nabla\theta)}_{\text{Poiseuille}} = \underbrace{\frac{\omega}{2}\,\rho_0\frac{\partial(\theta h)}{\partial\theta}}_{\text{Couette}}$$

where the **effective diffusion coefficient** is:

$$\Gamma = \frac{\rho_0\,\beta\,h^3}{12\mu}$$

Note $\Gamma = \beta\,\gamma$; for full-film, $\Gamma$ is independent of $\theta$, making the equation **linear** in $\theta$.

**Couette term as convection.** The RHS is the 1-D divergence of a Couette flux carrying $\theta$:

$$\rho_0\frac{\omega}{2}\frac{\partial(\theta h)}{\partial\theta} = \frac{\partial}{\partial\theta}\!\left[\rho_0\frac{\omega R h}{2}\,\theta\right] = \text{div}_\theta(\mathbf{F}\,\theta)$$

where the face flux at the east face of cell $(i,j)$ is:

$$F_e = \rho_0\,\frac{\omega}{2}\,R\,h_e\,\Delta z$$

with $h_e = \tfrac{1}{2}(h_i + h_{i+1})$. This flux is passed to `FVM::divergence` (upwind scheme). The assembled equation (FVM sign convention):

$$\underbrace{-\text{div}(\Gamma\,\nabla\theta)\cdot V}_{\text{laplacian}} + \underbrace{\text{div}(\mathbf{F}\,\theta)\cdot V}_{\text{divergence}} = 0$$

**Boundary conditions.** Since $\theta(p_{cav}) = 1$ by the EOS, the z-boundary condition $p = p_{cav}$ translates to $\theta = 1$.

**Solve sequence per timestep (Phase A.1):**

1. Compute $\Gamma = \rho_0\beta h^3/(12\mu)$ at physical + theta ghost cells.
2. Compute Couette face fluxes $F_e = \rho_0\,(\omega/2)\,R\,h_e\,\Delta z$ for $i \in [-1, N_\theta)$ (the $i=-1$ ghost is needed for the west face of the first cell in `FVM::divergence`).
3. Reset and assemble:
   - `FVM::laplacian(sys, Gamma, mesh)`
   - `FVM::divergence(sys, couette_flux, zero_flux_z, theta, UPWIND, mesh)`
   - Dirichlet z-BC: $\theta = 1$ at $z = 0,\,L$
4. Solve $A\,\theta^{n+1} = b$.
5. Clamp $\theta \leftarrow \max(\theta, 1)$ (full-film constraint; relaxed in Phase A.2).
6. Recover $p = p_{cav} + \beta\ln(\theta)$, update $\rho = \rho_0\,\theta$.

**Relationship to Gumbel solver.** Both formulations are equivalent in the full-film ($\theta \ge 1$) limit at small compressibility. Historically the two paths differed by $\approx 1.2\%$ because the Gumbel wedge term was centered while the Elrod path was upwind; since the WP-3 alignment both paths use the configured `theta_convection_scheme` and the matched full-film comparison agrees to $\approx 0.01\%$ (`tests/test_schemes.cpp`). The switch function $g(\theta)$ (Phase A.2) handles the cavitated region.

**Weighted time derivative.** The FVM operator `FVM::ddt_weighted(sys, theta, weight, dt, mesh)` discretises the liquid-content storage $\partial(\rho_0\theta h)/\partial t$ with weight $\rho_0 h$. It adds $\rho_0 h\,V/\Delta t$ to $a_P$ and $\rho_0 h\,V\,\theta^n/\Delta t$ to `source`. This capacity term is assembled for **every** transient solve (`solution_mode = TRANSIENT`), fixed or moving bearing — gating it on bearing motion would make a fixed-bearing transient (e.g. an `omega_ramp_time` startup) quasi-steady in $\theta$ with film-history physics silently absent. The explicit squeeze source $-\rho_0\theta^n\,\partial h/\partial t$ is added whenever a `dh_dt` field is present. `STEADY_STATE` mode omits both, giving the quasi-steady operating-point solve.

---

## 8. Linear System

`LinearSystem` wraps a PETSc `MPIAIJ` matrix. Each physical cell $(i,j)$ maps to global row $(\text{offset}_\theta + i)\cdot N_z + j$. Periodic $\theta$-wrapping is applied in the column indices. The system is solved with:
- **KSP**: BiCGStab (`KSPBCGS`)
- **PC**: Block Jacobi (`PCBJACOBI`)
- **Tolerance**: relative $10^{-8}$, max 1000 iterations

---

## 9. MPI Parallelism

The domain is decomposed in the $\theta$ direction: each MPI rank owns $N_\theta / N_{ranks}$ contiguous $\theta$-columns with the full $z$ extent. Ghost layers (default 2 cells) are exchanged each step via `Communicator::update_ghosts`, which packs/unpacks the non-contiguous ghost strips using `MPI_Sendrecv`.

For the Elrod-Adams solver, this exchange is also required inside the nonlinear
active-set loop. After each outer solve updates physical `theta`, the code
clamps and snaps near-full-film cells, then exchanges `theta` ghost cells before
the next iteration recomputes $g(\theta)$ on cross-rank faces. Without that
inner exchange, MPI partition interfaces can use stale cavitation flags and
produce rank-aligned pressure or `film_content` contours.

---

## 10. Notation Summary

| Symbol | Description | SI unit |
|--------|-------------|---------|
| $R$ | Shaft radius | m |
| $c$ | Radial clearance | m |
| $e$ | Eccentricity | m |
| $L$ | Bearing length | m |
| $\omega$ | Shaft angular velocity | rad/s |
| $\mu$ | Dynamic viscosity | Pa·s |
| $\rho_0$ | Reference density | kg/m³ |
| $\beta$ | Bulk modulus | Pa |
| $p_{cav}$ | Cavitation pressure | Pa |
| $h$ | Local film thickness | m |
| $p$ | Hydrodynamic pressure | Pa |
| $\gamma$ | Diffusion coefficient $\rho h^3/(12\mu)$ | m·s |
| $\theta_{film}$ | Fractional film content (EOS variable) | — |

---

## 11. Macroscopic Properties

After solving the Reynolds equation, the code integrates the local fields to compute the global bearing performance characteristics.

### Pressure and Viscous Forces

The force convention is now the force applied by the fluid to the moving outer
bearing surface. The shaft remains fixed. Gauge pressure $(p-p_{cav})$ is used
because uniform ambient pressure exerts zero net force on a closed body.

For bearing-surface radius $r_b = R+h$, the pressure area vector is:

$$d\mathbf{A}_b =
\left(r_b\mathbf{e}_r - \frac{\partial h}{\partial\theta}\mathbf{e}_\theta
      - r_b\frac{\partial h}{\partial z}\mathbf{e}_z\right)d\theta\,dz$$

so the pressure-only bearing force is:

$$\mathbf{F}^{p}_{bearing}
  = \int_0^L \int_0^{2\pi} (p-p_{cav})\,d\mathbf{A}_b$$

The solver writes this as `pressure_force_x`, `pressure_force_y`, and
`pressure_force_z`. Legacy `load_x`, `load_y`, and `load_z` output requests are
mapped to these pressure-force fields for older configs, but duplicate `load_*`
runtime fields are no longer allocated.

The shear force on the bearing surface uses the bearing-side wall shear:

$$\tau_{\theta,b} = \frac{\mu \omega R}{h}
  - g(\theta_{film}) \frac{h}{2R}\frac{\partial p}{\partial \theta}$$

$$\tau_{z,b} = -g(\theta_{film})\frac{h}{2}\frac{\partial p}{\partial z}$$

where $g(\theta_{film}) = 1$ in the full film and $g = 0$ in the cavitated
region. The viscous bearing force is:

$$F_x^\tau = \int_0^L \int_0^{2\pi} -\tau_{\theta,b}\sin\theta \cdot r_b\,d\theta\,dz$$

$$F_y^\tau = \int_0^L \int_0^{2\pi} \tau_{\theta,b}\cos\theta \cdot r_b\,d\theta\,dz$$

$$F_z^\tau = \int_0^L \int_0^{2\pi} \tau_{z,b} \cdot r_b\,d\theta\,dz$$

These are written as `viscous_force_x`, `viscous_force_y`, and
`viscous_force_z`.

### Friction Torque

The friction torque $T$ on the shaft arises from the circumferential shear:

$$T = \int_0^L \int_0^{2\pi} \tau_\theta \cdot R^2 \, d\theta \, dz$$

### Moving-bearing Fluid Force

For the moving-bearing model, the shaft is fixed and the outer bearing surface
moves. The fluid force used in the bearing equation of motion is the total
bearing-surface force:

$$\mathbf{F}^{fluid}_{bearing}
  = \mathbf{F}^{p}_{bearing} + \mathbf{F}^{\tau}_{bearing}$$

The solver writes this moving-bearing resultant as `fluid_force_x`,
`fluid_force_y`, and `fluid_force_z`. The independent configured applied load is
written separately as `external_load_x`, `external_load_y`, and
`external_load_z`.

---

## 12. Moving Outer-Bearing Motion

The implemented motion model fixes the shaft in space and moves the outer
bearing surface. Let $(x_b(t), y_b(t), z_b(t))$ be the bearing-center
displacement in the global frame. The axial displacement $z_b$ is tracked for
force and motion output, but it does not change the current 2D film geometry.

The moving-bearing film thickness is:

$$h(\theta,z,t) = c + x_b(t)\cos\theta + y_b(t)\sin\theta
                 - \alpha_x(z-L/2)\cos\theta
                 - \alpha_y(z-L/2)\sin\theta$$

This sign convention is the opposite of the static shaft-eccentricity
components. If `bearing_initial_from_attitude = true`, the initial bearing
position is:

$$x_b(0) = -e\cos\psi,\qquad y_b(0) = -e\sin\psi$$

which exactly reproduces the static gap $h = c - e\cos(\theta-\psi)$. During
motion, the equivalent attitude angle of the minimum gap is:

$$\psi_b = \operatorname{atan2}(-y_b, -x_b)$$

The moving-gap time derivative is:

$$\frac{\partial h}{\partial t}
  = \dot{x}_b\cos\theta + \dot{y}_b\sin\theta$$

For the Elrod variable $\theta_{film}$:

$$\frac{\partial(\rho_0\theta_{film}h)}{\partial t}
  = \rho_0 h\frac{\partial\theta_{film}}{\partial t}
  + \rho_0\theta_{film}\frac{\partial h}{\partial t}$$

The existing weighted time derivative covers the first term for static
geometry. Dynamic geometry also needs the explicit or Picard-lagged
$\theta_{film}\,\partial h/\partial t$ contribution.

The rigid bearing support model is:

$$m_b\ddot{x}_b + c_x\dot{x}_b + k_x x_b =
  F_x^{fluid} + F_x^{external}$$

$$m_b\ddot{y}_b + c_y\dot{y}_b + k_y y_b =
  F_y^{fluid} + F_y^{external}$$

$$m_b\ddot{z}_b + c_z\dot{z}_b + k_z z_b =
  F_z^{fluid} + F_z^{external}$$

The external load is configured directly and remains independent of bearing
motion. The visible implemented choices are `EULER_EXPLICIT`,
`EULER_IMPLICIT`, and `CRANK_NICOLSON`. `CRANK_NICOLSON` is the formal name for
the mixed implicit/explicit trapezoidal update; the parser still accepts older
`RK2`/`RK4` aliases for compatibility, but the GUI and default config no longer
advertise them as pressure or thermal methods.

---

## 13. Thermal and ETHD Model

The current solver implements the first thermal step: an optional
film-averaged lubricant energy equation with constant material properties,
viscous heat generation, and wall heat transfer. Temperature-dependent
viscosity, elastic deformation, and fully coupled ETHD remain planned follow-on
models.

### 13.1 Lubricant Energy Equation

A depth-averaged lubricant temperature model is:

$$\rho c_p h\left(\frac{\partial T}{\partial t}
  + \bar{u}_\theta\frac{1}{R}\frac{\partial T}{\partial\theta}
  + \bar{u}_z\frac{\partial T}{\partial z}\right)
  =
  \frac{1}{R^2}\frac{\partial}{\partial\theta}
  \left(k h\frac{\partial T}{\partial\theta}\right)
  + \frac{\partial}{\partial z}
  \left(k h\frac{\partial T}{\partial z}\right)
  + \Phi - Q_w$$

where:
- $T$ is the film-averaged lubricant temperature.
- $c_p$ is lubricant specific heat.
- $k$ is lubricant thermal conductivity.
- $\bar{u}_\theta$ and $\bar{u}_z$ are film-averaged velocities.
- $\Phi$ is viscous heat generation.
- $Q_w$ is net heat transfer to journal and bearing walls.

A practical initial dissipation model is:

$$\Phi \approx \frac{\mu U^2}{h}
       + \frac{h^3}{12\mu}\left[
         \left(\frac{1}{R}\frac{\partial p}{\partial\theta}\right)^2
         + \left(\frac{\partial p}{\partial z}\right)^2\right]$$

The first implementation stores $\Phi$ in `heat_generation` as heat generation
per bearing surface area [W/m^2]. The first term is Couette shear heating. The
second term is pressure-flow heating from Poiseuille shear and is active only in
full-film cells. Wall heat loss is represented by:

$$Q_w = h_j(T - T_j) + h_b(T - T_b)$$

where $h_j$ and $h_b$ are heat-transfer coefficients to the journal and bearing
surfaces.

For Elrod-Adams cavitation, pressure-flow heat uses the same active faces as
the pressure-driven velocity. A face contributes only when both adjacent cells
are full film:

$$g_f =
\begin{cases}
1, & \theta_{film,L}\ge 1 \text{ and } \theta_{film,R}\ge 1\\
0, & \text{otherwise}
\end{cases}$$

The cell pressure-gradient magnitude is evaluated as an average of squared
active face gradients, for example in the circumferential direction:

$$\left|\nabla_s p\right|^2_{\theta,P}
\approx
\frac{1}{2}g_e\left(\frac{p_E-p_P}{R\Delta\theta}\right)^2
+\frac{1}{2}g_w\left(\frac{p_P-p_W}{R\Delta\theta}\right)^2$$

with the same form in $z$. Cavitated faces and exterior axial boundary faces
therefore are not differenced against the boundary pressure. At the physical
axial ends, the boundary cell uses the adjacent interior active face as a
one-sided bounded gradient, so a smooth interior pressure field remains smooth
right up to the output boundary. This is a post-processing regularisation of the
hard JFO active set: cavitated cells still have $p=p_{cav}$, and the discrete
pressure plateau at the rupture/reformation boundary is expected for the current
first-order Elrod-Adams model.

The finite-volume equation solved for `temperature` is:

$$\frac{\rho c_p h V_s}{\Delta t}T_P^{n+1}
  + \mathcal{L}(T^{n+1}) + (h_j+h_b)V_s T_P^{n+1}
  =
  \frac{\rho c_p h V_s}{\Delta t}T_P^n
  + \Phi V_s + (h_jT_j+h_bT_b)V_s$$

for transient runs, where $V_s=R\,\Delta\theta\,\Delta z$ is the cell surface
area used by the thin-film FVM operators and $\mathcal{L}$ includes thermal
diffusion and upwind advection by film-averaged velocities. The advection term
is assembled in non-conservative form, equivalent to:

$$\nabla_s\cdot(\rho c_p h\bar{\mathbf{u}}T)
  - T\nabla_s\cdot(\rho c_p h\bar{\mathbf{u}})$$

This is important near supply grooves and pressure boundaries, where the
Reynolds pressure solve has local mass source/sink behavior. Using the
conservative form alone adds a spurious $T\nabla_s\cdot(\rho c_p h\bar{\mathbf{u}})$
term and can generate nonphysical hot/cold boundary cells.

Inlet features come in two flavors. An *open region* (default, e.g. a
large-clearance pocket of the same recirculating oil) only constrains pressure
or film content in the Reynolds solve and is thermally transparent to the energy
equation. An *actual fed inlet* (`inlet_* ... t_supply`) additionally pins those
cells to a fresh-oil supply temperature `t_supply`.

The axial boundary heat fluxes and pressure-gradient diagnostics use the
**same effective boundary pressure the flow solver imposed**. Under
Elrod-Adams, the solver derives boundary film content from `bc_z_*_val`,
`p_cav`, and `bulk_modulus`, so the effective pressure is
$p_{bc}=\max(p_{value},p_{cav})$. The Gumbel pressure-Dirichlet path uses
`bc_z_*_val` directly. This keeps pressure, velocity, force, and thermal
post-processing on one boundary contract and removes manual film-content
boundary inputs from the physical setup.

Axial boundaries use upwind thermal convection driven by that consistent
boundary pressure. Outflow is always zero-gradient (the leaving oil carries its
cell temperature). The inflow temperature is selected per side by
`bc_z_*_thermal`: `OPEN` (submerged / same-oil) treats inflow as zero-gradient
too — the entering oil is the same recirculating film oil, so no external
temperature is imposed; `RESERVOIR` (fed inlet/outlet) carries
`temperature_reference` in on inflow. `INLET_OUTLET` axial pressure boundaries
use the same pressure-gradient treatment for boundary heat-generation
diagnostics as the thermal boundary flux, so near-boundary viscous heating is
not silently dropped. With
`solution_mode = STEADY_STATE`, the transient capacity term is omitted and the
main solver runs one pressure/thermal step at $t=0$. With
`temperature_model = ISOTHERMAL`, `temperature` remains fixed at
`temperature_initial` while `heat_generation` is still refreshed for output.

Transient simulations start from a uniform outlet-compatible pressure field.
For pressure cavitation models, the initial pressure is the first non-Neumann
axial outlet pressure. For Elrod-Adams, the outlet boundary is a film-content
condition, so the initialized pressure is recovered from the EOS:

$$p^0 = p_{cav} + \beta\ln(\theta_{bc}), \qquad \theta^0 =
\max\left[\exp\left(\frac{p^0-p_{cav}}{\beta}\right),\theta_{min}\right].$$

The first transient solve uses the same transient terms as later solves. To
avoid an instantaneous jump from a static initialized state to full surface
speed, the configured shaft speed can be ramped:

$$\omega(t) = \omega_{target}\,
\min\left(\max\left(\frac{t}{t_{ramp}},0\right),1\right),$$

where `omega_ramp_time = t_ramp`. If `omega_ramp_time <= 0`, transient runs use
full `omega` immediately. `solution_mode = STEADY_STATE` ignores the ramp and
uses the target speed for its single operating-point solve.

For moving-bearing Elrod-Adams solves, the transient term
$\partial(\rho_l\theta h)/\partial t$ is expanded as
$\rho_l h\,\partial\theta/\partial t + \rho_l\theta\,\partial h/\partial t$.
The weighted capacity term $\rho_l h\,\partial\theta/\partial t$ is always
assembled for transient moving-bearing pressure solves, including when
`pressure_time_method = EULER_EXPLICIT`. Omitting that storage term while
retaining the lagged squeeze source would not conserve liquid content and can
create artificial startup cavitation. The explicit part is the lagged
$\theta^n\,\partial h/\partial t$ source.

### 13.2 Temperature-dependent Fluid Properties

Thermal hydrodynamic coupling enters Reynolds through field-valued liquid
solution density and viscosity. The Reynolds diffusion coefficient becomes:

$$\gamma(\theta,z) = \frac{\rho(T,p)h^3}{12\mu(T,p)}$$

and the Elrod-Adams diffusion coefficient becomes:

$$\Gamma_g(\theta,z) =
  \frac{\rho(T,p)\beta h^3}{12\mu(T,p)}g(\theta_{film})$$

The constant-property limit must reduce exactly to the current isothermal
solver.

### 13.2.1 Oil plus Dissolved-Gas Property Layer

The implemented fluid-property layer separates the liquid solution from free
gas cavitation:

- `rho_liquid_solution`, `mu_liquid_solution`, `cp_liquid_solution`, and
  `k_liquid_solution` describe the liquid oil plus dissolved gas.
- `theta` remains the Elrod-Adams/JFO liquid-content variable.
- `alpha_gas` represents released free gas bubbles transported with the film.
  It is not forcibly clipped to the JFO void fraction, so gas can persist into
  pressure-recovery/full-film cells and disappear only by finite-rate
  resorption or transport.

The default `fluid_property_model = CONSTANT` sets:

$$\rho_l = \rho_{oil}, \qquad \mu_l = \mu_{oil}, \qquad c_d = 0$$

and the Reynolds/JFO equations recover the previous pure-oil behavior.

For `OIL_DISSOLVED_GAS` or `GAS_CAVITATION_MIXTURE`, dissolved gas is stored as
a liquid-solution mass fraction:

$$c_d = \frac{m_{gas,dissolved}}{m_{liquid\ solution}}$$

The equilibrium concentration is selected by `oil_gas_solution_model`:

$$c_{d,eq} =
\begin{cases}
H_{ref}\exp\!\left[E_H\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right]p,
  & \text{HENRY}\\
B \frac{p}{p_{ref}}\frac{T_{ref}}{T}\frac{\rho_{g,ref}}{\rho_{oil}}, & \text{BUNSEN}\\
\text{linear table}(p), & \text{TABLE}
\end{cases}$$

with propane/R290 free-gas density evaluated by the ideal-gas relation:

$$\rho_g = \frac{\max(p,p_{floor})}{R_g T}$$

where $R_g = 188.55\ \mathrm{J/(kg\,K)}$ for propane and
$287.05\ \mathrm{J/(kg\,K)}$ for air.

The liquid-solution density and viscosity are selected independently:

$$\rho_l =
\begin{cases}
\rho_{oil}, & \text{PURE\_OIL}\\
\left(\frac{1-c_d}{\rho_{oil}} + \frac{c_d}{\rho_{dg,liq}}\right)^{-1}, & \text{MASS\_VOLUME\_MIXING}\\
\text{linear table}(c_d), & \text{TABLE}
\end{cases}$$

where `dissolved_gas_liquid_density` supplies $\rho_{dg,liq}$. If it is zero,
the code falls back to ideal-gas density for backward compatibility, but this
is not recommended for dissolved propane in oil.

For `EMPIRICAL_CORRELATION`, the liquid-solution viscosity is:

$$\mu_l =
\mu_{oil,ref}
\exp\!\left[E_\mu\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right]
\exp(a_c c_d)
\exp\!\left[\alpha_p(p-p_{ref})\right]$$

$$\mu_l =
\begin{cases}
\mu_{oil}, & \text{PURE\_OIL}\\
\mu_{oil}\exp(a_c c_d), & \text{LOG\_MIXING}\\
\mu_{oil,ref}\exp\!\left[E_\mu\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right]
  \exp(a_c c_d)\exp\!\left[\alpha_p(p-p_{ref})\right],
  & \text{EMPIRICAL\_CORRELATION}\\
\text{linear table}(c_d), & \text{TABLE}
\end{cases}$$

The field written as `mu_liquid_solution` is this liquid-solution viscosity.
For `GAS_CAVITATION_MIXTURE`, the effective film viscosity written as `mu` is
selected by `gas_mixture_viscosity_model` (all satisfy $\bar\mu(\alpha_g{=}0)=\mu_l$;
mass quality $x_g = m_g/(m_g + \rho_l\theta h)$):

$$\bar\mu =
\begin{cases}
\mu_l(1 + 2.5\,\alpha_g), & \text{EINSTEIN\_DILUTE}\\
\alpha_g\,\mu_g + (1-\alpha_g)\,\mu_l, & \text{DUKLER\_VOID}\\
\left(\dfrac{x_g}{\mu_g} + \dfrac{1-x_g}{\mu_l}\right)^{-1}, & \text{MCADAMS\_QUALITY}\\
\mu_l\,(1-\alpha_g/\alpha_{g,max})^{-2.5\,\alpha_{g,max}}, & \text{KRIEGER\_DOUGHERTY}\\
x_g\,\mu_g + (1-x_g)\,\mu_l, & \text{LINEAR\_QUALITY}
\end{cases}$$

`EINSTEIN_DILUTE` (the backward-compatible default) is a dilute-suspension
result valid only for $\alpha_g \lesssim 0.1$ and is directionally wrong as
$\alpha_g \to 1$ (it thickens; gas must thin the film toward $\mu_g$).
Validation therefore rejects it when `gas_alpha_max > 0.6` (Krieger-Dougherty
packing bound). `DUKLER_VOID` is the void-weighted linear analogue of the
Grando, Priest & Prata (2006) density closure; `MCADAMS_QUALITY` is the
homogeneous two-phase standard; `LINEAR_QUALITY` is the quality-weighted form
quoted by Grando 2006 ($\bar\mu = \chi\mu_g + (1-\chi)\mu_l$) for exact
replication of that reference. The quality/void models require the free-gas
viscosity `mu_gas`. The shipped R290/PZ68 cases select `MCADAMS_QUALITY` with
`gas_alpha_max = 0.6`.

For dissolved-gas (solution) viscosity beyond the calibration point, the
validated representation for R290/oil pairs is the Eyring-MTSM activity model
(measured 303-348 K isotherms); the implemented linear Henry law and
$\exp(a_c c_d)$ correction are local tangents around the (10 bar, 40 °C)
calibration point. Measured isotherms should be supplied through the existing
`solubility_table` / `density_table` / `viscosity_table` keys; the shipped
R290/PZ68 tables are still pending real data and the configs say so.

The density used by JFO is:

$$\rho = \rho_l \theta$$

so liquid mass conservation remains tied to the Elrod-Adams film-content
variable.

When `fluid_property_model = GAS_CAVITATION_MIXTURE`, finite-rate
release/resorption updates the dissolved and free gas diagnostics:

$$\Delta m_g =
k_{dg}(c_d - c_{d,eq})\,\rho_l h\,\alpha_l\,\Delta t$$

with sign convention `gas_mass_transfer > 0` for release from solution into
free gas and `gas_mass_transfer < 0` for resorption. The implementation caps
the reported/resolved free-gas volume fraction by `gas_alpha_max`:

$$0 \le \alpha_g \le \alpha_{g,max}$$

Dissolved-gas transport is a segregated advection-diffusion step on $c_d$ using
the same liquid film mass flux as the JFO equation. Free gas mass per bearing
area is transported conservatively with the same film-averaged velocity and no
slip. Both transport solves are skipped when `dt <= 0` (steady operating-point
runs) to avoid solving boundaryless steady advection problems.

Pressure supply/inlet cells are treated as reservoir-composition cells for the
gas model. They impose the configured incoming liquid mixture,

$$c_d = c_{d,initial}, \qquad m_{g,free}=0, \qquad \alpha_g=0,$$

and the diagnostic source is reset to

$$\dot{m}_{dg}=0.$$

This prevents a pressure inlet from deleting free gas silently or reporting
local release/resorption inside the imposed reservoir region. Axial pressure
boundaries are open gas-transport boundaries: outflow convects dissolved and
free gas out of the domain, while inflow carries `dissolved_gas_initial` and
zero free gas. A non-inlet cell can still show `gas_mass_transfer > 0` at high
absolute pressure when the local temperature or tabulated/correlation model
makes the liquid supersaturated, i.e. when $c_d > c_{d,eq}(p,T)$.

The provided `config_r290_pz68.txt` fills the 40 C, 10 bar R290/PZ68 starting
point explicitly:

| Config key | Value | Meaning |
|------------|-------|---------|
| `p_cav` | `3e4 Pa` | Absolute JFO cavitation plateau pressure. This is not a force reference pressure and must not be `0` gauge. See the cavity-pressure caveat in §6 — measured cavity pressures of gas-laden oils are not constant. |
| `bc_z_south_val`, `bc_z_north_val` | `1e6 Pa` | Absolute boundary pressure: submerged ends at the 10 bar housing pressure. Because $p_{sat}(T_{init}, c_{d,init}) = 1\,\mathrm{MPa}$, boundary oil enters exactly saturated and releases gas on any pressure drop — startup validation logs this ratio. The quick smoke case (`config_r290_pz68_quick.txt`) instead vents at `1e5 Pa` with `p_cav = 1e5 Pa`. |
| `property_reference_pressure` | `1e6 Pa` | Absolute 10 bar reference point for the Henry and Barus fits. |
| `dissolved_gas_henry_coeff` | `2.175e-7 1/Pa` | Chart-calibrated $H_{ref}=0.2175/10^6$ for 21.75% propane solubility. |
| `dissolved_gas_henry_temp_coeff` | `1800 K` | van't Hoff temperature slope; solubility decreases as temperature rises. |
| `dissolved_gas_liquid_density` | `500 kg/m^3` | Liquid-phase density used for dissolved propane volume mixing. |
| `solution_viscosity_gas_coeff` | `-11.40` | Log-mixing coefficient calibrated so the solution viscosity is about `6.843 cSt` at 40 C and 10 bar. |
| `viscosity_temperature_coeff` | `4000 K` | Andrade oil-viscosity temperature slope for ISO VG 68 trend. |
| `viscosity_pressure_coeff` | `1.5e-8 1/Pa` | Barus pressure-viscosity coefficient. |
| `gas_mass_transfer_rate` | `500 1/s` | Release/resorption rate; numerical/kinetic model input, not read from the chart. |
| `dissolved_gas_diffusivity` | `1e-9 m^2/s` | Effective dissolved-gas diffusivity; regularization/model input, not read from the chart. |

The CLI and GUI validate these dependencies before launch. For example, HENRY
requires `dissolved_gas_henry_coeff > 0`, MASS_VOLUME_MIXING requires
`dissolved_gas_liquid_density > 0`, empirical/log viscosity requires a non-zero
`solution_viscosity_gas_coeff`, and `GAS_CAVITATION_MIXTURE` currently requires
`ELROD_ADAMS` so gas release is coupled to the mass-conserving JFO
liquid-content solve.

This is a segregated nonlinear JFO/property/gas update, not a PISO method. A
true PISO label would require a Navier-Stokes pressure-velocity backend.

### 13.3 Elastic and Thermal Deformation

The effective film thickness for EHD/ETHD is:

$$h_{eff} = h_{kinematic} + \delta_p + \delta_T$$

where $\delta_p$ is pressure-induced elastic deformation and $\delta_T$ is
thermal expansion. A local-compliance first implementation is:

$$\delta_p(\theta,z) = C_p(\theta,z)[p(\theta,z)-p_{ref}]$$

and:

$$\delta_T(\theta,z) = \alpha_s t_s [T_s(\theta,z)-T_{ref}]$$

where $C_p$ is local structural compliance, $\alpha_s$ is solid thermal
expansion coefficient, and $t_s$ is an effective solid thickness. This can later
be replaced by an influence-matrix convolution or an external structural solve
without changing the Reynolds interface, as long as it returns a deflection
field on the same mesh.

### 13.4 Coupled ETHD Iteration

The planned ETHD timestep uses a Picard loop:

1. Update moving-bearing state and kinematic film thickness.
2. Compute elastic and thermal deformation.
3. Update $h_{eff}$.
4. Update $\mu(T,p)$ and any density or mixture properties.
5. Solve Reynolds/Elrod for pressure or film content.
6. Compute velocities and viscous heat generation.
7. Solve the energy equation.
8. Repeat until pressure, temperature, viscosity, and film thickness residuals
   satisfy the configured coupled tolerance.

Thermal-only, deformation-only, and fully coupled ETHD modes should remain
separately switchable so each coupling path can be validated in isolation.

---

## 14. GUI and Output Selection

The Windows GUI is an orchestration layer, not a separate numerical model. It
writes the same `SimulationConfig` values consumed by the CLI solver, launches
`pancake.exe`, and reads the solver's ASCII VTK output for preview. Therefore
the governing equations, discretisation, cavitation model, and macroscopic
property calculations above are unchanged by using the GUI.

The live summary reports a nominal circumferential Couette Courant number,

$$Co_\theta = \frac{|\omega|\,\Delta t}{2\,\Delta\theta}, \qquad
\Delta\theta = \frac{2\pi}{N_\theta}.$$

This is the cell-crossing ratio for the surface-driven Reynolds advection term.
It is an input-summary diagnostic only; it does not include pressure-driven
Poiseuille velocities, which are known only after the pressure field has been
solved.

The GUI's Motion / Loads section exposes the same moving-bearing controls that
can be written manually in `config.txt`. The applied load entry is a clean input
form for the vector $\mathbf{F}^{ext}$: the in-plane magnitude and direction are
converted to `external_load_x` and `external_load_y`, while the z entry writes
`external_load_z`. This conversion is purely an input convenience; the bearing
equation of motion still uses the Cartesian external-load vector defined in
Section 12.

The GUI's Energy section exposes `temperature_model`, thermal wall
temperatures, wall heat-transfer coefficients, volumetric heat capacity, and
thermal conductivity. These controls write the same config keys used by the CLI
solver. The preview unit labels report `temperature` in K and `heat_generation`
in W/m^2.

The output-selection parameters only affect persistence:

| Parameter | Effect |
|-----------|--------|
| `output_write_3d` | Writes the cylindrical structured-grid VTK/PVD result tree. |
| `output_write_flat` | Writes the unwrapped flat structured-grid VTK/PVD result tree. |
| `output_fields` | Selects which cell-data arrays are written. Coordinate helper arrays `theta_rad` and `z_m` are still written with each enabled VTK file. |

Disabling a field does not remove it from the simulation state. For example,
unchecking `velocity` skips the visualization vector, but the solver still
computes velocities before force and torque post-processing. When enabled,
`velocity` is written as the active 3-component VTK cell vector `U`. Curved
output stores Cartesian components on the cylinder; flat output stores
`(u_theta,u_z,0)` on the unwrapped plane. Force/load/bearing component fields
are grouped into vector arrays `Fp`, `Fv`, `F`, `Fext`, and `xB`.
Legacy `load_x`, `load_y`, and `load_z` output selections write
`pressure_force_x`, `pressure_force_y`, and `pressure_force_z` respectively.
Disabling both
`output_write_3d` and `output_write_flat` runs the numerical model without VTK
field persistence.

---

## 15. Conservation, Regime, and Convergence Diagnostics

`src/diagnostics.hpp/cpp` computes per-step global mass balances on the
physical domain ($dA = R\,d\theta\,dz$):

$$M_l = \int_\Omega \rho_l\,\theta\,h\,dA, \qquad
M_{gas,tot} = \int_\Omega \big(c_d\,\rho_l\,\theta\,h + m_g\big)\,dA$$

The liquid closure residual mirrors the discrete operators the Elrod solve
actually assembles (units kg/s per cell): the storage term is formed from
$\theta$ vs $\theta^n$ with the same $\rho_l h$ capacity weight and squeeze
source, boundary fluxes use the same Dirichlet half-cell formulas built on
$\Gamma_{base}$, and inlet penalty cells are re-evaluated as physical operator
residuals (the mass the penalty row injects). The normalized residual

$$r_l = \frac{\Delta t\,\big(S_{storage} - \Phi_{boundary} - S_{inlet}\big) - M_{clamp}}{\max(M_l, \epsilon)}$$

sits at the linear-solver tolerance (`linear_rtol`) in closed and open
domains — $O(10^{-12})$ with `linear_rtol = 1e-13` per `tests/test_diagnostics.cpp`,
meeting the PLAN A.2 target. Clamp and snap operations are **not mass-neutral**;
their per-cell mass deltas are accumulated into the `cavitation_clamp_mass` and
`gas_clamp_mass` fields by the Reynolds and gas-transport solves and appear as
explicit terms in the balance. The per-step results (masses, fluxes, inlet
source, clamp masses, residuals, Elrod outer iterations and flag flips) are
appended to `<output_dir>/diagnostics.csv` every `diagnostics_interval` steps
and summarized on the DETAILED log line.

Two deliberate visibility limits remain (they motivate WP-1/WP-2 of the audit
plan): variable-property runs refresh $\rho_l$ after the segregated solve, so
$r_l$ then measures the property lag of the sweep; and the dissolved-gas
transport is solved in material (non-conservative) form, so $r_{gas}$ in
advecting cases includes that discretization choice. The closed-cell
release/resorption exchange cancels exactly in $M_{gas,tot}$.

Startup regime guards in `SimulationConfig::validate()`:

- **Laminar validity**: $Re_c = \rho\,\omega R\,c/\mu$ and Taylor number
  $Ta = Re_c\sqrt{c/R}$; a warning is issued above the Taylor-vortex onset
  $Ta \approx 41$ (Hamrock; San Andrés), where the laminar Reynolds equation
  stops being valid.
- **Saturation context** (gas models): the saturation pressure
  $p_{sat}(T, c_d) = c_d / H(T)$ (Henry inversion; Bunsen and table models are
  inverted analogously via `SimulationConfig::saturation_pressure`) is logged
  together with $p_{cav}/p_{sat}$, and a warning fires when a boundary or inlet
  feed pressure is at or above $p_{sat}$ — that oil enters saturated and
  releases gas on any pressure drop. For the shipped R290/PZ68 constants
  $p_{sat} = 0.2175/2.175\times10^{-7} = 1.0\ \mathrm{MPa}$, exactly the 10 bar
  feed pressure.

---

## 16. Oil–Refrigerant Gaseous Cavitation (PZ68S + R290)

This section records the *physical model* for cavitation of an oil with a
dissolved refrigerant (the target case is PAG-class PZ68S oil saturated with
R290/propane). It is the canonical summary; the full survey, the
literature-verification ledger, and the cited sources live in
`docs/CAVITATION_OIL_REFRIGERANT.md`. The implementation status and the
work-package breakdown are in `PLAN.md` Phase E (WP-12 data backend + onset,
WP-1 two-phase pressure coupling).

### 16.1 Which species cavitates — the two-threshold rule

A lubricant film has **two** low-pressure release thresholds. As pressure falls
in the divergent wedge, the species with the **higher release pressure** comes
out first. For oil + dissolved R290:

- **Dissolved refrigerant outgasses first**, at the solution **bubble-point**
  $p_{sat}(T, c_d)$ — *gaseous cavitation*. For the shipped R290/PZ68 state
  $p_{sat}\approx 1.0$ MPa.
- **Oil vaporization** would require reaching the oil's own vapour pressure
  (essentially zero for nonvolatile PAG) — *vaporous cavitation*, effectively
  never reached; kept off by default.

So the operative mechanism across the whole realistic $(p,T)$ map is
**gaseous/outgassing cavitation governed by $p_{sat}$**, not boiling of the oil.
A single fixed $p_{cav}$ (the classical JFO threshold) cannot represent this:
the onset must track $p_{sat}$, which depends on local temperature and dissolved
fraction. This is exactly the defect in the shipped config, where the Elrod
threshold is pinned at the oil value $p_{cav}=3.0\times10^4$ Pa while propane
releases at $\sim$1 MPa.

### 16.2 Phase picture across the PTSV map

The working fluid is **gas dissolved in liquid oil**, becoming a two-phase film
once the solubility limit is crossed — neither a pure liquid–liquid mixture nor
a pure gas-in-liquid problem, but a transition keyed on the saturation surface:

| Local state | Phase | Closure |
|---|---|---|
| $p > p_{sat}$ (compressed) | single-phase liquid solution (R290 fully dissolved) | mixture $\mu_{sol}(T,p,c_d)$, $\rho_{sol}(T,p,c_d)$; elliptic Reynolds |
| $p \approx p_{sat}$ | bubble-point onset | $p_{sat}$ from solubility inversion |
| $p < p_{sat}$ (divergent) | two-phase: desorbed R290 gas + liquid | void fraction $\alpha = 1-\theta$; two-phase $\bar\mu(\alpha)$, $\bar\rho(\alpha)$; cavitated branch |

### 16.3 Generalized JFO complementarity

The classical JFO complementarity holds with $p_{cav}$ promoted to the local
bubble-point field:

$$p \ge p_{sat}(T,c_d),\qquad \theta \le 1,\qquad (p - p_{sat})\,(1-\theta)=0,$$

i.e. either full film ($p\ge p_{sat}$, $\theta=1$) or cavitated ($p=p_{sat}$,
$\theta<1$). The effective onset threshold used by the solver is the **larger**
of the two species' release pressures (general, not R290-specific):

$$p_{cav,\mathrm{eff}} = \max\!\big(p_{cav,\mathrm{oil}},\; p_{sat}(T,c_d)\big).$$

With WP-1 void coupling the cavitated plateau generalizes further to the
gas-filled value $p_{void}=\max(p_{cav,\mathrm{eff}},\, m_g R_g T/((1-\theta)h))$,
of which $p_{cav,\mathrm{eff}}$ is the gas-free floor.

### 16.4 Property closures and the supplier PTSV data

The supplier chart (`data/R290_PZ68S/r290_pz68s.csv`) is a $(p,T)$ grid of the
**saturated** solution's solubility (mass fraction R290) and **kinematic**
viscosity. It is ingested as the liquid-solution surface:

- **Solubility** $c_{sat}(p,T)$ sets the dissolved fraction *and*, inverted at
  fixed $T$, the onset threshold $p_{sat}(T,c_d)$.
- **Viscosity**: tabulated $\nu(p,T)$ → dynamic $\mu=\nu\,\rho_{sol}$ per cell
  (the chart has no density column; $\rho_{sol}$ comes from the density model).
- **Density** of the *solution* is mass–volume mixing of oil and dissolved-R290
  partial volumes; the **two-phase** density $\bar\rho=(1-\alpha)\rho_l+\alpha\rho_g$
  is carried by the Elrod $\theta$.
- **Two-phase viscosity** uses a homogeneous quality-based model (McAdams) with a
  Krieger–Dougherty packing bound; it reduces to $\mu_l$ at zero void.

The dissolved-gas (solution) viscosity below saturation is, as a local-tangent
fallback to the table, $\mu_l(T,c_d,p)=\mu_{ref}\,e^{E_\mu(1/T-1/T_{ref})}\,
e^{a_c c_d}\,e^{\alpha_p(p-p_{ref})}$ (Andrade temperature × dissolved-gas
thinning × Barus pressure), calibrated at (40 °C, 1.0 MPa).

### 16.5 Model provenance and honest validity

This year's model corrects three errors in the earlier (Kunz-vaporous) approach:
mass-conserving Elrod–Adams JFO replaces non-conservative volume-fraction
manipulation; the dissolved refrigerant is treated as **liquid** (its partial
molar volume is liquid-like) rather than as an ideal gas; and mixing moved from
volume-fraction to **mass-fraction (quality)**, which satisfies both single-phase
limits and stays physical at high void. Honest gaps (tracked as WP-12 and WP-1):
the linear Henry solubility is a one-point tangent unverified at high pressure
(→ the 2-D TABLE backend), the Barus $\alpha_p$ is a mineral-oil default not a
PZ68 measurement, the release-rate constant is provisional, and the desorbed gas
is currently **volumetrically inert to pressure** (→ WP-1). The supplier PTSV
data itself is a chart-digitized estimation (R² ≈ 0.9995 over 0.1–4 MPa,
−60…160 °C), not first-principles VLE, and carries **no density**.
