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

### 5.5 Convection (TVD)

`FVM::divergence` assembles $\text{div}(F\,\phi)$ for a known face flux field $F$. The upwind (first-order) scheme adds implicit contributions; TVD schemes (van Leer, MINMOD) add an explicit deferred correction using $\phi^{\text{old}}$.

---

## 6. Boundary Conditions

### Circumferential ($\theta$)
**Periodic**: ghost cells are exchanged via MPI `Sendrecv` between adjacent ranks (rank 0 ↔ rank $N_{ranks}-1$ for the wrap-around).

### Axial ($z$)
**Dirichlet** $p = p_{cav}$ at $z = 0$ and $z = L$ (open bearing ends). Implemented by adding the boundary face flux to the stencil using a half-cell distance:

$$a_P \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}, \qquad \text{source} \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}\,p_{cav}$$

---

## 7. Reynolds Solvers

Two solvers are available, selected by `cfg.cavitation_model`.

Before each solve, the flooded initial guess is set from the first configured
inlet supply pressure:

$$p^0 = p_{inlet},\qquad
\theta^0 = \max\!\left[\theta_{film}(p_{inlet}),\,\theta_{min}\right]$$

If no inlet is configured, the fallback is $p^0=p_{cav}$ and the corresponding
full-film value $\theta^0=1$.

### 7.1 Gumbel Solver (`Reynolds::solve`)

The solver follows a **lagged-coefficient** approach at each time step:

1. Compute $\gamma^n = \rho^n h^3 / (12\mu)$ and $(\rho h)^n$ at physical + theta ghost cells (using already-synced ghost values of $\rho$ and $h$).
2. Reset and assemble the linear system:
   - `FVM::laplacian` with $\gamma^n$
   - Dirichlet z-boundary corrections
   - Couette source: $-\frac{\omega}{2}\frac{\partial(\rho h)^n}{\partial\theta}$ (centered difference, sign convention above)
3. Solve $A\,p^{n+1} = b$ via PETSc KSP (BiCGStab + block Jacobi).
4. Apply Gumbel cavitation: $p^{n+1} \leftarrow \max(p^{n+1},\, p_{cav})$.
5. Update $\rho^{n+1} = \rho_0\,\theta_{film}(p^{n+1})$ via EOS.

The transient $\partial(\rho h)/\partial t$ term is omitted for the static-geometry steady-state solve.

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

**Relationship to Gumbel solver.** Both formulations are equivalent in the full-film ($\theta \ge 1$) limit at small compressibility. Phase A.1 validation shows peak pressures agree to within $\approx 1.2\%$ (residual comes from centered-difference vs upwind Couette discretisation). Phase A.2 will add the switch function $g(\theta)$ to handle the cavitated region.

**Weighted time derivative.** The new FVM operator `FVM::ddt_weighted(sys, theta, h_field, dt, mesh)` discretises $\partial(\theta h)/\partial t = h\,\partial\theta/\partial t$ for static geometry. It adds $h\,V/\Delta t$ to $a_P$ and $h\,V\,\theta^n/\Delta t$ to `source`. Not yet active in Phase A.1 (steady-state).

---

## 8. Linear System

`LinearSystem` wraps a PETSc `MPIAIJ` matrix. Each physical cell $(i,j)$ maps to global row $(\text{offset}_\theta + i)\cdot N_z + j$. Periodic $\theta$-wrapping is applied in the column indices. The system is solved with:
- **KSP**: BiCGStab (`KSPBCGS`)
- **PC**: Block Jacobi (`PCBJACOBI`)
- **Tolerance**: relative $10^{-8}$, max 1000 iterations

---

## 9. MPI Parallelism

The domain is decomposed in the $\theta$ direction: each MPI rank owns $N_\theta / N_{ranks}$ contiguous $\theta$-columns with the full $z$ extent. Ghost layers (default 2 cells) are exchanged each step via `Communicator::update_ghosts`, which packs/unpacks the non-contiguous ghost strips using `MPI_Sendrecv`.

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
`pressure_force_z`. The legacy fields `load_x`, `load_y`, and `load_z` remain as
aliases for the same pressure-only resultant.

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
motion. The motion integrator supports `EULER_EXPLICIT`, `EULER_IMPLICIT`,
`CRANK_NICOLSON`, `RK2`, and `RK4`. `CRANK_NICOLSON` is the formal name for the
mixed implicit/explicit trapezoidal update; the parser also accepts aliases such
as `SEMI_IMPLICIT`, `IMEX`, and `CRANK_NICHOLSON`.

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

The finite-volume equation solved for `temperature` is:

$$\frac{\rho c_p h V_s}{\Delta t}T_P^{n+1}
  + \mathcal{L}(T^{n+1}) + (h_j+h_b)V_s T_P^{n+1}
  =
  \frac{\rho c_p h V_s}{\Delta t}T_P^n
  + \Phi V_s + (h_jT_j+h_bT_b)V_s$$

for transient runs, where $V_s=R\,\Delta\theta\,\Delta z$ is the cell surface
area used by the thin-film FVM operators and $\mathcal{L}$ includes thermal
diffusion and upwind advection by film-averaged velocities. With
`solution_mode = STEADY_STATE`, the transient capacity term is omitted and the
main solver runs one pressure/thermal step at $t=0$. With
`temperature_model = ISOTHERMAL`, `temperature` remains fixed at
`temperature_initial` while `heat_generation` is still refreshed for output.

### 13.2 Temperature-dependent Fluid Properties

Thermal hydrodynamic coupling enters Reynolds through field-valued viscosity
and density. A configurable first model is:

$$\mu(T,p) = \mu_{ref}
  \exp[-a_T(T - T_{ref}) + a_p(p - p_{ref})]$$

For a purely thermal viscosity model, set $a_p=0$. The Reynolds diffusion
coefficient becomes:

$$\gamma(\theta,z) = \frac{\rho(T,p)h^3}{12\mu(T,p)}$$

and the Elrod-Adams diffusion coefficient becomes:

$$\Gamma_g(\theta,z) =
  \frac{\rho(T,p)\beta h^3}{12\mu(T,p)}g(\theta_{film})$$

The constant-property limit must reduce exactly to the current isothermal
solver.

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
solver. The preview unit labels report `temperature` in K and
`heat_generation` in W/m^2.

The output-selection parameters only affect persistence:

| Parameter | Effect |
|-----------|--------|
| `output_write_3d` | Writes the cylindrical structured-grid VTK/PVD result tree. |
| `output_write_flat` | Writes the unwrapped flat structured-grid VTK/PVD result tree. |
| `output_fields` | Selects which cell-data arrays are written. Coordinate helper arrays `theta_rad` and `z_m` are still written with each enabled VTK file. |

Disabling a field does not remove it from the simulation state. For example,
unchecking `velocity` skips the visualization vector array but the solver still
computes velocities before load and torque post-processing. Disabling both
`output_write_3d` and `output_write_flat` runs the numerical model without VTK
field persistence.
