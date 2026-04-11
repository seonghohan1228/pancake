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

with $a_P = a_E + a_W + a_N + a_S$ (plus any boundary contributions). This stencil represents **the negative divergence times cell area**:

$$-\text{div}(\gamma\,\nabla p)\cdot A_{cell} = \text{source}$$

**Sign convention (critical):** For a PDE $\text{div}(\gamma\,\nabla p) = f$, the source term passed to `FVM::add_source` must be $-f$, not $+f$. For the Reynolds equation with Couette source $S$:

$$\text{div}(\gamma\,\nabla p) = S \implies \texttt{source\_field} = -S$$

### 5.4 Time Derivative

The `FVM::ddt` operator discretises $\partial\phi/\partial t$ with a first-order backward Euler scheme, adding:

$$a_P \mathrel{+}= \frac{A_{cell}}{\Delta t}, \qquad \text{source} \mathrel{+}= \frac{A_{cell}}{\Delta t}\,\phi^n$$

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

## 7. Reynolds Solver (`Reynolds::solve`)

The solver follows a **lagged-coefficient** approach at each time step:

1. Compute $\gamma^n = \rho^n h^3 / (12\mu)$ and $(\rho h)^n$ at physical + theta ghost cells (using already-synced ghost values of $\rho$ and $h$).
2. Reset and assemble the linear system:
   - `FVM::laplacian` with $\gamma^n$
   - Dirichlet z-boundary corrections
   - Couette source: $-\frac{\omega}{2}\frac{\partial(\rho h)^n}{\partial\theta}$ (centered difference, sign convention above)
3. Solve $A\,p^{n+1} = b$ via PETSc KSP (BiCGStab + block Jacobi).
4. Apply Gumbel cavitation: $p^{n+1} \leftarrow \max(p^{n+1},\, p_{cav})$.
5. Update $\rho^{n+1} = \rho_0\,\theta_{film}(p^{n+1})$ via EOS.

The transient $\partial(\rho h)/\partial t$ term is omitted for the static-geometry steady-state solve (it is zero for constant $h$ and converged $\rho$).

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
