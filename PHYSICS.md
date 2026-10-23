# Physics and Numerical Methods

Current mixture formulation: [coupled mass and volume](docs/COUPLED_CAVITATION.md),
[thermal energy](docs/COUPLED_ENERGY.md), and [R410A EOS](docs/R410A_EOS.md).
These describe `cavitation_model = COUPLED`, including its phase and material
limits, and supersede the earlier split-mixture descriptions below. Historical
sections are retained for understanding previous results and other solvers.
Neither proposed FW68D mixture is experimentally validated.

For the current choices across both new and older solver paths, use the
[model-selection guide](docs/MODEL_SELECTION.md). It identifies active models,
incompatible options, branch-specific Ferron work and proposed models that
are not implemented. A smooth load is not a physical validation criterion.

## Current coupled formulation (2026-09-07)

The state is absolute pressure $p$, oil inventory $M_o$, dissolved inventory
$M_d$, free inventory $M_g$ (kg/m²), and optional common film temperature $T$.
Define $c_d=M_d/(M_o+M_d)$ and $\theta_l=(M_o+M_d)/(\rho_l h)$.
Oil and total dissolved/free species obey separate finite-volume conservation
equations. A third species equation supplies finite-rate release/absorption,
or local equilibrium complementarity when explicitly selected. The same
nonlinear system imposes

$$M_g=\rho_g(p-p_v(T),T)h(1-\theta_l),\quad p\ge p_v(T),\quad 0\le\theta_l\le1.$$

Gas-free cells recover JFO complementarity. A gas-containing cavity has
variable pressure from its inventory, temperature and available volume;
the initial dissolved charge's bubble point is not a universal cavity
pressure. `film_content` is the liquid volume fraction on this path, not
the compressible Elrod universal variable of the older solvers.

Thermal mode conserves $\mathcal U=\sum_i M_i u_i$ using the component
mass fluxes, compression work $-p(\partial_t h+\nabla_s\cdot q)$,
Couette/Poiseuille dissipation, wall heat transfer and in-plane conduction.
Supply mass carries enthalpy. The dissolved and free energy reference
difference accounts for desorption cooling/absorption heating once.
The linked derivations give the discrete equations, constitutive assumptions,
input conventions, tests and remaining thermal vapor limitations.

Axial pressure boundaries select reservoir temperature/composition only for
inflow. Outflow uses the interior state, including the mobility properties.
Internal pressure patches use their net mass source to select incoming or
outgoing composition and enthalpy. No patch is required for an end-fed bearing.
Grooves currently extend across the axial length; holes are circular film
patches. Neither resolves three-dimensional feed plumbing.

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

**Steady-state incompressible** (standard form): set $\rho = \text{const}$ and $\partial/\partial t = 0$. This gives the classical Reynolds equation; half-Sommerfeld pressure clipping is a separate, non-mass-conserving cavitation approximation.

---

## 4. Barotropic Equation of State and Cavitation

A **barotropic EOS** links density to pressure and handles sub-ambient (cavitated) regions via a **universal variable** $\theta_{film}$ (fractional film content):

| Region | Condition | Pressure | Density |
|--------|-----------|----------|---------|
| Full film | $\theta_{film} \ge 1$ | $p = p_{cav} + \beta\ln\theta_{film}$ | $\rho = \rho_0\,\theta_{film}$ |
| Cavitated | $\theta_{film} < 1$ | $p = p_{cav}$ | $\rho = \rho_0\,\theta_{film}$ |

Inverse relations:

$$\theta_{film}(p) = \begin{cases} \exp\!\left(\dfrac{p - p_{cav}}{\beta}\right) & p > p_{cav} \\ 1 & p \le p_{cav} \end{cases}$$

where $\beta$ is the bulk modulus. The **Gumbel (half-Sommerfeld) selection**
clamps pressure below $p_{cav}$ after the pressure solve. Other solver
selections impose different cavitation formulations as described in §7 and
the current model-selection guide; clipping is not the universal algorithm.

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
- `LINEAR` — second-order central interpolation $\phi_f = \tfrac12(\phi_U+\phi_D)$; sharpest in smooth regions but **unbounded** (can overshoot at fronts). A textbook reference scheme.
- `VANLEER` — van Leer TVD limited high-order deferred correction
  $\phi_f = \phi_U + \tfrac12\psi(r)(\phi_D-\phi_U)$ with
  $r = (\phi_U-\phi_{UU})/(\phi_D-\phi_U)$; sharper fronts, no new extrema (bounded).
  The triple $(\phi_{UU},\phi_U,\phi_D)$ is taken from the side the flux comes
  from, so it flips with the sign of $F_f$: on the east face of cell $i$ it is
  $(\phi_{i-1},\phi_i,\phi_{i+1})$ for $F_e \ge 0$ and
  $(\phi_{i+2},\phi_{i+1},\phi_i)$ for $F_e < 0$. This matters on the energy and
  gas-transport equations, where the Poiseuille term can reverse a face flux
  even with $\omega > 0$; on the Elrod $\theta$ equation the face flux is pure
  Couette and carries the sign of $\omega$ throughout.
  **Scope limit:** the deferred correction is applied on the circumferential
  faces only. Axial ($z$) faces remain first-order upwind under `VANLEER`.
- `TYPE_DIFFERENCING` (theta only) — Vijayaraghavan & Keith (1989): central
  interpolation across faces where both adjacent cells are full-film
  ($g = 1$), first-order upwind across and inside the cavitated region. This
  is the canonical Elrod-Adams discretization; full upwind everywhere is a
  conservative-but-diffusive default.

Validation (`tests/test_schemes.cpp`): step advection is bounded and sharper
under TVD than upwind; observed refinement orders are $\approx 1.09$ (upwind)
and $\approx 1.99$ (TVD van Leer) on a smooth profile with $\Delta t \propto \Delta x^2$.
Directional symmetry is pinned separately: the mirrored profile advected with
mirrored (negative) face flux reproduces the forward run's $L_1$ error to
$10^{-9}$ and stays inside the initial bounds.

---

## 6. Boundary Conditions

### Circumferential ($\theta$)
**Periodic**: ghost cells are exchanged via MPI `Sendrecv` between adjacent ranks (rank 0 ↔ rank $N_{ranks}-1$ for the wrap-around).

### Axial ($z$)
**Dirichlet** $p = p_{bc}$ at $z = 0$ and $z = L$, where $p_{bc}$ is the configured
absolute boundary pressure `bc_z_*_val`. Implemented by adding the boundary face
flux to the stencil using a half-cell distance:

$$a_P \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}, \qquad \text{source} \mathrel{+}= \frac{\gamma\,R\,\Delta\theta}{\Delta z/2}\,p_{bc}$$

A reservoir end uses the actual reservoir absolute pressure. A vented
reservoir is approximately at atmospheric pressure, not generally at the
oil vapor pressure. If $p_{bc}>p_{cav}$, the Elrod path converts it to boundary
content $\theta_{bc}=\exp[(p_{bc}-p_{cav})/\beta]>1$. `DIRICHLET` already allows
reverse flow; it does not prescribe flow direction. `NEUMANN` selects zero
pressure gradient. Legacy `INLET_OUTLET` has solver-dependent semantics and
is rejected by FBNS and COUPLED; see the current model-selection guide.

**Cavity-pressure caveat.** The model holds the cavitated-zone pressure at the
single constant `p_cav`. Measurements show the cavity pressure of gas-laden
oils is neither constant nor equal to the vapor pressure — it tracks the gas
release state and drops with speed (Etsion & Ludwig 1982; Braun & Hendricks
1984). A constant `p_cav` is therefore a model input to be chosen per case, and
the WP-1 void coupling (`PLAN.md` Phase E) is the planned mechanistic
replacement; see §16 for the oil–refrigerant gaseous-cavitation theory and
`docs/CAVITATION_OIL_REFRIGERANT.md` for the cited basis.

**Reformation reflood at a cavitated end (`consistent_boundary_flux`).** The
barotropic Dirichlet link above maps the boundary *film content* through the
$\theta_{bc}=\exp[(p_{bc}-p_{cav})/\beta]$ relation, which is only valid for a
full-film boundary cell. When a boundary cell is cavitated ($\theta<1$), the
physical process at a submerged end is **reformation**: reservoir liquid floods
back in at the Poiseuille rate set by the cavitation plateau,

$$\dot m_{\text{reflood}} = \rho_\ell\,\frac{h^3}{12\mu}\,\frac{p_{bc}-p_{cav}}{\Delta z/2}\,R\,\Delta\theta \;\ge\; 0 ,$$

an **inflow** (the bulk modulus cancels, so the rate is independent of $\beta$).
With `consistent_boundary_flux` the liquid balance, energy equation, and
dissolved-gas equation all use this single reflood rate
(`SimulationConfig::reformation_boundary_velocity`). The reflooding liquid is
reservoir liquid, so it convects the reservoir state into the cell: the
**reservoir temperature** `temperature_reference` (energy, RESERVOIR mode) and the
**reservoir composition** `dissolved_gas_initial` at the full liquid density
$\rho_\ell$ (gas). Without the flag, or for a full-film boundary cell, the
boundary face reduces to the barotropic Dirichlet link above (default, unchanged).

---

## 7. Reynolds Solvers

Five selections are available through `cfg.cavitation_model`: FULL_SOMMERFELD
and GUMBEL share the pressure-primary routine (§7.1), JFO uses Elrod–Adams
(§7.2), FBNS uses complementarity (§7.3), and COUPLED uses simultaneous
conservation, phase volume and optional energy (the current formulation above).
Gas-free JFO and FBNS impose the JFO free boundary with different numerical
algorithms and different full-film compressibility treatment.

Before the time loop, the pressure field is initialized from the first
non-Neumann axial outlet pressure, not from a supply groove. For
`JFO`, the initialized and boundary film content are derived from that
pressure and the bulk modulus:

$$p^0 = \max(p_{bc},p_{cav}),\qquad
\theta^0 =
\max\!\left[\exp\!\left(\frac{p^0-p_{cav}}{\beta}\right),\theta_{min}\right]$$

Legacy `bc_z_*_theta` config keys are parsed for old files but are not used by
the solver. This keeps boundary and initial conditions pressure-based:
changing pressure or bulk modulus changes the corresponding Elrod film content
through the EOS. Supply inlets are still imposed during the linear solve.

`FBNS` has no bulk modulus, so its film content starts at $\theta^0 = 1$; the
solver detects that as a cold start and builds its own uncavitated initial guess
(§7.3.4).

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

**The full-film test is one function (`FullFilm`, `src/full_film.hpp`).** A cell holds an unbroken liquid film when its film content reaches the ceiling $\theta_{full}$,

$$g = 1 \iff \theta \ge \theta_{full}, \qquad \theta_{full} = 1 - \alpha_{gas},$$

and a face carries Poiseuille flow only when *both* cells it joins do. Everything that gates on rupture asks this one question: the Poiseuille terms of the momentum, energy and dissolved-gas equations, the shear post-processing, the `consistent_boundary_flux` reflood branch, the Elrod rupture-front flip counter, the TYPE_DIFFERENCING face-centrality mask, and the `cavitated_fraction` diagnostic. `FullFilm::bound` is the single definition of $\theta_{full}$.

Currently $\alpha_{gas} = 0$ in that ceiling for every model: free gas is treated as volumetrically inert in the film-content equation, so the liquid may occupy the whole gap and $\theta_{full} = 1$ exactly, which is the historical test. The reason the ceiling is written as one function rather than fifteen literals is that it is not a local decision — moving it to $1 - \alpha_{gas}$ at a subset of the sites would switch off Poiseuille transport across the whole domain as soon as any gas came out of solution, silently and with a converged residual (`docs/FBNS_GAS_COUPLING.md` §3). The Elrod operator carries the same ceiling separately as `VoidParams::theta_full` (§16), which *is* gas-aware under `gas_pressure_coupling = VOID_COUPLED`; the two merge when the resolved-gas bound lands.

**Conservative time derivative (2026-09-06).** Every transient JFO solve adds
$\rho_l^{n+1} h^{n+1}\Delta A/\Delta t$ to the theta diagonal and
$(\rho_l h\theta)^n\Delta A/\Delta t$ to its source. Stored density and gap
include property and squeeze changes in the inventory difference; there is no
additional squeeze source. Steady mode omits storage. See §18.

### 7.3 FBNS Complementarity Solver (`FBNS::solve`, `cavitation_model = FBNS`)

The Elrod–Adams path of §7.2 reaches the JFO conditions indirectly: it introduces
an artificial bulk modulus $\beta$, writes $p = p_{cav} + \beta\ln\theta$, and lets
a switch function $g(\theta)$ turn Poiseuille diffusion off where $\theta < 1$. The
free boundary is then located by reading the sign of a quantity that is
$O(\beta^{-1})$ small near the front, which places the rupture line with
first-order error and a large prefactor (§7.3.5).

`cavitation_model = FBNS` imposes the same JFO conditions **directly**, following
Woloszynski, Podsiadlo & Stachowiak, *Efficient Solution to the Cavitation Problem
in Hydrodynamic Lubrication*, Tribology Letters **58**:18 (2015). There is no bulk
modulus and no switch function anywhere in the model.

#### 7.3.1 The complementarity form

The liquid is incompressible. The unknowns are the pressure $p$ and the film
content $\theta$ of every cell, subject to

$$p - p_{cav} \;\ge\; 0, \qquad \theta_{full} - \theta \;\ge\; 0, \qquad (p - p_{cav})\,(\theta_{full} - \theta) \;=\; 0 ,$$

where $\theta_{full}$ is the largest film content the cell can hold. Without a
resolved gas phase $\theta_{full} = 1$ and this is the classical statement;
§7.3.7 moves it to $1 - \alpha_{gas}$ and nothing else in this section changes.

A cell is either full ($\theta = \theta_{full}$, pressure free to rise) or
cavitated ($p = p_{cav}$, content free to fall) — never both, never neither.
Mass balance is
the ordinary Reynolds equation with the film content multiplying the Couette flux
and the storage term, and the pressure driving Poiseuille flow:

$$\nabla\!\cdot\!\left(\frac{\rho_l h^3}{12\mu}\,\nabla p\right)
\;-\; \frac{\partial}{R\,\partial\vartheta}\!\left(\rho_l\,\frac{\omega R}{2}\,h\,\theta\right)
\;-\; \frac{\partial(\rho_l h \theta)}{\partial t} \;=\; 0 ,$$

with $\vartheta$ the circumferential coordinate. Note that no $g(\theta)$ gates the
Poiseuille term: in a cavitated region complementarity already pins $p = p_{cav}$,
so $\nabla p = 0$ there and the term switches itself off. The reformation jump is
likewise automatic — at the reformation face the upstream cell delivers
$\rho_l\,\omega R\,h_f\,\theta/2$ and the downstream full-film cell absorbs the
difference through its pressure gradient.

#### 7.3.2 Fischer–Burmeister reformulation

The complementarity triple is equivalent to the single smooth equation

$$\mathrm{FB}(a,b) \;=\; a + b - \sqrt{a^2+b^2} \;=\; 0 ,
\qquad a = \frac{p - p_{cav}}{p_{ref}}, \qquad b = \theta_{full} - \theta ,$$

because $\mathrm{FB}(a,b) = 0 \iff a\ge0,\;b\ge0,\;ab=0$. The system is now
*unconstrained and continuously differentiable*, so it admits a Newton method — no
active-set logic, no flag array, no outer flag-update loop.

**Scaling matters.** $\mathrm{FB}$ weights its two arguments equally. Feeding raw
pascals against a film content in $[0,1]$ makes the pressure branch dominate by six
orders of magnitude, drives $\partial \mathrm{FB}/\partial p \to 0$ everywhere, and
degenerates the column swap of §7.3.3. The pressure is therefore nondimensionalized
by $p_{ref}$ (`fbns_pressure_scale`; $0$ = derived automatically as the Couette
pressure scale $6\mu|\omega|R^2/c^2$, raised to cover any supplied overpressure).
The mass balance is likewise divided by a flux scale $g_{ref}$ — the largest
assembled row sum — which only sets what `fbns_tol` means and cannot change the
solution.

#### 7.3.3 Newton with a column swap and a Schur reduction

Write $\alpha = \theta_{full} - \theta$ (the cavity fraction) so both
Fischer–Burmeister
derivatives come out non-negative, and let $G$ be the scaled mass-balance residual
and $F$ the vector of $\mathrm{FB}$ values. Newton's step solves

$$\begin{bmatrix} J_{F,p} & J_{F,\alpha} \\ J_{G,p} & J_{G,\alpha} \end{bmatrix}
\begin{bmatrix}\delta p^{*} \\ \delta\alpha\end{bmatrix} = -\begin{bmatrix}F \\ G\end{bmatrix},$$

where $J_{F,p}$ and $J_{F,\alpha}$ are **diagonal**,

$$J_{F,p}^{(j)} = 1 - \frac{a_j}{\sqrt{a_j^2+b_j^2}}, \qquad
J_{F,\alpha}^{(j)} = 1 - \frac{b_j}{\sqrt{a_j^2+b_j^2}},$$

both in $[0,2]$, and $J_{G,p} = A$ (Poiseuille), $J_{G,\alpha} = B$ (upwind Couette
plus backward-Euler storage). Under `fluid_property_model = CONSTANT` and
`ISOTHERMAL`, $A$ and $B$ are **constant** through the whole Newton iteration:
the mobility, the face fluxes and the storage weight depend only on the film
geometry, and the upwind direction follows the face flux sign, not the solution.
Every nonlinearity lives in $F$.

At the non-differentiable point $a_j = b_j = 0$ the Jacobian arguments — and only
the Jacobian arguments, never $F$ itself — are perturbed to $\pm\varepsilon$
($\varepsilon = 10^{-14}$ on the scaled variables), which gives
$J_{F,p} = J_{F,\alpha} = 1 - 1/\sqrt2$ there.

**Column swap.** For each cell the unknown whose FB derivative is larger is moved
into the first block:

$$r_j = j,\; r_{j+N} = j+N \quad\text{if } J_{F,p}^{(j)} \ge J_{F,\alpha}^{(j)},
\qquad r_j = j+N,\; r_{j+N} = j \quad\text{otherwise},$$

giving $\;A_F = \mathrm{diag}\big(\max(J_{F,p},J_{F,\alpha})\big)$ and
$B_F = \mathrm{diag}\big(\min(J_{F,p},J_{F,\alpha})\big)$. Because
$\max(J_{F,p}^{(j)},J_{F,\alpha}^{(j)}) \in [1-\tfrac{\sqrt2}{2},\,2]$,
$\operatorname{cond}(A_F) \le 4 + 2\sqrt2$ regardless of the state, and $B_F \to 0$
at the solution. $A_G$ and $B_G$ are the correspondingly column-mixed versions of
$A$ and $B$; since $A$ and $B$ share the 5-point stencil, so do $A_G$, $B_G$ and
the Schur complement, which is why the existing `LinearSystem` wrapper carries the
reduced system unchanged.

**Schur reduction.** One sparse solve per Newton step:

$$\big(B_G - A_G A_F^{-1} B_F\big)\,\delta b = -G + A_G A_F^{-1} F,
\qquad A_F\,\delta a = -F - B_F\,\delta b ,$$

the second being a diagonal scaling. $\delta p^{*}$ and $\delta\alpha$ follow by
undoing the swap. At the solution $B_F = 0$, so the reduced matrix is exactly
$B_G$. The linear solve therefore costs what an uncavitated pressure solve costs;
this is the result the paper's title refers to.

**Update.** $p^{*} \leftarrow \max(p^{*} + \omega_N\,\delta p^{*},\,0)$ and
$\alpha \leftarrow \mathrm{clamp}(\alpha + \omega_N\,\delta\alpha,\;0,\;\theta_{full}-\theta_{min})$
with $\omega_N =$ `fbns_relaxation` (default $1$, i.e. a full Newton step). The
truncation is a safeguard, not the model: it snaps each cell exactly onto the
complementary axis, so $\mathrm{FB} = 0$ holds identically from the first iterate
onward and the mass balance $G$ alone drives convergence.

#### 7.3.4 Discretization, boundaries and initial guess

Every face expression matches the finite-volume forms of §5–§7 so the two
mass-conserving solvers can be compared face by face:

| Term | Face expression (east face of cell $P$) |
|---|---|
| Poiseuille | $\mathrm{harm}\big(\gamma_P,\gamma_E\big)\,\dfrac{\Delta z}{R\,\Delta\vartheta}\,(p_P - p_E)$, $\;\gamma = \rho_l h^3/(12\mu)$ |
| Couette | $F_e\,\theta_{up}$, $\;F_e = \bar\rho_l\,\omega R\,\bar h\,\Delta z/2$, first-order upwind on $\mathrm{sign}(F_e)$ |
| Storage | $\big[(\rho_l h\theta)_P^{n+1}-(\rho_l h\theta)_P^n\big]\Delta A/\Delta t$ (complete inventory; §18) |

- **Circumferential**: periodic, handled by the ghost exchange and by
  `LinearSystem`'s wrapping column indices. Nothing special is needed — $\mathrm{FB}$
  is pointwise, so periodicity never touches $A_F$ or $B_F$.
- **Axial `DIRICHLET`**: a half-cell Poiseuille link
  $\gamma_P\,R\Delta\vartheta/(\tfrac12\Delta z)\,(p_P - p_{bc})$ in the mass balance.
  A cavitated boundary cell therefore re-floods at the *physical* Poiseuille rate
  driven by $p_{bc} - p_{cav}$; the Elrod path needs the `consistent_boundary_flux`
  option to approximate this, because its barotropic link overstates re-flooding by
  $\sim\beta/(p_{bc}-p_{cav})$.
- **Axial `NEUMANN`**: no flux, no term. No $\theta$ boundary condition is required
  on either axial end, because the Couette transport is purely circumferential.
- **`INLET_OUTLET`** is rejected: it switches on the sign of the solution and would
  break the differentiability the Newton step relies on.
- **Inlet groove / circular inlet**: both rows of the cell are replaced, which is
  the FBNS boundary treatment — the mass balance becomes $p^{*} = p^{*}_{supply}$
  and the complementarity row becomes $\alpha = 0$ with
  $(J_{F,p}, J_{F,\alpha}) = (0,1)$. Forcing that pair also forces the swap, which
  is what keeps the Schur complement non-singular on those rows. No penalty
  coefficient is used.
- **Initial guess**: on a cold start (no cavitation information in $\theta$ yet) the
  uncavitated problem $G(p,\alpha{=}0) = 0$ is solved once and clamped to
  $p \ge p_{cav}$. A transient step warm-starts from the previous step instead.
  The cold-start test is $\theta \ge 1$ everywhere — a *sentinel* comparison
  against the value `main.cpp` initialises $\theta$ to, deliberately **not** the
  full-film test $\theta \ge \theta_{full}$. Under a moved ceiling (§7.3.7)
  $\theta = 1 - \alpha_{gas}$ is an ordinary converged uncavitated state, so
  testing against $\theta_{full}$ would declare a cold start on every gassy step,
  discard the warm start and refill $\theta$ down to the new ceiling in one go — a
  spurious $\partial\theta/\partial t$ sink that drove the pressure to
  $O(10^3)$ bar on the R290/PZ68 case. The ceiling change belongs to Newton.

#### 7.3.5 Verified behaviour

Against the exact incompressible-JFO solution of the infinitely long grooved
bearing (`cases/00c_jfo_verification/exact_1d_jfo.py`: rupture $302.766608^\circ$,
reformation $79.006725^\circ$, $p_{max} = 7.597141$ MPa, $W' = 495463.6$ N/m),
with the front reported at the face between the last full-film and the first
cavitated cell:

| $N_\vartheta$ | cell | FBNS rupture err | in cells | JFO rupture err | in cells | FBNS $p_{max}$ | JFO $p_{max}$ | FBNS $W'$ | JFO $W'$ |
|---|---|---|---|---|---|---|---|---|---|
| 120 | 3.000° | $-2.767^\circ$ | $-0.92$ | $-23.77^\circ$ | $-7.9$ | $+0.54\%$ | $+50.4\%$ | $+0.70\%$ | $+42.0\%$ |
| 240 | 1.500° | $-1.267^\circ$ | $-0.84$ | $-17.77^\circ$ | $-11.8$ | $+0.28\%$ | $+39.9\%$ | $+0.27\%$ | $+41.0\%$ |
| 480 | 0.750° | $-0.517^\circ$ | $-0.69$ | $-10.27^\circ$ | $-13.7$ | $+0.14\%$ | $+24.0\%$ | $+0.07\%$ | $+31.1\%$ |
| 960 | 0.375° | $-0.142^\circ$ | $-0.38$ | $-5.02^\circ$ | $-13.4$ | $+0.07\%$ | $+12.2\%$ | $-0.04\%$ | $+17.8\%$ |

The reformation front lands within $0.66$ cells at every resolution
($+1.99^\circ, +0.49^\circ, +0.49^\circ, +0.12^\circ$). The JFO column keeps a
front error of $7$–$13$ cell widths — it shrinks in degrees but not relative to the
mesh, which is the defect that motivated this solver.

Convergence is Newton-like and grid-independent: **4 Newton iterations** at every
resolution from $N_\vartheta = 120$ to $960$. The residual history at
$N_\vartheta = 480$ is

$$6.40\times10^{-4} \to 6.02\times10^{-4} \to 3.93\times10^{-5} \to 1.23\times10^{-5} \to 3.49\times10^{-12},$$

the last step dropping seven orders at once: once the active set (which cells are
cavitated) stops changing, the remaining system is exactly linear and the semismooth
Newton step terminates finitely.

The per-cell continuity residual (§15.1) on the converged $N_\vartheta = 480$ run is
$1.5\times10^{-11}$ kg/s at worst, $3.5\times10^{-10}$ of the peak Couette face flux
— the Krylov tolerance — and $3\times10^{-14}$ of it on the rupture and reformation
cells themselves. The groove cells are exactly zero (their imbalance is the modelled
supply, reported separately as `inlet_mass_source`).

#### 7.3.6 Thermohydrodynamic coupling

The energy equation and the single-phase liquid property layer reach the FBNS
operator only through its per-cell coefficients — the Poiseuille mobility

$$\Gamma_{i,j} = \frac{\rho_{l,i,j}\, h_{i,j}^{3}}{12\,\mu_{i,j}},$$

which `FBNS::solve` rebuilds at the top of every call from the current `mu` and
`rho_liquid_solution` fields. They are then **frozen for the duration of that
Newton iteration**. This is deliberate, and it is the same split the Elrod path
uses:

- The Jacobian blocks $A$ and $B$ are differentiated from exactly the residual
  `FBNS::cell_outflow` evaluates, because both read the same frozen $\Gamma$. No
  $\partial\mu/\partial T$ or $\partial\mu/\partial p$ term appears in either. A
  Jacobian carrying those terms while the residual did not is what would make
  Newton limp; freezing both keeps the local quadratic rate and the finite
  termination of §7.3.5 intact.
- The $\mu(T,p)$ feedback is closed one level up, by the segregated outer (Picard)
  loop in `main.cpp`: properties $\to$ FBNS $\to$ properties $\to$ gas $\to$
  velocities $\to$ energy, repeated until the normalised pressure and temperature
  changes both fall below `coupling_tolerance`. FBNS warm-starts from the previous
  outer iterate, so the later Newton solves cost 0–3 iterations.

Consequences for the energy equation itself:

- **Liquid-fraction weighting.** `theta` carries identical semantics under both
  mass-conserving models, so `cavitated_film_model = STRIATED` weights heat
  generation, stored energy, conduction and convected enthalpy by
  $\theta \in [\theta_{min},1]$ for FBNS through exactly the same code path (§13.1).
  Full-film FBNS cells hold $\theta = 1.0$ *exactly* — the Newton update clamps the
  cavity fraction at zero — so the `theta >= 1` full-film test needs no tolerance.
- **Poiseuille faces are never switched off under FBNS.** Elrod carries a
  barotropic $p = p_{cav} + \beta\ln\theta$ through the ruptured zone, which is not
  the physical pressure there, so its enthalpy flux is gated by the full-film
  switch. FBNS pins a ruptured cell at $p = p_{cav}$ by complementarity, so every
  gradient it produces is physical: between two ruptured cells it vanishes on its
  own, and at a rupture front or an open axial end it is precisely the reformation
  flux the FBNS mass balance conserved. Gating it would leave the energy equation
  convecting enthalpy on a different velocity field from the one carrying mass —
  the defect `consistent_boundary_flux` exists to patch on the Elrod path, and
  which FBNS gets for free. The full-film switch is *kept* on the Poiseuille
  dissipation term, which is a striated-film modelling claim rather than a
  conservation requirement and can only touch the single front cells.

**Verified** (`tests/test_fbns.cpp`, tests 8–10, effectively-1D grid
$N_\vartheta = 720$, $N_z = 4$, Neumann ends, Andrade $E_\mu = 3700$ K):

| Check | Result |
|---|---|
| SINGLE_PHASE with all property coefficients zero vs CONSTANT | bit-identical $p$, $\theta$ |
| converged coupled $\mu$ field replayed through an isothermal FBNS solve | $\lvert\Delta p\rvert/p_{ref} = 1.1\times10^{-9}$, $\lvert\Delta\theta\rvert = 5.8\times10^{-9}$, against a $\mu$ that moved 41% |
| FBNS vs Elrod temperature, shared full-film arc, $\varepsilon = 0.3$ | $\max\lvert\Delta T\rvert = 0.44$ K over 465 cells |
| Picard iterations to `coupling_tolerance` $=10^{-8}$ | 10 (FBNS), 9 (Elrod) |
| load capacity, thermal vs isothermal | $1.653\times10^{4}$ N vs $1.807\times10^{4}$ N |

The **fronts** are deliberately not compared against Elrod. Elrod's rupture front
sits 12–18° from the exact 1-D solution on this bearing (§7.3.5) and that gap is a
property of its barotropic smearing, not of the grid: it stays near 13° from
$N_\vartheta = 240$ to $720$, and raising $\beta$ to sharpen the switch makes it
worse *and* stops the Elrod outer loop converging (200+ iterations at
$\beta = 10^{10}$ against FBNS's 13). The exact benchmark of §7.3.5 is the only
front reference that means anything. The residual $\Delta T$ above is Elrod's
misplaced front leaking into its pressure field and hence into the Poiseuille
share of the enthalpy transport, which is ~25% of the total here; it scales with
the thermal rise (0.24 K at $\varepsilon = 0.2$, 0.44 K at 0.3, 0.70 K at 0.4,
1.42 K at 0.6) and is not a property of the coupling.

**Known instability, both cavitation models.** A groove that pins pressure but is
left thermally transparent (`feeds_fresh_oil = false`) diverges under a strong
Andrade slope. The Poiseuille dissipation $h^{3}\lvert\nabla p\rvert^{2}/(12\mu)$
scales as $1/\mu$, and at the groove edge $\nabla p$ is held fixed by the pressure
pin, so hotter $\to$ thinner $\to$ hotter is an unbounded positive feedback: on the
`grid_G3_fbns` geometry with $E_\mu = 3700$ K it reaches $T \sim 10^{6}$ K under
FBNS *and* under Elrod (JFO), at every under-relaxation from 1.0 down to 0.15.
Giving the groove its supply temperature (`inlet_groove = 90 18 1.713e5 313.15`)
removes the artefact — a fed groove is a reservoir, not a sheared film — and the
same case then converges in 9 Picard iterations. This is a pre-existing property
of the energy model, unrelated to which cavitation solver is selected.

#### 7.3.7 Gaseous cavitation: the moved full-film ceiling

Under `fluid_property_model = TWO_PHASE` a resolved free-gas phase occupies a
volume fraction $\alpha_{gas}$ of the gap, so the liquid can only fill what is
left and the ceiling of §7.3.1 becomes

$$\theta_{full} = 1 - \alpha_{gas}, \qquad b = (1 - \alpha_{gas}) - \theta \ \ge\ 0 .$$

**Why this composes.** $\alpha_{gas}$ is *lagged*: the gas path is operator-split
and runs after the flow solve (`GasTransport::solve` then
`FluidProperties::update_gas_state` in the outer loop of `main.cpp`), so within one
Newton solve it is a coefficient, not an unknown. The new $b$ is therefore affine
in $\theta$ with the same derivative as before,

$$\frac{\partial b}{\partial \theta} = -1 ,$$

and consequently:

- $J_{G,\alpha} = B$ is **unchanged** — the mass balance only ever used
  $\partial(\text{outflow})/\partial\alpha = -\partial(\text{outflow})/\partial\theta$;
- $J_{F,p}$ and $J_{F,\alpha}$ keep their Fischer–Burmeister forms, so the column
  swap, the bound $\operatorname{cond}(A_F) \le 4 + 2\sqrt2$ and the Schur
  complement of §7.3.3 are untouched;
- the Poiseuille mobility needs nothing new: the free gas already reaches it
  through the mixture viscosity written into the `mu` field (§7.3.6, §13.2.1).

**One definition of the ceiling.** The literal test `theta >= 1.0` used to appear
at fifteen sites in `reynolds.cpp`, `energy.cpp`, `gas_transport.cpp` and
`diagnostics.cpp`. Moving a ceiling in fifteen places moves it wrong: the instant
any gas came out of solution every one of those switches would turn off at once
and the solver would report a fully cavitated bearing carrying no Poiseuille flow,
no Poiseuille enthalpy and no dissolved-gas transport — silently, with a converged
Newton residual. All fifteen therefore read `FullFilm::bound`
(`src/full_film.hpp`), which is the only place $\theta_{full}$ exists. It returns
exactly $1$ for every model except FBNS + `TWO_PHASE`, which is what makes the
Elrod and Gumbel paths byte-identical by construction rather than by inspection.

**Outer-loop damping.** The moved ceiling closes a feedback loop
$\alpha_{gas} \to \theta_{full} \to \theta \to p \to \alpha_{gas}$ that the
volumetrically-inert gas model never had. It is the same loop
`gas_pressure_coupling = VOID_COUPLED` opens on the Elrod path, and it takes the
same lever: `FluidProperties::update_gas_state` under-relaxes the free-gas state by
`void_coupling_relaxation` (default $0.3$) whenever the gas displaces liquid. The
fixed point is unchanged; only the path to it is damped.

**Verified** (`tests/test_fbns.cpp`, tests 11–13, 00c bearing at
$N_\vartheta = 240$, saturated propane feed, five transient steps from the
converged steady film):

| Check | Result |
|---|---|
| $\max\big(-(p-p_{cav})/p_{ref}\big)$ | $0$ exactly |
| $\max\big(\theta - (1-\alpha_{gas})\big)$ | $0$ exactly |
| $\max\big\lvert (p-p_{cav})\,((1-\alpha_{gas})-\theta)\big\rvert$ | $0$ exactly |
| per-cell $\lvert r\rvert$ / peak Couette flux, away from the groove | $5.1\times10^{-15}$ |
| free gas released, FBNS vs JFO | $\alpha_{gas}^{max} = 2.50\times10^{-2}$ vs $8.44\times10^{-2}$ |
| Newton iterations per step | $\le 7$ |
| $\alpha_{gas} \equiv 0$ vs `SINGLE_PHASE`, three steps | bit-identical $p$, $\theta$ |
| `FullFilm::bound` over 24 model combinations | exactly $1$ except FBNS + `TWO_PHASE` |

The **fronts are logged and asserted against nothing**, for the reason given in
§7.3.6: Elrod's rupture front sits 12–18° from the exact 1-D solution on this
bearing, so agreeing with it would be evidence of a defect. The complementarity
invariants above are what the solver owns. Note also that the *written output*
pairs $\theta$ with an $\alpha_{gas}$ that the gas update has already advanced, so
$\theta \le 1-\alpha_{gas}$ need not hold in a VTK file; it holds at the point the
FBNS solve returns, which is where the constraint lives.

**Historical stiffness result, prior to the §18 conservation repair.** These
transient results used a mass-inconsistent relaxation and must be rerun. The
relaxation setting now applies only to steady feed flash. On
`config_r290_pz68_quick.txt` the feed entered exactly at
its bubble point and the flash releases enough propane to fill $\sim$23% of the gap
in one $5\,\mu$s step. The ceiling drops by that much at once, the resulting
$\partial\theta/\partial t$ sink swings the load between $5.8\times10^{2}$ N and
$2.4\times10^{5}$ N on alternate steps, and the amplitude decays only $\sim$30% per
step pair. The swing is **timestep-independent** — identical at
$\Delta t = 5\times10^{-6}$, $5\times10^{-8}$ and $5\times10^{-9}$ s — because the
released fraction is itself proportional to $\Delta t$, so this is the model's
displacement rate rather than a temporal-resolution artefact. Every Newton solve
still converges (61/61 steps, 2–11 iterations, mean 6.3) and the gas mass stays
bounded, but the load history over the first tens of steps is transient noise, not
a converged operating point. `void_coupling_relaxation` is the lever; the Elrod
path with `VOID_COUPLED` does not resolve this case either (it drifts monotonically
instead of oscillating), so this is a property of a volumetrically-active gas model
on a saturated feed, not of FBNS.

#### 7.3.8 Scope and limits

Supported: `temperature_model` `ISOTHERMAL` and `ENERGY_EQUATION`,
`fluid_property_model` `CONSTANT`, `SINGLE_PHASE` and `TWO_PHASE` (all
`density_model` / `liquid_viscosity_model` options), `solution_mode`
`STEADY_STATE` and `TRANSIENT` (backward Euler), `motion_model` `STATIC` and
`MOVING_BEARING`, groove and circular inlets, `DIRICHLET`/`NEUMANN` axial
boundaries, `cavitated_film_model` `FULL_FILM` and `STRIATED`,
`cross_film_thermal_model` `POLYNOMIAL` (§7.3.9), and MPI decomposition in
$\vartheta$.

**Rejected at config validation** rather than silently ignored, because each
changes the complementarity datum or the operator's stencil rather than its
coefficients:

- `gas_pressure_coupling = VOID_COUPLED`. It replaces the scalar cavitation datum
  $p_{cav}$ with a per-cell plateau field `p_void`, and that field is advanced and
  under-relaxed **only** inside `Reynolds::solve_elrod`'s outer flag-update loop
  (`reynolds.cpp`). FBNS has no such loop, so it would run the plateau frozen at
  its initial value $p_{cav}$ — silently reducing to `NONE` while reporting
  `VOID_COUPLED`. Making it work needs a plateau update site of its own plus its
  own validation; see `docs/FBNS_GAS_COUPLING.md` §5.
- `INLET_OUTLET` axial boundaries (they switch on the sign of the solution and
  would break the differentiability the Newton step relies on).

**$p_{cav}$ is a scalar on this path, and correctly so.** Under
`cavitation_threshold = LOCAL_PSAT` the datum is
$\max\big(p_{cav},\,p_{sat}(c_{d,init},T_{init})\big)$ — evaluated at the *feed*
state, hence uniform over the domain by construction, so that it acts as a
pressure datum without injecting a net force (`SimulationConfig::effective_p_cav`).
The Fischer–Burmeister argument $a = (p - p_{cav})/p_{ref}$ therefore already reads
the correct local value and needs no per-cell array. The only genuinely per-cell
cavitation pressure in the codebase is `p_void` under `VOID_COUPLED`, which is
rejected above.

`theta_convection_scheme` is ignored (FBNS uses first-order upwind, and a
higher-order upwind stencil would add $\vartheta$-second-neighbour columns to $B$);
`bulk_modulus` is unused.

#### 7.3.9 Generalized (Dowson) coefficients under FBNS

With `cross_film_thermal_model = POLYNOMIAL` the FBNS mass balance is assembled
from the generalized Reynolds equation of §13.1.3 rather than the classical one:

$$\nabla\!\cdot\!\left(\rho_l F_2\,\nabla p\right)
\;-\; \frac{\partial}{R\,\partial\vartheta}\!\left(\rho_l\,\frac{\omega R}{2}\,\frac{2F_1}{F_0}\,\theta\right)
\;-\; \frac{\partial(\rho_l h \theta)}{\partial t} \;=\; 0 .$$

Only the two *flux* coefficients move: $\rho_l h^3/12\mu \to \rho_l F_2$ in the
Poiseuille term and $h \to 2F_1/F_0$ in the Couette term. The storage term keeps
the geometric $h$, because it measures a volume of liquid, not a flux. This is
what puts $\mu(T(y))$ **inside** the film integrals instead of evaluating
$h^3/12\mu$ at a depth-averaged viscosity; §13.1.3 explains why the two are not
the same thing and why the effective viscosity implied by $F_2$ is the lower one.

**The Newton structure is untouched.** `CrossFilm::update` rebuilds $F_0,F_1,F_2$
once per outer (Picard) iteration, *before* `FBNS::solve` is called, so inside one
Newton solve they are frozen coefficients in exactly the sense $\mu$ already was
(§7.3.3): $A$ and $B$ are still built once at the top of the solve, the mass
balance is still linear in $(p,\alpha)$, every nonlinearity still lives in $F$,
and the residual and the Jacobian read the same frozen mobilities — which is the
condition for the quadratic local rate. A Jacobian carrying
$\partial F_2/\partial T$ terms against a residual that did not would be the
limping case, and is not what is done.

The complementarity datum is unaffected: $F_2 > 0$ everywhere, so a cavitated cell
still pins itself at $p = p_{cav}$ and its Poiseuille flux still vanishes through
$\nabla p = 0$ rather than through a switch function.

`Reynolds::solve_elrod`, `FBNS::solve` and `Diagnostics` read these coefficients
from one definition (`src/reynolds_coefficients.hpp`), so the residual §15.1
reports is the residual of the operator that was actually assembled — see the note
at the end of §15.1.

### 7.4 Outer (Picard) coupling loop (`src/coupling.hpp`, `main.cpp`)

In `solution_mode = STEADY_STATE` the segregated sweep

$$\text{properties} \to \text{Reynolds/FBNS} \to \text{properties} \to \text{gas}
\to \text{velocities} \to \text{energy}$$

is repeated until the normalised max changes

$$R_p = \frac{\max_{i,j}|p^{k}-p^{k-1}|}{p_{ref}},\qquad
R_T = \frac{\max_{i,j}|T^{k}-T^{k-1}|}{T_{ref}}$$

both fall below `coupling_tolerance`, with $p_{ref} =$
`property_reference_pressure` and $T_{ref} =$ `temperature_reference`. TRANSIENT
runs exactly one sweep per timestep. Iterates are ranked by the single scalar
$\rho^{k} = \max(R_p^{k}, R_T^{k})$, which is exactly the quantity that has to
fall below the tolerance.

**Under-relaxation.** After the residual of sweep $k$ is formed,

$$p \leftarrow \omega\, p^{k} + (1-\omega)\, p^{k-1},\qquad
T \leftarrow \omega\, T^{k} + (1-\omega)\, T^{k-1}$$

with $\omega =$ `coupling_relaxation` $\in (0,1]$. Two schedules:

* `coupling_relaxation_mode = FIXED` (default) holds $\omega$ constant. With
  $\omega = 1$ this is a no-op and the loop is the plain segregated sweep.
* `coupling_relaxation_mode = ADAPTIVE` starts at $\omega =$
  `coupling_relaxation` and, after each sweep,
  $\omega \leftarrow \max(\omega_{min}, \beta\omega)$ if $\rho^{k} > \rho^{k-1}$,
  else $\omega \leftarrow \min(\omega_{max}, \gamma\omega)$, with
  $\beta =$ `coupling_relaxation_backoff` (0.5), $\gamma =$
  `coupling_relaxation_growth` (1.1) and $\omega_{min} =$
  `coupling_relaxation_min` (0.1). It targets the signature of a limit cycle —
  loop gain slightly above 1 with a reversing residual trend — without needing a
  per-case $\omega$. Measured performance on the Ferron conjugate case: §13.1.2.

Note the pre-existing limitation (HANDOFF **D8**): the relaxation acts on $p$ and
$T$, while under Elrod the primary state is $\theta$. Under FBNS the relaxed $p$
enters the next sweep only through the property update and the dissipation term,
because FBNS re-solves $(p,\theta)$ from scratch each sweep.

**Best-iterate rollback.** A Picard loop is not a descent method: nothing forces
$\rho^{k}$ to be monotone, and a loop that exhausts `coupling_max_iters` can end
on an iterate that is worse than one it already passed through. That is exactly
what the conjugate Ferron cases do (§13.1.2), where the delivered iterate 200 was
6.3x worse than iterate 56. The loop therefore snapshots the whole field set at
every new minimum of $\rho$ and, **on cap exhaustion without meeting the
tolerance**, restores that snapshot. Three properties make this safe:

1. A loop that *meets* the tolerance always delivers its last iterate, so every
   converging run is bit-identical to the pre-rollback solver.
2. The snapshot is taken *before* under-relaxation, so the delivered state is a
   genuine end-of-sweep solve rather than a blend of two sweeps.
3. The snapshot covers every registered field including ghost layers and
   `old_data`, so $p$, $\theta$, $T$, the wall state, the properties and the
   integrated forces are restored as one self-consistent iterate.

`coupling_keep_best_iterate = false` restores the historical last-iterate
behaviour.

**Non-convergence is reported, not implied.** A steady run that ends without
meeting `coupling_tolerance` prints `NOT CONVERGED`, the cap it hit, and the
best-vs-final residual comparison, in red, regardless of `log_outer_iters`;
`results/diagnostics.csv` carries `coupling_iters`, `coupling_R_p`,
`coupling_R_T`, `coupling_converged`, `coupling_best_iter`,
`coupling_best_residual`, `coupling_restored_best` and `coupling_relaxation` for
every step. `log_coupling_history = true` additionally prints $R_p$, $R_T$ and
$\omega$ for every sweep, which is how the traces in §13.1.2 were produced.

**What the rollback does not do.** It returns the best iterate the loop visited;
it does not make that iterate converged. On the Ferron conjugate case the best
iterate still sits at $\rho = 4.3\times10^{-5}$, 43x above the tolerance. The
rollback is an honesty and robustness measure; under-relaxation is the fix.

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

The older liquid-temperature path described here supports temperature-dependent
viscosity and optional polynomial cross-film/conjugate-wall closures. COUPLED
has a separate mixture-energy formulation described at the start of this file.
Elastic deformation and calculated thermal expansion of the gap remain planned.

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

Moving-bearing JFO and FBNS now difference the complete liquid inventory
$(\rho_l h\theta)^{n+1}-(\rho_l h\theta)^n$ (§18). This includes density as
well as geometry changes; the old product-rule approximation omitted density
storage and has been removed from these paths.

### 13.1.1 Groove Hot-Oil Carry-Over

A fed groove does not deliver pure supply oil. The film sweeping into it
recirculates and mixes with the fresh charge, so the oil leaving the groove is
warmer than `t_supply`. `groove_mixing_model` selects the closure:

$$T_{groove}(z) = \lambda\,T_{supply} + (1-\lambda)\,T_{upstream}(z)$$

- `PINNED` (default) is $\lambda = 1$: the historical hard pin, bit-identical to
  the model-free solver.
- `FIXED_LAMBDA` takes $\lambda$ from the inlet's own `mixing_lambda` token.
- `FLOW_WEIGHTED` computes it from the two flow rates,

$$\lambda = \frac{Q_{supply}}{Q_{supply} + Q_{up}}, \qquad
  Q_{up} = \sum_j \tfrac{1}{2}|\omega| R\,h\,\theta_{film}\,\Delta z ,$$

  with $Q_{supply}$ the inlet's configured `supply_flow_rate` and $Q_{up}$ the
  solved circumferential Couette carry-over crossing the groove's leading edge.

$T_{upstream}(z)$ is the flow-weighted temperature of the cell immediately
upstream of the groove's leading edge, taken per axial row so the axial structure
of the incoming film survives the mixing. The leading edge is defined in the
sense the shaft entrains ($\mathrm{sign}\,\omega$), and the sample cell may be
owned by any rank, so both the temperature and the flux are gathered with a sum
reduction over an array indexed by the (undecomposed) axial index. The mixed
temperature is imposed on groove cells by replacing their linear-system rows
with prescribed values at the existing diagonal scale.

### 13.1.2 Cross-Film Temperature Profile and Conjugate Walls

The depth-averaged energy equation of §13.1 represents one temperature across
the gap. Its adequacy depends on cross-film gradients from wall temperatures,
viscous dissipation and convection. Brinkman number and reduced Peclet number
$\rho c_p U h^2/(kR)$ help assess these effects; a converged depth-average
solution alone cannot establish that those gradients are negligible.
`cross_film_thermal_model = POLYNOMIAL` replaces the assumption with a
reconstructed profile. `NONE` (the default) removes the entire model and is
bit-identical to the previous solver.

Let $\eta = y/h \in [0,1]$ with $y = 0$ at the shaft (journal) surface, which
moves at $U = \omega R$, and $y = h$ at the bush surface. The profile is a
quartic in three zero-mean shape functions,

$$T(\eta) = \bar T + a\,\psi_1(\eta) + b\,\psi_2(\eta) + c\,\psi_4(\eta),$$

$$\psi_1 = \eta - \tfrac12, \qquad
  \psi_2 = \eta - \eta^2 - \tfrac16, \qquad
  \psi_4 = \eta^2(1-\eta)^2 - \tfrac1{30},$$

each of which integrates to zero over $[0,1]$. The depth average of $T(\eta)$ is
therefore $\bar T$ **by construction**: the profile is collocated with the
depth-averaged balance of §13.1 and can never contradict it.

Three conditions close the three coefficients.

**1. Wall values.** $T(0) = T_s$ (shaft) and $T(1) = T_b$ (bush) give

$$a = T_b - T_s, \qquad b + \tfrac15 c = 6\,\Delta T_w,
  \qquad \Delta T_w = \bar T - \tfrac12 (T_s + T_b).$$

**2. Midplane collocation of the cross-film energy balance** (order 4 only;
`cross_film_profile_order = 2` sets $c = 0$):

$$k\,\frac{\partial^2 T}{\partial y^2}
  = \rho c_p\!\left(u\frac{\partial T}{\partial s} + w\frac{\partial T}{\partial z}\right)
    - \Phi(y).$$

With $\psi_2'' = -2$ and $\psi_4''(\tfrac12) = -1$ this is $-2b - c = X$ where
$X = (h^2/k)\,\mathrm{rhs}(\tfrac12)$, so

$$c = -\tfrac53 X - 20\,\Delta T_w, \qquad
  b = \tfrac13 X + 10\,\Delta T_w.$$

In the pure conduction-plus-dissipation equilibrium ($\Delta T_w = S/6$,
$X = -2S$, $S = \mu U^2/2k$) this returns $c = 0$, $b = S$: the classical
parabola $T = T_w + S\,\eta(1-\eta)$ is recovered **exactly**.

**3. Conjugate walls** supply $T_s$ and $T_b$ (below).

**Two stated approximations in the collocation.** The right-hand side is
evaluated with the Couette part of the flow only ($\tau = -U/F_0$, uniform across
the gap), and the convective term is saturated against the dissipative one,

$$X = X_\Phi + X_{conv}\frac{|X_\Phi|}{|X_\Phi| + |X_{conv}|},
  \qquad X_\Phi = -\frac{h^2}{k}\Phi(\tfrac12),
  \qquad X_{conv} = \frac{h^2}{k}\rho c_p\,u\frac{\partial\bar T}{\partial s}.$$

Both are exact in the conduction-dominated limit the profile is derived for. The
saturation exists because one midplane collocation cannot represent the thin wall
layers that appear at reduced Peclet $\gg 1$; unbounded, the convective term
produced several-hundred-kelvin cores on the Ferron rig. The Couette-only
evaluation exists because a $\partial p/\partial s$ dependence makes the profile
shape jump whenever an Elrod cell flips between full film and cavitated, which
turns the segregated solve into a limit cycle (measured: the coupling residual
$R_p$ sat at 9.4 bar, identical at 40 and at 120 outer iterations, and
under-relaxation did not move it). Both are recorded with their measured evidence
in `src/cross_film.cpp`. The dissipation source, the $F$-integrals, the wall
fluxes and the Reynolds coefficients all still use the full pressure field.

**Conjugate walls.** The wall conduction fluxes follow from the profile,

$$q_{shaft} = \frac{k}{h}(a + b), \qquad q_{bush} = \frac{k}{h}(b - a),$$

positive meaning heat leaving the film. (The optional `heat_flux_journal` /
`heat_flux_bearing` output fields of §15.2 are these with the sign flipped, so
that a field named after a surface is positive *into* the fluid.) The bush closes locally through 1-D
radial conduction in series with outer-surface convection, referred to the bore
area,

$$q_{bush} = G\,(T_b - T_{b,out}), \qquad
  \frac1G = \frac{R\ln(r_o/R)}{k_{bush}} + \frac{R}{r_o\,h_{out}},$$

where $T_{b,out}$ is `bush_outer_temperature` when positive and otherwise
`thermal_ambient_temperature`. This is linear in $T_b$ and solves in closed form.
`bush_outer_heat_transfer = 0`
makes $G = 0$, an adiabatic bush; `bush_outer_boundary = FIXED_TEMPERATURE` is the
$h_{out}\to\infty$ limit of the same series resistance, $G = k_{bush}/(R\ln(r_o/R))$.
This is `bush_thermal_model = LOCAL`, the default. **Every watt that enters the
bore leaves through the wall directly underneath it**, so the bush surface
temperature can only peak where the film is hottest — on the Ferron rig, at the
groove. §13.1.2.1 replaces it. The shaft is treated as **circumferentially and
axially isothermal**, its single temperature $T_s$ set by the global film-to-shaft
balance. This reduced closure does not reproduce Ferron's axisymmetric $r$-$z$
shaft conduction or its prescribed axial-end temperature.

$$\oint q_{shaft}\,\mathrm{d}A
  = h_{shaft}\,A_{shaft}\,(T_s - T_{s,amb}),$$

where $T_{s,amb}$ is `shaft_ambient_temperature` when positive and otherwise
`thermal_ambient_temperature`.

With `shaft_ambient_heat_transfer = 0` the shaft is globally adiabatic: the
integral of the film-to-shaft flux is zero, so heat entering the shaft in the
loaded arc leaves it again upstream. **That is the transport path the
constant-$T_{wall}$ Newton closure has no representation of**, and it is why no
value of `journal_heat_transfer` could move the predicted temperature peak.
Substituting the bush relation $T_b = \alpha + \beta T_s$ into the shaft balance
closes $T_s$ in one step; only four scalars need reducing across ranks.

The net wall term entering the depth-averaged equation is
$Q_w = q_{shaft} + q_{bush} = 2(k/h)\,b$, applied with the consistent Newton
linearisation

$$Q_w(T) \approx Q_w(\bar T_{lag})
  + \frac{\partial Q_w}{\partial \bar T}\,(T - \bar T_{lag}),
  \qquad \frac{\partial Q_w}{\partial \bar T}
  = 2\frac{k}{h}\,b_{slope}\left(1 - \frac{b_{slope}}{2}\frac{\partial T_b}{\partial\bar T}\right),$$

with $b_{slope} = 10$ (order 4) or $6$ (order 2). Taking $T_b$ as frozen would
give $12k/h$, the fully-sunk-wall value, and over-damp the adiabatic-bush limit
by a factor of four. When this closure is active, `journal_heat_transfer`,
`bearing_heat_transfer`, `journal_wall_temperature` and
`bearing_wall_temperature` are ignored, and the config validator emits a warning
saying so.

#### 13.1.2.1 Conjugate Bush: 2-D (r, theta) Conduction

`bush_thermal_model = CONDUCTION_2D` replaces the local closure by steady heat
conduction in an axially collapsed bush annulus. This is a reduction of the
three-dimensional bush equation used by Ferron et al., not its full spatial
model. It supplies a transport path absent from the local closure:
**circumferential redistribution inside the bush**. Heat picked up under the
loaded arc can leave the bush at a different angular position, which is what
turns the measured bush temperature down before the groove. With the local
closure the bush surface necessarily peaks where the film is hottest — near the
groove, at about $348^\circ$ on the Ferron rig, against a measured $258^\circ$.

**Governing equation.** In the annulus $R \le r \le r_o$, $\theta$ periodic,

$$\frac1r\frac{\partial}{\partial r}\!\left(r\,k_{bush}\frac{\partial T}{\partial r}\right)
  + \frac{1}{r^2}\frac{\partial}{\partial\theta}\!\left(k_{bush}\frac{\partial T}{\partial\theta}\right) = 0 .$$

The implemented bush is **$z$-collapsed** with no axial dependence, while
Ferron's equation includes the axial second derivative. Axial end-face cooling
is also absent here. The paper's $h_b=80\ \mathrm{W/(m^2K)}$ applies to those
end faces, not the outer cylinder, whose temperature is prescribed. No universal
percentage bound on the neglected axial heat loss is established by this code
(Ferron et al., *[J. Lubr. Technol.](https://doi.org/10.1115/1.3254632)*, 1983).

**Inner boundary: flux continuity with the film.** The cross-film profile makes
the film-side flux *linear* in the bore temperature. Writing $b$ from §13.1.2 and
eliminating $T_m = (T_s + T_b)/2$,

$$q_{bush} = A_0 + A_1 T_s - B\,T_b, \qquad
  A_0 = \frac{k}{h}\!\left(b_{const} + b_{slope}\bar T\right), \quad
  A_1 = \frac{k}{h}\!\left(1 - \tfrac{b_{slope}}{2}\right), \quad
  B = \frac{k}{h}\!\left(1 + \tfrac{b_{slope}}{2}\right) > 0,$$

which is the same relation the local closure inverts against $q = G(T_b-T_{b,out})$.
Here it is imposed as a **Robin condition on the bore**, so the bush problem stays
linear and needs no iteration. Because the bush is $z$-collapsed, the coupling
uses the **$z$-average** of $A_0$, $A_1$ and $B$ at each $\theta$, and returns one
$T_b(\theta)$ that the cross-film profile then uses at every $z$ of that column.

**Outer boundary.** `bush_outer_boundary = CONVECTION` gives
$-k_{bush}\partial T/\partial r = h_{out}(T - T_{b,out})$ at $r = r_o$
(`bush_outer_heat_transfer`; $T_{b,out}$ is `bush_outer_temperature` when positive
and otherwise `thermal_ambient_temperature`; $h_{out} = 0$ is an adiabatic outer
surface, leaving the bush a pure circumferential redistributor).
`FIXED_TEMPERATURE` holds $T(r_o) = T_{b,out}$, the $h_{out}\to\infty$ limit of the
same expression.

**Discretisation.** Finite volume, uniform in $\theta$ on the solver's own theta
grid, uniform in $r$ with `bush_radial_cells` (default 16) cells. Every radial
face uses the **exact logarithmic shell conductance**
$G_r = k_{bush}\,\Delta\theta/\ln(r_2/r_1)$, including the two half cells at the
boundaries, so a purely radial problem is reproduced to round-off at any radial
resolution (asserted in `tests/test_bush_conduction.cpp`). Circumferential faces
use $G_\theta = k_{bush}\,\Delta r/(r_c\Delta\theta)$. The bore condition is put in
series with the bore-to-first-node half cell, which eliminates $T_b$ from the
linear system without losing accuracy:

$$Q_{\to\text{cell }0} = \lambda\,S\,(A - B\,T_0), \qquad
  \lambda = \frac{C_{1/2}}{C_{1/2} + S B}, \qquad
  T_b = T_0 + \frac{Q}{C_{1/2}},$$

with $S = R\Delta\theta$ the bore face length and $C_{1/2}$ the half-cell
conductance. Setting $B$ from the film and $A$ to zero recovers a pure Neumann
bore; setting $\Delta\theta$-uniform data recovers the local closure exactly.

**Linear solve.** Direct, not iterative. The system is block tridiagonal in
$\theta$ with $n_r \times n_r$ blocks plus one periodic wrap, so blocks
$0\ldots n_\theta-2$ are swept with block Thomas carrying $n_r$ extra right-hand
side columns ($X_i = Y_i + Z_i X_{n_\theta-1}$), and one dense $n_r \times n_r$
Schur solve closes the wrap. Cost is $O(n_\theta n_r^3)$, about $10^7$ flops at
$720\times16$ — negligible beside one film pressure solve, with no tolerance to
tune. Line relaxation was tried first and rejected: the bush grid is close to
isotropic, which is exactly where it stalls.

**Coupling to the shaft.** Because $q_{bush}$ is affine in $T_s$, the annulus is
solved for *both* right-hand sides on one factorisation, giving
$T_b(\theta) = \alpha(\theta) + \beta(\theta)\,T_s$ exactly. The global shaft
balance then closes $T_s$ in one step, unchanged from the local closure. The
solve is serial and deterministic: the three coupling arrays are all-reduced onto
the global $\theta$ grid and every rank runs the same arithmetic on the same
input, which removes the broadcast and any rank-0 asymmetry.

**Newton linearisation.** Under `CONDUCTION_2D` the exact
$\partial T_b/\partial\bar T$ is nonlocal — one film cell barely moves a bush its
neighbours also feed — so the bush is taken as **frozen**,
$\partial T_b/\partial\bar T = 0$, giving the fully-sunk-wall coefficient
$2(k/h)\,b_{slope}$. Two alternatives were implemented and rejected:

* the bore admittance with the neighbouring bush cells frozen reproduces the
  frozen-bush number to 5 % (measured $\partial T_b/\partial\bar T = 0.053$ on
  the Ferron rig) but grows without bound as `bush_radial_cells` is refined, so
  the damping of the film energy equation would depend on the bush mesh;
* the series-to-ambient conductance $G$ is mesh-free but far too small: on a
  nearly adiabatic bush it drives the bracket negative, the clamp zeroes the wall
  coefficient, and with the wall term fully explicit the temperature KSP diverged
  outright on the Ferron rig at 720 x 192.

This affects only the path to the fixed point; the converged answer always uses
the solved $T_b$.

**Coupled-loop convergence (measured, Ferron rig, 720 x 192).** With a
`FIXED_TEMPERATURE` outer surface the Picard loop converges normally (21–27 outer
iterations). With a weak `CONVECTION` outer surface — $h = 80$ or
$10\ \mathrm{W/(m^2K)}$, i.e. a bush that is almost globally adiabatic, so the
film's temperature level is set by a very weak global balance — the loop
**limit-cycles** instead of reaching `coupling_tolerance`. The residual
$\max(R_p,R_T)$ decays geometrically to $\approx 4.3\times10^{-5}$ at outer 56 and
then oscillates between $5\times10^{-5}$ and $3.6\times10^{-4}$ for as many sweeps
as it is given, never approaching the $10^{-6}$ tolerance. In physical units the
cycle is about 27 Pa of pressure ripple against a 1.4 MPa peak and 5 mK of
temperature ripple. At 360 x 96 both variants converge normally, so this is a
property of the segregated loop on this grid, not of the bush model.

**The cycle is a loop gain marginally above one, and under-relaxation removes it.**
Measured on the `rerun_cht_robin10` equilibrium (24 ranks, afml-r6, sweeps to
`coupling_tolerance` $=10^{-6}$):

| $\omega$ = `coupling_relaxation` | 1.0 | 0.95 | 0.9 | 0.8 | 0.7 | 0.5 | 0.3 | ADAPTIVE |
|---|---|---|---|---|---|---|---|---|
| outer sweeps | **200, never** | 48 | 53 | 64 | 76 | 113 | 197 | **52** |

Even $\omega = 0.95$ converges, which locates the loop gain just above 1. The
ADAPTIVE schedule (`coupling_relaxation_mode = ADAPTIVE`, §7.4) reaches the
tolerance in 52 sweeps without a tuned $\omega$: it starts at the ceiling, halves
the factor on the first residual rise and lets it recover, ending at
$\omega = 0.974$ with a minimum of 0.5.

**The limit cycle did not corrupt the archived answers.** Every converged variant
above and the non-converged $\omega=1$ run agree on the delivered fields: the
bush-side wall temperature profile differs by at most $10^{-3}$ K (the VTK write
precision), $W$ by 0.11 N in 4002 N, $p_{max}$ by 40 Pa in 1.405 MPa, $T_{max}$ by
1 mK, and $h_{min}$ and the cavitated fraction to every printed digit. The only
visibly mobile quantity is the *argmax cell* of a nearly flat maximum, which hops
$\pm1$ cell (0.5°). Details and the full residual traces:
`cases/01_ferron_wallbc_fields/PICARD_FIX_REPORT.md`.

The shaft stays globally isothermal in both implemented bush models. Ferron's
shaft is circumferentially uniform but varies in radius and axial position,
with axial-end temperature 43°C. A global scalar shaft is a further reduction,
not the same boundary-value problem.

### 13.1.3 Generalized (Dowson) Reynolds Equation

The same $T(\eta)$ gives $\mu(\eta)$, so the Poiseuille coefficient is no longer
$h^3/12\mu$. Integrating the thin-film momentum balance
$\partial p/\partial x = \partial(\mu\,\partial u/\partial y)/\partial y$ twice
with $u(0) = U$, $u(h) = 0$ gives the exact flux

$$q = -F_2\frac{\partial p}{\partial x} + U\frac{F_1}{F_0},$$

$$F_0 = \int_0^h \frac{\mathrm{d}y}{\mu}, \qquad
  F_1 = \int_0^h \frac{y\,\mathrm{d}y}{\mu}, \qquad
  F_2 = \int_0^h \frac{y\left(y - F_1/F_0\right)}{\mu}\,\mathrm{d}y .$$

$F_2$ **replaces** $h^3/12\mu$, and $2F_1/F_0$ replaces $h$ in the Couette flux.
Both reduce to their classical values identically when $\mu$ does not vary across
the gap, which `tests/test_cross_film.cpp` asserts to machine precision. This is
not a correction of $h^3/12\mu$ by an effective viscosity: because $\mu$ is
convex in $T$ and the core is the hottest part of the film, the effective
viscosity implied by $F_2$ is *lower* than $\mu(\bar T)$, and the load falls.

The dissipation source becomes exact under the same integrals:

$$\Phi = \frac{U^2}{F_0} + F_2\left|\nabla p\right|^2 .$$

The Couette-Poiseuille cross term integrates to zero identically, so the two
contributions simply add, and both reduce to $\mu U^2/h$ and
$h^3|\nabla p|^2/12\mu$ at uniform viscosity. The energy equation's convective
flux uses the same $q$ the pressure solve conserved, so enthalpy and mass ride on
one velocity field.

The $F$-integrals are evaluated by Gauss-Legendre quadrature
(`cross_film_quadrature_points`, default 6) on
$\mu(\eta) = \mu_{cell}\,\mu_l(T(\eta))/\mu_l(\bar T)$: the cell's own viscosity
field already carries the pressure, dissolved-gas and mixture corrections, so
only the temperature ratio is applied on top of it and every property model keeps
working through the reconstruction. The ratio is clamped to $[0.05, 20]$ as a
numerical guard; on the validation rigs it stays inside about $[0.6, 1.6]$.

`CrossFilm::update` runs once per outer (Picard) iteration, before the Reynolds
solve, and the caller exchanges ghosts for $F_0$, $F_1$, $F_2$ and $2F_1/F_0$
because the Reynolds assembly reads them on both theta ghost layers.

Both cavitation models use these coefficients. The Elrod–Adams path carries them
in $\Gamma_{base} = \rho_0\bar\beta F_2$ and in its Couette face flux; the FBNS
path carries them in $\rho_l F_2$ and the same Couette thickness (§7.3.9). The
single definition lives in `src/reynolds_coefficients.hpp`.

### 13.2 Temperature-dependent Fluid Properties

Thermal hydrodynamic coupling enters Reynolds through field-valued liquid
solution density and viscosity. The Reynolds diffusion coefficient becomes:

$$\gamma(\theta,z) = \frac{\rho(T,p)h^3}{12\mu(T,p)}$$

(replaced by $\rho(T,p)\,F_2$ when the cross-film model of §13.1.3 is active,
and the Elrod form below likewise by $\rho_0\beta F_2$ — see §13.1.3)

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

For `SINGLE_PHASE` or `TWO_PHASE`, dissolved gas is stored as
a liquid-solution mass fraction:

$$c_d = \frac{m_{gas,dissolved}}{m_{liquid\ solution}}$$

The equilibrium concentration is selected by `solubility_model`:

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

where `dissolved_gas_liquid_density` supplies the effective partial liquid-volume
density $\rho_{dg,liq}$. MIXTURE requires a positive value. A free-gas EOS is
not a dissolved-species partial-volume model; the former fallback was removed.

For `EMPIRICAL`, the liquid-solution viscosity is:

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

#### Polynomial (Ferron) viscosity — `liquid_viscosity_model = POLYNOMIAL`

Ferron, Frene & Boncompain (*ASME J. Lubr. Technol.* **105** (1983) 422) close
their THD journal-bearing model with a nondimensional quadratic viscosity law,

$$\bar\mu = \frac{\mu}{\mu_i} = k_0 - k_1\bar T + k_2\bar T^2 ,
\qquad \bar T = \frac{T}{T_i},$$

so that in this solver

$$\mu_\ell(T) = \mu\,\left(k_0 - k_1\bar T + k_2\bar T^2\right).$$

**Unit convention — $\bar T$ is a ratio of temperatures in degrees Celsius, not
Kelvin.** This is not a stylistic choice; the coefficients only work one way.
The paper's nomenclature reads "$T$ = temperature °C" and "$\bar T = T/T_0$
nondimensional temperature", with $T_0$ the inlet film temperature in °C, and
its Table 1 gives the ISO VG32 set

$$k_0 = 3.287,\quad k_1 = 3.064,\quad k_2 = 0.777,
\qquad \mu_i = 0.0277\ \mathrm{Pa\,s\ at\ } T_i = 40\ ^\circ\mathrm{C}.$$

Read as a Celsius ratio these reproduce the oil: at $40\ ^\circ$C, $\bar T = 1$
and $k_0 - k_1 + k_2 = 1.000$ exactly; at $80\ ^\circ$C, $\bar T = 2$ and
$\bar\mu = 0.267$, i.e. $\mu = 7.4$ mPa·s, which is the measured VG32 drop over
that span. Read as a Kelvin ratio the same coefficients give
$\bar T = 353.15/313.15 = 1.128$ and $\bar\mu = 0.82$ — a 18 % drop across
40 K, which is not an oil at all.

`property_reference_temperature` stays in **Kelvin**, like every other
temperature in the configuration; the conversion
$T[^\circ\mathrm{C}] = T[\mathrm{K}] - 273.15$ is applied inside the viscosity
law (`SimulationConfig::polynomial_reference_temperature_celsius`). A reference
at or below $0\ ^\circ$C is rejected by `validate()`, since the Celsius ratio
has no usable normalisation there.

Two further properties of the law follow from it being a *fit*, not a
constitutive relation:

- **It is a pure temperature law.** No Barus pressure term and no dissolved-gas
  term are applied, even if `viscosity_pressure_coeff` or
  `solution_viscosity_gas_coeff` are set — those belong to `EMPIRICAL`.
- **It is only valid inside its fitted range.** Ferron's own quadratic has a
  minimum at $\bar T = k_1/2k_2 = 1.97$ ($\approx 79\ ^\circ$C) and turns
  *upward* beyond it, and a set with $k_2 = 0$ crosses zero and goes negative.
  The implementation floors $\mu_\ell$ at $10^{-3}\mu$ so the Reynolds mobility
  stays positive-definite, and logs one warning per run when the floor is hit.

Because $\mu$ is defined *at* the reference temperature, $k_0 - k_1 + k_2$
should equal 1. `validate()` warns when it deviates by more than 5 %.

The field written as `mu_liquid_solution` is this liquid-solution viscosity.
For `TWO_PHASE`, the effective film viscosity written as `mu` is selected by
`mixture_viscosity_model` (all satisfy $\bar\mu(\alpha_g{=}0)=\mu_l$; mass quality
$x_g = m_g/(m_g + \rho_l\theta h)$):

$$\bar\mu =
\begin{cases}
\alpha_g\,\mu_g + (1-\alpha_g)\,\mu_l, & \text{DUKLER}\\
x_g\,\mu_g + (1-x_g)\,\mu_l, & \text{LINEAR}\\
\left(\dfrac{x_g}{\mu_g} + \dfrac{1-x_g}{\mu_l}\right)^{-1}, & \text{MCADAMS}
\end{cases}$$

`MCADAMS` (the default) is the homogeneous two-phase standard, with correct limits
at both $x_g\to 0$ and $x_g\to 1$. `LINEAR` is the quality-weighted form quoted by
Grando, Priest & Prata (2006) ($\bar\mu = \chi\mu_g + (1-\chi)\mu_l$),
refrigerant-oil validated. `DUKLER` is the legacy void-weighted linear blend, kept
for comparison with prior results. All three require the free-gas viscosity
`mu_gas`. The shipped R290/PZ68 cases select `MCADAMS` with `gas_alpha_max = 0.6`.
(The earlier `EINSTEIN_DILUTE` — a dilute-suspension result that thickens with gas,
directionally wrong as $\alpha_g\to 1$ — and the `KRIEGER_DOUGHERTY` packing law
were removed.)

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

When `fluid_property_model = TWO_PHASE`, finite-rate
release/resorption updates dissolved, free and liquid inventories (§18).
The former explicit approximation was:

$$\Delta m_g =
k_{dg}(c_d - c_{d,eq})\,\rho_l h\,\alpha_l\,\Delta t$$

with sign convention `gas_mass_transfer > 0` for release from solution into
free gas and `gas_mass_transfer < 0` for resorption. The implementation caps
the reported/resolved free-gas volume fraction by `gas_alpha_max`:

$$0 \le \alpha_g \le \alpha_{g,max}$$

Dissolved-gas transport is a conservative segregated advection-diffusion step on
$c_d$. FBNS uses its exact discrete liquid flux; JFO retains a reconstructed
carrier velocity. Free gas mass per bearing
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
requires `dissolved_gas_henry_coeff > 0`, MIXTURE requires
`dissolved_gas_liquid_density > 0`, empirical/log viscosity requires a non-zero
`solution_viscosity_gas_coeff`, and `TWO_PHASE` requires JFO or FBNS. These
configuration checks do not certify physical gas-pressure coupling.

This is a segregated nonlinear JFO/property/gas update, not a PISO method. A
true PISO label would require a Navier-Stokes pressure-velocity backend.

### 13.3 Elastic and Thermal Deformation

This subsection describes proposed equations. The current geometry routine
does not implement these deformation terms; a measured hot clearance can
instead be supplied as a fixed geometric input.

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
computes velocities before force and torque post-processing. The one exception is
`mass_residual` (§15.1), a pure diagnostic that is registered and evaluated only
when selected, so an unselected run is bit-identical to one built without it. When enabled,
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
$\Gamma_{base}$, and inlet reservoir cells are re-evaluated as physical operator
residuals (the mass the prescribed-value row injects). The normalized residual

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

### 15.1 Per-cell continuity residual (`mass_residual`)

A global balance that closes to round-off cannot distinguish a cell-by-cell
conservative scheme from large local errors that cancel on summation. The
optional output field `mass_residual` (select it in `output_fields`; JFO and FBNS
only) exposes the *per-cell* liquid-continuity residual of the converged state, in
**kg/s, cell-integrated — not area-normalized**, so it carries the same units as
every rate in the global balance. For cell $P$ with faces $f \in \{e,w,n,s\}$,

$$r_P \;=\; \underbrace{\sum_{f} \Big[ \Gamma_f\,(\theta_P - \theta_{F}) \frac{A_f}{\delta_f} \;+\; F_f\,\theta_f \Big]}_{\text{Poiseuille} \,+\, \text{Couette outflow}}
\;+\; \underbrace{\rho_l h\,\frac{\theta_P - \theta_P^{\,n}}{\Delta t}\,\Delta A \;+\; \rho_l\,\theta_P^{\,n}\,\frac{\partial h}{\partial t}\,\Delta A}_{\text{storage} \,+\, \text{squeeze (transient only)}}
\;-\; \frac{m_{clamp,P}}{\Delta t}$$

with exactly the quantities the assembly uses and no re-derived approximation:
$\Gamma_f = \mathrm{harm}\big(\Gamma_{base}g(\theta)\big)$ on interior faces,
$F_f = \bar\rho_0\,\omega R\,\bar h\,\Delta z/2$ for the Couette face flux,
$\theta_f$ from the configured `theta_convection_scheme` (upwind, central,
Type-differencing on full-film faces, or the van Leer limiter), and
`SimulationConfig::elrod_boundary_outflow` on an active axial Dirichlet or
`INLET_OUTLET` face. Interior faces telescope between neighbours, so summing the
field over the whole domain returns the global closure exactly:

$$\sum_P r_P \;=\; S_{storage} - \Phi_{boundary} - S_{inlet} - \frac{M_{clamp}}{\Delta t}
\;=\; \frac{r_l\,M_l}{\Delta t}$$

Inlet reservoir cells carry only the clamp term. Their flux imbalance *is* the
modelled supply source, which the global balance already reports separately as
$S_{inlet}$; leaving it in the field would double-count it in the sum above.
The inlet row is applied after the axial boundary rows. At a full-length groove,
it therefore overwrites the axial boundary link in the two overlapping end
cells. Both $\Phi_{boundary}$ and $S_{inlet}$ omit those discarded faces. Adding
them back only in diagnostics creates a contradictory pressure jump over
$\Delta z/2$ and a false groove source proportional to $1/\Delta z$.

**Under `cavitation_model = FBNS`** the field mirrors the FBNS operator instead
(`FBNS::cell_outflow`), because that solver assembles different face fluxes: the
Poiseuille term is $\mathrm{harm}(\rho_l h^3/12\mu)\,(p_P-p_F)\,A_f/\delta_f$ —
driven by pressure, with no bulk modulus and no switch function — the Couette face
value is always first-order upwind, and axial Dirichlet faces use the plain
half-cell pressure link of §7.3.4. There is no clamp term: the FBNS truncations act
on Newton iterates, not on the converged state, so `cavitation_clamp_mass` stays
zero and the sum identity reads $\sum_P r_P = -\Phi_{boundary} - S_{inlet}$ in
steady state. Inlet cells are zeroed for the same reason as above.

**One definition of the coefficients.** The residual above is a *measurement* of
the assembled operator, so every coefficient in it — $\Gamma_{base}$, the switch
function's full-film ceiling, the Couette face flux, the FBNS mobility — is taken
from `src/reynolds_coefficients.hpp`, the same header `Reynolds::solve_elrod` and
`FBNS::solve` assemble from. The transport equations draw the *same* face
velocities from that header: the energy equation's axial boundary faces and every
face of the dissolved-gas and free-gas transport now call
`ReynoldsCoeff::couette_face_velocity` and
`ReynoldsCoeff::poiseuille_face_velocity`, so under the cross-film model they
carry $-g F_2\,\partial p/\partial s$ and $U F_1/F_0$ rather than the classical
$-g h^3/(12\mu)\,\partial p/\partial s$ and $U h/2$. Before that, the energy
equation transported heat through its *interior* faces with $F_2$ and out through
its *axial ends* with $h^3/12\mu$ — on a case with Dirichlet axial ends, which is
the Ferron cross-film geometry, those mismatched faces are the oil outflow. Both
helpers return the classical expression verbatim (same floating-point
association) when the cross-film model is off, so non-cross-film runs are
bit-identical. Local copies had drifted: the diagnostics still used
$\rho_0\beta h^3/12\mu$ with the *scalar* `bulk_modulus` and $h$ in the Couette
flux after the assembly had moved to $\rho_0\bar\beta F_2$ and $2F_1/F_0$, so on
exactly the cross-film and `VOID_COUPLED` runs this field exists to police it
reported a fabricated imbalance of the size of that difference (measured at ~10 %
of the mobility on the cross-film test bearing). `tests/test_fbns.cpp` test 14
holds both paths to round-off with the cross-film model on, fronts included.

Reading the field: a cell that is not a groove cell, not clamped and not on the
rupture or reformation front must sit at round-off of its own row — the
barotropic conductance $\rho_0\beta h^3/(12\mu)\,A_f/\delta_f$, not the Couette
throughput, is what sets that scale, because it exceeds the throughput by
$O(\beta/p)$. `tests/test_mass_residual.cpp` measures $7\times10^{-15}$ of the
row scale on a steady grooved bearing. Anything larger in such a cell is a defect
in the assembled operator, not a physical source. The same test perturbs one
converged cell and confirms the global $r_l$ does not move at all while the field
rises to $2\times10^3$ times the Couette flux in the five affected cells — the
blindness this field exists to remove.

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

### 15.2 Film-surface outputs (`temperature_journal`, `heat_flux_journal`, …)

The written `temperature` field is $\bar T$, the depth average of §13.1 — under
the cross-film model it is the *mean* of a profile whose two ends can differ by
several kelvin, so it alone cannot be compared against a wall thermocouple. Four
optional output fields expose the ends of that profile. They require
`cross_film_thermal_model = POLYNOMIAL`; selecting one under `NONE` is a
validation **error**, because there is no profile to evaluate and a field of
zeros would read as a result.

| `output_fields` key | VTK name | definition | unit |
| --- | --- | --- | --- |
| `temperature_journal` | `T_journal` | $T(\eta = 0)$ | K |
| `temperature_bearing` | `T_bearing` | $T(\eta = 1)$ | K |
| `heat_flux_journal` | `q_journal` | $-\dfrac{k}{h}\left.\dfrac{\partial T}{\partial\eta}\right|_{\eta=0} = -q_{shaft}$ | W/m² |
| `heat_flux_bearing` | `q_bearing` | $+\dfrac{k}{h}\left.\dfrac{\partial T}{\partial\eta}\right|_{\eta=1} = -q_{bush}$ | W/m² |

The temperatures are the profile of §13.1.2 evaluated at the two walls, and
therefore reproduce the conjugate wall temperatures $T_s$ and $T_b$ exactly —
they are written from the profile rather than copied from the wall solve so the
two can never disagree. The temperatures are absolute kelvin, like every
temperature in the solver.

**Sign convention: the fluxes are positive INTO the fluid.** This is the
*opposite* sign to the internal $q_{shaft} = (k/h)(a+b)$ and
$q_{bush} = (k/h)(b-a)$ of §13.1.2, which are positive for heat *leaving* the
film; the output fields are their negatives. The opposite signs on the same
$\partial T/\partial\eta$ in the table are geometric: at $\eta = 0$ the fluid
lies at larger $\eta$, at $\eta = 1$ at smaller $\eta$. Consequences worth
checking a result against: a wall colder than the film gives a **negative**
into-fluid flux there, an adiabatic wall ($h_{out} = 0$ for the bush) gives an
**identically zero** field cell by cell, and an adiabatic *shaft*
($h_{shaft} = 0$) gives a locally nonzero `heat_flux_journal` whose area
integral is zero — that upstream heat path is the whole point of the conjugate
shaft. `tests/test_cross_film.cpp` asserts all of it, including the equality
with the closed-form profile and the fact that an unselected run is
byte-identical: the fields are allocated only when selected, so they cost
nothing and change nothing when they are not.

---

## 16. Oil–Refrigerant Gaseous Cavitation (PZ68S + R290)

This section records the historical proposed model for PAG-class PZ68S oil
saturated with R290/propane. It is not a general thermodynamic law or a validated
model for the current oil-R410A target. The max-threshold and local plateau
construction below remains experimental; §18 states current limits. The full survey, the
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
fraction. This is why the propane configs set `cavitation_threshold = LOCAL_PSAT`
(so the Elrod plateau is $\max(p_{cav}, p_{sat})\approx1$ MPa). With the code
default `SCALAR_PCAV` the threshold is instead pinned at the oil value
$p_{cav}=3.0\times10^4$ Pa while propane releases at $\sim$1 MPa — the cavitated
zone then sits at 30 kPa and is far too small (a verified 30 kPa vs 1 MPa plateau,
68 vs 506 cavitated cells on the shipped case). Both the threshold and the WP-1
`gas_pressure_coupling` are exposed in the GUI Cavitation panel.

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

With WP-1 void coupling the cavitated plateau generalizes further to a gas-filled
value $p_{void}$, of which $p_{cav,\mathrm{eff}}$ is the gas-free floor. The
plateau is obtained as the **fixed point** of the gas inventory rather than
evaluated explicitly from the current free-gas mass. Requiring the gas released
between the bubble point and the local pressure to exactly fill the void
$(1-\theta)$ at ideal-gas density $\rho_g=p/(R_sT)$ gives

$$H\,(p_{sat}-p)\,\rho_l\,\theta\,h\cdot\frac{R_sT}{p} = (1-\theta)\,h
\quad\Longrightarrow\quad
\boxed{\;p_{void}=\frac{L\,\theta\,p_{sat}}{(1-\theta)+L\,\theta}\;}$$

where

$$L \equiv H\,\rho_l\,R_s\,T$$

is the dimensionless **Ostwald solubility coefficient** (gas volume dissolved per
unit liquid volume) — the Henry slope $H$, the liquid density and the species gas
constant collapse into that single group, so the plateau depends on the dissolved
species *only* through $L$. Limits: $p_{void}\to p_{sat}$ as $\theta\to1$,
$p_{void}\to0$ as $\theta\to0$, monotone in between, and never above $p_{sat}$.

This replaces the earlier explicit form $p_{void}=m_gR_gT/((1-\theta)h)$, which is
the same balance but unbounded: the flash (§16.2) releases gas wherever
$p<p_{sat}$ regardless of film content, so a full-film cell could carry gas and
divide it by a floored, near-zero void, yielding $O(10^4)$ bar plateaus. Note the
tempting cap "free gas $\le$ the available void $(1-\theta)$" is **not** a valid
repair: it suppresses outgassing in a supersaturated *full* film, which is the
phenomenon the WP-1 pressure-cap gate asserts. Gas creates its own void by
displacing liquid.

Because $\alpha_g$, $\theta_{full}=1-\alpha_g$ and $\theta$ are mutually
dependent, the coupled iteration limit-cycles if taken undamped — the rupture
front oscillates with zero flag flips and the Elrod loop stalls at a finite mass
residual. `void_coupling_relaxation` $w$ (default $0.3$) under-relaxes $p_{void}$,
$\alpha_g$ and $m_g$ toward their new values each outer iteration, which fixes the
*inner* Elrod loop.

> **Status: `VOID_COUPLED` does not converge and must not be used for
> quantitative work.** The *outer* segregated coupling between the flash
> (`equilibrate_gas`) and the Elrod solve exhausts `coupling_max_iters` at
> $R_p\sim10^0$–$10^1$ against a $10^{-9}$ tolerance. The diagnostic is
> relaxation dependence: at a fixed mesh, $w=0.10/0.30/0.60$ give peak pressures
> of $12.1/6.29/6.42$ bar and cavity arcs of $64/100/108^\circ$. Since
> under-relaxation cannot move a fixed point, those states are not fixed points.
> Refining the mesh does not help. The cure is a numerical-methods change —
> an implicit/Newton coupling, or folding the flash into the Elrod switch so
> $\alpha_g$ and $\theta$ are solved simultaneously rather than alternately.
> Use `gas_pressure_coupling = NONE` meanwhile.

The closed form above remains exact and useful on its own, independent of the
solver: it can be applied *a posteriori* to a converged `NONE` solution to test
whether the dissolved-gas inventory could sustain the void that solution
predicts.

**Consequence for gaseous cavitation.** At a given void fraction the sustainable
plateau is set by solubility alone. With $p_{sat}=1$ bar and $\alpha=0.4$,
air ($L=0.098$) supports only $0.13$ bar — effectively a vaporous cavity — while
CO$_2$ ($L=1.08$) supports $0.62$ bar, a genuinely gaseous one. This is why a
weakly soluble gas cannot produce a gas-filled cavity at moderate void, however
much the Elrod solution says the film has ruptured.

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

## 17. Configurable Model Options (quick reference)

Each model selector below is a config key. Values are listed **simplest → most
physical / accurate**; the first-listed or **bold** value is the default. The GUI
hover text mirrors these; this section is the canonical explanation.

### 17.1 Cavitation / film rupture — `cavitation_model`

| Value | Treatment | Use |
|---|---|---|
| `FULL_SOMMERFELD` | No cavitation; sub-cavity (tensile) $p$ kept. Antisymmetric Sommerfeld solution. | **Reference baseline** for overlay plots only (not a real operating point). |
| `GUMBEL` | Half-Sommerfeld: solve, then clamp $p\!\ge\!p_{cav}$. Non-conservative. | Fast, robust load-capacity estimates. |
| **`JFO`** | Mass-conserving Jakobsson–Floberg–Olsson via Elrod–Adams film content $\theta$ (§4, §7.2): $p=p_{cav}+\beta\ln\theta$ for $\theta\!\ge\!1$, $\theta<1$ cavitated. | The production model for thermal / gas-coupled runs (physical rupture **and** reformation). |
| **`FBNS`** | The same JFO conditions as a complementarity problem in $(p,\theta)$, solved by Fischer–Burmeister–Newton–Schur (§7.3). No bulk modulus, no switch function. | Isothermal, constant-property runs where free-boundary accuracy matters: front error $\le\!1$ cell vs $7$–$13$ for `JFO`, at 4 Newton iterations. |

### 17.2 Lubricant property fidelity — `fluid_property_model`

| Value | What varies | Free gas? |
|---|---|---|
| **`CONSTANT`** | Fixed $\mu,\rho$. | no |
| `SINGLE_PHASE` | Liquid $\mu(c_d,T,p)$, $\rho(c_d)$ vary; gas stays in solution. | no (cavity is the classic massless JFO void) |
| `TWO_PHASE` | Adds finite-rate gas release/resorption and a free-gas phase in cavitated zones (§16). | yes — needs `JFO` + gas params |

Note: film rupture (the cavitation model) operates in **every** tier; the tier only
decides whether the released phase is *tracked*. The dominant oil–refrigerant
mechanism (gaseous outgassing at $p_{sat}$) needs `TWO_PHASE` + `gas_pressure_coupling`.

### 17.3 Solubility — `solubility_model`

| Value | Law |
|---|---|
| **`HENRY`** | $c_{sat}=H(T)\,p$, $H(T)=H_{ref}e^{E_H(1/T-1/T_{ref})}$ (van 't Hoff). |
| `TABLE` | Measured $c_{sat}(p,T)$ surface (1-D or 2-D PTSV data); bubble point by first-crossing inversion (§16.4). |

### 17.4 Liquid-solution density — `density_model`

| Value | $\rho_\ell$ |
|---|---|
| **`CONSTANT`** | $\rho_{oil}$ (dissolved gas does not change density). |
| `MIXTURE` | Mass/volume mixing of two **explicit** phase densities: oil `rho` and `dissolved_gas_liquid_density`. $1/\rho_\ell=(1-c_d)/\rho_{oil}+c_d/\rho_{g,\ell}$. |
| `TABLE` | Measured $\rho(p,T)$. |

### 17.5 Liquid-solution viscosity — `liquid_viscosity_model`

| Value | $\mu_\ell$ |
|---|---|
| **`CONSTANT`** | $\mu_{oil}$. |
| `EMPIRICAL` | $\mu_{oil}\,e^{E_\mu(1/T-1/T_{ref})}\,e^{a_c c_d}\,e^{\alpha_p(p-p_{ref})}$ (Andrade × dissolved-gas thinning × Barus). Set $E_\mu=\alpha_p=0$ for pure gas-thinning (the old log-mixing law). |
| `POLYNOMIAL` | $\mu_{oil}\,(k_0-k_1\bar T+k_2\bar T^2)$, Ferron, Frene & Boncompain (1983), with $\bar T=T[^\circ\mathrm{C}]/T_{ref}[^\circ\mathrm{C}]$ — a **Celsius** ratio. Pure temperature law: no pressure or dissolved-gas term. Coefficients `viscosity_poly_k0/k1/k2` (default 1, 0, 0 = identity); their VG32 set is 3.287, 3.064, 0.777 with $\mu=0.0277$ Pa·s and `property_reference_temperature = 313.15`. See §13.2.1. |
| `TABLE` | Measured kinematic $\nu(p,T)$, converted to dynamic per cell with $\rho_\ell$. |

### 17.6 Mixture (free-gas) viscosity — `mixture_viscosity_model`

All satisfy $\mu(\alpha\!=\!0)=\mu_\ell$; mass quality $x=m_g/(m_g+\rho_\ell\theta h)$.

| Value | $\mu_{mix}$ | Note |
|---|---|---|
| `DUKLER` | $\alpha\,\mu_g+(1-\alpha)\mu_\ell$ | Void-weighted linear. **Legacy/deprecated** (kept for comparison). |
| `LINEAR` | $x\,\mu_g+(1-x)\mu_\ell$ | Grando 2006; refrigerant-oil validated. |
| **`MCADAMS`** | $1/\mu=x/\mu_g+(1-x)/\mu_\ell$ | Homogeneous standard; correct asymptotes. Default. |

(The earlier `EINSTEIN_DILUTE` — which thickens with gas, wrong as $\alpha\!\to\!1$ — and
`KRIEGER_DOUGHERTY` suspension law were removed.)

### 17.7 Convection scheme — `theta_/thermal_/gas_convection_scheme`

| Value | Face interpolation | Property |
|---|---|---|
| **`UPWIND`** | 1st-order upwind. | Bounded, diffusive (smears fronts). |
| `LINEAR` | 2nd-order central. | Sharp but **unbounded** (oscillates at fronts) — reference/smooth-field use. |
| `VANLEER` | van Leer TVD limited. | 2nd-order **and** bounded — recommended for cavitation fronts. |
| `TYPE_DIFFERENCING` | Central on full-film faces, upwind on cavitated faces (Vijayaraghavan & Keith 1989). | JFO film-content ($\theta$) only. |

### 17.8 Bearing motion integrator — `motion_time_method`

Applies only to the rigid-body equation of motion (field solves are backward-Euler).
`EULER_EXPLICIT`, **`EULER_IMPLICIT`** (default, stable), `CRANK_NICOLSON`, `RK2`,
`RK4` (accurate orbits). The three middle/explicit options are kept as references.

### 17.9 Axial thermal inflow — `bc_z_*_thermal`

| Value | Incoming oil temperature | Situation |
|---|---|---|
| **`ZERO_GRADIENT`** | $\partial T/\partial z=0$ (local film temperature). | Submerged / open end of the same recirculating oil. |
| `CONSTANT` | Fixed `temperature_reference`. | Fed supply / reservoir end at a controlled temperature. |

Feature gates (default-off ≡ historical behavior): `cavitation_threshold`
(`SCALAR_PCAV`/`LOCAL_PSAT`, §16.1), `gas_pressure_coupling` (`NONE`/`VOID_COUPLED`,
§16.3/WP-1), `cavitated_film_model` (`FULL_FILM`/`STRIATED`), `steady_gas_model`
(`FROZEN`/`EQUILIBRIUM`), `consistent_boundary_flux` (§6).

## 18. Conservative gas update and experimental scope (2026-09-06)

### 18.1 Inventories and flow storage

Use bearing-area inventories in kg/m²:

$$M_l=\rho_l h\theta,\qquad M_d=c_dM_l,\qquad M_g=m_{g,free},
\qquad M_o=(1-c_d)M_l.$$

Here $M_l$ includes dissolved species, $M_o$ is base oil, and $M_g$ is the
resolved free species. The total mixture inventory is $M_l+M_g$; the total gas
species inventory is $M_d+M_g$. The `rho` field is $\rho_l\theta$, the liquid
inventory per gap volume, rather than the total two-phase density.

JFO and FBNS flow storage is now

$$\frac{M_l^{n+1}-M_l^n}{\Delta t}\Delta A.$$

The driver snapshots `rho_liquid_solution`, `h` and `theta` once per physical
step. Geometry change is already included, so adding $\rho_l\theta\dot h$
again would double-count squeeze. Low-level callers without a stored gap may
use the legacy prescribed-squeeze fallback $h^n=h-\Delta t\dot h$; new drivers
should store the actual previous geometry. Frozen-coefficient Newton/Picard
still leaves pressure/property coupling error until the outer state closes.

### 18.2 Conservative transport and phase exchange

Before local reaction, dissolved species solves

$$\frac{c_d^* M_l^*-c_d^n M_l^n}{\Delta t}\Delta A
+\sum_f F_{l,f}c_{d,f}^* - \sum_f D_f M_{l,f}\nabla c_d^*\cdot\mathbf n_f A_f=0.$$

FBNS supplies its exact harmonic-mobility Poiseuille and upwind Couette face
fluxes $F_{l,f}$. The former subtraction of $c_d\nabla\cdot F_l$ is gone.
Free inventory uses conservative advection at a reconstructed homogeneous
liquid velocity. This no-slip closure does not model bubble migration or
coalescence. Reservoir rows prescribe feed concentration and zero free gas
inside the solves; displaced-row residuals measure the imposed gas source.

For reaction at frozen $p,T$, let $c_e=c_{d,eq}(p,T)$ and $k$ be the transfer
rate. The implemented local kinetic law is

$$\dot M_l=\dot M_d=-k(M_d-c_eM_l),\qquad
\dot M_g=k(M_d-c_eM_l).$$

Thus $M_o=M_l-M_d$ and $M_d+M_g$ are invariants. Define
$q=M_d-c_eM_l$; then $\dot q=-k(1-c_e)q$. The exact unrestricted transfer is

$$\Delta M=M_l(c_d-c_e)
\frac{1-\exp[-k(1-c_e)\Delta t]}{1-c_e}.$$

The $c_e\to1$ limit of the quotient is $k\Delta t$. `expm1` and a small-argument
expansion retain accuracy for short steps. Bound release by the available
dissolved mass and configured free-gas capacity; bound absorption by existing
free mass. A capacity-limited release leaves supersaturation rather than
discarding species. Then reconstruct

$$M_l'=M_l-\Delta M,\quad M_d'=c_dM_l-\Delta M,\quad
M_g'=M_g+\Delta M,\quad c_d'=M_d'/M_l',\quad
\theta'=M_l'/[\rho_l(c_d',p,T)h].$$

`gas_mass_transfer` is $\Delta M/(h\Delta t)$, positive for release. It is a
sink in liquid continuity and a source in free gas; it cancels from total gas
species conservation. No numerical relaxation is applied to physical transient
exchange. Analytic local integration prevents kinetic overshoot but does not
make the split flow/transport/energy integration exact.

### 18.3 Steady flash and property consistency

The steady feed flash is a local composition closure. For feed mass fraction
$c_f=(c_dM_l+M_g)/(M_l+M_g)$,

$$M_g=\frac{c_f-c_d}{1-c_f}M_l.$$

At an unrestricted flash $c_d=\min(c_e,c_f)$. If the volume cap limits $M_g$,
or steady numerical relaxation changes it, reconstruct
$c_d=c_f-(1-c_f)M_g/M_l$. This preserves the specified local feed fraction.
It is not conservative steady transport through the bearing.

The two-dimensional viscosity table stores kinematic $\nu(p,T)$, so the shared
property helper returns $\mu=\rho_l\nu$. Both mean-film and cross-film sampling
now call that helper. A saturation table contains composition implicitly and
cannot establish off-saturation $\mu(p,T,c_d)$ without additional data.

### 18.4 Admissibility and physical scope

At the gas-update state define

$$\alpha_{EOS}=\frac{M_g}{\rho_g h},\qquad
E_V=\max_{cells}[0,\theta+\alpha_{EOS}-1],\qquad
E_{EOS}=\max_{cells}|\alpha_g-\alpha_{EOS}|.$$

`gas_volume_excess` and `gas_eos_defect` report these dimensionless diagnostics.
They use the stored gas-step `rho_gas`; a subsequent temperature update may
introduce additional thermodynamic lag. A cap may hide gas volume in the
reported alpha, but cannot hide it from these diagnostics. The mild synthetic
air test already gives $E_V\simeq2.7\times10^{-4}$. Conservation is therefore
verified separately from volume admissibility.

FBNS currently constrains pressure against a scalar equilibrium plateau using
lagged gas volume. Pressure, species masses and volume must be solved
consistently before claiming a physical gaseous cavity pressure. The transient
driver has one segregated sweep; steady convergence additionally checks
$R_g=\max(\|\Delta c_d\|_\infty/\max(c_{d,initial},10^{-6}),
\|\Delta\alpha_g\|_\infty)$. This detects an evolving gas state but is not an
EOS/volume closure test. Reported Newton convergence applies to its inner solve.

Classical JFO is an equilibrium liquid-content model. There is no independent
vapor inventory, multicomponent gas EOS, nucleation model or latent/dissolution
enthalpy term here. AIR and PROPANE/R290 are the only parsed free species;
R410A is rejected. For oil-R410A, component solubility/partitioning and mixture
properties are required. For oil-air, measured dissolved and entrained air
content and rate data are required. Neither experimental system is validated
by the present regression tests. See the [audit and priorities](docs/PHYSICS_AUDIT_2026-09.md).

### 18.5 Physical time and geometry

The transient driver writes the initial fields at $t=0$ without advancing
them. For each interval, $\Delta t_n=\min(\Delta t,t_{end}-t_n)$, it predicts
the bearing geometry at $t_{n+1}$ using the previous fluid force, updates gap
geometry, and advances flow, species and energy with that same interval.
The completed fields and bearing position are saved at $t_{n+1}$. Floating-point
roundoff near the final time is absorbed into the last interval rather than
creating a tiny extra step. `end_t = 0` writes only the initial state.

This replaces the old convention that performed an advance at label zero and
included an extra interval at the final label. It also aligns the gap change
with old/new liquid storage. The fluid force in the geometry predictor remains
explicit; choosing RK4 for the rigid-body subsystem does not make the coupled
flow/gas/thermal evolution fourth order. New transient comparisons require
timestep convergence of the complete coupled sequence.

The CLI regression checks the exact insulated uniform-film solution
$T(t)=T_0+\mu(\omega R/h)^2t/(\rho c_p)$ and a decoupled constant axial load
$z(t)=F_z t^2/(2m)$, including a half-sized final interval. Initial values,
integrated heating and saved geometry therefore test physical elapsed time,
not only the labels written to the PVD file.

### 18.6 Measured consequence and transient coupling failure

The [paired comparison](cases/90_diagnostics/physics_comparison/SUMMARY.md)
holds input data fixed between pre-audit `b9dcaf3` and repaired `2bc4d3c`.
Steady pure-oil Ferron and optical-rig pressure, temperature and liquid-content
fields are unchanged at saved precision. A closed fixed-p,T reaction test
loses 11.548% of total gas in the old implementation after 20 ms; the repaired
maximum relative inventory drift is 1.33e-15 and follows the analytic law in
section 18.2. This isolates mass loss from boundary transport.

The extended synthetic air bearing test is a failed acceptance case for the
present pressure/volume closure. At 1 ms, maximum gas volume fraction changes
from 0.4966% to 1.6558% after repair. Gas conservation improves, but pressure
and load alternate between physical timesteps. Refining from 10 to 5 to
2.5 microseconds makes the oscillation faster and slightly larger, while
maximum gas fraction changes by only 0.00088 percentage points. A converged
inner Newton solve therefore does not establish a converged transient.

The pressure solve uses lagged alpha in its full-film ceiling; the subsequent
species step changes alpha without re-solving pressure and gap occupancy.
This splitting is consistent with the observed instability, though the paired
case does not isolate every feedback term. At 10 microseconds, excess gas plus
liquid volume reaches 0.003358 of the gap even while the gas mass residual
stays below 1.4e-12. No stable physical load increment, oscillation frequency
or gas-pressure validation can be inferred from those histories. The coupled
pressure, EOS, volume and species acceptance problem remains open.

### 18.7 End reservoirs and asymmetric kinetics (2026-09-07)

The specified experimental target is now FW68D with air or R410A at 20–80 °C,
reservoirs up to 20 bar absolute, and refrigerant charging near 10 bar. Charge
temperature and FW68D mixture data remain required. The full physical
interpretation and literature comparison are in
[FW68D cavitation review](docs/FW68D_CAVITATION_REVIEW.md).

An axial pressure outlet is a bidirectional reservoir. With outward liquid
mass flux $F_l$, the dissolved-species flux is

$$
F_d=\max(F_l,0)c_{d,\mathrm{cell}}+
    \min(F_l,0)c_{d,\mathrm{reservoir}}.
$$

Each end now has its own `bc_z_*_dissolved_gas` (negative sentinel −1 inherits
the global initial concentration), `bc_z_*_temperature` (0 inherits reference
temperature) and `bc_z_*_free_gas_mass_fraction` (default 0). Temperature is in
kelvin. The temperature state is used on inflow only when `bc_z_*_thermal =
CONSTANT`; `ZERO_GRADIENT` extrapolates the interior temperature. An explicit
temperature with `ZERO_GRADIENT` is rejected. The outgoing state remains the
interior state under either thermal mode.

For incoming free-species mass fraction $q_r=M_g/(M_g+M_l)$ and common phase
transport velocity, the incoming free-species mass flux is
$F_g=F_l q_r/(1-q_r)$. This ratio is not a void fraction. Outgoing free-species
flux uses the solved free mass per area and the same carrier velocity. Boundary
diagnostics retain the actual terms used by the linear solves.

The species boundary inputs apply to transient `TWO_PHASE` transport; the
steady local feed flash has no conservative species transport operator. These
inputs do not make different steady reservoir compositions physically supported.

`LOCAL_PSAT` still denotes a global initial bubble-point pressure floor, not a
coupled local thermodynamic state. Previously the boundary helper's maximum
with that floor could alter a specified reservoir pressure even after config
validation. Both end and internal reservoir pressures are now validated
against the **effective** floor. A conflict is rejected. This prevents solving
a different reservoir problem; it does not provide a desorption-pressure EOS.

The local reaction of section 18.2 now uses

$$
k=\begin{cases}k_r,&c_d\ge c_{d,eq},\\k_a,&c_d<c_{d,eq},\end{cases}
\qquad \dot M_l=\dot M_d=-k(M_d-c_{d,eq}M_l),\quad \dot M_g=-\dot M_d.
$$

`gas_mass_transfer_rate` is $k_r$; `gas_absorption_rate` is $k_a$, with −1
inheriting $k_r$ for compatibility. Either rate can be zero independently.
For frozen $p,T,c_{d,eq}$, the transferred mass over a timestep remains

$$
\Delta M=(M_d-c_{d,eq}M_l)
\frac{1-\exp[-(1-c_{d,eq})k\Delta t]}{1-c_{d,eq}},
$$

limited by available donor inventory and the documented release cap. The limit
as $c_{d,eq}\to1$ is evaluated without cancellation. Updating $M_l$, $M_d$ and
$M_g$ with equal and opposite transfers conserves oil and total species mass.
Pressure-dependent gas volume is then evaluated from the EOS, even if both
rates vanish. Thus compression can shrink a bubble without absorbing gas.
Separate rates are needed to test short-cycle observations of release followed
by compression with negligible reabsorption; the literature's paraffin rates
are not assigned to FW68D (Reinke et al.,
*[Fluids](https://doi.org/10.3390/fluids7050159)*, 2022).

R410A remains unimplemented in the production solver. The independent
[EOS diagnostic](cases/90_diagnostics/cavitation_requirements/README.md) uses
CoolProp 7.2.0 to demonstrate that even an ideal R410A gas has material density
errors in the specified operating window, and that phase identity matters.
It provides no FW68D solubility or off-equilibrium solution properties.
The classical equilibrium JFO vapor limit remains available, but no general
joint gas/vapor component-and-energy model is claimed by this addition.

## Steady optical-rig screening added 2026-09-11

The initial case 02b geometry is superseded by the user-confirmed
**D = 43 mm, L = 40 mm and radial clearance c = 42.5 micrometres** in
[case 02c](cases/02c_oil_corrected_geometry/README.md). The equations below remain
the same, but all forces, positions and cavities are recalculated on that
geometry. No previous equilibrium field is rescaled or reused. The groove and
500 Pa floor are retained assumptions, not newly confirmed dimensions or oil
properties.

[Case 02b](cases/02b_oil_load_rpm/README.md) uses the existing isothermal FBNS
closure. With circumferential coordinate $x=R\theta$, liquid filling fraction
$\vartheta$, positive rotation $\omega$ and uniform viscosity $\mu$, the
depth-integrated fluxes are

\[
q_\theta=-\frac{h^3}{12\mu R}\frac{\partial p}{\partial\theta}
          +\frac{\omega R h\vartheta}{2},\qquad
q_z=-\frac{h^3}{12\mu}\frac{\partial p}{\partial z},\qquad
\frac{1}{R}\frac{\partial q_\theta}{\partial\theta}
 +\frac{\partial q_z}{\partial z}=0.
\]

The full-film pressure and cavity fraction obey
$p\ge p_c$, $0\le\vartheta\le1$ and
$(p-p_c)(1-\vartheta)=0$. Both ends and the full-length axial groove impose
$p=100000$ Pa. The assumed floor is $p_c=500$ Pa. The geometric law is
$h=c-e\cos(\theta-\psi)$, and an outer damped Newton search adjusts
$(e/c,\psi)$ to balance both pressure-force components against the applied load.
This outer position search is distinct from the FBNS nonlinear solve and does
not integrate a physical journal trajectory or use a physical timestep.

Scaling with $P_\mu=\mu\omega(R/c)^2$ yields two useful coordinates for the
fixed-geometry scan,

\[
S=\frac{\mu(N/60)}{W/(2RL)}(R/c)^2,\qquad
\Pi_c=\frac{p_{\rm end}-p_c}{P_\mu}.
\]

Thus equal $N/W$ does not alone guarantee an equal cavity, because the absolute
boundary-pressure margin is independently scaled by $N$. Constant density
multiplies the mass flux but adds no inertial mechanism to this Reynolds model.
The plotted $1-\vartheta$ is an equilibrium film-volume deficit. Assigning a low
pressure floor does not calculate phase-change rates, latent heat or dissolved
gas release. The constant 0.210 Pa s viscosity and rigid geometry also omit
thermal and elastic feedback; the computed peak pressure and minimum gap must
be checked before transferring high-load cases to the physical rig.

### Temperature sensitivity requested for the corrected geometry

Case 02c adds a steady energy solve while retaining constant viscosity and the
previously balanced position. It uses the STRIATED liquid-energy closure with
liquid depth $\vartheta h$, viscous dissipation and in-plane heat transport.
For the prescribed wall-cooling branch, the depth-integrated sink is
$H_j(T-T_j)+H_b(T-T_b)$, with $T_j=T_b=T_{\rm in}=293.15$ K and
$H_j=H_b=10000$ W/(m2 K). The other branch uses $H_j=H_b=0$.
The provisional $c_p=2000$ J/(kg K) and $k=0.13$ W/(m K) give
$\rho c_p=1.874\times10^6$ J/(m3 K) for the user's density. Reservoir inflow
and the pinned groove carry $T_{\rm in}$. These are cooling scenarios, not
identified experimental boundary conditions. Since $\mu(T)$ is not supplied,
this added temperature field has no viscosity feedback on pressure or position
and must not be described as a coupled THD result.


## FW68D / R410A digitized solution laws (2026-09-15)

The [source table and calibration](data/R410A_FW68D/README.md) preserve 125 chart
points, including 15 pure-oil points and 11 points at each refrigerant mass fraction
from 0.1 through 1.0. In these laws,

$$c_d = \frac{m_{\mathrm{R410A,dissolved}}}{m_{\mathrm{oil}}+m_{\mathrm{R410A,dissolved}}},
\qquad T_C=T_K-273.15,\qquad p=p_{\mathrm{absolute}}.$$

The CSV stores temperature in degrees Celsius, pressure in MPa absolute and kinematic
viscosity in mm²/s. Loading converts these to K, Pa and m²/s. No atmospheric pressure
is added. The user's explicit mass and pressure basis overrides the unspecified basis
metadata in the original digitization.

### Actual-composition kinematic viscosity

`liquid_viscosity_model = KINEMATIC_TABLE` reads constant-composition curves with
separate temperature nodes, without rectangular-grid resampling. Let $g(\nu)=\ln\nu$
for `kinematic_viscosity_log_interpolation = true`, or $g(\nu)=\nu$ for `false`.
On each enclosing concentration curve $i$, its enclosing temperatures define

$$f_i=\frac{T-T_{i,j}}{T_{i,j+1}-T_{i,j}},\qquad
v_i(T)=(1-f_i)g(\nu_{i,j})+f_i g(\nu_{i,j+1}).$$

With $w=(c_d-c_i)/(c_{i+1}-c_i)$, the model is

$$\nu(c_d,T)=g^{-1}\!\left((1-w)v_i(T)+wv_{i+1}(T)\right),\qquad
\mu_l(c_d,p,T)=\rho_l(c_d,p,T)\,\nu(c_d,T).$$

Exact nodes return their own curve value. The logarithmic option preserves positive
viscosity and interpolates ratios. The source chart supplies equilibrium points;
using these curves at off-equilibrium composition assumes negligible independent
pressure dependence at fixed composition and temperature. No Barus factor is applied.
`viscosity_pressure_coeff` must be zero for this model.

### Equilibrium solubility and its inverse

`solubility_model = CURVE_TABLE` groups the same nonzero-concentration rows by their
labelled isobar. Within an isobar, it linearly interpolates $c_{\mathrm{sat}}$ in
physical temperature. Between enclosing absolute pressures,

$$c_{\mathrm{sat}}(p,T)=(1-q)c_k(T)+q c_{k+1}(T),\qquad
q=\frac{p-p_k}{p_{k+1}-p_k}.$$

The bubble pressure helper inverts this same piecewise-linear isotherm and returns
its first crossing. Its zero-charge convention remains $p_{\mathrm{sat}}(0,T)=0$.
Pure-oil viscosity rows have no pressure data and create no artificial solubility
anchors. The equilibrium concentration is distinct from the transported actual
concentration used for viscosity.

`solution_table_clamp = true` holds each contributing curve's endpoint when its
query lies outside the digitized temperature range; it also holds the outer
pressure or concentration endpoint. This policy preserves the 100% saturation
limit but also creates a 10% lower concentration limit where the solubility chart
ends. It must not be read as measured behavior in the dilute region. With `false`,
queries outside either contributing curve fail. No linear extrapolation is used
for runtime properties. The pressure axis covers 0.1–4 MPa; valid temperature
intervals depend on composition or isobar and are listed in the data README.

### Requested density closures

The configurable density choices are

$$\rho_l=\rho_0,\quad \rho_0=940\;\mathrm{kg\,m^{-3}}\quad
\text{for CONSTANT},$$

$$\rho_l(T)=a+bT_C,\quad a=954.9\;\mathrm{kg\,m^{-3}},\quad
b=-0.7284\;\mathrm{kg\,m^{-3}\,K^{-1}}\quad
\text{for LINEAR\_TEMPERATURE}.$$

`rho`, `density_temperature_intercept` and `density_temperature_slope` expose these
parameters. At 40 °C the linear law gives 925.764 kg/m³. Nonpositive or nonfinite
linear density is rejected; it is not silently floored. The coupled solver uses
$1/\rho_l(T)$ for both liquid partial specific volumes under this
composition-independent density closure, including its pressure-work contribution
at a feed. Existing heat capacities and caloric assumptions are unchanged.

Both density options are user-specified liquid-solution approximations with no
composition or pressure dependence. They do not supply a measured pure-R410A
liquid density. The pure-refrigerant chart endpoint describes liquid kinematic
viscosity; the free-gas EOS and `mu_gas` remain separate models. The top-level
`fluid_property_model = CONSTANT` continues to bypass variable property laws.

## Experimental bubble branch integration (2026-09-15)

The separate `pancake_bubble` target retains the stereographic hemisphere and
finite-volume thickness-transport implementation from `feature/bubble`. It does
not change the Reynolds, cavitation, thermal or refrigerant/oil models above.
Only its geometry, mask, surface operators and gravity thickness transport are
implemented. Its original full-model specification is preserved in
[experimental/bubble/PHYSICS.md](experimental/bubble/PHYSICS.md), with an explicit
implementation-status notice. The remaining momentum, thermal, evaporation,
Marangoni, disjoining-pressure and rupture models are unfinished.

## Solver audit corrections (2026-09-15)

### Periodic inlet geometry and prescribed values

All bearing paths now use one inlet membership calculation. For a cell angle
$\theta_i$ and inlet angle $\theta_s$, the angular displacement is

$$\Delta\theta=\operatorname{remainder}(\theta_i-\theta_s,2\pi).$$

A groove includes a cell when $|\Delta\theta|$ is at most half its angular width.
A circular inlet uses $R^2\Delta\theta^2+(z_j-z_s)^2\le r_s^2$.
Equivalent full-revolution shifts therefore give the same pressure, thermal
and diagnostic masks, including neighbour cells across MPI partitions.
The last matching pressure reservoir retains precedence when inlets overlap.

For a prescribed pressure, Elrod content or supply temperature $\phi_s$, replace
the equation $a_P\phi_P=\sum_N a_N\phi_N+b$ by

$$s\phi_P=s\phi_s,\qquad s=\max(|a_P|,10^{-30}),\qquad a_N=0.$$

Only that row is replaced. Neighbour equations retain their coupling to the
prescribed cell. This is the existing gas-transport boundary method; the previous
$10^{20}$ pressure/temperature penalty is removed. Inlet mass supply remains the
imbalance of the displaced physical balance, as described in the diagnostics
sections. A prescribed reservoir row is not a closed-cell conservation equation.

### Liquid supply enthalpy with linear density

The coupled internal supply term remains
$\sum_k\Delta M_{k,s}[u_k(T_s)+p_s v_k(T_s)]$ for inflow and uses the local
film state for outflow. With composition-independent linear density,

$$v_o(T_s)=v_d(T_s)=\frac{1}{954.9-0.7284(T_s-273.15)}.$$

The configurable coefficients replace the numbers above when supplied. Previously,
the pressure-work contribution used the receiving cell's liquid volume even
when the injected fluid had a different temperature. It now uses $T_s$.
A uniform pressure-fed cooling test checks the independent energy balance and
the additional liquid mass drawn into the fixed gap.

Temperature-dependent volume combined with the existing constant-partial-heat-
capacity internal energies remains a caloric approximation. A complete
thermodynamically consistent thermal-expansion model is not established by this
boundary correction. That choice is recorded for discussion in
[the audit](docs/SOLVER_AUDIT_2026_09_15.md).

### Computation and stopping controls

The coupled solver refills its cell-state buffer for each residual evaluation.
Once the final residual satisfies convergence, the same states supply physical
acceptance checks and output fluxes. Reuse introduces no interpolation or lag
between Newton trials. EOS root bracketing also reuses the pressure derivative
already calculated with pressure; its branch selection and tolerance are unchanged.

PETSc runtime stopping controls override per-solve defaults. The coupled
accepted-state infinity norm must still be below $10^{-10}$, with the existing
positivity, table-domain and free-gas phase checks. Looser stopping options cannot
turn a rejected state into an accepted solution. Failed solves now report their
actual iteration count. Every MPI partition must own at least as many physical
theta cells as the requested ghost depth before an exchange begins.


## Steady experimental load matching (2026-09-27)

The [six FW68D cases](cases/02d_fw68d_steady_six/README.md) select
`STEADY_STATE` and `STATIC` for each fixed-geometry pressure/energy/CHT solve.
Physical time-storage terms and bearing-motion integration are inactive.
The equilibrium bearing offset is obtained by an external search over
$q=(\varepsilon,\phi)$, using the total fluid force on the bearing:

$$\mathbf r(q)=\mathbf F_{\mathrm{fluid}}(q)+\mathbf W,
\qquad \frac{\|\mathbf r\|_2}{\|\mathbf W\|_2}<0.002.$$

Every force evaluation must first converge its pressure and temperature
coupling; a capped or restored-best inner iterate is not an equilibrium.
The applied load vector and feed positions must share a specified frame.
The current `load_angle_deg` input does not itself set this external target.

TC0 defines feed/reservoir and external-fluid temperature references, not a
uniform film or wall temperature. For convection it supplies $T_\infty$ in
$q''=h_{\mathrm{ext}}(T_{\mathrm{wall}}-T_\infty)$; the heat-transfer
coefficient and solid conduction geometry remain independent inputs.
The supplied viscosity reference fixes the experimental Sommerfeld coordinate;
local film viscosity remains temperature dependent in the steady solve.

This is a campaign specification, not a solver change or a convergence result.
Feed and thermal boundaries remain unresolved. Quantitative torque comparison
also awaits correction of the mean-viscosity wall-stress post-processing to
match the existing cross-film flow integrals; see the case plan.


### FW68D groove and external natural convection update (2026-09-27)

The original interpretation in this section used load (0,-W) with groove
theta=90 degrees. It was superseded on September 28 after the user clarified
a push from below: applied force (0,+W), groove theta=90 degrees, and target
fluid-on-bearing force (0,-W). The groove is opposite the actuator contact
location, not opposite the force vector. Its 5 mm bore arc converts
to 13.2983125 degrees; the native groove spans all axial film cells. Its
reservoir pressure and temperature match the bearing ends, with TC0 supplying
the temperature. The requested cavity plateau is 0 Pa, pending its datum;
no conversion to absolute pressure has been assumed.

The housing is stainless steel, outer radius 42.5 mm, length 40 mm. The user
selected surrounding still oil at TC0. Adopt an estimated horizontal-cylinder
mean natural-convection coefficient at the oil temperature halfway between
the mean outer wall and TC0 (Churchill and Chu,
[*Int. J. Heat Mass Transfer*](https://doi.org/10.1016/0017-9310(75)90222-7), 1975).
See [the case boundary plan](cases/02d_fw68d_steady_six/README.md) for the
correlation, conditional values and unspecified material properties.

The r-theta solid model has no axial end-face heat loss. With L/D=0.47 and
annular end area about 79% of lateral area, the cylinder correlation is an
approximate lateral boundary, not a full finite-housing thermal model.
Existing CHT accepts fixed h; self-consistent updating of h during equilibrium
remains planned, not implemented or validated by these input records.

## Steady FW68D thermal closure and wall stress (2026-09-27)

The groove and both axial reservoirs are 101325 Pa absolute at the selected
case TC0. The requested cavity plateau is interpreted as 0 Pa absolute.
Validation permits this only for gas-free SINGLE_PHASE FBNS with SCALAR_PCAV,
no gas-pressure coupling, zero initial dissolved charge/Henry coefficient and
zero vapor-temperature coefficient. Negative pressure and zero pressure in
other gas/coupled configurations remain invalid. FBNS uses the pressure offset
p-p_cav and a separately positive reference scale; it does not divide by p_cav.

### Wall forces with cross-film viscosity

For shaft surface speed U=omega R and y measured from shaft to bush, define
$F_0=\int_0^h\mu(y)^{-1}\,dy$ and $F_1=\int_0^h y\mu(y)^{-1}\,dy$.
The positive resisting shaft stress and fluid-on-bush tangential stress are

$$\tau_{s}=\frac{U}{F_0}+\frac{F_1}{F_0}\frac{\partial p}{\partial s},
\qquad \tau_b=\frac{U}{F_0}-\left(h-\frac{F_1}{F_0}\right)
\frac{\partial p}{\partial s}.$$

The bush axial stress is
$\tau_{b,z}=-(h-F_1/F_0)\partial p/\partial z$.
Apply the full-film gate to pressure-driven stresses and the existing STRIATED
liquid-filling factor to Couette stress in cavities. With cross-film closure
inactive, retain the classical mu U/h and h/2 expressions. For uniform mu,
F0=h/mu and F1=h^2/(2mu), so the new formula has the same limit. The correction
makes reported torque and wall forces use the same viscosity integrals as flow.

### Adopted materials and exposed shaft ends

The user authorizes ordinary or similar thermophysical properties. Use BVC32
PVE refrigeration oil as the cp/k analogue: cp=2010 J/(kg K), k=0.13 W/(m K),
the manufacturer's values at 40 C, held constant
([BITZER technical data](https://www.bitzer.de/shared_media/html/est-500/en-GB/271261579271263371.html)).
Keep the supplied FW68D nu(T) and rho(T), so the experimental reference S does
not change. The film equation uses constant volumetric rho_cp=rho(TC0)*2010
within each case; external-film property evaluation uses cp=2010. This is a
caloric approximation, not a measured FW68D heat-capacity curve.

Use 304 stainless k=15 W/(m K), held constant, as the housing analogue
([Outokumpu physical properties](https://www.outokumpu.com/en/products/product-ranges/-/media/files/products/core/outokumpu-core-range-datasheet.pdf)).
Steady solids have no storage term, so solid heat capacity/density are not
required by this reduced CHT model.

Both shaft ends are exposed to surrounding oil at TC0. Their combined area is
$A_e=2\pi R^2=0.0029044024$ m2, giving native shaft_ambient_area_ratio=R/L=0.5375.
The shaft remains one equilibrated temperature; axial solid conduction is not
resolved. For rotating ends in otherwise still oil, use the laminar high-Pr
free-disk estimate

$$Re_R=\frac{|\omega|R^2}{\nu},\qquad
Nu_R\simeq0.62 Re_R^{1/2}Pr^{1/3},\qquad h_e=\frac{k_o Nu_R}{R},$$

with properties at (Tshaft+TC0)/2. This high-Pr limit follows rotating-disk
heat/mass-transfer theory (Liu and Stewart,
[*Int. J. Heat Mass Transfer*](https://doi.org/10.1016/0017-9310(72)90181-0), 1972).
The heat sink is Q_e=h_e A_e(Tshaft-TC0). This coefficient is distinct from
natural convection at the stationary housing. Finite end geometry, confinement,
shaft conduction and mixed buoyancy/rotation are model uncertainties, not
measured corrections. Record disk Re/Pr with the results.

### External coefficient iteration

The campaign driver wraps the existing steady solver with updates of housing
and shaft coefficients. Each accepted inner evaluation must meet pressure and
temperature coupling criteria. For the current constant-k r-theta housing with
uniform outer Robin coefficient, angular averaging gives

$$\overline T_{outer}=T_\infty+
\frac{\overline T_{bore}-T_\infty}
{1+h_b r_o\ln(r_o/R)/k_b}.$$

This reconstructs the mean outer surface temperature without treating the bore
as the outside wall. Evaluate the housing Churchill-Chu coefficient at
(overline(Touter)+TC0)/2 and the end-disk coefficient at (Tshaft+TC0)/2; iterate
with relaxation. Require relative coefficient change below 1e-3 and wall
changes below 0.01 K before accepting the external thermal closure. Force
balance then requires fluid-on-bearing (0,+W) for applied (0,-W).
The housing has no axial end-face cooling in this r-theta model; this omission
and the short L/D=0.47 remain explicit limits on the thermal comparison.

### FW68D preflight qualifications (2026-09-27)

The initial F convergence plot exposed a tolerance mismatch: FBNS could stop
at its initial iterate once the scaled residual passed its default 1e-9,
leaving pressure unchanged and making the outer pressure-change norm zero.
The prior pressure update still exceeded the campaign 1e-6 outer tolerance.
Set campaign `fbns_tol=1e-12` so the outer criterion is not hidden by that floor.
This is a numerical setting change, with no change to the physical equations.

The PINNED groove replaces mean-film energy rows by TC0, while cross-film
reconstruction still treats those cells as sheared film with conjugate walls.
The pinned row therefore acts as an implicit supply-energy sink. A warm shaft
and fixed mean can produce an algebraically sub-ambient bore temperature,
which bush conduction spreads into neighboring cells. The F diagnostic found
a 0.218 K groove undershoot and 0.069 K in two adjacent columns. Balanced
solid wall fluxes do not remove this local modeling inconsistency. Preserve
these values as baseline diagnostics; do not interpret them as physical
sub-ambient cooling or assume negligible effects on other outputs. The correction below derives reservoir exchange from the solved mass flux;
it does not prescribe an unmeasured mixing-flow rate.


## Conservative steady FW68D closure repair (2026-09-28)

The September 27 baseline failed for two separate reasons: undamped quartic
collocation updates stalled the pressure/temperature fixed point, and a
prescribed groove mean temperature could force an artificial cold part of the
cross-film profile against warmer walls. Fixed-position convergence is only
an inner calculation. The six experimental conditions require an outer search
for eccentricity and attitude satisfying the two fluid-force components.
The following changes supersede the PINNED/constant-volumetric-capacity
assumptions in the initial campaign record; the earlier failures are retained.

### One discrete liquid flux for mass and enthalpy

For FBNS, energy now uses the public Reynolds mass-flux operator on every
face, including its harmonic pressure mobility and upwind liquid fraction
on the Couette term. In schematic form,

\[
\dot m_{i+1/2}=G_{i+1/2}(p_i-p_{i+1})
                 +C_{i+1/2}\theta_{\rm upwind},\qquad
\dot H_f=c_p\dot m_f(T_{\rm upwind}-T_0).
\]

Pressure-driven flux is not multiplied by an averaged cavity fraction. FBNS
complementarity makes pressure constant inside a cavity, while a physical
pressure flux can remain at a rupture/reformation front. Pressure-work
heating therefore remains active at those fronts. The dissipation law is
\(\theta U^2/F_0+F_2|\nabla p|^2\) for a striated cavity, with the full-film
Couette fraction equal to one. The existing thin-film dissipation and
cross-film approximations remain; no vapor latent-heat model is added.

A positive `liquid_specific_heat` selects constant mass-specific heat and
updates `rho_cp_liquid_solution = rho_liquid_solution * liquid_specific_heat`.
FW68D uses 2010 J/(kg K). Zero retains the earlier volumetric-capacity input
for configurations that have not selected this option.

### Pressure-fed groove energy

`groove_mixing_model = FLUX_BALANCED` leaves the groove pressure constraint
in place and solves the groove film temperature. The reservoir exchange is
computed from the divergence of the same signed face mass fluxes:

\[
q_g=\sum_{f\in\partial g}\dot m_f,\qquad
Q_{g,\mathrm{in}}=c_p\left[q_g^+(T_s-T_0)+q_g^-(T_g-T_0)\right].
\]

Here \(q_g^+=\max(q_g,0)\), \(q_g^-=\min(q_g,0)\). Net supply enters at
TC0; drainage leaves at the solved groove temperature. The groove is an
ideal pressure-fed reservoir coupled to the existing film cell, not a
resolved deep-groove recirculation region. No independent fresh-oil flow is
invented. At a full-length groove/end overlap, the groove row owns the
reservoir exchange; the discarded axial pressure face is not counted again.
This closure is currently restricted to steady, gas-free, single-phase FBNS,
upwind thermal convection, and positive constant mass-specific heat.

### Conservation and admissibility

The global energy ledger sums the same discrete fluxes and sources:

\[
R_E=H_{\rm axial,out}+Q_{\rm walls,out}-Q_{\rm viscous}
       -Q_{\rm groove,in}-Q_{\rm mass\ correction}.
\]

Periodic/internal face enthalpy fluxes cancel. The final term records the
existing non-groove continuity-defect correction explicitly; it cannot hide
a mass error because mass closure is tested separately. Wall heat is evaluated
from the actual CHT flux, so its lag is included in the nonlinear residual.
The solver keeps iterating after small pressure/temperature changes until
\(|R_E|/\max(S_E,1\,\mathrm W)<\texttt{coupling_tolerance}\), where \(S_E\)
is the sum of absolute ledger contributions. The campaign requires 1e-6.

The liquid rate residual is storage minus net boundary/feed inflow and clamp
rate, in kg/s. Its scale includes absolute reservoir exchanges, storage and
clamp rates. This dimensionless check is independent of inactive steady
`dt` metadata; the campaign requires absolute relative residual at most 1e-8.
The legacy timestep-normalized diagnostic remains for compatibility.

The raw quartic profile is checked at both walls and every stationary point
inside the film. Its global minimum/maximum are written to diagnostics.
Nonfinite/nonpositive walls, profiles, viscosity integrals, or a failed bush
conduction solve terminate the run collectively. No temperature floor or
viscosity-ratio clamp repairs the result. For this common-ambient, positive-
heating campaign, every film/wall/profile temperature must be at least
TC0 minus 0.001 K numerical allowance. This is a rejection criterion, not a
clipping operation or a general constraint on other boundary conditions.

### Coupled equilibrium and reusable initial fields

The raw collocation fixed-point defect is checked separately as
max |X_target-X_previous| / max(|T_reference|,1 K) below 1e-6. Testing this
unrelaxed defect prevents source damping from hiding a changing internal
profile. Exact profile extrema must remain inside the FW68D viscosity table.

The quartic collocation update uses `cross_film_relaxation = 0.3` without
changing mesh, profile order, quadrature or residual tolerances. Each outer
force trial closes pressure, film/solid temperature and both external heat-
transfer coefficients. The position solver bounds Newton steps and backtracks
when a trial fails or increases the vector force error. Only
\(\|\mathbf F_{\rm fluid}+(0,+W)\|/W<0.002\), together with all thermal and
physical checks, is an accepted bearing position.

A successful gas-free steady FBNS solve writes `steady_state.dat` with
pressure, liquid content, mean temperature and quartic collocation state.
`steady_initial_state_file` loads those values only as an initial guess.
The current geometry, load and boundary configuration still define the new
problem. Grid dimensions, R, c, L, field order and values are validated; MPI
rank decomposition may differ. The campaign uses only same-case states with
matching source/executable/input provenance and rechecks their diagnostics.
Neither a saved state nor a cached fixed-position result establishes load
balance. Failed attempts and their files remain available for diagnosis.


### FW68D upward push and force-vector correction (2026-09-28)

With the journal centre fixed and the actuator below the bearing, the applied
force on the housing is W_app=(0,+W). Fluid-on-bearing equilibrium is
F_fluid+W_app=0, so the two-component residual is (Fx,Fy+W), normalized by W.
The feed groove remains at theta=90 degrees and the contact location is
theta=270 degrees. Setting load_angle_deg=90 alone has no effect on the
static inner solve; the explicit external-load components and outer force
residual must agree. This correction does not change journal RPM.

For the static film law h=c-e*cos(theta-psi), the bearing centre relative to
the journal is b=-e*(cos(psi),sin(psi)). Upward-bearing equilibria therefore
require searching the opposite attitude half-plane from the invalidated
downward-load campaign. The original E temperature peak near the top groove
belonged to the wrong load configuration; it is not evidence for the corrected
experiment. All previous case outputs are retained, but their numerical
convergence is insufficient for physical acceptance. Fresh inputs/results and
provenance prevent reuse as corrected equilibria.


### Steady pressure-update acceptance (2026-09-28)

Every steady FBNS call takes at least one Newton correction before accepting
its mass/complementarity residual. Previously a warm state already below
fbns_tol could return at iteration zero, manufacture a zero outer pressure
change, and pass the pressure-update gate without measuring a correction.
The minimum correction makes the outer pressure change observable, including
at very small residuals; all existing thermal, conservation and force gates
remain required. An exact zero-RHS state takes a safe zero correction.
Transient FBNS retains its existing initial-residual exit.


## 2026-09-28 - Remove the absolute-temperature linear accuracy floor

The corrected-load D run hit4000 sweeps with Rp8.62e-15 and RT1.51e-15,
but energy residual7.224642e-6W /6.406604W=1.12769e-6 still failed the
unchanged1e-6 gate. The linear energy equation used absolute Kelvin with
relative tolerance1e-8, so a large common300K offset obscured a small heat
budget. Solve the identical assembled equation for temperature rise instead:
A*(T-Tref)=b-A*Tref, including actual axial stencil and Dirichlet rows.
Restore absolute temperature after the solve; physics and tolerances are unchanged.

A low-heat regression reproduces the old four-rank residual1.33436e-6W;
the shifted formulation gives9.64e-14W, with the same rise when the reference
moves from300 to900K. One/four-rank energy tests and28 affected native CTests
pass. Log R_E alongside field residuals once the energy gate is evaluated,
so stationary field updates cannot conceal a failed energy balance.

This numerical correction is distinct from A's21.7mK profile undershoot at
the groove entrance; that model-consistency issue remains under investigation.


## 2026-09-28 - Position coupled into the steady bearing solve

The requested steady equilibrium has no physical time derivatives. Bearing
position remains an unknown; setting all storage terms to zero does not fix
that position. The earlier external force controller solved the same algebraic
force condition, but fully converged temperature for each trial geometry. That
execution path is superseded for this campaign by position corrections inside
the pressure/cavitation/energy coupling iteration.

Let the bearing-centre displacement from the fixed journal centre be
\(\mathbf b=(b_x,b_y)\). With the existing tilt convention,

\[
h(\theta,z)=c+b_x\cos\theta+b_y\sin\theta
 -(s_x\cos\theta+s_y\sin\theta)(z-L/2),\qquad \partial h/\partial t=0.
\]

At every steady coupling iteration, solve the steady Reynolds/cavitation
system with current temperature-dependent coefficients, integrate the total
fluid force on the bearing, and form

\[
\mathbf r_F(\mathbf b)=\mathbf F_{\rm fluid}(p,\phi,T,\mathbf b)
                     +\mathbf W_{\rm applied}=\mathbf0.
\]

For this rig \(\mathbf W_{\rm applied}=(0,+W)\). Both Cartesian force
components are required. The pressure force and viscous surface force use the
same geometry and pressure as the accepted pressure solve.

A two-variable finite-difference Newton correction acts on \(b_x/c,b_y/c\).
Each derivative probe freezes the current thermal state and solves only
pressure/cavitation. For a frozen dimensionless cross-film viscosity profile,
changing the local gap by the ratio \(r=h_{\rm trial}/h_{\rm base}\) scales
\(F_0,F_1,F_2\) by \(r,r^2,r^3\), respectively, and the effective Couette
gap by \(r\). Keeping those coefficients unchanged would give an incorrect
force Jacobian. Restore the same base fields before each probe, bound the
position correction, and backtrack until the trial has positive clearance,
converged pressure/cavitation and a reduced force norm.

Continue the energy/CHT update at the accepted geometry and repeat the complete
steady coupling iteration. This is a nonlinear algebraic equilibrium solve;
no bearing mass, acceleration, physical velocity, pseudo-time storage or
squeeze-film source is introduced by a position correction. The groove remains
at its prescribed angular location.

The same final iterate must satisfy force balance, pressure and temperature
change tolerances, cross-film profile convergence and the discrete global
energy balance. The campaign separately verifies mass balance, admissible
fields, and exterior heat-transfer-coefficient closure. A low pressure or
temperature change with a high force imbalance cannot declare convergence.
Final bearing position and its initial guess are recorded separately. Failed
or capped iterations remain diagnostic evidence, never accepted equilibria.

The earlier groove-profile undershoot is a separate limitation. Position
coupling does not justify clipping temperature or weakening its acceptance
criterion. The unsuccessful reservoir-collocation experiment is preserved
outside the source snapshot and excluded from this implementation.

### September 28 correction: resolved film energy and exposed housing ends

The FW68D campaign now selects `cross_film_thermal_model = FINITE_VOLUME`.
The old `POLYNOMIAL` option is retained for reproducibility. Its restricted
quartic profile and convection limiter can produce subambient interior
temperatures with positive dissipation and equal-temperature reservoirs.
Changing convergence tolerances cannot repair that approximation.

Divide the film into N equal layers, with cell-centred temperatures T_k and
piecewise constant viscosity mu_k. Integrate the inverse viscosity exactly:

    F0 = integral(1/mu dy), F1 = integral(y/mu dy),
    F2 = integral(y^2/mu dy) - F1^2/F0.

The Reynolds mobility remains F2 and the Couette transport thickness remains
2 F1/F0. Layer face fluxes follow the corresponding integrated velocity
profiles and sum to the same conservative FBNS mass flux used by pressure.
Internal transverse fluxes follow layer continuity. For each layer solve

    sum_faces(cp * mass_flux_out * T_upwind)
      + sum_neighbours(G * (T_k - T_neighbour)) = viscous_heat_k.

The groove supplies positive net inflow at TC0 and drains negative net flow
at the solved layer temperature. Axial inflow uses TC0. A continuity-defect
correction preserves constant temperature under the residual of the pressure
solve; it is included explicitly in the global energy ledger. Dissipation
is the exact layer integral of the squared Couette and pressure-driven shear,
not a clipped or fitted heat source. Its film integral is
fill * U^2/F0 + F2 * |grad(p)|^2. The layer temperatures, bore temperatures,
housing nodes and one shaft temperature are solved in one sparse system.
The unknown is T minus the reference temperature, avoiding cancellation
against a large absolute-Kelvin baseline. Viscosity is updated on the next
outer iteration. No physical time derivative is present.

`bush_thermal_model = CONDUCTION_3D` resolves radial, circumferential and axial
housing conduction. `bush_axial_cells` divides the fluid axial cell count;
film heat fluxes are conservatively summed onto each solid axial slab.
Radial links use logarithmic cylindrical resistances. Both annular ends have
convection to TC0 in series with the half-cell conduction resistance.
`bush_end_heat_transfer` specifies a fixed coefficient for verification;
`bush_end_natural_convection = true` instead updates it locally from the
solved end-surface temperature. For the horizontal shaft, the annular ends
are vertical. Use a vertical-plate analogue of height D_outer:

    Nu = [0.825 + 0.387 Ra^(1/6) /
          (1 + (0.492/Pr)^(9/16))^(8/27)]^2, h = k_oil Nu / D_outer.

Oil properties are evaluated at the mean of surface and ambient temperature.
This is an approximation to finite annular geometry, not measured cooling
(Churchill and Chu, *[Int. J. Heat Mass Transf.](https://doi.org/10.1016/0017-9310(75)90243-4)*, 1975).
The outer cylinder retains its separate horizontal-cylinder correlation.
Its coefficient iteration now uses `temperature_housing_outer`, the actual
resolved surface temperature; the old axially uniform annulus conversion
from mean bore temperature is invalid with exposed ends. Housing fields
`cross_film_solid_k` retain each radial node at each axial slab, allowing
TC1-3 interpolation at r=26 mm, z=20 mm and theta=240,270,300 degrees.
The solid inner radius retains the thin-film approximation R, 42.5 um inside
the physical bore. The shaft retains a global temperature and two rotating
end-disk convective sinks; its axial conduction is not spatially resolved.

`steady_pressure_inner_only = true` uses the requested inner FBNS residual
of 1e-6 without an additional pressure-field-change threshold. Force,
temperature, all layer/solid changes and global energy must still close.
The temperature floor remains TC0 minus 0.001 K; there is no temperature
clipping and no calibration to measured thermocouples or torque.
`steady_output_interval` saves complete, iteration-numbered checkpoints.
A STOP file requests export after the current iteration and returns failure,
never a falsely converged result. Version-2 initial-state files preserve
all film layers; version-1 files seed them uniformly from the saved mean.

Validation includes exact Couette conjugate temperatures, second-order
through-film refinement (8/16/32 layers), signed groove flow and reverse
rotation, exact outer-plus-end housing heat loss, natural-convection fixed
point closure, MPI rank comparisons and stopped-run recovery.

### September 29 - Pressure-work sensitivity not retained

A single Case B sensitivity added reversible liquid compression and
decompression work, beta*T*u.grad(p), to the steady energy equation.
It changed TC3 by -0.01336 K and did not reproduce the measured 1.334 K
TC2-to-TC3 drop. The user rejected the added complexity, so the option was
removed and the validated viscous-heating/CHT implementation restored.
The source snapshot, converged results and full-circle comparison remain
in `results/fw68d_pressure_work_20260929`. Reproducing that trial requires
its archived source; it is not part of the current solver model.

## Homogeneous thermal-expansion sensitivity (2026-10-17)

The optional external geometry iteration uses the solved isothermal shaft
temperature, $T_s$, and the volume-mean monolithic bearing/housing temperature,
$T_b$. Both bodies use the specified 304 analogue, $\alpha=16\times10^{-6}$
K$^{-1}$. The reference temperature is $T_0=20$ C, with
$D_{s0}=0.043$ m, $D_{b0}=0.043085$ m, $L_0=0.04$ m and $D_{o0}=0.085$ m:

$$
D_s=D_{s0}[1+\alpha(T_s-T_0)],\qquad
D_b=D_{b0}[1+\alpha(T_b-T_0)],\qquad
c=\frac{D_b-D_s}{2},
$$

$$
L=L_0[1+\alpha(T_b-T_0)],\qquad
D_o=D_{o0}[1+\alpha(T_b-T_0)].
$$

Temperature differences are in kelvin. Every update starts from these original
reference dimensions; it never multiplies an already-expanded geometry by
another expansion factor. Equal shaft and bearing warming therefore gives
$c=c_0[1+\alpha(T-T_0)]$ rather than a spurious change caused by accumulated strain.

For the actual solid mesh of a solve, radial edges run from the shaft radius
$R=D_s/2$ to $R_o=D_o/2$. With equal angular and axial cell sizes,

$$
w_k=\frac{r_{k+1}^2-r_k^2}{R_o^2-R^2},\qquad
T_b=\sum_{k=0}^{15}w_k\left(\frac{1}{N_\theta N_{z,s}}
\sum_{i=0}^{N_\theta-1}\sum_{j=0}^{N_{z,s}-1}T_{kij}\right).
$$

Here $N_\theta=720$ and $N_{z,s}=16$. Each native solid axial slab is repeated
into 12 of the 192 film output rows. Postprocessing verifies this repetition
and uses all 16 radial fields; bore temperature and TC0 do not replace the
volume mean. The conduction inner radius remains $R$, while the physical bore
is $R+c$. This preserves the existing thin-film approximation during the test.

Each geometry receives a fresh coupled force, cavitation, energy and exterior
cooling solve. Only cooling iterations at the same geometry reuse a validated
checkpoint. The native configuration updates $R,c,L,R_o$, so the surface speed
$\omega R$, gap field, torque integral, conduction mesh and exposed areas follow
the new dimensions. Exterior cooling uses $D_o$ for the housing and $R$ for
the rotating disks. Two shaft disks have area $2\pi R^2$, giving the existing
normalized cooling-area ratio $R/L$. The groove keeps its original angular
width, 90-degree centre and full-length convention. Sensor radii scale as
$r_{TC}=0.026(L/L_0)$ m, at $z=L/2$ and 240/270/300 degrees. The shaft centre
remains fixed; load-balanced bearing displacement adds no translational velocity.

After closure, the new $T_s,T_b$ propose the next dimensions. Feedback converges
when $|\Delta c|\leq10^{-9}$ m and each of
$|\Delta D_s|,|\Delta D_b|,|\Delta L|,|\Delta D_o|\leq10^{-8}$ m. At most five
geometry updates are allowed; reaching that cap without closure is reported
as unresolved. These are numerical stopping criteria, not measurement errors.
Fresh nominal and expanded repeats tighten force tolerance from 0.002 to
$2\times10^{-5}$, FBNS/thermal criteria from $10^{-6}$ to $10^{-8}$, exterior
coefficient closure from $10^{-3}$ to $10^{-5}$ and surface-temperature change
from 0.01 to 0.0001 K. Paired responses separate expansion from stopping error.

The model excludes stress, restraint, nonuniform thermoelastic distortion,
elastic load deformation, contact and different stainless grades. It does not
change measured temperatures, oil laws or the corrected experimental torque.
See [Thermal-expansion campaign](cases/02d_fw68d_steady_six/THERMAL_EXPANSION.md)
for the implementation trace, evidence and results.

## Conditional undersaturated FW68D–R410A scenario (18 October 2026)

The 17 selected mixture cases use the same nominal geometry and resolved steady
Reynolds/FBNS/film-energy/housing-conduction equations as the accepted pure-oil
campaign. The measured refrigerant content of approximately 30% is not assigned
a mass or volume basis. The digitized chart's 30% label identifies one
fixed viscosity curve, with no concentration evolution or phase-release model.

For adjacent decoded nodes $(T_i,\nu_i)$ and $(T_{i+1},\nu_{i+1})$, with Celsius
temperature and kinematic viscosity in mm²/s, evaluate

$$w=\frac{T-T_i}{T_{i+1}-T_i},\qquad
\nu(T)=\exp[(1-w)\ln\nu_i+w\ln\nu_{i+1}],\qquad
\mu(T)=950\,\nu(T)\,10^{-6}\ {\rm Pa\,s}.$$

The calibrated chart coordinates are decoded in logarithmic Kelvin on the
temperature axis and logarithmic viscosity on the ordinate before this Celsius
interpolation. Native Kelvin queries have identical interpolation weights.
Local solved cross-film temperatures determine viscosity. The accepted node
range is −26.0429574434 to 84.0172809983 °C. Evaluation outside this range fails.
A viscosity-only table contains one solver coordinate $c=0$ and
blank pressure. This denotes no tracked species, not the physical absence of
refrigerant. It leaves solubility empty and cannot support CURVE_TABLE release.

The user confirms common chamber/groove/end pressure and liquid around exposed
surfaces. Each case therefore applies $p_{\rm res}=10^5p_{\rm measured,bar}$ Pa
absolute at both ends and the groove. The selected surrogate separately uses
$p_{\rm cav}=0$ Pa absolute as an oil-vapor-limit idealization while suppressing
refrigerant release. This is not mixture bubble pressure, refrigerant partial
pressure or a conclusion that local degassing cannot occur. Native R410A EOS
data alone cannot replace missing mixture solubility and thermal inputs.

The conditional bulk density is constant at 950 kg/m³. Fluid $c_p=2010$ J/(kg K)
and $k=0.13$ W/(m K) use the BVC32-PVE analogue. Steel $k_s=15$ W/(m K) uses
the 304 stainless analogue. Neither fluid value is a measured mixture property.
Exterior natural convection separately uses the declared Boussinesq analogue
$\beta=0.7284/950$ K⁻¹. Film density remains 950 kg/m³, and volumetric heat
capacity is $950c_p$ J/(m³ K). Buoyancy affects only exterior convection
correlations:

$$\alpha=\frac{k}{\rho c_p},\qquad Pr=\frac{\nu}{\alpha},\qquad
Ra_\ell=\frac{g\beta|T_s-T_0|\ell^3}{\nu\alpha}.$$

The cylinder and annular-face correlations use the same $\beta$, while rotating
shaft disks use $Re_R=|\omega|R^2/\nu$. All properties are evaluated at the local
surface/bath mean temperature. The existing annular implementation obtains
$\beta=-b/\rho$ from `density_temperature_slope`, denoted by $b$. Under CONSTANT
density, that slope is explicitly a buoyancy-only input. Validation requires a finite
negative slope and emits this assumption. LINEAR_TEMPERATURE behavior remains
unchanged. Both input paths feed the same existing heat-transfer equations.

TC0 supplies feed and external liquid temperatures, not prescribed wall
temperatures. The 17 experimental plotting coordinates retain measured TC2:

$$S_{\rm ref}=\frac{\mu(TC2)(RPM/60)LD}{W}\left(\frac Rc\right)^2,
\qquad f_{\rm CFD}=\frac{M_{\rm CFD}}{WR}.$$

Native torque is cylindrical journal shear integrated over theta and axial
position, including Couette and pressure shear with the existing striated
cavity treatment. Project it onto rotation for comparison. Exposed disk drag,
seals, motor losses and resolved groove sidewall drag are not included. Compare
against the saved bearing residual and separately compare $M_{\rm CFD}+M_{\rm
reference}(RPM)$ against measured directional motor torque. This reference's
transfer to pressurized mixture is unverified. Sample housing TC1/2/3 at 26 mm
radius, 20 mm axial position and 240/270/300°. TC0 has no independent sensor-error
test because it is the boundary proxy. Tight tolerances and the lowest-S film
layer refinement quantify numerical sensitivity without validating the physical
composition, cavity threshold or thermal analogues.

Normal cases use an FBNS pressure/complementarity tolerance of $10^{-10}$;
tight checks use $10^{-12}$. Normal force and temperature/profile limits remain
0.002 and $10^{-6}$; tight checks retain $2\times10^{-5}$ and $10^{-8}$.
The mass-balance magnitude limit remains $10^{-8}$, with normal/tight exterior
coefficient closure of $10^{-3}/10^{-5}$ and wall-temperature change of
0.01/0.0001 K. Global energy limits remain $10^{-6}$ normally and $10^{-8}$
for tight checks.

The initial M02 pilot had signed relative mass residual $-4.553754639\times10^{-8}$,
so it failed the independent mass gate despite pressure-equation residual
$G=6.023062\times10^{-9}$. A small scaled pressure residual alone does not
establish the required global mass balance. All three initial pilot attempts
remain preserved in `lge-mixture-20261018`. The retry uses the isolated
`lge-mixture-20261018-corrected` directory and the tighter FBNS settings above;
the earlier short pressure-refined attempt is also retained. Final cooling
verification uses the recorded seed of that call. Physical inputs, equations
and acceptance gates are unchanged.
That campaign subsequently completed all 17 cases and three focused checks;
the original source, inputs and results remain preserved separately.

# 2026-10-19 conditional LGE gas comparison and steady thermal extension

The separate 02f follow-on has launched full-grid remote pilots. It compares paired
steady COUPLED calculations with gas release disabled/enabled for the same
17 experimental operating windows. Both branches use the full chart
`nu(c,T)`, density law, resolved cross-film and housing energy model, geometry,
cooling, boundaries and load balance. `coupling_relaxation = 0.3` is common.
The control uses `SINGLE_PHASE` and no phase equilibrium; the gas-on case uses
`TWO_PHASE` with instantaneous local equilibrium. The only branch switches
are fluid-property model and the equilibrium flag, so their difference
isolates release within this conditional numerical model. The earlier 02e
completed campaign remains the distinct fixed-viscosity surrogate and is not
used as a matched gas-off control.

The source chart's pressure and concentration axes are unspecified. The
explicitly authorized scenario treats pressure as absolute and the
concentration coordinate as refrigerant mass fraction. Initial dissolved
composition is uniformly 0.30 across the initial field and both axial feeds;
incoming free-gas mass fraction is zero. Composition and phase fractions are
uniform through the gap, and equilibrium uses the film-mean temperature.
Layer temperatures remain separate energy unknowns. The solver computes equilibrium bubble pressure by inverting
the chart's first crossing and requires strict support;
unsupported chart or temperature queries fail without clamp or extrapolation.
Free gas must also remain within the stable real-R410A EOS domain at all film
layers and the journal/bearing walls. No fitted release or absorption rate is
used: gas-off suppresses phase transfer and gas-on assumes instantaneous
equilibrium.

Both branches retain the pure-oil proxy
\(\rho(T_C)=954.9-0.7284T_C\;\mathrm{kg\,m^{-3}}\), which came from a
different PVE/R32 lubricant context and is not an FW68D/R410A mixture-density
measurement (Zeng et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2021.12.009)*, 2022). Direct PVE/R410A literature distinguishes
actual concentration from saturation solubility; its oil grade and reported
scope do not establish FW68D coefficients (Fukuta et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2009.09.015)*, 2010). No verified
exact-grade density or heat-of-solution table was identified in the review.

The steady thermal extension alternates the coupled pressure/species
equilibrium solve with finite-volume energy solves across the film layers and
through the resolved housing network. It iterates local temperatures into
full `nu(c,T)` and phase properties, and accounts for component face enthalpy,
pressure work, dissipation, wall heat flux, feed enthalpy and solid conduction.
It is not a cross-film transient Newton solve. Caloric inputs are explicit
conditional analogues: oil \(c_p=2010\) J kg⁻¹ K⁻¹, liquid conductivity
0.13 W m⁻¹ K⁻¹, steel conductivity 15 W m⁻¹ K⁻¹, and dissolved-refrigerant
heat capacity (1656.9291 J kg⁻¹ K⁻¹), desorption internal energy
(192486.2475 J kg⁻¹), gas conductivity (0.0129583 W m⁻¹ K⁻¹), and gas
viscosity (1.26525×10⁻⁵ Pa s) from the [conditional caloric inputs](results/lge_mixture_density_gas_20261019/property_review/conditional_caloric_inputs.json).
These values are referenced at 293.15 K and 100 kPa. The dissolved-phase
analogue is pure saturated-liquid R410A at 20 °C with constant heat capacity
and zero excess mixing energy; it is not measured FW68D heat of solution.
Inputs are not fitted to torque or temperature.

The driver rejects nonzero return codes, unconverged pressure/coupling,
iteration-cap exhaustion, restored-best iterates, pressure residual at or
above \(10^{-10}\), mass residual above the tighter of the requested tolerance
and \(10^{-8}\), or volume, thermal, species, profile, global energy and
independent load residuals outside configured limits. It also rejects
nonfinite/nonphysical fields and any unsupported viscosity/chart or free-gas
EOS state. These gates establish accepted numerical solutions under the
scenario assumptions; they do not establish experimental accuracy. The
gas-on 24×8 pilot and matched 48×16 M16 relaxation-0.3 checks passed strict
raw-field acceptance. Three full-grid attempts have been interrupted and
preserved while the M02 initialization and position solve are diagnosed.
The full campaign is not complete and
no campaign-level gas-on/off results are available.

## Cavitating initial guess and load feedback

The optional `coupled_fbns_initialization` first solves gas-free FBNS pressure
and film content at the unchanged starting geometry, composition and
temperature. It supplies an initial guess only. Conservative areal inventories
are then reconstructed as
\[
m_o=(1-c_d)\rho_l h\theta_l,\qquad
m_d=c_d\rho_l h\theta_l,\qquad m_g=0,
\]
where \(c_d\) is dissolved refrigerant mass fraction. Current and pseudo-time
old inventories are initialized consistently. Every subsequent pressure,
phase, thermal and load solution must pass the original COUPLED gates. The
option defaults off and rejects configurations with prescribed free gas that
the preliminary gas-free solve would discard. A preliminary FBNS result is
never a completed COUPLED solution.

With phase equilibrium disabled, the nonlinear species rows are the dissolved
and free balances, respectively. Subtracting the dissolved balance from total
species conservation gives
\[
\partial_t m_g+\nabla_s\cdot F_g=\Gamma,
\]
including the same incoming and outgoing boundary sources. Solving this row
directly avoids cancellation between nearly equal dissolved and total-species
equations when free gas is absent. This is an equation transformation, with no
change to release physics. The original normalized total-species residual is
still included in final acceptance; pressure, volume and mass bounds remain
unchanged.

The resolved steady COUPLED iteration also damps the inverse-viscosity
profile. With layer index \(k\), its physical target is
\[
v_k^*=\left[\theta_l\mu_l(c_d,p,T_k)+(1-\theta_l)\mu_g\right]^{-1},
\qquad
v_k^{n+1}=\omega v_k^*+(1-\omega)v_k^n.
\]
The first update uses \(v_k^*\) directly. The current outer-controller
relaxation \(\omega\) applies to all layers before computing the complete
Dowson integrals. Pressure mobility, velocity weights and shear therefore
use one consistent profile; the conservative phase inventories are not
blended. This damps the observed alternating fill/viscosity iteration without
changing its physical fixed point. After the pressure and energy solves,
acceptance additionally requires
\[
R_\mu=\max_{i,j,k}\frac{|v_k-v_k^*(\theta_l,c_d,p,T_k)|}{v_k^*(\theta_l,c_d,p,T_k)}
<\texttt{coupling_tolerance}.
\]
This comparison uses the unrelaxed current-state law and is independent of
\(\omega\). The native summary records `viscosity_residual`, and the campaign
checks it independently from the saved layer temperatures and inverse
viscosities. Transient and other cavitation paths retain their prior update.

The measured external load is prescribed, while the shaft position is solved:
\[
R_F=\frac{\lVert\mathbf F_{\rm film}+\mathbf W\rVert}{\lVert\mathbf W\rVert}
<\texttt{steady_force_tolerance}.
\]
For the aligned circular geometry, \(h_{\min}=c-\sqrt{x^2+y^2}\). A viscosity
reduction generally decreases capacity at fixed position; the load-controlled
solve generally responds with greater eccentricity and a smaller minimum
gap. Thus the Couette shear scale \(\mu U/h\) does not fall in direct
proportion to viscosity. The position, pressure, local temperature, composition
and void distribution must be solved together before assigning a torque change.

The position Jacobian and line search evaluate disposable steady pressure
states. Their failed direct solves are rejected so that the load search can
try the opposite finite difference or a smaller displacement; pseudo-time
continuation remains available for the main pressure state. Accepted trials
must satisfy the same equations and admissibility gates, and rejected trials
restore both current fields and old histories. A previous accepted CFD shaft
position is a numerical starting guess only, never an imposed eccentricity.

This load balance contains fluid pressure only. No surface-roughness or
asperity-contact load/friction law is included, so a small predicted gap does
not itself introduce contact friction. Reported journal torque remains the
integrated fluid shear.

The 02f campaign's wall-clock deadline is an operational control, not a
convergence criterion. Its default is disabled after the inherited two-hour
limit stopped progressing full-grid solves. A requested `--timeout SECONDS`
still produces a rejected, preserved result if reached. Iteration caps and
all residual, constitutive, conservation and admissibility gates remain in
force; removing the clock limit does not make an unconverged state acceptable.

### 19 October selected-case inner tolerance

The user reduced the mixture comparison to the exact ten selected windows,
P01--P10, and chose an inner pressure/species tolerance of 1e-8. The new
`coupled_nonlinear_tolerance` sets the default PETSc absolute stopping threshold
and the independent infinity-norm acceptance of the scaled conservative
residual vector. Its native default remains 1e-10 for other configurations.
The 02f driver writes 1e-8 and no longer forces a separate 1e-13 PETSc override.
The raw-field verifier reads the same frozen configuration value.

This is an inner algebraic residual tolerance, not a relative torque error or
the outer thermal criterion. Normal outer coupling remains 1e-6, force balance
remains 0.002, and oil/species conservation remains below
min(coupling_tolerance, 1e-8). Pressure positivity, phase bounds, bubble-point
equilibrium, property-table/EOS domains, positive gaps and energy checks remain
required. No result is promoted solely because its nonlinear norm is small.
The previous 800 rpm gas-on attempt repeatedly reached norms below 1e-10 but
never reached its imposed 1e-13 threshold before the first thermal update.
Local probes and subsequent remote checks must establish the effect of this
change; no full-grid accuracy improvement is asserted before those results.

The free-gas inventory roundoff guard permits a scaled negative inventory no
larger than max(1e-9, min(coupled_nonlinear_tolerance, 1e-8)); the default guard
is unchanged. This avoids rejecting a 1e-8-converged state merely for a few
1e-9 of negative scaled inventory. It does not clip fields or bypass the
independent [-1e-8, 1+1e-8] phase-fraction and conservation checks. The first
selected 1100 rpm coarse probe exposed this mismatch and is retained as a
failed check, not counted as accepted evidence.

The pressure bound p >= p_v(T) follows the same rule (20 October). The coupled
unknown and its residual are scaled by

$$p_\mathrm{scale}=\max\left(10^5,\ p_\mathrm{south},\ p_\mathrm{north},\ \frac{\mu|\omega|R^2}{c^2}\right),$$

so an accepted Newton state may satisfy the cavitation complementarity only to
within the inner tolerance in scaled units. A cell is admissible when

$$p-p_v(T)\ \ge\ -\max\left(10^{-9},\ \min(\epsilon_\mathrm{NL},10^{-8})\right)p_\mathrm{scale},$$

where $\epsilon_\mathrm{NL}$ is `coupled_nonlinear_tolerance`. At 12.7 bar and
1e-8 this is about 1.3e-2 Pa, nine orders below the supply pressure and
physically negligible. The previous fixed 1e-6 Pa bound rejected the gas-off
P01/P07 pilots at -7.0e-3 and -4.0e-3 Pa while they were still converging.
The independent verifier applies the identical bound. Fields are not clipped.

### 20 October film-fraction admissibility

In a full-film cell the volume complementarity residual

$$\mathcal{R}_V=\frac{2(ab-g)}{a+b+\sqrt{a^2+b^2+2g}},\quad a=\frac{\rho_g c}{g_\mathrm{scale}},\ b=1-\theta,\ g=\frac{m_g c}{g_\mathrm{scale}h},$$

reduces for small $b$ and $g$ to $\mathcal{R}_V\approx b-g/a=b-\alpha_g$, where
$\alpha_g=m_g/(\rho_g h)$ is the written gas fraction. A Newton state accepted
with $|\mathcal{R}_V|\le r$ and an admitted negative free-gas round-off
($m_g\ge-r\,g_\mathrm{scale}$) can therefore have $\theta-1\le r-\alpha_g$. The
solver and verifier accept

$$\theta\le1+r+\max(0,-\alpha_g),\qquad r=\max\left(10^{-9},\min(\epsilon_\mathrm{NL},10^{-8})\right),$$

instead of a fixed $1+10^{-8}$. The verifier also requires
$\alpha_g\ge-\min\left(r(c/h)/a_\mathrm{ideal},\sqrt{rc/(2h)}\right)$, the first
term being the free-gas guard with an ideal-gas $a$ (conservative for $Z\le1$)
and the second following from $a^2+b^2+2g\ge0$. Gas-on P01 was rejected at
$\theta-1=3.8\times10^{-8}$ with scaled free gas $-2.4\times10^{-9}$ while its
outer iteration was still converging.

### 20 October relaxed outer coupling tolerance

The remaining selected jobs use an outer coupling tolerance of
$\epsilon_\mathrm{c}=10^{-4}$ instead of $10^{-6}$, by explicit user choice.
The outer loop is an under-relaxed fixed point ($\omega=0.3$) between pressure,
temperature, gas inventory, viscosity and the cross-film profile; each listed
residual is a relative change per outer iteration and must fall below
$\epsilon_\mathrm{c}$. With $\omega=0.3$ the residual contracts by roughly 0.7
per iteration, so stopping at $10^{-4}$ removes about $\ln(100)/\ln(1/0.7)\approx13$
iterations per thermal pass. In the accepted pilots the fluid force at
$R_p<10^{-4}$ differed from its $10^{-6}$ value by at most 0.011%, below the
0.2% force-balance tolerance. The same value sets the energy-balance and
cross-film collocation/viscosity gates; inner Newton tolerance, conservation
(1e-8), force balance, cooling closure and admissibility are unchanged.

The independent bubble-point check uses the same outer tolerance. In free-gas
cells it requires $|p-p_b(c,T)|\le\max(1\ \mathrm{Pa},\ \epsilon_\mathrm{c}\,p_\mathrm{scale})$,
because the written $T$ and $c$ lag the last hydraulic solve by a relative
change of order $\epsilon_\mathrm{c}$. At $10^{-4}$ the measured mismatch was
1.8-4.2 Pa against an allowance of about 127 Pa; at $10^{-6}$ it was below
0.02 Pa.

## Mixed lubrication: asperity contact and average flow factors (23 October 2026)

Smooth rigid-surface Reynolds models cannot carry load or torque through
asperity contact. When the minimum film approaches a few times the composite
roughness, part of the load is supported by asperities and the shaft torque gains
a boundary-friction contribution. The implementation follows the standard
load-sharing approach used for mixed-lubricated journal bearings: a Patir-Cheng
average Reynolds equation for the fluid and a Greenwood-Tripp statistical
contact model for the asperities. Both are independent switches and default off.

### Film ratio

With shaft and bush RMS roughness $R_{q,s}$, $R_{q,b}$,

$$\sigma=\sqrt{R_{q,s}^2+R_{q,b}^2},\qquad H=\frac{h}{\sigma}.$$

Conventional regimes: $H>3$-4 full film, $1<H<3$ mixed, $H<1$ boundary.

### Greenwood-Tripp asperity contact (`asperity_contact_model = GREENWOOD_TRIPP`)

$$p_c=K E' F_{5/2}(H),\qquad K=\frac{16\sqrt2\,\pi}{15}(\eta\beta\sigma)^2\sqrt{\frac{\sigma}{\beta}},\qquad
\frac1{E'}=\frac{1-\nu_1^2}{E_1}+\frac{1-\nu_2^2}{E_2},$$

$$F_{5/2}(H)=4.4086\times10^{-5}(4-H)^{6.804}\ (H<4),\qquad 0\ (H\ge4).$$

$\eta$ is the summit density, $\beta$ the summit radius and $\sigma$ the RMS
roughness. $\eta\beta\sigma$ is nearly constant for engineering surfaces,
0.03-0.05 (Whitehouse and Archard); $\sigma/\beta$ is typically
$10^{-4}$-$10^{-2}$ and should come from a surface measurement. The contact
pressure acts on the bearing-surface normal like the fluid pressure, so the
steady equilibrium balances fluid plus contact reaction against the external
load. The asperity shear is $\tau_c=f_b p_c$ with boundary friction coefficient
$f_b$ (0.1 is the common choice for lubricated steel contacts). It drags the
bearing along the rotation and adds $f_b p_c R$ per unit shaft area to the shaft
torque. The Reynolds solve itself is unchanged by contact.

### Patir-Cheng average flow model (`flow_factor_model = PATIR_CHENG`)

For isotropic Gaussian roughness ($\gamma=1$) the averaged mass flux per unit width is

$$\dot m=-\rho\,\phi_x\frac{h^3}{12\mu}\nabla p+\rho\frac{U}{2}\theta\left(h+\sigma\phi_s\right),$$

$$\phi_x=1-0.90e^{-0.56H}\ (H\ge0.5),\qquad
\Phi_s=\begin{cases}1.899H^{0.98}e^{-0.92H+0.05H^2}&H\le5\\1.126e^{-0.25H}&H>5\end{cases},\qquad
\phi_s=\frac{R_{q,s}^2-R_{q,b}^2}{\sigma^2}\Phi_s,$$

with the shaft as the moving surface 1. The shear flow factor vanishes for
equally rough surfaces. The wall shear on the shaft is

$$\tau=\frac{\mu U}{h}\left(\phi_f+\frac{R_{q,s}^2-R_{q,b}^2}{\sigma^2}\Phi_{fs}\right)+\phi_{fp}\frac{h}{2}\frac{\partial p}{\partial x},$$

$$\phi_f=\tfrac{35}{32}z\left[(1-z^2)^3\ln\frac{z+1}{z-1}+\tfrac{z}{15}\left(66+z^2(30z^2-80)\right)\right]\ (H>3),$$

$$\phi_f=\tfrac{35}{32}z\left[(1-z^2)^3\ln\frac{z+1}{\varepsilon^*}+\tfrac1{60}\left(-55+z(132+z(345+z(-160+z(-405+z(60+147z))))))\right]\ (H\le3),$$

$z=H/3$, $\varepsilon^*=1/300$, $\Phi_{fs}=11.1H^{2.31}e^{-2.38H+0.11H^2}$ for
$H\le7$ (0 above), $\phi_{fp}=1-1.40e^{-0.66H}$. Fits below $H=0.5$ are held at
their $H=0.5$ value. The coefficients were checked by continuity: both $\Phi_s$
branches give 0.3226 at $H=5$, both $\phi_f$ branches give $7/6$ at $H=3$, and
$\phi_f\to1$ as $H\to\infty$ (the $z^5$, $z^3$ and $z$ terms cancel and the
$1/z$ term gives exactly 1).

Assumptions and limits:

- The expected local gap $E[h_T]$ is replaced by $h$ (within 0.5% for $H\ge2$),
  so the contact factor of Wu and Zheng is omitted; results for $H<2$ are
  indicative only.
- Flow factors apply to the JFO and FBNS fluxes (Poiseuille mobility and the
  mass-conserving Couette flux, which carries the shear-flow term with the film
  content). COUPLED and Gumbel do not apply them; validation rejects the
  combination. Asperity contact works with every cavitation model.
- The $\phi_s$ and $\Phi_{fs}$ sign convention (shaft as surface 1) follows
  Patir and Cheng (1979) but has not been verified here against an independent
  unequal-roughness solution; equal roughness makes both terms vanish.
- Rigid, aligned surfaces: no elastic deformation of the film, misalignment,
  edge loading, wear or running-in. These usually lower the real local minimum
  gap, so mixed lubrication in practice starts at larger computed $h_{min}$.

References: Greenwood and Tripp, Proc. IMechE 185 (1970-71) 625-633; Patir and
Cheng, J. Lubr. Technol. 100 (1978) 12-17 and 101 (1979) 220-229; Arcoumanis et
al. (1997) for the $F_{5/2}$ fit; Gelinck and Schipper, Tribol. Int. 33 (2000)
175-181; Lu, Khonsari and Gelinck, J. Tribol. 128 (2006) 789-794; Song, Wu and
Yuan, AIP Advances 12 (2022) 015213 (journal bearing, Patir-Cheng with
Greenwood-Tripp, $f_b=0.1$).
