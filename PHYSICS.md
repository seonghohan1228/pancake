# Physics and Numerical Methods

## 1. Problem Overview

**pancake** simulates the dynamics of a spherical soap bubble film of radius $R$. The film is a thin liquid shell of local thickness $h(\mathbf{x}, t) \ll R$. The governing physics are:

- **Gravitational drainage** — body force drives flow from the top towards the equator.
- **Marangoni convection** — surface-tension gradients $\nabla_s \sigma(T)$ drive tangential flow.
- **Evaporation** — diffusion-limited vapour loss thins the film (Sultan model).
- **Disjoining pressure** — van der Waals repulsion stabilises very thin films and triggers rupture.

The domain is the **upper hemisphere** (above the equatorial rim). The rim is treated as a fixed boundary: no-slip velocity, pinned film thickness $h_{rim}$, pinned temperature $T_{rim}$.

---

## 2. Coordinate System — Inverse Stereographic Projection

The hemisphere is parameterised by the **inverse stereographic projection** from a 2D square domain $[-L_{box}, L_{box}]^2$ (with $L_{box} > R$) onto the hemisphere of radius $R$:

$$x = \frac{2R^2\,u}{r_{2D}^2 + R^2}, \quad y = \frac{2R^2\,v}{r_{2D}^2 + R^2}, \quad z = R\,\frac{R^2 - r_{2D}^2}{r_{2D}^2 + R^2}$$

where $r_{2D}^2 = u^2 + v^2$. The projection is **conformal** with scale factor:

$$\lambda(u,v) = \frac{2R^2}{r_{2D}^2 + R^2}$$

so that all angles are preserved and $|\partial\mathbf{x}/\partial u| = |\partial\mathbf{x}/\partial v| = \lambda$.

Key mapping properties:
- $(u, v) = (0, 0)$ maps to the **north pole** $(0, 0, R)$.
- $r_{2D} = R$ maps to the **equator** (rim).
- $r_{2D} > R$ maps to the southern hemisphere and is **excluded** by the circular mask.

The structured 2D grid has $n_u \times n_v$ uniform cells with spacing $\Delta u = \Delta v = 2L_{box}/n_u$. Each cell is projected to 3D for geometry computation and VTK output.

---

## 3. Geometry Engine

Per-cell quantities are precomputed in `StereoMesh::compute_geometry()`:

**Cell area** — via cross product of 3D cell diagonals $\mathbf{d}_1 = \mathbf{p}_{NE} - \mathbf{p}_{SW}$, $\mathbf{d}_2 = \mathbf{p}_{NW} - \mathbf{p}_{SE}$:
$$A_{i,j} = \tfrac{1}{2}|\mathbf{d}_1 \times \mathbf{d}_2|$$

The sum $\sum_{i,j} A_{i,j}$ over active cells converges to the analytic hemisphere area $2\pi R^2$ as resolution increases (error $< 0.14\%$ at $64 \times 64$).

**Face chord lengths** — 3D distance between the two endpoints of each cell face:
$$L_{e} = |\mathbf{p}(u_{i+1/2}, v_{j+1/2}) - \mathbf{p}(u_{i+1/2}, v_{j-1/2})|, \quad L_{n} = |\mathbf{p}(u_{i+1/2}, v_{j+1/2}) - \mathbf{p}(u_{i-1/2}, v_{j+1/2})|$$

**Centre-to-centre distances** — 3D distance between adjacent cell centres:
$$d_u^{(i,j)} = |\mathbf{c}_{i+1,j} - \mathbf{c}_{i,j}|, \quad d_v^{(i,j)} = |\mathbf{c}_{i,j+1} - \mathbf{c}_{i,j}|$$

**Outward surface normal** — for a sphere, the outward normal equals the normalised position vector:
$$\hat{n}_{i,j} = \mathbf{c}_{i,j} / |\mathbf{c}_{i,j}|$$

**Surface gravity** — the component of gravity tangent to the surface:
$$\mathbf{g}_s = \mathbf{g} - (\mathbf{g} \cdot \hat{n})\hat{n}$$

Decomposed into $(u, v)$ tangent components by projecting onto the normalised partial derivatives of the stereographic map:
$$g_u = \mathbf{g}_s \cdot \hat{e}_u, \quad g_v = \mathbf{g}_s \cdot \hat{e}_v, \quad \hat{e}_u = \frac{\partial\mathbf{x}/\partial u}{|\partial\mathbf{x}/\partial u|}$$

Expected behaviour: $|\mathbf{g}_s| \to 0$ at the north pole; $|\mathbf{g}_s| \to |\mathbf{g}|$ at the equator.

---

## 4. Circular Domain Mask

Cells with $r_{2D} > R$ are **blanked** (inactive). Active cells adjacent to at least one blanked or out-of-domain cell are **rim cells** carrying Dirichlet boundary conditions.

Ghost cells outside the domain are filled with **zero-gradient** (Neumann) values: the boundary physical cell value is copied into all ghost layers. This is handled by `Communicator::fill_v_ghosts` (v-direction) and the open-boundary logic in `Communicator::update_ghosts` (u-direction).

---

## 5. Surface FVM Operators (`SurfaceFVM` namespace)

All operators act on the 2D stereographic grid, using the precomputed 3D geometry from `StereoMesh`. They accumulate contributions into a `BubbleLinearSystem` (5-point stencil: $a_P, a_E, a_W, a_N, a_S$, source).

Masked (inactive) cells are skipped entirely. Rim cells are populated normally and then overwritten by `BubbleLinearSystem::apply_mask()`.

### ddt — transient term

$$\frac{\partial\phi}{\partial t} \approx \frac{\phi - \phi^n}{\Delta t}$$

Integrated over cell area $A_{i,j}$:

$$a_P \mathrel{+}= \frac{A_{i,j}}{\Delta t}, \qquad b \mathrel{+}= \frac{A_{i,j}}{\Delta t}\,\phi^n_{i,j}$$

### ddt_weighted — weighted transient

Same as `ddt` with an additional cell weight $w_{i,j}$ (used for the energy equation where $w = \rho c_p h$):

$$a_P \mathrel{+}= \frac{w_{i,j}\,A_{i,j}}{\Delta t}, \qquad b \mathrel{+}= \frac{w_{i,j}\,A_{i,j}}{\Delta t}\,\phi^n_{i,j}$$

### laplacian — diffusion

For each face between cells $P$ and $N$ (neighbour):

$$c_f = \gamma_f \cdot \frac{L_{face}}{d_{PN}}$$

where $\gamma_f$ is the harmonic mean of $\gamma_P$ and $\gamma_N$, $L_{face}$ is the 3D face chord length, and $d_{PN}$ is the 3D centre-to-centre distance. Accumulated:

$$a_P \mathrel{+}= c_f, \qquad a_{dir} \mathrel{+}= c_f \quad (dir \in \{E, W, N, S\})$$

West faces at MPI rank boundaries use `StereoMesh::east_face_len(gu-1, j)` and `d_u_east(gu-1, j)`, which compute geometry analytically for any global cell index without MPI communication.

### divergence — convection

Pre-computed face fluxes $F_e, F_w, F_n, F_s$ (outward positive). Upwind scheme (implicit):

$$a_P \mathrel{+}= \max(F_e, 0) + \max(-F_w, 0) + \max(F_n, 0) + \max(-F_s, 0)$$
$$a_E \mathrel{+}= \max(-F_e, 0), \quad a_W \mathrel{+}= \max(F_w, 0), \quad a_N \mathrel{+}= \max(-F_n, 0), \quad a_S \mathrel{+}= \max(F_s, 0)$$

TVD deferred correction (Van Leer / Minmod) adds an explicit high-order correction to the source using $\phi^n$ from the previous timestep.

### add_source — volumetric source

$$b \mathrel{+}= s_{i,j} \cdot A_{i,j}$$

### gradient — cell-centred gradient

Central difference with 3D chord lengths as denominators:

$$(\nabla\phi)_u = \frac{\phi_{i+1,j} - \phi_{i-1,j}}{d_e + d_w}, \qquad (\nabla\phi)_v = \frac{\phi_{i,j+1} - \phi_{i,j-1}}{d_n + d_s}$$

---

## 6. Thickness Transport

The film thickness equation (Stokes thin-film limit):

$$\frac{\partial h}{\partial t} + \nabla_s \cdot (h\,\mathbf{u}_s) = \dot{J}_{evap}$$

### Gravity-driven Stokes velocity

In the creeping-flow (low-Re) limit, the depth-averaged surface velocity due to gravity is:

$$\mathbf{u}_s = \frac{\rho_l\,h^2}{3\mu}\,\mathbf{g}_s$$

where $\mathbf{g}_s$ is the surface-tangent component of gravity (precomputed in `StereoMesh`).

### Face flux interpolation

Cell-centred velocities are interpolated to face centres by simple averaging:

$$F_e = \frac{u_u(i,j) + u_u(i+1,j)}{2} \cdot L_{e}^{3D}$$

Ghost layers of $u_u, u_v$ must be updated via `Communicator::update_ghosts` before calling `compute_face_fluxes`.

### FVM assembly

`ThicknessTransport::solve` assembles:

$$\text{ddt}(h) + \text{divergence}(F_u, F_v, h) + \text{add\_source}(J_{evap}) + \text{apply\_mask}(h_{rim}) \to \text{solve}$$

Post-solve, $h$ is clamped to $h \ge h_{min}$ and blanked cells are zeroed.
