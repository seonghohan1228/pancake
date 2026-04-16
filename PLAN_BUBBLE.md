# Soap Bubble Thin-Film Multiphysics Solver — Master Implementation Plan

This document is the complete, staged implementation plan for adapting the `pancake` journal bearing FVM solver into a multiphysics solver for spherical soap bubble dynamics. It is written to be consumed directly by an implementing agent (Sonnet) and must be followed phase by phase, with user review between each phase.

**Branch**: `feature/bubble`

---

## Execution Protocol

> **CRITICAL**: Execute this plan in strict sequential phases. Complete **one phase at a time**, deliver the code, tests, and compilation proof, then **stop and wait for user review** before proceeding. Do not write the entire codebase at once.

> **REUSE RULE**: Before writing any new module, check `AUDIT.md` (generated in Phase 1). If the existing codebase already implements a needed feature (linear solvers, TVD schemes, PETSc wrappers, config parsing), use it directly. Do not rewrite.

> **CONVENTION RULE**: All new code must follow the conventions in `CLAUDE.md` exactly: snake_case, minimal comments (but complete header documentation), scientific terminology, modular swappable design.

---

# Part I — Architecture Overview

## 1.1 Physics Summary

We simulate a hemispherical soap bubble film clamped at its equatorial rim. The film is a thin viscous sheet whose thickness $h(\mathbf{x}, t)$ evolves under:

- **Gravitational drainage**: gravity projects onto the curved surface, driving downward flow that thins the top.
- **Thermal Marangoni convection**: non-uniform temperature creates surface tension gradients $\nabla_s \sigma(T)$ that drive flow from hot (low $\sigma$) to cold (high $\sigma$) regions.
- **Evaporation**: the Sultan diffusion-limited model removes mass from the film, coupling temperature to thickness loss.
- **Disjoining pressure**: van der Waals attraction across the thin film ($\Pi \propto h^{-3}$) accelerates thinning and triggers rupture.
- **Capillary pressure**: Laplace pressure $\Delta P = 2\sigma/R_{bubble}$ maintains the spherical shape (quasi-static assumption for the bubble geometry; the film surface itself is fixed as a hemisphere).

The computational domain is the **hemispherical surface** of the bubble, discretised via inverse stereographic projection of a uniform 2D Cartesian grid.

## 1.2 Key Architectural Differences from Journal Bearing Solver

| Aspect | Journal Bearing (`pancake`) | Soap Bubble (`pancake-bubble`) |
|--------|----------------------------|-------------------------------|
| Domain | Cylinder $(θ, z)$ | Hemisphere via inverse stereographic $(u, v)$ |
| Metric | Uniform: $R\,\Delta\theta\,\Delta z$ | Non-uniform: $A_{cell}(u,v)$ from projection |
| Primary PDEs | Reynolds eqn (single elliptic) | 2D Navier-Stokes + continuity + energy (coupled) |
| Film thickness | Static geometry input | Dynamic transported scalar |
| Velocity | Derived post-solve | Primary unknown (momentum eqn) |
| Pressure | Primary unknown | Correction variable (SIMPLE/PISO) |
| Boundary | Periodic $θ$, Dirichlet $z$ | Circular rim (no-slip, pinned $h$, $T$) |
| Invalid cells | None | Corners outside circular mask |

## 1.3 Module Map (New and Reused)

```
pancake/src/
├── bubble_config.hpp/cpp       [NEW]  Bubble-specific SimulationConfig
├── stereo_mesh.hpp/cpp         [NEW]  Inverse stereographic grid + geometry engine
├── mask.hpp/cpp                [NEW]  Circular domain mask + ghost blanking
├── surface_operators.hpp/cpp   [NEW]  Surface gradient, divergence, laplacian (metric-aware)
├── thickness_transport.hpp/cpp [NEW]  Continuity equation for h
├── momentum.hpp/cpp            [NEW]  2D N-S with SIMPLE/PISO
├── energy.hpp/cpp              [NEW]  Energy equation + Sultan evaporation
├── marangoni.hpp/cpp           [NEW]  σ(T) model + Marangoni body force
├── disjoining.hpp/cpp          [NEW]  van der Waals Π(h) + rupture detection
├── bubble_io.hpp/cpp           [NEW]  VTK output for hemispherical surface
├── bubble_main.cpp             [NEW]  Top-level time loop
│
├── field.hpp/cpp               [REUSE] Unchanged — ghost cells, old-time, staggered
├── linear_system.hpp/cpp       [REUSE] Adapted — per-cell metric replaces uniform volume
├── communicator.hpp/cpp        [REUSE] Adapted — open (non-periodic) boundaries
├── fvm.hpp/cpp                 [ADAPT] Generalize metric terms OR write surface_operators
├── utils.hpp                   [REUSE] Logger, BC helpers
├── config.hpp/cpp              [REUSE] Parser infrastructure (extend for bubble params)
```

## 1.4 Overall FVM Solver Loop (Per Timestep)

This is the target algorithm. Each numbered step corresponds to a function call in `bubble_main.cpp`:

```
for each timestep:
    1. Update geometry (if dynamic) & recompute 3D cell areas A(i,j)
       Apply circular mask; enforce rim BCs on all fields

    2. Compute evaporation flux J(i,j) from current T field (Sultan model)

    3. Solve Continuity FVM for new h:
       ∂h/∂t + ∇_s·(h u_s) = -J/ρ_l

    4. Solve Energy FVM for new T:
       ρ c_p h (∂T/∂t + u_s·∇_s T) = ∇_s·(k h ∇_s T) - J ΔH_v + h_conv(T_amb - T)

    5. Update σ(T), compute Marangoni body force F_M = -γ_T ∇_s T

    6. SIMPLE/PISO pressure-velocity coupling loop:
       a. Momentum predictor → u*, v* (face velocities via Rhie-Chow)
       b. Pressure correction equation → p'
       c. Correct velocities and pressure
       d. Check SIMPLE convergence, repeat if needed

    7. Compute disjoining pressure Π(h), add to pressure field

    8. Check residuals, adaptive Δt (CFL + Marangoni + disjoining stability)

    9. Check rupture: if min(h) < h_rupture → halt

   10. Output VTK if write interval reached

   11. Store old-time fields, advance t
```

---

# Part II — Mathematical Formulation

All equations below use the **surface operators** $\nabla_s$, $\nabla_s \cdot$, $\nabla_s^2$ which act on the 2D hemispherical manifold. In the FVM discretisation, these reduce to face fluxes weighted by 3D face areas and distances derived from the inverse stereographic mapping.

## 2.1 Inverse Stereographic Projection

Map a uniform Cartesian grid $(u, v)$ on $[-L_{box}, L_{box}]^2$ to the hemisphere of radius $R$:

$$r_{2D}^2 = u^2 + v^2$$

$$x = \frac{2R^2 u}{r_{2D}^2 + R^2}, \qquad y = \frac{2R^2 v}{r_{2D}^2 + R^2}, \qquad z = R\,\frac{R^2 - r_{2D}^2}{r_{2D}^2 + R^2}$$

The **south pole** $(0, 0, -R)$ maps to $r_{2D} \to \infty$; the **north pole** $(0, 0, R)$ maps to the origin $(u, v) = (0, 0)$.

**Domain validity**: Only cells with $r_{2D}^2 \le R_{mask}^2$ are physical. Cells in the square corners where $r_{2D}^2 > R_{mask}^2$ are masked (ghost-blanked). We choose $R_{mask}$ slightly larger than $R$ to capture the full hemisphere with some buffer, or exactly $R$ for a hemispherical cap pinned at the equator ($z = 0$).

**Metric tensor** on the surface (for the stereographic parameterisation):

$$g_{uu} = g_{vv} = \left(\frac{2R^2}{r_{2D}^2 + R^2}\right)^2, \qquad g_{uv} = 0$$

The conformal factor:

$$\lambda(u,v) = \frac{2R^2}{r_{2D}^2 + R^2}$$

so the metric is $ds^2 = \lambda^2 (du^2 + dv^2)$. This is a **conformal** mapping: angles are preserved, and the area element is $dA = \lambda^2\,du\,dv$.

## 2.2 3D Cell Area Computation

For a quadrilateral cell with 3D vertices $\mathbf{p}_1, \mathbf{p}_2, \mathbf{p}_3, \mathbf{p}_4$ (ordered), the area is:

$$A_{cell} = \frac{1}{2}\left|(\mathbf{p}_3 - \mathbf{p}_1) \times (\mathbf{p}_4 - \mathbf{p}_2)\right|$$

Similarly, each face has a 3D length:

$$L_{face} = |\mathbf{p}_{end} - \mathbf{p}_{start}|$$

and the distance between adjacent cell centres projected onto 3D:

$$d_{face} = |\mathbf{c}_{neighbour} - \mathbf{c}_{owner}|_{3D}$$

These replace the uniform $R\,\Delta\theta$, $\Delta z$ in the existing FVM operators.

## 2.3 Continuity Equation (Thickness Transport)

$$\frac{\partial h}{\partial t} + \nabla_s \cdot (h\,\mathbf{u}_s) = -\frac{J}{\rho_l}$$

FVM integral form over cell $\Omega$ with area $A$:

$$A\,\frac{h^{n+1} - h^n}{\Delta t} + \sum_{faces} (h\,\mathbf{u}_s \cdot \hat{n})\,L_{face} = -\frac{J}{\rho_l}\,A$$

where $h$ at faces uses TVD interpolation (van Leer or Superbee) for boundedness.

**Non-negativity**: After solve, clamp $h \ge h_{min}$ (e.g., $h_{min} = 10\text{ nm}$).

## 2.4 2D Momentum Equation (Surface Navier-Stokes)

$$\rho_l\left(\frac{\partial \mathbf{u}_s}{\partial t} + \mathbf{u}_s \cdot \nabla_s \mathbf{u}_s\right) = -\nabla_s p + \mu \nabla_s^2 \mathbf{u}_s + \mathbf{F}_{Marangoni} + \mathbf{F}_{drag} + \rho_l \mathbf{g}_s$$

where:
- $\mathbf{g}_s$ is gravity projected onto the local tangent plane of the hemisphere
- $\mathbf{F}_{Marangoni} = -\gamma_T \nabla_s T$ with $\gamma_T = d\sigma/dT$
- $\mathbf{F}_{drag}$ is air drag on the film (optional, can be disabled)

**SIMPLE algorithm** adapted for the surface:

1. **Momentum predictor**: Solve for $\mathbf{u}^*$ using lagged pressure $p^n$.
2. **Rhie-Chow face velocities**: Interpolate $\mathbf{u}^*$ to faces with pressure gradient correction to suppress checkerboard modes.
3. **Pressure correction**: Solve $\nabla_s \cdot \left(\frac{A_{face}}{\hat{a}_P} \nabla_s p'\right) = \nabla_s \cdot (h\,\mathbf{u}^*)$ where $\hat{a}_P$ is the momentum diagonal.
4. **Velocity correction**: $\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{1}{\hat{a}_P}\nabla_s p'$.
5. **Pressure update**: $p^{n+1} = p^n + \alpha_p\,p'$ (under-relaxation $\alpha_p \sim 0.3$).

**Velocity components**: Use a **covariant** representation in $(u, v)$ coordinates on the surface. The momentum equations are solved component-wise for $u_u$ and $u_v$ (the contravariant velocity components in the $(u, v)$ parameterisation), with cross-metric corrections handled explicitly.

## 2.5 Energy Equation

$$\rho_l c_p h\left(\frac{\partial T}{\partial t} + \mathbf{u}_s \cdot \nabla_s T\right) = \nabla_s \cdot (k_l h\,\nabla_s T) - J\,\Delta H_v + h_{conv}(T_{ambient} - T)$$

FVM form: this is a standard convection-diffusion equation for $T$ with:
- Diffusion coefficient: $\Gamma_T = k_l h$ (conductivity $\times$ thickness)
- Convective flux: $h\,\mathbf{u}_s$ (same as continuity fluxes)
- Source: $-J\,\Delta H_v + h_{conv}(T_{ambient} - T)$
- "Volume" weighting: $\rho_l c_p h$ multiplies the transient term (use `ddt_weighted` pattern)

## 2.6 Evaporation — Sultan Diffusion-Limited Model

$$J(T) = \frac{\rho_{gas}\,D_{AB}\,(\omega_{sat}(T) - \omega_\infty)}{R_c}$$

where:
- $\rho_{gas}$ = air density
- $D_{AB}$ = binary diffusion coefficient (water vapour in air)
- $\omega_{sat}(T)$ = saturation mass fraction at local temperature (Clausius-Clapeyron or Antoine)
- $\omega_\infty$ = ambient humidity mass fraction
- $R_c$ = characteristic diffusion length (typically bubble radius $R$)

This is evaluated cell-by-cell as an explicit source — no PDE solve needed.

## 2.7 Marangoni Stress

$$\mathbf{F}_{Marangoni} = -\gamma_T\,\nabla_s T, \qquad \gamma_T = \frac{d\sigma}{dT} < 0 \text{ (water)}$$

Surface tension model (linear):

$$\sigma(T) = \sigma_0 + \gamma_T (T - T_{ref})$$

The Marangoni force is an explicit body force in the momentum equation. It drives flow from hot (thin, low $\sigma$) regions toward cold (thick, high $\sigma$) regions.

## 2.8 Disjoining Pressure

$$\Pi(h) = \frac{A_{Hamaker}}{6\pi\,\max(h, h_{min})^3}$$

where $A_{Hamaker} \sim 10^{-20}\text{ J}$ for soap film (water-air-water).

This adds a pressure term: the effective pressure becomes $p_{eff} = p - \Pi(h)$, which modifies the pressure gradient in momentum:

$$-\nabla_s p_{eff} = -\nabla_s p + \nabla_s \Pi(h) = -\nabla_s p + \Pi'(h)\,\nabla_s h$$

## 2.9 Rupture Criterion

Halt the simulation when:

$$\min_{cells} h(i,j) < h_{rupture} = 10\text{ nm}$$

Report the location $(u, v)$ and time of rupture. This is the primary output of interest.

## 2.10 Boundary Conditions (Rim / Equator)

At the circular mask boundary ($r_{2D} = R_{mask}$, corresponding to the equator $z = 0$):

| Field | BC Type | Value |
|-------|---------|-------|
| $\mathbf{u}_s$ | No-slip Dirichlet | $\mathbf{u}_s = 0$ |
| $h$ | Dirichlet (pinned) | $h = h_{rim}$ |
| $T$ | Dirichlet | $T = T_{rim}$ |
| $p$ | Zero-gradient (Neumann) | $\partial p / \partial n = 0$ |

Cells outside the circular mask are blanked: they do not participate in the FVM solve. The mask is enforced via:
1. A boolean `mask(i,j)` field (1 = active, 0 = blanked).
2. In LinearSystem assembly: masked cells get $a_P = 1$, $\text{source} = \phi_{BC}$, all off-diagonals = 0.
3. In VTK output: a `vtkGhostType` array marks blanked cells for ParaView.

## 2.11 Time Integration & Stability

Adaptive timestep based on:

$$\Delta t = \min\left(\Delta t_{CFL},\;\Delta t_{Marangoni},\;\Delta t_{disjoining},\;\Delta t_{max}\right)$$

where:

$$\Delta t_{CFL} = C_{CFL} \cdot \min_{cells} \frac{\sqrt{A_{cell}}}{|\mathbf{u}_s|}, \qquad C_{CFL} \sim 0.5$$

$$\Delta t_{Marangoni} = C_M \cdot \min_{cells} \frac{\rho_l\,(\sqrt{A_{cell}})^2}{|\gamma_T|\,|\nabla_s T|\,\sqrt{A_{cell}}/h}$$

$$\Delta t_{disjoining} = C_\Pi \cdot \min_{cells} \frac{\rho_l\,(\sqrt{A_{cell}})^2}{|\Pi'(h)|\,h}$$

Use implicit Euler for diffusion (unconditionally stable) and explicit treatment of Marangoni and disjoining forces with adaptive $\Delta t$.

---

# Part III — Phase-by-Phase Implementation Checklist

## Phase 1: Codebase Audit & Project Setup

**Goal**: Create `feature/bubble` branch, audit existing code, set up configuration and build system for the bubble solver. Produce `AUDIT.md`.

### Tasks

- [ ] **1.1** Create and switch to branch `feature/bubble` from `main`.

- [ ] **1.2** Generate `AUDIT.md` at repo root. For each existing source file, document:
  - What it does.
  - Whether it is **reusable as-is**, **needs adaptation**, or **not applicable** for the bubble solver.
  - Specific functions/classes to reuse and what modifications (if any) are needed.

  Expected reuse assessment:
  | Component | Verdict |
  |-----------|---------|
  | `field.hpp/cpp` | Reuse as-is |
  | `linear_system.hpp/cpp` | Adapt (per-cell metric instead of uniform volume) |
  | `fvm.hpp/cpp` | Reference only — write new `surface_operators` with per-cell metrics |
  | `communicator.hpp/cpp` | Adapt (non-periodic boundaries) |
  | `config.hpp/cpp` | Reuse parser, extend with bubble parameters |
  | `utils.hpp` | Reuse as-is |
  | `io.hpp/cpp` | Replace (need unstructured VTK for hemisphere surface) |
  | `reynolds.hpp/cpp` | Not applicable — skip |
  | `film_thickness.hpp/cpp` | Not applicable — skip |
  | `equation_of_state.hpp` | Not applicable — skip |

- [ ] **1.3** Create `src/bubble_config.hpp/cpp` extending the configuration system.

  Required parameters (with sensible defaults):
  ```cpp
  struct BubbleConfig {
      // Geometry
      double R_bubble = 0.01;        // Bubble radius [m]
      double L_box = 0.012;          // Half-width of 2D square domain [m]
      int n_u = 64, n_v = 64;        // Grid resolution in (u,v)

      // Film properties
      double h_initial = 1.0e-6;     // Initial uniform thickness [m]
      double h_rim = 1.0e-6;         // Rim (equator) thickness BC [m]
      double h_min = 1.0e-8;         // Minimum allowed thickness [m]
      double h_rupture = 1.0e-8;     // Rupture threshold [m]

      // Fluid properties
      double rho_l = 1000.0;         // Liquid density [kg/m^3]
      double mu = 1.0e-3;            // Dynamic viscosity [Pa.s]
      double sigma_0 = 0.025;        // Surface tension at T_ref [N/m]
      double gamma_T = -1.7e-4;      // dσ/dT [N/(m·K)]
      double k_l = 0.6;              // Thermal conductivity [W/(m·K)]
      double c_p = 4186.0;           // Specific heat [J/(kg·K)]

      // Evaporation (Sultan model)
      double rho_gas = 1.2;          // Air density [kg/m^3]
      double D_AB = 2.5e-5;          // Diffusion coeff H2O in air [m^2/s]
      double omega_inf = 0.01;       // Ambient humidity mass fraction [-]
      double T_ref = 293.15;         // Reference temperature [K]
      double T_ambient = 293.15;     // Ambient temperature [K]
      double T_rim = 293.15;         // Rim temperature BC [K]
      double delta_H_v = 2.26e6;     // Latent heat of vaporisation [J/kg]
      double h_conv = 10.0;          // Convective heat transfer coeff [W/(m^2·K)]

      // Disjoining pressure
      double A_hamaker = 1.0e-20;    // Hamaker constant [J]

      // Gravity
      double g_x = 0.0, g_y = 0.0, g_z = -9.81;  // Gravity vector [m/s^2]

      // Numerics
      double dt = 1.0e-6;            // Initial timestep [s]
      double dt_max = 1.0e-3;        // Maximum timestep [s]
      double end_t = 10.0;           // End time [s]
      double cfl = 0.5;              // CFL number
      double alpha_p = 0.3;          // Pressure under-relaxation (SIMPLE)
      double alpha_u = 0.7;          // Velocity under-relaxation (SIMPLE)
      int max_simple_iters = 50;     // Max SIMPLE iterations per timestep
      double simple_tol = 1.0e-4;    // SIMPLE convergence tolerance
      ConvectionScheme advection_scheme = ConvectionScheme::TVD_VANLEER;

      // Toggles (modularity)
      bool enable_evaporation = true;
      bool enable_marangoni = true;
      bool enable_disjoining = true;
      bool enable_gravity = true;

      // Output
      std::string output_dir = "results_bubble";
      double write_interval = 0.01;

      void load_from_file(const std::string& path);
  };
  ```

- [ ] **1.4** Update `CMakeLists.txt`:
  - Add new source files (initially empty stubs).
  - Add a separate target `pancake_bubble` alongside the existing `pancake`.
  - Ensure both targets link MPI and PETSc.

- [ ] **1.5** Create `bubble_main.cpp` with a skeleton:
  - PETSc init, config load, mesh creation, field registration, empty time loop, PETSc finalize.
  - Must compile and run (producing no output beyond a log message).

### Verification (Phase 1)
- [ ] `AUDIT.md` exists and is comprehensive.
- [ ] `cmake .. && make pancake_bubble` compiles without errors.
- [ ] `./pancake_bubble` runs and exits cleanly (prints "Bubble solver initialized" or similar).
- [ ] Existing `pancake` target still compiles and all existing tests pass.

---

## Phase 2: Decoupled Dynamic Mesh & Topology

**Goal**: Implement the inverse stereographic mesh, circular mask, geometry engine, and VTK output for the hemispherical surface. No physics yet — just geometry and I/O.

### Tasks

- [ ] **2.1** Create `src/stereo_mesh.hpp/cpp` — the **Logical Grid Manager + Geometry Engine**.

  **Logical grid** (static, 2D Cartesian):
  ```cpp
  class StereoMesh {
  public:
      StereoMesh(const BubbleConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);

      int n_u_global, n_v_global;        // Global grid size
      int n_u_local, n_v_local;          // Local (this rank)
      int offset_u, offset_v;            // Local starting index
      double du, dv;                     // Uniform spacing in (u,v)
      double R;                          // Bubble radius

      // 2D logical coordinates of cell centres
      double u_center(int i) const;      // u = -L_box + (offset_u + i + 0.5)*du
      double v_center(int j) const;      // v = -L_box + (offset_v + j + 0.5)*dv

      // 2D logical coordinates of cell vertices (corners)
      double u_vertex(int i) const;      // u = -L_box + (offset_u + i)*du
      double v_vertex(int j) const;      // v = -L_box + (offset_v + j)*dv

      // 3D mapping: (u,v) → (x,y,z) on hemisphere
      void map_to_3d(double u, double v, double& x, double& y, double& z) const;

      // Precomputed per-cell geometric quantities
      std::vector<double> cell_area;     // A(i,j) — 3D area of projected quad
      std::vector<double> face_u_length; // |face| for u-direction faces (east/west)
      std::vector<double> face_v_length; // |face| for v-direction faces (north/south)
      std::vector<double> dist_u;        // 3D centre-to-centre distance in u-direction
      std::vector<double> dist_v;        // 3D centre-to-centre distance in v-direction

      // Surface normal and tangent vectors at cell centres
      std::vector<std::array<double,3>> cell_normal;   // Unit outward normal
      std::vector<std::array<double,3>> tangent_u;     // ∂r/∂u (unnormalised)
      std::vector<std::array<double,3>> tangent_v;     // ∂r/∂v (unnormalised)

      // Gravity projection onto surface tangent plane at each cell
      std::vector<std::array<double,2>> g_surface;     // (g_u, g_v) components

      // Domain decomposition: 1D in u-direction (matching existing pattern)
      int rank, size;

      // Index helper
      int idx(int i, int j) const;       // Flat index including ghosts
  };
  ```

  **Geometry engine** (compute all metrics):
  ```cpp
  void StereoMesh::compute_geometry() {
      // For each cell (i,j):
      //   1. Compute 4 corner vertices in 3D via map_to_3d
      //   2. Cell area via cross-product formula:
      //      A = 0.5 * |diag1 × diag2|
      //      where diag1 = p3 - p1, diag2 = p4 - p2
      //   3. Face lengths: |p2 - p1| for u-faces, |p4 - p1| for v-faces
      //   4. Centre-to-centre distances: |c_neighbour - c_owner| in 3D
      //   5. Surface tangent vectors: ∂r/∂u, ∂r/∂v (analytic from projection)
      //   6. Surface normal: (∂r/∂u × ∂r/∂v) / |∂r/∂u × ∂r/∂v|
      //   7. Gravity projection: g_s = g - (g·n)n, decompose into (u,v) components
  }
  ```

  **Inverse stereographic map** (inline):
  ```cpp
  void StereoMesh::map_to_3d(double u, double v, double& x, double& y, double& z) const {
      double r2 = u*u + v*v;
      double denom = r2 + R*R;
      x = 2.0 * R*R * u / denom;
      y = 2.0 * R*R * v / denom;
      z = R * (R*R - r2) / denom;
  }
  ```

- [ ] **2.2** Create `src/mask.hpp/cpp` — circular domain mask.

  ```cpp
  class Mask {
  public:
      Mask(const StereoMesh& mesh, double R_mask);

      bool is_active(int i, int j) const;     // true if r2D < R_mask
      bool is_rim(int i, int j) const;         // true if active but adjacent to masked cell
      int n_active_local() const;
      int n_active_global() const;             // MPI-reduced

      // Apply Dirichlet BC at rim for a scalar field
      void enforce_rim_bc(Field& field, double value) const;

      // Apply no-slip at rim for velocity components
      void enforce_rim_noslip(Field& u_field, Field& v_field) const;

      // Zero out masked cells in a field (for safety)
      void blank_masked(Field& field, double fill = 0.0) const;

  private:
      const StereoMesh& mesh;
      std::vector<bool> active;   // Flat boolean array
      std::vector<bool> rim;      // Flat boolean array
  };
  ```

  **Rim detection**: A cell is a rim cell if it is active AND at least one of its 4 neighbours is inactive (or is at the domain boundary).

- [ ] **2.3** Create `src/bubble_io.hpp/cpp` — VTK output for the hemispherical surface.

  Write **VTK Unstructured Grid** (`.vtu`) files since the masked domain is not rectangular:
  - Only output active cells.
  - Vertex coordinates are the 3D-projected positions.
  - Cell type: `VTK_QUAD` (type 9).
  - Cell data: all registered scalar fields.
  - Include a `vtkGhostType` array that blanks inactive cells (alternative: simply omit them).
  - Rank 0 writes a `.pvd` time-series collection file.
  - Parallel: each rank writes its own `.vtu`; rank 0 writes a `.pvtu` parallel descriptor.

- [ ] **2.4** Adapt `src/communicator.hpp/cpp` for the bubble mesh:
  - Support **non-periodic** boundaries in both $u$ and $v$ (no wrapping).
  - Ghost exchange still via `MPI_Sendrecv` for 1D decomposition in $u$.
  - External boundaries (domain edges) get zero-gradient ghost fill by default.

- [ ] **2.5** Register geometry fields and write initial output in `bubble_main.cpp`:
  ```cpp
  fields.add("cell_area", mesh);    // Filled from StereoMesh::cell_area
  fields.add("mask", mesh);         // 1.0 for active, 0.0 for blanked
  fields.add("g_u", mesh);          // Surface gravity u-component
  fields.add("g_v", mesh);          // Surface gravity v-component
  ```
  Write a single VTK timestep at $t=0$ showing the hemispherical grid with cell areas and mask.

### Verification (Phase 2)
- [ ] **V2.1**: Open the VTK output in ParaView. Visually confirm:
  - The grid forms a hemisphere (or hemispherical cap).
  - Corner cells outside the circle are blanked/absent.
  - Cell area colouring shows smooth variation (largest at north pole where the conformal factor peaks, smallest near equator).
- [ ] **V2.2**: Sum all `cell_area` values across ranks; compare to analytic hemisphere area $2\pi R^2$. Relative error should be $< 1\%$ for $64 \times 64$ grid.
- [ ] **V2.3**: Gravity projection test: at the north pole $(0,0,R)$, $g_s$ should be zero (gravity is purely normal). At the equator, $|g_s|$ should approach $|g|$. Print these values.
- [ ] **V2.4**: Rim detection: count rim cells, verify they form an annular ring. Visualise in ParaView.
- [ ] **V2.5**: Parallel: run with 2 and 4 ranks; verify identical VTK output (up to partitioning).

---

## Phase 3: Thickness Transport & Gravitational Drainage

**Goal**: Solve the continuity equation for $h$ on the hemispherical surface with gravitational drainage only (no evaporation, no Marangoni, no pressure coupling). This validates the surface FVM operators and TVD transport.

### Tasks

- [ ] **3.1** Create `src/surface_operators.hpp/cpp` — metric-aware FVM operators.

  These are the curved-surface analogues of the existing `fvm.hpp` operators. They replace the uniform cylindrical metrics ($R\,\Delta\theta$, $\Delta z$) with per-cell quantities from `StereoMesh`.

  ```cpp
  namespace SurfaceFVM {

      /// Time derivative: ∂φ/∂t over cell area A(i,j).
      /// Adds A(i,j)/dt to aP, A(i,j)/dt * φ^n to source.
      void ddt(LinearSystem& sys, const Field& phi, double dt,
               const StereoMesh& mesh, const Mask& mask);

      /// Weighted time derivative: ∂(w·φ)/∂t.
      /// Adds w(i,j)*A(i,j)/dt to aP and similar to source.
      void ddt_weighted(LinearSystem& sys, const Field& phi, const Field& weight,
                        double dt, const StereoMesh& mesh, const Mask& mask);

      /// Surface Laplacian: ∇_s·(γ ∇_s φ).
      /// Face coefficient = γ_face * face_length / centre_distance.
      /// γ_face via harmonic averaging.
      void laplacian(LinearSystem& sys, const Field& gamma,
                     const StereoMesh& mesh, const Mask& mask);

      /// Surface divergence of convective flux: ∇_s·(F φ).
      /// flux_u(i,j) = volumetric flux through east face = (u_s · n_face) * face_length * h_face
      /// TVD with deferred correction using phi.old.
      void divergence(LinearSystem& sys,
                      const Field& flux_u, const Field& flux_v, const Field& phi,
                      ConvectionScheme scheme,
                      const StereoMesh& mesh, const Mask& mask);

      /// Explicit source: adds source(i,j) * A(i,j) to RHS.
      void add_source(LinearSystem& sys, const Field& source_field,
                      const StereoMesh& mesh, const Mask& mask);
  }
  ```

  **Key implementation detail for `laplacian`**:
  ```cpp
  // East face coefficient (replaces d_z/(R*d_theta) in cylindrical)
  double gamma_e = harmonic_avg(gamma(i,j), gamma(i+1,j));
  double coeff_e = gamma_e * mesh.face_u_length[face_idx] / mesh.dist_u[face_idx];
  ```

  **Masked cells**: All operators must skip masked cells (or the mask is enforced by setting trivial equations for masked cells in LinearSystem before solve).

- [ ] **3.2** Adapt `src/linear_system.hpp/cpp` for the bubble mesh.

  Changes needed:
  - Constructor takes `StereoMesh` instead of `Mesh`.
  - Global row mapping: $(offset_u + i) \times n_v_{global} + j$ (replacing theta-based).
  - **Remove periodic wrapping** in $u$-direction (no periodicity).
  - Add masking support: before `solve()`, set masked cells to identity equation ($a_P = 1$, source = prescribed value, off-diagonals = 0).

  Approach: create a new class `BubbleLinearSystem` or parameterise the existing `LinearSystem`. Prefer the former to avoid breaking the journal bearing solver.

- [ ] **3.3** Create `src/thickness_transport.hpp/cpp`:

  ```cpp
  namespace ThicknessTransport {
      /// Solve ∂h/∂t + ∇_s·(h u_s) = -J/ρ_l for new h.
      /// For Phase 3, u_s is a prescribed gravity-driven Stokes flow (no momentum solve).
      void solve(Field& h, const Field& u_u, const Field& u_v,
                 const Field& evap_rate,
                 BubbleLinearSystem& sys, const StereoMesh& mesh,
                 const Mask& mask, const BubbleConfig& cfg, double dt);
  }
  ```

  Assembly:
  1. `sys.reset()`
  2. `SurfaceFVM::ddt(sys, h, dt, mesh, mask)` — transient
  3. Compute face fluxes: $F_e = (u_u)_e \cdot L_{face,e}$ (face velocity $\times$ face length, with $h$ interpolated via TVD)
  4. `SurfaceFVM::divergence(sys, flux_u, flux_v, h, scheme, mesh, mask)` — convection
  5. `SurfaceFVM::add_source(sys, evap_source, mesh, mask)` — evaporation (zero in Phase 3)
  6. Enforce rim BC: $h = h_{rim}$ at rim cells (large-diagonal penalty or direct injection)
  7. Enforce mask: identity equation for blanked cells
  8. `sys.solve(h)`
  9. Clamp: $h = \max(h, h_{min})$

- [ ] **3.4** Implement a simple gravity-driven velocity for testing:

  In the thin-film limit with no inertia and no pressure gradient, gravity drainage gives:

  $$\mathbf{u}_{gravity} = \frac{\rho_l h^2}{3\mu}\,\mathbf{g}_s$$

  This is the lubrication-theory drainage velocity. Compute $u_u(i,j)$ and $u_v(i,j)$ from the precomputed `g_surface` and current $h$.

- [ ] **3.5** Wire up in `bubble_main.cpp`:
  ```
  Time loop:
    1. Compute gravity-driven velocity from current h
    2. Compute face fluxes
    3. Solve thickness transport (no evaporation)
    4. Enforce BCs and clamping
    5. Output
    6. Store old, advance
  ```

### Verification (Phase 3)
- [ ] **V3.1 Mass conservation**: With $J = 0$ and $\mathbf{u}_s = 0$, thickness must remain at $h_{initial}$ everywhere. Check $\max|h - h_{initial}| < 10^{-12}$.
- [ ] **V3.2 Gravity drainage**: Enable gravity. The film should thin at the top (near north pole) and thicken at the bottom (near equator). Verify:
  - Total mass $M = \sum h(i,j) \cdot A(i,j)$ is conserved (no evaporation).
  - $h$ at north pole decreases monotonically.
  - $h$ near equator increases (pooling effect).
  Visualise $h$ field on hemisphere in ParaView at several timesteps.
- [ ] **V3.3 TVD boundedness**: $h$ must remain $\ge h_{min}$ everywhere. No overshoots beyond $h_{initial}$ in the first few timesteps.
- [ ] **V3.4 Grid convergence**: Run with $32^2$, $64^2$, $128^2$. The drainage pattern should converge. Total mass conservation error should decrease with resolution.

---

## Phase 4: 2D Navier-Stokes Solver (SIMPLE)

**Goal**: Implement the full 2D surface Navier-Stokes with SIMPLE pressure-velocity coupling. Replace the prescribed drainage velocity from Phase 3 with the solved velocity field.

### Tasks

- [ ] **4.1** Create `src/momentum.hpp/cpp`:

  ```cpp
  namespace Momentum {
      /// Solve one SIMPLE iteration for velocity components u_u, u_v and pressure p.
      /// Returns residuals for convergence check.
      struct SimpleResiduals {
          double u_residual, v_residual, p_residual, continuity_residual;
      };

      SimpleResiduals solve_simple_iteration(
          Field& u_u, Field& u_v, Field& p,
          const Field& h,
          const Field& F_marangoni_u, const Field& F_marangoni_v,  // Body forces
          const Field& disjoining_grad_u, const Field& disjoining_grad_v,
          BubbleLinearSystem& sys_u, BubbleLinearSystem& sys_v, BubbleLinearSystem& sys_p,
          const StereoMesh& mesh, const Mask& mask, const BubbleConfig& cfg, double dt);

      /// Full SIMPLE loop: iterate until convergence or max_iters.
      void solve(Field& u_u, Field& u_v, Field& p,
                 const Field& h,
                 const Field& F_marangoni_u, const Field& F_marangoni_v,
                 const Field& disjoining_grad_u, const Field& disjoining_grad_v,
                 BubbleLinearSystem& sys_u, BubbleLinearSystem& sys_v, BubbleLinearSystem& sys_p,
                 const StereoMesh& mesh, const Mask& mask, const BubbleConfig& cfg, double dt);
  }
  ```

  **SIMPLE iteration** (each call to `solve_simple_iteration`):

  **Step A — Momentum predictor** (solve for $u^*_u$, $u^*_v$ separately):

  For the $u$-momentum equation:
  $$\rho_l\frac{u_u - u_u^n}{\Delta t} + \rho_l\,\nabla_s \cdot (\mathbf{u}\,u_u) = -\frac{\partial p^n}{\partial u}\frac{1}{\lambda^2} + \frac{\mu}{h}\nabla_s^2(h\,u_u) + F_u$$

  FVM assembly:
  1. `sys_u.reset()`
  2. `SurfaceFVM::ddt(sys_u, u_u, dt, mesh, mask)` — transient (weighted by $\rho_l$)
  3. `SurfaceFVM::divergence(sys_u, mass_flux_u, mass_flux_v, u_u, scheme, mesh, mask)` — convection (weighted by $\rho_l$)
  4. `SurfaceFVM::laplacian(sys_u, mu_field, mesh, mask)` — diffusion (negative sign: goes to LHS)
  5. Add pressure gradient source: $\text{source}(i,j) += -(\partial p / \partial u)_{face} \cdot A(i,j)$
  6. Add body forces: gravity, Marangoni, disjoining, drag
  7. Under-relax: $a_P \to a_P / \alpha_u$, $\text{source} += (1-\alpha_u)/\alpha_u \cdot a_P \cdot u_u^n$
  8. Enforce rim no-slip (penalty)
  9. Enforce mask
  10. `sys_u.solve(u_u)`

  Repeat analogously for $u_v$.

  **Step B — Rhie-Chow face velocities**:

  At each face $f$ between cells $P$ and $N$:
  $$u_f = \overline{u^*}_f - \hat{d}_f\left(\left(\nabla p\right)_f - \overline{\nabla p}_f\right)$$

  where:
  - $\overline{u^*}_f$ = linear interpolation of cell-centre velocities to face
  - $\hat{d}_f = \overline{A(i,j) / a_P}_f$ = interpolated reciprocal diagonal coefficient
  - $(\nabla p)_f$ = compact face gradient: $(p_N - p_P) / d_{PN}$
  - $\overline{\nabla p}_f$ = interpolated cell-centre gradient

  Compute mass flux: $\dot{m}_f = \rho_l \cdot u_{f,\perp} \cdot L_{face}$

  **Step C — Pressure correction equation**:

  $$\nabla_s \cdot \left(\hat{d}_f \cdot L_{face}^2 / d_{PN}\,\nabla_s p'\right) = \nabla_s \cdot (\dot{m}^*_f)$$

  This is a Poisson-type equation. Assemble using `SurfaceFVM::laplacian` with $\gamma = \hat{d}_f$ and the mass imbalance as source.

  Fix pressure at one cell (e.g., the north pole) to remove the null space of the Poisson equation. Or use the Laplace pressure $p_{ref} = 2\sigma/R$ at the north pole.

  **Step D — Corrections**:
  - Face velocity: $\dot{m}_f^{n+1} = \dot{m}^*_f - \hat{d}_f \cdot L_{face} / d_{PN} \cdot (p'_N - p'_P)$
  - Cell velocity: $u_u^{n+1} = u^*_u - (A/a_P) \cdot (\partial p'/\partial u)$
  - Pressure: $p^{n+1} = p^n + \alpha_p\,p'$

- [ ] **4.2** Store momentum diagonal $a_P$ after predictor solve for use in Rhie-Chow and pressure correction.

  Add a method to `BubbleLinearSystem`:
  ```cpp
  void extract_diagonal(Field& a_P_field) const;
  ```

- [ ] **4.3** Implement surface gradient operator for pressure:
  ```cpp
  namespace SurfaceFVM {
      /// Compute cell-centre gradient of phi in surface coordinates.
      /// Returns (∂φ/∂u, ∂φ/∂v) scaled by 1/λ² (physical gradient).
      void gradient(const Field& phi, Field& grad_u, Field& grad_v,
                    const StereoMesh& mesh, const Mask& mask);
  }
  ```
  Use Green-Gauss gradient: $\nabla \phi|_P = \frac{1}{A_P} \sum_{faces} \phi_f \cdot \hat{n}_f \cdot L_f$

- [ ] **4.4** Wire up SIMPLE in `bubble_main.cpp`. Replace the prescribed gravity velocity with the SIMPLE solve. The time loop becomes:
  ```
  1. Compute J (zero for now)
  2. Solve thickness transport with current u_u, u_v
  3. Update σ(T), Marangoni (zero for now)
  4. SIMPLE loop for u_u, u_v, p
  5. Check residuals, adapt dt
  6. Output, advance
  ```

- [ ] **4.5** Add pressure reference point. Since the momentum equation is driven by $\nabla p$ (not $p$ itself), the pressure Poisson has a null space. Fix one cell (preferably the north pole, or use mean pressure = Laplace pressure $2\sigma/R$).

### Verification (Phase 4)
- [ ] **V4.1 Stokes flow**: Disable convection (set $\rho_l \to 0$ in convective term or use very low Re). With gravity only, verify the drainage velocity matches the thin-film lubrication solution from Phase 3 to within discretisation error.
- [ ] **V4.2 No-flow equilibrium**: With gravity off and uniform initial conditions, the velocity should remain zero and pressure uniform. Check $\max|u| < 10^{-12}$.
- [ ] **V4.3 SIMPLE convergence**: Monitor SIMPLE residuals per outer iteration. They should decrease monotonically. Log and plot.
- [ ] **V4.4 Mass conservation**: $\sum \dot{m}_f = 0$ for each cell after pressure correction (continuity). Check global mass imbalance $< 10^{-8}$.
- [ ] **V4.5 Parallel consistency**: Results with 1, 2, 4 ranks should be identical (within solver tolerance).

---

## Phase 5: Thermodynamics & Evaporation

**Goal**: Add the energy equation and Sultan evaporation model. Couple temperature to Marangoni forces in the momentum equation.

### Tasks

- [ ] **5.1** Create `src/energy.hpp/cpp`:

  ```cpp
  namespace Energy {
      /// Solve energy equation for new T.
      void solve(Field& T, const Field& h, const Field& u_u, const Field& u_v,
                 const Field& evap_rate,
                 BubbleLinearSystem& sys, const StereoMesh& mesh,
                 const Mask& mask, const BubbleConfig& cfg, double dt);
  }
  ```

  Assembly:
  1. `sys.reset()`
  2. `SurfaceFVM::ddt_weighted(sys, T, rho_cp_h, dt, mesh, mask)` — transient ($\rho c_p h \cdot \partial T/\partial t$)
     where `rho_cp_h(i,j) = rho_l * c_p * h(i,j)`
  3. Convection: `SurfaceFVM::divergence(sys, heat_flux_u, heat_flux_v, T, scheme, mesh, mask)`
     where `heat_flux_u(i,j) = rho_l * c_p * h_face * u_u(i,j) * L_face_u`
  4. Diffusion: `SurfaceFVM::laplacian(sys, k_h, mesh, mask)` where `k_h(i,j) = k_l * h(i,j)`
  5. Evaporative cooling source: $S_{evap}(i,j) = -J(i,j) \cdot \Delta H_v$
  6. Convective heat loss: $S_{conv}(i,j) = h_{conv}(T_{ambient} - T(i,j))$
     Split: implicit part $-h_{conv}$ to $a_P$, explicit part $h_{conv} \cdot T_{ambient}$ to source.
  7. `SurfaceFVM::add_source(sys, S_explicit, mesh, mask)`
  8. Enforce rim BC: $T = T_{rim}$ at rim cells
  9. Enforce mask
  10. `sys.solve(T)`

- [ ] **5.2** Implement the **Sultan evaporation model** in a standalone module:

  Create `src/evaporation.hpp/cpp`:
  ```cpp
  namespace Evaporation {
      /// Compute evaporation rate J(T) at each cell [kg/(m²·s)].
      void compute_sultan(Field& J, const Field& T,
                          const StereoMesh& mesh, const Mask& mask,
                          const BubbleConfig& cfg);

      /// Saturation mass fraction via Clausius-Clapeyron approximation.
      double omega_sat(double T);
  }
  ```

  Antoine equation for water saturation pressure:
  $$\log_{10}(p_{sat}/\text{Pa}) = A - B/(C + T/\text{K})$$
  with $A = 8.07131$, $B = 1730.63$, $C = 233.426$ (valid 1-100 C).

  Then $\omega_{sat} = 0.622 \cdot p_{sat} / (p_{atm} - p_{sat})$ (mixing ratio to mass fraction conversion).

  Sultan model:
  $$J(i,j) = \frac{\rho_{gas} \cdot D_{AB} \cdot (\omega_{sat}(T(i,j)) - \omega_\infty)}{R_c}$$

  where $R_c = R_{bubble}$ (characteristic diffusion length).

- [ ] **5.3** Create `src/marangoni.hpp/cpp`:

  ```cpp
  namespace Marangoni {
      /// Compute surface tension from temperature.
      void compute_sigma(Field& sigma, const Field& T, const BubbleConfig& cfg);
      // σ(i,j) = σ_0 + γ_T * (T(i,j) - T_ref)

      /// Compute Marangoni body force = -γ_T * ∇_s T.
      void compute_force(Field& F_u, Field& F_v, const Field& T,
                         const StereoMesh& mesh, const Mask& mask,
                         const BubbleConfig& cfg);
  }
  ```

  Uses `SurfaceFVM::gradient(T, dT_du, dT_dv, mesh, mask)` then:
  $$F_u(i,j) = -\gamma_T \cdot dT\_du(i,j), \qquad F_v(i,j) = -\gamma_T \cdot dT\_dv(i,j)$$

- [ ] **5.4** Wire up in `bubble_main.cpp`. The time loop now becomes:
  ```
  1. Compute J from current T (Sultan model) — if enabled
  2. Solve thickness transport with current u, J
  3. Solve energy equation for new T — if enabled
  4. Update σ(T), compute Marangoni force — if enabled
  5. SIMPLE loop for u_u, u_v, p (now includes Marangoni body force)
  6. Check residuals, adapt dt
  7. Check rupture
  8. Output, advance
  ```

- [ ] **5.5** Add adaptive timestep based on Marangoni stability:
  $$\Delta t_{Marangoni} = C_M \cdot \min_{cells} \frac{\rho_l \cdot A_{cell}}{|\gamma_T| \cdot |\nabla_s T|}$$

### Verification (Phase 5)
- [ ] **V5.1 Isothermal baseline**: With `enable_evaporation = false` and `enable_marangoni = false`, results should match Phase 4 exactly.
- [ ] **V5.2 Evaporation only**: Disable Marangoni, enable evaporation. With no flow ($\mathbf{u} = 0$), $h$ should decrease uniformly at rate $J/\rho_l$. Verify against analytic:
  $$h(t) = h_0 - \frac{J}{\rho_l} t$$
- [ ] **V5.3 Marangoni flow**: Apply a temperature gradient (e.g., hot at top, cold at bottom). Verify flow from hot to cold. The flow direction should be consistent with $\gamma_T < 0$ (water).
- [ ] **V5.4 Energy conservation**: With no evaporation and no convective loss, total thermal energy $\sum \rho c_p h T \cdot A$ should be conserved. Check to within solver tolerance.
- [ ] **V5.5 Coupled test**: Enable all physics. The top of the bubble should cool (evaporation), drive Marangoni flow downward, and thin. Visualise the coupled fields in ParaView.

---

## Phase 6: Disjoining Pressure & Rupture Detection

**Goal**: Add van der Waals disjoining pressure and simulation halting when the film ruptures.

### Tasks

- [ ] **6.1** Create `src/disjoining.hpp/cpp`:

  ```cpp
  namespace Disjoining {
      /// Compute disjoining pressure Π(h) at each cell.
      void compute(Field& Pi, const Field& h, const BubbleConfig& cfg);
      // Π(i,j) = A_hamaker / (6π max(h(i,j), h_min)³)

      /// Compute ∇_s Π for use as body force in momentum.
      /// This is equivalent to Π'(h) * ∇_s h.
      void compute_gradient(Field& dPi_du, Field& dPi_dv,
                            const Field& h, const StereoMesh& mesh,
                            const Mask& mask, const BubbleConfig& cfg);

      /// Check if any cell has h < h_rupture. Returns rupture info.
      struct RuptureInfo {
          bool ruptured = false;
          double h_min;             // Global minimum h
          double u_loc, v_loc;      // Location of minimum
          double x_loc, y_loc, z_loc; // 3D coordinates
      };
      RuptureInfo check_rupture(const Field& h, const StereoMesh& mesh,
                                const Mask& mask, const BubbleConfig& cfg);
  }
  ```

  **Disjoining pressure**:
  $$\Pi(h) = \frac{A_{Hamaker}}{6\pi\,[\max(h, h_{min})]^3}$$

  **Gradient** (for momentum source):
  $$\nabla_s \Pi = \Pi'(h)\,\nabla_s h = -\frac{A_{Hamaker}}{2\pi\,h^4}\,\nabla_s h$$

  Use `SurfaceFVM::gradient(h, dh_du, dh_dv, ...)` then multiply.

- [ ] **6.2** Add disjoining pressure to momentum predictor in `momentum.cpp`:
  - Disjoining gradient acts as a body force alongside Marangoni and gravity.
  - Alternatively, add $\Pi(h)$ directly to the pressure field: $p_{eff} = p - \Pi(h)$, then use $\nabla_s p_{eff}$ in momentum.

  **Recommended approach**: Add $\nabla_s \Pi$ as an explicit body force. This is simpler and avoids modifying the pressure Poisson.

- [ ] **6.3** Add disjoining stability constraint to adaptive timestep:
  $$\Delta t_{disjoining} = C_\Pi \cdot \min_{cells} \frac{\rho_l \cdot (\sqrt{A_{cell}})^2}{\left|\Pi'(h)\right| \cdot h}$$

- [ ] **6.4** Implement rupture detection in the time loop:
  ```cpp
  auto info = Disjoining::check_rupture(h, mesh, mask, cfg);
  if (info.ruptured) {
      log("RUPTURE detected at t = " + to_string(t), Color::RED);
      log("  h_min = " + to_string(info.h_min) + " at (" +
          to_string(info.x_loc) + ", " + to_string(info.y_loc) + ", " +
          to_string(info.z_loc) + ")", Color::RED);
      // Write final output
      IO::write_timestep(t, step, mesh, fields, cfg);
      break;  // Exit time loop
  }
  ```

  The `check_rupture` function uses `MPI_Allreduce` with `MPI_MINLOC` to find the global minimum $h$ and its location.

- [ ] **6.5** Final integration in `bubble_main.cpp`. The complete time loop is now:
  ```
  for each timestep:
      1. Compute J from T (Sultan model)
      2. Solve continuity for h
      3. Solve energy for T
      4. Compute σ(T), Marangoni force
      5. Compute disjoining Π(h), ∇_s Π
      6. SIMPLE loop (momentum + pressure with all body forces)
      7. Adaptive Δt = min(CFL, Marangoni, disjoining, dt_max)
      8. Check rupture → halt if triggered
      9. Output if interval reached
     10. Store old fields, advance t
  ```

### Verification (Phase 6)
- [ ] **V6.1 Disjoining only**: Set up a film with a thin spot ($h = 50\text{ nm}$ in a small patch, $h = 1\,\mu\text{m}$ elsewhere). Disable evaporation/Marangoni/gravity. The disjoining pressure should accelerate thinning of the thin spot. Verify the spot thins further while surrounding film barely changes.
- [ ] **V6.2 Rupture detection**: Continue V6.1 until $h < h_{rupture}$. Confirm the simulation halts with correct rupture location.
- [ ] **V6.3 Stability**: Verify the adaptive timestep correctly limits $\Delta t$ as the film thins. Plot $\Delta t$ vs. time — it should decrease as $h_{min}$ decreases.
- [ ] **V6.4 Full physics test**: Enable all physics. Run a hemispherical bubble with:
  - $R = 1\text{ cm}$, $h_0 = 1\,\mu\text{m}$
  - Gravity pointing down
  - Ambient temperature 20 C, 50% humidity
  - Default fluid properties (water + surfactant)

  Expected behaviour:
  - Gravitational drainage thins the top.
  - Evaporative cooling at the top creates a temperature gradient.
  - Marangoni flow (cold top → hot bottom) partially opposes drainage.
  - Disjoining pressure accelerates final thinning.
  - Rupture occurs at the top of the bubble after some time.

  Produce a ParaView visualisation of $h$, $T$, $|\mathbf{u}_s|$, and $p$ fields on the hemisphere at several timesteps leading up to rupture.
- [ ] **V6.5 Matrix condition check**: Monitor PETSc KSP iteration counts throughout the simulation. Flag if iterations spike (indicating ill-conditioning from thin-film disjoining forces). If problematic, document and propose preconditioner improvements.

---

# Part IV — Physical Constants Reference

For convenience, collect all physical constants and their default values:

| Symbol | Name | Default Value | Unit |
|--------|------|---------------|------|
| $R$ | Bubble radius | 0.01 | m |
| $h_0$ | Initial thickness | $10^{-6}$ | m |
| $\rho_l$ | Liquid density | 1000 | kg/m$^3$ |
| $\mu$ | Dynamic viscosity | $10^{-3}$ | Pa$\cdot$s |
| $\sigma_0$ | Surface tension (ref) | 0.025 | N/m |
| $\gamma_T$ | $d\sigma/dT$ | $-1.7 \times 10^{-4}$ | N/(m$\cdot$K) |
| $k_l$ | Thermal conductivity | 0.6 | W/(m$\cdot$K) |
| $c_p$ | Specific heat | 4186 | J/(kg$\cdot$K) |
| $\rho_{gas}$ | Air density | 1.2 | kg/m$^3$ |
| $D_{AB}$ | Diffusion coeff (H$_2$O in air) | $2.5 \times 10^{-5}$ | m$^2$/s |
| $\omega_\infty$ | Ambient humidity | 0.01 | — |
| $T_{ref}$ | Reference temperature | 293.15 | K |
| $\Delta H_v$ | Latent heat | $2.26 \times 10^6$ | J/kg |
| $h_{conv}$ | Convective HTC | 10 | W/(m$^2\cdot$K) |
| $A_{Hamaker}$ | Hamaker constant | $10^{-20}$ | J |
| $g$ | Gravitational acceleration | 9.81 | m/s$^2$ |
| $h_{rupture}$ | Rupture threshold | $10^{-8}$ | m |

---

# Part V — Documentation Requirements

After each phase, update the following documents:

1. **`CHANGES.md`**: Append a dated entry following the existing format (see current entries for style). List new files, modified files, and the rationale for each change.

2. **`PHYSICS_BUBBLE.md`**: Append the mathematical derivations and FVM discretisation for the equations implemented in that phase. Follow the style of the existing `PHYSICS.md` (numbered sections, LaTeX equations, stencil diagrams, sign convention notes).

3. **`PLAN_BUBBLE.md`** (this file): Mark completed tasks with `[x]`.

4. **`README.md`**: After Phase 2 (first runnable output), add a "Bubble Solver" section with build and run instructions.

---

# Part VI — Risk Register

| Risk | Mitigation |
|------|-----------|
| Stereographic projection singularity at $r_{2D} \to \infty$ | Circular mask excludes cells beyond $R_{mask}$; never evaluate near south pole |
| Disjoining pressure stiffness ($h^{-3}$) | Clamp $h \ge h_{min}$ in $\Pi$; adaptive timestep; implicit treatment if needed |
| SIMPLE non-convergence on curved surface | Start with low CFL and many SIMPLE iters; under-relax aggressively ($\alpha_p = 0.1$) initially |
| Checkerboard pressure oscillations | Rhie-Chow interpolation is mandatory; verify no odd-even decoupling |
| Metric singularity at conformal factor extremes | The conformal factor $\lambda$ is bounded for $r_{2D} \le R$; no singularity in the hemispherical domain |
| PETSc solver breakdown for ill-conditioned systems | Monitor iteration count; switch to GMRES with ILU if BiCGStab stalls |
| Evaporation driving $h$ negative | TVD + non-negative clamping + implicit evaporation term (linearise $J$ w.r.t. $h$ if needed) |
