# Development Plan

**pancake** is a multiphysics thin-film solver for a hemispherical soap bubble. The film thickness $h(\mathbf{x}, t)$ evolves on a stereographically projected hemisphere under gravitational drainage, thermal Marangoni convection, evaporation (Sultan diffusion-limited model), and disjoining pressure (van der Waals).

---

## Architecture Overview

**Domain**: Inverse stereographic projection $(u, v) \to (x, y, z)$ on a hemisphere of radius $R$. The 2D grid has uniform spacing in $(u, v)$; cell geometry is non-uniform and precomputed in 3D.

**Module map (current state):**
```
src/
  bubble_config.hpp/cpp         — BubbleConfig + ConvectionScheme
  stereo_mesh.hpp/cpp           — Inverse stereographic grid + geometry engine
  mask.hpp/cpp                  — Circular domain mask + ghost blanking
  communicator.hpp/cpp          — MPI ghost exchange (open, non-periodic)
  bubble_io.hpp/cpp             — VTK unstructured output (vtu/pvtu/pvd)
  field.hpp/cpp                 — Scalar fields with ghost layers
  bubble_linear_system.hpp/cpp  — PETSc 5-point sparse system
  surface_operators.hpp/cpp     — ddt, laplacian, divergence, add_source, gradient
  thickness_transport.hpp/cpp   — Gravity velocity, face fluxes, h transport
  utils.hpp                     — Rank-0 coloured logger
  bubble_main.cpp               — Entry point and time loop

  [Planned]
  momentum.hpp/cpp              — 2D N-S with SIMPLE
  energy.hpp/cpp                — Energy equation
  evaporation.hpp/cpp           — Sultan evaporation model
  marangoni.hpp/cpp             — σ(T) model + Marangoni body force
  disjoining.hpp/cpp            — van der Waals Π(h) + rupture detection
```

**Solver loop (per timestep):**
```
1.  Store old-time fields
2.  Compute evaporation flux J(i,j) from current T  [Phase 5]
3.  Solve continuity for h:   ∂h/∂t + ∇_s·(h u_s) = -J/ρ_l
4.  Solve energy for T:       ρ c_p h (∂T/∂t + u_s·∇_s T) = ∇_s·(k h ∇_s T) - J ΔH_v + h_conv(T_amb - T)  [Phase 5]
5.  Update σ(T), Marangoni force F_M = -γ_T ∇_s T  [Phase 5]
6.  Compute disjoining Π(h), ∇_s Π  [Phase 6]
7.  SIMPLE loop: momentum predictor → Rhie-Chow → pressure correction → correct  [Phase 4]
8.  Adaptive Δt = min(CFL, Marangoni, disjoining, dt_max)
9.  Check rupture: if min(h) < h_rupture → halt  [Phase 6]
10. Output VTK if write interval reached
```

---

## Phase 1 — Infrastructure & Configuration

**Status: complete** (2026-04-16).

- [x] **1.1** Create `src/bubble_config.hpp/cpp` with `BubbleConfig` struct and `load_from_file()`.
- [x] **1.2** Set up `CMakeLists.txt` with `pancake_bubble` as the sole target.
- [x] **1.3** Create `bubble_main.cpp` skeleton: PETSc init, config load, empty time loop.
- [x] **1.4** Create `tests/CMakeLists.txt` with `add_bubble_test()` helper.

---

## Phase 2 — Inverse Stereographic Mesh & VTK Output

**Goal**: Inverse stereographic mesh, circular mask, geometry engine, VTK unstructured output. No physics — validate geometry and I/O.

**Status: complete** (2026-04-16).

- [x] **2.1** Create `src/stereo_mesh.hpp/cpp`.

  Inverse stereographic map — $(u,v) \to (x,y,z)$ on hemisphere of radius $R$:
  $$x = \frac{2R^2 u}{r_{2D}^2 + R^2}, \quad y = \frac{2R^2 v}{r_{2D}^2 + R^2}, \quad z = R\,\frac{R^2 - r_{2D}^2}{r_{2D}^2 + R^2}$$

  Per-cell geometry precomputed: cell area (cross product), face chord lengths, centre-to-centre distances, surface normal, surface gravity components.

- [x] **2.2** Create `src/mask.hpp/cpp` — circular domain mask.

  Cells with $r_{2D} > R$ are blanked. Rim cells are active cells adjacent to at least one blanked neighbour.

- [x] **2.3** Create `src/bubble_io.hpp/cpp` — VTK unstructured (`.vtu/pvtu/pvd`) output.

- [x] **2.4** Adapt `src/communicator.hpp/cpp` — non-periodic ghost exchange; domain edges get zero-gradient ghost fill.

- [x] **2.5** Write initial VTK output and area verification in `bubble_main.cpp`.

- [x] **2.6** Tests: `test_stereo_mesh` (11 checks), `test_mask` (5 checks).

---

## Phase 3 — Thickness Transport & Gravitational Drainage

**Goal**: Solve $\partial h/\partial t + \nabla_s \cdot (h \mathbf{u}_s) = 0$ with prescribed gravity-driven Stokes velocity.

**Status: complete** (2026-04-16).

- [x] **3.1** Create `src/bubble_linear_system.hpp/cpp` — PETSc 5-point sparse system with `apply_mask()` (row replacement for masked and rim cells).

- [x] **3.2** Create `src/surface_operators.hpp/cpp` — metric-aware FVM operators.

  Laplacian face coefficient: $c_e = \gamma_e \cdot L_{face,e} / d_{PN,e}$ (harmonic $\gamma$, 3D geometry).

  `ConvectionScheme` moved to `bubble_config.hpp`; no dependency on bearing-era headers.

- [x] **3.3** Implement gravity-driven velocity (thin-film Stokes):
  $$\mathbf{u}_s = \frac{\rho_l h^2}{3\mu}\,\mathbf{g}_s$$

- [x] **3.4** Create `src/thickness_transport.hpp/cpp`: `compute_gravity_velocity`, `compute_face_fluxes`, `solve` (ddt + divergence + optional evap source + rim BC + mask + PETSc solve + clamp $h \ge h_{min}$).

- [x] **3.5** Wire up in `bubble_main.cpp`.

- [x] **3.6** Tests: `test_surface_operators` (6 checks), `test_thickness_transport` (6 checks).

---

## Phase 4 — 2D Navier-Stokes (SIMPLE)

**Goal**: Full 2D surface Navier-Stokes with SIMPLE pressure-velocity coupling. Replaces prescribed drainage velocity with a solved velocity field.

**Status: pending**

- [ ] **4.1** Create `src/momentum.hpp/cpp`.

  **SIMPLE iteration:**

  *Step A — Momentum predictor* (solved separately for $u_u$ and $u_v$):
  $$\rho_l\frac{u_u - u_u^n}{\Delta t} + \rho_l\,\nabla_s \cdot (\mathbf{u}\,u_u) = -\frac{1}{\lambda^2}\frac{\partial p^n}{\partial u} + \frac{\mu}{h}\nabla_s^2(h\,u_u) + F_u$$

  Under-relaxation: $a_P \to a_P/\alpha_u$, source $\mathrel{+}= (1-\alpha_u)/\alpha_u \cdot a_P \cdot u_u^n$.

  *Step B — Rhie-Chow face velocities*:
  $$u_f = \bar{u}^*_f - \hat{d}_f\left[(\nabla p)_f - \overline{(\nabla p)}_f\right]$$

  *Step C — Pressure correction*:
  $$\nabla_s \cdot \!\left(\hat{d}_f \cdot \frac{L_{face}^2}{d_{PN}}\nabla_s p'\right) = \nabla_s \cdot \dot{m}^*_f$$

  *Step D — Corrections*: face fluxes, cell velocities, $p \mathrel{+}= \alpha_p p'$.

- [ ] **4.2** Add `extract_diagonal(Field& a_P_field)` to `BubbleLinearSystem` for Rhie-Chow.

- [ ] **4.3** Pin pressure at north pole to Laplace pressure $p_{ref} = 2\sigma/R$ (remove null space).

- [ ] **4.4** Wire SIMPLE loop in `bubble_main.cpp`.

**Verification:**
- No-flow equilibrium: uniform IC + no gravity → $|\mathbf{u}| < 10^{-12}$.
- Stokes limit: low-Re gravity drainage matches Phase 3 Stokes velocity.

---

## Phase 5 — Thermodynamics & Evaporation

**Status: pending**

- [ ] **5.1** Create `src/evaporation.hpp/cpp` — Sultan diffusion-limited model:
  $$J(T) = \frac{\rho_{gas}\,D_{AB}\,(\omega_{sat}(T) - \omega_\infty)}{R_c}$$

- [ ] **5.2** Create `src/energy.hpp/cpp`:
  $$\rho_l c_p h\left(\frac{\partial T}{\partial t} + \mathbf{u}_s \cdot \nabla_s T\right) = \nabla_s \cdot (k_l h\,\nabla_s T) - J\,\Delta H_v + h_{conv}(T_{amb} - T)$$

- [ ] **5.3** Create `src/marangoni.hpp/cpp`: $\mathbf{F}_M = -\gamma_T\,\nabla_s T$.

- [ ] **5.4** Adaptive Marangoni timestep.

- [ ] **5.5** Wire up in `bubble_main.cpp`.

---

## Phase 6 — Disjoining Pressure & Rupture Detection

**Status: pending**

- [ ] **6.1** Create `src/disjoining.hpp/cpp` — van der Waals:
  $$\Pi(h) = \frac{A_{Hamaker}}{6\pi\,h^3}, \qquad \nabla_s \Pi = -\frac{A_{Hamaker}}{2\pi\,h^4}\,\nabla_s h$$

- [ ] **6.2** Add $\nabla_s \Pi$ as explicit body force in momentum predictor.

- [ ] **6.3** Disjoining stability timestep constraint.

- [ ] **6.4** Rupture detection: `MPI_Allreduce` with `MPI_MINLOC` for global minimum $h$.

- [ ] **6.5** Full physics integration in `bubble_main.cpp`.

---

## Phase Summary

| Phase | Description | Status |
|-------|-------------|--------|
| B1 | Infrastructure & configuration | Complete |
| B2 | Stereographic mesh & VTK I/O | Complete |
| B3 | Thickness transport + gravity drainage | Complete |
| B4 | 2D Navier-Stokes (SIMPLE) | Pending |
| B5 | Thermodynamics & evaporation | Pending |
| B6 | Disjoining pressure & rupture | Pending |

---

## Physical Constants Reference

| Symbol | Name | Default | Unit |
|--------|------|---------|------|
| $R$ | Bubble radius | 0.01 | m |
| $h_0$ | Initial thickness | $10^{-6}$ | m |
| $\rho_l$ | Liquid density | 1000 | kg/m$^3$ |
| $\mu$ | Dynamic viscosity | $10^{-3}$ | Pa·s |
| $\sigma_0$ | Surface tension (ref) | 0.025 | N/m |
| $\gamma_T$ | $d\sigma/dT$ | $-1.7 \times 10^{-4}$ | N/(m·K) |
| $k_l$ | Thermal conductivity | 0.6 | W/(m·K) |
| $c_p$ | Specific heat | 4186 | J/(kg·K) |
| $\rho_{gas}$ | Air density | 1.2 | kg/m$^3$ |
| $D_{AB}$ | Diffusion coeff (H$_2$O in air) | $2.5 \times 10^{-5}$ | m$^2$/s |
| $\omega_\infty$ | Ambient humidity | 0.01 | — |
| $T_{ref}$ | Reference temperature | 293.15 | K |
| $\Delta H_v$ | Latent heat | $2.26 \times 10^6$ | J/kg |
| $h_{conv}$ | Convective HTC | 10 | W/(m$^2$·K) |
| $A_{Hamaker}$ | Hamaker constant | $10^{-20}$ | J |
| $g$ | Gravitational acceleration | 9.81 | m/s$^2$ |
| $h_{rupture}$ | Rupture threshold | $10^{-8}$ | m |

---

## Risk Register

| Risk | Mitigation |
|------|-----------|
| Stereographic singularity at $r_{2D} \to \infty$ | Circular mask excludes cells beyond $R$ |
| Disjoining pressure stiffness ($h^{-3}$) | Clamp $h \ge h_{min}$ in $\Pi$; adaptive timestep |
| SIMPLE non-convergence on curved surface | Low CFL; aggressive under-relaxation ($\alpha_p = 0.1$) initially |
| Checkerboard pressure oscillations | Rhie-Chow face interpolation mandatory |
| Pressure null space (closed surface) | Pin pressure at north pole to Laplace pressure $2\sigma/R$ |
| PETSc solver breakdown (ill-conditioned) | Monitor KSP count; switch to GMRES+ILU if BiCGStab stalls |
| Evaporation driving $h$ negative | TVD + non-negative clamping; implicit linearisation of $J(h)$ if needed |
