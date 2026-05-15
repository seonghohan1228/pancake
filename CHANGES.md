# Change Log

---

## Phase B3 — Surface FVM operators, thickness transport, gravity drainage (2026-04-16)

### New files

**`src/bubble_linear_system.hpp/cpp`**

PETSc-backed 5-point sparse linear system for the stereographic surface mesh.
- Takes `StereoMesh` instead of any bearing-era mesh type.
- Global row index: `(offset_u + i) * n_v_global + j` (no periodic wrap).
- `apply_mask()`: overwrites masked rows with identity (a_P=1, rhs=0) and rim rows with identity (a_P=1, rhs=bc).
- `solve()`: BiCGStab + block Jacobi; E/W connectivity includes bounds check for non-periodic domain.

**`src/surface_operators.hpp/cpp`**

Metric-aware FVM operators for the curved surface. All use 3D chord lengths and areas from `StereoMesh`. Implements: `ddt`, `ddt_weighted`, `laplacian`, `divergence` (upwind + TVD deferred correction), `add_source`, `gradient`. The `active_at(mesh, gu, j)` helper analytically checks mask membership for any global cell index, enabling correct masked-neighbour detection without MPI communication.

`ConvectionScheme` moved from the deleted `fvm.hpp` into `bubble_config.hpp`.

**`src/thickness_transport.hpp/cpp`**

- `compute_gravity_velocity`: thin-film Stokes velocity $\mathbf{u}_s = \rho_l h^2 / (3\mu) \cdot \mathbf{g}_s$.
- `compute_face_fluxes`: face-centred flux by simple averaging of adjacent cell velocities.
- `solve`: full FVM assembly (ddt + divergence + optional evap source) + rim BC + mask + PETSc solve + h_min clamp.

**`tests/test_surface_operators.cpp`** — 6 checks: ddt a_P, ddt source, ddt_weighted a_P, laplacian interior conservation, add_source, divergence with zero flux.

**`tests/test_thickness_transport.cpp`** — 6 checks: zero-gravity → zero velocity, pole velocity zero, off-pole velocity nonzero, face flux for uniform velocity, zero-flux solve preserves h, h clamped to h_min.

### Modified files

**`src/communicator.cpp`** — bug fix: swapped send buffers in `MPI_Sendrecv` calls.

The original code sent the left physical strip to the left neighbor and the right strip to the right neighbor in each Sendrecv — meaning rank 0 sent to rank -1 and each rank sent the wrong boundary. Correct pattern: send right strip to `right_rank`, send left strip to `left_rank`, with matching tags. Ghost values are now correct under MPI.

**`src/bubble_main.cpp`** — Phase 3 wiring: `store_old_time` → `compute_gravity_velocity` → `update_ghosts` → `compute_face_fluxes` → `ThicknessTransport::solve` → `update_ghosts(h)`.

**`src/bubble_config.hpp`** — `ConvectionScheme` enum moved here; output dir changed from `results_bubble` to `results`.

**`CMakeLists.txt`** / **`tests/CMakeLists.txt`** — all bearing-era source files and tests removed; Phase 3 sources added; 4 bubble tests registered.

**`CLAUDE.md`**, **`README.md`**, **`PLAN.md`**, **`PHYSICS.md`**, **`AUDIT.md`**, **`CHANGES.md`** — rewritten to describe the bubble solver only; all journal-bearing content removed.

### Test results (2 MPI ranks, 32×32 grid)

```
test_stereo_mesh:         11 passed, 0 failed
test_mask:                 5 passed, 0 failed
test_surface_operators:    6 passed, 0 failed
test_thickness_transport:  6 passed, 0 failed
```

---

## Phase B2 — Inverse Stereographic Mesh & VTK Output (2026-04-16)

### New files

**`src/stereo_mesh.hpp/cpp`** — `StereoMesh` class.
- Inverse stereographic projection $(u,v) \to (x,y,z)$ on hemisphere of radius $R$.
- 1D MPI domain decomposition in $u$.
- Precomputed per-cell geometry: `cell_area`, `face_len_e/n`, `d_u/v`, outward unit normals, surface gravity components $(g_u, g_v)$.
- `Vec3` helper struct for 3D vector arithmetic.

**`src/mask.hpp/cpp`** — `Mask` class.
- Circular domain mask: cells with $r_{2D} > R$ are blanked.
- Rim detection: active cells adjacent to at least one blanked or domain-edge neighbour.
- `blank_masked`, `enforce_rim_bc`, `enforce_rim_noslip` helpers.

**`src/bubble_io.hpp/cpp`** — `BubbleIO` namespace.
- Writes VTK unstructured grid (`.vtu`) per rank, `.pvtu` parallel descriptor (rank 0), and `.pvd` time collection (rank 0).
- Each active cell emitted as `VTK_QUAD` (type 9) with 4 projected 3D vertices (SW, SE, NE, NW order).

**`tests/test_stereo_mesh.cpp`** — 11 checks: north-pole projection, equator projection, area sum within 1% of $2\pi R^2$, unit normals, gravity at pole/equator, positive geometry, `east_face_len` helper.

**`tests/test_mask.cpp`** — 5 checks: active flags, center cell active, rim adjacency, `enforce_rim_bc`, `blank_masked`.

### Modified files

**`src/communicator.hpp/cpp`** — replaced with open-boundary (non-periodic) ghost exchange. Takes `StereoMesh`. Domain edges get zero-gradient ghost fill.

**`src/field.hpp/cpp`** — removed `Mesh` dependency. Primary constructor now takes raw dimensions `(name, n_u, n_v, n_ghost, GridLocation)`.

**`src/bubble_main.cpp`** — Phase 2 wiring: `StereoMesh`, `Mask`, geometry diagnostic fields, area verification, initial VTK output.

---

## Phase B1 — Infrastructure & Configuration (2026-04-16)

### New files

**`src/bubble_config.hpp/cpp`** — `BubbleConfig` struct.
- Covers geometry, film thickness, liquid/gas/evaporation properties, disjoining pressure, gravity, time-stepping, SIMPLE solver, physics toggles, and output parameters.
- `load_from_file` uses key=value text parser.

**`src/bubble_main.cpp`** — entry point skeleton with PETSc init, config load, empty time loop.

**`CMakeLists.txt`** — single `pancake_bubble` target; no bearing target on this branch.

**`tests/CMakeLists.txt`** — `add_bubble_test()` helper for registering test executables.
