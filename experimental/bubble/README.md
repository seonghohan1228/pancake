# Experimental bubble solver

Imported from `origin/feature/bubble` (72736c6). This is a separate, optional
executable; the default `pancake` executable remains the journal-bearing solver.

Implemented: stereographic hemisphere geometry, circular mask, open MPI halo
exchange, surface operators and gravity-driven thickness transport. Momentum
SIMPLE coupling, energy, evaporation, Marangoni flow, disjoining pressure and
rupture detection remain unfinished. Their configuration fields are placeholders
and do not activate those models. Adaptive stepping is also unfinished; the
current driver uses `dt`. The archived physics/plan documents describe the intended
model and must not be read as a claim that it is complete or validated.

## Build and test

Use the root project's MPI/PETSc discovery and Windows runtime deployment:

```bash
cmake -S . -B build-linux -DCMAKE_BUILD_TYPE=Release -DPANCAKE_BUILD_BUBBLE=ON
cmake --build build-linux -j
ctest --test-dir build-linux --output-on-failure
./build-linux/experimental/bubble/pancake_bubble -c bubble_config.txt
```

On Windows, add `-DPANCAKE_BUILD_BUBBLE=ON` to
`cmake --preset windows-native-mingw`, then use the usual build preset. The
executable is `build-windows-mingw/experimental/bubble/pancake_bubble.exe`.
Set a short `end_t` explicitly for smoke checks. Use a new or empty `output_dir`;
the executable rejects occupied directories to preserve previous results.

Integration also replaces non-standard pi constants with C++20 `std::numbers`,
destroys PETSc objects before finalization, and rejects failed linear solves.
The original four unit executables are registered with CTest at 1, 2 and 4 MPI
ranks, plus a CLI test covering output and preservation of existing files.

## Configuration

Parameters are loaded from `bubble_config.txt` at runtime. Missing file → hardcoded defaults.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `R_bubble` | 0.01 m | Bubble radius |
| `L_box` | 0.012 m | Half-width of stereographic 2D domain |
| `n_u`, `n_v` | 64 | Grid cells in each direction |
| `h_initial` | 1×10⁻⁶ m | Initial film thickness |
| `h_rim` | 1×10⁻⁶ m | Rim Dirichlet BC on h |
| `h_min` | 1×10⁻⁸ m | Minimum allowed thickness (clamping floor) |
| `h_rupture` | 1×10⁻⁸ m | Rupture threshold |
| `rho_l` | 1000 kg/m³ | Liquid density |
| `mu` | 1×10⁻³ Pa·s | Dynamic viscosity |
| `sigma_0` | 0.025 N/m | Surface tension at T_ref |
| `gamma_T` | −1.7×10⁻⁴ N/(m·K) | dσ/dT |
| `end_t` | 10.0 s | Simulation end time |
| `dt` | 1×10⁻⁶ s | Initial time step |
| `enable_gravity` | true | Toggle gravitational drainage |
| `enable_evaporation` | true | Placeholder; evaporation unfinished |
| `enable_marangoni` | true | Placeholder; Marangoni unfinished |
| `enable_disjoining` | true | Placeholder; disjoining pressure unfinished |
| `output_dir` | results | Output directory |
| `write_interval` | 0.01 s | Time between VTK snapshots |

---

## Project Structure

```
src/
  bubble_config.hpp/cpp         — BubbleConfig (all parameters) + ConvectionScheme
  stereo_mesh.hpp/cpp           — Inverse stereographic grid + geometry engine
  mask.hpp/cpp                  — Circular domain mask + rim BC helpers
  communicator.hpp/cpp          — MPI ghost exchange (open, non-periodic)
  bubble_io.hpp/cpp             — VTK unstructured (vtu/pvtu/pvd) output
  field.hpp/cpp                 — Scalar fields with ghost layers
  bubble_linear_system.hpp/cpp  — PETSc 5-point sparse system
  surface_operators.hpp/cpp     — Surface FVM operators
  thickness_transport.hpp/cpp   — Film thickness transport solver
  utils.hpp                     — Rank-0 coloured logger
  bubble_main.cpp               — Entry point and time loop

tests/
  test_stereo_mesh.cpp          — Geometry: projection, area, normals, gravity
  test_mask.cpp                 — Mask: active flags, rim detection, BCs
  test_surface_operators.cpp    — FVM operators: ddt, laplacian, divergence, source
  test_thickness_transport.cpp  — Transport: velocity, fluxes, solve, clamping
```

---

## Implementation Status

| Phase | Description | Status |
|-------|-------------|--------|
| B1 | Infrastructure & configuration | Complete |
| B2 | Stereographic mesh & VTK I/O | Complete |
| B3 | Thickness transport + gravity drainage | Complete |
| B4 | 2D Navier-Stokes (SIMPLE) | Pending |
| B5 | Thermodynamics & evaporation | Pending |
| B6 | Disjoining pressure & rupture | Pending |
