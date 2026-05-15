# pancake

Multiphysics thin-film solver for a hemispherical soap bubble, written in C++20 with MPI and PETSc. The film thickness $h(\mathbf{x}, t)$ evolves on a stereographically projected hemisphere under gravitational drainage, Marangoni convection, evaporation, and disjoining pressure.

See [PHYSICS.md](PHYSICS.md) for governing equations and numerical methods.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| C++20 compiler (GCC ≥ 11 or Clang ≥ 14) | Language standard |
| CMake ≥ 3.20 | Build system |
| MPI (OpenMPI or MPICH) | Domain decomposition and ghost-cell exchange |
| PETSc ≥ 3.15 (real scalar build) | Sparse linear algebra (KSP/PC) |

**Installing PETSc on Debian/Ubuntu:**
```bash
sudo apt install libpetsc-real-dev
```
The default CMake configuration expects PETSc at `/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`. Override with:
```bash
cmake .. -DPETSC_DIR=/path/to/petsc
```

---

## Build & Run

```bash
git clone https://github.com/seonghohan1228/pancake.git
cd pancake
mkdir build && cd build
cmake ..
make -j$(nproc)

# Single process
./pancake_bubble

# Parallel
mpirun -n <N> ./pancake_bubble

# Custom config
./pancake_bubble --config bubble_config.txt
```

Results are written to `results/`. Open `results/results.pvd` in ParaView.

---

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
| `enable_evaporation` | true | Toggle Sultan evaporation model |
| `enable_marangoni` | true | Toggle Marangoni convection |
| `enable_disjoining` | true | Toggle van der Waals disjoining pressure |
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
