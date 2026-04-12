# pancake

2D finite-volume solver for compressible thin-film lubrication in a journal bearing. Solves the Reynolds equation on a cylindrical $(\theta, z)$ mesh with MPI parallelism and outputs VTK files for ParaView.

See [PHYSICS.md](PHYSICS.md) for the governing equations, FVM discretisation, and solver details.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| C++20 compiler (GCC ≥ 11 or Clang ≥ 14) | Language standard |
| CMake ≥ 3.20 | Build system |
| MPI (e.g., OpenMPI, MPICH) | Domain decomposition and ghost-cell exchange |
| PETSc ≥ 3.15 (real scalar build) | Sparse linear algebra (KSP/PC solvers) |

**Installing PETSc on Debian/Ubuntu:**
```bash
sudo apt install libpetsc-real-dev
```
The default CMake configuration expects PETSc at `/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`. Override with:
```bash
cmake .. -DPETSC_DIR=/path/to/petsc
```

---

## Build

```bash
git clone https://github.com/yourusername/pancake.git
cd pancake
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## Run

```bash
# Single process
./build/pancake

# Parallel (recommended)
mpirun -n <N> ./build/pancake
```

Results are written to `results/`. Open `results/results.pvd` in ParaView for visualisation.

---

## Tests

```bash
# Unit tests (run from build/)
./tests/test_eos
mpirun -n 2 ./tests/test_linear_system
mpirun -n 2 ./tests/test_fvm
```

| Test | Validates |
|------|-----------|
| `test_eos` | Barotropic EOS round-trip, cavitation branch, density |
| `test_linear_system` | PETSc diagonal solve, 1D Helmholtz (periodic θ, 2 ranks) |
| `test_fvm` | `ddt` coefficients, 2D Poisson (2nd-order accuracy), upwind divergence |

---

## Configuration

Simulation parameters are loaded from `config.txt` at runtime. If the file is missing, hardcoded defaults in `src/config.hpp` are used.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cavitation_model` | ELROD_ADAMS | GUMBEL or ELROD_ADAMS |
| `R` | 0.01 m | Shaft radius |
| `c` | 0.001 m | Radial clearance |
| `e` | 0.0008 m | Eccentricity |
| `L` | 0.02 m | Bearing length |
| `omega` | 100 rad/s | Shaft angular velocity |
| `mu` | 0.01 Pa·s | Dynamic viscosity |
| `rho` | 900 kg/m³ | Reference density |
| `p_cav` | 0 Pa | Cavitation pressure |
| `bulk_modulus` | 1×10⁵ Pa | Compressibility |
| `n_theta_global` | 120 | Circumferential cells |
| `n_z_global` | 40 | Axial cells |
| `end_t` | 1.0 s | Simulation end time |
| `dt` | 0.01 s | Time step |

Example `config.txt`:
```ini
R = 0.01
c = 0.001
omega = 100.0
cavitation_model = ELROD_ADAMS
```

---

## Project Structure

```
src/
  config.hpp           — SimulationConfig (all parameters)
  mesh.hpp/cpp         — Structured cylindrical grid, MPI decomposition
  field.hpp/cpp        — Cell-centred scalar fields with ghost layers
  communicator.hpp/cpp — MPI ghost-cell exchange (periodic θ)
  fvm.hpp/cpp          — FVM operators: ddt, laplacian, divergence, add_source
  linear_system.hpp/cpp— PETSc sparse system assembly and KSP solve
  film_thickness.hpp/cpp— Static eccentric film geometry
  equation_of_state.hpp— Barotropic EOS (inline)
  reynolds.hpp/cpp     — Reynolds equation assembler and solver
  io.hpp/cpp           — VTK (.vtr / .pvd) output
  utils.hpp            — Rank-0 logger, z-boundary condition helpers
  main.cpp             — Entry point, time loop
tests/
  test_eos.cpp
  test_linear_system.cpp
  test_fvm.cpp
```
