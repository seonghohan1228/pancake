# pancake

`pancake` is a 2D finite-volume solver for compressible thin-film
lubrication in a journal bearing. It solves the Reynolds equation on a
cylindrical `(theta, z)` mesh with MPI parallelism and writes VTK/PVD files
for ParaView.

See [PHYSICS.md](PHYSICS.md) for the governing equations, FVM discretisation,
and solver details.

## Platforms

The project is intended to remain cross-platform from one CMake source tree.
Linux and Windows are both first-class targets:

| Platform | Supported toolchains |
|----------|----------------------|
| Linux | GCC >= 11 or Clang >= 14 |
| Windows | MSVC 2022, v143 toolset |

Use separate build directories for each platform, normally `build-linux/` from
WSL/Linux and `build-windows/` from Windows. See [BUILDING.md](BUILDING.md) for
the current build recipes and dependency notes.

The checked-in source currently uses MPI and PETSc. The Windows-port assessment
in [ASSESSMENT.md](ASSESSMENT.md) records the dependency mismatch with the
strategic Trilinos target and lists the work needed before a native Windows
build is expected to succeed.

## Dependencies

| Package | Purpose |
|---------|---------|
| C++20 compiler | Language standard |
| CMake >= 3.21 | Build system |
| MPI, such as OpenMPI, MPICH, or MS-MPI | Domain decomposition and ghost-cell exchange |
| PETSc >= 3.15, real scalar build | Current sparse linear algebra backend |

On Debian/Ubuntu, PETSc can be installed with:

```bash
sudo apt install libpetsc-real-dev
```

The current CMake files default to the Debian PETSc path
`/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real`. Override it with:

```bash
cmake -S . -B build-linux -DPETSC_DIR=/path/to/petsc
```

## Build

```bash
git clone https://github.com/seonghohan1228/pancake.git
cd pancake
cmake -S . -B build-linux
cmake --build build-linux -j
```

## Run

```bash
# Single process
./build-linux/pancake

# Parallel
mpirun -n <N> ./build-linux/pancake
```

Results are written to `results/`. Open `results/results.pvd` in ParaView for
visualisation.

## Tests

Tests are built under the selected build directory. Until CTest registration is
added, run the executables directly:

```bash
./build-linux/tests/test_eos
mpirun -n 2 ./build-linux/tests/test_linear_system
mpirun -n 2 ./build-linux/tests/test_fvm
mpirun -n 2 ./build-linux/tests/test_elrod_a1
mpirun -n 2 ./build-linux/tests/test_elrod_a2
mpirun -n 2 ./build-linux/tests/test_inlets
mpirun -n 2 ./build-linux/tests/test_velocity
mpirun -n 2 ./build-linux/tests/test_bc
```

## Configuration

Simulation parameters are loaded from `config.txt` at runtime. If the file is
missing, hardcoded defaults in `src/config.hpp` are used.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cavitation_model` | `ELROD_ADAMS` | `GUMBEL`, `ELROD_ADAMS`, or `FULL_SOMMERFELD` |
| `R` | `0.01` m | Shaft radius |
| `c` | `0.0001` m | Radial clearance |
| `e` | `0.00008` m | Eccentricity |
| `L` | `0.02` m | Bearing length |
| `omega` | `100` rad/s | Shaft angular velocity |
| `mu` | `0.01` Pa s | Dynamic viscosity |
| `rho` | `900` kg/m3 | Reference density |
| `p_cav` | `0` Pa | Cavitation pressure |
| `bulk_modulus` | `1e5` Pa | Compressibility |
| `n_theta_global` | `120` | Circumferential cells |
| `n_z_global` | `40` | Axial cells |
| `end_t` | `1.0` s | Simulation end time |
| `dt` | `0.01` s | Time step |

Example:

```ini
R = 0.01
c = 0.001
omega = 100.0
cavitation_model = ELROD_ADAMS
```

## Project Structure

```text
src/
  config.hpp            SimulationConfig and text config parsing
  mesh.hpp/cpp          Structured cylindrical grid, MPI decomposition
  field.hpp/cpp         Scalar fields with ghost layers
  communicator.hpp/cpp  MPI ghost-cell exchange
  fvm.hpp/cpp           FVM operators
  linear_system.hpp/cpp Current PETSc sparse system wrapper
  film_thickness.hpp/cpp Film geometry and inlet indicator
  equation_of_state.hpp Barotropic EOS
  reynolds.hpp/cpp      Reynolds equation assembly and post-processing
  io.hpp/cpp            VTK/PVD output
  utils.hpp             Rank-0 logger and boundary helpers
  main.cpp              CLI entry point and time loop
tests/
  test_*.cpp            Focused solver and physics checks
```
