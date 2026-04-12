# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**pancake** is a 2D numerical solver for thin-film flow (journal bearing simulation), written in C++20 with MPI parallelism. It solves the Reynolds equation on a cylindrical mesh (theta × z) and outputs VTK files for ParaView visualization.

## Build & Run

```bash
# Build (from repo root)
cd build && cmake .. && make

# Run single process
./build/pancake

# Run parallel
mpirun -n <N> ./build/pancake
```

Dependencies: C++20 compiler, MPI (found via CMake's `find_package(MPI)`).

Output goes to `results/` directory; open `results/results.pvd` in ParaView.

## Architecture

The solver follows a standard CFD pattern: mesh → fields → time loop → output.

- **SimulationConfig** (`config.hpp`): Plain struct holding all simulation parameters (geometry, grid, time, physics, output). Currently hardcoded defaults, no file parsing.
- **Mesh** (`mesh.hpp/cpp`): Structured cylindrical grid. 1D domain decomposition in the theta direction across MPI ranks. Each rank owns a contiguous slice of theta cells with the full z extent.
- **Field / Fields** (`field.hpp/cpp`): Cell-centered scalar fields with ghost layers (default 2). `Field` stores data in a flat vector with ghost padding; `operator()(i,j)` handles ghost-cell offset arithmetic. `Fields` is a named container (`std::map<string, unique_ptr<Field>>`).
- **Communicator** (`communicator.hpp/cpp`): MPI ghost-cell exchange. Periodic in theta (rank 0 ↔ rank N-1 wrap). Packs/unpacks non-contiguous ghost strips via buffers and `MPI_Sendrecv`.
- **IO** (`io.hpp/cpp`): Writes per-rank `.vtr` (VTK rectilinear grid) files and a single `.pvd` collection file (rank 0 only).
- **Utils** (`utils.hpp`): Rank-0-only colored logger and z-boundary condition helpers (`FIXED` / `ZERO_GRADIENT`).

### Key design details

- Ghost cell indexing: `field(i,j)` with `i` in `[-n_ghost, n_theta_phys + n_ghost)` and `j` in `[-n_ghost, n_z_phys + n_ghost)`. Physical cells start at index 0.
- Domain decomposition is theta-only (1D); z is not decomposed (`n_z_local == n_z_global`).
- `Field::old_data` stores the previous timestep for time derivatives (`store_old_time()`).
- The solver loop in `main.cpp` is still being developed — equation assembly and solve steps are marked with TODOs.

### Conventions

The codebase follows the following conventions

- Use snakecase for all variables, functions, and file names. Exceptions include Class, names that require capitalization.
- Use minimal comments. Only when others might be confused with the shorthand. However, do not skip comments for function, class, etc. Users should be able to know how to use it just by looking at the headers.
- Respect scientific notations and terminology.
- Do not unnecessarily shorten variable names unless the name becomes absurdly long.
- Keep functionality and solvers modular and human readable. Users should be able to swap functionalities easily (e.g., equations, boundary conditions, etc.)
- Make sure it is difficult to accidentally corrupt data fields.
- Complex models (cavitation model, solvers) must be easy to disable at the highest level.
- When installing packages and a significant implementation, always update README.md.
- After edits, write detailed documentation (append edits and reasons for updates) so past edits can be reviewed later.
- Write a detailed mechanism of the current algorithm and physical model using detailed explanation and equations in PHYSICS.md
- To keep track of development plan, refer to PLAN.md and update accordingly. However, do not delete tasks that are done.
