# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**pancake** is a multiphysics thin-film solver for a hemispherical soap bubble, written in C++20 with MPI parallelism and PETSc for sparse linear algebra. The film thickness $h(\mathbf{x}, t)$ evolves on a stereographically projected hemisphere under gravitational drainage, Marangoni convection, evaporation, and disjoining pressure.

## Build & Run

```bash
# From repo root
cd build && cmake .. && make -j$(nproc)

# Single process
./pancake_bubble

# Parallel
mpirun -n <N> ./pancake_bubble

# Custom config (defaults to bubble_config.txt)
./pancake_bubble --config bubble_config.txt
```

Dependencies: C++20 compiler, MPI, PETSc ≥ 3.15 (real scalar build).

Output goes to `results/`; open `results/results.pvd` in ParaView.

## Architecture

```
src/
  bubble_config.hpp/cpp    — BubbleConfig: all solver parameters + ConvectionScheme
  stereo_mesh.hpp/cpp      — Inverse stereographic mesh + precomputed 3D geometry
  mask.hpp/cpp             — Circular domain mask + rim BC helpers
  communicator.hpp/cpp     — MPI ghost exchange (open, non-periodic boundaries)
  bubble_io.hpp/cpp        — VTK unstructured (vtu/pvtu/pvd) output
  field.hpp/cpp            — Field/Fields: scalar fields with ghost layers
  bubble_linear_system.hpp/cpp  — PETSc 5-point sparse system for surface FVM
  surface_operators.hpp/cpp     — ddt, laplacian, divergence, add_source, gradient
  thickness_transport.hpp/cpp   — Gravity velocity, face fluxes, h transport solve
  utils.hpp                — Rank-0 coloured logger
  bubble_main.cpp          — Entry point and time loop
```

### Key design details

- **Inverse stereographic projection**: $(u, v) \to (x, y, z)$ on hemisphere of radius $R$. The $(u, v)$ plane is the working domain; all geometry is precomputed in 3D in `StereoMesh::compute_geometry()`.
- **1D MPI decomposition in u**: each rank owns a contiguous slab of $u$-cells with the full $v$-extent. Ghost exchange is non-periodic (zero-gradient at domain edges).
- **Circular mask**: cells with $r_{2D} > R$ are blanked. Rim cells carry Dirichlet BCs applied via `BubbleLinearSystem::apply_mask()` (row replacement).
- **Surface FVM**: all operators in `SurfaceFVM` use 3D chord lengths and areas from `StereoMesh`. The `active_at(mesh, gu, j)` helper checks mask membership analytically for any global cell index.
- **Ghost indexing**: `field(i, j)` with $i \in [-n\_ghost, n\_u\_local + n\_ghost)$ and $j \in [-n\_ghost, n\_v\_global + n\_ghost)$. Physical cells start at index 0.
- **PETSc KSP**: BiCGStab + block Jacobi. Global row = $(offset\_u + i) \times n\_v\_global + j$ (no periodic wrap).

## Conventions

- Snake case for all variables, functions, and file names; exceptions for class names.
- Minimal comments — only where logic is non-obvious. Headers must be self-documenting.
- Respect scientific notation and terminology.
- Keep solvers modular: users should be able to swap equations, BCs, and physics toggles easily.
- Physics toggles (`enable_gravity`, `enable_evaporation`, …) must be honoured at the highest level.
- After any significant implementation: update README.md, append to CHANGES.md, update PHYSICS.md with equations and mechanisms, update PLAN.md task checkboxes.
