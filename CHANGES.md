# Change Log

---

## 2026-04-11 — Reynolds solver, FVM fix, and documentation

### New files

**`src/reynolds.hpp` / `src/reynolds.cpp`**

Implements `Reynolds::solve`, which assembles and solves the steady-state compressible Reynolds equation at each time step:

$$\frac{1}{R^2}\frac{\partial}{\partial\theta}\!\left[\gamma\frac{\partial p}{\partial\theta}\right] + \frac{\partial}{\partial z}\!\left[\gamma\frac{\partial p}{\partial z}\right] = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta}$$

with $\gamma = \rho h^3/(12\mu)$.

Steps inside the solver:
1. Compute $\gamma$ and $\rho h$ over physical + theta ghost cells (avoids a separate MPI sync).
2. Assemble via `FVM::laplacian` (diffusion) + manual Dirichlet z-BC correction + `FVM::add_source` (negated Couette term, per FVM sign convention).
3. Solve with PETSc KSP.
4. Apply Gumbel cavitation clamp: $p \leftarrow \max(p, p_{cav})$.
5. Update density via EOS: $\rho = \rho_0\,\theta_{film}(p)$.

The transient $\partial(\rho h)/\partial t$ term is left out for the static-geometry case (zero for constant $h$ at steady state). It is a natural extension for Phase 4+ dynamic bearing work.

**`PHYSICS.md`**

Created. Documents the governing equations, EOS, FVM discretisation (including the sign convention), boundary conditions, solver algorithm, and notation.

**`CHANGES.md`**

This file. Created to satisfy the "append edits and reasons" requirement from `CLAUDE.md`.

---

### Modified files

**`src/fvm.cpp` — bug fix: incorrect `ew_scale`**

*Before:* `ew_scale = d_z / (R * R * d_theta)`  
*After:* `ew_scale = d_z / (R * d_theta)`

Root cause: the east/west FVM face coefficient for the θ-direction is derived from the surface divergence theorem:

$$\oint \gamma(\nabla p \cdot \hat{n})\,dl \;\ni\; \gamma_e \cdot \frac{\Delta z}{R\,\Delta\theta} \cdot (p_E - p_P)$$

The arc-length gradient contributes one factor of $1/R$ (not $1/R^2$). The previous code applied an additional spurious $1/R$ from mistakenly omitting the area-element factor $R$ when integrating $(1/R^2)\partial_\theta[\gamma\partial_\theta p]$ over the cell. This caused all Laplacian solutions to be off by $1/R$ (a factor of 100 for $R = 0.01\,\text{m}$), failing the 2D Poisson test.

**`tests/test_fvm.cpp` — bug fix: wrong source sign in `test_laplacian_2d`**

*Before:* `p_rhs(i,j) = -lambda * sin(...)  // passed as source_field`  
*After:* `p_rhs(i,j) = lambda * sin(...)`

The FVM stencil assembles $-\text{div}(\gamma\nabla p)\cdot V = \text{source}$, so for PDE $\text{div}(\gamma\nabla p) = f$ one must pass `source_field = -f`. With $f = -\lambda\,p_{exact}$, the correct source field is $+\lambda\,p_{exact}$. The original test passed $f$ directly, giving the wrong sign on the RHS and a solution of $-p_{exact}$ instead of $p_{exact}$.

Also added a `max_err` print before the assertion to aid future debugging.

**`src/reynolds.cpp` — Couette source sign**

The Couette term $S = \frac{\omega}{2}\frac{\partial(\rho h)}{\partial\theta}$ is the RHS of $\text{div}(\gamma\nabla p) = S$, so `source_field = -S`. The Couette contribution is stored with a leading minus sign accordingly.

**`src/main.cpp` — wire in Reynolds solver; fix post-finalize MPI abort**

- Added `#include "reynolds.hpp"` and `#include "linear_system.hpp"`.
- `LinearSystem sys(mesh)` created once before the time loop (PETSc setup cost paid once).
- `Reynolds::solve(fields, sys, mesh, cfg)` replaces the `// TODO` in the time loop.
- All PETSc objects (`LinearSystem`, `Fields`) now live inside a nested scope block that closes before `PetscFinalize()`. Previously, their destructors (which call `KSPDestroy`, `MatDestroy`, etc.) ran after `PetscFinalize()`, triggering an Open MPI post-finalize abort.

**`README.md`**

- Corrected dependency list: removed OpenMP (not used), added PETSc with install instructions.
- Added test table, configuration reference, and project-structure overview.
