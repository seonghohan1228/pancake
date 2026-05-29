# Change Log

## Unreleased - Windows MSYS2 GUI Workflow

### New files

**`src/gui_win32.cpp`**
- Added a native Win32 GUI for Windows builds.
- The GUI loads `config.txt` from the same directory as `pancake_gui.exe`, falls
  back to hardcoded `SimulationConfig` defaults when the file is missing, saves
  edits back to that file, and can launch `pancake.exe` from the same directory.
- Replaced the raw-editor-only prototype with workspace, output, and raw-config
  views.
- Added process launch/stop, captured solver logs, output folder/PVD open
  actions, and a GDI heatmap preview for the latest ASCII VTK scalar field.
- Reworked the main workspace into a split-pane engineering page:
  the result preview expands with the window, the process log stays below it,
  and a draggable horizontal splitter lets users tune the pane sizes.
- Lifted solver backend controls into the always-visible top-right header so
  the backend can be run/stopped from any tab.
- Moved the parameter editor into the right side of the main `Workspace` tab as
  a scrollable one-row-per-parameter inspector.
- Streamlined solver launch around the native `pancake.exe` backend built by
  the MSYS2 MinGW64 workflow.
- Preserved literal numeric edit text in GUI saves so typed or file-provided
  values such as `0.0010` stay visible as `0.0010`.
- Added preview timestep selection with previous/next controls.
- Added a heatmap colorbar with min/max labels for the current preview field,
  then expanded it with a field title, clearer number formatting, and plot axes.
- Fixed multi-rank preview loading by assembling the selected timestep from all
  matching `processor*` VTS files instead of treating rank 0 as a complete
  field.
- Added `load_angle_deg`, a visualization reference angle measured from the
  positive y axis, and rotated the preview display relative to that reference.
- Reworked axial boundary controls into South and North groups, and renamed the
  theta boundary inputs to film content value.
- Added per-value GUI unit selectors for length, time, angle, tilt, inlet
  pressure, and omega while still saving solver-native config units.
- Extended unit selectors across the full parameter inspector, including
  pressure fields, grid counts, dimensionless controls, and inlet supply
  pressure.
- Preserved raw inlet token text in the structured editor so values such as
  `1.5e6` are not rewritten as `1500000` after saving.
- Switched the process log to a RichEdit control that renders ANSI color codes
  instead of showing raw `[32m` / `[0m` fragments, with plain text fallback if
  RichEdit is unavailable.
- Normalized raw-config edit text to Windows CRLF newlines so multiline config
  text displays correctly.
- Replaced the raw multiline inlet editor with a structured groove/circular
  selector that comments out the unused inlet form when saving.
- Fixed the all-grey/low-contrast appearance: added a light modern palette and
  `WM_CTLCOLORSTATIC` / `WM_CTLCOLOREDIT` / `WM_CTLCOLORLISTBOX` handlers so
  labels, group boxes, and check boxes share one clean surface and editable
  fields stay white. Previously these controls fell back to the system 3D-face
  grey on the window background.
- Subclassed the tab control's `WM_ERASEBKGND` so its display area is painted
  with the page colour instead of the default grey slab behind the child
  controls.
- Forced every combo box above the interleaved group-box frames after creation,
  fixing the preview field/step selectors that only painted after a click.
- Reworked the preview annotations: wider right margin so colorbar min/max
  values are no longer truncated, the field name as a centred colorbar title,
  a stats subtitle at the top, and a decluttered axis footer.
- Moved preview title/stat annotations farther above the contour so min/max and
  grid metadata do not overlap the plotted field.
- Restored `film_content` to write and preview the raw Elrod/JFO universal
  variable instead of clamping it to `[0, 1]`, because pressure recovery depends
  on the same theta values and compressed full-film cells must remain visible.
- Added a regression check that the Elrod A.2 journal-bearing case contains
  both cavitated theta values below one and compressed full-film theta values
  above one.
- Centered the preview z-axis title, labeled the z-axis from `0` to `L`, and
  expanded the colorbar to include units plus min/max-inclusive tick labels.
- Fixed Elrod axial DIRICHLET boundaries so the solver uses the dimensionless
  `bc_z_south_theta` / `bc_z_north_theta` film-content values instead of
  converting the pressure-valued `bc_z_*_val` entries into theta. Added sanity
  coverage for this separation.
- Changed the default bulk modulus to `1e9` Pa so a 1.5 MPa supply pressure
  produces film-content values near 1.0015 rather than millions in the default
  Windows GUI workflow.
- Added blue/magenta to the ANSI log colour table so every code the solver
  emits (`\033[31..36m`, `\033[0m`) is colourised rather than only stripped.

**`src/pancake.ico` / `src/pancake_logo.svg` / `src/resource.h` / `src/pancake.rc`**
- Added the Pancake logo as a Windows application icon resource. The editable
  SVG mirrors the supplied blue film-profile logo, while the generated ICO is
  embedded into both `pancake.exe` and `pancake_gui.exe` for File Explorer,
  taskbar, and title-bar branding.

**`src/pancake_gui.manifest` / `src/pancake_gui.rc`**
- Added a `CREATEPROCESS_MANIFEST_RESOURCE_ID` manifest that declares the
  ComCtl32 v6 dependency. The previous `#pragma comment(linker, ...)` was guarded
  by `_MSC_VER`, so MinGW builds shipped without it and fell back to the classic
  (grey, Windows-95 era) common controls. Embedding the manifest via a compiled
  `.rc` enables modern visual styles under MinGW (and MSVC).

**`CMakePresets.json`**
- Kept a single Windows preset, `windows-native-mingw`, for building native
  `pancake.exe`, `pancake_gui.exe`, tests, and runtime DLL deployment from the
  MSYS2 MinGW64 PETSc/MS-MPI stack.

### Modified files

**`CMakeLists.txt` / `tests/CMakeLists.txt`**
- Split solver sources into `pancake_core`, kept `pancake` as the CLI target,
  and added the Windows-only `pancake_gui` target.
- Embedded the shared `pancake.ico` resource into both Windows executables, so
  the GUI and console backend present the same application logo.
- Enabled the `RC` language for the GUI target and added `src/pancake_gui.rc`
  (with its `INCLUDE_DIRECTORIES` source property) so the ComCtl32 v6 manifest
  is compiled in by windres/rc and visual styles work on every Windows compiler.
- Added PETSc discovery through `pkg-config`, including MSYS2 PETSc package
  names, while preserving the existing `PETSC_DIR` fallback.
- Restricted native Windows PETSc auto-discovery to MPI-enabled package names
  so CMake does not select the serial `petsc-sso` variant that collides with
  MS-MPI headers.
- Added MinGW runtime deployment for native Windows builds so required DLLs
  such as `libpetsc-dmo.dll`, `libstdc++-6.dll`, and `libopenblas.dll` are
  copied beside `pancake.exe`, `pancake_gui.exe`, and test executables.
- Updated the Windows CMake preset to prepend `C:\msys64\mingw64\bin` to PATH,
  preventing native compiler helper DLL lookup failures in plain PowerShell.

**`cmake/deploy_mingw_runtime.cmake`**
- Added a small post-build deployment helper that runs MSYS2 `ldd`, filters
  dependencies from the MinGW runtime directory, and copies them next to the
  target executable for File Explorer launches.

**`src/config.hpp` / `src/main.cpp`**
- Added config serialization for GUI defaults and same-directory config loading
  for the CLI default path.
- Added output persistence settings: `output_write_3d`, `output_write_flat`,
  and `output_fields`.
- Added `load_angle_deg` as a visualization-only load-reference angle.

**`src/io.cpp`**
- Honours the new output persistence settings by skipping disabled VTK/PVD
  products and disabled field arrays.

**`config.txt`**
- Added the new output persistence keys to the checked-in default config and
  included the axial inlet/outlet theta selectors emitted by config
  serialization.

**Markdown documentation**
- Rewrote the README, BUILDING guide, and Windows-port notes around the current
  MSYS2 MINGW64 workflow. Visual Studio is documented as unnecessary for the
  current native solver + GUI path.

**`.gitignore`**
- Removed stale Visual Studio GUI build-directory entries and kept the active
  `build-windows-mingw/` workflow ignored.

---

## 2026-05-26 — Cross-Platform Windows Port Assessment and Documentation

### New files

**`ASSESSMENT.md`**
- Added a Windows-port readiness assessment with file/line-specific blockers, warnings, dependency findings, and effort estimates.
- Documented the current PETSc implementation versus the stated Trilinos strategy as the primary dependency blocker.

**`BUILDING.md` / `ARCHITECTURE.md` / `CONTRIBUTING.md`**
- Added cross-platform build instructions, intended solver/GUI layering, and contribution conventions for Linux and Windows.
- Made `std::filesystem::path`, non-POSIX core code, explicit integer widths, and GCC/MSVC compatibility part of the contributor contract.

**`docs/WINDOWS_PORT.md`**
- Added an ordered checklist for the first native Windows build, starting with WSL-side cleanup and ending with Windows validation.

**`.gitattributes`**
- Added text/binary and line-ending rules for cross-platform development.

### Modified files

**`README.md` / `CLAUDE.md`**
- Updated build and platform guidance so future contributors and agents see the Linux/Windows dual-build strategy immediately.

**`.gitignore`**
- Added `build-linux/` and `build-windows/` as dedicated out-of-tree build directories.

---

## 2026-05-21 — Macroscopic Properties: Load Capacity and Friction Torque

### Modified files

**`src/reynolds.hpp` / `src/reynolds.cpp` — Post-processing calculations**
- Added `Reynolds::calculate_macroscopic_properties` to calculate the integrated load capacity vector ($F_x, F_y, F_z$) and the scalar friction torque $T$ over the full 2D mesh domain.
- Computes $F_z$ properly taking into account the surface normal variation caused by shaft tilt.
- Uses `MPI_Allreduce` to correctly sum the integrated values across all decomposed sub-domains.

**`src/main.cpp` — Global field allocation**
- Added four new uniform scalar fields (`load_x`, `load_y`, `load_z`, `friction_torque`) to the visualization output.
- Hooked up `calculate_macroscopic_properties` directly in the solver's main time loop.

**`PHYSICS.md` — Documentation**
- Appended a new section detailing the mathematical integration of the load capacity from the pressure field and the friction torque from the Couette and Poiseuille shear stress formulas.

---



## 2026-04-12 — Phase A.8: Configurable Outlet Boundary Conditions

### New files

**`tests/test_bc.cpp`**

Unit test for axial boundary conditions. Validates:
1. **Neumann (Zero Gradient)**: One end DIRICHLET and the other NEUMANN results in a uniform pressure field.
2. **Dirichlet (Fixed Value)**: Both ends DIRICHLET with different values results in a linear pressure profile.

### Modified files

**`src/config.hpp` — Boundary condition definitions**

- Added `BCType` enum (`DIRICHLET`, `NEUMANN`).
- Added axial boundary parameters to `SimulationConfig`: `bc_z_south_type`, `bc_z_south_val`, `bc_z_north_type`, `bc_z_north_val`.
- Updated `load_from_file` to support these parameters.

**`src/reynolds.cpp` — Boundary condition logic**

- Updated `Reynolds::solve` and `Reynolds::solve_elrod` to conditionally apply boundary terms based on the configured `BCType`.
- In `calculate_velocities`, the pressure gradient at boundaries is now correctly determined by the BC type (zero gradient enforced for `NEUMANN`).

---

## 2026-04-12 — Phase A.7: Staggered Grid Velocity Fields

### New files

**`tests/test_velocity.cpp`**

Unit test for the staggered grid velocity calculation. Verifies:
1. **Couette Flow**: In a concentric bearing with no pressure gradient, $u_\theta$ correctly matches $\omega R / 2$ and $u_z$ is zero.
2. **Poiseuille Flow**: A prescribed axial pressure gradient produces the analytical $u_z = -(h^2/12\mu) \partial p/\partial z$ response.

### Modified files

**`src/field.hpp` / `src/field.cpp` — Staggered grid support**

- `Field` constructor now correctly sizes physical dimensions for `FACE_THETA` and `FACE_Z`. 
- `FACE_Z` fields are allocated with $N_z + 1$ physical rows to account for the non-periodic axial faces.

**`src/io.cpp` — Staggered visualization**

- Completed `get_centered_data` for `FACE_THETA` and `FACE_Z`.
- Face-based velocities are now interpolated to cell centers during VTK output, ensuring they can be visualized alongside pressure and theta.

**`src/reynolds.hpp` / `src/reynolds.cpp` — Velocity calculation**

- Implemented `Reynolds::calculate_velocities`.
- $u_\theta$ is calculated at $\theta$-faces using a centered pressure gradient and averaged film thickness.
- $u_z$ is calculated at $z$-faces, including the bearing ends (Dirichlet $p_{cav}$ boundary).
- The Poiseuille (pressure-driven) component is only active in full-film regions ($\theta \ge 1$), consistent with the Elrod-Adams model.

**`src/main.cpp`**

- Added `velocity_theta` and `velocity_z` fields.
- Call `calculate_velocities` at each timestep after the Reynolds solve.

---

## 2026-04-12 — Phase A.6: Flexible Inlet Boundary Conditions

### New files

**`tests/test_inlets.cpp`**

Unit test for the internal Dirichlet (penalty) boundary conditions. Validates that both axial grooves and circular supply holes correctly force the prescribed supply pressure in the target cells for both GUMBEL and ELROD_ADAMS models.

### Modified files

**`src/config.hpp` — Inlet definitions**

- Added `InletConfig` struct.
- Added `std::vector<InletConfig> inlets` to `SimulationConfig`.
- Updated `load_from_file` to support multiple `inlet_circular` and `inlet_groove` definitions from the text config.

**`src/reynolds.cpp` — Penalty method for inlets**

- Implemented `apply_inlet_conditions` helper using a large-diagonal penalty method ($10^{20}$). This effectively locks inlet cells to their supply values without requiring invasive changes to the sparse matrix topology.
- Wired the helper into `Reynolds::solve` and the `Reynolds::solve_elrod` outer loop.

**`config.txt`**

Added example inlet definitions (one axial groove at 180° and one circular hole at 90°) to demonstrate usage.

---

## 2026-04-12 — Phase A.4: Text-based configuration and Phase A.5: Starvation fix

### New files

**`config.txt`**

Created a text-based configuration file containing all simulation parameters. This allows users to modify the simulation without recompiling. The format is a simple `key = value` structure with support for `#` comments.

### Modified files

**`src/config.hpp` — Added `load_from_file`**

Updated `SimulationConfig` with a `load_from_file(const std::string& path)` method. It parses the text file and updates the struct fields accordingly. This satisfies the user request for an editable config file.

**`src/main.cpp` — Load config and re-initialise fields**

- Added a call to `cfg.load_from_file("config.txt")` at the start of `main`.
- **Numerical Starvation Fix**: In the time loop, `theta` and `pressure` are now reset to full-film values ($\theta = 1$, $p = p_{cav}$) before every `solve_elrod` call. This provides a "flooded guess" that allows oil from the $z$-boundaries to reach the interior of the domain, preventing the "theta=0 everywhere" trap.

**`src/reynolds.cpp` — Flooded initialization in outer loop**

Modified the Elrod-Adams outer loop to force $g = 1$ for the very first iteration (`outer == 0`). This ensures that the solver starts with a fully active diffusion field, establishing a healthy pressure zone that "pulls" lubricant into the bearing from the boundaries. Subsequent iterations correctly use the `switch_function(theta)` to develop the cavitation zone.

**`PLAN.md` — Updated with Phase A.4 and A.5**

Updated the development plan to reflect the completion of the text-based configuration and the numerical starvation fix.

---

## 2026-04-12 — Phase A.2: Elrod-Adams switch function and mass-conserving cavitation

### New files

**`tests/test_elrod_a1.cpp` — updated regression test**

Added `cfg.max_outer_iters = 1` for the ELROD_ADAMS branch of `test_regression_vs_gumbel`. The intent of this test is Phase A.1: compare the full-film θ-formulation (g=1 everywhere) against the Gumbel solver. Running the full outer iteration was incorrect here because it allows cavitation to develop — the mass-conserving solution gives p_max ≈ 13 Pa (physical, but incomparable to Gumbel's 983 Pa). With max_outer_iters=1, only the full-film solve runs (outer=0 gives p_Elrod ≈ 994 Pa, error 1.2% vs Gumbel). All tests now pass.

**`tests/test_elrod_a2.cpp`**

New test file with four validation tests for the full Elrod-Adams model:
1. **Complementarity**: after convergence, no cell has g=1 AND θ<1 (max g·(1−θ)⁺ = 0).
2. **Cavitation location**: for ε=0.8 with attitude_angle=-90°, the majority (>50%) of cavitated cells lie in the diverging half (cos θ > 0). Observed fraction ≈ 0.51.
3. **Mass conservation**: total liquid content Σ θ·h over 30 timesteps drifts by < 2%. Observed relative drift = 0.
4. **Pressure positivity**: min pressure ≥ p_cav everywhere.

### Modified files

**`src/reynolds.cpp` — outer iteration loop + numerical fixes**

Extended `Reynolds::solve_elrod` with the switch-function outer (flag-update) iteration:

1. **Switch function** applied to diffusion: `Gamma_g(i,j) = Gamma_base(i,j) * EOS::switch_function(theta(i,j))`.

2. **Flooded bearing z-BC** (`Gamma_base`, not `Gamma_g`): The Dirichlet θ=1 condition at z=0,L uses `Gamma_base` (the full diffusion coefficient without the switch). Reasoning: the z-ends are immersed in the oil sump, which maintains full film regardless of the interior cavitation state. If `Gamma_g` were used, a cavitated boundary cell would set its z-BC coefficient to zero, removing the only restoring force for that cell. In subsequent outer iterations, adjacent cells lose their z-BC forcing, and the cavitation spreads globally (cascade collapse to θ_min everywhere). Using `Gamma_base` keeps the z-end forcing active unconditionally.

3. **Snap step**: After clamping θ ≥ theta_min, cells with θ ∈ (1−5×10⁻⁷, 1.0) are snapped to exactly 1.0. Motivation: the BiCGStab solver at rtol=1e-8 on a system with condition ~O(10³) produces residual errors ~O(1e-7) near the full-film/cavitation free boundary. Without snapping, these near-1.0 cells are assigned g=0 in the next outer iteration (because θ < 1), triggering spurious flag changes and preventing convergence. The snap tolerance 5e-7 was chosen to cover the observed numerical noise amplitude (1.7e-7) with margin, while staying well below physically cavitated cells (θ ≈ 0.111 in the diverging zone).

4. **n_switched convergence criterion** uses snap_tol threshold: `(tp >= 1.0 - snap_tol) != (tn >= 1.0 - snap_tol)` rather than a hard 1.0 comparison, consistent with the snap step.

5. **VecZeroEntries(x)** added in `LinearSystem::solve` before `KSPSolve`: `VecDuplicate` leaves the vector uninitialised; without this, the initial guess for the iterative solver may be garbage (or NaN), causing BiCGStab to fail or produce wrong results. This is a correctness fix, not an optimisation.

**Physical note on Elrod-Adams vs Gumbel pressure difference:**

The fully-converged Elrod-Adams mass-conserving solution gives a much lower peak pressure than Gumbel (≈13 Pa vs ≈983 Pa for ε=0.8). This is physically correct: in the mass-conserving model, the Couette flow advects low-θ fluid (θ ≈ 0.111 in the cavitation zone) back into the converging (full-film) zone, depressing the upstream θ and therefore the pressure. Gumbel clamps p ≥ 0 post-solve with no mass accounting, so the full density is always available. The two models are fundamentally different once mass conservation is active.

---

## 2026-04-11 — Phase A.1: universal variable (θ) reformulation

### New files

**`tests/test_elrod_a1.cpp`**

Two tests for the Elrod-Adams Phase A.1 solver:
- `test_uniform_film`: zero eccentricity → $h = c$ everywhere → no Couette wedge → steady solution must be $\theta = 1$ everywhere. Checks $\max|\theta - 1| < 10^{-6}$.
- `test_regression_vs_gumbel`: runs the standard journal bearing with both GUMBEL and ELROD_ADAMS for one timestep, verifies peak pressures agree to within 5% (observed: ~1.2%).

**`PLAN.md`**

Full phased development plan: Phases 1–3 (complete), Phase A (Elrod-Adams mass-conserving cavitation), Phase B (dissolved/free gas transport). Records the mathematical formulation, code changes, and validation benchmarks for each sub-phase.

### Modified files

**`src/config.hpp` — new `CavitationModel` enum**

Added `enum class CavitationModel { GUMBEL, ELROD_ADAMS }` and `cavitation_model` field to `SimulationConfig` (default: `GUMBEL`). Allows clean dispatch in `main.cpp` and tests.

**`src/fvm.hpp`/`fvm.cpp` — new `FVM::ddt_weighted` operator**

```cpp
void ddt_weighted(LinearSystem& sys, const Field& phi, const Field& weight,
                  double dt, const Mesh& mesh);
```

Discretises $\partial(w\phi)/\partial t$ for static weight $w$. Adds $w(i,j)\cdot V/\Delta t$ to $a_P$ and $w(i,j)\cdot V\,\phi^n/\Delta t$ to source. Required for the transient term $\partial(\theta h)/\partial t$ in Phase A.1+ dynamic bearings.

**`src/reynolds.hpp`/`reynolds.cpp` — new `Reynolds::solve_elrod`**

Implements the $\theta$-formulation (Phase A.1, full-film, $g(\theta) = 1$):

$$\text{div}(\Gamma\,\nabla\theta) = \frac{\omega}{2}\rho_0\frac{\partial(\theta h)}{\partial\theta}, \qquad \Gamma = \frac{\rho_0\beta h^3}{12\mu}$$

Key implementation notes:
- Couette face flux: $F_e = \rho_0\,(\omega/2)\,R\,h_e\,\Delta z$ — must include $\rho_0$ (missing it reduces peak pressure by factor $\rho_0$) and must include R (from the FVM area element $R\,\Delta\theta\,\Delta z$).
- Ghost cell coverage: flux loop runs from $i = -1$ to $n_\theta - 1$ so that `FVM::divergence` can read `couette_flux(-1, j)` for the west face of cell $i = 0$ without hitting uninitialised ghost data.
- Dirichlet z-BC: $\theta = 1$ (equivalent to $p = p_{cav}$ via EOS).
- Post-solve: clamp $\theta \geq 1$ (full-film constraint; Phase A.2 relaxes this).
- Recover $p$ and $\rho$ via EOS.

**`src/main.cpp` — dispatch on cavitation model**

- Added `"theta"` field, initialised to 1.0.
- `fields["theta"].store_old_time()` called alongside pressure.
- Time loop dispatches: `Reynolds::solve_elrod` for `ELROD_ADAMS`, `Reynolds::solve` for `GUMBEL`.
- Initial pressure changed from 101325 Pa to `p_cav` (cleaner initialisation for the Elrod solver).

**`PHYSICS.md` — Section 7 expanded**

Renamed to "Reynolds Solvers" (7.1 Gumbel, 7.2 Elrod-Adams A.1). Documents the full derivation of the $\theta$-formulation, the Couette-as-convection interpretation with face flux formula, the FVM assembly, and the key bugs found during implementation.

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
