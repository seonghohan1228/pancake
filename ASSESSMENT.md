# Cross-Platform Assessment

Audit date: 2026-05-26

Scope: source code, tests, CMake files, root Markdown docs, and repository
metadata in the checked-in tree. This pass made documentation and metadata
changes only; no solver code was changed.

## Blockers

### [Blocker] Stated Trilinos stack does not match the checked-in code

- Evidence: `src/linear_system.hpp:3` includes `<petsc.h>`, and
  `src/linear_system.hpp:40` through `src/linear_system.hpp:42` store PETSc
  `Mat`, `Vec`, and `KSP` handles directly.
- Evidence: `src/linear_system.cpp:16` through `src/linear_system.cpp:36`
  constructs PETSc objects and configures PETSc KSP/PC solvers.
- Evidence: `src/main.cpp:1` includes `<petsc.h>` and `src/main.cpp:14`
  initializes PETSc.
- Evidence: no `Trilinos`, `Tpetra`, `Epetra`, `Belos`, `Amesos2`, `Ifpack2`,
  `Kokkos`, `MueLu`, `ML`, or `AztecOO` include or CMake link entry was found.

Impact: the project prompt says Trilinos is fixed, but the repository currently
uses PETSc. A Windows port cannot be planned accurately until this dependency
target is aligned. If the current PETSc implementation is kept for the first
native Windows executable, PETSc must be built and found on Windows. If Trilinos
is the required backend, `src/linear_system.*` must be migrated before the
Windows dependency plan is executable.

Effort: 1-3 days to make the existing PETSc backend discoverable and buildable
on Windows; several days to weeks if this is also the Trilinos migration pass.

### [Blocker] PETSc discovery is hardcoded to a Debian/Ubuntu path

- Evidence: `CMakeLists.txt:13` through `CMakeLists.txt:16` hardcode
  `/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real` and search only below
  that root.
- Evidence: `tests/CMakeLists.txt:1` through `tests/CMakeLists.txt:4` repeats
  the same hardcoded PETSc path instead of using a root target or inherited
  package result.

Impact: native Windows CMake configuration will fail immediately. It also makes
Linux builds brittle across PETSc versions and package managers.

Fix: centralize dependency discovery in the root CMake file. For the current
PETSc backend, provide a portable imported target or package-config path. For
the intended Trilinos backend, use `find_package(Trilinos REQUIRED)` and link
the actual packages used by the code.

Effort: 2-4 hours after the dependency target is decided.

### [Blocker] Non-standard `M_PI` is used in production code and tests

- Evidence: `src/mesh.hpp:3` includes `<math.h>` and `src/mesh.hpp:18` uses
  `M_PI`.
- Evidence: `src/film_thickness.cpp:8`, `src/film_thickness.cpp:43`,
  `src/film_thickness.cpp:44`, and `src/film_thickness.cpp:49` use `M_PI`.
- Evidence: `src/reynolds.cpp:25`, `src/reynolds.cpp:26`, and
  `src/reynolds.cpp:33` use `M_PI`.
- Evidence: tests also use `M_PI` at `tests/test_fvm.cpp:76` and
  `tests/test_inlets.cpp:68` through `tests/test_inlets.cpp:82`.

Impact: MSVC does not expose `M_PI` from `<cmath>` by default. These files are
likely to fail the first MSVC compile.

Fix: include `<numbers>` and use `std::numbers::pi_v<double>`, or define a
project-local `constexpr double pi` in a shared numeric constants header.

Effort: 20-40 minutes.

### [Blocker] Several files rely on transitive standard-library includes

- Evidence: `src/field.cpp:29` uses `std::fill` but `src/field.cpp:1` only
  includes `field.hpp`; neither `src/field.cpp` nor `src/field.hpp` includes
  `<algorithm>`.
- Evidence: `src/field.cpp:32` throws `std::runtime_error`, but no included
  project header includes `<stdexcept>` explicitly.
- Evidence: `src/linear_system.cpp:47` through `src/linear_system.cpp:52` use
  `std::fill` without an explicit `<algorithm>` include in that translation
  unit.
- Evidence: `src/reynolds.cpp:289` uses `std::fprintf` without `<cstdio>`.
- Evidence: tests such as `tests/test_bc.cpp:58` and
  `tests/test_velocity.cpp:50` use standard library functions through
  transitive includes.

Impact: GCC may compile this because another header happens to include the
needed declaration. MSVC's standard library headers are not guaranteed to expose
the same transitive declarations, so this can become compile failures.

Fix: add explicit standard headers to each translation unit or header that uses
the corresponding symbol.

Effort: 30-60 minutes.

## Warnings

### [Warning] Filesystem paths are mostly string concatenations

- Evidence: `src/config.hpp:67` stores `output_dir` as `std::string` and
  `src/config.hpp:70` accepts config paths as `const std::string&`.
- Evidence: `src/io.cpp:51`, `src/io.cpp:58`, and `src/io.cpp:69` build
  `results.pvd` paths with `cfg.output_dir + "/results.pvd"`.
- Evidence: `src/io.cpp:79` builds `cfg.output_dir + "/flat/results.pvd"`.
- Evidence: `src/io.cpp:101`, `src/io.cpp:161`, `src/io.cpp:196`,
  `src/io.cpp:281`, `src/io.cpp:327`, `src/io.cpp:330`, `src/io.cpp:338`,
  `src/io.cpp:348`, and `src/io.cpp:349` build output paths manually.

Impact: Windows accepts `/` in many APIs, but raw string path handling makes
non-ASCII paths, drive roots, UNC paths, and future GUI file pickers harder to
support consistently.

Fix: make config/input/output path APIs use `std::filesystem::path`; build
paths with `/` path composition; pass paths to streams directly.

Effort: 1-2 hours.

### [Warning] Integer-size choices are adequate for small cases but not robust

- Evidence: mesh dimensions and offsets are `int` in `src/mesh.hpp:10` through
  `src/mesh.hpp:13`.
- Evidence: global PETSc sizes are computed as `int` in
  `src/linear_system.cpp:6` and `src/linear_system.cpp:14`.
- Evidence: `src/linear_system.cpp:55` returns global row numbers as `int`,
  then converts them to `PetscInt` at `src/linear_system.cpp:80`,
  `src/linear_system.cpp:88`, `src/linear_system.cpp:94`,
  `src/linear_system.cpp:100`, and `src/linear_system.cpp:106`.
- Evidence: the only explicit wide counter found is `long long` in
  `tests/test_elrod_a2.cpp:122`, reduced with `MPI_LONG_LONG` at
  `tests/test_elrod_a2.cpp:135` and `tests/test_elrod_a2.cpp:136`.

Impact: this is not a Windows x64 `long` bug in the current source, but large
grids can overflow 32-bit `int` before reaching PETSc. Windows and Linux both
need the same explicit size policy.

Fix: use `std::int32_t` for intentionally small counts, `std::size_t` for
container indexing, and `PetscInt` or `std::int64_t` for global matrix rows and
global cell counts. Prefer `std::int64_t` plus `MPI_INT64_T` for portable test
reductions.

Effort: 1-2 hours.

### [Warning] CMake does not define reusable core/library targets

- Evidence: `CMakeLists.txt:21` through `CMakeLists.txt:25` glob all
  `src/*.cpp` directly into the `pancake` executable.
- Evidence: `tests/CMakeLists.txt:7` through `tests/CMakeLists.txt:15` glob the
  same source files again for every test executable.

Impact: the planned Qt/VTK GUI should link a core solver library. The current
shape makes the future GUI and tests duplicate source selection and dependency
link logic.

Fix: create a `pancake_core` library target from all solver sources except
`main.cpp`. Link the CLI executable and tests to that target.

Effort: 2-4 hours.

### [Warning] Tests are executable-only, not registered with CTest

- Evidence: `tests/CMakeLists.txt:19` through `tests/CMakeLists.txt:26` create
  test executables but no `enable_testing()` or `add_test()` appears in the root
  or test CMake files.

Impact: Windows and future CI validation will be manual and easy to skip.

Fix: add `enable_testing()` at the root and register serial and MPI tests with
`add_test()`, using `MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, and per-test
process counts.

Effort: 1-2 hours.

### [Warning] CMake minimum is slightly below the stated modern baseline

- Evidence: `CMakeLists.txt:1` requires CMake 3.20.

Impact: the requested cross-platform baseline is CMake >= 3.21 for better
modern multi-config/toolchain behavior. This is unlikely to be the first
failure, but it should be aligned before Windows work starts.

Fix: raise the root minimum to 3.21 after confirming the Linux environment has
that version.

Effort: 10-20 minutes.

### [Warning] VTK/PVD writer declares little endian while writing ASCII

- Evidence: `src/io.cpp:105`, `src/io.cpp:165`, `src/io.cpp:201`,
  `src/io.cpp:285`, `src/io.cpp:332`, and `src/io.cpp:340` write
  `byte_order="LittleEndian"`.
- Evidence: scalar/vector arrays are written with `format="ascii"` at
  `src/io.cpp:111`, `src/io.cpp:127`, `src/io.cpp:149`, `src/io.cpp:207`,
  `src/io.cpp:231`, and `src/io.cpp:258`.

Impact: there is no current binary endianness/alignment issue because the VTK
payloads are ASCII. If binary output is added later, endianness and VTK data
array encoding must be handled deliberately.

Fix: no immediate change needed for ASCII output. Document binary-output rules
before changing VTK output format.

Effort: 15 minutes for documentation; 1-2 hours if binary output is added.

## Style Findings

### [Style] `file(GLOB ...)` hides source changes from CMake regeneration

- Evidence: `CMakeLists.txt:21`, `CMakeLists.txt:22`, and
  `tests/CMakeLists.txt:7` use `file(GLOB ...)` for source/header lists.

Impact: not Windows-specific, but explicit target source lists make future
library splitting and review cleaner.

Fix: replace globs with explicit target source lists, or use
`target_sources()` on `pancake_core`.

Effort: 30-60 minutes.

### [Style] CLI config path is still a raw string

- Evidence: `src/main.cpp:17` defaults to `"config.txt"`, `src/main.cpp:22`
  stores the argument in `std::string`, and `src/main.cpp:28` passes it to
  `load_from_file`.

Impact: works for simple paths, but future GUI/open-file flows should use
`std::filesystem::path` end-to-end.

Fix: change config APIs to `std::filesystem::path`.

Effort: included in the filesystem warning above.

## Explicit Negative Findings

- POSIX-only headers: none found in source or tests. The audit found MPI and
  PETSc headers only, not `<unistd.h>`, `<dirent.h>`, `<sys/stat.h>`,
  `<sys/types.h>`, `<pwd.h>`, or `<fcntl.h>`.
- Direct shell/process calls: none found. No `system()`, `popen()`, `fork()`,
  or `exec*()` uses were found.
- GCC-only headers/extensions: none found. No `<bits/stdc++.h>`,
  `__attribute__`, `typeof`, GCC statement expressions, or `#pragma GCC` uses
  were found.
- Header case mismatches: none found. All project includes match lowercase file
  names present under `src/`.
- Custom CMake helpers: the current tree contains `cmake/deploy_mingw_runtime.cmake`
  to copy MinGW/PETSc runtime DLLs beside native Windows executables.

## Dependency Audit

### Current Code Dependencies

| Dependency | Evidence | Current Windows note |
|------------|----------|----------------------|
| MPI | `CMakeLists.txt`; `<mpi.h>` in communicator, mesh, IO, Reynolds, and tests | Use Microsoft MPI plus the MSYS2 `mingw-w64-x86_64-msmpi` package. |
| PETSc | `CMakeLists.txt`; `src/linear_system.hpp`; `src/main.cpp` | Use the MSYS2 `mingw-w64-x86_64-petsc` package. CMake should select the MPI-enabled `petsc-dmo` pkg-config module, not the serial `petsc-sso` package. |

### Trilinos Package Usage

No Trilinos packages are actually used in the checked-in source. The audit found
zero includes or link entries for:

- Epetra
- Tpetra
- Belos
- Amesos2
- Ifpack2
- ML
- AztecOO
- Kokkos
- MueLu
- Teuchos

This means there is no current source-derived Trilinos package command to run.
The package list must be regenerated after the code is migrated to Trilinos.
The current native Windows workflow proceeds with the existing PETSc backend
through MSYS2 MinGW64.

### Other Libraries

- BLAS/LAPACK: no direct code or CMake dependency found. They are transitive
  through PETSc if PETSc is used.
- Qt: no current code or CMake dependency. Planned GUI dependency.
- VTK: no current library dependency. Current output is handwritten VTK XML/PVD
  text in `src/io.cpp`. Planned GUI/embedded visualization dependency.
- HDF5, Boost, Eigen, fmt, spdlog: no current direct dependency found.

## Executive Summary

- Blockers: 4
- Warnings: 6
- Style findings: 2
- POSIX-only source usage: 0 findings
- Direct shell/process usage: 0 findings
- GCC-only header/extension usage: 0 findings
- Actual Trilinos package usage: 0 packages found

Estimated effort to make the checked-in code Windows-buildable depends on the
linear algebra decision. If the first Windows build keeps the current PETSc
backend, expect roughly 2-4 days: half a day for source/CMake portability, plus
1-3 days for PETSc/MS-MPI build and CMake discovery. If the project must move to
the stated Trilinos backend before Windows work, add the Trilinos migration
effort first; that is likely several days to weeks depending on package choice
and solver parity requirements.
