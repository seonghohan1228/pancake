# Windows Port Checklist

This checklist is ordered for the first native Windows build. Do the
source-level cleanup in WSL/Linux first, then move to Windows.

## 1. Pre-Port Cleanup in WSL

1. Create a branch.

   ```bash
   git switch -c feature/windows-port-prep
   ```

2. Resolve the linear algebra dependency mismatch.

   Current source uses PETSc:

   - `src/linear_system.hpp:3`
   - `src/linear_system.cpp:16` through `src/linear_system.cpp:36`
   - `src/main.cpp:1`
   - `CMakeLists.txt:13` through `CMakeLists.txt:16`

   The strategic target says Trilinos. Before opening Visual Studio, decide
   which backend this branch is porting:

   - For the shortest path to a current-code Windows build, keep PETSc and make
     PETSc discovery portable.
   - For the strategic Trilinos path, migrate `src/linear_system.hpp` and
     `src/linear_system.cpp` to the selected Trilinos packages first, then
     update CMake.

   Do not leave the source half-PETSc and half-Trilinos.

3. Replace non-standard `M_PI`.

   Files:

   - `src/mesh.hpp:18`
   - `src/film_thickness.cpp:8`
   - `src/film_thickness.cpp:43`
   - `src/film_thickness.cpp:44`
   - `src/film_thickness.cpp:49`
   - `src/reynolds.cpp:25`
   - `src/reynolds.cpp:26`
   - `src/reynolds.cpp:33`
   - `tests/test_fvm.cpp:76`
   - `tests/test_inlets.cpp:68` through `tests/test_inlets.cpp:82`

   Recommended fix: include `<numbers>` and use
   `std::numbers::pi_v<double>`.

4. Add explicit standard-library includes.

   Files:

   - `src/field.cpp`: add `<algorithm>` for `std::fill` and `<stdexcept>` for
     `std::runtime_error`.
   - `src/linear_system.cpp`: add `<algorithm>` for `std::fill`.
   - `src/reynolds.cpp`: add `<cstdio>` for `std::fprintf`.
   - Tests using `std::max` or `std::min`: add `<algorithm>` where missing,
     especially `tests/test_bc.cpp` and `tests/test_velocity.cpp`.

5. Convert path handling to `std::filesystem::path`.

   Files and changes:

   - `src/config.hpp:67`: change `output_dir` from `std::string` to
     `std::filesystem::path`.
   - `src/config.hpp:70`: change `load_from_file(const std::string& path)` to
     `load_from_file(const std::filesystem::path& path)`.
   - `src/main.cpp:17` through `src/main.cpp:28`: store the config path as
     `std::filesystem::path`.
   - `src/io.cpp:51`, `src/io.cpp:58`, `src/io.cpp:69`, `src/io.cpp:79`,
     `src/io.cpp:101`, `src/io.cpp:161`, `src/io.cpp:196`, `src/io.cpp:281`,
     `src/io.cpp:327`, `src/io.cpp:330`, `src/io.cpp:338`, `src/io.cpp:348`,
     and `src/io.cpp:349`: replace raw string concatenation with path
     composition, for example `cfg.output_dir / "results.pvd"`.

6. Tighten global index types before large Windows tests.

   Files:

   - `src/mesh.hpp:10` through `src/mesh.hpp:13`
   - `src/linear_system.cpp:6`
   - `src/linear_system.cpp:14`
   - `src/linear_system.cpp:55`
   - `src/linear_system.cpp:80` through `src/linear_system.cpp:106`

   Use `PetscInt` for PETSc rows/sizes or the Trilinos equivalent for global
   indices. Keep local loop counters as `int` only where dimensions are known
   to fit.

7. Build and run the Linux checks before moving to Windows.

   ```bash
   cmake -S . -B build-linux -G Ninja -DCMAKE_BUILD_TYPE=Debug
   cmake --build build-linux -j
   ./build-linux/tests/test_eos
   mpirun -n 2 ./build-linux/tests/test_linear_system
   mpirun -n 2 ./build-linux/tests/test_fvm
   mpirun -n 2 ./build-linux/tests/test_elrod_a1
   mpirun -n 2 ./build-linux/tests/test_elrod_a2
   mpirun -n 2 ./build-linux/tests/test_inlets
   mpirun -n 2 ./build-linux/tests/test_velocity
   mpirun -n 2 ./build-linux/tests/test_bc
   ```

## 2. CMake Hygiene in WSL

1. Raise the root CMake minimum to 3.21 or newer.

   File: `CMakeLists.txt:1`

2. Split a reusable core target.

   Current files:

   - `CMakeLists.txt:21` through `CMakeLists.txt:25`
   - `tests/CMakeLists.txt:7` through `tests/CMakeLists.txt:15`

   Target shape:

   ```cmake
   add_library(pancake_core ...)
   target_include_directories(pancake_core PUBLIC src)
   add_executable(pancake src/main.cpp)
   target_link_libraries(pancake PRIVATE pancake_core)
   ```

3. Replace hardcoded dependency paths.

   Current files:

   - `CMakeLists.txt:13` through `CMakeLists.txt:16`
   - `tests/CMakeLists.txt:1` through `tests/CMakeLists.txt:4`

   The root CMake file should own dependency discovery and tests should link
   `pancake_core`.

4. Register tests with CTest.

   Add `enable_testing()` at the root and use `add_test()` for each executable.
   For MPI tests, use CMake's MPI launcher variables instead of hardcoding
   `mpirun`.

5. Guard compiler-specific warnings and options.

   Use target-scoped options:

   ```cmake
   if(MSVC)
     target_compile_options(pancake_core PRIVATE /W4)
   else()
     target_compile_options(pancake_core PRIVATE -Wall -Wextra)
   endif()
   ```

6. Rebuild on Linux after every CMake change.

## 3. Repository Setup for Dual Builds

1. Confirm `.gitignore` includes:

   ```gitignore
   build-linux/
   build-windows/
   ```

2. Confirm `.gitattributes` is committed with:

   ```gitattributes
   * text=auto
   *.sh   text eol=lf
   *.bat  text eol=crlf
   *.ps1  text eol=crlf
   *.png  binary
   *.pdf  binary
   *.vtu  binary
   *.pvd  text eol=lf
   ```

3. Commit the WSL cleanup before starting the native Windows build.

## 4. Windows Toolchain Installation

1. Install Visual Studio 2022 Community.

2. Select the "Desktop development with C++" workload.

3. Include these components:

   - MSVC v143 toolset
   - Windows 11 SDK
   - CMake tools for Windows
   - Ninja, optional but useful

4. Install Git for Windows.

5. Optional: install Windows Terminal.

## 5. vcpkg Setup

Open a Developer PowerShell or Developer Command Prompt:

```bat
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg integrate install
setx VCPKG_ROOT C:\vcpkg
```

Open a new terminal so `%VCPKG_ROOT%` is available.

## 6. Dependency Installation

Current checked-in source:

```bat
%VCPKG_ROOT%\vcpkg install mpi:x64-windows
```

This installs the vcpkg MPI meta-port, which resolves to MS-MPI on Windows.

PETSc blocker:

- The audit did not find an official vcpkg PETSc port.
- If keeping PETSc, build PETSc from source against MSVC and MS-MPI, then pass
  the installation to CMake through `PETSC_DIR`.
- Keep the PETSc build and `pancake` build on normal Windows paths, not
  `\\wsl$`.

Trilinos blocker:

- The current source uses zero Trilinos packages.
- If moving to Trilinos before Windows, first migrate the source and then list
  the exact packages found from includes and CMake links.
- At audit time no official vcpkg `trilinos` port was found. If still true,
  build Trilinos from source on Windows with the selected packages enabled.

Future GUI dependencies:

```bat
%VCPKG_ROOT%\vcpkg install qtbase[widgets,opengl]:x64-windows vtk[qt,utf8]:x64-windows
```

Add `vtk[mpi]` only if the GUI truly needs VTK's MPI feature.

## 7. First Windows Build

1. Clone to a Windows filesystem path.

   ```bat
   mkdir C:\src
   cd C:\src
   git clone https://github.com/seonghohan1228/pancake.git
   cd pancake
   ```

2. Configure after the linear algebra dependency blocker is resolved.

   PETSc current-backend example:

   ```bat
   cmake -S . -B build-windows -G "Visual Studio 17 2022" -A x64 ^
     -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake ^
     -DVCPKG_TARGET_TRIPLET=x64-windows ^
     -DPETSC_DIR=C:\src\petsc-install
   ```

   Trilinos future-backend example:

   ```bat
   cmake -S . -B build-windows -G "Visual Studio 17 2022" -A x64 ^
     -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake ^
     -DVCPKG_TARGET_TRIPLET=x64-windows ^
     -DTrilinos_DIR=C:\src\trilinos-install\lib\cmake\Trilinos
   ```

3. Build.

   ```bat
   cmake --build build-windows --config Debug
   ```

4. Run.

   ```bat
   build-windows\Debug\pancake.exe -c config.txt
   mpiexec -n 2 build-windows\Debug\pancake.exe -c config.txt
   ```

5. Visual Studio workflow:

   - Open the repo folder as a CMake folder.
   - Select an `x64-Debug` configuration.
   - Set the vcpkg toolchain file in CMake settings.
   - Configure, build, then run `pancake.exe`.

## 8. Validation

1. Run all tests on Windows.

   Until CTest is added, run direct commands from `build-windows\Debug` or the
   generated test output directory.

2. Run the same reference simulation on Linux and Windows.

   ```bash
   mpirun -n 2 ./build-linux/pancake -c config.txt
   ```

   ```bat
   mpiexec -n 2 build-windows\Debug\pancake.exe -c config.txt
   ```

3. Compare the output structure.

   - Same number of timesteps in `results/results.pvd`.
   - Same number of `.pvts` and per-rank `.vts` files.
   - Same mesh extents and field names.

4. Compare field values with numerical tolerance.

   Recommended first tolerance for the same backend and same MPI rank count:

   - Dimensionless fields such as film content: relative L2 <= `1e-8`, max
     absolute <= `1e-7`.
   - Pressure and velocity fields: relative L2 <= `1e-8`, max relative <=
     `1e-6`, with an absolute floor of `1e-8` for near-zero values.
   - Integrated load and torque fields: relative difference <= `1e-6`.

   If the Windows port also changes the linear algebra backend, first require
   all physical invariants and tests to pass, then accept relative L2 <= `1e-6`
   for field comparisons until solver tolerance parity is established.

## 9. Daily Workflow Going Forward

- Commit from either OS; Git history is shared.
- Build Linux in `build-linux/` and Windows in `build-windows/`.
- Re-test both platforms after any change to:
  - file I/O or paths
  - CMake dependency discovery
  - compiler flags
  - MPI communication
  - global indexing
  - linear solver setup
  - new third-party dependencies
- Keep Linux sanitizers and debugging available in WSL/Linux.
- Use MSVC regularly because it catches different portability bugs than GCC or
  Clang.

## 10. Known Issues / FAQ

### CMake cannot find PETSc on Windows

Expected until PETSc is built and exported for Windows or the code is migrated
to Trilinos. Do not work around this with a hardcoded local path in committed
CMake.

### `M_PI` is undefined on MSVC

Replace all uses with `std::numbers::pi_v<double>` during pre-port cleanup.

### Build is slow under `\\wsl$`

Clone and build under `C:\src\pancake` or another normal Windows path.

### Output differs slightly between Linux and Windows

Small numerical differences are expected across compilers and math libraries.
Use the tolerance policy in the validation section and investigate larger
differences by checking solver iteration counts, MPI rank count, and dependency
versions.

### GUI dependencies are not installed yet

That is expected. The current CLI port should come first. Install Qt/VTK only
when the GUI branch starts.
