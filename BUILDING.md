# Building

`pancake` is intended to build from one CMake source tree on Linux and Windows.
Keep platform builds in separate directories:

```text
build-linux/
build-windows/
```

The current checked-in source uses MPI and PETSc. The strategic Windows-port
target in [ASSESSMENT.md](ASSESSMENT.md) calls out that the repository does not
yet contain Trilinos code even though Trilinos is the intended long-term solver
stack.

## Linux

### Prerequisites

- GCC >= 11 or Clang >= 14
- CMake >= 3.21
- Ninja or Make
- MPI, such as OpenMPI or MPICH
- PETSc real scalar build

### Dependency Installation

Debian/Ubuntu:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build openmpi-bin libopenmpi-dev libpetsc-real-dev
```

Spack alternative:

```bash
spack install gcc cmake ninja openmpi petsc
spack load gcc cmake ninja openmpi petsc
```

### Configure

```bash
cmake -S . -B build-linux -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DPETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

Adjust `PETSC_DIR` to match your PETSc installation.

### Build

```bash
cmake --build build-linux -j
```

### Run

```bash
./build-linux/pancake -c config.txt
mpirun -n 2 ./build-linux/pancake -c config.txt
```

### Tests

CTest registration is not present yet, so run tests directly:

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

## Windows

Native Windows builds are planned alongside the WSL/Linux build. Use a normal
Windows filesystem path such as `C:\src\pancake`; avoid building through
`\\wsl$` paths.

### Prerequisites

- Visual Studio 2022 Community or newer
- Workload: Desktop development with C++
- Components: MSVC v143, Windows 11 SDK, CMake tools for Windows
- Git for Windows
- vcpkg at `C:\vcpkg`
- Optional: Windows Terminal

### Dependency Installation

Install vcpkg:

```bat
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg integrate install
setx VCPKG_ROOT C:\vcpkg
```

Install the dependency that vcpkg can provide for the current source:

```bat
C:\vcpkg\vcpkg install mpi:x64-windows
```

The current code also requires PETSc, but the audit did not find an official
vcpkg PETSc port. A first native Windows build is blocked until either:

- PETSc is built from source for MSVC and MS-MPI, then passed to CMake through
  `PETSC_DIR`, or
- the planned Trilinos backend is implemented and the CMake files use
  `find_package(Trilinos)` for the exact packages used by the source.

Do not introduce a Windows-only dependency path; Linux must remain first-class.

### Configure

After the linear algebra dependency blocker is resolved:

```bat
cmake -S . -B build-windows -G "Visual Studio 17 2022" -A x64 ^
  -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake ^
  -DVCPKG_TARGET_TRIPLET=x64-windows ^
  -DPETSC_DIR=C:\path\to\petsc
```

If the solver has been migrated to Trilinos by then, replace `PETSC_DIR` with
the relevant `Trilinos_DIR` or vcpkg/source-build toolchain settings.

### Build

```bat
cmake --build build-windows --config Debug
```

### Run

```bat
build-windows\Debug\pancake.exe -c config.txt
mpiexec -n 2 build-windows\Debug\pancake.exe -c config.txt
```

## Future Qt/VTK GUI

The planned GUI should link a cross-platform core solver library instead of
duplicating solver code in the executable. vcpkg can provide Qt and VTK for the
GUI phase:

```bat
C:\vcpkg\vcpkg install qtbase[widgets,opengl]:x64-windows vtk[qt,utf8]:x64-windows
```

On Linux, use distro, Spack, or vcpkg packages consistently with the selected
developer environment. Keep GUI code in a separate executable target that links
the solver library.
