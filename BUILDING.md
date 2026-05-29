# Building

`pancake` builds from one CMake source tree. The current native Windows
development workflow is MSYS2 MINGW64. Visual Studio is not required for the
current solver or GUI build.

## Windows: MSYS2 MINGW64

### Prerequisites

Install Microsoft MPI Runtime and SDK:

```text
https://www.microsoft.com/en-us/download/details.aspx?id=105289
```

Install both `msmpisetup.exe` and `msmpisdk.msi`.

Install MSYS2:

```text
https://www.msys2.org/
```

Open **MSYS2 MINGW64** and update packages:

```bash
pacman -Syu
```

If prompted, close and reopen **MSYS2 MINGW64**, then run:

```bash
pacman -Syu
```

Install build dependencies:

```bash
pacman -S --needed git base-devel mingw-w64-x86_64-toolchain \
  mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja \
  mingw-w64-x86_64-pkgconf mingw-w64-x86_64-msmpi \
  mingw-w64-x86_64-petsc
```

### Clone

```bash
mkdir -p /c/Dev
cd /c/Dev
git clone https://github.com/seonghohan1228/pancake.git
cd /c/Dev/pancake
```

For an existing checkout:

```bash
cd /c/Dev/pancake
```

### Configure And Build

```bash
cmake --fresh --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

The same commands can be run from PowerShell in `C:\Dev\pancake`. The preset
prepends `C:\msys64\mingw64\bin` to `PATH`, so native compiler helper programs
and MinGW DLLs are discoverable without launching a separate MSYS2 shell.

Expected configure output:

```text
PETSc pkg-config module: petsc-dmo
```

If configure output mentions `petsc-sso`, CMake found the serial PETSc package.
Run the fresh configure command above again from **MSYS2 MINGW64**.

Build outputs:

```text
build-windows-mingw/pancake.exe
build-windows-mingw/pancake_gui.exe
```

The build copies required runtime DLLs, including `libpetsc-dmo.dll`, beside
`pancake.exe`, `pancake_gui.exe`, and the test executables so File Explorer
launches and plain `ctest` work.

### Run

```bash
cd /c/Dev/pancake/build-windows-mingw
./pancake.exe -c config.txt
mpiexec -n 2 ./pancake.exe -c config.txt
./pancake_gui.exe
```

### Test

```bash
cd /c/Dev/pancake
ctest --test-dir build-windows-mingw --output-on-failure
```

Known current status: `test_elrod_a1` fails its pressure regression threshold;
the other native tests pass after the MinGW/PETSc runtime is available.

## Linux

Linux remains supported. Keep it in a separate build directory:

```bash
cmake -S . -B build-linux -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build-linux -j
./build-linux/pancake -c config.txt
mpirun -n 2 ./build-linux/pancake -c config.txt
ctest --test-dir build-linux --output-on-failure
```

Install typical Debian/Ubuntu dependencies with:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build openmpi-bin libopenmpi-dev libpetsc-real-dev
```

Use `-DPETSC_DIR=/path/to/petsc` if your PETSc installation is not discoverable
by the default CMake search.
