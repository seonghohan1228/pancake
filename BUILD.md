# Building the Pancake GUI

`pancake_gui.exe` is a single statically linked Windows executable (Dear ImGui
docking + ImPlot on Win32/Direct3D 11). It ships together with the solver
`pancake.exe` in one folder; nothing else is required at runtime. Settings and
window layout are written to `%LOCALAPPDATA%\PancakeGui\`, never next to the
exe, so it runs from read-only locations.

## Prerequisites

Either toolchain works; the build is the same CMake tree.

**MSVC (the `release` preset)**
- Visual Studio 2022 Build Tools (or newer) with the "Desktop development
  with C++" workload (provides `cl.exe`, x64)
- CMake >= 3.21 and Ninja (both included in the VS Build Tools installer)
- Run the commands below from an **x64 Native Tools Command Prompt** so
  `cl.exe` is on PATH.

**MinGW (the `release-mingw` preset)**
- MSYS2 with the `mingw-w64-x86_64-toolchain` group installed at `C:\msys64`
- CMake >= 3.21 and Ninja on PATH

Internet access is needed at **configure** time only: CMake FetchContent
downloads pinned releases of Dear ImGui (v1.91.9b, docking branch) and ImPlot
(v0.17). Fonts and stb are vendored in `third_party/`.

## Build

```
cmake --preset release          # or: cmake --preset release-mingw
cmake --build --preset release  # or: cmake --build --preset release-mingw
```

The executable lands at:

```
build/release/pancake_gui.exe         (MSVC)
build/release-mingw/pancake_gui.exe   (MinGW)
```

Both presets enforce a static C/C++ runtime (`/MT` for MSVC, `-static` for
MinGW): the import table contains only DLLs that ship with Windows (kernel32,
user32, gdi32, shell32, ole32, comdlg32, dwmapi, d3d11, d3dcompiler_47,
msvcrt). Verify with `dumpbin /DEPENDENTS pancake_gui.exe` or
`objdump -p pancake_gui.exe | findstr "DLL Name"`.

## Distribution

### Two-file build (no DLLs)

With the static libraries of MSYS2 PETSc and the Fortran runtime
(`pacman -S mingw-w64-x86_64-gcc-fortran`, matching the installed GCC), the
solver links statically:

```
cmake -S . -B build-static -G Ninja -DCMAKE_BUILD_TYPE=Release -DPANCAKE_STATIC_SOLVER=ON
cmake --build build-static
```

`pancake.exe` (about 60 MB) then imports only Windows DLLs and `msmpi.dll`,
which the Microsoft MPI runtime installs into Windows. A working folder is
`pancake.exe`, `pancake_gui.exe`, `fluids\` and `examples\`.

### Release folder

Build the full tree (solver + GUI), then collect the release folder:

```
powershell -File tools/make_release.ps1            # -> release/
```

The folder holds everything an engineer needs besides Microsoft MPI:

| Item | Why |
|---|---|
| `pancake_gui.exe` | GUI (static; Windows DLLs only) |
| `pancake.exe` | solver |
| `libpetsc-dmo.dll`, `libopenblas.dll`, `libmetis.dll`, `libparmetis.dll`, `libgfortran-5.dll`, `libgomp-1.dll`, `libquadmath-0.dll`, `libhwloc-15.dll`, `libltdl-7.dll`, `libstdc++-6.dll`, `libgcc_s_seh-1.dll`, `libwinpthread-1.dll` | solver runtime (MinGW, PETSc, BLAS, METIS), deployed by the build |
| `fluids/` | oil-refrigerant pairs; the shipped `example_vg68_rx` is synthetic |
| `examples/` | example cases and the example load table |
| `README.txt` | what to install and how to start |

Install separately: the Microsoft MPI runtime (`msmpisetup.exe`). It provides
`msmpi.dll`, which `pancake.exe` needs even on one core, and `mpiexec` for more
than one core. Not shipped: `data/` and `tests/data/lge/` (LGE test data).

The GUI looks for `pancake.exe` next to itself first; another path can be set
in *Tools > Settings*, together with the fluid directory.

## Solver + tests (full tree)

The existing full build is unchanged:

```
cmake --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

This produces `build-windows-mingw/pancake.exe`, `pancake_gui.exe`, and the
test suite (see `docs/WINDOWS_PORT.md`).
