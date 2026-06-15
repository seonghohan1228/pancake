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

Copy these two files into one folder:

```
pancake_gui.exe
pancake.exe        (plus the solver's own runtime DLLs - PETSc/MS-MPI/MinGW -
                    as deployed by the solver build)
```

The GUI looks for `pancake.exe` next to itself first; a different path can be
set in *Tools > Settings* inside the app.

## Solver + tests (full tree)

The existing full build is unchanged:

```
cmake --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
```

This produces `build-windows-mingw/pancake.exe`, `pancake_gui.exe`, and the
test suite (see `docs/WINDOWS_PORT.md`).
