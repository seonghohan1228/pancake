# Architecture

`pancake` should remain a cross-platform C++ solver with a clear boundary
between numerical core code, user interfaces, and file output.

## Intended Layers

```text
Cross-platform core solver library
  - mesh and domain decomposition
  - fields and ghost-cell exchange
  - finite-volume operators
  - Reynolds solvers and post-processing
  - linear algebra backend wrapper

Command-line executable
  - parses configuration
  - owns the time loop
  - writes VTK/PVD output

Native Windows GUI executable
  - exposes tabbed parameter controls backed by SimulationConfig
  - edits config.txt next to pancake_gui.exe
  - writes output-selection settings for VTK/PVD products and field names
  - can launch or stop pancake.exe without closing the GUI
  - streams solver logs and previews ASCII VTK scalar fields by timestep
  - rotates the preview around a load-reference angle measured from +y

Future Qt + VTK result-viewer GUI executable
  - owns richer result viewers and field inspection tools
  - links the core solver library
  - embeds VTK widgets for contours and field inspection

External visualization
  - PVD/VTK output remains consumable by ParaView
```

The CMake tree exposes the solver as `pancake_core`, keeps `src/main.cpp` as a
thin CLI executable, and builds the Windows-only `pancake_gui` target separately
from the MSYS2 MinGW64 workflow. The GUI remains an external orchestrator: it
writes config, launches the solver process, and reads VTK files rather than
duplicating Reynolds solver logic.

## Portability Rules

- `std::filesystem::path` is the canonical type for paths.
- Core solver code must not use POSIX-only headers such as `<unistd.h>`,
  `<dirent.h>`, `<sys/stat.h>`, or `<pwd.h>`.
- Platform-specific behavior belongs behind a small abstraction with Linux and
  Windows implementations.
- New code must compile with GCC, Clang, and MSVC.
- The Linux build must not be weakened to make Windows easier.

## Linear Algebra

The checked-in source currently wraps PETSc in `src/linear_system.hpp` and
`src/linear_system.cpp`. The strategic porting target identifies Trilinos as
the fixed long-term numerical library. Treat that mismatch as an architectural
blocker: the source tree should expose one solver-facing linear-system wrapper,
with backend details isolated behind that wrapper and configured portably by
CMake.

## File Output

PVD/VTK output is part of the stable interoperability surface. The CLI and
future result-viewer GUI may both produce output, but path creation and file
naming should flow through `std::filesystem::path` so Windows, Linux, and
non-ASCII paths are handled consistently. Runtime configuration defaults to
`config.txt` beside the executable unless the CLI is given `-c` or `--config`.
