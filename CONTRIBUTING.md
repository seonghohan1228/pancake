# Contributing

`pancake` is cross-platform by policy. Changes should preserve Linux and
Windows builds from the same CMake source tree.

## Git Workflow

Use one task, one branch, one verified merge, then delete the branch.

1. Begin with a clean, current `main` and create `codex/<short-scope>`.
2. Keep the branch limited to one coherent change. Never commit local agent
   settings, build directories, generated results, or unrelated edits.
3. Preserve unrelated work already in the checkout; do not silently absorb or
   discard it.
4. Run the relevant build/tests and update the required documentation before
   committing with an imperative, scoped message.
5. Push or open a pull request only when the repository owner requests it.
   Prefer squash or fast-forward integration; never rewrite shared `main`.
6. Delete the task branch locally and remotely after it is merged.

## Code Conventions

- Use `std::filesystem::path` for filesystem paths.
- Do not build paths by concatenating raw strings with `/` or `\`.
- Do not use POSIX-only headers in core code. If a feature requires a platform
  API, isolate it behind a platform abstraction and provide a Windows
  implementation.
- Prefer `std::int64_t`, `std::int32_t`, and `std::size_t` for values with
  meaningful width or size semantics.
- Avoid `long` for 64-bit data. It is 8 bytes on typical Linux x64 builds and
  4 bytes on Windows x64.
- All new code must compile with GCC and MSVC. Clang support should be kept
  clean as well.
- Warnings-as-errors are recommended on both GCC/Clang and MSVC once the
  current portability cleanup is complete.
- Tests must pass on both platforms before merge.

## Build Discipline

- Use out-of-tree builds: `build-linux/` for Linux and `build-windows-mingw/`
  for the current native Windows workflow.
- Keep dependencies discoverable through CMake targets. Avoid hardcoded
  `/usr`, `/opt`, or local Windows paths in committed CMake files.
- Guard compiler-specific flags with compiler checks, for example `if(MSVC)`
  or `if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")`.
- Prefer target-scoped CMake commands such as `target_compile_options`,
  `target_include_directories`, and `target_link_libraries`.

## Documentation

- Update [README.md](README.md) when user-facing build, run, or dependency
  behavior changes.
- Update [PHYSICS.md](PHYSICS.md) when the mathematical model or numerical
  method changes.
- Update [PLAN.md](PLAN.md) when completing or changing planned development
  phases. Do not delete completed tasks.
- Append meaningful entries to [CHANGES.md](CHANGES.md) so prior edits and
  their reasons remain reviewable.

## Review Checklist

- Linux build still configures and builds from a clean `build-linux/`.
- Windows build impact has been considered, especially for file I/O,
  dependencies, compiler flags, and integer widths.
- New output files use stable relative paths and do not depend on the current
  working directory unexpectedly.
- Reference simulations produce comparable PVD/VTK output after numerical
  tolerance is applied.
