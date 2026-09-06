# Physics repair comparison

Compare pre-audit `b9dcaf3` with repaired `2bc4d3c` using the same numerical
inputs. [SUMMARY.md](SUMMARY.md) explains the measured changes and the limits
of the repaired model. This case adds measurements and documentation; it does
not change either solver revision.

| Input | Purpose | Resolution / duration |
|---|---|---|
| `ferron.txt` | Archived conjugate steady iterate, fixed geometry and supply flow | 180 × 48 |
| `optical.txt` | Optical-rig P1 condition, prescribed 0.8 bar absolute plateau | 720 × 192 |
| `air.txt` | Synthetic finite-rate oil-air release with flow and energy | 48 × 8, 10 µs steps to 1 ms |
| `timing.txt` | Exact insulated Couette heating and decoupled axial motion | 24 × 4, 10 ms steps to 15 ms |
| `phase_exchange_probe.cpp` | Closed local reaction at fixed pressure and temperature | 100 µs steps to 20 ms |

`input_provenance.json` records the original Ferron input and the read-only
optical template hashes, plus the archived optical summary row. The optical
input freezes R = 24.65 mm, c = 350 µm, L = 50 mm, e/c = 0.5, 1500 rpm,
mu = 0.068 Pa s and end pressure = 101325 Pa. The viscosity is the existing
campaign estimate, not a newly measured oil property. Thermal options omitted
from that isothermal template are inert under this configuration.

The local reaction probe calls the same public functions in both revisions.
It has no transport or reservoir and measures liquid, oil, dissolved and free
inventories directly. Initial liquid fraction is 0.5, leaving room for the
released gas. It is a numerical conservation check, not an experimental bearing.

## Reproduce

Use two clean source checkouts of the stated commits, and build them with the
same compiler and PETSc/MPI libraries. The commands below use shell variables
as placeholders for absolute paths; use the corresponding PowerShell syntax
on Windows. Follow [BUILDING.md](../../../BUILDING.md) for platform dependencies.

```sh
git worktree add --detach "$previous_source" b9dcaf3
git worktree add --detach "$current_source" 2bc4d3c
cmake -S "$previous_source" -B "$previous_build" -DCMAKE_BUILD_TYPE=Release -DPANCAKE_BUILD_GUI=OFF -DPANCAKE_BUILD_TESTS=OFF -DCMAKE_PROJECT_INCLUDE="$case_dir/probe_targets.cmake"
cmake -S "$current_source" -B "$current_build" -DCMAKE_BUILD_TYPE=Release -DPANCAKE_BUILD_GUI=OFF -DPANCAKE_BUILD_TESTS=OFF -DCMAKE_PROJECT_INCLUDE="$case_dir/probe_targets.cmake"
cmake --build "$previous_build" --target pancake phase_exchange_probe -j 4
cmake --build "$current_build" --target pancake phase_exchange_probe -j 4
python "$case_dir/run_comparison.py" --previous "$previous_build/pancake" --current "$current_build/pancake" --work "$new_work_dir"
python "$case_dir/run_refinement.py" --solver "$current_build/pancake" --work "$new_work_dir"
python "$case_dir/reduce_comparison.py" --work "$new_work_dir"
python "$case_dir/plot_comparison.py"
```

Use `.exe` for Windows executables. The runner fixes OMP/OpenBLAS threads to
one and launches one MPI rank. Supply a **new isolated work directory** because
the historical executable removes its configured output directory. Production
results and the dirty Ferron/optical checkout were not modified. No remote run,
mirror, push or fetch was used.

Plotting requires NumPy and Matplotlib. For this task they were installed under
the ignored `build/plot-packages` directory with `pip --target`; setting
`PYTHONPATH` to that directory makes them available without a global install.
The saved figures used NumPy 2.5.2 and Matplotlib 3.11.1. `plot_style.py` is a
copy of the user's publication style. Every figure is a single 300 dpi PNG.
Plotting the archived CSVs needs no solver run. The CMake project-include hook
adds only the common probe; it does not patch either core.

## Evidence and checks

- `evidence/run_manifest.json` records input and executable hashes, return
  codes, step counts and observed elapsed times. A null outer-coupling verdict
  means no transient fixed-point iteration was performed.
- `evidence/summary.json` contains field differences and direct inventory
  checks. CSVs retain both the old time labels and corrected comparison times.
- Only the static oil-air comparison shifts old labels by +10 µs. Both curves
  then cover the same 100 physical advances. The old extra step at 1.01 ms is
  retained in the CSV and excluded from aligned metrics and figures.
- Steady profiles use the center-adjacent axial row on an even grid. Ferron
  angle is measured from its 90° groove. Optical angle uses the mesh origin.
- `pressure_load_N` is the norm of `Fp`, as in the archived campaign. `load_N`
  is the norm of total fluid force `F`, including shear. These vectors are
  repeated across cells; they must not be summed over the mesh.
- Liquid-deficit arc thresholds and axial selection are explicit in metric
  names. In gas-active cases this deficit does not by itself distinguish vapor,
  released gas and numerical volume error.
- The reducer checks finite VTK fields, matching inputs, inner convergence,
  aligned times, direct oil/gas conservation and independent exact reaction
  and timing solutions. Saved VTK values have about six significant digits;
  equal fields mean equal at that output precision. Direct probe CSVs use
  17-digit output.
- `air_refinement` preserves every step at 5 and 2.5 µs. The timestep check
  exposes a failure, so it is not an assertion that the load curves converge.

Raw VTK remains in the ignored isolated build directory. Portable profiles,
histories, diagnostics, inputs, manifests and logs are archived here. This
comparison does not replace the historical experimental validation scores.
