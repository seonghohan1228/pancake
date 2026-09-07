# Boundary, kinetic and material checks

2026-09-07. Interpretation: [FW68D cavitation review](../../../docs/FW68D_CAVITATION_REVIEW.md).
These checks do **not** establish a working R410A bearing solver or repair the
pressure/volume instability recorded in [the paired comparison](../physics_comparison/SUMMARY.md).

## Evidence and figures

* `kinetics_probe.cpp` calls the current C++ phase-exchange routine on a
  prescribed-pressure parcel. `evidence/kinetics.csv` records the full histories
  for equal release/absorption rates and for zero absorption. Both use the
  conservative repaired code; “equal rates” is not the old mass-losing revision.
  Pressure is 1 bar to 20 ms, ramps linearly to 4 bar at 30 ms, then stays there
  to 40 ms. Temperature is 313.15 K, initial liquid content 0.5, oil density
  950 kg/m³, initial dissolved mass fraction 0.002, Henry coefficient 1e-8 Pa⁻¹,
  release 20 s⁻¹ and timestep 0.1 ms. Rates and properties are illustrative.
  Gap volume is not constrained by a bearing pressure solve in this probe.
* `evidence/r410a_eos.csv` is generated with CoolProp 7.2.0 `HEOS::R410A`.
  It samples 20/40/60/80 °C at 1–20 bar absolute and records phase identity,
  real density, ideal R410A density and the incorrect propane substitution.
  It contains pure refrigerant data, **not** FW68D solubility or solution
  properties. The CSV is not wired into the solver.
* `figures/absorption_vs_compression.png` shows free-gas mass and volume for
  those two kinetic choices. Gas disappears through absorption in one curve;
  compression alone retains its mass in the other.
* `figures/r410a_eos.png` shows gaseous states only, with lighter blue for
  increasing temperature. Relative error is `100*(rho_ideal/rho_EOS-1)`.
  The 20 °C curve stops before the nongaseous states. It is not the solution
  phase boundary.
* `figures/etsion_reformation.png` replots the **archived** `05_etsion_extent`
  FBNS results and measured angles, with the original ±5° systematic angular
  uncertainty. A measured angle of 355° is shown as −5°. No new bearing runs
  or new gas-validation scores are implied by this figure.

`evidence/checks.json` stores EOS provenance and kinetic inventory drift;
`evidence/provenance.json` identifies the exact tested source and binaries.

## Tests

The source changes add analytical reversal checks to `test_gas_transport` and
`test_energy`, serialization/invalid-state checks to `test_config_keys`, and
independent absorption/compression checks to `test_physics_audit`.

* Native Windows: all 40 registered tests pass; full log in `ctest_windows.log`.
  The new gas-transport reversal test also passes at one and four MPI ranks.
* Linux/WSL: 37 tests pass on the first run; three six-rank processes cannot
  launch because the default advertised core slots are fewer than six.
  `ctest_linux_mpi6_retry.txt` records all three passing with
  `OMPI_MCA_rmaps_base_oversubscribe=1`. Both logs are retained. This was a
  launcher resource limit, not a numerical failure. All 40 tests have therefore
  executed successfully on Linux.

The reversal test imposes a linear axial pressure field between distinct
reservoirs, swaps the pressure ordering, and compares every cell with the
backward-Euler upwind recurrence. Its integrated dissolved-plus-free inventory
change must equal the applied boundary flux. The thermal test includes the
saved viscous heating source in the same recurrence. These are operator
verification tests, not empirical bearing calibration.

Species reservoir inputs operate in transient `TWO_PHASE` transport. They do
not add conservative species transport to the steady local feed-flash model.

## Reproduce

From the repository root, configure the optional diagnostic executable:

```text
cmake -S . -B <build> -DCMAKE_PROJECT_INCLUDE=<absolute-path>/cases/90_diagnostics/cavitation_requirements/probe_targets.cmake
cmake --build <build> --target kinetics_probe
<build>/kinetics_probe cases/90_diagnostics/cavitation_requirements/evidence/kinetics.csv
python cases/90_diagnostics/cavitation_requirements/plot_checks.py
```

Use `.exe` on Windows and make MPI/PETSc runtime libraries discoverable as
described in `BUILDING.md`. The plot script requires NumPy, Matplotlib and
CoolProp 7.2.0. For this run they were installed in ignored task-local
`build/plot-packages` and `build/thermo-packages`, exposed via `PYTHONPATH`.
They are not dependencies of the C++ bearing solver. Figures use the lab's
shared publication style and are individually inspected after rendering.

The optional CMake hook adds only the diagnostic executable; it does not change
the production solver. Remove `CMAKE_PROJECT_INCLUDE` from the build cache if
the probe is no longer needed.
