# September 2026 physics verification

This case records numerical verification after the conservative inventory,
phase-exchange and gas-transport repairs. It contains **no experimental
calibration** for oil-air or oil-R410A. See the [audit](../../../docs/PHYSICS_AUDIT_2026-09.md)
and [physical derivation](../../../PHYSICS.md#18-conservative-gas-update-and-experimental-scope-2026-09-06).

## Reproduce

Run from the repository root using a built solver and Python with NumPy. Each
runner writes configs, solver logs and a summary below the requested output
directory; solver products are replaced on rerun.

```sh
python cases/90_diagnostics/physics_audit/run_verification.py \
  --solver build-linux/pancake --output build-linux/audit-verification

python cases/90_diagnostics/physics_audit/run_air_smoke.py \
  --solver build-linux/pancake --output build-linux/audit-air
```

On native Windows, use `build-windows-mingw/pancake.exe` and, if needed,
`--mpiexec "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"`. The air runner
launches 1, 2 and 4 ranks sequentially. The full CTest suite also includes
6-rank tests; OpenMPI hosts with fewer slots may need
`OMPI_MCA_rmaps_base_oversubscribe=1`. Both runners set BLAS/OpenMP threads to one.

Run the complete regression suite separately:

```sh
ctest --test-dir build-linux --output-on-failure
```

## Criteria and limits

`run_verification.py` reuses the four existing exact 1-D FBNS configs and their
extractor from `cases/00c_jfo_verification`. Pass requires rupture/reformation
error at most one cell and absolute peak/load error below 1%. The exact reference
is rupture 302.766608 degrees, reformation 79.006725 degrees, peak absolute
pressure 7.597141 MPa and load per length 495463.6 N/m. The extraction uses the
cell-face front estimator and midpoint pressure integration already registered
for that benchmark. Wall time includes launch, solution and VTK output.

`air_smoke.txt` is a synthetic 48x8 transient at 200 rad/s with a 200 kPa groove,
100 kPa ends, 40 C supply/walls, 0.2% initial dissolved air, a chosen Henry slope
and a 20/s exchange rate. It exercises the CLI flow/transport/reaction/energy
sequence. Numerical pass requires all inner solves to converge, normalized gas
mass residual below 1e-8, and cross-rank final-field differences below 1e-5.
The VTK output has limited decimal precision, so identical saved fields do not
imply bit-identical internal fields. Unit tests check the full-precision
inventories independently.

**Volume errors are reported, not accepted as physically valid.** The air case
exhibits positive `gas_volume_excess`; its numerical pass does not certify the
segregated pressure-volume closure. The corrected driver records the initial
VTK state and ten physical advances through 0.1 ms. Diagnostics contain the ten
completed steps. These synthetic coefficients must not be used to identify
experimental time constants. The separate CLI time regression checks uniform
Couette heating and rigid-body acceleration through a shortened final interval.

## Recorded evidence

- [Exact 1-D numerical results](verification_windows.json): solver SHA-256,
  reference errors, residuals and timings for the native Windows rerun.
- [Synthetic air numerical results](air_windows.json): solver SHA-256, MPI
  comparisons, mass and volume errors, and final field ranges.
- [Final Windows CTest output](ctest_windows.txt) and
  [final Linux CTest output](ctest_linux.txt): all 40 tests and outcomes.

The evidence is tied to the executable hashes. The source commit containing
this directory is the audit revision; rebuilding on another toolchain need not
reproduce the binary hash. Full build directories and VTK fields are local
artifacts and are not committed. Historical Ferron or Etsion results are
reviewed in the audit; they were not silently rerun or revalidated by these tests.
