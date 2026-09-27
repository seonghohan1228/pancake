# FW68D steady execution record — 2026-09-27

The six user-selected conditions and two E temperature variants are specified
in [the case plan](README.md). No load/thermal equilibrium has yet been accepted.

## Reproducible snapshots and checks

- Native GCC/MinGW build and focused one/four-rank stress/cavitation checks pass.
- Remote `f754147` and `5163417` snapshots: fresh Release builds and focused tests pass.
  The final snapshot also passes explicit four-rank force and FBNS tests.
- Source archive SHA256: `8804dc0aa8bcbce2f7c167edaf5c239885dddb1d0226f002c608f1b0f3aec08d`.
- Remote executable SHA256: `a3141fce7bdab6aaa2f4c98abbb1a400ff32e183a23fe64c8b7c5d596f6f3e1e`.
- Workstation: afml-r6, GCC 11.4, OpenMPI 4.1.2, PETSc 3.15.5.
- Corrected driver/configuration snapshot `5163417` has identical C++ source.
  Its archive SHA256 is `304432e258a7a3a6a192d551dd087ac71fed99dfcbb67c8fffebf86537cf74ad`.

## Brief multicore test

Two repetitions per rank count, fixed-position F, 720×192, 16 bush radial
cells, quartic cross-film profile, 12 outer sweeps. Every rank occupies a
separate physical core; OMP/OpenBLAS/MKL threads are one. All 10 trials
returned the expected nonconverged code 4 at the deliberate cap. These timings
measure startup work, not time to a converged equilibrium. Final pressure
residuals agree within 7e-8 and temperature residuals agree to displayed
precision across rank counts.

| MPI ranks | Median wall time (s) | Median core-seconds |
|---:|---:|---:|
| 2 | 75.20 | 150.39 |
| 4 | 39.55 | 158.19 |
| 6 | 29.22 | 175.32 |
| 12 | 17.99 | 215.87 |
| 24 | 11.26 | 270.34 |

Use four ranks per case for batch throughput: nearly half the two-rank wall
time for about 5% more core-seconds. Six such cases fit on 24 physical cores.
Use 24 ranks when reducing the latency of one diagnostic is more useful.
The optimum can differ in the later nonlinear iterations and with concurrent
memory traffic; this is the supported choice from the requested brief test.

Full inputs, logs, rankfiles, diagnostics, fields and hashes are retained in
`C:/Dev/pancake/results/fw68d_steady_20260927/` and remotely under
`/home/user/pancake-runs/fw68d-steady-20260927-f754147/benchmarks/benchmark_20260927T100147Z_3904602`.

## Invalid and capped checks

Two native 200-sweep attempts returned code 4 and inherited an incorrect
313.15 K housing reservoir. They are invalid as requested physical cases and
remain saved. A corrected native 12-sweep check also returned code 4 at its
intended cap, exercised the output parser, and was correctly rejected by the
driver. None supplies a temperature/friction/cavity comparison result.

The current cross-film path does not expose a global heat-balance residual.
Pressure/temperature convergence, wall fluxes and mass diagnostics are checked
separately; energy-balance verification remains unavailable.

## MPI binding correction

The workstation OpenMPI 4.1.2 left all ranks unbound when given `--cpu-list`.
Explicit rankfiles were verified to bind ranks to separate physical cores,
including disjoint groups 0–3 and 4–7. The driver uses `--rankfile` and records
its hash and binding output. This changes process placement only.

## Fixed-position preflight and inner-tolerance correction

Snapshot875a58b completed a24-rank F trial at epsilon0.2 and attitude165 degrees
in262 outer iterations,101.71s, return code0. It remained a fixed-coefficient,
fixed-load-position trial: housing-h mismatch51.3%, force mismatch6.35%. Two
linear-solver fallbacks recovered successfully. Solid wall flux integrals agree
with their imposed external cooling to within5e-10W.

The convergence curve revealed that the final pressure residual dropped from
6.49e-5 to exactly zero because FBNS accepted its initial iterate at its default
1e-9 tolerance. This does not demonstrate the requested1e-6 outer pressure
tolerance. Retain that trial as a tolerance diagnostic. Tighten inner FBNS to
1e-12 and repeat before accepting pressure convergence or launching production.

The trial also shows a localized bore undershoot of0.218K below TC0 in the
pinned groove, extending0.069K below TC0 into two neighboring columns. This
requires a model review before interpreting local groove-region temperature.
No clamping or adjustment of those saved temperatures has been applied.
