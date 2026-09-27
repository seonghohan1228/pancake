# FW68D steady execution record — 2026-09-27

The six user-selected conditions and two E temperature variants are specified
in [the case plan](README.md). The production batch started at 19:33 KST. No load/thermal equilibrium had
been accepted at launch; the fixed-position numerical check below passed.

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

## Fixed-position checks and current execution

The first long F trial used snapshot `875a58b`, eccentricity ratio 0.2 and
attitude 165 degrees. It returned code 0 in 262 iterations, but its pressure
residual dropped from 6.49e-5 to exactly zero because FBNS accepted its initial
iterate at the default 1e-9 tolerance. Retain this as a tolerance diagnostic;
it does not demonstrate the requested 1e-6 outer pressure tolerance.

The corrected `ae0efd2` trial used inner FBNS tolerance 1e-12. It converged in
484 iterations and 163.84 s on 24 pinned ranks. The final pressure change was
9.86777e-7, down smoothly from 1.00547e-6, after an actual FBNS Newton update.
The final temperature change was 3.56775e-9 and FBNS residual 5.38e-18. Two
linear-solver fallback recoveries succeeded; return code was 0.

The film maximum was 28.9183 °C, shaft 28.8635 °C and mean bore 28.8762 °C,
with TC0 = 28.3716 °C. These are fixed-position, frozen-coefficient diagnostic
values, not the final Case F prediction. Housing h remained 51.4% from its
updated estimate, and load balance remained open. Both solid wall flux
integrals match imposed external convection to within 5e-10 W. This does not
supply the missing global fluid energy-balance diagnostic.

The pinned groove still produces a local bore undershoot, about 0.216 K below
TC0 in the corrected trial. The sheared quartic reconstruction and fixed mean
inlet temperature are inconsistent locally. Keep this baseline limitation
explicit; do not clamp the field or interpret groove-region wall temperature
as physical sub-ambient cooling. Influence on other results requires a
separate sensitivity check.

Complete corrected preflight evidence is local under
`C:/Dev/pancake/results/fw68d_steady_20260927/remote_ae0efd2/`.
The archive SHA256 is
`7967e832a82e1e446d0ee94b42dd4dc2b0b355ee8e5424c21fb0e4fe3550abee`.
The inspected convergence plot is `preflight_F/F/0001_trial_thermal1_dfbc2e12/F_fixed_position_convergence.png`.

## Production launch

The eight-case batch started on 2026-09-27 at 10:33:16 UTC, 19:33:16 KST,
with launcher PID 3922893. Source snapshot `ae0efd2` is isolated at
`/home/user/pancake-runs/fw68d-steady-20260927-ae0efd2/`.
Its source archive SHA256 is
`72ad3ed23fd72039e6bf5a7559c89d9a161b51f01076d3165433d4b245163875`.
The executable remains `a3141fce7bdab6aaa2f4c98abbb1a400ff32e183a23fe64c8b7c5d596f6f3e1e`.

A-F started on six disjoint four-core slots, and E_low/E_high are queued.
The verified groups are 0-3, 4-7, 8-11, 14-17, 18-21 and 22-25. All six
initial logs report four bound ranks. Each case independently closes external
h and load balance; failure preserves its evidence and does not stop the other
cases. The batch is a numerical baseline with the thermal model limitations
above, not an experimentally validated result.

Live status and evidence are in `production/<case>/batch_status.json`,
`production/<case>/driver.log` and each numbered CFD evaluation directory.
`production/batch_manifest.json` records completed cases, while per-case status
files distinguish running cases from the waiting queue. `batch_launch.json`
records the launch command and PID. The prepared launcher is retained locally
at `C:/Dev/pancake/results/fw68d_steady_20260927/launch_batch.py`.
