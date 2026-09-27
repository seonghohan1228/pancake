# FW68D transient CHT campaign and MPI benchmark

The user replaced this 24-point transient campaign plan with the
[six-case steady A–F campaign](../../02d_fw68d_steady_six/README.md).
The benchmark and failed-startup evidence below remain historical records.

This case plans the FW68D-only cavitation comparison and measures workstation
MPI performance on a short synthetic case. See [measured results](RESULTS.md) for
the current G6 test outcome and the archived coarse-grid timing comparison.
The G6 startup hit the minimum-film guard, so its MPI optimum is unresolved.
The benchmark is not an experimental
prediction. Source physics is unchanged from the recorded source commit.

## Comparison matrix

| Speed (rpm) | Target Sommerfeld numbers | Base cases |
|---|---|---:|
| 200 | 0.030, 0.032, 0.035, 0.070, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90, 1.10, 1.30, 1.45 | 15 |
| 1500 | 0.063, 0.075, 0.090, 0.113, 0.70, 0.80, 0.90, 1.10, 1.21 | 9 |

Start with the four outer endpoints: (200 rpm, 0.030), (200 rpm, 1.45),
(1500 rpm, 0.063), and (1500 rpm, 1.21). They count toward the 24 base cases.
Repeat each on a finer spatial mesh and separately at half the timestep.
The initial budget is 32 campaign executions, plus short numerical checks.
Further refinements are conditional on failing the acceptance targets.

Use the previously accepted G6 film mesh (720 by 192) and dt = 2e-5 s (20 us).
Retain the earlier CHT settings: 16 radial bush cells, polynomial order 4 and
six quadrature points; their independent convergence was not established.
The endpoint spatial repeats use G7, 1080 by 288, with the same radial and
cross-film settings; the temporal repeats use dt = 1e-5 s (10 us).
These inherit the [prior convergence protocol](../../01e_ferron_terra_repro/README.md).
Use the endpoint checks to assess transfer to the present FW68D conditions.
Compare friction and mean force
within 1%, peak/sensor temperature within 0.5 K, cavity area within one percentage
point, front locations within 2 degrees, and minimum gap within 2%.

Require finite fields, positive gap, converged field solves and satisfactory
conservation. The bearing equation includes inertia and any supports, so do not
mistake an instantaneous dynamic force imbalance for failure of static force
balance. Compare the full equation of motion; apply a 0.2% static force-closure
target only after settling. Compare statistics across consecutive full-revolution
windows (at least 0.3 s at 200 rpm and 0.04 s at 1500 rpm), extending the windows
for slower orbital dynamics. Target less than 1% mean friction/cavity-area drift
and 0.2 K temperature drift before reporting settled averages.

The physical inputs for each target remain the experimental load, reference
viscosity convention, inlet/bath temperature, supply pressures and geometry,
thermal boundaries, and bearing/support properties. The image alone does not
supply them. With the established geometry R=0.0215 m, L=0.040 m and radial
clearance c=42.5e-6 m, the existing convention gives

    S = mu_ref * (rpm / 60) / (W / (2 * R * L)) * (R / c)**2
    W [N] = 7.336286043829293 * mu_ref [Pa s] * rpm / S

Match measured load and thermal conditions before plotting predictions against
the experiment. Preserve the separate experimental series and their histories.

## Model and pre-production checks

The G6/20-us synthetic startup failed the minimum-film guard before step 3
on both 24 and 28 MPI ranks. Before production, establish compatible initial
bearing/flow states and calibrated bearing mass/support/load data, and resolve
the lagged-force motion update at the requested timestep. Re-run a short
admissibility check before further MPI timing. Keep 720 by 192 and 20 us as
the selected protocol; no smaller timestep has been substituted.

Use TRANSIENT, MOVING_BEARING, FBNS, STRIATED, SINGLE_PHASE with zero gas charge,
KINEMATIC_TABLE at the pure-oil concentration, ENERGY_EQUATION, POLYNOMIAL
cross-film reconstruction and CONDUCTION_2D bush conduction. The journal centre
is fixed while the journal rotates; MOVING_BEARING describes bearing-centre
translation. True bearing rotation with a nonrotating journal is not that model.

The current CHT closure solves steady radial/circumferential bush conduction and
one global shaft temperature at each fluid step. It has no solid heat capacity
or bush axial conduction. Report settled-temperature comparisons under this
closure, not solid warm-up predictions.

Before quantitative friction comparison, make wall forces and torque consistent
with the cross-film F0/F1 integrals. The present torque still uses the mean oil
viscosity. If a speed ramp is used, pass the time-dependent speed to the cross-film
update; the present call uses the final configured speed. The timing benchmark
uses zero ramp time and does not apply either repair. Time-split transient
pressure, temperature and force evaluation also requires timestep refinement.

Record torque and its friction-coefficient convention, wall and film temperature
maps and sensor coordinates, liquid-filling maps, cavity area/fronts, bearing
position and minimum gap. A liquid-filling threshold is a model-based cavity
boundary, not a resolved optical bubble edge.

## Benchmark purpose and controls

The synthetic benchmark uses the selected G6 mesh, FW68D pure-oil viscosity,
cavitation, film energy, reduced CHT and bearing motion. Its imposed temperatures,
cavitation floor, support mass and load are provisional performance-test inputs.
A short benchmark interval does not establish thermal or mechanical settling.

Build and check locally, transfer a source/configuration snapshot to a new
workstation directory, hash-check it, and build there. Preserve all previous
runs and failed attempts. Do not reuse a historical executable for this test.

Measure MPI ranks 1, 2, 4, 7, 14, 24 and 28 without SMT. Set threaded numerical
libraries to one thread per MPI rank and verify physical-core binding. Compare
short/long runs with otherwise identical inputs to estimate time per additional
step using paired total runtimes. This reduces fixed launch-cost effects but
still measures the startup problem, not a long-time settled workload. Repeat leading candidates and compare
actual concurrent throughput on disjoint core sets. Report wall time,
core-seconds, numerical agreement and measured throughput separately.

Archive the source/configuration/executable hashes, machine and library versions,
MPI ranks and bindings, command lines, return codes, complete logs, diagnostics
and final field summaries. Synchronize evidence back locally before reporting
an optimum. A measured optimum applies to this mesh/model/workload only.

## Reproduction

The current G6 template reproduces the failed synthetic startup. The driver
stops on its first invalid run; the commands below are not a passing G6
performance recipe until that startup is resolved.

The timing driver targets Open MPI on the Linux workstation. Its pure-Python
self-check also runs on Windows; native Windows solver checks are recorded
separately. It uses only the Python standard library.

From the source root after building the solver:

```sh
python3 cases/90_diagnostics/fw68d_mpi_benchmark/run_benchmark.py --self-check
python3 cases/90_diagnostics/fw68d_mpi_benchmark/run_benchmark.py \
  --solver build-linux/pancake --source-dir . \
  --template cases/90_diagnostics/fw68d_mpi_benchmark/config_template.txt \
  --work results/fw68d_mpi_scaling --source-sha "$(git rev-parse HEAD)" \
  --ranks 1,2,4,7,14,24,28 --short-steps 2 --long-steps 6
```

Use a new work directory for every invocation; the driver refuses overwrites.
For throughput, select `--mode throughput --groups 7,7,7,7` (or the measured
candidate layout), still within 28 physical cores. CPU IDs and binding must be
verified for this host before reusing the script on another machine.

The two/six-step subtraction estimates the incremental cost of steps 3-6 of this
startup problem. It is not evidence of long-time steady per-step cost. Recheck
rank choice after changing grid, solver, physics or nonlinear workload.
