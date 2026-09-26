# FW68D MPI test with the prior convergence settings — 2026-09-27

The campaign and benchmark now use the previously selected **G6 720×192 film
mesh and dt = 2e-5 s (20 us)**. The current synthetic moving-bearing startup
fails the positive-film guard before step 3 at these settings. **No accepted
multicore optimum is available for this G6 case.** The earlier coarse-grid
rank recommendation must not be treated as a G6 result.

## Settings carried forward

The [frozen campaign protocol](../../VALIDATION_PLAN.md) and
[completed Ferron D-series convergence tests](../../01e_ferron_terra_repro/README.md)
select G6 and 20 us. The D-series checked 20/10/5 us through 0.10 s. Its film
mesh ladder reached 1620×432 while retaining G6 as the operational resolution.
Retain its CHT implementation settings: 16 radial bush cells, polynomial order
4 and six quadrature points. Their separate radial/profile independence was
not established by that film-plane mesh ladder.

The present FW68D plan uses G7 1080×288 and 10 us for the four endpoint checks,
keeping the existing 24 base cases plus eight verification repeats. It does
not restart the entire old convergence ladder. See [the simulation plan](README.md).

## What the new test found

The numerical C++ source was unchanged. A native Windows two-rank/two-step
preflight passed, followed by a fresh Release build in an isolated workstation
directory. The 24-rank runs used explicit physical-core rankfiles, one numerical
library thread per rank and the same output cadence.

| G6 test | Requested steps | Result | Wall time |
|---|---:|---|---:|
| Native Windows, 2 ranks | 2 | Exit 0; both steps converged; positive gap | Preflight only |
| afml-r6, 24 ranks | 2 | Exit 0; required fields valid | 14.102 s |
| afml-r6, 24 ranks | 6 | Exit 4; minimum-film guard before step 3 | 13.903 s |
| afml-r6, 28 ranks | 6 | Exit 4; same guard and trajectory | 15.196 s |

The two completed workstation steps required 57 and 73 FBNS Newton iterations.
Both converged; absolute liquid residuals were below 5e-16 at 24 ranks and
1e-15 at 28 ranks. At step 2
(t=40 us), minimum film thickness was 37.7364 um, maximum film temperature
358.246 K, and fluid force approximately (-40.42, -439.85) kN. The following
motion update produced a logged minimum film thickness of approximately
**-60 um**, triggering the configured stop before step 3.

This is an inadmissible moving-bearing trajectory, not an FBNS convergence
failure in the completed steps. Negative model gap is not evidence of contact
in the real apparatus. The fixture's 1.64 kg bearing mass, 1 kN load, initial
eccentricity and thermal conditions are synthetic. The prior Ferron protocol
used a different geometry and dynamic setup; its accepted timestep does not
by itself establish stability of this FW68D startup.

The code path is explicit about the ordering: [main.cpp](../../../src/main.cpp)
advances bearing motion with the previous step's fluid force, then checks
film thickness before the new flow solve. With zero stiffness/damping,
[journal_motion.cpp](../../../src/journal_motion.cpp) gives a force contribution
to displacement of dt²(F_y + W_y)/m. The step-2 force contributes about -107.5 um
in y at 20 us, exceeding the 42.5-um clearance. Both rank counts reproduce
the guard, so changing MPI decomposition does not resolve this trajectory.

A failed six-step run cannot be paired with a successful two-step run to
estimate seconds per step. The shorter failed wall time is therefore not a
performance improvement. Remaining scaling/throughput trials were stopped;
the timestep and physical inputs were not adjusted to obtain a passing timing.

Before selecting production ranks, establish admissible initialization and
bearing dynamics for the calibrated physical case at the requested settings,
then repeat the short MPI comparison. Keep the selected 720×192/20-us protocol
in the plan while that issue is investigated.

## Earlier coarse-grid probe

The initial 360×96, 1-us probe completed 64 accepted workstation executions.
It found lowest isolated core cost at one rank, highest tested queue throughput
with 14 concurrent two-rank jobs, and lowest single-case latency at 28 ranks.
The measured concurrent batches were 24 one-rank cases in 73.634 s and 14
two-rank cases in 39.029 s. These results belong to that smaller, shorter-step
fixture; they do not establish the current G6 rank choice.

Its [full initial report](../../../results/fw68d_mpi_benchmark_20260927/initial_probe_report.md)
and numerical audit remain archived. The complete coarse-grid source/build/run
archive is `remote/evidence-complete.tar.gz` below the local evidence root;
SHA-256 `828b4b2651d0b42cbb88774bd3830c64ca93fe04f4426486fa83ef76400d28ee`.
The archive includes all raw frames and the earlier MPI placement failure,
which is excluded from accepted timings.

## Provenance

Numerical source commit: `ea337f9024ba8032a280d8e0c26d66b9615259a3`.

| Artifact | SHA-256 |
|---|---|
| Source archive | `be3126386171e80447a68f1183f33c14ec87e9882627449c052078deb654f679` |
| G6 configuration template | `db43b8ac7ec089374baa3cb07cb3475b44230d198ae7ed72ac31349cba369d08` |
| Timing driver | `f7a08ac332c4d8f9e3b324a296042a3c2a9008b5693268608bd120cc6bb4c959` |
| FW68D table | `90057c6b2497506b1b6d4c371f6b938764910aa6cba305919eb5447fc43ed45c` |
| Workstation executable | `ad79ee366575e036b5509f13a203acbab5fdc34b7a64b263976c34eff8bc2ddb` |
| G6 complete evidence archive | `ebe93be96efe5e3781acb786ce3462b94fa92e11d2fa36991e1c71fc337ae84e` |
| Native executable | `76decbcaaa96a17949bf4729df35d1ea9fb27ff7c416b5a372d4fc703d70dd28` |

Workstation: afml-r6, two Xeon E5-2660 v4 CPUs, 28 physical cores. Open MPI
rankfiles pin one rank per distinct physical core; SMT is not used to add ranks.
The fresh Linux Release build uses GCC 11.4, CMake 3.22.1, Open MPI 4.1.2 and
PETSc 3.15.5. The native checks use GCC 16.1, PETSc 3.25.1 and MS-MPI
10.1.12498.52. The timing-driver self-check passes on both platforms.

Remote G6 evidence:
`/home/user/pancake-runs/fw68d-cht-g6-20260927-ea337f9-bench02`.
The previous probe remains untouched at
`/home/user/pancake-runs/fw68d-cht-20260927-ea337f9-bench01`.

Local evidence root:
`C:/Dev/pancake/results/fw68d_mpi_benchmark_20260927/`.
`local`, `remote` and `analysis` contain the initial probe;
`g6_local` contains the native G6 check; `g6_remote` contains the verified
complete G6 archive and extracted logs, diagnostics and all available raw fields.
`analysis/analysis_g6.json` records the independent G6 audit. Workstation
manifests preserve exact
input/executable hashes, commands, binding, return codes, logs, diagnostic rows
and available field summaries, including failed runs. Generated data stays
outside Git.
