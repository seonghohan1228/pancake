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


## September 28 repair and validation

The original A/B/C first-position thermal solves reached their 2000-sweep
caps. E later failed during its position search; D/F and queued variants were
stopped after diagnosis. Only verified processes belonging to the obsolete
campaign were stopped, and all files plus a process/stop record were retained.
No original case is relabeled as a converged equilibrium.

The original C position at epsilon=0.6, attitude 140 degrees converged after
reducing quartic-source relaxation to 0.3: 513 sweeps at both 8 and 4 MPI ranks,
531.8 and 1000.5 seconds respectively. Its fluid force was only about 469 N
against 1973 N applied load, and the PINNED thermal defect remained. This is
numerical damping evidence, not an accepted position or temperature result.

The repair uses FLUX_BALANCED groove energy, FBNS-consistent mass/enthalpy
faces, constant mass-specific cp, exact profile extrema and raw collocation
residual checks, a global energy ledger, and a timestep-independent mass-rate
check. Saved steady fields accelerate subsequent same-case position/h trials.
The native Windows full build succeeds. All 52 native CTest checks pass after
four FBNS variants were rebuilt/rerun with a tighter test solve tolerance:
1e-13 instead of 1e-11, retaining the original 1e-9 physical-flux assertion.
The changed energy model exposed that mismatch between the test's inner
Jacobian-scaled norm and its physical-flow residual requirement.

A 720x192, four-rank, three-sweep startup check returned the expected cap
failure 4 and wrote the new energy/profile columns. Its executable SHA256 is
77a75f8f0cd0b0ef14385ab3840b2dfc47ce610fb8ce99e5e647239197ae8cdb.
This check is not convergence evidence. Native validation is retained at
`results/fw68d_physical_v3_short_native/` in the managed worktree.

The corrected remote source is isolated under
`/home/user/pancake-runs/fw68d-physical-20260928-a2ba963-c6d9373b/`.
Source archive SHA256:
`c6d9373b8e16f82dd87f1e59c1680e66b444fb06e5d6a504330439dfedcc9f30`.
The subsequent test-only tolerance patch is recorded separately; it does
not change production source or executable physics. Remote trials and
accepted equilibria will be appended when their checks complete.


The fresh v3 workstation Release build and 14 focused CTest checks passed,
as did explicit four-rank I/O and diagnostics checks and the case-driver
self-check. Production executable SHA256 is
`d9d19941151107d308ca484a1a6a5c212366c041553e83cbb67a032f8c55396a`.
Compiler/MPI/PETSc are GCC 11.4 / OpenMPI 4.1.2 / PETSc 3.15.5.


### Corrected E inner trial and force-search check

At 2026-09-27 15:30:31 UTC (September 28 KST), E's first fixed-position
trial completed in 616.160 s on eight ranks after 641 coupling sweeps.
R_p=1.59038e-7, R_T=7.91947e-10 and raw R_X=4.50689e-9; the relative
energy defect is 9.91995e-7 and liquid-rate defect is 4.62e-16. The exact
profile range is 328.658--333.54 K, above TC0=324.71526 K. The minimum gap
is 36.13 micrometres; pressure is 0--229.8 kPa. These are inner-solve checks.

The force is (-52.333,191.208) N against the required (0,461.252) N,
so its vector load error remains 59.635%. Its initial external-coefficient
mismatch is 15.657%. The record is `trial_only_not_thermally_closed`, not an
accepted bearing equilibrium. The next cooling evaluation reduced the h
mismatch to 1.106% while again passing the inner energy/physical gates.

A full E cooling/position search was launched with the existing driver at
15:31:38 UTC, PID 4057589, reusing this verified checkpoint. It uses physical
cores 16--23 and the same `diagnostic_trials/E` work directory. A/C remain
in their first trials. B has not launched: automatic approval review blocked
its launcher upload; the destination-specific user approval is pending.

The standalone check `python tests/test_fw68d_force_balance.py` passes. It
uses a known two-component force equilibrium from two different positions,
rejects the first full candidate, and verifies retained history plus recovery.
It tests the search algorithm, not the CFD force law or uniqueness of the
physical equilibrium. The current motion model assumes fixed journal centre
with a rotating journal surface; confirmation of that interpretation of
"journal stationary" is pending.


### Thermal closure and iteration-cap follow-up

E reached external thermal closure on its fourth evaluation: relative h
mismatch 6.20e-5, wall-temperature change 0.00134 K, and inner energy
defect 9.97e-7. Its force error was still about 59.66%, so the position
search continued. For the downloaded E evaluations, shaft torque times
angular speed and integrated viscous heat agree within 0.003%; this is a
consistency check, not an independently complete mechanical energy ledger.

A and C exhausted 2000 sweeps at 15:50:38/43 UTC, returning solver code 4.
A had R_p=1.81038e-6 and energy residual -0.000808094 W on a 9.11394 W
scale (relative 8.87e-5). C's retained best iterate had energy residual
-2.63664e-5 W on a 6.49650 W scale (relative 4.06e-6); it is explicitly
marked restored and unconverged. Their exact profile minima, 303.398 K
and 302.558 K, remained above their respective TC0 values. These failures
are insufficient iteration budgets for the slowly decaying coupled thermal
mode; neither is accepted.

The existing driver/executable were restarted for full A/C cooling and
position searches with `--coupling-max-iters 4000 --timeout 7200`, leaving
all mesh, physical and tolerance settings unchanged. A driver PID 4068349
uses cores 0--7 and log `A_force4000_driver.log`; C PID 4068483 uses cores
8--15 and `C_force4000_driver.log`. Existing work directories preserve the
failed attempts; no failed state is reused. These searches and the E search
remain in progress. No final experimental comparison or verified merge has
been claimed.


### Confirmed motion and execution scope

The user confirmed on September 28 KST that the journal rotates at the
previously supplied RPM while its centre remains fixed. A/C/D/F use 200 rpm;
B/E use 1500 rpm, including E's TC0 variants. The bearing centre is an
equilibrium unknown: the outer search updates eccentricity and attitude
until fluid force balances the applied load. The bearing housing has no
imposed rotation. This confirms the motion convention already used in the
current runs; no physical configuration change or RPM substitution is needed.
The user also answered yes to the explicit B workstation-transfer/run
question. The prepared B launcher is being retried through normal approval.

Automatic approval review rejected the B launcher upload a second time,
despite the preceding yes response. It requires explicit authorization of
the named local launcher transfer to user@143.248.99.6. A full-sentence
approval question is pending; no B upload or launch succeeded. Existing
A/C/E computations continue.


## September 28 11:34 KST: first accepted equilibrium

Case E completed the full thermal and position search. The final result is
`diagnostic_trials/E/E/result_5bf30b10.json`, supported by inner evaluation
`0050_force6_step1_thermal2_e9e1d1ec`. The executable and production-source
hashes are unchanged from v3. Epsilon=0.5915322966040504 and attitude
131.08340747690076 degrees; fluid force is (-0.0000174,461.652433) N
against the target (0,461.252102) N. Relative vector load error is
0.0008679225 (0.0867923%), below 0.2%.

The final inner solve returned 0, with coupling converged and no restored
iterate. Absolute relative mass defect is 1.27215e-9, energy defect
9.68616e-7, and raw profile defect 1.74979e-8. External-coefficient mismatch
is 1.66634e-5 and wall-temperature change 0.000297556 K. All configured
acceptance gates pass. Mean-film temperatures span 53.07397--58.93606 C;
the exact reconstructed profile spans 51.87065--59.32023 C, above TC0.
Minimum gap is 17.35998 micrometres, maximum pressure 882.0366 kPa,
shaft torque 0.29517068 Nm and cavitated area fraction 0.5067274 using
film content below 0.999. This establishes numerical/model admissibility;
experimental agreement is still to be assessed.

The complete E final fields, result and force history, config, logs, source
snapshot and build/test provenance were returned under
`C:/Dev/pancake/results/fw68d_physical_20260928/final/E/`.
At the same status check A/C searches were active with best thermally
closed vector-load errors 35.5025% and 8.83798%; they were not accepted.
B remained unlaunched pending the explicit transfer authorization.

At 11:34:56 KST the already-deployed driver started D, followed sequentially
by F, E_low and E_high, on E's freed cores 16--23 (eight ranks). PID 151952
uses work directory `remaining_cases`, log `remaining_cases_driver.log`,
4000 coupling iterations maximum and 7200 s per inner solve. D's first
config and MPI startup were verified. No new source was uploaded for this
launch. Later cases remain queued; an earlier failure stops this driver and
requires inspection before the remaining queue can continue.


### Case B launched without a new source transfer

Read-only inspection confirmed that the workstation already held Case B's
experimental inputs, the v3 driver/executable, and its cores 24--27 rankfile.
A separate normal approval review permitted direct invocation of those
existing files. The previously rejected launcher upload was unnecessary and
was not retried; no local source/input file was transferred for this launch.

B started at 2026-09-28 11:59:29 KST, driver PID 162977, work directory
`diagnostic_trials/B`, log `B_force4000_driver.log`. All four MPI ranks were
verified on physical cores 24--27. Its first generated config uses omega
157.07963267948966 rad/s (1500 rpm), applied load 3237.813830374732 N,
TC0=334.63691239412935 K and coupling cap 4000; geometry, mesh, physical
models and acceptance tolerances remain v3. First config SHA256 is
94c12bba983db41fc463aac992b63ec05380cc70639826d067e410b36987ca7c.
B is running a full thermal/position search, and its former upload-approval
question is no longer a prerequisite for this execution.


## September 28: prior campaign invalidated by confirmed upward push

The user confirmed the actuator below the bearing pushes it upward, with the
groove at the top. Previous runs explicitly used external_load_y=-W and the
outer search targeted Fy=+W. Every old case, including the numerically closed
E above, is invalid for the confirmed experiment. Its 66.02% friction mismatch
and near-groove temperature interpretation are results of the wrong physical
configuration. They are preserved for audit and are not final comparisons.

Corrected requirements are external_load_y=+W, fluid-force target Fy=-W,
residual (Fx,Fy+W), and groove theta=90 degrees. The new campaign requires a
fresh work directory and verified inputs; old equilibria are not reusable.
Stopping, tests, hashes and replacement launches will be recorded below.

### Corrected-direction local validation

The driver now shares one applied-load vector (0,+W) across configuration,
force-error checks and the Newton residual. Its startup angles are A300,
B/C320, D335 and E/F/variants345 degrees; these are only initial guesses.
Angle is unrestricted and eccentricity remains bounded. Saved config hashes
and the explicit upward-load components must pass before a cached result,
position or checkpoint can be reused. Python self-check and the two-start
force-search/old-cache rejection regression pass, as does git diff --check.

Native Release rebuild and the existing test_forces/test_journal_motion pass.
A four-rank C diagnostic at eps0.6, psi300 used the full720x192 grid, exact
TC0 and pressure/geometry/viscosity inputs, with energy disabled solely to
isolate the force sign. Return code0 in7.110s; fluid force is
(-657.835884,-1334.774681)N, minimum gap17.00024um at300.25degrees, and
the bearing centre lies above the journal. Pressure/filling/gap checks pass.
This does not establish a load or thermal equilibrium. Evidence is under
C:/Dev/pancake/results/fw68d_physical_20260928/load_orientation_correction/.

Native executable SHA256:
ffdf96db9829b075b480c7d9ad99e0ba4e9ae7c221fc82530abe45efe72a3640.
Production C++ aggregate hash remains
fc32d207b6bad1bb8e7736a16a95629d38c5525fc4328d428058eef8293e9163.
Driver SHA256 after correction:
eb8bdae2efbb2b905de36da10fe7620dcdd18bdeb9c9e9fe75d652284e8edeb1.

The attempted targeted stop of the old workstation campaign was rejected by
automatic approval review, which required explicit authorization to interrupt
the jobs. No termination signals were sent. The exact host/campaign stop
question is pending; no corrected remote campaign has been launched yet.


### Authorized stop and additional pressure-convergence fix

The user's explicit stop instruction superseded the earlier approval blocker.
At 2026-09-28 04:36 UTC all 36 old campaign processes were stopped using targeted
SIGTERM; no solver processes remained and all files were preserved. Evidence:
load_orientation_correction/termination_authorized.json and resource_readiness.json.

The final audit found a second convergence issue: an initial FBNS residual exit
could skip the actual pressure update and report R_p=0. Steady calls now require
at least one Newton correction; transient behavior and tolerances are unchanged.
Exact zero-RHS and sub-tolerance perturbed warm-state tests pass at 1/2/4/6 native
MPI ranks. The prior a46dbb7 remote build passed its 14 focused tests, four-rank
I/O/diagnostics and Python checks, but is staging evidence only: it must not
launch production because it predates this additional correction.

Replacement scheduling uses seven independent four-rank workers on the 28
physical cores. A-F and E_low start first; E_high uses the first released slot.
Six slots are socket-local and one spans the remaining two cores per socket.
Explicit rankfiles and one numerical thread per rank avoid SMT oversubscription.
Prior timing evidence supports this resource-efficient layout but does not prove
it is the absolute fastest makespan for the corrected coupled cases. A case
failure releases its slot and does not cancel the other cases or the eighth run.
All production acceptance gates and the original converged resolution remain.

Final native Release build and all 52 CTests pass (46.45 s). Driver self-check
and force-balance/cache-rejection regression also pass. Native executable
SHA256: 4a9f462ed7ddd3504f94deb2a8544449b58cf9f00db16187660170552d05cbb7.
The executable includes the minimum steady FBNS correction. Full evidence is
load_orientation_correction/native_final_validation.json and native_final_ctest.log.
