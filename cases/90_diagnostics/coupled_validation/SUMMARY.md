# Workstation results and their interpretation

The final-source campaign runs on `afml-r6`, a dual Xeon E5-2660 v4 machine
with 28 physical cores, using GCC, OpenMPI 4.1.2 and PETSc 3.15.5. The
transferred source was checked by SHA-256 before building. All **51/51**
CTest checks pass there, as on native Windows and local WSL. The original
remote and local campaign checkouts are preserved.

## Paired oil-air timestep test

This is a controlled synthetic case, not an FW68D prediction. Both selections
run from the same executable, geometry and material inputs on eight ranks.
The segregated selection is the repaired legacy FBNS/TWO_PHASE combination;
it is not the gas-free FBNS model that passed the earlier JFO verification.

| Timestep (µs) | Segregated load at 1 ms (N) | Coupled load at 1 ms (N) | Segregated maximum excess gap volume | Coupled maximum volume defect |
|---:|---:|---:|---:|---:|
| 10 | 1258.067 | 598.653598 | 6.681 | 2.22e-11 |
| 5 | 3259.939 | 598.648229 | 2.856 | 1.46e-11 |
| 2.5 | 1401.515 | 598.645544 | 1.548 | 1.87e-12 |

The defect is normalized by the geometrical gap. A value of one means one
whole extra gap volume, so the old result is physically inadmissible even
if its pressure solve converged. The coupled terminal-load difference
between 5 and 2.5 µs is 0.002685 N, about 0.00045%. Its maximum normalized
species step residual is below 2.4e-17 in these runs.

![Load histories](figures/air_load_timestep.png)

The two panels deliberately use different load units (kN and N). The very
large early segregated spike is retained. Both omit only the unadvanced
initial frame. The coupled curves are nearly coincident without temporal
averaging. This establishes timestep consistency for this case; it does not
calibrate the synthetic solubility or phase-transfer constants.

![Volume consistency](figures/air_volume_closure.png)

The original pre-audit gas-inventory loss and unchanged pure-oil steady
fields are documented in the separate
[historical before/after comparison](../physics_comparison/SUMMARY.md).
They should not be confused with this later pressure/volume repair.

## Thermal effect

These cases use synthetic viscosity, solubility, partial heat capacities and
desorption energy. They run for 20 ms at 600 rad/s on 48 × 8 cells, with
reservoirs initially at 40°C. The oil-air and oil-R410A charges and supply
pressures differ, so this is not an isolated experiment comparing species.

| Mixture | Isothermal load (N) | Thermal load (N) | Load change | Maximum film temperature (°C) |
|---|---:|---:|---:|---:|
| Oil–air | 1120.600 | 1040.764 | −7.12% | 45.804 |
| Oil–R410A | 695.342 | 683.260 | −1.74% | 43.655 |

![Thermal effect](figures/thermal_effect.png)

The differences show why an isothermal result cannot answer the user's
heating question. They do not predict the magnitude of heating in FW68D.
For the oil-air case, timestep refinement from 100 to 50 to 25 µs gives
1040.7644, 1040.7837 and 1040.7955 N. The corresponding maximum temperatures
are 45.8041, 45.8014 and 45.8001°C. Energy conservation is evaluated using
the actual coupled fluxes; the maximum global residual is below 6e-9 W in
the two 48 × 8 thermal comparison runs.

Spatial refinement is recorded separately in `refinement_evidence`.
A small timestep error does not remove the spatial error of the 48 × 8 grid.

| Mesh | Load at 20 ms (N) | Maximum temperature (°C) | Eight-rank elapsed time (s) |
|---|---:|---:|---:|
| 24 × 4 | 1119.2444 | 45.7144 | 1.42 |
| 48 × 8 | 1040.7644 | 45.8041 | 2.72 |
| 96 × 16 | 1014.8616 | 45.8526 | 8.14 |
| 192 × 32 | 1007.5049 | 45.9037 | 36.42 |
| 384 × 64 | 1005.1555 | 45.9615 | 235.77 |

The last mesh doubling changes load by 0.234% and maximum temperature by
0.058°C. Thus 192 × 32 is much closer to the finest result than the quick
48 × 8 case, whose load is 3.54% above it. This is an observed refinement
difference for this case, not a universal error bound or validation against
an experiment. The large cost increase on the finest grid also shows why
runtime should be compared at a specified spatial accuracy.

![Spatial and temporal refinement](figures/thermal_resolution.png)

## Independent JFO and gas-cavity comparisons

The pure-vapor benchmark compares with an independently integrated
one-dimensional JFO free-boundary solution. Exact rupture and reformation
are 302.766608° and 79.006725°. The 360- and 720-cell coupled pressure
errors are 0.3490% and 0.1737% in the stated relative L2 norm. At 1440 and
2880 cells these decrease to 0.0865% and 0.0432%. On the finest grid the
detected rupture/reformation are 302.75°/79.125°, both within one cell
(0.125°) of the independent values. Full results are in
`jfo_refinement_evidence`.

The 180-cell error is unusually small (0.0299%) because the cell-centered
groove mask changes its represented edge alignment on that mesh. It is an
accidental cancellation, not evidence that 180 cells are more accurate.
Convergence must be assessed across meshes representing the same groove
edges. The face-gap transport repair improves the finer-grid result compared
with the earlier fully upwind gap discretization.

![JFO solution and refinement](figures/jfo_verification.png)

The photographed gas cavity is **not reproduced**. At a liquid-deficit
threshold of 0.001, the equilibrium prediction is 165° on both 120 × 32
and 240 × 64 grids, versus a measured 235° extent. At threshold 0.01 the
predictions are 156° and 154.5°. These discrepancies are much larger than
the change under refinement or the endpoint-derived ±10° extent bound.

Zero absorption predicts gas over the entire circumference at thresholds
up to 0.01; raising the threshold to 0.05 gives 210°. That threshold
sensitivity does not establish agreement with a photograph. The model's
homogeneous phase assumption and unmeasured charge/kinetics remain possible
causes; the comparison does not uniquely identify one cause
(Etsion and Ludwig, *[J. Lubr. Technol.](https://doi.org/10.1115/1.3253174)*, 1982).

![Gas-cavity extent](figures/etsion_gas_comparison.png)

## Runtime and MPI check

A 96 × 16 thermal pilot covering 2 ms took 5.429, 2.070 and 1.468 seconds
on 1, 4 and 8 ranks. Its load was 1056.071207217 N in all three runs
to the shown precision. Eight ranks gave a 3.70× elapsed-time speedup over
one rank, including process startup and output. These are single timing
samples for this workload, not a scaling guarantee on larger meshes.

The 48 × 8, 20-ms thermal runs took 2.72 seconds for air and 8.94 seconds
for R410A on eight ranks. The 240 × 64 steady gas-cavity calculation took
105.1 seconds. Timings from the older local campaigns must not be substituted
for these workstation measurements.

## R410A phase-domain and numerical failures

With 20-bar ends, finite-rate 60°C and 80°C cases complete. The 40°C case
is rejected when resolved free refrigerant leaves its stable gas domain.
The initial 20°C run failed at the first step in both finite-rate and
equilibrium modes. Serial and direct-LU retries and timestep halving retained
that failure. Their diagnostics identified a zero-pivot matrix before the
first Newton iteration. The actual defect was a phase-domain check on a
tiny outgoing trial gas perturbation, even though the physical solution
was gas-free. This prevented forming the numerical Jacobian. The corrected
outgoing trial uses EOS continuation, while accepted physical free phases
retain strict checks. A new regression fails before that repair and passes
after it on all three tested platforms. `r410a_20C_diagnostics` retains
the original failures; corrected bearing cases use separate `*_trialfix_evidence`
manifests.

After this repair, all four 20-bar local-equilibrium runs complete:

| Initial/reservoir temperature (°C) | Maximum film temperature (°C) | Load (N) | Maximum film pressure (bar absolute) |
|---:|---:|---:|---:|
| 20 | 27.378 | 2095.227 | 33.629 |
| 40 | 44.908 | 1454.317 | 27.793 |
| 60 | 63.351 | 1126.379 | 24.433 |
| 80 | 82.298 | 973.714 | 22.472 |

All four use the same synthetic 10 wt% initial dissolved charge, not a
separately equilibrated 10-bar charge at every temperature. The internal
pressure patch is 10 bar and both ends are 20 bar; flow direction follows
the calculated field. Peak film pressure exceeds the reservoir pressure.
The finite-rate 20°C and 40°C reruns still retain free refrigerant outside
the accepted stable gas domain and are rejected. The correction therefore
removes a numerical failure without concealing the remaining phase-model limit.

Do not turn a last available time slice from any failed case into a final
prediction. The manifests preserve every failed return code. Instantaneous
dissolution changes the physical assumption and must be reported explicitly.
The fixed-blend EOS does not supply the missing FW68D solubility/caloric
data or a separate refrigerant-rich liquid phase.

## Provenance and reproduction

`final_evidence`, `final_thermal_evidence`, `final_equilibrium_evidence` and
`final_benchmark_evidence` contain configurations, logs, reduced histories,
diagnostics and source/binary hashes. `pilot_mpi*` records the rank comparison.
The R410A trial-state repair follows this primary snapshot; its two changed
reservoir guards only affect the R410A branch. The air/JFO equations and
recorded refinement results are unchanged. R410A cases are rerun against
the repaired binary and retained separately rather than overwriting failures.
The initial benchmark dispatch lacked a transferred config template and
stopped before running its cases; after transferring the template the full
benchmark batch completed. This packaging failure is retained in the remote
campaign log and is not a physics result.

Generated evidence and figures are separate from the implementation commit.
The [reproduction commands](README.md) regenerate them. The detailed
[implementation review](../../../docs/COUPLED_IMPLEMENTATION_REVIEW.md)
explains why each correction was needed and the remaining physical scope.
