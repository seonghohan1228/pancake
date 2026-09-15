# Solver audit and repairs, 2026-09-15

## Scope and decisions

Reviewed the bearing solver's property evaluation, coupled conservation and
energy, inlet boundaries, linear solves, MPI exchange and existing regression
tests. The user requested discussion before removing any model. All model
selectors remain available. No new solver backend or dependency was introduced.
This is a bounded repair pass, not certification of every physical model.
The production `src/` changes remove 79 net lines; added regression checks and
this record are separate from that count.

## Repairs

| Issue | Change and affected paths | Verification |
|---|---|---|
| Repeated R410A/property work after a converged solve | Reuse the cell states from the final residual for acceptance checks, field publication and face-flux output. Refill the buffer on every residual evaluation, including perturbed Newton states. | Coupled conservation, energy and MPI tests; unchanged benchmark output. |
| Duplicate EOS work during root bracketing | Obtain pressure and its derivative from one Helmholtz evaluation. Keep the same bracket, Newton iteration, tolerance and phase limits. | Existing independent CoolProp density/caloric references and coupled R410A tests. |
| Inlet angles outside one revolution gave incorrect masks | Replace six geometry implementations with `Mesh::inside_inlet` and `Mesh::inlet_at`. Use periodic angular distance and retain last-matching-pressure-reservoir precedence. | Groove/hole masks and pressure are invariant under -720 and +1440 degree shifts; overlap and existing thermal-boundary tests. |
| Pressure and temperature pins used a diagonal penalty of 1e20 | Reuse the gas-transport method: replace the entire row, retain its diagonal scale, and leave neighbouring rows coupled to the prescribed cell. | Analytic axial diffusion at coefficient scales 1e-12 and 1e12; inlet, Elrod, thermal and mass-balance tests. |
| PETSc stopping controls were overwritten at each solve | Honor runtime SNES tolerances/limits and KSP relative tolerance/iteration limits, including fallback prefixes. Record iterations used by failed coupled attempts. | Limits of zero and one iteration stop coupled solves without overwriting accepted fields; a runtime KSP tolerance overrides a deliberately loose caller tolerance. |
| Too many MPI ranks could forward ghost values as physical cells | Reject partitions with fewer physical theta cells than the exchanged ghost depth, on every rank before communication. | An undersized partition is rejected; existing multi-rank tests still pass. |
| Linear-density supply enthalpy used the receiving temperature for flow work | Evaluate incoming liquid specific volume at the supply temperature. Outflow continues to use the cell state. | Independently checked cold-feed energy balance during uniform liquid cooling. |
| PETSc initialization failure did not stop the time loop | Return the initialization error before reading config or producing results. | An invalid PETSc boolean must fail before `Starting Simulation`. |
| Constant-viscosity warning incorrectly ruled out all thermal feedback | Say that viscosity feedback is disabled; density, gas compressibility and phase equilibrium can still depend on temperature. | Configuration checks and source inspection. |

The cold-feed repair makes the implemented enthalpy source consistent with its
selected density law. It does not supply a complete thermal-expansion caloric
model. That distinction remains a model decision below.

## Workstation measurements

Baseline source: `8f5cfad796b62f998dcac4796ed7e989c62a8309`.
Workstation: afml-r6, GCC 11.4, Open MPI 4.1.2, PETSc 3.15.5, Release build,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`.

Both cases use 96 x 24 cells, 20 steps of 0.1 ms, omega=1000 rad/s, fixed
eccentricity, 1 MPa absolute end pressures and identical output cadence. The
FW68D case uses COUPLED, SINGLE_PHASE, the actual-concentration kinematic table,
w=0.1 and 40 degrees C. The oil case uses FBNS, temperature-dependent empirical
viscosity and the energy equation. Neither case has an internal supply patch.
The isolated cold-feed test covers the changed internal enthalpy source.

Timing uses three repeats per case and rank count. Whole-process wall times
include MPI startup and output; PETSc elapsed time excludes most launcher cost.
The first smaller trial was retained but excluded because startup dominated it.

| Case | MPI ranks | Baseline wall median | Repaired wall median | Baseline / repaired PETSc elapsed median |
|---|---:|---:|---:|---:|
| FW68D | 1 | 1.669 s | 1.168 s | 1.225 / 0.712 s |
| FW68D | 4 | 0.967 s | 0.817 s | 0.470 / 0.326 s |
| Oil thermal | 1 | 1.268 s | 1.268 s | 0.808 / 0.811 s |
| Oil thermal | 4 | 0.767 s | 0.767 s | 0.299 / 0.310 s |

The final FW68D run takes 30% less wall time on one rank and 16% less on four.
An earlier candidate repeat gave a 26% four-rank reduction, so this short-case
parallel timing has noticeable variation. The final table uses only the final
executable. Oil thermal timings show no meaningful improvement.

| Verification | Native Windows | Remote afml-r6 |
|---|---|---|
| Release build | GCC 16.1, MS-MPI, PETSc 3.25.1; solver and GUI | GCC 11.4, Open MPI 4.1.2, PETSc 3.15.5 |
| Full CTest | 65/65, 48.42 s | 65/65, 80.01 s |
| Final benchmark cases | Four successful runs, 1/4 ranks | Twelve successful runs, three repeats at 1/4 ranks |

The 65 checks include 52 bearing tests and 13 optional bubble tests. The first
new wrapped-inlet pressure assertion initially compared different density
iterates; resetting the initial density made the comparison equivalent. The
failure is preserved, and the corrected regression passes. Cold-feed and startup
failure checks are included in the final full suites.
The final nonfinite-angle assertion was then rerun through `test_inlets` on
both machines and passed. All 157 transferred source/build/data fixture files
match byte for byte between the tested local and workstation snapshots.

Saved pressure, liquid density, dynamic viscosity and temperature are compared
at the VTK output precision. This is numerical regression evidence, not
experimental FW68D validation or a mesh/time convergence study.
Windows and Linux also agree: the largest saved pressure difference is
3.5e-9 Pa and the largest temperature difference is 5.7e-13 K in these cases.

The baseline one-rank FW68D profile spent about 0.34 s building two numerical
Jacobians versus 0.044 s in linear solves. There were 64 colored residual calls
for those Jacobians. These timings overlap with the parent SNES events and must
not be added to them. Reusing accepted states removes substantial work outside
the nonlinear solve as well. PETSc's [profiling guide](https://petsc.org/release/manual/profiling/)
explains the event hierarchy; [SNES options](https://petsc.org/release/manualpages/SNES/SNESSetFromOptions/)
document the runtime stopping controls.

## Remaining physics and feature decisions

These models have different limits. A converged residual cannot establish their
experimental accuracy. No removal or new restriction was made in this pass.

| Priority | Item | Assessment and proposed next step |
|---|---|---|
| High | Legacy gas transport and `steady_gas_model=EQUILIBRIUM` | The legacy local feed flash does not solve conservative steady species transport. JFO/FBNS gas updates do not impose the same simultaneous EOS-volume closure as COUPLED. Decide whether to keep these as explicitly approximate comparisons or restrict production gas cases to COUPLED. |
| High | LINEAR_TEMPERATURE with coupled energy | Density varies with T while liquid internal energies retain the incompressible caloric approximation. Decide between documenting that approximation, restricting the combination, or deriving and validating a consistent thermal-expansion model with identified caloric inputs. |
| High | Property-table endpoints | Endpoint clamping extends data beyond their measured range. In particular, the FW68D solubility curves start at w=0.1; the low-pressure clamped value is not a measured dilute Henry limit. Choose strict bounds for a data-domain study or explicitly justify the extension. Neither constant nor linear liquid density supplies a measured pure-R410A liquid density. |
| Medium | GUMBEL and FULL_SOMMERFELD | Gumbel clips pressure after solving and does not conserve cavitation mass. Full Sommerfeld retains tensile pressure. They remain useful limiting comparisons, but should not be confused with mass-conserving rupture/reformation models. Discuss retention before removing either. |
| Medium | First-order phase/thermal transport | COUPLED uses upwind spatial transport and backward Euler in time. These dissipate resolved gradients. Assess load, cavity extent and temperature under mesh/time refinement before selecting a higher-order conservative, bounded scheme. A different bearing-motion integrator does not raise the order of these field equations. |
| Medium | Reduced thermal geometry | COUPLED has one film temperature; the legacy cross-film/solid closure resolves different physics. The reduced bush/shaft model is not full three-dimensional conjugate heat transfer. Match the chosen closure to the measured temperatures before comparing experiments. |
| Low | Experimental bubble target | It remains optional and unfinished. No bearing solver simplification requires deleting it. |

## Next performance work

1. Profile an agreed production workload, including released refrigerant and
   thermal transients, at a fixed accuracy requirement. The present timings do
   not measure all those regimes.
2. For COUPLED, the colored finite-difference Jacobian still repeats property and
   flux evaluations. An analytic/local derivative implementation could reduce
   this cost, but adds maintenance and phase-boundary derivative work. Do not
   introduce it before measuring a production bottleneck and checking derivatives.
3. For pure-oil thermal runs, this pass did not show a runtime improvement.
   Profile field-map lookup, repeated scratch-field allocation, ghost exchange
   and thermal/cross-film updates before changing their layout or communication.
4. Retain PETSc and the existing solver interface. A backend migration is not
   justified by these profiles. Do not change overlap, tolerance or MPI count
   on the assumption that more work or more ranks must be faster.

## Reproducibility

Evidence is preserved locally at `C:/Dev/pancake-audit-20260915` and remotely at
`/home/user/pancake-runs/audit-20260915`. The old workstation checkout and its
baseline executable were preserved. The directory contains source/fixture
archives, build and test logs, benchmark configs, executable/input hashes,
PETSc profiles, output comparisons and the small stdlib Python benchmark scripts.
`source-manifest.json` and `verification.json` record the source and executable
identities. The tested source aggregate SHA-256 is
`1a8f0b5c3ccec6d7e2f87835cdbc47bda7229981f4a7b0e5b1a6df2a2550d883`.
The initial Windows profiling command exposed a PETSc viewer-syntax issue with
drive-letter colons. It is retained as a failed run. The portable benchmark uses
`:petsc.log` in each run's own working directory.
