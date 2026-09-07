# Coupled cavitation and thermal verification

This campaign distinguishes mathematical verification, material-property
checks and comparison with a gas-cavity experiment. The bearing constants
used in synthetic oil–air and oil–R410A runs are **not FW68D data**.

## Current execution state

The [workstation campaign](SUMMARY.md) ran on `afml-r6` in the isolated
directory `/home/user/pancake_coupled_20260907`. It includes paired solvers,
thermal and spatial refinement, R410A phase-domain tests and the Etsion gas
cavity. The current C++ source passes all 51 tests on native Windows, local
WSL and that remote Linux workstation. Manifests record hostname, operating
system, MPI rank count, source and binary hashes. The original worktrees
and existing remote campaigns are preserved.

Existing `evidence`, `thermal_evidence`, `equilibrium_evidence` and
`benchmark_evidence` folders preserve **intermediate** local implementations.
Their source hashes precede the final face-gap Couette transport and outgoing
reservoir-property repairs. Do not use them as final-source or remote timing
evidence. Workstation results use `final_*`, `refinement_evidence`,
`jfo_refinement_evidence`, `pilot_mpi*` and `*_trialfix_evidence`; the latter
rerun R410A after repairing the outgoing Newton-trial domain gate. Raw VTK
output remains under ignored build directories on the execution machine.

## Reproduction

Use Python with NumPy and Matplotlib, and a Release solver built from the
same source snapshot. No CoolProp library is required to run the solver.
From the repository root, a Linux workstation can run:

```bash
python3 cases/90_diagnostics/coupled_validation/run_validation.py \
  --binary build-linux/pancake --work build/coupled-comparison \
  --evidence cases/90_diagnostics/coupled_validation/final_evidence --ranks 8
python3 cases/90_diagnostics/coupled_validation/run_thermal.py \
  --binary build-linux/pancake --work build/coupled-thermal \
  --evidence cases/90_diagnostics/coupled_validation/final_thermal_evidence --ranks 8
python3 cases/90_diagnostics/coupled_validation/run_thermal.py \
  --binary build-linux/pancake --work build/coupled-equilibrium \
  --evidence cases/90_diagnostics/coupled_validation/final_equilibrium_evidence \
  --ranks 8 --only 20bar --r410a-equilibrium
python3 cases/90_diagnostics/coupled_validation/run_benchmarks.py \
  --binary build-linux/pancake --work build/coupled-benchmarks \
  --evidence cases/90_diagnostics/coupled_validation/final_benchmark_evidence --ranks 8
```

`--ranks 1` runs without an MPI launcher. `--mpi-launcher` selects a different
launcher. Run cases sequentially when measuring elapsed time. Eight ranks
were fastest among 1/4/8 in the recorded 96 × 16 pilot, not an established
optimum for every mesh. Use a separate evidence
folder for each rank count and filtered campaign: the manifest describes
that invocation, not an automatic merge of previous invocations.

`--only` selects a name substring for thermal and benchmark campaigns.
`--end-time` overrides thermal campaign duration for pilot checks.
`--timeout` sets a per-case wall-clock limit; a timeout is recorded as failure
124. Drivers return nonzero if any case fails. A failed run supplies a log and return code, never a successful final
result extracted from its last partial output. The solver's exit 4 must also
be retained. Finite-rate R410A cases may fail if they carry a resolved gas
beyond the subcritical dew line; instantaneous dissolution is a different
physical assumption, not an interchangeable numerical repair.

## What each comparison establishes

* `run_validation.py`: paired segregated/coupled oil-air histories at
  10, 5 and 2.5 µs, plus actual mass and gas-volume diagnostics. Same material,
  geometry and boundary conditions across the pair.
* `run_thermal.py`: isothermal versus heating, timestep/mesh refinement,
  and 20-bar reservoir R410A cases at 20, 40, 60 and 80°C. Synthetic caloric
  constants verify the implementation; they do not identify experimental
  heating or desorption rates.
* `run_benchmarks.py`: independent one-dimensional JFO pressure/front
  solution on 180/360/720/1440/2880 circumferential cells, and the first Etsion/Ludwig
  photographed gas-cavity condition on 120 × 32 and 240 × 64 meshes. An
  additional zero-absorption case probes gas retention without fitting a rate.
* `plot_campaigns.py`: thermal, resolution, JFO and cavity-extent plots from
  a selected evidence folder, plus the separately labeled FVC68B viscosity
  reference. Failed cases are excluded explicitly from successful curves.

For Etsion/Ludwig, the selected experimental cavity runs from 200° to 75°
through the angular wrap, giving 235° extent. Its two endpoint uncertainties
of ±5° give a conservative ±10° extent bound; that is a derived bound,
not a separately reported standard deviation. The gas charge is assumed,
not measured in the paper. Compare a range of liquid-deficit detection
thresholds rather than fitting one threshold to force agreement
(Etsion and Ludwig, *[J. Lubr. Technol.](https://doi.org/10.1115/1.3253174)*, 1982).

The source/property scope is documented in
[COUPLED_CAVITATION](../../../docs/COUPLED_CAVITATION.md),
[COUPLED_ENERGY](../../../docs/COUPLED_ENERGY.md) and
[the FVC68B table note](../../ref_data/fvc68b_r410a_1999_viscosity.md).
