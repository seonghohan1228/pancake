# Case 01e: Terra Ferron reproduction and boundary-condition ladder

## Purpose

This is a new, isolated campaign.  It never reads prior case JSON, VTK, or
solver logs.  It re-runs selected archived configurations with the corrected
inlet-source accounting and records all solver calls below this directory.

The labels are deliberate:

- **MODEL-CLOSURE** means the groove supply is inferred from this model's
  PINNED-pass steady mass balance.  It is not Ferron's measured fresh flow.
- **SURROGATE** means a spatially uniform scalar shaft temperature.  It is not
  Ferron's axial shaft-temperature boundary condition.
- No temperature agreement is an experimental validation claim until the
  fresh supply and all thermal boundaries are independently established.

The present solver still omits Ferron's bush axial-face convection and full
`(r, theta, z)` bush conduction. R2 clamps the existing scalar shaft near 43 C;
it does not implement Ferron's `(r, z)` shaft with 43 C axial ends. Consequently
R0--R2 can receive numerical and grid verdicts, but not a Ferron temperature-
validation PASS.

For this campaign, axial bush cooling is deliberately neglected by user
decision. R1/R2/D instead enforce Ferron's outer-radius boundary *form*,
`T(R2)=Tb`. Ferron reports `Ta=40 C` and `Ti=40 C` but does not publish a
numerical `Tb`; the configured `Tb=40 C` is therefore an explicit `Tb=Ta`
assumption, not a measured paper value.

## Pre-registered scope

| mode | runs | grid | boundary condition | purpose |
|---|---|---:|---|---|
| R0 | `r0_g6`, `r0_g7` | 720x192, 1080x288 | legacy 2-D bush, outer Robin `h=10`, shaft adiabatic | corrected reproduction of case 00d |
| R1 | `r1_g6`, `r1_g7` | 720x192, 1080x288 | outer cylinder fixed at 40 C, shaft adiabatic | Ferron outer-cylinder condition |
| R2 | `r2_g6`, `r2_g7` | 720x192, 1080x288 | R1 plus scalar shaft held near 43 C | explicit surrogate sensitivity |
| C | Fig. 4 THD, Fig. 5 six-load, cross-film A/B/C, conjugate fixed/Robin80 | G6 | archived Claude matrix, corrected source | reproduction comparison |
| D-grid | `d_grid_g1`--`g7` | 60x16 through 1080x288 | R1 thermal closure | fresh steady CHT grid check |
| D-dt | `d_dt_t0..t4` | G6 | R1 closure, moving-bearing 4 kN step | fresh physical-time timestep check |

Every equilibrium executes two independent force-balance passes: `PINNED`, then
`FLOW_WEIGHTED` using the first pass's corrected source.  Coupling must finish
with both residuals below `1e-6`; a best iterate is not accepted.

R2 requires the solver keys `shaft_ambient_temperature` and
`bush_outer_temperature`.  The former is intentionally a uniform scalar
surrogate, not a replacement for Ferron's r-z shaft condition.

## Commands

Develop and smoke-test on Windows.  Run one mode at a time on afml-r6, with no
concurrent solver jobs:

```bash
PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_campaign.py R0

python3 cases/01e_ferron_terra_repro/analyse_campaign.py
python3 cases/01e_ferron_terra_repro/export_visual_data.py
python3 cases/01e_ferron_terra_repro/plot_campaign.py
```

For the fresh discretisation sequence, run the steady grid set first, then the
physical-time ladder, still one process group at a time:

```bash
PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_campaign.py D_GRID

PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_timestep.py

python3 cases/01e_ferron_terra_repro/analyse_campaign.py
python3 cases/01e_ferron_terra_repro/analyse_discretization.py
python3 cases/01e_ferron_terra_repro/export_visual_data.py
python3 cases/01e_ferron_terra_repro/plot_campaign.py
python3 cases/01e_ferron_terra_repro/plot_discretization.py
```

An infrastructure/launch failure may be replaced explicitly with
`run_timestep.py --rerun-failed`; completed or physically incomplete runs with
`rc=0` remain immutable.

Use `--smoke` only for a 180x48 harness check; it is never scored.  Before a
remote batch: verify the remote commit hash, run the Windows/Linux gate, and
checksum the returned archive. Solver calls are isolated below
`runs/production/`; coarse harness calls stay below `runs/smoke/`. The campaign
also writes `campaign_equilibrium.json`, `campaign_analysis.json`, and plots.

## Deliverables and verdicts

`analyse_campaign.py` writes an explicit per-case verdict table.  It reports
completed/coupling/force-balance status, R0 G6-to-G7 grid checks, and measured
Ferron residuals where an output `T_bearing` field is available.  Archive values
are copied into `baseline.json` as labelled comparison data, not imported at run
time.  `plot_campaign.py` writes 300 dpi, no-title PNGs; `export_visual_data.py`
writes the compact local field archive used for profiles and contours.

## D: discretisation add-on

This is deliberately distinct from Claude's Case 00.

- **D-grid** is steady and force-balanced on 60x16, 120x32, 240x64, 360x96,
  480x128, 720x192, and 1080x288 grids. It measures spatial discretisation
  of the current R1 conjugate thermal model, at the Fig. 4 condition.  It is
  not a pseudo-time test: `dt` is ignored by `STEADY_STATE`.
- **D-dt** is an actual physical-time calculation: `TRANSIENT` plus
  `MOVING_BEARING`, released from `eps=0.1` under a 4 kN load at 2000 rpm.  It
  uses the corrected source inferred by D-grid G6 and thus cannot run before
  that case has closed.  It tests the discretisation of the transient path,
  not a replacement for the steady force-balanced solution.

D-dt's ladder is `3.2e-4`, `1.6e-4`, `8e-5`, `4e-5`, `2e-5`, `1e-5`, and
`5e-6` s to `0.06 s`, written every `5e-4 s`. The three largest steps are
intentional stability probes.
Failure of a coarse probe is archived but does not prevent the finer timesteps
from running.
The proposed acceptance comparison is T2-to-T3: final eccentricity <=0.5%,
final force magnitude <=1.0%, final attitude <=0.25 deg, and final maximum
film temperature <=0.10 K; T3-to-T4 must be no worse.  A failed film guard,
non-converged FBNS step, or missing final frame fails that ladder leg. Transient
mode advances one segregated flow/energy sweep per physical timestep, so the
steady Picard-coupling residual is intentionally not defined for D-dt.

After the remote ladder, extract its lightweight eccentricity histories before
transferring the compact JSON and plotting on Windows:

```bash
python3 cases/01e_ferron_terra_repro/analyse_discretization.py
python3 cases/01e_ferron_terra_repro/export_timestep_history.py
python3 cases/01e_ferron_terra_repro/plot_report.py
```

`plot_report.py` produces the requested total-cell grid convergence, full
transient/timestep convergence, Ferron pressure validation, and CHT temperature
comparison figures. Convergence panels reserve 10% of their data range above
the maximum and 40% below the minimum.

Claude's original Case 00 instead used: (i) fixed-eccentricity, steady,
isothermal JFO grids 60x16 to 480x128 and a later FBNS timestep tail; and
(ii) a **physical** 4000 rpm / 6 kN moving-bearing transient with a 10 kg
bearing. Its nominal THD pair is not used in D-grid; it is reproduced only in
the separately labelled, minimally repaired legacy mode below because its
archived viscosity was thermally inert.

## REPAIRED_LEGACY: exact Claude GCT reproduction

`run_legacy_gct.py` is the separate historical reproduction, not D-grid.  It
retains the original 00a grids: isothermal JFO G1--G5 = 60x16, 120x32,
240x64, 360x96, 480x128; fixed `eps=0.6`, `psi=-90 deg`, 4000 rpm, and
`STEADY_STATE`.  The legacy `dt=1e-5 s` is a solver guard only and is ignored.
There is no imposed external load because eccentricity is fixed.

It also carries the original nominal energy-equation G3/G4 pair.  Those inputs
are `REPAIRED_LEGACY`, not exact: the current `SINGLE_PHASE` validator rejects
two zero gas coefficients.  The two coefficients are changed to `1e-12` only;
the machine-readable before/after/reason record is
`legacy_repair_ledger.json`.  The legacy Andrade coefficient remains `0.034 K`,
so that pair is energy-only and must never be described as a thermoviscous grid
check.  These figures are report-only numerical reproductions, not Ferron
validation or a replacement for D-grid.

```bash
PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_legacy_gct.py ISOTHERMAL

PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_legacy_gct.py THD

python3 cases/01e_ferron_terra_repro/analyse_legacy_gct.py
python3 cases/01e_ferron_terra_repro/plot_legacy_gct.py
```
