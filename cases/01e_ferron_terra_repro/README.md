# Case 01e: Terra Ferron reproduction and boundary-condition ladder

## Discretisation status (closed 2026-08-11)

Both discretisation ladders are complete; the frozen campaign protocol is
unchanged by them.

- **D-grid, G1--G8.** The 1620 x 432 rung closed the ladder at the fine end.
  G7-to-G8 passes every registered gate with the worst check at 20 % of its
  limit: eccentricity 0.050 % (limit 0.5 %), attitude 0.030 deg (0.25 deg),
  `p_max` 0.007 % (1.0 %), peak bearing temperature 0.012 K (0.25 K), peak
  location 0.055 deg (2.0 deg). **G6 remains the campaign resolution.**
- **D-dt, extended to 0.10 s.** The 40/20/10/5 microsecond ladder now runs to
  0.10 s of physical time. T2-to-T3 and T3-to-T4 both pass: final eccentricity
  within 0.0012 %, attitude within 0.0004 deg, load to machine precision, film
  temperature identical. **`dt = 2e-5 s` therefore holds over the longer
  settling window, not only over 0.06 s.**
- The 40 microsecond rung stays a report-only stability probe. Its final-frame
  load reads 11 kN against 4016 N for the finer rungs; that is the frame landing
  on a transient oscillation, visible directly in
  `terra_report_timestep_convergence.png`, not a different equilibrium. Its
  50 ms mean eccentricity differs from the converged rungs by 0.15 %.

Ledgers: `d_grid_campaign_{equilibrium,analysis}.json`,
`d_timestep_100ms_runs.json`, `d_discretization_100ms_analysis.json`,
`d_timestep_100ms_history.json`.

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
| D-grid | `d_grid_g1`--`g8` | 60x16 through 1620x432 | R1 thermal closure | fresh steady CHT grid check |
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

PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_timestep.py --extended-100ms

python3 cases/01e_ferron_terra_repro/analyse_campaign.py
python3 cases/01e_ferron_terra_repro/analyse_discretization.py
python3 cases/01e_ferron_terra_repro/export_visual_data.py
python3 cases/01e_ferron_terra_repro/plot_campaign.py
python3 cases/01e_ferron_terra_repro/plot_discretization.py
```

For one unattended sequential batch, launch
`run_g8_dt100_extension.py` while G8 is active; it waits for the G8 ledger entry
before starting the 0.10 s ladder, analysis, and history export.

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
  480x128, 720x192, 1080x288, and 1620x432 grids. It measures spatial discretisation
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
The separately archived `--extended-100ms` mode reruns only 40/20/10/5
microseconds to 0.10 s for the final convergence figure.
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
python3 cases/01e_ferron_terra_repro/analyse_discretization.py --extended-100ms
python3 cases/01e_ferron_terra_repro/export_timestep_history.py --extended-100ms
python3 cases/01e_ferron_terra_repro/plot_report.py
```

`plot_report.py` produces the requested total-cell grid convergence, full
transient/timestep convergence, Ferron pressure and R1-G7 bearing-interface
temperature validation, and CHT temperature comparison figures. It also plots
the time reached by every coarse stability
probe, including runs that cannot enter the common endpoint comparison.
The CHT comparison includes fixed-wall, depth-averaged, outer-Robin,
fixed-outer-temperature CHT, and fully adiabatic cross-film treatments. The
outer-Robin curve is a heat-transfer boundary, not a prescribed wall-flux mode.
Convergence panels reserve 10% of their data range above the maximum and 40%
below the minimum.

Claude's original Case 00 instead used: (i) fixed-eccentricity, steady,
isothermal JFO grids 60x16 to 480x128 and a later FBNS timestep tail; and
(ii) a **physical** 4000 rpm / 6 kN moving-bearing transient with a 10 kg
bearing. Its nominal THD pair is not used in D-grid; it is reproduced only in
the separately labelled, minimally repaired legacy mode below because its
archived viscosity was thermally inert.

## F6: 4000 rpm / 6000 N with digitised Ferron supply flow

This mode targets Ferron's Figure 6 temperature condition independently of the
older model-derived inlet closure. The force-balance search starts from
`eps=0.58`, uses G6 (`720x192`), and holds the digitised Figure 10 fresh supply
fixed throughout the search:

- `Q*=0.5128803764`, where `Q*=Qi/(L*c*U)`;
- `L=0.08 m`, hot clearance `c=152 um`, `R=0.05 m`, and 4000 rpm;
- `Qi=1.306195764e-4 m3/s = 7.83717 L/min`.

Because `Qi` is supplied independently, F6 runs one direct `FLOW_WEIGHTED`
force balance. It does not perform the old PINNED pass or replace the measured
flow with the solver's groove-source diagnostic.

The four registered cases are depth-averaged THD, cross-film adiabatic,
cross-film with both walls fixed at 40 C, and cross-film plus CHT with the bush
outer cylinder fixed at `Tb=40 C`. There is no `h=80` case: Ferron's `h_b=80`
W/(m2 K) belongs to the axial bush faces, and axial cooling is deliberately
deferred. The CHT result is therefore a reduced-model comparison, not a full
reproduction of Ferron's three-dimensional bush conduction.

Run on the isolated Ubuntu checkout with one 24-rank process group:

```bash
PANCAKE_EXE=/path/to/pancake PANCAKE_MPI=mpirun PANCAKE_NPROC=24 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/01e_ferron_terra_repro/run_campaign.py F6

python3 cases/01e_ferron_terra_repro/export_fig6_visual_data.py
```

After checksum-verified transfer, plot on Windows:

```bash
python cases/01e_ferron_terra_repro/plot_fig6.py
```

Ferron's Figure 6 triangles are thermocouple readings approximately 0.5 mm
below the bearing bore. The depth-averaged curve is a film-average temperature;
the three cross-film curves report the bearing-side interface temperature.

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
