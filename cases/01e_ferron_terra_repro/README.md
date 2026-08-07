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

## Pre-registered scope

| mode | runs | grid | boundary condition | purpose |
|---|---|---:|---|---|
| R0 | `r0_g6`, `r0_g7` | 720x192, 1080x288 | legacy 2-D bush, outer Robin `h=10`, shaft adiabatic | corrected reproduction of case 00d |
| R1 | `r1_g6`, `r1_g7` | 720x192, 1080x288 | outer cylinder fixed at 40 C, shaft adiabatic | Ferron outer-cylinder condition |
| R2 | `r2_g6`, `r2_g7` | 720x192, 1080x288 | R1 plus scalar shaft held near 43 C | explicit surrogate sensitivity |
| C | Fig. 4 THD, Fig. 5 six-load, cross-film A/B/C, conjugate fixed/Robin80 | G6 | archived Claude matrix, corrected source | reproduction comparison |

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
