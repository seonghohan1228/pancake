# Case 00d: Ferron conjugate-THD grid sensitivity

## Question

Are the quantities quoted from the final Ferron 1983 conjugate comparison
insensitive to film-plane grid refinement around the campaign grid G6?

This case closes a specific limitation in `cases/VALIDATION_STATUS.md`. It does
not repeat the generic case-00 study and does not claim universal grid
independence.

## Frozen physics and operating point

- Ferron et al. (1983), Fig. 4: 2000 rpm and 4 kN.
- FBNS cavitation, energy equation, Ferron polynomial viscosity.
- Polynomial cross-film temperature and 2-D bush conduction.
- Bush: `k = 250 W/(m K)`, outer radius 0.10 m, 16 radial cells, outer
  convection `h = 10 W/(m^2 K)` to 313.15 K.
- Globally adiabatic shaft.
- Two-pass groove closure: `PINNED`, then `FLOW_WEIGHTED` using the first pass's
  converged supply flow.
- Force balance is re-closed independently on every grid. Frozen-geometry
  comparisons are not accepted.
- Steady pressure-temperature coupling uses the untuned `ADAPTIVE` relaxation
  schedule and must meet `coupling_tolerance = 1e-6`; best-iterate rollback is
  not accepted as convergence.

The thermal block and force-balance algorithm are imported from
`cases/01_ferron_wallbc_fields/force_balance.py` so this case cannot silently
drift from the production comparison.

## Grids

The refinement ratio is 1.5 in both in-plane directions.

| tag | theta x z | role |
|---|---:|---|
| G5 | 480 x 128 | coarse |
| G6 | 720 x 192 | campaign |
| G7 | 1080 x 288 | fine |

Bush radial cells and cross-film quadrature points stay fixed. This is a
film-plane grid check; it does not establish radial-bush or cross-film-order
independence.

## Pre-registered G6-to-G7 acceptance checks

The experiment has about 1 K local scatter and about 20 deg disagreement between
open digitisations of the plateau maximum. Refinement effects must be small
relative to those limits:

| quantity | acceptance |
|---|---:|
| equilibrium eccentricity ratio | <= 0.5% relative |
| attitude angle | <= 0.25 deg |
| cell-centred peak pressure | <= 1.0% relative |
| mean error over digitised bush temperatures | <= 0.10 K change |
| RMS error over digitised bush temperatures | <= 0.10 K change |
| predicted bush-temperature peak | <= 0.25 K change |
| predicted peak location | <= 2.0 deg circular change |

A non-monotone triplet is reported as such; no Richardson/GCI number is forced
onto oscillatory data. Passing this table supports the quoted Ferron condition
only.

## Execution

Develop and smoke-test on Windows, then run the production matrix on the Ubuntu
workstation from the exact committed hash:

```bash
PANCAKE_EXE=/path/to/pancake \
PANCAKE_NPROC=24 PANCAKE_MPI=mpirun \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python3 cases/00d_ferron_cht_grid/run_grid.py

python3 cases/00d_ferron_cht_grid/analyse_grid.py
```

`python run_grid.py --smoke` uses 180 x 48 and is not scored. Solver output stays
ignored inside each run directory; configs, logs, equilibrium JSON, and the
analysis are the provenance record.
