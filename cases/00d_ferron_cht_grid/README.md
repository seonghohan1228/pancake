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

## Result

Run on `afml-r6` from hash-verified solver commit `df0802b`, 24 MPI ranks. Both
the PINNED and FLOW_WEIGHTED passes re-closed force balance independently and
met the adaptive coupling tolerance on every grid.

| grid | Q_supply [L/min] | eps | psi [deg] | p_max [MPa abs.] | mean [K] | RMS [K] | peak [C] | peak x [deg] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| G5 | 4.937 | 0.638338 | -40.1890 | 1.45508 | +3.944 | 4.135 | 50.495 | 261.375 |
| G6 | 5.058 | 0.637330 | -40.1132 | 1.45395 | +3.821 | 4.016 | 50.379 | 261.250 |
| G7 | 5.176 | 0.636401 | -40.0442 | 1.45282 | +3.710 | 3.909 | 50.274 | 261.166 |

G6 to G7 passes eccentricity (0.146%), attitude (0.069 deg), peak pressure
(0.078%), peak temperature (0.105 K), and peak location (0.084 deg). It fails
the pre-registered mean-error and RMS-error changes: 0.111 K and 0.107 K against
0.100 K limits. **Overall verdict: FAIL.** G6 is not certified for the quoted
conjugate temperature metrics; G7 is the best available result, not a certified
grid-independent replacement.

### Defect exposed and corrected

The first matrix reported Q_supply = 16.26, 22.04, and 30.65 L/min. That nearly
linear `n_z` growth was not physical grid sensitivity. The inlet row overwrites
the axial boundary link where the full-length groove reaches each bearing end,
but `mass_balance()` restored those discarded faces. Their fixed pressure jump
over `dz/2` generated a false source proportional to `1/dz`.

Commit `df0802b` makes inlet-source and boundary accounting follow the assembled
row precedence for both FBNS and Elrod. A manufactured two-grid regression is
independent of `n_z` to 3.5e-16 relative. The corrected supply flow has a much
smaller 4.94 to 5.18 L/min trend. Because the old FLOW_WEIGHTED pass used the
contaminated flow, its favorable +0.02 K mean / 1.17 K RMS / 46.71 C peak is
superseded. The corrected G6 comparison is +3.82 K mean / 4.02 K RMS / 50.38 C
peak. Its midplane pressure peak is 1.353 MPa gauge at 199.75 deg, +7.31% and
1.25 deg from the digitised experimental peak.
