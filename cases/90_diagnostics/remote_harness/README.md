# 90-REMOTE-HARNESS - the loose afml-r6 driver scripts

Reference copies of the scripts that lived directly in `afml-r6:~/` and drove the
remote batches. They are archived because several numbers quoted in case
summaries were produced by them and by nothing else; they are **not** runnable
from this directory, because most carry absolute `/home/user/...` paths.

| file | drove | belongs to |
|---|---|---|
| `sweep.sh` | Case 00b dt tail T2/T3/T4, n = 24, sequential | `00_grid_dt_independence` |
| `sweep_t3t4.sh` | the T3/T4 rerun with threads pinned (1 thread/rank) | `00_grid_dt_independence` |
| `scaling.sh` | strong-scaling probe on `dt_T2_fbns`, n = 4..56, unpinned | `00_grid_dt_independence` (`scaling_probe_n*`, `scaling_results.txt`) |
| `scaling2.sh` | same probe with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` | `00_grid_dt_independence` (`scaling_results2.txt`) |
| `analyse_remote.py` | Case 00b tail analysis: post-settling means over t in [0.03, 0.06] s, per-step limit-cycle amplitude from `diagnostics.csv` | `00_grid_dt_independence` |
| `spectra.py` | detrended DFT of per-step `cavitated_fraction` - the script that caught the fake 33 Hz "half-whirl" (window bin-1 artifact) | `00_grid_dt_independence` |
| `check_w.py` | raw per-frame `Fp` / `xB` dump; written to resolve `W_std = 0` against a moving eccentricity | `00_grid_dt_independence` |
| `gate.py` | standalone one-run extractor (W, p_max, T_max, h_min, torque) used as the cross-platform gate | stage 01 (all cases) |
| `extract.py` | the wall-BC field extraction (`extract.json`) incl. post-derived Newton-cooling fluxes for the `favg` group | `01_ferron_wallbc_fields` |
| `smoke.py` | coarse-grid (180 x 24) config smoke test for all four wall-BC groups | `01_ferron_wallbc_fields` |
| `diag_h10_diag.py` | the frozen-geometry h = 10 W/(m^2 K) probe | `01_ferron_conjugate` (see its README's CORRECTION note) |

## The frozen-geometry warning attached to `diag_h10_diag.py`

That probe produced the claim, later shown wrong, that h = 10 would not close the
cold bias of the h = 80 conjugate result. Holding `eps` fixed is not the
comparison that matters: re-closing the force balance at h = 10 moves the
equilibrium from eps = 0.5934 to 0.6033 and the added eccentricity supplies the
heat. Full correction in `../../01_ferron_conjugate/README.md`.

## Where the run directories are

The runs these scripts produced live with their cases. The three loose remote run
directories were folded in on 2026-08-05:

- `~/diag_h10` -> `~/pancake/cases/01_ferron_conjugate/diag_h10_frozen_geometry/`
- `~/gate_fixed313` -> `~/pancake/cases/01_ferron_conjugate/gate_fixed313/`
- `~/smoke` -> `~/pancake/cases/01_ferron_wallbc_fields/smoke_coarse_180x24/`
