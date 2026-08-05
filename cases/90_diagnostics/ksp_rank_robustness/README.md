# 90-KSP - rank-dependent linear-solve divergence, minimal reproducer

Evidence behind `CHANGES.md` 2026-08-05 "Rank-robust linear solves: ASM
preconditioning and a fallback ladder" (commit `f0e772d`). Originally run in
`afml-r6:/tmp/ksp_repro`; relocated here 2026-08-05 so the evidence outlives
`/tmp`.

## What failed, and where it hurt

The Ferron conjugate-bush campaign (`../../01_ferron_conjugate`) lost its
`outer_robin80` variant twice, at different force-balance iterates on different
rank counts:

- 4 ranks (local): died at `outer_robin80/fb_iter_06`.
- 24 ranks (afml-r6): died at `outer_robin80/fb_iter_03` - the aborted archive is
  `~/pancake/cases/01_ferron_conjugate/r6_24rank_rerun/` on the remote host.

The Python driver reported it as `KeyError: 'force_angle_deg'`, which is just the
driver reading an empty metrics dict after the solver called `MPI_ABORT`. The real
line is in the solver log:

    [FATAL] linear solve 'temperature' failed after 25 iterations:
            KSP diverged (DIVERGED_DTOL, reason -4)

Deterministic per decomposition, not a lottery: the same config converges at 1, 2
and 6 ranks and aborts at 4 and 24.

## The reproducer

Each directory is one solver call at a *frozen* geometry taken straight from the
failing force-balance iterate, so no Newton search is involved and the failure is
one solve deep. Both operating points are the 720 x 192 Ferron Fig. 4 point
(2000 rpm / 4 kN, FBNS + `ENERGY_EQUATION` + cross-film + `CONDUCTION_2D` bush,
`bush_outer_heat_transfer = 80`):

| dir | binary | geometry | outcome |
|---|---|---|---|
| `base_i03` | pre-fix | `fb_iter_03` (eps psi = -37.013497) | **aborted**, 2 x `DIVERGED_DTOL` on `temperature` |
| `base_i06` | pre-fix | `fb_iter_06` (psi = -36.950026) | finished, 47374 Krylov iters |
| `i03` | post-fix (`f0e772d`) | same as `base_i03` | finished; hit the same 2 divergences, **recovered on fallback solver 1** |
| `i06` | post-fix | same as `base_i06` | finished, 47374 iters (unchanged) |
| `i03_bcgsl` | post-fix + BCGSL | `fb_iter_03` | finished clean, 46618 iters |
| `i03_bj_bcgsl` | post-fix + block-Jacobi + BCGSL | `fb_iter_03` | finished clean, 49910 iters |

The two `bj_bcgsl` / `bcgsl` variants were the option-space probe that separated
"the preconditioner truncates the theta coupling" from "BiCGStab's non-monotone
residual turns a poor preconditioner into a divergence". Both effects are real:
ASM fixes the first, the fallback ladder covers the second.

`i03` is the load-bearing row. The post-fix binary still *meets* the bad solve -
ASM does not make the operator well conditioned - but the fallback ladder turns
an `MPI_ABORT` into a warning and a converged step, at the cost of one extra
solve. That is why the fix is described as a ladder and not a cure; the residual
conditioning issue is listed as still open in `CHANGES.md`.

## Files

- `run.sh` - the harness: `run.sh <exe> <dir> <ranks> [petsc opts...]`. Wipes
  `results/`, pins `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, binds to core,
  tees to `run.log`, then prints OK/FAIL, wall time, solve count, fallback
  recoveries and total Krylov iterations.
- `compare_runs.py` - field-by-field max-abs difference between the final flat
  frames of two run directories. Used to confirm the fix changes no converged
  field. Needs `solver_io.py` on `sys.path`; its hard-coded `/tmp/ksp_repro`
  path is stale here, pass the directory instead.
- `solver_io.py` - the shared pvts reader (same file as the 01 case harnesses).
- `<run>/config.txt`, `<run>/run.log` - inputs and full solver output, kept.

## Where the bulk data is

The `results/` trees (137 MB, 24 ranks x 6 runs) stay on afml-r6 at
`~/pancake/cases/90_diagnostics/ksp_rank_robustness/`. Only configs, logs and
scripts are mirrored locally - the run logs already contain every number quoted
above and in `CHANGES.md`.

## Reproducing

Needs a Linux build and 24 ranks; not runnable on the Windows side as archived
(`run.sh` is bash and the rank count is the point). On afml-r6:

    bash run.sh ~/pancake/build-linux/pancake <dir> 24
