# Stage 90 - diagnostics

Not validation cases. Everything here exists to explain or defend a *solver*
behaviour that a validation case tripped over, so the validation case itself can
stay a clean physics comparison. Nothing in this stage is compared against
experiment, and nothing here carries a pass/fail verdict from the user.

| dir | what it defends | verdict |
|---|---|---|
| `ksp_rank_robustness/` | the rank-dependent linear-solve divergence that cost `01_ferron_conjugate/outer_robin80` two force-balance runs; evidence behind commit `f0e772d` | closed by `f0e772d` (ASM + fallback ladder) |
| `remote_harness/` | the loose driver/analysis scripts that lived in `afml-r6:~/` and drove remote batches for stages 00 and 01 | reference copies, not runnable in place |

## Why this stage exists at all

Before this consolidation the two things above lived outside `cases/`, on the
remote host, in `/tmp` and the home directory - i.e. in the two places most
likely to be cleaned. Both are evidence for claims that are made in `CHANGES.md`
and in case summaries, so they are archived here. See `cases/INDEX.md` for where
the bulk data for each of them lives.
