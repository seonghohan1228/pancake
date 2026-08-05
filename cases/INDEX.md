# cases/ - validation campaign index

Front door for the validation campaign. Every case, test and diagnostic run for
the campaign lives under `cases/`; nothing is kept in a home directory, `/tmp` or
a scratch tree any more. Read `VALIDATION_PLAN.md` for the ladder and the decision
gates, `REFERENCE_DATA.md` for what the digitised reference curves are and what is
still missing, then this file for where anything actually is.

Last consolidated: 2026-08-05 (branch `codex/cases-cleanup`, from `f0e772d`).

## Stage numbering

The directory prefix is the campaign stage, and it deliberately matches the case
numbers used in `VALIDATION_PLAN.md`, in `PLAN.md` and in the commit history:

| stage | meaning |
|---|---|
| `00*` | discretisation and solver verification - no experiment involved |
| `01*` | Ferron et al. (1983) validation family - the load-bearing comparison |
| `02*` | static locus (load capacity vs eccentricity) |
| `90*` | diagnostics: solver behaviour, not physics; no verdict |

Stages 03 (transient to equilibrium), 04 (second THD dataset) and 05 (gaseous
cavitation) are planned but have no directory yet.

**No directory was renamed during the consolidation.** The existing prefixes
already encode the stage, they are the vocabulary that `VALIDATION_PLAN.md`,
`PLAN.md`, `CHANGES.md` and every case summary already speak, and roughly sixty
archived `config.txt` and `run.log` files name their own case in a header comment
- renaming would leave those archives self-inconsistent for no functional gain.
The only structural addition is stage `90_diagnostics/`, which is new and
therefore breaks nothing. If you are looking for a `10_*` Ferron family or a
`20_*` locus from an earlier plan, they are `01_*` and `02_*` here.

## The cases

Status is the *user's* verdict where one exists; the agent-side state is in
parentheses. "Bulk data" is where the `results/` VTK trees live - local means the
Windows worktree, r6 means `afml-r6:~/pancake/cases/`.

| stage | dir | purpose | status | key figures | bulk data |
|---|---|---|---|---|---|
| 00 | `00_grid_dt_independence/` | grid + timestep independence on the Ferron geometry; freezes G6 720x192 and dt = 2e-5 s for every later case | **PASSED** (user, 08-05) | `gci_00a_fbns.png`, `dt_00b_fbns.png`, `e_vs_t_00b.png`, `gci_00a.png` (superseded JFO generation) | both (13489 local / 22570 r6; T2-T4 tails are r6-only, 5.4 GB each) |
| 00c | `00c_jfo_verification/` | exact 1-D JFO benchmark - solver correctness, no experiment | **CLOSED** by replacing Elrod with FBNS (front <= 1 cell) | none (tabular, `extract_front.py` output) | both (131 / 131) |
| 01 | `01_ferron_prelim/` | Fig. 6 point (4000 rpm / 6 kN), isothermal, two bracketing viscosities | done (agent); superseded as a *result* by 01-THD, kept as the isothermal baseline and the sign-convention derivation | `ferron_midplane_pressure.png` | both (117 / 117) |
| 01 | `01_ferron_thd/` | Fig. 4 point (2000 rpm / 4 kN) + Fig. 5 load line, energy equation + Ferron polynomial viscosity, `cross_film = NONE` | done (agent); under-predicts T by 5-7 K, which motivated the cross-film case | `ferron_thd_pressure.png`, `ferron_thd_temperature.png`, `ferron_thd_eccentricity.png` | both (872 / 872) |
| 01 | `01_ferron_thd_crossfilm/` | same point with `cross_film_thermal_model = POLYNOMIAL`, three wall-BC variants | done (agent); scored the wall BC and found the peak-location failure of any local closure | `ferron_crossfilm_temperature.png`, `ferron_crossfilm_pressure.png`, `ferron_Tbearing_comparison.png` | both (347 / 347) |
| 01 | `01_ferron_conjugate/` | same point with `bush_thermal_model = CONDUCTION_2D` (2-D bush conduction), two outer BCs | **verdict pending** (user); this is the CHT verdict campaign | `ferron_conjugate_temperature.png`, `ferron_conjugate_pressure.png` | both (173 local / 656 r6 - r6 also holds the 24-rank rerun and two probes, see below) |
| 01 | `01_ferron_wallbc_fields/` | final four-way wall-BC comparison with the film-surface T and wall-flux fields exposed; first production use of the KSP fix | **verdict pending** (user) | `ferron_wallbc_temperature.png`, `ferron_wallbc_heatflux.png` | r6 (95 local / 1170 r6 - 983 MB of VTK stays remote by decision) |
| 01 | `01_ferron_vs_reynsim/` | three-way overlay: experiment vs the previous-generation MATLAB/ReynSim solver vs pancake | done (agent); "This year" inside the band, "Last year" -17.6 % on the peak | `reynsim_midplane_pressure.png`, `reynsim_eccentricity.png`, `fig4_three_way.png` | both (98 / 98) |
| 02 | `02_static_locus/` | eps vs Sommerfeld and attitude locus at L/D = 1 and 0.5 | **POSTPONED** (user, 08-05) - no reference table in the library; solver side complete and kept as supporting material | `static_locus.png` | both (268 / 268) |
| 90 | `90_diagnostics/ksp_rank_robustness/` | minimal reproducer for the rank-dependent `DIVERGED_DTOL`; evidence behind `f0e772d` | closed by `f0e772d` | none | r6 (configs + logs local, 137 MB of results r6) |
| 90 | `90_diagnostics/remote_harness/` | reference copies of the loose `afml-r6:~/` driver and analysis scripts | reference only | none | both (text only) |

Shared, not a case: `ref_data/` (digitised Ferron Figs. 4, 5, 6 - CSV + source
PNG, format defined in `REFERENCE_DATA.md`), `sync_to_r6.sh` (manifest mirror to
afml-r6), `sync_*.log` (its transcripts).

## What was folded in, and from where (2026-08-05)

| was | now | why it is not a duplicate |
|---|---|---|
| `afml-r6:~/cases_run/01_ferron_conjugate` | `r6:~/pancake/cases/01_ferron_conjugate/r6_24rank_rerun/` | a **different run**: 24 ranks, whose `outer_robin80` leg aborted at `fb_iter_03`. Its `outer_fixed313` leg reproduces the 4-rank archive to every printed digit (eps = 0.593397, psi = -37.5036), which is the rank-independence evidence. Kept whole, 421 files. |
| `afml-r6:/tmp/ksp_repro` | `cases/90_diagnostics/ksp_rank_robustness/` (configs, logs, scripts) and `r6:~/pancake/cases/90_diagnostics/ksp_rank_robustness/` (bulk) | new case, nothing to duplicate |
| `afml-r6:~/diag_h10` | `r6:.../01_ferron_conjugate/diag_h10_frozen_geometry/` | the frozen-geometry h = 10 probe whose conclusion was later corrected; kept as the record of the wrong turn |
| `afml-r6:~/gate_fixed313` | `r6:.../01_ferron_conjugate/gate_fixed313/` | the cross-platform gate run for that case |
| `afml-r6:~/smoke` | `r6:.../01_ferron_wallbc_fields/smoke_coarse_180x24/` | coarse 180 x 24 config smoke for the four wall-BC groups |
| loose scripts in `afml-r6:~/` | copied to `cases/90_diagnostics/remote_harness/` (originals left in place) | copies, deliberately - the originals still work from `~` |

`afml-r6:~/cases_run/` is now gone (it was empty after the move). Nothing was
deleted anywhere: every move was `mv`, every script pull was `cp`.

## Deliberately left where it is

- `01_ferron_vs_reynsim/bin/` - a pinned Windows `pancake.exe` plus its MinGW
  DLLs. It is the binary that produced this case's numbers, so it belongs to the
  case; it is untracked (`*.exe`, `*.dll` are gitignored at the repo root).
- `01_ferron_vs_reynsim/smoke/`, `01_ferron_vs_reynsim/diag/` - the coarse smoke
  run and the four-run KSP diagnostic set (`serial`, `repeat`, `relax_0.5`,
  `relax_0.7`) that `CHANGES.md` quotes for the ASM iteration table. In-case,
  because they are that case's own provenance.
- `validation/` at the repo root - the *pre-campaign* analytical benchmark
  harness (Christopherson, Brito, cross-solver, grid convergence, `case_archive/`).
  Tracked source, referenced by `PLAN.md` WP-7, and not part of this campaign's
  ladder. Its `_work*/` scratch trees are gitignored. Do not merge it into
  `cases/`.
- `afml-r6:~/dt_T{3,4}_fbns_results.tar.gz`, `~/T2_results.tar.gz`,
  `/tmp/remote_only_pull.tar.gz` (330 MB total) - compressed **transfer**
  snapshots of trees that are already unpacked under
  `00_grid_dt_independence/`. Redundant, kept because nothing is deleted without
  the user saying so. Safe to remove on request.
- `%TEMP%\claude\C--Dev-pancake\<session>\scratchpad\` on Windows - agent scratch
  (ad-hoc probe runs, PDFs, throwaway builds, and a detached-HEAD git worktree
  `base_b74db94` used for an A/B binary comparison). Ephemeral by construction and
  **not** campaign data; nothing there is referenced by any case. It is not
  mirrored and will not survive a temp clean.

## Known provenance inconsistencies

- `00_grid_dt_independence/dt_fbns_analysis.json` records `end_t = 0.15` for
  `dt_T2_fbns`, `dt_T3_fbns` and `dt_T4_fbns`, while the `config.txt` now on disk
  in those directories says `end_t = 0.06`. The JSON is the record of the run as
  executed and has been left alone; re-running `analyse_dt_fbns.py` will rewrite
  those three fields to 0.06. Do not treat the regenerated file as authoritative
  without checking the run logs.
- `01_ferron_conjugate/conjugate_analysis.json` was missing from the local archive
  and was regenerated on 2026-08-05 from the local 4-rank run directories
  (`analyse_conjugate.py`, no arguments, deterministic).

## Housekeeping rules for this directory

1. Never delete or overwrite a `results/` tree. Superseded output is renamed
   (`results_*stale*`, `results_partial_*`) with its `run_*.log` kept alongside.
2. One case = one directory, with `README.md` (what and why, before the run) and
   `SUMMARY.md` (what came out, after). Every solver call archived as its own
   subdirectory with the `config.txt` that produced it.
3. `.gitignore` in this directory keeps `results*/` and `run*.log` untracked. The
   configs, scripts, JSON metrics, CSVs, figures and prose are tracked. Bulk VTK
   is never committed; it is mirrored to afml-r6 instead.
4. Mirror with `bash cases/sync_to_r6.sh` after each batch. It is
   manifest-based, append-only and idempotent, ships only what differs by
   (path, size), and reports remote-only files at the end - which is how remote
   run output gets noticed.
5. Analysis and plot scripts must run standalone: `python <script>` from inside
   the case directory, no solver run needed. All 20 of them were verified to do
   so on 2026-08-05.
