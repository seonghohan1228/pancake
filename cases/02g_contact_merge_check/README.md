# 02g - Contact-model merge check (October 2026)

Gate for merging `codex/integrate-2026-10` into `main`. The branch keeps the
Greenwood-Tripp + Patir-Cheng model (`75299d8`) and ports from
`codex/gt-contact` (`86a45a1`) the exact $F_{5/2}$ table, the direct
`asperity_gt_k` input and the legacy key aliases. See PHYSICS.md, "Mixed
lubrication", and CHANGES.md, 2026-10-28.

Status (2026-10-28): prepared, **not run** (Seongho: have the runs ready, do not rerun yet).

## Check sets

| Set | Cases | What it shows | Pass criterion |
|---|---|---|---|
| `af-smooth` | A-F, contact off, full load + thermal closure from `experiment_recalibrated_20261023.json` | the merge changed nothing with contact off | torque equals `references/smooth_recalibrated_AF.csv` (`M_cfd_new`) to the printed digits; all `load_and_thermal_closed` |
| `af-contact` | A, B: smooth aligned, smooth tilt, mixed tilt sourced, mixed tilt upper | sourced parameters with the exact $F_{5/2}$ | smooth variants equal `references/mixed_lubrication_AB.csv`; mixed variants may differ only through contact at $H<8$ (the 10-23 fit was zero for $H\ge4$) |
| `mixture` | M01, M03, M05, M11, M16: smooth, sourced, upper | contact where films are thin (M01/M03, $h_{min}$ 1.9-2.5 µm) | smooth equals `references/mixed_sensitivity_summary_20261023.json`; sourced/upper reported, no pass value (10-23 used generic parameters) |

Parameters: `cases/contact_parameters_sourced_20261022.json` (Seongho's decision,
2026-10-28). Setup: A-E use `steady_setup_fbns1e-8.json`, F the default
`steady_setup.json`; both choices were checked against the 10-23 record hashes
(`aeac5d95...`, `6a6670b1...`).

## Commands

Build the branch first (Linux and native Windows; see BUILDING.md), then run
from this folder:

```bash
S=../../build-linux/pancake
W=/path/to/runs/merge_check_2026_10
python3 check.py af-smooth  --solver $S --work $W/af_smooth --ranks 2
python3 check.py af-contact --solver $S --work $W/af_contact --closed $W/af_smooth --cases A,B --ranks 2
python3 check.py mixture    --solver $S --work $W/mixture --ranks 4 \
    --base /home/user/pancake-runs/lge-mixture-density-off-20261019/production
python3 check.py compare    --work $W
```

Add `--dry-run` to any run command to write only the configs and manifests.
The `mixture` base is the accepted density-off production of 2026-10-19 on the
workstation (the same base as the 10-23 sensitivity), so the smooth controls are
directly comparable. For the recalibrated loads use the 2026-10-23
density-recal production instead and compare smooth to smooth only.

Run cost (2 ranks, 720x192): `af-smooth` about 6 x 3-4 thermal passes x 3-5 min;
`af-contact` 8 runs x 4-5 min; `mixture` 15 runs x 2-3 min at 4 ranks.

## Provenance

Every run directory gets `manifest.json`: commit, branch, dirty flag, solver
SHA-256, parameter-file SHA-256, source run, command, return code. A dirty tree
is refused unless `--allow-dirty` (the diff is then saved as `source.diff`).
