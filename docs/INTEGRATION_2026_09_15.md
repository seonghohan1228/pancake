# Integration verification — 2026-09-15

## Result

The journal-bearing, Ferron, FW68D/R410A and pending optical-rig code is combined
on main. The unfinished bubble branch is retained as a separate optional target.
Verified solver source commit: `c7e29b0745087b271aef88ff8c07eab909190dd4`.
Subsequent documentation commits do not change the tested solver source.

| Check | Native Windows (ATLAS) | Remote workstation (afml-r6) |
|---|---|---|
| Release solver build | GCC 16.1, MS-MPI, PETSc 3.25.1 | GCC 11.4, Open MPI 4.1.2, PETSc 3.15.5 |
| GUI build | Passed | Windows-only |
| Journal-bearing CTest | 52 passed | 52 passed |
| Optional bubble CTest | 13 passed | 13 passed |
| FW68D importer unit check | Passed | Passed |
| FW68D density smoke runs | 4 passed: 2 models × 1/2 ranks | 4 passed: 2 models × 1/2 ranks |
| Bubble transient smoke runs | 3 passed: 1/2/4 ranks | 3 passed: 1/2/4 ranks |
| Optical-rig template smoke | Passed, 2 ranks | Passed, 2 ranks |

The initial workstation transfer omitted the composition-viscosity CSV fixture.
After supplying the unchanged file from the same commit, that test passed on
rerun; the other 64 had already passed. Native bubble tests initially exposed
out-of-rank mask indexing and a false near-pole velocity assertion. Both were
fixed and all 13 bubble checks rerun successfully before transfer.

At 40 °C, FW68D density is 940 kg/m³ with CONSTANT and 925.764 kg/m³ with
LINEAR_TEMPERATURE. The 128-cell outputs give dynamic viscosities
0.033438291 and 0.032931879 Pa·s, respectively (rounded). The 593 active bubble
cells produce identical saved thickness values at 1/2/4 ranks, ranging from
9.99804e-6 to 1.00089e-5 m after 0.003 s. Every smoke process returns zero;
the optical-rig steady solve reports converged coupling. Recorded scalar ranges
also match between Windows and Linux. These are regression/smoke checks, not
new validation of the full experimental bubble physics or a rerun of the long
Ferron/optical operating-window campaigns.

## Workstation and evidence

Fresh build and run directory:
`/home/user/pancake-runs/integration-c7e29b0-20260915`.
The pre-existing `/home/user/pancake` checkout and previous runs were preserved.

The directory contains source/fixture archives, `build-linux`, `logs` and `smoke`.
Logs record the hostname, build versions, source archive and executable SHA-256
hashes. `smoke/summary.json` records each input hash, MPI count, return code,
wall time and checked field ranges. Results and logs were copied back to
`C:/Dev/pancake-integration-backup-20260915` along with the reproducible smoke
script and native records.

| Artifact | SHA-256 |
|---|---|
| Source archive | `630c7ab0b13f2f00c5a6490590889502de0c2365fe1b60b9e47b22c0cf7c702e` |
| Fixture archive | `4a2023b751f9194f06b2a0bea14f8f542a59f0c6536b2ecc8b65517079bb24c5` |
| Linux pancake | `74c74194ac622a7596d0697936c17db502d39fa4f620ca1a80dcc2cc709b7eb7` |
| Linux bubble | `f622e560c39af0d5e4b477c0970d71709f719cbdf742137704b6c2ea3e447639` |

## Branch and archive preservation

- Integrate the Ferron branch (446327c) and pending product documentation and
  optical-rig scripts. Preserve local agent settings separately.
- Workstation commit a7a992e has identical source/tests to f0e772d already in
  main; merge its history while retaining the more complete local test record.
  The workstation's other four branch tips were already ancestors of main.
- Preserve feature/bubble (72736c6) under `experimental/bubble`; retain the main
  journal-bearing solver. feature/journal_bearing was already an ancestor.
- Three generated tar archives exceed the push size limit. Remove only their
  blobs from unpublished commits in an isolated mirror. Verify that every
  published origin ref is unchanged and that the resulting tip fast-forwards
  origin/main. The before/after tip trees differ only by those three archives;
  the code and data fixtures are identical. No force push is needed.
- Preserve original history in `original-branches.bundle` and
  `integrated-unfiltered.bundle`, the original workstation branches in
  `C:/Dev/pancake-workstation-branches-20260915.bundle`, pending tracked edits in
  patch files, and the commit mapping in `commit-map.txt` under the backup
  directory. Historical commit IDs in campaign records remain meaningful through
  these original-history bundles and mapping.

Excluded generated archives:

- `cases/00_grid_dt_independence/dt_T2_fbns/T2_results.tar.gz`
- `cases/00_grid_dt_independence/dt_T3_fbns/dt_T3_fbns_results.tar.gz`
- `cases/00_grid_dt_independence/dt_T4_fbns/dt_T4_fbns_results.tar.gz`
