# Case 01 wall-BC field batch - final four-case comparison (Ferron Fig. 4)

2000 rpm / 4 kN, FBNS, G6, frozen protocol. Binary: codex/thermal-ksp-robustness
@ f0e772d (KSP rank-robustness fix - first production use; rerun_cht_robin10
closed its force balance where it previously diverged on both hosts). Run on
afml-r6, 24 pinned ranks; remote HEAD hash-verified; cross-platform gate
identical to all printed digits (W differs 1e-5 N relative). Reruns reproduce
their archived predecessors (eps exact to 6 digits where conditions identical).

## T_bearing vs the 16 digitised Fig. 4 points (favg is T-bar, labelled)

| group | quantity | peak C | location | mean K | rms K | 0-140 | 140-270 | >270 |
|---|---|---|---|---|---|---|---|---|
| experiment | - | 49.04 | 258.2 | - | - | - | - | - |
| rerun_xf_adiabatic | T_bearing | 50.25 | 347.7 | +0.18 | 2.07 | -0.96 | -0.47 | +3.80 |
| rerun_xf_bush40_J40 | T_bearing | 40.00 | 172.3 | -4.05 | 5.01 | -0.75 | -5.43 | -5.90 |
| rerun_cht_robin10 | T_bearing | 46.71 | 265.2 | +0.02 | 1.17 | +1.06 | -0.47 | -0.43 |
| favg_walls40 | T-bar | 42.14 | 260.7 | -2.83 | 3.69 | -0.26 | -3.97 | -4.07 |

CHT (2-D bush, outer Robin h = 10 to 40 C, shaft globally adiabatic) is the best
variant on every aggregate measure: peak location within 7 deg, mean +0.02 K,
rms 1.17 K, post-peak turndown reproduced. Peak value remains -2.3 K under the
measurement (early arc runs ~+1 K hot instead).

CORRECTION recorded: an earlier frozen-geometry probe claimed h=10 would not
lift the temperature level vs h=80. Wrong - re-closing the force balance at
h=10 moves eps 0.5934 -> 0.6033 and the added eccentricity supplies the heat.
Frozen-geometry probes are not comparable to re-closed equilibria
(cases/01_ferron_conjugate/README.md).

## Net wall heat (W, positive into the fluid)

| group | bush | shaft | total |
|---|---|---|---|
| rerun_xf_adiabatic | 0.0 (exact) | 0.0 (exact) | 0.0 |
| rerun_xf_bush40_J40 | -211.5 | -211.5 | -423.1 |
| rerun_cht_robin10 | -1.9 | -0.0 (exact) | -1.9 |
| favg_walls40 | -208.6 | -289.8 | -498.4 |

Mechanism (heat-flux figure): the conjugate bush injects up to +17 kW/m2 into
the film through the early arc and extracts ~-10 kW/m2 through the loaded arc -
circumferential redistribution made visible. Invariants hold on the solved
fields: adiabatic fluxes and the conjugate shaft net are exactly zero.

## Equilibria

rerun_xf_adiabatic eps 0.608601, rerun_xf_bush40_J40 0.564609,
favg_walls40 0.567497, rerun_cht_robin10 0.603314 (W +2.37 N, angle -90.013).

## Picard cap defect: checked, and the numbers above are unaffected

Every `rerun_cht_robin10` solver call hit `coupling_max_iters = 200` without
meeting `coupling_tolerance = 1e-6`, and the loop was returning the last iterate
rather than the best one it had passed through (6.3x worse in residual). Fixed on
`codex/picard-best-iterate`: the loop now rolls back to its lowest-residual
iterate on cap exhaustion and reports NON-CONVERGENCE loudly.

Re-measured at the same equilibrium with the fix, and against nine converged
under-relaxed settings: the bush-side `T_bearing` midplane profile moves by at
most 1e-3 K (the VTK write precision), W by 0.11 N in 4002 N, p_max by 40 Pa in
1.405 MPa and T_max by 1 mK. **Peak 46.71 C @ 265.2 deg, mean +0.02 K, rms 1.17 K
and the -2.3 K peak residual are unchanged**, and the temperature figure did not
need regenerating.

The loop does converge under under-relaxation - `coupling_relaxation = 0.95` in 48
sweeps, `coupling_relaxation_mode = ADAPTIVE` in 52, against 200-and-never at 1.0
- so the frozen protocol's `coupling_relaxation = 1.0` is now a known-poor choice
for the conjugate variants even though it costs nothing in the answer. Full
report, traces and the run-by-run table: `PICARD_FIX_REPORT.md`, figure
`ferron_picard_convergence.png`.

## Figures / files

ferron_wallbc_temperature.png, ferron_wallbc_heatflux.png, equilibrium.json,
extract.json, per-group fb_iter_NN archives (provenance + extracted profiles
local; bulk VTK archives, 983 MB, remain on afml-r6 pending user decision on
mirroring - see progress.log).
