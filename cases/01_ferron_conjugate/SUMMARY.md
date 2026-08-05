# Case 01-CONJUGATE - results

Ferron et al. (1983) Fig. 4 point (2000 rpm, 4 kN), `bush_thermal_model =
CONDUCTION_2D`, G6 720 x 192, FBNS + `ENERGY_EQUATION` + cross-film polynomial,
globally adiabatic shaft. Protocol and every stated approximation are in
`README.md`. Nothing here is tuned to the experiment.

Numbers below are read out of `conjugate_analysis.json` (regenerated 2026-08-05
by `python analyse_conjugate.py`); `equilibrium.json` carries the converged
points and `<variant>/fb_log.json` every iterate.

## What converged and what did not

| variant | outer BC | state |
|---|---|---|
| `outer_fixed313` | `FIXED_TEMPERATURE` at 313.15 K | **converged**, both the PINNED and the FLOW_WEIGHTED pass; 11 solver calls, 593 s |
| `outer_robin80` | `CONVECTION`, h = 80 W/(m^2 K) | **aborted** after 7 solver calls (4 ranks: at `fb_iter_06`; 24 ranks: at `fb_iter_03`) |

The abort is not a physics result. It is the rank-dependent linear-solve
divergence, `[FATAL] linear solve 'temperature' failed ... DIVERGED_DTOL`, which
the Python driver surfaces as `KeyError: 'force_angle_deg'` when it reads the
empty metrics dict of an aborted call. Diagnosis, reproducer and fix: commit
`f0e772d` and `cases/90_diagnostics/ksp_rank_robustness/`. `outer_robin80` has
**not** been re-run since the fix; the post-fix production batch in
`cases/01_ferron_wallbc_fields/` supersedes it for the verdict.

## Converged equilibrium, `outer_fixed313`

| quantity | PINNED pass (`fb_iter_06`) | FLOW_WEIGHTED pass (`fb_iter_10`) |
|---|---|---|
| eps | 0.586599 | **0.593397** |
| psi (solver input) | -36.9839 deg | **-37.5036 deg** |
| classical attitude phi | 53.02 deg | 52.51 deg |
| W | 4000.15 N (+0.0038 %) | 4000.48 N (+0.0120 %) |
| force angle | -89.99944 deg | -90.0104 deg |
| p peak, gauge | 1.2838 MPa @ 200.25 deg | 1.2936 MPa @ 199.75 deg |
| T_bearing peak | 42.28 C @ 250.75 deg | 42.49 C @ 247.75 deg |
| T_shaft | 43.53 C | 44.45 C |
| h_min | 62.84 um | 61.80 um |
| cavitated area fraction | 0.2611 | 0.2629 |
| friction torque | 2.778 N m | 2.729 N m |
| Q_supply from the pass-1 mass balance | 3.6163e-4 m^3/s (21.70 L/min) | carried in |
| Picard iterations | 26 | 31 |

Worst liquid residual over the group: 1.06e-10. Experiment for the same point:
p peak 1.2602 MPa @ 201.0 deg, T peak 49.04 C @ 258.2 deg, inlet 40.49 C, 16
digitised temperature points.

Pressure agrees to +1.9 % on the peak and 1.25 deg on its location (flow pass),
with no fitted shift.

## Temperature vs the 16 Fig. 4 points

`T_bearing` (bush-side film surface) against the thermocouples, mapping
`solver_theta = (90 + x) mod 360`, groove band +-12 deg excluded from peaks:

| pass | mean K | rms K | max abs K | 0-140 | 140-270 | >270 |
|---|---|---|---|---|---|---|
| PINNED | -3.19 | 3.75 | 6.78 | -1.10 | -3.99 | -4.53 |
| FLOW_WEIGHTED | -2.97 | 3.58 | 6.59 | -0.87 | -3.76 | -4.38 |
| reference: local-bush adiabatic (`01_ferron_thd_crossfilm` variant A) | +0.18 | 2.07 | 6.93 | - | - | - |

**The peak moved, which is what this case was built to test.** The local-closure
adiabatic run put the bush peak at 347.7 deg against a measured 258.2 deg; the
2-D conducting bush puts it at 247.8-250.8 deg, within 8-11 deg. The measured
mid-arc turndown appears, and the wall heat flux shows the mechanism: the bush
injects up to +20 kW/m^2 into the film through the early arc and extracts about
-20 kW/m^2 under the loaded arc.

**The level is too cold**, by about 3 K in the mean and 6.6 K at the measured
peak. That is this variant's outer boundary, not the conduction model: a fixed
313.15 K outer surface is a much stronger sink than the rig. Re-closing the force
balance with a weak outer Robin coefficient (h = 10 W/(m^2 K)) moves the
equilibrium to eps = 0.6033 and brings the mean error to +0.02 K, rms 1.17 K -
see `cases/01_ferron_wallbc_fields/SUMMARY.md`.

## The frozen-geometry trap, recorded

An earlier diagnostic held the geometry fixed at the converged `outer_fixed313`
point, swapped h_out from 80 to 10, measured +0.2 K and -0.5 deg, and concluded
the outer coefficient could not close the cold bias. **That conclusion was
wrong.** Holding eps fixed removes the channel that matters: re-closing the force
balance lets eps rise and the added eccentricity supplies the heat. Full
correction in `README.md`; the probe is archived on afml-r6 at
`~/pancake/cases/01_ferron_conjugate/diag_h10_frozen_geometry/`.

## Archives

- `outer_fixed313/fb_iter_00..10`, `outer_robin80/fb_iter_00..06` - one directory
  per solver call, each with the `config.txt` that produced it (local, 4 ranks).
- `r6_24rank_rerun/` (afml-r6 only) - the 24-rank reproduction. Its
  `outer_fixed313` leg lands on eps = 0.593397, psi = -37.5036, identical to the
  4-rank result to every printed digit; its `outer_robin80` leg hit the same KSP
  abort three iterates earlier. That pair is the rank-independence evidence.
- `gate_fixed313/` (afml-r6 only) - the cross-platform gate run for this case.
- `conjugate_analysis.json`, `equilibrium.json`, `<variant>/fb_log.json`,
  `driver.log`, `progress.log`.

Figures: `ferron_conjugate_temperature.png`, `ferron_conjugate_pressure.png`.

## Verdict

User verdict pending. The peak-location question this case was built to answer is
answered; the temperature-level question is answered by the
`01_ferron_wallbc_fields` batch, not here.
