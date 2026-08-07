# Case 01-CONJUGATE - results

> **Historical thermal comparison, superseded by case 00d (`df0802b`).** Its
> FLOW_WEIGHTED supply rate came from the same grid-divergent inlet diagnostic.
> Use case 00d for current pressure, temperature, equilibrium, and grid evidence.

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
| `outer_robin80` | `CONVECTION`, h = 80 W/(m^2 K) | originally **aborted** after 7 solver calls (4 ranks: at `fb_iter_06`; 24 ranks: at `fb_iter_03`); **closed on 2026-08-06** after the KSP fix, 11 solver calls - see below |

The abort is not a physics result. It is the rank-dependent linear-solve
divergence, `[FATAL] linear solve 'temperature' failed ... DIVERGED_DTOL`, which
the Python driver surfaces as `KeyError: 'force_angle_deg'` when it reads the
empty metrics dict of an aborted call. Diagnosis, reproducer and fix: commit
`f0e772d` and `cases/90_diagnostics/ksp_rank_robustness/`. `outer_robin80` **has now been re-run**
(2026-08-06, `codex/picard-best-iterate`, 24 ranks) and closes; its numbers are
below. The post-fix production batch in `cases/01_ferron_wallbc_fields/` remains
the reference for the verdict, but h = 80 is no longer an open hole.

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

## `outer_robin80` now closes, and it is nearly as good as h = 10

Re-run on `codex/picard-best-iterate` (24 ranks, afml-r6): `fb_iter_06`, the call
that used to abort on the rank-dependent KSP divergence, runs, and the force
balance closes for the first time.

| | eps | psi [deg] | W [N] | p peak gauge | h_min | cav. frac. | torque |
|---|---|---|---|---|---|---|---|
| `outer_robin80` PINNED | 0.592143 | -36.9500 | 4000.04 | - | - | - | - |
| `outer_robin80` FLOW_WEIGHTED | **0.602770** | **-37.7046** | 4002.17 | 1.3034 MPa | 60.38 um | 0.2645 | 2.627 N m |

Scored against the same 16 Fig. 4 points:

| variant | mean K | rms K | max abs K | peak | peak residual |
|---|---|---|---|---|---|
| `outer_fixed313` FLOW (control, reproduces the archive) | -2.971 | 3.581 | 6.586 | 42.488 C @ 247.75 deg | -6.55 K |
| **`outer_robin80` FLOW (new)** | **-0.201** | **1.214** | 2.615 | 46.450 C @ 264.75 deg | -2.59 K |
| `rerun_cht_robin10` h = 10 (the production group) | +0.017 | 1.173 | 2.358 | 46.712 C @ 265.25 deg | -2.33 K |

This refines the CORRECTION in `README.md` rather than contradicting it. The
-1.31 K / 1.76 K recorded there for h = 80 came from a *near-equilibrium* state,
because the run could not close. Force-balanced, h = 80 gives -0.20 K / 1.21 K -
within 0.26 K of h = 10, and 2.8 K better than `FIXED_TEMPERATURE`. **What
supplies the missing heat is closing the force balance, not the value of h.**
Re-closed independently under `coupling_relaxation_mode = ADAPTIVE`, where every
call converges, the equilibrium is eps = 0.602769 / psi = -37.7045 - the same to
six digits, from 501 outer sweeps instead of 2178.

## The Picard cap, and why it does not change the table above

`outer_fixed313` converges normally (26-31 outer sweeps, `R_p = 0` exactly).
`outer_robin80` does not: four of the six calls that ran to completion hit
`coupling_max_iters = 200` with `R_p` around 3e-4, i.e. 300x above
`coupling_tolerance` (the other two scraped in at 187 and 190 sweeps), and the loop was
returning the last iterate of a limit cycle rather than the best one it had
passed. Fixed on `codex/picard-best-iterate` (best-iterate rollback on cap
exhaustion + loud non-convergence reporting), and the cycle itself is removed by
under-relaxation - `coupling_relaxation = 0.95` converges in 48 sweeps where 1.0
never does.

Measured consequence for the delivered fields: none beyond the printed digits.
On the sibling h = 10 equilibrium the bush-side wall temperature profile moves by
at most 1e-3 K between the iterate-200 state, the rolled-back best iterate and
nine converged relaxation settings. Full report and traces:
`../01_ferron_wallbc_fields/PICARD_FIX_REPORT.md`.

## The frozen-geometry trap, recorded

An earlier diagnostic held the geometry fixed at the converged `outer_fixed313`
point, swapped h_out from 80 to 10, measured +0.2 K and -0.5 deg, and concluded
the outer coefficient could not close the cold bias. **That conclusion was
wrong.** Holding eps fixed removes the channel that matters: re-closing the force
balance lets eps rise and the added eccentricity supplies the heat. Full
correction in `README.md`; the probe is archived on afml-r6 at
`~/pancake/cases/01_ferron_conjugate/diag_h10_frozen_geometry/`.

**Refinement, now that both variants are force-balanced (2026-08-06).** The
probe's *magnitude* was in fact about right: re-closed, h = 10 sits +0.26 K above
h = 80 on the peak, against the probe's +0.2 K. What was wrong was the reference
it was compared against - the h = 80 number then on the books (45.38 C) came from
a never-closed near-equilibrium state and understated force-balanced h = 80 by
about 1.1 K. So the probe's error was not the +0.2 K it measured; it was reading
that 0.2 K as "cannot close a 1.3 K bias" when 1.1 K of that bias was an artefact
of the unconverged reference. The lesson is unchanged and now has two independent
supports. (The +-0.5 deg on peak *location* is one mesh cell and is argmax
quantisation - see `../01_ferron_wallbc_fields/PICARD_FIX_REPORT.md` §4.3.)

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
