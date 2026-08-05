# Case 01-CONJUGATE

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Ferron, Frene & Boncompain (1983) Fig. 4 point (2000 rpm, 4 kN) re-run with
`bush_thermal_model = CONDUCTION_2D`: steady 2-D `(r, theta)` heat conduction in
the bush annulus, Ferron's own bush model, in place of the per-cell 1-D radial
closure.

**Why.** Case 01-THD-crossfilm scored the wall boundary condition and found the
one thing no local closure can supply. Its variant A (both walls adiabatic)
reproduced the *magnitude* of the measured bush temperature rise (mean +0.18 K,
rms 2.07 K over the 16 digitised points) but put the peak at **347.7 deg**
against a measured **258.2 deg**, because an adiabatic film keeps heating all the
way round to the groove and a local bush closure sends every watt straight out
through the wall underneath it. The measured mid-arc turndown is bush conduction
physics: heat entering the bush under the loaded arc leaves it somewhere else.
This case tests whether putting that physics in moves the peak.

Nothing here is tuned to the experiment.

## Files

- `config_template.txt` - parametrised config. Frozen from
  `../01_ferron_thd_crossfilm/config_template.txt`; the only changes are the bush
  block and the four extra film-surface output fields.
- `solver_io.py`, `plot_style.py` - copies of the 01-THD harness, unchanged.
- `force_balance.py` - the same two-pass Newton driver, two outer-BC variants.
- `analyse_conjugate.py`, `conjugate_analysis.json` - every reported number.
- `plot_conjugate.py` - the two deliverable figures.
- `equilibrium.json`, `<variant>/fb_log.json` - converged points and every iterate.
- `outer_fixed313/`, `outer_robin80/` - one archived run per solver call.
- `progress.log` - running record.

## Variants

Both run the same bush: `bush_thermal_model = CONDUCTION_2D`,
`bush_thermal_conductivity = 250` W/(m K), `bush_outer_radius = 0.10` m,
`bush_radial_cells = 16`, and a **globally adiabatic shaft**
(`shaft_ambient_heat_transfer = 0`, Ferron's assumption and the one variant A of
the LOCAL case was run under). They differ only in the outer surface:

| variant | `bush_outer_boundary` | `bush_outer_heat_transfer` |
|---|---|---|
| `outer_fixed313` | `FIXED_TEMPERATURE` | ignored (carried as 0) |
| `outer_robin80` | `CONVECTION` | 80 W/(m^2 K) |

with `thermal_ambient_temperature = 313.15` K (40 C) as the datum in both.

## Stated approximations

- The bush equation is **z-collapsed**, as Ferron's is: it consumes the
  z-averaged film-side flux at each theta and returns one bore temperature per
  theta. The axial end faces of the bush are therefore not modelled; at the
  paper's h = 80 W/(m^2 K) over the two end rings that is of order 1 % of the
  heat, and it is stated rather than absorbed into a coefficient.
- The shaft stays circumferentially isothermal (globally adiabatic here), which
  is also Ferron's assumption.
- Everything else - grid (720 x 192), force-balance protocol and tolerances,
  the two-pass groove closure, the Ferron polynomial viscosity, FBNS cavitation,
  the cross-film order-4 collocation and its two documented approximations -
  is frozen from case 01-THD-crossfilm.

## Everything else is frozen from case 01-THD-crossfilm

Same geometry, same operating point, same 720 x 192 grid, same damped-Newton
force balance to 0.1 % in |W| and 0.05 deg in force direction, same PINNED ->
FLOW_WEIGHTED two-pass groove closure with Q_s taken from the first pass's own
steady mass balance. Every solver call, including the Jacobian probes, is
archived.

## Experiment mapping

`cases/ref_data/ferron1983_fig4_2000rpm_4kN_temperature.csv` and
`..._pressure.csv`. The CSV abscissa is degrees from the supply-groove centre in
the direction of rotation; the solver groove centre is at theta = 90 deg and
rotation is +theta, so `solver_theta = (90 + x) mod 360`. Both are plotted
against `(theta - 90) mod 360`. No fitted shift.

The model quantity compared to the thermocouples is `T_bearing`
(`temperature_bearing`), the bush-side film surface temperature - the closest
available quantity to Ferron's sensors 0.5 mm below the bush inner surface.
Peaks are quoted with a +-12 deg band around the groove removed, because the
groove pins the film temperature over a single cell and the conjugate closure
reacts to that step with a one-cell spike. Same convention as the LOCAL case.

## Relation to the later production plot (traceability)

This case - the **verdict** campaign - uses the paper's own outer coefficient,
**h = 80 W/(m^2 K)**, for the Robin variant. The later production plot batch
(`rerun_cht_robin10`) uses a user-selected **h = 10 W/(m^2 K)** instead. Both
must stay quotable, so the difference was measured directly rather than assumed:

| outer BC | G referred to the bore | T_bearing peak | peak location |
|---|---|---|---|
| h = 80 (this case) | 156.5 W/(m^2 K) | 45.38 C | 266.3 deg |
| h = 10 (production) | 19.9 W/(m^2 K) | 45.58 C | 265.8 deg |

An 7.8x weaker outer sink moves the bush surface by **+0.2 K and -0.5 deg**. The
reason is that the outer path is not what sets the bore temperature: the
film-side coupling is B ~ 5200 W/(m^2 K), so G is 3.0 % of the local balance at
h = 80 and 0.38 % at h = 10, and the circumferential fin length that sets the
shape, sqrt(k_b t/(B+G)), changes from 55.4 deg to 56.1 deg of arc. The shape is
a k_bush effect, not an h_out effect.

(Measured at the converged `outer_fixed313` geometry, eps = 0.593400,
psi = -37.5036, PINNED, 720 x 192, 24 ranks on afml-r6 - a single diagnostic
point, not a force-balanced run.)

**CORRECTION, from the force-balanced production run.** The diagnostic above was
run at *frozen geometry*, and on that basis it was stated here that h = 10 would
not close the ~1.3 K cold bias of the h = 80 result. **That was wrong**, because
holding eps fixed is not the comparison that matters: re-closing the force
balance at h = 10 moves the equilibrium from eps = 0.5934 to **eps = 0.6033**,
and the extra eccentricity adds most of the missing temperature. The production
group `rerun_cht_robin10` (see `cases/01_ferron_wallbc_fields/`) lands at

| | peak | location | mean err | rms |
|---|---|---|---|---|
| h = 80, near-eq (this case) | 45.38 C | 266.3 deg | -1.31 K | 1.76 K |
| h = 10, force-balanced | 46.71 C | 265.2 deg | **+0.02 K** | **1.17 K** |

so the cold bias does close, to a mean error of +0.02 K over the 16 points. The
lesson is about the diagnostic, not the physics: a frozen-geometry probe cannot
be compared against a re-closed equilibrium.

## Results

`SUMMARY.md` in this directory. Short version: `outer_fixed313` converged and
moved the bush-temperature peak from 347.7 deg (local closure, adiabatic) to
247.8 deg against a measured 258.2 deg, at a level about 3 K too cold;
`outer_robin80` was aborted by the rank-dependent linear-solve divergence
(reproducer and fix: `cases/90_diagnostics/ksp_rank_robustness/`, commit
`f0e772d`) and has not been re-run since.

## Remote-only archives for this case

On afml-r6 under `~/pancake/cases/01_ferron_conjugate/`:

- `r6_24rank_rerun/` - the 24-rank reproduction (folded in from `~/cases_run/`
  on 2026-08-05). `outer_fixed313` matches the 4-rank result to every printed
  digit; `outer_robin80` hit the same KSP abort at `fb_iter_03`.
- `diag_h10_frozen_geometry/` - the frozen-geometry h = 10 probe whose conclusion
  is corrected above.
- `gate_fixed313/` - the cross-platform gate run.
