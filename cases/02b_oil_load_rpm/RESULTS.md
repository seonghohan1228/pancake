# Oil load and rpm scan results

**Superseded.** These results used the wrong geometry for the user's rig.
Use [the corrected 43 mm / 40 mm / 42.5 micrometre scan](../02c_oil_corrected_geometry/RESULTS.md).
The numerical checks below apply only to the archived assumed geometry.

The remote scan completed 25 force-balanced steady points and six selected
mesh repeats, plus a 100 rpm / 1 kN probe. All 25 matrix points cavitate at the
assumed 500 Pa absolute floor. This operating window compares established
cavity extent and filling deficit; it does not supply a noncavitating control.

The [setup and model assumptions](README.md) are part of these results. The
geometry is the inherited 49.3 mm diameter, 50 mm length and 350 micrometre
radial clearance. The assumed axial groove is 10 degrees wide at 90 degrees.
Both ends and the groove are at 1 bar absolute. The liquid has constant
viscosity 0.210 Pa s and density 937 kg/m3.

## Recommended comparison cuts

1. At **1000 rpm**, compare **1, 2 and 3 kN**. The cavity responds strongly over
   this range. At 5 kN, pressure nearly doubles from the 3 kN value while cavity
   arc grows only another 4.5 degrees, so 5 kN is a weak additional optical
   contrast for the greater loading.
2. At **3 kN**, compare **500, 1000 and 2000 rpm**. Increasing speed broadens the
   cavity while reducing its maximum filling deficit and increasing the minimum
   film. The 1500 rpm point adds little mid-span arc relative to 2000 rpm.

These cuts share 3 kN / 1000 rpm, requiring **five distinct comparison points**.
They are candidate conditions within the model. They are not certified rig
operating limits. In particular, the 3 kN / 500 rpm film is only about 15
micrometres, and its peak pressure is about 69 bar. The PMMA wall and thermal
response need checking before transferring these loads to the experiment.

| Load (kN) | Speed (rpm) | Mid-span cavity arc (deg) | Cavity area (%) | Minimum film (micrometres) | Peak absolute pressure (bar) |
|---:|---:|---:|---:|---:|---:|
| 1 | 1000 | 90.0 | 18.62 | 69.4 | 13.3 |
| 2 | 1000 | 103.5 | 22.24 | 39.2 | 32.1 |
| 3 | 1000 | 109.5 | 23.71 | 27.8 | 54.2 |
| 3 | 500 | 91.5 | 20.20 | 15.0 | 69.0 |
| 3 | 2000 | 114.0 | 25.78 | 50.0 | 43.4 |
| 5 | 1000 | 114.0 | 24.96 | 17.8 | 106.0 |
| 5 | 250 | 75.0 | 17.20 | 5.0 | 177.7 |

Values in the table and mosaics come from the common 240 by 64 grid. Refinement
at 3 kN changes the 500/1000/2000 rpm arcs to 90.75/108/112.5 degrees.

## What lower speed adds

At 1 kN, 250 rpm gives a 66-degree arc and 20.7 micrometre minimum film. It is
a useful optional extension if a smaller cavity is needed. At 100 rpm, the arc
falls to 51 degrees but minimum film falls to 8.9 micrometres. This probe still
does not remove cavitation. At 5 kN / 250 rpm, the 480 by 128 check gives 4.95
micrometres and 178.4 bar, so this is an extreme of the rigid calculation and
not a sensible first experimental comparison.

At 1 kN the speed trend is non-monotonic: arcs at 500/1000/1500/2000 rpm are
81/90/88.5/87 degrees. Load balance changes eccentricity with speed, so more rpm
does not guarantee a longer cavity. At 3 kN the calculated friction power rises
from 36 W at 500 rpm to 295 W at 2000 rpm. The scan holds viscosity fixed and
does not predict the resulting steady oil temperature.

## Numerical checks

- The existing executable's FBNS, inlet and force tests all passed before the
  sweep. No source changes or package installations were required.
- All 25 coarse equilibria pass vector load error below 0.2 percent, FBNS
  residual below 1e-7, finite fields, positive gap and the prescribed pressure
  and filling-fraction bounds. Largest load error is 0.146 percent.
- The largest relative mismatch of net end flux and groove supply is
  1.64e-7. The maximum coarse coupling pressure residual passes 1e-6.
- Six fresh force balances at 480 by 128 change cavity arc by at most 1.5
  degrees, cavity area by at most 0.388 percentage points, peak pressure by less
  than 0.4 percent, and minimum gap by less than 1 percent. These checks support
  ranking the tested conditions; they do not establish a converged onset curve.

## Contours and evidence

The housing angle is displayed from -180 to 180 degrees to keep the cavity
continuous across the plot. It is not re-centred independently for each case.
Grey marks the groove; orange marks filling fraction 0.999. Each mosaic uses
one shared color scale. The pressure plot uses a square-root color normalization
to show low-pressure structure alongside the high-load peaks; tick values
remain physical absolute pressures in bar.

![Load and speed maps](load_rpm_map.png)

![Cavity filling deficit](cavity_contours.png)

![Absolute pressure contours](pressure_contours.png)

- `sweep_summary.csv` and `sweep_summary.json`: all 25 equilibrium metrics,
  Sommerfeld number and scaled cavitation margin.
- `validation.json`: point-by-point checks and six grid comparisons.
- `provenance.json`: host, executable SHA-256, source commit, settings and trial
  count. The compact archive includes the final configurations, logs, diagnostics
  and field planes. Every trial's raw VTK output stays on `afml-r6` under
  `/home/user/pancake_oil_scan_20260911/runs`.

The two independent nondimensional axes are $S\propto N/W$ and
$\Pi_c\propto1/N$, with exact definitions in the [case description](README.md).
The 500 Pa floor is an assumed equilibrium limit, not a measured vapor-pressure
law or a kinetic vaporization calculation.

Next step: confirm the current rig dimensions and groove geometry against this
case before using the five-point cross as the experimental sweep.
