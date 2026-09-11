# Corrected oil sweep results

This is the replacement scan for **43 mm diameter, 40 mm length and 42.5
micrometre radial clearance**. The former 49.3 mm / 50 mm / 350 micrometre
calculation and its recommended cuts are superseded.

All 25 main points are newly force-balanced steady solutions using viscosity
0.210 Pa s, density 937 kg/m3, and 1 bar absolute pressure at the ends and groove.
Six selected points were repeated at 480 by 128 and two 100 rpm probes were
added. The full [setup and retained assumptions](README.md) are part of these
results. No old equilibrium field is scaled or reused.

## Two useful sweep axes

| Sweep | Hold fixed | Conditions | Maximum filling deficit, 1 - vartheta |
|---|---|---|---|
| Load | 500 rpm | 1, 3, 5 kN | 15.1, 46.6, 62.9 percent |
| Speed | 3 kN | 250, 500, 1000, 2000 rpm | 66.8, 46.6, 27.1, 14.6 percent |

The shared 3 kN / 500 rpm point gives **six distinct conditions**. If only
500--2000 rpm is convenient, omit the 250 rpm extension and use five points.
The 1500 rpm cases lie between the 1000 and 2000 rpm values and do not add a
separate regime for a first comparison.

**More rpm gives a broader, more liquid-filled cavity at fixed load.** The
maximum filling deficit changes much more clearly than the arc length. At
500 rpm, load strongly increases the deficit even when cavity area is nearly
unchanged. This is why the new recommendation uses 500 rpm for the load axis,
and does not carry forward the first scan's 1000 rpm recommendation.

| Load (kN) | Speed (rpm) | Cavity area (%) | Summed mid-span cavity arc (deg) | Maximum deficit (%) | Minimum film (micrometres) | Peak absolute pressure (bar) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 500 | 20.30 | 115.5 | 15.1 | 35.5 | 10.6 |
| 3 | 500 | 24.65 | 108.0 | 46.6 | 24.6 | 36.6 |
| 5 | 500 | 25.04 | 96.0 | 62.9 | 18.7 | 68.6 |
| 3 | 250 | 24.36 | 97.5 | 66.8 | 16.9 | 43.0 |
| 3 | 1000 | 30.35 | 133.5 | 27.1 | 31.7 | 33.4 |
| 3 | 2000 | 33.10 | 144.0 | 14.6 | 36.5 | 32.8 |

Rounded values are from the common 240 by 64 matrix. The groove separates
the upstream and downstream cavity patches in several cases. The reported arc
sums all cavitated mid-span cells; it is not necessarily one connected arc.
Area and arc use the threshold vartheta < 0.999, and should not be read as
measured optical visibility or an onset criterion.

## Lower speed and thermal interpretation

At 100 rpm / 1 kN, maximum filling deficit is 58.1 percent, cavity arc 87 degrees
and minimum film 19.2 micrometres. At 100 rpm / 3 kN, the deficit reaches 83.7
percent, but minimum film falls to 9.1 micrometres. That condition offers a
stronger deficit at a much smaller gap; 250 rpm is the more moderate extension
of the main matrix. Both low-speed probes still cavitate under the assumed
500 Pa floor. No noncavitating control was found in the scanned conditions.

At 3 kN, calculated friction power is 11/39/141/548 W at
250/500/1000/2000 rpm. Thus 2000 rpm is a useful weak-deficit comparison, but
holding viscosity at 0.210 Pa s requires temperature control or a subsequent
thermal calculation. These runs are steady pressure solutions, not predictions
of steady oil temperature. The former large-clearance minimum-film and pressure
limits do not apply to this geometry. No PMMA material or wall design is
assumed solely from the old rig dimensions.

## Numerical checks

- All 25 equilibrium configurations were checked against R = 0.0215 m,
  L = 0.040 m, c = 42.5e-6 m, the requested oil properties and both end pressures.
  All match. Their largest vector force error is 0.155 percent.
- Every main point passes the FBNS residual threshold 1e-7, finite physical
  fields, positive film and the prescribed pressure/filling bounds. Largest
  relative end-flux versus groove-supply imbalance is 7.16e-9.
- Six 480 by 128 repeats change maximum filling deficit by at most **0.240
  percentage points**. The contrast used to select the cuts is robust at this
  screening level.
- Arc changes are larger, up to **6 degrees**, and cavity area changes by up to
  **1.348 percentage points**. These are not grid-independent cavity boundaries.
  Cell-centre groove masking also changes the represented width between meshes;
  a grid aligned to the groove edges is appropriate for later front-location
  work. Use the present sweep to rank conditions by deficit, not to validate
  exact rupture/reformation angles.

## Requested temperature plots

An additional 18 steady energy-equation calculations cover the nine displayed
load/speed points under two explicit cooling assumptions. Both use provisional
20 C inlet oil, cp = 2000 J/(kg K), conductivity = 0.13 W/(m K) and constant
0.210 Pa s viscosity. The cooled branch uses 20 C walls with H = 10000
W/(m2 K) per wall; the adiabatic branch sets both wall coefficients to zero.
The groove is an ideal 20 C thermal reservoir. These are sensitivity scenarios,
not measured thermal conditions or a self-consistent viscosity-coupled solution.

| Speed (rpm), at 3 kN | Cooled-wall peak (C) | Adiabatic-wall peak (C) |
|---:|---:|---:|
| 500 | 20.6 | 45.9 |
| 1000 | 21.7 | 61.6 |
| 2000 | 25.8 | 96.5 |

The temperature contours are film-averaged liquid temperatures. STRIATED
weights thermal mass, convection, conduction and shear by liquid filling.
It does not compute a separate vapor temperature. All thermal runs converged
below 1e-6 in the pressure and temperature coupling checks, retained the
0.2-percent force tolerance, and preserved the requested constant viscosity.
The existing solver energy-equation test also passed.

At 3 kN / 2000 rpm, repeating both thermal scenarios on the 480 by 128
force-balanced field changes the cooled peak by 0.001 K and the adiabatic peak
by -0.249 K. `temperature_grid_check.json` records this selected check.

The comparison image below holds W = 3 kN. It plots temperature rise relative
to 20 C, with **different color ranges for the two rows**, explicitly labelled
on the colorbars. The separate full mosaics use one scale for each cooling
scenario. The high adiabatic temperatures show why actual cooling and mu(T)
are needed before treating a thermal value as a rig prediction.

![Temperature rise at 3 kN for two cooling assumptions](temperature_comparison_contours.png)

[All nine cooled-wall temperature contours](temperature_cooled_contours.png)
and [all nine adiabatic-wall temperature contours](temperature_adiabatic_contours.png)
are also saved. `temperature_summary.json` contains every thermal result and
its input assumptions.

## Pressure, cavity and archived evidence

Grey marks the prescribed groove and orange the vartheta = 0.999 contour.
Fields use the same housing-angle coordinate and color bounds across panels.
The cavity color scale spans zero to 0.7 for the selected mosaic, rather than
the previous case's zero-to-one scale. Pressure is absolute, with a square-root
color normalization to show structure at both light and heavy loads. The
parameter maps interpolate within the sampled matrix; black dots are actual
calculations.

![Load and rpm maps](load_rpm_map.png)

![Corrected cavity contours](cavity_contours.png)

![Corrected pressure contours](pressure_contours.png)

`sweep_summary.csv`, `sweep_summary.json`, `validation.json` and
`oil_scan_compact.tar.gz` contain the numerical results and compact fields.
`provenance.json` identifies the executable, source commit and every fixed input.
All raw trial VTK files remain on `afml-r6` in
`/home/user/pancake_oil_corrected_20260911/runs`.

For a dimensionless plot, the two independent axes remain the Sommerfeld
number S and scaled pressure margin Pi_c, now recomputed with the corrected
geometry. Their exact definitions are in the [case description](README.md).

The retained 500 Pa floor and 10-degree unloaded-side groove are assumptions.
The filling deficit is an equilibrium cavity proxy and does not establish that
an experimental bubble contains pure oil vapor.
