# Oil load and speed screening, 2026-09-11

**Superseded for the user's rig.** The user corrected the journal diameter to
43 mm, bearing length to 40 mm and radial clearance to 42.5 micrometres.
[Case 02c](../02c_oil_corrected_geometry/README.md) reruns the scan with those
dimensions. This case remains an archive of the original assumed geometry;
its contour values and recommended sweep must not be used for the corrected rig.

This is a steady, isothermal operating-window scan for the optical journal-bearing
rig. It uses the user's oil properties, 0.210 Pa s and 937 kg/m3, and equal 100000 Pa
absolute pressure at both ends and the full-length axial groove.

## Assumptions and scope

- Geometry inherited from case 02 Rev B: radius 24.65 mm, length 50 mm, radial
  clearance 350 micrometres. The user did not provide replacement dimensions.
- A 10-degree groove is centred at housing angle 90 degrees. The supported
  external load points at +90 degrees, and the target pressure force is -90
  degrees. Thus the groove is on the unloaded side. Groove location and width
  remain assumptions pending confirmation of the rig drawing.
- The groove pressure is prescribed everywhere along its length. This represents
  an ideal connection to the end reservoirs and does not calculate pressure loss
  along a physical groove or inlet passage.
- FBNS with SCALAR_PCAV and no dissolved gas uses a prescribed 500 Pa absolute
  cavitation floor, inherited from the previous rig screening. This is not a
  measured vapor pressure for the newly supplied oil. Viscosity and density alone
  do not determine its vapor pressure.
- The liquid is isothermal with constant properties. No thermal viscosity change,
  shaft or PMMA deformation, surface roughness, contact, vaporization kinetics,
  or gas release is included. The cavity variable is a depth-averaged liquid
  filling fraction, not a resolved bubble population.

These assumptions make the scan useful for ranking comparison points. They do
not establish that a cavity seen experimentally consists of oil vapor. The
mass-conserving filling-fraction interpretation is consistent with the general
[Elrod-Adams model description](https://doc.comsol.com/6.4/doc/com.comsol.help.models.cfd.journal_bearing_cavitation/journal_bearing_cavitation.html);
the actual implementation used here is pancake FBNS.

## Reproduction

The executable is the already-built remote solver at
`/home/user/pancake_fig6_446327c/build-linux/pancake`, source commit
`446327ca4643e8806104a167ae2075759e38b4e6`. Its SHA-256 is
`d2431b788661af31fd29b4e07cfac4f4d02b6cf66171efd0b6ae81447977ec64`.
The source tree has no modified files under `src/`. No solver changes were made.

Run directory on `afml-r6`: `/home/user/pancake_oil_scan_20260911`.
`settings.json` contains all sweep and geometry parameters. Each trial position
has a hash-keyed directory with its exact configuration, solver log, diagnostics,
raw VTK files, compact field planes and metrics. Successful equilibrium points
are recorded separately so a resumed sweep skips completed work.

```sh
python3 run_scan.py pilot
python3 run_scan.py scan
python3 run_scan.py refine
python3 validate_scan.py
python3 plot_scan.py
```

The main scan is 25 points, 1/2/3/4/5 kN at 250/500/1000/1500/2000 rpm on a
240 by 64 grid. Each point solves both eccentricity and attitude by damped
Newton iteration, with a finite-difference force Jacobian. The vector force
error must be below 0.2 percent of target load. Each inner solution must pass
the FBNS and outer-coupling residual checks. Continuation uses the preceding
load at the same speed; speeds are independent. The MPI rank count is two.

The `refine` mode repeats selected points at 480 by 128 with a fresh force
balance. `validate_scan.py` verifies positive films, finite pressure and filling
fraction, cavitation floor, force closure, FBNS residual and the steady balance
of end flux against groove supply. This is a coarse screening with selected
refinement checks, not a complete grid-independence study.

Plots use housing angle, without shifting each case to its own line of centres.
The grey strip marks the groove and the orange contour marks filling fraction
0.999. Cavity area and mid-span arc use that same threshold. Void fraction is
`max(0, 1 - film_content)`. All panels of a field mosaic share one color scale.
The load-speed maps interpolate within the sampled matrix; black dots identify
the actual calculations. No onset boundary is inferred outside those samples.

## Two dimensionless axes

For fixed geometry, define mean load pressure and the viscous pressure scale as

\[
\bar p=W/(2RL),\qquad P_\mu=\mu\omega(R/c)^2,\qquad \omega=2\pi N/60.
\]

Two useful independent coordinates are the Sommerfeld number and the scaled
pressure margin to the cavitation floor:

\[
S=\frac{\mu(N/60)}{\bar p}(R/c)^2,
\qquad \Pi_c=\frac{p_{\rm end}-p_c}{P_\mu}.
\]

These follow by scaling the steady Reynolds flux and load integral. At fixed
oil and geometry, `S` is proportional to `N/W`, while `Pi_c` is proportional
to `1/N`. A constant `N/W` comparison therefore does not keep the cavitation
pressure margin fixed. Density is used in the mass flux but does not provide
an independent inertial control in this constant-property Reynolds model.

The directly adjustable plotting axes remain load and rpm. The dimensionless
coordinates are included in the summary data to distinguish loading from the
available pressure margin without claiming an empirical data collapse.

## Change record

- 2026-09-11: Create an isolated task worktree because the original worktree has
  unrelated Ferron and optical-rig changes. Reuse the existing checked remote
  executable instead of rebuilding or modifying the solver for a quick scan.
- Replace the earlier rig oil properties and add a pressure-pinned axial groove.
  Close the vector load at every point rather than prescribing eccentricity.
- Save raw evidence remotely, compact results and 300 dpi contours locally.
  Generated outputs are ignored by Git so this change adds only the reproducible
  driver, plotting, numerical audit and documentation.
