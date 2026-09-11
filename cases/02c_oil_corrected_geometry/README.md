# Corrected-geometry oil scan, 2026-09-11

This case supersedes the geometry and operating recommendations in
[case 02b](../02b_oil_load_rpm/README.md). The user confirmed a **43 mm journal
diameter, 40 mm bearing length and 42.5 micrometre radial clearance** after the
first scan. Every trial here is a fresh calculation with that geometry.

## Inputs

| Input | Value | Status |
|---|---:|---|
| Journal radius R | 21.5 mm | User diameter divided by two |
| Bearing length L | 40 mm | User supplied |
| Radial clearance c | 42.5 micrometres | User supplied, not diametral |
| Dynamic viscosity | 0.210 Pa s | User supplied, constant |
| Density | 937 kg/m3 | User supplied, constant |
| Both end pressures | 100000 Pa absolute | 1 bar, absolute interpretation retained |
| Axial groove pressure | 100000 Pa absolute | Ideal end-reservoir connection |
| Groove centre / width | 90 / 10 degrees | Retained assumption, unloaded side |
| Cavitation floor | 500 Pa absolute | Retained assumption, not measured oil vapor pressure |

The groove is full length, with pressure prescribed everywhere. Pressure loss
along an actual feed passage is not solved. The imposed external load points
at +90 degrees and the target fluid pressure force points at -90 degrees in
the fixed housing frame. The shaft rotates in the positive theta direction.

The calculation uses the existing isothermal, constant-property FBNS solver
with SCALAR_PCAV, no dissolved gas, and a static bearing position. For each
applied load, an external damped Newton search adjusts eccentricity and attitude
until the vector force error is below 0.2 percent. These steady nonlinear
solutions do not integrate a physical transient or a steady thermal balance.

## Execution and evidence

Remote host: `afml-r6`, accessed through `user@143.248.99.6`.
Run directory: `/home/user/pancake_oil_corrected_20260911`.
Executable: `/home/user/pancake_fig6_446327c/build-linux/pancake`.
Source commit: `446327ca4643e8806104a167ae2075759e38b4e6`.
Executable SHA-256:
`d2431b788661af31fd29b4e07cfac4f4d02b6cf66171efd0b6ae81447977ec64`.

The source and binary are unchanged from the checked first scan. Separate
case and output directories prevent cached equilibria from the previous
geometry being used. Configurations, residual logs and raw VTK files for every
trial stay on the workstation. The compact archive contains final field planes,
configurations, diagnostics, logs, metrics and the executable provenance.

```sh
python3 run_scan.py pilot
python3 run_scan.py scan
python3 run_scan.py refine
python3 run_scan.py low
python3 validate_scan.py
python3 plot_scan.py
python3 run_temperature.py
python3 run_temperature.py refine
python3 plot_temperature.py
python3 bundle_results.py
```

The main matrix contains 1/2/3/4/5 kN and 250/500/1000/1500/2000 rpm on a
240 by 64 grid. Two ranks serve each solver call; five independent speed
continuations can run concurrently. Selected 480 by 128 force-balanced
repeats check mesh sensitivity. Additional 100 rpm probes assess whether
lower speed offers more useful contrast.

`validate_scan.py` audits FBNS residual, positive film, finite fields,
pressure/filling bounds, force closure and net end-flux versus groove supply.
`plot_scan.py` reads final compact arrays using their physical theta/z
coordinates. Plots show the same housing frame, from -180 to 180 degrees, for
every case. Grey marks the groove and orange the filling-fraction 0.999 contour.
The output labelled `1 - vartheta` is the depth-averaged liquid-volume deficit.
The peak deficit and cavity area are different observables: a broad cavity may
still contain mostly liquid. The load/rpm map therefore shows area, maximum
deficit and minimum film thickness.

## Two physical axes and limits

The adjustable axes are load W and rpm N. For fixed geometry and liquid,
independent dimensionless coordinates follow from the Reynolds pressure scale:

\[
P_\mu=\mu\omega(R/c)^2,\quad \omega=2\pi N/60,\quad
S=\frac{\mu(N/60)}{W/(2RL)}(R/c)^2,\quad
\Pi_c=\frac{p_{\rm end}-p_c}{P_\mu}.
\]

Thus S is proportional to N/W and Pi_c to 1/N. The numerical values are
recomputed using this case's confirmed dimensions and included in the summary.
Density multiplies mass flux but supplies no inertial term in the constant-
property Reynolds model used here.

The prescribed 500 Pa floor is an equilibrium cavitation proxy. Neither these
oil properties nor the pressure/filling contours identify actual vapor pressure,
bubble nucleation, dissolved gas release or phase-change kinetics. Viscosity
is held fixed; thermal changes, structural deformation, roughness and contact
are absent. The former 50 mm PMMA-bore geometry is not imported as a confirmed
material or structural design for the corrected bearing.

## Change record

- 2026-09-11: Correct all three geometry inputs explicitly requested by the
  user. Preserve the first scan as superseded evidence and rerun in isolation.
- Reuse the prior driver and checks, but initialize the much lower expected
  eccentricity for the smaller clearance. All accepted forces still come from
  Newton-balanced solver fields, not from that initial estimate.
- Recompute contour ranges and add maximum filling deficit to the parameter
  map because it carries contrast that cavity arc alone misses.

See [RESULTS.md](RESULTS.md) for the corrected sweep recommendation and values.

## Requested temperature plots

After the isothermal pressure sweep, the user requested temperature plots.
`run_temperature.py` therefore solves the existing steady film-averaged energy
equation at nine displayed equilibria under two provisional cooling assumptions.
The inlet is 20 C, cp is 2000 J/(kg K), conductivity is 0.13 W/(m K), and both
walls are either adiabatic or cooled toward 20 C with H = 10000 W/(m2 K) on each
wall. These values are assumed because the user has not supplied thermal data.
The corresponding volumetric heat capacity is rho cp = 1.874e6 J/(m3 K).

STRIATED weights cavity liquid heat capacity, convection, conductivity and shear
by filling fraction. The temperature is that of the liquid, not a separately
resolved vapor phase. The groove pins entering/recirculating oil to 20 C in
this ideal-reservoir closure; actual inlet mixing or feed resistance is not
identified. Viscosity and density stay fixed, so pressure and force-balanced
positions do not respond to heating. This is a temperature-sensitivity
calculation with frozen properties, not a self-consistent THD prediction.

Each thermal solve checks finite temperatures, residuals below 1e-6, the original
force tolerance and unchanged viscosity. Two-dimensional temperature planes and
the provisional thermal inputs are archived alongside each final result.
