# Journal-bearing analytical benchmarks

Verification of the isothermal, rigid, prescribed-eccentricity Reynolds solve
against the two exact limits that bracket any finite bearing.

```bash
python run_benchmarks.py ../../build-windows-mingw/pancake.exe
```

`analytic.py` holds the closed forms; `run_benchmarks.py` generates configs, runs
`pancake.exe`, and reports the error normalised by the pressure amplitude.

## Results (2026-07-31, 720 x 40/60 mesh, eps = 0.6, 1000 rpm)

Case A -- infinitely-long (Sommerfeld) limit, `FULL_SOMMERFELD`, mid-plane
p(theta). Residual is the axial end effect and must vanish as L/D grows.

| L/D | max err | rms err |
|-----|---------|---------|
| 1   | 22.86%  | 14.94%  |
| 2   | 7.32%   | 4.94%   |
| 4   | 0.76%   | 0.51%   |
| 8   | 0.10%   | 0.04%   |

At L/D = 8 the peak pressure is 14.5944 bar against an exact 14.5923 bar (0.014%).

Case B -- short-bearing (Ocvirk) limit, `GUMBEL`, full p(theta, z). Residual is
the neglected circumferential pressure flow, which is O((L/D)^2).

| L/D    | max err | rms err | error ratio |
|--------|---------|---------|-------------|
| 0.5    | 37.12%  | 7.87%   | --          |
| 0.25   | 14.63%  | 2.95%   | 2.54        |
| 0.125  | 4.40%   | 0.86%   | 3.32        |
| 0.0625 | 1.17%   | 0.22%   | 3.77        |

The error ratio approaching 4 under successive halving of L/D confirms the
expected second-order convergence to the Ocvirk limit.

## Scope

These are verification against exact solutions of the same PDE the solver
discretizes. They bound discretization error in both operator terms; they say
nothing about whether the Reynolds equation models a given experiment, and they
do not exercise the cavitation or gas models (Case A has no cavitation, Case B
uses Gumbel clipping).

---

# Cavitation and energy-equation validation (2026-07-31)

The benchmarks above validate the Reynolds discretisation only. This section
covers the mass-conserving Elrod-Adams (JFO) free boundary and the energy
equation, which they do not touch.

```bash
python run_raimondi_boyd.py ../../build-windows-mingw/pancake.exe   # Target A
python run_ferron.py        ../../build-windows-mingw/pancake.exe   # Target B
python plot_validation.py                                          # figures
python make_tables.py                                              # these tables
```

| module | contents |
|---|---|
| `raimondi_boyd.py` | the Raimondi & Boyd chart data, its provenance audit, and the exact infinite-length Swift-Stieber solution |
| `ferron.py` | Ferron's rig geometry, oil properties and digitised measurements, with the source of each number |
| `run_raimondi_boyd.py` | Target A sweep and the free-boundary resolution audit |
| `run_ferron.py` | Target B thermohydrodynamic comparison |
| `plot_validation.py` | the figures |
| `make_tables.py` | generates the result tables below from the result JSONs |

Run outputs land in `_work/` (gitignored); the committed artefacts are the result
JSONs, the figures and this file.

## Reference data and its provenance

`RB_TABLE` reproduces "Table (1): Dimensionless performance parameters for full
journal bearing" from p.192 of the University of Babylon lecture notes *Sliding
Bearings*, the standard textbook reprint of Raimondi & Boyd, ASLE Trans. **1**,
159-209 (1958). Textbook reprints are known to disagree in the third digit, so
the transcription is not taken on trust: `swift_stieber_infinite` solves the
infinitely-long bearing with the Reynolds rupture condition to machine precision,
and it reproduces the table's L/D = inf block to **0.6% in S and 0.4 deg in the
attitude angle** across eps = 0.1-0.9. A second, independent consistency check
comes from Al-Hamood et al. (2015), SSRG IJME **2**(10), 17-20, whose Table (3)
reports eps = 0.28 and phi = 70 deg from the charts at conditions that give
S = 1.39; log-interpolating the L/D = 1/2 block at S = 1.39 returns eps = 0.28
and phi = 69.5 deg.

Ferron's parameters come from the open-access reproduction Alakhramsing, van
Ostayen & Eling, *Lubricants* **3**(2), 256-280 (2015), Table 1, retrieved from
the TU Delft institutional repository. **Note that this contradicts the commonly
circulated secondary values**: the rig is L/D = 0.80 (L = 80 mm) with
C = 145 um (C/R = 2.9e-3), not L/D = 0.55 with C/R = 1.5e-3.

## The free boundary has a resolution requirement, and it is not obvious

Elrod carries pressure inside the film content, `p = p_cav + beta ln(theta)`, so
the entire full-film pressure signal lives in `theta - 1 = p/beta`. How well the
mesh resolves that excursion controls the accuracy of the rupture boundary, and
the accuracy degrades sharply as `p_max/beta` shrinks -- that is, as the film
becomes *more* incompressible, which is the physically correct limit.

Measured at L/D = 8, eps = 0.6 against the exact Swift-Stieber rupture angle of
213.08 deg (`run_raimondi_boyd.convergence_study`):

| beta [Pa] | n_theta | rupture [deg] | error [deg] | S vs exact |
|---|---|---|---|---|
| 1e9 | 720  | 192.25 | -20.83 | -19.7% |
| 1e9 | 2880 | 206.06 | -7.02  | -8.2%  |
| 1e8 | 720  | 211.25 | -1.83  | +2.1%  |
| 1e8 | 1440 | 212.62 | -0.46  | +6.6%  |
| 1e8 | 2880 | 213.31 | +0.23  | +9.0%  |
| 1e8 | 5760 | 213.59 | +0.51  | +10.2% |

The rupture angle converges cleanly to the exact value, so the cavitation
physics is right. But **the beta = 1e9 / 720 combination is wrong by 21 degrees
of arc while reporting itself converged** -- `mc = 0`, `ns = 0`, `R_p = 0`, and a
zero exit code. The residual +10% in S at L/D = 8 is end leakage, which does not
vanish at that aspect ratio for a cavitated film.

Production runs therefore use `bulk_modulus = 1e8` on a 2880-cell mesh, and every
swept case is repeated at a second viscosity: S and the attitude angle are
invariant under viscosity in the incompressible limit, so any change is a direct
measure of the combined compressibility and under-resolution error.

## Target A -- Raimondi & Boyd, finite-width cavitated film

`cavitation_model = JFO`, `cavitation_threshold = SCALAR_PCAV` at ambient,
prescribed eccentricity, isothermal. Two inlet conditions are run at every point
because they are different physics:

- **fed film** -- a 4 deg full-film supply groove at the maximum-gap position. In
  JFO the groove is a Dirichlet `theta = 1`, which is exactly the "complete film
  at the maximum gap" that Raimondi & Boyd assume. This is the like-for-like
  comparison.
- **ungrooved** -- no supply, the configuration of the optical rig.

### Fed film -- the like-for-like comparison

L/D = 0.5:

| eps | mu [Pa.s] | S sim | S chart | dS | phi sim | phi chart | dphi | rupture [deg] | p_max/beta |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 0.010 | 1.99139 | 2.0300 | -1.9% | 72.57 | 74.94 | -2.37 | 191.94 | 0.0005 |
| 0.4 | 0.010 | 0.76067 | 0.7790 | -2.4% | 59.12 | 61.45 | -2.33 | 193.44 | 0.0015 |
| 0.4 | 0.003 | 0.74528 | 0.7790 | -4.3% | 58.21 | 61.45 | -3.24 | 188.06 | 0.0004 |
| 0.6 | 0.010 | 0.30698 | 0.3190 | -3.8% | 45.78 | 48.14 | -2.36 | 193.94 | 0.0043 |
| 0.8 | 0.010 | 0.08816 | 0.0923 | -4.5% | 31.49 | 33.31 | -1.82 | 193.31 | 0.0203 |
| 0.8 | 0.003 | 0.08465 | 0.0923 | -8.3% | 30.26 | 33.31 | -3.05 | 191.56 | 0.0062 |
| 0.9 | 0.010 | 0.02957 | 0.0313 | -5.5% | 22.11 | 23.66 | -1.55 | 192.19 | 0.0794 |

L/D = 1:

| eps | mu [Pa.s] | S sim | S chart | dS | phi sim | phi chart | dphi | rupture [deg] | p_max/beta |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 0.010 | 0.61711 | 0.6310 | -2.2% | 69.94 | 74.02 | -4.08 | 208.19 | 0.0015 |
| 0.4 | 0.010 | 0.25180 | 0.2640 | -4.6% | 59.12 | 63.10 | -3.98 | 206.44 | 0.0041 |
| 0.4 | 0.003 | 0.23860 | 0.2640 | -9.6% | 55.86 | 63.10 | -7.24 | 201.19 | 0.0013 |
| 0.6 | 0.010 | 0.11640 | 0.1210 | -3.8% | 47.91 | 50.58 | -2.67 | 203.19 | 0.0102 |
| 0.8 | 0.010 | 0.04284 | 0.0446 | -3.9% | 34.52 | 36.24 | -1.72 | 198.19 | 0.0360 |
| 0.8 | 0.003 | 0.04137 | 0.0446 | -7.2% | 33.31 | 36.24 | -2.93 | 196.69 | 0.0110 |
| 0.9 | 0.010 | 0.01771 | 0.0188 | -5.8% | 24.62 | 26.45 | -1.83 | 194.69 | 0.1120 |

Every case converged on the Elrod tolerance with `R_p < 1e-9`; none exhausted its
iteration budget.

At the baseline viscosity, **max |dS| = 5.8% against a 10% criterion and
max |dphi| = 4.08 deg against 5 deg, across eps = 0.2-0.9 at both aspect ratios.**
The errors are one-signed -- the solver carries slightly more load than the chart
everywhere -- and they shrink as `p_max/beta` grows, which identifies them as
residual free-boundary under-resolution rather than a modelling difference.

The mu = 0.003 rows quantify that. Dropping the viscosity by 3.3x drops
`p_max/beta` by the same factor and moves L/D = 1, eps = 0.4 from dS = -4.6% to
-9.6% and dphi from -4.0 to -7.2 deg -- which would have **failed** the attitude
criterion. Roughly half the residual disagreement at the baseline is numerical.

### Ungrooved -- total collapse

| L/D | eps | load [N] | p_max gauge [Pa] | theta_min | S sim | S chart |
|---|---|---|---|---|---|---|
| 0.5 | 0.2 | 0.001 | 74 | 1.000e-06 | 4.99e+04 | 2.0300 |
| 0.5 | 0.4 | 0.002 | 114 | 1.000e-06 | 2.558e+04 | 0.7790 |
| 0.5 | 0.6 | 0.006 | 399 | 1.000e-06 | 1.025e+04 | 0.3190 |
| 0.5 | 0.8 | 0.028 | 2658 | 1.111e-01 | 2219 | 0.0923 |
| 0.5 | 0.9 | 0.144 | 1.603e+04 | 5.264e-02 | 433.4 | 0.0313 |
| 1   | 0.2 | 0.007 | 140 | 1.000e-06 | 1.811e+04 | 0.6310 |
| 1   | 0.4 | 0.018 | 454 | 1.000e-06 | 6860 | 0.2640 |
| 1   | 0.6 | 0.067 | 1593 | 2.500e-01 | 1873 | 0.1210 |
| 1   | 0.8 | 0.311 | 1.059e+04 | 1.111e-01 | 401.7 | 0.0446 |
| 1   | 0.9 | 1.763 | 6.36e+04 | 5.269e-02 | 70.9 | 0.0188 |

Load capacity falls by three to six orders of magnitude. The Sommerfeld number is
1e2-1e5 times the chart value, which is another way of saying there is no bearing
action at all. Two distinct end states appear:

- **eps <= 0.4-0.6:** `theta_min` sits on the solver clamp at 1e-6. The film is
  driven to total loss and the result is the clamp, not physics.
- **eps >= 0.8:** the mid-plane void matches the analytic starved state
  `theta h = const` to rms 6.2e-4 at L/D = 0.5, eps = 0.9 (`theta_min` 0.052673
  against 0.052632 analytic).

Both follow directly from the Elrod switch function being zero in cavitated
cells: nothing can re-flood the film, so a cavity that forms in the divergent
region persists into the convergent region and shortens the wedge to nothing.
This is the same starved branch the optical-rig sweep brief warns about, here
reproduced in a *steady* solve at every eccentricity tested.

## Target B -- Ferron, Frene & Boncompain (1983), thermohydrodynamic

`temperature_model = ENERGY_EQUATION` with Ferron's own oil, not the optical
rig's. Ferron specifies viscosity by an exponential law
`mu = mu_0 exp(-beta (T - T_s))` with `beta` = 3.3e-2 1/K; pancake uses the
Andrade form `mu = mu_ref exp(E_mu (1/T - 1/T_ref))`, so a least-squares refit
over 40-55 C gives **E_mu = 3353 K**, accurate to 0.56% in viscosity across that
range. The optical-rig value of 4803 K is not used anywhere here.

Alakhramsing Figure 9 specifies its operating point by **load** (4 kN at
2000 rpm), not by eccentricity, so each case is run at the eccentricity that
carries 4 kN. The eccentricity that comes out is then itself a prediction,
testable against the measured eps-load curve.

The rig has a single axial groove feeding the unloaded region, so the groove is
placed diametrically opposite the load line; since the load direction is an
output, the groove angle is found by fixed-point iteration. That is the one
convention-dependent input in the comparison, which is why the peak pressure
magnitude and the peak-to-rupture arc -- neither of which depends on it -- are
reported alongside the position, which does.

### Reference point, W = 4 kN at 2000 rpm

| case | eps for 4 kN | p_max [MPa] | vs meas | p_max X/pi | vs meas [deg] | rupture X/pi | peak-to-rupture [deg] | T_max [C] | vs meas [K] | Q [l/min] | max dp per coupling iter [Pa] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| THD, adiabatic | 0.534 | 1.162 | -9.9% | +0.104 | +4.2 | +0.420 | 57.0 | 47.48 | **-1.62** | 8.52 | 159 |
| THD, metal bush | 0.529 | 1.176 | -8.8% | +0.101 | +3.8 | +0.401 | 54.0 | 43.00 | **-6.10** | 7.09 | 59.1 |
| isothermal | 0.517 | 1.168 | -9.5% | +0.100 | +3.5 | +0.405 | 55.0 | 40.00 | -9.10 | 6.93 | 0 |

Measured: p_max = 1.29 MPa at X/pi = +0.08; T_max = 49.1 C at X/pi = +0.55
(a 9.1 K rise); rupture bracketed by X/pi = +0.30 to +0.42, i.e. 40-61 deg past
the peak. Measured side flow is 2.6-4.7 l/min at 1500 rpm and 3.6-8.6 at
4000 rpm (Figure 10b), so 8.5 l/min at 2000 rpm is at the high end.

Against the acceptance criteria: **peak pressure within 15% -- pass** (-9.9%);
**peak angular position within 10 deg -- pass** (+4.2 deg); **rupture location
within 10 deg -- pass**, it lands inside the measured bracket and the
peak-to-rupture arc of 57 deg sits inside the measured 40-61 deg; **bearing-wall
temperature rise within 5 K -- pass adiabatic** (-1.62 K), **fail with a metal
bush** (-6.10 K); **temperature profile shape qualitatively correct -- fail**,
the peak is predicted at X/pi = +0.95 against a measured +0.55.

The last column is the convergence evidence in physical units. `R_p` is
normalised by 1e5 Pa, so the coupling loop has stopped moving the pressure by
more than 159 Pa on a 1.16 MPa peak (0.014%) in the worst case. The configured
`coupling_tolerance` of 1e-8 corresponds to 1e-3 Pa, which the segregated Picard
loop does not reach on the finer mesh within any sane budget and which nothing
here can use.

### Mesh resolution

| mesh | p_max [MPa] | p_max X/pi | rupture X/pi | T_max [C] | load [N] |
|---|---|---|---|---|---|
| 720x40  | 1.1775 | +0.1005 | +0.4005 | 47.491 | 3979 |
| 1440x60 | 1.1617 | +0.1036 | +0.4202 | 47.477 | 3870 |

1.4% in peak pressure, 0.6 deg in peak position, 3.5 deg in rupture position,
0.01 K in temperature. All well inside the acceptance tolerances, and the
measured rupture is itself only bracketed to 22 deg by the spacing of Ferron's
pressure taps.

### Thermal bracket, at fixed eccentricity

Running the three headline cases at the eccentricity that carries 4 kN makes the
isothermal-vs-THD bracket vacuous: matching the load forces the peak pressures
together (1.162 / 1.176 / 1.168 MPa above). The bracket only means anything at a
**common eccentricity**, where the isothermal film keeps the 40 C inlet viscosity
everywhere and must carry more:

| at eps = 0.5345 | p_max [MPa] | W [N] | T_max [C] | mu_min/mu_0 |
|---|---|---|---|---|
| THD | 1.1777 | 3979 | 47.49 | 0.771 |
| isothermal | 1.2751 | 4272 | 40.00 | 1.000 |

The energy equation removes **8.3%** of the peak pressure through a viscosity
knock-down to 0.771 of the inlet value. That is the effect size that justifies
enabling it.

The stated criterion was that isothermal must sit *above* the measurement. It
gives 1.2751 MPa against a measured 1.29 +/- 0.03 MPa, so it falls 1.1% short --
inside the digitising uncertainty, but the strict test is not met.

## Target C -- Cole & Hughes (1956), cavity photographs: NOT SOURCED

Cole, J. A. and Hughes, C. J., "Oil Flow and Film Extent in Complete Journal
Bearings", Proc. IMechE **170**, 499-510 (1956) is paywalled at SAGE, and no
reproduction of its quantitative film-extent data could be found. The sources
checked -- San Andres' cavitation lecture notes, a WPI doctoral thesis
specifically on rupture-point movement, an IntechOpen chapter on cavitation
evolution, Bayada & Chupin (2013), and Unpaywall/Semantic Scholar for repository
copies -- all cite it qualitatively only, for the streamer structure of the
cavitated film. Jakobsson & Floberg (1957) was attempted as a substitute and is
also unavailable (a Chalmers technical report, not online).

**Independently of the data, the reformation half of this target is not testable
with this solver.** There are only two available inlet conditions and neither
produces a physically determined reformation angle: ungrooved, the film never
reforms anywhere (see the collapse below); grooved, reformation is pinned at the
groove by construction, so its location is an input rather than a prediction.
The rupture angle *is* validated, by the exact Swift-Stieber comparison above.

## Verdict 1 -- is the cavitation model validated?

**Yes for a fed film, no for the ungrooved optical-rig configuration.** These are
not shades of the same answer; they are different regimes and only one of them
was validated.

Trustworthy, with a fed film, on a mesh where `p_max/beta >= 4e-3`:

- **Load capacity / Sommerfeld number.** Within 5.8% of the Raimondi & Boyd
  charts over eps = 0.2-0.9 at L/D = 1/2 and 1, one-signed high.
- **Attitude angle.** Within 4.08 deg over the same range.
- **Rupture angle.** Converges to the exact Swift-Stieber value to 0.23 deg at
  2880 cells and 0.51 deg at 5760. This is the single best-validated quantity
  here.
- **Cavity extent measured from the rupture point**, insofar as it follows from
  the rupture angle and mass conservation in the cavity.

Not trustworthy:

- **Anything computed without a supply.** In the ungrooved configuration the
  solver has no load capacity at any eccentricity tested. Predicted pressures,
  loads and void fractions from that configuration describe the starved branch or
  the `theta_min` clamp, not a bearing.
- **The reformation angle**, in any configuration. Ungrooved it does not exist;
  grooved it is an input. Elrod-Adams sets the switch function to zero in
  cavitated cells, so pressure-driven re-flooding is absent by construction --
  yet Etsion & Ludwig (1982) photographed strong reverse-flow re-flooding in a
  real submerged bearing.
- **Cavity pressure structure.** The cavity is a single uniform pressure by
  construction, against ~50 kPa of measured variation.
- **Any cavitated result on an unaudited mesh.** The resolution table above is
  the reason: an under-resolved free boundary passes every convergence check the
  solver reports.

## Verdict 2 -- is the energy equation validated?

**Partially: the film temperature magnitude is right, its distribution is not,
and the underlying pressure field carries too much load.** Target B is the only
target in this set that tests temperature at all, so everything below rests on a
single experiment.

Trustworthy:

- **Peak film temperature, under near-adiabatic walls.** 47.48 C against a
  measured 49.1 C, a 1.62 K error on a 9.1 K measured rise, against a 5 K
  criterion.
- **The size and direction of the thermal effect on pressure.** At a fixed
  eccentricity the energy equation knocks 8.3% off the peak pressure (1.2751 ->
  1.1777 MPa) through a viscosity drop to 0.771 of the inlet value. This is the
  quantity that justifies enabling the energy equation at all, and it behaves
  correctly.
- **Circumferential pressure profile shape.** Peak within 4.2 deg of the measured
  position and rupture inside the measured bracket, on the same solve.

Not trustworthy:

- **The location of the temperature peak.** Predicted at X/pi = +0.95, measured
  at +0.55 -- roughly 71 deg out. The predicted profile climbs monotonically to
  the groove; the measurement peaks and then falls. Alakhramsing et al. report
  exactly the same failure for their own model and attribute it to how heat
  transport in the cavitated region is treated, so it is a limitation of 2-D THD
  film models rather than something peculiar to pancake -- but it fails the
  stated criterion that the peak be in the right place.
- **The wall heat-transfer coefficient.** Not a prediction: it is an input the
  answer is very sensitive to. Going from adiabatic to h = 5000 W/(m2.K) moves
  the film rise from 7.48 K to 3.00 K and the temperature error from -1.6 K to
  -6.1 K, while moving the peak pressure by only 1.2%. Temperature results
  should be quoted with the assumed wall condition attached.
- **The isothermal bracket does not enclose the measurement.** The criterion was
  that isothermal must overpredict; it gives 1.2751 MPa against a measured
  1.29 +/- 0.03 MPa, falling 1.1% short. Within digitising uncertainty, but the
  strict test is not met, which says the pressure field is slightly
  under-predicted even with no thermal softening at all.

### The load-capacity error, and why it is probably a scale error

Load is over-predicted at every measured point, and the error grows with
eccentricity:

| eps (measured) | W measured [N] | W simulated [N] | error |
|---|---|---|---|
| 0.47 | 2000 | 3048  | +52.4% |
| 0.70 | 5000 | 8402  | +68.0% |
| 0.78 | 8000 | 12996 | +62.5% |

Equivalently, all three thermal treatments need eps = 0.52-0.53 to carry the
4 kN of Figure 9, where Ferron measured 0.67. It is not a thermal effect --
adiabatic, metal-bush and isothermal agree on the eccentricity to within 0.02.

It is also **not a groove-placement artefact**. Sweeping the groove angle at
fixed eccentricity (0.5345, 2000 rpm, THD adiabatic) over a 90 deg range moves
the load by under 5%:

| groove [deg] | 270 | 290 | 306 | 320 | 340 | 0 |
|---|---|---|---|---|---|---|
| W [N] | 3966 | 4016 | 3980 | 3995 | 4067 | 4155 |
| p_max [MPa] | 1.1534 | 1.1711 | 1.1772 | 1.1814 | 1.1863 | 1.1902 |

306 deg is the self-consistent placement (groove opposite the load line). The
total spread, 3966 to 4155 N, is 4.8% -- an order of magnitude smaller than the
discrepancy being explained.

The two apparently contradictory results reconcile under a single scale factor:

- at Ferron's measured eccentricity (0.67) the model gives W = 7238 N against
  4000 N and p_max = 2.40 MPa against 1.29 -- both about 80% high;
- at Ferron's measured *load* (4 kN) it gives p_max = 1.16 MPa against 1.29 --
  10% low.

Scaling the whole pressure field by `k` and matching the load at eps = 0.67 gives
`k` = 4000/7238 = 0.553, which then predicts p_max = 2.400 x 0.553 = **1.33 MPa
against a measured 1.29** -- a 2.8% match. One scalar fixes both.

Since `p ~ mu omega R^2 / c^2`, that scale factor corresponds to a running
clearance about 35% above the nominal 145 um, or an effective viscosity of
0.0153 rather than 0.0277 Pa.s. Alakhramsing et al. name thermal growth of the
clearance as an unmodelled source of discrepancy in exactly this comparison, and
Ferron's bush is bronze on a steel shaft.

**This is a hypothesis consistent with the data, arrived at after seeing it, not
an independent confirmation, and it is not a licence to tune the clearance until
the benchmark passes.** What it does establish is that the *shape* of the
predicted field -- peak location, rupture location, profile form and thermal
response -- is right, and that the disagreement is concentrated in one overall
scale. Settling it needs a number from the primary Ferron paper (hot running
clearance, or the measured eccentricity datum), which is paywalled.

## Limitations that bound what any of this can show

- **No axial re-flooding of a cavitated film.** The Elrod switch function is zero
  in cavitated cells, so liquid cannot re-enter by pressure flow. This is the
  direct cause of the ungrooved collapse.
- **Uniform cavity pressure.** The cavity sits at a single `p_cav` (or `p_sat`).
  Etsion & Ludwig (1982) measured ~50 kPa of variation inside a real cavitation
  zone.
- **Zero gas release under `LOCAL_PSAT`.** The threshold clamps `p >= p_sat`, so
  the flash releases nothing and `alpha_gas` is identically zero. Not exercised
  here (Target A uses `SCALAR_PCAV`), but it applies to the optical-rig sweep.
- **`VOID_COUPLED` does not converge** and was not enabled.
- **The Andrade viscosity law is a local tangent.** Ferron's oil is specified by
  an exponential-in-temperature law; the refit to pancake's Andrade form is
  accurate to 0.56% in viscosity over 40-55 C and should not be extrapolated.
- **Exit code is not convergence.** Every case in this section records the Elrod
  `mc`/`ns` and the coupling `R_p`, and is accepted only on those.
