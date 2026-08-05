# Case 01-THD - Ferron et al. (1983), FBNS + energy equation at 2000 rpm

Purpose: place the solver at the Ferron Fig. 4 operating point (2000 rpm, 4 kN)
with the energy equation and the Ferron polynomial viscosity law switched on, and
sweep the Fig. 5 load line. Case 01-prelim showed an isothermal film brackets the
measured pressure between two hand-picked viscosities; this case removes the hand
picking and lets mu(T) come out of the solve.

Nothing in this case is tuned to the experiment. Every deviation below is
reported as found.

Branch `codex/viscosity-polynomial`, HEAD `6dacec0`. Build preset
`windows-native-mingw` (Release), `ctest` 27/27 green before the batch.

## Parameter table

Common to every run. Base config `cases/00_grid_dt_independence/grid_G6_fbns/config.txt`
(G6 grid, 720 x 192); the parametrised copy is `config_template.txt`.

| item | value | provenance |
|---|---|---|
| geometry | R = 0.05 m, c = 1.52e-4 m, L = 0.08 m | user-verified (Case 00/01-prelim frozen set) |
| grid | n_theta = 720, n_z = 192 (G6) | grid-independent per Case 00a |
| speed | omega = 209.43951 rad/s (2000 rpm) | user-verified |
| mode | `STEADY_STATE` + `STATIC` + `ENERGY_EQUATION`, `cavitation_model = FBNS` | protocol |
| cavitation datum | p_cav = 1.0e5 Pa (absolute) | protocol |
| fluid model | `fluid_property_model = SINGLE_PHASE`, no dissolved gas | protocol |
| viscosity | `liquid_viscosity_model = POLYNOMIAL`, k0 = 3.287, k1 = 3.064, k2 = 0.777, mu = 0.0277 Pa s at T_ref | user-verified (Ferron eq. for their oil) |
| reference T | `property_reference_temperature` = `temperature_initial` = `temperature_reference` = 313.15 K | user-verified |
| oil thermal | rho = 860 kg/m^3, rho c_p = 1.72e6 J/(m^3 K), k = 0.13 W/(m K) | user-verified |
| feed | `inlet_groove = 90 18 1.713e5 313.15 ...` - 18 deg axial groove at theta = 90 deg, 0.70 bar gauge, fed at 40 C | user-verified; thermally fed groove is mandatory (a transparent groove runs away thermally) |
| groove mixing | `groove_mixing_model = FLOW_WEIGHTED` | protocol; see the Q_supply note below |
| bush wall | `bearing_heat_transfer` = 7200 W/(m^2 K), `bearing_wall_temperature` = 313.15 K | Tier 1, user-signed: h_b = K_b / (R1 ln(R2/R1)) with K_b = 250 W/(m K), R1 = 0.05 m, R2 = 0.10 m -> 7213, rounded to 7200 |
| journal wall | TWO variants bracketing the shaft: **J-adiabatic** `journal_heat_transfer = 0`; **J-iso40** `journal_heat_transfer = 10000`, `journal_wall_temperature = 313.15 K` | Tier 1, user-signed |
| cross-film | `cross_film_thermal_model = NONE` (default) | protocol did not enable it; consequences stated under the T figure |
| axial ends | Dirichlet 101325 Pa, thermal `ZERO_GRADIENT` | protocol |
| Picard coupling | `coupling_max_iters = 200`, `coupling_tolerance = 1e-6`, `coupling_relaxation = 1.0` | numerics only |
| ranks | `mpiexec -n 4` | local machine |

### Deviation 1 - the FLOW_WEIGHTED groove needs a supply flow rate that was not signed off

`groove_mixing_model = FLOW_WEIGHTED` sets the groove temperature from
lambda = Q_supply / (Q_supply + Q_up), and the solver **requires** a positive
`supply_flow_rate` token on the inlet (`src/config.hpp`, the FLOW_WEIGHTED
validation branch). The signed-off parameter set does not give one and the
digitised Ferron material in `cases/ref_data/` does not contain a measured feed
flow. Rather than invent a number, every operating point is solved twice:

1. **pass `pinned`** - `groove_mixing_model = PINNED` (lambda = 1, the historical
   hard pin, needs no Q_supply). Converged by the Newton search.
2. **pass `flow`** - `FLOW_WEIGHTED` with Q_supply taken from pass 1's own steady
   mass balance, Q_supply = `inlet_mass_source` / rho. In steady state the net
   liquid injected by the groove penalty cells *is* the fresh-oil charge
   (everything else crossing the groove is recirculated film), so this is the
   solver's self-consistent value, not a fitted one. The Newton search is then
   re-converged from the pass-1 point.

Q_supply is held fixed inside pass 2 rather than re-derived each iterate; the
residual inconsistency is small (the pass-2 mass balance returns 21.23 L/min
against the 20.87 L/min it was given at 4 kN, i.e. 1.7 %). All reported results
are the `flow` pass. The `pinned` pass is kept as the hot-oil carry-over
sensitivity and is tabulated alongside.

Derived Q_supply per point: 14.0 (2 kN), 17.7 (3 kN), 20.9 (4 kN), 26.1 (6 kN),
30.1 (8 kN), 33.3 (10 kN) L/min. These are large but they are dominated by
pressure-driven side leakage off the hydrodynamic pressure hill, not by the
0.70 bar groove overpressure, and they scale with load as that reading predicts.

## Force-balance search

Damped Newton on (eps, `attitude_angle_deg`) with a forward-difference Jacobian
(`force_balance.py`), identical in method and sign convention to Case 01-prelim:
the pressure-force target is Fp = (0, -W) N, i.e. magnitude W at -90 deg.
Tolerance |W| to 0.1 % and force direction to 0.05 deg. Every solver call,
including the Jacobian probes, is archived as `<group>/fb_iter_NN/` with its
config, log and `results/`; `<group>/fb_log.json` maps each call to its metrics.

## Part A - Fig. 4 point (2000 rpm, 4 kN)

### Equilibrium

| quantity | J-adiabatic | J-iso40 |
|---|---|---|
| eccentricity ratio eps = e/c | **0.573803** | **0.567497** |
| `attitude_angle_deg` (min-film line, solver input) | **-36.3237 deg** | **-36.1187 deg** |
| classical attitude angle phi = psi - angle(Fp) | 53.68 deg | 53.88 deg |
| eccentricity e | 8.7218e-5 m | 8.6259e-5 m |
| Fp | (0.37, -3999.70) N | (0.06, -3999.88) N |
| **W closure** | 3999.70 N (-0.0075 % of 4000) | 3999.88 N (-0.0030 %) |
| force direction | -90.0053 deg (0.0053 deg off) | -90.0009 deg (0.0009 deg off) |
| h_min | 64.78 um | 65.74 um |
| p_max (gauge, whole field) | 1.2703 MPa | 1.2653 MPa |
| T_max (whole field) | 43.92 C at 279.8 deg | 42.20 C at 256.8 deg |
| T_min (groove) | 40.26 C | 40.30 C |
| mu range | 0.0238 - 0.0274 Pa s | 0.0255 - 0.0276 Pa s |
| cavitated area fraction | 0.258 | 0.259 |
| friction torque | 2.808 N m | 2.900 N m |
| derived Q_supply | 3.4788e-4 m^3/s (20.87 L/min) | 3.4581e-4 m^3/s (20.75 L/min) |
| Newton steps (pass `flow`) | 1 | 1 |
| Picard iterations at the solution | 8 (R_p = 0, R_T = 8.2e-8) | 6 |

Groove-mixing sensitivity (`pinned` vs `flow`): eps moves 0.5709 -> 0.5738
(+0.5 %) for J-adiabatic and 0.5670 -> 0.5675 (+0.1 %) for J-iso40. Hot-oil
carry-over is a small effect here because lambda comes out near 0.8.

### Figure 1 - mid-plane pressure

`ferron_thd_pressure.png` (300 dpi PNG, lab style, `plot_thd.py`).

**Caption / mapping assumption.** Mid-span film pressure at the force-balanced
equilibrium for both journal-wall variants, against the digitised Ferron et al.
(1983) Fig. 4 experiment at 2000 rpm and 4 kN. The experimental abscissa is
degrees from the supply-groove CENTRE in the direction of rotation; the solver
groove centre is at theta = 90 deg and rotation is +theta, so the experiment is
placed at solver_theta = (90 + x) mod 360 and both curves are drawn against
(theta - 90) mod 360. This is the stated mapping applied as-is: no shift was
fitted. Pressure is gauge (solver absolute minus 101325 Pa). The theta origin of
the digitisation is flagged as unconfirmed against the PDF in the CSV header.

| | peak gauge pressure | peak location (deg from groove) |
|---|---|---|
| experiment | 1.2602 MPa | 201.0 |
| J-adiabatic | 1.2702 MPa (+0.8 %) | 199.75 (-1.25) |
| J-iso40 | 1.2653 MPa (+0.4 %) | 199.75 (-1.25) |

Reading: the pressure field is reproduced essentially exactly. Peak magnitude is
inside 1 % and peak location inside 1.3 deg for both wall variants, and the
rising flank, the crest and the rupture shoulder all lie on the data. The two
variants are visually indistinguishable, so the pressure solution carries no
information about the journal-side wall condition at this operating point.

This is the headline improvement over Case 01-prelim, where the isothermal film
needed a hand-chosen viscosity to bracket the peak and still ruptured too early
on the falling flank. With mu(T) solved, the rupture shoulder now sits on the
data through 240 deg. Beyond 250 deg the solver holds p_cav (-1.3 kPa gauge)
while the experiment reads a few tens of kPa gauge; that offset is the chosen
p_cav plus the real sub-ambient/gaseous behaviour the FBNS datum does not model.

### Figure 2 - mid-plane temperature

`ferron_thd_temperature.png` (300 dpi PNG, lab style, `plot_thd.py`).

**Caption.** Mid-span film temperature at the same equilibria, against the
digitised Ferron Fig. 4 bush temperature. Two quantities are being compared and
they are not the same thing:

- The experiment is a **bush temperature measured 0.5 mm below the bush inner
  surface**, not a film temperature.
- With `cross_film_thermal_model = NONE` the solver has **no cross-film profile at all**:
  the plotted curve is the depth-averaged film temperature. There is no
  "bush-side of the film" value to plot in this configuration. Producing one
  would require `cross_film_thermal_model = POLYNOMIAL`, which is a model change outside
  the signed-off protocol and brings its own unverified conjugate-wall
  parameters, so it was not enabled.

The expected sign of the systematic offset is therefore **model >= experiment**:
a depth-averaged film should sit above a sub-surface bush reading by a few K.

**The observed sign is the opposite, and that is the main negative result of this
case.**

| | peak T | peak location (deg from groove) | rise above the 40.5 C inlet |
|---|---|---|---|
| experiment | 49.04 C | 258.2 | 8.6 K |
| J-adiabatic (mid-plane) | 43.78 C | 281.8 | 3.3 K |
| J-iso40 (mid-plane) | 42.14 C | 260.7 | 1.6 K |

- The model under-predicts the peak by 5.3 K (J-adiabatic) to 6.9 K (J-iso40),
  and under-predicts the *rise* by a factor of 2.6 to 5.4. The journal-wall
  bracket does not contain the experiment: both bounds are too cold, so the
  bracket fails as a bracket.
- Peak LOCATION is not a discriminator here. This digitisation puts the measured
  maximum at 258 deg; the independent digitisation noted in the CSV header
  (Alakhramsing 2015 Fig. 9c, same condition) puts it at 279 deg. J-iso40 lands
  on the first (260.7), J-adiabatic on the second (281.8). A ~20 deg
  reproduction spread in the reference means neither variant is falsified or
  confirmed by location.
- The inlet level is right: the solved groove temperature is 40.26-40.30 C
  against a measured 40.5 C, so the fed-groove boundary condition and the
  flow-weighted carry-over are behaving.

### Deviation 2 - where the missing heat goes (diagnostic, not a deliverable)

The film energy balance is closed almost entirely by the Tier-1 bush term. At
h_b = 7200 W/(m^2 K) the bush is a very strong sink, so the film is pushed back
towards the 313.15 K wall almost as fast as shear generates heat. Adding the
J-iso40 journal term nearly doubles the total wall coefficient and roughly halves
the rise, which is exactly the behaviour of a wall-dominated balance.

`bush_htc_sensitivity.py` re-runs the converged J-adiabatic geometry (eps and psi
frozen, load therefore no longer exactly 4 kN) at lower bush conductivities.
Nothing is re-fitted; this is a diagnostic only and none of it enters the figures.

| K_b [W/(m K)] | h_b [W/(m^2 K)] | T_max [C] | rise [K] | W [N] | p_max gauge [MPa] |
|---|---|---|---|---|---|
| 250 (signed off) | 7213 | 43.92 | 3.66 | 4000 | 1.270 |
| 120 (leaded bronze) | 3462 | 45.48 | 5.13 | 3906 | 1.240 |
| 50 (steel) | 1443 | 47.13 | 6.71 | 3821 | 1.212 |

The temperature deficit is owned by the bush wall closure, not by the energy
equation, the viscosity law or the cavitation model. **Flag for the user:
K_b = 250 W/(m K) is far above steel (~45) or leaded bronze (~60); it is closer
to aluminium or copper. If the intended bush material is bronze or steel, the
signed-off h_b is 2-5x too large and the temperature comparison should be
re-scored after that value is re-checked.** This is raised as a parameter query,
not applied - the deliverable figures use the signed-off 7200.

[COORDINATOR NOTE, added on filing: K_b = 250 W/(m K) is the value printed in
Ferron's Table 1 (verified against the PDF text), so it stays. The suspect
Tier-1 input is instead the 40 C outer-bush datum: the paper fixes T = T_b at
R2 where T_b is the MEASURED bush surface temperature (nomenclature: "Tb = bush
surface temperature degC"), a per-test value not printed in Table 1 - and the
rig runs ~45 C. The Tier-1 datum assumption, not K_b, is the leading suspect
for the 5-7 K deficit. See the verdict discussion in chat.]

## Part B - Fig. 5 load sweep (2000 rpm, J-adiabatic)

Force balance at each load, `flow` pass. The 4 kN point is the Part A
J-adiabatic result, reused rather than re-run.

| W [N] | eps | psi [deg] | phi [deg] | h_min [um] | p_max gauge [MPa] | T_max [C] | cav. frac. | torque [N m] | Q_supply [L/min] | Picard | Newton |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2000 | 0.4061 | -25.258 | 64.75 | 90.28 | 0.571 | 42.93 | 0.218 | 2.498 | 14.03 | 7 | 1 |
| 3000 | 0.5051 | -31.675 | 58.33 | 75.22 | 0.907 | 43.42 | 0.241 | 2.658 | 17.66 | 7 | 1 |
| 4000 | 0.5738 | -36.324 | 53.68 | 64.78 | 1.270 | 43.92 | 0.258 | 2.808 | 20.87 | 8 | 1 |
| 6000 | 0.6636 | -42.665 | 47.34 | 51.13 | 2.064 | 44.92 | 0.281 | 3.075 | 26.10 | 9 | 1 |
| 8000 | 0.7207 | -46.883 | 43.12 | 42.46 | 2.929 | 45.88 | 0.296 | 3.305 | 30.12 | 10 | 1 |
| 10000 | 0.7607 | -49.962 | 40.04 | 36.38 | 3.855 | 46.81 | 0.306 | 3.508 | 33.31 | 10 | 1 |

Every point closed |W| to better than 0.1 % and the force direction to better
than 0.05 deg. phi decreases monotonically with load as a plain bearing must.

### Figure 3 - eccentricity vs load

`ferron_thd_eccentricity.png` (300 dpi PNG, lab style, `plot_thd.py`).

**Caption / rescale.** Equilibrium eccentricity ratio at 2000 rpm against the
digitised Ferron Fig. 5 experiment. The digitised trace is the C = 145 um curve;
this case runs c = 152 um, so the experimental points are drawn **rescaled by
145/152 = 0.9539**, which converts a fixed measured displacement to the present
clearance basis. This is the rescale the protocol specifies; no other adjustment
is applied.

Solver vs experiment interpolated onto the solver loads:

| W [N] | eps (solver) | eps (exp, C = 145 um) | exp x 145/152 | delta | exp - 0.10 (C_eff ~ 171 um) | delta |
|---|---|---|---|---|---|---|
| 2000 | 0.4061 | 0.4901 | 0.4676 | -0.0615 | 0.3901 | +0.0159 |
| 3000 | 0.5051 | 0.5961 | 0.5687 | -0.0636 | 0.4961 | +0.0090 |
| 4000 | 0.5738 | 0.6628 | 0.6322 | -0.0584 | 0.5628 | +0.0111 |
| 6000 | 0.6636 | 0.7494 | 0.7149 | -0.0513 | 0.6494 | +0.0142 |
| 8000 | 0.7207 | 0.8125 | 0.7751 | -0.0544 | 0.7125 | +0.0082 |
| 10000 | 0.7607 | 0.8491 | 0.8100 | -0.0494 | 0.7491 | +0.0115 |

Reading:

- The SHAPE is right. The solver curve is parallel to the experiment over the
  whole 2-10 kN range: the residual is -0.055 +- 0.005, i.e. a constant offset,
  not a growing or load-dependent error. A physics error in load capacity would
  not be flat in load.
- A constant eps offset is the signature of a clearance-basis disagreement, which
  the reference data itself flags. Ferron's eps is a proximity-probe displacement
  divided by an ASSUMED clearance, and the same figure carries a second
  C = 152 um legend entry. Applying instead the later re-reduction noted in the
  CSV header (Lubricants 11(7) 308, 2023: reduce the published eps by 0.10,
  C_eff ~ 171 um) leaves +0.011 +- 0.003, i.e. agreement to about one digitising
  scatter width across every load.
- Both readings cannot be right, and this case cannot resolve which clearance the
  experiment actually ran at. What it does establish is that the solver's eps(W)
  curve is consistent with the Ferron data to within the published clearance
  ambiguity, and that the ambiguity (0.10 in eps) is an order of magnitude larger
  than anything the solver contributes.
- The measured points are also softer than the solver at low load and the offset
  narrows slightly with load (-0.062 at 2 kN to -0.049 at 10 kN). That residual
  tilt is consistent with the under-predicted temperature: too little viscosity
  drop means too much load capacity, most visible where the film is thickest.

## Numerical health

- **Force balance.** Every group converged. The `pinned` pass took 2 Newton
  steps from the cold guess; the `flow` pass took 1 from the `pinned` solution.
  No Jacobian probe failed and no run hit the eps clamp.
- **Picard (outer pressure/temperature) coupling.** 6-16 iterations, rising with
  load as the viscosity feedback strengthens. Every reported solution satisfied
  `coupling_tolerance = 1e-6`; the worst final residual across the batch is
  R_T = 8.2e-8, and R_p is exactly 0 at most solutions (once the FBNS
  complementarity flags stop flipping, the pressure solve reproduces bitwise).
  `coupling_max_iters = 200` was never approached.
- **FBNS mass conservation.** `liquid_residual` is at round-off (|.| < 1e-10) in
  every archived final frame.
- **Cavitated fraction** grows smoothly with load, 0.218 to 0.306. No flag
  thrashing, no stalled rupture front.

## Runtime

89 archived solver calls at 720 x 192 on `mpiexec -n 4`, 12-22 s each,
1434 s of solver wall time in total (about 24 min). Plus 3 sensitivity runs
(47 s) and 2 setup probes.

| group | solver calls | wall [s] | max Picard |
|---|---|---|---|
| fig4_Jadiabatic | 11 | 164.4 | 8 |
| fig4_Jiso40 | 11 | 137.9 | 6 |
| fig5_W2000 | 14 | 189.4 | 7 |
| fig5_W3000 | 11 | 154.7 | 8 |
| fig5_W6000 | 14 | 224.9 | 10 |
| fig5_W8000 | 14 | 259.4 | 13 |
| fig5_W10000 | 14 | 303.3 | 16 |

## Deviations and caveats, collected

1. **Q_supply for FLOW_WEIGHTED was not in the signed-off set.** Derived
   self-consistently from each point's own steady mass balance (see Deviation 1).
   Not fitted, but not independently verified against the paper either.
2. **No bush-side film temperature exists in this configuration.**
   `cross_film_thermal_model = NONE`, so Figure 2 plots the depth-averaged film
   temperature. The comparison to a 0.5 mm sub-surface bush thermocouple is
   therefore doubly indirect.
3. **The temperature comparison fails in the wrong direction.** Model is 5-7 K
   COLDER than the experiment where physics says it should be a few K warmer.
   Traced to the Tier-1 bush closure; see Deviation 2 and the coordinator note.
4. **The journal-wall bracket does not bracket.** J-adiabatic and J-iso40 differ
   by only 1.7 K in peak temperature and are indistinguishable in pressure,
   because the bush term dominates the wall balance. The bracket as specified
   does not test what it was meant to test.
5. **Fig. 5 clearance basis is unresolved.** The prescribed 145/152 rescale
   leaves a flat -0.055 offset; the alternative published re-reduction leaves
   +0.011. Both are reported; neither is selected.
6. **Fig. 4 theta origin is unconfirmed against the PDF** (inherited caveat, see
   the CSV headers). If wrong, every experimental point in Figures 1 and 2 shifts
   by the attitude angle. The 1.25 deg pressure-peak agreement is strong
   circumstantial support for the origin used.
7. **Ferron's measured eps and attitude at this exact condition were not
   available** as a separate digitised quantity, so the Part A equilibrium is
   scored against the Fig. 5 curve interpolated to 4 kN rather than against a
   stated Fig. 4 value.

## Files

- `config_template.txt` - parametrised config (e, attitude, grid, journal HTC,
  groove mixing model, Q_supply).
- `solver_io.py` - run driver and VTK/pvts readers (01-prelim copy, 4 ranks,
  plus thermal metrics, Picard-health parsing and a generic midplane profile).
- `force_balance.py` - the two-pass Newton driver.
- `analyse_thd.py`, `thd_analysis.json` - every reported number.
- `plot_thd.py`, `plot_style.py` - the three deliverable figures.
- `equilibrium.json`, `<group>/fb_log.json` - converged points and every iterate.
- `fig4_Jadiabatic/`, `fig4_Jiso40/`, `fig5_W*/` - one archived run per solver
  call (`fb_iter_NN/` with config, run.log and results).
- `bush_htc_sensitivity.py`, `sensitivity_bush/`, `sensitivity_bush.json` -
  the diagnostic in Deviation 2. Not a deliverable.
- `probe_cfg/` (180 x 48 config validation), `probe_full/` (720 x 192 timing
  probe) - setup runs, kept for the record.
- `progress.log` - running record.
