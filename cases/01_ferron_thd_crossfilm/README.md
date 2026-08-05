# Case 01-THD-crossfilm

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Ferron, Frene & Boncompain (1983) Fig. 4 point (2000 rpm, 4 kN) re-run with
`cross_film_thermal_model = POLYNOMIAL` under `cavitation_model = FBNS`, in three
thermal wall variants.

Why: case 01-THD (`../01_ferron_thd`) under-predicted the measured temperature by
5-7 K, and one of the two reasons it could not be scored properly was that with
`cross_film_thermal_model = NONE` the solver has no cross-film profile, so the
only temperature it can plot is the depth-averaged film value, while the
experiment is a bush temperature measured 0.5 mm below the bush inner surface.
The cross-film model produces a bush-side film surface temperature, which is the
closest model quantity to that thermocouple.

Nothing here is tuned to the experiment.

## Files

- `config_template.txt` - parametrised config (e, attitude, grid, the three wall
  HTCs, groove mixing model, Q_supply).
- `solver_io.py`, `plot_style.py` - copies of the case 01-THD harness, unchanged.
- `force_balance.py` - the two-pass Newton driver, three wall variants.
- `probe.py`, `probe_cfg/` - coarse-grid config validation, and the measurement
  that the large-HTC limit really is a fixed-wall boundary condition.
- `analyse_thd.py`, `thd_analysis.json` - every reported number.
- `plot_thd.py` - the two deliverable figures.
- `equilibrium.json`, `<variant>/fb_log.json` - converged points and every iterate.
- `walls_adiabatic/`, `bush40_Jadiabatic/`, `bush40_J40/` - one archived run per
  solver call (`fb_iter_NN/` with config, run.log and results).
- `progress.log` - running record.

## Wall variants, and what they map onto in the config

The conjugate cross-film closure has **no prescribed-wall-temperature option**.
The shaft is circumferentially isothermal at a temperature solved from the global
film-to-shaft heat balance; the bush conducts radially to an outer convection
boundary. A "fixed 40 C wall" is therefore expressed as the large-conductance
limit, in which the surface is driven onto `thermal_ambient_temperature`:

| variant | `shaft_ambient_heat_transfer` | `bush_outer_heat_transfer` | `bush_thermal_conductivity` |
|---|---|---|---|
| A `walls_adiabatic` | 0 | 0 | 50 (inert) |
| B `bush40_Jadiabatic` | 0 | 1e7 | 1e7 |
| C `bush40_J40` | 1e7 | 1e7 | 1e7 |

with `thermal_ambient_temperature = 313.15` and `bush_outer_radius = 0.10` in all
three. `shaft_ambient_heat_transfer = 0` is a *globally* adiabatic shaft: the
integral of the film-to-shaft flux is zero, so heat entering the shaft in the
loaded arc leaves it again upstream. That is a different boundary condition from
a per-cell adiabatic shaft, and it is the one the model implements.

`probe_cfg/probe.json` measures what the limits actually realise on a coarse grid:
every "fixed" surface lands on 40.000 +/- 0.001 C, so the mapping is exact to
1 mK. Under the cross-film model `journal_heat_transfer`,
`bearing_heat_transfer`, `journal_wall_temperature` and `bearing_wall_temperature`
are ignored by the solver; they are pinned to inert values in the template so
nothing in the file looks effective when it is not.

## Everything else is frozen from case 01-THD

Geometry (R = 0.05 m, c = 1.52e-4 m, L = 0.08 m), grid (720 x 192, G6),
omega = 209.43951 rad/s, `cavitation_model = FBNS`, p_cav = 1.0e5 Pa,
`fluid_property_model = SINGLE_PHASE`, `liquid_viscosity_model = POLYNOMIAL`
(k0 = 3.287, k1 = 3.064, k2 = 0.777, mu = 0.0277 Pa s at 313.15 K), oil
rho = 860, rho c_p = 1.72e6, k = 0.13, the 18 deg fed groove at theta = 90 deg
(0.70 bar gauge, 40 C), `groove_mixing_model = FLOW_WEIGHTED` with the same
self-consistent two-pass Q_supply closure, Dirichlet ambient axial ends with
`ZERO_GRADIENT` thermal, `coupling_tolerance = 1e-6`, `mpiexec -n 4`.

## Model caveats that travel with every number here

PHYSICS.md 13.1.2 states two approximations inside the order-4 cross-film
collocation, and they apply to every temperature in this case:

1. The midplane collocation is **Couette-only**: the pressure-gradient
   contribution to the shear profile is dropped, because including it makes the
   profile shape a function of the local pressure gradient, which jumps whenever
   a cell flips between full film and cavitated, and the segregated solve then
   limit-cycles.
2. The convective term is **saturated against the dissipative one**. A single
   midplane collocation cannot represent the wall layers that appear once the
   reduced Peclet number leaves order one, and on this rig it is 10-30. Where the
   saturation binds, the quartic degenerates towards a fixed flattened shape
   rather than resolving the wall layers.

So the cross-film profile is a reconstruction collocated with the depth-averaged
energy balance, not a solution of the full 3-D energy equation.
