# Six steady FW68D comparisons: A–F

This campaign replaces the earlier 24-point transient plan with the six cases
selected by the user. It uses a steady pressure–energy–CHT solve and an outer
load-balance search. The package now includes the groove, housing and still-oil
cooling. The user supplied 1 atm reservoirs and exposed shaft ends and
authorized ordinary/similar material properties. Local and remote focused checks pass; production equilibria have not yet
been accepted.

## Accepted cases

Values displayed here are rounded. `experiment.json` retains full precision,
raw ranges, timestamps, pressure observations and comparison targets from the
supplied JSON. `provenance.json` records the original JSON/CSV paths and hashes.

| Case | Experimental run | RPM | Inlet/ambient reference TC0 (°C) | Measured load (N) | Experimental S | Comparison |
|---|---|---:|---:|---:|---:|---|
| A | 03 | 200 | 29.67355 | 5276.271 | 0.03284994 | Low S |
| B | 02 | 1500 | 61.48691 | 3237.814 | 0.09531855 | Pair B/C |
| C | 04 | 200 | 28.64013 | 1973.003 | 0.09476705 | Pair B/C |
| D | 04 | 200 | 28.30434 | 382.524 | 0.50098435 | Flatter region |
| E | 01 | 1500 | 51.56526 | 461.252 | 0.99599165 | Pair E/F |
| F | 04 | 200 | 28.37163 | 190.634 | 1.00031890 | Pair E/F |

For E, add inlet/ambient-temperature sensitivity cases at 51.23781185626001 °C
and 51.92767223620965 °C, holding its measured mean load and nominal speed
fixed. Recompute properties and re-close load balance in each variant. Their
reference S values are 1.00915308 and 0.98162479. These two additional cases
bring the initial matrix to eight equilibria; each equilibrium may require
multiple CFD evaluations. No physical timestep refinement applies to steady mode.

## Solver and load control

- Geometry: R=0.0215 m, L=0.040 m, radial clearance c=42.5e-6 m; bore 43.085 mm.
- Film mesh: 720×192; bush radial cells: 16; cross-film order 4, quadrature 6.
- `solution_mode=STEADY_STATE`, `motion_model=STATIC`, `omega_ramp_time=0`.
- Pure-FW68D `SINGLE_PHASE`, FBNS cavitation, STRIATED cavity shear, film energy,
  polynomial cross-film temperature and `CONDUCTION_2D` bush conduction.
- Retain dt=20 us as configuration metadata. Steady Reynolds/energy storage
  terms are disabled and bearing motion is not integrated.
- Inner pressure/temperature coupling: both normalized changes below 1e-6,
  maximum 2000 sweeps, fixed relaxation 1. An iteration cap or restored best
  iterate without convergence is a failed trial.
- Outer equilibrium: vary eccentricity and attitude until the total fluid force
  on the bearing balances the specified applied load vector to 0.2%.

Eccentricity and attitude are solved outputs, so measured values are not needed
for this load-controlled calculation. Choose applied load along -y (270 degrees)
and the opposite groove at +y (90 degrees). The fluid-on-bearing target is
(0,+W). `load_angle_deg` currently does not rotate the solver force target;
the outer search must use explicit force-vector components in the chosen frame.
The static geometry in each inner solve represents an equilibrium position of
the loaded bearing, with the journal centre fixed and the journal rotating.

The reusable pattern is the damped two-variable Newton search in
[the corrected oil scan](../02c_oil_corrected_geometry/run_scan.py), with both
pressure and temperature convergence checks from
[the Ferron CHT workflow](../01e_ferron_terra_repro/force_balance.py).
Neither script may be run unchanged: both contain case-specific geometry,
feed and thermal assumptions. The new setup is in `steady_setup.json`, which
is a preparation record, not a runnable solver configuration.

## Temperature and property assignment

Use TC0 for feed/reservoir temperatures and surrounding fluid. The numerical
initial film guess is TC0 + 1 K. Compute the film and wall temperatures.
Do not impose TC0 on either wall without a separately specified wall boundary.
Logged pressure is retained as an observation; no boundary uses it until its
sensor location is known. Compare TC2 only at the corresponding sensor position.

`fw68d_experiment.csv` matches the supplied experimental viscosity relation:
log-linear kinematic viscosity between its rounded 10 °C temperature nodes,
with rho=954.9−0.7284 T_C and mu=rho nu. The checked-in general table uses
slightly shifted digitized temperature nodes, giving reference viscosities up
to 1.36% different for these cases. The case table replaces only the 15 pure-oil
rows; its other 110 rows remain solely to satisfy the existing table loader.
Zero refrigerant concentration is used for this campaign. No refrigerant
solubility or release calculation is enabled.

`viscosity_checks.json` records all six independent viscosity, S and friction
recalculations, CSV/JSON agreement and E sensitivity values. No extrapolation
of the kinematic table is allowed. The supplied linear density law is an
adopted PVE approximation; cases B/E exceed its stated 15–50 °C fit range.
The present energy model uses constant volumetric heat capacity `rho_cp`;
its value must be specified consistently with the selected oil heat capacity.

## Specified groove and housing

The user specifies a 5 mm groove opposite the applied load, at the same
reservoir conditions as both bearing ends. Treat the width as a circumferential
arc at the 43.085 mm bore and use the native full-length axial groove:
width = (0.005/0.0215425) * (180/pi) = 13.2983125 degrees; center = 90 degrees.
At 720 theta cells the native cell-center mask selects 26 columns (13 degrees,
4.88784 mm bore arc). Retain the nominal 5 mm input and report this discretization.
Groove pressure equals the common end pressure; its supplied temperature is
case TC0. The native fed-groove PINNED closure fixes reservoir cells to TC0.

The stainless housing has OD 85 mm, outer radius 42.5 mm and length 40 mm.
Its physical radial wall is 20.9575 mm thick. The r-theta solid solver starts
at journal R=21.5 mm, a thin-film approximation to bore radius 21.5425 mm.
Use 304 stainless as the authorized ordinary-material analogue, with constant
conductivity 15 W/(m K), from the
[Outokumpu Core datasheet](https://www.outokumpu.com/en/products/product-ranges/-/media/files/products/core/outokumpu-core-range-datasheet.pdf).

Groove and end reservoirs are 101325 Pa absolute. Retain the requested cavity
plateau as 0 Pa absolute, distinct from the reservoirs. The scoped validator
change permits zero only for this gas-free SINGLE_PHASE/FBNS/SCALAR_PCAV path;
negative pressures and zero pressure in other gas/coupled paths remain invalid.

## Adopted oil properties and exposed shaft ends

Use constant cp=2010 J/(kg K) and k=0.13 W/(m K), the published 40 C thermal
properties of BITZER BVC32 PVE, as the authorized analogous oil values
([manufacturer data](https://www.bitzer.de/shared_media/html/est-500/en-GB/271261579271263371.html)).
Keep the supplied FW68D viscosity and density laws. Film energy uses constant
rho_cp=rho(TC0)*cp within each case. These cp/k values are not FW68D measurements.

Both exposed shaft ends have radius 21.5 mm and total area 0.0029044024 m2.
Use shaft_ambient_area_ratio=R/L=0.5375 and ambient temperature TC0. The existing
CHT model solves one global shaft temperature, not axial solid conduction.
For the rotating disks use Nu_R=0.62 Re_R^(1/2) Pr^(1/3), Re_R=|omega| R^2/nu,
and h_end=Nu_R*k/R at mean shaft/oil temperature. This is a high-Pr, laminar,
unconfined-disk estimate (Liu and Stewart,
[*Int. J. Heat Mass Transfer*](https://doi.org/10.1016/0017-9310(72)90181-0), 1972).
Record Re/Pr and distinguish this rotating-end coefficient from housing cooling.

## Natural convection to surrounding oil

The user confirmed still surrounding FW68D at TC0 and requested estimated h.
Use convection with a solved housing-surface temperature. For an assumed
horizontal cylinder, estimate mean lateral-surface h with the correlation
(Churchill and Chu, [*Int. J. Heat Mass Transfer*](https://doi.org/10.1016/0017-9310(75)90222-7), 1975):

$$Nu_D=\left[0.60+\frac{0.387 Ra_D^{1/6}}
{[1+(0.559/Pr)^{9/16}]^{8/27}}\right]^2,\qquad h=Nu_D k_o/D.$$

Use D=0.085 m and properties at (mean outer-surface temperature+TC0)/2:
alpha=k_o/(rho cp), Pr=nu/alpha, Ra_D=g beta |T_surface-TC0| D^3/(nu alpha),
and beta=0.7284/rho from the supplied density approximation. Use the supplied
FW68D viscosity table. The boundary is q''=h(T_surface-TC0).

The driver iterates the housing and shaft coefficients within every force
evaluation, using solved mean surface temperatures and relaxation. The mean
housing outer temperature is reconstructed from the bore mean and exact mean
radial/Robin resistance. Accept coefficient changes below 1e-3 relative and
wall-temperature changes below 0.01 K, alongside both inner coupling residuals
below 1e-6. A fixed-coefficient timing trial is explicitly labelled unclosed.

L/D=0.47 is short, so this long-cylinder lateral-surface estimate is approximate.
The two annular end faces have about 79% of the outside cylindrical area, but
axial solid conduction/end-face cooling is absent from the r-theta CHT model.
Confinement and mounting contacts are also unresolved. Report these limits
and test cooling sensitivity; this is not a full finite-housing heat-loss model.

## Remaining comparison information

TC2 and pressure sensor positions are needed to compare measured temperatures
and pressure at corresponding points. They do not block the numerical solves.
Full-length groove, horizontal housing, isothermal shaft and analogous
materials are explicit modeling assumptions. Revisit cooling sensitivity when
assessing experimental agreement; chosen coefficient/property variations are
sensitivity scenarios, not measured confidence bounds.

## Execution and acceptance

`run_steady.py` renders the exact cases and iterates eccentricity/attitude to
fluid-on-bearing target (0,+W). The F0/F1 wall-stress correction is implemented
and checked against variable-viscosity and uniform-viscosity limits. Every
solver attempt retains its config, log, hashes, diagnostic table and fields;
failed or capped inner/thermal/force solves remain failed.

```sh
python run_steady.py --self-check
python run_steady.py --solver /absolute/path/pancake --ranks 2 \
  --work /absolute/path/new-results --cases F,E,D,C,B,A,E_low,E_high
```

A frozen-position timing call is `--trial F .15 165 --thermal-iterations 1`.
It is not load-balanced or externally thermally converged. On Open MPI only,
`--rankfile` assigns one rank to each chosen physical core. The workstation
left ranks unbound with `--cpu-list`, so that option is not used. The user-requested resolution remains 720 by 192 throughout.

Follow the repository workflow: local build and short checks, isolated remote
source/configuration snapshot, fresh remote build, converged equilibrium runs,
and return of evidence. Preserve every failed evaluation and record source,
input and executable hashes, rank layout, residuals and return codes. Existing
coarse transient timings do not establish an optimum for these steady solves.

A completed comparison must report inner thermal/pressure convergence and outer
vector load closure separately; finite fields, positive gap and liquid filling;
minimum gap, eccentricity and attitude; film and wall temperature maps; torque
and f=M/(W R); cavity area/fronts and filling deficit. Keep mass/heat-balance
diagnostics with the results and identify any diagnostic the solver does not
provide. Numerical convergence does not itself establish experimental agreement.

The old failed transient runs remain in
[the MPI benchmark record](../90_diagnostics/fw68d_mpi_benchmark/RESULTS.md).
Their startup failure is not a result for this steady six-case campaign.

### Preflight correction and numerical settings (2026-09-27)

Two native diagnostic trials reached the 200-iteration cap (return code 4).
They also inherited a 313.15 K housing reservoir override from the old timing
fixture. Both are invalid for the requested TC0 housing boundary and remain
saved as failed attempts. The case driver explicitly sets the housing, shaft,
groove and axial-end reservoir references to each case TC0.

Use fixed coupling relaxation 1, best-iterate retention and a cap of 2000;
pressure and temperature tolerances remain 1e-6. Initial film temperature is
TC0 + 1 K, a numerical guess only. Reaching the cap or restoring a best iterate
is still a failure. The initial trials varied several numerical settings and
cannot identify which individual change improved their convergence.

The ENERGY_EQUATION/cross-film path leaves the global energy-residual CSV field
blank. The driver records that diagnostic as unavailable, never zero or
verified. Wall heat fluxes remain available for inspection.

See [the execution record](EXECUTION.md) for verified builds, short MPI timings,
failed attempts and current physical-run status.
