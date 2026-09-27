# Six steady FW68D comparisons: A–F

This campaign replaces the earlier 24-point transient plan with the six cases
selected by the user. It uses a steady pressure–energy–CHT solve and an outer
load-balance search. The package now includes the groove, housing and still-oil
cooling. Pressure datum and material/shaft thermal inputs remain unresolved;
no physical cases have been run.

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
  maximum 200 sweeps, adaptive relaxation. An iteration cap or restored best
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

Use TC0 for confirmed feed/reservoir temperatures, surrounding fluid and the
numerical initial temperature guess. Compute the film and wall temperatures.
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
Stainless grade/conductivity remains unspecified; the generic solver default
50 W/(m K) is not assigned. As a reference only, the
[Outokumpu Core datasheet](https://www.outokumpu.com/en/products/product-ranges/-/media/files/products/core/outokumpu-core-range-datasheet.pdf)
lists 15 W/(m K) at 20 degrees C for 304. This does not identify the rig material.

The requested cavitation pressure is 0 Pa; its gauge/absolute datum is pending.
The current validator requires positive absolute cavity pressure, though exact
zero is mathematically compatible with this pure-oil FBNS pressure-offset
formulation. Zero absolute would need a narrow validation change and regression
check. No positive floor or atmospheric conversion has been substituted. The
common groove/end reservoir pressure remains a separate unresolved input.

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

Preliminary estimates below use **unconfirmed timing-fixture assumptions**
cp=2000 J/(kg K) and k_o=0.13 W/(m K). They are not production boundary inputs
or simulated temperatures. Wall excesses are illustrative.

| Case | h at 5 K excess | h at 10 K excess | h at 20 K excess |
|---|---:|---:|---:|
| A | 42.06 | 53.66 | 71.34 |
| B | 63.02 | 80.10 | 104.99 |
| C | 41.43 | 52.83 | 70.22 |
| D | 41.22 | 52.57 | 69.87 |
| E | 56.19 | 71.32 | 93.51 |
| F | 41.26 | 52.62 | 69.94 |

All h values are W/(m2 K). Full-precision estimates and Ra/Pr are recorded in
`steady_setup.json`. Recalculate when oil thermal properties are established.

The solver currently takes fixed h per solve. The planned thermal iteration
updates mean h from the converged mean outer-surface temperature inside each
force evaluation; this coupling is not implemented. A pressure/temperature
solve with frozen trial h does not verify the natural-convection closure.

L/D=0.47 is short, so this long-cylinder lateral-surface estimate is approximate.
The two annular end faces have about 79% of the outside cylindrical area, but
axial solid conduction/end-face cooling is absent from the r-theta CHT model.
Confinement and mounting contacts are also unresolved. Report these limits
and test cooling sensitivity; this is not a full finite-housing heat-loss model.

## Inputs still required

- Common groove/end absolute pressure and whether cavity zero is gauge or
  absolute. Logged pressure remains an observation until assigned.
- Stainless grade or conductivity, shaft cooling condition, oil heat capacity
  and conductivity. These are not determined by still-oil cooling.
- TC2 and pressure sensor positions for comparisons, not for solving the PDEs.

Full-length groove and horizontal housing are explicit model interpretations.
Unresolved numeric production inputs remain null in the preparation record.

## Execution and acceptance

Once the physical boundaries are established, make cross-film wall shear
consistent with the existing F0/F1 flow integrals before quantitative friction
comparison. The current force/torque routine still uses mean-film viscosity.
Verify the variable-viscosity wall stress and its constant-viscosity limit with
the existing force tests. This source correction has not been made in this
input-preparation task.

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
