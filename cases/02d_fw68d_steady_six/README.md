# Six steady FW68D comparisons: A–F

This campaign replaces the earlier 24-point transient plan with the six cases
selected by the user. It uses a steady pressure–energy–CHT solve and an outer
load-balance search. The input package is prepared; simulations await the feed,
load-direction and wall/solid thermal specifications listed below.

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
for this load-controlled calculation. The load direction relative to feed ports
is needed. `load_angle_deg` currently does not rotate the solver force target;
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

## Inputs still required

1. Load direction relative to the groove/ports and the chosen rotation frame.
2. Feed/groove geometry, supply pressure or flow, axial-end pressure and the
   logged pressure sensor location. The earlier 90°/10° groove and 1-bar feed
   are not newly confirmed by the six-case attachment.
3. Bush outer radius and conductivity, outer convection coefficient or measured
   wall-temperature boundary, and shaft cooling condition. For a convective
   boundary TC0 supplies the external fluid temperature; the coefficient is
   still required. The bore diameter does not determine the bush outer radius.
4. Oil heat capacity/volumetric heat capacity and thermal conductivity.
5. Cavitation pressure/plateau assumption. The earlier 500-Pa value was a timing
   assumption; it has not been assigned to these six experimental cases.
6. TC2 location for a like-for-like thermal comparison.

The package deliberately records unresolved inputs as null. No assumed pressure,
wall temperature, heat-transfer coefficient or initial eccentricity has been
substituted for an experimental input.

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
