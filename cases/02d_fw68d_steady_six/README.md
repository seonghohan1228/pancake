# Six steady FW68D comparisons: A–F

This campaign replaces the earlier 24-point transient plan with the six cases
selected by the user. It solves steady pressure, cavitation, energy and bearing
position together in C++. Python updates only the exterior cooling coefficients.
The package now includes the groove, housing and still-oil
cooling. The user supplied 1 atm reservoirs and exposed shaft ends and
authorized ordinary/similar material properties. The original batch is preserved as a failed/stopped baseline. The September 28
repair adds conservative feed energy and physical acceptance checks; the new
full-resolution campaign is being validated before accepting equilibria.

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
- Film mesh: 720×192; housing radial/axial cells: 16/16; film layers: 16.
- `solution_mode=STEADY_STATE`, `motion_model=STEADY_EQUILIBRIUM`, `omega_ramp_time=0`.
- Pure-FW68D `SINGLE_PHASE`, FBNS cavitation, STRIATED cavity shear, film energy,
  finite-volume cross-film temperature and `CONDUCTION_3D` housing conduction.
- Retain dt=20 us as configuration metadata. Steady Reynolds/energy storage
  terms are disabled and bearing motion is not integrated.
- Inner pressure/temperature coupling: converged FBNS pressure/complementarity
  residuals below 1e-6 and normalized temperature change below 1e-6,
  maximum 2000 sweeps, outer field relaxation 1,
  and normalized global energy defect below 1e-6. Pressure field change is
  recorded without another stopping threshold. An iteration cap or restored best
  iterate without convergence is a failed trial.
- Native equilibrium: update bearing position within the pressure/energy loop
  until total fluid force balances the applied load vector to 0.2%. Numerical
  position steps are limited to 0.05 clearance, with pressure-only position
  probes separated by 0.001 clearance. No physical bearing motion is integrated.

Eccentricity and attitude are solved outputs, so measured values are not needed
for this load-controlled calculation. The user confirmed an upward push from
below: applied load (0,+W), actuator contact at the bottom (270 degrees), and
the groove at the top (90 degrees). The groove is opposite the actuator contact,
not opposite the applied force vector. The fluid-on-bearing target is (0,-W). `load_angle_deg` currently does not rotate the solver force target;
the native solve uses explicit force-vector components in the chosen frame.
Bearing position changes are numerical equilibrium corrections. The journal
centre stays fixed, and the journal rotates at the specified speed.

`steady_setup.json` records inputs and acceptance. `run_steady.py` generates
native configurations and retains each cooling-coefficient evaluation. It
uses the final solved position and validated field checkpoint as the next
coefficient evaluation's initial guess. There is no external load search.

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
The campaign sets `liquid_specific_heat=2010` J/(kg K) and updates volumetric
capacity with the local liquid density. The legacy `rho_cp` value remains a
reference/default for other configurations.

## Specified groove and housing

The user specifies a 5 mm groove opposite the applied load, at the same
reservoir conditions as both bearing ends. Treat the width as a circumferential
arc at the 43.085 mm bore and use the native full-length axial groove:
width = (0.005/0.0215425) * (180/pi) = 13.2983125 degrees; center = 90 degrees.
At 720 theta cells the native cell-center mask selects 26 columns (13 degrees,
4.88784 mm bore arc). Retain the nominal 5 mm input and report this discretization.
Groove pressure equals the common end pressure; its supplied temperature is
case TC0. `FLUX_BALANCED` solves groove film temperature, using the signed
net FBNS mass flux to supply oil at TC0 or drain oil at local temperature.

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
Keep the supplied FW68D viscosity and density laws. Film energy uses
rho_cp(T)=rho(T)*cp. These cp/k values are not FW68D measurements.

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

The driver updates the housing and shaft coefficients after each native
load/thermal solve, using solved mean surface temperatures and relaxation. The mean
housing outer temperature is read from the resolved outer-surface field.
Accept coefficient changes below 1e-3 relative and
wall-temperature changes below 0.01 K, alongside both inner coupling residuals
below 1e-6. A fixed-coefficient timing trial is explicitly labelled unclosed.

L/D=0.47 is short, so this long-cylinder lateral-surface estimate is approximate.
The two annular end faces have about 79% of the outside cylindrical area.
The user confirmed both are exposed to surrounding oil at TC0. They now use
local natural convection coupled to three-dimensional housing conduction.
A vertical-plate analogue with characteristic height 85 mm estimates their
heat transfer (Churchill and Chu, *[Int. J. Heat Mass Transf.](https://doi.org/10.1016/0017-9310(75)90243-4)*, 1975).
Confinement and mounting contacts remain unresolved; these correlations
estimate cooling rather than supplying measured coefficients.

## Remaining comparison information

TC1-3 are at r=26 mm and the axial midplane, at 240/270/300 degrees.
Only the logged pressure sensor location remains unspecified for comparison.
Full-length groove, horizontal housing, isothermal shaft and analogous
materials are explicit modeling assumptions. Revisit cooling sensitivity when
assessing experimental agreement; chosen coefficient/property variations are
sensitivity scenarios, not measured confidence bounds.

## Execution and acceptance

`run_steady.py` renders the exact cases with native equilibrium for the
fluid-on-bearing target (0,-W). The F0/F1 wall-stress correction is implemented
and checked against variable-viscosity and uniform-viscosity limits. Every
solver attempt retains its config, log, hashes, diagnostic table and fields;
failed or capped inner/thermal/force solves remain failed.

```sh
python run_steady.py --self-check
python run_steady.py --solver /absolute/path/pancake --ranks 2 \
  --work /absolute/path/new-results --cases F,E,D,C,B,A,E_low,E_high
```

A frozen-position timing call is `--trial F .15 345 --thermal-iterations 1`.
It is not load-balanced or externally thermally converged. On Open MPI only,
`--rankfile` assigns one rank to each chosen physical core. The workstation
left ranks unbound with `--cpu-list`, so that option is not used. The user-requested resolution remains 720 by 192 throughout.

Follow the repository workflow: local build and short checks, isolated remote
source/configuration snapshot, fresh remote build, converged equilibrium runs,
and return of evidence. Preserve every failed evaluation and record source,
input and executable hashes, rank layout, residuals and return codes. Existing
coarse transient timings do not establish an optimum for these steady solves.

A completed comparison must report pressure/thermal convergence, native
vector load closure and exterior cooling-coefficient closure separately;
finite fields, positive gap and liquid filling;
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

### Pinned-groove thermal limitation

The current model pins the groove depth-average temperature to TC0 while
applying the sheared quartic film profile and conjugate wall closure there.
This can create a sub-ambient bore patch in the groove and neighboring cells,
even with nonnegative bulk temperature rise and balanced solid wall fluxes.
The fixed-F diagnostic showed a 0.218 K groove undershoot and 0.069 K in two
adjacent columns. These wall values are a model artifact, not physical cooling
below the reservoir. Keep the raw fields, report their minima, and avoid using
groove-region wall temperatures as experimental predictions. Its influence on
friction/cavitation remains a sensitivity question; it is not assumed negligible.


## September 28 physical acceptance and continuation

The original PINNED groove produced subambient bore/profile temperatures and
undamped quartic updates stalled several inner solves. Those results are
retained but rejected. The replacement uses pressure-derived mass/enthalpy
exchange, with no imposed mixing flow, and computes the exact internal profile
minimum/maximum. No clipping is used to pass acceptance.

Each trial must have a clean solver return, mean-temperature change
and FBNS residuals below 1e-6, relative mass-rate residual
at most 1e-8, and normalized global energy defect at most 1e-6. Require finite
positive gap/viscosity/density, pressure at or above the cavity floor within
numerical tolerance, film content in [0,1], and film/wall/profile temperatures
at least TC0 minus 0.001 K. External coefficient closure remains 1e-3 relative
and 0.01 K. A final case additionally needs vector force error below 0.2%.
Fixed-position trials also require pressure field change below 1e-6.

A fixed-position `--trial` uses `motion_model=STATIC` and cannot pass as a
native equilibrium. Production runs require the final `Equilibrium: CONVERGED`
marker and independently verify the VTK total-fluid-force vector against the
applied load. Native position is read from VTK `xB`; reported eccentricity and
attitude describe that final position. `seed_eps` and `seed_psi_deg` preserve
the configuration's initial guess.

The driver accepts cached runs only with matching source, executable, driver,
input, rank count, rankfile and native-mode provenance. It rechecks physical
and convergence evidence before reusing a same-case checkpoint. Legacy
external-search results remain preserved but cannot masquerade as native results.

### September 28 native equilibrium driver update

Removed the external position Newton loop and its full thermal solves for
base, finite-difference and backtracking geometry evaluations. Each exterior
coefficient update now runs one integrated native steady solve. The existing
housing and shaft correlations, coefficient tolerances, case inputs and full
resolution are retained. The driver self-check exercises final-position signs,
invalid position rejection, independent/native force failures, cache provenance,
trial isolation and reuse of the solved position for the next coefficient update.


## September 28 load-direction correction

The earlier downward-load campaign, including numerically converged Case E,
was invalidated when the user confirmed a push upward from below. Those
outputs remain preserved as historical evidence of the wrong physical setup.
The corrected driver must use residual (Fx,Fy+W), keep the groove at 90 degrees,
and solve new bearing positions. Previous friction and temperature plots are
not finalized experimental predictions. No physical field may simply be
rotated: the groove remains fixed while the force reverses.


### User-requested inner-pressure stop, 28 September 2026

The user explicitly selected the inner pressure-equation residual, rather than
the pressure change between coupled CFD iterations, with threshold 1e-6.
All seven active cases in revision 23fcbe7 were terminated after verifying that
they had already met this criterion. B had previously failed its thermal-profile
physical check. Saved fields, checkpoints and logs were preserved; the solver
has no stop-and-export handler, so unsaved intermediate fields were not written.
The local campaign now renders fbns_tol=1e-6. This does not change the separate
force, thermal, global conservation or physicality requirements for accepting
a final experimental prediction, and it does not automatically restart cases.

### September 28 recovery driver changes

The recovered campaign uses `FINITE_VOLUME` with 16 layers and native
`steady_pressure_inner_only=true`. It continues after the inner FBNS residual
reaches 1e-6 until temperature, force, profile, global mass/energy and exterior
cooling closure pass. Pressure field change remains in the diagnostic record.

`steady_output_interval=50` retains intermediate fields and a restart file in
`results/iterations/`. At a driver runtime deadline, `results/STOP` requests a
completed-sweep export; a 120-second grace period precedes process-tree cleanup.
Stopped outputs keep their failed status and cannot become accepted results.

`--initial-state-root /path/to/prior/production` permits a preserved same-case
state to initialize a fresh solve after changes to the numerical algorithm.
The driver checks the original configuration hash, experiment/property-table
hashes, geometry, load direction and physical settings, reads the saved bearing
position, and records the prior status and state hash. This option imports an
initial guess only; accepted-result caching still requires matching current
source, executable, driver, inputs and MPI placement. Native I/O validates all
state values and initializes new layer fields from legacy mean-temperature seeds.

TC1/TC2/TC3 are comparison probes at radius 26 mm, axial midplane and angles
240/270/300 degrees. These positions do not set the thermal boundary conditions.

The corrected source snapshot was rebuilt on afml-r6 after 59 local tests
passed. All 26 selected remote tests and both driver checks passed. A full-mesh
Case B iteration took 61.96/38.04/27.26 s on 2/4/8 ranks. The campaign therefore
runs six primary cases on four cores each and the two E temperature bounds
on two cores each, covering all 28 physical cores without a queued case.
These short timing pilots were intentionally capped and are not final results.

### Completed resolved campaign, September 28

Remote snapshot: `/home/user/pancake-runs/fw68d-resolved-20260928-e414b50f`.
All A-F and E-low/high equilibria now pass inner pressure, force, temperature,
layer/solid change, mass, energy, exterior coefficient and temperature-floor
checks. Primary cases took 18, 22, 23, 24, 32 and 24 total coupled iterations
respectively, including all cooling-coefficient passes. E-low/high were moved
from two to eight freed cores using saved-state recovery; their earlier stop
records are deliberate reallocations, not accepted CFD solutions.

The maximum load residual is 0.1445%, maximum relative energy residual
5.10e-10 and maximum relative mass residual 7.43e-10. The fully coupled B
check with 32 film layers changes torque by 0.0048%, peak pressure by 0.0226%
and sensor temperatures by at most 0.0021 K. Independent solid conduction
reproduces saved A/B fields within 5e-6 K; 16-to-32 radial/axial refinement
at fixed interface heat flux changes sensor temperatures by <=0.0024 K.

Results and plots are in
`C:/Dev/pancake/results/fw68d_physical_20260928/resolved_results/RESULTS.md`.
The report includes conditions, all experimental points, sensor errors,
pressure/temperature/filling contours, convergence histories and the full
raw-output archive. Numerical closure is established, but friction errors
range from -38.4% to +103.2% and E is 4.1-4.8 K hotter than experiment.
F has no cavity: its minimum absolute pressure is 24.60 kPa, above p_cav=0.
The nominal material/cooling and ideal groove assumptions remain model
limitations; no coefficient or torque offset was tuned to the measurements.

## Independent thermal boundary scenarios

See [the sensitivity campaign](SENSITIVITY_CAMPAIGN.md) for the explicit five-temperature setup interface, eight-run proposal, archived baseline match, commands and acceptance. With no scenario selected, historical baseline behavior is preserved. Older torque-error prose above uses superseded experimental subtraction; the campaign uses `new_exp_torque_nm`.
