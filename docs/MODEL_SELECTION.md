# Model selection, unsteadiness and experimental comparisons

Reviewed 2026-09-08 against the source at `f964288`. This is the current
model-selection reference; older sections of PHYSICS.md also contain historical
implementations and proposed work. The current main checkout is
`C:/Dev/pancake-physics-audit`. The original `C:/Dev/pancake` checkout remains
on `codex/ferron-fig6-flow` at `446327c`, with separate uncommitted results and
documentation. Its Ferron campaign is not a run of the new COUPLED solver.

## A steady load is not a physical acceptance criterion

Gas release, bubble compression and reabsorption can be time dependent even
under nominally constant operating conditions. Local intermittent release does
not require a similarly erratic integrated bearing force: pressure acts over
the whole bearing, and local fluctuations can cancel in the force integral.
The present model has no random forcing, bubble-number population or resolved
interface. It can predict deterministic transients from its transport,
kinetics, thermal feedback and bearing motion, but not individual nucleation
events or collapse pressure waves.

The reported improvement was timestep consistency and simultaneous mass/volume
admissibility. The former split model could have a converged pressure solve
while the subsequently updated phases occupied more than the available gap.
Its alternating frequency followed the timestep and its load did not settle
under refinement. These are reasons to reject those numerical fluctuations;
smoothness alone is not one. Backward Euler and upwind transport also damp
resolved fluctuations, so a smooth coupled curve does not prove that a physical
instability is absent.

For a proposed physical oscillation, compare time histories over equal physical
windows after startup, refine space/time and nonlinear tolerances, and check
mass, volume and energy at every step. Require convergence of frequency,
amplitude, phase and growth/decay rate for periodic motion. For irregular or
chaotic motion compare statistical quantities such as mean, RMS, spectrum and
distributions; exact pointwise agreement is not generally the appropriate test.
Record the actual pressure/speed/load forcing and the temporal bandwidth of the
measurement. Do not interpret a cycle in steady nonlinear iterations as a
cycle in physical time.

High-speed narrow-gap observations distinguish release from compression and
delayed absorption. They justify separate kinetics, not a universal stochastic
load model or transferable FW68D rate constants
(Reinke et al., *[Fluids](https://doi.org/10.3390/fluids7050159)*, 2022).

## Shared mechanical approximation

The fluid domain is the unwrapped bearing surface `(s=R theta,z)`. Pressure is
assumed uniform across the clearance. Lubrication momentum balance eliminates
the cross-clearance velocity solve; for uniform viscosity,

\[
q_s={Uh\over2}-{h^3\over12\mu}\partial_s p,\qquad
q_z=-{h^3\over12\mu}\partial_z p,\qquad U=\omega R.
\]

These are volume fluxes per unit width in full film. Cavitation changes the
mass/film-content transport described below. The physical assumptions are
Newtonian laminar flow, no slip, a thin slowly varying gap, negligible fluid
inertia and pressure equilibration across the gap. The inertia criterion is
the reduced lubrication Reynolds number, approximately
`rho U h^2/(mu ell)`, for the relevant in-plane length `ell`; merely checking
`rho U h/mu` against one is not the general criterion. Three-dimensional VTK
output does not make this a three-dimensional Navier–Stokes calculation.

Static geometry is `h=c-ex(z) cos(theta)-ey(z) sin(theta)`, with optional linear
misalignment. Moving-bearing mode integrates a rigid outer housing's position
and velocity under film force, external load, springs and damping; the shaft
axis is fixed. Fluid force is lagged between motion advances. Choosing RK4 for
housing motion does not make the fluid or the coupled fluid/structure scheme
fourth order. The fluid storage equations use backward Euler.

Elastic deformation, conjugate thermal expansion of the clearance, contact,
roughness and a general shaft/rotor structural model are not implemented by
the geometry routine. A measured hot clearance can be prescribed. Proposed
ETHD equations in PHYSICS.md are not active implementation.

## Pressure and cavitation choices

| `cavitation_model` | Actual formulation | Appropriate use and limit |
|---|---|---|
| `FULL_SOMMERFELD` | Reynolds pressure without cavitation clipping | Analytic/full-film verification, or an operating state known to remain full film. Does not predict rupture or metastable tensile-liquid physics. |
| `GUMBEL` | Solve pressure and then clip below `p_cav` | Historical half-Sommerfeld baseline. Does not conserve mass through rupture/reformation; unsuitable for cavity, leakage or phase-inventory validation. |
| `JFO` | Elrod–Adams universal content, `rho=rho0 theta`, `p=pcav+beta ln(theta)` in full film and a cavity plateau below `theta=1` | Mass-conserving rupture/reformation with barotropic full-film compressibility. Requires bulk modulus and mesh/time studies. |
| `FBNS` | Incompressible JFO complementarity, `p>=pcav`, `theta<=1`, `(p-pcav)(1-theta)=0` in the gas-free case | Useful pure-oil/JFO reference. Different numerical solution of the JFO conditions, not inherently better cavitation physics than JFO. No liquid bulk-modulus storage. |
| `COUPLED` | Pressure plus conserved oil, dissolved and free-species mass, and optional energy, solved simultaneously with gas EOS and gap volume | Current path for oil–air and oil–R410A release/compression/absorption. Homogeneous reduced model with the phase and thermal restrictions below. |

For JFO, `film_content>1` means compressed liquid under the barotropic EOS.
For gas-free FBNS and COUPLED it is a liquid volume fraction bounded by one.
Identically named output fields therefore need model-specific interpretation.

The older TWO_PHASE additions are still available for history reproduction:
JFO `VOID_COUPLED` shifts the plateau and modifies bulk modulus and occupied
volume using a lagged gas state; FBNS displaces its liquid ceiling to
`1-alpha_gas`, also lagged. They are not the simultaneous COUPLED equations.
The documented segregated pressure/volume failure makes them unsuitable as
the validated gas-bearing reference. FBNS rejects `VOID_COUPLED`.

`SCALAR_PCAV` is a prescribed plateau for the older JFO paths. A plateau fitted
to a gas-bearing experiment is an empirical input, not an oil vapor-pressure
measurement. Legacy `LOCAL_PSAT` raises that plateau to the feed bubble-point
estimate. It does not predict a spatially varying gas-cavity pressure and must
not silently raise a requested end pressure. COUPLED rejects this switch and
requires `SCALAR_PCAV` and `gas_pressure_coupling=NONE`: there, `p_cav` denotes
the oil vapor-pressure reference, and gas pressure follows its own inventory.

## What the coupled mixture model represents

There is one oil component and one dissolved/free species: AIR, R410A, or the
retained propane option. AIR and R410A cannot both be independently transported
in one run. Air is lumped into one ideal-gas species; nitrogen/oxygen
fractionation is not represented. R410A is one fixed-blend pseudo-pure species,
not separate dissolved R32 and R125 inventories.

Let `Mo,Md,Mg` be kg/m² of oil, dissolved species and free species. Then

\[
c_d={M_d\over M_o+M_d},\quad
\theta_l={M_o+M_d\over\rho_lh},\quad
M_g=\rho_g(p-p_v(T),T)h(1-\theta_l).
\]

Oil and total dissolved-plus-free species are conserved separately. Free gas
and liquid share the continuum transport velocity. Dissolved Fickian diffusion
has an equal opposite oil mass flux. Mobility uses mixture density and
`mu_m=theta_l mu_l+(1-theta_l) mu_g`, the DUKLER volume-weighted law. This is an
assumed effective rheology, not a validated formula for every bubble or oil
stripe geometry. Gas viscosity and conductivity are supplied constants.

For finite-rate exchange the dissolved balance contains `-Gamma`, where

\[
\Gamma=k(1-c_{eq})(M_d-M_{d,eq}),\qquad
M_{d,eq}=\min\left(M_d+M_g,{c_{eq}\over1-c_{eq}}M_o\right).
\]

`gas_mass_transfer_rate` is the release rate; `gas_absorption_rate` controls
absorption independently, with `-1` inheriting the release rate. Zero
absorption retains free mass but still allows EOS compression. These are
phenomenological rates without explicit interfacial area, bubble radius or
diffusion boundary layers. With release enabled, supersaturation can release
gas even from a cell with no pre-existing bubble: no nucleation barrier is
modeled. A smooth rate law need not generate erratic nucleation histories.

`coupled_phase_equilibrium=true` instead imposes
`ceq-cd>=0`, `Mg>=0`, `(ceq-cd)Mg=0`. This is the instantaneous equilibrium
limit. A useful choice criterion is the transfer time relative to residence
and forcing times: rapid exchange supports an equilibrium approximation;
comparable times require finite kinetics; negligible exchange approaches
frozen inventories. Steady flow can still have finite-rate spatial
disequilibrium. `STEADY_STATE` is not synonymous with phase equilibrium.

The vapor limit is `p>=pv(T)` with negligible oil-vapor mass. When free-species
mass vanishes it recovers JFO complementarity. With dissolved gas released,
the bubble contains free species and equilibrium oil vapor and can have
variable pressure above `pv`. They do not occupy two separate additive voids.
There is no transported oil-vapor inventory, evaporation/condensation rate,
oil-vapor latent heat, bubble surface tension, slip, coalescence/breakup,
Rayleigh–Plesset inertia, collapse shock or erosion model. The model can be
useful for averaged film mass and pressure while failing cavity topology or
high-frequency bubble dynamics. The Etsion extent mismatch remains unresolved.

R410A density and internal energy use the published Helmholtz EOS
(Lemmon, *[Int. J. Thermophys.](https://doi.org/10.1023/A:1025048800563)*, 2003).
The implemented API spans 200–500 K and up to 50 MPa, with subcritical free
fluid restricted to the stable gas branch below its dew pressure. Liquid oil
containing dissolved refrigerant may exist above that pure-fluid dew pressure.
A resolved free phase beyond the stable gas domain is rejected; metastable
bubbles and a separate refrigerant-rich liquid are not modeled. Above the
critical temperature the free phase may be dense. This EOS is substantially
better than ideal gas for R410A compression, but it supplies no oil-mixture
equilibrium or caloric calibration and no blend fractionation.

## Material options

| Setting | Meaning and caveat |
|---|---|
| `fluid_property_model=CONSTANT` | Fixed liquid viscosity/density; no dissolved species. An energy field can heat without changing load through viscosity. |
| `SINGLE_PHASE` | Actual dissolved concentration changes selected liquid properties; phase exchange disabled. Use for pure-oil thermal viscosity with zero dissolved charge, or a dissolved-only control. |
| `TWO_PHASE` | Enables release/absorption and free inventory. Use with COUPLED for the present gas-bearing target. |
| `solubility_model=HENRY` | `ceq=H(T) pg`, `H=Href exp[EH(1/T-1/Tref)]` on COUPLED. Dilute empirical law, not a nonideal oil/refrigerant flash. |
| `solubility_model=TABLE` | Measured equilibrium fraction against pressure, optionally temperature. Existing table interpolation clamps at edges; ensure all states lie within the measured domain rather than treating that clamp as extrapolation. |
| `density_model=CONSTANT` / `MIXTURE` | Fixed solution density or additive oil/dissolved partial volumes. COUPLED thermal mode assumes incompressible partial liquid volumes and rejects the extra concentration correction. No liquid thermal expansion/compressibility is derived from the R410A free-fluid EOS. |
| `density_model=TABLE` | Older concentration-only or pressure/temperature table path; rejected by COUPLED. |
| `liquid_viscosity_model=EMPIRICAL` | Andrade temperature, exponential concentration and Barus pressure factors. Requires fitted coefficients. |
| `POLYNOMIAL` | Ferron-form quadratic viscosity in the ratio of Celsius temperatures; pure temperature dependence. The Ferron oil coefficients are not FW68D coefficients. |
| `TABLE` | Older kinematic-viscosity table converted to dynamic viscosity. Saturated pressure/temperature form rejected by COUPLED. |
| `COMPOSITION_TABLE` | Dynamic viscosity against actual dissolved mass fraction and temperature, logarithmic interpolation, strict range rejection; no pressure axis. |
| Mixture `DUKLER` / `LINEAR` / `MCADAMS` | Volume-weighted arithmetic / free-mass-quality arithmetic / free-mass-quality harmonic viscosity. COUPLED currently accepts DUKLER only. None is a universal experimental validation certificate. |

FW68D lacks an identified joint calibration of solubility, viscosity, partial
volume and dissolved/free caloric data here. The reference FVC68B table is a
different oil. In addition, independently fitted Henry temperature dependence
and heat of desorption must be thermodynamically consistent; conservation by
itself does not enforce that consistency.

## Boundaries and supply geometry

* Circumference: periodic pressure, transported fields and temperature, with
  an eccentric gap that is fixed in the bearing frame for static geometry.
* Axial ends: `DIRICHLET` prescribes absolute pressure, not velocity or the
  sign of flow. Both ends can simultaneously contain inflow and outflow at
  different circumferential positions. Half-cell distances define the end
  fluxes. This is the appropriate setting for the user's two reservoirs.
* COUPLED end inflow uses that side's dissolved and free mass fractions;
  `bc_z_*_thermal=CONSTANT` supplies its measured reservoir temperature.
  Outflow transports local state. `ZERO_GRADIENT` extrapolates temperature
  even on inflow, an assumed thermal recirculation condition rather than a
  known-temperature reservoir. Open ends have zero diffusive heat flux.
* `NEUMANN` is zero pressure gradient and no axial pressure-driven flow; with
  no axial Couette transport this is a sealed/symmetry end in COUPLED.
* `INLET_OUTLET` is not needed to allow reverse flow. Legacy JFO conditionally
  applies its pressure link on inflow, while the pressure-primary legacy
  routine treats it as Dirichlet. FBNS and COUPLED reject it. Avoid it as a
  common cross-model validation boundary.
* No internal inlet: periodic film supplied/drained through its ends. A fully
  closed domain requires transient inventory evolution on COUPLED.
* Groove: a pressure patch of specified angular width over the full axial
  length. Circular hole: circular patch in unwrapped `(s,z)`, size is radius
  in metres. These constrain film cells without adding resolved groove depth,
  feed-pipe resistance or entrance jets.
* COUPLED patch: net required supply mass follows the local balance. Incoming
  mass carries configured composition and enthalpy; negative supply drains
  local mixture. Its temperature is solved, not pinned. Optional fresh-oil
  temperature is used on inflow; otherwise local temperature is extrapolated.
* Legacy thermal patch: `PINNED` fixes fresh-oil temperature; `FIXED_LAMBDA`
  imposes an empirical fresh/recirculated mixture; `FLOW_WEIGHTED` uses supplied
  fresh flow and estimated upstream Couette flow. Pressure determines mass
  source separately. The prescribed fresh flow used for heat mixing need not
  equal the hydraulic source diagnostic. COUPLED rejects the two mixing
  overrides; its required `PINNED` enum value does not pin temperature there.

## Exactly what changed in the energy equation

The older `Energy::solve` advances a liquid-temperature equation of the form

\[
\rho_l c_{p,l}h_l(\partial_tT+\bar{\boldsymbol u}\cdot\nabla_sT)
=\nabla_s\cdot(k_lh_l\nabla_sT)+\Phi-Q_w.
\]

`h_l=h` under FULL_FILM, or liquid-weighted under STRIATED in cavitated cells.
A divergence correction implements material-temperature advection. That is
not intrinsically wrong for an incompressible, constant-composition liquid
with consistent continuity and sources. It is not a general mixture-energy
balance when phase inventories, composition and specific energies change.
Pressure, properties and temperature are segregated; transient mode takes one
such sweep, while steady mode iterates the fields.

COUPLED advances `Ucal=Mo uo+Md ud+Mg ug` in J/m². In continuous balance form,
with internal reservoir sources made explicit,

\[
\partial_t\mathcal U+\nabla_s\cdot F_u
=-p(\partial_th+\nabla_s\cdot q)+\Phi-Q_w
+\nabla_s\cdot(k_mh\nabla_sT)
+\sum_i S_i(u_{i,in}+pv_{i,in}).
\]

The [discrete derivation](COUPLED_ENERGY.md) specifies the face pressure work,
diffusion enthalpy and signed draining sources. The practical differences are:

1. Store the actual component energies instead of a changing liquid capacity
   multiplied by temperature alone. Conserve old/new inventories at the same
   physical times.
2. Transport energy with the same discrete component mass fluxes used by the
   mass balance, including independent pressure-flow and Couette donors.
3. Include compression/expansion work. It matters for bubbles even when their
   mass does not change; incompressible single-phase flow has a different limit.
4. Distinguish dissolved and free energy. `uo=cp,o (T-Tref)`,
   `ud=cp,d (T-Tref)-Delta udes`; free internal energy comes from ideal-gas cv
   or R410A EOS. Positive specified desorption energy gives release cooling
   and absorption heating without a second, double-counted latent source.
5. Use exact discrete pressure-flow power `-qp Delta p`, including the axial
   half-cell and pressure-patch neighbourhood, rather than the older gated
   cell-gradient estimate. Couette heating remains `mu_m U²/h`.
6. Solve temperature in the same nonlinear system as pressure and species;
   solubility, gas density and viscosity respond to the accepted temperature.
7. Apply directional supply enthalpy and prescribed wall heat removal,
   `Hj(T-Tj)+Hb(T-Tb)`. This path has no conjugate solid calculation.

The earlier timing repair also removed an extra physical advance: an old
15-ms request with 10-ms steps heated for 20 ms. This was separate from
changing energy physics. The current thermal/air comparison is a synthetic
fixed-eccentricity case; its 7.12% load reduction cannot be transferred to a
load-controlled Ferron experiment. At prescribed load, thermal feedback instead
changes equilibrium eccentricity, attitude, minimum gap, pressure and flow.

That 7.12% compares thermal versus isothermal runs of the new solver. It is
not an old-versus-new thermal-solver comparison and does not measure improved
agreement with an experiment.

## The older cross-film model and Ferron

`cross_film_thermal_model=POLYNOMIAL` reconstructs a quadratic or quartic
`T(y)` whose mean equals the solved film-average temperature. The quartic uses
one midplane collocation, a Couette-only approximation for its shape equation
and a bounded convective contribution. These are modeling approximations;
the profile is not a resolved cross-film energy PDE. It evaluates viscosity
integrals F0/F1/F2 and generalized Dowson transport. For the reconstructed
viscosity, `q=-F2 grad p+U F1/F0` and dissipation `U²/F0+F2|grad p|²` replace
the uniform-viscosity expressions. Resolving different wall/core viscosities
can affect force even at the same film-average temperature.

`bush_thermal_model=LOCAL` supplies a radial thermal resistance per cell.
`CONDUCTION_2D` solves steady `(r,theta)` bush conduction using axially averaged
film coupling. It redistributes heat around the circumference but has no bush
axial conduction or end-face cooling. The shaft has one global temperature
from a net heat balance. Outer bush cooling can be convective, fixed-temperature
or adiabatic. In main, the ambient reference is shared; the separate Ferron
branch adds independent bush-outer and shaft-ambient temperature inputs.
Those extra scalar inputs do not add missing spatial solid physics.

This legacy cross-film path works with JFO/FBNS; COUPLED currently rejects it.
FULL_FILM/STRIATED are older cavity heat-capacity, conduction, shear and
transport choices. STRIATED weights liquid-bearing contributions by liquid
content and is useful for a pure-oil cavity sensitivity study. It still does
not resolve individual oil streamers or supply gas thermodynamics. COUPLED
computes its own phase inventories and volume-weighted properties and rejects
STRIATED, so these choices cannot be freely combined.

The paper is *A Study of the Thermohydrodynamic Performance of a Plain Journal
Bearing Comparison Between Theory and Experiments* by Ferron, Frene and
Boncompain. It uses a three-dimensional film thermal calculation, a
three-dimensional bush conduction equation, and an axisymmetric `(r,z)` shaft
with 43°C axial ends. The outer bush cylinder has prescribed temperature;
80 W/(m² K) belongs to the bush axial faces, not that outer cylinder. Groove
recirculation enters through a mixing balance. Measurements are steady after
thermal equilibration, not high-speed bubble/load measurements
(Ferron et al., *[J. Lubr. Technol.](https://doi.org/10.1115/1.3254632)*, 1983,
pp. 423–426).

The thermocouple junctions were typically within 0.5 mm of the bore, with
reported temperature accuracy ±1°C. Compare predicted temperature at the
sensor location if available; bearing-interface temperature is a closer
observable than film-average temperature but still differs by sensor depth.
The paper's Fig. 6 theoretical peak is itself about 4°C below experiment.
Matching that theoretical curve is distinct from matching the measurement.
The paper also distinguishes cold clearance near 145 µm from measured hot
clearance near 152 µm; eccentricity must be interpreted with its clearance
(Ferron et al., *[J. Lubr. Technol.](https://doi.org/10.1115/1.3254632)*, 1983,
pp. 425–427).

The previous audit's 180×48 Ferron smoke case was unchanged across the
inventory/time repair, at saved precision. That was not a new COUPLED
Ferron validation. No new paired Ferron run was conducted for this explanatory
review. The separate branch's Fig. 6 archive illustrates model sensitivity:

| Archived thermal choice, 4000 rpm / 6000 N | Equilibrium eccentricity | Maximum film-average temperature (°C) | Hydraulic source (L/min) |
|---|---:|---:|---:|
| Depth averaged | 0.5453 | 49.990 | 8.454 |
| Cross-film adiabatic | 0.6575 | 63.086 | 10.323 |
| Cross-film walls at 40°C | 0.5097 | 45.959 | 7.670 |
| Cross-film / 2D bush, fixed outer temperature | 0.6127 | 54.800 | 8.959 |

These are read-only archived values from
`C:/Dev/pancake/cases/01e_ferron_terra_repro/fig6_analysis.json`, not fresh
measurements, new runs, or thermocouple temperatures. All four used
7.837 L/min as the digitized fresh-flow input for the thermal mixing model.
That differs from the pressure-derived source in the last column. A model can
therefore match a temperature with a fresh-flow mixing input while its
hydraulics predict a different supply. That disagreement needs explicit
assessment of groove/plumbing approximations, not a concealed retuning.

## Validation choices for discussion

| Question | Initial model choice | Evidence required before a validation claim |
|---|---|---|
| Is the hydrodynamic pressure discretization correct? | FULL_SOMMERFELD in a full-film analytic case, then gas-free FBNS and JFO | Pressure, force and flux error under mesh refinement; correct squeeze storage for transients. |
| Does the vapor/JFO free boundary conserve mass? | Gas-free FBNS/JFO; COUPLED with zero dissolved/free charge as a limit check | Independent rupture/reformation and mass-flow solution. JFO does not validate oil-vapor evaporation kinetics or collapse. |
| What explains Ferron pressure and temperatures? | SINGLE_PHASE with zero dissolved charge, measured oil mu(T), FBNS/JFO, ENERGY_EQUATION; compare depth-averaged and polynomial/conjugate paths | Match load, hot clearance, supply/mixing, end conditions, thermal boundaries and sensor observable independently. The reduced solid models limit the claim. |
| Do pressure reservoirs handle reverse flow? | COUPLED with DIRICHLET at both ends and independently specified reservoir state | Signed end fluxes, analytic transport/energy cases, conserved inventories; test groove, hole and absent feed separately. |
| Is air-cavity compression represented? | COUPLED, ideal AIR, fixed free mass and zero absorption control | Pressure-volume response and energy under compression before fitting phase-transfer rates. |
| Is gaseous cavitation captured in a flooded bearing? | COUPLED, AIR, actual charge, equilibrium and finite-rate comparisons | Etsion-type cavity pressure/extent, independent charge/rates, mesh and detection-threshold sensitivity. Current extent comparison fails. |
| Is release/absorption timing captured? | COUPLED finite-rate transient with measured pressure forcing | Independently characterized charge and kinetics; time-resolved void/pressure data and converged temporal statistics. No claim about resolved bubble topology. |
| Is FW68D–R410A heating predicted? | COUPLED R410A with energy after mixture-property calibration | Consistent solubility, concentration-viscosity and caloric data; measured thermal boundaries and temperature observable; full phase-domain checks. |
| Does a resolved refrigerant phase condense or fractionate? | No current complete model | Add condensed refrigerant and multicomponent equilibrium/kinetics before claiming support. |

The choice should follow the observable and experimental regime. The present
code has more developed pure-oil wall-temperature modeling on one path and
more consistent mixture mass/energy modeling on another. Combining those
capabilities requires a derived mixture cross-film/conjugate closure and its
verification; selecting both existing switches is not that derivation.

## Simulation workflow

For all subsequent campaigns: edit source locally, compile locally and run
short relevant checks, synchronize the selected source/configuration snapshot
to an isolated remote directory, recompile there, execute the cases remotely,
then copy results and provenance back locally for analysis and plots.
Record source revision plus dirty-source hashes if applicable, input and
executable hashes, compiler/PETSc/MPI versions, rank count, exit status,
convergence diagnostics and timings. Keep old results immutable and retain
failed cases. Do not use a failed case's last frame as its completed result.
This is source/result synchronization, not authorization to push Git history
or overwrite another checkout. This review changes documentation only and
does not launch another CFD campaign.
