# Conditional LGE gas release comparison (2026-10-19)

The active request now selects P01--P10, one exact measurement window at each
speed from 200 to 1100 rpm. The 200, 500 and 1100 rpm selections differ from
the original seventeen-case matrix. The new campaign freezes their measured
load, chamber pressure, TC0 and comparison data directly, with no interpolation
of results from nearby earlier windows. Its twenty gas-off/on jobs are followed
by four tighter outer/load checks at P01/P10 and two 32-layer checks at P01.

The user requested `coupled_nonlinear_tolerance = 1e-8`. This now controls the
inner nonlinear stop and native/raw verification consistently. The separate
1e-13 override described in the historical notes below is superseded for this
campaign. Native defaults outside this campaign remain 1e-10. Outer thermal,
load, mass, energy and physical admissibility gates remain enabled. Earlier
attempts are retained; local coarse probes are not production results.

The original follow-on prepared matched steady COUPLED calculations for the 17
selected M01–M17 operating windows. It is a conditional model comparison,
not an experimental validation. Input preparation and all 34 configuration
checks are complete. Matched 48×16 pilots using relaxation 0.3 passed strict
raw-field acceptance. The isolated remote build passed all 65 native tests,
and the full 720×192 campaign is running. These local pilot checks do not
provide campaign-level gas-effect results; remote cases must pass the gates
below before they enter the comparison.

The fifth remote attempt exposed an inherited two-hour wall-clock limit per
native solve. M02 gas-off was still reducing its coupled residual when the
driver stopped it; the three gas-on pilots also reached that limit. Their
failed outputs are retained. The gas worker now defaults to no wall-clock
deadline. `--timeout SECONDS` remains available for explicitly bounded checks;
`--timeout 0` selects the default. Each native record stores the selected
deadline. Iteration caps and every physical, load and conservation gate remain
unchanged. The two progressing gas-off pilots from that attempt continue in
their original snapshot while the retry uses disjoint free cores.

## Matched branches

Both branches use the same geometry, measured operating conditions, absolute
reservoir/groove pressure, full digitized viscosity chart `nu(c,T)`, density
law, `COUPLED` steady finite-volume thermal model, mesh, wall/feed conditions,
solid conduction, cooling and load-balance procedure. Coupling relaxation is
0.3 in both branches. Dissolved composition is initialized uniformly at 0.30
across the field and supplied at both axial ends; feed free-gas mass fraction
is zero. Composition and phase fractions are uniform through the gap;
equilibrium uses film-mean temperature while layer energy and viscosity are
resolved separately. The branch
differences are limited to:

| Branch | Fluid property model | Phase equilibrium | Interpretation |
|---|---|---|---|
| `gas_off` | `SINGLE_PHASE` | disabled | Composition transports without release |
| `gas_on` | `TWO_PHASE` | instantaneous local equilibrium | Release follows local bubble-pressure complementarity |

The gas-on minus gas-off result can attribute release effects within this
shared conditional model. It cannot be compared as a controlled gas effect
against the earlier 02e FBNS/STRIATED fixed-chart-30 cases, which use a
different solver and cavitation structure.

## Conditional material assumptions

The chart metadata leaves pressure and concentration bases unspecified. Under
the user-authorized scenario, the pressure coordinate is interpreted as
absolute pressure and the concentration coordinate as refrigerant mass
fraction. The approximately 30% charge observation does not independently
establish that basis; no fitted or inferred transfer rate is used. The full
viscosity surface is evaluated without clamping or extrapolation. The
pressure-temperature chart is inverted at each local composition and
temperature for the first supported bubble-pressure crossing. Unsupported
bubble points fail the run. A free R410A phase must also remain within the
stable EOS domain at every film layer and journal/bearing wall.

The same pure-oil density proxy is retained in both branches:

\[
\rho(T_C)=954.9-0.7284T_C\quad\mathrm{kg/m^3}.
\]

This law derives from advanced PVE VG68 formulated for R32, not an identified
FW68D/R410A mixture density (Zeng et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2021.12.009)*, 2022). The reviewed PVE/R410A study distinguishes actual concentration from saturation solubility, but does not establish FW68D coefficients (Fukuta et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2009.09.015)*, 2010). No verified exact-grade mixture density or heat-of-solution data were found in the bounded source review.

The caloric input record is
[`conditional_caloric_inputs.json`](conditional_caloric_inputs.json). The full
source review and checked manuscript excerpts are retained with the campaign
evidence in `results/lge_mixture_density_gas_20261019/property_review/` in the
parent working checkout.
It uses pure saturated-liquid R410A at 20 °C as a dissolved-phase caloric
analogue with constant heat capacity 1656.9291 J/(kg K), desorption internal
energy 192486.2475 J/kg, and zero excess mixing energy; free-gas properties
use the native real-gas EOS, with gas conductivity 0.0129583 W/(m K) and
viscosity 1.26525×10⁻⁵ Pa s at the 293.15 K, 100 kPa reference. Oil heat
capacity is 2010 J/(kg K), liquid conductivity 0.13 W/(m K), and steel
conductivity 15 W/(m K). These are explicit proxies. The R410A analogue is
not measured FW68D heat of solution, and none of the inputs is fitted to
torque or temperature.

A future density comparison could use additive liquid specific volumes,
`1/rho_mix = (1-w)/rho_oil + w/rho_R410A_liquid + v_excess(T,p,w)`, with
refrigerant mass fraction `w` and excess specific volume `v_excess`. Setting
`v_excess = 0` is a conditional ideal-volume model. It needs a confirmed
concentration basis and exact-grade oil/mixture density measurements before
it can be treated as a FW68D correlation. The R410A gas density cannot serve
as the dissolved-liquid density. Fukuta et al.'s ideal-volume conversion
uses saturated-liquid refrigerant density and does not establish an
actual-pressure partial molar volume model.

## Thermal method and acceptance

The steady thermal extension alternates the coupled pressure/species solve
with a finite-volume cross-film energy solve and the resolved housing
conduction network. It updates layer temperatures, local `nu(c,T)`, phase
properties, wall temperatures and heat fluxes until the coupled solution and
thermal profile converge. The discrete energy balance includes component
face enthalpy, pressure work, viscous and pressure-flow dissipation, wall heat
flux, feed enthalpy, and solid conduction. This is a steady outer iteration;
it does not add cross-film temperature layers to the transient Newton state.
See [the coupled energy method](../../docs/COUPLED_ENERGY.md) and
[the coupled conservation model](../../docs/COUPLED_CAVITATION.md).

The driver requires zero solver return code, native pressure convergence,
no exhausted cap or restored-best iterate, pressure residual below
\(10^{-10}\), and coupling, volume, pressure/species/thermal, profile and
global energy residuals within their requested gates. Oil and species balance
must satisfy the tighter of configured tolerance and \(10^{-8}\). Independent
integrated load, finite positive physical fields, strict table support and
stable free-gas EOS checks must also pass. Failed attempts and outputs are
preserved; no solver result can be accepted merely because it reached a
nominal iteration cap.

Both branches use GMRES with overlapping additive Schwarz and local LU
subdomain solves (`-coupled_sub_pc_type lu`). The default local ILU(0)
reached its 1000-iteration cap on the production grid. Controlled local
comparisons reproduced the same first-step equilibrium with LU, and the
full-grid remote probe reached the unchanged steady pressure gate directly.
The solver command and its nonlinear/linear progress logs are retained.
See the [PETSc ASM options](https://petsc.org/release/manualpages/PC/PCASM/).

The inner SNES absolute tolerance is `1e-13` in both normal and focused runs.
At the approximately 1.4 MPa pressure scale, the former `1e-10` stopping
tolerance permitted a negative-pressure error larger than the unchanged
`1e-6 Pa` admissibility bound. Focused runs additionally tighten force,
thermal and exterior-cooling closure. Gas-off uses the algebraically
equivalent explicit free-species balance, avoiding subtraction of nearly
equal dissolved and total-species equations. The original total-species
residual remains in final acceptance.

The position Jacobian uses `steady_position_fd_step = 0.0001`, a 4.25 nm
displacement at the nominal clearance. The former 42.5 nm difference was
several percent of M02's minimum gap and produced a non-descent Newton
direction on the reduced grid. The smaller step reached the unchanged
0.2% force criterion without changing geometry or measured load.

Disposable load-position probes and line-search candidates use direct steady
pressure solves. A failed candidate is rejected without pseudo-time recovery;
the base pressure state retains that recovery. The verified density-only CFD
shaft positions can seed both matched branches identically. These are numerical
initial guesses, not fixed geometry or reused thermal/species checkpoints;
each branch must independently close its prescribed measured load.

For the steady COUPLED path, the outer controller also damps each layer's
inverse viscosity before computing all Dowson integrals. This addresses the
observed alternating liquid-fill/viscosity iteration. The first profile uses
the physical law directly; later updates blend positive old and current
profiles. Pressure mobility, velocity weights and shear use those same
integrals. Conservative phase inventories are not relaxed. A separate
`viscosity_residual` compares the used profile with the unrelaxed physical
law at the current fill, composition and layer temperatures. It must satisfy
the coupling tolerance, and the campaign independently recomputes it from
saved layer fields. Small relaxation cannot hide a constitutive mismatch.

Both branches enable `coupled_fbns_initialization`. One preliminary gas-free
FBNS pressure/film-content solve at the prescribed starting position,
composition and temperature initializes the conservative COUPLED inventories.
It supplies a cavitating initial guess; the subsequent COUPLED equations and
all acceptance gates are unchanged. A fully flooded starting guess stalled
M02's Newton iteration even with converged linear solves, a different shaft
position and trust-region globalization. The isolated cavitating-seed check
completed three strict coupled load/thermal updates in both branches. Those
capped diagnostics are preserved and are not accepted campaign results.

## Inputs and next campaign step

Evidence provenance and numerical assumptions are collected under
`results/lge_mixture_density_gas_20261019/property_review/`. The case input
builder derives from the preserved 02e evidence/preflight and writes 34
matched configurations, along with strict-domain and provenance manifests.
The latest native Windows suite passes 66 checks. The prior remote Linux
snapshot also passed 66; each new snapshot must pass its complete suite before launch.
The 40-job campaign includes six paired pilots, 28 remaining cases, four
tighter-tolerance checks (M02/M16), and two 32-layer checks (M02). Four MPI
ranks are bound to separate physical cores per case, with at most 24 cores
across concurrent campaigns. Source/input/executable hashes, runtime versions,
return codes and convergence diagnostics are retained with all attempts.

The mass-fraction/absolute-pressure interpretation, pure-oil density proxy,
caloric values and instantaneous-equilibrium assumption remain conditional.
Condition stability and the measured approximately 30% charge do not prove
an in-run steady composition or validate the mixture model.
