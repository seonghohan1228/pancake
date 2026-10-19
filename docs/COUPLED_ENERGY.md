# Energy, compression and phase transfer

`cavitation_model = COUPLED` with `temperature_model = ENERGY_EQUATION`
solves temperature in the same Newton system as pressure and the oil,
dissolved-species and free-species inventories. The gas EOS, solubility and
viscosity therefore use the accepted film temperature. The legacy liquid-only
energy solver is not called after this solve.

## Conserved energy

Let \(M_i\) be mass per bearing surface area, \(u_i\) specific internal
energy, \(h\) the gap, and \(q_f\) signed volume flow through a face. The
stored energy is \(\mathcal U=\sum_i M_i u_i\), in J/m². The discrete
backward-Euler balance is

\[
\mathcal U^{n+1}-\mathcal U^n+
{\Delta t\over A}\sum_f F_{u,f}^{n+1}
=-p^{n+1}\left[h^{n+1}-h^n+
{\Delta t\over A}\sum_f q_f^{n+1}\right]
+\Delta t(\Phi-Q_w)-{\Delta t\over A}\sum_f Q_{k,f}
+\sum_i\Delta M_{i,\rm supply}(u_{i,\rm supply}+p v_{i,\rm supply}).
\]

Face sums are outward from the cell. \(F_u\) transports the component
internal energies with exactly the pressure and Couette component mass
fluxes used by conservation. The two contributions may have opposite donors.
Dissolved diffusion has an opposite oil diffusion flux and transports the
difference of their partial enthalpies. There is no independent guessed
velocity field in this equation.

The term in square brackets gives compression/expansion work. Without it,
a sealed bubble changes volume while remaining artificially isothermal.
An injected mass carries enthalpy because the reservoir performs flow work.
For a draining pressure patch, the outgoing enthalpy uses the local film
state. A specified supply temperature is used only for actual inflow; the
temperature of a groove cell is not pinned to the supply temperature.

All old inventories, energy and gap remain fixed during a physical time
step. Iterating or calling the solve twice does not heat the fluid twice.
Steady mode removes storage from the species and energy equations together.

## Caloric properties and heat of solution

The liquid components are treated as incompressible, with additive partial
volumes and constant partial heat capacities. A concentration-only density
correction that breaks those additive volumes is rejected in thermal mode.
`LINEAR_TEMPERATURE` changes liquid volume while retaining this caloric
approximation; it is not a complete thermal-expansion thermodynamic closure.
Its incoming pressure-work term now uses the density at the supply temperature,
as verified by the cold-feed balance in the
[September 15 audit](SOLVER_AUDIT_2026_09_15.md). The validity policy for using
this density law in thermal calculations remains an explicit model discussion.
The reference convention is

\[
u_o=c_{p,o}(T-T_*),\qquad
u_d=c_{p,d}(T-T_*)-\Delta u_{\rm des,*},
\qquad u_g=u_{g,\rm EOS}(p_g,T)-u_{g,\rm EOS}(p_*,T_*).
\]

For an ideal gas, the last expression is \(c_{v,g}(T-T_*)\). For R410A,
both density and internal energy come from the fixed-composition Helmholtz
EOS. The oil heat capacity is `rho_cp/rho`. Required mixture inputs are
`dissolved_heat_capacity` and `desorption_internal_energy`; ideal gas also
requires `gas_heat_capacity_cv`. No air or propane value is silently assigned
to dissolved R410A. A positive desorption energy cools the mixture as gas
is released, and the reverse transfer heats it. This occurs through the
changing \(M_i u_i\) inventories. Adding a second \(\Gamma L\) source would
count the same heat twice.

The desorption input is an **internal-energy difference**, not latent heat
of the pure refrigerant. If measured data provide the partial desorption
enthalpy at the reference state, convert with

\[
\Delta u_{\rm des,*}=\Delta h_{\rm des,*}-p_*(v_{g,*}-v_{d,*}).
\]

Solubility and caloric fits must describe the same mixture. In the dilute,
ideal-gas van't Hoff limit the solubility temperature slope constrains the
desorption enthalpy. Arbitrary Henry slopes and caloric constants can pass
a conservation test without representing a thermodynamically calibrated oil.
The synthetic tests deliberately verify the numerical terms; they are not
FW68D thermochemical data. Nonideal R410A/oil equilibrium and partial enthalpy
data are still required for an experimental prediction.

## Dissipation, conduction and walls

The volume-weighted mixture viscosity is
\(\mu_m=\theta_l\mu_l+(1-\theta_l)\mu_g\). Couette dissipation is
\(\mu_m U^2/h\), in W/m². Pressure-flow dissipation on a face is
\(-q_{p,f}(p_R-p_L)\), which is nonnegative. Each internal face divides that
power equally between its neighboring cells; an axial half-cell pressure
drop belongs wholly to the boundary cell. This includes pressure heating at
the end faces and at resolved pressure gradients near a supply.

In-plane conduction uses \(h[\theta_l k_l+(1-\theta_l)k_g]\), harmonically
interpolated. Open ends prescribe advected inflow temperature and have zero
diffusive heat flux. Wall heat removal is
\(Q_w=H_j(T-T_j)+H_b(T-T_b)\). The wall coefficients and wall temperatures
must be specified or measured; zero coefficients mean adiabatic walls.

This is a common, depth-averaged phase temperature. The steady COUPLED
comparison also supports the separate finite-volume `resolved_cross_film`
energy path described below. It is an outer steady thermal iteration around
the conservative coupled pressure/species solve, not an additional unknown in
the transient Newton system.

## Steady resolved cross-film extension (2026-10-19)

For steady COUPLED calculations with `cross_film_thermal_model = FINITE_VOLUME`,
the solver alternates the coupled pressure/species equilibrium solve with a
finite-volume energy solve across the film layers. Layer temperatures feed
the local full `nu(c,T)` viscosity and phase properties on the next coupled
solve. The energy balance uses the same component mass fluxes, phase internal
energies and pressure-work terms as COUPLED conservation. It includes viscous
and pressure-flow dissipation, axial enthalpy transport, wall fluxes, and
groove/feed enthalpy. Film layers exchange heat with resolved journal and
bearing surfaces; the optional radial/axial solid network adds conduction
through the housing and its outer boundary. Layer and solid temperatures are
iterated with the fluid/geometry solution until the profile and global energy
residual gates pass.

This extension is steady only. It does not make the transient five-unknown
thermal Newton formulation cross-film resolved. The conditional LGE scenario
uses the same `cross_film_layers`, solid network, wall/feed temperatures,
cooling closures and common coupling relaxation (0.3) in both gas-off and
gas-on branches. Oil heat capacity (2010 J kg⁻¹ K⁻¹), liquid conductivity
(0.13 W m⁻¹ K⁻¹), steel conductivity (15 W m⁻¹ K⁻¹), and the conditional
R410A caloric values are proxies, not measured FW68D/R410A solution data.
In particular, the pure-R410A liquid internal-energy difference is an
analogue; no FW68D heat of solution or excess mixing energy has been measured.

The 02f driver rejects nonzero solver exit codes, nonconverged pressure or
coupling, exhausted iteration caps and restored-best states. It checks native
pressure residual below 1e-10; volume, pressure/species/thermal and profile
residuals against the requested coupling tolerance; oil/species balance at
the tighter of that tolerance and 1e-8; global energy residual against the
requested tolerance; independent load balance; finite positive admissible
fields; strict viscosity/chart support; and stable R410A EOS support at every
free-gas layer and wall. These are numerical/model-domain gates, not
experimental validation criteria.

For phase/component \(a\in\{o,d,g\}\) and film layer \(k\), the steady
discrete balance has the form

\[
\sum_{f,a}\dot m_{a,k,f}u_a(T_{\mathrm{donor},k,f},p_{\mathrm{donor},f})
+\sum_{\eta,a}\dot m_{a,k,\eta}u_a(T_{\mathrm{donor},\eta},p)
+W_k+Q_{\mathrm{cond},k}=D_k+E_{\mathrm{feed},k}.
\]

Face fluxes retain the accepted hydraulic pressure-driven, Couette and
diffusive contributions with their own donor states. Their layer weights
sum to the accepted depth-integrated phase flux. The internal transverse
fluxes enforce each layer's component balance and cancel when summed across
the gap. Phase reaction is distributed uniformly through the layers. Its
energy enters through the conserved component fluxes and internal energies,
so no separate latent-heat source is added.

With the common reference \((T_*,p_*)\),

\[
u_o=c_{p,o}(T-T_*),\qquad
u_d=c_{p,d}(T-T_*)-\Delta u,\qquad
u_g=u_{\mathrm{EOS}}(T,p)-u_{\mathrm{EOS}}(T_*,p_*).
\]

The liquid \(c_p\)-based internal-energy approximation follows the existing
incompressible caloric model. The gas energy is evaluated at actual local
pressure. The pressure-work contribution, summed over the layers, is
\(p\,\mathrm{div}\,Q+\mathrm{div}\,F_{\mathrm{diff},pV}-W_{\mathrm{feed}}\).
The accepted hydraulic dissipation is distributed with the resolved shear
profile while preserving its cell sum. Conduction links between film,
shaft and housing are symmetric, so interface exchanges cancel in the
combined energy ledger. A common relaxation factor is applied to film,
shaft and housing temperatures before the next hydraulic iteration.

Composition and phase fractions remain depth uniform, and phase equilibrium
uses film-mean temperature. This is a modeling limitation when cross-film
thermal gradients materially change local equilibrium composition.

## Vapor and refrigerant phase limits

Oil vapor is treated in the JFO limit of negligible vapor inventory. The
model does not transport oil vapor mass or its evaporation enthalpy. It is
appropriate only when those contributions are negligible relative to the
liquid and dissolved/free refrigerant energy. A temperature-dependent vapor
pressure alone does not supply a general thermal evaporation model.

R410A may dissolve in the oil or occupy the stable free-fluid branch. Above
the pure free-fluid dew line, a subcritical R410A bubble cannot remain a
stable gas. Equilibrium dissolution can remove that free mass when the
liquid solution has capacity. A finite-rate calculation that transports a
resolved bubble beyond that domain stops; the code does not mislabel
condensed refrigerant as gas. A separate refrigerant-rich liquid phase and
R32/R125 fractionation require a multicomponent flash model and mixture
property data. No universal absorption rate is imposed to conceal this limit.

Thermodynamic domain checks distinguish a supplied reservoir state from a
Newton trial. An outgoing trial state uses the same EOS continuation as
an interior trial. Treating a finite-difference perturbation of zero gas
mass as a real outgoing bubble previously invalidated the Jacobian for a
20°C, 20-bar gas-free oil solution. The corrected gate permits that
derivative evaluation; the accepted interior state still undergoes the
strict free-phase domain check. A dedicated regression verifies the linear
pressure profile of a dissolved, gas-free solution above the pure R410A
dew pressure, on one, two and four MPI ranks.

## Checks and outputs

`tests/test_coupled_energy.cpp` checks exact uniform Couette heating, implicit
wall cooling, an independently calculated adiabatic bubble compression,
closed desorption energy, a heated pure-vapor limit, and inlet/outlet energy
transport in both directions. It also checks grooves and circular holes,
outlet-temperature independence with variable viscosity, repeated nonlinear
calls at fixed old time, and real-fluid R410A release at 20–80°C.

`internal_energy` is J/m² and `energy_residual` is the cell-integrated
balance in W. `heat_generation` is W/m² and `gas_mass_transfer` is kg/(m³ s),
positive for release. CSV columns `energy_residual_W` and `internal_energy_J`
report global thermal quantities; they are blank for other energy paths.
The face mass/volume flux fields are also available for auditing. VTK
interpolates face fields to cell centers for display, so use the solver's
balance diagnostics when checking exact finite-volume conservation.
