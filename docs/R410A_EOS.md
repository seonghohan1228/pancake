# R410A free-phase equation of state

`src/r410a.cpp` independently evaluates the 21-term residual Helmholtz
equation for fixed-composition, pseudo-pure R410A
(Lemmon, *[Int. J. Thermophys.](https://doi.org/10.1023/A:1025048800563)*, 2003).
Coefficients were checked against the JSON delivered by CoolProp 7.2.0;
the [CoolProp fluid documentation](https://coolprop.org/fluid_properties/fluids/R410A.html)
identifies that EOS and its reference constants.

\[
\alpha^r(\delta,\tau)=\sum_i n_i\delta^{d_i}\tau^{t_i}
\exp(-\delta^{l_i}),\qquad
p=\rho R_sT(1+\delta\alpha^r_\delta).
\]

Terms with \(l_i=0\) omit the exponential. The reducing temperature is
344.494 K, density 459.0300696 kg/m³, and specific gas constant
8.314472 / 0.0725854 J/(kg K). An analytic isothermal pressure derivative
supports safeguarded Newton inversion on the dilute, stable gas branch.
Temperature terms are cached per thread and integer density powers are
reused, avoiding repeated general-power evaluations during field solves.

The current free-fluid API accepts 200–500 K and pressures up to 50 MPa.
The pressure ceiling must exceed the reservoir pressure: hydrodynamic film
pressure can exceed 20 bar even with reservoirs at or below 20 bar.
Below the critical
temperature it rejects pressures above the R410A dew curve, rather than
silently returning a liquid density labeled as gas. At 20 °C and 20 bar,
pure R410A is outside this gas domain. A gas-free oil solution can still
be pressurized to that condition; a cell retaining a separate refrigerant
phase there needs an additional condensation/flash treatment.
Above the critical temperature, the stable free-fluid branch may be dense;
calling it an ideal gas would be inappropriate.

Thermal calculations also evaluate the ideal and residual Helmholtz terms:

\[
u=R_sT\tau(\alpha^0_\tau+\alpha^r_\tau),\qquad
c_v=-R_s\tau^2(\alpha^0_{\tau\tau}+\alpha^r_{\tau\tau}).
\]

The mixture energy solver subtracts a gas reference internal energy before
combining the free and dissolved inventories. That reference must be paired
with the specified desorption internal-energy difference; see
[the energy formulation](COUPLED_ENERGY.md). The real-fluid EOS provides the
pressure dependence of gas internal energy, which a constant ideal-gas heat
capacity cannot represent.

The coupled Newton residual has a continuation outside the gas domain so
it can visit intermediate guesses and gas-free pressurized cells. This is
not accepted as a physical free-phase state: cells containing a resolved
free refrigerant inventory must pass the strict EOS domain check.

The implementation does not resolve R32/R125 fractionation. It is a
fixed-blend approximation and must not be described as a full binary
refrigerant/oil equilibrium calculation. Nor does a pure-refrigerant EOS
provide oil solubility, dissolved partial volume, mixture viscosity, or
enthalpy of dissolution.

`test_r410a` compares density against 15 independently generated CoolProp
states at 20, 40, 60 and 80 °C, including film pressures above 20 bar and
100 bar at 80 °C. Four caloric reference states check internal energy and
heat capacity. Derivative checks and rejected-domain cases supplement them.
`test_coupled_cavitation` checks real-fluid volume closure and species
conservation during finite-rate release at all four temperatures.
`test_coupled_energy` uses those thermodynamic functions in conservative
thermal release at the same four temperatures.

This is implementation verification against another evaluation of the same
published EOS, not independent experimental validation of a refrigerant/oil
mixture. In particular, no FW68D–R410A solubility curve has been inferred
from propane, FVS50, FVC68D, or ISO-68 POE data.
