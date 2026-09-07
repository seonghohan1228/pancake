#pragma once

/// Fixed-composition R410A pseudo-pure EOS (Lemmon, 2003). SI throughout.
/// This describes the free refrigerant, not its solubility in a particular oil.
namespace R410A {
constexpr double gas_constant = 8.314472 / 0.0725854;
constexpr double critical_temperature = 344.494;
constexpr double critical_density = 459.0300696;
constexpr double maximum_pressure = 50e6; // EOS range; bearing pressures may exceed reservoir pressures.

/// Dew pressure below the critical temperature; infinity above it.
double dew_pressure(double temperature);
/// Helmholtz pressure and isothermal derivative on a specified density branch.
double pressure(double density, double temperature);
double pressure_derivative(double density, double temperature);
/// Specific internal energy and isochoric heat capacity [J/kg, J/(kg K)].
/// Same Helmholtz reference convention as CoolProp; density branch supplied by caller.
double internal_energy(double density, double temperature);
double heat_capacity_cv(double density, double temperature);
/// Stable gas density. Rejects condensed/two-phase states and extrapolation.
/// Supported envelope: 200--500 K, up to 50 MPa, below dew pressure when subcritical.
double gas_density(double pressure, double temperature);
}
