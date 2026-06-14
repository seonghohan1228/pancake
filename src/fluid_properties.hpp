#pragma once

#include <vector>

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

namespace FluidProperties {

/// Specific gas constant for the configured dissolved/free gas species [J/(kg K)].
double gas_specific_constant(DissolvedGasSpecies species);

/// Ideal-gas density with pressure and temperature floors for numerical safety.
double gas_density(double pressure, double temperature, const SimulationConfig& cfg);

/// Equilibrium dissolved-gas mass fraction for the configured solution model.
double equilibrium_dissolved_gas(double pressure, double temperature, const SimulationConfig& cfg);

/// Liquid oil plus dissolved-gas solution density [kg/m^3].
double liquid_solution_density(double dissolved_gas,
                               double pressure,
                               double temperature,
                               const SimulationConfig& cfg);

/// Liquid oil plus dissolved-gas solution dynamic viscosity [Pa s].
double liquid_solution_viscosity(double dissolved_gas,
                                 double pressure,
                                 double temperature,
                                 const SimulationConfig& cfg);

/// Effective film viscosity with free gas present, per the configured
/// gas_mixture_viscosity_model. mass_quality is m_g / (m_g + rho_l theta h)
/// (used by the quality-based models only). All models return mu_liquid at
/// alpha_gas = 0.
double gas_mixture_viscosity(double mu_liquid,
                             double alpha_gas,
                             double mass_quality,
                             const SimulationConfig& cfg);

/// Piecewise-linear interpolation for empirical property tables.
double interpolate_property_table(const std::vector<PropertyTablePoint>& table,
                                  double x,
                                  double fallback);

/// Refresh rho/mu/thermal property fields from pressure, temperature, theta, and dissolved gas.
void update_solution_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

/// Finite-rate release/resorption between dissolved gas and free gas in JFO void space.
void update_gas_state(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg, double dt);

/// Steady-state equilibrium gas flash: set the dissolved fraction to the local
/// saturation value and release the excess (relative to the feed composition) as
/// free gas. Used in STEADY_STATE with steady_gas_model = EQUILIBRIUM, where time
/// advance (and therefore the finite-rate transport) is not meaningful.
void equilibrate_gas(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace FluidProperties
