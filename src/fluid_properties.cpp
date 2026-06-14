#include "fluid_properties.hpp"

#include <algorithm>
#include <cmath>

namespace FluidProperties {
namespace {

constexpr double kMinPositive = 1.0e-30;

double clamp_mass_fraction(double value, const SimulationConfig& cfg) {
    const double upper = std::max(0.0, cfg.dissolved_gas_max);
    return std::clamp(value, 0.0, upper);
}

double center_field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double bounded_temperature(double temperature, const SimulationConfig& cfg) {
    const double fallback = cfg.property_reference_temperature > 0.0
        ? cfg.property_reference_temperature
        : 300.0;
    return std::max(temperature > 0.0 ? temperature : fallback, 1.0);
}

bool property_model_enabled(const SimulationConfig& cfg) {
    return cfg.fluid_property_model != FluidPropertyModel::CONSTANT;
}

bool pressure_inlet_cell(const Fields& fields, int i, int j) {
    return fields.has("inlet_indicator") && fields["inlet_indicator"](i, j) > 0.5;
}

}  // namespace

double gas_specific_constant(DissolvedGasSpecies species) {
    switch (species) {
        case DissolvedGasSpecies::AIR: return 287.05;
        case DissolvedGasSpecies::PROPANE: return 188.55;
    }
    return 188.55;
}

double gas_density(double pressure, double temperature, const SimulationConfig& cfg) {
    const double p_abs = std::max(pressure, cfg.gas_pressure_floor);
    const double T = bounded_temperature(temperature, cfg);
    return p_abs / (gas_specific_constant(cfg.dissolved_gas_species) * T);
}

double interpolate_property_table(const std::vector<PropertyTablePoint>& table,
                                  double x,
                                  double fallback) {
    if (table.empty()) return fallback;
    if (x <= table.front().x) return table.front().value;
    if (x >= table.back().x) return table.back().value;

    for (size_t idx = 1; idx < table.size(); ++idx) {
        const auto& right = table[idx];
        if (x > right.x) continue;
        const auto& left = table[idx - 1];
        const double span = right.x - left.x;
        if (span <= 0.0) return right.value;
        const double t = (x - left.x) / span;
        return left.value + t * (right.value - left.value);
    }
    return fallback;
}

double equilibrium_dissolved_gas(double pressure, double temperature, const SimulationConfig& cfg) {
    if (!property_model_enabled(cfg)) return 0.0;

    const double p_abs = std::max(pressure, 0.0);
    double dissolved = 0.0;
    switch (cfg.oil_gas_solution_model) {
        case OilGasSolutionModel::HENRY: {
            // c_sat = H(T) p. van't Hoff temperature dependence: solubility falls with T.
            double henry = cfg.dissolved_gas_henry_coeff;
            if (cfg.dissolved_gas_henry_temp_coeff != 0.0) {
                const double T = bounded_temperature(temperature, cfg);
                const double T_ref = bounded_temperature(cfg.property_reference_temperature, cfg);
                henry *= std::exp(cfg.dissolved_gas_henry_temp_coeff * (1.0 / T - 1.0 / T_ref));
            }
            dissolved = henry * p_abs;
            break;
        }
        case OilGasSolutionModel::BUNSEN: {
            const double p_ref = std::max(cfg.property_reference_pressure, cfg.gas_pressure_floor);
            const double T_ref = bounded_temperature(cfg.property_reference_temperature, cfg);
            const double rho_g_ref = gas_density(p_ref, T_ref, cfg);
            dissolved = cfg.dissolved_gas_bunsen_coeff * (p_abs / p_ref) *
                (T_ref / bounded_temperature(temperature, cfg)) * rho_g_ref / std::max(cfg.rho, kMinPositive);
            break;
        }
        case OilGasSolutionModel::TABLE:
            dissolved = cfg.solubility_table_2d.empty()
                ? interpolate_property_table(cfg.solubility_table, p_abs, 0.0)
                : cfg.solubility_table_2d.interpolate(p_abs, bounded_temperature(temperature, cfg), 0.0);
            break;
    }
    return clamp_mass_fraction(dissolved, cfg);
}

double liquid_solution_density(double dissolved_gas,
                               double pressure,
                               double temperature,
                               const SimulationConfig& cfg) {
    if (!property_model_enabled(cfg) || cfg.density_model == DensityModel::PURE_OIL) {
        return std::max(cfg.rho, kMinPositive);
    }

    const double c_g = clamp_mass_fraction(dissolved_gas, cfg);
    double rho_solution = cfg.rho;
    switch (cfg.density_model) {
        case DensityModel::PURE_OIL:
            rho_solution = cfg.rho;
            break;
        case DensityModel::MASS_VOLUME_MIXING: {
            // Dissolved gas occupies ~its liquid-phase molar volume, not the gas phase.
            // Use the configured liquid-phase density when given; else fall back to ideal gas.
            const double dg_rho = cfg.dissolved_gas_liquid_density > 0.0
                ? cfg.dissolved_gas_liquid_density
                : std::max(gas_density(pressure, temperature, cfg), kMinPositive);
            const double specific_volume =
                (1.0 - c_g) / std::max(cfg.rho, kMinPositive) + c_g / std::max(dg_rho, kMinPositive);
            rho_solution = 1.0 / std::max(specific_volume, kMinPositive);
            if (cfg.solution_density_gas_coeff != 0.0) {
                rho_solution *= std::max(0.0, 1.0 + cfg.solution_density_gas_coeff * c_g);
            }
            break;
        }
        case DensityModel::TABLE:
            rho_solution = cfg.density_table_2d.empty()
                ? interpolate_property_table(cfg.density_table, c_g, cfg.rho)
                : cfg.density_table_2d.interpolate(pressure, bounded_temperature(temperature, cfg), cfg.rho);
            break;
    }
    return std::max(rho_solution, kMinPositive);
}

double liquid_solution_viscosity(double dissolved_gas,
                                 double pressure,
                                 double temperature,
                                 const SimulationConfig& cfg) {
    if (!property_model_enabled(cfg) || cfg.viscosity_model == ViscosityModel::PURE_OIL) {
        return std::max(cfg.mu, kMinPositive);
    }

    const double c_g = clamp_mass_fraction(dissolved_gas, cfg);
    double mu_solution = cfg.mu;
    switch (cfg.viscosity_model) {
        case ViscosityModel::PURE_OIL:
            mu_solution = cfg.mu;
            break;
        case ViscosityModel::LOG_MIXING:
            mu_solution = cfg.mu * std::exp(cfg.solution_viscosity_gas_coeff * c_g);
            break;
        case ViscosityModel::EMPIRICAL_CORRELATION: {
            // mu_l(T, p, c_d) = mu_oil * exp[E_mu (1/T - 1/T_ref)]   (Andrade temperature)
            //                          * exp[a_c c_d]                (dissolved-gas thinning)
            //                          * exp[alpha_p (p - p_ref)].   (Barus piezo-viscosity)
            const double p_ref = cfg.property_reference_pressure;
            const double T_ref = bounded_temperature(cfg.property_reference_temperature, cfg);
            const double T = bounded_temperature(temperature, cfg);
            const double mu_oil_T = cfg.mu * std::exp(cfg.viscosity_temperature_coeff * (1.0 / T - 1.0 / T_ref));
            mu_solution = mu_oil_T
                * std::exp(cfg.solution_viscosity_gas_coeff * c_g)
                * std::exp(cfg.viscosity_pressure_coeff * (pressure - p_ref));
            break;
        }
        case ViscosityModel::TABLE:
            mu_solution = interpolate_property_table(cfg.viscosity_table, c_g, cfg.mu);
            break;
    }
    return std::max(mu_solution, kMinPositive);
}

double gas_mixture_viscosity(double mu_liquid,
                             double alpha_gas,
                             double mass_quality,
                             const SimulationConfig& cfg) {
    const double mu_l = std::max(mu_liquid, kMinPositive);
    const double alpha_max = std::clamp(cfg.gas_alpha_max, kMinPositive, 1.0);
    const double alpha = std::clamp(alpha_gas, 0.0, alpha_max);
    const double quality = std::clamp(mass_quality, 0.0, 1.0);
    const double mu_g = std::max(cfg.mu_gas, kMinPositive);

    switch (cfg.gas_mixture_viscosity_model) {
        case GasMixtureViscosityModel::EINSTEIN_DILUTE:
            return mu_l * (1.0 + 2.5 * alpha);
        case GasMixtureViscosityModel::DUKLER_VOID:
            return alpha * mu_g + (1.0 - alpha) * mu_l;
        case GasMixtureViscosityModel::MCADAMS_QUALITY:
            return 1.0 / (quality / mu_g + (1.0 - quality) / mu_l);
        case GasMixtureViscosityModel::KRIEGER_DOUGHERTY: {
            // Cap the packing argument to keep the divergence finite at alpha_max.
            const double packing = std::min(alpha / alpha_max, 0.99);
            return mu_l * std::pow(1.0 - packing, -2.5 * alpha_max);
        }
        case GasMixtureViscosityModel::LINEAR_QUALITY:
            return quality * mu_g + (1.0 - quality) * mu_l;
    }
    return mu_l;
}

void update_solution_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    const bool has_pressure = fields.has("pressure");
    const bool has_temperature = fields.has("temperature");
    const bool has_theta = fields.has("theta");

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double pressure = has_pressure ? fields["pressure"](i, j) : cfg.initial_pressure();
            const double temperature = has_temperature ? fields["temperature"](i, j) : cfg.temperature_initial;
            const double theta = has_theta ? fields["theta"](i, j) : 1.0;

            if (fields.has("dissolved_gas") && cfg.fluid_property_model == FluidPropertyModel::CONSTANT) {
                fields["dissolved_gas"](i, j) = 0.0;
            } else if (fields.has("dissolved_gas")) {
                fields["dissolved_gas"](i, j) = clamp_mass_fraction(fields["dissolved_gas"](i, j), cfg);
            }
            const double dissolved = center_field_value(
                fields, "dissolved_gas",
                property_model_enabled(cfg) ? cfg.dissolved_gas_initial : 0.0,
                i, j);

            const double rho_liquid = liquid_solution_density(dissolved, pressure, temperature, cfg);
            double mu_liquid = liquid_solution_viscosity(dissolved, pressure, temperature, cfg);
            // 2-D table viscosity is the KINEMATIC nu(p,T) [m^2/s]; convert per cell to
            // dynamic mu = nu * rho_solution (the supplier chart carries no density column).
            if (cfg.viscosity_model == ViscosityModel::TABLE && !cfg.viscosity_table_2d.empty()) {
                const double nu = cfg.viscosity_table_2d.interpolate(
                    pressure, bounded_temperature(temperature, cfg),
                    cfg.mu / std::max(rho_liquid, kMinPositive));
                mu_liquid = std::max(nu * rho_liquid, kMinPositive);
            }
            const double cp_liquid = cfg.rho_cp / std::max(rho_liquid, kMinPositive);
            const double rho_g = gas_density(pressure, temperature, cfg);

            // Effective film viscosity. The dissolved gas already thinned mu_liquid;
            // free gas modifies the mixture per the configured model (Einstein dilute
            // thickening, void/quality-weighted thinning toward mu_gas, or
            // Krieger-Dougherty packing).
            double mu_mixture = mu_liquid;
            if (cfg.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE &&
                fields.has("alpha_gas")) {
                const double alpha_g = fields["alpha_gas"](i, j);
                double quality = 0.0;
                if (fields.has("free_gas_mass") && fields.has("h")) {
                    const double m_g = std::max(fields["free_gas_mass"](i, j), 0.0);
                    const double h_val = std::max(fields["h"](i, j), cfg.min_film_thickness);
                    const double liquid_mass = rho_liquid * std::max(theta, 0.0) * h_val;
                    quality = m_g / std::max(m_g + liquid_mass, kMinPositive);
                }
                mu_mixture = gas_mixture_viscosity(mu_liquid, alpha_g, quality, cfg);
            }

            if (fields.has("rho_liquid_solution")) fields["rho_liquid_solution"](i, j) = rho_liquid;
            if (fields.has("mu_liquid_solution")) fields["mu_liquid_solution"](i, j) = mu_liquid;
            if (fields.has("cp_liquid_solution")) fields["cp_liquid_solution"](i, j) = cp_liquid;
            if (fields.has("rho_cp_liquid_solution")) fields["rho_cp_liquid_solution"](i, j) = rho_liquid * cp_liquid;
            if (fields.has("k_liquid_solution")) fields["k_liquid_solution"](i, j) = cfg.thermal_conductivity;
            if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
            if (fields.has("mu")) fields["mu"](i, j) = mu_mixture;
            if (fields.has("rho")) fields["rho"](i, j) = rho_liquid * theta;
        }
    }
}

void update_gas_state(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg, double dt) {
    const bool active = cfg.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE &&
        dt > 0.0 && cfg.gas_mass_transfer_rate > 0.0;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double pressure = center_field_value(fields, "pressure", cfg.initial_pressure(), i, j);
            const double temperature = center_field_value(fields, "temperature", cfg.temperature_initial, i, j);
            const double h = std::max(center_field_value(fields, "h", cfg.c, i, j), cfg.min_film_thickness);
            const double theta = center_field_value(fields, "theta", 1.0, i, j);
            const double rho_liquid = center_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
            const double rho_g = gas_density(pressure, temperature, cfg);
            const double alpha_max = std::clamp(cfg.gas_alpha_max, 0.0, 1.0);
            const double release_mass_cap = alpha_max * rho_g * h;

            double dissolved = center_field_value(fields, "dissolved_gas", 0.0, i, j);
            double free_mass = center_field_value(fields, "free_gas_mass", 0.0, i, j);
            double source = 0.0;

            if (pressure_inlet_cell(fields, i, j)) {
                dissolved = property_model_enabled(cfg) ? cfg.dissolved_gas_initial : 0.0;
                free_mass = 0.0;
            } else if (active) {
                dissolved = clamp_mass_fraction(dissolved, cfg);
                free_mass = std::max(free_mass, 0.0);
                const double equilibrium = equilibrium_dissolved_gas(pressure, temperature, cfg);
                const double liquid_mass_area =
                    std::max(rho_liquid * h * std::max(theta, cfg.theta_min), kMinPositive);

                if (dissolved > equilibrium && release_mass_cap > free_mass) {
                    const double fraction_release =
                        std::min((dissolved - equilibrium) * cfg.gas_mass_transfer_rate * dt,
                                 dissolved);
                    const double mass_release = std::min(fraction_release * liquid_mass_area,
                                                         release_mass_cap - free_mass);
                    dissolved -= mass_release / liquid_mass_area;
                    free_mass += mass_release;
                    source = mass_release / (h * dt);
                } else if (dissolved < equilibrium && free_mass > 0.0) {
                    const double available_capacity =
                        std::max(0.0, (std::min(equilibrium, cfg.dissolved_gas_max) - dissolved) *
                                      liquid_mass_area);
                    const double kinetic_limit = available_capacity * cfg.gas_mass_transfer_rate * dt;
                    const double mass_resorb = std::min({kinetic_limit, available_capacity, free_mass});
                    dissolved += mass_resorb / liquid_mass_area;
                    free_mass -= mass_resorb;
                    source = -mass_resorb / (h * dt);
                }
            } else if (cfg.fluid_property_model != FluidPropertyModel::GAS_CAVITATION_MIXTURE) {
                free_mass = 0.0;
            }

            free_mass = std::max(free_mass, 0.0);
            const double alpha = (rho_g > 0.0 && h > 0.0)
                ? std::clamp(free_mass / (rho_g * h), 0.0, alpha_max)
                : 0.0;

            if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = clamp_mass_fraction(dissolved, cfg);
            if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = free_mass;
            if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = alpha;
            if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = source;
            if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
        }
    }
}

void equilibrate_gas(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.fluid_property_model != FluidPropertyModel::GAS_CAVITATION_MIXTURE) return;
    const double alpha_max = std::clamp(cfg.gas_alpha_max, 0.0, 1.0);
    const double c_feed = clamp_mass_fraction(cfg.dissolved_gas_initial, cfg);

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double pressure = center_field_value(fields, "pressure", cfg.initial_pressure(), i, j);
            const double temperature = center_field_value(fields, "temperature", cfg.temperature_initial, i, j);
            const double h = std::max(center_field_value(fields, "h", cfg.c, i, j), cfg.min_film_thickness);
            const double theta = center_field_value(fields, "theta", 1.0, i, j);
            const double rho_liquid = center_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
            const double rho_g = gas_density(pressure, temperature, cfg);

            if (pressure_inlet_cell(fields, i, j)) {
                if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = c_feed;
                if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = 0.0;
                if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = 0.0;
                if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = 0.0;
                if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
                continue;
            }

            // Local flash from the feed composition: gas above the saturation limit comes out.
            const double c_eq = equilibrium_dissolved_gas(pressure, temperature, cfg);
            const double dissolved = std::min(c_eq, c_feed);
            const double released = std::max(0.0, c_feed - dissolved);
            const double liquid_mass_area = rho_liquid * std::max(theta, 0.0) * h;
            double free_mass = std::min(released * liquid_mass_area, alpha_max * rho_g * h);
            free_mass = std::max(free_mass, 0.0);
            const double alpha = (rho_g > 0.0 && h > 0.0)
                ? std::clamp(free_mass / (rho_g * h), 0.0, alpha_max) : 0.0;

            if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = dissolved;
            if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = free_mass;
            if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = alpha;
            if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = 0.0;
            if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
        }
    }
}

}  // namespace FluidProperties
