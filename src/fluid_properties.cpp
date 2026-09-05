#include "fluid_properties.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "full_film.hpp"
#include "utils.hpp"

namespace FluidProperties {
namespace {

constexpr double kMinPositive = 1.0e-30;
constexpr double kKelvinOffset = 273.15;
/// Lower bound on the POLYNOMIAL viscosity, as a fraction of the reference mu.
/// The polynomial is a fit, not a physical law: extrapolated far enough it
/// crosses zero and turns negative, which would make the Reynolds operator
/// indefinite. 1e-3 mu is far below any oil yet still strictly positive.
constexpr double kPolynomialViscosityFloor = 1.0e-3;

/// One-shot notice that the POLYNOMIAL law has been clamped. Per rank-0 process,
/// not per cell: a clamped run touches thousands of cells per timestep and the
/// message is about the model being extrapolated, not about any one cell.
void warn_polynomial_viscosity_floor(double T_C, double mu_bar) {
    static bool warned = false;
    if (warned) return;
    warned = true;
    Utils::log("WARNING: POLYNOMIAL viscosity fell to mubar = " + std::to_string(mu_bar) +
                   " at T = " + std::to_string(T_C) + " degC and was floored at " +
                   std::to_string(kPolynomialViscosityFloor) +
                   " mu. The polynomial is being extrapolated outside its fitted "
                   "temperature range; check viscosity_poly_k0/k1/k2 and the film "
                   "temperatures. Reported once per run.",
               Utils::Color::YELLOW);
}

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

/// Numerical relaxation for the steady feed-flash iteration only. Transient
/// phase exchange is physical time integration and must not use this factor.
double gas_state_relaxation(const SimulationConfig& cfg) {
    const bool gas_is_volumetric =
        cfg.gas_pressure_coupling == GasPressureCoupling::VOID_COUPLED ||
        FullFilm::gas_displaces_liquid(cfg);
    if (!gas_is_volumetric) return 1.0;
    return std::clamp(cfg.void_coupling_relaxation, 0.01, 1.0);
}

double equilibrium_dissolved_gas(double pressure, double temperature, const SimulationConfig& cfg) {
    if (!property_model_enabled(cfg)) return 0.0;

    const double p_abs = std::max(pressure, 0.0);
    double dissolved = 0.0;
    switch (cfg.solubility_model) {
        case SolubilityModel::HENRY: {
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
        case SolubilityModel::TABLE:
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
    if (!property_model_enabled(cfg) || cfg.density_model == DensityModel::CONSTANT) {
        return std::max(cfg.rho, kMinPositive);
    }

    const double c_g = clamp_mass_fraction(dissolved_gas, cfg);
    double rho_solution = cfg.rho;
    switch (cfg.density_model) {
        case DensityModel::CONSTANT:
            rho_solution = cfg.rho;
            break;
        case DensityModel::MIXTURE: {
            if (c_g == 0.0) break;
            // The dissolved species occupies its partial liquid volume. The free
            // gas EOS cannot supply this property.
            if (!(cfg.dissolved_gas_liquid_density > 0.0))
                throw std::invalid_argument("MIXTURE density requires dissolved_gas_liquid_density > 0.");
            const double dg_rho = cfg.dissolved_gas_liquid_density;
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
    if (!property_model_enabled(cfg) || cfg.liquid_viscosity_model == LiquidViscosityModel::CONSTANT) {
        return std::max(cfg.mu, kMinPositive);
    }

    const double c_g = clamp_mass_fraction(dissolved_gas, cfg);
    double mu_solution = cfg.mu;
    switch (cfg.liquid_viscosity_model) {
        case LiquidViscosityModel::CONSTANT:
            mu_solution = cfg.mu;
            break;
        case LiquidViscosityModel::EMPIRICAL: {
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
        case LiquidViscosityModel::POLYNOMIAL: {
            // Ferron, Frene & Boncompain (1983), nondimensional viscosity:
            //   mubar = mu/mu_i = k0 - k1 Tbar + k2 Tbar^2,   Tbar = T/T_i,
            // where the paper's T and T_i are in DEGREES CELSIUS (its nomenclature
            // reads "T = temperature degC" and "Tbar = T/T_0"). Their VG32 set
            // (3.287, 3.064, 0.777, mu_i = 0.0277 Pa s at T_i = 40 degC) returns
            // exactly 1 at 40 degC and 0.267 at 80 degC; the same coefficients fed a
            // Kelvin ratio return 0.82 at 80 degC, which is not the measured oil.
            // Solver temperatures are Kelvin, so the conversion is explicit here.
            // This is a pure temperature law: no Barus and no dissolved-gas factor.
            const double T_ref_C = cfg.polynomial_reference_temperature_celsius();
            if (!(T_ref_C > 0.0)) {
                mu_solution = cfg.mu;  // validate() rejects this; stay finite anyway
                break;
            }
            const double T_C = bounded_temperature(temperature, cfg) - kKelvinOffset;
            const double T_bar = T_C / T_ref_C;
            const double mu_bar = cfg.viscosity_poly_k0
                - cfg.viscosity_poly_k1 * T_bar
                + cfg.viscosity_poly_k2 * T_bar * T_bar;
            mu_solution = cfg.mu * mu_bar;
            // A quadratic with a real root goes negative past it, and a linear-only
            // set always does. Floor the film at a small positive fraction of mu so
            // the Reynolds mobility stays positive-definite, and say so once.
            const double mu_floor = kPolynomialViscosityFloor * std::max(cfg.mu, kMinPositive);
            if (mu_solution < mu_floor) {
                warn_polynomial_viscosity_floor(T_C, mu_bar);
                mu_solution = mu_floor;
            }
            break;
        }
        case LiquidViscosityModel::TABLE:
            if (!cfg.viscosity_table_2d.empty()) {
                const double rho_l = liquid_solution_density(c_g, pressure, temperature, cfg);
                mu_solution = rho_l * cfg.viscosity_table_2d.interpolate(
                    pressure, bounded_temperature(temperature, cfg), cfg.mu / rho_l);
            } else {
                mu_solution = interpolate_property_table(cfg.viscosity_table, c_g, cfg.mu);
            }
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

    switch (cfg.mixture_viscosity_model) {
        case MixtureViscosityModel::DUKLER:  // legacy: void-fraction-weighted linear blend
            return alpha * mu_g + (1.0 - alpha) * mu_l;
        case MixtureViscosityModel::LINEAR:  // mass-quality-weighted linear (Grando 2006)
            return quality * mu_g + (1.0 - quality) * mu_l;
        case MixtureViscosityModel::MCADAMS:  // homogeneous reciprocal blend (default)
            return 1.0 / (quality / mu_g + (1.0 - quality) / mu_l);
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
            const double cp_liquid = cfg.rho_cp / std::max(rho_liquid, kMinPositive);
            const double rho_g = gas_density(pressure, temperature, cfg);

            // Dissolved gas modifies the liquid viscosity; free gas enters through
            // the selected void-fraction or mass-quality mixture law.
            double mu_mixture = mu_liquid;
            if (cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE &&
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
    const bool active = cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE &&
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
                const double liquid_mass_area = std::max(rho_liquid * h * theta, kMinPositive);
                // dM_l/dt = dM_d/dt = -k (M_d - c_eq M_l), with oil M_l-M_d
                // invariant. Integrate exactly at frozen p,T; expm1 retains small steps.
                const double rate_dt = cfg.gas_mass_transfer_rate * dt;
                const double exponent = (1.0 - equilibrium) * rate_dt;
                const double kinetic_time = exponent > 1e-8
                    ? -std::expm1(-exponent) / (1.0 - equilibrium)
                    : rate_dt * (1.0 - 0.5 * exponent);
                double transfer = liquid_mass_area * (dissolved - equilibrium) * kinetic_time;
                transfer = std::clamp(transfer, -free_mass,
                    std::min(dissolved * liquid_mass_area, std::max(release_mass_cap - free_mass, 0.0)));
                const double remaining = liquid_mass_area - transfer;
                if (remaining > kMinPositive) {
                    dissolved = (dissolved * liquid_mass_area - transfer) / remaining;
                    free_mass += transfer;
                    const double rho_new = liquid_solution_density(dissolved, pressure, temperature, cfg);
                    if (fields.has("theta")) fields["theta"](i, j) = remaining / (rho_new * h);
                    if (fields.has("rho_liquid_solution")) fields["rho_liquid_solution"](i, j) = rho_new;
                    if (fields.has("rho")) fields["rho"](i, j) = remaining / h;
                    source = transfer / (h * dt);
                }
            } else if (cfg.fluid_property_model != FluidPropertyModel::TWO_PHASE) {
                free_mass = 0.0;
            }

            free_mass = std::max(free_mass, 0.0);
            const double alpha = (rho_g > 0.0 && h > 0.0)
                ? std::clamp(free_mass / (rho_g * h), 0.0, alpha_max)
                : 0.0;

            if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = clamp_mass_fraction(dissolved, cfg);
            // Physical timesteps have no iterative relaxation: applying it to only
            // one phase destroys mass and invents a pressure-dependent delay.
            if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = free_mass;
            if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = alpha;
            if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = source;
            if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
        }
    }
}

void equilibrate_gas(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.fluid_property_model != FluidPropertyModel::TWO_PHASE) return;
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
            const double alpha_max = std::clamp(cfg.gas_alpha_max, 0.0, 1.0);
            const double c_eq = equilibrium_dissolved_gas(pressure, temperature, cfg);
            double dissolved = std::min(c_eq, c_feed);
            const double liquid_mass_area = rho_liquid * std::max(theta, 0.0) * h;
            // Feed mass fraction is defined relative to solution + free gas.
            // c_feed = (c_d M_l + M_g)/(M_l + M_g). A volume cap must leave
            // unreleased gas dissolved; it cannot discard the excess inventory.
            double free_mass = std::min((c_feed - dissolved) * liquid_mass_area /
                std::max(1.0 - c_feed, kMinPositive), alpha_max * rho_g * h);
            free_mass = std::max(free_mass, 0.0);

            const double w = gas_state_relaxation(cfg);
            free_mass = (1.0 - w) * center_field_value(fields, "free_gas_mass", 0.0, i, j) + w * free_mass;
            free_mass = std::clamp(free_mass, 0.0, std::min(alpha_max * rho_g * h,
                c_feed * liquid_mass_area / std::max(1.0 - c_feed, kMinPositive)));
            dissolved = c_feed - (1.0 - c_feed) * free_mass / std::max(liquid_mass_area, kMinPositive);
            const double alpha = std::clamp(free_mass / (rho_g * h), 0.0, alpha_max);
            if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = dissolved;
            if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = free_mass;
            if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = alpha;
            if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = 0.0;
            if (fields.has("rho_gas")) fields["rho_gas"](i, j) = rho_g;
        }
    }
}

}  // namespace FluidProperties
