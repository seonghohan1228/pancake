#include "fluid_properties.hpp"
#include "property_table_io.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, const char* message) {
    check(std::abs(actual - expected) <= 1e-11 * std::max(std::abs(expected), 1e-6), message);
}
template<class Function> void rejects(Function function, const char* message) {
    bool rejected = false;
    try { function(); } catch (const std::domain_error&) { rejected = true; }
    check(rejected, message);
}
}

int main() {
    const auto source = std::filesystem::path(PANCAKE_SOURCE_DIR);
    SimulationConfig cfg;
    cfg.load_from_file(source / "config_r410a_fw68d.txt");
    std::vector<std::string> errors, warnings;
    FluidProperties::load_config_property_tables(cfg, errors);
    cfg.validate(errors, warnings);
    for (const auto& error : errors) std::cerr << error << '\n';
    check(errors.empty(), "FW68D example must load and validate");
    const auto& table = cfg.kinematic_composition_table;
    check(table.axis.size() == 11 && cfg.solubility_curve_table.axis.size() == 11,
          "11 composition curves and 11 absolute-pressure curves");
    std::size_t count = 0;
    cfg.solution_table_clamp = false;
    for (std::size_t i = 0; i < table.axis.size(); ++i) {
        for (auto point : table.curves[i]) {
            ++count;
            near(FluidProperties::liquid_solution_viscosity(table.axis[i], 1e6, point.x, cfg),
                 940 * point.value, "Every original nu node, including c=0 and c=1, must be reproduced");
        }
    }
    check(count == 125, "Preserve all 125 chart points");
    for (std::size_t i = 0; i < cfg.solubility_curve_table.axis.size(); ++i) {
        for (auto point : cfg.solubility_curve_table.curves[i]) {
            near(FluidProperties::equilibrium_dissolved_gas(cfg.solubility_curve_table.axis[i], point.x, cfg),
                 point.value, "Solubility source nodes must reproduce mass fraction");
        }
    }
    // Calibration-derived source checks catch MPa/Pa, Celsius/K, cSt/SI and %/fraction mistakes.
    near(table.curves.front()[6].x, 313.16773902022654, "Pure-oil calibrated temperature");
    near(table.curves.front()[6].value, 74.42707818477315e-6, "Pure-oil calibrated kinematic viscosity");
    near(cfg.solubility_curve_table.axis.front(), 100000, "Pressure is absolute, with no atmospheric offset");

    // Ragged temperature grids must be interpolated separately before concentration blending.
    PropertyCurveTable ragged{{0, 1}, {{{300, 1}, {400, 9}}, {{320, 4}, {380, 16}}}};
    near(ragged.interpolate(.5, 350, false, true), std::sqrt(24.), "Log interpolation in T and c");
    near(ragged.interpolate(.5, 350, false, false), 7.5, "Linear interpolation in T and c");
    near(ragged.interpolate(0, 300, false, true), 1, "Unused neighbour must not cause range failure");
    rejects([&] { ragged.interpolate(.5, 310, false); }, "Reject outside either contributing curve");
    near(ragged.interpolate(.5, 310, true), 2.9, "Explicit endpoint clamp");
    rejects([&] { ragged.interpolate(2, 350, false); }, "Reject concentration range");
    rejects([&] { ragged.interpolate(0, std::numeric_limits<double>::quiet_NaN(), true); }, "Reject NaN");

    cfg.solution_table_clamp = true;
    const double t = 313.15, p = 1e6;
    const double concentration = FluidProperties::equilibrium_dissolved_gas(p, t, cfg);
    near(cfg.saturation_pressure(concentration, t), p, "Bubble-pressure inverse uses the same isotherm");
    cfg.solution_table_clamp = false;
    rejects([&] { FluidProperties::equilibrium_dissolved_gas(1, t, cfg); }, "Reject pressure below data");
    rejects([&] { cfg.saturation_pressure(.01, t); }, "Reject unsupported bubble-pressure inversion");
    cfg.solution_table_clamp = true;
    const double viscosity_constant = FluidProperties::liquid_solution_viscosity(.1, p, t, cfg);
    cfg.density_model = DensityModel::LINEAR_TEMPERATURE;
    near(FluidProperties::liquid_solution_density(.1, p, t, cfg), 925.764, "Celsius linear density at 40 C");
    near(FluidProperties::liquid_solution_viscosity(.1, p, t, cfg),
         viscosity_constant * 925.764 / 940, "mu conversion follows selected density");
    cfg.density_temperature_intercept = 1000;
    cfg.density_temperature_slope = -1;
    near(FluidProperties::liquid_solution_density(.9, p, t, cfg), 960, "Density coefficients are configurable");
    rejects([&] { FluidProperties::liquid_solution_density(0, p, 2000, cfg); }, "Reject nonpositive density");
    cfg.kinematic_viscosity_log_interpolation = false;
    cfg.solution_table_clamp = false;
    SimulationConfig copy;
    copy.load_from_text(cfg.to_config_text());
    check(copy.parse_errors.empty() && copy.solubility_model == SolubilityModel::CURVE_TABLE &&
          copy.liquid_viscosity_model == LiquidViscosityModel::KINEMATIC_TABLE &&
          copy.density_model == DensityModel::LINEAR_TEMPERATURE &&
          copy.solution_table_file == cfg.solution_table_file && !copy.solution_table_clamp &&
          !copy.kinematic_viscosity_log_interpolation, "Round-trip every model and option");
    near(copy.density_temperature_intercept, 1000, "Round-trip intercept");
    near(copy.density_temperature_slope, -1, "Round-trip slope");

    const auto bad = std::filesystem::temp_directory_path() / "pancake_fw68d_invalid.csv";
    // c=0 is no-species bookkeeping for a fixed-property surrogate, not its measured composition.
    {
        std::ofstream out(bad);
        out << "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n0,30,,8\n0,40,,2\n";
    }
    SimulationConfig fixed = cfg;
    fixed.solution_table_file = bad;
    fixed.cavitation_model = CavitationModel::FBNS;
    fixed.dissolved_gas_species = DissolvedGasSpecies::AIR;
    fixed.dissolved_gas_initial = 0;
    fixed.dissolved_gas_henry_coeff = 0;
    fixed.solubility_model = SolubilityModel::HENRY;
    fixed.density_model = DensityModel::CONSTANT;
    fixed.rho = 950;
    fixed.kinematic_viscosity_log_interpolation = true;
    errors.clear();
    FluidProperties::load_config_property_tables(fixed, errors);
    fixed.validate(errors, warnings);
    check(errors.empty(), "Viscosity-only fixed-property configuration must load and validate");
    check(fixed.solubility_curve_table.empty() && fixed.kinematic_composition_table.axis.size() == 1,
          "Viscosity-only input must not invent a solubility curve");
    near(FluidProperties::liquid_solution_viscosity(0, p, 303.15, fixed), 950 * 8e-6,
         "Singleton table preserves the first Celsius/SI node");
    near(FluidProperties::liquid_solution_viscosity(0, p, 308.15, fixed), 950 * 4e-6,
         "Singleton table interpolates log nu linearly in Celsius");
    near(FluidProperties::liquid_solution_viscosity(0, p, 313.15, fixed), 950 * 2e-6,
         "Singleton table preserves the last Celsius/SI node");
    rejects([&] { FluidProperties::liquid_solution_viscosity(0, p, 303.14, fixed); },
            "Singleton table must reject temperature below its range");
    rejects([&] { FluidProperties::liquid_solution_viscosity(0, p, 313.16, fixed); },
            "Singleton table must reject temperature above its range");
    rejects([&] { fixed.kinematic_composition_table.interpolate(.1, 308.15, false, true); },
            "Singleton table must reject a tracked composition outside its zero node");
    auto thermal = fixed;
    thermal.solution_mode = SolutionMode::STEADY_STATE;
    thermal.temperature_model = TemperatureModel::ENERGY_EQUATION;
    thermal.cross_film_thermal_model = CrossFilmModel::FINITE_VOLUME;
    thermal.groove_mixing_model = GrooveMixingModel::FLUX_BALANCED;
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90;
    groove.size = 15;
    groove.p_supply = 1500000;
    groove.feeds_fresh_oil = true;
    groove.t_supply = 303.15;
    thermal.inlets = {groove};
    thermal.thermal_convection_scheme = ConvectionScheme::UPWIND;
    thermal.liquid_specific_heat = 2010;
    thermal.bush_thermal_model = BushThermalModel::CONDUCTION_3D;
    thermal.bush_outer_radius = thermal.R + .02;
    thermal.bush_axial_cells = thermal.n_z_global;
    thermal.bush_end_natural_convection = true;
    thermal.density_temperature_slope = -.7284;
    errors.clear();
    warnings.clear();
    thermal.validate(errors, warnings);
    for (const auto& error : errors) std::cerr << error << '\n';
    check(errors.empty(), "Constant-density conjugate housing accepts explicit Boussinesq buoyancy");
    check(std::any_of(warnings.begin(), warnings.end(), [](const auto& warning) {
        return warning.find("only for exterior housing-end buoyancy") != std::string::npos;
    }), "Constant-density buoyancy assumption must be reported");
    near(FluidProperties::liquid_solution_density(0, p, 350, thermal), 950,
         "Buoyancy slope must not change the constant liquid density");
    for (double slope : {0., 1., std::numeric_limits<double>::quiet_NaN(),
                         -std::numeric_limits<double>::infinity()}) {
        thermal.density_temperature_slope = slope;
        errors.clear();
        thermal.validate(errors, warnings);
        check(!errors.empty(), "Constant-density end convection rejects invalid buoyancy slopes");
    }
    thermal.density_temperature_slope = -.7284;
    thermal.density_model = DensityModel::MIXTURE;
    thermal.dissolved_gas_liquid_density = 1000;
    errors.clear();
    thermal.validate(errors, warnings);
    check(std::any_of(errors.begin(), errors.end(), [](const auto& error) {
        return error.find("Housing end natural convection requires") != std::string::npos;
    }), "End convection must still reject unsupported density models");
    fixed.solubility_model = SolubilityModel::CURVE_TABLE;
    errors.clear();
    fixed.validate(errors, warnings);
    check(!errors.empty(), "Viscosity-only input cannot support CURVE_TABLE solubility");

    // Failure must clear both tables, including when reloading a previously valid file.
    for (const char* row : {".1,40,0.1,NaN", ".1,40,-1,10", "1.1,40,1,10", ".1,40,1,10junk"}) {
        { std::ofstream out(bad); out << "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n" << row << '\n'; }
        std::string error;
        PropertyCurveTable solubility = cfg.solubility_curve_table, viscosity = table;
        FluidProperties::load_solution_table_csv(bad, solubility, viscosity, error);
        check(!error.empty() && solubility.empty() && viscosity.empty(), "Bad CSV must fail without stale data");
    }
    for (const char* rows : {"0,30,,8\n", "0,30,,8\n0,30,,2\n",
                            ".3,30,,8\n.3,40,,2\n",
                            "0,30,,8\n0,40,,2\n.3,30,1,4\n.3,40,1,1\n",
                            ".3,30,1,8\n.3,40,2,2\n"}) {
        { std::ofstream out(bad); out << "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n" << rows; }
        std::string error;
        PropertyCurveTable solubility = cfg.solubility_curve_table, viscosity = table;
        FluidProperties::load_solution_table_csv(bad, solubility, viscosity, error);
        check(!error.empty() && solubility.empty() && viscosity.empty(),
              "Incomplete viscosity or partial solubility input must retain strict rejection");
    }
    std::filesystem::remove(bad);
    std::cout << "PASS: FW68D source nodes, interpolation, density, SI conversion, bounds and config\n";
}
