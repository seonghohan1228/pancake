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

    // Failure must clear both tables, including when reloading a previously valid file.
    const auto bad = std::filesystem::temp_directory_path() / "pancake_fw68d_invalid.csv";
    for (const char* row : {".1,40,0.1,NaN", ".1,40,-1,10", "1.1,40,1,10", ".1,40,1,10junk"}) {
        { std::ofstream out(bad); out << "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n" << row << '\n'; }
        std::string error;
        PropertyCurveTable solubility = cfg.solubility_curve_table, viscosity = table;
        FluidProperties::load_solution_table_csv(bad, solubility, viscosity, error);
        check(!error.empty() && solubility.empty() && viscosity.empty(), "Bad CSV must fail without stale data");
    }
    std::filesystem::remove(bad);
    std::cout << "PASS: FW68D source nodes, interpolation, density, SI conversion, bounds and config\n";
}
