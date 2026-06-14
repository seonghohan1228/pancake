#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "config.hpp"
#include "field.hpp"
#include "fluid_properties.hpp"
#include "mesh.hpp"

namespace {

void check(bool condition, const char* message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

bool contains_text(const std::vector<std::string>& values, const std::string& needle) {
    return std::any_of(values.begin(), values.end(), [&](const std::string& value) {
        return value.find(needle) != std::string::npos;
    });
}

Fields make_property_fields(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields fields;
    fields.add("pressure", mesh).fill(cfg.initial_pressure());
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);
    fields.add("mu", mesh).fill(cfg.mu);
    fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
    fields.add("mu_liquid_solution", mesh).fill(cfg.mu);
    fields.add("cp_liquid_solution", mesh).fill(cfg.rho_cp / cfg.rho);
    fields.add("rho_cp_liquid_solution", mesh).fill(cfg.rho_cp);
    fields.add("k_liquid_solution", mesh).fill(cfg.thermal_conductivity);
    fields.add("temperature", mesh).fill(cfg.temperature_initial);
    fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
    fields.add("free_gas_mass", mesh).fill(0.0);
    fields.add("alpha_gas", mesh).fill(0.0);
    fields.add("gas_mass_transfer", mesh).fill(0.0);
    fields.add("rho_gas", mesh).fill(0.0);
    fields.add("inlet_indicator", mesh).fill(0.0);
    return fields;
}

SimulationConfig make_valid_r290_validation_config() {
    SimulationConfig cfg;
    cfg.p_cav = 1.0e5;
    cfg.bc_z_south_val = 1.0e6;
    cfg.bc_z_north_val = 1.0e6;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.dissolved_gas_species = DissolvedGasSpecies::PROPANE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.density_model = DensityModel::MASS_VOLUME_MIXING;
    cfg.viscosity_model = ViscosityModel::EMPIRICAL_CORRELATION;
    cfg.dissolved_gas_initial = 0.2175;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.175e-7;
    cfg.dissolved_gas_henry_temp_coeff = 1800.0;
    cfg.dissolved_gas_liquid_density = 500.0;
    cfg.gas_mass_transfer_rate = 500.0;
    cfg.dissolved_gas_diffusivity = 1.0e-9;
    cfg.solution_viscosity_gas_coeff = -11.40;
    cfg.viscosity_temperature_coeff = 4000.0;
    cfg.viscosity_pressure_coeff = 1.5e-8;
    cfg.gas_pressure_floor = 1.0e5;
    cfg.gas_alpha_max = 0.6;
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::MCADAMS_QUALITY;
    cfg.mu_gas = 8.7e-6;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;
    return cfg;
}

void test_gas_config_validation() {
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    SimulationConfig cfg = make_valid_r290_validation_config();
    cfg.validate(errors, warnings);
    check(errors.empty(), "valid R290/PZ68 gaseous cavitation config should pass validation");

    SimulationConfig zero_pressure = cfg;
    zero_pressure.p_cav = 0.0;
    zero_pressure.validate(errors, warnings);
    check(contains_text(errors, "p_cav"), "validation should reject gauge-zero p_cav");

    SimulationConfig missing_henry = cfg;
    missing_henry.dissolved_gas_henry_coeff = 0.0;
    missing_henry.validate(errors, warnings);
    check(contains_text(errors, "HENRY"), "validation should reject missing Henry coefficient");

    SimulationConfig non_jfo = cfg;
    non_jfo.cavitation_model = CavitationModel::GUMBEL;
    non_jfo.validate(errors, warnings);
    check(contains_text(errors, "ELROD_ADAMS"), "gas cavitation should require Elrod-Adams/JFO");
    check(contains_text(warnings, "GUMBEL"), "Gumbel should be visibly warned as non-mass-conserving");

    SimulationConfig no_initial_gas = cfg;
    no_initial_gas.dissolved_gas_initial = 0.0;
    no_initial_gas.validate(errors, warnings);
    check(errors.empty(), "zero initial gas should be a warning, not a hard error");
    check(contains_text(warnings, "zero initial dissolved gas"),
          "zero initial gas warning should explain why no gaseous cavitation occurs");
}

void test_config_parsing_and_tables() {
    SimulationConfig cfg;
    cfg.load_from_text(
        "fluid_property_model = gas cavitation mixture\n"
        "dissolved_gas_species = R290\n"
        "oil_gas_solution_model = table\n"
        "density_model = mass volume mixing\n"
        "viscosity_model = empirical\n"
        "solubility_table = 100000:0.01, 200000:0.02\n");

    check(cfg.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE,
          "fluid_property_model alias should parse");
    check(cfg.dissolved_gas_species == DissolvedGasSpecies::PROPANE,
          "R290 alias should parse as propane");
    check(cfg.oil_gas_solution_model == OilGasSolutionModel::TABLE,
          "table solution model should parse");
    check(cfg.density_model == DensityModel::MASS_VOLUME_MIXING,
          "density model alias should parse");
    check(cfg.viscosity_model == ViscosityModel::EMPIRICAL_CORRELATION,
          "viscosity model alias should parse");
    check(cfg.solubility_table.size() == 2,
          "property table should parse x:value pairs");
    check(std::abs(FluidProperties::interpolate_property_table(cfg.solubility_table, 150000.0, 0.0) - 0.015) < 1.0e-14,
          "property table interpolation should be linear");
}

void test_pure_oil_fallback() {
    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 2;
    cfg.rho = 960.0;
    cfg.mu = 0.12;
    cfg.dissolved_gas_initial = 0.2;

    Mesh mesh(cfg);
    Fields fields = make_property_fields(mesh, cfg);
    FluidProperties::update_solution_fields(fields, mesh, cfg);

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            check(std::abs(fields["rho"](i, j) - cfg.rho) < 1.0e-14,
                  "constant model should keep pure-oil density");
            check(std::abs(fields["mu"](i, j) - cfg.mu) < 1.0e-14,
                  "constant model should keep pure-oil viscosity");
            check(std::abs(fields["dissolved_gas"](i, j)) < 1.0e-14,
                  "constant model should disable dissolved gas");
        }
    }
}

void test_solution_property_models() {
    SimulationConfig cfg;
    cfg.fluid_property_model = FluidPropertyModel::OIL_DISSOLVED_GAS;
    cfg.rho = 960.0;
    cfg.mu = 0.10;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 1.0e-8;

    const double c_low = FluidProperties::equilibrium_dissolved_gas(1.0e5, 300.0, cfg);
    const double c_high = FluidProperties::equilibrium_dissolved_gas(2.0e5, 300.0, cfg);
    check(c_high > c_low, "Henry dissolved gas concentration should increase with pressure");

    cfg.density_model = DensityModel::MASS_VOLUME_MIXING;
    const double gas_rho = FluidProperties::gas_density(1.0e5, 300.0, cfg);
    const double mixture_rho = FluidProperties::liquid_solution_density(0.03, 1.0e5, 300.0, cfg);
    check(mixture_rho > gas_rho && mixture_rho < cfg.rho,
          "mass-volume density should remain between gas and pure oil density");

    cfg.viscosity_model = ViscosityModel::LOG_MIXING;
    cfg.solution_viscosity_gas_coeff = -4.0;
    const double mixture_mu = FluidProperties::liquid_solution_viscosity(0.05, 1.0e5, 300.0, cfg);
    check(mixture_mu > 0.0 && mixture_mu < cfg.mu,
          "log viscosity model should apply gas-concentration coefficient");

    cfg.viscosity_model = ViscosityModel::TABLE;
    cfg.viscosity_table = {{0.0, 0.10}, {0.10, 0.05}};
    check(std::abs(FluidProperties::liquid_solution_viscosity(0.05, 1.0e5, 300.0, cfg) - 0.075) < 1.0e-14,
          "table viscosity should interpolate by dissolved gas concentration");
}

void test_r290_pz68_reference_point() {
    SimulationConfig cfg;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.dissolved_gas_species = DissolvedGasSpecies::PROPANE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.density_model = DensityModel::MASS_VOLUME_MIXING;
    cfg.viscosity_model = ViscosityModel::EMPIRICAL_CORRELATION;
    cfg.rho = 960.0;
    cfg.mu = 0.0653;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.175e-7;
    cfg.dissolved_gas_henry_temp_coeff = 1800.0;
    cfg.dissolved_gas_liquid_density = 500.0;
    cfg.solution_viscosity_gas_coeff = -11.40;
    cfg.viscosity_temperature_coeff = 4000.0;
    cfg.viscosity_pressure_coeff = 1.5e-8;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;

    const double pressure = 1.0e6;
    const double temperature = 313.15;
    const double c_ref = FluidProperties::equilibrium_dissolved_gas(pressure, temperature, cfg);
    check(std::abs(c_ref - 0.2175) < 1.0e-12,
          "R290/PZ68 Henry coefficient should reproduce 21.75% at 40 C and 10 bar");

    const double rho_l = FluidProperties::liquid_solution_density(c_ref, pressure, temperature, cfg);
    const double mu_l = FluidProperties::liquid_solution_viscosity(c_ref, pressure, temperature, cfg);
    const double nu_cst = mu_l / rho_l * 1.0e6;
    check(std::abs(nu_cst - 6.843) < 0.05,
          "R290/PZ68 viscosity calibration should reproduce the supplied 40 C, 10 bar chart point");

    const double c_hot = FluidProperties::equilibrium_dissolved_gas(pressure, 333.15, cfg);
    check(c_hot < c_ref,
          "positive Henry temperature coefficient should reduce propane solubility at higher temperature");
}

void test_release_and_resorption() {
    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 2;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.dissolved_gas_henry_coeff = 1.0e-8;
    cfg.dissolved_gas_max = 0.5;
    cfg.gas_mass_transfer_rate = 20.0;
    cfg.c = 1.0e-4;
    cfg.rho = 960.0;
    cfg.gas_alpha_max = 1.0;

    Mesh mesh(cfg);
    Fields fields = make_property_fields(mesh, cfg);
    fields["pressure"].fill(1.0e5);
    fields["theta"].fill(0.5);
    fields["dissolved_gas"].fill(FluidProperties::equilibrium_dissolved_gas(1.0e5, 300.0, cfg));
    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);
    check(std::abs(fields["gas_mass_transfer"](0, 0)) < 1.0e-14,
          "saturated state should have zero gas mass transfer");

    fields["dissolved_gas"].fill(0.02);
    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);
    check(fields["gas_mass_transfer"](0, 0) > 0.0,
          "supersaturated liquid should release gas at lower pressure");
    check(fields["alpha_gas"](0, 0) >= 0.0 && fields["alpha_gas"](0, 0) <= cfg.gas_alpha_max + 1.0e-14,
          "released free gas should be bounded by gas_alpha_max");
    const double alpha_after_release = fields["alpha_gas"](0, 0);

    fields["pressure"].fill(3.0e6);
    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);
    check(fields["gas_mass_transfer"](0, 0) < 0.0,
          "pressure recovery should resorb free gas");
    check(fields["alpha_gas"](0, 0) < alpha_after_release,
          "resorption should reduce free gas volume fraction");
}

void test_free_gas_is_not_deleted_in_full_film() {
    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 2;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.dissolved_gas_henry_coeff = 2.0e-7;
    cfg.dissolved_gas_max = 0.5;
    cfg.gas_mass_transfer_rate = 20.0;
    cfg.c = 1.0e-4;
    cfg.rho = 960.0;
    cfg.gas_alpha_max = 0.5;
    cfg.gas_pressure_floor = 1.0e5;

    Mesh mesh(cfg);
    Fields fields = make_property_fields(mesh, cfg);
    fields["pressure"].fill(1.0e6);
    fields["theta"].fill(1.001);
    fields["dissolved_gas"].fill(FluidProperties::equilibrium_dissolved_gas(1.0e6, 300.0, cfg));
    const double rho_g = FluidProperties::gas_density(1.0e6, 300.0, cfg);
    const double initial_free_mass = 0.1 * rho_g * cfg.c;
    fields["free_gas_mass"].fill(initial_free_mass);

    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);
    check(std::abs(fields["free_gas_mass"](0, 0) - initial_free_mass) / initial_free_mass < 1.0e-12,
          "saturated full-film cells must not delete free gas");
    check(std::abs(fields["gas_mass_transfer"](0, 0)) < 1.0e-12,
          "saturated full-film free gas should have zero transfer, not silent deletion");
    check(fields["alpha_gas"](0, 0) > 0.0,
          "full-film free gas should remain visible as alpha_gas");

    fields["dissolved_gas"].fill(0.5 * FluidProperties::equilibrium_dissolved_gas(1.0e6, 300.0, cfg));
    fields["free_gas_mass"].fill(initial_free_mass);
    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);
    check(fields["gas_mass_transfer"](0, 0) < 0.0,
          "undersaturated full-film liquid should resorb free gas with a negative source");
    check(fields["free_gas_mass"](0, 0) < initial_free_mass,
          "finite-rate resorption should reduce free gas mass");
}

void test_pressure_inlet_imposes_initial_gas_state() {
    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 2;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.dissolved_gas_initial = 0.20;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.0e-7;
    cfg.dissolved_gas_henry_temp_coeff = 2000.0;
    cfg.property_reference_temperature = 300.0;
    cfg.gas_mass_transfer_rate = 20.0;
    cfg.c = 1.0e-4;
    cfg.rho = 960.0;
    cfg.gas_alpha_max = 1.0;
    cfg.gas_pressure_floor = 1.0e5;

    Mesh mesh(cfg);
    Fields fields = make_property_fields(mesh, cfg);
    fields["pressure"].fill(1.0e6);
    fields["temperature"].fill(330.0);
    fields["theta"].fill(1.0);
    fields["dissolved_gas"].fill(cfg.dissolved_gas_initial);
    fields["free_gas_mass"].fill(1.0e-6);
    fields["inlet_indicator"](0, 0) = 1.0;

    FluidProperties::update_solution_fields(fields, mesh, cfg);
    FluidProperties::update_gas_state(fields, mesh, cfg, 0.01);

    check(std::abs(fields["dissolved_gas"](0, 0) - cfg.dissolved_gas_initial) < 1.0e-14,
          "pressure inlet cell should keep the configured initial dissolved gas fraction");
    check(std::abs(fields["free_gas_mass"](0, 0)) < 1.0e-14,
          "pressure inlet cell should impose zero incoming free gas");
    check(std::abs(fields["alpha_gas"](0, 0)) < 1.0e-14,
          "pressure inlet cell should impose zero incoming free-gas volume fraction");
    check(std::abs(fields["gas_mass_transfer"](0, 0)) < 1.0e-14,
          "pressure inlet cell should not report local gas release or resorption");
    check(fields["gas_mass_transfer"](0, 1) > 0.0,
          "non-inlet supersaturated cell should still release gas");
}

void test_gas_mixture_viscosity_models() {
    SimulationConfig cfg;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.mu = 0.05;
    cfg.mu_gas = 1.0e-5;
    cfg.gas_alpha_max = 1.0;
    const double mu_l = cfg.mu;

    const GasMixtureViscosityModel all_models[] = {
        GasMixtureViscosityModel::EINSTEIN_DILUTE,
        GasMixtureViscosityModel::DUKLER_VOID,
        GasMixtureViscosityModel::MCADAMS_QUALITY,
        GasMixtureViscosityModel::KRIEGER_DOUGHERTY,
        GasMixtureViscosityModel::LINEAR_QUALITY};
    for (GasMixtureViscosityModel model : all_models) {
        cfg.gas_mixture_viscosity_model = model;
        check(std::abs(FluidProperties::gas_mixture_viscosity(mu_l, 0.0, 0.0, cfg) - mu_l) <
                  1.0e-15 * mu_l,
              "every mixture viscosity model must return mu_liquid at alpha_gas = 0");
    }

    // Gas-limit asymptotes: void/quality-weighted models thin toward mu_gas.
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::DUKLER_VOID;
    check(std::abs(FluidProperties::gas_mixture_viscosity(mu_l, 1.0, 0.0, cfg) - cfg.mu_gas) <
              1.0e-12,
          "DUKLER_VOID must reach mu_gas at alpha_gas = 1");
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::MCADAMS_QUALITY;
    check(std::abs(FluidProperties::gas_mixture_viscosity(mu_l, 1.0, 1.0, cfg) - cfg.mu_gas) <
              1.0e-12,
          "MCADAMS_QUALITY must reach mu_gas at quality = 1");
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::LINEAR_QUALITY;
    check(std::abs(FluidProperties::gas_mixture_viscosity(mu_l, 1.0, 1.0, cfg) - cfg.mu_gas) <
              1.0e-12,
          "LINEAR_QUALITY must reach mu_gas at quality = 1");
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::KRIEGER_DOUGHERTY;
    const double kd_limit = FluidProperties::gas_mixture_viscosity(mu_l, 1.0, 0.0, cfg);
    check(std::isfinite(kd_limit) && kd_limit > 10.0 * mu_l,
          "KRIEGER_DOUGHERTY must stay finite and diverge upward near packing");

    // Monotonicity: thinning models decrease with gas content, Einstein is the
    // documented dilute exception (increases).
    auto monotone_decreasing = [&](GasMixtureViscosityModel model, bool use_quality) {
        cfg.gas_mixture_viscosity_model = model;
        double previous = mu_l + 1.0;
        for (double x = 0.0; x <= 1.0 + 1.0e-12; x += 0.1) {
            const double mu_mix = use_quality
                ? FluidProperties::gas_mixture_viscosity(mu_l, 0.5, x, cfg)
                : FluidProperties::gas_mixture_viscosity(mu_l, x, 0.0, cfg);
            if (mu_mix > previous + 1.0e-15) return false;
            previous = mu_mix;
        }
        return true;
    };
    check(monotone_decreasing(GasMixtureViscosityModel::DUKLER_VOID, false),
          "DUKLER_VOID must decrease monotonically with alpha_gas");
    check(monotone_decreasing(GasMixtureViscosityModel::MCADAMS_QUALITY, true),
          "MCADAMS_QUALITY must decrease monotonically with quality");
    check(monotone_decreasing(GasMixtureViscosityModel::LINEAR_QUALITY, true),
          "LINEAR_QUALITY must decrease monotonically with quality");
    cfg.gas_mixture_viscosity_model = GasMixtureViscosityModel::EINSTEIN_DILUTE;
    check(FluidProperties::gas_mixture_viscosity(mu_l, 0.5, 0.0, cfg) > mu_l,
          "EINSTEIN_DILUTE thickens with alpha_gas (documented dilute behavior)");

    // Validation guards.
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    SimulationConfig guard = make_valid_r290_validation_config();
    guard.gas_mixture_viscosity_model = GasMixtureViscosityModel::EINSTEIN_DILUTE;
    guard.gas_alpha_max = 1.0;
    guard.validate(errors, warnings);
    check(contains_text(errors, "EINSTEIN_DILUTE"),
          "EINSTEIN_DILUTE with gas_alpha_max > 0.6 must fail validation");

    guard = make_valid_r290_validation_config();
    guard.mu_gas = 0.0;
    guard.validate(errors, warnings);
    check(contains_text(errors, "mu_gas"),
          "quality/void-weighted mixture viscosity must require mu_gas > 0");

    guard = make_valid_r290_validation_config();
    guard.validate(errors, warnings);
    check(errors.empty(), "shipped MCADAMS_QUALITY configuration must validate cleanly");
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    test_config_parsing_and_tables();
    test_gas_config_validation();
    test_pure_oil_fallback();
    test_solution_property_models();
    test_r290_pz68_reference_point();
    test_release_and_resorption();
    test_free_gas_is_not_deleted_in_full_film();
    test_pressure_inlet_imposes_initial_gas_state();
    test_gas_mixture_viscosity_models();

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_fluid_properties\n";

    MPI_Finalize();
    return 0;
}
