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

bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

bool has_substr(const std::vector<std::string>& items, const std::string& needle) {
    return std::any_of(items.begin(), items.end(),
                       [&](const std::string& s) { return s.find(needle) != std::string::npos; });
}

SimulationConfig gas_steady_config() {
    SimulationConfig cfg;
    cfg.p_cav = 1.0e5;
    cfg.bc_z_south_val = 1.0e6;
    cfg.bc_z_north_val = 1.0e6;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    cfg.solubility_model = SolubilityModel::HENRY;
    cfg.density_model = DensityModel::CONSTANT;
    cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
    cfg.dissolved_gas_initial = 0.2175;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.175e-7;
    cfg.gas_mass_transfer_rate = 500.0;
    cfg.gas_alpha_max = 0.6;
    cfg.mu_gas = 8.7e-6;
    cfg.mixture_viscosity_model = MixtureViscosityModel::MCADAMS;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;
    cfg.temperature_initial = 313.15;
    return cfg;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // --- steady + gas + FROZEN honesty guard --------------------------------
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig frozen = gas_steady_config();
        frozen.steady_gas_model = SteadyGasModel::FROZEN;
        frozen.allow_frozen_gas_steady = false;
        frozen.validate(errors, warnings);
        check(has_substr(errors, "freezes the gas model"),
              "guard: STEADY_STATE + gas + FROZEN must error");

        SimulationConfig equil = gas_steady_config();
        equil.steady_gas_model = SteadyGasModel::EQUILIBRIUM;
        errors.clear(); warnings.clear();
        equil.validate(errors, warnings);
        check(!has_substr(errors, "freezes the gas model"), "guard: EQUILIBRIUM is allowed");
        check(errors.empty(), "guard: EQUILIBRIUM steady gas config is otherwise valid");

        SimulationConfig allowed = gas_steady_config();
        allowed.steady_gas_model = SteadyGasModel::FROZEN;
        allowed.allow_frozen_gas_steady = true;
        errors.clear(); warnings.clear();
        allowed.validate(errors, warnings);
        check(!has_substr(errors, "freezes the gas model"),
              "guard: allow_frozen_gas_steady opts out of the error");
    }

    // --- WP-6 cavitated_film_model = STRIATED requires Elrod-Adams ----------
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig striated = gas_steady_config();
        striated.steady_gas_model = SteadyGasModel::EQUILIBRIUM;  // avoid the unrelated steady guard
        striated.cavitated_film_model = CavitatedFilmModel::STRIATED;
        striated.cavitation_model = CavitationModel::GUMBEL;
        striated.validate(errors, warnings);
        check(has_substr(errors, "STRIATED requires cavitation_model = JFO"),
              "WP-6: STRIATED + GUMBEL must error");

        striated.cavitation_model = CavitationModel::JFO;
        errors.clear(); warnings.clear();
        striated.validate(errors, warnings);
        check(!has_substr(errors, "STRIATED requires"), "WP-6: STRIATED + JFO is valid");
    }

    // --- coupling config round-trip -----------------------------------------
    {
        SimulationConfig cfg;
        cfg.load_from_text(
            "coupling_max_iters = 12\n"
            "coupling_tolerance = 1e-5\n"
            "coupling_relaxation = 0.7\n"
            "steady_gas_model = EQUILIBRIUM\n"
            "allow_frozen_gas_steady = true\n");
        check(cfg.coupling_max_iters == 12, "parse coupling_max_iters");
        check(close(cfg.coupling_tolerance, 1e-5, 1e-12), "parse coupling_tolerance");
        check(close(cfg.coupling_relaxation, 0.7, 1e-12), "parse coupling_relaxation");
        check(cfg.steady_gas_model == SteadyGasModel::EQUILIBRIUM, "parse steady_gas_model");
        check(cfg.allow_frozen_gas_steady, "parse allow_frozen_gas_steady");
    }

    // --- equilibrate_gas local flash ----------------------------------------
    {
        SimulationConfig cfg = gas_steady_config();
        cfg.steady_gas_model = SteadyGasModel::EQUILIBRIUM;
        cfg.n_theta_global = 8;
        cfg.n_z_global = 4;
        cfg.R = 0.02;
        cfg.c = 4.25e-5;
        cfg.L = 0.04;
        Mesh mesh(cfg);

        Fields fields;
        fields.add("pressure", mesh).fill(5.0e5);  // below p_sat (1.0 MPa) -> outgassing
        fields.add("temperature", mesh).fill(313.15);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("theta", mesh).fill(1.0);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("gas_mass_transfer", mesh).fill(0.0);
        fields.add("rho_gas", mesh).fill(0.0);
        fields.add("inlet_indicator", mesh).fill(0.0);

        FluidProperties::equilibrate_gas(fields, mesh, cfg);
        // c_eq(5e5) = 2.175e-7 * 5e5 = 0.108750 < feed -> some R290 comes out of solution.
        check(close(fields["dissolved_gas"](0, 0), 0.10875, 1e-5),
              "equilibrate_gas: dissolved = local saturation below p_sat");
        check(fields["free_gas_mass"](0, 0) > 1e-9, "equilibrate_gas: free gas released");
        check(fields["alpha_gas"](0, 0) > 0.0 && fields["alpha_gas"](0, 0) <= cfg.gas_alpha_max + 1e-12,
              "equilibrate_gas: void fraction within [0, alpha_max]");

        // At/above p_sat the feed stays fully dissolved (no free gas).
        fields["pressure"].fill(1.0e6);
        fields["dissolved_gas"].fill(cfg.dissolved_gas_initial);
        fields["free_gas_mass"].fill(0.0);
        fields["alpha_gas"].fill(0.0);
        FluidProperties::equilibrate_gas(fields, mesh, cfg);
        check(close(fields["dissolved_gas"](0, 0), cfg.dissolved_gas_initial, 1e-6),
              "equilibrate_gas: fully dissolved at p_sat");
        check(close(fields["free_gas_mass"](0, 0), 0.0, 1e-12), "equilibrate_gas: no free gas at p_sat");
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_coupling: OK\n";
    MPI_Finalize();
    return 0;
}
