#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static double global_max_abs(double local) {
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

static void add_energy_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("temperature", mesh).fill(cfg.temperature_initial);
    fields.add("heat_generation", mesh).fill(0.0);
    fields["temperature"].store_old_time();
}

static void test_energy_config_aliases() {
    SimulationConfig cfg;
    cfg.load_from_text(
        "solution_mode = single step\n"
        "temperature_model = thd\n"
        "temperature_initial = 315\n"
        "journal_heat_transfer = 1000\n");

    check(cfg.solution_mode == SolutionMode::STEADY_STATE,
          "single step alias should map to STEADY_STATE");
    check(cfg.temperature_model == TemperatureModel::ENERGY_EQUATION,
          "THD alias should map to ENERGY_EQUATION");
    check(std::abs(cfg.temperature_initial - 315.0) < 1e-14,
          "temperature_initial parse failed");
    check(std::abs(cfg.journal_heat_transfer - 1000.0) < 1e-14,
          "journal_heat_transfer parse failed");
}

static void test_isothermal_refreshes_heat_only() {
    SimulationConfig cfg;
    cfg.n_theta_global = 16;
    cfg.n_z_global = 4;
    cfg.omega = 500.0;
    cfg.temperature_model = TemperatureModel::ISOTHERMAL;
    cfg.temperature_initial = 321.0;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    const double expected_heat = cfg.mu * std::pow(cfg.omega * cfg.R, 2.0) / cfg.c;
    double local_heat_err = 0.0;
    double local_temp_err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local_heat_err = std::max(local_heat_err, std::abs(fields["heat_generation"](i, j) - expected_heat));
            local_temp_err = std::max(local_temp_err, std::abs(fields["temperature"](i, j) - cfg.temperature_initial));
        }
    }

    check(global_max_abs(local_heat_err) < 1e-9,
          "isothermal solve should compute Couette heat generation");
    check(global_max_abs(local_temp_err) < 1e-12,
          "isothermal solve should not change temperature");
}

static void test_steady_wall_balance() {
    SimulationConfig cfg;
    cfg.n_theta_global = 16;
    cfg.n_z_global = 4;
    cfg.omega = 1000.0;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.temperature_initial = 300.0;
    cfg.journal_wall_temperature = 300.0;
    cfg.bearing_wall_temperature = 300.0;
    cfg.journal_heat_transfer = 1000.0;
    cfg.bearing_heat_transfer = 1000.0;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    const double heat = cfg.mu * std::pow(cfg.omega * cfg.R, 2.0) / cfg.c;
    const double expected_temperature =
        cfg.journal_wall_temperature + heat / (cfg.journal_heat_transfer + cfg.bearing_heat_transfer);

    double local_temp_err = 0.0;
    double local_heat_err = 0.0;
    double local_temp_min = 1.0e300;
    double local_temp_max = -1.0e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local_temp_err = std::max(local_temp_err, std::abs(fields["temperature"](i, j) - expected_temperature));
            local_heat_err = std::max(local_heat_err, std::abs(fields["heat_generation"](i, j) - heat));
            local_temp_min = std::min(local_temp_min, fields["temperature"](i, j));
            local_temp_max = std::max(local_temp_max, fields["temperature"](i, j));
        }
    }

    const double global_temp_err = global_max_abs(local_temp_err);
    if (global_temp_err >= 1e-6) {
        double global_temp_min = 0.0;
        double global_temp_max = 0.0;
        MPI_Allreduce(&local_temp_min, &global_temp_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_temp_max, &global_temp_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cerr << std::setprecision(16)
                      << "steady wall balance expected " << expected_temperature
                      << " got range [" << global_temp_min << ", " << global_temp_max
                      << "] max_err " << global_temp_err
                      << " heat " << heat << "\n";
        }
    }
    check(global_temp_err < 1e-6, "steady energy solve should match uniform wall heat balance");
    check(global_max_abs(local_heat_err) < 1e-9,
          "steady energy solve should preserve Couette heat generation");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_energy_config_aliases();
    test_isothermal_refreshes_heat_only();
    test_steady_wall_balance();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_energy\n";

    PetscFinalize();
    return 0;
}
