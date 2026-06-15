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
#include "film_thickness.hpp"
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

static double global_min_value(double local) {
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global;
}

static double global_max_value(double local) {
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

static void test_steady_initial_solve_ignores_old_temperature() {
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
    fields["temperature"].fill(1.0);
    fields["temperature"].store_old_time();
    fields["temperature"].fill(900.0);
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    const double heat = cfg.mu * std::pow(cfg.omega * cfg.R, 2.0) / cfg.c;
    const double expected_temperature =
        cfg.journal_wall_temperature + heat / (cfg.journal_heat_transfer + cfg.bearing_heat_transfer);

    double local_temp_err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            local_temp_err = std::max(local_temp_err, std::abs(fields["temperature"](i, j) - expected_temperature));

    check(global_max_abs(local_temp_err) < 1e-6,
          "steady initial energy solve should ignore previous temperature guesses");
}

static void test_pressure_inlet_does_not_fix_temperature() {
    SimulationConfig cfg;
    cfg.n_theta_global = 36;
    cfg.n_z_global = 6;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.temperature_reference = 250.0;
    cfg.journal_wall_temperature = 300.0;
    cfg.bearing_wall_temperature = 300.0;
    cfg.journal_heat_transfer = 1000.0;
    cfg.bearing_heat_transfer = 1000.0;
    cfg.omega = 200.0;

    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0;
    groove.size = 20.0;
    groove.p_supply = 2.0e5;
    cfg.inlets = {groove};

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    Field& inlet_indicator = fields.add("inlet_indicator", mesh);
    FilmThickness::compute_inlet_indicator(inlet_indicator, mesh, cfg);
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double local_min_reference_delta = 1.0e300;
    int local_inlet_cells = 0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            if (inlet_indicator(i, j) > 0.5) {
                local_min_reference_delta = std::min(
                    local_min_reference_delta,
                    std::abs(fields["temperature"](i, j) - cfg.temperature_reference));
                ++local_inlet_cells;
            }
        }
    }

    int global_inlet_cells = 0;
    MPI_Allreduce(&local_inlet_cells, &global_inlet_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(global_inlet_cells > 0, "test setup should include inlet cells");
    check(global_min_value(local_min_reference_delta) > 10.0,
          "pressure inlet cells should not be fixed to a thermal Dirichlet value");
}

static void test_axial_inflow_temperature_is_bounded() {
    SimulationConfig cfg;
    cfg.n_theta_global = 12;
    cfg.n_z_global = 12;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.omega = 0.0;
    cfg.temperature_initial = 300.0;
    cfg.temperature_reference = 330.0;
    cfg.journal_wall_temperature = 300.0;
    cfg.bearing_wall_temperature = 300.0;
    cfg.journal_heat_transfer = 500.0;
    cfg.bearing_heat_transfer = 500.0;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 2.0e5;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 0.0;
    // ELROD derives its boundary film content from bc_z_*_val and bulk modulus;
    // request a fed reservoir so inflow carries temperature_reference.
    cfg.bc_z_south_thermal = ThermalInflowMode::CONSTANT;
    cfg.bc_z_north_thermal = ThermalInflowMode::CONSTANT;
    cfg.inlets.clear();

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double z_c = (j + 0.5) * mesh.get_d_z();
            fields["pressure"](i, j) =
                cfg.bc_z_south_val + (cfg.bc_z_north_val - cfg.bc_z_south_val) * z_c / cfg.L;
        }
    }
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double local_min = 1.0e300;
    double local_max = -1.0e300;
    double local_south_max = -1.0e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double value = fields["temperature"](i, j);
            local_min = std::min(local_min, value);
            local_max = std::max(local_max, value);
            if (j == 0) local_south_max = std::max(local_south_max, value);
        }
    }

    const double global_min = global_min_value(local_min);
    const double global_max = global_max_value(local_max);
    const double global_south_max = global_max_value(local_south_max);
    check(global_min >= cfg.journal_wall_temperature - 1e-8,
          "axial thermal inflow should not undershoot wall temperature");
    check(global_max <= cfg.temperature_reference + 1e-8,
          "axial thermal inflow should not overshoot reference inflow temperature without heat source");
    check(global_south_max > cfg.journal_wall_temperature + 1e-4,
          "south pressure inflow should carry temperature_reference into the energy solve");
}

static void test_north_axial_inflow_temperature_is_bounded() {
    SimulationConfig cfg;
    cfg.n_theta_global = 12;
    cfg.n_z_global = 12;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.omega = 0.0;
    cfg.temperature_initial = 300.0;
    cfg.temperature_reference = 330.0;
    cfg.journal_wall_temperature = 300.0;
    cfg.bearing_wall_temperature = 300.0;
    cfg.journal_heat_transfer = 500.0;
    cfg.bearing_heat_transfer = 500.0;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 0.0;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 2.0e5;
    cfg.bc_z_south_thermal = ThermalInflowMode::CONSTANT;
    cfg.bc_z_north_thermal = ThermalInflowMode::CONSTANT;
    cfg.inlets.clear();

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double z_c = (j + 0.5) * mesh.get_d_z();
            fields["pressure"](i, j) =
                cfg.bc_z_south_val + (cfg.bc_z_north_val - cfg.bc_z_south_val) * z_c / cfg.L;
        }
    }
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double local_min = 1.0e300;
    double local_max = -1.0e300;
    double local_north_max = -1.0e300;
    const int north_j = mesh.n_z_local - 1;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double value = fields["temperature"](i, j);
            local_min = std::min(local_min, value);
            local_max = std::max(local_max, value);
            if (j == north_j) local_north_max = std::max(local_north_max, value);
        }
    }

    check(global_min_value(local_min) >= cfg.journal_wall_temperature - 1e-8,
          "north axial thermal inflow should not undershoot wall temperature");
    check(global_max_value(local_max) <= cfg.temperature_reference + 1e-8,
          "north axial thermal inflow should not overshoot reference temperature without heat source");
    check(global_max_value(local_north_max) > cfg.journal_wall_temperature + 1e-4,
          "north pressure inflow should carry temperature_reference into the energy solve");
}

static void test_outflow_boundaries_do_not_impose_reference_temperature() {
    SimulationConfig cfg;
    cfg.n_theta_global = 12;
    cfg.n_z_global = 12;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.omega = 0.0;
    cfg.temperature_initial = 300.0;
    cfg.temperature_reference = 1000.0;
    cfg.journal_wall_temperature = 300.0;
    cfg.bearing_wall_temperature = 300.0;
    cfg.journal_heat_transfer = 1000.0;
    cfg.bearing_heat_transfer = 1000.0;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 0.0;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 0.0;
    // Request fed reservoirs so this exercises the RESERVOIR path; the interior
    // pressure stays above the boundary pressure, so both ends are outflow and
    // must NOT pull in temperature_reference.
    cfg.bc_z_south_thermal = ThermalInflowMode::CONSTANT;
    cfg.bc_z_north_thermal = ThermalInflowMode::CONSTANT;
    cfg.inlets.clear();

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double z_c = (j + 0.5) * mesh.get_d_z();
            fields["pressure"](i, j) = 5.0e3 * std::sin(3.14159265358979323846 * z_c / cfg.L);
        }
    }
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double local_max = -1.0e300;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            local_max = std::max(local_max, fields["temperature"](i, j));

    check(global_max_value(local_max) < 310.0,
          "axial outflow boundaries should not impose temperature_reference");
}

static void test_boundary_heat_is_bounded_and_uses_field_viscosity() {
    // The Poiseuille heat term uses a bounded, full-film-gated gradient: at the axial
    // ends it differences against the interior cell (NOT the half-cell-to-boundary
    // value), so a boundary cell does not get a spurious dissipation spike from the BC.
    // With a pressure field linear in z, dp/dz is the same constant in the interior and
    // at the boundary, so heat is uniform and equals c^3 (dp/dz)^2 / (12 mu_field).
    SimulationConfig cfg;
    cfg.n_theta_global = 12;
    cfg.n_z_global = 6;
    cfg.temperature_model = TemperatureModel::ISOTHERMAL;
    cfg.omega = 0.0;
    cfg.mu = 0.01;  // NOT used: the field mu below must take precedence
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 2.0e5;   // deliberately far from the interior field -> would spike if used
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 2.0e5;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    Field& mu = fields.add("mu", mesh);
    mu.fill(0.02);
    const double slope = 1.0e8;  // dp/dz [Pa/m]
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double z_c = (j + 0.5) * mesh.get_d_z();
            fields["pressure"](i, j) = 5.0e4 + slope * z_c;  // linear in z, uniform in theta
        }
    }
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    const double expected_heat = cfg.c * cfg.c * cfg.c * slope * slope / (12.0 * 0.02);
    double local_err = 0.0;
    double local_boundary_err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double err = std::abs(fields["heat_generation"](i, j) - expected_heat);
            local_err = std::max(local_err, err);
            if (j == 0 || j == mesh.n_z_local - 1) local_boundary_err = std::max(local_boundary_err, err);
        }
    }

    check(global_max_abs(local_err) / expected_heat < 1.0e-9,
          "Poiseuille heat should use the field viscosity and the bounded interior gradient");
    check(global_max_abs(local_boundary_err) / expected_heat < 1.0e-9,
          "boundary cells must match the interior (no half-cell-to-BC dissipation spike)");
}

static void test_cavitation_front_heat_uses_active_face_average() {
    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 1;
    cfg.temperature_model = TemperatureModel::ISOTHERMAL;
    cfg.omega = 0.0;
    cfg.c = 5.0e-5;
    cfg.mu = 0.02;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);

    fields["theta"].fill(0.5);
    fields["pressure"].fill(cfg.p_cav);
    const double p_left = 1.0e6;
    const double p_right = 2.0e6;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const int global_i = mesh.offset_theta + i;
        if (global_i == 1 || global_i == 2) fields["theta"](i, 0) = 1.1;
        if (global_i == 1) fields["pressure"](i, 0) = p_left;
        if (global_i == 2) fields["pressure"](i, 0) = p_right;
    }
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    const double dp_ds = (p_right - p_left) / (cfg.R * mesh.get_d_theta());
    const double expected_heat = cfg.c * cfg.c * cfg.c * 0.5 * dp_ds * dp_ds / (12.0 * cfg.mu);
    double local_err = 0.0;
    double local_cavitated_heat = 0.0;
    int local_front_cells = 0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const int global_i = mesh.offset_theta + i;
        const double heat = fields["heat_generation"](i, 0);
        if (global_i == 1 || global_i == 2) {
            local_err = std::max(local_err, std::abs(heat - expected_heat));
            ++local_front_cells;
        } else {
            local_cavitated_heat = std::max(local_cavitated_heat, std::abs(heat));
        }
    }

    int global_front_cells = 0;
    MPI_Allreduce(&local_front_cells, &global_front_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(global_front_cells == 2, "test setup should include two full-film front cells");
    check(global_max_abs(local_err) / expected_heat < 1.0e-12,
          "cavitation-front heat should average active face gradients");
    check(global_max_abs(local_cavitated_heat) < 1.0e-12,
          "cavitated cells should not receive pressure-flow heat");
}

static void test_elrod_axial_boundary_uses_pressure_value() {
    // Regression guard for the pressure-derived boundary contract. Under
    // ELROD_ADAMS the axial DIRICHLET boundary derives film content from
    // bc_z_*_val, p_cav, and bulk_modulus. Legacy bc_z_*_theta values are ignored.
    // A full-film interior below bc_z_*_val must therefore see physical axial
    // inflow from the reservoir.
    SimulationConfig cfg;
    cfg.n_theta_global = 12;
    cfg.n_z_global = 12;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.omega = 0.0;
    cfg.temperature_initial = 300.0;
    cfg.temperature_reference = 1000.0;  // hot reservoir, to expose any spurious inflow
    cfg.journal_wall_temperature = 360.0;
    cfg.bearing_wall_temperature = 360.0;
    cfg.journal_heat_transfer = 1000.0;
    cfg.bearing_heat_transfer = 1000.0;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 1.5e6;
    cfg.bc_z_north_val = 1.5e6;
    cfg.bc_z_south_theta = 1.0;          // deprecated: must not override bc_z_south_val
    cfg.bc_z_north_theta = 1.0;
    cfg.bc_z_south_thermal = ThermalInflowMode::CONSTANT;
    cfg.bc_z_north_thermal = ThermalInflowMode::CONSTANT;
    cfg.inlets.clear();

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_energy_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            fields["pressure"](i, j) = 5.0e4;  // > p_cav (full film) and < bc_z_*_val
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double local_max = -1.0e300;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            local_max = std::max(local_max, fields["temperature"](i, j));

    check(global_max_value(local_max) > 400.0,
          "ELROD axial boundary should use bc_z_val-derived pressure; legacy theta must not suppress inflow");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_energy_config_aliases();
    test_isothermal_refreshes_heat_only();
    test_steady_wall_balance();
    test_steady_initial_solve_ignores_old_temperature();
    test_pressure_inlet_does_not_fix_temperature();
    test_axial_inflow_temperature_is_bounded();
    test_north_axial_inflow_temperature_is_bounded();
    test_outflow_boundaries_do_not_impose_reference_temperature();
    test_boundary_heat_is_bounded_and_uses_field_viscosity();
    test_cavitation_front_heat_uses_active_face_average();
    test_elrod_axial_boundary_uses_pressure_value();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_energy\n";

    PetscFinalize();
    return 0;
}
