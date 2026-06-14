#include <algorithm>
#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "gas_transport.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static double global_max_abs(double local) {
    double g = 0.0; MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return g;
}
static double global_sum(double local) {
    double g = 0.0; MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g;
}

static SimulationConfig mixture_config() {
    SimulationConfig cfg;
    cfg.n_theta_global = 24;
    cfg.n_z_global = 8;
    cfg.omega = 200.0;
    cfg.mu = 0.05;
    cfg.rho = 950.0;
    cfg.p_cav = 1.0e5;
    cfg.bc_z_south_val = 1.0e5;
    cfg.bc_z_north_val = 1.0e5;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.dt = 1.0e-4;
    cfg.dissolved_gas_initial = 0.12;
    cfg.dissolved_gas_max = 1.0;
    cfg.dissolved_gas_diffusivity = 0.0;
    return cfg;
}

static void add_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    fields.add("pressure", mesh).fill(cfg.p_cav);   // uniform -> no Poiseuille, pure Couette
    fields.add("h", mesh).fill(cfg.c);
    fields.add("theta", mesh).fill(1.0);            // full film -> Couette flux active
    fields.add("rho", mesh).fill(cfg.rho);
    fields.add("mu", mesh).fill(cfg.mu);
    fields.add("dissolved_gas", mesh);
    fields.add("free_gas_mass", mesh).fill(0.0);
    fields.add("alpha_gas", mesh).fill(0.0);
    fields.add("gas_mass_transfer", mesh).fill(0.0);
    fields.add("inlet_indicator", mesh).fill(0.0);
}

// Uniform dissolved gas must be preserved by transport (constant-preserving scheme).
static void test_preserves_uniform() {
    SimulationConfig cfg = mixture_config();
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    fields["dissolved_gas"].fill(0.12);
    fields["dissolved_gas"].store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            err = std::max(err, std::abs(fields["dissolved_gas"](i, j) - 0.12));
    check(global_max_abs(err) < 1e-6, "uniform dissolved gas must be preserved by transport");
}

// No transport for non-mixture models.
static void test_noop_for_constant_model() {
    SimulationConfig cfg = mixture_config();
    cfg.fluid_property_model = FluidPropertyModel::CONSTANT;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            fields["dissolved_gas"](i, j) = 0.05 + 0.01 * i;
    fields["dissolved_gas"].store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            err = std::max(err, std::abs(fields["dissolved_gas"](i, j) - (0.05 + 0.01 * i)));
    check(global_max_abs(err) < 1e-14, "CONSTANT model must leave dissolved gas untouched");
}

// No transport without a physical timestep; this avoids assembling a singular
// steady advection equation for steady-state operating-point runs.
static void test_noop_without_timestep() {
    SimulationConfig cfg = mixture_config();
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            fields["dissolved_gas"](i, j) = 0.05 + 0.01 * i;
    fields["dissolved_gas"].store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, 0.0);

    double err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            err = std::max(err, std::abs(fields["dissolved_gas"](i, j) - (0.05 + 0.01 * i)));
    check(global_max_abs(err) < 1e-14, "zero timestep must leave dissolved gas untouched");
}

// A non-uniform field is advected (changes), stays bounded, and the mass-weighted
// total is conserved within numerics (neutral pressure ends, periodic theta, steady mass).
static void test_advects_and_conserves() {
    SimulationConfig cfg = mixture_config();
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    const double kPi = 3.14159265358979323846;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double th = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
        for (int j = 0; j < mesh.n_z_local; ++j)
            fields["dissolved_gas"](i, j) = 0.1 + 0.05 * std::sin(th);  // periodic, smooth
    }
    fields["dissolved_gas"].store_old_time();
    comm.update_ghosts(fields);

    const double vol = mesh.R * mesh.get_d_theta() * mesh.get_d_z();
    double m0_local = 0.0, change_local = 0.0, min_v = 1e9, max_v = -1e9;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            m0_local += fields["rho"](i, j) * fields["h"](i, j) * fields["dissolved_gas"](i, j) * vol;

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double m1_local = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double v = fields["dissolved_gas"](i, j);
            m1_local += fields["rho"](i, j) * fields["h"](i, j) * v * vol;
            change_local = std::max(change_local, std::abs(v - fields["dissolved_gas"].old(i, j)));
            min_v = std::min(min_v, v);
            max_v = std::max(max_v, v);
        }

    const double m0 = global_sum(m0_local);
    const double m1 = global_sum(m1_local);
    check(global_max_abs(change_local) > 1e-6, "advection should move a non-uniform dissolved-gas field");
    check(std::abs(m1 - m0) / std::abs(m0) < 1e-3, "transport should conserve mass-weighted dissolved gas");
    double gmin = 0.0, gmax = 0.0;
    MPI_Allreduce(&min_v, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&max_v, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(gmin >= -1e-9 && gmax <= cfg.dissolved_gas_max + 1e-9, "dissolved gas must stay bounded");
}

static void test_free_gas_mass_advects_and_conserves() {
    SimulationConfig cfg = mixture_config();
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    fields["dissolved_gas"].fill(0.1);

    const double kPi = 3.14159265358979323846;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double th = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
        for (int j = 0; j < mesh.n_z_local; ++j)
            fields["free_gas_mass"](i, j) = 2.0e-6 + 1.0e-6 * std::sin(th);
    }
    fields["dissolved_gas"].store_old_time();
    fields["free_gas_mass"].store_old_time();
    comm.update_ghosts(fields);

    const double area = mesh.R * mesh.get_d_theta() * mesh.get_d_z();
    double m0_local = 0.0, change_local = 0.0, min_v = 1e9;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            m0_local += fields["free_gas_mass"](i, j) * area;

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double m1_local = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double v = fields["free_gas_mass"](i, j);
            m1_local += v * area;
            change_local = std::max(change_local, std::abs(v - fields["free_gas_mass"].old(i, j)));
            min_v = std::min(min_v, v);
        }

    const double m0 = global_sum(m0_local);
    const double m1 = global_sum(m1_local);
    check(global_max_abs(change_local) > 1e-10, "advection should move a non-uniform free-gas field");
    check(std::abs(m1 - m0) / std::abs(m0) < 1e-3, "transport should conserve free gas mass");
    double gmin = 0.0;
    MPI_Allreduce(&min_v, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    check(gmin >= -1e-14, "free gas mass must stay non-negative");
}

static void test_supply_cells_impose_initial_composition() {
    SimulationConfig cfg = mixture_config();
    cfg.omega = 0.0;
    cfg.dissolved_gas_initial = 0.23;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    fields["dissolved_gas"].fill(0.04);
    fields["free_gas_mass"].fill(2.0e-6);
    fields["alpha_gas"].fill(0.35);
    fields["gas_mass_transfer"].fill(1.0);
    fields["inlet_indicator"](0, 0) = 1.0;
    fields["dissolved_gas"].store_old_time();
    fields["free_gas_mass"].store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    check(std::abs(fields["dissolved_gas"](0, 0) - cfg.dissolved_gas_initial) < 1.0e-14,
          "supply cell should impose initial dissolved gas composition");
    check(std::abs(fields["free_gas_mass"](0, 0)) < 1.0e-14,
          "supply cell should impose zero free gas in the entering mixture");
    check(std::abs(fields["alpha_gas"](0, 0)) < 1.0e-14,
          "supply cell should reset free-gas volume fraction");
    check(std::abs(fields["gas_mass_transfer"](0, 0)) < 1.0e-14,
          "supply cell should reset gas mass-transfer diagnostic");
}

static void test_axial_inflow_uses_initial_dissolved_gas() {
    SimulationConfig cfg = mixture_config();
    cfg.omega = 0.0;
    cfg.dissolved_gas_initial = 0.25;
    cfg.bc_z_south_val = 2.0e5;
    cfg.bc_z_north_val = 2.0e5;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    fields["pressure"].fill(1.0e5);
    fields["dissolved_gas"].fill(0.05);
    fields["dissolved_gas"].store_old_time();
    fields["free_gas_mass"].store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double south_change = 0.0;
    double north_change = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        south_change = std::max(south_change, fields["dissolved_gas"](i, 0) - 0.05);
        north_change = std::max(north_change, fields["dissolved_gas"](i, mesh.n_z_local - 1) - 0.05);
    }
    check(global_max_abs(south_change) > 1.0e-8,
          "south axial inflow should carry dissolved_gas_initial into the domain");
    check(global_max_abs(north_change) > 1.0e-8,
          "north axial inflow should carry dissolved_gas_initial into the domain");
}

static void test_axial_outflow_removes_free_gas() {
    SimulationConfig cfg = mixture_config();
    cfg.omega = 0.0;
    cfg.bc_z_south_val = 1.0e5;
    cfg.bc_z_north_val = 1.0e5;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    add_fields(fields, mesh, cfg);
    fields["pressure"].fill(2.0e5);
    fields["dissolved_gas"].fill(0.1);
    fields["free_gas_mass"].fill(2.0e-6);
    fields["dissolved_gas"].store_old_time();
    fields["free_gas_mass"].store_old_time();
    comm.update_ghosts(fields);

    const double area = mesh.R * mesh.get_d_theta() * mesh.get_d_z();
    double m0_local = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            m0_local += fields["free_gas_mass"](i, j) * area;

    LinearSystem sys(mesh);
    GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);

    double m1_local = 0.0;
    double min_v = 1e9;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            m1_local += fields["free_gas_mass"](i, j) * area;
            min_v = std::min(min_v, fields["free_gas_mass"](i, j));
        }

    const double m0 = global_sum(m0_local);
    const double m1 = global_sum(m1_local);
    check(m1 < m0, "axial pressure outflow should remove transported free gas");
    double gmin = 0.0;
    MPI_Allreduce(&min_v, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    check(gmin >= -1e-14, "axial free gas outflow should keep free_gas_mass non-negative");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    test_preserves_uniform();
    test_noop_for_constant_model();
    test_noop_without_timestep();
    test_advects_and_conserves();
    test_free_gas_mass_advects_and_conserves();
    test_supply_cells_impose_initial_composition();
    test_axial_inflow_uses_initial_dissolved_gas();
    test_axial_outflow_removes_free_gas();
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_gas_transport\n";
    PetscFinalize();
    return 0;
}
