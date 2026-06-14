#include <algorithm>
#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "equation_of_state.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

struct Range {
    double min_value = 1e300;
    double max_value = -1e300;
};

static Range global_range(const Field& field, const Mesh& mesh) {
    Range local;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local.min_value = std::min(local.min_value, field(i, j));
            local.max_value = std::max(local.max_value, field(i, j));
        }
    }

    Range global;
    MPI_Allreduce(&local.min_value, &global.min_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local.max_value, &global.max_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

static Fields make_fields(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta",    mesh).fill(1.0);
    fields.add("h",        mesh).fill(cfg.c);
    fields.add("rho",      mesh).fill(cfg.rho);
    FilmThickness::compute_static(fields["h"], mesh, cfg);
    return fields;
}

static void solve_elrod_once(Fields& fields, LinearSystem& sys, const Mesh& mesh, Communicator& comm, const SimulationConfig& cfg) {
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();
    fields["pressure"].store_old_time();
    Reynolds::solve_elrod(fields, sys, mesh, cfg);
    comm.update_ghosts(fields);
}

static SimulationConfig base_config() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
    cfg.n_theta_global = 32;
    cfg.n_z_global = 12;
    cfg.e = 0.0;
    cfg.omega = 0.0;
    cfg.bulk_modulus = 1e5;
    cfg.log_outer_iters = false;
    cfg.outer_tol = 1e-12;
    cfg.max_outer_iters = 20;
    cfg.inlets.clear();
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    return cfg;
}

static void test_elrod_dirichlet_derives_film_content_from_pressure() {
    SimulationConfig cfg = base_config();
    cfg.bulk_modulus = 1.0e9;
    cfg.bc_z_south_val = 1.5e6;
    cfg.bc_z_north_val = 1.5e6;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields = make_fields(mesh, cfg);
    LinearSystem sys(mesh);
    solve_elrod_once(fields, sys, mesh, comm, cfg);

    const Range theta_range = global_range(fields["theta"], mesh);
    const Range pressure_range = global_range(fields["pressure"], mesh);
    const double expected_theta = cfg.elrod_boundary_theta(cfg.bc_z_south_val);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  pressure-derived theta for pressure=1.5e6: expected="
                  << expected_theta << " range=["
                  << theta_range.min_value << ", " << theta_range.max_value << "]\n";
        std::cout << "  pressure range: [" << pressure_range.min_value
                  << ", " << pressure_range.max_value << "] Pa\n";
    }

    const double theta_tol = 1e-6;
    const double pressure_tol = cfg.bulk_modulus * 1e-6;
    check(std::abs(theta_range.min_value - expected_theta) < theta_tol,
          "south/north theta min should be derived from pressure and bulk modulus");
    check(std::abs(theta_range.max_value - expected_theta) < theta_tol,
          "south/north theta max should be derived from pressure and bulk modulus");
    check(std::abs(pressure_range.max_value - cfg.bc_z_south_val) < pressure_tol,
          "pressure should recover from pressure-derived theta");
}

static void test_legacy_boundary_film_content_is_ignored() {
    SimulationConfig cfg = base_config();
    cfg.bulk_modulus = 1.0e9;
    cfg.bc_z_south_val = 1.5e6;
    cfg.bc_z_north_val = 1.5e6;
    cfg.bc_z_south_theta = 1.2;  // deprecated input: must not drive the solve
    cfg.bc_z_north_theta = 1.2;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields = make_fields(mesh, cfg);
    LinearSystem sys(mesh);
    solve_elrod_once(fields, sys, mesh, comm, cfg);

    const Range theta_range = global_range(fields["theta"], mesh);
    const Range pressure_range = global_range(fields["pressure"], mesh);
    const double expected_theta = cfg.elrod_boundary_theta(cfg.bc_z_south_val);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  legacy theta ignored: expected pressure-derived theta="
                  << expected_theta << " range=["
                  << theta_range.min_value << ", " << theta_range.max_value << "]\n";
        std::cout << "  recovered pressure expected=" << cfg.bc_z_south_val
                  << " actual max=" << pressure_range.max_value << " Pa\n";
    }

    const double theta_tol = 1e-6;
    const double pressure_tol = cfg.bulk_modulus * 1e-6;
    check(std::abs(theta_range.min_value - expected_theta) < theta_tol,
          "legacy south theta input should not override pressure-derived theta");
    check(std::abs(theta_range.max_value - expected_theta) < theta_tol,
          "legacy north theta input should not override pressure-derived theta");
    check(std::abs(pressure_range.max_value - cfg.bc_z_south_val) < pressure_tol,
          "recovered pressure should follow bc_z_*_val, not legacy theta BC");
}

static void test_pressure_supply_with_realistic_bulk_modulus() {
    SimulationConfig cfg = base_config();
    cfg.n_theta_global = 48;
    cfg.n_z_global = 16;
    cfg.bulk_modulus = 1e9;
    cfg.bc_z_south_val = 1.5e6;
    cfg.bc_z_north_val = 1.5e6;
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0;
    groove.size = 20.0;
    groove.p_supply = 1.5e6;
    cfg.inlets.push_back(groove);

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields = make_fields(mesh, cfg);
    LinearSystem sys(mesh);
    solve_elrod_once(fields, sys, mesh, comm, cfg);

    const Range theta_range = global_range(fields["theta"], mesh);
    const double supply_theta = EOS::theta_from_pressure(groove.p_supply, cfg.p_cav, cfg.bulk_modulus);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  realistic-beta supply theta expected=" << supply_theta
                  << " field range=[" << theta_range.min_value << ", "
                  << theta_range.max_value << "]\n";
    }

    check(theta_range.min_value >= cfg.elrod_boundary_theta(cfg.bc_z_south_val) - 1e-6,
          "theta should not drop below pressure-derived outlet film content in this no-source case");
    check(theta_range.max_value > 1.0001, "pressure inlet did not raise theta");
    check(theta_range.max_value < 1.01, "realistic bulk modulus should keep theta near one for 1.5 MPa");
}

static void test_elrod_initial_pressure_uses_outlet_pressure() {
    SimulationConfig cfg = base_config();
    cfg.bulk_modulus = 1.0e9;
    cfg.p_cav = 2.0e4;
    cfg.bc_z_south_val = 1.5e6;
    cfg.bc_z_north_val = 1.0e5;
    cfg.bc_z_south_theta = 1.0;  // deprecated input should be ignored

    check(std::abs(cfg.initial_pressure() - cfg.bc_z_south_val) < 1.0e-6,
          "Elrod initial pressure should use outlet pressure, not legacy film-content BC");
}

static void test_explicit_pressure_keeps_elrod_capacity_term() {
    SimulationConfig cfg = base_config();
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.dt = 1e-13;
    cfg.omega = 0.0;
    cfg.bc_z_south_val = cfg.p_cav;
    cfg.bc_z_north_val = cfg.p_cav;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields = make_fields(mesh, cfg);
    fields.add("dh_dt", mesh).fill(0.0);
    fields["theta"].fill(0.5);
    fields["pressure"].fill(cfg.p_cav);

    LinearSystem sys(mesh);
    solve_elrod_once(fields, sys, mesh, comm, cfg);

    const Range theta_range = global_range(fields["theta"], mesh);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  explicit pressure capacity theta range: ["
                  << theta_range.min_value << ", " << theta_range.max_value << "]\n";
    }

    check(theta_range.min_value > 0.45, "transient Elrod capacity should retain old theta lower bound");
    check(theta_range.max_value < 0.75, "transient solve skipped Elrod capacity term");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_elrod_bc...\n";

    test_elrod_dirichlet_derives_film_content_from_pressure();
    test_legacy_boundary_film_content_is_ignored();
    test_pressure_supply_with_realistic_bulk_modulus();
    test_elrod_initial_pressure_uses_outlet_pressure();
    test_explicit_pressure_keeps_elrod_capacity_term();

    if (rank == 0) std::cout << "PASS: test_elrod_bc\n";
    PetscFinalize();
    return 0;
}
