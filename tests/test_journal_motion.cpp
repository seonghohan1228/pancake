#include <algorithm>
#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "journal_motion.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static double global_max_abs_diff(const Field& a, const Field& b, const Mesh& mesh) {
    double local = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local = std::max(local, std::abs(a(i, j) - b(i, j)));
        }
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

static void test_moving_bearing_matches_static_attitude() {
    SimulationConfig cfg;
    cfg.e = 2.5e-5;
    cfg.attitude_angle_deg = 35.0;
    cfg.tilt_slope_x = 2.0e-4;
    cfg.tilt_slope_y = -1.0e-4;

    Mesh mesh(cfg);
    Field h_static("h_static", mesh);
    Field h_moving("h_moving", mesh);
    Field dh_dt("dh_dt", mesh);

    FilmThickness::compute_static(h_static, mesh, cfg);
    const JournalMotion::BearingState state = JournalMotion::initial_state(cfg);
    FilmThickness::compute_moving_bearing(h_moving, dh_dt, mesh, cfg, state);

    const double h_err = global_max_abs_diff(h_static, h_moving, mesh);
    check(h_err < 1e-15, "moving-bearing initial state must reproduce static h");

    double max_dhdt = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            max_dhdt = std::max(max_dhdt, std::abs(dh_dt(i, j)));
    MPI_Allreduce(MPI_IN_PLACE, &max_dhdt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(max_dhdt < 1e-15, "initial zero bearing velocity must give zero dh_dt");

    const double attitude = JournalMotion::equivalent_attitude_angle_deg(state, cfg);
    check(std::abs(attitude - cfg.attitude_angle_deg) < 1e-12,
          "bearing displacement must recover equivalent attitude angle");
}

static void test_config_aliases() {
    SimulationConfig cfg;
    cfg.load_from_text(
        "motion_model = moving-bearing\n"
        "pressure_time_method = imex\n"
        "motion_time_method = semi implicit\n"
        "temperature_time_method = crank-nicholson\n");

    check(cfg.motion_model == MotionModel::MOVING_BEARING, "motion_model alias parse failed");
    check(cfg.pressure_time_method == TimeSteppingMethod::CRANK_NICOLSON,
          "IMEX alias should map to CRANK_NICOLSON");
    check(cfg.motion_time_method == TimeSteppingMethod::CRANK_NICOLSON,
          "semi implicit alias should map to CRANK_NICOLSON");
    check(cfg.temperature_time_method == TimeSteppingMethod::CRANK_NICOLSON,
          "crank-nicholson alias should map to CRANK_NICOLSON");
}

static void test_motion_integrators_constant_load() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.bearing_initial_from_attitude = false;
    cfg.bearing_mass = 2.0;
    cfg.external_load_x = 4.0;

    JournalMotion::BearingState state;
    state.vx = 1.0;

    cfg.motion_time_method = TimeSteppingMethod::EULER_EXPLICIT;
    JournalMotion::advance(state, {}, cfg, 0.5);
    check(std::abs(state.x - 0.5) < 1e-14, "explicit Euler position mismatch");
    check(std::abs(state.vx - 2.0) < 1e-14, "explicit Euler velocity mismatch");

    state = {};
    cfg.motion_time_method = TimeSteppingMethod::RK4;
    JournalMotion::advance(state, {}, cfg, 0.5);
    check(std::abs(state.x - 0.25) < 1e-14, "RK4 position mismatch for constant acceleration");
    check(std::abs(state.vx - 1.0) < 1e-14, "RK4 velocity mismatch for constant acceleration");

    state = {};
    cfg.motion_time_method = TimeSteppingMethod::CRANK_NICOLSON;
    JournalMotion::advance(state, {}, cfg, 0.5);
    check(std::abs(state.x - 0.25) < 1e-14, "Crank-Nicolson position mismatch for constant acceleration");
    check(std::abs(state.vx - 1.0) < 1e-14, "Crank-Nicolson velocity mismatch for constant acceleration");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_moving_bearing_matches_static_attitude();
    test_config_aliases();
    test_motion_integrators_constant_load();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_journal_motion\n";

    PetscFinalize();
    return 0;
}
