#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

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
    check(cfg.motion_time_method == TimeSteppingMethod::CRANK_NICOLSON,
          "semi implicit alias should map to CRANK_NICOLSON");
    // Field solves are backward-Euler; the legacy keys must parse with a warning.
    check(cfg.parse_warnings.size() == 2,
          "legacy pressure/temperature_time_method keys should each warn once");
    check(cfg.parse_warnings[0].find("ignored") != std::string::npos,
          "legacy time-method warning should say the key is ignored");
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


static std::filesystem::path write_rank_file(const std::string& stem, const std::string& text) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const auto path = std::filesystem::temp_directory_path() /
        ("pancake_" + stem + "_" + std::to_string(rank) + ".csv");
    std::ofstream(path) << text;
    return path;
}

static bool near(double a, double b, double tol) { return std::abs(a - b) <= tol; }

static void test_load_table_interpolation_and_wrap() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.load_model = LoadModel::TABLE;
    cfg.load_table_abscissa = LoadTableAbscissa::PHASE_DEG;
    cfg.load_table_format = LoadTableFormat::POLAR;
    cfg.load_frequency = 10.0;
    // Header, comment, a jump between 179 and 180 deg, and no 360 deg row.
    const auto path = write_rank_file("load_polar",
        "phase_deg,F_N,angle_deg\n# comment\n0,100,90\n179,100,90\n180,300,90\n350,200,90\n");
    std::string error;
    cfg.load_table = JournalMotion::read_load_table(path, cfg, error);
    check(error.empty() && cfg.load_table.size() == 4, "polar load table must parse with header and comment");

    auto load_at_phase = [&](double phase) {
        return JournalMotion::external_load(cfg, phase / (360.0 * cfg.load_frequency));
    };
    check(near(load_at_phase(89.5).y, 100.0, 1e-9) && near(load_at_phase(89.5).x, 0.0, 1e-9),
          "constant segment must interpolate to its value");
    check(near(load_at_phase(179.5).y, 200.0, 1e-9), "jump segment interpolates linearly over its 1 deg");
    check(near(load_at_phase(180.0).y, 300.0, 1e-9), "jump row must be reproduced without smoothing");
    check(near(load_at_phase(355.0).y, 150.0, 1e-9), "periodic wrap joins the last row to the first + 360");
    check(near(load_at_phase(720.0 + 355.0).y, 150.0, 1e-9), "phase must wrap over later cycles");
    check(near(JournalMotion::load_phase_deg(cfg, 1.25 / cfg.load_frequency), 90.0, 1e-9),
          "load phase must advance at load_frequency and wrap at 360");

    cfg.load_phase_offset_deg = 90.0;
    check(near(JournalMotion::external_load(cfg, 0.0).y, 100.0, 1e-9), "phase offset shifts the cycle");
    cfg.load_phase_offset_deg = 0.0;

    const auto bad = write_rank_file("load_bad", "0,1,0\n360,2,0\n");
    std::string bad_error;
    check(JournalMotion::read_load_table(bad, cfg, bad_error).empty() && !bad_error.empty(),
          "PHASE_DEG table with a 360 deg row must be rejected");
    const auto unsorted = write_rank_file("load_unsorted", "0,1,0\n20,2,0\n10,3,0\n");
    std::string unsorted_error;
    check(JournalMotion::read_load_table(unsorted, cfg, unsorted_error).empty() && !unsorted_error.empty(),
          "non-increasing abscissa must be rejected");
    std::filesystem::remove(path);
    std::filesystem::remove(bad);
    std::filesystem::remove(unsorted);
}

static void test_load_direction_across_360() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.load_model = LoadModel::TABLE;
    cfg.load_frequency = 1.0;
    // Direction 350 deg -> 10 deg: the midpoint must point along +x, not -x.
    const auto path = write_rank_file("load_wrap", "0,50,350\n10,50,10\n");
    std::string error;
    cfg.load_table = JournalMotion::read_load_table(path, cfg, error);
    check(error.empty(), "wrap table must parse");
    const auto mid = JournalMotion::external_load(cfg, 5.0 / 360.0);
    check(mid.x > 49.0 && near(mid.y, 0.0, 1e-9), "direction interpolation must not cross through 180 deg");
    cfg.load_table_format = LoadTableFormat::CARTESIAN;
    cfg.load_table_scale = 2.0;
    cfg.load_table_angle_offset_deg = 90.0;
    const auto cart = write_rank_file("load_cart", "0,1,0\n");
    cfg.load_table = JournalMotion::read_load_table(cart, cfg, error);
    const auto rotated = JournalMotion::external_load(cfg, 0.123);
    check(near(rotated.x, 0.0, 1e-12) && near(rotated.y, 2.0, 1e-12),
          "Cartesian row must be scaled and rotated by the angle offset");
    std::filesystem::remove(path);
    std::filesystem::remove(cart);
}

static void test_time_table_and_constant_equivalence() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.load_model = LoadModel::TABLE;
    cfg.load_table_abscissa = LoadTableAbscissa::TIME_S;
    cfg.load_table_format = LoadTableFormat::CARTESIAN;
    const auto path = write_rank_file("load_time", "0.0,0,0\n1.0,10,-20\n");
    std::string error;
    cfg.load_table = JournalMotion::read_load_table(path, cfg, error);
    const auto half = JournalMotion::external_load(cfg, 0.5);
    check(near(half.x, 5.0, 1e-12) && near(half.y, -10.0, 1e-12), "TIME_S table must interpolate in time");
    bool threw = false;
    try { (void)JournalMotion::external_load(cfg, 1.5); } catch (const std::runtime_error&) { threw = true; }
    check(threw, "TIME_S table must refuse times outside its range");

    // One-row table, zero-speed rotating load and the constant keys agree.
    SimulationConfig constant;
    constant.external_load_x = 120.0;
    constant.external_load_y = -50.0;
    constant.external_load_z = 3.0;
    SimulationConfig one_row = constant;
    one_row.motion_model = MotionModel::MOVING_BEARING;
    one_row.load_model = LoadModel::TABLE;
    one_row.load_table_format = LoadTableFormat::CARTESIAN;
    one_row.load_frequency = 1.0;
    one_row.load_table = {{0.0, 120.0, -50.0}};
    SimulationConfig still = constant;
    still.load_model = LoadModel::ROTATING;
    still.load_magnitude = std::hypot(120.0, -50.0);
    still.load_rotation_start_deg = std::atan2(-50.0, 120.0) * 180.0 / M_PI;
    still.load_rotation_speed = 0.0;
    for (double t : {0.0, 0.37, 2.5}) {
        const auto a = JournalMotion::external_load(constant, t);
        const auto b = JournalMotion::external_load(one_row, t);
        const auto c = JournalMotion::external_load(still, t);
        check(a.x == b.x && a.y == b.y && a.z == b.z, "one-row table must equal external_load_x|y|z");
        check(near(a.x, c.x, 1e-12) && near(a.y, c.y, 1e-12) && a.z == c.z,
              "zero-speed rotating load must equal external_load_x|y|z");
    }
    std::filesystem::remove(path);
}

static void test_rotating_load_direction() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.load_model = LoadModel::ROTATING;
    cfg.load_magnitude = 10.0;
    cfg.load_rotation_start_deg = 30.0;
    cfg.load_rotation_speed = M_PI;
    const auto a = JournalMotion::external_load(cfg, 0.5);
    check(near(std::hypot(a.x, a.y), 10.0, 1e-12), "rotating load keeps its magnitude");
    check(near(std::atan2(a.y, a.x) * 180.0 / M_PI, 120.0, 1e-9), "rotating load turns at load_rotation_speed");
}

// Linear stand-in for the film force: F = -k_f x - c_f v (v = backward-Euler velocity).
static JournalMotion::Vector3 linear_film(const JournalMotion::BearingState& s, double k, double c) {
    return {-k * s.x - c * s.vx, -k * s.y - c * s.vy, 0.0};
}

static void test_implicit_coupling_linear_film() {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.motion_coupling = MotionCoupling::IMPLICIT;
    cfg.c = 1.0;
    cfg.motion_coupling_tolerance = 1e-9;  // m/dt^2 x carries ~1e-10 relative rounding
    cfg.motion_coupling_position_tolerance = 1e-13;
    cfg.motion_coupling_fd_step = 1e-6;
    cfg.motion_coupling_step_limit = 1.0;
    cfg.bearing_mass = 1.0;
    cfg.external_load_x = 100.0;
    cfg.external_load_y = -40.0;
    const double k = 100.0, c = 4.0;

    // (a) Each step must land on the closed-form backward-Euler update.
    JournalMotion::ImplicitCoupling coupling(cfg);
    JournalMotion::BearingState state;
    const double dt = 1e-4;
    double t = 0.0;
    int max_solves = 0;
    for (int n = 0; n < 10000; ++n) {
        const JournalMotion::BearingState start = state;
        t += dt;
        JournalMotion::BearingState trial = coupling.begin_step(state, dt, t);
        while (!coupling.update(trial, linear_film(trial, k, c))) {}
        max_solves = std::max(max_solves, coupling.evaluations());
        const double a = cfg.bearing_mass / (dt * dt) + c / dt + k;
        const double x_be = (cfg.external_load_x + cfg.bearing_mass * (start.x / (dt * dt) + start.vx / dt)
                             + c * start.x / dt) / a;
        check(near(trial.x, x_be, 1e-9 * std::abs(x_be) + 1e-14), "implicit step must equal the backward-Euler solution");
        state = trial;
    }
    check(max_solves <= 4, "a linear film force must converge within the probe solves plus one Newton step");
    // (b) Against the analytic underdamped response to a step load (O(dt) scheme).
    const double wn = std::sqrt(k / cfg.bearing_mass), zeta = c / (2.0 * std::sqrt(k * cfg.bearing_mass));
    const double wd = wn * std::sqrt(1.0 - zeta * zeta);
    const double x_exact = cfg.external_load_x / k *
        (1.0 - std::exp(-zeta * wn * t) * (std::cos(wd * t) + zeta / std::sqrt(1.0 - zeta * zeta) * std::sin(wd * t)));
    check(near(state.x, x_exact, 0.01 * cfg.external_load_x / k), "implicit coupling must follow the analytic response");
    check(near(state.y / state.x, cfg.external_load_y / cfg.external_load_x, 1e-9), "axes must stay independent");
}

static void test_implicit_coupling_stiff_film_is_stable() {
    // Film-like squeeze damping: the explicit (lagged-force) limit is dt < 2 m / c = 2e-6 s.
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::MOVING_BEARING;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.c = 1e-4;
    cfg.bearing_mass = 1e-3;
    cfg.external_load_x = 1.0;
    cfg.motion_coupling_tolerance = 1e-10;
    cfg.motion_coupling_position_tolerance = 1e-9;
    const double k = 1e6, c = 1e3, dt = 1e-4;

    JournalMotion::ImplicitCoupling coupling(cfg);
    JournalMotion::BearingState implicit_state, explicit_state;
    JournalMotion::Vector3 lagged{};
    cfg.motion_time_method = TimeSteppingMethod::EULER_IMPLICIT;
    double t = 0.0;
    for (int n = 0; n < 200; ++n) {
        t += dt;
        JournalMotion::BearingState trial = coupling.begin_step(implicit_state, dt, t);
        while (!coupling.update(trial, linear_film(trial, k, c))) {}
        implicit_state = trial;
        JournalMotion::advance(explicit_state, lagged, cfg, dt, t);
        lagged = linear_film(explicit_state, k, c);
    }
    const double static_x = cfg.external_load_x / k;
    check(near(implicit_state.x, static_x, 1e-6 * static_x), "implicit coupling must settle at W / k");
    check(!std::isfinite(explicit_state.x) || std::abs(explicit_state.x) > 1e3 * static_x,
          "the lagged explicit coupling is expected to diverge at this dt (the reason for IMPLICIT)");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_moving_bearing_matches_static_attitude();
    test_config_aliases();
    test_motion_integrators_constant_load();
    test_load_table_interpolation_and_wrap();
    test_load_direction_across_360();
    test_time_table_and_constant_equivalence();
    test_rotating_load_direction();
    test_implicit_coupling_linear_film();
    test_implicit_coupling_stiff_film_is_stable();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_journal_motion\n";

    PetscFinalize();
    return 0;
}
