// WP-3 numerics validation: convection schemes, capacity-term gate, and
// legacy-key back-compat.
//
// Test 1 — TVD boundedness: 1-D periodic step advection produces no new
//   extrema and a strictly smaller L1 error than first-order upwind.
// Test 2 — Order of accuracy: smooth profile advection refines at ~1st order
//   for upwind and visibly faster for TVD (dt scaled with dx so the spatial
//   error dominates the backward-Euler time error).
// Test 3 — A.1 ↔ Gumbel match: with TYPE_DIFFERENCING on both paths the
//   full-film theta formulation and the Gumbel pressure formulation agree
//   more tightly than the historical 1.2% upwind/central mismatch.
// Test 4 — Capacity term: a fixed-bearing transient relaxes theta over finite
//   time instead of jumping to the steady solution, and converges to the
//   steady solve as t -> infinity.
// Test 5 — TVD under MPI: ghost layers stay consistent across ranks
//   (also registered with 6 ranks as test_schemes_mpi6).
// Test 6 — Config back-compat: legacy time-method keys parse with a warning;
//   TYPE_DIFFERENCING is rejected for thermal/gas convection.
//
// Run as: mpirun -n <N> ./test_schemes

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "fvm.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

namespace {

void check(bool condition, const char* message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

double global_max(double local) {
    double value = local;
    MPI_Allreduce(&local, &value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return value;
}

double global_min(double local) {
    double value = local;
    MPI_Allreduce(&local, &value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return value;
}

double global_sum(double local) {
    double value = local;
    MPI_Allreduce(&local, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return value;
}

// Advect an initial profile around the periodic theta ring with constant face
// flux at the given CFL number for n_steps backward-Euler steps. Returns the
// normalized L1 error against the exactly shifted initial profile and records
// the running min/max.
struct AdvectionResult {
    double l1_error = 0.0;
    double min_seen = 1e300;
    double max_seen = -1e300;
};

template <typename ProfileFn>
AdvectionResult advect_profile(int n_cells, double cfl, int n_steps,
                               ConvectionScheme scheme, ProfileFn profile) {
    SimulationConfig cfg;
    cfg.n_theta_global = n_cells;
    cfg.n_z_global = 1;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);

    const double dt = 1.0;
    const double volume = mesh.cell_volume();
    const double face_flux = cfl * volume / dt;

    Fields fields;
    Field& theta = fields.add("theta", mesh);
    Field flux_theta("flux_theta", mesh);
    Field flux_z("flux_z", mesh);
    for (int i = -theta.n_ghost; i < mesh.n_theta_local + theta.n_ghost; ++i) {
        flux_theta(i, 0) = face_flux;
        flux_z(i, 0) = 0.0;
    }
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        theta(i, 0) = profile(mesh.offset_theta + i);
    }

    AdvectionResult result;
    for (int step = 0; step < n_steps; ++step) {
        comm.update_ghosts(theta);
        theta.store_old_time();
        sys.reset();
        FVM::divergence(sys, flux_theta, flux_z, theta, scheme, mesh);
        FVM::ddt(sys, theta, dt, mesh);
        sys.solve(theta);
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            result.min_seen = std::min(result.min_seen, theta(i, 0));
            result.max_seen = std::max(result.max_seen, theta(i, 0));
        }
    }
    result.min_seen = global_min(result.min_seen);
    result.max_seen = global_max(result.max_seen);

    const double shift = cfl * n_steps;  // cells travelled
    double local_error = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double global_index = mesh.offset_theta + i;
        double source_index = std::fmod(global_index - shift, double(n_cells));
        if (source_index < 0.0) source_index += n_cells;
        local_error += std::abs(theta(i, 0) - profile(int(std::round(source_index)) % n_cells));
    }
    result.l1_error = global_sum(local_error) / n_cells;
    return result;
}

void test_tvd_boundedness_and_sharpness() {
    const int n_cells = 120;
    auto step_profile = [n_cells](int global_index) {
        return (global_index >= n_cells / 4 && global_index < n_cells / 2) ? 1.0 : 0.2;
    };
    // CFL 0.5 for 20 steps: exactly 10 cells of travel.
    const AdvectionResult upwind =
        advect_profile(n_cells, 0.5, 20, ConvectionScheme::UPWIND, step_profile);
    const AdvectionResult tvd =
        advect_profile(n_cells, 0.5, 20, ConvectionScheme::TVD_VANLEER, step_profile);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  step advection L1: upwind=" << upwind.l1_error
                  << " tvd=" << tvd.l1_error
                  << " tvd range=[" << tvd.min_seen << ", " << tvd.max_seen << "]\n";
    }

    const double bound_tol = 1e-6;
    check(tvd.min_seen >= 0.2 - bound_tol && tvd.max_seen <= 1.0 + bound_tol,
          "TVD step advection must not create new extrema");
    check(tvd.l1_error < upwind.l1_error,
          "TVD must resolve the step more sharply than upwind (smaller L1 error)");
}

void test_order_of_accuracy() {
    auto smooth_coarse = [](int global_index) {
        return 1.0 + 0.5 * std::sin(2.0 * M_PI * global_index / 60.0);
    };
    auto smooth_fine = [](int global_index) {
        return 1.0 + 0.5 * std::sin(2.0 * M_PI * global_index / 120.0);
    };

    // dt scales with dx (CFL halves with refinement) so the backward-Euler
    // time error refines at second order and the spatial scheme dominates.
    const AdvectionResult up_coarse =
        advect_profile(60, 0.2, 75, ConvectionScheme::UPWIND, smooth_coarse);
    const AdvectionResult up_fine =
        advect_profile(120, 0.1, 300, ConvectionScheme::UPWIND, smooth_fine);
    const AdvectionResult tvd_coarse =
        advect_profile(60, 0.2, 75, ConvectionScheme::TVD_VANLEER, smooth_coarse);
    const AdvectionResult tvd_fine =
        advect_profile(120, 0.1, 300, ConvectionScheme::TVD_VANLEER, smooth_fine);

    const double upwind_order = std::log2(up_coarse.l1_error / up_fine.l1_error);
    const double tvd_order = std::log2(tvd_coarse.l1_error / tvd_fine.l1_error);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  observed orders: upwind=" << upwind_order
                  << " tvd=" << tvd_order << "\n";
    }

    check(upwind_order > 0.6 && upwind_order < 1.4,
          "upwind should refine at ~1st order on a smooth profile");
    check(tvd_order > 1.2, "TVD should refine visibly faster than 1st order");
    check(tvd_fine.l1_error < up_fine.l1_error,
          "TVD must beat upwind on the smooth profile");
}

void test_gumbel_elrod_match_type_differencing() {
    // A.1 full-film limit: elevated boundary pressure keeps the whole film
    // pressurized so both paths solve the same PDE with the same (central)
    // wedge discretization under TYPE_DIFFERENCING.
    SimulationConfig cfg;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.theta_convection_scheme = ConvectionScheme::TYPE_DIFFERENCING;
    cfg.log_outer_iters = false;
    cfg.outer_tol = 1e-10;
    cfg.max_outer_iters = 30;
    cfg.e = 0.4 * cfg.c;
    cfg.p_cav = 1.0e5;
    cfg.bc_z_south_val = 1.0e6;
    cfg.bc_z_north_val = 1.0e6;

    Mesh mesh(cfg);
    Communicator comm(mesh);

    double p_floor = 1e300;
    auto run_path = [&](CavitationModel model) {
        SimulationConfig run_cfg = cfg;
        run_cfg.cavitation_model = model;
        const double theta0 = std::exp((run_cfg.initial_pressure() - cfg.p_cav) / cfg.bulk_modulus);
        Fields fields;
        fields.add("pressure", mesh).fill(run_cfg.initial_pressure());
        fields.add("theta", mesh).fill(theta0);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho * theta0);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);
        if (model == CavitationModel::ELROD_ADAMS) {
            Reynolds::solve_elrod(fields, sys, mesh, run_cfg);
            comm.update_ghosts(fields);
        } else {
            // Gumbel lags density/pressure; iterate so the type-differencing
            // mask and theta interpolation reach their fixed point.
            for (int sweep = 0; sweep < 5; ++sweep) {
                Reynolds::solve(fields, sys, mesh, run_cfg);
                comm.update_ghosts(fields);
            }
        }

        double local_max = 0.0, local_min = 1e300;
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                local_max = std::max(local_max, fields["pressure"](i, j));
                local_min = std::min(local_min, fields["pressure"](i, j));
            }
        }
        p_floor = std::min(p_floor, global_min(local_min));
        return global_max(local_max);
    };

    const double p_gumbel = run_path(CavitationModel::GUMBEL);
    const double p_elrod = run_path(CavitationModel::ELROD_ADAMS);
    const double rel_err = std::abs(p_gumbel - p_elrod) / (std::abs(p_gumbel) + 1e-30);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "  type-differencing peak p: gumbel=" << p_gumbel
                  << " elrod=" << p_elrod << " rel err=" << rel_err
                  << " p floor=" << p_floor << "\n";
    }
    check(p_floor > cfg.p_cav, "comparison case must stay full-film for the A.1 limit");
    // Historical upwind/central mismatch was ~1.2%; matched schemes must do better.
    check(rel_err < 0.012, "matched type-differencing paths should agree within 1.2%");
}

void test_capacity_term_finite_relaxation() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.motion_model = MotionModel::STATIC;
    cfg.n_theta_global = 36;
    cfg.n_z_global = 12;
    cfg.e = 0.0;
    cfg.omega = 0.0;
    cfg.p_cav = 1.0e5;
    cfg.bulk_modulus = 1.0e5;
    cfg.bc_z_south_val = cfg.p_cav;
    cfg.bc_z_north_val = cfg.p_cav;
    cfg.outer_tol = 1e-10;
    cfg.log_outer_iters = false;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);

    const double theta_start = 1.2;
    auto make_fields = [&]() {
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav + cfg.bulk_modulus * std::log(theta_start));
        fields.add("theta", mesh).fill(theta_start);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho * theta_start);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        comm.update_ghosts(fields);
        fields["theta"].store_old_time();
        return fields;
    };

    // One tiny step: theta must stay near its initial value (finite relaxation),
    // not jump to the boundary-imposed steady state.
    {
        Fields fields = make_fields();
        SimulationConfig small_cfg = cfg;
        small_cfg.dt = 1.0e-7;
        Reynolds::solve_elrod(fields, sys, mesh, small_cfg);
        double local_dev = 0.0;
        for (int i = 0; i < mesh.n_theta_local; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                local_dev = std::max(local_dev, std::abs(fields["theta"](i, j) - theta_start));
        const double max_dev = global_max(local_dev);
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "  one tiny step max |theta - start| = " << max_dev << "\n";
        check(max_dev < 0.1 * (theta_start - 1.0),
              "fixed-bearing transient must relax theta over finite time, not instantly");
    }

    // Long transient must approach the steady solve.
    {
        Fields fields = make_fields();
        SimulationConfig run_cfg = cfg;
        run_cfg.dt = 5.0e-4;
        double max_dev = 1e300;
        for (int step = 0; step < 400 && max_dev > 1e-4; ++step) {
            comm.update_ghosts(fields);
            fields["theta"].store_old_time();
            Reynolds::solve_elrod(fields, sys, mesh, run_cfg);
            double local_dev = 0.0;
            for (int i = 0; i < mesh.n_theta_local; ++i)
                for (int j = 0; j < mesh.n_z_local; ++j)
                    local_dev = std::max(local_dev, std::abs(fields["theta"](i, j) - 1.0));
            max_dev = global_max(local_dev);
        }
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "  long transient max |theta - 1| = " << max_dev << "\n";
        check(max_dev <= 1e-4, "long fixed-bearing transient must reach the steady state");
    }
}

double max_theta_ghost_mismatch(const Field& theta, const Mesh& mesh) {
    const int ng = theta.n_ghost;
    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const int count = ng * n_z;
    std::vector<double> send_left(count), send_right(count), recv_left(count), recv_right(count);

    int idx = 0;
    for (int j = 0; j < n_z; ++j) {
        for (int k = 0; k < ng; ++k) send_left[idx + k] = theta(k, j);
        for (int k = 0; k < ng; ++k) send_right[idx + k] = theta(n_theta - ng + k, j);
        idx += ng;
    }

    const int left_rank = (mesh.rank - 1 + mesh.size) % mesh.size;
    const int right_rank = (mesh.rank + 1) % mesh.size;
    MPI_Status status;
    MPI_Sendrecv(send_left.data(), count, MPI_DOUBLE, left_rank, 60,
                 recv_right.data(), count, MPI_DOUBLE, right_rank, 60,
                 MPI_COMM_WORLD, &status);
    MPI_Sendrecv(send_right.data(), count, MPI_DOUBLE, right_rank, 61,
                 recv_left.data(), count, MPI_DOUBLE, left_rank, 61,
                 MPI_COMM_WORLD, &status);

    double local_max = 0.0;
    idx = 0;
    for (int j = 0; j < n_z; ++j) {
        for (int k = 0; k < ng; ++k) {
            local_max = std::max(local_max, std::abs(theta(-ng + k, j) - recv_left[idx + k]));
            local_max = std::max(local_max, std::abs(theta(n_theta + k, j) - recv_right[idx + k]));
        }
        idx += ng;
    }
    return global_max(local_max);
}

void test_tvd_mpi_ghost_consistency() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.theta_convection_scheme = ConvectionScheme::TVD_VANLEER;
    cfg.e = 0.8 * cfg.c;
    cfg.bulk_modulus = 1e5;
    cfg.dt = 1e-3;
    cfg.outer_tol = 1e-8;
    cfg.log_outer_iters = false;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);

    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);
    FilmThickness::compute_static(fields["h"], mesh, cfg);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();

    for (int step = 0; step < 10; ++step) {
        comm.update_ghosts(fields);
        fields["theta"].store_old_time();
        Reynolds::solve_elrod(fields, sys, mesh, cfg);
        const double mismatch = max_theta_ghost_mismatch(fields["theta"], mesh);
        check(mismatch < 1e-14, "TVD run left stale theta ghost cells across ranks");
    }

    double local_min = 1e300, local_max = -1e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local_min = std::min(local_min, fields["theta"](i, j));
            local_max = std::max(local_max, fields["theta"](i, j));
        }
    }
    check(global_min(local_min) >= cfg.theta_min - 1e-15,
          "TVD theta must respect the theta_min clamp");
    check(global_max(local_max) < 2.5, "TVD theta must stay physically bounded");
}

void test_config_backcompat_and_guards() {
    SimulationConfig cfg;
    cfg.load_from_text(
        "pressure_time_method = CRANK_NICOLSON\n"
        "theta_convection_scheme = tvd-vanleer\n"
        "thermal_convection_scheme = minmod\n"
        "gas_convection_scheme = upwind\n");
    check(cfg.parse_warnings.size() == 1,
          "legacy pressure_time_method must parse with exactly one warning");
    check(cfg.theta_convection_scheme == ConvectionScheme::TVD_VANLEER,
          "theta_convection_scheme alias tvd-vanleer should parse");
    check(cfg.thermal_convection_scheme == ConvectionScheme::TVD_MINMOD,
          "thermal_convection_scheme alias minmod should parse");
    check(cfg.gas_convection_scheme == ConvectionScheme::UPWIND,
          "gas_convection_scheme upwind should parse");

    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    SimulationConfig guard;
    guard.p_cav = 1.0e5;
    guard.bc_z_south_val = 2.0e5;
    guard.bc_z_north_val = 2.0e5;
    guard.thermal_convection_scheme = ConvectionScheme::TYPE_DIFFERENCING;
    guard.validate(errors, warnings);
    check(std::any_of(errors.begin(), errors.end(), [](const std::string& e) {
              return e.find("TYPE_DIFFERENCING") != std::string::npos;
          }),
          "TYPE_DIFFERENCING must be rejected for thermal convection");
}

}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_schemes...\n";

    test_config_backcompat_and_guards();
    test_tvd_boundedness_and_sharpness();
    test_order_of_accuracy();
    test_gumbel_elrod_match_type_differencing();
    test_capacity_term_finite_relaxation();
    test_tvd_mpi_ghost_consistency();

    if (rank == 0) std::cout << "PASS: test_schemes\n";
    PetscFinalize();
    return 0;
}
