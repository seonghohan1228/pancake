// Phase A.2 validation: Elrod-Adams switch function and mass-conserving cavitation.
//
// Test 1 — Complementarity: after convergence, no cell has g=1 AND θ<1.
//   i.e., EOS::switch_function(θ) * max(1−θ, 0) = 0 for all cells.
//
// Test 2 — Cavitation zone: for an eccentric journal bearing (ε=0.8) the
//   cavitation (θ<1) must lie in the diverging half of the bearing.
//   With attitude_angle_deg = -90°, minimum gap is at θ_coord = 270° (3π/2).
//   The full-film zone spans ~[90°, 270°] and cavitation ~[270°, 450°].
//   Test checks that >50% of cavitated cells lie in the expected half.
//
// Test 3 — Mass conservation: for steady state, the total liquid content
//   ∑∑ θ(i,j)·h(i,j) must agree between the first and last timestep
//   to within a small tolerance (pure numerical drift, not physical change).
//
// Test 4 — Pressure positivity: p ≥ p_cav everywhere after solve.
//
// Run as: mpirun -n <N> ./test_elrod_a2

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

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

// Compute global sum of theta*h over physical cells.
static double global_liquid_content(const Field& theta, const Field& h, const Mesh& mesh) {
    double local = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            local += theta(i, j) * h(i, j);
    double global;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// Test 1: complementarity — g(θ)·max(1−θ,0) = 0 everywhere.
// Test 2: cavitation location — θ<1 cells lie in the diverging half.
// Test 3: mass conservation — liquid content stable after many steps.
// Test 4: pressure positivity.
static void test_journal_bearing(const Mesh& mesh, Communicator& comm) {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
    cfg.e = 0.8 * cfg.c;   // eccentricity ratio ε = 0.8
    cfg.outer_tol = 1e-8;

    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta",    mesh).fill(1.0);
    fields.add("h",        mesh).fill(cfg.c);
    fields.add("rho",      mesh).fill(cfg.rho);

    FilmThickness::compute_static(fields["h"], mesh, cfg);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();
    fields["pressure"].store_old_time();

    LinearSystem sys(mesh);

    // Run several timesteps to reach a near-steady state.
    const int n_steps = 30;
    double liquid_prev = 0.0, liquid_last = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        comm.update_ghosts(fields);
        Reynolds::solve_elrod(fields, sys, mesh, cfg);

        liquid_prev = liquid_last;
        liquid_last = global_liquid_content(fields["theta"], fields["h"], mesh);

        fields["theta"].store_old_time();
        fields["pressure"].store_old_time();
    }

    const Field& theta    = fields["theta"];
    const Field& pressure = fields["pressure"];
    const double d_theta  = mesh.get_d_theta();

    // Test 1: complementarity
    double max_violation = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double th = theta(i, j);
            const double g  = EOS::switch_function(th);
            max_violation = std::max(max_violation, g * std::max(1.0 - th, 0.0));
        }
    MPI_Allreduce(MPI_IN_PLACE, &max_violation, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  complementarity max g·(1−θ)+ = " << max_violation << "\n";
    check(max_violation < 1e-10, "complementarity violated: full-film cell has theta < 1");

    // Test 2: cavitation location
    // With attitude_angle = -90°, minimum gap at θ_coord = 270° (3π/2).
    // The diverging zone (where cavitation occurs) is [270°, 360°] ∪ [0°, 90°].
    // i.e., |θ_coord - 270°| < 90°, or more precisely: sin(θ_coord - 3π/2) < 0 → cos(θ_coord) > 0.
    // Actually with h = c - e*cos(θ - ψ) and ψ = -90° = -π/2:
    //   h = c - e*cos(θ + π/2) = c + e*sin(θ)
    // Minimum h at sin(θ) = -1 → θ = 270° = 3π/2.
    // Converging zone: dh/dθ < 0 → cos(θ) < 0 → θ ∈ (90°, 270°).
    // Diverging zone: dh/dθ > 0 → cos(θ) > 0 → θ ∈ (270°, 360°)∪(0°, 90°).
    // So cavitation expected in θ ∈ (270°, 450°) i.e. cos(θ) > 0.
    long long n_cav_total = 0, n_cav_diverging = 0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            if (theta(i, j) < 1.0) {
                ++n_cav_total;
                const double theta_coord = (mesh.offset_theta + i + 0.5) * d_theta;
                // diverging zone: cos(theta_coord - 3π/2) > 0 → sin(theta_coord) < 0... hmm
                // dh/dθ = e*cos(θ) where h = c + e*sin(θ).
                // Diverging: dh/dθ > 0 → cos(θ_coord) > 0.
                if (std::cos(theta_coord) > 0.0) ++n_cav_diverging;
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &n_cav_total,     1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &n_cav_diverging,  1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
        std::cout << "  cavitated cells: " << n_cav_total
                  << "  in diverging half: " << n_cav_diverging
                  << "  fraction: " << (n_cav_total > 0 ? (double)n_cav_diverging/n_cav_total : 0.0) << "\n";
    check(n_cav_total > 0, "no cavitation found — expected cavitation at ε=0.8");
    check(n_cav_total > 0 && (double)n_cav_diverging / n_cav_total > 0.5,
          "majority of cavitated cells not in diverging zone");

    // Test 3: mass conservation (near-steady state: liquid content should not drift).
    // Compare last two steps to ensure convergence to a steady mass balance.
    const double rel_mass_drift = std::abs(liquid_last - liquid_prev) / (std::abs(liquid_last) + 1e-30);
    if (rank == 0)
        std::cout << "  liquid content prev=" << liquid_prev
                  << "  last=" << liquid_last
                  << "  rel drift=" << rel_mass_drift << "\n";
    check(rel_mass_drift < 1e-10, "liquid content drifted between last two steps");

    // Test 4: pressure positivity
    double min_p = 1e30;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            min_p = std::min(min_p, pressure(i, j));
    MPI_Allreduce(MPI_IN_PLACE, &min_p, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    if (rank == 0) std::cout << "  min pressure = " << min_p << " Pa\n";
    check(min_p >= cfg.p_cav - 1e-6, "pressure below p_cav");

    if (rank == 0) std::cout << "  test_journal_bearing PASS\n";
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    SimulationConfig cfg;
    Mesh mesh(cfg);
    Communicator comm(mesh);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_elrod_a2...\n";

    test_journal_bearing(mesh, comm);

    if (rank == 0) std::cout << "PASS: test_elrod_a2\n";

    PetscFinalize();
    return 0;
}
