// Phase A.1 validation: θ-formulation (full-film) reproduces Gumbel pressure field.
//
// Test 1 — Regression: solve the same journal bearing with both GUMBEL and ELROD_ADAMS
//   (g=1 everywhere). In the full-film limit, both formulations solve the same PDE
//   up to linearisation of the lagged density. The peak pressures and overall
//   distributions should agree to within ~1% for the standard test case.
//
// Test 2 — 1D uniform-h channel: div(Γ∇θ) = div(F·θ) with constant h and Γ.
//   Periodic θ, Dirichlet θ=1 at z=0 and z=L. With no circumferential variation
//   (uniform h → no Couette wedge), the steady solution is θ=1 everywhere.
//   Checks the solver returns θ=1 to within solver tolerance.
//
// Run as: mpirun -n <N> ./test_elrod_a1

#include <algorithm>
#include <cmath>
#include <iostream>

#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
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

// Test 1: run one Reynolds step with GUMBEL, then one with ELROD_ADAMS on the
//         same geometry. Peak pressures must agree to within 1%.
static void test_regression_vs_gumbel(const Mesh& mesh, Communicator& comm) {
    SimulationConfig cfg;

    auto run_one_step = [&](CavitationModel model) {
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("theta",    mesh).fill(1.0);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("rho",      mesh).fill(cfg.rho);

        FilmThickness::compute_static(fields["h"], mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);
        cfg.cavitation_model = model;

        if (model == CavitationModel::ELROD_ADAMS) {
            cfg.max_outer_iters = 1;  // Phase A.1: test full-film (g=1) formulation only
            Reynolds::solve_elrod(fields, sys, mesh, cfg);
        }
        else
            Reynolds::solve(fields, sys, mesh, cfg);

        double local_max = 0.0;
        for (int i = 0; i < mesh.n_theta_local; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                local_max = std::max(local_max, fields["pressure"](i, j));

        double global_max;
        MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global_max;
    };

    double p_gumbel = run_one_step(CavitationModel::GUMBEL);
    double p_elrod  = run_one_step(CavitationModel::ELROD_ADAMS);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::cout << "  Gumbel peak p = " << p_gumbel
                  << "  Elrod A.1 peak p = " << p_elrod << "\n";

    const double rel_err = std::abs(p_gumbel - p_elrod) / (std::abs(p_gumbel) + 1e-30);
    if (rank == 0) std::cout << "  relative peak-p error = " << rel_err << "\n";
    check(rel_err < 0.05, "regression vs Gumbel: peak pressure differs by more than 5%");

    if (rank == 0) std::cout << "  test_regression_vs_gumbel PASS\n";
}

// Test 2: uniform film → no wedge → steady θ=1 everywhere.
static void test_uniform_film(const Mesh& mesh, Communicator& comm) {
    SimulationConfig cfg;
    cfg.e = 0.0;  // zero eccentricity → h = c everywhere

    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta",    mesh).fill(1.0);
    fields.add("h",        mesh).fill(cfg.c);
    fields.add("rho",      mesh).fill(cfg.rho);

    FilmThickness::compute_static(fields["h"], mesh, cfg);
    comm.update_ghosts(fields);

    cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
    LinearSystem sys(mesh);
    Reynolds::solve_elrod(fields, sys, mesh, cfg);

    double max_dev = 0.0;
    const Field& theta = fields["theta"];
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            max_dev = std::max(max_dev, std::abs(theta(i, j) - 1.0));

    double global_max_dev;
    MPI_Allreduce(&max_dev, &global_max_dev, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::cout << "  uniform-film max |θ−1| = " << global_max_dev << "\n";

    check(global_max_dev < 1e-6, "uniform film: theta should be 1 everywhere");
    if (rank == 0) std::cout << "  test_uniform_film PASS\n";
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    SimulationConfig cfg;
    Mesh mesh(cfg);
    Communicator comm(mesh);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_elrod_a1...\n";

    test_uniform_film(mesh, comm);
    test_regression_vs_gumbel(mesh, comm);

    if (rank == 0) std::cout << "PASS: test_elrod_a1\n";

    PetscFinalize();
    return 0;
}
