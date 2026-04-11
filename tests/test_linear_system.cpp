// Phase 0 validation: parallel sparse linear system solve.
//
// Test 1 – Diagonal system: a_P = 1, source = val → solution = val.
//   Verifies basic assembly and solve on all ranks.
//
// Test 2 – 1D Helmholtz in theta (n_z = 1):
//   d²p/dθ² - alpha² * p = f   (periodic in theta, theta ∈ [0, 2π])
//   p_exact = sin(2θ),  f = -(4 + alpha²) * sin(2θ) * cell_vol
//   where theta_coeff = 1/(R² * dθ²) and alpha² dampens the null space.
//
// Run as: mpirun -n <N> ./test_linear_system

#include <cassert>
#include <cmath>
#include <iostream>

#include <petsc.h>

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Test 1: trivial diagonal system
static void test_diagonal(const Mesh& mesh) {
    LinearSystem sys(mesh);
    Fields fields;
    Field& sol = fields.add("sol", mesh);

    const double val = 42.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            sys.a_p(i, j)    = 1.0;
            sys.source(i, j) = val;
        }

    sys.solve(sol);

    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            check(std::abs(sol(i, j) - val) < 1e-7, "diagonal system: sol != val");

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  test_diagonal PASS\n";
}

// Test 2: 1D discrete Helmholtz in theta on a UNIT CIRCLE (R=1) for good conditioning.
// Source uses the exact discrete eigenvalue so p = sin(2θ_i) is the exact discrete solution.
// Any error is purely from the linear solver, not discretization.
static void test_helmholtz_1d(const Mesh& mesh) {
    // Override R to 1.0 for this test so lap_coeff ~ O(1) → condition number ~ O(n²)
    const double R_test    = 1.0;
    const double alpha2    = 1.0;
    const double d_theta   = mesh.get_d_theta();
    const double lap_coeff = 1.0 / (R_test * R_test * d_theta * d_theta);

    // Exact discrete eigenvalue for mode k=2: λ = -4*lap_coeff*sin²(d_theta)
    const double disc_eigval  = 4.0 * lap_coeff * std::sin(d_theta) * std::sin(d_theta);
    const double source_coeff = disc_eigval + alpha2;

    LinearSystem sys(mesh);
    Fields fields;
    Field& sol = fields.add("sol", mesh);

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double theta_i = (mesh.offset_theta + i + 0.5) * d_theta;
        const double p_exact = std::sin(2.0 * theta_i);
        sys.a_p(i, 0)    = 2.0 * lap_coeff + alpha2;
        sys.a_e(i, 0)    = lap_coeff;
        sys.a_w(i, 0)    = lap_coeff;
        sys.source(i, 0) = source_coeff * p_exact;
    }

    sys.solve(sol, 1e-12);

    double max_err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double theta_i = (mesh.offset_theta + i + 0.5) * d_theta;
        const double p_exact = std::sin(2.0 * theta_i);
        max_err = std::max(max_err, std::abs(sol(i, 0) - p_exact));
    }

    double global_max_err;
    MPI_Allreduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  helmholtz_1d max_err = " << global_max_err << "\n";
    check(global_max_err < 1e-8, "1D Helmholtz: max error too large");
    if (rank == 0) std::cout << "  test_helmholtz_1d PASS\n";
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // Use n_z = 1 for the 1D Helmholtz test
    SimulationConfig cfg;
    cfg.n_z_global = 1;

    Mesh mesh(cfg);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_linear_system...\n";

    test_diagonal(mesh);
    test_helmholtz_1d(mesh);

    if (rank == 0) std::cout << "PASS: test_linear_system\n";

    PetscFinalize();
    return 0;
}
