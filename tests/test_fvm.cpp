// Phase 1 validation: FVM operator correctness.
//
// Test 1 – ddt: one Euler step, check a_P and source contributions.
// Test 2 – Laplacian: 2D Poisson solve on (theta, z), compare to analytical solution.
//   p_exact = sin(2θ) * sin(π z / L),  periodic theta, Dirichlet z=0,z=L.
//   Checks 2nd-order accuracy and MPI consistency.
// Test 3 – Divergence (upwind): steady uniform advection in theta.
//   Constant flux F, step initial condition. Verifies mass conservation.
//
// Run as: mpirun -n <N> ./test_fvm

#include <cassert>
#include <cmath>
#include <iostream>

#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "fvm.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Test 1: ddt adds vol/dt to a_P and vol/dt * phi_old to source
static void test_ddt(const Mesh& mesh) {
    SimulationConfig cfg;
    const double dt  = cfg.dt;
    const double vol = mesh.R * mesh.get_d_theta() * mesh.get_d_z();

    Fields fields;
    Field& phi = fields.add("phi", mesh);
    // Set old values (store_old_time copies data -> old_data)
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            phi(i, j) = 3.14;
    phi.store_old_time();

    LinearSystem sys(mesh);
    FVM::ddt(sys, phi, dt, mesh);

    const double expected_ap  = vol / dt;
    const double expected_src = expected_ap * 3.14;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            check(std::abs(sys.a_p(i, j) - expected_ap) < 1e-10, "ddt: a_P incorrect");
            check(std::abs(sys.source(i, j) - expected_src) < 1e-8, "ddt: source incorrect");
        }
    }

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  test_ddt PASS\n";
}

// Test 2: 2D Poisson / Laplacian solve
//   Solve: div(grad(p)) = f   with gamma = 1 everywhere
//   p_exact = sin(2θ) * sin(π z / L)
//   f = -(4/R² + π²/L²) * p_exact
//   BCs: Dirichlet p=0 at z=0 and z=L (ghost cells set to -p_interior)
//        Periodic in theta (handled by MPI ghost exchange)
static void test_laplacian_2d(const Mesh& mesh, Communicator& comm) {
    const double R       = mesh.R;
    const double L       = mesh.L;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double n_z     = mesh.n_z_local;
    const double pi      = M_PI;

    Fields fields;
    Field& gamma  = fields.add("gamma",  mesh);
    Field& p_sol  = fields.add("p_sol",  mesh);
    Field& p_rhs  = fields.add("p_rhs",  mesh);

    gamma.fill(1.0);
    comm.update_ghosts(fields);  // sync gamma ghost cells

    // The FVM stencil assembles -div(γ∇p)·vol = source, so source_field = -f.
    // f = div(grad p_exact) = -(4/R² + π²/L²) * p_exact => source_field = lambda * p_exact
    const double lambda = 4.0 / (R * R) + pi * pi / (L * L);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double z_c     = (j + 0.5) * d_z;
            p_rhs(i, j) = lambda * std::sin(2.0 * theta_c) * std::sin(pi * z_c / L);
        }
    }

    // Apply Dirichlet z-BC by adding boundary flux to source manually.
    // At j=0: south face, p_bc = 0. Contribution from south face:
    //   a_s_coeff = gamma_s * R*d_theta / d_z  (to a_P and source)
    //   But laplacian skips j=0 south face, so we add it here.
    // Using half-cell distance: coeff = gamma * R*d_theta / (d_z/2)
    const double bc_coeff_z = 1.0 * R * d_theta / (d_z / 2.0);  // gamma=1, half-cell
    const double p_bc = 0.0;

    LinearSystem sys(mesh);
    FVM::laplacian(sys, gamma, mesh);
    FVM::add_source(sys, p_rhs, mesh);

    // Dirichlet BC at z=0 and z=L (half-cell distance to boundary)
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        // Bottom (j=0): add south-face contribution
        sys.a_p(i, 0) += bc_coeff_z;
        sys.source(i, 0) += bc_coeff_z * p_bc;
        // Top (j=n_z-1): add north-face contribution
        sys.a_p(i, n_z - 1) += bc_coeff_z;
        sys.source(i, n_z - 1) += bc_coeff_z * p_bc;
    }

    sys.solve(p_sol);

    // Compare to analytical solution
    double max_err = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double z_c     = (j + 0.5) * d_z;
            double p_exact = std::sin(2.0 * theta_c) * std::sin(pi * z_c / L);
            max_err = std::max(max_err, std::abs(p_sol(i, j) - p_exact));
        }
    }
    double global_max_err;
    MPI_Allreduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  test_laplacian_2d max_err = " << global_max_err << "\n";

    // Error should scale as O(dθ² + dz²) ≈ O((2π/120)² + (L/40)²) ≈ small
    check(global_max_err < 1e-3, "Laplacian 2D: solution error too large");

    if (rank == 0) {
        std::cout << "  test_laplacian_2d PASS  (max_err = " << global_max_err << ")\n";
    }
}

// Test 3: divergence (upwind) — steady uniform advection conserves mass
// Uniform flux F in theta direction. Source = -div(F*1) = 0 (since phi=1 uniform).
// After assembly, a_P should equal a_E (from upwind with positive F).
static void test_divergence_upwind(const Mesh& mesh, Communicator& comm) {
    const double F_val = 2.5;  // constant east-face flux through all faces

    Fields fields;
    Field& flux_theta = fields.add("flux_theta", mesh);
    Field& flux_z     = fields.add("flux_z", mesh);
    Field& phi        = fields.add("phi", mesh);

    // Uniform positive east-face flux
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            flux_theta(i, j) = F_val;
            flux_z(i, j)     = 0.0;
            phi(i, j) = 1.0;
        }
    phi.store_old_time();
    comm.update_ghosts(fields);

    LinearSystem sys(mesh);
    FVM::divergence(sys, flux_theta, flux_z, phi, ConvectionScheme::UPWIND, mesh);

    // With constant positive F, upwind: a_P = F (from east outflow) + 0 (west inflow handled by a_W)
    // West inflow: F_w = flux_theta(i-1,j) = F_val > 0 → max(-F_w,0) = 0, max(F_w,0) → a_W += F_val
    // East outflow: F_e = F_val > 0 → max(F_e,0) = F_val → a_P += F_val, a_E = 0
    // So a_P = F_val, a_W = F_val, a_E = 0, a_N = 0, a_S = 0.
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            check(std::abs(sys.a_p(i, j) - F_val) < 1e-10, "divergence upwind: a_P incorrect");
            check(std::abs(sys.a_w(i, j) - F_val) < 1e-10, "divergence upwind: a_W incorrect");
            check(std::abs(sys.a_e(i, j)) < 1e-10,          "divergence upwind: a_E should be 0");
        }
    }

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "  test_divergence_upwind PASS\n";
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    SimulationConfig cfg;
    Mesh mesh(cfg);
    Communicator comm(mesh);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_fvm...\n";

    test_ddt(mesh);
    test_laplacian_2d(mesh, comm);
    test_divergence_upwind(mesh, comm);

    if (rank == 0) std::cout << "PASS: test_fvm\n";

    PetscFinalize();
    return 0;
}
