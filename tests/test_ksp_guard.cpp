// Regression test for the KSP convergence guard in LinearSystem::solve.
//
// Before the guard existed, a diverged or NaN-producing linear solve was copied
// straight into the field and the run went on to report success. This test feeds
// the solver a right-hand side containing a NaN and requires that the guard says
// so and aborts.
//
// The test is registered with PASS_REGULAR_EXPRESSION on the guard's own message
// rather than WILL_FAIL: any nonzero exit would satisfy WILL_FAIL, including a
// PETSc-internal error that never reached our check. Matching the message is what
// makes this a test of the guard rather than of PETSc.
//
// It therefore aborts on success. That is intentional; do not "fix" the abort.

#include <cmath>
#include <iostream>
#include <limits>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_ksp_guard...\n";

    SimulationConfig cfg;
    cfg.n_theta_global = 16;
    cfg.n_z_global = 4;
    Mesh mesh(cfg);
    LinearSystem sys(mesh);
    Field solution("solution", mesh);

    // Well-posed diagonal system except for a single poisoned source term, so the
    // failure comes from the solve rather than from a singular preconditioner.
    sys.reset();
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            sys.a_p(i, j) = 1.0;
            sys.source(i, j) = 1.0;
        }
    }
    if (mesh.offset_theta == 0) {
        sys.source(0, 0) = std::numeric_limits<double>::quiet_NaN();
    }

    // Expected to abort inside solve(). Reaching the next line means the guard
    // let a non-finite solution through.
    sys.solve(solution, 1e-8, 1000, "ksp guard probe");

    if (rank == 0) {
        std::cerr << "FAIL: solve() returned on a NaN right-hand side; "
                     "the convergence guard did not fire\n";
    }
    PetscFinalize();
    return 1;
}
