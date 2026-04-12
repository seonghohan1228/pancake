#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

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

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    {
        SimulationConfig cfg;
        cfg.e = 0.0; // Concentric
        cfg.omega = 100.0;
        cfg.R = 0.01;
        cfg.mu = 0.01;
        cfg.p_cav = 0.0;

        Mesh mesh(cfg);
        Fields fields;
        fields.add("pressure", mesh).fill(0.0);
        fields.add("theta",    mesh).fill(1.0);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z",     mesh, 2, GridLocation::FACE_Z);

        // Test 1: Uniform Couette flow (no pressure gradient)
        Reynolds::calculate_velocities(fields, mesh, cfg);

        double u_expected = 0.5 * cfg.omega * cfg.R;
        const Field& u_theta = fields["velocity_theta"];
        const Field& u_z     = fields["velocity_z"];

        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                check(std::abs(u_theta(i, j) - u_expected) < 1e-10, "Uniform u_theta mismatch");
            }
        }
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j <= mesh.n_z_local; ++j) {
                check(std::abs(u_z(i, j)) < 1e-10, "Uniform u_z mismatch");
            }
        }

        // Test 2: Pressure gradient in z
        // Set linear pressure profile p(z) = p0 * (1 - z/L)
        // dp/dz = -p0/L
        // uz = - (h^2 / 12mu) * (-p0/L) = h^2 * p0 / (12 * mu * L)
        double p0 = 1.0e5;
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                double z_c = (j + 0.5) * mesh.get_d_z();
                fields["pressure"](i, j) = p0 * (1.0 - z_c / cfg.L);
            }
        }
        
        Reynolds::calculate_velocities(fields, mesh, cfg);
        
        double uz_expected = (cfg.c * cfg.c * p0) / (12.0 * cfg.mu * cfg.L);
        // Note: our calculate_velocities handles boundaries too.
        // For interior faces (j=1 to n_z-1), dp/dz should be exactly -p0/L
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 1; j < mesh.n_z_local; ++j) {
                check(std::abs(u_z(i, j) - uz_expected) < 1e-7, "Poiseuille u_z mismatch");
            }
        }

        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "PASS: test_velocity\n";
    }

    PetscFinalize();
    return 0;
}
