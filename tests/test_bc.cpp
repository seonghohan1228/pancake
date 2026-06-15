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
        cfg.cavitation_model = CavitationModel::GUMBEL;
        cfg.omega = 0.0; // Pure diffusion
        cfg.bc_z_south_type = BCType::DIRICHLET;
        cfg.bc_z_south_val = 1.0;
        cfg.bc_z_north_type = BCType::NEUMANN; // Zero gradient at north

        Mesh mesh(cfg);
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("theta",    mesh).fill(1.0);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("rho",      mesh).fill(cfg.rho);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z",     mesh, 2, GridLocation::FACE_Z);

        FilmThickness::compute_static(fields["h"], mesh, cfg);
        LinearSystem sys(mesh);

        Reynolds::solve(fields, sys, mesh, cfg);

        const Field& pressure = fields["pressure"];
        
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // For pure diffusion with DIRICHLET on one side and NEUMANN on the other,
        // and NO source terms (since h is constant and omega is 0),
        // the solution should be UNIFORM pressure equal to the Dirichlet value.
        
        double max_err = 0.0;
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                max_err = std::max(max_err, std::abs(pressure(i, j) - cfg.bc_z_south_val));
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (rank == 0) std::cout << "  max error vs uniform p0 = " << max_err << "\n";
        check(max_err < 1e-6, "Neumann boundary did not result in uniform pressure");

        // Now test DIRICHLET on both sides with different values
        cfg.bc_z_north_type = BCType::DIRICHLET;
        cfg.bc_z_north_val = 0.0;
        
        Reynolds::solve(fields, sys, mesh, cfg);
        
        // Should be linear: p(z) = p_south * (1 - z/L) + p_north * (z/L)
        // Here p(z) = 1.0 * (1 - z/L)
        max_err = 0.0;
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                double z_c = (j + 0.5) * mesh.get_d_z();
                double p_expected = cfg.bc_z_south_val * (1.0 - z_c / cfg.L);
                max_err = std::max(max_err, std::abs(pressure(i, j) - p_expected));
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank == 0) std::cout << "  max error vs linear profile = " << max_err << "\n";
        check(max_err < 1e-6, "Dirichlet boundaries did not result in linear pressure profile");

        if (rank == 0) std::cout << "PASS: test_bc\n";
    }

    PetscFinalize();
    return 0;
}
