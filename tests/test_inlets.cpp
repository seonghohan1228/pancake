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
        cfg.omega = 0.0; // Pure diffusion test
        
        // Add an axial groove at theta = 180 deg
        InletConfig groove;
        groove.type = InletConfig::Type::GROOVE;
        groove.theta = 180.0;
        groove.size = 10.0; // 10 degrees wide
        groove.p_supply = 2.0e5;
        cfg.inlets.push_back(groove);

        // Add a circular inlet at theta = 90 deg, z = L/2
        InletConfig hole;
        hole.type = InletConfig::Type::CIRCULAR;
        hole.theta = 90.0;
        hole.z = 0.5 * cfg.L;
        hole.size = 0.002; // 2mm radius
        hole.p_supply = 3.0e5;
        cfg.inlets.push_back(hole);

        Mesh mesh(cfg);
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("rho",      mesh).fill(cfg.rho);

        FilmThickness::compute_static(fields["h"], mesh, cfg);
        LinearSystem sys(mesh);

        Reynolds::solve(fields, sys, mesh, cfg);

        const Field& pressure = fields["pressure"];
        const double d_theta = mesh.get_d_theta();
        const double d_z = mesh.get_d_z();
        const double R = mesh.R;

        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Check cells inside the groove
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - M_PI);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            if (dth <= 0.5 * (groove.size * M_PI / 180.0)) {
                for (int j = 0; j < mesh.n_z_local; ++j) {
                    check(std::abs(pressure(i, j) - groove.p_supply) < 1e-6, "Groove pressure mismatch");
                }
            }
        }

        // Check cells inside the hole
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - 0.5 * M_PI);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            for (int j = 0; j < mesh.n_z_local; ++j) {
                const double z_c = (j + 0.5) * d_z;
                const double dist_sq = (R * dth) * (R * dth) + (z_c - hole.z) * (z_c - hole.z);
                if (dist_sq <= hole.size * hole.size) {
                    check(std::abs(pressure(i, j) - hole.p_supply) < 1e-6, "Circular hole pressure mismatch");
                }
            }
        }

        if (rank == 0) std::cout << "PASS: test_inlets\n";
    }

    PetscFinalize();
    return 0;
}
