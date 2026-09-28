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
#include "fbns.hpp"
#include "communicator.hpp"

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

        cfg.cavitation_model = CavitationModel::FBNS;
        cfg.bc_z_south_val = p0;
        cfg.bc_z_north_val = 0.0;
        Communicator comm(mesh);
        comm.update_ghosts(fields);
        Reynolds::calculate_velocities(fields, mesh, cfg);
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j)
                check(std::abs(u_theta(i, j) - u_expected) < 1e-10,
                      "FBNS full-film Couette limit mismatch");
            for (int j = 0; j <= mesh.n_z_local; ++j)
                check(std::abs(u_z(i, j) - uz_expected) < 1e-7,
                      "FBNS full-film Poiseuille limit mismatch");
        }
    }

    {
        SimulationConfig cfg;
        cfg.n_theta_global = 8;
        cfg.n_z_global = 3;
        cfg.cavitation_model = CavitationModel::FBNS;
        cfg.p_cav = 0.0;
        cfg.bc_z_south_val = 1.0e5;
        cfg.bc_z_north_val = 1.0e5;
        Mesh mesh(cfg);
        Communicator comm(mesh);
        Fields fields;
        for (const char* name : {"pressure", "theta", "h", "mu", "rho_liquid_solution",
                                "cross_film_f2", "cross_film_couette_h"}) fields.add(name, mesh);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z", mesh, 2, GridLocation::FACE_Z);
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const int global_i = mesh.offset_theta + i;
            for (int j = 0; j < mesh.n_z_local; ++j) {
                const bool cavitated = global_i % 2 == 0;
                fields["pressure"](i, j) = cavitated ? 0.0 : 2.0e5 * (j + 1);
                fields["theta"](i, j) = cavitated ? 0.4 : 1.0;
                fields["h"](i, j) = cfg.c * (1 + 0.1 * global_i);
                fields["mu"](i, j) = cfg.mu * (1 + 0.2 * global_i);
                fields["rho_liquid_solution"](i, j) = cfg.rho * (1 + 0.01 * global_i);
                fields["cross_film_f2"](i, j) = 2 * std::pow(fields["h"](i, j), 3) /
                    (12 * fields["mu"](i, j));
                fields["cross_film_couette_h"](i, j) = 0.8 * fields["h"](i, j);
            }
        }
        comm.update_ghosts(fields);
        const auto matches_flux = [](double actual, double expected) {
            return std::abs(actual - expected) < 1e-12 * std::max(std::abs(expected), 1e-20);
        };
        for (const auto model : {CrossFilmModel::NONE, CrossFilmModel::POLYNOMIAL}) {
            cfg.cross_film_thermal_model = model;
            for (const double omega : {100.0, -100.0}) {
                cfg.omega = omega;
                Reynolds::calculate_velocities(fields, mesh, cfg);
                for (int i = 0; i < mesh.n_theta_local; ++i) {
                    for (int j = 0; j < mesh.n_z_local; ++j) {
                        const double flux = FBNS::mass_flux_theta(fields, mesh, cfg, i - 1, j);
                        const int upwind = flux >= 0 ? i - 1 : i;
                        const double rho_face = 0.5 * (fields["rho_liquid_solution"](i - 1, j) +
                                                      fields["rho_liquid_solution"](i, j));
                        const double h_face = 0.5 * (fields["h"](i - 1, j) + fields["h"](i, j));
                        const double recovered = fields["velocity_theta"](i, j) * rho_face * h_face *
                            fields["theta"](upwind, j) * mesh.get_d_z();
                        check(matches_flux(recovered, flux),
                              "FBNS theta velocity must recover conserved rupture-front mass flux");
                    }
                    for (int j = 0; j <= mesh.n_z_local; ++j) {
                        const int south = std::max(j - 1, 0);
                        const int north = std::min(j, mesh.n_z_local - 1);
                        const double rho_face = 0.5 * (
                            fields["rho_liquid_solution"](i, south) * fields["theta"](i, south) +
                            fields["rho_liquid_solution"](i, north) * fields["theta"](i, north));
                        const double h_face = 0.5 * (fields["h"](i, south) + fields["h"](i, north));
                        const double flux = FBNS::mass_flux_z(fields, mesh, cfg, i, j - 1);
                        const double recovered = fields["velocity_z"](i, j) * rho_face * h_face *
                            mesh.R * mesh.get_d_theta();
                        check(matches_flux(recovered, flux),
                              "FBNS z velocity must recover interior and reformation boundary mass flux");
                    }
                }
            }
        }
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "PASS: test_velocity\n";
    }

    PetscFinalize();
    return 0;
}
