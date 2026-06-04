#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

namespace {
constexpr double pi = 3.14159265358979323846264338327950288;

void check(bool condition, const char* message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cerr << "rank " << rank << " FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

Fields make_force_fields(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);
    fields.add("pressure_force_x", mesh).fill(0.0);
    fields.add("pressure_force_y", mesh).fill(0.0);
    fields.add("pressure_force_z", mesh).fill(0.0);
    fields.add("load_x", mesh).fill(0.0);
    fields.add("load_y", mesh).fill(0.0);
    fields.add("load_z", mesh).fill(0.0);
    fields.add("viscous_force_x", mesh).fill(0.0);
    fields.add("viscous_force_y", mesh).fill(0.0);
    fields.add("viscous_force_z", mesh).fill(0.0);
    fields.add("fluid_force_x", mesh).fill(0.0);
    fields.add("fluid_force_y", mesh).fill(0.0);
    fields.add("fluid_force_z", mesh).fill(0.0);
    fields.add("friction_torque", mesh).fill(0.0);
    return fields;
}
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    {
        SimulationConfig cfg;
        cfg.n_theta_global = 240;
        cfg.n_z_global = 20;
        cfg.R = 0.01;
        cfg.c = 1.0e-4;
        cfg.e = 0.0;
        cfg.L = 0.02;
        cfg.mu = 0.01;
        cfg.omega = 0.0;
        cfg.p_cav = 0.0;

        Mesh mesh(cfg);
        Fields fields = make_force_fields(mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);

        const double p0 = 2.0e5;
        const double d_theta = mesh.get_d_theta();
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta = (mesh.offset_theta + i + 0.5) * d_theta;
            for (int j = 0; j < mesh.n_z_local; ++j) {
                fields["pressure"](i, j) = p0 * (1.0 + std::cos(theta));
            }
        }

        const Reynolds::ForceComponents forces =
            Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        const double expected_fx = p0 * pi * (cfg.R + cfg.c) * cfg.L;
        check(std::abs(forces.pressure_x - expected_fx) / expected_fx < 1.0e-3,
              "pressure force x should include bearing surface radius");
        check(std::abs(forces.pressure_y) < 1.0e-8 * expected_fx,
              "pressure force y should cancel by symmetry");
        check(std::abs(fields["pressure_force_x"](0, 0) - forces.pressure_x) < 1.0e-8,
              "pressure_force_x field should store pressure force resultant");
        check(std::abs(fields["load_x"](0, 0) - forces.pressure_x) < 1.0e-8,
              "legacy load_x field should alias pressure force resultant");
    }

    {
        SimulationConfig cfg;
        cfg.n_theta_global = 120;
        cfg.n_z_global = 40;
        cfg.R = 0.01;
        cfg.c = 1.0e-4;
        cfg.e = 0.0;
        cfg.L = 0.02;
        cfg.mu = 0.01;
        cfg.omega = 0.0;
        cfg.p_cav = 0.0;
        const double p0 = 1.0e5;
        cfg.bc_z_south_val = p0;
        cfg.bc_z_north_val = 0.0;

        Mesh mesh(cfg);
        Fields fields = make_force_fields(mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);

        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                const double z = (j + 0.5) * mesh.get_d_z();
                fields["pressure"](i, j) = p0 * (1.0 - z / cfg.L);
            }
        }

        const Reynolds::ForceComponents forces =
            Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        const double expected_fz = pi * (cfg.R + cfg.c) * cfg.c * p0;
        check(std::abs(forces.viscous_z - expected_fz) / expected_fz < 1.0e-3,
              "viscous z force should use bearing-side pressure-gradient shear");
    }

    {
        SimulationConfig cfg;
        cfg.n_theta_global = 120;
        cfg.n_z_global = 20;
        cfg.R = 0.01;
        cfg.c = 1.0e-4;
        cfg.e = 0.0;
        cfg.L = 0.02;
        cfg.mu = 0.01;
        cfg.omega = 100.0;
        cfg.p_cav = 0.0;

        Mesh mesh(cfg);
        Fields fields = make_force_fields(mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);

        const Reynolds::ForceComponents forces =
            Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        const double expected_torque =
            2.0 * pi * cfg.mu * cfg.omega * cfg.R * cfg.R * cfg.R * cfg.L / cfg.c;
        check(std::abs(forces.friction_torque - expected_torque) / expected_torque < 1.0e-12,
              "friction torque should include shaft surface geometry");
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "PASS: test_forces\n";

    PetscFinalize();
    return 0;
}
