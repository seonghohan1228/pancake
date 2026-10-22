#include <cmath>
#include <iostream>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"
#include "contact.hpp"
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

Fields make_fields(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields fields;
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);
    fields.add("mu", mesh).fill(cfg.mu);
    return fields;
}

SimulationConfig base_config() {
    SimulationConfig cfg;
    cfg.n_theta_global = 240;
    cfg.n_z_global = 20;
    cfg.R = 0.0215;
    cfg.L = 0.04;
    cfg.c = 1.0e-6;
    cfg.e = 0.0;
    cfg.mu = 0.1;
    cfg.omega = 20.0;
    cfg.p_cav = 0.0;
    return cfg;
}
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // F_5/2(0) = 2^(3/4) Gamma(7/4) / sqrt(2 pi) for a Gaussian height distribution.
    const double f0_exact = std::pow(2.0, 0.75) * std::tgamma(1.75) / std::sqrt(2.0 * pi);
    check(std::abs(Contact::f52(0.0) / f0_exact - 1.0) < 1.0e-6, "F_5/2(0) analytic value");
    for (double lambda : {0.37, 1.0, 2.5, 4.0, 6.0}) {
        const double exact = Contact::f52_exact(lambda);
        check(std::abs(Contact::f52(lambda) / exact - 1.0) < 2.0e-3, "F_5/2 table interpolation");
        check(Contact::f52(lambda) < Contact::f52(lambda - 0.1), "F_5/2 decreases with lambda");
    }
    check(Contact::f52(9.0) == 0.0, "F_5/2 is zero beyond the table");

    // Concentric uniform gap: closed-form contact load and torque, zero net force.
    {
        SimulationConfig cfg = base_config();
        cfg.contact_model = ContactModel::GREENWOOD_TRIPP;
        cfg.contact_roughness_sigma = 0.5e-6;
        cfg.contact_gt_k = 1.0e-3;
        cfg.contact_modulus = 1.1e11;
        cfg.contact_boundary_friction = 0.07;
        Mesh mesh(cfg);
        Fields fields = make_fields(mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        const auto forces = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        const double p_c = cfg.contact_gt_k * cfg.contact_modulus * Contact::f52_exact(2.0);
        const double load = p_c * 2.0 * pi * cfg.R * cfg.L;
        check(std::abs(forces.contact_load / load - 1.0) < 2.0e-3, "uniform contact load");
        check(std::abs(forces.contact_torque / (cfg.contact_boundary_friction * load * cfg.R) - 1.0) < 2.0e-3,
              "boundary torque = mu_b W_c R");
        check(std::hypot(forces.contact_x, forces.contact_y) < 1.0e-9 * load, "concentric net contact force vanishes");
    }

    // Eccentric gap: net contact force pushes the bearing toward the narrow gap side
    // reaction, i.e. along +x for minimum gap at theta = 0; NONE gives zero contact.
    {
        SimulationConfig cfg = base_config();
        cfg.c = 2.0e-6;
        cfg.e = 1.5e-6;
        cfg.attitude_angle_deg = 0.0;
        cfg.contact_model = ContactModel::GREENWOOD_TRIPP;
        Mesh mesh(cfg);
        Fields fields = make_fields(mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        const auto on = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        check(on.contact_x > 0.0 && std::abs(on.contact_x) > 10.0 * std::abs(on.contact_z),
              "contact force points toward the minimum gap");
        check(std::abs(on.fluid_x - (on.pressure_x + on.viscous_x + on.contact_x)) <= 1e-9 * std::abs(on.contact_x),
              "fluid force includes contact force");

        cfg.contact_model = ContactModel::NONE;
        const auto off = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        check(off.contact_load == 0.0 && off.contact_torque == 0.0 && off.contact_x == 0.0, "NONE has no contact");
        check(std::abs((on.friction_torque - on.contact_torque) - off.friction_torque)
                  <= 1e-12 * std::abs(off.friction_torque), "film torque unchanged by contact model");
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_contact passed\n";
    PetscFinalize();
    return 0;
}
