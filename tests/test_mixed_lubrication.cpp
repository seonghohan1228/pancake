// Mixed lubrication closures: Greenwood-Tripp asperity contact and
// Patir-Cheng flow / shear-stress factors.
//
// Test 1 - Closed-form checks: exact F5/2 against SciPy values, branch continuity of the published fits,
//   smooth-surface limits, and a hand-computed contact pressure.
// Test 2 - Models off: every factor is exactly 1 and contact exactly 0.
// Test 3 - Uniform gap: contact force cancels by symmetry and the contact
//   torque equals f_b p_c 2 pi R^2 L.
// Test 4 - Eccentric gap: contact pushes the bearing toward the minimum gap and
//   matches an independent fine quadrature; it is included in fluid_force.
// Test 5 - Flux factors: Poiseuille mobility scales by phi_x and the Couette
//   flux gains (U/2) sigma phi_s for unequal roughness.
// Test 6 - Validation rejects flow factors with COUPLED and missing modulus.
//
// Run as: mpirun -n <N> ./test_mixed_lubrication

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"
#include "mixed_lubrication.hpp"
#include "reynolds.hpp"
#include "reynolds_coefficients.hpp"

namespace {
constexpr double pi = 3.14159265358979323846264338327950288;

void check(bool condition, const std::string& message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cerr << "rank " << rank << " FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

bool close(double a, double b, double rel) {
    return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

SimulationConfig rough_config() {
    SimulationConfig cfg;
    cfg.n_theta_global = 360;
    cfg.n_z_global = 8;
    cfg.R = 0.0215;
    cfg.c = 4.25e-5;
    cfg.e = 0.0;
    cfg.L = 0.04;
    cfg.mu = 0.008;
    cfg.rho = 930.0;
    cfg.omega = 20.0;
    cfg.p_cav = 0.0;
    cfg.cavitation_model = CavitationModel::FBNS;
    cfg.asperity_contact_model = AsperityContactModel::GREENWOOD_TRIPP;
    cfg.flow_factor_model = FlowFactorModel::PATIR_CHENG;
    cfg.roughness_rq_shaft = 0.6e-6;
    cfg.roughness_rq_bush = 0.8e-6;   // sigma = 1.0 um
    cfg.asperity_eta_beta_sigma = 0.04;
    cfg.asperity_sigma_over_beta = 1e-3;
    cfg.contact_modulus = 107e9;
    cfg.boundary_friction_coefficient = 0.1;
    return cfg;
}

Fields force_fields(const Mesh& mesh) {
    Fields fields;
    fields.add("pressure", mesh).fill(0.0);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(0.0);
    fields.add("rho", mesh).fill(930.0);
    for (const char* name : {"pressure_force_x", "pressure_force_y", "pressure_force_z",
                             "viscous_force_x", "viscous_force_y", "viscous_force_z",
                             "fluid_force_x", "fluid_force_y", "fluid_force_z", "friction_torque",
                             "contact_pressure", "contact_force_x", "contact_force_y",
                             "contact_force_z", "contact_friction_torque"})
        fields.add(name, mesh).fill(0.0);
    return fields;
}
}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    using namespace MixedLubrication;

    // Test 1
    // Exact F5/2 reference values: SciPy quad, relative tolerance 1e-13.
    {
        const double h_ref[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0};
        const double f_ref[] = {0.616634218996844, 0.24040157067407217, 0.08056233557081148,
                                0.005423705197484413, 0.00017087299621442443,
                                2.3533810525876063e-06, 3.142968531919884e-11};
        for (int k = 0; k < 7; ++k) {
            check(close(greenwood_tripp_f52_exact(h_ref[k]), f_ref[k], 1e-6), "exact F5/2 quadrature");
            check(close(greenwood_tripp_f52(h_ref[k]), f_ref[k], 1e-4), "tabulated F5/2");
        }
        check(greenwood_tripp_f52(-1.0) == greenwood_tripp_f52(0.0), "F5/2 clamps H < 0 to 0");
        check(greenwood_tripp_f52(8.0) == 0.0 && greenwood_tripp_f52(9.0) == 0.0, "F5/2 is zero for H >= 8");
        double previous = greenwood_tripp_f52(0.0);
        for (double h = 0.01; h < 8.0; h += 0.01) {
            const double value = greenwood_tripp_f52(h);
            check(value < previous && value > 0.0, "F5/2 decreases monotonically on [0, 8)");
            previous = value;
        }
    }
    check(close(1.899 * std::pow(5.0, 0.98) * std::exp(-0.92 * 5 + 0.05 * 25), 1.126 * std::exp(-1.25), 1e-3),
          "Patir-Cheng shear flow factor branches meet at H = 5");
    check(close(shear_stress_factor_f(3.0 - 1e-9), shear_stress_factor_f(3.0 + 1e-9), 1e-6) &&
          close(shear_stress_factor_f(3.0), 35.0 / 32.0 * 16.0 / 15.0, 1e-9),
          "phi_f branches meet at H = 3 with value 7/6");
    check(std::abs(shear_stress_factor_f(30.0) - 1.0) < 2e-3, "phi_f -> 1 for a thick film");
    check(std::abs(pressure_flow_factor_isotropic(40.0) - 1.0) < 1e-9 &&
          std::abs(shear_stress_factor_fp(40.0) - 1.0) < 1e-9, "pressure factors -> 1 for a thick film");
    check(shear_stress_factor_fs_isotropic(7.5) == 0.0 && shear_stress_factor_fs_isotropic(3.0) > 0.0,
          "Phi_fs is positive and cut off above H = 7");
    {
        SimulationConfig cfg = rough_config();
        const double k = 16.0 * std::sqrt(2.0) * pi / 15.0 * 0.04 * 0.04 * std::sqrt(1e-3);
        check(close(contact_pressure(2.0e-6, cfg), k * 107e9 * greenwood_tripp_f52(2.0), 1e-12),
              "Greenwood-Tripp contact pressure at H = 2");
        check(contact_pressure(8.0e-6, cfg) == 0.0, "no contact at H >= 8");
        check(contact_pressure(4.0e-6, cfg) > 0.0, "contact onset tail at H = 4 is kept");
        SimulationConfig direct = cfg;
        direct.asperity_gt_k = 1e-3;
        check(close(contact_pressure(2.0e-6, direct), 1e-3 * 107e9 * greenwood_tripp_f52(2.0), 1e-12),
              "asperity_gt_k overrides the eta beta sigma form");
    }

    // Test 1b - legacy codex/gt-contact keys map onto this model with warnings.
    {
        SimulationConfig cfg;
        cfg.load_from_text("contact_model = GREENWOOD_TRIPP\ncontact_roughness_sigma = 0.4e-6\n"
                           "contact_gt_k = 0.001\ncontact_modulus = 1.099e11\n"
                           "contact_boundary_friction = 0.07\n");
        check(cfg.parse_errors.empty(), "legacy contact keys parse");
        check(cfg.parse_warnings.size() == 4, "each legacy contact key warns once");
        check(cfg.asperity_contact_model == AsperityContactModel::GREENWOOD_TRIPP &&
              cfg.asperity_gt_k == 0.001 && cfg.boundary_friction_coefficient == 0.07,
              "legacy contact keys set the merged fields");
        check(close(composite_roughness(cfg), 0.4e-6, 1e-12) &&
              roughness_variance_difference(cfg) == 0.0, "legacy sigma splits equally");
        check(close(contact_pressure(1.6e-6, cfg), 1e-3 * 1.099e11 * greenwood_tripp_f52(4.0), 1e-12),
              "legacy inputs give the codex/gt-contact pressure law");
    }

    // Test 2
    {
        SimulationConfig cfg;
        check(contact_pressure(1e-9, cfg) == 0.0 && pressure_flow_factor(1e-9, cfg) == 1.0 &&
              shear_flow_thickness(1e-9, cfg) == 0.0 && couette_shear_factor(1e-9, cfg) == 1.0 &&
              pressure_shear_factor(1e-9, cfg) == 1.0, "models off reproduce smooth surfaces");
    }

    // Test 3
    {
        SimulationConfig cfg = rough_config();
        cfg.flow_factor_model = FlowFactorModel::NONE;
        Mesh mesh(cfg);
        Fields fields = force_fields(mesh);
        const double h0 = 2.0e-6;
        fields["h"].fill(h0);
        const auto forces = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        const double pc = contact_pressure(h0, cfg);
        const double expected_torque = 0.1 * pc * 2.0 * pi * cfg.R * cfg.R * cfg.L;
        check(close(forces.contact_torque, expected_torque, 1e-9), "uniform-gap contact torque");
        check(std::abs(forces.contact_x) < 1e-9 * pc * cfg.R * cfg.L &&
              std::abs(forces.contact_y) < 1e-9 * pc * cfg.R * cfg.L, "uniform-gap contact force cancels");
        check(close(forces.friction_torque - forces.contact_torque,
                    cfg.mu * cfg.omega * cfg.R / h0 * 2.0 * pi * cfg.R * cfg.R * cfg.L, 1e-9),
              "hydrodynamic torque unchanged by contact");
        check(fields["contact_pressure"](0, 0) == pc && fields["contact_friction_torque"](0, 0) == forces.contact_torque,
              "contact fields store the cell pressure and the torque");
    }

    // Test 4
    {
        SimulationConfig cfg = rough_config();
        cfg.flow_factor_model = FlowFactorModel::NONE;
        cfg.omega = 0.0;
        Mesh mesh(cfg);
        Fields fields = force_fields(mesh);
        const double eps = 0.94;
        const double d_theta = mesh.get_d_theta();
        for (int i = -1; i <= mesh.n_theta_local; ++i) {  // h_theta reads one ghost layer
            const double theta = (mesh.offset_theta + i + 0.5) * d_theta;
            for (int j = 0; j < mesh.n_z_local; ++j) fields["h"](i, j) = cfg.c * (1.0 + eps * std::cos(theta));
        }
        const auto forces = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        // Independent quadrature: F_x = L int p_c(h) [(R+h) cos t + h' sin t] dt with h' = -c eps sin t.
        double expected_x = 0.0;
        const int n = 200000;
        for (int k = 0; k < n; ++k) {
            const double t = (k + 0.5) * 2.0 * pi / n;
            const double h = cfg.c * (1.0 + eps * std::cos(t));
            const double hp = -cfg.c * eps * std::sin(t);
            expected_x += contact_pressure(h, cfg) * ((cfg.R + h) * std::cos(t) + hp * std::sin(t)) * 2.0 * pi / n;
        }
        expected_x *= cfg.L;
        check(forces.contact_x < 0.0, "contact pushes the bearing toward the minimum gap at theta = pi");
        check(close(forces.contact_x, expected_x, 2e-2), "eccentric contact force matches quadrature");
        check(std::abs(forces.contact_y) < 1e-6 * std::abs(forces.contact_x), "eccentric contact force is symmetric");
        check(close(forces.fluid_x, forces.pressure_x + forces.viscous_x + forces.contact_x, 1e-12),
              "fluid_force includes the contact reaction for the load balance");
    }

    // Test 5
    {
        SimulationConfig cfg = rough_config();
        Mesh mesh(cfg);
        Fields fields;
        fields.add("h", mesh).fill(3.0e-6);
        fields.add("rho", mesh).fill(cfg.rho);
        const double smooth = 3.0e-6 * 3.0e-6 * 3.0e-6 / (12.0 * cfg.mu);
        check(close(ReynoldsCoeff::poiseuille(fields, cfg, 1.0, 0, 0),
                    smooth * pressure_flow_factor_isotropic(3.0), 1e-12), "Poiseuille mobility scaled by phi_x");
        const double rough_flux = ReynoldsCoeff::couette_face_flux(fields, mesh, cfg, 0, 0);
        cfg.flow_factor_model = FlowFactorModel::NONE;
        const double smooth_flux = ReynoldsCoeff::couette_face_flux(fields, mesh, cfg, 0, 0);
        const double vr = (0.36 - 0.64) / 1.0;
        const double expected = cfg.rho * cfg.omega * cfg.R * 0.5 * mesh.get_d_z() *
                                1.0e-6 * vr * shear_flow_factor_isotropic(3.0);
        check(close(rough_flux - smooth_flux, expected, 1e-9) && rough_flux < smooth_flux,
              "rough stationary bush reduces Couette transport by (U/2) sigma phi_s");
    }

    // Test 6
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig cfg = rough_config();
        cfg.cavitation_model = CavitationModel::COUPLED;
        cfg.contact_modulus = 0.0;
        cfg.validate(errors, warnings);
        bool coupled = false, modulus = false;
        for (const auto& e : errors) {
            coupled |= e.find("flow_factor_model = PATIR_CHENG requires") != std::string::npos;
            modulus |= e.find("contact_modulus > 0") != std::string::npos;
        }
        check(coupled && modulus, "validation rejects COUPLED flow factors and a missing modulus");
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_mixed_lubrication passed\n";
    PetscFinalize();
    return 0;
}
