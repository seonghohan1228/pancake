#include "energy.hpp"

#include <algorithm>
#include <cmath>

#include "fvm.hpp"

namespace Energy {
namespace {

double full_film_switch(const Field& theta, int i, int j) {
    return theta(i, j) >= 1.0 ? 1.0 : 0.0;
}

double dp_dz_at_cell(const Field& p, int i, int j, int n_z, double d_z,
                     const SimulationConfig& cfg) {
    if (n_z == 1) return 0.0;
    if (j == 0) {
        const double p_s = (cfg.bc_z_south_type == BCType::DIRICHLET) ? cfg.bc_z_south_val : p(i, 0);
        return (p(i, 0) - p_s) / (0.5 * d_z);
    }
    if (j == n_z - 1) {
        const double p_n = (cfg.bc_z_north_type == BCType::DIRICHLET) ? cfg.bc_z_north_val : p(i, n_z - 1);
        return (p_n - p(i, n_z - 1)) / (0.5 * d_z);
    }
    return (p(i, j + 1) - p(i, j - 1)) / (2.0 * d_z);
}

void compute_heat_generation(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    Field& heat = fields["heat_generation"];

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double U = cfg.omega * R;

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double g = full_film_switch(theta, i, j);
            const double dp_dtheta = (p(i + 1, j) - p(i - 1, j)) / (2.0 * d_theta);
            const double dp_ds = dp_dtheta / R;
            const double dp_dz = dp_dz_at_cell(p, i, j, n_z, d_z, cfg);
            const double couette = cfg.mu * U * U / h_val;
            const double poiseuille = g * h_val * h_val * h_val *
                (dp_ds * dp_ds + dp_dz * dp_dz) / (12.0 * cfg.mu);
            heat(i, j) = couette + poiseuille;
        }
    }
}

void add_wall_heat_loss(LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg) {
    const double h_wall = cfg.journal_heat_transfer + cfg.bearing_heat_transfer;
    if (h_wall <= 0.0) return;

    const double wall_temperature =
        (cfg.journal_heat_transfer * cfg.journal_wall_temperature +
         cfg.bearing_heat_transfer * cfg.bearing_wall_temperature) / h_wall;
    const double area = mesh.R * mesh.get_d_theta() * mesh.get_d_z();
    const double coeff = h_wall * area;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            sys.a_p(i, j) += coeff;
            sys.source(i, j) += coeff * wall_temperature;
        }
    }
}

}  // namespace

void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg) {
    compute_heat_generation(fields, mesh, cfg);
    if (cfg.temperature_model == TemperatureModel::ISOTHERMAL) return;

    Field& temperature = fields["temperature"];
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const Field& heat = fields["heat_generation"];

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;

    Field conductivity("thermal_conductivity_h", mesh);
    Field capacity("thermal_capacity", mesh);
    for (int i = -h.n_ghost; i < n_theta + h.n_ghost; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            conductivity(i, j) = cfg.thermal_conductivity * h_val;
            if (i >= 0 && i < n_theta) capacity(i, j) = cfg.rho_cp * h_val;
        }
    }

    Field flux_theta("thermal_flux_theta", mesh);
    Field flux_z("thermal_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i + 1, j)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i + 1, j) >= 1.0) ? 1.0 : 0.0;
            const double dp_dtheta = (p(i + 1, j) - p(i, j)) / d_theta;
            const double u_theta = 0.5 * cfg.omega * R
                - g_face * h_face * h_face * dp_dtheta / (12.0 * cfg.mu * R);
            flux_theta(i, j) = cfg.rho_cp * h_face * u_theta * d_z;
        }
    }

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) flux_z(i, j) = 0.0;
        for (int j = 0; j + 1 < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i, j + 1)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i, j + 1) >= 1.0) ? 1.0 : 0.0;
            const double dp_dz = (p(i, j + 1) - p(i, j)) / d_z;
            const double u_z = -g_face * h_face * h_face * dp_dz / (12.0 * cfg.mu);
            flux_z(i, j) = cfg.rho_cp * h_face * u_z * R * d_theta;
        }
    }

    sys.reset();
    FVM::laplacian(sys, conductivity, mesh);
    FVM::divergence(sys, flux_theta, flux_z, temperature, ConvectionScheme::UPWIND, mesh);
    if (cfg.solution_mode != SolutionMode::STEADY_STATE) {
        FVM::ddt_weighted(sys, temperature, capacity, cfg.dt, mesh);
    }
    FVM::add_source(sys, heat, mesh);
    add_wall_heat_loss(sys, mesh, cfg);

    sys.solve(temperature);
}

}  // namespace Energy
