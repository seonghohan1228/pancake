#include "gas_transport.hpp"

#include <algorithm>
#include <cmath>

#include "fvm.hpp"

namespace GasTransport {
namespace {

double field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return std::max(field_value(fields, "mu", cfg.mu, i, j), 1.0e-30);
}

bool uses_pressure_boundary(BCType type) {
    return type == BCType::DIRICHLET || type == BCType::INLET_OUTLET;
}

double full_film_switch(const Field& theta, int i, int j) {
    return theta(i, j) >= 1.0 ? 1.0 : 0.0;
}

// Liquid-solution mass per unit bearing area, rho*h. rho = rho_l * theta is the
// Elrod density, so rho*h is exactly the liquid mass the JFO solve conserves.
double liquid_mass_area(const Fields& fields, const SimulationConfig& cfg, int i, int j,
                        double h_val) {
    return std::max(field_value(fields, "rho", cfg.rho, i, j) * h_val, 1.0e-30);
}

double south_boundary_velocity(const Field& p, const Field& h, const Field& theta,
                               const Fields& fields, const Mesh& mesh,
                               const SimulationConfig& cfg, int i) {
    if (!uses_pressure_boundary(cfg.bc_z_south_type)) return 0.0;
    const double h_face = std::max(h(i, 0), cfg.min_film_thickness);
    const double p_bc = cfg.cavitation_model == CavitationModel::ELROD_ADAMS
        ? cfg.elrod_boundary_pressure(cfg.bc_z_south_val)
        : cfg.bc_z_south_val;
    const double dp_dz = (p(i, 0) - p_bc) / (0.5 * mesh.get_d_z());
    const double mu_face = dynamic_viscosity(fields, cfg, i, 0);
    return -full_film_switch(theta, i, 0) * h_face * h_face * dp_dz / (12.0 * mu_face);
}

double north_boundary_velocity(const Field& p, const Field& h, const Field& theta,
                               const Fields& fields, const Mesh& mesh,
                               const SimulationConfig& cfg, int i) {
    if (!uses_pressure_boundary(cfg.bc_z_north_type)) return 0.0;
    const int j = mesh.n_z_local - 1;
    const double h_face = std::max(h(i, j), cfg.min_film_thickness);
    const double p_bc = cfg.cavitation_model == CavitationModel::ELROD_ADAMS
        ? cfg.elrod_boundary_pressure(cfg.bc_z_north_val)
        : cfg.bc_z_north_val;
    const double dp_dz = (p_bc - p(i, j)) / (0.5 * mesh.get_d_z());
    const double mu_face = dynamic_viscosity(fields, cfg, i, j);
    return -full_film_switch(theta, i, j) * h_face * h_face * dp_dz / (12.0 * mu_face);
}

double south_boundary_mass_flux(const Field& p, const Field& h, const Field& theta,
                                const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, int i) {
    const double u_z = south_boundary_velocity(p, h, theta, fields, mesh, cfg, i);
    const double h_face = std::max(h(i, 0), cfg.min_film_thickness);
    const double rho_face = std::max(field_value(fields, "rho", cfg.rho, i, 0), 1.0e-30);
    return rho_face * h_face * u_z * mesh.R * mesh.get_d_theta();
}

double north_boundary_mass_flux(const Field& p, const Field& h, const Field& theta,
                                const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, int i) {
    const int j = mesh.n_z_local - 1;
    const double u_z = north_boundary_velocity(p, h, theta, fields, mesh, cfg, i);
    const double h_face = std::max(h(i, j), cfg.min_film_thickness);
    const double rho_face = std::max(field_value(fields, "rho", cfg.rho, i, j), 1.0e-30);
    return rho_face * h_face * u_z * mesh.R * mesh.get_d_theta();
}

double south_boundary_volume_flux(const Field& p, const Field& h, const Field& theta,
                                  const Fields& fields, const Mesh& mesh,
                                  const SimulationConfig& cfg, int i) {
    return south_boundary_velocity(p, h, theta, fields, mesh, cfg, i) *
        mesh.R * mesh.get_d_theta();
}

double north_boundary_volume_flux(const Field& p, const Field& h, const Field& theta,
                                  const Fields& fields, const Mesh& mesh,
                                  const SimulationConfig& cfg, int i) {
    return north_boundary_velocity(p, h, theta, fields, mesh, cfg, i) *
        mesh.R * mesh.get_d_theta();
}

// Non-conservative correction: convert the conservative div(F c) assembled by
// FVM::divergence into the material form F.grad(c) by subtracting c_P * div(F).
// Axial boundary faces are open: outflow carries gas out, and inflow carries the
// configured initial dissolved composition.
void add_nonconservative_correction(LinearSystem& sys, const Field& flux_theta,
                                    const Field& flux_z, const Field& p,
                                    const Field& h, const Field& theta,
                                    const Fields& fields,
                                    const SimulationConfig& cfg,
                                    const Mesh& mesh) {
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double flux_e = flux_theta(i, j);
            const double flux_w = flux_theta(i - 1, j);
            const double flux_n = (j + 1 < mesh.n_z_local)
                ? flux_z(i, j)
                : north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
            const double flux_s = (j > 0)
                ? flux_z(i, j - 1)
                : south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
            sys.a_p(i, j) -= (flux_e - flux_w + flux_n - flux_s);
        }
    }
}

void add_dissolved_axial_boundary_convection(LinearSystem& sys, const Field& p,
                                             const Field& h, const Field& theta,
                                             const Fields& fields,
                                             const Mesh& mesh,
                                             const SimulationConfig& cfg) {
    const double c_in = std::clamp(cfg.dissolved_gas_initial, 0.0, cfg.dissolved_gas_max);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double flux_s = south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_s >= 0.0) sys.source(i, 0) += flux_s * c_in;  // inflow
        else sys.a_p(i, 0) += -flux_s;                         // outflow

        const int j = mesh.n_z_local - 1;
        const double flux_n = north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_n < 0.0) sys.source(i, j) += -flux_n * c_in;  // inflow
        else sys.a_p(i, j) += flux_n;                          // outflow
    }
}

void add_free_gas_axial_boundary_convection(LinearSystem& sys, const Field& p,
                                            const Field& h, const Field& theta,
                                            const Fields& fields,
                                            const Mesh& mesh,
                                            const SimulationConfig& cfg) {
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double flux_s = south_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_s < 0.0) sys.a_p(i, 0) += -flux_s;  // south outflow; inflow has zero free gas

        const int j = mesh.n_z_local - 1;
        const double flux_n = north_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_n > 0.0) sys.a_p(i, j) += flux_n;   // north outflow; inflow has zero free gas
    }
}

// Imposes the configured reservoir composition on supply cells. The implied
// total-gas mass change [kg per cell] is accumulated into "gas_clamp_mass" when
// present so the diagnostics balance can account for the overwrite.
void apply_supply_cell_gas_composition(Fields& fields, const Mesh& mesh,
                                       const SimulationConfig& cfg) {
    if (!fields.has("inlet_indicator")) return;
    const double c_in = std::clamp(cfg.dissolved_gas_initial, 0.0, cfg.dissolved_gas_max);
    Field* gas_clamp = fields.has("gas_clamp_mass") ? &fields["gas_clamp_mass"] : nullptr;
    const double cell_area = mesh.cell_volume();
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            if (fields["inlet_indicator"](i, j) <= 0.5) continue;
            if (gas_clamp != nullptr) {
                const double h_val = std::max(field_value(fields, "h", cfg.c, i, j),
                                              cfg.min_film_thickness);
                const double rho_h = liquid_mass_area(fields, cfg, i, j, h_val);
                const double c_old = field_value(fields, "dissolved_gas", 0.0, i, j);
                const double m_old = field_value(fields, "free_gas_mass", 0.0, i, j);
                (*gas_clamp)(i, j) += ((c_in - c_old) * rho_h - m_old) * cell_area;
            }
            if (fields.has("dissolved_gas")) fields["dissolved_gas"](i, j) = c_in;
            if (fields.has("free_gas_mass")) fields["free_gas_mass"](i, j) = 0.0;
            if (fields.has("alpha_gas")) fields["alpha_gas"](i, j) = 0.0;
            if (fields.has("gas_mass_transfer")) fields["gas_mass_transfer"](i, j) = 0.0;
        }
    }
}

}  // namespace

void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh,
           const SimulationConfig& cfg, double dt) {
    if (cfg.fluid_property_model != FluidPropertyModel::GAS_CAVITATION_MIXTURE) return;
    if (!fields.has("dissolved_gas")) return;
    if (dt <= 0.0) return;

    Field& c = fields["dissolved_gas"];
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    Field* gas_clamp = fields.has("gas_clamp_mass") ? &fields["gas_clamp_mass"] : nullptr;
    if (gas_clamp != nullptr) gas_clamp->fill(0.0);
    const double cell_area = mesh.cell_volume();

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double D_eff = std::max(cfg.dissolved_gas_diffusivity, 0.0);

    // Cell coefficients: capacity = rho*h (conserved liquid mass), diffusion = D*rho*h.
    Field capacity("gas_capacity", mesh);
    Field diffusion("gas_diffusion_h", mesh);
    for (int i = -h.n_ghost; i < n_theta + h.n_ghost; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double m = liquid_mass_area(fields, cfg, i, j, h_val);
            diffusion(i, j) = D_eff * m;
            if (i >= 0 && i < n_theta) capacity(i, j) = m;
        }
    }

    // Face mass fluxes F = rho_face * h_face * u_face * (face length), with the same
    // film-averaged velocity (Couette + cavitation-gated Poiseuille) as Reynolds/energy.
    Field flux_theta("gas_flux_theta", mesh);
    Field flux_z("gas_flux_z", mesh);
    Field free_gas_flux_theta("free_gas_flux_theta", mesh);
    Field free_gas_flux_z("free_gas_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i + 1, j)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i + 1, j) >= 1.0) ? 1.0 : 0.0;
            const double dp_dtheta = (p(i + 1, j) - p(i, j)) / d_theta;
            const double mu_face = 0.5 * (dynamic_viscosity(fields, cfg, i, j) +
                                          dynamic_viscosity(fields, cfg, i + 1, j));
            const double rho_face = 0.5 * (field_value(fields, "rho", cfg.rho, i, j) +
                                           field_value(fields, "rho", cfg.rho, i + 1, j));
            const double u_theta = 0.5 * cfg.omega * R
                - g_face * h_face * h_face * dp_dtheta / (12.0 * mu_face * R);
            flux_theta(i, j) = rho_face * h_face * u_theta * d_z;
            free_gas_flux_theta(i, j) = u_theta * d_z;
        }
    }
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            flux_z(i, j) = 0.0;
            free_gas_flux_z(i, j) = 0.0;
        }
        for (int j = 0; j + 1 < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i, j + 1)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i, j + 1) >= 1.0) ? 1.0 : 0.0;
            const double dp_dz = (p(i, j + 1) - p(i, j)) / d_z;
            const double mu_face = 0.5 * (dynamic_viscosity(fields, cfg, i, j) +
                                          dynamic_viscosity(fields, cfg, i, j + 1));
            const double rho_face = 0.5 * (field_value(fields, "rho", cfg.rho, i, j) +
                                           field_value(fields, "rho", cfg.rho, i, j + 1));
            const double u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
            flux_z(i, j) = rho_face * h_face * u_z * R * d_theta;
            free_gas_flux_z(i, j) = u_z * R * d_theta;
        }
    }

    sys.reset();
    if (D_eff > 0.0) FVM::laplacian(sys, diffusion, mesh);
    // TYPE_DIFFERENCING is theta-only (validated); fall back to upwind defensively.
    const ConvectionScheme gas_scheme =
        cfg.gas_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING
            ? ConvectionScheme::UPWIND : cfg.gas_convection_scheme;
    FVM::divergence(sys, flux_theta, flux_z, c, gas_scheme, mesh);
    add_dissolved_axial_boundary_convection(sys, p, h, theta, fields, mesh, cfg);
    if (cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0) {
        FVM::ddt_weighted(sys, c, capacity, dt, mesh);
    }
    add_nonconservative_correction(sys, flux_theta, flux_z, p, h, theta, fields, cfg, mesh);

    sys.solve(c, cfg.linear_rtol);

    const double c_max = std::max(0.0, cfg.dissolved_gas_max);
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double c_clamped = std::clamp(c(i, j), 0.0, c_max);
            if (gas_clamp != nullptr && c_clamped != c(i, j)) {
                const double h_val = std::max(h(i, j), cfg.min_film_thickness);
                (*gas_clamp)(i, j) += (c_clamped - c(i, j)) *
                    liquid_mass_area(fields, cfg, i, j, h_val) * cell_area;
            }
            c(i, j) = c_clamped;
        }
    }
    apply_supply_cell_gas_composition(fields, mesh, cfg);

    if (fields.has("free_gas_mass")) {
        Field& free_gas_mass = fields["free_gas_mass"];
        sys.reset();
        FVM::divergence(sys, free_gas_flux_theta, free_gas_flux_z, free_gas_mass,
                        gas_scheme, mesh);
        add_free_gas_axial_boundary_convection(sys, p, h, theta, fields, mesh, cfg);
        FVM::ddt(sys, free_gas_mass, dt, mesh);
        sys.solve(free_gas_mass, cfg.linear_rtol);

        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double m_clamped = std::max(0.0, free_gas_mass(i, j));
                if (gas_clamp != nullptr && m_clamped != free_gas_mass(i, j)) {
                    (*gas_clamp)(i, j) += (m_clamped - free_gas_mass(i, j)) * cell_area;
                }
                free_gas_mass(i, j) = m_clamped;
            }
        }
        apply_supply_cell_gas_composition(fields, mesh, cfg);
    }
}

}  // namespace GasTransport
