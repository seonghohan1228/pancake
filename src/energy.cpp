#include "energy.hpp"

#include <algorithm>
#include <cmath>

#include "fvm.hpp"

namespace Energy {
namespace {

double full_film_switch(const Field& theta, int i, int j) {
    return theta(i, j) >= 1.0 ? 1.0 : 0.0;
}

double full_film_face_switch(const Field& theta, int i0, int j0, int i1, int j1) {
    return (theta(i0, j0) >= 1.0 && theta(i1, j1) >= 1.0) ? 1.0 : 0.0;
}

bool uses_pressure_boundary(BCType type) {
    return type == BCType::DIRICHLET || type == BCType::INLET_OUTLET;
}

// Pressure the flow solver actually imposed at an axial boundary. Elrod derives
// its boundary film content from the configured pressure and bulk modulus, with
// pressure clamped at p_cav, so thermal fluxes use the same effective pressure.
double solved_boundary_pressure(const SimulationConfig& cfg, double bc_val) {
    if (cfg.cavitation_model == CavitationModel::JFO) {
        return cfg.elrod_boundary_pressure(bc_val);
    }
    return bc_val;
}

double field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double positive_field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return std::max(field_value(fields, name, fallback, i, j), 1.0e-30);
}

double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "mu", cfg.mu, i, j);
}

double volumetric_heat_capacity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "rho_cp_liquid_solution", cfg.rho_cp, i, j);
}

double thermal_conductivity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "k_liquid_solution", cfg.thermal_conductivity, i, j);
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
            const double mu_val = dynamic_viscosity(fields, cfg, i, j);
            // STRIATED: in a cavitated (striated) cell only the liquid fraction theta
            // carries Couette shear; full-film cells (g=1) are unchanged.
            const double shear_w = (cfg.cavitated_film_model == CavitatedFilmModel::STRIATED && g <= 0.0)
                ? std::max(theta(i, j), 0.0) : 1.0;
            const double couette = mu_val * U * U / h_val * shear_w;

            // Pressure-flow dissipation must follow the same active faces as the
            // Reynolds/velocity solve. Cavitated faces carry no Poiseuille flow, so
            // their pressure jump is not a cell-centered gradient. Average squared
            // active face gradients; inactive and exterior faces contribute zero.
            double poiseuille = 0.0;
            if (g > 0.0) {
                double grad_sq = 0.0;
                const double g_e = full_film_face_switch(theta, i, j, i + 1, j);
                const double g_w = full_film_face_switch(theta, i - 1, j, i, j);
                if (g_e > 0.0) {
                    const double dp_ds_e = (p(i + 1, j) - p(i, j)) / (d_theta * R);
                    grad_sq += 0.5 * dp_ds_e * dp_ds_e;
                }
                if (g_w > 0.0) {
                    const double dp_ds_w = (p(i, j) - p(i - 1, j)) / (d_theta * R);
                    grad_sq += 0.5 * dp_ds_w * dp_ds_w;
                }

                if (n_z > 1) {
                    if (j == 0) {
                        const double g_n = full_film_face_switch(theta, i, j, i, j + 1);
                        if (g_n > 0.0) {
                            const double dp_dz_n = (p(i, j + 1) - p(i, j)) / d_z;
                            grad_sq += dp_dz_n * dp_dz_n;
                        }
                    } else if (j == n_z - 1) {
                        const double g_s = full_film_face_switch(theta, i, j - 1, i, j);
                        if (g_s > 0.0) {
                            const double dp_dz_s = (p(i, j) - p(i, j - 1)) / d_z;
                            grad_sq += dp_dz_s * dp_dz_s;
                        }
                    } else {
                        const double g_n = full_film_face_switch(theta, i, j, i, j + 1);
                        const double g_s = full_film_face_switch(theta, i, j - 1, i, j);
                        if (g_n > 0.0) {
                            const double dp_dz_n = (p(i, j + 1) - p(i, j)) / d_z;
                            grad_sq += 0.5 * dp_dz_n * dp_dz_n;
                        }
                        if (g_s > 0.0) {
                            const double dp_dz_s = (p(i, j) - p(i, j - 1)) / d_z;
                            grad_sq += 0.5 * dp_dz_s * dp_dz_s;
                        }
                    }
                }
                poiseuille = h_val * h_val * h_val * grad_sq / (12.0 * mu_val);
            }
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

double south_boundary_flux(const Field& p, const Field& h, const Field& theta,
                           const Fields& fields, const Mesh& mesh,
                           const SimulationConfig& cfg, int i) {
    if (!uses_pressure_boundary(cfg.bc_z_south_type)) return 0.0;
    const double face_length = mesh.R * mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double h_face = std::max(h(i, 0), cfg.min_film_thickness);
    const double mu_face = dynamic_viscosity(fields, cfg, i, 0);
    const double rho_cp_face = volumetric_heat_capacity(fields, cfg, i, 0);
    // Cavitated boundary cell under consistent_boundary_flux re-floods at the Elrod
    // reformation rate; the inflowing reservoir liquid carries temperature_reference
    // (RESERVOIR mode) through add_axial_boundary_convection below.
    double u_z;
    if (cfg.is_reformation_boundary(theta(i, 0))) {
        u_z = cfg.reformation_boundary_velocity(cfg.bc_z_south_val, h_face, mu_face, 0.5 * d_z, true);
    } else {
        const double g_face = full_film_switch(theta, i, 0);
        const double p_bc = solved_boundary_pressure(cfg, cfg.bc_z_south_val);
        const double dp_dz = (p(i, 0) - p_bc) / (0.5 * d_z);
        u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
    }
    return rho_cp_face * h_face * u_z * face_length;
}

double north_boundary_flux(const Field& p, const Field& h, const Field& theta,
                           const Fields& fields, const Mesh& mesh,
                           const SimulationConfig& cfg, int i) {
    if (!uses_pressure_boundary(cfg.bc_z_north_type)) return 0.0;
    const int j = mesh.n_z_local - 1;
    const double face_length = mesh.R * mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double h_face = std::max(h(i, j), cfg.min_film_thickness);
    const double mu_face = dynamic_viscosity(fields, cfg, i, j);
    const double rho_cp_face = volumetric_heat_capacity(fields, cfg, i, j);
    double u_z;
    if (cfg.is_reformation_boundary(theta(i, j))) {
        u_z = cfg.reformation_boundary_velocity(cfg.bc_z_north_val, h_face, mu_face, 0.5 * d_z, false);
    } else {
        const double g_face = full_film_switch(theta, i, j);
        const double p_bc = solved_boundary_pressure(cfg, cfg.bc_z_north_val);
        const double dp_dz = (p_bc - p(i, j)) / (0.5 * d_z);
        u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
    }
    return rho_cp_face * h_face * u_z * face_length;
}

void add_axial_boundary_convection(LinearSystem& sys, const Field& p, const Field& h,
                                   const Field& theta, const Mesh& mesh,
                                   const Fields& fields, const SimulationConfig& cfg) {
    // RESERVOIR: fluid entering through the boundary carries temperature_reference
    //   (an actual fed inlet/outlet to a reservoir).
    // OPEN: the entering fluid is the same recirculating film oil, so the boundary
    //   face is zero-gradient for temperature regardless of flow direction. The
    //   non-conservative correction below already accounts for div(F) at the cell,
    //   so the open boundary face contributes no net source — no spurious heating
    //   or cooling. flux_south is the +z-flux at the south face (>0 means inflow at
    //   the south end); flux_north is the +z-flux at the north face (>0 means outflow).
    const bool south_reservoir = cfg.bc_z_south_thermal == ThermalInflowMode::CONSTANT;
    const bool north_reservoir = cfg.bc_z_north_thermal == ThermalInflowMode::CONSTANT;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double flux_south = south_boundary_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_south >= 0.0 && south_reservoir) {
            sys.source(i, 0) += flux_south * cfg.temperature_reference;  // fed inflow
        } else {
            sys.a_p(i, 0) += -flux_south;  // outflow, or open zero-gradient inflow
        }

        const int j = mesh.n_z_local - 1;
        const double flux_north = north_boundary_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_north < 0.0 && north_reservoir) {
            sys.source(i, j) += -flux_north * cfg.temperature_reference;  // fed inflow
        } else {
            sys.a_p(i, j) += flux_north;  // outflow, or open zero-gradient inflow
        }
    }
}

// Actual fed inlets (feeds_fresh_oil) pump fresh oil in at a fixed supply
// temperature, so their cells are pinned to inlet.t_supply via the penalty method.
// Open large-clearance regions (feeds_fresh_oil == false) impose no thermal
// condition and are left thermally transparent to the energy solve.
void apply_inlet_thermal_conditions(LinearSystem& sys, const Mesh& mesh,
                                    const SimulationConfig& cfg) {
    if (cfg.inlets.empty()) return;

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;

    for (const auto& inlet : cfg.inlets) {
        if (!inlet.feeds_fresh_oil) continue;  // open region: thermally transparent
        const double theta_in_rad = inlet.theta * M_PI / 180.0;
        const double size_rad = inlet.size * M_PI / 180.0;

        for (int i = 0; i < n_theta; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - theta_in_rad);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            for (int j = 0; j < n_z; ++j) {
                const double z_c = (j + 0.5) * d_z;
                bool inside = false;
                if (inlet.type == InletConfig::Type::CIRCULAR) {
                    const double dist_sq = (R * dth) * (R * dth) + (z_c - inlet.z) * (z_c - inlet.z);
                    if (dist_sq <= inlet.size * inlet.size) inside = true;
                } else {  // GROOVE
                    if (dth <= 0.5 * size_rad) inside = true;
                }

                if (inside) {
                    const double penalty = 1e20;
                    sys.a_p(i, j) = penalty;
                    sys.source(i, j) = penalty * inlet.t_supply;
                }
            }
        }
    }
}

void add_nonconservative_convection_correction(LinearSystem& sys,
                                               const Field& flux_theta,
                                               const Field& flux_z,
                                               const Field& p,
                                               const Field& h,
                                               const Field& theta,
                                               const Fields& fields,
                                               const Mesh& mesh,
                                               const SimulationConfig& cfg) {
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double flux_e = flux_theta(i, j);
            const double flux_w = flux_theta(i - 1, j);
            const double flux_n = (j + 1 < mesh.n_z_local)
                ? flux_z(i, j)
                : north_boundary_flux(p, h, theta, fields, mesh, cfg, i);
            const double flux_s = (j > 0)
                ? flux_z(i, j - 1)
                : south_boundary_flux(p, h, theta, fields, mesh, cfg, i);
            const double flux_divergence = flux_e - flux_w + flux_n - flux_s;
            sys.a_p(i, j) -= flux_divergence;
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
            conductivity(i, j) = thermal_conductivity(fields, cfg, i, j) * h_val;
            if (i >= 0 && i < n_theta) {
                capacity(i, j) = volumetric_heat_capacity(fields, cfg, i, j) * h_val;
            }
        }
    }

    Field flux_theta("thermal_flux_theta", mesh);
    Field flux_z("thermal_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i + 1, j)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i + 1, j) >= 1.0) ? 1.0 : 0.0;
            const double dp_dtheta = (p(i + 1, j) - p(i, j)) / d_theta;
            const double mu_face = 0.5 * (
                dynamic_viscosity(fields, cfg, i, j) +
                dynamic_viscosity(fields, cfg, i + 1, j));
            const double rho_cp_face = 0.5 * (
                volumetric_heat_capacity(fields, cfg, i, j) +
                volumetric_heat_capacity(fields, cfg, i + 1, j));
            const double u_theta = 0.5 * cfg.omega * R
                - g_face * h_face * h_face * dp_dtheta / (12.0 * mu_face * R);
            flux_theta(i, j) = rho_cp_face * h_face * u_theta * d_z;
        }
    }

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) flux_z(i, j) = 0.0;
        for (int j = 0; j + 1 < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i, j + 1)), cfg.min_film_thickness);
            const double g_face = (theta(i, j) >= 1.0 && theta(i, j + 1) >= 1.0) ? 1.0 : 0.0;
            const double dp_dz = (p(i, j + 1) - p(i, j)) / d_z;
            const double mu_face = 0.5 * (
                dynamic_viscosity(fields, cfg, i, j) +
                dynamic_viscosity(fields, cfg, i, j + 1));
            const double rho_cp_face = 0.5 * (
                volumetric_heat_capacity(fields, cfg, i, j) +
                volumetric_heat_capacity(fields, cfg, i, j + 1));
            const double u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
            flux_z(i, j) = rho_cp_face * h_face * u_z * R * d_theta;
        }
    }

    sys.reset();
    FVM::laplacian(sys, conductivity, mesh);
    // TYPE_DIFFERENCING is theta-only (validated); fall back to upwind defensively.
    const ConvectionScheme thermal_scheme =
        cfg.thermal_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING
            ? ConvectionScheme::UPWIND : cfg.thermal_convection_scheme;
    FVM::divergence(sys, flux_theta, flux_z, temperature, thermal_scheme, mesh);
    if (cfg.solution_mode != SolutionMode::STEADY_STATE) {
        FVM::ddt_weighted(sys, temperature, capacity, cfg.dt, mesh);
    }
    FVM::add_source(sys, heat, mesh);
    add_wall_heat_loss(sys, mesh, cfg);
    add_axial_boundary_convection(sys, p, h, theta, mesh, fields, cfg);
    add_nonconservative_convection_correction(sys, flux_theta, flux_z, p, h, theta, fields, mesh, cfg);
    apply_inlet_thermal_conditions(sys, mesh, cfg);

    sys.solve(temperature, cfg.linear_rtol);
}

}  // namespace Energy
