#include "gas_transport.hpp"

#include <algorithm>
#include <cmath>

#include "full_film.hpp"
#include "fbns.hpp"
#include "communicator.hpp"
#include "fvm.hpp"
#include "reynolds_coefficients.hpp"

namespace GasTransport {
namespace {

struct ReservoirRow {
    int i, j;
    double a_p, a_e, a_w, a_n, a_s, source;
};

// Impose reservoir values during the solve so downstream cells see the boundary
// composition. Retain the displaced equations to measure the reservoir source.
std::vector<ReservoirRow> impose_reservoir(LinearSystem& sys, const Fields& fields,
                                         const Mesh& mesh, double value) {
    std::vector<ReservoirRow> rows;
    if (!fields.has("inlet_indicator")) return rows;
    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            if (fields["inlet_indicator"](i, j) <= 0.5) continue;
            rows.push_back({i, j, sys.a_p(i, j), sys.a_e(i, j), sys.a_w(i, j),
                            sys.a_n(i, j), sys.a_s(i, j), sys.source(i, j)});
            sys.a_e(i, j) = sys.a_w(i, j) = sys.a_n(i, j) = sys.a_s(i, j) = 0;
            sys.a_p(i, j) = std::max(std::abs(sys.a_p(i, j)), 1e-30);
            sys.source(i, j) = sys.a_p(i, j) * value;
        }
    }
    return rows;
}

void account_reservoir(const std::vector<ReservoirRow>& rows, Field& concentration,
                       Field* gas_clamp, const Mesh& mesh, double dt) {
    Communicator comm(mesh);
    comm.update_ghosts(concentration);
    if (!gas_clamp) return;
    for (const auto& row : rows) {
        const int i = row.i, j = row.j;
        const double source = row.a_p * concentration(i, j) - row.a_e * concentration(i + 1, j)
            - row.a_w * concentration(i - 1, j) - row.a_n * concentration(i, j + 1)
            - row.a_s * concentration(i, j - 1) - row.source;
        (*gas_clamp)(i, j) += source * dt;
    }
}

double field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return std::max(field_value(fields, "mu", cfg.mu, i, j), 1.0e-30);
}

bool uses_pressure_boundary(BCType type) {
    return type == BCType::DIRICHLET || type == BCType::INLET_OUTLET;
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
    const double mu_face = dynamic_viscosity(fields, cfg, i, 0);
    // Cavitated boundary cell under consistent_boundary_flux: the Elrod solve
    // re-floods it at the physical reformation rate, so convect that same inflow.
    if (cfg.is_reformation_boundary(theta(i, 0), FullFilm::bound(fields, cfg, i, 0))) {
        return cfg.reformation_boundary_velocity(cfg.bc_z_south_val, h_face, mu_face,
                                                 0.5 * mesh.get_d_z(), true);
    }
    const double p_bc = cfg.cavitation_model == CavitationModel::JFO
        ? cfg.elrod_boundary_pressure(cfg.bc_z_south_val)
        : cfg.bc_z_south_val;
    const double dp_dz = (p(i, 0) - p_bc) / (0.5 * mesh.get_d_z());
    return ReynoldsCoeff::poiseuille_face_velocity(
        fields, cfg, FullFilm::gate(fields, cfg, i, 0), dp_dz, h_face, mu_face, 1.0,
        i, 0, i, 0);
}

double north_boundary_velocity(const Field& p, const Field& h, const Field& theta,
                               const Fields& fields, const Mesh& mesh,
                               const SimulationConfig& cfg, int i) {
    if (!uses_pressure_boundary(cfg.bc_z_north_type)) return 0.0;
    const int j = mesh.n_z_local - 1;
    const double h_face = std::max(h(i, j), cfg.min_film_thickness);
    const double mu_face = dynamic_viscosity(fields, cfg, i, j);
    if (cfg.is_reformation_boundary(theta(i, j), FullFilm::bound(fields, cfg, i, j))) {
        return cfg.reformation_boundary_velocity(cfg.bc_z_north_val, h_face, mu_face,
                                                 0.5 * mesh.get_d_z(), false);
    }
    const double p_bc = cfg.cavitation_model == CavitationModel::JFO
        ? cfg.elrod_boundary_pressure(cfg.bc_z_north_val)
        : cfg.bc_z_north_val;
    const double dp_dz = (p_bc - p(i, j)) / (0.5 * mesh.get_d_z());
    return ReynoldsCoeff::poiseuille_face_velocity(
        fields, cfg, FullFilm::gate(fields, cfg, i, j), dp_dz, h_face, mu_face, 1.0,
        i, j, i, j);
}

// Density convected through an axial boundary face. Full-film faces carry the
// Elrod density rho = rho_l*theta; a cavitated reflood instead brings full-density
// reservoir liquid rho_l, matching the Elrod liquid reflood mass exactly (the bulk
// modulus cancels, leaving rho_l*h^3/(12 mu) * dp/dz).
double boundary_liquid_density(const Fields& fields, const SimulationConfig& cfg,
                               const Field& theta, int i, int j) {
    const char* name = cfg.is_reformation_boundary(theta(i, j), FullFilm::bound(fields, cfg, i, j)) ? "rho_liquid_solution" : "rho";
    return std::max(field_value(fields, name, cfg.rho, i, j), 1.0e-30);
}

double south_boundary_mass_flux(const Field& p, const Field& h, const Field& theta,
                                const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, int i) {
    if (FBNS::active(cfg)) return FBNS::mass_flux_z(fields, mesh, cfg, i, -1);
    const double u_z = south_boundary_velocity(p, h, theta, fields, mesh, cfg, i);
    const double h_face = std::max(h(i, 0), cfg.min_film_thickness);
    const double rho_face = boundary_liquid_density(fields, cfg, theta, i, 0);
    return rho_face * h_face * u_z * mesh.R * mesh.get_d_theta();
}

double north_boundary_mass_flux(const Field& p, const Field& h, const Field& theta,
                                const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, int i) {
    const int j = mesh.n_z_local - 1;
    if (FBNS::active(cfg)) return FBNS::mass_flux_z(fields, mesh, cfg, i, j);
    const double u_z = north_boundary_velocity(p, h, theta, fields, mesh, cfg, i);
    const double h_face = std::max(h(i, j), cfg.min_film_thickness);
    const double rho_face = boundary_liquid_density(fields, cfg, theta, i, j);
    return rho_face * h_face * u_z * mesh.R * mesh.get_d_theta();
}

double south_boundary_volume_flux(const Field& p, const Field& h, const Field& theta,
                                  const Fields& fields, const Mesh& mesh,
                                  const SimulationConfig& cfg, int i) {
    if (FBNS::active(cfg)) return FBNS::mass_flux_z(fields, mesh, cfg, i, -1) /
        std::max(ReynoldsCoeff::liquid_mass_area(fields, cfg, i, 0), 1e-30);
    return south_boundary_velocity(p, h, theta, fields, mesh, cfg, i) *
        mesh.R * mesh.get_d_theta();
}

double north_boundary_volume_flux(const Field& p, const Field& h, const Field& theta,
                                  const Fields& fields, const Mesh& mesh,
                                  const SimulationConfig& cfg, int i) {
    if (FBNS::active(cfg)) return FBNS::mass_flux_z(fields, mesh, cfg, i, mesh.n_z_local - 1) /
        std::max(ReynoldsCoeff::liquid_mass_area(fields, cfg, i, mesh.n_z_local - 1), 1e-30);
    return north_boundary_velocity(p, h, theta, fields, mesh, cfg, i) *
        mesh.R * mesh.get_d_theta();
}

void add_dissolved_axial_boundary_convection(LinearSystem& sys, const Field& p,
                                             const Field& h, const Field& theta,
                                             const Fields& fields,
                                             const Mesh& mesh,
                                             const SimulationConfig& cfg) {
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double flux_s = south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_s >= 0.0) sys.source(i, 0) += flux_s * cfg.axial_inflow_dissolved_gas(true);
        else sys.a_p(i, 0) += -flux_s;                         // outflow

        const int j = mesh.n_z_local - 1;
        const double flux_n = north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_n < 0.0) sys.source(i, j) += -flux_n * cfg.axial_inflow_dissolved_gas(false);
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
        if (flux_s < 0.0) sys.a_p(i, 0) += -flux_s;
        else sys.source(i, 0) += south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i) *
                                    cfg.axial_inflow_free_gas_ratio(true);

        const int j = mesh.n_z_local - 1;
        const double flux_n = north_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_n > 0.0) sys.a_p(i, j) += flux_n;
        else sys.source(i, j) -= north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i) *
                                    cfg.axial_inflow_free_gas_ratio(false);
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

double axial_gas_inflow(const Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.fluid_property_model != FluidPropertyModel::TWO_PHASE || !fields.has("dissolved_gas")) return 0;
    double inflow = 0;
    if (fields.has("gas_boundary_mass_flux")) {
        for (int j = 0; j < mesh.n_z_local; ++j)
            for (int i = 0; i < mesh.n_theta_local; ++i) inflow += fields["gas_boundary_mass_flux"](i, j);
        return inflow;
    }
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const Field& c = fields["dissolved_gas"];
    const double c_in = std::clamp(cfg.dissolved_gas_initial, 0.0, cfg.dissolved_gas_max);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double fs = south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        inflow += fs * (fs >= 0 ? c_in : c(i, 0));
        if (fs < 0) inflow += south_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i) *
            field_value(fields, "free_gas_mass", 0, i, 0);
        const int j = mesh.n_z_local - 1;
        const double fn = north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        inflow -= fn * (fn < 0 ? c_in : c(i, j));
        if (fn > 0) inflow -= north_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i) *
            field_value(fields, "free_gas_mass", 0, i, j);
    }
    return inflow;
}

void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh,
           const SimulationConfig& cfg, double dt) {
    if (cfg.fluid_property_model != FluidPropertyModel::TWO_PHASE) return;
    if (!fields.has("dissolved_gas")) return;
    if (dt <= 0.0) return;

    Field& c = fields["dissolved_gas"];
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    Field* gas_clamp = fields.has("gas_clamp_mass") ? &fields["gas_clamp_mass"] : nullptr;
    if (gas_clamp != nullptr) gas_clamp->fill(0.0);
    const double cell_area = mesh.cell_volume();
    Field& saved_boundary = fields.has("gas_boundary_mass_flux") ? fields["gas_boundary_mass_flux"]
        : fields.add("gas_boundary_mass_flux", mesh);
    saved_boundary.fill(0);

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
            const double g_face = FullFilm::face_gate(fields, cfg, i, j, i + 1, j);
            const double dp_dtheta = (p(i + 1, j) - p(i, j)) / d_theta;
            const double mu_face = 0.5 * (dynamic_viscosity(fields, cfg, i, j) +
                                          dynamic_viscosity(fields, cfg, i + 1, j));
            const double rho_face = 0.5 * (field_value(fields, "rho", cfg.rho, i, j) +
                                           field_value(fields, "rho", cfg.rho, i + 1, j));
            // Both terms come from the same header the Reynolds and FBNS
            // assemblies use, so gas rides on the velocity field that produced the
            // pressure rather than on a second, classically-approximated one.
            const double u_theta =
                ReynoldsCoeff::couette_face_velocity(fields, cfg, R, h_face,
                                                     i, j, i + 1, j) +
                ReynoldsCoeff::poiseuille_face_velocity(fields, cfg, g_face, dp_dtheta,
                                                        h_face, mu_face, R,
                                                        i, j, i + 1, j);
            flux_theta(i, j) = rho_face * h_face * u_theta * d_z;
            if (FBNS::active(cfg)) flux_theta(i, j) = FBNS::mass_flux_theta(fields, mesh, cfg, i, j);
            free_gas_flux_theta(i, j) = u_theta * d_z;
            if (FBNS::active(cfg)) {
                const int upwind = flux_theta(i, j) >= 0 ? i : i + 1;
                const double rho_l_face = 0.5 * (ReynoldsCoeff::base_density(fields, cfg, i, j) +
                                                  ReynoldsCoeff::base_density(fields, cfg, i + 1, j));
                free_gas_flux_theta(i, j) = flux_theta(i, j) /
                    std::max(rho_l_face * h_face * theta(upwind, j), 1e-30);
            }
        }
    }
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            flux_z(i, j) = 0.0;
            free_gas_flux_z(i, j) = 0.0;
        }
        for (int j = 0; j + 1 < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i, j + 1)), cfg.min_film_thickness);
            const double g_face = FullFilm::face_gate(fields, cfg, i, j, i, j + 1);
            const double dp_dz = (p(i, j + 1) - p(i, j)) / d_z;
            const double mu_face = 0.5 * (dynamic_viscosity(fields, cfg, i, j) +
                                          dynamic_viscosity(fields, cfg, i, j + 1));
            const double rho_face = 0.5 * (field_value(fields, "rho", cfg.rho, i, j) +
                                           field_value(fields, "rho", cfg.rho, i, j + 1));
            const double u_z =
                ReynoldsCoeff::poiseuille_face_velocity(fields, cfg, g_face, dp_dz,
                                                        h_face, mu_face, 1.0,
                                                        i, j, i, j + 1);
            flux_z(i, j) = rho_face * h_face * u_z * R * d_theta;
            if (FBNS::active(cfg)) flux_z(i, j) = FBNS::mass_flux_z(fields, mesh, cfg, i, j);
            free_gas_flux_z(i, j) = u_z * R * d_theta;
            if (FBNS::active(cfg)) free_gas_flux_z(i, j) = flux_z(i, j) /
                std::max(rho_face * h_face, 1e-30);
        }
    }

    sys.reset();
    if (D_eff > 0.0) FVM::laplacian(sys, diffusion, mesh);
    // TYPE_DIFFERENCING is theta-only (validated); fall back to upwind defensively.
    const ConvectionScheme gas_scheme =
        cfg.gas_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING
            ? ConvectionScheme::UPWIND : cfg.gas_convection_scheme;
    // See the note in Energy::solve: phi.old on the live field lags a whole
    // timestep, so the deferred correction needs the current iterate instead.
    Field c_lag("dissolved_gas_lag", mesh);
    c_lag.data = c.data;
    c_lag.store_old_time();
    FVM::divergence(sys, flux_theta, flux_z, c_lag, gas_scheme, mesh);
    add_dissolved_axial_boundary_convection(sys, p, h, theta, fields, mesh, cfg);
    if (cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0) {
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                sys.a_p(i, j) += capacity(i, j) * cell_area / dt;
                sys.source(i, j) += ReynoldsCoeff::old_liquid_mass_area(fields, cfg, dt, i, j) *
                    c.old(i, j) * cell_area / dt;
            }
        }
    }

    const double c_in = std::clamp(cfg.dissolved_gas_initial, 0.0, cfg.dissolved_gas_max);
    const auto dissolved_reservoir = impose_reservoir(sys, fields, mesh, c_in);
    sys.solve(c, cfg.linear_rtol, 1000, "dissolved gas concentration");
    account_reservoir(dissolved_reservoir, c, gas_clamp, mesh, dt);
    for (int i = 0; i < n_theta; ++i) {
        const double fs = south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        const int j = n_z - 1;
        const double fn = north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i);
        saved_boundary(i, 0) += fs * (fs >= 0 ? cfg.axial_inflow_dissolved_gas(true) : c(i, 0));
        saved_boundary(i, j) -= fn * (fn < 0 ? cfg.axial_inflow_dissolved_gas(false) : c(i, j));
    }

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
    if (fields.has("free_gas_mass")) {
        Field& free_gas_mass = fields["free_gas_mass"];
        sys.reset();
        FVM::divergence(sys, free_gas_flux_theta, free_gas_flux_z, free_gas_mass,
                        gas_scheme, mesh);
        add_free_gas_axial_boundary_convection(sys, p, h, theta, fields, mesh, cfg);
        FVM::ddt(sys, free_gas_mass, dt, mesh);
        const auto free_reservoir = impose_reservoir(sys, fields, mesh, 0);
        sys.solve(free_gas_mass, cfg.linear_rtol, 1000, "free gas mass");
        account_reservoir(free_reservoir, free_gas_mass, gas_clamp, mesh, dt);

        for (int i = 0; i < n_theta; ++i) {
            const double fs = south_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i);
            const double fn = north_boundary_volume_flux(p, h, theta, fields, mesh, cfg, i);
            if (fs < 0) saved_boundary(i, 0) += fs * free_gas_mass(i, 0);
            else saved_boundary(i, 0) += south_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i) *
                                            cfg.axial_inflow_free_gas_ratio(true);
            if (fn > 0) saved_boundary(i, n_z - 1) -= fn * free_gas_mass(i, n_z - 1);
            else saved_boundary(i, n_z - 1) -= north_boundary_mass_flux(p, h, theta, fields, mesh, cfg, i) *
                                                 cfg.axial_inflow_free_gas_ratio(false);
        }

        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double m_clamped = std::max(0.0, free_gas_mass(i, j));
                if (gas_clamp != nullptr && m_clamped != free_gas_mass(i, j)) {
                    (*gas_clamp)(i, j) += (m_clamped - free_gas_mass(i, j)) * cell_area;
                }
                free_gas_mass(i, j) = m_clamped;
            }
        }
    }
    apply_supply_cell_gas_composition(fields, mesh, cfg);
}

}  // namespace GasTransport
