#include "diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include <mpi.h>

#include "equation_of_state.hpp"

namespace Diagnostics {
namespace {

double field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double positive_field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return std::max(field_value(fields, name, fallback, i, j), 1.0e-30);
}

double base_density(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
}

double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "mu", cfg.mu, i, j);
}

double harmonic_avg(double a, double b) {
    if (a + b == 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

double gamma_base(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    const double h_val = fields["h"](i, j);
    return base_density(fields, cfg, i, j) * cfg.bulk_modulus *
        h_val * h_val * h_val / (12.0 * dynamic_viscosity(fields, cfg, i, j));
}

// Geometric inlet membership, mirroring Reynolds::apply_inlet_conditions.
bool inlet_cell(const Mesh& mesh, const SimulationConfig& cfg, int i, int j) {
    if (cfg.inlets.empty()) return false;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
    const double z_c = (j + 0.5) * d_z;
    for (const auto& inlet : cfg.inlets) {
        const double theta_in_rad = inlet.theta * M_PI / 180.0;
        const double size_rad = inlet.size * M_PI / 180.0;
        double dth = std::abs(theta_c - theta_in_rad);
        if (dth > M_PI) dth = 2.0 * M_PI - dth;
        if (inlet.type == InletConfig::Type::CIRCULAR) {
            const double dist_sq = (mesh.R * dth) * (mesh.R * dth) +
                (z_c - inlet.z) * (z_c - inlet.z);
            if (dist_sq <= inlet.size * inlet.size) return true;
        } else {
            if (dth <= 0.5 * size_rad) return true;
        }
    }
    return false;
}

// Couette mass face flux at the east face of cell (i,j), mirroring solve_elrod.
double couette_face_flux(const Fields& fields, const Mesh& mesh,
                         const SimulationConfig& cfg, int i, int j) {
    const Field& h = fields["h"];
    const double h_e = 0.5 * (h(i, j) + h(i + 1, j));
    const double rho_e = 0.5 * (base_density(fields, cfg, i, j) +
                                base_density(fields, cfg, i + 1, j));
    return rho_e * cfg.omega * mesh.R * h_e * 0.5 * mesh.get_d_z();
}

// Face value of theta carried by the Couette flux (upwind; central on full-film
// faces under TYPE_DIFFERENCING, matching the assembled deferred correction).
double couette_face_theta(const Field& theta, const SimulationConfig& cfg,
                          double flux, int i, int j) {
    const double theta_up = (flux >= 0.0) ? theta(i, j) : theta(i + 1, j);
    if (cfg.theta_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING &&
        theta(i, j) >= 1.0 && theta(i + 1, j) >= 1.0) {
        return 0.5 * (theta(i, j) + theta(i + 1, j));
    }
    return theta_up;
}

// Net liquid mass outflow rate [kg/s] of one Elrod inlet penalty cell,
// re-evaluating the physical operators the penalty row replaced.
double elrod_inlet_cell_outflow(const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, double dt, int i, int j) {
    const Field& theta = fields["theta"];
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double ew_scale = d_z / (R * d_theta);
    const double ns_scale = R * d_theta / d_z;

    auto gamma_g = [&](int ii, int jj) {
        return gamma_base(fields, cfg, ii, jj) * EOS::switch_function(theta(ii, jj));
    };

    double outflow = 0.0;

    // Poiseuille diffusion through the four faces (z-boundary faces use the
    // Dirichlet half-cell formula the assembly adds; Neumann faces carry none).
    outflow += harmonic_avg(gamma_g(i, j), gamma_g(i + 1, j)) * ew_scale *
        (theta(i, j) - theta(i + 1, j));
    outflow += harmonic_avg(gamma_g(i - 1, j), gamma_g(i, j)) * ew_scale *
        (theta(i, j) - theta(i - 1, j));
    if (j + 1 < n_z) {
        outflow += harmonic_avg(gamma_g(i, j), gamma_g(i, j + 1)) * ns_scale *
            (theta(i, j) - theta(i, j + 1));
    } else if (cfg.bc_z_north_type == BCType::DIRICHLET) {
        outflow += cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, j), theta(i, j),
                                              cfg.bc_z_north_val, R, d_theta, d_z);
    }
    if (j > 0) {
        outflow += harmonic_avg(gamma_g(i, j - 1), gamma_g(i, j)) * ns_scale *
            (theta(i, j) - theta(i, j - 1));
    } else if (cfg.bc_z_south_type == BCType::DIRICHLET) {
        outflow += cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, j), theta(i, j),
                                              cfg.bc_z_south_val, R, d_theta, d_z);
    }

    // Couette convection (theta direction only).
    const double flux_e = couette_face_flux(fields, mesh, cfg, i, j);
    const double flux_w = couette_face_flux(fields, mesh, cfg, i - 1, j);
    outflow += flux_e * couette_face_theta(theta, cfg, flux_e, i, j);
    outflow -= flux_w * couette_face_theta(theta, cfg, flux_w, i - 1, j);

    // Liquid-content storage plus the explicit squeeze term, mirroring assembly.
    if (cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0) {
        const double rho_l = base_density(fields, cfg, i, j);
        outflow += rho_l * fields["h"](i, j) * mesh.cell_volume() *
            (theta(i, j) - theta.old(i, j)) / dt;
        if (fields.has("dh_dt")) {
            outflow += rho_l * theta.old(i, j) * fields["dh_dt"](i, j) * mesh.cell_volume();
        }
    }
    return outflow;
}

// Net liquid inflow [kg/s] through the axial boundaries, mirroring the
// boundary terms the active cavitation model assembles.
double liquid_boundary_inflow(const Fields& fields, const Mesh& mesh,
                              const SimulationConfig& cfg) {
    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const Field& pressure = fields["pressure"];

    double inflow = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS) {
            const Field& theta = fields["theta"];
            const bool south_active = cfg.bc_z_south_type == BCType::DIRICHLET ||
                (cfg.bc_z_south_type == BCType::INLET_OUTLET &&
                 pressure(i, 0) < cfg.bc_z_south_val);
            if (south_active) {
                inflow -= cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, 0), theta(i, 0),
                                                     cfg.bc_z_south_val, R, d_theta, d_z);
            }
            const bool north_active = cfg.bc_z_north_type == BCType::DIRICHLET ||
                (cfg.bc_z_north_type == BCType::INLET_OUTLET &&
                 pressure(i, n_z - 1) < cfg.bc_z_north_val);
            if (north_active) {
                inflow -= cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, n_z - 1), theta(i, n_z - 1),
                                                     cfg.bc_z_north_val, R, d_theta, d_z);
            }
        } else {
            const Field& h = fields["h"];
            auto gamma_p = [&](int ii, int jj) {
                const double h_val = h(ii, jj);
                return field_value(fields, "rho", cfg.rho, ii, jj) *
                    h_val * h_val * h_val / (12.0 * dynamic_viscosity(fields, cfg, ii, jj));
            };
            if (cfg.bc_z_south_type != BCType::NEUMANN) {
                const double cs = gamma_p(i, 0) * R * d_theta / (0.5 * d_z);
                inflow += cs * (cfg.bc_z_south_val - pressure(i, 0));
            }
            if (cfg.bc_z_north_type != BCType::NEUMANN) {
                const double cn = gamma_p(i, n_z - 1) * R * d_theta / (0.5 * d_z);
                inflow += cn * (cfg.bc_z_north_val - pressure(i, n_z - 1));
            }
        }
    }
    return inflow;
}

// Net total-gas inflow [kg/s] through the axial boundaries, mirroring the
// gas-transport boundary contracts (inflow carries the configured reservoir
// composition with zero free gas; outflow convects local values out).
double gas_boundary_inflow(const Fields& fields, const Mesh& mesh,
                           const SimulationConfig& cfg) {
    if (cfg.fluid_property_model != FluidPropertyModel::GAS_CAVITATION_MIXTURE) return 0.0;
    if (!fields.has("dissolved_gas")) return 0.0;

    const Field& pressure = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const Field& c = fields["dissolved_gas"];
    const int n_z = mesh.n_z_local;
    const double d_z = mesh.get_d_z();
    const double c_in = std::clamp(cfg.dissolved_gas_initial, 0.0, cfg.dissolved_gas_max);

    auto pressure_boundary = [](BCType type) {
        return type == BCType::DIRICHLET || type == BCType::INLET_OUTLET;
    };

    double inflow = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        // South face: u_z > 0 is inflow.
        if (pressure_boundary(cfg.bc_z_south_type)) {
            const double h_face = std::max(h(i, 0), cfg.min_film_thickness);
            const double mu_face = dynamic_viscosity(fields, cfg, i, 0);
            double u_z;
            double rho_face;
            if (cfg.is_reformation_boundary(theta(i, 0))) {
                // Mirror the gas-transport reflood: full-density reservoir liquid.
                u_z = cfg.reformation_boundary_velocity(cfg.bc_z_south_val, h_face, mu_face,
                                                        0.5 * d_z, true);
                rho_face = std::max(field_value(fields, "rho_liquid_solution", cfg.rho, i, 0), 1.0e-30);
            } else {
                const double p_bc = cfg.elrod_boundary_pressure(cfg.bc_z_south_val);
                const double dp_dz = (pressure(i, 0) - p_bc) / (0.5 * d_z);
                const double g_face = theta(i, 0) >= 1.0 ? 1.0 : 0.0;
                u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
                rho_face = field_value(fields, "rho", cfg.rho, i, 0);
            }
            const double volume_flux = u_z * mesh.R * mesh.get_d_theta();
            const double mass_flux = rho_face * h_face * volume_flux;
            if (mass_flux >= 0.0) {
                inflow += mass_flux * c_in;  // reservoir composition, zero free gas
            } else {
                inflow += mass_flux * c(i, 0) +
                    volume_flux * field_value(fields, "free_gas_mass", 0.0, i, 0);
            }
        }
        // North face: u_z > 0 is outflow.
        if (pressure_boundary(cfg.bc_z_north_type)) {
            const int j = n_z - 1;
            const double h_face = std::max(h(i, j), cfg.min_film_thickness);
            const double mu_face = dynamic_viscosity(fields, cfg, i, j);
            double u_z;
            double rho_face;
            if (cfg.is_reformation_boundary(theta(i, j))) {
                u_z = cfg.reformation_boundary_velocity(cfg.bc_z_north_val, h_face, mu_face,
                                                        0.5 * d_z, false);
                rho_face = std::max(field_value(fields, "rho_liquid_solution", cfg.rho, i, j), 1.0e-30);
            } else {
                const double p_bc = cfg.elrod_boundary_pressure(cfg.bc_z_north_val);
                const double dp_dz = (p_bc - pressure(i, j)) / (0.5 * d_z);
                const double g_face = theta(i, j) >= 1.0 ? 1.0 : 0.0;
                u_z = -g_face * h_face * h_face * dp_dz / (12.0 * mu_face);
                rho_face = field_value(fields, "rho", cfg.rho, i, j);
            }
            const double volume_flux = u_z * mesh.R * mesh.get_d_theta();
            const double mass_flux = rho_face * h_face * volume_flux;
            if (mass_flux < 0.0) {
                inflow += -mass_flux * c_in;
            } else {
                inflow -= mass_flux * c(i, j) +
                    volume_flux * field_value(fields, "free_gas_mass", 0.0, i, j);
            }
        }
    }
    return inflow;
}

}  // namespace

MassBalance mass_balance(const Fields& fields, const Mesh& mesh,
                         const SimulationConfig& cfg, double dt) {
    MassBalance balance;
    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double cell_area = mesh.cell_volume();
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const bool transient = cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0;
    const bool gas_active =
        cfg.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE &&
        fields.has("dissolved_gas") && fields.has("free_gas_mass");

    // Storage terms mirror the assembled capacity (current rho_l h times the
    // theta change) plus the explicit squeeze source when dynamic h is present.
    double liquid_storage_rate = 0.0;  // kg/s
    double gas_storage_rate = 0.0;
    double cavitated_cells = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double rho_l = base_density(fields, cfg, i, j);
            const double rho_h = rho_l * theta(i, j) * h(i, j);
            balance.liquid_mass += rho_h * cell_area;
            if (theta(i, j) < 1.0) cavitated_cells += 1.0;

            if (transient) {
                liquid_storage_rate += rho_l * h(i, j) * cell_area *
                    (theta(i, j) - theta.old(i, j)) / dt;
                if (fields.has("dh_dt")) {
                    liquid_storage_rate += rho_l * theta.old(i, j) *
                        fields["dh_dt"](i, j) * cell_area;
                }
            }

            if (gas_active) {
                const double c_val = fields["dissolved_gas"](i, j);
                const double m_g = fields["free_gas_mass"](i, j);
                balance.gas_mass += (c_val * rho_h + m_g) * cell_area;
                if (transient) {
                    gas_storage_rate += (rho_l * theta(i, j) * h(i, j) *
                                             (c_val - fields["dissolved_gas"].old(i, j)) +
                                         (m_g - fields["free_gas_mass"].old(i, j))) *
                        cell_area / dt;
                }
            }

            balance.clamp_mass_liquid += field_value(fields, "cavitation_clamp_mass", 0.0, i, j);
            balance.clamp_mass_gas += field_value(fields, "gas_clamp_mass", 0.0, i, j);

            if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS &&
                inlet_cell(mesh, cfg, i, j)) {
                balance.inlet_mass_source +=
                    elrod_inlet_cell_outflow(fields, mesh, cfg, dt, i, j);
            }
        }
    }

    balance.liquid_boundary_flux = liquid_boundary_inflow(fields, mesh, cfg);
    balance.gas_boundary_flux = gas_boundary_inflow(fields, mesh, cfg);

    double sums[10] = {
        balance.liquid_mass, balance.gas_mass,
        balance.liquid_boundary_flux, balance.gas_boundary_flux,
        balance.inlet_mass_source,
        balance.clamp_mass_liquid, balance.clamp_mass_gas,
        liquid_storage_rate, gas_storage_rate, cavitated_cells};
    MPI_Allreduce(MPI_IN_PLACE, sums, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    balance.liquid_mass = sums[0];
    balance.gas_mass = sums[1];
    balance.liquid_boundary_flux = sums[2];
    balance.gas_boundary_flux = sums[3];
    balance.inlet_mass_source = sums[4];
    balance.clamp_mass_liquid = sums[5];
    balance.clamp_mass_gas = sums[6];
    liquid_storage_rate = sums[7];
    gas_storage_rate = sums[8];
    balance.cavitated_fraction =
        sums[9] / (double(mesh.n_theta_global) * double(mesh.n_z_global));

    // Step closure: storage - inflow - inlet source - clamp/dt should be the
    // linear-solver residual. Normalized by the resident mass.
    const double step = transient ? dt : std::max(cfg.dt, 1.0e-30);
    balance.liquid_residual =
        (liquid_storage_rate - balance.liquid_boundary_flux - balance.inlet_mass_source -
         balance.clamp_mass_liquid / step) * step / std::max(balance.liquid_mass, 1.0e-30);
    if (gas_active) {
        balance.gas_residual =
            (gas_storage_rate - balance.gas_boundary_flux - balance.clamp_mass_gas / step) *
            step / std::max(balance.gas_mass, 1.0e-30);
    }
    return balance;
}

void append_csv(const std::filesystem::path& path, double time, int step,
                const MassBalance& balance, const Reynolds::ElrodStats& stats) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) return;

    const bool write_header = !std::filesystem::exists(path);
    std::ofstream csv(path, std::ios::app);
    if (!csv.is_open()) return;
    if (write_header) {
        csv << "step,time,outer_iters,flag_flips,max_theta_change,converged,"
               "liquid_mass,liquid_boundary_flux,inlet_mass_source,clamp_mass_liquid,"
               "liquid_residual,gas_mass,gas_boundary_flux,clamp_mass_gas,gas_residual,"
               "cavitated_fraction\n";
    }
    csv.precision(15);
    csv << step << ',' << time << ',' << stats.outer_iters << ',' << stats.flag_flips << ','
        << stats.max_theta_change << ',' << (stats.converged ? 1 : 0) << ','
        << balance.liquid_mass << ',' << balance.liquid_boundary_flux << ','
        << balance.inlet_mass_source << ',' << balance.clamp_mass_liquid << ','
        << balance.liquid_residual << ',' << balance.gas_mass << ','
        << balance.gas_boundary_flux << ',' << balance.clamp_mass_gas << ','
        << balance.gas_residual << ',' << balance.cavitated_fraction << '\n';
}

}  // namespace Diagnostics
