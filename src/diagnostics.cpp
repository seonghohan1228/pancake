#include "diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include <mpi.h>

#include "equation_of_state.hpp"
#include "fbns.hpp"
#include "full_film.hpp"
#include "gas_transport.hpp"
#include "reynolds_coefficients.hpp"

namespace Diagnostics {
namespace {

// Every coefficient below comes from ReynoldsCoeff, the single definition the
// Elrod and FBNS assemblies use. This module *measures* those assemblies, so a
// re-derived approximation here is not a shortcut but a wrong answer: the
// previous local copies used h^3/(12 mu) and cfg.bulk_modulus where the assembly
// had already moved to the Dowson F2 and the mixture beta_bar, which made
// mass_residual report a fabricated imbalance on exactly the cross-film and
// void-coupled runs it exists to police.
using ReynoldsCoeff::base_density;
using ReynoldsCoeff::couette_face_flux;
using ReynoldsCoeff::dynamic_viscosity;
using ReynoldsCoeff::field_value;
using ReynoldsCoeff::void_params;

double harmonic_avg(double a, double b) {
    if (a + b == 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

/// Elrod diffusion coefficient before the switch function, exactly as
/// Reynolds::solve_elrod builds it: rho_0 beta_bar F2 under the generalized
/// (cross-film) Reynolds equation, rho_0 beta h^3/(12 mu) classically.
double gamma_base(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return ReynoldsCoeff::elrod_gamma_base(fields, cfg, i, j);
}

double van_leer(double r) {
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

// Face value of theta carried by the Couette flux at the east face of (i,j),
// formed exactly as FVM::divergence assembles it: implicit upwind plus the
// scheme's deferred correction. LINEAR is central on every face,
// TYPE_DIFFERENCING only on full-film faces, VANLEER reproduces the limiter
// stencil (including its use of the far-upwind cell for both flux directions).
double couette_face_theta(const Fields& fields, const SimulationConfig& cfg,
                          double flux, int i, int j) {
    const Field& theta = fields["theta"];
    const double theta_up = (flux >= 0.0) ? theta(i, j) : theta(i + 1, j);
    switch (cfg.theta_convection_scheme) {
        case ConvectionScheme::UPWIND:
            return theta_up;
        case ConvectionScheme::LINEAR:
            return 0.5 * (theta(i, j) + theta(i + 1, j));
        case ConvectionScheme::TYPE_DIFFERENCING:
            return FullFilm::is_full_film_face(fields, cfg, i, j, i + 1, j)
                ? 0.5 * (theta(i, j) + theta(i + 1, j)) : theta_up;
        case ConvectionScheme::VANLEER: {
            const double denom = theta(i, j) - theta(i - 1, j);
            const double num = theta(i + 1, j) - theta(i, j);
            double r = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
            if (flux < 0.0) r = 1.0 / (r + 1e-14);
            return theta_up + 0.5 * van_leer(r) * (theta_up - theta(i - 1, j));
        }
    }
    return theta_up;
}

// True when solve_elrod adds an axial Dirichlet link on this boundary face:
// always under DIRICHLET, and under INLET_OUTLET only while the cell sits below
// the reservoir pressure (inflow). Single predicate so the per-cell outflow and
// the domain boundary flux cannot disagree about which faces are open.
bool axial_link_active(const Fields& fields, BCType type, double bc_val, int i, int j) {
    if (type == BCType::DIRICHLET) return true;
    if (type != BCType::INLET_OUTLET) return false;
    return fields["pressure"](i, j) < bc_val;
}

// Net liquid mass outflow rate [kg/s] of one Elrod cell, evaluating the physical
// operators solve_elrod assembles. Interior faces telescope between neighbours,
// so summing this over the domain leaves only the axial boundary faces and the
// storage terms. Used both for the replaced inlet rows and for mass_residual.
double elrod_cell_outflow(const Fields& fields, const Mesh& mesh,
                          const SimulationConfig& cfg, double dt, int i, int j,
                          bool include_axial_boundary_faces = true) {
    const Field& theta = fields["theta"];
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double ew_scale = d_z / (R * d_theta);
    const double ns_scale = R * d_theta / d_z;

    auto gamma_g = [&](int ii, int jj) {
        return gamma_base(fields, cfg, ii, jj) *
            EOS::switch_function(theta(ii, jj), void_params(fields, cfg, ii, jj).theta_full);
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
    } else if (include_axial_boundary_faces &&
               axial_link_active(fields, cfg.bc_z_north_type, cfg.bc_z_north_val, i, j)) {
        outflow += cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, j), theta(i, j),
                                              FullFilm::bound(fields, cfg, i, j),
                                              cfg.bc_z_north_val, R, d_theta, d_z);
    }
    if (j > 0) {
        outflow += harmonic_avg(gamma_g(i, j - 1), gamma_g(i, j)) * ns_scale *
            (theta(i, j) - theta(i, j - 1));
    } else if (include_axial_boundary_faces &&
               axial_link_active(fields, cfg.bc_z_south_type, cfg.bc_z_south_val, i, j)) {
        outflow += cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, j), theta(i, j),
                                              FullFilm::bound(fields, cfg, i, j),
                                              cfg.bc_z_south_val, R, d_theta, d_z);
    }

    // Couette convection (theta direction only).
    const double flux_e = couette_face_flux(fields, mesh, cfg, i, j);
    const double flux_w = couette_face_flux(fields, mesh, cfg, i - 1, j);
    outflow += flux_e * couette_face_theta(fields, cfg, flux_e, i, j);
    outflow -= flux_w * couette_face_theta(fields, cfg, flux_w, i - 1, j);

    // Conservative storage, using the same old inventory as assembly.
    if (cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0) {
        outflow += (ReynoldsCoeff::liquid_mass_area(fields, cfg, i, j) -
                    ReynoldsCoeff::old_liquid_mass_area(fields, cfg, dt, i, j)) * mesh.cell_volume() / dt;
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
    if (FBNS::active(cfg)) {
        // FBNS carries the axial ends as a plain Poiseuille Dirichlet link on the
        // pressure, with no bulk-modulus conversion and no switch function, so
        // the boundary flux is read straight off the pressure drop.
        auto gamma_p = [&](int ii, int jj) {
            return ReynoldsCoeff::fbns_mobility(fields, cfg, ii, jj);
        };
        for (int i = 0; i < n_theta; ++i) {
            if (cfg.bc_z_south_type == BCType::DIRICHLET &&
                !mesh.inlet_at(cfg, i, 0)) {
                const double cs = gamma_p(i, 0) * R * d_theta / (0.5 * d_z);
                inflow += cs * (std::max(cfg.bc_z_south_val, cfg.effective_p_cav()) - pressure(i, 0));
            }
            if (cfg.bc_z_north_type == BCType::DIRICHLET &&
                !mesh.inlet_at(cfg, i, n_z - 1)) {
                const double cn = gamma_p(i, n_z - 1) * R * d_theta / (0.5 * d_z);
                inflow += cn * (std::max(cfg.bc_z_north_val, cfg.effective_p_cav()) - pressure(i, n_z - 1));
            }
        }
        return inflow;
    }
    for (int i = 0; i < n_theta; ++i) {
        if (cfg.cavitation_model == CavitationModel::JFO) {
            const Field& theta = fields["theta"];
            if (!mesh.inlet_at(cfg, i, 0) &&
                axial_link_active(fields, cfg.bc_z_south_type, cfg.bc_z_south_val, i, 0)) {
                inflow -= cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, 0), theta(i, 0),
                                                     FullFilm::bound(fields, cfg, i, 0),
                                                     cfg.bc_z_south_val, R, d_theta, d_z);
            }
            if (!mesh.inlet_at(cfg, i, n_z - 1) &&
                axial_link_active(fields, cfg.bc_z_north_type, cfg.bc_z_north_val, i, n_z - 1)) {
                inflow -= cfg.elrod_boundary_outflow(gamma_base(fields, cfg, i, n_z - 1), theta(i, n_z - 1),
                                                     FullFilm::bound(fields, cfg, i, n_z - 1),
                                                     cfg.bc_z_north_val, R, d_theta, d_z);
            }
        } else {
            // Gumbel / full-Sommerfeld: Reynolds::solve assembles rho (not rho_l)
            // times the same Poiseuille coefficient.
            auto gamma_p = [&](int ii, int jj) {
                return ReynoldsCoeff::poiseuille(
                    fields, cfg, field_value(fields, "rho", cfg.rho, ii, jj), ii, jj);
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
    return GasTransport::axial_gas_inflow(fields, mesh, cfg);
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
        cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE &&
        fields.has("dissolved_gas") && fields.has("free_gas_mass");

    // Difference the complete inventories. Phase transfer is an internal source.
    double liquid_storage_rate = 0.0;  // kg/s
    double gas_storage_rate = 0.0;
    double cavitated_cells = 0.0;
    double liquid_exchange_scale = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double rho_l = base_density(fields, cfg, i, j);
            const double rho_h = rho_l * theta(i, j) * h(i, j);
            balance.liquid_mass += rho_h * cell_area;
            if (!FullFilm::is_full_film(fields, cfg, i, j)) cavitated_cells += 1.0;

            if (transient) {
                liquid_storage_rate += (rho_h -
                    ReynoldsCoeff::old_liquid_mass_area(fields, cfg, dt, i, j)) * cell_area / dt;
                if (gas_active) liquid_storage_rate +=
                    field_value(fields, "gas_mass_transfer", 0.0, i, j) * h(i, j) * cell_area;
            }

            if (gas_active) {
                const double c_val = fields["dissolved_gas"](i, j);
                const double m_g = fields["free_gas_mass"](i, j);
                balance.gas_mass += (c_val * rho_h + m_g) * cell_area;
                if (fields.has("rho_gas") && fields["rho_gas"](i, j) > 0 && h(i, j) > 0) {
                    const double alpha_eos = m_g / (fields["rho_gas"](i, j) * h(i, j));
                    balance.gas_volume_excess = std::max(balance.gas_volume_excess,
                        theta(i, j) + alpha_eos - 1.0);
                    balance.gas_eos_defect = std::max(balance.gas_eos_defect,
                        std::abs(field_value(fields, "alpha_gas", 0, i, j) - alpha_eos));
                }
                if (transient) {
                    gas_storage_rate += (rho_h * c_val -
                                         ReynoldsCoeff::old_liquid_mass_area(fields, cfg, dt, i, j) *
                                             fields["dissolved_gas"].old(i, j) +
                                         (m_g - fields["free_gas_mass"].old(i, j))) *
                        cell_area / dt;
                }
            }

            balance.clamp_mass_liquid += field_value(fields, "cavitation_clamp_mass", 0.0, i, j);
            balance.clamp_mass_gas += field_value(fields, "gas_clamp_mass", 0.0, i, j);

            if (mesh.inlet_at(cfg, i, j)) {
                // Inlet row replacement occurs after axial-BC assembly, so an
                // overlapping groove cell no longer owns an axial boundary
                // link. Restoring that discarded face here makes a full-length
                // groove source grow as 1/dz.
                if (FBNS::active(cfg)) {
                    const double feed = FBNS::cell_outflow(fields, mesh, cfg, dt, i, j, false);
                    balance.inlet_mass_source += feed;
                    liquid_exchange_scale += std::abs(feed);
                } else if (cfg.cavitation_model == CavitationModel::JFO) {
                    balance.inlet_mass_source +=
                        elrod_cell_outflow(fields, mesh, cfg, dt, i, j, false);
                }
            }
        }
    }

    balance.liquid_boundary_flux = liquid_boundary_inflow(fields, mesh, cfg);
    balance.gas_boundary_flux = gas_boundary_inflow(fields, mesh, cfg);
    if (FBNS::active(cfg)) {
        for (int i = 0; i < n_theta; ++i) {
            if (!mesh.inlet_at(cfg, i, 0))
                liquid_exchange_scale += std::abs(FBNS::mass_flux_z(fields, mesh, cfg, i, -1));
            if (!mesh.inlet_at(cfg, i, n_z - 1))
                liquid_exchange_scale += std::abs(FBNS::mass_flux_z(fields, mesh, cfg, i, n_z - 1));
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &liquid_exchange_scale, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double sums[10] = {
        balance.liquid_mass, balance.gas_mass,
        balance.liquid_boundary_flux, balance.gas_boundary_flux,
        balance.inlet_mass_source,
        balance.clamp_mass_liquid, balance.clamp_mass_gas,
        liquid_storage_rate, gas_storage_rate, cavitated_cells};
    MPI_Allreduce(MPI_IN_PLACE, sums, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double volume_errors[2] = {balance.gas_volume_excess, balance.gas_eos_defect};
    MPI_Allreduce(MPI_IN_PLACE, volume_errors, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    balance.gas_volume_excess = volume_errors[0];
    balance.gas_eos_defect = volume_errors[1];
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
    balance.liquid_residual_kg_s = liquid_storage_rate - balance.liquid_boundary_flux -
        balance.inlet_mass_source - balance.clamp_mass_liquid / step;
    balance.liquid_balance_scale_kg_s = std::max({liquid_exchange_scale,
        std::abs(liquid_storage_rate), std::abs(balance.liquid_boundary_flux),
        std::abs(balance.inlet_mass_source), std::abs(balance.clamp_mass_liquid / step), 1.0e-30});
    balance.liquid_residual_rel = balance.liquid_residual_kg_s / balance.liquid_balance_scale_kg_s;
    if (gas_active) {
        balance.gas_residual =
            (gas_storage_rate - balance.gas_boundary_flux - balance.clamp_mass_gas / step) *
            step / std::max(balance.gas_mass, 1.0e-30);
    }
    return balance;
}

void mass_residual(Fields& fields, const Mesh& mesh,
                   const SimulationConfig& cfg, double dt) {
    if (!fields.has("mass_residual")) return;  // not selected: nothing computed
    Field& residual = fields["mass_residual"];
    residual.fill(0.0);
    if (!cfg.mass_conserving_cavitation()) return;

    // FBNS assembles a different operator from solve_elrod (pressure-driven
    // Poiseuille without the bulk modulus or the switch function), so the
    // residual mirrors that operator instead. No clamp correction is needed:
    // the Newton clamps act on iterates, not on the converged state.
    if (FBNS::active(cfg)) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            for (int j = 0; j < mesh.n_z_local; ++j) {
                residual(i, j) = mesh.inlet_at(cfg, i, j)
                    ? 0.0 : FBNS::cell_outflow(fields, mesh, cfg, dt, i, j);
                if (dt > 0 && !mesh.inlet_at(cfg, i, j) &&
                    cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE)
                    residual(i, j) += field_value(fields, "gas_mass_transfer", 0, i, j) *
                        fields["h"](i, j) * mesh.cell_volume();
            }
        }
        return;
    }

    // Same step normalization mass_balance() uses for the clamp rate, so the two
    // agree term by term.
    const bool transient = cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0;
    const double step = transient ? dt : std::max(cfg.dt, 1.0e-30);

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            // Inlet reservoir cells: their flux imbalance IS the modelled supply
            // source (mass_balance reports it as inlet_mass_source), so only the
            // clamp correction is left here.
            const double outflow = mesh.inlet_at(cfg, i, j)
                ? 0.0 : elrod_cell_outflow(fields, mesh, cfg, dt, i, j);
            residual(i, j) = outflow -
                field_value(fields, "cavitation_clamp_mass", 0.0, i, j) / step;
            if (transient && !mesh.inlet_at(cfg, i, j) &&
                cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE)
                residual(i, j) += field_value(fields, "gas_mass_transfer", 0, i, j) *
                    fields["h"](i, j) * mesh.cell_volume();
        }
    }
}

void append_csv(const std::filesystem::path& path, double time, int step,
                const MassBalance& balance, const Reynolds::ElrodStats& stats,
                const Coupling::Status* coupling) {
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
               "cavitated_fraction,"
               "coupling_iters,coupling_R_p,coupling_R_T,coupling_converged,"
               "coupling_best_iter,coupling_best_residual,coupling_restored_best,"
               "coupling_relaxation,coupling_R_g,gas_volume_excess,gas_eos_defect,energy_residual_W,internal_energy_J,"
               "liquid_residual_kg_s,liquid_balance_scale_kg_s,liquid_residual_rel,"
               "energy_balance_residual_W,energy_balance_scale_W,energy_viscous_heat_W,"
               "energy_wall_heat_W,energy_axial_enthalpy_out_W,energy_groove_enthalpy_in_W,"
               "energy_mass_imbalance_correction_W,energy_feed_mass_kg_s,energy_pressure_work_W,"
               "cross_film_profile_t_min_K,cross_film_profile_t_max_K,cross_film_profile_residual\n";
    }
    csv.precision(15);
    csv << step << ',' << time << ',' << stats.outer_iters << ',' << stats.flag_flips << ','
        << stats.max_theta_change << ',' << (stats.converged ? 1 : 0) << ','
        << balance.liquid_mass << ',' << balance.liquid_boundary_flux << ','
        << balance.inlet_mass_source << ',' << balance.clamp_mass_liquid << ','
        << balance.liquid_residual << ',' << balance.gas_mass << ','
        << balance.gas_boundary_flux << ',' << balance.clamp_mass_gas << ','
        << balance.gas_residual << ',' << balance.cavitated_fraction;
    if (coupling) {
        csv << ',' << coupling->iters << ',' << coupling->R_p << ',' << coupling->R_T
            << ',' << (coupling->converged ? 1 : 0)
            << ',' << coupling->best_iter << ',' << coupling->best_residual
            << ',' << (coupling->restored_best ? 1 : 0)
            << ',' << coupling->relaxation << ',' << coupling->R_g;
    } else {
        csv << ",,,,,,,,,";
    }
    csv << ',' << balance.gas_volume_excess << ',' << balance.gas_eos_defect << ',';
    if (std::isfinite(balance.energy_residual_w)) csv << balance.energy_residual_w;
    csv << ',';
    if (std::isfinite(balance.internal_energy_j)) csv << balance.internal_energy_j;
    csv << ',';
    if (std::isfinite(balance.liquid_residual_kg_s)) csv << balance.liquid_residual_kg_s;
    csv << ',';
    if (std::isfinite(balance.liquid_balance_scale_kg_s)) csv << balance.liquid_balance_scale_kg_s;
    csv << ',';
    if (std::isfinite(balance.liquid_residual_rel)) csv << balance.liquid_residual_rel;
    if (balance.energy) {
        const auto& energy = *balance.energy;
        csv << ',' << energy.residual_W << ',' << energy.scale_W
            << ',' << energy.viscous_heat_W << ',' << energy.wall_heat_W
            << ',' << energy.axial_enthalpy_out_W << ',' << energy.groove_enthalpy_in_W
            << ',' << energy.mass_imbalance_correction_W << ',' << energy.feed_mass_kg_s << ',' << energy.pressure_work_W;
    } else csv << ",,,,,,,,,";
    csv << ',';
    if (std::isfinite(balance.cross_film_profile_t_min_K)) csv << balance.cross_film_profile_t_min_K;
    csv << ',';
    if (std::isfinite(balance.cross_film_profile_t_max_K)) csv << balance.cross_film_profile_t_max_K;
    csv << ',';
    if (std::isfinite(balance.cross_film_profile_residual)) csv << balance.cross_film_profile_residual;
    csv << '\n';
}

}  // namespace Diagnostics
