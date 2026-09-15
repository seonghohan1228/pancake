#include "energy.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <mpi.h>

#include "cross_film.hpp"
#include "full_film.hpp"
#include "fvm.hpp"
#include "reynolds_coefficients.hpp"

namespace Energy {
namespace {

/// Fraction of the gap occupied by liquid, for weighting the thermal terms.
///
/// A ruptured cell holds liquid theta*h per unit area, not h, in both the
/// striated and the homogeneous picture. Under STRIATED every thermal term that
/// scales with liquid content -- heat generation, stored energy, conduction and
/// convected enthalpy -- is weighted by this. Under FULL_FILM (the default) the
/// cavitated zone is not thermally modelled and this returns 1, which is the
/// historical behaviour.
///
/// `g` is the caller's full-film switch, passed in so cell and face users share
/// one definition of "cavitated".
double liquid_fraction(const Field& theta, const SimulationConfig& cfg,
                       double g, int i, int j) {
    if (cfg.cavitated_film_model != CavitatedFilmModel::STRIATED || g > 0.0) {
        return 1.0;
    }
    return std::clamp(theta(i, j), 0.0, 1.0);
}

/// Liquid fraction on a face, as the mean of the two cells it separates.
double liquid_fraction_face(const Field& theta, const SimulationConfig& cfg,
                            double g_face, int i0, int j0, int i1, int j1) {
    if (cfg.cavitated_film_model != CavitatedFilmModel::STRIATED || g_face > 0.0) {
        return 1.0;
    }
    return 0.5 * (std::clamp(theta(i0, j0), 0.0, 1.0) +
                  std::clamp(theta(i1, j1), 0.0, 1.0));
}

bool uses_pressure_boundary(BCType type) {
    return type == BCType::DIRICHLET || type == BCType::INLET_OUTLET;
}

// Pressure the flow solver actually imposed at an axial boundary. Both
// mass-conserving models clamp the configured value at p_cav (Elrod through
// elrod_boundary_pressure, FBNS through its axial_boundary_pressure), so the
// thermal fluxes must use the same effective pressure or they will transport
// enthalpy against a boundary gradient the flow solve never saw.
double solved_boundary_pressure(const SimulationConfig& cfg, double bc_val) {
    if (cfg.mass_conserving_cavitation()) {
        return cfg.elrod_boundary_pressure(bc_val);
    }
    return bc_val;
}

/// Poiseuille activity of a face for the convected enthalpy.
///
/// Elrod carries a barotropic pressure p = p_cav + beta ln(theta) through the
/// ruptured zone, which is not the physical pressure there, so its Poiseuille
/// flux is switched off wherever the film is not full (`g_face`).
///
/// FBNS pins a ruptured cell at p = p_cav by complementarity, so every pressure
/// gradient it produces is the physical one and no switch belongs in the flux:
/// between two ruptured cells the gradient vanishes on its own, and at a rupture
/// front or an open axial end the flux is exactly the reformation flux the FBNS
/// mass balance conserved. Switching it off would leave the energy equation
/// convecting enthalpy on a different velocity field from the one carrying mass
/// -- the same defect `consistent_boundary_flux` exists to repair on the Elrod
/// path, which FBNS gets for free.
double poiseuille_face_activity(const SimulationConfig& cfg, double g_face) {
    return cfg.cavitation_model == CavitationModel::FBNS ? 1.0 : g_face;
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
            const double g = FullFilm::gate(fields, cfg, i, j);
            const double mu_val = dynamic_viscosity(fields, cfg, i, j);
            // Only the liquid fraction theta carries Couette shear in a cavitated
            // cell; full-film cells (g=1) are unchanged. Shares one helper with the
            // capacity, conduction and convective flux weights so the four cannot
            // drift apart: weighting shear alone conserves the leading-order
            // temperature rise by accident (source and mass flux both scale with
            // theta and the error cancels) while leaving conduction and storage
            // over-weighted by 1/theta.
            const double liquid = liquid_fraction(theta, cfg, g, i, j);
            // Depth-integrated Couette dissipation. With uniform viscosity the shear
            // stress is mu U/h and the integral is mu U^2/h. With a cross-film
            // temperature profile the stress is still uniform across the gap (the
            // thin-film momentum balance with no pressure gradient) but its value is
            // U/F0, and int tau^2/mu dy = U^2/F0 exactly. The hot core raises F0
            // above h/mu(T_mean), so this is strictly less than the uniform-viscosity
            // estimate -- the correct, not an approximated, replacement.
            const double couette_specific = (CrossFilm::active(cfg) && fields.has("cross_film_f0") &&
                                             fields["cross_film_f0"](i, j) > 0.0)
                ? U * U / fields["cross_film_f0"](i, j)
                : mu_val * U * U / h_val;
            const double couette = couette_specific * liquid;

            // Pressure-flow dissipation must follow the same active faces as the
            // Reynolds/velocity solve. Cavitated faces carry no Poiseuille flow, so
            // their pressure jump is not a cell-centered gradient. Average squared
            // active face gradients; inactive and exterior faces contribute zero.
            //
            // The full-film switch is kept here for every cavitation model,
            // including FBNS. Attributing a depth-integrated h^3|grad p|^2/(12 mu)
            // to a ruptured cell is a modelling claim about a striated film, not a
            // conservation requirement the way the convected enthalpy above is; the
            // only cells it can touch are the single rupture/reformation front cells
            // (deeper in the cavity FBNS holds p = p_cav, so grad p vanishes), where
            // Couette shear outweighs it by orders of magnitude in a journal bearing.
            double poiseuille = 0.0;
            if (g > 0.0) {
                double grad_sq = 0.0;
                const double g_e = FullFilm::face_gate(fields, cfg, i, j, i + 1, j);
                const double g_w = FullFilm::face_gate(fields, cfg, i - 1, j, i, j);
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
                        const double g_n = FullFilm::face_gate(fields, cfg, i, j, i, j + 1);
                        if (g_n > 0.0) {
                            const double dp_dz_n = (p(i, j + 1) - p(i, j)) / d_z;
                            grad_sq += dp_dz_n * dp_dz_n;
                        }
                    } else if (j == n_z - 1) {
                        const double g_s = FullFilm::face_gate(fields, cfg, i, j - 1, i, j);
                        if (g_s > 0.0) {
                            const double dp_dz_s = (p(i, j) - p(i, j - 1)) / d_z;
                            grad_sq += dp_dz_s * dp_dz_s;
                        }
                    } else {
                        const double g_n = FullFilm::face_gate(fields, cfg, i, j, i, j + 1);
                        const double g_s = FullFilm::face_gate(fields, cfg, i, j - 1, i, j);
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
                // Pressure-flow dissipation equals the pressure work -q.grad p, and
                // q = -F2 grad p for the pressure-driven part, so the exact
                // depth-integrated value is F2 |grad p|^2. F2 -> h^3/(12 mu) under
                // uniform viscosity, which is the historical expression. The
                // Couette-Poiseuille cross term integrates to zero identically, so
                // the two contributions simply add.
                poiseuille = (CrossFilm::active(cfg) && fields.has("cross_film_f2"))
                    ? fields["cross_film_f2"](i, j) * grad_sq
                    : h_val * h_val * h_val * grad_sq / (12.0 * mu_val);
            }
            heat(i, j) = couette + poiseuille;
        }
    }
}

// Conjugate wall closure. CrossFilm::update has already solved the shaft and bush
// surface temperatures and reduced them to a per-cell flux and its derivative with
// respect to the depth-averaged temperature, so the energy equation only has to
// apply the Newton linearisation
//     Q_w(T) ~ Q_w(T_lag) + dQ/dT (T - T_lag)
// around the lagged film temperature. This replaces the constant-T_wall cooling
// entirely, which is why the two h_wall coefficients are ignored when it is active.
void add_conjugate_wall_heat_loss(LinearSystem& sys, const Fields& fields,
                                  const Field& temperature, const Mesh& mesh) {
    const Field& flux = fields["cross_film_wall_flux"];
    const Field& coeff = fields["cross_film_wall_coeff"];
    const double area = mesh.R * mesh.get_d_theta() * mesh.get_d_z();

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double c = coeff(i, j) * area;
            sys.a_p(i, j) += c;
            sys.source(i, j) += c * temperature(i, j) - flux(i, j) * area;
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
    if (cfg.is_reformation_boundary(theta(i, 0), FullFilm::bound(fields, cfg, i, 0))) {
        u_z = cfg.reformation_boundary_velocity(cfg.bc_z_south_val, h_face, mu_face, 0.5 * d_z, true);
    } else {
        const double g_face =
            poiseuille_face_activity(cfg, FullFilm::gate(fields, cfg, i, 0));
        const double p_bc = solved_boundary_pressure(cfg, cfg.bc_z_south_val);
        const double dp_dz = (p(i, 0) - p_bc) / (0.5 * d_z);
        // Same generalized mobility as the interior faces and the Reynolds solve.
        // Hand-rolling h^3/(12 mu) here left the energy equation transporting heat
        // through its interior faces with Dowson's F2 and out through its axial
        // ends with the classical law -- on the Ferron cross-film case, which has
        // Dirichlet axial ends, those are the oil outflow faces.
        u_z = ReynoldsCoeff::poiseuille_face_velocity(fields, cfg, g_face, dp_dz,
                                                      h_face, mu_face, 1.0, i, 0, i, 0);
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
    if (cfg.is_reformation_boundary(theta(i, j), FullFilm::bound(fields, cfg, i, j))) {
        u_z = cfg.reformation_boundary_velocity(cfg.bc_z_north_val, h_face, mu_face, 0.5 * d_z, false);
    } else {
        const double g_face =
            poiseuille_face_activity(cfg, FullFilm::gate(fields, cfg, i, j));
        const double p_bc = solved_boundary_pressure(cfg, cfg.bc_z_north_val);
        const double dp_dz = (p_bc - p(i, j)) / (0.5 * d_z);
        u_z = ReynoldsCoeff::poiseuille_face_velocity(fields, cfg, g_face, dp_dz,
                                                      h_face, mu_face, 1.0, i, j, i, j);
    }
    return rho_cp_face * h_face * u_z * face_length;
}

void add_axial_boundary_convection(LinearSystem& sys, const Field& p, const Field& h,
                                   const Field& theta, const Mesh& mesh,
                                   const Fields& fields, const SimulationConfig& cfg) {
    // CONSTANT: entering fluid carries that end reservoir's temperature.
    // ZERO_GRADIENT: entering fluid extrapolates the adjacent film temperature.
    //   This assumes thermal recirculation, rather than a specified reservoir. The
    //   non-conservative correction below already accounts for div(F) at the cell,
    //   so the open boundary face contributes no net source — no spurious heating
    //   or cooling. flux_south is the +z-flux at the south face (>0 means inflow at
    //   the south end); flux_north is the +z-flux at the north face (>0 means outflow).
    const bool south_reservoir = cfg.bc_z_south_thermal == ThermalInflowMode::CONSTANT;
    const bool north_reservoir = cfg.bc_z_north_thermal == ThermalInflowMode::CONSTANT;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double flux_south = south_boundary_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_south >= 0.0 && south_reservoir) {
            sys.source(i, 0) += flux_south * cfg.axial_inflow_temperature(true);
        } else {
            sys.a_p(i, 0) += -flux_south;  // outflow, or open zero-gradient inflow
        }

        const int j = mesh.n_z_local - 1;
        const double flux_north = north_boundary_flux(p, h, theta, fields, mesh, cfg, i);
        if (flux_north < 0.0 && north_reservoir) {
            sys.source(i, j) += -flux_north * cfg.axial_inflow_temperature(false);
        } else {
            sys.a_p(i, j) += flux_north;  // outflow, or open zero-gradient inflow
        }
    }
}

/// Hot-oil carry-over at one fed inlet.
///
/// A real groove does not deliver pure supply oil: the film sweeping into it
/// recirculates and mixes with the fresh charge, so the oil leaving the groove is
///
///     T_groove = lambda T_supply + (1 - lambda) T_upstream,
///
/// with lambda the fresh-oil share of the mixture. `t_mix` is filled per axial row
/// (the groove spans the full width and the upstream film is not uniform in z);
/// `lambda_out` receives the row-independent mixing fraction actually used.
///
/// The upstream film state is read one cell circumferentially ahead of the
/// feature's leading edge, in the sense the shaft entrains. That cell can live on
/// any rank, so its temperature and its Couette carry-over flux are gathered with
/// a sum reduction over an array indexed by the (undecomposed) axial index.
void groove_mixing_temperatures(const InletConfig& inlet, const Fields& fields,
                                const Mesh& mesh, const SimulationConfig& cfg,
                                std::vector<double>& t_mix, double& lambda_out) {
    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    t_mix.assign(n_z, inlet.t_supply);
    lambda_out = 1.0;
    if (cfg.groove_mixing_model == GrooveMixingModel::PINNED) return;

    const Field& temperature = fields["temperature"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const double d_z = mesh.get_d_z();
    const int downstream = cfg.omega >= 0.0 ? 1 : -1;

    // Leading edge in global theta indices: the first inlet cell met when walking
    // downstream, so the upstream sample sits one cell before it.
    std::vector<double> gathered(static_cast<size_t>(2) * n_z, 0.0);
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            if (mesh.inside_inlet(inlet, i, j)) continue;
            // (i) is outside; it is the upstream sample only if the next cell
            // downstream is inside the feature.
            if (!mesh.inside_inlet(inlet, i + downstream, j)) continue;
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double liquid = std::clamp(theta(i, j), 0.0, 1.0);
            // Volumetric Couette carry-over through this row, m^3/s.
            const double q_up = 0.5 * std::abs(cfg.omega) * mesh.R * h_val * liquid * d_z;
            gathered[j] += temperature(i, j) * q_up;
            gathered[n_z + j] += q_up;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, gathered.data(), static_cast<int>(gathered.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double q_up_total = 0.0;
    for (int j = 0; j < n_z; ++j) q_up_total += gathered[n_z + j];

    double lambda = 1.0;
    if (cfg.groove_mixing_model == GrooveMixingModel::FIXED_LAMBDA) {
        lambda = std::clamp(inlet.mixing_lambda, 0.0, 1.0);
    } else if (q_up_total > 0.0) {
        // Flow-weighted closure: the groove delivers the supply stream and the
        // recirculating film in proportion to the two flow rates.
        lambda = inlet.supply_flow_rate / (inlet.supply_flow_rate + q_up_total);
    }
    lambda_out = lambda;

    for (int j = 0; j < n_z; ++j) {
        const double q = gathered[n_z + j];
        const double t_up = q > 0.0 ? gathered[j] / q : inlet.t_supply;
        t_mix[j] = lambda * inlet.t_supply + (1.0 - lambda) * t_up;
    }
}

// Actual fed inlets (feeds_fresh_oil) pump fresh oil in at a fixed supply
// temperature, so their rows impose that value directly. Under the default
// groove_mixing_model = PINNED the pinned value is inlet.t_supply exactly; the
// mixing models pin the groove to the hot-oil carry-over mixture instead.
// Open large-clearance regions (feeds_fresh_oil == false) impose no thermal
// condition and are left thermally transparent to the energy solve.
void apply_inlet_thermal_conditions(LinearSystem& sys, const Fields& fields,
                                    const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.inlets.empty()) return;

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    std::vector<double> t_mix;
    double lambda = 1.0;

    for (const auto& inlet : cfg.inlets) {
        if (!inlet.feeds_fresh_oil) continue;  // open region: thermally transparent
        groove_mixing_temperatures(inlet, fields, mesh, cfg, t_mix, lambda);

        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                if (!mesh.inside_inlet(inlet, i, j)) continue;
                sys.set_dirichlet(i, j, t_mix[j]);
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
            // Liquid depth theta*h, not the full gap: a ruptured cell neither
            // stores nor conducts as though it were solid oil.
            const double liquid_h =
                h_val * liquid_fraction(theta, cfg, FullFilm::gate(fields, cfg, i, j), i, j);
            conductivity(i, j) = thermal_conductivity(fields, cfg, i, j) * liquid_h;
            if (i >= 0 && i < n_theta) {
                capacity(i, j) = volumetric_heat_capacity(fields, cfg, i, j) * liquid_h;
            }
        }
    }

    Field flux_theta("thermal_flux_theta", mesh);
    Field flux_z("thermal_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i + 1, j)), cfg.min_film_thickness);
            const double g_face = FullFilm::face_gate(fields, cfg, i, j, i + 1, j);
            // Separate switches on purpose: `g_flow` decides whether the face
            // carries Poiseuille flow at all (FBNS always does, see
            // poiseuille_face_activity), while `g_face` stays the "is this face
            // full film" test the STRIATED liquid weighting keys on.
            const double g_flow = poiseuille_face_activity(cfg, g_face);
            const double dp_dtheta = (p(i + 1, j) - p(i, j)) / d_theta;
            const double mu_face = 0.5 * (
                dynamic_viscosity(fields, cfg, i, j) +
                dynamic_viscosity(fields, cfg, i + 1, j));
            const double rho_cp_face = 0.5 * (
                volumetric_heat_capacity(fields, cfg, i, j) +
                volumetric_heat_capacity(fields, cfg, i + 1, j));
            const double u_theta = 0.5 * cfg.omega * R
                - g_flow * h_face * h_face * dp_dtheta / (12.0 * mu_face * R);
            // Enthalpy rides on the liquid that is actually there. Carrying the
            // full film through the cavitated arc overstates the mass flux by
            // 1/theta -- a factor of ~2.5 at the measured theta_min ~ 0.39 --
            // which damps the temperature rise reaching the next inlet.
            const double liquid_face =
                liquid_fraction_face(theta, cfg, g_face, i, j, i + 1, j);
            if (CrossFilm::active(cfg) && fields.has("cross_film_f2")) {
                // Volumetric flux per unit width from the same generalized law the
                // Reynolds solve conserved, q = U F1/F0 - F2 grad p, so enthalpy and
                // mass ride on one velocity field rather than two.
                const double couette_face = 0.5 * (fields["cross_film_couette_h"](i, j) +
                                                   fields["cross_film_couette_h"](i + 1, j));
                const double f2_face = 0.5 * (fields["cross_film_f2"](i, j) +
                                              fields["cross_film_f2"](i + 1, j));
                const double q_face =
                    0.5 * cfg.omega * R * couette_face - g_flow * f2_face * dp_dtheta / R;
                flux_theta(i, j) = rho_cp_face * liquid_face * q_face * d_z;
            } else {
                flux_theta(i, j) = rho_cp_face * liquid_face * h_face * u_theta * d_z;
            }
        }
    }

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) flux_z(i, j) = 0.0;
        for (int j = 0; j + 1 < n_z; ++j) {
            const double h_face = std::max(0.5 * (h(i, j) + h(i, j + 1)), cfg.min_film_thickness);
            const double g_face = FullFilm::face_gate(fields, cfg, i, j, i, j + 1);
            const double g_flow = poiseuille_face_activity(cfg, g_face);
            const double dp_dz = (p(i, j + 1) - p(i, j)) / d_z;
            const double mu_face = 0.5 * (
                dynamic_viscosity(fields, cfg, i, j) +
                dynamic_viscosity(fields, cfg, i, j + 1));
            const double rho_cp_face = 0.5 * (
                volumetric_heat_capacity(fields, cfg, i, j) +
                volumetric_heat_capacity(fields, cfg, i, j + 1));
            const double u_z = -g_flow * h_face * h_face * dp_dz / (12.0 * mu_face);
            const double liquid_face =
                liquid_fraction_face(theta, cfg, g_face, i, j, i, j + 1);
            if (CrossFilm::active(cfg) && fields.has("cross_film_f2")) {
                const double f2_face = 0.5 * (fields["cross_film_f2"](i, j) +
                                              fields["cross_film_f2"](i, j + 1));
                const double q_face = -g_flow * f2_face * dp_dz;
                flux_z(i, j) = rho_cp_face * liquid_face * q_face * R * d_theta;
            } else {
                flux_z(i, j) = rho_cp_face * liquid_face * h_face * u_z * R * d_theta;
            }
        }
    }

    sys.reset();
    FVM::laplacian(sys, conductivity, mesh);
    // TYPE_DIFFERENCING is theta-only (validated); fall back to upwind defensively.
    const ConvectionScheme thermal_scheme =
        cfg.thermal_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING
            ? ConvectionScheme::UPWIND : cfg.thermal_convection_scheme;
    // Deferred corrections read phi.old, which on the live field is the previous
    // timestep: stale across outer iterations, and in STEADY_STATE (one timestep,
    // many outer iterations) never refreshed at all, which silently degrades every
    // high-order scheme to first-order upwind. Carry the current iterate through a
    // lagged copy, as solve_elrod does for theta.
    Field temperature_lag("temperature_lag", mesh);
    temperature_lag.data = temperature.data;
    temperature_lag.store_old_time();
    FVM::divergence(sys, flux_theta, flux_z, temperature_lag, thermal_scheme, mesh);
    if (cfg.solution_mode != SolutionMode::STEADY_STATE) {
        FVM::ddt_weighted(sys, temperature, capacity, cfg.dt, mesh);
    }
    FVM::add_source(sys, heat, mesh);
    // The conjugate closure replaces the Newton cooling outright: with the shaft
    // and bush surface temperatures solved, a prescribed-wall-temperature term on
    // top of them would double-count the same heat path.
    if (CrossFilm::active(cfg) && fields.has("cross_film_wall_flux")) {
        add_conjugate_wall_heat_loss(sys, fields, temperature, mesh);
    } else {
        add_wall_heat_loss(sys, mesh, cfg);
    }
    add_axial_boundary_convection(sys, p, h, theta, mesh, fields, cfg);
    add_nonconservative_convection_correction(sys, flux_theta, flux_z, p, h, theta, fields, mesh, cfg);
    apply_inlet_thermal_conditions(sys, fields, mesh, cfg);

    sys.solve(temperature, cfg.linear_rtol, 1000, "temperature");
}

}  // namespace Energy
