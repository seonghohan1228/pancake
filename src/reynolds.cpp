#include "reynolds.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <mpi.h>

#include "cross_film.hpp"
#include "equation_of_state.hpp"
#include "full_film.hpp"
#include "fvm.hpp"
#include "communicator.hpp"
#include "fluid_properties.hpp"
#include "reynolds_coefficients.hpp"

namespace Reynolds {

namespace {
double bounded_film_content(double theta, const SimulationConfig& cfg) {
    return std::max(theta, cfg.theta_min);
}

double field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

double positive_field_value(const Fields& fields, const char* name, double fallback, int i, int j) {
    return std::max(field_value(fields, name, fallback, i, j), 1.0e-30);
}

double theta_base_density(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
}

double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "mu", cfg.mu, i, j);
}

// The Poiseuille mobility, the generalized Couette thickness and the void/EOS
// parameters are shared with FBNS::solve and with Diagnostics, which must
// re-evaluate exactly the operator assembled here. See reynolds_coefficients.hpp.
using ReynoldsCoeff::VoidParams;
using ReynoldsCoeff::couette_thickness;
using ReynoldsCoeff::gas_inventory_plateau;
using ReynoldsCoeff::poiseuille_coefficient;
using ReynoldsCoeff::void_params;

double reynolds_diffusion(const Fields& fields, const SimulationConfig& cfg,
                          double prefactor, int i, int j) {
    return ReynoldsCoeff::poiseuille(fields, cfg, prefactor, i, j);
}

// Pressure the solver actually imposed at an axial boundary. Elrod derives its
// film-content boundary from the configured pressure and bulk modulus, clamped
// at p_cav. Velocity and force post-processing must match the solve.
double solved_boundary_pressure(const SimulationConfig& cfg, double bc_val) {
    if (cfg.cavitation_model == CavitationModel::JFO) {
        return cfg.elrod_boundary_pressure(bc_val);
    }
    return bc_val;
}

double van_leer_limiter(double r) {
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

// Explicit face interpolation of phi for the Gumbel Couette flux, matched to the
// configured theta scheme so the Gumbel and Elrod paths discretize the wedge term
// identically. Values phi_uu, phi_u, phi_d follow the upwind direction; central_face
// enables central interpolation under TYPE_DIFFERENCING.
double couette_face_value(double phi_uu, double phi_u, double phi_d,
                          ConvectionScheme scheme, bool central_face) {
    switch (scheme) {
        case ConvectionScheme::UPWIND:
            return phi_u;
        case ConvectionScheme::LINEAR:
            return 0.5 * (phi_u + phi_d);  // central everywhere (2nd-order, unbounded)
        case ConvectionScheme::VANLEER: {
            const double denom = phi_u - phi_uu;
            const double num = phi_d - phi_u;
            const double r = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
            return phi_u + 0.5 * van_leer_limiter(r) * (phi_u - phi_uu);
        }
        case ConvectionScheme::TYPE_DIFFERENCING:
            return central_face ? 0.5 * (phi_u + phi_d) : phi_u;
    }
    return phi_u;
}
}

/// Replace inlet rows with their prescribed pressure or Elrod film content.
static void apply_inlet_conditions(LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.inlets.empty()) return;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const auto* inlet = mesh.inlet_at(cfg,i,j);
            if (!inlet) continue;
            const double value = cfg.cavitation_model == CavitationModel::JFO
                ? EOS::theta_from_pressure(inlet->p_supply, cfg.effective_p_cav(), cfg.bulk_modulus)
                : inlet->p_supply;
            sys.set_dirichlet(i,j,value);
        }
    }
}

void solve(Fields& fields, LinearSystem& sys,
           const Mesh& mesh, const SimulationConfig& cfg)
{
    Field& pressure  = fields["pressure"];
    const Field& h   = fields["h"];
    Field& rho       = fields["rho"];

    const int n_theta    = mesh.n_theta_local;
    const int n_z        = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;
    const int ng         = h.n_ghost;

    // Compute gamma = ρ·(Poiseuille coefficient) over physical cells plus theta
    // ghost layers: ρh³/(12μ) normally, ρF2 under the generalized (cross-film)
    // Reynolds equation. h, rho and the cross-film theta-ghost cells must already
    // be synced before calling this function.
    Field gamma_f("gamma",  mesh);
    for (int i = -ng; i < n_theta + ng; ++i) {
        for (int j = 0; j < n_z; ++j) {
            gamma_f(i, j) = reynolds_diffusion(fields, cfg, rho(i, j), i, j);
        }
    }

    // Assemble diffusion term: ∇·(γ ∇p)
    sys.reset();
    FVM::laplacian(sys, gamma_f, mesh);

    // Dirichlet BC: pressure = bc_z_val at z = 0 (south) and z = L (north).
    // The laplacian skips these boundary faces; add them here using half-cell distance.
    for (int i = 0; i < n_theta; ++i) {
        if (cfg.bc_z_south_type == BCType::DIRICHLET || cfg.bc_z_south_type == BCType::INLET_OUTLET) {
            const double coeff_s = gamma_f(i, 0)       * R * d_theta / (0.5 * d_z);
            sys.a_p(i, 0)        += coeff_s;
            sys.source(i, 0)     += coeff_s * cfg.bc_z_south_val;
        }
        if (cfg.bc_z_north_type == BCType::DIRICHLET || cfg.bc_z_north_type == BCType::INLET_OUTLET) {
            const double coeff_n = gamma_f(i, n_z - 1) * R * d_theta / (0.5 * d_z);
            sys.a_p(i, n_z - 1)  += coeff_n;
            sys.source(i, n_z - 1) += coeff_n * cfg.bc_z_north_val;
        }
    }

    // Couette source: stencil assembles -div(γ∇p)·vol = source_field·vol,
    // so source_field = -(ω/2 · ∂(ρh)/∂θ), evaluated as a face-flux divergence
    // matched to the Elrod path: ρh factorizes into ρ₀h (central face value,
    // like the Elrod Couette flux) times the transported content θ = ρ/ρ₀,
    // which is interpolated with the configured theta scheme.
    Field couette_src("couette_src", mesh);
    Field rho0_h("rho0_h", mesh);
    Field theta_lag_f("theta_lag", mesh);
    for (int i = -ng; i < n_theta + ng; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double rho0 = theta_base_density(fields, cfg, i, j);
            rho0_h(i, j) = rho0 * couette_thickness(fields, cfg, i, j);
            theta_lag_f(i, j) = rho(i, j) / rho0;
        }
    }
    const bool downwind_is_east = cfg.omega >= 0.0;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const bool central_e = pressure(i, j) > cfg.effective_p_cav() && pressure(i + 1, j) > cfg.effective_p_cav();
            const bool central_w = pressure(i - 1, j) > cfg.effective_p_cav() && pressure(i, j) > cfg.effective_p_cav();
            double theta_e, theta_w;
            if (downwind_is_east) {
                theta_e = couette_face_value(theta_lag_f(i - 1, j), theta_lag_f(i, j),
                                             theta_lag_f(i + 1, j),
                                             cfg.theta_convection_scheme, central_e);
                theta_w = couette_face_value(theta_lag_f(i - 2, j), theta_lag_f(i - 1, j),
                                             theta_lag_f(i, j),
                                             cfg.theta_convection_scheme, central_w);
            } else {
                theta_e = couette_face_value(theta_lag_f(i + 2, j), theta_lag_f(i + 1, j),
                                             theta_lag_f(i, j),
                                             cfg.theta_convection_scheme, central_e);
                theta_w = couette_face_value(theta_lag_f(i + 1, j), theta_lag_f(i, j),
                                             theta_lag_f(i - 1, j),
                                             cfg.theta_convection_scheme, central_w);
            }
            const double phi_e = 0.5 * (rho0_h(i, j) + rho0_h(i + 1, j)) * theta_e;
            const double phi_w = 0.5 * (rho0_h(i - 1, j) + rho0_h(i, j)) * theta_w;
            couette_src(i, j) = -0.5 * cfg.omega * (phi_e - phi_w) / d_theta;
        }
    }
    FVM::add_source(sys, couette_src, mesh);

    if (cfg.solution_mode != SolutionMode::STEADY_STATE && fields.has("dh_dt")) {
        const Field& dh_dt = fields["dh_dt"];
        Field squeeze_src("squeeze_src", mesh);
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                squeeze_src(i, j) = -rho(i, j) * dh_dt(i, j);
            }
        }
        FVM::add_source(sys, squeeze_src, mesh);
    }

    apply_inlet_conditions(sys, mesh, cfg);

    sys.solve(pressure, cfg.linear_rtol, 1000, "pressure");

    // Gumbel cavitation condition: enforce p >= p_cav
    if (cfg.cavitation_model == CavitationModel::GUMBEL) {
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                pressure(i, j) = std::max(pressure(i, j), cfg.effective_p_cav());
    }

    // Update density via barotropic EOS (rho = rho_0 * theta(p))
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j) {
            double theta_p = EOS::theta_from_pressure(pressure(i, j), cfg.effective_p_cav(), cfg.bulk_modulus);
            // In Full Sommerfeld, we should ensure density doesn't drop due to negative pressure 
            // if we want to follow classical liquid theory, but we'll follow EOS here.
            rho(i, j) = theta_base_density(fields, cfg, i, j) * theta_p;
            if (fields.has("theta")) fields["theta"](i, j) = theta_p;
        }
}

ElrodStats solve_elrod(Fields& fields, LinearSystem& sys,
                       const Mesh& mesh, const SimulationConfig& cfg)
{
    ElrodStats stats;
    Field& theta    = fields["theta"];
    Field& pressure = fields["pressure"];
    const Field& h  = fields["h"];
    Field& rho      = fields["rho"];

    const int n_theta    = mesh.n_theta_local;
    const int n_z        = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;
    const int ng         = h.n_ghost;

    // Γ_base = ρ₀β·(Poiseuille coefficient) — effective diffusion without switch
    // function: ρ₀βh³/(12μ) normally, ρ₀βF2 under the generalized Reynolds
    // equation. h and the cross-film theta-ghost cells must already be synced.
    Field Gamma_base("Gamma_base", mesh);
    for (int i = -ng; i < n_theta + ng; ++i) {
        for (int j = 0; j < n_z; ++j) {
            Gamma_base(i, j) = ReynoldsCoeff::elrod_gamma_base(fields, cfg, i, j);
        }
    }

    // Couette convection face flux: F_e(i,j) = ρ₀·(ω/2)·R·h_e·Δz at the east face.
    // Extends to i = -1 so FVM::divergence can access the west-face ghost at i=0.
    Field couette_flux("couette_flux", mesh);
    Field zero_flux_z("zero_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            couette_flux(i, j) = ReynoldsCoeff::couette_face_flux(fields, mesh, cfg, i, j);
            zero_flux_z(i, j)  = 0.0;
        }
    }

    // Γ_g = Γ_base · g(θ): diffusion disabled in cavitated cells.
    Field Gamma_g("Gamma_g", mesh);

    // Lagged iterate carrier for deferred-correction schemes (FVM::divergence reads
    // phi.old). theta.old itself must stay at the previous timestep for ddt_weighted.
    const bool deferred_scheme =
        cfg.theta_convection_scheme != ConvectionScheme::UPWIND;
    Field theta_lag("theta_lag", mesh);
    Field central_mask_theta("central_mask_theta", mesh);
    Field central_mask_z("central_mask_z", mesh);

    // Buffer for outer-iteration convergence check (flat array, physical cells only).
    const int n_local = n_theta * n_z;
    std::vector<double> theta_prev(n_local);
    Communicator theta_comm(mesh);

    Field* clamp_mass = fields.has("cavitation_clamp_mass")
        ? &fields["cavitation_clamp_mass"] : nullptr;
    const double cell_area = mesh.cell_volume();

    // Outer (flag-update) iteration: iterate until the switch pattern converges.
    for (int outer = 0; outer < cfg.max_outer_iters; ++outer) {

        // Save θ for convergence check.
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                theta_prev[i * n_z + j] = theta(i, j);

        // VOID_COUPLED: advance the lagged plateau once per outer iteration under
        // relaxation. Updating it here (rather than inside every void_params call)
        // breaks the p_void <-> θ limit cycle; w = 1 reproduces the unrelaxed,
        // non-convergent behaviour.
        if (fields.has("p_void") && outer > 0) {
            Field& p_void = fields["p_void"];
            const double w = std::clamp(cfg.void_coupling_relaxation, 0.01, 1.0);
            const double floor_p = cfg.effective_p_cav();
            for (int i = -ng; i < n_theta + ng; ++i)
                for (int j = 0; j < n_z; ++j) {
                    const double target =
                        std::max(floor_p, gas_inventory_plateau(fields, cfg, theta(i, j), i, j));
                    p_void(i, j) = (1.0 - w) * p_void(i, j) + w * target;
                }
        }

        // Apply switch function g(θ) to get the effective diffusion coefficient.
        // For the first outer iteration (outer=0), we force g=1 everywhere.
        // This is a "flooded initialization" that prevents numerical starvation
        // by establishing a pressure zone from the oil supply boundaries.
        for (int i = -ng; i < n_theta + ng; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double g = (outer == 0) ? 1.0
                    : EOS::switch_function(theta(i, j), void_params(fields, cfg, i, j).theta_full);
                Gamma_g(i, j) = Gamma_base(i, j) * g;
            }
        }

        // Assemble: -div(Γ_g·∇θ)·V + div(F·θ)·V = 0
        // i.e., div(Γ_g·∇θ) = div(F·θ)  (Elrod-Adams with cavitation switch)
        sys.reset();
        FVM::laplacian(sys, Gamma_g, mesh);
        if (deferred_scheme) {
            // Deferred corrections lag one outer iteration behind the implicit
            // upwind part; carry the current iterate through theta_lag.old.
            theta_lag.data = theta.data;
            theta_lag.store_old_time();
            if (cfg.theta_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING) {
                for (int i = -1; i < n_theta; ++i) {
                    for (int j = 0; j < n_z; ++j) {
                        central_mask_theta(i, j) =
                            (theta(i, j) >= void_params(fields, cfg, i, j).theta_full &&
                             theta(i + 1, j) >= void_params(fields, cfg, i + 1, j).theta_full) ? 1.0 : 0.0;
                        central_mask_z(i, j) = 0.0;
                    }
                }
            }
            FVM::divergence(sys, couette_flux, zero_flux_z, theta_lag,
                            cfg.theta_convection_scheme, mesh,
                            &central_mask_theta, &central_mask_z);
        } else {
            FVM::divergence(sys, couette_flux, zero_flux_z, theta,
                            ConvectionScheme::UPWIND, mesh);
        }

        // Conservative storage of rho_l h theta, including density and gap changes.
        if (cfg.solution_mode != SolutionMode::STEADY_STATE) {
            for (int i = 0; i < n_theta; ++i) {
                for (int j = 0; j < n_z; ++j) {
                    const double area_dt = mesh.cell_volume() / cfg.dt;
                    sys.a_p(i, j) += theta_base_density(fields, cfg, i, j) * h(i, j) * area_dt;
                    sys.source(i, j) += ReynoldsCoeff::old_liquid_mass_area(fields, cfg, cfg.dt, i, j) * area_dt;
                }
            }
        }

        // Dirichlet z-BC: θ = θ_boundary at z=0 and z=L.
        // Uses Gamma_base (not Gamma_g) to represent the flooded-bearing condition.
        // One link per axial boundary face. Full-film cells use the barotropic
        // Gamma_base Dirichlet link; under consistent_boundary_flux a cavitated
        // cell instead re-floods at the physical Poiseuille reformation rate
        // cs*(p_bc - p_cav)/beta (an explicit inflow, no implicit theta link),
        // which the Gamma_base link overstates by ~beta/(p_bc-p_cav).
        auto add_boundary_link = [&](double bc_val, int i, int j) {
            if (cfg.consistent_boundary_flux && !FullFilm::is_full_film(fields, cfg, i, j)) {
                // Cavitated: explicit reformation inflow = -outflow (shared formula).
                sys.source(i, j) -= cfg.elrod_boundary_outflow(
                    Gamma_base(i, j), theta(i, j), FullFilm::bound(fields, cfg, i, j),
                    bc_val, R, d_theta, d_z);
            } else {
                // Full-film: implicit barotropic link (effective outflow matches
                // cfg.elrod_boundary_outflow's full-film branch, used by diagnostics).
                const double cs = Gamma_base(i, j) * R * d_theta / (0.5 * d_z);
                sys.a_p(i, j) += cs;
                sys.source(i, j) += cs * bounded_film_content(cfg.elrod_boundary_theta(bc_val), cfg);
            }
        };
        for (int i = 0; i < n_theta; ++i) {
            // South boundary (z=0). INLET_OUTLET adds the link only on inflow.
            if (cfg.bc_z_south_type == BCType::DIRICHLET) {
                add_boundary_link(cfg.bc_z_south_val, i, 0);
            } else if (cfg.bc_z_south_type == BCType::INLET_OUTLET) {
                if (fields["pressure"](i, 0) < cfg.bc_z_south_val) add_boundary_link(cfg.bc_z_south_val, i, 0);
            }

            // North boundary (z=L).
            if (cfg.bc_z_north_type == BCType::DIRICHLET) {
                add_boundary_link(cfg.bc_z_north_val, i, n_z - 1);
            } else if (cfg.bc_z_north_type == BCType::INLET_OUTLET) {
                if (fields["pressure"](i, n_z - 1) < cfg.bc_z_north_val) add_boundary_link(cfg.bc_z_north_val, i, n_z - 1);
            }
        }

        apply_inlet_conditions(sys, mesh, cfg);

        sys.solve(theta, cfg.linear_rtol, 1000, "elrod film content");

        // Clamp θ ≥ θ_min and snap near-full-film cells to exactly 1.0, accounting
        // the implied mass so the diagnostics balance still closes (clamps are not
        // mass-neutral). Only the final outer iteration's deltas persist in theta,
        // so the accumulator restarts every iteration.
        //
        // Snap rationale: solver precision (rtol=1e-8, condition ~O(10^3)) produces
        // errors ~O(1e-5) near the free boundary, causing spurious g=0 assignments
        // in the next outer iteration. Physically cavitated cells have θ << 1 and
        // are not affected.
        const double snap_tol = 5e-7;
        if (clamp_mass != nullptr) clamp_mass->fill(0.0);
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double theta_solved = theta(i, j);
                double theta_new = std::max(theta_solved, cfg.theta_min);
                const double tf = void_params(fields, cfg, i, j).theta_full;
                if (theta_new > tf - snap_tol && theta_new < tf) theta_new = tf;
                if (theta_new != theta_solved) {
                    if (clamp_mass != nullptr) {
                        (*clamp_mass)(i, j) += theta_base_density(fields, cfg, i, j) *
                            h(i, j) * cell_area * (theta_new - theta_solved);
                    }
                    theta(i, j) = theta_new;
                }
            }
        }

        // Convergence: max |θ_new − θ_prev| and count of flag (g) changes.
        // Use snap_tol as threshold so near-interface cells do not trigger restarts.
        theta_comm.update_ghosts(theta);

        double max_change = 0.0;
        int    n_switched = 0;
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double tp = theta_prev[i * n_z + j];
                const double tn = theta(i, j);
                max_change = std::max(max_change, std::abs(tn - tp));
                // Same ceiling the flux switches use, widened by snap_tol so a
                // near-interface cell does not count as a flip on round-off.
                const double tf = FullFilm::bound(fields, cfg, i, j);
                if (FullFilm::holds(tp, tf, snap_tol) != FullFilm::holds(tn, tf, snap_tol)) {
                    ++n_switched;
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &max_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &n_switched, 1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD);

        if (cfg.log_outer_iters) {
            int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) {
                double mn = 1e30, mx = -1e30;
                for (int i = 0; i < n_theta; ++i)
                    for (int j = 0; j < n_z; ++j) { mn = std::min(mn, theta(i,j)); mx = std::max(mx, theta(i,j)); }
                std::fprintf(stderr, "outer=%d mc=%.3e ns=%d th=[%.6f,%.6f]\n", outer, max_change, n_switched, mn, mx);
            }
        }
        stats.outer_iters = outer + 1;
        stats.flag_flips += n_switched;
        stats.max_theta_change = max_change;
        if (max_change < cfg.outer_tol && n_switched == 0) {
            stats.converged = true;
            break;
        }
    }

    // Recover pressure and density from EOS.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const VoidParams vp = void_params(fields, cfg, i, j);
            pressure(i, j) = EOS::pressure_from_theta(theta(i, j), vp.theta_full, vp.p_void, vp.beta_bar);
            rho(i, j)      = theta_base_density(fields, cfg, i, j) * theta(i, j);
        }
    }
    return stats;
}

void calculate_velocities(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    Field& u_theta = fields["velocity_theta"];
    Field& u_z     = fields["velocity_z"];

    const int n_theta = mesh.n_theta_local;
    const int n_z     = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;
    const double omega   = cfg.omega;

    // u_theta at FACE_THETA: interface between (i-1, j) and (i, j)
    // u_theta(i, j) is the west face of cell (i, j)
    // We need to calculate for i=0 to n_theta-1. Face i=n_theta is the same as i=0 for periodic.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            // Face values at west face of cell (i, j)
            // Poiseuille term active only if both cells are full film
            double g_f = FullFilm::face_gate(fields, cfg, i - 1, j, i, j);
            double dp_dtheta = (p(i, j) - p(i-1, j)) / d_theta;

            if (CrossFilm::active(cfg)) {
                // Depth-averaged velocity q/h from the generalized Reynolds flux
                // q = U F1/F0 - F2 dp/ds, the same flux the pressure solve conserved.
                const double h_f = std::max(0.5 * (h(i - 1, j) + h(i, j)), cfg.min_film_thickness);
                const double couette_f = 0.5 * (couette_thickness(fields, cfg, i - 1, j) +
                                                couette_thickness(fields, cfg, i, j));
                const double poise_f = 0.5 * (poiseuille_coefficient(fields, cfg, i - 1, j) +
                                              poiseuille_coefficient(fields, cfg, i, j));
                u_theta(i, j) = 0.5 * omega * R * (couette_f / h_f)
                    - g_f * poise_f / (h_f * R) * dp_dtheta;
            } else {
                const double h_f = 0.5 * (h(i-1, j) + h(i, j));
                const double mu_f = 0.5 * (
                    dynamic_viscosity(fields, cfg, i - 1, j) +
                    dynamic_viscosity(fields, cfg, i, j));
                u_theta(i, j) = 0.5 * omega * R - g_f * (h_f * h_f) / (12.0 * mu_f * R) * dp_dtheta;
            }
        }
    }

    // u_z at FACE_Z: interface between (i, j-1) and (i, j)
    // u_z(i, j) is the south face of cell (i, j). 
    // Field u_z has n_z_local + 1 faces in z direction.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j <= n_z; ++j) {
            double p_prev, p_next, h_f, g_f, mu_f, poise_f;
            double dz_eff = d_z;

            if (j == 0) { // South boundary face
                p_next = p(i, 0);
                if (cfg.bc_z_south_type == BCType::DIRICHLET) {
                    p_prev = solved_boundary_pressure(cfg, cfg.bc_z_south_val);
                    dz_eff = 0.5 * d_z;
                } else {
                    p_prev = p_next; // Neumann: dp/dz = 0
                }
                h_f = h(i, 0);
                mu_f = dynamic_viscosity(fields, cfg, i, 0);
                poise_f = poiseuille_coefficient(fields, cfg, i, 0);
                g_f = FullFilm::gate(fields, cfg, i, 0);
            } else if (j == n_z) { // North boundary face
                p_prev = p(i, n_z - 1);
                if (cfg.bc_z_north_type == BCType::DIRICHLET) {
                    p_next = solved_boundary_pressure(cfg, cfg.bc_z_north_val);
                    dz_eff = 0.5 * d_z;
                } else {
                    p_next = p_prev; // Neumann: dp/dz = 0
                }
                h_f = h(i, n_z - 1);
                mu_f = dynamic_viscosity(fields, cfg, i, n_z - 1);
                poise_f = poiseuille_coefficient(fields, cfg, i, n_z - 1);
                g_f = FullFilm::gate(fields, cfg, i, n_z - 1);
            } else { // Interior face
                p_prev = p(i, j - 1);
                p_next = p(i, j);
                h_f = 0.5 * (h(i, j-1) + h(i, j));
                mu_f = 0.5 * (
                    dynamic_viscosity(fields, cfg, i, j - 1) +
                    dynamic_viscosity(fields, cfg, i, j));
                poise_f = 0.5 * (poiseuille_coefficient(fields, cfg, i, j - 1) +
                                 poiseuille_coefficient(fields, cfg, i, j));
                g_f = FullFilm::face_gate(fields, cfg, i, j - 1, i, j);
            }

            double dp_dz = (p_next - p_prev) / dz_eff;
            if (CrossFilm::active(cfg)) {
                u_z(i, j) = - g_f * poise_f / std::max(h_f, cfg.min_film_thickness) * dp_dz;
            } else {
                u_z(i, j) = - g_f * (h_f * h_f) / (12.0 * mu_f) * dp_dz;
            }
        }
    }
}

ForceComponents calculate_macroscopic_properties(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];

    const int n_theta = mesh.n_theta_local;
    const int n_z     = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;
    const double omega   = cfg.omega;

    double local_pressure_x = 0.0;
    double local_pressure_y = 0.0;
    double local_pressure_z = 0.0;
    double local_viscous_x = 0.0;
    double local_viscous_y = 0.0;
    double local_viscous_z = 0.0;
    double local_torque = 0.0;

    for (int i = 0; i < n_theta; ++i) {
        const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
        const double cos_t = std::cos(theta_c);
        const double sin_t = std::sin(theta_c);

        for (int j = 0; j < n_z; ++j) {
            // Gauge pressure is used because uniform ambient p_cav has zero net force.
            // The force convention is the force applied by the fluid to the moving
            // bearing surface. The bearing pressure area vector is
            // ((R+h)e_r - h_theta e_theta - (R+h)h_z e_z) dtheta dz.
            const double p_gauge = p(i, j) - cfg.effective_p_cav();
            const double h_val = h(i, j);
            const double r_surface = R + h_val;
            const double h_theta = (h(i + 1, j) - h(i - 1, j)) / (2.0 * d_theta);
            double h_z = 0.0;
            if (n_z == 1) {
                h_z = 0.0;
            } else if (j == 0) {
                h_z = (h(i, 1) - h(i, 0)) / d_z;
            } else if (j == n_z - 1) {
                h_z = (h(i, n_z - 1) - h(i, n_z - 2)) / d_z;
            } else {
                h_z = (h(i, j + 1) - h(i, j - 1)) / (2.0 * d_z);
            }
            const double area_scale = d_theta * d_z;
            const double area_x = (r_surface * cos_t + h_theta * sin_t) * area_scale;
            const double area_y = (r_surface * sin_t - h_theta * cos_t) * area_scale;
            const double area_z = -r_surface * h_z * area_scale;

            local_pressure_x += p_gauge * area_x;
            local_pressure_y += p_gauge * area_y;
            local_pressure_z += p_gauge * area_z;

            const double dA_bearing = r_surface * d_theta * d_z;

            // Shear force on the bearing surface. The Couette contribution drags
            // the bearing in +e_theta; pressure-gradient shear opposes dp/ds.
            const double dp_dtheta = (p(i + 1, j) - p(i - 1, j)) / (2.0 * d_theta);
            
            // Poiseuille part is only active when full film
            double g_f = FullFilm::gate(fields, cfg, i, j);
            const double mu_val = dynamic_viscosity(fields, cfg, i, j);
            // STRIATED: only the liquid fraction theta carries Couette shear in cavitated cells.
            const double shear_w = (cfg.cavitated_film_model == CavitatedFilmModel::STRIATED && g_f <= 0.0)
                ? std::max(theta(i, j), 0.0) : 1.0;

            double couette_shear = mu_val * omega * R / h_val;
            double shaft_pressure_arm = 0.5 * h_val;
            if (CrossFilm::active(cfg) && fields.has("cross_film_f0") &&
                fields.has("cross_film_f1") && fields["cross_film_f0"](i, j) > 0.0) {
                const double f0 = fields["cross_film_f0"](i, j);
                const double f1 = fields["cross_film_f1"](i, j);
                couette_shear = omega * R / f0;
                shaft_pressure_arm = f1 / f0;
            }
            const double tau_theta_bearing = couette_shear * shear_w
                - g_f * ((h_val - shaft_pressure_arm) / R) * dp_dtheta;

            double dp_dz = 0.0;
            if (j == 0) {
                const double p_s = (cfg.bc_z_south_type == BCType::DIRICHLET)
                    ? solved_boundary_pressure(cfg, cfg.bc_z_south_val)
                    : p(i, 0);
                dp_dz = (p(i, 0) - p_s) / (0.5 * d_z);
            } else if (j == n_z - 1) {
                const double p_n = (cfg.bc_z_north_type == BCType::DIRICHLET)
                    ? solved_boundary_pressure(cfg, cfg.bc_z_north_val)
                    : p(i, n_z - 1);
                dp_dz = (p_n - p(i, n_z - 1)) / (0.5 * d_z);
            } else {
                dp_dz = (p(i, j + 1) - p(i, j - 1)) / (2.0 * d_z);
            }
            const double tau_z_bearing = -g_f * (h_val - shaft_pressure_arm) * dp_dz;

            local_viscous_x += tau_theta_bearing * (-sin_t) * dA_bearing;
            local_viscous_y += tau_theta_bearing * ( cos_t) * dA_bearing;
            local_viscous_z += tau_z_bearing * dA_bearing;

            const double dA_shaft = R * d_theta * d_z;
            const double tau_torque = couette_shear * shear_w
                + g_f * (shaft_pressure_arm / R) * dp_dtheta;
            local_torque += tau_torque * R * dA_shaft;
        }
    }

    // Reduce over all MPI ranks
    double local_vals[7] = {
        local_pressure_x, local_pressure_y, local_pressure_z,
        local_viscous_x, local_viscous_y, local_viscous_z,
        local_torque};
    double global_vals[7] = {};
    MPI_Allreduce(local_vals, global_vals, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    ForceComponents result;
    result.pressure_x = global_vals[0];
    result.pressure_y = global_vals[1];
    result.pressure_z = global_vals[2];
    result.viscous_x = global_vals[3];
    result.viscous_y = global_vals[4];
    result.viscous_z = global_vals[5];
    result.fluid_x = result.pressure_x + result.viscous_x;
    result.fluid_y = result.pressure_y + result.viscous_y;
    result.fluid_z = result.pressure_z + result.viscous_z;
    result.friction_torque = global_vals[6];

    // Fill the scalar fields
    if (fields.has("pressure_force_x")) fields["pressure_force_x"].fill(result.pressure_x);
    if (fields.has("pressure_force_y")) fields["pressure_force_y"].fill(result.pressure_y);
    if (fields.has("pressure_force_z")) fields["pressure_force_z"].fill(result.pressure_z);
    if (fields.has("viscous_force_x")) fields["viscous_force_x"].fill(result.viscous_x);
    if (fields.has("viscous_force_y")) fields["viscous_force_y"].fill(result.viscous_y);
    if (fields.has("viscous_force_z")) fields["viscous_force_z"].fill(result.viscous_z);
    if (fields.has("fluid_force_x")) fields["fluid_force_x"].fill(result.fluid_x);
    if (fields.has("fluid_force_y")) fields["fluid_force_y"].fill(result.fluid_y);
    if (fields.has("fluid_force_z")) fields["fluid_force_z"].fill(result.fluid_z);
    if (fields.has("friction_torque")) fields["friction_torque"].fill(result.friction_torque);

    return result;
}

}  // namespace Reynolds
