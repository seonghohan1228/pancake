#include "fbns.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <mpi.h>

#include "full_film.hpp"

namespace FBNS {

namespace {

// Stencil slot order used by every 5-point block in this file.
constexpr int kP = 0, kE = 1, kW = 2, kN = 3, kS = 4;

// Perturbation that keeps the Fischer-Burmeister derivatives defined at the
// non-differentiable point a = b = 0 (Woloszynski et al. 2015, sec. 2.1). It is
// applied to the Jacobian arguments only, never to the residual itself, so the
// converged solution is unperturbed. Both arguments are O(1) after scaling.
constexpr double kFisherBurmeisterEpsilon = 1.0e-14;

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

/// Liquid Poiseuille mobility rho_l h^3 / (12 mu) [kg m / (Pa s)] per unit
/// pressure gradient. Unlike the Elrod path there is no switch function and no
/// bulk modulus: cavitated cells sit at p = p_cav by complementarity, so their
/// pressure gradient — and hence their Poiseuille flux — vanishes on its own.
double poiseuille_gamma(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    const double h_val = fields["h"](i, j);
    return base_density(fields, cfg, i, j) * h_val * h_val * h_val /
        (12.0 * dynamic_viscosity(fields, cfg, i, j));
}

/// Couette mass flux [kg/s] through the east face of (i,j), carrying theta = 1.
double couette_face_flux(const Fields& fields, const Mesh& mesh,
                         const SimulationConfig& cfg, int i, int j) {
    const Field& h = fields["h"];
    const double h_e = 0.5 * (h(i, j) + h(i + 1, j));
    const double rho_e = 0.5 * (base_density(fields, cfg, i, j) +
                                base_density(fields, cfg, i + 1, j));
    return rho_e * cfg.omega * mesh.R * h_e * 0.5 * mesh.get_d_z();
}

/// Inlet covering cell (i,j), or nullptr. Mirrors the geometric membership test
/// in Reynolds::apply_inlet_conditions, including its "last matching inlet wins"
/// overwrite order, so the two solvers pin identical cells.
const InletConfig* inlet_at(const Mesh& mesh, const SimulationConfig& cfg, int i, int j) {
    const InletConfig* found = nullptr;
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
            if (dist_sq <= inlet.size * inlet.size) found = &inlet;
        } else if (dth <= 0.5 * size_rad) {
            found = &inlet;
        }
    }
    return found;
}

/// Pressure imposed on an axial Dirichlet boundary. Config validation already
/// requires bc values >= p_cav; the max is a guard, not a model choice.
double axial_boundary_pressure(const SimulationConfig& cfg, double bc_val) {
    return std::max(bc_val, cfg.effective_p_cav());
}

/// Signed perturbation of a Fischer-Burmeister argument for the Jacobian only.
double fb_aux(double value) {
    if (std::abs(value) >= kFisherBurmeisterEpsilon) return value;
    return value >= 0.0 ? kFisherBurmeisterEpsilon : -kFisherBurmeisterEpsilon;
}

double global_max(double local) {
    double global = local;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

double global_min(double local) {
    double global = local;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global;
}

}  // namespace

bool active(const SimulationConfig& cfg) {
    return cfg.cavitation_model == CavitationModel::FBNS;
}

double pressure_scale(const SimulationConfig& cfg) {
    if (cfg.fbns_pressure_scale > 0.0) return cfg.fbns_pressure_scale;

    // Couette pressure scale 6 mu U R / h_ref^2 with U = omega R and h_ref = c.
    // This is the peak hydrodynamic pressure to within an O(1) factor, which is
    // all the Fischer-Burmeister conditioning needs.
    const double mu_ref = std::max(cfg.mu, 1.0e-30);
    const double c_ref = std::max(cfg.c, 1.0e-30);
    double p_ref = 6.0 * mu_ref * std::abs(cfg.omega) * cfg.R * cfg.R / (c_ref * c_ref);

    const double p_cav = cfg.effective_p_cav();
    for (const auto& inlet : cfg.inlets) {
        p_ref = std::max(p_ref, inlet.p_supply - p_cav);
    }
    if (cfg.bc_z_south_type == BCType::DIRICHLET) {
        p_ref = std::max(p_ref, cfg.bc_z_south_val - p_cav);
    }
    if (cfg.bc_z_north_type == BCType::DIRICHLET) {
        p_ref = std::max(p_ref, cfg.bc_z_north_val - p_cav);
    }
    return std::max(p_ref, 1.0);
}

double cell_outflow(const Fields& fields, const Mesh& mesh,
                    const SimulationConfig& cfg, double dt, int i, int j) {
    const Field& pressure = fields["pressure"];
    const Field& theta = fields["theta"];
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double ew_scale = d_z / (R * d_theta);
    const double ns_scale = R * d_theta / d_z;
    const double bnd_scale = R * d_theta / (0.5 * d_z);

    auto gamma_p = [&](int ii, int jj) { return poiseuille_gamma(fields, cfg, ii, jj); };

    double outflow = 0.0;

    // Poiseuille diffusion through the four faces. Axial Dirichlet faces use the
    // half-cell distance; NEUMANN faces carry no flux.
    outflow += harmonic_avg(gamma_p(i, j), gamma_p(i + 1, j)) * ew_scale *
        (pressure(i, j) - pressure(i + 1, j));
    outflow += harmonic_avg(gamma_p(i - 1, j), gamma_p(i, j)) * ew_scale *
        (pressure(i, j) - pressure(i - 1, j));
    if (j + 1 < n_z) {
        outflow += harmonic_avg(gamma_p(i, j), gamma_p(i, j + 1)) * ns_scale *
            (pressure(i, j) - pressure(i, j + 1));
    } else if (cfg.bc_z_north_type == BCType::DIRICHLET) {
        outflow += gamma_p(i, j) * bnd_scale *
            (pressure(i, j) - axial_boundary_pressure(cfg, cfg.bc_z_north_val));
    }
    if (j > 0) {
        outflow += harmonic_avg(gamma_p(i, j - 1), gamma_p(i, j)) * ns_scale *
            (pressure(i, j) - pressure(i, j - 1));
    } else if (cfg.bc_z_south_type == BCType::DIRICHLET) {
        outflow += gamma_p(i, j) * bnd_scale *
            (pressure(i, j) - axial_boundary_pressure(cfg, cfg.bc_z_south_val));
    }

    // Couette transport of the film content, first-order upwind on the face flux
    // sign (the sliding direction), as in the FBNS discretization.
    const double flux_e = couette_face_flux(fields, mesh, cfg, i, j);
    const double flux_w = couette_face_flux(fields, mesh, cfg, i - 1, j);
    outflow += flux_e * ((flux_e >= 0.0) ? theta(i, j) : theta(i + 1, j));
    outflow -= flux_w * ((flux_w >= 0.0) ? theta(i - 1, j) : theta(i, j));

    // Backward-Euler storage of the liquid content plus the explicit squeeze term.
    if (dt > 0.0) {
        const double rho_l = base_density(fields, cfg, i, j);
        const double vol = mesh.cell_volume();
        outflow += rho_l * fields["h"](i, j) * vol * (theta(i, j) - theta.old(i, j)) / dt;
        if (fields.has("dh_dt")) {
            outflow += rho_l * theta.old(i, j) * fields["dh_dt"](i, j) * vol;
        }
    }
    return outflow;
}

Stats solve(Fields& fields, LinearSystem& sys, Communicator& comm,
            const Mesh& mesh, const SimulationConfig& cfg)
{
    Stats stats;
    Field& pressure = fields["pressure"];
    Field& theta = fields["theta"];
    Field& rho = fields["rho"];

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const int n_local = n_theta * n_z;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double ew_scale = d_z / (R * d_theta);
    const double ns_scale = R * d_theta / d_z;
    const double bnd_scale = R * d_theta / (0.5 * d_z);
    const double vol = mesh.cell_volume();

    const double p_cav = cfg.effective_p_cav();
    const double p_ref = pressure_scale(cfg);
    const double dt = (cfg.solution_mode == SolutionMode::STEADY_STATE) ? 0.0 : cfg.dt;
    const double theta_floor = std::clamp(cfg.theta_min, 0.0, 1.0);
    const double relaxation = std::clamp(cfg.fbns_relaxation, 0.0, 1.0);

    auto index = [n_z](int i, int j) { return i * n_z + j; };
    auto gamma_p = [&](int ii, int jj) { return poiseuille_gamma(fields, cfg, ii, jj); };

    // Full-film ceiling of a cell: 1 without a resolved gas phase, 1 - alpha_gas
    // with one. The complementarity variable is the cavity fraction measured from
    // that ceiling, alpha = theta_full - theta, so the second JFO condition reads
    // theta_full - theta >= 0 rather than 1 - theta >= 0.
    //
    // alpha_gas is a lagged coefficient here: the gas path runs after the flow
    // solve in the outer loop, so theta_full is constant through this Newton
    // solve, d(alpha)/d(theta) = -1 as before, and the Jacobian blocks A and B
    // below are unchanged. With no gas phase every expression collapses to the
    // literal 1 it replaced, bit for bit.
    auto theta_full = [&](int ii, int jj) { return FullFilm::bound(fields, cfg, ii, jj); };
    auto alpha_ceiling = [&](int ii, int jj) {
        return std::max(theta_full(ii, jj) - theta_floor, 0.0);
    };

    // --- Jacobian blocks -----------------------------------------------------
    // J_G,p = A (Poiseuille) and J_G,alpha = B (Couette + storage), with the
    // cavity fraction alpha = 1 - theta as the second unknown so both
    // Fischer-Burmeister diagonals come out non-negative. Both blocks are
    // constant through the Newton iteration: the mobility, the face fluxes and
    // the storage weight are built once, here, from the film geometry and the
    // per-cell "mu" / "rho_liquid_solution" fields as the caller last left them,
    // and the upwind direction follows the face flux sign rather than the
    // solution. The mass balance is therefore linear in (p, alpha) and every
    // Newton nonlinearity lives in the FB equation.
    //
    // Under a thermal or piezo-viscous run mu varies cell to cell and changes
    // between calls; freezing it inside one Newton solve is the deliberate
    // Picard split (see the header). It keeps this Jacobian exactly consistent
    // with the residual cell_outflow evaluates -- both read the same frozen
    // gamma -- which is what preserves the quadratic Newton rate. A Jacobian
    // that carried d(mu)/dT terms while the residual did not would be the
    // limping case.
    std::array<std::vector<double>, 5> A, B;
    for (auto& slot : A) slot.assign(n_local, 0.0);
    for (auto& slot : B) slot.assign(n_local, 0.0);
    std::vector<char> is_inlet(n_local, 0);
    std::vector<double> inlet_pstar(n_local, 0.0);

    double row_scale_max = 0.0;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const int c = index(i, j);

            const double coeff_e = harmonic_avg(gamma_p(i, j), gamma_p(i + 1, j)) * ew_scale;
            const double coeff_w = harmonic_avg(gamma_p(i - 1, j), gamma_p(i, j)) * ew_scale;
            A[kE][c] = -coeff_e;
            A[kW][c] = -coeff_w;
            double diag = coeff_e + coeff_w;
            if (j + 1 < n_z) {
                const double coeff_n = harmonic_avg(gamma_p(i, j), gamma_p(i, j + 1)) * ns_scale;
                A[kN][c] = -coeff_n;
                diag += coeff_n;
            } else if (cfg.bc_z_north_type == BCType::DIRICHLET) {
                diag += gamma_p(i, j) * bnd_scale;
            }
            if (j > 0) {
                const double coeff_s = harmonic_avg(gamma_p(i, j - 1), gamma_p(i, j)) * ns_scale;
                A[kS][c] = -coeff_s;
                diag += coeff_s;
            } else if (cfg.bc_z_south_type == BCType::DIRICHLET) {
                diag += gamma_p(i, j) * bnd_scale;
            }
            A[kP][c] = diag;

            // d(outflow)/d(alpha) = -d(outflow)/d(theta) for the upwind Couette
            // stencil (and for the storage term).
            const double flux_e = couette_face_flux(fields, mesh, cfg, i, j);
            const double flux_w = couette_face_flux(fields, mesh, cfg, i - 1, j);
            B[kP][c] = -(std::max(flux_e, 0.0) + std::max(-flux_w, 0.0));
            B[kE][c] = std::max(-flux_e, 0.0);
            B[kW][c] = std::max(flux_w, 0.0);
            if (dt > 0.0) {
                B[kP][c] -= base_density(fields, cfg, i, j) * fields["h"](i, j) * vol / dt;
            }

            double row_scale = 0.0;
            for (int k = 0; k < 5; ++k) {
                row_scale += std::abs(A[k][c]) * p_ref + std::abs(B[k][c]);
            }
            row_scale_max = std::max(row_scale_max, row_scale);

            const InletConfig* inlet = inlet_at(mesh, cfg, i, j);
            if (inlet != nullptr) {
                is_inlet[c] = 1;
                inlet_pstar[c] = (inlet->p_supply - p_cav) / p_ref;
            }
        }
    }

    // Row scale of the mass balance [kg/s]. Only normalizes the residual (and
    // therefore fbns_tol); it cannot change the Newton solution.
    const double g_ref = std::max(global_max(row_scale_max), 1.0e-30);

    // Nondimensionalize: unknowns become p* = (p - p_cav)/p_ref and alpha, and
    // the mass balance is measured in units of g_ref. Feeding raw pascals into
    // FB(a, b) = a + b - sqrt(a^2 + b^2) against a film content in [0, 1] would
    // make the pressure branch dominate by orders of magnitude and collapse the
    // column swap below.
    for (int c = 0; c < n_local; ++c) {
        for (int k = 0; k < 5; ++k) {
            A[k][c] *= p_ref / g_ref;
            B[k][c] /= g_ref;
        }
    }

    // Inlet cells carry an internal Dirichlet pressure. Both rows are replaced,
    // exactly as the FBNS boundary treatment prescribes: the mass balance becomes
    // p* = p*_supply, and the complementarity row becomes alpha = 0 with
    // (J_F,p, J_F,alpha) = (0, 1). Forcing that pair also forces the column swap,
    // which is what keeps the Schur complement non-singular on these rows.
    for (int c = 0; c < n_local; ++c) {
        if (!is_inlet[c]) continue;
        for (int k = 0; k < 5; ++k) { A[k][c] = 0.0; B[k][c] = 0.0; }
        A[kP][c] = 1.0;
    }

    // --- Initial guess -------------------------------------------------------
    // Cold start (no cavitation information in theta yet): solve the full-film
    // problem G(p, alpha = 0) = 0 and clamp to p >= p_cav, per Woloszynski et al.
    // A warm start from the previous timestep is kept when it exists.
    //
    // The literal 1 here is deliberate and must NOT become FullFilm::bound. It is
    // a sentinel test against the value main.cpp initializes theta to, not a
    // full-film test. Under a moved ceiling theta = 1 - alpha_gas is an ordinary
    // *converged* uncavitated state, so testing theta >= theta_full would declare
    // a cold start on every gassy step, throw the warm start away, and refill
    // theta down to the new ceiling in one go -- a spurious d(theta)/dt sink that
    // drove the pressure to O(10^3) bar on the R290/PZ68 case. Testing against
    // the sentinel keeps the warm start and leaves the ceiling change to Newton,
    // where it belongs.
    double theta_min_local = 1.0e300;
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j)
            theta_min_local = std::min(theta_min_local, theta(i, j));
    const bool cold_start = global_min(theta_min_local) >= 1.0;

    Field delta_b_field("fbns_delta_b", mesh);
    Field swap_flag("fbns_swap", mesh);
    Field schur_d("fbns_schur_d", mesh);
    Field schur_f("fbns_schur_f", mesh);
    Field pstar_field("fbns_pstar", mesh);

    std::vector<double> residual_g(n_local, 0.0);
    std::vector<double> residual_fb(n_local, 0.0);
    std::vector<double> diag_af(n_local, 1.0);
    std::vector<double> diag_bf(n_local, 0.0);

    if (cold_start) {
        theta.fill(1.0);
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                theta(i, j) = theta_full(i, j);
        comm.update_ghosts(theta);
        comm.update_ghosts(pressure);
        for (int i = -pressure.n_ghost; i < n_theta + pressure.n_ghost; ++i)
            for (int j = 0; j < n_z; ++j)
                pstar_field(i, j) = (pressure(i, j) - p_cav) / p_ref;

        // A p*_new = A p*_cur - G(p*_cur, alpha = 0): one linear solve, exact
        // because the mass balance is linear in p*.
        sys.reset();
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const int c = index(i, j);
                double lhs = A[kP][c] * pstar_field(i, j) +
                             A[kE][c] * pstar_field(i + 1, j) +
                             A[kW][c] * pstar_field(i - 1, j);
                if (j + 1 < n_z) lhs += A[kN][c] * pstar_field(i, j + 1);
                if (j > 0)       lhs += A[kS][c] * pstar_field(i, j - 1);
                const double g_value = is_inlet[c]
                    ? (pstar_field(i, j) - inlet_pstar[c])
                    : cell_outflow(fields, mesh, cfg, dt, i, j) / g_ref;

                sys.a_p(i, j) = A[kP][c];
                sys.a_e(i, j) = -A[kE][c];
                sys.a_w(i, j) = -A[kW][c];
                sys.a_n(i, j) = -A[kN][c];
                sys.a_s(i, j) = -A[kS][c];
                sys.source(i, j) = lhs - g_value;
            }
        }
        sys.solve(pstar_field, cfg.linear_rtol, 1000, "fbns full-film initial guess");
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                pressure(i, j) = p_cav + std::max(pstar_field(i, j), 0.0) * p_ref;
    }

    // --- Newton iteration ----------------------------------------------------
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int iter = 0; iter <= cfg.fbns_max_iters; ++iter) {
        comm.update_ghosts(pressure);
        comm.update_ghosts(theta);

        double max_g = 0.0;
        double max_fb = 0.0;
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const int c = index(i, j);
                const double p_star = (pressure(i, j) - p_cav) / p_ref;
                const double alpha = theta_full(i, j) - theta(i, j);

                residual_g[c] = is_inlet[c]
                    ? (p_star - inlet_pstar[c])
                    : cell_outflow(fields, mesh, cfg, dt, i, j) / g_ref;

                double j_fb_p, j_fb_alpha;
                if (is_inlet[c]) {
                    residual_fb[c] = alpha;
                    j_fb_p = 0.0;
                    j_fb_alpha = 1.0;
                } else {
                    residual_fb[c] = p_star + alpha -
                        std::sqrt(p_star * p_star + alpha * alpha);
                    const double a_aux = fb_aux(p_star);
                    const double b_aux = fb_aux(alpha);
                    const double norm_aux = std::sqrt(a_aux * a_aux + b_aux * b_aux);
                    j_fb_p = 1.0 - a_aux / norm_aux;
                    j_fb_alpha = 1.0 - b_aux / norm_aux;
                }

                // Column swap: the larger FB derivative becomes the diagonal of
                // A_F, which bounds cond(A_F) <= 4 + 2 sqrt(2) (Woloszynski et
                // al., Theorem 1) and leaves B_F -> 0 at the solution.
                const bool swapped = j_fb_p < j_fb_alpha;
                diag_af[c] = std::max(swapped ? j_fb_alpha : j_fb_p, 1.0e-12);
                diag_bf[c] = swapped ? j_fb_p : j_fb_alpha;
                swap_flag(i, j) = swapped ? 1.0 : 0.0;
                schur_d(i, j) = diag_bf[c] / diag_af[c];
                schur_f(i, j) = residual_fb[c] / diag_af[c];

                max_g = std::max(max_g, std::abs(residual_g[c]));
                max_fb = std::max(max_fb, std::abs(residual_fb[c]));
            }
        }

        double norms[2] = {max_g, max_fb};
        MPI_Allreduce(MPI_IN_PLACE, norms, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        stats.residual_mass = norms[0];
        stats.residual_complementarity = norms[1];
        stats.residual = std::max(norms[0], norms[1]);
        stats.residual_history.push_back(stats.residual);

        if (cfg.log_outer_iters && rank == 0) {
            std::fprintf(stderr, "fbns iter=%d |G|=%.6e |FB|=%.6e\n",
                         iter, stats.residual_mass, stats.residual_complementarity);
        }

        if (stats.residual < cfg.fbns_tol) {
            stats.converged = true;
            break;
        }
        if (iter == cfg.fbns_max_iters) break;

        comm.update_ghosts(swap_flag);
        comm.update_ghosts(schur_d);
        comm.update_ghosts(schur_f);

        // --- Schur complement ------------------------------------------------
        // With the swap applied column by column, A_G and B_G are column
        // mixtures of A and B that keep the same 5-point sparsity. The reduced
        // system is (B_G - A_G A_F^-1 B_F) db = -G + A_G A_F^-1 F.
        sys.reset();
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const int c = index(i, j);
                const int nb_i[5] = {i, i + 1, i - 1, i, i};
                const int nb_j[5] = {j, j, j, j + 1, j - 1};
                double schur[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
                double rhs = -residual_g[c];

                for (int k = 0; k < 5; ++k) {
                    if (k == kN && j + 1 >= n_z) continue;  // no north cell to couple to
                    if (k == kS && j <= 0) continue;
                    const int ii = nb_i[k];
                    const int jj = nb_j[k];
                    const bool swapped_col = swap_flag(ii, jj) > 0.5;
                    const double a_g = swapped_col ? B[k][c] : A[k][c];
                    const double b_g = swapped_col ? A[k][c] : B[k][c];
                    schur[k] = b_g - a_g * schur_d(ii, jj);
                    rhs += a_g * schur_f(ii, jj);
                }

                sys.a_p(i, j) = schur[kP];
                sys.a_e(i, j) = -schur[kE];
                sys.a_w(i, j) = -schur[kW];
                sys.a_n(i, j) = -schur[kN];
                sys.a_s(i, j) = -schur[kS];
                sys.source(i, j) = rhs;
            }
        }
        sys.solve(delta_b_field, cfg.linear_rtol, 1000, "fbns schur complement");

        // --- Update ----------------------------------------------------------
        double max_theta_change = 0.0;
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const int c = index(i, j);
                const double delta_b = delta_b_field(i, j);
                const double delta_a =
                    (-residual_fb[c] - diag_bf[c] * delta_b) / diag_af[c];
                const bool swapped = swap_flag(i, j) > 0.5;
                const double delta_p_star = swapped ? delta_b : delta_a;
                const double delta_alpha = swapped ? delta_a : delta_b;

                const double tf = theta_full(i, j);
                const double p_star = (pressure(i, j) - p_cav) / p_ref;
                const double alpha = tf - theta(i, j);
                const double p_star_new = std::max(p_star + relaxation * delta_p_star, 0.0);
                const double alpha_new =
                    std::clamp(alpha + relaxation * delta_alpha, 0.0, alpha_ceiling(i, j));

                pressure(i, j) = p_cav + p_star_new * p_ref;
                const double theta_new = tf - alpha_new;
                max_theta_change = std::max(max_theta_change, std::abs(theta_new - theta(i, j)));
                theta(i, j) = theta_new;
            }
        }
        stats.max_theta_change = global_max(max_theta_change);
        stats.newton_iters = iter + 1;
    }

    comm.update_ghosts(pressure);
    comm.update_ghosts(theta);

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            rho(i, j) = base_density(fields, cfg, i, j) * theta(i, j);
        }
    }
    return stats;
}

}  // namespace FBNS
