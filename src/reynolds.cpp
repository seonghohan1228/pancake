#include "reynolds.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <mpi.h>

#include "equation_of_state.hpp"
#include "fvm.hpp"

namespace Reynolds {

/// Apply internal Dirichlet boundary conditions for oil supply inlets using the penalty method.
static void apply_inlet_conditions(LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg) {
    if (cfg.inlets.empty()) return;

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;

    for (const auto& inlet : cfg.inlets) {
        const double theta_in_rad = inlet.theta * M_PI / 180.0;
        const double size_rad = inlet.size * M_PI / 180.0;
        
        for (int i = 0; i < n_theta; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            
            // Handle periodic wrap-around for theta distance
            double dth = std::abs(theta_c - theta_in_rad);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            for (int j = 0; j < n_z; ++j) {
                const double z_c = (j + 0.5) * d_z;
                bool inside = false;

                if (inlet.type == InletConfig::Type::CIRCULAR) {
                    const double dist_sq = (R * dth) * (R * dth) + (z_c - inlet.z) * (z_c - inlet.z);
                    if (dist_sq <= inlet.size * inlet.size) inside = true;
                } else if (inlet.type == InletConfig::Type::GROOVE) {
                    if (dth <= 0.5 * size_rad) inside = true;
                }

                if (inside) {
                    double val_sup = inlet.p_supply;
                    if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS) {
                        val_sup = EOS::theta_from_pressure(inlet.p_supply, cfg.p_cav, cfg.bulk_modulus);
                    }
                    
                    const double penalty = 1e20;
                    sys.a_p(i, j) = penalty;
                    sys.source(i, j) = penalty * val_sup;
                }
            }
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

    // Compute gamma = ρh³/(12μ) and ρh over physical cells plus theta ghost layers.
    // h and rho theta-ghost cells must already be synced before calling this function.
    Field gamma_f("gamma",  mesh);
    Field rho_h_f("rho_h",  mesh);
    for (int i = -ng; i < n_theta + ng; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val   = h(i, j);
            const double rho_val = rho(i, j);
            gamma_f(i, j) = rho_val * h_val * h_val * h_val / (12.0 * cfg.mu);
            rho_h_f(i, j) = rho_val * h_val;
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
    // so source_field = -(ω/2 · ∂(ρh)/∂θ).  Centered difference uses ±1 ghost cells.
    Field couette_src("couette_src", mesh);
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            couette_src(i, j) = -0.5 * cfg.omega *
                (rho_h_f(i + 1, j) - rho_h_f(i - 1, j)) / (2.0 * d_theta);
        }
    }
    FVM::add_source(sys, couette_src, mesh);

    apply_inlet_conditions(sys, mesh, cfg);

    sys.solve(pressure);

    // Gumbel cavitation condition: enforce p >= p_cav
    if (cfg.cavitation_model == CavitationModel::GUMBEL) {
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                pressure(i, j) = std::max(pressure(i, j), cfg.p_cav);
    }

    // Update density via barotropic EOS (rho = rho_0 * theta(p))
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j) {
            double theta_p = EOS::theta_from_pressure(pressure(i, j), cfg.p_cav, cfg.bulk_modulus);
            // In Full Sommerfeld, we should ensure density doesn't drop due to negative pressure 
            // if we want to follow classical liquid theory, but we'll follow EOS here.
            rho(i, j) = cfg.rho * theta_p;
            fields["theta"](i, j) = theta_p;
        }
}

void solve_elrod(Fields& fields, LinearSystem& sys,
                 const Mesh& mesh, const SimulationConfig& cfg)
{
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

    // Γ_base = ρ₀βh³/(12μ) — effective diffusion without switch function.
    // h theta-ghost cells must already be synced before calling this function.
    Field Gamma_base("Gamma_base", mesh);
    for (int i = -ng; i < n_theta + ng; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double hv = h(i, j);
            Gamma_base(i, j) = cfg.rho * cfg.bulk_modulus * hv * hv * hv / (12.0 * cfg.mu);
        }
    }

    // Couette convection face flux: F_e(i,j) = ρ₀·(ω/2)·R·h_e·Δz at the east face.
    // Extends to i = -1 so FVM::divergence can access the west-face ghost at i=0.
    Field couette_flux("couette_flux", mesh);
    Field zero_flux_z("zero_flux_z", mesh);
    for (int i = -1; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_e = 0.5 * (h(i, j) + h(i + 1, j));
            couette_flux(i, j) = cfg.rho * cfg.omega * R * h_e * 0.5 * d_z;
            zero_flux_z(i, j)  = 0.0;
        }
    }

    // Γ_g = Γ_base · g(θ): diffusion disabled in cavitated cells.
    Field Gamma_g("Gamma_g", mesh);

    // Buffer for outer-iteration convergence check (flat array, physical cells only).
    const int n_local = n_theta * n_z;
    std::vector<double> theta_prev(n_local);

    // Outer (flag-update) iteration: iterate until the switch pattern converges.
    for (int outer = 0; outer < cfg.max_outer_iters; ++outer) {

        // Save θ for convergence check.
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                theta_prev[i * n_z + j] = theta(i, j);

        // Apply switch function g(θ) to get the effective diffusion coefficient.
        // For the first outer iteration (outer=0), we force g=1 everywhere.
        // This is a "flooded initialization" that prevents numerical starvation
        // by establishing a pressure zone from the oil supply boundaries.
        for (int i = -ng; i < n_theta + ng; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double g = (outer == 0) ? 1.0 : EOS::switch_function(theta(i, j));
                Gamma_g(i, j) = Gamma_base(i, j) * g;
            }
        }

        // Assemble: -div(Γ_g·∇θ)·V + div(F·θ)·V = 0
        // i.e., div(Γ_g·∇θ) = div(F·θ)  (Elrod-Adams with cavitation switch)
        sys.reset();
        FVM::laplacian(sys, Gamma_g, mesh);
        FVM::divergence(sys, couette_flux, zero_flux_z, theta, ConvectionScheme::UPWIND, mesh);

        // Dirichlet z-BC: θ = θ_boundary at z=0 and z=L.
        // Uses Gamma_base (not Gamma_g) to represent the flooded-bearing condition.
        for (int i = 0; i < n_theta; ++i) {
            // South boundary (z=0)
            if (cfg.bc_z_south_type == BCType::DIRICHLET) {
                const double cs = Gamma_base(i, 0) * R * d_theta / (0.5 * d_z);
                sys.a_p(i, 0) += cs;
                double th_s = EOS::theta_from_pressure(cfg.bc_z_south_val, cfg.p_cav, cfg.bulk_modulus);
                sys.source(i, 0) += cs * th_s;
            } else if (cfg.bc_z_south_type == BCType::INLET_OUTLET) {
                // If inflow (u_z > 0 at south face), use fixed theta
                // In solve_elrod, u_z at south face (j=0) is roughly proportional to -(p_0 - p_south)/0.5dz
                // But we can use the velocity field from the previous step/iteration if available.
                // For simplicity, we can assume Dirichlet for now or check current pressure.
                // A better way: check current pressure relative to bc_z_south_val.
                if (fields["pressure"](i, 0) < cfg.bc_z_south_val) { // Inflow
                    const double cs = Gamma_base(i, 0) * R * d_theta / (0.5 * d_z);
                    sys.a_p(i, 0) += cs;
                    sys.source(i, 0) += cs * cfg.bc_z_south_theta;
                }
                // Else: outflow -> zeroGradient (default in FVM::laplacian)
            }

            // North boundary (z=L)
            if (cfg.bc_z_north_type == BCType::DIRICHLET) {
                const double cn = Gamma_base(i, n_z - 1) * R * d_theta / (0.5 * d_z);
                sys.a_p(i, n_z - 1) += cn;
                double th_n = EOS::theta_from_pressure(cfg.bc_z_north_val, cfg.p_cav, cfg.bulk_modulus);
                sys.source(i, n_z - 1) += cn * th_n;
            } else if (cfg.bc_z_north_type == BCType::INLET_OUTLET) {
                if (fields["pressure"](i, n_z - 1) < cfg.bc_z_north_val) { // Inflow
                    const double cn = Gamma_base(i, n_z - 1) * R * d_theta / (0.5 * d_z);
                    sys.a_p(i, n_z - 1) += cn;
                    sys.source(i, n_z - 1) += cn * cfg.bc_z_north_theta;
                }
                // Else: outflow -> zeroGradient
            }
        }

        apply_inlet_conditions(sys, mesh, cfg);

        sys.solve(theta);

        // Clamp θ ≥ θ_min to prevent unphysical negative film content.
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                theta(i, j) = std::max(theta(i, j), cfg.theta_min);

        // Snap near-full-film cells to exactly 1.0. Solver precision (rtol=1e-8,
        // condition ~O(10^3)) produces errors ~O(1e-5) near the free boundary,
        // causing spurious g=0 assignments in the next outer iteration. Snap any
        // cell within snap_tol of 1.0 upward to suppress this numerical chatter.
        // Physically cavitated cells have θ << 1 and are not affected.
        const double snap_tol = 5e-7;
        for (int i = 0; i < n_theta; ++i)
            for (int j = 0; j < n_z; ++j)
                if (theta(i, j) > 1.0 - snap_tol && theta(i, j) < 1.0)
                    theta(i, j) = 1.0;

        // Convergence: max |θ_new − θ_prev| and count of flag (g) changes.
        // Use snap_tol as threshold so near-interface cells do not trigger restarts.
        double max_change = 0.0;
        int    n_switched = 0;
        for (int i = 0; i < n_theta; ++i) {
            for (int j = 0; j < n_z; ++j) {
                const double tp = theta_prev[i * n_z + j];
                const double tn = theta(i, j);
                max_change = std::max(max_change, std::abs(tn - tp));
                if ((tp >= 1.0 - snap_tol) != (tn >= 1.0 - snap_tol)) ++n_switched;
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
        if (max_change < cfg.outer_tol && n_switched == 0) break;
    }

    // Recover pressure and density from EOS.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            pressure(i, j) = EOS::pressure_from_theta(theta(i, j), cfg.p_cav, cfg.bulk_modulus);
            rho(i, j)      = cfg.rho * theta(i, j);
        }
    }
}

void calculate_velocities(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    const Field& p = fields["pressure"];
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    Field& u_theta = fields["velocity_theta"];
    Field& u_z     = fields["velocity_z"];

    const int n_theta = mesh.n_theta_local;
    const int n_z     = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;
    const double omega   = cfg.omega;
    const double mu      = cfg.mu;

    // u_theta at FACE_THETA: interface between (i-1, j) and (i, j)
    // u_theta(i, j) is the west face of cell (i, j)
    // We need to calculate for i=0 to n_theta-1. Face i=n_theta is the same as i=0 for periodic.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            // Face values at west face of cell (i, j)
            double h_f = 0.5 * (h(i-1, j) + h(i, j));
            
            // Poiseuille term active only if both cells are full film
            double g_f = (theta(i-1, j) >= 1.0 && theta(i, j) >= 1.0) ? 1.0 : 0.0;
            
            double dp_dtheta = (p(i, j) - p(i-1, j)) / d_theta;
            u_theta(i, j) = 0.5 * omega * R - g_f * (h_f * h_f) / (12.0 * mu * R) * dp_dtheta;
        }
    }

    // u_z at FACE_Z: interface between (i, j-1) and (i, j)
    // u_z(i, j) is the south face of cell (i, j). 
    // Field u_z has n_z_local + 1 faces in z direction.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j <= n_z; ++j) {
            double p_prev, p_next, h_f, g_f;
            double dz_eff = d_z;

            if (j == 0) { // South boundary face
                p_next = p(i, 0);
                if (cfg.bc_z_south_type == BCType::DIRICHLET) {
                    p_prev = cfg.bc_z_south_val;
                    dz_eff = 0.5 * d_z;
                } else {
                    p_prev = p_next; // Neumann: dp/dz = 0
                }
                h_f = h(i, 0);
                g_f = (theta(i, 0) >= 1.0) ? 1.0 : 0.0;
            } else if (j == n_z) { // North boundary face
                p_prev = p(i, n_z - 1);
                if (cfg.bc_z_north_type == BCType::DIRICHLET) {
                    p_next = cfg.bc_z_north_val;
                    dz_eff = 0.5 * d_z;
                } else {
                    p_next = p_prev; // Neumann: dp/dz = 0
                }
                h_f = h(i, n_z - 1);
                g_f = (theta(i, n_z - 1) >= 1.0) ? 1.0 : 0.0;
            } else { // Interior face
                p_prev = p(i, j - 1);
                p_next = p(i, j);
                h_f = 0.5 * (h(i, j-1) + h(i, j));
                g_f = (theta(i, j-1) >= 1.0 && theta(i, j) >= 1.0) ? 1.0 : 0.0;
            }

            double dp_dz = (p_next - p_prev) / dz_eff;
            u_z(i, j) = - g_f * (h_f * h_f) / (12.0 * mu) * dp_dz;
        }
    }
}

}  // namespace Reynolds
