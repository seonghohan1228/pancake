#include "reynolds.hpp"

#include <algorithm>
#include <cmath>

#include "equation_of_state.hpp"
#include "fvm.hpp"

namespace Reynolds {

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

    // Dirichlet BC: p = p_cav at z = 0 (south) and z = L (north).
    // The laplacian skips these boundary faces; add them here using half-cell distance.
    for (int i = 0; i < n_theta; ++i) {
        const double coeff_s = gamma_f(i, 0)       * R * d_theta / (0.5 * d_z);
        const double coeff_n = gamma_f(i, n_z - 1) * R * d_theta / (0.5 * d_z);
        sys.a_p(i, 0)        += coeff_s;
        sys.source(i, 0)     += coeff_s * cfg.p_cav;
        sys.a_p(i, n_z - 1)  += coeff_n;
        sys.source(i, n_z - 1) += coeff_n * cfg.p_cav;
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

    sys.solve(pressure);

    // Gumbel cavitation condition: enforce p >= p_cav
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j)
            pressure(i, j) = std::max(pressure(i, j), cfg.p_cav);

    // Update density via barotropic EOS (rho = rho_0 * theta(p))
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j)
            rho(i, j) = cfg.rho * EOS::theta_from_pressure(
                pressure(i, j), cfg.p_cav, cfg.bulk_modulus);
}

}  // namespace Reynolds
