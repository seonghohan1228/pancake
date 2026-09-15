#include "thickness_transport.hpp"

#include <algorithm>

#include "surface_operators.hpp"

namespace ThicknessTransport {

void compute_gravity_velocity(Field& u_u, Field& u_v,
                              const Field& h,
                              const StereoMesh& mesh, const Mask& mask,
                              const BubbleConfig& cfg)
{
    const double coeff_base = cfg.rho_l / (3.0 * cfg.mu);

    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j) || mask.is_rim(i, j)) {
                u_u(i, j) = 0.0;
                u_v(i, j) = 0.0;
                continue;
            }
            const double h_val  = h(i, j);
            const double factor = coeff_base * h_val * h_val;
            u_u(i, j) = factor * mesh.g_u[mesh.idx(i, j)];
            u_v(i, j) = factor * mesh.g_v[mesh.idx(i, j)];
        }
    }
}

void compute_face_fluxes(Field& flux_u, Field& flux_v,
                          const Field& u_u, const Field& u_v,
                          const StereoMesh& mesh, const Mask& mask)
{
    // flux_u(i,j) = east-face flux of cell (i,j)
    // Velocity at east face: average of cell (i,j) and ghost/neighbour (i+1,j)
    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            const double u_e = 0.5 * (u_u(i, j) + u_u(i + 1, j));
            flux_u(i, j) = u_e * mesh.face_len_e[mesh.idx(i, j)];

            const double v_n = 0.5 * (u_v(i, j) + u_v(i, j + 1));
            flux_v(i, j) = v_n * mesh.face_len_n[mesh.idx(i, j)];
        }
    }
}

void solve(Field& h,
           const Field& flux_u, const Field& flux_v,
           const Field* evap_rate,
           BubbleLinearSystem& sys,
           const StereoMesh& mesh, const Mask& mask,
           const BubbleConfig& cfg, double dt)
{
    sys.reset();

    // Transient term: d(h)/dt
    SurfaceFVM::ddt(sys, h, dt, mesh, mask);

    // Convection: div(h * u_s) — face fluxes already contain the velocity;
    // divergence operator multiplies by phi (which is h here)
    SurfaceFVM::divergence(sys, flux_u, flux_v, h,
                           cfg.advection_scheme, mesh, mask);

    // Optional evaporation source (negative = film thinning)
    if (evap_rate)
        SurfaceFVM::add_source(sys, *evap_rate, mesh, mask);

    // Rim Dirichlet BC: h = h_rim
    // Build a field containing the rim value everywhere
    Field bc(h.name + "_bc", mesh.n_u_local, mesh.n_v_global, 0, GridLocation::CENTER);
    bc.fill(cfg.h_rim);
    sys.apply_mask(mask, bc);

    // Solve
    sys.solve(h);

    // Post-solve: clamp to h_min and blank masked cells
    for (int i = 0; i < mesh.n_u_local; ++i)
        for (int j = 0; j < mesh.n_v_global; ++j)
            h(i, j) = mask.is_active(i, j)
                    ? std::max(h(i, j), cfg.h_min)
                    : 0.0;
}

}  // namespace ThicknessTransport
