#include "surface_operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

double harmonic(double a, double b) {
    if (a + b == 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

double van_leer(double r) { return (r + std::abs(r)) / (1.0 + std::abs(r)); }
double minmod(double r)   { return std::max(0.0, std::min(1.0, r)); }

double limiter(double r, ConvectionScheme scheme) {
    switch (scheme) {
        case ConvectionScheme::TVD_VANLEER: return van_leer(r);
        case ConvectionScheme::TVD_MINMOD:  return minmod(r);
        default: return 0.0;
    }
}

// Returns true if a cell at global position (gu, j) lies within the circular mask.
bool active_at(const StereoMesh& mesh, int gu, int j) {
    double uc = -mesh.L_box + (gu + 0.5) * mesh.du;
    double vc = -mesh.L_box + (j  + 0.5) * mesh.dv;
    return uc*uc + vc*vc <= mesh.R * mesh.R;
}

}  // namespace

namespace SurfaceFVM {

void ddt(BubbleLinearSystem& sys, const Field& phi, double dt,
         const StereoMesh& mesh, const Mask& mask)
{
    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;
            const double coeff = mesh.cell_area[mesh.idx(i, j)] / dt;
            sys.a_p(i, j)    += coeff;
            sys.source(i, j) += coeff * phi.old(i, j);
        }
    }
}

void ddt_weighted(BubbleLinearSystem& sys, const Field& phi, const Field& weight,
                  double dt, const StereoMesh& mesh, const Mask& mask)
{
    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;
            const double coeff = weight(i, j) * mesh.cell_area[mesh.idx(i, j)] / dt;
            sys.a_p(i, j)    += coeff;
            sys.source(i, j) += coeff * phi.old(i, j);
        }
    }
}

void laplacian(BubbleLinearSystem& sys, const Field& gamma,
               const StereoMesh& mesh, const Mask& mask)
{
    for (int i = 0; i < mesh.n_u_local; ++i) {
        const int gu = mesh.global_u(i);
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;

            // East face: connects (gu,j) → (gu+1,j)
            if (gu + 1 < mesh.n_u_global && active_at(mesh, gu + 1, j)) {
                const double g_f = harmonic(gamma(i, j), gamma(i + 1, j));
                const double L   = mesh.face_len_e[mesh.idx(i, j)];
                const double d   = mesh.d_u[mesh.idx(i, j)];
                const double c   = g_f * L / d;
                sys.a_e(i, j) += c;
                sys.a_p(i, j) += c;
            }

            // West face: connects (gu-1,j) → (gu,j)
            if (gu - 1 >= 0 && active_at(mesh, gu - 1, j)) {
                const double g_f = harmonic(gamma(i - 1, j), gamma(i, j));
                const double L   = mesh.east_face_len(gu - 1, j);   // east face of left cell
                const double d   = mesh.d_u_east(gu - 1, j);
                const double c   = g_f * L / d;
                sys.a_w(i, j) += c;
                sys.a_p(i, j) += c;
            }

            // North face: connects (gu,j) → (gu,j+1)
            if (j + 1 < mesh.n_v_global && active_at(mesh, gu, j + 1)) {
                const double g_f = harmonic(gamma(i, j), gamma(i, j + 1));
                const double L   = mesh.face_len_n[mesh.idx(i, j)];
                const double d   = mesh.d_v[mesh.idx(i, j)];
                const double c   = g_f * L / d;
                sys.a_n(i, j) += c;
                sys.a_p(i, j) += c;
            }

            // South face: connects (gu,j-1) → (gu,j)
            if (j - 1 >= 0 && active_at(mesh, gu, j - 1)) {
                const double g_f = harmonic(gamma(i, j - 1), gamma(i, j));
                // North face of (i, j-1) = south face of (i, j)
                const double L   = mesh.face_len_n[mesh.idx(i, j - 1)];
                const double d   = mesh.d_v[mesh.idx(i, j - 1)];
                const double c   = g_f * L / d;
                sys.a_s(i, j) += c;
                sys.a_p(i, j) += c;
            }
        }
    }
}

void divergence(BubbleLinearSystem& sys,
                const Field& flux_u, const Field& flux_v, const Field& phi,
                ConvectionScheme scheme,
                const StereoMesh& mesh, const Mask& mask)
{
    const bool is_tvd = (scheme == ConvectionScheme::TVD_VANLEER ||
                         scheme == ConvectionScheme::TVD_MINMOD);

    for (int i = 0; i < mesh.n_u_local; ++i) {
        const int gu = mesh.global_u(i);
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;

            // Face fluxes: outward positive.
            // flux_u(i,j)   = east outward flux of (i,j)
            // flux_u(i-1,j) = east outward flux of left neighbour = west inward flux of (i,j)
            const double F_e = flux_u(i,     j);
            const double F_w = flux_u(i - 1, j);  // ghost access for i==0
            const double F_n = (j + 1 < mesh.n_v_global) ? flux_v(i, j)     : 0.0;
            const double F_s = (j > 0)                   ? flux_v(i, j - 1) : 0.0;

            // Upwind implicit
            sys.a_p(i, j) += std::max( F_e, 0.0) + std::max(-F_w, 0.0)
                           +  std::max( F_n, 0.0) + std::max(-F_s, 0.0);
            sys.a_e(i, j) += std::max(-F_e, 0.0);
            sys.a_w(i, j) += std::max( F_w, 0.0);
            if (j + 1 < mesh.n_v_global) sys.a_n(i, j) += std::max(-F_n, 0.0);
            if (j > 0)                   sys.a_s(i, j) += std::max( F_s, 0.0);

            if (!is_tvd) continue;

            // TVD deferred correction using phi.old (same pattern as fvm.cpp)

            // East face
            {
                const double phi_up = (F_e >= 0.0) ? phi.old(i, j) : phi.old(i + 1, j);
                const double denom  = phi.old(i, j) - phi.old(i - 1, j);
                const double num    = phi.old(i + 1, j) - phi.old(i, j);
                double r = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
                if (F_e < 0.0) r = 1.0 / (r + 1e-14);
                const double phi_ho = phi_up + 0.5 * limiter(r, scheme) * (phi_up - phi.old(i - 1, j));
                sys.source(i, j) -= F_e * (phi_ho - phi_up);
            }

            // West face
            {
                const double phi_up = (F_w >= 0.0) ? phi.old(i - 1, j) : phi.old(i, j);
                const double denom  = phi.old(i - 1, j) - phi.old(i - 2, j);
                const double num    = phi.old(i,     j) - phi.old(i - 1, j);
                double r = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
                if (F_w < 0.0) r = 1.0 / (r + 1e-14);
                const double phi_ho = phi_up + 0.5 * limiter(r, scheme) * (phi_up - phi.old(i - 2, j));
                sys.source(i, j) += F_w * (phi_ho - phi_up);
            }
        }
    }
}

void add_source(BubbleLinearSystem& sys, const Field& source_field,
                const StereoMesh& mesh, const Mask& mask)
{
    for (int i = 0; i < mesh.n_u_local; ++i)
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;
            sys.source(i, j) += source_field(i, j) * mesh.cell_area[mesh.idx(i, j)];
        }
}

void gradient(const Field& phi, Field& grad_u, Field& grad_v,
              const StereoMesh& mesh, const Mask& mask)
{
    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) {
                grad_u(i, j) = 0.0;
                grad_v(i, j) = 0.0;
                continue;
            }
            // Central difference; use half-distances for denominator
            const double d_e = mesh.d_u[mesh.idx(i, j)];
            const double d_w = mesh.d_u_east(mesh.global_u(i) - 1, j);
            grad_u(i, j) = (phi(i + 1, j) - phi(i - 1, j)) / (d_e + d_w);

            const double d_n = mesh.d_v[mesh.idx(i, j)];
            const double d_s = (j > 0) ? mesh.d_v[mesh.idx(i, j - 1)] : d_n;
            grad_v(i, j) = (phi(i, j + 1) - phi(i, j - 1)) / (d_n + d_s);
        }
    }
}

}  // namespace SurfaceFVM
