#include "fvm.hpp"

#include <algorithm>
#include <cmath>

namespace {

double harmonic_avg(double a, double b) {
    if (a + b == 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

double van_leer(double r) {
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

double limiter(double r, ConvectionScheme scheme) {
    return scheme == ConvectionScheme::VANLEER ? van_leer(r) : 0.0;
}

}  // namespace

namespace FVM {

void ddt(LinearSystem& sys, const Field& phi, double dt, const Mesh& mesh) {
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double vol     = mesh.R * d_theta * d_z;
    const double coeff   = vol / dt;

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            sys.a_p(i, j)    += coeff;
            sys.source(i, j) += coeff * phi.old(i, j);
        }
    }
}

void ddt_weighted(LinearSystem& sys, const Field& phi, const Field& weight,
                  double dt, const Mesh& mesh)
{
    const double base_vol = mesh.R * mesh.get_d_theta() * mesh.get_d_z();

    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double coeff   = weight(i, j) * base_vol / dt;
            sys.a_p(i, j)    += coeff;
            sys.source(i, j) += coeff * phi.old(i, j);
        }
    }
}

void laplacian(LinearSystem& sys, const Field& gamma, const Mesh& mesh) {
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double R       = mesh.R;

    // coeff for east/west faces: gamma_f * face_length / arc_dist
    //   face_length = d_z, arc dist between cell centres = R * d_theta
    //   → gamma_f * d_z / (R * d_theta)
    const double ew_scale = d_z / (R * d_theta);

    // coeff for north/south faces: face_area_z / d_z
    //   face_area = R * d_theta, dist = d_z
    //   → gamma_f * R * d_theta / d_z
    const double ns_scale = R * d_theta / d_z;

    const int n_theta = mesh.n_theta_local;
    const int n_z     = mesh.n_z_local;

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            // East face
            double ge = harmonic_avg(gamma(i, j), gamma(i + 1, j));
            double ce = ge * ew_scale;
            sys.a_e(i, j) += ce;
            sys.a_p(i, j) += ce;

            // West face
            double gw = harmonic_avg(gamma(i - 1, j), gamma(i, j));
            double cw = gw * ew_scale;
            sys.a_w(i, j) += cw;
            sys.a_p(i, j) += cw;

            // North face (skip at top z-boundary)
            if (j + 1 < n_z) {
                double gn = harmonic_avg(gamma(i, j), gamma(i, j + 1));
                double cn = gn * ns_scale;
                sys.a_n(i, j) += cn;
                sys.a_p(i, j) += cn;
            }

            // South face (skip at bottom z-boundary)
            if (j > 0) {
                double gs = harmonic_avg(gamma(i, j - 1), gamma(i, j));
                double cs = gs * ns_scale;
                sys.a_s(i, j) += cs;
                sys.a_p(i, j) += cs;
            }
        }
    }
}

void divergence(LinearSystem& sys,
                const Field& flux_theta, const Field& flux_z, const Field& phi,
                ConvectionScheme scheme, const Mesh& mesh,
                const Field* central_mask_theta, const Field* central_mask_z)
{
    const int n_theta = mesh.n_theta_local;
    const int n_z     = mesh.n_z_local;
    const bool is_tvd = (scheme == ConvectionScheme::VANLEER);
    const bool is_linear = (scheme == ConvectionScheme::LINEAR);
    const bool is_type = (scheme == ConvectionScheme::TYPE_DIFFERENCING &&
                          central_mask_theta != nullptr);

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            // Fluxes at cell faces (outward positive)
            // flux_theta(i,j) is at the east face of (i,j)
            double F_e = flux_theta(i, j);
            double F_w = flux_theta(i - 1, j);  // east face of western neighbor
            double F_n = (j + 1 < n_z) ? flux_z(i, j)     : 0.0;
            double F_s = (j > 0)        ? flux_z(i, j - 1) : 0.0;

            // Upwind implicit stencil
            sys.a_p(i, j) += std::max(F_e, 0.0) + std::max(-F_w, 0.0)
                           +  std::max(F_n, 0.0) + std::max(-F_s, 0.0);
            sys.a_e(i, j) += std::max(-F_e, 0.0);
            sys.a_w(i, j) += std::max( F_w, 0.0);
            if (j + 1 < n_z) sys.a_n(i, j) += std::max(-F_n, 0.0);
            if (j > 0)        sys.a_s(i, j) += std::max( F_s, 0.0);

            if (is_type || is_linear) {
                // Deferred central correction. LINEAR applies pure central (2nd-order,
                // unbounded) on every face; TYPE_DIFFERENCING (Vijayaraghavan & Keith
                // 1989) applies it only on full-film faces, leaving cavitated faces
                // pure upwind.
                auto theta_central = [&](int ii) {
                    return is_linear ||
                           (central_mask_theta != nullptr && (*central_mask_theta)(ii, j) > 0.5);
                };
                auto z_central = [&](int jj) {
                    return is_linear ||
                           (central_mask_z != nullptr && (*central_mask_z)(i, jj) > 0.5);
                };
                if (theta_central(i)) {
                    const double phi_e_up = (F_e >= 0.0) ? phi.old(i, j) : phi.old(i + 1, j);
                    const double phi_e_c  = 0.5 * (phi.old(i, j) + phi.old(i + 1, j));
                    sys.source(i, j) -= F_e * (phi_e_c - phi_e_up);
                }
                if (theta_central(i - 1)) {
                    const double phi_w_up = (F_w >= 0.0) ? phi.old(i - 1, j) : phi.old(i, j);
                    const double phi_w_c  = 0.5 * (phi.old(i - 1, j) + phi.old(i, j));
                    sys.source(i, j) += F_w * (phi_w_c - phi_w_up);
                }
                if (j + 1 < n_z && z_central(j)) {
                    const double phi_n_up = (F_n >= 0.0) ? phi.old(i, j) : phi.old(i, j + 1);
                    const double phi_n_c  = 0.5 * (phi.old(i, j) + phi.old(i, j + 1));
                    sys.source(i, j) -= F_n * (phi_n_c - phi_n_up);
                }
                if (j > 0 && z_central(j - 1)) {
                    const double phi_s_up = (F_s >= 0.0) ? phi.old(i, j - 1) : phi.old(i, j);
                    const double phi_s_c  = 0.5 * (phi.old(i, j - 1) + phi.old(i, j));
                    sys.source(i, j) += F_s * (phi_s_c - phi_s_up);
                }
                continue;
            }

            if (!is_tvd) continue;

            // TVD deferred correction: explicit high-order fix added to source.
            // Uses phi.old (previous timestep / last outer iteration).
            // East face
            double phi_e_up = (F_e >= 0.0) ? phi.old(i, j) : phi.old(i + 1, j);
            {
                double denom = phi.old(i, j) - phi.old(i - 1, j);
                double num   = phi.old(i + 1, j) - phi.old(i, j);
                double r     = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
                if (F_e < 0.0) r = 1.0 / (r + 1e-14);
                double phi_e_ho = phi_e_up + 0.5 * limiter(r, scheme) * (phi_e_up - phi.old(i - 1, j));
                sys.source(i, j) -= F_e * (phi_e_ho - phi_e_up);
            }

            // West face
            double phi_w_up = (F_w >= 0.0) ? phi.old(i - 1, j) : phi.old(i, j);
            {
                double denom = phi.old(i - 1, j) - phi.old(i - 2, j);
                double num   = phi.old(i, j) - phi.old(i - 1, j);
                double r     = (std::abs(denom) > 1e-14) ? num / denom : 0.0;
                if (F_w < 0.0) r = 1.0 / (r + 1e-14);
                double phi_w_ho = phi_w_up + 0.5 * limiter(r, scheme) * (phi_w_up - phi.old(i - 2, j));
                sys.source(i, j) += F_w * (phi_w_ho - phi_w_up);
            }
        }
    }
}

void add_source(LinearSystem& sys, const Field& source_field, const Mesh& mesh) {
    const double d_theta = mesh.get_d_theta();
    const double d_z     = mesh.get_d_z();
    const double vol     = mesh.R * d_theta * d_z;

    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            sys.source(i, j) += source_field(i, j) * vol;
}

}  // namespace FVM
