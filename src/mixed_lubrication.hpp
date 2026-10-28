#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "config.hpp"

/// Mixed lubrication closures for rough journal bearings (load-sharing model).
///
/// Both closures depend on the local film ratio H = h / sigma, with the composite
/// roughness sigma = sqrt(Rq_shaft^2 + Rq_bush^2). Each is switched independently:
///
///   asperity_contact_model = GREENWOOD_TRIPP
///     Asperity contact pressure (Greenwood & Tripp 1970-71)
///       p_c = K E' F_5/2(H),  K = (16 sqrt(2) pi / 15) (eta beta sigma)^2 sqrt(sigma/beta),
///     with 1/E' = (1-nu1^2)/E1 + (1-nu2^2)/E2 and the exact Gaussian integral
///     F_5/2(H) = (2 pi)^-1/2 int_H^inf (s - H)^5/2 exp(-s^2/2) ds, tabulated on
///     0 <= H <= 8 (zero above; F_5/2(8) ~ 1e-17). K may be given directly
///     (asperity_gt_k > 0, as measured surface data report it); otherwise it is
///     built from asperity_eta_beta_sigma and asperity_sigma_over_beta. Asperity shear on the
///     shaft is boundary_friction_coefficient * p_c. Contact load adds to the
///     hydrodynamic load in the force balance; the Reynolds solve is unchanged.
///
///   flow_factor_model = PATIR_CHENG
///     Average Reynolds flow factors for isotropic Gaussian roughness (gamma = 1)
///     from Patir & Cheng (1978, 1979):
///       phi_x = 1 - 0.90 exp(-0.56 H)                         (H >= 0.5)
///       Phi_s = 1.899 H^0.98 exp(-0.92 H + 0.05 H^2)          (H <= 5)
///             = 1.126 exp(-0.25 H)                            (H > 5)
///       phi_s = (Rq_shaft^2 - Rq_bush^2)/sigma^2 Phi_s, with the shaft as moving surface 1
///     and the shear-stress factors phi_f, phi_fs (A3 = 11.1, a4 = 2.31, a5 = 2.38,
///     a6 = 0.11, zero for H > 7) and phi_fp = 1 - 1.40 exp(-0.66 H).
///     The expected gap E[h_T] is taken equal to h (within 0.5% for H >= 2).
///     Valid with cavitation_model = JFO or FBNS only.
///
/// Every function returns the smooth-surface value (1 for factors, 0 for contact)
/// when its model is NONE.
namespace MixedLubrication {

inline bool contact_active(const SimulationConfig& cfg) {
    return cfg.asperity_contact_model == AsperityContactModel::GREENWOOD_TRIPP;
}

inline bool flow_factors_active(const SimulationConfig& cfg) {
    return cfg.flow_factor_model == FlowFactorModel::PATIR_CHENG;
}

/// Composite RMS roughness sigma [m].
inline double composite_roughness(const SimulationConfig& cfg) {
    return std::hypot(cfg.roughness_rq_shaft, cfg.roughness_rq_bush);
}

/// Film ratio H = h / sigma.
inline double film_ratio(double h, const SimulationConfig& cfg) {
    return h / composite_roughness(cfg);
}

/// Exact Greenwood-Tripp integral F_5/2(H) by composite Simpson quadrature
/// (used to build the table; ported from codex/gt-contact 86a45a1).
inline double greenwood_tripp_f52_exact(double H) {
    constexpr double pi = 3.14159265358979323846;
    const double lower = H;
    const double upper = std::max(H, 0.0) + 12.0;
    const int n = 4000;  // even
    const double ds = (upper - lower) / n;
    auto integrand = [H](double s) {
        const double x = s - H;
        return x <= 0.0 ? 0.0 : std::pow(x, 2.5) * std::exp(-0.5 * s * s);
    };
    double sum = integrand(lower) + integrand(upper);
    for (int k = 1; k < n; ++k) sum += (k % 2 ? 4.0 : 2.0) * integrand(lower + k * ds);
    return sum * ds / 3.0 / std::sqrt(2.0 * pi);
}

/// Greenwood-Tripp integral F_5/2(H) for Gaussian summit heights: exact values
/// tabulated on [0, 8] with log-linear interpolation, zero for H >= 8, and H < 0
/// clamped to 0. The common fit 4.4086e-5 (4 - H)^6.804 (Hu et al. 1994) is 11 %
/// low at H = 0, 4x low at H = 3 and zero for H >= 4, i.e. wrong at contact onset.
inline double greenwood_tripp_f52(double H) {
    constexpr int n_table = 1601;
    constexpr double h_max = 8.0;
    static const std::array<double, n_table> log_table = [] {
        std::array<double, n_table> table{};
        for (int k = 0; k < n_table; ++k)
            table[k] = std::log(greenwood_tripp_f52_exact(h_max * k / (n_table - 1)));
        return table;
    }();
    if (!(H < h_max)) return 0.0;
    const double x = std::max(H, 0.0) / h_max * (n_table - 1);
    const int k = std::min(static_cast<int>(x), n_table - 2);
    const double w = x - k;
    return std::exp((1.0 - w) * log_table[k] + w * log_table[k + 1]);
}

/// Greenwood-Tripp roughness parameter K [-]: asperity_gt_k when set (> 0),
/// otherwise (16 sqrt(2) pi / 15) (eta beta sigma)^2 sqrt(sigma / beta).
inline double greenwood_tripp_k(const SimulationConfig& cfg) {
    if (cfg.asperity_gt_k > 0.0) return cfg.asperity_gt_k;
    constexpr double pi = 3.14159265358979323846;
    return 16.0 * std::sqrt(2.0) * pi / 15.0 *
           cfg.asperity_eta_beta_sigma * cfg.asperity_eta_beta_sigma *
           std::sqrt(cfg.asperity_sigma_over_beta);
}

/// Asperity contact pressure p_c [Pa] at local gap h [m].
inline double contact_pressure(double h, const SimulationConfig& cfg) {
    if (!contact_active(cfg)) return 0.0;
    const double sigma = composite_roughness(cfg);
    return greenwood_tripp_k(cfg) * cfg.contact_modulus * greenwood_tripp_f52(std::max(h, 0.0) / sigma);
}

/// Isotropic pressure flow factor phi_x = phi_z; held at its H = 0.5 value below the fit range.
inline double pressure_flow_factor_isotropic(double H) {
    return 1.0 - 0.90 * std::exp(-0.56 * std::max(H, 0.5));
}

/// Isotropic single-surface shear flow factor Phi_s(H).
inline double shear_flow_factor_isotropic(double H) {
    H = std::max(H, 0.5);
    return H <= 5.0 ? 1.899 * std::pow(H, 0.98) * std::exp(-0.92 * H + 0.05 * H * H)
                    : 1.126 * std::exp(-0.25 * H);
}

/// Mean shear stress factor phi_f(H) for the truncated polynomial height
/// distribution, with the contact cut-off eps* = 1/300 for H <= 3.
inline double shear_stress_factor_f(double H) {
    const double z = std::max(H, 1e-6) / 3.0;
    const double a = (1.0 - z * z) * (1.0 - z * z) * (1.0 - z * z);
    if (H > 3.0) {
        return 35.0 / 32.0 * z * (a * std::log((z + 1.0) / (z - 1.0)) +
                                  z / 15.0 * (66.0 + z * z * (30.0 * z * z - 80.0)));
    }
    const double polynomial = -55.0 + z * (132.0 + z * (345.0 + z * (-160.0 +
                              z * (-405.0 + z * (60.0 + 147.0 * z)))));
    return 35.0 / 32.0 * z * (a * std::log((z + 1.0) * 300.0) + polynomial / 60.0);
}

/// Isotropic single-surface roughness shear-stress factor Phi_fs(H).
inline double shear_stress_factor_fs_isotropic(double H) {
    if (H > 7.0) return 0.0;
    H = std::max(H, 0.5);
    return 11.1 * std::pow(H, 2.31) * std::exp(-2.38 * H + 0.11 * H * H);
}

/// Pressure-driven shear stress factor phi_fp(H).
inline double shear_stress_factor_fp(double H) {
    return 1.0 - 1.40 * std::exp(-0.66 * std::max(H, 0.5));
}

/// Variance-ratio difference V_r1 - V_r2 between shaft (moving) and bush roughness.
inline double roughness_variance_difference(const SimulationConfig& cfg) {
    const double sigma2 = cfg.roughness_rq_shaft * cfg.roughness_rq_shaft +
                          cfg.roughness_rq_bush * cfg.roughness_rq_bush;
    return (cfg.roughness_rq_shaft * cfg.roughness_rq_shaft -
            cfg.roughness_rq_bush * cfg.roughness_rq_bush) / sigma2;
}

/// Pressure flow factor at gap h; 1 when flow factors are off.
inline double pressure_flow_factor(double h, const SimulationConfig& cfg) {
    return flow_factors_active(cfg) ? pressure_flow_factor_isotropic(film_ratio(h, cfg)) : 1.0;
}

/// Extra Couette transport thickness sigma * phi_s [m]; 0 when flow factors are off.
inline double shear_flow_thickness(double h, const SimulationConfig& cfg) {
    if (!flow_factors_active(cfg)) return 0.0;
    return composite_roughness(cfg) * roughness_variance_difference(cfg) *
           shear_flow_factor_isotropic(film_ratio(h, cfg));
}

/// Multiplier of the shaft Couette shear mu U / h: phi_f + (V_r1 - V_r2) Phi_fs.
inline double couette_shear_factor(double h, const SimulationConfig& cfg) {
    if (!flow_factors_active(cfg)) return 1.0;
    const double H = film_ratio(h, cfg);
    return shear_stress_factor_f(H) + roughness_variance_difference(cfg) * shear_stress_factor_fs_isotropic(H);
}

/// Multiplier of the pressure-driven wall shear; 1 when flow factors are off.
inline double pressure_shear_factor(double h, const SimulationConfig& cfg) {
    return flow_factors_active(cfg) ? shear_stress_factor_fp(film_ratio(h, cfg)) : 1.0;
}

}  // namespace MixedLubrication
