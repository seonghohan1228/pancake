#pragma once

#include <array>
#include <cmath>

#include "config.hpp"

/// Asperity contact for mixed lubrication (Greenwood-Tripp).
///
/// Nominal asperity contact pressure on the local film thickness h:
///
///     p_c(h) = K E* F_{5/2}(h / sigma),
///     K      = (16 sqrt(2) pi / 15) (eta beta sigma)^2 sqrt(sigma / beta),
///     1/E*   = (1 - nu_1^2)/E_1 + (1 - nu_2^2)/E_2,
///
/// with the exact Gaussian integral
///
///     F_n(lambda) = (2 pi)^{-1/2} int_lambda^inf (s - lambda)^n exp(-s^2/2) ds.
///
/// K is entered directly (`contact_gt_k`) because measured journal-bearing data
/// report the product, not eta and beta separately. The 16 sqrt(2) pi/15 prefactor
/// pairs only with the Hertz modulus E* above; the 8 sqrt(2) pi/15 form pairs
/// with E' = 2 E*. Boundary friction on the asperity load is mu_b p_c. See
/// PHYSICS.md section "Asperity contact" for sources and validity limits.
namespace Contact {

inline bool active(const SimulationConfig& cfg) {
    return cfg.contact_model == ContactModel::GREENWOOD_TRIPP;
}

/// Exact F_{5/2}(lambda) by composite Simpson quadrature (used to build the table).
inline double f52_exact(double lambda) {
    const double upper = std::max(lambda, 0.0) + 12.0;
    const double lower = lambda;
    const int n = 4000;  // even
    const double ds = (upper - lower) / n;
    auto g = [lambda](double s) {
        const double x = s - lambda;
        return x <= 0.0 ? 0.0 : std::pow(x, 2.5) * std::exp(-0.5 * s * s);
    };
    double sum = g(lower) + g(upper);
    for (int k = 1; k < n; ++k) sum += (k % 2 ? 4.0 : 2.0) * g(lower + k * ds);
    return sum * ds / 3.0 / std::sqrt(2.0 * M_PI);
}

/// Tabulated F_{5/2} on lambda in [0, 8]; log-linear interpolation; zero beyond 8
/// (F_{5/2}(8) ~ 1e-17). lambda < 0 is clamped to 0: the film-thickness guard
/// rejects nonpositive gaps before this is reached.
inline double f52(double lambda) {
    constexpr int N = 1601;
    constexpr double LMAX = 8.0;
    static const std::array<double, N> table = [] {
        std::array<double, N> t{};
        for (int k = 0; k < N; ++k) t[k] = std::log(f52_exact(LMAX * k / (N - 1)));
        return t;
    }();
    if (!(lambda < LMAX)) return 0.0;
    const double x = std::max(lambda, 0.0) / LMAX * (N - 1);
    const int k = std::min(static_cast<int>(x), N - 2);
    const double w = x - k;
    return std::exp((1.0 - w) * table[k] + w * table[k + 1]);
}

/// Asperity contact pressure [Pa] at local film thickness h [m].
inline double pressure(double h, const SimulationConfig& cfg) {
    if (!active(cfg)) return 0.0;
    return cfg.contact_gt_k * cfg.contact_modulus * f52(h / cfg.contact_roughness_sigma);
}

}  // namespace Contact
