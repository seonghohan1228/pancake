#pragma once

#include <cmath>

/// Barotropic equation of state for the compressible thin-film model.
/// Relates the universal variable theta (fractional film content) to pressure and density.
///
/// Full-film region (theta >= 1):  p = p_cav + beta * ln(theta)
/// Cavitated region (theta < 1):   p = p_cav, theta is the fraction of gap filled with liquid
///
/// Density: rho = rho_0 * theta (both regions)
namespace EOS {

    /// Generalized barotropic EOS for the two-phase (WP-1) void coupling. The
    /// full-film ceiling is theta_full = 1 - alpha_g (gas occupies the rest of the
    /// gap); the cavitated-zone plateau is the local p_void; beta_bar is the
    /// mixture bulk modulus. Pressurized (theta >= theta_full):
    ///   p = p_void + beta_bar * ln(theta / theta_full); cavitated: p = p_void.
    inline double pressure_from_theta(double theta, double theta_full, double p_void, double beta_bar) {
        const double tf = (theta_full > 0.0) ? theta_full : 1.0;
        return (theta >= tf) ? p_void + beta_bar * std::log(theta / tf) : p_void;
    }

    inline double theta_from_pressure(double p, double theta_full, double p_void, double beta_bar) {
        const double tf = (theta_full > 0.0) ? theta_full : 1.0;
        return (p > p_void) ? tf * std::exp((p - p_void) / beta_bar) : tf;
    }

    /// Single-phase back-compat overloads (theta_full = 1, p_void = p_cav, beta_bar = beta).
    inline double pressure_from_theta(double theta, double p_cav, double beta) {
        return pressure_from_theta(theta, 1.0, p_cav, beta);
    }

    inline double theta_from_pressure(double p, double p_cav, double beta) {
        return theta_from_pressure(p, 1.0, p_cav, beta);
    }

    inline double density(double theta, double rho_0) {
        return rho_0 * theta;
    }

    /// dp/dtheta used for the effective diffusion coefficient in the theta-formulation.
    inline double dp_dtheta(double theta, double beta) {
        return (theta >= 1.0) ? beta / theta : 0.0;
    }

    /// Elrod-Adams switch function: 1 in full-film (theta >= theta_full), 0 in
    /// cavitated. Multiplied into the diffusion coefficient Γ to disable
    /// Poiseuille flow in cavitation. theta_full = 1 recovers the single-phase case.
    inline double switch_function(double theta, double theta_full) {
        return (theta >= theta_full) ? 1.0 : 0.0;
    }

    inline double switch_function(double theta) {
        return switch_function(theta, 1.0);
    }

    /// Smooth approximation of the switch function for initialisation or diagnostics.
    /// eps controls the transition width; approaches the hard switch as eps -> 0.
    inline double switch_function_smooth(double theta, double eps) {
        return 0.5 * (1.0 + std::tanh((theta - 1.0) / eps));
    }

}  // namespace EOS
