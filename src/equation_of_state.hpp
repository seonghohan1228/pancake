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

    inline double pressure_from_theta(double theta, double p_cav, double beta) {
        return (theta >= 1.0) ? p_cav + beta * std::log(theta) : p_cav;
    }

    inline double theta_from_pressure(double p, double p_cav, double beta) {
        return (p > p_cav) ? std::exp((p - p_cav) / beta) : 1.0;
    }

    inline double density(double theta, double rho_0) {
        return rho_0 * theta;
    }

    /// dp/dtheta used for the effective diffusion coefficient in the theta-formulation.
    inline double dp_dtheta(double theta, double beta) {
        return (theta >= 1.0) ? beta / theta : 0.0;
    }

    /// Elrod-Adams switch function: 1 in full-film (theta >= 1), 0 in cavitated (theta < 1).
    /// Multiplied into the diffusion coefficient Γ to disable Poiseuille flow in cavitation.
    inline double switch_function(double theta) {
        return (theta >= 1.0) ? 1.0 : 0.0;
    }

    /// Smooth approximation of the switch function for initialisation or diagnostics.
    /// eps controls the transition width; approaches the hard switch as eps -> 0.
    inline double switch_function_smooth(double theta, double eps) {
        return 0.5 * (1.0 + std::tanh((theta - 1.0) / eps));
    }

}  // namespace EOS
