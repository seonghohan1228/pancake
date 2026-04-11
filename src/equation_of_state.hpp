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

}  // namespace EOS
