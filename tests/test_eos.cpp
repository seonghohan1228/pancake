// Phase 2 validation: barotropic EOS round-trip and property checks.
// Run as: ./test_eos
#include <cassert>
#include <cmath>
#include <iostream>

#include "equation_of_state.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        std::exit(1);
    }
}

int main() {
    const double p_cav = 0.0;
    const double beta  = 1e5;
    const double rho_0 = 900.0;
    const double tol   = 1e-10;

    // Round-trip only valid for theta >= 1 (full-film): theta -> p -> theta
    // Cavitated (theta < 1) maps to p = p_cav, which maps back to theta = 1 (boundary).
    // This is intentional: pressure alone cannot recover film content in cavitated region.
    for (double theta : {1.0, 1.5, 2.0}) {
        double p      = EOS::pressure_from_theta(theta, p_cav, beta);
        double theta2 = EOS::theta_from_pressure(p, p_cav, beta);
        check(std::abs(theta2 - theta) < tol, "round-trip theta->p->theta (full-film)");
    }

    // Cavitated region: p = p_cav for theta < 1
    check(EOS::pressure_from_theta(0.5, p_cav, beta) == p_cav,
          "cavitated region pressure == p_cav");
    check(EOS::pressure_from_theta(0.9, p_cav, beta) == p_cav,
          "cavitated region pressure == p_cav (0.9)");

    // Full-film: p > p_cav for theta > 1
    check(EOS::pressure_from_theta(1.5, p_cav, beta) > p_cav,
          "full-film pressure > p_cav");

    // Density: rho = rho_0 * theta
    for (double theta : {0.5, 1.0, 1.5}) {
        check(std::abs(EOS::density(theta, rho_0) - rho_0 * theta) < tol,
              "density = rho_0 * theta");
    }

    // dp/dtheta = beta/theta for theta >= 1, else 0
    check(std::abs(EOS::dp_dtheta(1.5, beta) - beta / 1.5) < tol,
          "dp_dtheta full-film");
    check(EOS::dp_dtheta(0.5, beta) == 0.0,
          "dp_dtheta cavitated == 0");

    // theta_from_pressure: p <= p_cav returns theta = 1 (not < 1)
    check(EOS::theta_from_pressure(p_cav - 1.0, p_cav, beta) == 1.0,
          "sub-cavitation pressure returns theta=1");

    // --- WP-1 generalized two-phase EOS (theta_full, p_void, beta_bar) ------
    {
        const double pv = 1.0e5;   // local void plateau
        const double bb = 5.0e4;   // mixture bulk modulus
        const double tf = 0.6;     // full-film ceiling (alpha_g = 0.4)

        // Back-compat: theta_full = 1 reproduces the single-phase overloads exactly.
        for (double theta : {0.5, 1.0, 1.5, 2.0}) {
            check(EOS::pressure_from_theta(theta, 1.0, pv, bb) ==
                  EOS::pressure_from_theta(theta, pv, bb), "4-arg theta_full=1 == 3-arg pressure");
            check(EOS::switch_function(theta, 1.0) == EOS::switch_function(theta),
                  "switch theta_full=1 == 1-arg");
        }

        // The switch and onset move to theta_full, not 1.
        check(EOS::switch_function(0.7, tf) == 1.0, "full-film above theta_full");
        check(EOS::switch_function(0.5, tf) == 0.0, "cavitated below theta_full");
        check(std::abs(EOS::pressure_from_theta(tf, tf, pv, bb) - pv) < tol,
              "pressure at theta_full == p_void");
        check(EOS::pressure_from_theta(0.5, tf, pv, bb) == pv, "cavitated plateau == p_void");

        // Pressurized round-trip for theta >= theta_full.
        for (double theta : {0.6, 0.9, 1.2}) {
            const double p = EOS::pressure_from_theta(theta, tf, pv, bb);
            const double t2 = EOS::theta_from_pressure(p, tf, pv, bb);
            check(std::abs(t2 - theta) < 1e-9, "generalized round-trip theta->p->theta");
            check(p >= pv, "pressurized p >= p_void");
        }
    }

    std::cout << "PASS: test_eos\n";
    return 0;
}
