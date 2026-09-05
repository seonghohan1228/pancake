#pragma once

#include <algorithm>
#include <cmath>

#include "config.hpp"
#include "cross_film.hpp"
#include "field.hpp"
#include "mesh.hpp"

/// Single definition of the coefficients that appear in every assembled form of
/// the Reynolds equation: the Poiseuille mobility, the film thickness that
/// carries the Couette flux, the Couette face mass flux, and the generalized-EOS
/// parameters (full-film ceiling, cavitation plateau, mixture bulk modulus).
///
/// Why one header rather than one copy per solver. Three call sites must agree
/// exactly or the solver lies about itself:
///
///  1. `Reynolds::solve_elrod` / `Reynolds::solve` assemble the operator.
///  2. `FBNS::solve` assembles a different operator over the same coefficients.
///  3. `Diagnostics` re-evaluates the assembled operator to report the per-cell
///     continuity residual (`mass_residual`) and the global mass balance.
///
/// (3) is a *measurement* of (1) and (2). When the three drifted apart, the
/// diagnostics kept using h^3/(12 mu) and `cfg.bulk_modulus` after the assembly
/// had moved to the Dowson F2 and the mixture beta_bar, so `mass_residual`
/// reported a fabricated imbalance on exactly the runs it exists to police.
/// Routing all three through this header makes that class of defect a compile
/// error rather than a silent wrong number.
///
/// Every function reduces identically to the classical, uniform-viscosity,
/// single-phase expression when the cross-film model is off and no gas phase is
/// resolved, and does so with the original floating-point association, so those
/// runs stay bit-identical.
namespace ReynoldsCoeff {

inline double field_value(const Fields& fields, const char* name, double fallback,
                          int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

inline double positive_field_value(const Fields& fields, const char* name, double fallback,
                                   int i, int j) {
    return std::max(field_value(fields, name, fallback, i, j), 1.0e-30);
}

/// Liquid-solution density [kg/m^3] used as the base density of every flux.
inline double base_density(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
}

/// Liquid-solution inventory per bearing area [kg/m^2].
inline double liquid_mass_area(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return base_density(fields, cfg, i, j) * fields["h"](i, j) * fields["theta"](i, j);
}

/// Previous inventory. Drivers should snapshot h, theta and rho_liquid_solution.
/// Legacy prescribed-squeeze callers without h history use h_old = h - dt dh_dt.
inline double old_liquid_mass_area(const Fields& fields, const SimulationConfig& cfg,
                                   double dt, int i, int j) {
    const Field& h = fields["h"];
    const Field& theta = fields["theta"];
    const double h_old = h.has_old_time() ? h.old(i, j)
        : h(i, j) - dt * field_value(fields, "dh_dt", 0.0, i, j);
    double rho_old = base_density(fields, cfg, i, j);
    if (fields.has("rho_liquid_solution") && fields["rho_liquid_solution"].has_old_time()) {
        rho_old = fields["rho_liquid_solution"].old(i, j);
    }
    const double theta_old = theta.has_old_time() ? theta.old(i, j) : theta(i, j);
    return rho_old * h_old * theta_old;
}

/// Cell dynamic viscosity [Pa s], already carrying the pressure, dissolved-gas
/// and mixture corrections of the selected property model.
inline double dynamic_viscosity(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return positive_field_value(fields, "mu", cfg.mu, i, j);
}

/// Poiseuille coefficient of the (generalized) Reynolds equation, per unit density.
///
/// Uniform across-the-gap viscosity gives h^3/(12 mu). With the cross-film model
/// active it is Dowson's F2 = int_0^h y (y - F1/F0) dy/mu, which reduces to
/// h^3/(12 mu) exactly when mu does not vary through the film. The two are
/// alternatives, not a correction of one by the other (PHYSICS.md 13.1.3).
///
/// `prefactor` multiplies the coefficient from the left; it is an argument rather
/// than applied by the caller so the classical branch keeps the exact original
/// floating-point association and stays bit-identical when the model is off.
inline double poiseuille(const Fields& fields, const SimulationConfig& cfg,
                         double prefactor, int i, int j) {
    if (CrossFilm::active(cfg) && fields.has("cross_film_f2")) {
        return prefactor * positive_field_value(fields, "cross_film_f2", 0.0, i, j);
    }
    const double h_val = fields["h"](i, j);
    return prefactor * h_val * h_val * h_val / (12.0 * dynamic_viscosity(fields, cfg, i, j));
}

/// The same coefficient with no density weighting, for velocity reconstruction.
inline double poiseuille_coefficient(const Fields& fields, const SimulationConfig& cfg,
                                     int i, int j) {
    return poiseuille(fields, cfg, 1.0, i, j);
}

/// True when the assembled fluxes must use the generalized cross-film
/// coefficients (Dowson F0/F1/F2) instead of the classical uniform-viscosity
/// h^3/(12 mu) and h/2.
inline bool generalized(const Fields& fields, const SimulationConfig& cfg) {
    return CrossFilm::active(cfg) && fields.has("cross_film_f2");
}

/// Face average of a per-cell cross-film coefficient. Passing the same cell
/// twice returns that cell's value exactly, which is what an axial boundary
/// face wants (there is no second cell to average with).
inline double cross_film_face_average(const Fields& fields, const char* name,
                                      int i0, int j0, int i1, int j1) {
    return 0.5 * (positive_field_value(fields, name, 0.0, i0, j0) +
                  positive_field_value(fields, name, 0.0, i1, j1));
}

/// Depth-averaged Couette velocity through a circumferential face [m/s].
///
/// Classical: U/2. Generalized: U (F1/F0) / h, i.e. the same 2 F1/F0 Couette
/// thickness the Reynolds operator transports, divided by the film thickness the
/// caller multiplies back in. Returns the classical expression verbatim when the
/// model is off, so those runs stay bit-identical.
inline double couette_face_velocity(const Fields& fields, const SimulationConfig& cfg,
                                    double R, double h_face,
                                    int i0, int j0, int i1, int j1) {
    if (generalized(fields, cfg)) {
        return 0.5 * cfg.omega * R *
               cross_film_face_average(fields, "cross_film_couette_h", i0, j0, i1, j1) / h_face;
    }
    return 0.5 * cfg.omega * R;
}

/// Depth-averaged Poiseuille velocity through a face [m/s].
///
/// Classical: -g h^2 (dp/ds) / (12 mu L). Generalized: -g F2 (dp/ds) / (h L), so
/// that h times this velocity is the generalized volumetric flux -g F2 grad p
/// that the Reynolds solve conserved. `length_scale` is the metric factor the
/// caller carries in the denominator (R for a theta face, 1 for a z face); it is
/// an argument rather than applied outside so the classical branch keeps its
/// original floating-point association.
///
/// `gate` is the caller's Poiseuille activity switch, which differs by equation
/// (the Elrod full-film gate, or 1 everywhere under FBNS), so it is passed in
/// rather than decided here.
inline double poiseuille_face_velocity(const Fields& fields, const SimulationConfig& cfg,
                                       double gate, double dp_ds, double h_face,
                                       double mu_face, double length_scale,
                                       int i0, int j0, int i1, int j1) {
    if (generalized(fields, cfg)) {
        return -gate * cross_film_face_average(fields, "cross_film_f2", i0, j0, i1, j1) *
               dp_ds / (h_face * length_scale);
    }
    return -gate * h_face * h_face * dp_ds / (12.0 * mu_face * length_scale);
}

/// Film thickness that carries the Couette flux. The generalized Reynolds
/// equation transports U * F1/F0 rather than U * h/2, so this returns 2 F1/F0,
/// which equals h under uniform viscosity and lets every Couette expression keep
/// its explicit factor of one half.
inline double couette_thickness(const Fields& fields, const SimulationConfig& cfg,
                                int i, int j) {
    if (CrossFilm::active(cfg) && fields.has("cross_film_couette_h")) {
        return positive_field_value(fields, "cross_film_couette_h", 0.0, i, j);
    }
    return fields["h"](i, j);
}

/// Couette mass flux [kg/s] through the east face of (i,j), carrying a film
/// content of 1. Both the Elrod and the FBNS assemblies build this face flux,
/// and both must see the same generalized Couette thickness.
inline double couette_face_flux(const Fields& fields, const Mesh& mesh,
                                const SimulationConfig& cfg, int i, int j) {
    const double h_e = 0.5 * (couette_thickness(fields, cfg, i, j) +
                              couette_thickness(fields, cfg, i + 1, j));
    const double rho_e = 0.5 * (base_density(fields, cfg, i, j) +
                                base_density(fields, cfg, i + 1, j));
    return rho_e * cfg.omega * mesh.R * h_e * 0.5 * mesh.get_d_z();
}

/// Liquid Poiseuille mobility of the FBNS operator [kg m / (Pa s)]: rho_l F2, or
/// rho_l h^3/(12 mu) classically. Unlike the Elrod path there is no switch
/// function and no bulk modulus -- cavitated cells sit at p = p_cav by
/// complementarity, so their pressure gradient, and hence their Poiseuille flux,
/// vanishes on its own.
inline double fbns_mobility(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return poiseuille(fields, cfg, base_density(fields, cfg, i, j), i, j);
}

// WP-1 two-phase void coupling: per-cell generalized EOS parameters. Under NONE
// (or when no free-gas field exists) these default to the single-phase values
// (theta_full = 1, p_void = effective p_cav, beta_bar = bulk_modulus), so the
// Elrod solve is bit-identical to the volumetrically-inert gas model.
struct VoidParams {
    double theta_full = 1.0;  // full-film ceiling = 1 - alpha_g
    double p_void = 0.0;      // local cavitated-zone plateau pressure
    double beta_bar = 0.0;    // mixture bulk modulus
};

/// Gas-inventory-consistent cavitated-zone plateau pressure [Pa].
///
/// Requiring the gas released between the bubble point p_sat and the local
/// pressure p to exactly fill the void (1 - theta) at the ideal-gas density
/// p / (R_s T) gives
///     H (p_sat - p) rho_l theta h . R_s T / p = (1 - theta) h
///  => p = L theta p_sat / [ (1 - theta) + L theta ],      L = H rho_l R_s T,
/// the Ostwald coefficient. This is the fixed point of the explicit form
/// p = m_g R_s T / ((1 - theta) h), but bounded: it tends to p_sat as theta -> 1
/// rather than diverging like m_g / (1 - theta), and never exceeds p_sat. The
/// explicit form divided a full-film cell's released gas by a ~zero void and
/// produced O(1e4) bar plateaus whenever the flash left free gas in an
/// uncavitated cell. Returns 0 when no Henry solubility model applies.
inline double gas_inventory_plateau(const Fields& fields, const SimulationConfig& cfg,
                                    double theta, int i, int j) {
    const double T = std::max(field_value(fields, "temperature", cfg.temperature_initial, i, j), 1.0);
    const double rho_l = positive_field_value(fields, "rho_liquid_solution", cfg.rho, i, j);
    const double L_ostwald = cfg.ostwald_coefficient(T, rho_l);
    const double p_sat = cfg.saturation_pressure(cfg.dissolved_gas_initial, T);
    if (!(L_ostwald > 0.0) || !(p_sat > 0.0)) return 0.0;
    const double theta_c = std::clamp(theta, 0.0, 1.0);
    const double denom = (1.0 - theta_c) + L_ostwald * theta_c;
    if (!(denom > 0.0)) return p_sat;
    return std::min(L_ostwald * theta_c * p_sat / denom, p_sat);
}

inline VoidParams void_params(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    VoidParams vp;
    vp.p_void = cfg.effective_p_cav();
    vp.beta_bar = cfg.bulk_modulus;
    if (cfg.gas_pressure_coupling != GasPressureCoupling::VOID_COUPLED || !fields.has("alpha_gas")) {
        return vp;  // single-phase defaults -> identical to the inert-gas model
    }
    const double alpha_cap = std::min(std::max(cfg.gas_alpha_max, 0.0), 0.999);
    const double alpha_g = std::clamp(fields["alpha_gas"](i, j), 0.0, alpha_cap);
    if (alpha_g <= 0.0) return vp;  // gas-free cell -> single-phase, bit-identical to NONE
    vp.theta_full = 1.0 - alpha_g;
    const double theta = std::max(field_value(fields, "theta", 1.0, i, j), cfg.theta_min);
    const double p = std::max(field_value(fields, "pressure", vp.p_void, i, j), cfg.gas_pressure_floor);
    // Mixture bulk modulus weighted by volume fraction: liquid (1 - alpha_g =
    // theta_full) at the liquid beta, gas alpha_g at the isothermal value 1/p.
    // Reduces to beta exactly when alpha_g = 0, so the no-gas limit matches the
    // single-phase model (the vaporous-limit acceptance gate).
    const double inv_beta = vp.theta_full / cfg.bulk_modulus + alpha_g / p;
    vp.beta_bar = inv_beta > 0.0 ? 1.0 / inv_beta : cfg.bulk_modulus;

    // When the lagged plateau field is present (VOID_COUPLED runs), read the
    // under-relaxed value the outer loop maintains. Recomputing the closure on the
    // instantaneous theta at every assembly makes the p_void <-> theta feedback
    // limit-cycle: the rupture front oscillates without ever flipping a flag, so
    // the Elrod loop stalls at a finite mass residual and never converges.
    if (fields.has("p_void")) {
        vp.p_void = std::max(vp.p_void, field_value(fields, "p_void", vp.p_void, i, j));
        return vp;
    }
    vp.p_void = std::max(vp.p_void, gas_inventory_plateau(fields, cfg, theta, i, j));
    return vp;
}

/// Elrod diffusion coefficient before the switch function:
/// Gamma_base = rho_0 beta_bar (F2 | h^3/12mu). The barotropic theta-formulation
/// carries the bulk modulus here, which is exactly what makes it different from
/// the FBNS mobility above; both are needed and neither may stand in for the other.
inline double elrod_gamma_base(const Fields& fields, const SimulationConfig& cfg, int i, int j) {
    return poiseuille(fields, cfg,
                      base_density(fields, cfg, i, j) * void_params(fields, cfg, i, j).beta_bar,
                      i, j);
}

}  // namespace ReynoldsCoeff
