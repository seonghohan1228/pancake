#include "cross_film.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <mpi.h>
#include <stdexcept>
#include <vector>

#include "bush_conduction.hpp"
#include "fluid_properties.hpp"
#include "output_naming.hpp"
#include "resolved_cross_film.hpp"

namespace CrossFilm {
namespace {

/// Zero-mean shape functions on eta = y/h in [0, 1], y measured from the shaft
/// surface (eta = 0) to the bush surface (eta = 1).
///
///   psi1 = eta - 1/2                      antisymmetric, carries the wall difference
///   psi2 = eta - eta^2 - 1/6              parabolic core, the conduction shape
///   psi4 = eta^2 (1-eta)^2 - 1/30         quartic core flattening, the convection shape
///
/// All three integrate to zero over [0, 1], so the depth average of
/// T = T_bar + A psi1 + B psi2 + C psi4 is T_bar by construction: the profile is
/// collocated with, and never contradicts, the depth-averaged energy balance.
double psi1(double e) { return e - 0.5; }
double psi2(double e) { return e - e * e - 1.0 / 6.0; }
double psi4(double e) { const double q = e * (1.0 - e); return q * q - 1.0 / 30.0; }

void require_globally_valid(bool invalid, const char* message) {
    const int local = invalid ? 1 : 0;
    int global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global) throw std::runtime_error(message);
}

double profile_temperature(double t_bar, double a, double b, double c, double eta) {
    return t_bar + a * psi1(eta) + b * psi2(eta) + c * psi4(eta);
}

std::array<double, 2> profile_extrema(double t_bar, double a, double b, double c) {
    const auto derivative = [=](double eta) {
        const double x = eta - 0.5;
        return a - (2.0 * b + c) * x + 4.0 * c * x * x * x;
    };
    std::vector<double> limits{0.0, 1.0};
    if (c != 0.0) {
        const double root_squared = (c + 2.0 * b) / (12.0 * c);
        if (root_squared > 0.0 && root_squared < 0.25) {
            const double root = std::sqrt(root_squared);
            limits.push_back(0.5 - root);
            limits.push_back(0.5 + root);
        }
    }
    std::sort(limits.begin(), limits.end());
    std::array<double, 2> extrema{
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
    const auto include = [&](double eta) {
        const double value = profile_temperature(t_bar, a, b, c, eta);
        if (!std::isfinite(value)) {
            extrema.fill(std::numeric_limits<double>::quiet_NaN());
            return;
        }
        extrema[0] = std::min(extrema[0], value);
        extrema[1] = std::max(extrema[1], value);
    };
    for (double eta : limits) include(eta);
    for (size_t k = 1; k < limits.size(); ++k) {
        double left = limits[k - 1], right = limits[k];
        const double d_left = derivative(left), d_right = derivative(right);
        if (!std::isfinite(d_left) || !std::isfinite(d_right) ||
            !((d_left < 0.0 && d_right > 0.0) || (d_left > 0.0 && d_right < 0.0))) continue;
        for (int iteration = 0; iteration < 40; ++iteration) {
            const double mid = 0.5 * (left + right);
            if ((derivative(mid) < 0.0) == (d_left < 0.0)) left = mid;
            else right = mid;
        }
        include(0.5 * (left + right));
    }
    return extrema;
}

/// Gauss-Legendre nodes and weights on [0, 1], by Newton iteration on P_n.
void gauss_legendre_unit(int n, std::vector<double>& nodes, std::vector<double>& weights) {
    nodes.assign(n, 0.0);
    weights.assign(n, 0.0);
    for (int k = 0; k < n; ++k) {
        double x = std::cos(M_PI * (k + 0.75) / (n + 0.5));  // Chebyshev starting guess
        double dp = 1.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = 0.0;
            for (int m = 0; m < n; ++m) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * m + 1.0) * x * p1 - m * p2) / (m + 1.0);
            }
            dp = n * (x * p0 - p1) / (x * x - 1.0);
            const double dx = -p0 / dp;
            x += dx;
            if (std::abs(dx) < 1.0e-15) break;
        }
        nodes[k] = 0.5 * (1.0 + x);
        weights[k] = 1.0 / ((1.0 - x * x) * dp * dp);  // 2/((1-x^2)p'^2) halved for [0,1]
    }
}

double field_or(const Fields& fields, const char* name, double fallback, int i, int j) {
    return fields.has(name) ? fields[name](i, j) : fallback;
}

/// Bush wall thermal resistance referred to the film-side surface area
/// [m^2 K / W]: radial conduction through the bush in series with outer-surface
/// convection. Returns the conductance 1/R; zero means an adiabatic bush.
double bush_conductance(const SimulationConfig& cfg) {
    const double r_i = cfg.R;
    const double r_o = cfg.bush_outer_radius > r_i ? cfg.bush_outer_radius : 1.1 * r_i;
    const double wall = r_i * std::log(r_o / r_i) / cfg.bush_thermal_conductivity;
    // FIXED_TEMPERATURE is the h_out -> infinity limit of the same series
    // resistance, so the two outer conditions share one expression.
    if (cfg.bush_outer_boundary == BushOuterBoundary::FIXED_TEMPERATURE) {
        return wall > 0.0 ? 1.0 / wall : 0.0;
    }
    if (cfg.bush_outer_heat_transfer <= 0.0) return 0.0;
    const double resistance = wall + r_i / (r_o * cfg.bush_outer_heat_transfer);
    return resistance > 0.0 ? 1.0 / resistance : 0.0;
}

bool bush_is_2d(const SimulationConfig& cfg) {
    return cfg.bush_thermal_model == BushThermalModel::CONDUCTION_2D;
}

}  // namespace

void register_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    // Profile coefficients of T(eta) = T_bar + a psi1 + b psi2 + c psi4 [K].
    fields.add("cross_film_a", mesh).fill(0.0);
    fields.add("cross_film_b", mesh).fill(0.0);
    fields.add("cross_film_c", mesh).fill(0.0);
    // Relaxation carrier for the midplane collocation source X = (h^2/k) rhs [K].
    fields.add("cross_film_collocation", mesh).fill(0.0);
    // Dowson viscosity integrals. f2 replaces h^3/(12 mu) as the Poiseuille
    // coefficient [m^3/(Pa s)]; couette_h is 2*F1/F0, the effective film thickness
    // of the Couette flux [m], equal to h whenever mu is uniform across the gap.
    fields.add("cross_film_f0", mesh).fill(0.0);
    fields.add("cross_film_f1", mesh).fill(0.0);
    fields.add("cross_film_f2", mesh).fill(0.0);
    fields.add("cross_film_couette_h", mesh).fill(0.0);
    // Wall closure: net film->wall heat flux [W/m^2] and its derivative with
    // respect to the depth-averaged temperature [W/(m^2 K)].
    fields.add("cross_film_wall_flux", mesh).fill(0.0);
    fields.add("cross_film_wall_coeff", mesh).fill(0.0);
    // Conjugate wall temperatures [K], for output and inspection.
    fields.add("cross_film_t_shaft", mesh).fill(0.0);
    fields.add("cross_film_t_bush", mesh).fill(0.0);
    // Diagnostic: midplane temperature minus the depth-averaged temperature [K].
    fields.add("cross_film_core_rise", mesh).fill(0.0);

    if (cfg.cross_film_thermal_model == CrossFilmModel::FINITE_VOLUME)
        ResolvedCrossFilm::register_fields(fields, mesh, cfg);

    // Optional film-surface outputs. Allocated only when the user selected them
    // AND the model that defines them is on, so an unselected run allocates
    // nothing, `update()` skips the stores, and the written files are unchanged.
    // The config validator rejects a selection made under NONE, so an absent
    // field here is always a deliberate non-selection and never a silent zero.
    if (active(cfg)) {
        for (const char* name : OutputNaming::surface_fields()) {
            if (cfg.output_field_enabled(name)) fields.add(name, mesh).fill(0.0);
        }
    }
}

WallState update(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg, double relaxation) {
    if (cfg.cross_film_thermal_model == CrossFilmModel::FINITE_VOLUME)
        return ResolvedCrossFilm::update(fields, mesh, cfg, relaxation);
    WallState state;
    if (!active(cfg)) return state;

    const Field& h = fields["h"];
    const Field& p = fields["pressure"];
    const Field& temperature = fields["temperature"];
    Field& a_f = fields["cross_film_a"];
    Field& b_f = fields["cross_film_b"];
    Field& c_f = fields["cross_film_c"];
    Field& x_prev_f = fields["cross_film_collocation"];
    Field& f0_f = fields["cross_film_f0"];
    Field& f1_f = fields["cross_film_f1"];
    Field& f2_f = fields["cross_film_f2"];
    Field& couette_h_f = fields["cross_film_couette_h"];
    Field& wall_flux_f = fields["cross_film_wall_flux"];
    Field& wall_coeff_f = fields["cross_film_wall_coeff"];
    Field& t_shaft_f = fields["cross_film_t_shaft"];
    Field& t_bush_f = fields["cross_film_t_bush"];
    Field& core_rise_f = fields["cross_film_core_rise"];
    // Optional film-surface outputs; null unless selected (see register_fields).
    Field* t_journal_f = fields.has("temperature_journal") ? &fields["temperature_journal"] : nullptr;
    Field* t_bearing_f = fields.has("temperature_bearing") ? &fields["temperature_bearing"] : nullptr;
    Field* q_journal_f = fields.has("heat_flux_journal") ? &fields["heat_flux_journal"] : nullptr;
    Field* q_bearing_f = fields.has("heat_flux_bearing") ? &fields["heat_flux_bearing"] : nullptr;

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;
    const double U = cfg.omega * R;
    const double area = R * d_theta * d_z;
    const bool quartic = cfg.cross_film_profile_order >= 4;

    bool invalid_input = false;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = h(i, j);
            const double t_val = temperature(i, j);
            const double mu_val = field_or(fields, "mu", cfg.mu, i, j);
            const double k_val = field_or(fields, "k_liquid_solution", cfg.thermal_conductivity, i, j);
            const double rho_cp = field_or(fields, "rho_cp_liquid_solution", cfg.rho_cp, i, j);
            invalid_input |= !(std::isfinite(h_val) && h_val > 0.0 &&
                               std::isfinite(t_val) && t_val > 0.0 &&
                               std::isfinite(p(i, j)) &&
                               std::isfinite(mu_val) && mu_val > 0.0 &&
                               std::isfinite(k_val) && k_val > 0.0 &&
                               std::isfinite(rho_cp) && rho_cp > 0.0);
        }
    }
    require_globally_valid(invalid_input, "Cross-film input contains a nonphysical or nonfinite value");

    std::vector<double> nodes, weights;
    gauss_legendre_unit(std::max(2, cfg.cross_film_quadrature_points), nodes, weights);

    const double bush_g = bush_conductance(cfg);
    const double t_bush_outer = cfg.bush_outer_temperature > 0.0
        ? cfg.bush_outer_temperature : cfg.thermal_ambient_temperature;
    const double t_shaft_ambient = cfg.shaft_ambient_temperature > 0.0
        ? cfg.shaft_ambient_temperature : cfg.thermal_ambient_temperature;

    // ---------------------------------------------------------------------
    // Step 1: the quartic amplitude c, from one collocation of the cross-film
    // energy balance at the midplane.
    //
    //   k d2T/dy2 = rho cp (u dT/ds + w dT/dz) - Phi(y)
    //
    // Only the shape functions psi2 and psi4 have curvature; psi2'' = -2 and
    // psi4''(1/2) = -1 in eta, so with T'' = (1/h^2) d2T/deta2,
    //
    //   (k/h^2) (-2 b - c) = rhs_mid.
    //
    // Combined with the zero-mean basis and the two wall values (which give
    // b + c/5 = 6 (T_bar - T_wall_mean)) this closes the profile. `c` does not
    // depend on the wall temperatures, so it is solved first and then used by
    // the conjugate wall balance below.
    //
    // At order 2 the collocation is dropped and c = 0, leaving the classical
    // parabola whose amplitude follows from the mean and the two wall values.
    // ---------------------------------------------------------------------
    std::vector<double> x_collocation(static_cast<size_t>(n_theta) * n_z, 0.0);
    double profile_residual = 0.0;
    bool invalid_collocation = false;
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            if (!quartic) continue;
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double k_th = std::max(field_or(fields, "k_liquid_solution",
                                                  cfg.thermal_conductivity, i, j), 1.0e-30);
            const double rho_cp = std::max(field_or(fields, "rho_cp_liquid_solution",
                                                    cfg.rho_cp, i, j), 1.0e-30);
            // Lagged F-integrals: the profile they came from is one Picard
            // iteration old, which is the same lag the segregated loop already
            // carries on viscosity.
            double F0 = f0_f(i, j);
            if (!(F0 > 0.0)) {  // first call: uniform mu across the gap
                F0 = h_val / std::max(field_or(fields, "mu", cfg.mu, i, j), 1.0e-30);
            }
            const double mu_mid = std::max(field_or(fields, "mu", cfg.mu, i, j), 1.0e-30);

            // Couette-only collocation. The full stress is mu du/dy = (dp/ds) y + c1;
            // this drops the (dp/ds) part and keeps c1 = -U/F0, the shear that drags
            // the film with the shaft.
            //
            // The omission is deliberate and was measured, not assumed. Including
            // dp/ds makes the profile shape a function of the local pressure
            // gradient, which jumps discontinuously whenever an Elrod cell flips
            // between full film and cavitated. That closes a feedback loop
            // profile -> viscosity -> Couette flux -> rupture position -> profile,
            // and the segregated solve limit-cycles: on the Ferron rig the coupling
            // residual R_p sat at 9.4 bar, identical at 40 and at 120 outer
            // iterations, and under-relaxing the profile amplitude did not move it
            // (9.4 / 2.3 / 7.6 / 2.2 bar at w = 1.0 / 0.5 / 0.3 / 0.15). With the
            // pressure terms removed the same case converges to R_p = 3.6e-6 and
            // the load changes by 0.7%.
            //
            // What is given up is small and bounded: at Ferron's operating point the
            // Poiseuille share of the dissipation is about 6% of the Couette share,
            // and the Poiseuille share of the midplane velocity about 14%, and both
            // feed a shape correction that the saturation below already caps.
            // Everything else -- the F-integrals, the dissipation source
            // U^2/F0 + F2 |grad p|^2, the wall fluxes and the generalized Reynolds
            // coefficients -- still uses the full pressure field.
            const double c1_s = -U / F0;
            const double y_mid = 0.5 * h_val;
            const double tau_s = c1_s;
            const double phi_mid = tau_s * tau_s / mu_mid;
            const double u_mid = U + c1_s * y_mid / mu_mid;
            const double w_mid = 0.0;

            const double dT_ds = (temperature(i + 1, j) - temperature(i - 1, j)) / (2.0 * d_theta * R);
            const double dT_dz = (n_z > 1)
                ? (temperature(i, std::min(j + 1, n_z - 1)) - temperature(i, std::max(j - 1, 0))) /
                  ((std::min(j + 1, n_z - 1) - std::max(j - 1, 0)) * d_z)
                : 0.0;

            const double scale = h_val * h_val / k_th;
            const double x_dissipation = -scale * phi_mid;
            const double x_convection = scale * rho_cp * (u_mid * dT_ds + w_mid * dT_dz);
            // Saturate the convective term against the dissipative one.
            //
            // A single midplane collocation cannot represent the thin wall layers
            // that appear once the reduced Peclet rho c U h^2/(k R) leaves order
            // one, and on these rigs it is 10-30. Left unbounded the balance is
            // dominated by rho c u dT/ds -- measured at up to 5.6e3 K of equivalent
            // (h^2/k) rhs near a groove against 23 K of dissipation -- and produces
            // a physically impossible several-hundred-kelvin core. The saturation
            //     x_conv_eff = x_conv |x_phi| / (|x_phi| + |x_conv|)
            // is exact (x_conv_eff -> x_conv) in the conduction-dominated limit the
            // profile is derived for, so the classical parabola is recovered
            // exactly there, and it caps the convective flattening at parity with
            // dissipation instead of letting it run. Where it binds, the quartic
            // degenerates towards a fixed flattened shape rather than resolving the
            // wall layers: that is a stated limitation of the one-point collocation,
            // not a converged boundary-layer solution.
            const double x_conv_eff = (std::abs(x_dissipation) > 0.0)
                ? x_convection * std::abs(x_dissipation) /
                      (std::abs(x_dissipation) + std::abs(x_convection))
                : 0.0;
            // -2b - c = (h^2/k) rhs_mid  with  b = 6 dT_wall - c/5  gives
            //   c = -(5/3) [ (h^2/k) rhs_mid + 12 (T_bar - T_wall_mean) ].
            // T_wall_mean is not known yet, so the identity is applied after the
            // wall solve; store the part that does not depend on it.
            // Under-relax the collocation source itself, not the coefficient it
            // produces: relaxing `c` downstream would leave it inconsistent with
            // the wall balance solved from the same X, and the adiabatic-shaft
            // identity would stop holding exactly. Relaxing X keeps every algebraic
            // relation below exact at every call and only changes the path to the
            // fixed point.
            const double w = std::clamp(cfg.cross_film_relaxation, 0.01, 1.0);
            const double x_raw = x_dissipation + x_conv_eff;
            const double x_prev = x_prev_f(i, j);
            if (!std::isfinite(x_raw) || !std::isfinite(x_prev)) {
                invalid_collocation = true;
                continue;
            }
            profile_residual = std::max(profile_residual,
                std::abs(x_raw - x_prev) /
                std::max(std::abs(cfg.temperature_reference), 1.0));
            const double x_new = w * x_raw + (1.0 - w) * x_prev;
            x_prev_f(i, j) = x_new;
            x_collocation[static_cast<size_t>(i) * n_z + j] = x_new;
        }
    }
    require_globally_valid(invalid_collocation, "Cross-film collocation source is nonfinite");
    MPI_Allreduce(MPI_IN_PLACE, &profile_residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    state.profile_residual = profile_residual;

    // ---------------------------------------------------------------------
    // Step 2: conjugate walls.
    //
    // Shaft: circumferentially and axially isothermal at T_s, its temperature set
    // by a global film-to-shaft balance. This is a reduced closure: Ferron solves
    // an r-z shaft problem with prescribed 43 C axial-end temperatures. With
    // shaft_ambient_heat_transfer = 0 the shaft is globally adiabatic: the
    // integral of the film-to-shaft flux is zero, so heat entering the shaft in
    // the loaded arc leaves it again upstream.
    //
    // Bush (LOCAL, the default): per-cell 1-D radial conduction to the outer
    // boundary, q_bush = (T_b - T_bush_outer) * G, G referred to the bore. With
    // bush_thermal_model = CONDUCTION_2D the same T_b = alpha + beta T_s comes
    // instead from the annulus solve of step 2a; everything downstream is
    // unchanged, which is why only the alpha/beta assignment below is branched.
    //
    // With the profile fluxes q_shaft = (k/h)(a + b) and q_bush = (k/h)(b - a),
    // a = T_b - T_s and b = 6(T_bar - (T_s+T_b)/2) - c/5, the bush relation is
    // linear in T_b for a given T_s:
    //
    //   T_b = alpha + beta T_s,
    //   alpha = [kh (6 T_bar - c/5) + G T_amb] / (4 kh + G),  beta = -2 kh / (4 kh + G),
    //
    // and substituting into the shaft balance closes T_s in one step. Only four
    // scalars need reducing across ranks.
    // ---------------------------------------------------------------------
    const double shaft_sink = cfg.shaft_ambient_heat_transfer * cfg.shaft_ambient_area_ratio *
                              (2.0 * M_PI * R * mesh.L);
    double sums[4] = {0.0, 0.0, 0.0, 0.0};  // K, P, Q_alpha, Q_beta
    std::vector<double> alpha(static_cast<size_t>(n_theta) * n_z, 0.0);
    std::vector<double> beta(static_cast<size_t>(n_theta) * n_z, 0.0);

    // ---------------------------------------------------------------------
    // Step 2a (CONDUCTION_2D only): the bush annulus.
    //
    // The film-side flux is linear in the bore temperature,
    //   q_bush = A0 + A1 T_s - B T_b,
    //   A0 = kh (b_const + b_slope T_bar),  A1 = kh (1 - b_slope/2),
    //   B  = kh (1 + b_slope/2),
    // which is exactly the relation the LOCAL closure inverts against
    // q_bush = G (T_b - T_amb). Here it becomes a Robin condition on the bore of
    // a steady 2-D (r, theta) conduction problem instead (PHYSICS.md 13.1.2).
    //
    // Two modelling choices, both stated rather than hidden:
    //  * the bush is z-collapsed, so the coupling uses the z-AVERAGE of the
    //    film-side flux law at each theta and returns one bore temperature per
    //    theta. Ferron solves the r-theta-z bush equation instead;
    //  * Ferron's lateral (axial end-face) free convection is absent. A
    //    z-collapsed sink could approximate its integrated loss only, not the
    //    local r-theta-z solution.
    //
    // Both right-hand sides are solved, so T_b(theta) = alpha + beta T_s stays
    // exact and the shaft balance below closes in one step exactly as before.
    // The solve is serial and deterministic: the three coupling arrays are
    // all-reduced onto the global theta grid and every rank then runs the same
    // arithmetic on the same input, which removes the broadcast and any rank-0
    // asymmetry for a few hundred thousand flops per outer iteration.
    // ---------------------------------------------------------------------
    if (bush_is_2d(cfg)) {
        const int n_theta_global = mesh.n_theta_global;
        BushConduction::Coupling coupling;
        coupling.a0.assign(n_theta_global, 0.0);
        coupling.a1.assign(n_theta_global, 0.0);
        coupling.slope.assign(n_theta_global, 0.0);
        const double inv_nz = 1.0 / static_cast<double>(n_z);
        for (int i = 0; i < n_theta; ++i) {
            const int gi = mesh.offset_theta + i;
            double a0 = 0.0, a1 = 0.0, slope = 0.0;
            for (int j = 0; j < n_z; ++j) {
                const double h_val = std::max(h(i, j), cfg.min_film_thickness);
                const double k_th = std::max(field_or(fields, "k_liquid_solution",
                                                      cfg.thermal_conductivity, i, j), 1.0e-30);
                const double kh = k_th / h_val;
                const double X = x_collocation[static_cast<size_t>(i) * n_z + j];
                const double b_const = quartic ? (X / 3.0) : 0.0;
                const double b_slope = quartic ? 10.0 : 6.0;
                a0 += kh * (b_const + b_slope * temperature(i, j));
                a1 += kh * (1.0 - 0.5 * b_slope);
                slope += kh * (1.0 + 0.5 * b_slope);
            }
            coupling.a0[gi] = a0 * inv_nz;
            coupling.a1[gi] = a1 * inv_nz;
            coupling.slope[gi] = slope * inv_nz;
        }
        MPI_Allreduce(MPI_IN_PLACE, coupling.a0.data(), n_theta_global, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, coupling.a1.data(), n_theta_global, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, coupling.slope.data(), n_theta_global, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        const BushConduction::Solution bush =
            BushConduction::solve(coupling, n_theta_global, R, cfg);
        require_globally_valid(!bush.converged || !std::isfinite(bush.residual),
                               "Cross-film 2-D bush conduction did not converge");
        state.bush_2d_solved = true;
        state.bush_residual = bush.residual;
        for (int i = 0; i < n_theta; ++i) {
            const int gi = mesh.offset_theta + i;
            for (int j = 0; j < n_z; ++j) {
                const size_t idx = static_cast<size_t>(i) * n_z + j;
                alpha[idx] = bush.alpha[gi];
                beta[idx] = bush.beta[gi];
            }
        }
    }

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double k_th = std::max(field_or(fields, "k_liquid_solution",
                                                  cfg.thermal_conductivity, i, j), 1.0e-30);
            const double kh = k_th / h_val;
            const double t_bar = temperature(i, j);
            // X is the collocation source (h^2/k) rhs(1/2) from step 1. Eliminating
            // c between the collocation and the wall conditions,
            //   c = -(5/3) X - 20 (T_bar - T_m),
            //   b = 6 (T_bar - T_m) - c/5 = (1/3) X + 10 (T_bar - T_m),
            // leaves both coefficients affine in T_s and T_b, so the wall balances
            // below stay linear and close without iteration. At order 2, c = 0 and
            // b = 6 (T_bar - T_m).
            const double X = x_collocation[static_cast<size_t>(i) * n_z + j];
            // b = b_const + b_slope * (T_bar - T_m)
            const double b_const = quartic ? (X / 3.0) : 0.0;
            const double b_slope = quartic ? 10.0 : 6.0;
            // q_bush = kh (b - a) = kh (b_const + b_slope (T_bar - T_m) - T_b + T_s)
            //        = (T_b - T_amb) G
            // with T_m = (T_s + T_b)/2:
            //   kh [b_const + b_slope T_bar - (b_slope/2)(T_s + T_b) - T_b + T_s] = (T_b - T_amb) G
            //   T_b [kh (b_slope/2 + 1) + G] = kh [b_const + b_slope T_bar + T_s (1 - b_slope/2)] + G T_amb
            const size_t idx = static_cast<size_t>(i) * n_z + j;
            if (!bush_is_2d(cfg)) {
                const double denom = kh * (0.5 * b_slope + 1.0) + bush_g;
                alpha[idx] = (kh * (b_const + b_slope * t_bar) + bush_g * t_bush_outer) / denom;
                beta[idx] = kh * (1.0 - 0.5 * b_slope) / denom;
            }  // CONDUCTION_2D filled alpha/beta from the annulus solve in step 2a
            // q_shaft = kh (a + b) = kh [T_b - T_s + b_const + b_slope (T_bar - T_m)]
            //         = kh [b_const + b_slope T_bar + T_b (1 - b_slope/2) - T_s (1 + b_slope/2)]
            sums[0] += kh * area * (1.0 + 0.5 * b_slope);                    // coefficient of T_s
            sums[1] += kh * area * (b_const + b_slope * t_bar);              // constant part
            sums[2] += kh * area * (1.0 - 0.5 * b_slope) * alpha[idx];       // T_b constant part
            sums[3] += kh * area * (1.0 - 0.5 * b_slope) * beta[idx];        // T_b slope in T_s
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // sum q_shaft dA = sums[1] + sums[2] + T_s (sums[3] - sums[0])
    //                 = shaft_sink (T_s - T_shaft_ambient)
    const double lhs_slope = sums[3] - sums[0] - shaft_sink;
    double t_shaft = std::numeric_limits<double>::quiet_NaN();
    if (std::abs(lhs_slope) > 1.0e-30) {
        t_shaft = -(sums[1] + sums[2] + shaft_sink * t_shaft_ambient) / lhs_slope;
    }
    require_globally_valid(!std::isfinite(t_shaft) || t_shaft <= 0.0,
                           "Cross-film shaft temperature is nonphysical or nonfinite");
    state.t_shaft = t_shaft;

    // ---------------------------------------------------------------------
    // Step 3: profile coefficients, F-integrals, and the wall closure.
    // ---------------------------------------------------------------------
    double bush_area_sum = 0.0, bush_temp_sum = 0.0;
    double q_shaft_net = 0.0, q_bush_net = 0.0, core_rise_max = 0.0;
    double profile_t_min = std::numeric_limits<double>::infinity();
    double profile_t_max = -std::numeric_limits<double>::infinity();
    bool invalid_profile = false;
    // Physical cells only. The Reynolds assembly reads f2 and couette_h on both
    // theta ghost layers, which the caller supplies by ghost exchange rather than
    // by re-deriving a stencil that would need three ghost layers of p and T.
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const double h_val = std::max(h(i, j), cfg.min_film_thickness);
            const double k_th = std::max(field_or(fields, "k_liquid_solution",
                                                  cfg.thermal_conductivity, i, j), 1.0e-30);
            const double kh = k_th / h_val;
            const double t_bar = temperature(i, j);

            const size_t idx = static_cast<size_t>(i) * n_z + j;
            const double t_bush = alpha[idx] + beta[idx] * t_shaft;
            if (!std::isfinite(t_bush) || t_bush <= 0.0) {
                invalid_profile = true;
                continue;
            }

            const double dT_wall = t_bar - 0.5 * (t_shaft + t_bush);
            const double X = x_collocation[static_cast<size_t>(i) * n_z + j];
            const double a_val = t_bush - t_shaft;
            const double c_val = quartic ? -(5.0 / 3.0) * X - 20.0 * dT_wall : 0.0;
            const double b_val = 6.0 * dT_wall - 0.2 * c_val;
            if (!std::isfinite(a_val) || !std::isfinite(b_val) || !std::isfinite(c_val)) {
                invalid_profile = true;
                continue;
            }
            const auto extrema = profile_extrema(t_bar, a_val, b_val, c_val);
            if (!std::isfinite(extrema[0]) || !std::isfinite(extrema[1]) || extrema[0] <= 0.0) {
                invalid_profile = true;
                continue;
            }
            profile_t_min = std::min(profile_t_min, extrema[0]);
            profile_t_max = std::max(profile_t_max, extrema[1]);
            a_f(i, j) = a_val;
            b_f(i, j) = b_val;
            c_f(i, j) = c_val;
            t_bush_f(i, j) = t_bush;
            t_shaft_f(i, j) = t_shaft;
            core_rise_f(i, j) = a_val * psi1(0.5) + b_val * psi2(0.5) + c_val * psi4(0.5);

            // Film-surface outputs, evaluated from the same profile the
            // F-integrals below integrate. T(0) and T(1) reproduce the conjugate
            // wall temperatures because the wall values are what closes a and b;
            // they are written from the profile rather than copied so a future
            // change to the shape functions cannot leave the two disagreeing.
            //
            // The fluxes flip the module's internal sign convention: q_shaft and
            // q_bush are positive for heat LEAVING the film (PHYSICS.md 13.1.2),
            // while the output fields are positive INTO the fluid, which is the
            // convention a reader expects of a field named after the surface.
            if (t_journal_f) {
                (*t_journal_f)(i, j) =
                    t_bar + a_val * psi1(0.0) + b_val * psi2(0.0) + c_val * psi4(0.0);
            }
            if (t_bearing_f) {
                (*t_bearing_f)(i, j) =
                    t_bar + a_val * psi1(1.0) + b_val * psi2(1.0) + c_val * psi4(1.0);
            }
            if (q_journal_f) (*q_journal_f)(i, j) = -kh * (a_val + b_val);
            if (q_bearing_f) (*q_bearing_f)(i, j) = -kh * (b_val - a_val);

            // F-integrals with mu(T(eta)). The cell's own mu field already carries
            // the pressure, dissolved-gas and mixture corrections, so only the
            // temperature ratio is applied on top of it; that keeps every property
            // model working unchanged through the cross-film reconstruction.
            const double mu_cell = field_or(fields, "mu", cfg.mu, i, j);
            const double dissolved = field_or(fields, "dissolved_gas", 0.0, i, j);
            const double p_cell = p(i, j);
            double mu_ref = std::numeric_limits<double>::quiet_NaN();
            try {
                mu_ref = FluidProperties::liquid_solution_viscosity(dissolved, p_cell, t_bar, cfg);
            } catch (const std::exception&) {
                invalid_profile = true;
                continue;
            }
            if (!std::isfinite(mu_ref) || mu_ref <= 0.0) {
                invalid_profile = true;
                continue;
            }
            double int0 = 0.0, int1 = 0.0, int2 = 0.0;
            bool invalid_viscosity = false;
            for (size_t q = 0; q < nodes.size(); ++q) {
                const double e = nodes[q];
                const double t_local = profile_temperature(t_bar, a_val, b_val, c_val, e);
                double ratio = std::numeric_limits<double>::quiet_NaN();
                try {
                    ratio = FluidProperties::liquid_solution_viscosity(dissolved, p_cell, t_local, cfg) / mu_ref;
                } catch (const std::exception&) {
                    invalid_viscosity = true;
                    break;
                }
                if (!std::isfinite(ratio) || ratio <= 0.0) {
                    invalid_viscosity = true;
                    break;
                }
                const double inv_mu = 1.0 / (mu_cell * ratio);
                int0 += weights[q] * inv_mu;
                int1 += weights[q] * e * inv_mu;
                int2 += weights[q] * e * e * inv_mu;
            }
            if (invalid_viscosity) {
                invalid_profile = true;
                continue;
            }
            const double F0 = h_val * int0;
            const double F1 = h_val * h_val * int1;
            const double G2 = h_val * h_val * h_val * int2;
            const double F2 = G2 - F1 * F1 / F0;
            if (!std::isfinite(F0) || F0 <= 0.0 || !std::isfinite(F1) ||
                !std::isfinite(F2) || F2 <= 0.0) {
                invalid_profile = true;
                continue;
            }
            f0_f(i, j) = F0;
            f1_f(i, j) = F1;
            f2_f(i, j) = F2;
            couette_h_f(i, j) = 2.0 * F1 / F0;

            // Wall closure for the depth-averaged energy equation. The total loss
            // is q_shaft + q_bush = 2 kh b. Linearising in T_bar through the bush
            // relation (T_s is global and stays lagged) gives the derivative
            //   d(2 kh b)/dT_bar = 2 kh b_slope [1 - (b_slope/2) dT_b/dT_bar],
            // which is well below the 2 kh b_slope an explicit-wall linearisation
            // would use, and is what keeps the adiabatic-bush limit from being
            // over-damped.
            const double b_slope = quartic ? 10.0 : 6.0;
            // dT_b/dT_bar for the Newton linearisation only; the converged answer
            // always uses the solved T_b.
            //
            // LOCAL: the closed form of the per-cell bush relation.
            //
            // CONDUCTION_2D: the bush is taken as FROZEN, dT_b/dT_bar = 0, which
            // is the fully-sunk-wall value 2 (k/h) b_slope. That is very nearly
            // true here -- one film cell barely moves a bush its neighbours also
            // feed -- and it is the only choice that is both stable and free of
            // the bush mesh. Two alternatives were implemented and rejected:
            //   * the bore admittance with the neighbouring bush cells frozen
            //     reproduces this number to 5% (measured dT_b/dT_bar = 0.053 on
            //     the Ferron rig) but diverges as bush_radial_cells is refined,
            //     so the damping of the film energy equation would depend on the
            //     bush mesh;
            //   * the series-to-ambient conductance bush_g is mesh-free but far
            //     too small: on a nearly adiabatic bush it drives the whole
            //     bracket negative, the clamp below zeroes the wall coefficient,
            //     and with the wall term fully explicit the temperature KSP
            //     diverged outright on the Ferron rig at 720 x 192.
            const double denom = kh * (0.5 * b_slope + 1.0) + bush_g;
            const double dTb_dTbar = bush_is_2d(cfg) ? 0.0 : kh * b_slope / denom;
            const double q_wall = 2.0 * kh * b_val;
            wall_flux_f(i, j) = q_wall;
            wall_coeff_f(i, j) = std::max(
                2.0 * kh * b_slope * (1.0 - 0.5 * dTb_dTbar), 0.0);

            bush_area_sum += area;
            bush_temp_sum += area * t_bush;
            q_shaft_net += area * kh * (a_val + b_val);
            q_bush_net += area * kh * (b_val - a_val);
            core_rise_max = std::max(core_rise_max, core_rise_f(i, j));
        }
    }

    require_globally_valid(invalid_profile,
                           "Cross-film bush, profile, or viscosity integral is nonphysical");

    double reduce[4] = {bush_area_sum, bush_temp_sum, q_shaft_net, q_bush_net};
    MPI_Allreduce(MPI_IN_PLACE, reduce, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &core_rise_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &profile_t_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &profile_t_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    state.t_bush_mean = reduce[0] > 0.0 ? reduce[1] / reduce[0] : t_shaft;
    state.q_shaft_net = reduce[2];
    state.q_bush_net = reduce[3];
    state.core_rise_max = core_rise_max;
    state.profile_t_min = profile_t_min;
    state.profile_t_max = profile_t_max;
    return state;
}

}  // namespace CrossFilm
