#pragma once

#include <vector>

#include "config.hpp"
#include "communicator.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Mass-conserving cavitation as a complementarity problem, solved by the
/// Fischer-Burmeister-Newton-Schur (FBNS) algorithm of Woloszynski, Podsiadlo &
/// Stachowiak, "Efficient Solution to the Cavitation Problem in Hydrodynamic
/// Lubrication", Tribology Letters 58:18 (2015).
///
/// The unknowns are the pressure p and the liquid film content theta of every
/// cell. Instead of Elrod's barotropic switch, the JFO conditions are imposed
/// directly and pointwise:
///
///     p - p_cav >= 0,   theta_full - theta >= 0,   (p - p_cav)(theta_full - theta) = 0
///
/// which is equivalent to the smooth, unconstrained Fischer-Burmeister equation
/// FB(a, b) = a + b - sqrt(a^2 + b^2) = 0 with a = (p - p_cav)/p_ref and
/// b = theta_full - theta. There is no artificial bulk modulus anywhere in the
/// model, so `bulk_modulus` is ignored under this path.
///
/// theta_full is `FullFilm::bound`: 1 when the liquid may fill the whole gap, and
/// 1 - alpha_gas when a resolved free-gas phase (TWO_PHASE) occupies part of it.
/// alpha_gas is lagged -- the operator-split gas path updates it after the flow
/// solve -- so it is a coefficient of this Newton solve, not an unknown, and both
/// Jacobian blocks keep the forms they have without gas. See
/// `docs/FBNS_GAS_COUPLING.md`.
///
/// Required fields: "pressure", "theta", "h", "rho"; optional "mu",
/// "rho_liquid_solution", "dh_dt". Theta ghost cells of h, mu and the density
/// fields must be synced by Communicator before calling; the solver syncs the
/// unknowns itself.
///
/// Thermal / variable-property coupling is by frozen-coefficient Picard around
/// the Newton solve, exactly as the Elrod path does it. Every call reads the
/// current per-cell "mu" and "rho_liquid_solution" fields when it builds the
/// face mobilities, and holds them fixed for the duration of its own Newton
/// iteration; d(mu)/dT and d(mu)/dp never enter the Jacobian, so the Jacobian
/// stays exactly consistent with the residual it is differentiating and Newton
/// keeps its quadratic local rate. The caller's outer loop (see main.cpp) closes
/// the mu(T, p) feedback by re-solving until pressure and temperature stop
/// moving.
///
/// The cross-film model (`cross_film_thermal_model = POLYNOMIAL`) rides on the
/// same split. The face mobilities then come from the generalized (Dowson)
/// Reynolds equation of PHYSICS.md 13.1.3 -- Poiseuille rho_l F2 and Couette
/// thickness 2 F1/F0, i.e. mu(T(y)) integrated across the gap rather than a
/// depth-averaged mu -- and `CrossFilm::update` rebuilds those integrals once per
/// outer iteration, before this call. Inside one Newton solve they are frozen
/// coefficients like mu, so both Jacobian blocks and the residual read the same
/// values and the quadratic rate is preserved. `ReynoldsCoeff` holds the single
/// definition the Elrod path, this path and `Diagnostics` all use.
///
/// Scope (enforced in SimulationConfig::validate): fluid_property_model =
/// CONSTANT, SINGLE_PHASE or TWO_PHASE, gas_pressure_coupling = NONE,
/// DIRICHLET/NEUMANN axial boundaries. gas_pressure_coupling = VOID_COUPLED is
/// still refused: it makes the cavitation plateau a per-cell field that only
/// `Reynolds::solve_elrod`'s outer loop knows how to advance, so FBNS would run
/// it frozen at its initial value.
namespace FBNS {

    /// Convergence summary of the Newton iteration.
    /// `residual_history` holds max(|G|, |FB|) evaluated at every iterate,
    /// starting with the initial guess, so quadratic convergence is visible.
    struct Stats {
        int newton_iters = 0;
        double residual = 0.0;            // final max(|G|, |FB|), scaled
        double residual_mass = 0.0;       // final max |G| (mass balance, scaled)
        double residual_complementarity = 0.0;  // final max |FB|
        double max_theta_change = 0.0;    // max |delta theta| of the last accepted step
        bool converged = false;
        std::vector<double> residual_history;
    };

    /// True when FBNS is the selected cavitation model.
    bool active(const SimulationConfig& cfg);

    /// Pressure scale p_ref [Pa] used to nondimensionalize (p - p_cav) inside the
    /// Fischer-Burmeister function. Returns cfg.fbns_pressure_scale when it is
    /// positive, otherwise the Couette pressure scale 6 mu |omega| R^2 / c^2
    /// raised to cover any supplied overpressure.
    double pressure_scale(const SimulationConfig& cfg);

    /// Signed liquid mass flux [kg/s], positive toward increasing theta / z.
    /// theta index i denotes its east face; z index j denotes its north face,
    /// with j=-1 and j=n_z-1 the south and north physical boundaries.
    double mass_flux_theta(const Fields& fields, const Mesh& mesh,
                           const SimulationConfig& cfg, int i, int j);
    double mass_flux_z(const Fields& fields, const Mesh& mesh,
                       const SimulationConfig& cfg, int i, int j);

    /// Net liquid mass outflow [kg/s] of one cell under the FBNS operator:
    /// Poiseuille diffusion of p, upwind Couette transport of theta, and the
    /// backward-Euler storage of (rho h theta). Interior faces telescope between
    /// neighbours, so summing over the domain leaves only the axial boundary
    /// faces and the storage terms. `dt <= 0` drops the storage terms.
    /// `include_axial_boundary_faces = false` mirrors an inlet row that has
    /// overwritten an axial boundary link.
    ///
    /// This is the residual the Newton iteration drives to zero, and the same
    /// expression Diagnostics::mass_residual reports per cell.
    double cell_outflow(const Fields& fields, const Mesh& mesh,
                        const SimulationConfig& cfg, double dt, int i, int j,
                        bool include_axial_boundary_faces = true);

    /// Solve the complementarity problem. Sets "pressure", "theta" and "rho"
    /// (= rho_liquid * theta) and leaves both unknowns with valid ghost cells.
    Stats solve(Fields& fields, LinearSystem& sys, Communicator& comm,
                const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace FBNS
