#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "output_naming.hpp"

// Cavitation / film-rupture model, ordered simplest -> most physical:
//   FULL_SOMMERFELD - no cavitation; sub-cavity (tensile) pressures are kept.
//                     A reference baseline for verification / overlay plots only.
//   GUMBEL          - half-Sommerfeld: solve, then clamp p >= p_cav. Fast and
//                     robust but not mass-conserving (cavitation extent approximate).
//   JFO             - mass-conserving Jakobsson-Floberg-Olsson cavitation via the
//                     Elrod-Adams film-content algorithm. Physical rupture AND
//                     reformation; the production model.
//   FBNS            - the same JFO conditions written as a complementarity problem
//                     in (p, theta) and solved by the Fischer-Burmeister-Newton-Schur
//                     algorithm of Woloszynski, Podsiadlo & Stachowiak (2015).
//                     No artificial bulk modulus and no switch function: the free
//                     boundary is placed by p >= p_cav, theta <= theta_full and
//                     (p-p_cav)(theta_full-theta) = 0, with theta_full = 1 without a
//                     resolved gas phase and 1 - alpha_gas with one (FullFilm::bound).
//                     Scope: ISOTHERMAL or ENERGY_EQUATION, fluid_property_model
//                     CONSTANT / SINGLE_PHASE / TWO_PHASE. See PHYSICS.md 7.3.
enum class CavitationModel { FULL_SOMMERFELD, GUMBEL, JFO, FBNS, COUPLED };
enum class BCType { DIRICHLET, NEUMANN, INLET_OUTLET };
/// Thermal treatment of fluid entering through a pressure boundary.
///   ZERO_GRADIENT - submerged / open end: the entering oil is the same
///                   recirculating film oil, modelled as zero-gradient temperature
///                   (no externally fixed temperature is carried in).
///   CONSTANT      - fed supply / reservoir end: entering oil carries a fixed
///                   (constant) supply temperature temperature_reference.
enum class ThermalInflowMode { ZERO_GRADIENT, CONSTANT };
enum class MotionModel { STATIC, MOVING_BEARING, STEADY_EQUILIBRIUM };
// Integrator for the bearing equation of motion only (field solves are
// backward-Euler), ordered by increasing order/accuracy. EULER_IMPLICIT (stable)
// and RK4 (accurate) are maintained; the others are kept as references.
enum class TimeSteppingMethod { EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICOLSON, RK2, RK4 };
enum class SolutionMode { TRANSIENT, STEADY_STATE };
// Gas treatment in STEADY_STATE: FROZEN holds the gas state (gas_dt = 0, the
// historical behaviour); EQUILIBRIUM sets the dissolved fraction to the local
// saturation value each outer iteration (an equilibrium flash).
enum class SteadyGasModel { FROZEN, EQUILIBRIUM };
// Under-relaxation schedule of the outer (Picard) coupling loop:
//   FIXED    - coupling_relaxation is applied unchanged every sweep (historical).
//   ADAPTIVE - the factor backs off whenever the coupling residual RISES and
//              recovers while it falls. Targets the limit cycle the conjugate
//              (CONDUCTION_2D bush) cases show, where the loop gain sits just
//              above 1 and a fixed factor either over- or under-damps.
enum class CouplingRelaxationMode { FIXED, ADAPTIVE };
enum class TemperatureModel { ISOTHERMAL, ENERGY_EQUATION };
// Lubricant property fidelity, ordered simplest -> fullest:
//   CONSTANT     - fixed density and viscosity (plain oil).
//   SINGLE_PHASE - oil + dissolved gas: liquid mu/rho vary with the dissolved
//                  fraction (and T, p); gas stays in solution (no free phase).
//   TWO_PHASE    - adds finite-rate gas release/resorption and a free-gas phase
//                  in cavitated zones (the full gaseous-cavitation model).
enum class FluidPropertyModel { CONSTANT, SINGLE_PHASE, TWO_PHASE };
enum class DissolvedGasSpecies { AIR, PROPANE, R410A };
// Saturation / solubility law c_sat(p,T), simplest -> measured:
//   HENRY - c_sat = H(T) p (linear; optional van't Hoff T-dependence).
//   TABLE - measured c_sat(p,T) surface (1-D or 2-D PTSV data).
enum class SolubilityModel { HENRY, TABLE, CURVE_TABLE };
// Liquid-solution density, simplest -> measured:
//   CONSTANT - dissolved gas does not change liquid density (rho = rho_oil).
//   MIXTURE  - mass/volume mixing of oil and the dissolved-gas liquid phase;
//              uses both phase densities (rho and dissolved_gas_liquid_density).
//   TABLE    - measured rho(p,T).
enum class DensityModel { CONSTANT, MIXTURE, TABLE, LINEAR_TEMPERATURE };
// Liquid-solution viscosity, simplest -> measured:
//   CONSTANT  - dissolved gas does not change liquid viscosity (mu = mu_oil).
//   EMPIRICAL - mu_oil*exp[E_mu(1/T-1/T_ref)]*exp[a_c c_d]*exp[alpha_p(p-p_ref)]
//               (Andrade temperature x dissolved-gas thinning x Barus piezo-viscosity).
//   POLYNOMIAL- mu_oil*(k0 - k1 Tbar + k2 Tbar^2), Ferron, Frene & Boncompain
//               (1983). Tbar is a ratio of CELSIUS temperatures, T[C]/T_ref[C];
//               see viscosity_poly_k0 and PHYSICS.md. Pure temperature law: it
//               carries no pressure or dissolved-gas dependence.
//   TABLE     - measured kinematic nu(p,T), converted to dynamic per cell.
enum class LiquidViscosityModel { CONSTANT, EMPIRICAL, POLYNOMIAL, TABLE, COMPOSITION_TABLE, KINEMATIC_TABLE };
// Per-timestep terminal/log verbosity. Both modes print every timestep; DETAILED
// adds per-step diagnostics (film thickness, temperature range, fluid force, torque).
enum class OutputVerbosity { SIMPLIFIED, DETAILED };
// Face interpolation for convective terms, ordered by accuracy / complexity:
//   UPWIND            - first-order upwind (most diffusive, unconditionally bounded).
//   LINEAR            - second-order central interpolation (unbounded; can oscillate
//                       at sharp fronts). A textbook reference scheme.
//   VANLEER           - van Leer TVD limited high-order (sharp AND bounded).
//   TYPE_DIFFERENCING - Vijayaraghavan & Keith (1989): central interpolation across
//                       full-film faces, upwind across/inside cavitated faces.
//                       Only meaningful for the Elrod film-content (JFO) equation.
enum class ConvectionScheme { UPWIND, LINEAR, VANLEER, TYPE_DIFFERENCING };
// Krylov method for every inner linear solve (pressure, film content, temperature,
// gas transport). BCGS is the default; BCGSL is the more stable (and more
// expensive) BiCGStab(ell) variant, GMRES/FGMRES are restarted alternatives.
enum class KrylovMethod { BCGS, BCGSL, GMRES, FGMRES };
// Preconditioner for every inner linear solve:
//   ASM     - additive Schwarz, subdomain ILU(0), overlap set by linear_pc_overlap.
//             Keeps the theta coupling across rank boundaries; the default because
//             it is the only choice here whose convergence does not degrade with
//             rank count on the convection-dominated operators.
//   BJACOBI - block Jacobi, subdomain ILU(0). Identical to ASM at overlap 0 and
//             cheaper per iteration; the pre-2026-08 behaviour. Known to diverge
//             at some rank counts on the steady energy and Elrod equations.
//   JACOBI  - point Jacobi. Diagnostic only.
//   NONE    - unpreconditioned. Diagnostic only.
enum class LinearPreconditioner { ASM, BJACOBI, JACOBI, NONE };
// Effective film viscosity once free gas is present (all satisfy mu(alpha=0)=mu_l),
// ordered legacy -> recommended:
//   DUKLER  - void-fraction-weighted linear blend: alpha mu_g + (1-alpha) mu_l.
//             Legacy/deprecated; kept for comparison with prior results.
//   LINEAR  - mass-quality-weighted linear blend: x mu_g + (1-x) mu_l (Grando 2006);
//             refrigerant-oil validated.
//   MCADAMS - homogeneous two-phase standard: 1/mu = x/mu_g + (1-x)/mu_l with mass
//             quality x = m_g/(m_g + rho_l theta h). Correct asymptotes; the default.
enum class MixtureViscosityModel { DUKLER, LINEAR, MCADAMS };
// Cavitation onset/plateau pressure threshold:
//   SCALAR_PCAV - the configured oil cavitation pressure p_cav everywhere (default;
//                 reproduces prior behaviour exactly).
//   LOCAL_PSAT  - max(p_cav, p_sat) where p_sat is the dissolved-species bubble
//                 point at the feed state, so whichever species releases first
//                 (the higher release pressure) governs onset. General: any second
//                 fluid with a higher saturation pressure than the oil takes over,
//                 not specific to the R290/PZ68S pair. A uniform datum, so it shifts
//                 the cavitation extent without injecting a spurious net force.
enum class CavitationThreshold { SCALAR_PCAV, LOCAL_PSAT };
// Coupling of released free gas into the film pressure problem (WP-1):
//   NONE         - released gas is volumetrically inert (default; rho = rho_l*theta)
//                  on the Elrod (JFO) path. Note that under cavitation_model = FBNS
//                  the free gas always displaces liquid: the complementarity bound
//                  is theta <= 1 - alpha_gas whatever this key says, because that is
//                  the constraint itself rather than a pressure closure on top of it
//                  (PHYSICS.md 7.3.7). This key only selects the Elrod closure.
//   VOID_COUPLED - gas displaces liquid (full-film ceiling theta_full = 1 - alpha_g),
//                  sets a local plateau p_void, softens the film via a mixture bulk
//                  modulus beta_bar, and enters the mixture density rho = rho_l*theta
//                  + rho_g*alpha_g. Gaseous rupture at p_sat then emerges mechanistically.
enum class GasPressureCoupling { NONE, VOID_COUPLED };
// Thermal/shear treatment of the cavitated zone:
//   FULL_FILM - the whole gap shears and carries thermal mass (default; the
//               historical behaviour, overstates power loss as theta -> theta_min).
//   STRIATED  - in cavitated cells only the liquid fraction theta carries Couette
//               shear, thermal capacity, and enthalpy; full-film cells unchanged.
enum class CavitatedFilmModel { FULL_FILM, STRIATED };
// Asperity contact (mixed lubrication): NONE keeps smooth surfaces;
// GREENWOOD_TRIPP adds p_c = K E* F_5/2(h/sigma) and boundary shear mu_b p_c (contact.hpp).
enum class ContactModel { NONE, GREENWOOD_TRIPP };
// Cross-film (through-thickness) temperature treatment:
//   NONE       - the film is isothermal across the gap: viscosity is taken at the
//                depth-averaged temperature and wall heat loss is Newton cooling
//                to prescribed journal/bearing temperatures (historical default).
//   POLYNOMIAL - a low-order polynomial T(y) is reconstructed per cell (Lund &
//                Tonnesen style), supplying the wall conduction fluxes, the
//                conjugate wall temperatures, and the Dowson viscosity integrals
//                F0/F1/F2 that replace h^3/(12 mu) in the Reynolds equation.
// This is the single highest-level switch for the whole cross-film model: with
// NONE nothing in the cross-film module runs and the results are bit-identical
// to the pre-model solver.
enum class CrossFilmModel { NONE, POLYNOMIAL, FINITE_VOLUME };
// How the bush closes the film-to-bush heat flux (only read when the cross-film
// model is active):
//   LOCAL        - per-cell 1-D radial conduction through the wall in series with
//                  outer-surface convection (historical default). Every watt that
//                  enters the bore leaves through the wall directly underneath it,
//                  so the bush surface can only peak where the film is hottest.
//   CONDUCTION_2D- steady 2-D (r, theta) conduction in the bush annulus, Ferron's
//                  own bush model. Heat picked up in the loaded arc can leave at a
//                  different angular position, which is the circumferential
//                  redistribution a local closure has no representation of.
enum class BushThermalModel { LOCAL, CONDUCTION_2D, CONDUCTION_3D };
// Outer-surface boundary condition of the 2-D bush (also selects the outer limit
// of the LOCAL closure's series resistance):
//   CONVECTION        - Robin, bush_outer_heat_transfer to thermal_ambient_temperature
//                       (h = 0 is an adiabatic outer surface). Historical default.
//   FIXED_TEMPERATURE - the h -> infinity limit, outer surface held at
//                       thermal_ambient_temperature.
enum class BushOuterBoundary { CONVECTION, FIXED_TEMPERATURE };
// Thermal condition imposed on the cells of a fed inlet groove:
//   PINNED        - groove cells are held at t_supply (historical default).
//   FIXED_LAMBDA  - T = lambda*t_supply + (1-lambda)*T_upstream, lambda per inlet.
//   FLOW_WEIGHTED - the same mix, with lambda = Q_supply/(Q_supply + Q_upstream)
//                   from the inlet's supply flow rate and the solved upstream
//                   circumferential film flow.
//   FLUX_BALANCED - keep the groove film temperature unknown and exchange mass
//                   with the supply at the net FBNS face-flux imbalance.
enum class GrooveMixingModel { PINNED, FIXED_LAMBDA, FLOW_WEIGHTED, FLUX_BALANCED };

struct PropertyTablePoint {
    double x = 0.0;
    double value = 0.0;
};

/// 2-D property table on a rectangular (pressure, temperature) grid, used to
/// ingest measured PTSV data (e.g. data/R290_PZ68S/r290_pz68s.csv). Axes are
/// strictly increasing; values are row-major: value(ip,it) = values[ip*t.size()+it].
/// Everything is stored in SI (Pa, K, and the property's SI unit). The grid
/// resolution is whatever the source file contains (editable post-compile).
struct PropertyTable2D {
    std::vector<double> p;       // pressure axis [Pa], strictly increasing
    std::vector<double> t;       // temperature axis [K], strictly increasing
    std::vector<double> values;  // row-major, size p.size()*t.size()

    bool empty() const { return p.empty() || t.empty() || values.empty(); }
    double at(std::size_t ip, std::size_t it) const { return values[ip * t.size() + it]; }

    /// Bilinear interpolation with edge clamping; returns fallback when empty.
    double interpolate(double pq, double tq, double fallback) const {
        if (empty()) return fallback;
        std::size_t ip0, ip1, it0, it1;
        double fp, ft;
        bracket(p, pq, ip0, ip1, fp);
        bracket(t, tq, it0, it1, ft);
        const double v0 = at(ip0, it0) + ft * (at(ip0, it1) - at(ip0, it0));
        const double v1 = at(ip1, it0) + ft * (at(ip1, it1) - at(ip1, it0));
        return v0 + fp * (v1 - v0);
    }

    /// Property value along the tq-isotherm at pressure node ip (interp in T only).
    double isotherm_at(std::size_t ip, double tq) const {
        std::size_t it0, it1;
        double ft;
        bracket(t, tq, it0, it1, ft);
        return at(ip, it0) + ft * (at(ip, it1) - at(ip, it0));
    }

    /// Bracket q on a strictly-increasing axis: outputs the lo/hi node indices
    /// and the interpolation fraction; clamps to the end nodes outside the range.
    static void bracket(const std::vector<double>& axis, double q,
                        std::size_t& lo, std::size_t& hi, double& frac) {
        if (q <= axis.front()) { lo = hi = 0; frac = 0.0; return; }
        if (q >= axis.back()) { lo = hi = axis.size() - 1; frac = 0.0; return; }
        const auto upper = std::upper_bound(axis.begin(), axis.end(), q);
        hi = static_cast<std::size_t>(upper - axis.begin());
        lo = hi - 1;
        const double span = axis[hi] - axis[lo];
        frac = span > 0.0 ? (q - axis[lo]) / span : 0.0;
    }
};

/// Property curves with a separate temperature grid at each pressure or mass fraction.
/// Values and coordinates are SI. Interpolate in T, then the outer axis; optionally
/// interpolate log(value). CLAMP holds endpoints, otherwise out-of-range queries fail.
struct PropertyCurveTable {
    std::vector<double> axis;
    std::vector<std::vector<PropertyTablePoint>> curves;

    bool empty() const { return axis.empty() || curves.empty(); }

    double isotherm_at(std::size_t index, double temperature, bool clamp, bool logarithmic = false) const {
        const auto& curve = curves.at(index);
        if (curve.empty() || !std::isfinite(temperature) || temperature <= 0 ||
            (!clamp && (temperature < curve.front().x || temperature > curve.back().x)))
            throw std::domain_error("Solution curve table: temperature outside digitized range");
        const auto value = [&](double v) { return logarithmic ? std::log(v) : v; };
        auto right = std::lower_bound(curve.begin(), curve.end(), temperature,
            [](const auto& point, double t) { return point.x < t; });
        if (right == curve.end()) return value(curve.back().value);
        if (right == curve.begin() || right->x == temperature) return value(right->value);
        const auto& left = *(right - 1);
        const double weight = (temperature - left.x) / (right->x - left.x);
        return std::lerp(value(left.value), value(right->value), weight);
    }

    double interpolate(double x, double temperature, bool clamp, bool logarithmic = false) const {
        if (empty() || !std::isfinite(x) ||
            (!clamp && (x < axis.front() || x > axis.back())))
            throw std::domain_error("Solution curve table: pressure/concentration outside digitized range");
        std::size_t left, right;
        double weight;
        PropertyTable2D::bracket(axis, x, left, right, weight);
        // An exact node must not query an unused neighbouring temperature curve.
        if (weight == 0.0) right = left;
        if (weight == 1.0) left = right;
        const double a = isotherm_at(left, temperature, clamp, logarithmic);
        const double result = left == right ? a : std::lerp(a,
            isotherm_at(right, temperature, clamp, logarithmic), weight);
        return logarithmic ? std::exp(result) : result;
    }
};

/// Dynamic liquid viscosity as a function of actual composition and temperature.
/// Log interpolation preserves positivity. No extrapolation is permitted.
struct CompositionViscosityTable {
    std::vector<double> concentration, temperature, log_viscosity;
    bool empty() const { return concentration.empty() || temperature.empty() || log_viscosity.empty(); }
    bool contains(double c,double t) const {
        return !empty() && c>=concentration.front() && c<=concentration.back() &&
            t>=temperature.front() && t<=temperature.back();
    }
    double interpolate(double c,double t) const {
        if (!contains(c,t)) throw std::domain_error("Composition viscosity table queried outside its measured c,T range");
        std::size_t c0,c1,t0,t1; double fc,ft;
        PropertyTable2D::bracket(concentration,c,c0,c1,fc);
        PropertyTable2D::bracket(temperature,t,t0,t1,ft);
        const auto at=[&](std::size_t i,std::size_t j) { return log_viscosity[i*temperature.size()+j]; };
        return std::exp((1-fc)*((1-ft)*at(c0,t0)+ft*at(c0,t1))+fc*((1-ft)*at(c1,t0)+ft*at(c1,t1)));
    }
};

/// Invert a monotone-rising-then-flat saturation curve c_sat(p): return the
/// smallest p at which c_sat reaches c_d (the bubble point). A first-crossing
/// scan from low p makes it robust to the flat saturated tail and to small
/// non-monotone noise in measured VLE data. p and c are equal length with p
/// strictly increasing; returns 0 for a degenerate table.
inline double invert_csat(const std::vector<double>& p,
                          const std::vector<double>& c, double c_d) {
    if (p.size() < 2 || p.size() != c.size()) return 0.0;
    if (c_d <= c.front()) return p.front();
    if (c_d > c.back()) return p.back();  // strictly above range; '==' falls to the scan
    // (so c_d at the saturated plateau returns the onset p, not the tail end)
    for (std::size_t idx = 1; idx < p.size(); ++idx) {
        if (c[idx] < c_d) continue;            // c_d not yet reached
        const double span = c[idx] - c[idx - 1];
        if (span <= 0.0) return p[idx];
        const double frac = (c_d - c[idx - 1]) / span;
        return p[idx - 1] + frac * (p[idx] - p[idx - 1]);
    }
    return p.back();
}

inline std::string normalise_config_token(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    text = text.substr(first, last - first + 1);
    std::replace(text.begin(), text.end(), '-', '_');
    std::replace(text.begin(), text.end(), ' ', '_');
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return text;
}

inline const char* to_config_value(BCType type) {
    switch (type) {
        case BCType::DIRICHLET: return "DIRICHLET";
        case BCType::NEUMANN: return "NEUMANN";
        case BCType::INLET_OUTLET: return "INLET_OUTLET";
    }
    return "DIRICHLET";
}

inline const char* to_config_value(OutputVerbosity mode) {
    switch (mode) {
        case OutputVerbosity::SIMPLIFIED: return "SIMPLIFIED";
        case OutputVerbosity::DETAILED: return "DETAILED";
    }
    return "SIMPLIFIED";
}

inline const char* to_config_value(ThermalInflowMode mode) {
    switch (mode) {
        case ThermalInflowMode::ZERO_GRADIENT: return "ZERO_GRADIENT";
        case ThermalInflowMode::CONSTANT: return "CONSTANT";
    }
    return "ZERO_GRADIENT";
}

inline const char* to_config_value(CavitationModel model) {
    switch (model) {
        case CavitationModel::FULL_SOMMERFELD: return "FULL_SOMMERFELD";
        case CavitationModel::GUMBEL: return "GUMBEL";
        case CavitationModel::JFO: return "JFO";
        case CavitationModel::FBNS: return "FBNS";
        case CavitationModel::COUPLED: return "COUPLED";
    }
    return "JFO";
}

inline const char* to_config_value(ConvectionScheme scheme) {
    switch (scheme) {
        case ConvectionScheme::UPWIND: return "UPWIND";
        case ConvectionScheme::LINEAR: return "LINEAR";
        case ConvectionScheme::VANLEER: return "VANLEER";
        case ConvectionScheme::TYPE_DIFFERENCING: return "TYPE_DIFFERENCING";
    }
    return "UPWIND";
}

inline const char* to_config_value(KrylovMethod method) {
    switch (method) {
        case KrylovMethod::BCGS:   return "BCGS";
        case KrylovMethod::BCGSL:  return "BCGSL";
        case KrylovMethod::GMRES:  return "GMRES";
        case KrylovMethod::FGMRES: return "FGMRES";
    }
    return "BCGS";
}

inline const char* to_config_value(LinearPreconditioner pc) {
    switch (pc) {
        case LinearPreconditioner::ASM:     return "ASM";
        case LinearPreconditioner::BJACOBI: return "BJACOBI";
        case LinearPreconditioner::JACOBI:  return "JACOBI";
        case LinearPreconditioner::NONE:    return "NONE";
    }
    return "ASM";
}

inline const char* to_config_value(MixtureViscosityModel model) {
    switch (model) {
        case MixtureViscosityModel::DUKLER: return "DUKLER";
        case MixtureViscosityModel::LINEAR: return "LINEAR";
        case MixtureViscosityModel::MCADAMS: return "MCADAMS";
    }
    return "MCADAMS";
}

inline const char* to_config_value(MotionModel model) {
    switch (model) {
        case MotionModel::STATIC: return "STATIC";
        case MotionModel::MOVING_BEARING: return "MOVING_BEARING";
        case MotionModel::STEADY_EQUILIBRIUM: return "STEADY_EQUILIBRIUM";
    }
    return "STATIC";
}

inline const char* to_config_value(TimeSteppingMethod method) {
    switch (method) {
        case TimeSteppingMethod::EULER_EXPLICIT: return "EULER_EXPLICIT";
        case TimeSteppingMethod::EULER_IMPLICIT: return "EULER_IMPLICIT";
        case TimeSteppingMethod::CRANK_NICOLSON: return "CRANK_NICOLSON";
        case TimeSteppingMethod::RK2: return "RK2";
        case TimeSteppingMethod::RK4: return "RK4";
    }
    return "EULER_IMPLICIT";
}

inline const char* to_config_value(SolutionMode mode) {
    switch (mode) {
        case SolutionMode::TRANSIENT: return "TRANSIENT";
        case SolutionMode::STEADY_STATE: return "STEADY_STATE";
    }
    return "TRANSIENT";
}

inline const char* to_config_value(TemperatureModel model) {
    switch (model) {
        case TemperatureModel::ISOTHERMAL: return "ISOTHERMAL";
        case TemperatureModel::ENERGY_EQUATION: return "ENERGY_EQUATION";
    }
    return "ISOTHERMAL";
}

inline const char* to_config_value(CrossFilmModel model) {
    switch (model) {
        case CrossFilmModel::NONE: return "NONE";
        case CrossFilmModel::POLYNOMIAL: return "POLYNOMIAL";
        case CrossFilmModel::FINITE_VOLUME: return "FINITE_VOLUME";
    }
    return "NONE";
}

inline const char* to_config_value(BushThermalModel model) {
    switch (model) {
        case BushThermalModel::LOCAL: return "LOCAL";
        case BushThermalModel::CONDUCTION_2D: return "CONDUCTION_2D";
        case BushThermalModel::CONDUCTION_3D: return "CONDUCTION_3D";
    }
    return "LOCAL";
}

inline const char* to_config_value(BushOuterBoundary boundary) {
    switch (boundary) {
        case BushOuterBoundary::CONVECTION: return "CONVECTION";
        case BushOuterBoundary::FIXED_TEMPERATURE: return "FIXED_TEMPERATURE";
    }
    return "CONVECTION";
}

inline const char* to_config_value(GrooveMixingModel model) {
    switch (model) {
        case GrooveMixingModel::PINNED: return "PINNED";
        case GrooveMixingModel::FIXED_LAMBDA: return "FIXED_LAMBDA";
        case GrooveMixingModel::FLOW_WEIGHTED: return "FLOW_WEIGHTED";
        case GrooveMixingModel::FLUX_BALANCED: return "FLUX_BALANCED";
    }
    return "PINNED";
}

inline const char* to_config_value(FluidPropertyModel model) {
    switch (model) {
        case FluidPropertyModel::CONSTANT: return "CONSTANT";
        case FluidPropertyModel::SINGLE_PHASE: return "SINGLE_PHASE";
        case FluidPropertyModel::TWO_PHASE: return "TWO_PHASE";
    }
    return "CONSTANT";
}

inline const char* to_config_value(DissolvedGasSpecies species) {
    switch (species) {
        case DissolvedGasSpecies::AIR: return "AIR";
        case DissolvedGasSpecies::PROPANE: return "PROPANE";
        case DissolvedGasSpecies::R410A: return "R410A";
    }
    return "PROPANE";
}

inline const char* to_config_value(SolubilityModel model) {
    switch (model) {
        case SolubilityModel::HENRY: return "HENRY";
        case SolubilityModel::TABLE: return "TABLE";
        case SolubilityModel::CURVE_TABLE: return "CURVE_TABLE";
    }
    return "HENRY";
}

inline const char* to_config_value(DensityModel model) {
    switch (model) {
        case DensityModel::CONSTANT: return "CONSTANT";
        case DensityModel::MIXTURE: return "MIXTURE";
        case DensityModel::TABLE: return "TABLE";
        case DensityModel::LINEAR_TEMPERATURE: return "LINEAR_TEMPERATURE";
    }
    return "CONSTANT";
}

inline const char* to_config_value(LiquidViscosityModel model) {
    switch (model) {
        case LiquidViscosityModel::CONSTANT: return "CONSTANT";
        case LiquidViscosityModel::EMPIRICAL: return "EMPIRICAL";
        case LiquidViscosityModel::POLYNOMIAL: return "POLYNOMIAL";
        case LiquidViscosityModel::TABLE: return "TABLE";
        case LiquidViscosityModel::COMPOSITION_TABLE: return "COMPOSITION_TABLE";
        case LiquidViscosityModel::KINEMATIC_TABLE: return "KINEMATIC_TABLE";
    }
    return "CONSTANT";
}

inline KrylovMethod krylov_method_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "BCGSL" || value == "BICGSTABL" || value == "BCGS_L") return KrylovMethod::BCGSL;
    if (value == "GMRES") return KrylovMethod::GMRES;
    if (value == "FGMRES") return KrylovMethod::FGMRES;
    return KrylovMethod::BCGS;
}

inline LinearPreconditioner linear_preconditioner_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "BJACOBI" || value == "BLOCK_JACOBI") return LinearPreconditioner::BJACOBI;
    if (value == "JACOBI" || value == "POINT_JACOBI") return LinearPreconditioner::JACOBI;
    if (value == "NONE") return LinearPreconditioner::NONE;
    return LinearPreconditioner::ASM;
}

inline ConvectionScheme convection_scheme_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "LINEAR" || value == "CENTRAL" || value == "CD") return ConvectionScheme::LINEAR;
    // TVD_MINMOD was removed; the legacy token maps to the surviving van Leer limiter.
    if (value == "VANLEER" || value == "VAN_LEER" || value == "TVD_VANLEER" || value == "TVD" ||
        value == "TVD_MINMOD" || value == "MINMOD") {
        return ConvectionScheme::VANLEER;
    }
    if (value == "TYPE_DIFFERENCING" || value == "TYPE" || value == "HYBRID_TYPE" ||
        value == "VK1989") {
        return ConvectionScheme::TYPE_DIFFERENCING;
    }
    return ConvectionScheme::UPWIND;
}

inline MixtureViscosityModel mixture_viscosity_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "DUKLER" || value == "DUKLER_VOID") return MixtureViscosityModel::DUKLER;
    if (value == "LINEAR" || value == "LINEAR_QUALITY" || value == "GRANDO") {
        return MixtureViscosityModel::LINEAR;
    }
    // EINSTEIN_DILUTE and KRIEGER_DOUGHERTY were removed; MCADAMS is the default
    // and the survivor for any legacy/unknown token.
    return MixtureViscosityModel::MCADAMS;
}

inline MotionModel motion_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "STEADY_EQUILIBRIUM") return MotionModel::STEADY_EQUILIBRIUM;
    if (value == "MOVING_BEARING" || value == "DYNAMIC" || value == "BEARING") {
        return MotionModel::MOVING_BEARING;
    }
    return MotionModel::STATIC;
}

inline TimeSteppingMethod time_stepping_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "EULER_EXPLICIT" || value == "EXPLICIT_EULER" || value == "EXPLICIT") {
        return TimeSteppingMethod::EULER_EXPLICIT;
    }
    if (value == "EULER_IMPLICIT" || value == "IMPLICIT_EULER" || value == "IMPLICIT") {
        return TimeSteppingMethod::EULER_IMPLICIT;
    }
    if (value == "CRANK_NICOLSON" || value == "CRANK_NICHOLSON" ||
        value == "SEMI_IMPLICIT" || value == "SEMIIMPLICIT" || value == "IMEX") {
        return TimeSteppingMethod::CRANK_NICOLSON;
    }
    if (value == "RK2" || value == "RUNGE_KUTTA_2" || value == "MIDPOINT") {
        return TimeSteppingMethod::RK2;
    }
    if (value == "RK4" || value == "RUNGE_KUTTA_4") {
        return TimeSteppingMethod::RK4;
    }
    return TimeSteppingMethod::EULER_IMPLICIT;
}

inline SolutionMode solution_mode_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "STEADY_STATE" || value == "STEADY" || value == "SINGLE_STEP") {
        return SolutionMode::STEADY_STATE;
    }
    return SolutionMode::TRANSIENT;
}

inline OutputVerbosity output_verbosity_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "DETAILED" || value == "DETAIL" || value == "VERBOSE" || value == "FULL") {
        return OutputVerbosity::DETAILED;
    }
    return OutputVerbosity::SIMPLIFIED;  // SIMPLIFIED, SIMPLE, BRIEF, COMPACT, QUIET
}

inline ThermalInflowMode thermal_inflow_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "CONSTANT" || value == "RESERVOIR" || value == "FIXED" || value == "SUPPLY" ||
        value == "INFLOW_TEMPERATURE" || value == "DIRICHLET") {
        return ThermalInflowMode::CONSTANT;
    }
    // ZERO_GRADIENT, OPEN, SAME_OIL, ADIABATIC, NEUMANN all map to ZERO_GRADIENT.
    return ThermalInflowMode::ZERO_GRADIENT;
}

inline TemperatureModel temperature_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "ENERGY_EQUATION" || value == "ENERGY" || value == "THERMAL" || value == "THD") {
        return TemperatureModel::ENERGY_EQUATION;
    }
    return TemperatureModel::ISOTHERMAL;
}

inline CrossFilmModel cross_film_model_from_config(std::string value) {
    if (value == "FINITE_VOLUME" || value == "finite_volume") return CrossFilmModel::FINITE_VOLUME;
    value = normalise_config_token(value);
    if (value == "POLYNOMIAL" || value == "LUND_TONNESEN" || value == "PROFILE" ||
        value == "CROSS_FILM" || value == "ON" || value == "TRUE") {
        return CrossFilmModel::POLYNOMIAL;
    }
    return CrossFilmModel::NONE;  // NONE, OFF, FALSE, UNIFORM, DEPTH_AVERAGED
}

inline BushThermalModel bush_thermal_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "CONDUCTION_3D") return BushThermalModel::CONDUCTION_3D;
    if (value == "CONDUCTION_2D" || value == "CONDUCTION2D" || value == "2D" ||
        value == "FERRON" || value == "ANNULUS") {
        return BushThermalModel::CONDUCTION_2D;
    }
    return BushThermalModel::LOCAL;  // LOCAL, 1D, RADIAL
}

inline BushOuterBoundary bush_outer_boundary_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "FIXED_TEMPERATURE" || value == "FIXED" || value == "DIRICHLET" ||
        value == "ISOTHERMAL") {
        return BushOuterBoundary::FIXED_TEMPERATURE;
    }
    return BushOuterBoundary::CONVECTION;  // CONVECTION, ROBIN, HTC
}

inline GrooveMixingModel groove_mixing_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "FIXED_LAMBDA" || value == "FIXED" || value == "LAMBDA") {
        return GrooveMixingModel::FIXED_LAMBDA;
    }
    if (value == "FLOW_WEIGHTED" || value == "FLOW" || value == "MIXING") {
        return GrooveMixingModel::FLOW_WEIGHTED;
    }
    if (value == "FLUX_BALANCED") return GrooveMixingModel::FLUX_BALANCED;
    return GrooveMixingModel::PINNED;  // PINNED, SUPPLY, NONE
}

inline FluidPropertyModel fluid_property_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "SINGLE_PHASE" || value == "OIL_DISSOLVED_GAS" || value == "DISSOLVED_GAS" ||
        value == "OIL_PROPANE" || value == "PROPANE_OIL") {
        return FluidPropertyModel::SINGLE_PHASE;
    }
    if (value == "TWO_PHASE" || value == "GAS_CAVITATION_MIXTURE" ||
        value == "CAVITATION_MIXTURE" || value == "FREE_GAS") {
        return FluidPropertyModel::TWO_PHASE;
    }
    return FluidPropertyModel::CONSTANT;
}

inline DissolvedGasSpecies dissolved_gas_species_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "AIR") return DissolvedGasSpecies::AIR;
    if (value == "PROPANE" || value == "R290") return DissolvedGasSpecies::PROPANE;
    if (value == "R410A") return DissolvedGasSpecies::R410A;
    throw std::invalid_argument("Unsupported dissolved_gas_species '" + value +
        "'. Supported species: AIR, PROPANE (R290), R410A.");
}

inline SolubilityModel solubility_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "CURVE_TABLE") return SolubilityModel::CURVE_TABLE;
    if (value == "TABLE" || value == "TABULATED") return SolubilityModel::TABLE;
    // BUNSEN was removed; it is the same linear-in-p law as Henry, so map to HENRY.
    return SolubilityModel::HENRY;
}

inline DensityModel density_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "LINEAR_TEMPERATURE") return DensityModel::LINEAR_TEMPERATURE;
    if (value == "MIXTURE" || value == "MASS_VOLUME_MIXING" || value == "MASS_VOLUME" ||
        value == "MIXING") {
        return DensityModel::MIXTURE;
    }
    if (value == "TABLE" || value == "TABULATED") return DensityModel::TABLE;
    return DensityModel::CONSTANT;  // CONSTANT, PURE_OIL
}

inline LiquidViscosityModel liquid_viscosity_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "KINEMATIC_TABLE") return LiquidViscosityModel::KINEMATIC_TABLE;
    if (value == "COMPOSITION_TABLE") return LiquidViscosityModel::COMPOSITION_TABLE;
    if (value == "TABLE" || value == "TABULATED") return LiquidViscosityModel::TABLE;
    if (value == "POLYNOMIAL" || value == "POLY" || value == "FERRON") {
        return LiquidViscosityModel::POLYNOMIAL;
    }
    // LOG_MIXING was removed (it is EMPIRICAL with the T and p coefficients set to 0).
    if (value == "EMPIRICAL" || value == "EMPIRICAL_CORRELATION" || value == "CORRELATION" ||
        value == "LOG_MIXING" || value == "LOG") {
        return LiquidViscosityModel::EMPIRICAL;
    }
    return LiquidViscosityModel::CONSTANT;  // CONSTANT, PURE_OIL
}

struct InletConfig {
    enum class Type { CIRCULAR, GROOVE } type = Type::CIRCULAR;
    double theta = 0.0;     // Center angular position (deg)
    double z = 0.0;         // Center axial position (m)
    double size = 0.0;      // Radius for CIRCULAR (m), angular width for GROOVE (deg)
    double p_supply = 0.0;  // Supply pressure (Pa)
    // Thermal behavior. When feeds_fresh_oil is false (default) the feature is an
    // "open region" (e.g. a large-clearance pocket of the same recirculating oil):
    // it constrains pressure only and is thermally transparent to the energy solve.
    // When true it is an actual fed inlet that pumps fresh oil in at t_supply.
    bool feeds_fresh_oil = false;
    double t_supply = 0.0;  // Fresh-oil supply temperature (K), used when feeds_fresh_oil.
    // Hot-oil carry-over. Fraction of the groove-cell temperature supplied as fresh
    // oil: T = mixing_lambda*t_supply + (1-mixing_lambda)*T_upstream_film. 1.0 is the
    // historical hard pin. Used only when groove_mixing_model = FIXED_LAMBDA.
    double mixing_lambda = 1.0;
    // Volumetric fresh-oil supply rate through this feature [m^3/s]. Used only when
    // groove_mixing_model = FLOW_WEIGHTED, where lambda = Q_supply/(Q_supply + Q_up).
    double supply_flow_rate = 0.0;
    double dissolved_gas = -1.0; // Incoming liquid composition; -1 inherits initial charge.
    double free_gas_mass_fraction = 0.0; // Incoming free species / total mixture mass.
};

struct SimulationConfig {
    // Geometry
    double R = 0.01;
    double c = 0.0001;
    double e = 0.00008;
    double L = 0.02;
    double attitude_angle_deg = -90.0;
    double load_angle_deg = 0.0;

    // Grid
    int n_theta_global = 120;
    int n_z_global = 40;

    // Time
    double end_t = 1.0;
    double dt = 0.01;
    double write_interval = 0.1;
    SolutionMode solution_mode = SolutionMode::TRANSIENT;
    // Field solves (pressure/theta, temperature, gas) are backward-Euler in time.
    // Only the bearing equation of motion has a selectable integrator.
    TimeSteppingMethod motion_time_method = TimeSteppingMethod::EULER_IMPLICIT;

    // Convection face interpolation per transported variable. UPWIND is the
    // backward-compatible default; TYPE_DIFFERENCING applies to theta only.
    ConvectionScheme theta_convection_scheme = ConvectionScheme::UPWIND;
    ConvectionScheme thermal_convection_scheme = ConvectionScheme::UPWIND;
    ConvectionScheme gas_convection_scheme = ConvectionScheme::UPWIND;

    // Relative tolerance for the inner Krylov solves. The global mass-balance
    // residual reported by the diagnostics is bounded by this tolerance.
    double linear_rtol = 1e-8;

    // Inner linear solver. Defaults are chosen so that convergence does not
    // depend on the MPI rank count; see LinearSolverOptions in linear_system.hpp.
    KrylovMethod linear_ksp_type = KrylovMethod::BCGS;
    LinearPreconditioner linear_pc_type = LinearPreconditioner::ASM;
    // Additive-Schwarz subdomain overlap in cells; ignored unless linear_pc_type = ASM.
    // 0 makes ASM equivalent to block Jacobi and reinstates the rank dependence.
    int linear_pc_overlap = 1;
    // Number of escalation rungs (0-2) tried when a linear solve fails before the
    // run is aborted. Costs nothing on a solve that converges. 0 restores the
    // abort-on-first-failure behaviour.
    int linear_solver_fallbacks = 2;

    // Physics
    double omega = 100.0;
    double omega_ramp_time = 0.0;  // Transient startup ramp duration [s]; <=0 means full speed immediately.
    double mu = 0.01;
    double rho = 900.0;
    double p_cav = 0.0;
    double oil_vapor_temperature_coeff = 0.0; // Clausius-Clapeyron slope [K], zero holds p_cav fixed.
    double oil_vapor_pressure(double temperature) const {
        return p_cav*std::exp(oil_vapor_temperature_coeff*(1/property_reference_temperature-1/temperature));
    }
    double bulk_modulus = 1e9;

    // Misalignment (tilt)
    double tilt_slope_x = 0.0; // Slope in x-direction (m/m)
    double tilt_slope_y = 0.0; // Slope in y-direction (m/m)

    // Asperity contact (Greenwood-Tripp); see contact.hpp and PHYSICS.md
    ContactModel contact_model = ContactModel::NONE;
    double contact_roughness_sigma = 0.4e-6;  // composite RMS roughness sigma [m]
    double contact_gt_k = 1.0e-3;             // K = (16 sqrt2 pi/15)(eta beta sigma)^2 sqrt(sigma/beta) [-]
    double contact_modulus = 1.099e11;        // Hertz composite modulus E* [Pa]
    double contact_boundary_friction = 0.07;  // boundary friction coefficient mu_b [-]

    // Energy equation / thermal model
    TemperatureModel temperature_model = TemperatureModel::ISOTHERMAL;
    double temperature_initial = 300.0;
    double temperature_reference = 300.0;
    double journal_wall_temperature = 300.0;
    double bearing_wall_temperature = 300.0;
    double rho_cp = 1.92e6;              // Volumetric heat capacity [J/(m^3 K)]
    double liquid_specific_heat = 0.0;   // [J/(kg K)]; positive overrides legacy rho_cp/rho
    double thermal_conductivity = 0.15;  // Lubricant thermal conductivity [W/(m K)]
    double gas_heat_capacity_cv = -1.0; // Ideal free gas [J/(kg K)]; R410A uses its EOS.
    double dissolved_heat_capacity = -1.0; // Dissolved species [J/(kg K)].
    double desorption_internal_energy = std::numeric_limits<double>::quiet_NaN(); // u_g-u_d at property reference T,p [J/kg].
    double gas_thermal_conductivity = 0.0; // Explicit dilute-phase conductivity [W/(m K)].
    double journal_heat_transfer = 0.0;  // Wall heat-transfer coefficient [W/(m^2 K)]
    double bearing_heat_transfer = 0.0;  // Wall heat-transfer coefficient [W/(m^2 K)]

    // Cross-film temperature profile, conjugate walls, generalized Reynolds.
    // cross_film_thermal_model is the single switch for all three: NONE (default)
    // leaves the solver exactly as it was; POLYNOMIAL turns on the reconstructed
    // T(y), the conjugate shaft/bush wall temperatures (which then replace the
    // journal_heat_transfer / bearing_heat_transfer Newton cooling entirely), and
    // the Dowson F0/F1/F2 viscosity integrals in the Reynolds coefficient.
    CrossFilmModel cross_film_thermal_model = CrossFilmModel::NONE;
    // Polynomial order of T(y): 2 (parabola from the two wall values and the
    // depth-averaged mean) or 4 (adds one midplane collocation of the cross-film
    // energy balance, which flattens the core where convection dominates).
    int cross_film_profile_order = 4;
    // Number of Gauss-Legendre points used for the F0/F1/F2 integrals across the gap.
    int cross_film_quadrature_points = 6;
    /// Number of equal finite-volume layers across the gap.
    int cross_film_layers = 16;
    // Under-relaxation of the order-4 midplane collocation source between outer
    // iterations. 1.0 (no relaxation) converges on the validation rigs; the lever
    // exists because the collocation is the one part of the model that reads a
    // gradient of the solution. Only affects the path to the fixed point.
    double cross_film_relaxation = 1.0;
    // Shaft: circumferentially and axially isothermal, its temperature solved from the global
    // film-to-shaft heat balance. This is a reduced model, not Ferron's r-z shaft
    // conduction problem. shaft_ambient_heat_transfer = 0 (default) makes the
    // shaft globally adiabatic, so the integral of the film-to-shaft flux is zero
    // and the shaft carries heat from the loaded arc back upstream.
    double shaft_ambient_heat_transfer = 0.0;   // [W/(m^2 K)] shaft to ambient
    double shaft_ambient_area_ratio = 1.0;      // exposed shaft area / bearing surface area
    double thermal_ambient_temperature = 300.0; // Ambient/sump temperature [K]
    // Optional shaft sink temperature [K]. <= 0 uses thermal_ambient_temperature.
    double shaft_ambient_temperature = 0.0;
    // Bush: 1-D radial conduction from the film surface through the bush wall to an
    // outer-surface convection boundary. bush_outer_heat_transfer = 0 (default)
    // makes the bush adiabatic at its outer surface, i.e. no net heat leaves the film
    // through the bush.
    double bush_thermal_conductivity = 50.0;    // [W/(m K)] bush material
    double bush_outer_radius = 0.0;             // [m]; <=0 means R + 0.1*R
    double bush_outer_heat_transfer = 0.0;      // [W/(m^2 K)] bush outer surface to ambient
    // Optional outer-bush boundary temperature [K]. <= 0 uses
    // thermal_ambient_temperature for both convection and fixed-temperature BCs.
    double bush_outer_temperature = 0.0;
    // CONDUCTION_2D replaces the local radial closure above by a steady 2-D
    // (r, theta) conduction solve in the annulus [R, bush_outer_radius], which is
    // what lets the bush move heat circumferentially. LOCAL (default) is
    // bit-identical to the solver before the option existed.
    BushThermalModel bush_thermal_model = BushThermalModel::LOCAL;
    // Radial cells of the 2-D bush mesh. The radial profile is logarithmic and the
    // discretisation is exact for it, so this only resolves the (weak) radial
    // variation of the circumferential spreading; 16 is ample.
    int bush_radial_cells = 16;
    int bush_axial_cells = 16; // Solid axial slabs, each coupled to an integral number of film cells.
    double bush_end_heat_transfer = 0.0; // Both annular faces; zero is insulated.
    bool bush_end_natural_convection = false; // Vertical-plate analogue in the same surrounding oil.
    // Outer-surface condition of the bush, for both closures.
    BushOuterBoundary bush_outer_boundary = BushOuterBoundary::CONVECTION;

    // Hot-oil carry-over at fed grooves. PINNED (default) reproduces the historical
    // hard pin; FLUX_BALANCED uses the solved FBNS mass flux as reservoir exchange.
    GrooveMixingModel groove_mixing_model = GrooveMixingModel::PINNED;

    // Liquid oil plus dissolved gas property model. Defaults keep pure-oil behavior.
    FluidPropertyModel fluid_property_model = FluidPropertyModel::CONSTANT;
    DissolvedGasSpecies dissolved_gas_species = DissolvedGasSpecies::PROPANE;
    SolubilityModel solubility_model = SolubilityModel::HENRY;
    DensityModel density_model = DensityModel::CONSTANT;
    double density_temperature_intercept = 954.9; // rho(T_C = 0) [kg/m^3]
    double density_temperature_slope = -0.7284;   // d(rho)/d(T_C) [kg/(m^3 K)]
    // CONSTANT density: end natural convection uses this only as beta = -slope/rho.
    LiquidViscosityModel liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
    double dissolved_gas_initial = 0.0;          // Dissolved gas mass fraction [kg gas / kg liquid solution]
    double dissolved_gas_max = 1.0;              // Upper bound for dissolved gas mass fraction
    double dissolved_gas_henry_coeff = 0.0;      // Henry slope H at the reference T [1/Pa]: c_sat = H(T) p
    double dissolved_gas_henry_temp_coeff = 0.0; // van't Hoff slope E_H [K]: H(T)=H_ref exp(E_H(1/T-1/T_ref))
    double dissolved_gas_liquid_density = 0.0;   // Partial liquid-volume density [kg/m^3]; required by MIXTURE
    double gas_mass_transfer_rate = 0.0;         // Release rate [1/s]; legacy absorption fallback.
    double gas_absorption_rate = -1.0;          // -1 inherits release rate; 0 disables absorption.
    bool coupled_phase_equilibrium = false;    // COUPLED: equilibrium instead of finite-rate exchange.
    /// Separate absorption kinetics; bubble compression alone exchanges no mass.
    double effective_gas_absorption_rate() const {
        return gas_absorption_rate >= 0 ? gas_absorption_rate : gas_mass_transfer_rate;
    }
    double dissolved_gas_diffusivity = 0.0;      // Effective Fickian diffusivity of dissolved gas [m^2/s]
    double solution_density_gas_coeff = 0.0;     // Optional density correction by concentration
    double solution_viscosity_gas_coeff = 0.0;   // Log-mixing viscosity exponent a_c [-]: mu_l = mu_oil exp(a_c c_d)
    double viscosity_temperature_coeff = 0.0;    // Andrade slope E_mu [K]: mu_oil(T)=mu exp(E_mu(1/T-1/T_ref))
    double viscosity_pressure_coeff = 0.0;       // Barus piezo-viscosity alpha_p [1/Pa]: mu *= exp(alpha_p(p-p_ref))
    // POLYNOMIAL (Ferron, Frene & Boncompain, ASME J. Lubr. Technol. 105 (1983) 422):
    //   mu(T) = mu * (k0 - k1 Tbar + k2 Tbar^2),  Tbar = T[degC] / T_ref[degC].
    // Tbar is a ratio of CELSIUS temperatures, not Kelvin. The paper's nomenclature
    // reads "T = temperature degC" and "Tbar = T/T0", and its own coefficients
    // (k0=3.287, k1=3.064, k2=0.777, mu_i=0.0277 Pa s at T_i=40 degC) sum to
    // k0-k1+k2 = 1.000 at Tbar=1 and to 0.267 at 80 degC, which is the measured
    // VG32 drop; the same coefficients on a Kelvin ratio give 0.82 and are wrong.
    // property_reference_temperature stays in KELVIN like every other temperature
    // in this config; the Celsius conversion happens inside the viscosity law.
    // Defaults (1, 0, 0) are the identity mu(T) = mu.
    double viscosity_poly_k0 = 1.0;              // Ferron k0 [-] (constant term)
    double viscosity_poly_k1 = 0.0;              // Ferron k1 [-] (linear term, subtracted)
    double viscosity_poly_k2 = 0.0;              // Ferron k2 [-] (quadratic term)
    double gas_alpha_max = 1.0;                  // Maximum resolved free-gas volume fraction
    double gas_pressure_floor = 1.0;             // Ideal-gas pressure floor [Pa]
    MixtureViscosityModel mixture_viscosity_model = MixtureViscosityModel::MCADAMS;
    GasPressureCoupling gas_pressure_coupling = GasPressureCoupling::NONE;  // released gas into pressure (WP-1)
    // Numerical relaxation of the Elrod p_void iteration and steady feed flash.
    // Physical transient phase transfer and EOS gas volume are never relaxed.
    double void_coupling_relaxation = 0.3;
    double mu_gas = 0.0;                         // Free-gas dynamic viscosity [Pa s]; required by DUKLER/MCADAMS/LINEAR models
    double property_reference_pressure = 101325.0;
    double property_reference_temperature = 300.0;
    std::vector<PropertyTablePoint> solubility_table;
    std::vector<PropertyTablePoint> density_table;
    std::vector<PropertyTablePoint> viscosity_table;
    // 2-D (pressure, temperature) property tables loaded from external CSV files.
    // When non-empty they take precedence over the 1-D tables above for the same
    // property. Solubility is the saturated mass fraction c_sat(p,T); viscosity is
    // the KINEMATIC viscosity nu(p,T) [m^2/s] (converted to dynamic per cell using
    // the local solution density); density is rho(p,T) [kg/m^3].
    PropertyTable2D solubility_table_2d;
    PropertyTable2D viscosity_table_2d;
    PropertyTable2D density_table_2d;
    std::filesystem::path solubility_table_file;
    std::filesystem::path viscosity_table_file;
    std::filesystem::path viscosity_composition_file;
    CompositionViscosityTable viscosity_composition_table;
    // Digitized c_mass_fraction,T_C,P_MPa,viscosity_mm2s data, including c=0 and c=1.
    std::filesystem::path solution_table_file;
    PropertyCurveTable solubility_curve_table; // c_sat(p_abs,T)
    PropertyCurveTable kinematic_composition_table; // nu(c,T) [m^2/s]
    bool solution_table_clamp = true;
    bool kinematic_viscosity_log_interpolation = true;
    std::filesystem::path density_table_file;
    // Directory of the loaded config file; relative *_table_file paths resolve
    // against it first, then the CWD. Set by load_from_file().
    std::filesystem::path config_dir;

    // Bearing motion: outer bearing moves while the shaft remains fixed.
    MotionModel motion_model = MotionModel::STATIC;
    double steady_force_tolerance = 0.002;
    double steady_position_step_limit = 0.05; // maximum correction / clearance
    double steady_position_fd_step = 0.001; // Jacobian perturbation / clearance
    bool bearing_initial_from_attitude = true;
    double bearing_initial_x = 0.0;
    double bearing_initial_y = 0.0;
    double bearing_initial_z = 0.0;
    double bearing_initial_vx = 0.0;
    double bearing_initial_vy = 0.0;
    double bearing_initial_vz = 0.0;
    double bearing_mass = 1.0;
    double bearing_stiffness_x = 0.0;
    double bearing_stiffness_y = 0.0;
    double bearing_stiffness_z = 0.0;
    double bearing_damping_x = 0.0;
    double bearing_damping_y = 0.0;
    double bearing_damping_z = 0.0;
    double external_load_x = 0.0;
    double external_load_y = 0.0;
    double external_load_z = 0.0;
    double min_film_thickness = 1e-9;
    bool stop_on_nonpositive_film = true;

    // Axial Boundary Conditions
    BCType bc_z_south_type = BCType::DIRICHLET;
    double bc_z_south_val  = 0.0;
    double bc_z_south_theta = 1.0; // Deprecated compatibility input; Elrod derives theta from bc_z_south_val.
    ThermalInflowMode bc_z_south_thermal = ThermalInflowMode::ZERO_GRADIENT; // Energy inflow temperature treatment
    BCType bc_z_north_type = BCType::DIRICHLET;
    double bc_z_north_val  = 0.0;
    double bc_z_north_theta = 1.0; // Deprecated compatibility input; Elrod derives theta from bc_z_north_val.
    ThermalInflowMode bc_z_north_thermal = ThermalInflowMode::ZERO_GRADIENT; // Energy inflow temperature treatment
    // Pressure reservoirs allow both flow directions. These states apply only
    // to inflow; outflow transports the adjacent interior state.
    double bc_z_south_temperature = 0.0; // 0 inherits temperature_reference [K].
    double bc_z_north_temperature = 0.0;
    double bc_z_south_dissolved_gas = -1.0; // -1 inherits dissolved_gas_initial.
    double bc_z_north_dissolved_gas = -1.0;
    double bc_z_south_free_gas_mass_fraction = 0.0; // Free gas / total mixture mass.
    double bc_z_north_free_gas_mass_fraction = 0.0;

    /// Temperature carried by fluid entering from the selected end reservoir.
    double axial_inflow_temperature(bool south) const {
        const double value = south ? bc_z_south_temperature : bc_z_north_temperature;
        return value > 0 ? value : temperature_reference;
    }
    /// Dissolved-species mass fraction within incoming liquid solution.
    double axial_inflow_dissolved_gas(bool south) const {
        const double value = south ? bc_z_south_dissolved_gas : bc_z_north_dissolved_gas;
        return value >= 0 ? value : dissolved_gas_initial;
    }
    /// Incoming free-gas mass per liquid-solution mass, for a common velocity.
    double axial_inflow_free_gas_ratio(bool south) const {
        const double value = south ? bc_z_south_free_gas_mass_fraction : bc_z_north_free_gas_mass_fraction;
        return value / (1.0 - value);
    }

    // Inlets
    std::vector<InletConfig> inlets;

    // Cavitation
    CavitationModel cavitation_model = CavitationModel::JFO;
    CavitationThreshold cavitation_threshold = CavitationThreshold::SCALAR_PCAV;
    CavitatedFilmModel cavitated_film_model = CavitatedFilmModel::FULL_FILM;
    // WP-11: at an axial boundary, a cavitated cell re-floods at the physical
    // Poiseuille reformation rate rather than the barotropic Gamma_base link
    // (which overstates re-flooding by ~beta/(p_bc-p_cav) since p = p_cav + beta
    // ln(theta) is invalid where the cell is cavitated). Default off = prior link.
    bool   consistent_boundary_flux = false;
    int    max_outer_iters = 50;    // Elrod-Adams outer flag-update iteration limit
    double outer_tol       = 1e-6;  // Convergence tolerance on max change between outer iters
    double theta_min       = 1e-6;  // Minimum film content
    bool   log_outer_iters = true;  // Whether to print outer iteration progress
    // Print the coupling residual of EVERY outer (Picard) sweep, not just the
    // last. Off by default (it is one line per sweep, up to coupling_max_iters
    // of them); on, it is the only way to see whether a loop is converging,
    // stagnating or limit-cycling, which is what the conjugate cases needed.
    bool   log_coupling_history = false;

    // FBNS (cavitation_model = FBNS). The complementarity solve has no artificial
    // bulk modulus, so bulk_modulus is unused; these replace max_outer_iters /
    // outer_tol for that path.
    int    fbns_max_iters = 100;      // Newton iteration limit per solve
    double fbns_tol       = 1e-9;     // max-norm tolerance on the scaled (G, FB) residuals
    double fbns_relaxation = 1.0;     // Newton step damping in (0, 1]; 1 = full step
    // Pressure scale p_ref [Pa] that nondimensionalizes (p - p_cav) inside the
    // Fischer-Burmeister function. The FB function weights its two arguments
    // equally, so feeding raw pascals against a film content in [0, 1] destroys
    // the conditioning. 0 = derive it from the operating point (the Couette
    // pressure scale 6 mu |omega| R^2 / c^2 and any supplied overpressure).
    double fbns_pressure_scale = 0.0;

    // Outer (Picard) coupling between Reynolds, properties, and energy. In
    // STEADY_STATE the coupled sweep is repeated until the normalized pressure
    // and temperature residuals fall below coupling_tolerance (a real fixed
    // point, not a single segregated sweep). TRANSIENT advances one sweep per
    // step as before.
    int    coupling_max_iters = 30;
    double coupling_tolerance = 1e-6;        // on normalized max pressure/temperature change
    // FBNS equation convergence replaces the pressure-change stopping gate only.
    bool steady_pressure_inner_only = false;
    double coupling_relaxation = 1.0;        // under-relaxation in (0, 1]; 1 = none
    // Best-iterate rollback. A loop that exhausts coupling_max_iters without
    // meeting the tolerance may end on an iterate that is WORSE than one it
    // already passed through (measured on the conjugate Ferron cases: the
    // residual bottoms out near outer 51 and then oscillates upward for the
    // remaining 150 sweeps). With this on, such a run delivers the lowest-residual
    // iterate instead of the last one and says so in the log. A run that MEETS
    // the tolerance always delivers its last iterate, bit-identically to before.
    bool   coupling_keep_best_iterate = true;
    // Adaptive under-relaxation of the outer loop. Defaults are inert: FIXED
    // reproduces the historical schedule exactly.
    CouplingRelaxationMode coupling_relaxation_mode = CouplingRelaxationMode::FIXED;
    double coupling_relaxation_min = 0.1;     // floor for the adaptive factor
    double coupling_relaxation_backoff = 0.5; // multiplier applied when the residual rises
    double coupling_relaxation_growth = 1.1;  // multiplier applied while it falls
    SteadyGasModel steady_gas_model = SteadyGasModel::FROZEN;
    bool   allow_frozen_gas_steady = false;  // opt in to running gas cavitation frozen in steady

    // Output
    std::string output_dir = "results";
    // Converged steady pure-oil FBNS fields used only as a new solve's initial guess.
    std::string steady_initial_state_file;
    std::string filename_prefix = "solution";
    OutputVerbosity output_verbosity = OutputVerbosity::SIMPLIFIED;
    bool output_write_3d = true;
    bool output_write_flat = true;
    std::vector<std::string> output_fields = default_output_fields();

    // Diagnostics: per-step global mass balances and solver-convergence logging
    // (results/<run>/diagnostics.csv). Interval in steps; <= 0 disables.
    int diagnostics_interval = 1;
    // Full steady iterate output cadence; zero disables intermediate snapshots.
    int steady_output_interval = 0;

    // Non-fatal notes produced while parsing (e.g. ignored legacy keys).
    std::vector<std::string> parse_warnings;

    // Fatal parse problems (currently: unrecognized keys). `validate()` promotes
    // these into its error list, so every caller that already refuses to run on a
    // validation error refuses to run on a misspelled key too.
    std::vector<std::string> parse_errors;

    void load_from_file(const std::filesystem::path& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open config file '" << path
                      << "'. Using defaults.\n";
            return;
        }

        config_dir = path.parent_path();
        load_from_stream(file);
    }

    void load_from_text(const std::string& text) {
        std::istringstream input(text);
        load_from_stream(input);
    }

    void load_from_stream(std::istream& input) {
        inlets.clear();
        parse_warnings.clear();
        parse_errors.clear();

        std::string line;
        while (std::getline(input, line)) {
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);

            trim(line);
            if (line.empty()) continue;

            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            std::string key = line.substr(0, eq_pos);
            std::string value = line.substr(eq_pos + 1);

            trim(key);
            trim(value);

            if (key.empty()) continue;
            const bool empty_value_allowed =
                key == "output_fields" || key == "solubility_table" ||
                key == "density_table" || key == "viscosity_table";
            if (value.empty() && !empty_value_allowed) continue;

            if (key == "R") R = std::stod(value);
            else if (key == "c") c = std::stod(value);
            else if (key == "e") e = std::stod(value);
            else if (key == "L") L = std::stod(value);
            else if (key == "attitude_angle_deg") attitude_angle_deg = std::stod(value);
            else if (key == "load_angle_deg") load_angle_deg = std::stod(value);
            else if (key == "n_theta_global") n_theta_global = std::stoi(value);
            else if (key == "n_z_global") n_z_global = std::stoi(value);
            else if (key == "end_t") end_t = std::stod(value);
            else if (key == "dt") dt = std::stod(value);
            else if (key == "write_interval") write_interval = std::stod(value);
            else if (key == "solution_mode") solution_mode = solution_mode_from_config(value);
            else if (key == "pressure_time_method" || key == "temperature_time_method") {
                parse_warnings.push_back(
                    key + " is ignored: field solves are backward-Euler in time.");
            }
            else if (key == "motion_time_method") motion_time_method = time_stepping_from_config(value);
            else if (key == "theta_convection_scheme") theta_convection_scheme = convection_scheme_from_config(value);
            else if (key == "thermal_convection_scheme") thermal_convection_scheme = convection_scheme_from_config(value);
            else if (key == "gas_convection_scheme") gas_convection_scheme = convection_scheme_from_config(value);
            else if (key == "linear_rtol") linear_rtol = std::stod(value);
            else if (key == "linear_ksp_type") linear_ksp_type = krylov_method_from_config(value);
            else if (key == "linear_pc_type") linear_pc_type = linear_preconditioner_from_config(value);
            else if (key == "linear_pc_overlap") linear_pc_overlap = std::stoi(value);
            else if (key == "linear_solver_fallbacks") linear_solver_fallbacks = std::stoi(value);
            else if (key == "omega") omega = std::stod(value);
            else if (key == "omega_ramp_time") omega_ramp_time = std::stod(value);
            else if (key == "mu") mu = std::stod(value);
            else if (key == "rho") rho = std::stod(value);
            else if (key == "p_cav") p_cav = std::stod(value);
            else if (key == "oil_vapor_temperature_coeff") oil_vapor_temperature_coeff = std::stod(value);
            else if (key == "bulk_modulus") bulk_modulus = std::stod(value);
            else if (key == "tilt_slope_x") tilt_slope_x = std::stod(value);
            else if (key == "tilt_slope_y") tilt_slope_y = std::stod(value);
            else if (key == "contact_model") {
                if (value == "NONE") contact_model = ContactModel::NONE;
                else if (value == "GREENWOOD_TRIPP") contact_model = ContactModel::GREENWOOD_TRIPP;
            }
            else if (key == "contact_roughness_sigma") contact_roughness_sigma = std::stod(value);
            else if (key == "contact_gt_k") contact_gt_k = std::stod(value);
            else if (key == "contact_modulus") contact_modulus = std::stod(value);
            else if (key == "contact_boundary_friction") contact_boundary_friction = std::stod(value);
            else if (key == "temperature_model") temperature_model = temperature_model_from_config(value);
            else if (key == "temperature_initial") temperature_initial = std::stod(value);
            else if (key == "temperature_reference") temperature_reference = std::stod(value);
            else if (key == "journal_wall_temperature") journal_wall_temperature = std::stod(value);
            else if (key == "bearing_wall_temperature") bearing_wall_temperature = std::stod(value);
            else if (key == "rho_cp") rho_cp = std::stod(value);
            else if (key == "liquid_specific_heat") liquid_specific_heat = std::stod(value);
            else if (key == "thermal_conductivity") thermal_conductivity = std::stod(value);
            else if (key == "gas_heat_capacity_cv") gas_heat_capacity_cv = std::stod(value);
            else if (key == "dissolved_heat_capacity") dissolved_heat_capacity = std::stod(value);
            else if (key == "desorption_internal_energy") desorption_internal_energy = std::stod(value);
            else if (key == "gas_thermal_conductivity") gas_thermal_conductivity = std::stod(value);
            else if (key == "journal_heat_transfer") journal_heat_transfer = std::stod(value);
            else if (key == "bearing_heat_transfer") bearing_heat_transfer = std::stod(value);
            else if (key == "cross_film_thermal_model") cross_film_thermal_model = cross_film_model_from_config(value);
            else if (key == "cross_film_layers") cross_film_layers = std::stoi(value);
            else if (key == "cross_film_profile_order") cross_film_profile_order = std::stoi(value);
            else if (key == "cross_film_quadrature_points") cross_film_quadrature_points = std::stoi(value);
            else if (key == "cross_film_relaxation") cross_film_relaxation = std::stod(value);
            else if (key == "shaft_ambient_heat_transfer") shaft_ambient_heat_transfer = std::stod(value);
            else if (key == "shaft_ambient_area_ratio") shaft_ambient_area_ratio = std::stod(value);
            else if (key == "thermal_ambient_temperature") thermal_ambient_temperature = std::stod(value);
            else if (key == "shaft_ambient_temperature") shaft_ambient_temperature = std::stod(value);
            else if (key == "bush_thermal_conductivity") bush_thermal_conductivity = std::stod(value);
            else if (key == "bush_outer_radius") bush_outer_radius = std::stod(value);
            else if (key == "bush_outer_heat_transfer") bush_outer_heat_transfer = std::stod(value);
            else if (key == "bush_outer_temperature") bush_outer_temperature = std::stod(value);
            else if (key == "bush_thermal_model") bush_thermal_model = bush_thermal_model_from_config(value);
            else if (key == "bush_radial_cells") bush_radial_cells = std::stoi(value);
            else if (key == "bush_axial_cells") bush_axial_cells = std::stoi(value);
            else if (key == "bush_end_heat_transfer") bush_end_heat_transfer = std::stod(value);
            else if (key == "bush_end_natural_convection") {
                if (is_true(value)) bush_end_natural_convection = true;
                else if (is_false(value)) bush_end_natural_convection = false;
                else parse_errors.push_back("bush_end_natural_convection must be true or false.");
            }
            else if (key == "bush_outer_boundary") bush_outer_boundary = bush_outer_boundary_from_config(value);
            else if (key == "groove_mixing_model") groove_mixing_model = groove_mixing_model_from_config(value);
            else if (key == "fluid_property_model") fluid_property_model = fluid_property_model_from_config(value);
            else if (key == "dissolved_gas_species") {
                try { dissolved_gas_species = dissolved_gas_species_from_config(value); }
                catch (const std::invalid_argument& ex) { parse_errors.push_back(ex.what()); }
            }
            else if (key == "solubility_model" || key == "oil_gas_solution_model") solubility_model = solubility_model_from_config(value);
            else if (key == "density_model") density_model = density_model_from_config(value);
            else if (key == "density_temperature_intercept") density_temperature_intercept = std::stod(value);
            else if (key == "density_temperature_slope") density_temperature_slope = std::stod(value);
            else if (key == "liquid_viscosity_model" || key == "viscosity_model") liquid_viscosity_model = liquid_viscosity_model_from_config(value);
            else if (key == "dissolved_gas_initial") dissolved_gas_initial = std::stod(value);
            else if (key == "dissolved_gas_max") dissolved_gas_max = std::stod(value);
            else if (key == "dissolved_gas_henry_coeff") dissolved_gas_henry_coeff = std::stod(value);
            else if (key == "dissolved_gas_henry_temp_coeff") dissolved_gas_henry_temp_coeff = std::stod(value);
            else if (key == "dissolved_gas_bunsen_coeff") {
                parse_warnings.push_back(
                    "dissolved_gas_bunsen_coeff is ignored: the BUNSEN solubility model was "
                    "removed (use solubility_model = HENRY or TABLE).");
            }
            else if (key == "dissolved_gas_liquid_density") dissolved_gas_liquid_density = std::stod(value);
            else if (key == "gas_mass_transfer_rate") gas_mass_transfer_rate = std::stod(value);
            else if (key == "gas_absorption_rate") gas_absorption_rate = std::stod(value);
            else if (key == "coupled_phase_equilibrium") coupled_phase_equilibrium = is_true(value);
            else if (key == "dissolved_gas_diffusivity") dissolved_gas_diffusivity = std::stod(value);
            else if (key == "solution_density_gas_coeff") solution_density_gas_coeff = std::stod(value);
            else if (key == "solution_viscosity_gas_coeff") solution_viscosity_gas_coeff = std::stod(value);
            else if (key == "viscosity_temperature_coeff") viscosity_temperature_coeff = std::stod(value);
            else if (key == "viscosity_pressure_coeff") viscosity_pressure_coeff = std::stod(value);
            else if (key == "viscosity_poly_k0") viscosity_poly_k0 = std::stod(value);
            else if (key == "viscosity_poly_k1") viscosity_poly_k1 = std::stod(value);
            else if (key == "viscosity_poly_k2") viscosity_poly_k2 = std::stod(value);
            else if (key == "gas_alpha_max") gas_alpha_max = std::stod(value);
            else if (key == "gas_pressure_floor") gas_pressure_floor = std::stod(value);
            else if (key == "mixture_viscosity_model" || key == "gas_mixture_viscosity_model") mixture_viscosity_model = mixture_viscosity_model_from_config(value);
            else if (key == "gas_pressure_coupling") {
                if (value == "VOID_COUPLED" || value == "VOID") gas_pressure_coupling = GasPressureCoupling::VOID_COUPLED;
                else if (value == "NONE") gas_pressure_coupling = GasPressureCoupling::NONE;
            }
            else if (key == "void_coupling_relaxation") void_coupling_relaxation = std::stod(value);
            else if (key == "mu_gas") mu_gas = std::stod(value);
            else if (key == "diagnostics_interval") diagnostics_interval = std::stoi(value);
            else if (key == "steady_output_interval") steady_output_interval = std::stoi(value);
            else if (key == "property_reference_pressure") property_reference_pressure = std::stod(value);
            else if (key == "property_reference_temperature") property_reference_temperature = std::stod(value);
            else if (key == "solubility_table") solubility_table = parse_property_table(value);
            else if (key == "density_table") density_table = parse_property_table(value);
            else if (key == "viscosity_table") viscosity_table = parse_property_table(value);
            else if (key == "solubility_table_file") solubility_table_file = value;
            else if (key == "viscosity_table_file") viscosity_table_file = value;
            else if (key == "viscosity_composition_file") viscosity_composition_file = value;
            else if (key == "solution_table_file") solution_table_file = value;
            else if (key == "solution_table_clamp") solution_table_clamp = is_true(value);
            else if (key == "kinematic_viscosity_log_interpolation") kinematic_viscosity_log_interpolation = is_true(value);
            else if (key == "density_table_file") density_table_file = value;
            else if (key == "motion_model") motion_model = motion_model_from_config(value);
            else if (key == "steady_force_tolerance") steady_force_tolerance = std::stod(value);
            else if (key == "steady_position_step_limit") steady_position_step_limit = std::stod(value);
            else if (key == "steady_position_fd_step") steady_position_fd_step = std::stod(value);
            else if (key == "bearing_initial_from_attitude") {
                if (is_true(value)) bearing_initial_from_attitude = true;
                else if (is_false(value)) bearing_initial_from_attitude = false;
            }
            else if (key == "bearing_initial_x") bearing_initial_x = std::stod(value);
            else if (key == "bearing_initial_y") bearing_initial_y = std::stod(value);
            else if (key == "bearing_initial_z") bearing_initial_z = std::stod(value);
            else if (key == "bearing_initial_vx") bearing_initial_vx = std::stod(value);
            else if (key == "bearing_initial_vy") bearing_initial_vy = std::stod(value);
            else if (key == "bearing_initial_vz") bearing_initial_vz = std::stod(value);
            else if (key == "bearing_mass") bearing_mass = std::stod(value);
            else if (key == "bearing_stiffness_x") bearing_stiffness_x = std::stod(value);
            else if (key == "bearing_stiffness_y") bearing_stiffness_y = std::stod(value);
            else if (key == "bearing_stiffness_z") bearing_stiffness_z = std::stod(value);
            else if (key == "bearing_damping_x") bearing_damping_x = std::stod(value);
            else if (key == "bearing_damping_y") bearing_damping_y = std::stod(value);
            else if (key == "bearing_damping_z") bearing_damping_z = std::stod(value);
            else if (key == "external_load_x") external_load_x = std::stod(value);
            else if (key == "external_load_y") external_load_y = std::stod(value);
            else if (key == "external_load_z") external_load_z = std::stod(value);
            else if (key == "min_film_thickness") min_film_thickness = std::stod(value);
            else if (key == "stop_on_nonpositive_film") {
                if (is_true(value)) stop_on_nonpositive_film = true;
                else if (is_false(value)) stop_on_nonpositive_film = false;
            }
            else if (key == "bc_z_south_type") {
                if (value == "DIRICHLET") bc_z_south_type = BCType::DIRICHLET;
                else if (value == "NEUMANN") bc_z_south_type = BCType::NEUMANN;
                else if (value == "INLET_OUTLET") bc_z_south_type = BCType::INLET_OUTLET;
            }
            else if (key == "bc_z_south_val") bc_z_south_val = std::stod(value);
            else if (key == "bc_z_south_theta") bc_z_south_theta = std::stod(value);
            else if (key == "bc_z_south_thermal") bc_z_south_thermal = thermal_inflow_from_config(value);
            else if (key == "bc_z_north_type") {
                if (value == "DIRICHLET") bc_z_north_type = BCType::DIRICHLET;
                else if (value == "NEUMANN") bc_z_north_type = BCType::NEUMANN;
                else if (value == "INLET_OUTLET") bc_z_north_type = BCType::INLET_OUTLET;
            }
            else if (key == "bc_z_north_val") bc_z_north_val = std::stod(value);
            else if (key == "bc_z_north_theta") bc_z_north_theta = std::stod(value);
            else if (key == "bc_z_north_thermal") bc_z_north_thermal = thermal_inflow_from_config(value);
            else if (key == "bc_z_south_temperature") bc_z_south_temperature = std::stod(value);
            else if (key == "bc_z_north_temperature") bc_z_north_temperature = std::stod(value);
            else if (key == "bc_z_south_dissolved_gas") bc_z_south_dissolved_gas = std::stod(value);
            else if (key == "bc_z_north_dissolved_gas") bc_z_north_dissolved_gas = std::stod(value);
            else if (key == "bc_z_south_free_gas_mass_fraction") bc_z_south_free_gas_mass_fraction = std::stod(value);
            else if (key == "bc_z_north_free_gas_mass_fraction") bc_z_north_free_gas_mass_fraction = std::stod(value);
            else if (key == "cavitation_model") {
                if (value == "FULL_SOMMERFELD" || value == "SOMMERFELD")
                    cavitation_model = CavitationModel::FULL_SOMMERFELD;
                else if (value == "GUMBEL") cavitation_model = CavitationModel::GUMBEL;
                else if (value == "JFO" || value == "ELROD_ADAMS" || value == "ELROD")
                    cavitation_model = CavitationModel::JFO;
                else if (value == "FBNS" || value == "COMPLEMENTARITY")
                    cavitation_model = CavitationModel::FBNS;
                else if (value == "COUPLED") cavitation_model = CavitationModel::COUPLED;
            }
            else if (key == "cavitation_threshold") {
                if (value == "LOCAL_PSAT" || value == "PSAT")
                    cavitation_threshold = CavitationThreshold::LOCAL_PSAT;
                else if (value == "SCALAR_PCAV" || value == "PCAV" || value == "SCALAR")
                    cavitation_threshold = CavitationThreshold::SCALAR_PCAV;
            }
            else if (key == "cavitated_film_model") {
                if (value == "STRIATED") cavitated_film_model = CavitatedFilmModel::STRIATED;
                else if (value == "FULL_FILM") cavitated_film_model = CavitatedFilmModel::FULL_FILM;
            }
            else if (key == "consistent_boundary_flux") {
                if (is_true(value)) consistent_boundary_flux = true;
                else if (is_false(value)) consistent_boundary_flux = false;
            }
            else if (key == "coupling_max_iters") coupling_max_iters = std::stoi(value);
            else if (key == "coupling_tolerance") coupling_tolerance = std::stod(value);
            else if (key == "steady_pressure_inner_only") {
                if (is_true(value)) steady_pressure_inner_only = true;
                else if (is_false(value)) steady_pressure_inner_only = false;
                else parse_errors.push_back("steady_pressure_inner_only must be true or false.");
            }
            else if (key == "coupling_relaxation") coupling_relaxation = std::stod(value);
            else if (key == "coupling_keep_best_iterate") {
                if (is_true(value)) coupling_keep_best_iterate = true;
                else if (is_false(value)) coupling_keep_best_iterate = false;
            }
            else if (key == "coupling_relaxation_mode") {
                if (value == "ADAPTIVE") coupling_relaxation_mode = CouplingRelaxationMode::ADAPTIVE;
                else if (value == "FIXED") coupling_relaxation_mode = CouplingRelaxationMode::FIXED;
            }
            else if (key == "coupling_relaxation_min") coupling_relaxation_min = std::stod(value);
            else if (key == "coupling_relaxation_backoff") coupling_relaxation_backoff = std::stod(value);
            else if (key == "coupling_relaxation_growth") coupling_relaxation_growth = std::stod(value);
            else if (key == "steady_gas_model") {
                if (value == "EQUILIBRIUM") steady_gas_model = SteadyGasModel::EQUILIBRIUM;
                else if (value == "FROZEN") steady_gas_model = SteadyGasModel::FROZEN;
            }
            else if (key == "allow_frozen_gas_steady") {
                if (is_true(value)) allow_frozen_gas_steady = true;
                else if (is_false(value)) allow_frozen_gas_steady = false;
            }
            else if (key == "max_outer_iters") max_outer_iters = std::stoi(value);
            else if (key == "outer_tol") outer_tol = std::stod(value);
            else if (key == "fbns_max_iters") fbns_max_iters = std::stoi(value);
            else if (key == "fbns_tol") fbns_tol = std::stod(value);
            else if (key == "fbns_relaxation") fbns_relaxation = std::stod(value);
            else if (key == "fbns_pressure_scale") fbns_pressure_scale = std::stod(value);
            else if (key == "theta_min") theta_min = std::stod(value);
            else if (key == "log_outer_iters") {
                if (is_true(value)) log_outer_iters = true;
                else if (is_false(value)) log_outer_iters = false;
            }
            else if (key == "log_coupling_history") {
                if (is_true(value)) log_coupling_history = true;
                else if (is_false(value)) log_coupling_history = false;
            }
            else if (key == "output_dir") output_dir = value;
            else if (key == "steady_initial_state_file") steady_initial_state_file = value;
            else if (key == "filename_prefix") filename_prefix = value;
            else if (key == "output_verbosity") output_verbosity = output_verbosity_from_config(value);
            else if (key == "output_write_3d") {
                if (is_true(value)) output_write_3d = true;
                else if (is_false(value)) output_write_3d = false;
            }
            else if (key == "output_write_flat") {
                if (is_true(value)) output_write_flat = true;
                else if (is_false(value)) output_write_flat = false;
            }
            else if (key == "output_fields") output_fields = parse_output_fields(value);
            else if (key == "inlet_circular") {
                InletConfig inlet;
                inlet.type = InletConfig::Type::CIRCULAR;
                std::stringstream ss(value);
                ss >> inlet.theta >> inlet.z >> inlet.size >> inlet.p_supply;
                // Optional 5th token = fresh-oil supply temperature -> actual fed inlet.
                // Optional 6th/7th = hot-oil carry-over lambda and supply flow rate.
                if (ss >> inlet.t_supply) inlet.feeds_fresh_oil = true;
                // Read into temporaries: a failed extraction zeroes its target in
                // C++11 and later, which would silently overwrite mixing_lambda = 1.
                double token = 0.0;
                if (ss >> token) inlet.mixing_lambda = token;
                if (ss >> token) inlet.supply_flow_rate = token;
                inlets.push_back(inlet);
            }
            else if (key == "inlet_groove") {
                InletConfig inlet;
                inlet.type = InletConfig::Type::GROOVE;
                std::stringstream ss(value);
                ss >> inlet.theta >> inlet.size >> inlet.p_supply;
                // Optional 4th token = fresh-oil supply temperature -> actual fed inlet.
                // Optional 5th/6th = hot-oil carry-over lambda and supply flow rate.
                if (ss >> inlet.t_supply) inlet.feeds_fresh_oil = true;
                double token = 0.0;
                if (ss >> token) inlet.mixing_lambda = token;
                if (ss >> token) inlet.supply_flow_rate = token;
                inlets.push_back(inlet);
            }
            else if (key == "inlet_species") {
                std::stringstream ss(value);
                std::size_t index; double concentration, quality;
                if (!(ss >> index >> concentration >> quality) || index >= inlets.size())
                    parse_errors.push_back("inlet_species requires an existing zero-based inlet index, dissolved fraction and free mass fraction");
                else {
                    inlets[index].dissolved_gas = concentration;
                    inlets[index].free_gas_mass_fraction = quality;
                }
            }
            // Terminal branch: an unrecognized key is a hard error, not silence.
            // Without it a misspelling parses as a no-op and the run proceeds with
            // the default -- which is how `cross_film_model = NONE` sat in the
            // Ferron THD case templates through a whole validation campaign while
            // the parser only ever knew `cross_film_thermal_model`.
            else {
                parse_errors.push_back(unknown_key_message(key));
            }
        }
    }

    void save_to_file(const std::filesystem::path& path) const {
        std::ofstream file(path, std::ios::trunc);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open config file for writing: " + path.string());
        }
        file << to_config_text();
    }

    void validate(std::vector<std::string>& errors, std::vector<std::string>& warnings) const {
        errors.clear();
        warnings.clear();

        // A key the parser did not recognize was never applied, so everything
        // below it is validating defaults the user never asked for. Report it
        // first and let the caller's existing error path stop the run.
        errors.insert(errors.end(), parse_errors.begin(), parse_errors.end());

        auto finite = [](double value) { return std::isfinite(value); };
        auto error = [&](const std::string& message) { errors.push_back(message); };
        auto warn = [&](const std::string& message) { warnings.push_back(message); };
        auto require_positive = [&](double value, const char* label) {
            if (!finite(value) || value <= 0.0) {
                error(std::string(label) + " must be a finite positive value.");
            }
        };
        auto require_nonnegative = [&](double value, const char* label) {
            if (!finite(value) || value < 0.0) {
                error(std::string(label) + " must be a finite non-negative value.");
            }
        };
        auto require_pressure = [&](BCType type, double pressure, const char* label) {
            if (type == BCType::NEUMANN) return;
            if (!finite(pressure) || pressure <= 0.0) {
                error(std::string(label) + " must be a positive absolute pressure in Pa.");
            } else if (finite(p_cav) && p_cav > 0.0 && pressure < effective_p_cav()) {
                error(std::string(label) + " is below the effective cavitation pressure. "
                      "The selected pressure-floor model cannot impose this reservoir; "
                      "LOCAL_PSAT uses the initial charge bubble point as a global floor. "
                      "Do not raise the reservoir pressure to hide this conflict.");
            }
        };
        auto require_table = [&](const std::vector<PropertyTablePoint>& table,
                                 const char* label) {
            if (table.size() < 2) {
                error(std::string(label) + " needs at least two x:value entries when TABLE is selected.");
            }
        };

        require_positive(R, "R");
        require_positive(c, "c");
        require_positive(L, "L");
        require_nonnegative(e, "e");
        if (finite(c) && finite(e) && c > 0.0 && e >= c) {
            error("e must be smaller than c so the film thickness stays positive.");
        }
        if (n_theta_global <= 0 || n_z_global <= 0) {
            error("n_theta_global and n_z_global must be positive.");
        }
        require_positive(dt, "dt");
        require_positive(write_interval, "write_interval");
        require_positive(linear_rtol, "linear_rtol");
        if (linear_pc_overlap < 0) {
            error("linear_pc_overlap must be non-negative.");
        }
        if (linear_solver_fallbacks < 0 || linear_solver_fallbacks > 2) {
            error("linear_solver_fallbacks must be between 0 and 2.");
        }
        if (linear_solver_fallbacks == 0) {
            warn("linear_solver_fallbacks = 0: a linear solve that fails once aborts the run "
                 "instead of being retried on a different Krylov path.");
        }
        if (linear_pc_type == LinearPreconditioner::ASM && linear_pc_overlap == 0) {
            warn("linear_pc_overlap = 0 makes ASM equivalent to block Jacobi: the linear "
                 "solves become rank-dependent and can diverge on the energy and Elrod "
                 "equations. Use 1 or more.");
        }
        if (linear_pc_type != LinearPreconditioner::ASM) {
            warn(std::string("linear_pc_type = ") + to_config_value(linear_pc_type) +
                 " drops the coupling across MPI rank boundaries; convergence of the "
                 "convection-dominated solves then depends on the rank count. ASM is the "
                 "rank-robust choice.");
        }
        if (solution_mode == SolutionMode::TRANSIENT) require_nonnegative(end_t, "end_t");
        require_positive(mu, "mu");
        require_positive(rho, "rho");
        require_positive(bulk_modulus, "bulk_modulus");
        require_positive(theta_min, "theta_min");
        require_positive(min_film_thickness, "min_film_thickness");
        if (contact_model == ContactModel::GREENWOOD_TRIPP) {
            require_positive(contact_roughness_sigma, "contact_roughness_sigma");
            require_positive(contact_gt_k, "contact_gt_k");
            require_positive(contact_modulus, "contact_modulus");
            require_nonnegative(contact_boundary_friction, "contact_boundary_friction");
        }
        const bool zero_pressure_pure_oil = cavitation_model == CavitationModel::FBNS &&
            cavitation_threshold == CavitationThreshold::SCALAR_PCAV &&
            fluid_property_model == FluidPropertyModel::SINGLE_PHASE &&
            gas_pressure_coupling == GasPressureCoupling::NONE &&
            dissolved_gas_initial == 0.0 && dissolved_gas_henry_coeff == 0.0 &&
            oil_vapor_temperature_coeff == 0.0;
        if (!finite(p_cav) || p_cav < 0.0 || (p_cav == 0.0 && !zero_pressure_pure_oil)) {
            error("p_cav must be a positive absolute pressure, except 0 Pa for gas-free "
                  "SINGLE_PHASE FBNS with SCALAR_PCAV. Reservoir pressure belongs in "
                  "bc_z_south_val / bc_z_north_val.");
        }
        if (!finite(omega_ramp_time) || omega_ramp_time < 0.0) {
            error("omega_ramp_time must be finite and non-negative.");
        }
        if (temperature_model == TemperatureModel::ENERGY_EQUATION) {
            require_positive(rho_cp, "rho_cp");
            require_nonnegative(liquid_specific_heat, "liquid_specific_heat");
            require_positive(thermal_conductivity, "thermal_conductivity");
            require_nonnegative(journal_heat_transfer, "journal_heat_transfer");
            require_nonnegative(bearing_heat_transfer, "bearing_heat_transfer");

            // A constant viscosity disables viscous thermal feedback. Density,
            // gas compressibility and phase equilibrium may still depend on T.
            if (fluid_property_model == FluidPropertyModel::CONSTANT) {
                warn("temperature_model = ENERGY_EQUATION with fluid_property_model = CONSTANT: "
                     "viscosity is frozen at mu, so viscosity feedback from temperature is disabled. "
                     "Set fluid_property_model to a variable-property model to enable it.");
            } else if (liquid_viscosity_model == LiquidViscosityModel::CONSTANT) {
                warn("temperature_model = ENERGY_EQUATION with liquid_viscosity_model = CONSTANT: "
                     "viscosity is frozen at mu, so viscosity feedback from temperature is disabled. "
                     "Density and phase properties can still affect pressure or load. "
                     "Use a temperature-dependent liquid_viscosity_model for viscosity feedback.");
            }

            if (journal_heat_transfer == 0.0 && bearing_heat_transfer == 0.0 &&
                cross_film_thermal_model == CrossFilmModel::NONE) {
                warn("journal_heat_transfer = bearing_heat_transfer = 0: the film is fully "
                     "adiabatic and the configured wall temperatures are inert. Quote any "
                     "temperature result together with this wall condition.");
            }
            if (cross_film_thermal_model != CrossFilmModel::NONE &&
                (journal_heat_transfer > 0.0 || bearing_heat_transfer > 0.0)) {
                warn("cross_film_thermal_model is active: the conjugate shaft/bush model "
                     "replaces the Newton-cooling wall closure entirely, so "
                     "journal_heat_transfer, bearing_heat_transfer, journal_wall_temperature "
                     "and bearing_wall_temperature are IGNORED. Use "
                     "shaft_ambient_heat_transfer, bush_outer_heat_transfer and "
                     "thermal_ambient_temperature instead.");
            }
        }
        if (cross_film_thermal_model != CrossFilmModel::NONE) {
            if (cross_film_thermal_model == CrossFilmModel::FINITE_VOLUME) {
                if (cross_film_layers < 2 || cross_film_layers > 128)
                    error("cross_film_layers must be between 2 and 128.");
                if (solution_mode != SolutionMode::STEADY_STATE ||
                    cavitation_model != CavitationModel::FBNS ||
                    fluid_property_model != FluidPropertyModel::SINGLE_PHASE ||
                    groove_mixing_model != GrooveMixingModel::FLUX_BALANCED ||
                    thermal_convection_scheme != ConvectionScheme::UPWIND ||
                    liquid_specific_heat <= 0.0)
                    error("FINITE_VOLUME requires steady single-phase FBNS, FLUX_BALANCED, upwind energy and liquid_specific_heat > 0.");
            }
            if (temperature_model != TemperatureModel::ENERGY_EQUATION) {
                error("cross_film_thermal_model requires temperature_model = ENERGY_EQUATION; "
                      "there is no temperature field for it to build a profile from.");
            }
            if (cross_film_profile_order != 2 && cross_film_profile_order != 4) {
                error("cross_film_profile_order must be 2 or 4.");
            }
            if (cross_film_quadrature_points < 2 || cross_film_quadrature_points > 8) {
                error("cross_film_quadrature_points must be between 2 and 8.");
            }
            if (!finite(cross_film_relaxation) || cross_film_relaxation <= 0.0 ||
                cross_film_relaxation > 1.0) {
                error("cross_film_relaxation must be in (0, 1].");
            }
            require_positive(bush_thermal_conductivity, "bush_thermal_conductivity");
            require_nonnegative(shaft_ambient_heat_transfer, "shaft_ambient_heat_transfer");
            require_nonnegative(bush_outer_heat_transfer, "bush_outer_heat_transfer");
            require_positive(shaft_ambient_area_ratio, "shaft_ambient_area_ratio");
            require_positive(thermal_ambient_temperature, "thermal_ambient_temperature");
            if (!finite(shaft_ambient_temperature)) {
                error("shaft_ambient_temperature must be finite; <= 0 uses thermal_ambient_temperature.");
            }
            if (!finite(bush_outer_temperature)) {
                error("bush_outer_temperature must be finite; <= 0 uses thermal_ambient_temperature.");
            }
            if (bush_outer_radius > 0.0 && bush_outer_radius <= R) {
                error("bush_outer_radius must exceed R.");
            }
            if (bush_thermal_model == BushThermalModel::CONDUCTION_3D) {
                if (cross_film_thermal_model != CrossFilmModel::FINITE_VOLUME)
                    error("CONDUCTION_3D housing requires FINITE_VOLUME film energy.");
                if (bush_axial_cells < 1 || bush_axial_cells > n_z_global || n_z_global % bush_axial_cells != 0)
                    error("bush_axial_cells must be a positive divisor of n_z_global.");
                require_nonnegative(bush_end_heat_transfer, "bush_end_heat_transfer");
                if (bush_end_natural_convection) {
                    if (density_model == DensityModel::CONSTANT) {
                        if (!finite(density_temperature_slope) || density_temperature_slope >= 0.0)
                            error("CONSTANT density with housing end natural convection requires a finite negative density_temperature_slope for buoyancy.");
                        else
                            warn("CONSTANT density: density_temperature_slope is used only for exterior housing-end buoyancy; liquid density remains rho.");
                    } else if (density_model != DensityModel::LINEAR_TEMPERATURE)
                        error("Housing end natural convection requires LINEAR_TEMPERATURE or CONSTANT density with an explicit buoyancy slope.");
                }
            }
            else if (bush_end_heat_transfer != 0 || bush_end_natural_convection)
                error("Housing end convection requires CONDUCTION_3D.");
            if (bush_thermal_model != BushThermalModel::LOCAL) {
                if (bush_radial_cells < 1 || bush_radial_cells > 512) {
                    error("bush_radial_cells must be between 1 and 512.");
                }
                if (bush_outer_radius <= 0.0) {
                    warn("bush_thermal_model = CONDUCTION_2D with no bush_outer_radius: "
                         "the annulus falls back to 1.1*R, which is almost certainly not "
                         "the bush you meant. Set bush_outer_radius explicitly.");
                }
                if (bush_outer_boundary == BushOuterBoundary::CONVECTION &&
                    bush_outer_heat_transfer <= 0.0) {
                    warn("bush_thermal_model = CONDUCTION_2D with bush_outer_heat_transfer = 0 "
                         "and a CONVECTION outer boundary: the bush is globally adiabatic and "
                         "only redistributes heat circumferentially.");
                }
            }
            if (bush_outer_boundary == BushOuterBoundary::FIXED_TEMPERATURE &&
                bush_outer_heat_transfer > 0.0) {
                warn("bush_outer_boundary = FIXED_TEMPERATURE: the outer surface is held at "
                     "bush_outer_temperature when positive, otherwise thermal_ambient_temperature; "
                     "bush_outer_heat_transfer is IGNORED.");
            }
        }
        if (bush_thermal_model != BushThermalModel::LOCAL &&
            cross_film_thermal_model == CrossFilmModel::NONE) {
            error("bush_thermal_model = CONDUCTION_2D requires "
                  "cross_film_thermal_model = POLYNOMIAL; the conjugate bush only exists "
                  "inside the cross-film closure.");
        }
        if (groove_mixing_model != GrooveMixingModel::PINNED) {
            if (temperature_model != TemperatureModel::ENERGY_EQUATION) {
                error("groove_mixing_model requires temperature_model = ENERGY_EQUATION.");
            }
            bool any_fed = false;
            for (const auto& inlet : inlets) any_fed = any_fed || inlet.feeds_fresh_oil;
            if (!any_fed) {
                error("groove_mixing_model is set but no inlet feeds fresh oil, so the key "
                      "would have no effect. Give an inlet a t_supply token or leave the "
                      "model at PINNED.");
            }
            if (groove_mixing_model == GrooveMixingModel::FLUX_BALANCED &&
                (cavitation_model != CavitationModel::FBNS ||
                 fluid_property_model != FluidPropertyModel::SINGLE_PHASE ||
                 dissolved_gas_initial != 0.0 || dissolved_gas_henry_coeff != 0.0 ||
                 solution_mode != SolutionMode::STEADY_STATE ||
                 thermal_convection_scheme != ConvectionScheme::UPWIND ||
                 liquid_specific_heat <= 0.0)) {
                error("FLUX_BALANCED requires steady gas-free SINGLE_PHASE FBNS, UPWIND "
                      "thermal transport, and positive liquid_specific_heat in J/(kg K).");
            }
        }
        if (motion_model == MotionModel::MOVING_BEARING) {
            require_positive(bearing_mass, "bearing_mass");
        }
        if (motion_model == MotionModel::STEADY_EQUILIBRIUM) {
            if (solution_mode != SolutionMode::STEADY_STATE || cavitation_model != CavitationModel::FBNS ||
                fluid_property_model != FluidPropertyModel::SINGLE_PHASE || dissolved_gas_initial != 0.0 ||
                dissolved_gas_henry_coeff != 0.0 || gas_pressure_coupling != GasPressureCoupling::NONE)
                error("STEADY_EQUILIBRIUM requires steady gas-free SINGLE_PHASE FBNS with gas_pressure_coupling = NONE.");
            require_positive(steady_force_tolerance, "steady_force_tolerance");
            require_positive(steady_position_step_limit, "steady_position_step_limit");
            require_positive(steady_position_fd_step, "steady_position_fd_step");
            if (!finite(external_load_x) || !finite(external_load_y) ||
                std::hypot(external_load_x, external_load_y) <= 0.0)
                error("STEADY_EQUILIBRIUM requires a nonzero finite in-plane external load.");
            if (bearing_initial_z != 0.0 || (!bearing_initial_from_attitude &&
                (!finite(bearing_initial_x) || !finite(bearing_initial_y))))
                error("STEADY_EQUILIBRIUM requires a finite in-plane initial bearing position and bearing_initial_z = 0.");
            const double initial_x = bearing_initial_from_attitude ? -e*std::cos(attitude_angle_deg*3.14159265358979323846/180.0) : bearing_initial_x;
            const double initial_y = bearing_initial_from_attitude ? -e*std::sin(attitude_angle_deg*3.14159265358979323846/180.0) : bearing_initial_y;
            const double end_offset = 0.5*L;
            const double max_offset = std::max(std::hypot(initial_x-tilt_slope_x*end_offset,initial_y-tilt_slope_y*end_offset),
                std::hypot(initial_x+tilt_slope_x*end_offset,initial_y+tilt_slope_y*end_offset));
            if (!finite(max_offset) || c-max_offset <= min_film_thickness)
                error("STEADY_EQUILIBRIUM initial position violates continuous minimum film thickness.");
            if (external_load_z != 0.0 || bearing_stiffness_x != 0.0 || bearing_stiffness_y != 0.0)
                error("STEADY_EQUILIBRIUM supports two in-plane load components without bearing springs.");
        }

        require_pressure(bc_z_south_type, bc_z_south_val, "bc_z_south_val");
        require_pressure(bc_z_north_type, bc_z_north_val, "bc_z_north_val");
        for (bool south : {true, false}) {
            const std::string prefix = south ? "bc_z_south_" : "bc_z_north_";
            const double temperature = south ? bc_z_south_temperature : bc_z_north_temperature;
            const double concentration = south ? bc_z_south_dissolved_gas : bc_z_north_dissolved_gas;
            const double quality = south ? bc_z_south_free_gas_mass_fraction : bc_z_north_free_gas_mass_fraction;
            const auto thermal = south ? bc_z_south_thermal : bc_z_north_thermal;
            if (!finite(temperature) || temperature < 0)
                error(prefix + "temperature must be positive kelvin, or 0 to inherit temperature_reference.");
            if (temperature > 0 && thermal != ThermalInflowMode::CONSTANT)
                error(prefix + "temperature requires " + prefix + "thermal = CONSTANT.");
            if (!finite(concentration) || (concentration != -1.0 &&
                (concentration < 0 || concentration > dissolved_gas_max || concentration >= 1)))
                error(prefix + "dissolved_gas must be -1 (inherit), or a liquid mass fraction in [0,1).");
            if (!finite(quality) || quality < 0 || quality >= 1)
                error(prefix + "free_gas_mass_fraction must lie in [0,1).");
        }
        for (size_t idx = 0; idx < inlets.size(); ++idx) {
            const InletConfig& inlet = inlets[idx];
            const std::string prefix = "inlet[" + std::to_string(idx) + "]";
            if (!finite(inlet.theta)) error(prefix + ".theta must be a finite angle in degrees.");
            require_positive(inlet.size, (prefix + ".size").c_str());
            if (inlet.type == InletConfig::Type::CIRCULAR) {
                require_nonnegative(inlet.z, (prefix + ".z").c_str());
            }
            if (!finite(inlet.p_supply) || inlet.p_supply <= 0.0) {
                error(prefix + ".p_supply must be a positive absolute pressure in Pa.");
            } else if (finite(p_cav) && p_cav > 0.0 && inlet.p_supply < effective_p_cav()) {
                error(prefix + ".p_supply is below the effective cavitation pressure; "
                      "the pressure-floor model cannot impose this inlet reservoir.");
            }
            // t_supply defaults to 0 K and feeds_fresh_oil is set by the mere presence
            // of the optional token, so a stray or mistyped value would silently pin the
            // groove to absolute zero.
            if (inlet.feeds_fresh_oil && (!finite(inlet.t_supply) || inlet.t_supply <= 0.0)) {
                error(prefix + ".t_supply must be a positive absolute temperature in K "
                               "when the inlet feeds fresh oil.");
            }
            if (inlet.feeds_fresh_oil &&
                groove_mixing_model == GrooveMixingModel::FIXED_LAMBDA &&
                (!finite(inlet.mixing_lambda) || inlet.mixing_lambda < 0.0 ||
                 inlet.mixing_lambda > 1.0)) {
                error(prefix + ".mixing_lambda must be in [0, 1] under "
                               "groove_mixing_model = FIXED_LAMBDA.");
            }
            if (inlet.feeds_fresh_oil &&
                groove_mixing_model == GrooveMixingModel::FLOW_WEIGHTED &&
                (!finite(inlet.supply_flow_rate) || inlet.supply_flow_rate <= 0.0)) {
                error(prefix + ".supply_flow_rate must be a positive volumetric flow rate "
                               "[m^3/s] under groove_mixing_model = FLOW_WEIGHTED; without it "
                               "the mixing fraction is undefined.");
            }
        }

        if (cavitation_model == CavitationModel::GUMBEL) {
            warn("GUMBEL is a pressure clamp and is not mass-conserving JFO cavitation.");
        }

        if (output_field_enabled("mass_residual") &&
            cavitation_model != CavitationModel::JFO &&
            cavitation_model != CavitationModel::FBNS && cavitation_model != CavitationModel::COUPLED) {
            warn("output field mass_residual is the residual of the mass-conserving "
                 "cavitation balance (JFO or FBNS); under the selected cavitation model "
                 "it is written as zero.");
        }

        // The film-surface outputs are evaluations of the cross-film profile.
        // Without that model there is no profile to evaluate, and writing zeros
        // (or the depth-averaged value twice) would look like a physical result.
        if (cross_film_thermal_model == CrossFilmModel::NONE) {
            for (const char* name : OutputNaming::surface_fields()) {
                if (!output_field_enabled(name)) continue;
                error(std::string("output field ") + name + " is an evaluation of the "
                      "cross-film temperature profile at a film surface and only exists "
                      "under cross_film_thermal_model = POLYNOMIAL. Enable the cross-film "
                      "model, or remove the field from output_fields (the depth-averaged "
                      "field `temperature` is unaffected either way).");
            }
        }

        if (cavitation_model == CavitationModel::FBNS) {
            // Phase 3 scope. Temperature, the liquid solution properties and the
            // free-gas phase are in: the first two enter the FBNS operator only
            // through its per-cell coefficients (mu, rho_l), and the third through
            // the lagged full-film ceiling theta <= 1 - alpha_gas, which is affine
            // in theta and leaves both Jacobian blocks unchanged. What is still
            // refused makes the complementarity datum or the operator's stencil
            // something FBNS has no way to advance.
            if (gas_pressure_coupling != GasPressureCoupling::NONE) {
                error("cavitation_model = FBNS requires gas_pressure_coupling = NONE. "
                      "VOID_COUPLED turns the cavitation plateau into a per-cell field "
                      "(\"p_void\") that only Reynolds::solve_elrod's outer loop advances "
                      "and under-relaxes; FBNS has no such loop, so it would run the "
                      "plateau frozen at its initial value. Use cavitation_model = JFO for "
                      "the void-coupled plateau closure.");
            }
            // cross_film_thermal_model is supported: the FBNS face mobilities read
            // the same generalized-Dowson F-integrals the Elrod path uses
            // (rho_l F2 and the Couette thickness 2 F1/F0, ReynoldsCoeff), and
            // CrossFilm::update rebuilds them once per outer iteration, so inside
            // one Newton solve they are frozen coefficients exactly as mu is.
            if (bc_z_south_type == BCType::INLET_OUTLET || bc_z_north_type == BCType::INLET_OUTLET) {
                error("cavitation_model = FBNS supports DIRICHLET and NEUMANN axial boundaries; "
                      "INLET_OUTLET switches on the sign of the solution and would break the "
                      "differentiability the Newton solve relies on.");
            }
            if (theta_convection_scheme != ConvectionScheme::UPWIND) {
                warn("cavitation_model = FBNS transports the film content with first-order "
                     "upwind (the FBNS discretization); theta_convection_scheme is ignored.");
            }
            if (fbns_max_iters <= 0) error("fbns_max_iters must be positive.");
            if (!finite(fbns_tol) || fbns_tol <= 0.0) error("fbns_tol must be positive.");
            if (!finite(fbns_relaxation) || fbns_relaxation <= 0.0 || fbns_relaxation > 1.0) {
                error("fbns_relaxation must lie in (0, 1].");
            }
            if (!finite(fbns_pressure_scale) || fbns_pressure_scale < 0.0) {
                error("fbns_pressure_scale must be >= 0 (0 = derive it from the operating point).");
            }
        }

        if (cavitation_model == CavitationModel::COUPLED) {
            require_nonnegative(oil_vapor_temperature_coeff,"oil_vapor_temperature_coeff");
            if (cavitation_threshold != CavitationThreshold::SCALAR_PCAV || gas_pressure_coupling != GasPressureCoupling::NONE)
                error("COUPLED requires SCALAR_PCAV and gas_pressure_coupling = NONE; pressure follows the gas EOS and gap volume.");
            if (cross_film_thermal_model != CrossFilmModel::NONE)
                error("COUPLED uses a common film temperature with prescribed wall heat transfer; cross_film_thermal_model must be NONE.");
            if (mixture_viscosity_model != MixtureViscosityModel::DUKLER)
                error("COUPLED currently requires DUKLER volume-weighted mixture viscosity.");
            if (temperature_model == TemperatureModel::ENERGY_EQUATION) {
                if (thermal_convection_scheme != ConvectionScheme::UPWIND)
                    error("COUPLED energy transport requires UPWIND.");
                require_nonnegative(gas_thermal_conductivity,"gas_thermal_conductivity");
                if (solution_density_gas_coeff != 0)
                    error("COUPLED energy requires additive, incompressible partial liquid volumes: solution_density_gas_coeff = 0.");
                if (fluid_property_model != FluidPropertyModel::CONSTANT) {
                    require_positive(dissolved_heat_capacity,"dissolved_heat_capacity");
                    if (!finite(desorption_internal_energy))
                        error("COUPLED energy requires measured or explicitly assumed desorption_internal_energy; no oil-specific heat of solution is supplied.");
                    if (dissolved_gas_species != DissolvedGasSpecies::R410A)
                        require_positive(gas_heat_capacity_cv,"gas_heat_capacity_cv");
                }
                if (groove_mixing_model != GrooveMixingModel::PINNED)
                    error("COUPLED conserves inlet enthalpy using the resolved supply mass; empirical groove mixing factors are not supported.");
            }
            if (solution_mode == SolutionMode::STEADY_STATE && inlets.empty() &&
                bc_z_south_type == BCType::NEUMANN && bc_z_north_type == BCType::NEUMANN)
                error("A closed COUPLED domain requires transient inventory evolution; steady pressure/composition has no reservoir datum.");
            if (bc_z_south_type == BCType::INLET_OUTLET || bc_z_north_type == BCType::INLET_OUTLET)
                error("COUPLED uses DIRICHLET for a pressure reservoir with reverse flow, or NEUMANN for a sealed end.");
            if (theta_convection_scheme != ConvectionScheme::UPWIND || gas_convection_scheme != ConvectionScheme::UPWIND)
                error("COUPLED currently requires UPWIND species and film convection.");
            if (density_model == DensityModel::TABLE || liquid_viscosity_model == LiquidViscosityModel::TABLE)
                error("COUPLED cannot use saturated p,T density/viscosity tables for off-equilibrium composition; use a concentration-dependent constitutive law.");
            for (const auto& inlet : inlets) {
                if (!finite(inlet.dissolved_gas) || (inlet.dissolved_gas < 0 && inlet.dissolved_gas != -1) || inlet.dissolved_gas >= 1)
                    error("inlet_species dissolved fraction must be -1 (inherit) or in [0,1).");
                if (!finite(inlet.free_gas_mass_fraction) || inlet.free_gas_mass_fraction < 0 || inlet.free_gas_mass_fraction >= 1)
                    error("inlet_species free mass fraction must be in [0,1).");
            }
        }
        if (dissolved_gas_species == DissolvedGasSpecies::R410A && cavitation_model != CavitationModel::COUPLED)
            error("R410A requires COUPLED; the legacy segregated gas closures assume an ideal gas.");
        if (oil_vapor_temperature_coeff != 0 && cavitation_model != CavitationModel::COUPLED)
            error("Temperature-dependent oil vapor pressure currently requires COUPLED.");

        // Outer (Picard) coupling schedule.
        if (!steady_initial_state_file.empty() &&
            (solution_mode != SolutionMode::STEADY_STATE || cavitation_model != CavitationModel::FBNS ||
             fluid_property_model != FluidPropertyModel::SINGLE_PHASE ||
             gas_pressure_coupling != GasPressureCoupling::NONE ||
             cavitation_threshold != CavitationThreshold::SCALAR_PCAV || dissolved_gas_initial != 0.0)) {
            error("steady_initial_state_file requires steady gas-free single-phase FBNS with scalar p_cav.");
        }
        if (coupling_max_iters <= 0) error("coupling_max_iters must be positive.");
        if (!finite(coupling_tolerance) || coupling_tolerance <= 0.0) {
            error("coupling_tolerance must be positive.");
        }
        if (steady_output_interval < 0) error("steady_output_interval must be nonnegative.");
        if (steady_pressure_inner_only &&
            (solution_mode != SolutionMode::STEADY_STATE || cavitation_model != CavitationModel::FBNS)) {
            error("steady_pressure_inner_only requires steady FBNS.");
        }
        if (!finite(coupling_relaxation) || coupling_relaxation <= 0.0 || coupling_relaxation > 1.0) {
            error("coupling_relaxation must lie in (0, 1].");
        }
        if (coupling_relaxation_mode == CouplingRelaxationMode::ADAPTIVE) {
            if (solution_mode != SolutionMode::STEADY_STATE) {
                warn("coupling_relaxation_mode = ADAPTIVE only acts in STEADY_STATE; "
                     "TRANSIENT runs one sweep per step and never relaxes.");
            }
            if (!finite(coupling_relaxation_min) || coupling_relaxation_min <= 0.0 ||
                coupling_relaxation_min > coupling_relaxation) {
                error("coupling_relaxation_min must lie in (0, coupling_relaxation].");
            }
            if (!finite(coupling_relaxation_backoff) || coupling_relaxation_backoff <= 0.0 ||
                coupling_relaxation_backoff >= 1.0) {
                error("coupling_relaxation_backoff must lie in (0, 1).");
            }
            if (!finite(coupling_relaxation_growth) || coupling_relaxation_growth < 1.0) {
                error("coupling_relaxation_growth must be >= 1.");
            }
        }

        if (thermal_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING) {
            error("TYPE_DIFFERENCING is only available for theta_convection_scheme (it switches on the cavitation flag).");
        }
        if (gas_convection_scheme == ConvectionScheme::TYPE_DIFFERENCING) {
            error("TYPE_DIFFERENCING is only available for theta_convection_scheme (it switches on the cavitation flag).");
        }

        // Laminar-regime guard: Taylor-vortex onset for a rotating journal,
        // Ta = Re_c sqrt(c/R) with Re_c = rho omega R c / mu (Hamrock; San Andres).
        if (finite(rho) && finite(omega) && finite(R) && finite(c) && finite(mu) &&
            mu > 0.0 && c > 0.0 && R > 0.0) {
            const double couette_reynolds = rho * std::abs(omega) * R * c / mu;
            const double taylor = couette_reynolds * std::sqrt(c / R);
            if (taylor > 41.0) {
                std::ostringstream msg;
                msg << "Couette Reynolds number " << couette_reynolds
                    << " gives Taylor number " << taylor
                    << " > 41: flow is beyond the laminar Taylor-vortex limit; "
                       "the laminar Reynolds equation is not valid at this operating point.";
                warn(msg.str());
            }
        }

        const bool variable_properties = fluid_property_model != FluidPropertyModel::CONSTANT;
        const bool gas_cavitation = fluid_property_model == FluidPropertyModel::TWO_PHASE;

        // Does a dissolved gas species actually take part in this run?
        //
        // The solubility law and the log-mixing viscosity exponent are consulted
        // only when free gas is released or resorbed (TWO_PHASE), when the
        // cavitation datum is the local bubble point (LOCAL_PSAT), or when the oil
        // carries a nonzero dissolved fraction that thins and lightens it. A
        // SINGLE_PHASE run of pure oil -- the natural configuration for a THD or a
        // piezo-viscous study -- consults none of them, so zero is the correct
        // value for both coefficients there and must be legal.
        const bool gas_species_participates =
            gas_cavitation ||
            cavitation_threshold == CavitationThreshold::LOCAL_PSAT ||
            (finite(dissolved_gas_initial) && dissolved_gas_initial > 0.0);

        if (variable_properties) {
            if (!finite(dissolved_gas_initial) || !finite(dissolved_gas_max) ||
                dissolved_gas_initial < 0.0 || dissolved_gas_max <= 0.0 || dissolved_gas_max > 1.0 ||
                dissolved_gas_initial > dissolved_gas_max) {
                error("Dissolved gas values must satisfy 0 <= dissolved_gas_initial <= dissolved_gas_max <= 1 and dissolved_gas_max > 0.");
            }
            require_nonnegative(dissolved_gas_henry_temp_coeff, "dissolved_gas_henry_temp_coeff");
            require_nonnegative(dissolved_gas_diffusivity, "dissolved_gas_diffusivity");
            require_positive(gas_pressure_floor, "gas_pressure_floor");
            require_positive(property_reference_pressure, "property_reference_pressure");
            require_positive(property_reference_temperature, "property_reference_temperature");

            switch (solubility_model) {
                case SolubilityModel::HENRY:
                    if (!finite(dissolved_gas_henry_coeff) || dissolved_gas_henry_coeff < 0.0) {
                        error("dissolved_gas_henry_coeff must be finite and >= 0.");
                    } else if (dissolved_gas_henry_coeff == 0.0 && gas_species_participates) {
                        error("HENRY solubility requires dissolved_gas_henry_coeff > 0 when a "
                              "dissolved gas species participates (TWO_PHASE, LOCAL_PSAT, or "
                              "dissolved_gas_initial > 0); with H = 0 the saturation limit is "
                              "identically zero and all the gas would flash at once.");
                    }
                    break;
                case SolubilityModel::TABLE:
                    if (solubility_table_2d.empty()) require_table(solubility_table, "solubility_table");
                    break;
                case SolubilityModel::CURVE_TABLE:
                    if (solubility_curve_table.empty()) error("CURVE_TABLE requires solution_table_file.");
                    break;
            }

            switch (density_model) {
                case DensityModel::LINEAR_TEMPERATURE:
                    if (!finite(density_temperature_intercept) || !finite(density_temperature_slope))
                        error("Linear density coefficients must be finite.");
                    for (double t : {temperature_initial, property_reference_temperature}) {
                        if (!(density_temperature_intercept + density_temperature_slope * (t - 273.15) > 0))
                            error("LINEAR_TEMPERATURE density must be positive at initial/reference temperatures.");
                    }
                    break;
                case DensityModel::CONSTANT:
                    break;
                case DensityModel::MIXTURE:
                    // Both phase densities are explicit inputs: rho (oil liquid) and
                    // dissolved_gas_liquid_density (the dissolved-gas liquid phase).
                    if (!finite(dissolved_gas_liquid_density) || dissolved_gas_liquid_density <= 0.0) {
                        error("MIXTURE density requires dissolved_gas_liquid_density > 0 (the dissolved-gas liquid-phase density).");
                    }
                    break;
                case DensityModel::TABLE:
                    if (density_table_2d.empty()) require_table(density_table, "density_table");
                    break;
            }

            switch (liquid_viscosity_model) {
                case LiquidViscosityModel::CONSTANT:
                    break;
                case LiquidViscosityModel::EMPIRICAL:
                    // mu = mu_oil exp[E_mu(1/T - 1/T_ref)] exp[a_c c_d] exp[alpha_p(p - p_ref)].
                    // Any subset of the three exponents may legitimately be zero: an
                    // Andrade-only thermal run carries no dissolved gas, so a_c = 0 is
                    // the correct input, not a missing one. Only all three at zero is
                    // certainly a mistake, and only then is the model degenerate.
                    if (!finite(solution_viscosity_gas_coeff)) {
                        error("solution_viscosity_gas_coeff must be finite.");
                    } else if (solution_viscosity_gas_coeff == 0.0 && gas_species_participates) {
                        warn("EMPIRICAL viscosity with solution_viscosity_gas_coeff = 0 ignores "
                             "the dissolved gas that is present: the oil will not thin as gas "
                             "goes into solution.");
                    }
                    if (viscosity_temperature_coeff == 0.0 && viscosity_pressure_coeff == 0.0 &&
                        solution_viscosity_gas_coeff == 0.0) {
                        warn("EMPIRICAL viscosity with viscosity_temperature_coeff, "
                             "viscosity_pressure_coeff and solution_viscosity_gas_coeff all zero "
                             "is identical to liquid_viscosity_model = CONSTANT.");
                    }
                    require_nonnegative(viscosity_temperature_coeff, "viscosity_temperature_coeff");
                    require_nonnegative(viscosity_pressure_coeff, "viscosity_pressure_coeff");
                    break;
                case LiquidViscosityModel::POLYNOMIAL: {
                    // mu(T) = mu (k0 - k1 Tbar + k2 Tbar^2), Tbar = T[degC]/T_ref[degC].
                    if (!finite(viscosity_poly_k0) || !finite(viscosity_poly_k1) ||
                        !finite(viscosity_poly_k2)) {
                        error("viscosity_poly_k0, viscosity_poly_k1 and viscosity_poly_k2 must be finite.");
                    }
                    // Tbar is a Celsius ratio, so a reference at or below 0 degC has no
                    // usable normalisation (and 0 degC divides by zero).
                    const double T_ref_C = polynomial_reference_temperature_celsius();
                    if (!(T_ref_C > 0.0)) {
                        error("liquid_viscosity_model = POLYNOMIAL needs property_reference_temperature "
                              "above 273.15 K: the Ferron law normalises by the reference temperature "
                              "in DEGREES CELSIUS, which must be positive.");
                    }
                    // mu_ref is defined AT the reference temperature, so mubar(T_ref) =
                    // k0 - k1 + k2 must be 1. Ferron's own set gives exactly 1.000.
                    // This is a warning, not an error: a user may deliberately fit a set
                    // whose reference viscosity is carried by mu instead.
                    const double mu_bar_ref = viscosity_poly_k0 - viscosity_poly_k1 + viscosity_poly_k2;
                    if (finite(mu_bar_ref) && std::abs(mu_bar_ref - 1.0) > 0.05) {
                        warn("POLYNOMIAL viscosity: k0 - k1 + k2 = " + std::to_string(mu_bar_ref) +
                             " but mu is defined AT property_reference_temperature, where the "
                             "polynomial should return 1. The run will use mu * " +
                             std::to_string(mu_bar_ref) + " at the reference temperature.");
                    }
                    if (viscosity_poly_k1 == 0.0 && viscosity_poly_k2 == 0.0) {
                        warn("POLYNOMIAL viscosity with viscosity_poly_k1 = viscosity_poly_k2 = 0 "
                             "has no temperature dependence and is identical to "
                             "liquid_viscosity_model = CONSTANT.");
                    }
                    break;
                }
                case LiquidViscosityModel::TABLE:
                    if (viscosity_table_2d.empty()) require_table(viscosity_table, "viscosity_table");
                    break;
                case LiquidViscosityModel::COMPOSITION_TABLE:
                    if (viscosity_composition_table.empty()) error("COMPOSITION_TABLE requires viscosity_composition_file (c_mass_fraction,T_C,mu_Pa_s).");
                    if (viscosity_pressure_coeff != 0) error("COMPOSITION_TABLE has no pressure axis; supply a justified pressure-independent viscosity table, with viscosity_pressure_coeff = 0.");
                    break;
                case LiquidViscosityModel::KINEMATIC_TABLE:
                    if (kinematic_composition_table.empty()) error("KINEMATIC_TABLE requires solution_table_file.");
                    if (viscosity_pressure_coeff != 0) error("KINEMATIC_TABLE uses actual c,T and requires viscosity_pressure_coeff = 0.");
                    break;
            }
        }

        // 2-D table-file precedence/usage warnings (a loaded 2-D file takes
        // precedence over an inline 1-D table, and is only consumed under TABLE).
        auto warn_table_use = [&](const PropertyTable2D& t2d, bool model_is_table,
                                  const std::vector<PropertyTablePoint>& t1d, const char* file_key) {
            if (t2d.empty()) return;
            if (!model_is_table)
                warn(std::string(file_key) + " is loaded but the corresponding model is not TABLE; the 2-D table is unused.");
            else if (!t1d.empty())
                warn(std::string("both an inline 1-D table and ") + file_key +
                     " are set; the 2-D file takes precedence.");
        };
        warn_table_use(solubility_table_2d, solubility_model == SolubilityModel::TABLE,
                       solubility_table, "solubility_table_file");
        warn_table_use(viscosity_table_2d, liquid_viscosity_model == LiquidViscosityModel::TABLE,
                       viscosity_table, "viscosity_table_file");
        warn_table_use(density_table_2d, density_model == DensityModel::TABLE,
                       density_table, "density_table_file");

        if (cavitated_film_model == CavitatedFilmModel::STRIATED &&
            cavitation_model != CavitationModel::JFO &&
            cavitation_model != CavitationModel::FBNS) {
            error("cavitated_film_model = STRIATED requires cavitation_model = JFO or FBNS "
                  "(it weights the cavitated zone by the mass-conserving film content theta).");
        }

        if (gas_pressure_coupling == GasPressureCoupling::VOID_COUPLED) {
            if (fluid_property_model != FluidPropertyModel::TWO_PHASE) {
                warn("gas_pressure_coupling = VOID_COUPLED has no effect without "
                     "fluid_property_model = TWO_PHASE (no free-gas field).");
            }
            if (!finite(mu_gas) || mu_gas <= 0.0) {
                warn("gas_pressure_coupling = VOID_COUPLED uses the free-gas density; set mu_gas > 0 "
                     "and a consistent gas EOS for the released phase.");
            }
        }

        if (gas_cavitation) {
            if (!mass_conserving_cavitation()) {
                error("TWO_PHASE requires a mass-conserving cavitation model (JFO or FBNS) "
                      "so gas release is coupled to the liquid-content solve.");
            }
            if (finite(dissolved_gas_initial) && dissolved_gas_initial <= 0.0) {
                warn("TWO_PHASE is selected with zero initial dissolved gas; no gaseous cavitation will occur until gas enters the film.");
            }
            if (!finite(gas_mass_transfer_rate) || gas_mass_transfer_rate < 0.0) {
                error("TWO_PHASE requires gas_mass_transfer_rate >= 0 (0 disables release).");
            }
            if (!finite(gas_absorption_rate) || (gas_absorption_rate < 0 && gas_absorption_rate != -1))
                error("gas_absorption_rate must be -1 (inherit release rate) or >= 0 [1/s].");
            if (solution_mode == SolutionMode::STEADY_STATE &&
                steady_gas_model == SteadyGasModel::FROZEN && !allow_frozen_gas_steady && cavitation_model != CavitationModel::COUPLED) {
                error("TWO_PHASE in STEADY_STATE with steady_gas_model = FROZEN "
                      "silently freezes the gas model. Set steady_gas_model = EQUILIBRIUM, or "
                      "allow_frozen_gas_steady = true to accept the frozen behaviour.");
            }
            if (!finite(gas_alpha_max) || gas_alpha_max <= 0.0 || gas_alpha_max > 1.0) {
                error("gas_alpha_max must be in (0, 1] for TWO_PHASE.");
            }
            if (finite(dt) && finite(gas_mass_transfer_rate) && dt * gas_mass_transfer_rate > 1.0) {
                warn(cavitation_model==CavitationModel::COUPLED ?
                     "dt * gas_mass_transfer_rate > 1; backward Euler is implicit, but exchange accuracy still requires a timestep study." :
                     "dt * gas_mass_transfer_rate > 1; local exchange is integrated analytically, but transport/pressure splitting still requires a timestep study.");
            }

            switch (mixture_viscosity_model) {
                case MixtureViscosityModel::DUKLER:
                case MixtureViscosityModel::LINEAR:
                case MixtureViscosityModel::MCADAMS:
                    if (!finite(mu_gas) || mu_gas <= 0.0) {
                        error("The selected mixture_viscosity_model requires mu_gas > 0 (free-gas dynamic viscosity in Pa s).");
                    }
                    break;
            }

            // The charge bubble point is not the current reservoir saturation state.
            const double p_sat = saturation_pressure(dissolved_gas_initial, temperature_initial);
            if (finite(p_sat) && p_sat > 0.0) {
                std::ostringstream context;
                context << "Gaseous-cavitation context: p_sat(T_init, c_d_init) = " << p_sat
                        << " Pa, p_cav / p_sat = " << (p_cav / p_sat) << ".";
                warn(context.str());
                auto saturation_ratio_check = [&](BCType type, double value, const char* label) {
                    if (type == BCType::NEUMANN || !finite(value) || value <= 0.0) return;
                    const double ratio = value / p_sat;
                    if (ratio >= 1.0) {
                        std::ostringstream msg;
                        msg << label << " / p_sat = " << ratio
                            << " >= 1: the initial charge is at or below equilibrium solubility at this pressure. "
                               "Release requires crossing its temperature-dependent bubble point; reservoir composition may differ.";
                        warn(msg.str());
                    }
                };
                saturation_ratio_check(bc_z_south_type, bc_z_south_val, "bc_z_south_val");
                saturation_ratio_check(bc_z_north_type, bc_z_north_val, "bc_z_north_val");
                for (size_t idx = 0; idx < inlets.size(); ++idx) {
                    const std::string label = "inlet[" + std::to_string(idx) + "].p_supply";
                    saturation_ratio_check(BCType::DIRICHLET, inlets[idx].p_supply, label.c_str());
                }
            }
        }
    }

    /// Ostwald solubility coefficient L = V_gas(dissolved) / V_liquid, dimensionless.
    ///
    /// Derived from the configured Henry slope: with c_sat = H(T) p (mass fraction)
    /// and an ideal free-gas phase rho_g = p / (R_s T),
    ///     L = c_sat rho_liquid / rho_g = H rho_liquid R_s T,
    /// which is independent of pressure, as the Ostwald coefficient must be.
    /// Used by the VOID_COUPLED plateau closure; returns 0 when no Henry model
    /// applies, which makes callers fall back to the uncoupled threshold.
    double ostwald_coefficient(double temperature, double rho_liquid) const {
        if (fluid_property_model == FluidPropertyModel::CONSTANT) return 0.0;
        if (solubility_model != SolubilityModel::HENRY) return 0.0;
        const double T = temperature > 0.0 ? temperature : property_reference_temperature;
        double henry = dissolved_gas_henry_coeff;
        if (dissolved_gas_henry_temp_coeff != 0.0 && T > 0.0 &&
            property_reference_temperature > 0.0) {
            henry *= std::exp(dissolved_gas_henry_temp_coeff *
                              (1.0 / T - 1.0 / property_reference_temperature));
        }
        if (!(henry > 0.0)) return 0.0;
        const double rho_l = rho_liquid > 0.0 ? rho_liquid : rho;
        return henry * rho_l * gas_specific_constant_for_species() * T;
    }

    /// Reference temperature of the POLYNOMIAL viscosity law in DEGREES CELSIUS.
    ///
    /// Ferron, Frene & Boncompain (1983) nondimensionalise temperature as
    /// Tbar = T/T_0 with both in degC (their nomenclature: "T = temperature degC",
    /// "Tbar = T/T_0 nondimensional temperature"). Everything else in this config
    /// is Kelvin, so the single conversion point lives here and both the validator
    /// and FluidProperties call it rather than open-coding 273.15.
    /// Returns a non-positive value for references at or below 0 degC; validate()
    /// rejects those, and the viscosity law falls back to mu.
    double polynomial_reference_temperature_celsius() const {
        return property_reference_temperature - 273.15;
    }

    /// Specific gas constant R_s of the dissolved species [J/(kg K)].
    /// Mirrors FluidProperties::gas_specific_constant; duplicated here so the
    /// config header stays free of a dependency on the fluid-properties unit.
    double gas_specific_constant_for_species() const {
        switch (dissolved_gas_species) {
            case DissolvedGasSpecies::AIR: return 287.05;
            case DissolvedGasSpecies::PROPANE: return 188.55;
            case DissolvedGasSpecies::R410A: return 8.314472 / 0.0725854;
        }
        return 188.55;
    }

    /// Pressure at which oil with dissolved mass fraction c_d is saturated:
    /// inversion of the configured solubility law c_eq(p, T) = c_d.
    /// Returns 0 when no solubility model applies or c_d <= 0.
    double saturation_pressure(double c_d, double temperature) const {
        if (fluid_property_model == FluidPropertyModel::CONSTANT || c_d <= 0.0) return 0.0;
        const double T = temperature > 0.0 ? temperature : property_reference_temperature;
        switch (solubility_model) {
            case SolubilityModel::HENRY: {
                double henry = dissolved_gas_henry_coeff;
                if (dissolved_gas_henry_temp_coeff != 0.0 && T > 0.0 &&
                    property_reference_temperature > 0.0) {
                    henry *= std::exp(dissolved_gas_henry_temp_coeff *
                                      (1.0 / T - 1.0 / property_reference_temperature));
                }
                return henry > 0.0 ? c_d / henry : 0.0;
            }
            case SolubilityModel::CURVE_TABLE: {
                if (solubility_curve_table.empty()) throw std::domain_error("Missing solution curve table");
                std::vector<double> concentrations;
                for (std::size_t i = 0; i < solubility_curve_table.axis.size(); ++i)
                    concentrations.push_back(solubility_curve_table.isotherm_at(i, T, true));
                const double p = invert_csat(solubility_curve_table.axis, concentrations, c_d);
                if (!solution_table_clamp && std::abs(
                    solubility_curve_table.interpolate(p, T, false) - c_d) > 1e-10)
                    throw std::domain_error("Bubble pressure outside digitized solution table range");
                return p;
            }
            case SolubilityModel::TABLE: {
                // Invert c_sat = c_d for the smallest pressure (the bubble point).
                if (!solubility_table_2d.empty()) {
                    // 2-D surface: build the c_sat(p) curve along the local isotherm, then invert.
                    std::vector<double> c_iso(solubility_table_2d.p.size());
                    for (std::size_t ip = 0; ip < c_iso.size(); ++ip) {
                        c_iso[ip] = solubility_table_2d.isotherm_at(ip, T);
                    }
                    return invert_csat(solubility_table_2d.p, c_iso, c_d);
                }
                std::vector<double> pv, cv;
                pv.reserve(solubility_table.size());
                cv.reserve(solubility_table.size());
                for (const auto& pt : solubility_table) {
                    pv.push_back(pt.x);
                    cv.push_back(pt.value);
                }
                return invert_csat(pv, cv, c_d);
            }
        }
        return 0.0;
    }

    bool output_field_enabled(const std::string& name) const {
        const std::string key = normalise_key(name);
        if (std::find(output_fields.begin(), output_fields.end(), key) != output_fields.end()) return true;
        if (key == "pressure_force_x") {
            return std::find(output_fields.begin(), output_fields.end(), "load_x") != output_fields.end();
        }
        if (key == "pressure_force_y") {
            return std::find(output_fields.begin(), output_fields.end(), "load_y") != output_fields.end();
        }
        if (key == "pressure_force_z") {
            return std::find(output_fields.begin(), output_fields.end(), "load_z") != output_fields.end();
        }
        return false;
    }

    /// Effective cavitation/onset pressure. Under LOCAL_PSAT it is the larger of
    /// the oil cavitation pressure and the dissolved-species bubble point at the
    /// feed state, so whichever species releases first governs onset. Uniform, so
    /// it acts as a pressure datum without injecting a net force. Reduces to p_cav
    /// under SCALAR_PCAV or when no dissolved gas / solubility model is present.
    double effective_p_cav() const {
        if (cavitation_threshold != CavitationThreshold::LOCAL_PSAT) return p_cav;
        const double p_sat = saturation_pressure(dissolved_gas_initial, temperature_initial);
        return (std::isfinite(p_sat) && p_sat > p_cav) ? p_sat : p_cav;
    }

    /// True when the selected model resolves the JFO free boundary and therefore
    /// maintains a physical film content theta in [theta_min, 1].
    bool mass_conserving_cavitation() const {
        return cavitation_model == CavitationModel::JFO ||
               cavitation_model == CavitationModel::FBNS || cavitation_model == CavitationModel::COUPLED;
    }

    double initial_pressure() const {
        if (mass_conserving_cavitation()) {
            if (bc_z_south_type != BCType::NEUMANN) return elrod_boundary_pressure(bc_z_south_val);
            if (bc_z_north_type != BCType::NEUMANN) return elrod_boundary_pressure(bc_z_north_val);
            return effective_p_cav();
        }
        if (bc_z_south_type != BCType::NEUMANN) return bc_z_south_val;
        if (bc_z_north_type != BCType::NEUMANN) return bc_z_north_val;
        return effective_p_cav();
    }

    double elrod_boundary_pressure(double pressure_value) const {
        return std::max(pressure_value, effective_p_cav());
    }

    double elrod_boundary_theta(double pressure_value) const {
        const double beta = std::max(bulk_modulus, 1.0e-30);
        const double pressure = elrod_boundary_pressure(pressure_value);
        return std::max(std::exp((pressure - effective_p_cav()) / beta), theta_min);
    }

    /// Net liquid outflow through one axial Dirichlet/inflow boundary face of a
    /// boundary cell (positive = leaving the domain), in the Elrod theta-equation
    /// units. Full-film cells use the barotropic Gamma_base link; under
    /// consistent_boundary_flux a cavitated cell re-floods at the physical
    /// Poiseuille rate (an inflow, hence negative outflow). Single source of
    /// truth shared by solve_elrod's assembly and the diagnostics mass balance.
    ///
    /// `theta_full` is the cell's full-film ceiling, `FullFilm::bound(...)`. It is
    /// passed in rather than looked up because `FullFilm` includes this header;
    /// the ceiling still has exactly one definition, callers just carry the value
    /// across the include boundary. There is no default: a silent 1.0 here is the
    /// stale-ceiling bug the shared predicate exists to prevent.
    double elrod_boundary_outflow(double gamma_base_value, double theta_cell, double theta_full,
                                  double bc_val, double R, double d_theta, double d_z) const {
        const double cs = gamma_base_value * R * d_theta / (0.5 * d_z);
        if (consistent_boundary_flux && theta_cell < theta_full) {
            const double p_bc = elrod_boundary_pressure(bc_val);
            return -std::max(0.0, cs * (p_bc - effective_p_cav()) / std::max(bulk_modulus, 1.0e-30));
        }
        return cs * (theta_cell - std::max(elrod_boundary_theta(bc_val), theta_min));
    }

    /// True when a cavitated axial-boundary cell re-floods at the physical
    /// reformation rate (consistent_boundary_flux + Elrod) instead of the
    /// full-film Dirichlet link. The energy and dissolved-gas convection use this
    /// to switch their boundary face from the full-film Poiseuille gate to the
    /// reflood inflow, so the reformation liquid carries reservoir temperature /
    /// composition consistently with the Elrod theta-equation.
    ///
    /// `theta_full` is `FullFilm::bound(...)`; see elrod_boundary_outflow for why
    /// it is a parameter.
    bool is_reformation_boundary(double theta_cell, double theta_full) const {
        return consistent_boundary_flux &&
               cavitation_model == CavitationModel::JFO && theta_cell < theta_full;
    }

    /// +z Poiseuille velocity [m/s] of the reformation reflood at one open axial
    /// end, driven from the cavitation plateau (p_bc - p_cav) across the half-cell
    /// (half_dz = 0.5*d_z). is_south selects the end (south reflood is +z inflow,
    /// north reflood is -z inflow). The magnitude matches the inflow
    /// elrod_boundary_outflow imposes on the liquid balance, so the liquid, energy,
    /// and gas equations all see one reformation rate.
    double reformation_boundary_velocity(double bc_val, double h_face, double mu_face,
                                         double half_dz, bool is_south) const {
        const double drop = elrod_boundary_pressure(bc_val) - effective_p_cav();  // >= 0
        const double u_mag = h_face * h_face * drop /
            (12.0 * std::max(mu_face, 1.0e-30) * std::max(half_dz, 1.0e-30));
        return is_south ? u_mag : -u_mag;
    }

    double omega_at_time(double time) const {
        if (solution_mode == SolutionMode::STEADY_STATE || omega_ramp_time <= 0.0) return omega;
        const double fraction = std::clamp(time / omega_ramp_time, 0.0, 1.0);
        return omega * fraction;
    }

    std::string output_fields_text() const {
        std::ostringstream out;
        for (size_t i = 0; i < output_fields.size(); ++i) {
            if (i > 0) out << ", ";
            out << output_fields[i];
        }
        return out.str();
    }

    static std::string property_table_text(const std::vector<PropertyTablePoint>& table) {
        std::ostringstream out;
        out << std::setprecision(17);
        for (size_t i = 0; i < table.size(); ++i) {
            if (i > 0) out << ", ";
            out << table[i].x << ":" << table[i].value;
        }
        return out.str();
    }

    static std::vector<PropertyTablePoint> parse_property_table_text(std::string value) {
        return parse_property_table(std::move(value));
    }

    static std::vector<std::string> default_output_fields() {
        return {
            "pressure",
            "film_content",
            "h",
            "rho",
            "mu",
            "rho_liquid_solution",
            "mu_liquid_solution",
            "temperature",
            "heat_generation",
            "dissolved_gas",
            "alpha_gas",
            "gas_mass_transfer",
            "inlet_indicator",
            "velocity",
            "pressure_force_x",
            "pressure_force_y",
            "pressure_force_z",
            "viscous_force_x",
            "viscous_force_y",
            "viscous_force_z",
            "fluid_force_x",
            "fluid_force_y",
            "fluid_force_z",
            "external_load_x",
            "external_load_y",
            "external_load_z",
            "bearing_x",
            "bearing_y",
            "bearing_z",
            "bearing_attitude_angle",
            "friction_torque"};
    }

    std::string to_config_text() const {
        std::ostringstream out;
        out << std::setprecision(17);

        out << "# Geometry\n"
            << "R = " << R << "\n"
            << "c = " << c << "\n"
            << "e = " << e << "\n"
            << "L = " << L << "\n"
            << "attitude_angle_deg = " << attitude_angle_deg << "\n"
            << "load_angle_deg = " << load_angle_deg << "\n"
            << "tilt_slope_x = " << tilt_slope_x << "\n"
            << "tilt_slope_y = " << tilt_slope_y << "\n"
            << "contact_model = " << (contact_model == ContactModel::GREENWOOD_TRIPP ? "GREENWOOD_TRIPP" : "NONE") << "\n"
            << "contact_roughness_sigma = " << contact_roughness_sigma << "\n"
            << "contact_gt_k = " << contact_gt_k << "\n"
            << "contact_modulus = " << contact_modulus << "\n"
            << "contact_boundary_friction = " << contact_boundary_friction << "\n\n";

        out << "# Grid\n"
            << "n_theta_global = " << n_theta_global << "\n"
            << "n_z_global = " << n_z_global << "\n\n";

        out << "# Time\n"
            << "end_t = " << end_t << "\n"
            << "dt = " << dt << "\n"
            << "write_interval = " << write_interval << "\n"
            << "# solution_mode options: TRANSIENT, STEADY_STATE\n"
            << "solution_mode = " << to_config_value(solution_mode) << "\n"
            << "# Field solves are backward-Euler; only bearing motion has a selectable integrator.\n"
            << "# motion_time_method options: EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICOLSON, RK2, RK4\n"
            << "motion_time_method = " << to_config_value(motion_time_method) << "\n\n";

        out << "# Numerics\n"
            << "# Convection scheme options: UPWIND, LINEAR, VANLEER;\n"
            << "# theta additionally accepts TYPE_DIFFERENCING (central in full film, upwind in cavitation).\n"
            << "theta_convection_scheme = " << to_config_value(theta_convection_scheme) << "\n"
            << "thermal_convection_scheme = " << to_config_value(thermal_convection_scheme) << "\n"
            << "gas_convection_scheme = " << to_config_value(gas_convection_scheme) << "\n"
            << "linear_rtol = " << linear_rtol << "\n"
            << "# Inner Krylov solver: BCGS, BCGSL, GMRES, FGMRES.\n"
            << "# Preconditioner: ASM (default, rank-robust), BJACOBI, JACOBI, NONE.\n"
            << "linear_ksp_type = " << to_config_value(linear_ksp_type) << "\n"
            << "linear_pc_type = " << to_config_value(linear_pc_type) << "\n"
            << "linear_pc_overlap = " << linear_pc_overlap << "\n"
            << "# Escalation rungs (0-2) tried before a failed linear solve aborts the run.\n"
            << "linear_solver_fallbacks = " << linear_solver_fallbacks << "\n"
            << "# Per-step mass-balance / convergence CSV cadence in steps; <= 0 disables.\n"
            << "diagnostics_interval = " << diagnostics_interval << "\n"
            << "steady_output_interval = " << steady_output_interval << "\n"
            << "steady_pressure_inner_only = " << (steady_pressure_inner_only ? "true" : "false") << "\n\n";

        out << "# Physics\n"
            << "omega = " << omega << "\n"
            << "omega_ramp_time = " << omega_ramp_time << "\n"
            << "mu = " << mu << "\n"
            << "rho = " << rho << "\n"
            << "p_cav = " << p_cav << "\n"
            << "oil_vapor_temperature_coeff = " << oil_vapor_temperature_coeff << "\n"
            << "bulk_modulus = " << bulk_modulus << "\n\n";

        out << "# Energy Equation\n"
            << "# temperature_model options: ISOTHERMAL, ENERGY_EQUATION\n"
            << "temperature_model = " << to_config_value(temperature_model) << "\n"
            << "temperature_initial = " << temperature_initial << "\n"
            << "temperature_reference = " << temperature_reference << "\n"
            << "journal_wall_temperature = " << journal_wall_temperature << "\n"
            << "bearing_wall_temperature = " << bearing_wall_temperature << "\n"
            << "rho_cp = " << rho_cp << "\n"
            << "liquid_specific_heat = " << liquid_specific_heat << "\n"
            << "thermal_conductivity = " << thermal_conductivity << "\n"
            << "gas_heat_capacity_cv = " << gas_heat_capacity_cv << "\n"
            << "dissolved_heat_capacity = " << dissolved_heat_capacity << "\n"
            << "desorption_internal_energy = " << desorption_internal_energy << "\n"
            << "gas_thermal_conductivity = " << gas_thermal_conductivity << "\n"
            << "journal_heat_transfer = " << journal_heat_transfer << "\n"
            << "bearing_heat_transfer = " << bearing_heat_transfer << "\n"
            << "# cross_film_thermal_model options: NONE, POLYNOMIAL. POLYNOMIAL turns on the\n"
            << "# cross-film T(y) profile, the conjugate shaft/bush walls (which replace the\n"
            << "# two Newton-cooling coefficients above) and the Dowson F0/F1/F2 viscosity\n"
            << "# integrals in the Reynolds equation. NONE is bit-identical to the model-free solver.\n"
            << "cross_film_thermal_model = " << to_config_value(cross_film_thermal_model) << "\n"
            << "cross_film_layers = " << cross_film_layers << "\n"
            << "cross_film_profile_order = " << cross_film_profile_order << "\n"
            << "cross_film_quadrature_points = " << cross_film_quadrature_points << "\n"
            << "cross_film_relaxation = " << cross_film_relaxation << "\n"
            << "shaft_ambient_heat_transfer = " << shaft_ambient_heat_transfer << "\n"
            << "shaft_ambient_area_ratio = " << shaft_ambient_area_ratio << "\n"
            << "thermal_ambient_temperature = " << thermal_ambient_temperature << "\n"
            << "# <= 0 uses thermal_ambient_temperature\n"
            << "shaft_ambient_temperature = " << shaft_ambient_temperature << "\n"
            << "bush_thermal_conductivity = " << bush_thermal_conductivity << "\n"
            << "bush_outer_radius = " << bush_outer_radius << "\n"
            << "bush_outer_heat_transfer = " << bush_outer_heat_transfer << "\n"
            << "# <= 0 uses thermal_ambient_temperature\n"
            << "bush_outer_temperature = " << bush_outer_temperature << "\n"
            << "# bush_thermal_model options: LOCAL (per-cell 1-D radial closure),\n"
            << "# CONDUCTION_2D (steady (r, theta) conduction in the bush annulus).\n"
            << "bush_thermal_model = " << to_config_value(bush_thermal_model) << "\n"
            << "bush_radial_cells = " << bush_radial_cells << "\n"
            << "bush_axial_cells = " << bush_axial_cells << "\n"
            << "bush_end_heat_transfer = " << bush_end_heat_transfer << "\n"
            << "bush_end_natural_convection = " << (bush_end_natural_convection ? "true" : "false") << "\n"
            << "# bush_outer_boundary options: CONVECTION (Robin), FIXED_TEMPERATURE\n"
            << "bush_outer_boundary = " << to_config_value(bush_outer_boundary) << "\n"
            << "# groove_mixing_model options: PINNED, FIXED_LAMBDA, FLOW_WEIGHTED, FLUX_BALANCED\n"
            << "groove_mixing_model = " << to_config_value(groove_mixing_model) << "\n\n";

        out << "# Fluid Properties\n"
            << "# fluid_property_model options: CONSTANT, SINGLE_PHASE, TWO_PHASE\n"
            << "fluid_property_model = " << to_config_value(fluid_property_model) << "\n"
            << "# dissolved_gas_species options: AIR, PROPANE, R410A\n"
            << "dissolved_gas_species = " << to_config_value(dissolved_gas_species) << "\n"
            << "# solubility_model options: HENRY, TABLE, CURVE_TABLE\n"
            << "solubility_model = " << to_config_value(solubility_model) << "\n"
            << "# density_model options: CONSTANT, MIXTURE, TABLE, LINEAR_TEMPERATURE\n"
            << "density_model = " << to_config_value(density_model) << "\n"
            << "# Linear density [kg/m^3] = intercept + slope * (T[K] - 273.15).\n"
            << "density_temperature_intercept = " << density_temperature_intercept << "\n"
            << "density_temperature_slope = " << density_temperature_slope << "\n"
            << "# liquid_viscosity_model options: CONSTANT, EMPIRICAL, POLYNOMIAL, TABLE, COMPOSITION_TABLE, KINEMATIC_TABLE\n"
            << "liquid_viscosity_model = " << to_config_value(liquid_viscosity_model) << "\n"
            << "dissolved_gas_initial = " << dissolved_gas_initial << "\n"
            << "dissolved_gas_max = " << dissolved_gas_max << "\n"
            << "dissolved_gas_henry_coeff = " << dissolved_gas_henry_coeff << "\n"
            << "dissolved_gas_henry_temp_coeff = " << dissolved_gas_henry_temp_coeff << "\n"
            << "dissolved_gas_liquid_density = " << dissolved_gas_liquid_density << "\n"
            << "gas_mass_transfer_rate = " << gas_mass_transfer_rate << "\n"
            << "# gas_absorption_rate: -1 inherits release rate; 0 disables absorption\n"
            << "gas_absorption_rate = " << gas_absorption_rate << "\n"
            << "coupled_phase_equilibrium = " << (coupled_phase_equilibrium ? "true" : "false") << "\n"
            << "dissolved_gas_diffusivity = " << dissolved_gas_diffusivity << "\n"
            << "solution_density_gas_coeff = " << solution_density_gas_coeff << "\n"
            << "solution_viscosity_gas_coeff = " << solution_viscosity_gas_coeff << "\n"
            << "viscosity_temperature_coeff = " << viscosity_temperature_coeff << "\n"
            << "viscosity_pressure_coeff = " << viscosity_pressure_coeff << "\n"
            << "# POLYNOMIAL (Ferron 1983): mu(T) = mu * (k0 - k1 Tbar + k2 Tbar^2) with\n"
            << "# Tbar = T[degC] / property_reference_temperature[degC]. The ratio is in\n"
            << "# CELSIUS (a Kelvin ratio does not reproduce the paper's numbers); the config\n"
            << "# key itself stays in Kelvin. Ferron's VG32 set: 3.287, 3.064, 0.777 with\n"
            << "# mu = 0.0277 Pa.s and property_reference_temperature = 313.15 (40 degC).\n"
            << "viscosity_poly_k0 = " << viscosity_poly_k0 << "\n"
            << "viscosity_poly_k1 = " << viscosity_poly_k1 << "\n"
            << "viscosity_poly_k2 = " << viscosity_poly_k2 << "\n"
            << "gas_alpha_max = " << gas_alpha_max << "\n"
            << "gas_pressure_floor = " << gas_pressure_floor << "\n"
            << "# mixture_viscosity_model options: DUKLER, LINEAR, MCADAMS\n"
            << "mixture_viscosity_model = " << to_config_value(mixture_viscosity_model) << "\n"
            << "mu_gas = " << mu_gas << "\n"
            << "property_reference_pressure = " << property_reference_pressure << "\n"
            << "property_reference_temperature = " << property_reference_temperature << "\n"
            << "# Tables use x:value pairs. Solubility is keyed by pressure; density/viscosity by dissolved mass fraction.\n"
            << "solubility_table =" << (solubility_table.empty() ? "" : std::string(" ") + property_table_text(solubility_table)) << "\n"
            << "density_table =" << (density_table.empty() ? "" : std::string(" ") + property_table_text(density_table)) << "\n"
            << "viscosity_table =" << (viscosity_table.empty() ? "" : std::string(" ") + property_table_text(viscosity_table)) << "\n"
            << "# 2-D (P_MPa, T_C, value) CSV files; when set they take precedence over the inline\n"
            << "# tables above. Relative paths resolve against this config file's directory.\n"
            << "solubility_table_file = " << solubility_table_file.generic_string() << "\n"
            << "viscosity_table_file = " << viscosity_table_file.generic_string() << "\n"
            << "viscosity_composition_file = " << viscosity_composition_file.generic_string() << "\n"
            << "# CURVE_TABLE solubility and KINEMATIC_TABLE nu(c,T), using absolute P and refrigerant mass fraction.\n"
            << "solution_table_file = " << solution_table_file.generic_string() << "\n"
            << "solution_table_clamp = " << (solution_table_clamp ? "true" : "false") << "\n"
            << "kinematic_viscosity_log_interpolation = " << (kinematic_viscosity_log_interpolation ? "true" : "false") << "\n"
            << "density_table_file = " << density_table_file.generic_string() << "\n\n";

        out << "# Bearing Motion\n"
            << "# motion_model options: STATIC, MOVING_BEARING, STEADY_EQUILIBRIUM\n"
            << "motion_model = " << to_config_value(motion_model) << "\n"
            << "steady_force_tolerance = " << steady_force_tolerance << "\n"
            << "steady_position_step_limit = " << steady_position_step_limit << "\n"
            << "steady_position_fd_step = " << steady_position_fd_step << "\n"
            << "bearing_initial_from_attitude = " << (bearing_initial_from_attitude ? "true" : "false") << "\n"
            << "bearing_initial_x = " << bearing_initial_x << "\n"
            << "bearing_initial_y = " << bearing_initial_y << "\n"
            << "bearing_initial_z = " << bearing_initial_z << "\n"
            << "bearing_initial_vx = " << bearing_initial_vx << "\n"
            << "bearing_initial_vy = " << bearing_initial_vy << "\n"
            << "bearing_initial_vz = " << bearing_initial_vz << "\n"
            << "bearing_mass = " << bearing_mass << "\n"
            << "bearing_stiffness_x = " << bearing_stiffness_x << "\n"
            << "bearing_stiffness_y = " << bearing_stiffness_y << "\n"
            << "bearing_stiffness_z = " << bearing_stiffness_z << "\n"
            << "bearing_damping_x = " << bearing_damping_x << "\n"
            << "bearing_damping_y = " << bearing_damping_y << "\n"
            << "bearing_damping_z = " << bearing_damping_z << "\n"
            << "external_load_x = " << external_load_x << "\n"
            << "external_load_y = " << external_load_y << "\n"
            << "external_load_z = " << external_load_z << "\n"
            << "min_film_thickness = " << min_film_thickness << "\n"
            << "stop_on_nonpositive_film = " << (stop_on_nonpositive_film ? "true" : "false") << "\n\n";

        out << "# Axial pressure reservoirs (inflow and outflow are both permitted)\n"
            << "# Options: DIRICHLET, NEUMANN, INLET_OUTLET\n"
            << "# For JFO, boundary film content is derived from bc_z_*_val,\n"
            << "# p_cav, and bulk_modulus; legacy bc_z_*_theta inputs are ignored.\n"
            << "# bc_z_*_thermal: CONSTANT uses the reservoir temperature on inflow only;\n"
            << "#   ZERO_GRADIENT extrapolates interior temperature even during reverse flow.\n"
            << "# Temperature 0 and dissolved_gas -1 inherit their global reference/initial values.\n"
            << "# Free-gas mass fraction is free gas / (liquid solution + free gas), not void fraction.\n"
            << "bc_z_south_type = " << to_config_value(bc_z_south_type) << "\n"
            << "bc_z_south_val = " << bc_z_south_val << "\n"
            << "bc_z_south_thermal = " << to_config_value(bc_z_south_thermal) << "\n"
            << "bc_z_south_temperature = " << bc_z_south_temperature << "\n"
            << "bc_z_south_dissolved_gas = " << bc_z_south_dissolved_gas << "\n"
            << "bc_z_south_free_gas_mass_fraction = " << bc_z_south_free_gas_mass_fraction << "\n"
            << "bc_z_north_type = " << to_config_value(bc_z_north_type) << "\n"
            << "bc_z_north_val = " << bc_z_north_val << "\n"
            << "bc_z_north_thermal = " << to_config_value(bc_z_north_thermal) << "\n"
            << "bc_z_north_temperature = " << bc_z_north_temperature << "\n"
            << "bc_z_north_dissolved_gas = " << bc_z_north_dissolved_gas << "\n"
            << "bc_z_north_free_gas_mass_fraction = " << bc_z_north_free_gas_mass_fraction << "\n\n";

        out << "# Cavitation\n"
            << "# Options: FULL_SOMMERFELD, GUMBEL, JFO, FBNS\n"
            << "cavitation_model = " << to_config_value(cavitation_model) << "\n"
            << "# cavitation_threshold: SCALAR_PCAV (oil p_cav) or LOCAL_PSAT (max of p_cav\n"
            << "# and the dissolved-species bubble point, so the first-released species governs onset)\n"
            << "cavitation_threshold = "
            << (cavitation_threshold == CavitationThreshold::LOCAL_PSAT ? "LOCAL_PSAT" : "SCALAR_PCAV") << "\n"
            << "# gas_pressure_coupling: NONE (released gas volumetrically inert) or VOID_COUPLED\n"
            << "# (released gas displaces liquid and sets a local plateau p_void; needs TWO_PHASE + mu_gas)\n"
            << "gas_pressure_coupling = "
            << (gas_pressure_coupling == GasPressureCoupling::VOID_COUPLED ? "VOID_COUPLED" : "NONE") << "\n"
            << "max_outer_iters = " << max_outer_iters << "\n"
            << "outer_tol = " << outer_tol << "\n"
            << "theta_min = " << theta_min << "\n"
            << "log_outer_iters = " << (log_outer_iters ? "true" : "false") << "\n"
            << "# FBNS only (cavitation_model = FBNS). bulk_modulus is unused there.\n"
            << "# fbns_pressure_scale = 0 derives p_ref from the operating point.\n"
            << "fbns_max_iters = " << fbns_max_iters << "\n"
            << "fbns_tol = " << fbns_tol << "\n"
            << "fbns_relaxation = " << fbns_relaxation << "\n"
            << "fbns_pressure_scale = " << fbns_pressure_scale << "\n\n";

        out << "# Inlets\n"
            << "# inlet_circular = theta(deg) z(m) radius(m) p_supply(Pa) [t_supply(K)] [lambda] [q_supply(m3/s)]\n"
            << "# inlet_groove   = theta(deg) width(deg) p_supply(Pa) [t_supply(K)] [lambda] [q_supply(m3/s)]\n"
            << "# lambda and q_supply feed groove_mixing_model (FIXED_LAMBDA / FLOW_WEIGHTED)\n"
            << "# and are ignored under the default PINNED.\n"
            << "# Omit the optional t_supply to model an open large-clearance region of the\n"
            << "# same recirculating oil (pressure-only, thermally transparent). Provide it to\n"
            << "# model an actual fed inlet that pumps fresh oil in at t_supply.\n";
        if (inlets.empty()) {
            out << "# inlet_groove = 90.0 10.0 2e5\n";
        } else {
            for (const auto& inlet : inlets) {
                if (inlet.type == InletConfig::Type::CIRCULAR) {
                    out << "inlet_circular = " << inlet.theta << " " << inlet.z << " "
                        << inlet.size << " " << inlet.p_supply;
                } else {
                    out << "inlet_groove = " << inlet.theta << " "
                        << inlet.size << " " << inlet.p_supply;
                }
                if (inlet.feeds_fresh_oil) {
                    out << " " << inlet.t_supply;
                    if (groove_mixing_model != GrooveMixingModel::PINNED) {
                        out << " " << inlet.mixing_lambda << " " << inlet.supply_flow_rate;
                    }
                }
                out << "\n";
            }
        }
        out << "\n";

        for (std::size_t index = 0; index < inlets.size(); ++index)
            if (inlets[index].dissolved_gas != -1 || inlets[index].free_gas_mass_fraction != 0)
                out << "inlet_species = " << index << " " << inlets[index].dissolved_gas << " " << inlets[index].free_gas_mass_fraction << "\n";

        out << "# Output\n"
            << "# output_fields accepts comma-separated field names.\n"
            << "# Options: pressure, film_content, h, rho, mu,\n"
            << "#          rho_liquid_solution, mu_liquid_solution,\n"
            << "#          temperature, heat_generation, dissolved_gas, alpha_gas,\n"
            << "#          gas_mass_transfer,\n"
            << "#          inlet_indicator, velocity (active vector),\n"
            << "#          pressure_force_x, pressure_force_y, pressure_force_z,\n"
            << "#          viscous_force_x, viscous_force_y, viscous_force_z,\n"
            << "#          fluid_force_x, fluid_force_y, fluid_force_z,\n"
            << "#          external_load_x, external_load_y, external_load_z,\n"
            << "#          bearing_x, bearing_y, bearing_z, bearing_attitude_angle,\n"
            << "#          friction_torque,\n"
            << "#          mass_residual (per-cell continuity residual [kg/s], JFO/FBNS only),\n"
            << "#          temperature_journal, temperature_bearing (film surface temperatures\n"
            << "#          [K] at eta = 0 / 1), heat_flux_journal, heat_flux_bearing (wall heat\n"
            << "#          flux [W/m^2], positive INTO the fluid). The last four require\n"
            << "#          cross_film_thermal_model = POLYNOMIAL.\n"
            << "#          Legacy load_x/load_y/load_z requests map to pressure_force_*.\n"
            << "output_dir = " << output_dir << "\n"
            << "steady_initial_state_file = " << steady_initial_state_file << "\n"
            << "filename_prefix = " << filename_prefix << "\n"
            << "# output_verbosity options: SIMPLIFIED (compact per-step line), DETAILED (per-step diagnostics)\n"
            << "output_verbosity = " << to_config_value(output_verbosity) << "\n"
            << "output_write_3d = " << (output_write_3d ? "true" : "false") << "\n"
            << "output_write_flat = " << (output_write_flat ? "true" : "false") << "\n"
            << "output_fields = " << output_fields_text() << "\n";

        return out.str();
    }

    /// Every key `load_from_stream` accepts, including the deliberate legacy
    /// aliases (`viscosity_model`, `oil_gas_solution_model`, ...) and the
    /// warn-and-drop branches for removed models (`pressure_time_method`,
    /// `dissolved_gas_bunsen_coeff`). This is the manifest the unknown-key error
    /// suggests from; a new parser branch must be added here too, and a missing
    /// entry is loud rather than silent because the key then fails to parse.
    static const std::vector<std::string>& known_config_keys() {
        static const std::vector<std::string> keys = {
            "R", "c", "e", "L", "attitude_angle_deg", "load_angle_deg",
            "n_theta_global", "n_z_global", "end_t", "dt", "write_interval",
            "solution_mode", "pressure_time_method", "temperature_time_method",
            "motion_time_method", "theta_convection_scheme",
            "thermal_convection_scheme", "gas_convection_scheme", "linear_rtol",
            "linear_ksp_type", "linear_pc_type", "linear_pc_overlap",
            "linear_solver_fallbacks", "omega", "omega_ramp_time", "mu", "rho",
            "p_cav", "oil_vapor_temperature_coeff", "bulk_modulus", "tilt_slope_x", "tilt_slope_y",
            "contact_model", "contact_roughness_sigma", "contact_gt_k", "contact_modulus",
            "contact_boundary_friction",
            "temperature_model", "temperature_initial", "temperature_reference",
            "journal_wall_temperature", "bearing_wall_temperature", "rho_cp", "liquid_specific_heat",
            "thermal_conductivity", "gas_heat_capacity_cv", "dissolved_heat_capacity", "desorption_internal_energy", "gas_thermal_conductivity", "journal_heat_transfer",
            "bearing_heat_transfer", "cross_film_thermal_model",
            "cross_film_layers", "cross_film_profile_order", "cross_film_quadrature_points",
            "cross_film_relaxation", "shaft_ambient_heat_transfer",
            "shaft_ambient_area_ratio", "thermal_ambient_temperature",
            "shaft_ambient_temperature",
            "bush_thermal_conductivity", "bush_outer_radius",
            "bush_outer_heat_transfer", "bush_outer_temperature",
            "bush_thermal_model", "bush_radial_cells", "bush_axial_cells", "bush_end_heat_transfer", "bush_end_natural_convection",
            "bush_outer_boundary", "groove_mixing_model", "fluid_property_model",
            "dissolved_gas_species", "solubility_model", "oil_gas_solution_model",
            "density_model", "liquid_viscosity_model", "viscosity_model",
            "density_temperature_intercept", "density_temperature_slope",
            "solution_table_file", "solution_table_clamp", "kinematic_viscosity_log_interpolation",
            "dissolved_gas_initial", "dissolved_gas_max",
            "dissolved_gas_henry_coeff", "dissolved_gas_henry_temp_coeff",
            "dissolved_gas_bunsen_coeff", "dissolved_gas_liquid_density",
            "gas_mass_transfer_rate", "gas_absorption_rate", "coupled_phase_equilibrium", "inlet_species", "dissolved_gas_diffusivity",
            "solution_density_gas_coeff", "solution_viscosity_gas_coeff",
            "viscosity_temperature_coeff", "viscosity_pressure_coeff",
            "viscosity_poly_k0", "viscosity_poly_k1", "viscosity_poly_k2",
            "gas_alpha_max", "gas_pressure_floor", "mixture_viscosity_model",
            "gas_mixture_viscosity_model", "gas_pressure_coupling",
            "void_coupling_relaxation", "mu_gas", "diagnostics_interval", "steady_output_interval", "steady_pressure_inner_only",
            "property_reference_pressure", "property_reference_temperature",
            "solubility_table", "density_table", "viscosity_table",
            "solubility_table_file", "viscosity_table_file", "viscosity_composition_file", "density_table_file",
            "motion_model", "steady_force_tolerance", "steady_position_step_limit", "steady_position_fd_step", "bearing_initial_from_attitude", "bearing_initial_x",
            "bearing_initial_y", "bearing_initial_z", "bearing_initial_vx",
            "bearing_initial_vy", "bearing_initial_vz", "bearing_mass",
            "bearing_stiffness_x", "bearing_stiffness_y", "bearing_stiffness_z",
            "bearing_damping_x", "bearing_damping_y", "bearing_damping_z",
            "external_load_x", "external_load_y", "external_load_z",
            "min_film_thickness", "stop_on_nonpositive_film", "bc_z_south_type",
            "bc_z_south_val", "bc_z_south_theta", "bc_z_south_thermal",
            "bc_z_north_type", "bc_z_north_val", "bc_z_north_theta",
            "bc_z_north_thermal", "cavitation_model", "cavitation_threshold",
            "bc_z_south_temperature", "bc_z_north_temperature",
            "bc_z_south_dissolved_gas", "bc_z_north_dissolved_gas",
            "bc_z_south_free_gas_mass_fraction", "bc_z_north_free_gas_mass_fraction",
            "cavitated_film_model", "consistent_boundary_flux",
            "coupling_max_iters", "coupling_tolerance", "coupling_relaxation",
            "steady_gas_model",
            "allow_frozen_gas_steady", "max_outer_iters", "outer_tol",
            "fbns_max_iters", "fbns_tol", "fbns_relaxation",
            "fbns_pressure_scale", "theta_min", "log_outer_iters", "output_dir", "steady_initial_state_file",
            "filename_prefix", "output_verbosity", "output_write_3d",
            "output_write_flat", "output_fields", "inlet_circular", "inlet_groove",
        };
        return keys;
    }

    /// Levenshtein edit distance, used only to name the likely intended key.
    static int key_edit_distance(const std::string& a, const std::string& b) {
        const size_t n = b.size();
        std::vector<int> previous(n + 1), current(n + 1);
        for (size_t j = 0; j <= n; ++j) previous[j] = static_cast<int>(j);
        for (size_t i = 1; i <= a.size(); ++i) {
            current[0] = static_cast<int>(i);
            for (size_t j = 1; j <= n; ++j) {
                const int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
                current[j] = std::min({previous[j] + 1, current[j - 1] + 1,
                                       previous[j - 1] + cost});
            }
            previous = current;
        }
        return previous[n];
    }

    /// Underscore-separated words of a config key, for the suggestion heuristic.
    static std::vector<std::string> key_tokens(const std::string& key) {
        std::vector<std::string> tokens;
        std::string current;
        for (char ch : key) {
            if (ch == '_') {
                if (!current.empty()) tokens.push_back(current);
                current.clear();
            } else {
                current.push_back(ch);
            }
        }
        if (!current.empty()) tokens.push_back(current);
        return tokens;
    }

    /// Known key closest to `key`, or "" when nothing is close enough to be a
    /// plausible mistake.
    ///
    /// Two signals, because character distance alone is not enough. A typo
    /// (`bush_outer_radus`) is one or two edits away, so a length-scaled edit
    /// distance catches it. A *wrong name* for the right knob (`cross_film_model`
    /// for `cross_film_thermal_model`) is eight edits away but shares every one
    /// of its underscore-separated words, which is the case that actually shipped.
    static std::string closest_config_key(const std::string& key) {
        const std::vector<std::string> tokens = key_tokens(key);
        std::string best;
        int best_distance = std::numeric_limits<int>::max();
        size_t best_shared = 0;
        for (const std::string& candidate : known_config_keys()) {
            const int distance = key_edit_distance(key, candidate);
            const std::vector<std::string> candidate_tokens = key_tokens(candidate);
            size_t shared = 0;
            for (const std::string& token : tokens) {
                if (std::find(candidate_tokens.begin(), candidate_tokens.end(), token) !=
                    candidate_tokens.end()) {
                    ++shared;
                }
            }
            // More shared words wins; equal words fall back to fewer edits.
            if (shared > best_shared || (shared == best_shared && distance < best_distance)) {
                best_shared = shared;
                best_distance = distance;
                best = candidate;
            }
        }
        const int distance_limit = std::max(2, static_cast<int>(key.size()) / 3);
        const size_t shared_limit = std::max<size_t>(2, (tokens.size() + 1) / 2);
        if (best_distance <= distance_limit || best_shared >= shared_limit) return best;
        return std::string();
    }

    static std::string unknown_key_message(const std::string& key) {
        std::string message = "unrecognized config key '" + key + "'";
        const std::string suggestion = closest_config_key(key);
        if (!suggestion.empty()) message += "; did you mean '" + suggestion + "'?";
        return message;
    }

private:
    static std::string normalise_key(std::string text) {
        trim(text);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });
        return text;
    }

    static std::vector<std::string> parse_output_fields(std::string value) {
        std::replace(value.begin(), value.end(), ',', ' ');
        std::replace(value.begin(), value.end(), ';', ' ');
        std::replace(value.begin(), value.end(), '|', ' ');

        std::vector<std::string> fields;
        std::stringstream input(value);
        std::string token;
        while (input >> token) {
            const std::string key = normalise_key(token);
            if (!key.empty() && std::find(fields.begin(), fields.end(), key) == fields.end()) {
                fields.push_back(key);
            }
        }
        return fields;
    }

    static std::vector<PropertyTablePoint> parse_property_table(std::string value) {
        std::replace(value.begin(), value.end(), ',', ' ');
        std::replace(value.begin(), value.end(), ';', ' ');
        std::replace(value.begin(), value.end(), '|', ' ');

        std::vector<PropertyTablePoint> table;
        std::stringstream input(value);
        std::string token;
        while (input >> token) {
            const size_t delimiter = token.find(':');
            if (delimiter == std::string::npos) continue;
            std::string x_text = token.substr(0, delimiter);
            std::string value_text = token.substr(delimiter + 1);
            trim(x_text);
            trim(value_text);
            if (x_text.empty() || value_text.empty()) continue;
            table.push_back({std::stod(x_text), std::stod(value_text)});
        }
        std::sort(table.begin(), table.end(), [](const auto& a, const auto& b) {
            return a.x < b.x;
        });
        return table;
    }

    static void trim(std::string& text) {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            text.clear();
            return;
        }
        const auto last = text.find_last_not_of(" \t\r\n");
        text = text.substr(first, last - first + 1);
    }

    static bool is_true(const std::string& value) {
        return value == "true" || value == "TRUE" || value == "True" || value == "1";
    }

    static bool is_false(const std::string& value) {
        return value == "false" || value == "FALSE" || value == "False" || value == "0";
    }
};
