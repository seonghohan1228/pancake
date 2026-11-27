#pragma once

// Key results of one solved step, computed from the solver's VTK fields, plus
// the design checks shown with them. No toolkit dependency: unit-tested.
//
// Definitions (docs/GUI_GLOSSARY.md):
// - cavitated area: fraction of film cells with film content below 1 (the deck
//   metric; 57.09 % for the compressor example);
// - eccentricity ratio: |bearing displacement| / c;
// - attitude angle phi: angle from the load line to the line of centres
//   (journal displacement relative to the bearing), deck p18;
// - friction power: |friction torque x angular speed|;
// - shear share: |viscous film force| / |total film force|;
// - Lambda: h_min / sqrt(Rq_shaft^2 + Rq_bush^2).

#include "config.hpp"
#include "results.hpp"

#include <limits>
#include <string>
#include <vector>

namespace gui {

enum class Risk { None, Low, Medium, High };  // green / yellow / red

struct KeyResults {
    static constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    double h_min = kNaN;               // m
    double h_min_theta_deg = kNaN;     // fixed frame, where h is smallest (mid-plane)
    double h_max_theta_deg = kNaN;     // fixed frame, where h is largest (mid-plane)
    double eccentricity_ratio = kNaN;
    double attitude_deg = kNaN;        // deck definition, see above
    double friction_torque = kNaN;     // N m
    double friction_power = kNaN;      // W
    double p_max = kNaN;               // Pa
    double p_min = kNaN;               // Pa
    double cavitated_fraction = kNaN;  // 0..1
    double film_force = kNaN;          // N
    double shear_share = kNaN;         // 0..1
    double lambda = kNaN;              // film thickness ratio
    double t_max = kNaN;               // K (thermal runs)
    double t_mean = kNaN;              // K, mean film temperature (thermal runs)
    double time = 0.0;
    int step = 0;
};

/// Computes the key results of one grid. `cfg` gives c, omega, load direction
/// and roughness.
KeyResults compute_key_results(const ResultGrid& grid, const SimulationConfig& cfg);

/// Design checks on a result: film thickness ratio and gas release.
struct DesignChecks {
    Risk lambda_risk = Risk::None;       // Lambda < 3 red, < 5 yellow
    Risk gas_release_risk = Risk::None;  // p_min below p_sat (gas-off run) or cavitated at p_sat
    std::string lambda_text;
    std::string gas_release_text;
};

/// `p_sat` [Pa] is the saturation pressure of the case (0 if unknown);
/// `gas_model_on` says whether the run itself modelled gas release.
DesignChecks design_checks(const KeyResults& results, double p_sat, bool gas_model_on);

/// Field maps offered for a result, by solver form (plan 0.5).
///
/// FBNS: the cavity is liquid streaks at fraction theta plus a massless void at
/// the cavity pressure. No gas mass is tracked, so FBNS offers pressure, film
/// thickness, void fraction 1 - theta, effective density theta x rho_liquid and
/// temperature; never alpha_gas, dissolved_gas, rho or gas_mass_transfer.
/// COUPLED: pressure, film thickness, free-gas volume fraction, dissolved
/// content, gas mass transfer, mixture density and viscosity, temperature.
struct FieldView {
    std::string id;        // "p", "h", "void", "effective_density", or a grid field name
    std::string label;     // UI label
    std::string unit;      // SI unit of the values
    std::string tooltip;   // one line; model caveats
};

enum class ResultForm { Fbns, Coupled, Other };
ResultForm result_form(const SimulationConfig& cfg);

/// Views available for `grid` under the run's model (fields missing from the
/// grid are left out).
std::vector<FieldView> result_fields(const SimulationConfig& cfg, const ResultGrid& grid);

/// Values of a view, row-major like the grid fields (derived views computed).
/// Empty when the view is not available.
std::vector<double> field_view_values(const FieldView& view, const ResultGrid& grid,
                                      const SimulationConfig& cfg);

/// Percent difference (b - a) / |a| x 100; NaN when a is 0 or either is NaN.
double percent_difference(double a, double b);

}  // namespace gui
