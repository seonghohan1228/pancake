#pragma once

// Instant engineering estimates shown before a run (no solver call):
// Sommerfeld number, saturation pressure, pressure margin, the gaseous-
// cavitation onset verdict from the onset chart, and a rough run time.
//
// Onset chart (deck p43): pressure drop ratio Phi(S) = (p_amb - p_min) / P with
// specific load P = W / (L D), from gas-free solutions. Gas release starts when
// the pressure margin (p_amb - p_sat) / P falls below Phi(S). The curve comes
// from the tester bearing (L/D = 0.93, axial groove); for other bearings it is
// an estimate, and the run itself gives the direct check (film p_min vs p_sat).
// No toolkit dependency: unit-tested.

#include "case_model.hpp"
#include "fluid_library.hpp"

#include <string>
#include <vector>

namespace gui {

/// Design-risk verdict. Colours follow the risk: Expected = red, Borderline =
/// yellow, NotExpected = green (docs/GUI_LGE_PLAN.md section 6).
enum class OnsetVerdict { Unknown, NotExpected, Borderline, Expected };
const char* onset_verdict_text(OnsetVerdict verdict);

struct OnsetPoint {
    double sommerfeld = 0.0;
    double pressure_drop_ratio = 0.0;  // Phi
};

/// The onset curve (Phi vs S), sorted by S, plus where it came from.
struct OnsetCurve {
    std::vector<OnsetPoint> points;
    std::string source;
    /// Phi at S, log-log interpolation; NaN outside the curve's S range.
    double phi(double sommerfeld) const;
    bool covers(double sommerfeld) const;
};
/// The curve shipped with the GUI (tester bearing, gas-free runs G1-G8).
const OnsetCurve& tester_onset_curve();

struct Estimates {
    double viscosity = 0.0;          // Pa s at ambient T and content
    double specific_load = 0.0;      // P = W / (L D) [Pa]
    double sommerfeld = 0.0;         // S = mu N L D (R/c)^2 / W, N in rev/s
    double content = 0.0;            // kg/kg
    double saturation_pressure = 0.0;  // Pa
    double pressure_margin = 0.0;    // p_amb - p_sat [Pa]
    double margin_ratio = 0.0;       // (p_amb - p_sat) / P
    double onset_phi = 0.0;          // Phi(S) from the chart
    OnsetVerdict verdict = OnsetVerdict::Unknown;
    bool fluid_in_range = true;      // ambient state inside the fluid table
    double run_time_s = 0.0;         // rough wall-clock estimate
    std::string note;                // why something is unknown
};

/// Borderline band around the onset curve: |margin/P - Phi| <= band * Phi.
constexpr double kOnsetBorderlineBand = 0.2;

Estimates compute_estimates(const CaseModel& model, const FluidPair* pair);

/// Sommerfeld number S = mu N L D (R/c)^2 / W with N in rev/s.
double sommerfeld_number(double viscosity, double speed_rad_s, double length, double diameter,
                         double radial_clearance, double load);

/// Verdict from the margin ratio and Phi(S).
OnsetVerdict onset_verdict(double margin_ratio, double phi);

/// Rough run time [s]: cells x iterations x cost per cell-iteration / cores.
double estimate_run_time(const CaseModel& model);

}  // namespace gui
