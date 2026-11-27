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

/// Wall-time estimate for the chosen core count (docs/GUI_LGE_PLAN.md 0.2).
struct RunTime {
    double typical_s = 0.0;     // headline value
    double low_s = 0.0;
    double high_s = 0.0;
    double core_hours = 0.0;    // typical core-time (wall x processes)
    double per_revolution_s = 0.0;  // transient: typical wall time of one revolution
    int cores = 1;
    double effective_cores = 1.0;  // processes counted as physical cores
    std::string basis;          // which measured form the estimate scales
};

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
    RunTime run_time;                // wall-time estimate at the case's core count
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

/// Physical cores and hardware threads of this computer (0 = unknown). Set once
/// by the GUI; processes beyond the physical cores count as 0.25 of a core
/// each (hyperthreads), beyond the threads as nothing.
void set_machine_cores(int physical, int logical);
int machine_physical_cores();
int machine_logical_cores();

/// Run-time estimate from the measured core-times at 720 x 192 (plan 0.2),
/// scaled linearly with the cell count, by steps per revolution (transient)
/// and by the core count with parallel efficiency 1 / (1 + 0.02 (p - 1)).
RunTime estimate_run_time(const CaseModel& model, const FluidPair* pair);
/// "about 25 min (12-40 min) on 4 cores".
std::string format_run_time(const RunTime& run_time);

}  // namespace gui
