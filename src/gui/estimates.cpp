#include "estimates.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace gui {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

const char* onset_verdict_text(OnsetVerdict verdict) {
    switch (verdict) {
        case OnsetVerdict::Expected: return "Gas release expected";
        case OnsetVerdict::Borderline: return "Gas release borderline";
        case OnsetVerdict::NotExpected: return "Gas release not expected";
        case OnsetVerdict::Unknown: break;
    }
    return "Gas release: no estimate";
}

double OnsetCurve::phi(double sommerfeld) const {
    if (!covers(sommerfeld)) return std::numeric_limits<double>::quiet_NaN();
    const auto upper = std::upper_bound(points.begin(), points.end(), sommerfeld,
                                        [](double s, const OnsetPoint& p) { return s < p.sommerfeld; });
    if (upper == points.end()) return points.back().pressure_drop_ratio;
    if (upper == points.begin()) return points.front().pressure_drop_ratio;
    const OnsetPoint& a = *(upper - 1);
    const OnsetPoint& b = *upper;
    const double w = (std::log(sommerfeld) - std::log(a.sommerfeld)) /
                     (std::log(b.sommerfeld) - std::log(a.sommerfeld));
    return std::exp(std::lerp(std::log(a.pressure_drop_ratio), std::log(b.pressure_drop_ratio), w));
}

bool OnsetCurve::covers(double sommerfeld) const {
    return !points.empty() && std::isfinite(sommerfeld) &&
           sommerfeld >= points.front().sommerfeld && sommerfeld <= points.back().sommerfeld;
}

const OnsetCurve& tester_onset_curve() {
    // Phi = (p_amb - p_min) / P of the converged gas-free COUPLED runs G1-G8,
    // lge-onset-chart-20261116/gas_off/G*/accepted.json (isothermal 24.4 degC,
    // c = 0.274, chamber 200 bar, tester R = 21.5 mm, L = 40 mm, groove feed).
    // G4, G5 and G9 (S = 0.2 at two loads and at 16.6 bar) give the same Phi,
    // so the curve is dimensionless for this bearing.
    static const OnsetCurve curve = {
        {{0.02, 5.7280}, {0.05, 2.7481}, {0.1, 1.7644}, {0.2, 1.2137},
         {0.5, 0.8171}, {1.0, 0.6746}, {3.0, 0.5822}},
        "Tester bearing (L/D = 0.93, axial groove), gas-free runs, lge-onset-chart-20261116"};
    return curve;
}

double sommerfeld_number(double viscosity, double speed_rad_s, double length, double diameter,
                         double radial_clearance, double load) {
    if (load <= 0.0 || radial_clearance <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    const double revolutions_per_s = std::abs(speed_rad_s) / (2.0 * kPi);
    const double ratio = 0.5 * diameter / radial_clearance;
    return viscosity * revolutions_per_s * length * diameter * ratio * ratio / load;
}

OnsetVerdict onset_verdict(double margin_ratio, double phi) {
    if (!std::isfinite(margin_ratio) || !std::isfinite(phi) || phi <= 0.0) return OnsetVerdict::Unknown;
    if (std::abs(margin_ratio - phi) <= kOnsetBorderlineBand * phi) return OnsetVerdict::Borderline;
    return margin_ratio < phi ? OnsetVerdict::Expected : OnsetVerdict::NotExpected;
}

RunTime estimate_run_time(const CaseModel& m, const FluidPair* pair) {
    // Measured core-times at 720 x 192 on afml-r6 (rev. 3; core-time = wall x
    // processes): {low, typical, high} in core-hours.
    struct Range { double low, typical, high; };
    const bool coupled = uses_coupled(m, pair);
    RunTime r;
    r.cores = std::max(1, m.cores);
    Range core_h{};
    if (m.analysis == AnalysisType::Steady) {
        if (!m.thermal && !coupled) {
            core_h = {5.0 / 60.0, 10.0 / 60.0, 15.0 / 60.0};
            r.basis = "isothermal steady FBNS: 5-15 core-min";
        } else if (!m.thermal) {
            core_h = {0.8, 9.6, 29.6};  // 2-74 min wall, median 24 min, on 24 processes
            r.basis = "isothermal steady COUPLED: 2-74 min on 24 processes";
        } else if (!coupled) {
            core_h = {0.8, 2.5, 8.0};
            r.basis = "thermal steady FBNS: 0.8-8 core-h";
        } else {
            core_h = {10.0, 22.0, 50.0};
            r.basis = "thermal steady COUPLED: 10-50 core-h";
        }
    } else {
        // FBNS transient per revolution: 0.5 core-h at 360 steps, 3.7 core-h at
        // 1440 steps -> C(n) = 0.5 (n / 360)^1.44.
        const double steps = std::max(1, m.steps_per_revolution);
        const double per_rev = 0.5 * std::pow(steps / 360.0, 1.44);
        Range factor = {0.7, 1.0, 1.4};
        r.basis = "isothermal transient FBNS: 0.5 core-h/rev at 1 deg, 3.7 at 0.25 deg";
        if (m.thermal) {
            // Simplified thermal transient, measured 1.0-1.5 x isothermal per step.
            factor = {1.0, 1.5, 2.0};
            r.basis += "; thermal x1.0-2.0 (measured per step)";
        }
        const double revolutions = std::max(0.0, m.revolutions);
        core_h = {per_rev * factor.low * revolutions, per_rev * factor.typical * revolutions,
                  per_rev * factor.high * revolutions};
        r.per_revolution_s = per_rev * factor.typical * 3600.0;
    }
    int n_theta = m.n_theta, n_z = m.n_z;
    if (m.mesh_level != MeshLevel::Custom) mesh_level_cells(m.mesh_level, n_theta, n_z);
    const double cell_ratio = static_cast<double>(n_theta) * n_z / (720.0 * 192.0);
    const double p = r.cores;
    const double processes = p / (1.0 + 0.02 * (p - 1.0));  // p x efficiency
    const auto wall = [&](double core_hours) { return core_hours * cell_ratio * 3600.0 / processes; };
    r.low_s = wall(core_h.low);
    r.typical_s = wall(core_h.typical);
    r.high_s = wall(core_h.high);
    r.core_hours = core_h.typical * cell_ratio;
    r.per_revolution_s *= cell_ratio / processes;
    return r;
}

std::string format_run_time(const RunTime& r) {
    const auto text = [](double seconds) {
        char buffer[32];
        if (seconds < 90.0) std::snprintf(buffer, sizeof(buffer), "%.0f s", seconds);
        else if (seconds < 5400.0) std::snprintf(buffer, sizeof(buffer), "%.0f min", seconds / 60.0);
        else std::snprintf(buffer, sizeof(buffer), "%.1f h", seconds / 3600.0);
        return std::string(buffer);
    };
    return "about " + text(r.typical_s) + " (" + text(r.low_s) + " to " + text(r.high_s) + ") on " +
           std::to_string(r.cores) + (r.cores == 1 ? " core" : " cores");
}

Estimates compute_estimates(const CaseModel& m, const FluidPair* pair) {
    Estimates e;
    e.run_time = estimate_run_time(m, pair);
    e.specific_load = m.bearing_length > 0.0 && m.bearing_diameter > 0.0
                          ? m.load / (m.bearing_length * m.bearing_diameter)
                          : 0.0;
    if (!pair || !pair->loaded()) {
        e.note = "Choose a fluid pair for the estimates.";
        return e;
    }
    const PropertyValue content = m.saturated
        ? pair->saturated_content(m.ambient_pressure, m.ambient_temperature)
        : PropertyValue{m.content, m.content >= pair->content_min() && m.content <= pair->content_max()};
    e.content = content.value;
    const PropertyValue mu = pair->dynamic_viscosity(e.content, m.ambient_temperature);
    const PropertyValue p_sat = pair->saturation_pressure(e.content, m.ambient_temperature);
    e.viscosity = mu.value;
    e.saturation_pressure = p_sat.value;
    e.fluid_in_range = content.in_range && mu.in_range && p_sat.in_range;
    e.pressure_margin = m.ambient_pressure - e.saturation_pressure;
    if (m.load <= 0.0) {
        e.note = "Enter the load for S and the onset check.";
        return e;
    }
    e.sommerfeld = sommerfeld_number(e.viscosity, m.speed, m.bearing_length, m.bearing_diameter,
                                     m.radial_clearance, m.load);
    e.margin_ratio = e.pressure_margin / e.specific_load;
    const OnsetCurve& curve = tester_onset_curve();
    if (!curve.covers(e.sommerfeld)) {
        e.note = "S is outside the onset chart (0.02 to 3); run the case for the direct check.";
        e.onset_phi = std::numeric_limits<double>::quiet_NaN();
        // Below the chart's lowest S the film pressure drop only grows, so a
        // margin below the last Phi still means gas release.
        if (e.sommerfeld < curve.points.front().sommerfeld &&
            e.margin_ratio < curve.points.front().pressure_drop_ratio)
            e.verdict = OnsetVerdict::Expected;
        return e;
    }
    e.onset_phi = curve.phi(e.sommerfeld);
    e.verdict = onset_verdict(e.margin_ratio, e.onset_phi);
    return e;
}

}  // namespace gui
