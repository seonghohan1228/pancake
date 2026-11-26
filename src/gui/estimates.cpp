#include "estimates.hpp"

#include <algorithm>
#include <cmath>
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

double estimate_run_time(const CaseModel& m) {
    int n_theta = m.n_theta, n_z = m.n_z;
    if (m.mesh_level != MeshLevel::Custom) mesh_level_cells(m.mesh_level, n_theta, n_z);
    const double cells = static_cast<double>(n_theta) * n_z;
    // Calibration: compressor B0 mean (720 x 192 FBNS, load balance, 6 outer
    // iterations of ~4 film solves) took 101 s on 8 ranks of this workstation.
    const double seconds_per_cell_solve = 101.0 * 8.0 / (720.0 * 192.0 * 24.0);
    double solves = 24.0;  // steady load balance
    if (m.analysis == AnalysisType::Transient)
        solves = m.revolutions * m.steps_per_revolution * 3.0;  // ~3 coupled solves per step
    double factor = 1.0;
    if (m.thermal) factor *= 8.0;
    if (m.dissolved_gas && m.gas_model == GasModel::Transport) factor *= 10.0;
    const double parallel = std::max(1, m.cores) * 0.85 + 0.15;
    return cells * solves * seconds_per_cell_solve * factor / parallel;
}

Estimates compute_estimates(const CaseModel& m, const FluidPair* pair) {
    Estimates e;
    e.run_time_s = estimate_run_time(m);
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
