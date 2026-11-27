// Case report: one self-contained, printable HTML page (inline CSS and SVG, no
// external files) with the inputs, solver form, estimates, key results, design
// checks and figures drawn from the data of the run shown in Results.

#include "gui_app.hpp"

#include "result_data.hpp"
#include "win_util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string esc(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (char ch : text) {
        switch (ch) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += ch;
        }
    }
    return out;
}

std::string num(double value, const char* format = "%.4g") {
    if (!std::isfinite(value)) return "-";
    char text[64];
    std::snprintf(text, sizeof(text), format, value);
    return text;
}

/// Value in the engineering default unit of its family (first option), e.g. "20 mm".
std::string eng(double si, gui::UnitFamily family, const char* format = "%.4g") {
    const gui::UnitOption& unit = gui::unit_options(family).front();
    std::string text = num(gui::to_display(unit, si), format);
    if (std::isfinite(si) && unit.label[0]) text += std::string(" ") + unit.label;
    return text;
}

std::string hex(const ImVec4& c) {
    char text[16];
    std::snprintf(text, sizeof(text), "#%02x%02x%02x", static_cast<int>(c.x * 255.0f + 0.5f),
                  static_cast<int>(c.y * 255.0f + 0.5f), static_cast<int>(c.z * 255.0f + 0.5f));
    return text;
}

std::string risk_hex(gui::Risk risk) { return hex(ui::risk_color(risk)); }

std::string verdict_hex(gui::OnsetVerdict verdict) {
    switch (verdict) {
        case gui::OnsetVerdict::Expected: return hex(theme::kRiskHigh);
        case gui::OnsetVerdict::Borderline: return hex(theme::kRiskMedium);
        case gui::OnsetVerdict::NotExpected: return hex(theme::kRiskLow);
        case gui::OnsetVerdict::Unknown: break;
    }
    return hex(theme::kMuted);
}

/// Viridis-like palette, t in [0, 1].
std::string viridis(double t) {
    static const unsigned char anchors[9][3] = {{68, 1, 84},    {72, 40, 120},  {62, 73, 137},
                                                {49, 104, 142}, {38, 130, 142}, {31, 158, 137},
                                                {53, 183, 121}, {110, 206, 88}, {253, 231, 37}};
    t = std::clamp(std::isfinite(t) ? t : 0.0, 0.0, 1.0) * 8.0;
    const int i = std::min(static_cast<int>(t), 7);
    const double f = t - i;
    char text[16];
    std::snprintf(text, sizeof(text), "#%02x%02x%02x",
                  static_cast<int>(anchors[i][0] + f * (anchors[i + 1][0] - anchors[i][0])),
                  static_cast<int>(anchors[i][1] + f * (anchors[i + 1][1] - anchors[i][1])),
                  static_cast<int>(anchors[i][2] + f * (anchors[i + 1][2] - anchors[i][2])));
    return text;
}

double nice_step(double range) {
    if (!(range > 0.0)) return 1.0;
    const double raw = range / 5.0;
    const double base = std::pow(10.0, std::floor(std::log10(raw)));
    const double ratio = raw / base;
    return base * (ratio < 1.5 ? 1.0 : ratio < 3.5 ? 2.0 : ratio < 7.5 ? 5.0 : 10.0);
}

struct Series {
    std::vector<double> y;
    std::string colour;
    std::string label;
};

/// Line plot with axes and ticks as inline SVG; x over [x0, x0 + 360].
std::string svg_profile(const std::vector<double>& x, const std::vector<Series>& series, double x0,
                        const std::string& x_label, const std::string& y_label) {
    const double width = 700, height = 230, left = 64, right = 16, top = 14, bottom = 42;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const Series& s : series)
        for (double v : s.y)
            if (std::isfinite(v)) lo = std::min(lo, v), hi = std::max(hi, v);
    if (!std::isfinite(lo)) return "<p class=\"muted\">No data.</p>";
    if (hi <= lo) hi = lo + std::max(std::abs(lo) * 0.1, 1e-9);
    const double step = nice_step(hi - lo);
    lo = std::floor(lo / step) * step;
    hi = std::ceil(hi / step) * step;
    const auto px = [&](double v) { return left + (v - x0) / 360.0 * (width - left - right); };
    const auto py = [&](double v) { return top + (hi - v) / (hi - lo) * (height - top - bottom); };
    std::ostringstream svg;
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height
        << "\" class=\"fig\"><rect x=\"" << left << "\" y=\"" << top << "\" width=\"" << width - left - right
        << "\" height=\"" << height - top - bottom << "\" fill=\"none\" stroke=\"#444\"/>";
    for (double v = lo; v <= hi + step * 0.01; v += step)
        svg << "<line x1=\"" << left - 4 << "\" x2=\"" << left << "\" y1=\"" << py(v) << "\" y2=\"" << py(v)
            << "\" stroke=\"#444\"/><text x=\"" << left - 6 << "\" y=\"" << py(v) + 4
            << "\" text-anchor=\"end\">" << num(v, "%.4g") << "</text>";
    for (double v = x0; v <= x0 + 360.0 + 0.01; v += 60.0)
        svg << "<line x1=\"" << px(v) << "\" x2=\"" << px(v) << "\" y1=\"" << height - bottom << "\" y2=\""
            << height - bottom + 4 << "\" stroke=\"#444\"/><text x=\"" << px(v) << "\" y=\"" << height - bottom + 16
            << "\" text-anchor=\"middle\">" << num(v, "%.0f") << "</text>";
    svg << "<text x=\"" << (left + width - right) / 2 << "\" y=\"" << height - 6 << "\" text-anchor=\"middle\">"
        << esc(x_label) << "</text><text transform=\"translate(14," << (top + height - bottom) / 2
        << ") rotate(-90)\" text-anchor=\"middle\">" << esc(y_label) << "</text>";
    double legend_x = left + 8;
    for (const Series& s : series) {
        svg << "<polyline fill=\"none\" stroke=\"" << s.colour << "\" stroke-width=\"1.8\" points=\"";
        for (std::size_t i = 0; i < x.size() && i < s.y.size(); ++i)
            if (std::isfinite(s.y[i])) svg << px(x[i]) << ',' << py(s.y[i]) << ' ';
        svg << "\"/>";
        if (series.size() > 1) {
            svg << "<line x1=\"" << legend_x << "\" x2=\"" << legend_x + 18 << "\" y1=\"" << top + 12 << "\" y2=\""
                << top + 12 << "\" stroke=\"" << s.colour << "\" stroke-width=\"2\"/><text x=\"" << legend_x + 22
                << "\" y=\"" << top + 16 << "\">" << esc(s.label) << "</text>";
            legend_x += 30 + 7.0 * s.label.size();
        }
    }
    svg << "</svg>";
    return svg.str();
}

/// Unwrapped field map (block-averaged to at most 180 x 48 cells) as inline SVG.
/// `field` is row-major like the grid fields, in SI; `unit` converts it for display.
std::string svg_map(const ResultGrid& grid, const std::vector<double>& field, double zero, double start,
                    double aspect, const gui::UnitOption& unit, const std::string& x_label) {
    const auto* p = &field;
    if (grid.nx <= 0 || grid.ny <= 0 ||
        p->size() < static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny))
        return "<p class=\"muted\">Field not in the output.</p>";
    const auto order = result_data::angle_order(grid, zero, start);
    const int cols = std::min(grid.nx, 180), rows = std::min(grid.ny, 48);
    std::vector<double> cells(static_cast<std::size_t>(cols) * rows, 0.0);
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (int r = 0; r < rows; ++r) {
        const int j0 = r * grid.ny / rows, j1 = std::max(j0 + 1, (r + 1) * grid.ny / rows);
        for (int b = 0; b < cols; ++b) {
            const int k0 = b * grid.nx / cols, k1 = std::max(k0 + 1, (b + 1) * grid.nx / cols);
            double sum = 0.0;
            int count = 0;
            for (int j = j0; j < j1; ++j)
                for (int k = k0; k < k1; ++k) {
                    const double v = (*p)[static_cast<std::size_t>(j) * grid.nx + static_cast<std::size_t>(order.columns[k])];
                    if (std::isfinite(v)) sum += v, ++count;
                }
            const double v = count ? gui::to_display(unit, sum / count) : std::nan("");
            cells[static_cast<std::size_t>(r) * cols + b] = v;
            if (std::isfinite(v)) lo = std::min(lo, v), hi = std::max(hi, v);
        }
    }
    if (!std::isfinite(lo)) return "<p class=\"muted\">No values.</p>";
    if (hi <= lo) hi = lo + 1e-12;
    const double map_w = 600, map_h = std::clamp(600.0 * aspect, 50.0, 360.0);
    const double left = 52, top = 10, bar_x = left + map_w + 18;
    const double width = bar_x + 80, height = top + map_h + 44;
    const double cw = map_w / cols, ch = map_h / rows;
    std::ostringstream svg;
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << ' ' << height
        << "\" class=\"fig\" shape-rendering=\"crispEdges\">";
    for (int r = 0; r < rows; ++r)
        for (int b = 0; b < cols; ++b) {
            const double v = cells[static_cast<std::size_t>(r) * cols + b];
            // Row r = 0 is z = 0 (bottom of the map).
            svg << "<rect x=\"" << num(left + b * cw, "%.2f") << "\" y=\"" << num(top + map_h - (r + 1) * ch, "%.2f")
                << "\" width=\"" << num(cw + 0.05, "%.2f") << "\" height=\"" << num(ch + 0.05, "%.2f") << "\" fill=\""
                << viridis((v - lo) / (hi - lo)) << "\"/>";
        }
    svg << "<g shape-rendering=\"auto\"><rect x=\"" << left << "\" y=\"" << top << "\" width=\"" << map_w
        << "\" height=\"" << map_h << "\" fill=\"none\" stroke=\"#444\"/>";
    for (double v = start; v <= start + 360.0 + 0.01; v += 60.0) {
        const double x = left + (v - start) / 360.0 * map_w;
        svg << "<text x=\"" << x << "\" y=\"" << top + map_h + 15 << "\" text-anchor=\"middle\">" << num(v, "%.0f")
            << "</text>";
    }
    for (double v : {0.0, 0.5, 1.0})
        svg << "<text x=\"" << left - 6 << "\" y=\"" << top + map_h * (1.0 - v) + 4 << "\" text-anchor=\"end\">"
            << num(v, "%.1f") << "</text>";
    svg << "<text x=\"" << left + map_w / 2 << "\" y=\"" << height - 6 << "\" text-anchor=\"middle\">" << esc(x_label)
        << "</text><text transform=\"translate(12," << top + map_h / 2
        << ") rotate(-90)\" text-anchor=\"middle\">z / L</text>";
    for (int i = 0; i < 32; ++i)
        svg << "<rect x=\"" << bar_x << "\" y=\"" << num(top + map_h * (31 - i) / 32.0, "%.2f") << "\" width=\"14\" height=\""
            << num(map_h / 32.0 + 0.1, "%.2f") << "\" fill=\"" << viridis(i / 31.0) << "\"/>";
    svg << "<text x=\"" << bar_x + 18 << "\" y=\"" << top + 10 << "\">" << num(hi) << "</text><text x=\"" << bar_x + 18
        << "\" y=\"" << top + map_h << "\">" << num(lo) << "</text><text x=\"" << bar_x + 18 << "\" y=\""
        << top + map_h / 2 + 4 << "\">" << esc(unit.label) << "</text></g></svg>";
    return svg.str();
}

struct Table {
    std::ostringstream html;
    explicit Table(const char* title) { html << "<h2>" << title << "</h2><table>"; }
    void row(const std::string& name, const std::string& value, const std::string& colour = {}) {
        html << "<tr><th>" << esc(name) << "</th><td" << (colour.empty() ? "" : " style=\"color:" + colour + "\"")
             << ">" << esc(value) << "</td></tr>";
    }
    std::string str() { return html.str() + "</table>"; }
};

bool contains_word(std::string text, const char* word) {
    for (char& ch : text) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return text.find(word) != std::string::npos;
}

}  // namespace

bool GuiApp::export_case_report(const fs::path& html) {
    const RunRecord* record = shown_record_ >= 0 && shown_record_ < static_cast<int>(records_.size())
                                  ? &records_[static_cast<std::size_t>(shown_record_)] : nullptr;
    const gui::CaseModel& m = record ? record->job.model : case_;
    const SimulationConfig& cfg = record ? record->job.config : config_;
    const fs::path case_file = record ? record->job.config_file : case_path_;
    const std::string command = record ? record->command : command_line_for(case_path_, std::max(1, case_.cores));
    const gui::FluidPair* fluid = pair() && pair_.id == m.fluid_pair ? pair() : nullptr;
    const gui::Estimates est = gui::compute_estimates(m, fluid);
    const bool latest = selected_step_ < 0 || selected_step_ >= static_cast<int>(steps_.size()) - 1;
    gui::KeyResults results;
    bool have_results = false;
    if (record && record->results_ready && latest) results = record->results, have_results = true;
    else if (grid_) results = gui::compute_key_results(*grid_, cfg), have_results = true;
    const double p_sat = record ? record->job.saturation_pressure : est.saturation_pressure;
    const bool transient = m.analysis == gui::AnalysisType::Transient;

    SYSTEMTIME now{};
    GetLocalTime(&now);
    char date[32];
    std::snprintf(date, sizeof(date), "%04u-%02u-%02u %02u:%02u", now.wYear, now.wMonth, now.wDay, now.wHour, now.wMinute);
    const std::string title = case_file.empty() ? std::string("Unsaved case") : winutil::path_to_utf8(case_file.stem());

    std::ostringstream out;
    out << "<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>" << esc(title)
        << " - bearing case report</title><style>"
           "body{font-family:Segoe UI,Arial,sans-serif;font-size:10.5pt;color:#222;max-width:820px;margin:18px auto;padding:0 14px}"
           "h1{font-size:16pt;margin:0 0 2px}h2{font-size:11.5pt;margin:14px 0 4px;border-bottom:1px solid #bbb}"
           ".meta{color:#555;margin:0}.mono{font-family:Consolas,monospace;font-size:9pt;word-break:break-all}"
           ".cols{display:flex;gap:18px;flex-wrap:wrap}.cols>div{flex:1 1 380px}"
           "table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:2px 6px;border-bottom:1px solid #e3e3e3;"
           "vertical-align:top}th{font-weight:normal;color:#444;width:52%}td{font-weight:600}"
           ".muted{color:#777}.fig{width:100%;height:auto;font-size:10px;font-family:Segoe UI,Arial,sans-serif}"
           ".note{font-size:9pt;color:#555;margin-top:14px;border-top:1px solid #bbb;padding-top:6px}"
           "@media print{body{margin:0;max-width:none}h2{break-after:avoid}.fig{break-inside:avoid}}"
           "</style></head><body>";
    out << "<h1>" << esc(title) << "</h1><p class=\"meta\">Journal bearing case report, " << date << "</p>";
    out << "<p class=\"meta\">Case file: " << esc(case_file.empty() ? "(not saved)" : winutil::path_to_utf8(case_file))
        << "</p><p class=\"meta\">Command: <span class=\"mono\">" << esc(command) << "</span></p>";

    // Inputs
    Table geometry("Geometry");
    geometry.row("Bearing diameter D", eng(m.bearing_diameter, gui::UnitFamily::Length));
    geometry.row("Bearing length L", eng(m.bearing_length, gui::UnitFamily::Length));
    geometry.row("Radial clearance c", eng(m.radial_clearance, gui::UnitFamily::SmallLength));
    geometry.row("Oil feed", std::string(m.feed.type == gui::FeedType::Hole ? "hole" : "groove") + " at " +
                                 num(m.feed.angle_deg, "%.1f") + " deg (fixed frame)");
    if (m.feed.type == gui::FeedType::Hole) geometry.row("Hole diameter", eng(m.feed.hole_diameter, gui::UnitFamily::Length));
    else geometry.row("Groove width", num(m.feed.groove_width_deg, "%.1f") + " deg");
    geometry.row("Supply pressure above ambient", eng(m.feed.supply_overpressure, gui::UnitFamily::Pressure));
    if (m.tilt_angle != 0.0)
        geometry.row("Tilt angle, direction", eng(m.tilt_angle, gui::UnitFamily::TiltAngle) + ", " +
                                                  num(m.tilt_direction_deg, "%.1f") + " deg");

    Table fluid_table("Fluid");
    if (fluid) {
        fluid_table.row("Oil-refrigerant pair", fluid->name.empty() ? fluid->id : fluid->name);
        fluid_table.row("Oil / refrigerant", fluid->oil + " / " + fluid->refrigerant);
        fluid_table.row("Source", fluid->source.empty() ? "-" : fluid->source);
        if (!fluid->note.empty()) fluid_table.row("Note", fluid->note);
    } else {
        fluid_table.row("Oil-refrigerant pair", m.fluid_pair.empty() ? "(none)" : m.fluid_pair + " (not loaded)");
    }
    fluid_table.row("Refrigerant content", m.saturated ? "saturated at ambient pressure (" +
                                                             eng(est.content, gui::UnitFamily::Content, "%.3g") + ")"
                                                       : eng(m.content, gui::UnitFamily::Content, "%.3g"));
    if (m.roughness_shaft > 0.0 || m.roughness_bush > 0.0)
        fluid_table.row("Surface roughness Rq shaft / bush", eng(m.roughness_shaft, gui::UnitFamily::SmallLength) + " / " +
                                                                 eng(m.roughness_bush, gui::UnitFamily::SmallLength));

    Table operating("Operating conditions");
    operating.row("Analysis", transient ? "transient" : "steady");
    operating.row("Speed", eng(m.speed, gui::UnitFamily::Speed));
    if (transient) {
        operating.row("Load table", winutil::path_to_utf8(m.load_table.file.filename()));
        operating.row("Moving mass", eng(m.moving_mass, gui::UnitFamily::Mass));
        operating.row("Revolutions", num(m.revolutions, "%.3g"));
    } else if (m.position == gui::ShaftPosition::LoadBalance) {
        operating.row("Load", eng(m.load, gui::UnitFamily::Force) + " toward " + num(m.load_direction_deg, "%.0f") + " deg");
        operating.row("Specific load P = W/(LD)", eng(est.specific_load, gui::UnitFamily::Pressure));
    } else {
        operating.row("Fixed shaft position", "eccentricity ratio " + num(m.fixed_eccentricity_ratio, "%.3f") +
                                                  ", attitude " + num(m.fixed_attitude_deg, "%.1f") + " deg");
    }
    operating.row("Ambient pressure", eng(m.ambient_pressure, gui::UnitFamily::Pressure));
    operating.row("Ambient temperature", eng(m.ambient_temperature, gui::UnitFamily::Temperature, "%.1f"));

    Table model_table("Model and mesh");
    model_table.row("Dissolved-gas release", m.dissolved_gas ? "on (GC+VC)" : "off (VC only)");
    if (m.dissolved_gas)
        model_table.row("Gas model", m.gas_model == gui::GasModel::Plateau ? "pressure plateau at p_sat" : "gas transport");
    model_table.row("Thermal analysis", m.thermal ? "on" : "off (isothermal)");
    model_table.row("Solver form", gui::solver_form_text(m, fluid));
    model_table.row("Mesh", std::string(gui::mesh_level_name(m.mesh_level)) + ", " + std::to_string(m.n_theta) +
                                " x " + std::to_string(m.n_z) + " cells");
    if (transient) model_table.row("Time steps per revolution", std::to_string(m.steps_per_revolution));

    Table estimates("Estimates (before the run)");
    estimates.row("Viscosity at ambient state", eng(est.viscosity, gui::UnitFamily::Viscosity));
    estimates.row("Sommerfeld number S", num(est.sommerfeld, "%.4g"));
    estimates.row("Saturation pressure p_sat", eng(est.saturation_pressure, gui::UnitFamily::Pressure));
    estimates.row("Pressure margin p_amb - p_sat", eng(est.pressure_margin, gui::UnitFamily::Pressure));
    estimates.row("Margin / specific load", num(est.margin_ratio, "%.3g"));
    estimates.row("Onset chart Phi(S)", num(est.onset_phi, "%.3g"));
    estimates.row("Gaseous cavitation", gui::onset_verdict_text(est.verdict), verdict_hex(est.verdict));
    if (!est.fluid_in_range)
        estimates.row("Fluid table", "case state outside the table (clamped)", hex(theme::kRiskMedium));

    out << "<div class=\"cols\"><div>" << geometry.str() << fluid_table.str() << "</div><div>" << operating.str()
        << model_table.str() << estimates.str() << "</div></div>";

    // Results
    ui::UnitChoices& unit_map = units();
    const auto rv = [&](double si, gui::UnitFamily family, const char* id) {
        return std::isfinite(si) ? ui::format_value(si, family, unit_map, id) : std::string("-");
    };
    Table key("Key results");
    if (record) {
        const bool failed = record->outcome == SolverState::Failed;
        key.row("Run", std::string(to_string(record->outcome)) + (record->converged ? ", converged" : ", not converged"),
                failed || !record->converged ? hex(theme::kRiskHigh) : std::string());
    }
    if (have_results) {
        key.row("Minimum film thickness h_min", rv(results.h_min, gui::UnitFamily::SmallLength, result_data::kUnitFilm));
        key.row("Eccentricity ratio", num(results.eccentricity_ratio, "%.3f"));
        key.row("Attitude angle", num(results.attitude_deg, "%.1f") + " deg");
        key.row("Friction torque", rv(results.friction_torque, gui::UnitFamily::Torque, result_data::kUnitTorque));
        key.row("Friction power loss", rv(results.friction_power, gui::UnitFamily::Power, result_data::kUnitPower));
        key.row("Peak pressure", rv(results.p_max, gui::UnitFamily::Pressure, result_data::kUnitPressure));
        key.row("Cavitated area (void fraction > 0)", num(results.cavitated_fraction * 100.0, "%.1f") + " %");
        key.row("Shear share of film force", num(results.shear_share * 100.0, "%.1f") + " %");
        if (std::isfinite(results.lambda)) key.row("Film thickness ratio Lambda", num(results.lambda, "%.2f"));
        if (transient) key.row("Time of the step shown", num(results.time, "%.6g") + " s");
        const gui::DesignChecks checks = gui::design_checks(results, p_sat, m.dissolved_gas);
        if (!checks.gas_release_text.empty())
            key.row("Gas release check", checks.gas_release_text, risk_hex(checks.gas_release_risk));
        if (!checks.lambda_text.empty()) key.row("Film thickness check", checks.lambda_text, risk_hex(checks.lambda_risk));
    } else {
        key.row("Results", "no run results available");
    }
    out << key.str();

    if (grid_) {
        gui::KeyResults grid_results = gui::compute_key_results(*grid_, cfg);
        const double zero = result_data::zero_position(angle_zero_, m.feed.angle_deg, grid_results);
        const double start = gui::angle_axis_start(angle_zero_);
        const std::string x_label = gui::angle_axis_label(angle_zero_);
        const auto& pu = result_data::unit_of(unit_map, result_data::kUnitPressure, gui::UnitFamily::Pressure);
        const auto& hu = result_data::unit_of(unit_map, result_data::kUnitFilm, gui::UnitFamily::SmallLength);
        const double aspect = m.bearing_diameter > 0.0 ? m.bearing_length / (kPi * m.bearing_diameter) : 0.3;
        if (const auto* p = grid_->field("p"))
            out << "<h2>Film pressure, unwrapped (to scale, pi D : L)</h2>"
                << svg_map(*grid_, *p, zero, start, aspect, pu, x_label);
        // Cavitated region by model (plan 0.5): void fraction 1 - theta for the pressure
        // plateau, free-gas volume fraction for gas transport; never a gas field for FBNS.
        const gui::FieldView cavity = result_data::cavity_view(cfg);
        const std::vector<double> cavity_values = result_data::cavity_fraction(*grid_, cfg);
        if (!cavity_values.empty()) {
            static const gui::UnitOption fraction{"-", 1.0};
            out << "<h2>" << esc(cavity.label) << ", unwrapped</h2>"
                << svg_map(*grid_, cavity_values, zero, start, aspect, fraction, x_label);
            if (gui::result_form(cfg) == gui::ResultForm::Fbns)
                out << "<p class=\"muted\">Pressure-plateau model: the cavity is liquid streaks plus a massless "
                       "void; gas is not tracked.</p>";
        }
        const result_data::Profile profile = result_data::make_profile(*grid_, cfg, zero, start, 0.5);
        out << "<h2>Mid-plane pressure</h2>"
            << svg_profile(profile.angle, {{result_data::scaled(profile.p, pu), hex(theme::kAccent), "pressure"}}, start,
                           x_label, std::string("Pressure [") + pu.label + "]");
        out << "<h2>Mid-plane film thickness</h2>"
            << svg_profile(profile.angle, {{result_data::scaled(profile.h, hu), hex(theme::kAccent), "film"}}, start,
                           x_label, std::string("Film thickness [") + hu.label + "]");
    }

    // Notes
    out << "<div class=\"note\">";
    if (fluid && (contains_word(fluid->source, "synthetic") || contains_word(fluid->note, "example") ||
                  contains_word(fluid->note, "not for design")))
        out << "<p style=\"color:" << hex(theme::kRiskHigh) << "\"><b>Synthetic fluid data:</b> the pair "
            << esc(fluid->name.empty() ? fluid->id : fluid->name)
            << " is an example, not measured data. Do not use these results for design.</p>";
    const bool current_case = !case_path_.empty() && case_file == case_path_;
    if (current_case && !hidden_keys_.empty())
        out << "<p>The case file sets " << hidden_keys_.size()
            << " settings that the GUI does not show; they were used as written in the file.</p>";
    out << "<p>Estimates come from the onset chart of a tester bearing and are indicative; the run gives the direct "
           "check (film p_min against p_sat). Contact and wear are not modelled.</p></div></body></html>";

    std::ofstream file(html, std::ios::binary | std::ios::trunc);
    file << out.str();
    if (!file) {
        set_status("Report not written: cannot write " + winutil::path_to_utf8(html));
        return false;
    }
    set_status("Saved case report " + winutil::path_to_utf8(html));
    return true;
}
