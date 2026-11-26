// Results window: run picker, key-results card with the design checks, and the
// tabs Field maps | Polar | Profiles | Transient | Comparison, plus PNG, CSV and
// case-report export. Grids are parsed off the UI thread (ResultsLoader, and
// std::async for motion.csv and the two comparison runs).

#include "gui_app.hpp"

#include "exports.hpp"
#include "imgui.h"
#include "implot.h"
#include "result_data.hpp"
#include "win_util.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <future>

namespace fs = std::filesystem;
using result_data::CsvTable;
using result_data::MapCache;
using result_data::MapField;
using result_data::Profile;
using result_data::build_map;
using result_data::make_profile;
using result_data::scaled;
using result_data::MotionCache;
using result_data::CompareRun;
using result_data::CompareResult;
using result_data::ensure_motion;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr const char* kGasOnLabel = "gas on (GC+VC)";
constexpr const char* kGasOffLabel = "gas off (VC only)";



struct ResultsUi {
    MapField field = MapField::Pressure;
    float profile_z = 0.5f;
    int transient_x = 0;  // 0 time, 1 crank angle
    double cursor = 0.0;
    bool playing = false;
    double last_tick = 0.0;
    ImVec2 plot_min{0, 0}, plot_max{0, 0};
    std::string plot_name = "plot";
    fs::path pending_scan;
    result_data::ReloadGuard reload;
    const ResultGrid* results_grid = nullptr;
    int results_record = -2;
    gui::KeyResults grid_results;
    int estimate_record = -2;
    bool estimate_in_range = true;
    MapCache map;
    MotionCache motion;
    CompareRun compare[2];
};

ResultsUi& state() {
    static ResultsUi s;
    return s;
}

std::string value_text(double si, gui::UnitFamily family, ui::UnitChoices& units, const char* id,
                       const char* format = "%.4g") {
    return std::isfinite(si) ? ui::format_value(si, family, units, id, format) : std::string("-");
}

std::string plain_text(double value, const char* format) {
    if (!std::isfinite(value)) return "-";
    char text[64];
    std::snprintf(text, sizeof(text), format, value);
    return text;
}

void remember_plot(ResultsUi& s, const char* name) {
    s.plot_min = ImGui::GetItemRectMin();
    s.plot_max = ImGui::GetItemRectMax();
    s.plot_name = name;
}

std::string record_label(const RunRecord& r, int index) {
    char text[320];
    std::snprintf(text, sizeof(text), "#%d  %s  [%s]  %s  %s", index + 1, r.job.label.c_str(),
                  r.job.group.c_str(), r.started_at.c_str(), to_string(r.outcome));
    return text;
}

gui::Light risk_light(gui::Risk risk) {
    switch (risk) {
        case gui::Risk::Low: return gui::Light::Green;
        case gui::Risk::Medium: return gui::Light::Yellow;
        case gui::Risk::High: return gui::Light::Red;
        case gui::Risk::None: break;
    }
    return gui::Light::Idle;
}

/// Coloured status mark with one sentence.
void mark(gui::Light light, const std::string& text) {
    ui::status_dot(light);
    ImGui::SameLine();
    ImGui::PushStyleColor(ImGuiCol_Text, ui::light_color(light));
    ImGui::TextWrapped("%s", text.c_str());
    ImGui::PopStyleColor();
}

void unit_choice_row(const char* label, const char* id, gui::UnitFamily family,
                     ui::UnitChoices& choices) {
    const auto& options = gui::unit_options(family);
    result_data::unit_of(choices, id, family);
    int& choice = choices[id];
    ui::form_label(label, nullptr);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##unit", options[static_cast<std::size_t>(choice)].label)) {
        for (int i = 0; i < static_cast<int>(options.size()); ++i)
            if (ImGui::Selectable(options[static_cast<std::size_t>(i)].label, i == choice)) choice = i;
        ImGui::EndCombo();
    }
    ImGui::PopID();
}


void draw_field_map(ResultsUi& s, const MapCache& m, const std::string& bar_label,
                    const std::string& x_label, double start, double aspect, double feed_x,
                    double feed_y, double hmin_x) {
    const float font = ImGui::GetFontSize();
    const float bar_w = font * 6.5f;
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float pad_x = font * 4.0f, pad_y = font * 3.6f;
    float inner_w = std::max(avail.x - bar_w - ImGui::GetStyle().ItemSpacing.x - pad_x, 120.0f);
    float inner_h = inner_w * static_cast<float>(aspect);
    const float max_h = std::max(avail.y - ImGui::GetTextLineHeightWithSpacing() * 2.5f - pad_y, 80.0f);
    if (inner_h > max_h) {
        inner_h = max_h;
        inner_w = inner_h / static_cast<float>(aspect);
    }
    const ImVec2 size(inner_w + pad_x, inner_h + pad_y);
    ImGui::BeginGroup();
    ImPlot::PushColormap(ImPlotColormap_Viridis);
    if (ImPlot::BeginPlot("##field_map", size, ImPlotFlags_NoLegend | ImPlotFlags_NoMenus)) {
        ImPlot::SetupAxes(x_label.c_str(), "z / L", ImPlotAxisFlags_NoGridLines,
                          ImPlotAxisFlags_NoGridLines);
        ImPlot::SetupAxisLimits(ImAxis_X1, start, start + 360.0, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
        ImPlot::PlotHeatmap("##values", m.values.data(), m.rows, m.cols, m.vmin, m.vmax, nullptr,
                            ImPlotPoint(m.x_left, 0.0), ImPlotPoint(m.x_left + 360.0, 1.0));
        ImDrawList* draw = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        const ImU32 outline = ImGui::GetColorU32(theme::kHighlight);
        for (std::size_t i = 0; i + 3 < m.outline.size(); i += 4)
            draw->AddLine(ImPlot::PlotToPixels(m.outline[i], m.outline[i + 1]),
                          ImPlot::PlotToPixels(m.outline[i + 2], m.outline[i + 3]), outline, 2.0f);
        ImPlot::PopPlotClipRect();
        ImPlot::PushStyleColor(ImPlotCol_InlayText, ImVec4(1, 1, 1, 1));
        if (std::isfinite(hmin_x)) {
            ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.9f), 1.5f);
            ImPlot::PlotInfLines("##hmin", &hmin_x, 1);
            ImPlot::PlotText("h_min", hmin_x, 0.85, ImVec2(ImGui::GetFontSize() * 1.6f, 0));
        }
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Down, 7.0f, ImVec4(1, 1, 1, 1), 1.0f, ImVec4(0, 0, 0, 1));
        ImPlot::PlotScatter("##feed", &feed_x, &feed_y, 1);
        ImPlot::PlotText("feed", feed_x, feed_y, ImVec2(0, -ImGui::GetFontSize()));
        ImPlot::PopStyleColor();
        ImPlot::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale(bar_label.c_str(), m.vmin, m.vmax, ImVec2(bar_w, size.y), "%.3g");
    ImPlot::PopColormap();
    ImGui::EndGroup();
    remember_plot(s, "field_map");
    ImGui::TextDisabled("Purple outline: cavitated region (gas fraction 0.5). White line: minimum film. "
                        "Width and height to scale (pi D : L).");
}

void draw_polar(ResultsUi& s, const ResultGrid& grid, const gui::CaseModel& model, double zero,
                double start, const gui::KeyResults& r, double p_sat, ui::UnitChoices& units) {
    const auto* p = grid.field("p");
    if (!p) {
        ImGui::TextDisabled("This output has no pressure field.");
        return;
    }
    const auto order = result_data::angle_order(grid, zero, start);
    const auto values = result_data::row_values(grid, *p, grid.ny / 2, order);
    const double p_amb = model.ambient_pressure;
    double hi = 0.0, lo = 0.0;
    for (double v : values)
        if (std::isfinite(v)) hi = std::max(hi, v - p_amb), lo = std::min(lo, v - p_amb);
    if (p_sat > 0.0) hi = std::max(hi, p_sat - p_amb), lo = std::min(lo, p_sat - p_amb);

    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float side = std::max(220.0f, std::min(avail.x, avail.y - ImGui::GetTextLineHeightWithSpacing() * 3.5f));
    ImGui::BeginGroup();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##polar", ImVec2(side, side));
    const bool hovered = ImGui::IsItemHovered();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 centre(origin.x + side * 0.5f, origin.y + side * 0.5f);
    const float r_max = side * 0.40f, r_0 = r_max * 0.5f, r_min = r_max * 0.08f;
    double gain = hi > 0.0 ? (r_max - r_0) / hi : 0.0;
    if (lo < 0.0) {
        const double inward = (r_0 - r_min) / -lo;
        gain = gain > 0.0 ? std::min(gain, inward) : inward;
    }
    const auto at = [&](double display_deg, double radius) {
        const double a = display_deg * kPi / 180.0;
        return ImVec2(centre.x + static_cast<float>(radius * std::cos(a)),
                      centre.y - static_cast<float>(radius * std::sin(a)));
    };
    const ImU32 grid_colour = ImGui::GetColorU32(ImGuiCol_Border);
    const ImU32 text_colour = ImGui::GetColorU32(ImGuiCol_Text);
    draw->AddCircleFilled(centre, r_max * 1.05f, ImGui::GetColorU32(ImGuiCol_FrameBg), 128);
    for (int q = 0; q < 12; ++q) {
        const double angle = start + 30.0 * q;
        draw->AddLine(centre, at(angle, r_max * 1.05), grid_colour);
        if (q % 3 != 0) continue;
        char label[16];
        std::snprintf(label, sizeof(label), "%.0f", angle);
        const ImVec2 size = ImGui::CalcTextSize(label);
        const ImVec2 pos = at(angle, r_max * 1.16);
        draw->AddText(ImVec2(pos.x - size.x * 0.5f, pos.y - size.y * 0.5f), text_colour, label);
    }
    draw->AddCircle(centre, r_0, ImGui::GetColorU32(theme::kMuted), 128, 1.5f);
    if (p_sat > 0.0)
        draw->AddCircle(centre, static_cast<float>(r_0 + gain * (p_sat - p_amb)),
                        ImGui::GetColorU32(theme::kHighlight), 128, 1.5f);
    std::vector<ImVec2> points;
    points.reserve(values.size());
    for (std::size_t k = 0; k < values.size(); ++k)
        if (std::isfinite(values[k])) points.push_back(at(order.angles[k], r_0 + gain * (values[k] - p_amb)));
    if (points.size() > 2)
        draw->AddPolyline(points.data(), static_cast<int>(points.size()), ImGui::GetColorU32(theme::kAccent),
                          ImDrawFlags_Closed, 2.0f);
    const auto marker = [&](double fixed_deg, const char* label) {
        const double display = gui::display_angle(fixed_deg, zero, start);
        const ImVec2 pos = at(display, r_max * 1.05);
        draw->AddCircleFilled(pos, 4.5f, text_colour);
        draw->AddText(ImVec2(pos.x + 6.0f, pos.y - ImGui::GetFontSize()), text_colour, label);
    };
    marker(model.feed.angle_deg, "feed");
    if (std::isfinite(r.h_min_theta_deg)) marker(r.h_min_theta_deg, "h_min");
    // Direction of rotation: short arc with an arrow head (angles grow counter-clockwise).
    {
        std::vector<ImVec2> arc;
        for (int i = 0; i <= 16; ++i) arc.push_back(at(start + 30.0 + 2.5 * i, r_max * 0.3));
        const ImU32 muted = ImGui::GetColorU32(theme::kMuted);
        draw->AddPolyline(arc.data(), static_cast<int>(arc.size()), muted, 0, 2.0f);
        const ImVec2 tip = arc.back();
        const ImVec2 back = at(start + 60.0, r_max * 0.3);
        const float length = std::max(std::hypot(tip.x - back.x, tip.y - back.y), 1.0f);
        const float head = ImGui::GetFontSize() * 0.45f;
        const ImVec2 dir((tip.x - back.x) / length * head, (tip.y - back.y) / length * head);
        draw->AddTriangleFilled(ImVec2(tip.x + dir.x, tip.y + dir.y), ImVec2(tip.x - dir.y, tip.y + dir.x),
                                ImVec2(tip.x + dir.y, tip.y - dir.x), muted);
    }
    if (hovered && !values.empty()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        double a = std::atan2(-(mouse.y - centre.y), mouse.x - centre.x) * 180.0 / kPi;
        a = gui::display_angle(a, 0.0, start);
        std::size_t best = 0;
        for (std::size_t k = 1; k < order.angles.size(); ++k)
            if (std::abs(order.angles[k] - a) < std::abs(order.angles[best] - a)) best = k;
        ImGui::SetTooltip("%.1f deg: %s", order.angles[best],
                          value_text(values[best], gui::UnitFamily::Pressure, units,
                                     result_data::kUnitPressure).c_str());
    }
    ImGui::EndGroup();
    remember_plot(s, "polar_pressure");
    ImGui::TextWrapped("Mid-plane pressure (blue): radius grows with the pressure above ambient. Grey circle: "
                       "ambient pressure %s. Purple circle: saturation pressure %s.",
                       value_text(p_amb, gui::UnitFamily::Pressure, units, result_data::kUnitPressure).c_str(),
                       p_sat > 0.0 ? value_text(p_sat, gui::UnitFamily::Pressure, units,
                                                result_data::kUnitPressure).c_str()
                                   : "unknown");
    ImGui::TextDisabled("Angles grow counter-clockwise, in the direction of rotation.");
}

}  // namespace

void GuiApp::draw_results_window() {
    ResultsUi& s = state();
    s.reload.poll(loader_.busy(), loader_.revision(), steps_, selected_step_, grid_.get(),
                  [&](const StepInfo& info) { loader_.load_step(info); });
    if (!s.pending_scan.empty() && !loader_.busy()) {
        s.reload.arm(loader_.revision(), grid_.get());
        request_scan(s.pending_scan);
        s.pending_scan.clear();
    }
    if (!ImGui::Begin("Results")) {
        ImGui::End();
        return;
    }
    ui::UnitChoices& unit_map = units();
    const int record_count = static_cast<int>(records_.size());
    if (shown_record_ >= record_count) shown_record_ = -1;

    // --- Run picker, angle zero, units ---------------------------------------
    if (record_count > 0) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Run");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x * 0.5f, ImGui::GetFontSize() * 30.0f));
        const std::string preview = shown_record_ >= 0
            ? record_label(records_[static_cast<std::size_t>(shown_record_)], shown_record_)
            : std::string("Choose a run");
        if (ImGui::BeginCombo("##run", preview.c_str())) {
            for (int i = record_count - 1; i >= 0; --i) {
                ImGui::PushID(i);
                const RunRecord& r = records_[static_cast<std::size_t>(i)];
                if (ImGui::Selectable(record_label(r, i).c_str(), i == shown_record_) && i != shown_record_) {
                    shown_record_ = i;
                    selected_step_ = -1;
                    s.playing = false;
                    s.pending_scan = r.job.output_dir;
                    refresh_derived();
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        ui::item_tooltip("Runs of this session, newest first: label, group, start time and outcome.");
        ImGui::SameLine();
    }
    if (ImGui::Button("Open results folder...")) {
        const fs::path pvd = exports::open_file_dialog(
            hwnd, L"Result collection (results.pvd)\0*.pvd\0All files\0*.*\0");
        if (!pvd.empty()) {
            shown_record_ = -1;
            selected_step_ = -1;
            s.pending_scan = pvd.parent_path();
        }
    }
    ui::item_tooltip("Show the results of an earlier run: pick the results.pvd file in its results folder.");
    ImGui::SameLine();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Angle zero");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
    if (ImGui::BeginCombo("##angle_zero", gui::angle_zero_name(angle_zero_))) {
        for (int z = 0; z <= static_cast<int>(gui::AngleZero::OppositeFeed); ++z) {
            const auto zero = static_cast<gui::AngleZero>(z);
            if (ImGui::Selectable(gui::angle_zero_name(zero), zero == angle_zero_)) {
                angle_zero_ = zero;
                settings.angle_zero = z;
            }
        }
        ImGui::EndCombo();
    }
    ui::item_tooltip("Zero of every angle axis and of the polar plot. The maximum and minimum film zeros "
                     "follow the field shown.");
    ImGui::SameLine();
    if (ImGui::Button("Units")) ImGui::OpenPopup("result_units");
    ui::item_tooltip("Display units of the results.");
    if (ImGui::BeginPopup("result_units")) {
        if (ui::begin_form("units_form")) {
            unit_choice_row("Pressure", result_data::kUnitPressure, gui::UnitFamily::Pressure, unit_map);
            unit_choice_row("Film thickness", result_data::kUnitFilm, gui::UnitFamily::SmallLength, unit_map);
            unit_choice_row("Friction torque", result_data::kUnitTorque, gui::UnitFamily::Torque, unit_map);
            unit_choice_row("Friction power loss", result_data::kUnitPower, gui::UnitFamily::Power, unit_map);
            unit_choice_row("Temperature", result_data::kUnitTemperature, gui::UnitFamily::Temperature, unit_map);
            ui::end_form();
        }
        ImGui::EndPopup();
    }

    if (!grid_ && steps_.empty() && record_count == 0) {
        ImGui::Spacing();
        ImGui::TextDisabled("%s", results_error_.empty()
                                      ? "No results yet. Run the case (F5) or open a results folder."
                                      : results_error_.c_str());
        ImGui::End();
        return;
    }

    // --- Context of the shown run ---------------------------------------------
    const RunRecord* record = shown_record_ >= 0 ? &records_[static_cast<std::size_t>(shown_record_)] : nullptr;
    const SimulationConfig& cfg = record ? record->job.config : config_;
    const gui::CaseModel& model = record ? record->job.model : case_;
    if (grid_ && (s.results_grid != grid_.get() || s.results_record != shown_record_)) {
        s.results_grid = grid_.get();
        s.results_record = shown_record_;
        s.grid_results = gui::compute_key_results(*grid_, cfg);
    }
    if (s.estimate_record != shown_record_) {
        s.estimate_record = shown_record_;
        s.estimate_in_range = gui::compute_estimates(model, pair()).fluid_in_range;
    }
    const int last_step = static_cast<int>(steps_.size()) - 1;
    const int shown_step = selected_step_ < 0 ? last_step : std::min(selected_step_, last_step);
    const bool latest = shown_step == last_step;
    gui::KeyResults shown;
    bool have_results = false;
    if (record && record->results_ready && latest) {
        shown = record->results;
        have_results = true;
    } else if (grid_) {
        shown = s.grid_results;
        have_results = true;
    }
    const gui::KeyResults& grid_r = grid_ ? s.grid_results : shown;
    const double zero = result_data::zero_position(angle_zero_, model.feed.angle_deg, grid_r);
    const double start = gui::angle_axis_start(angle_zero_);
    const std::string x_label = gui::angle_axis_label(angle_zero_);
    const double p_sat = record ? record->job.saturation_pressure : estimates_.saturation_pressure;

    // --- Key-results card -------------------------------------------------------
    ui::section("Key results");
    if (record && record->outcome == SolverState::Failed)
        mark(gui::Light::Red, "The run failed (exit " + std::to_string(record->exit_code) +
                                  "): the solver log names the cause.");
    else if (record && record->outcome == SolverState::Cancelled)
        mark(gui::Light::Yellow, "The run was cancelled: the last saved step is shown.");
    else if (record && !record->converged)
        mark(gui::Light::Red, "The run did not converge: the results are not reliable.");
    if (!s.estimate_in_range)
        mark(gui::Light::Yellow, "The case state is outside the fluid table: properties were taken at the table edge.");
    if (record && !record->results_ready && !record->results_error.empty())
        ImGui::TextDisabled("Key results not read: %s", record->results_error.c_str());
    if (!have_results) {
        ImGui::TextDisabled(loader_.busy() ? "Reading the results..." : "No result step loaded.");
    } else {
        if (!steps_.empty())
            ImGui::TextDisabled("Output %d of %d, t = %.6g s%s", shown_step + 1, last_step + 1,
                                shown.time, record ? (record->converged ? ", converged" : ", not converged") : "");
        if (ImGui::BeginTable("key_results", 4, ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_RowBg |
                                                     ImGuiTableFlags_SizingStretchProp)) {
            const auto cell = [](const char* name, const char* help, const std::string& value) {
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(name);
                ui::item_tooltip(help);
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(value.c_str());
            };
            cell("Minimum film thickness h_min", "Smallest gap in the film.",
                 value_text(shown.h_min, gui::UnitFamily::SmallLength, unit_map, result_data::kUnitFilm, "%.3g"));
            cell("Eccentricity ratio", "Shaft centre offset per radial clearance, e / c.",
                 plain_text(shown.eccentricity_ratio, "%.3f"));
            cell("Attitude angle", "Angle from the load direction to the line of centres.",
                 plain_text(shown.attitude_deg, "%.1f deg"));
            cell("Friction torque", "Viscous torque on the shaft.",
                 value_text(shown.friction_torque, gui::UnitFamily::Torque, unit_map, result_data::kUnitTorque));
            cell("Friction power loss", "Friction torque times angular speed.",
                 value_text(shown.friction_power, gui::UnitFamily::Power, unit_map, result_data::kUnitPower));
            cell("Peak pressure", "Highest film pressure (absolute).",
                 value_text(shown.p_max, gui::UnitFamily::Pressure, unit_map, result_data::kUnitPressure));
            cell("Cavitated area", "Share of the film area filled partly by gas or vapour.",
                 plain_text(shown.cavitated_fraction * 100.0, "%.1f %%"));
            if (std::isfinite(shown.lambda))
                cell("Film thickness ratio Lambda", "Minimum film per combined surface roughness.",
                     plain_text(shown.lambda, "%.2f"));
            else
                cell("Lowest film pressure", "Lowest film pressure (absolute).",
                     value_text(shown.p_min, gui::UnitFamily::Pressure, unit_map, result_data::kUnitPressure));
            ImGui::EndTable();
        }
        ImGui::Text("Shear share of the film force: %s", plain_text(shown.shear_share * 100.0, "%.1f %%").c_str());
        ui::item_tooltip("Viscous (shear) part of the total film force; the rest is pressure.");
        const gui::DesignChecks checks = gui::design_checks(shown, p_sat, model.dissolved_gas);
        if (!checks.gas_release_text.empty()) mark(risk_light(checks.gas_release_risk), checks.gas_release_text);
        else ImGui::TextDisabled("Gas release check: the saturation pressure is unknown (choose a fluid pair).");
        if (!checks.lambda_text.empty()) mark(risk_light(checks.lambda_risk), checks.lambda_text);
        else ImGui::TextDisabled("Film thickness ratio: enter the surface roughness to check it.");
    }
    if (record) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Command");
        ImGui::SameLine();
        if (fonts.mono) ImGui::PushFont(fonts.mono);
        ImGui::TextWrapped("%s", record->command.c_str());
        if (fonts.mono) ImGui::PopFont();
        if (ImGui::SmallButton("Copy command")) {
            ImGui::SetClipboardText(record->command.c_str());
            set_status("Command copied.");
        }
        ui::item_tooltip("Runs the same case from a terminal.");
        ImGui::SameLine();
    }

    // --- Export ----------------------------------------------------------------
    if (ImGui::SmallButton("PNG")) {
        if (s.plot_max.x > s.plot_min.x) {
            const fs::path file = exports::save_file_dialog(hwnd, L"PNG image (*.png)\0*.png\0", L"png",
                                                            winutil::utf8_to_wide(s.plot_name + ".png"));
            if (!file.empty())
                pending_capture = CaptureRequest{s.plot_min.x, s.plot_min.y, s.plot_max.x, s.plot_max.y, file};
        } else {
            set_status("Show a plot first.");
        }
    }
    ui::item_tooltip("Save the plot shown in the tab below as a PNG image.");
    ImGui::SameLine();
    if (ImGui::SmallButton("CSV: key results") && have_results) {
        const fs::path file = exports::save_file_dialog(hwnd, L"CSV file (*.csv)\0*.csv\0", L"csv", L"key_results.csv");
        if (!file.empty()) {
            const auto row = [](const char* name, double value, const char* unit) {
                return std::vector<std::string>{name, plain_text(value, "%.8g"), unit};
            };
            const std::vector<std::vector<std::string>> rows = {
                row("minimum film thickness", shown.h_min, "m"),
                row("eccentricity ratio", shown.eccentricity_ratio, "-"),
                row("attitude angle", shown.attitude_deg, "deg"),
                row("friction torque", shown.friction_torque, "N m"),
                row("friction power loss", shown.friction_power, "W"),
                row("peak pressure", shown.p_max, "Pa"),
                row("lowest film pressure", shown.p_min, "Pa"),
                row("cavitated area", shown.cavitated_fraction, "-"),
                row("shear share of film force", shown.shear_share, "-"),
                row("film thickness ratio", shown.lambda, "-"),
                row("maximum temperature", shown.t_max, "K"),
                row("time", shown.time, "s")};
            std::string error;
            if (exports::write_csv_rows(file, {"quantity", "value", "unit"}, rows, error))
                set_status("Saved " + winutil::path_to_utf8(file));
            else set_status("CSV export failed: " + error);
        }
    }
    ui::item_tooltip("Key results of the step shown, in SI units.");
    ImGui::SameLine();
    if (ImGui::SmallButton("CSV: profiles") && grid_) {
        const fs::path file = exports::save_file_dialog(hwnd, L"CSV file (*.csv)\0*.csv\0", L"csv", L"profiles.csv");
        if (!file.empty()) {
            const Profile profile = make_profile(*grid_, zero, start, s.profile_z);
            std::string error;
            if (exports::write_csv_columns(file, {"angle_deg (" + x_label + ")", "p_Pa", "h_m", "gas_fraction", "T_K"},
                                           {profile.angle, profile.p, profile.h, profile.gas, profile.t}, error))
                set_status("Saved " + winutil::path_to_utf8(file));
            else set_status("CSV export failed: " + error);
        }
    }
    ui::item_tooltip("Pressure, film thickness, gas fraction and temperature along the angle at the axial "
                     "position of the Profiles tab, in SI units.");
    ImGui::SameLine();
    if (ImGui::SmallButton("Case report...")) {
        const std::wstring name = (case_path_.empty() ? std::wstring(L"case") : case_path_.stem().wstring()) + L"_report.html";
        const fs::path file = exports::save_file_dialog(hwnd, L"HTML page (*.html)\0*.html\0", L"html", name);
        if (!file.empty()) export_case_report(file);
    }
    ui::item_tooltip("One printable page: inputs, estimates, key results, checks and figures.");

    // --- Step selector -------------------------------------------------------------
    if (steps_.size() > 1) {
        int index = shown_step;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16.0f);
        if (ImGui::SliderInt("Output step", &index, 0, last_step, "%d")) {
            selected_step_ = index == last_step ? -1 : index;
            s.playing = false;
        }
        ui::item_tooltip("Saved output shown in the field maps, polar plot and profiles.");
        ImGui::SameLine();
        ImGui::TextDisabled("t = %.6g s", steps_[static_cast<std::size_t>(std::max(shown_step, 0))].time);
    }

    // --- Tabs ------------------------------------------------------------------------
    const bool transient = model.analysis == gui::AnalysisType::Transient ||
                           cfg.solution_mode == SolutionMode::TRANSIENT;
    int on = -1, off = -1;
    for (int i = record_count - 1; i >= 0; --i) {
        const RunRecord& r = records_[static_cast<std::size_t>(i)];
        if (r.job.group != "comparison") continue;
        if (on < 0 && r.job.label == kGasOnLabel) on = i;
        if (off < 0 && r.job.label == kGasOffLabel) off = i;
    }
    if (ImGui::BeginTabBar("results_tabs")) {
        if (ImGui::BeginTabItem("Field maps")) {
            const bool has_gas = grid_ && (grid_->field("alpha_gas") || grid_->field("film_content"));
            const bool has_t = grid_ && grid_->field("T");
            const auto radio = [&](const char* label, MapField field, bool enabled, const char* help) {
                if (!enabled) return;
                if (ImGui::RadioButton(label, s.field == field)) s.field = field;
                ui::item_tooltip(help);
                ImGui::SameLine();
            };
            radio("Pressure", MapField::Pressure, true, "Film pressure (absolute).");
            radio("Film thickness", MapField::Film, true, "Gap between shaft and bearing.");
            radio("Gas fraction", MapField::Gas, has_gas, "Volume fraction of gas or vapour in the gap.");
            radio("Temperature", MapField::Temperature, has_t, "Film temperature.");
            ImGui::NewLine();
            if ((s.field == MapField::Gas && !has_gas) || (s.field == MapField::Temperature && !has_t))
                s.field = MapField::Pressure;
            if (!grid_) {
                ImGui::TextDisabled(loader_.busy() ? "Reading the results..." : "No result step loaded.");
            } else {
                const gui::UnitOption* unit = nullptr;
                std::string bar = "Gas fraction [-]";
                if (s.field == MapField::Pressure) {
                    unit = &result_data::unit_of(unit_map, result_data::kUnitPressure, gui::UnitFamily::Pressure);
                    bar = std::string("Pressure [") + unit->label + "]";
                } else if (s.field == MapField::Film) {
                    unit = &result_data::unit_of(unit_map, result_data::kUnitFilm, gui::UnitFamily::SmallLength);
                    bar = std::string("Film thickness [") + unit->label + "]";
                } else if (s.field == MapField::Temperature) {
                    unit = &result_data::unit_of(unit_map, result_data::kUnitTemperature, gui::UnitFamily::Temperature);
                    bar = std::string("Temperature [") + unit->label + "]";
                }
                MapCache& m = s.map;
                if (m.grid != grid_.get() || m.field != s.field || m.zero != zero || m.start != start ||
                    m.scale != (unit ? unit->scale : 1.0) || m.offset != (unit ? unit->offset : 0.0))
                    build_map(m, *grid_, s.field, zero, start, unit);
                if (!m.available) {
                    ImGui::TextDisabled("This field is not in the output.");
                } else {
                    const double length = model.bearing_length > 0.0 ? model.bearing_length : 1.0;
                    const double aspect = model.bearing_diameter > 0.0 ? length / (kPi * model.bearing_diameter) : 0.3;
                    const double feed_x = gui::display_angle(model.feed.angle_deg, zero, start);
                    const double feed_y = model.feed.type == gui::FeedType::Hole && model.feed.hole_axial_position >= 0.0
                                              ? std::clamp(model.feed.hole_axial_position / length, 0.0, 1.0)
                                              : 0.5;
                    const double hmin_x = std::isfinite(grid_r.h_min_theta_deg)
                                              ? gui::display_angle(grid_r.h_min_theta_deg, zero, start) : kNaN;
                    draw_field_map(s, m, bar, x_label, start, std::clamp(aspect, 0.05, 2.0), feed_x, feed_y, hmin_x);
                }
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Polar")) {
            if (grid_) draw_polar(s, *grid_, model, zero, start, grid_r, p_sat, unit_map);
            else ImGui::TextDisabled("No result step loaded.");
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Profiles")) {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
            ImGui::SliderFloat("Axial position z / L", &s.profile_z, 0.0f, 1.0f, "%.2f");
            ui::item_tooltip("Where along the bearing length the profiles are taken (0.5 = mid-plane).");
            if (grid_) {
                const Profile profile = make_profile(*grid_, zero, start, s.profile_z);
                const auto& pu = result_data::unit_of(unit_map, result_data::kUnitPressure, gui::UnitFamily::Pressure);
                const auto& hu = result_data::unit_of(unit_map, result_data::kUnitFilm, gui::UnitFamily::SmallLength);
                const auto p = scaled(profile.p, pu), h = scaled(profile.h, hu);
                const std::string p_axis = std::string("Pressure [") + pu.label + "]";
                const std::string h_axis = std::string("Film thickness [") + hu.label + "]";
                const float height = std::max(260.0f, ImGui::GetContentRegionAvail().y);
                ImGui::BeginGroup();
                if (ImPlot::BeginPlot("##profiles", ImVec2(-1, height))) {
                    ImPlot::SetupAxis(ImAxis_X1, x_label.c_str());
                    ImPlot::SetupAxis(ImAxis_Y1, p_axis.c_str(), ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxis(ImAxis_Y2, h_axis.c_str(), ImPlotAxisFlags_AuxDefault | ImPlotAxisFlags_AutoFit);
                    ImPlot::SetupAxisLimits(ImAxis_X1, start, start + 360.0, ImPlotCond_Always);
                    ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
                    ImPlot::SetNextLineStyle(theme::kAccent, 2.0f);
                    if (!p.empty()) ImPlot::PlotLine("Pressure", profile.angle.data(), p.data(), static_cast<int>(p.size()));
                    const double ambient = gui::to_display(pu, model.ambient_pressure);
                    ImPlot::SetNextLineStyle(theme::kMuted, 1.0f);
                    ImPlot::PlotInfLines("Ambient pressure", &ambient, 1, ImPlotInfLinesFlags_Horizontal);
                    if (p_sat > 0.0) {
                        const double sat = gui::to_display(pu, p_sat);
                        ImPlot::SetNextLineStyle(theme::kHighlight, 1.0f);
                        ImPlot::PlotInfLines("Saturation pressure", &sat, 1, ImPlotInfLinesFlags_Horizontal);
                    }
                    ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
                    ImPlot::SetNextLineStyle(theme::kWarn, 2.0f);
                    if (!h.empty()) ImPlot::PlotLine("Film thickness", profile.angle.data(), h.data(), static_cast<int>(h.size()));
                    ImPlot::EndPlot();
                }
                ImGui::EndGroup();
                remember_plot(s, "profiles");
            } else {
                ImGui::TextDisabled("No result step loaded.");
            }
            ImGui::EndTabItem();
        }
        const fs::path base = record ? record->job.output_dir
                                     : (results_dir_.filename() == "flat" ? results_dir_.parent_path() : results_dir_);
        if (transient && ImGui::BeginTabItem("Transient")) {
            ensure_motion(s.motion, base / "motion.csv");
            if (ImGui::SmallButton("Reload")) s.motion.path.clear();
            ImGui::SameLine();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9.0f);
            ImGui::Combo("x axis", &s.transient_x, "Time\0Crank angle\0");
            ui::item_tooltip("Plot against time or against the shaft angle turned since the start.");
            ImGui::SameLine();
            if (ImGui::SmallButton(s.playing ? "Pause" : "Play")) {
                s.playing = !s.playing;
                s.last_tick = ImGui::GetTime();
            }
            ui::item_tooltip("Step the field maps through the saved outputs (about 5 per second).");
            if (s.playing && steps_.size() > 1 && !loader_.busy() && ImGui::GetTime() - s.last_tick >= 0.2) {
                s.last_tick = ImGui::GetTime();
                selected_step_ = (shown_step + 1) % static_cast<int>(steps_.size());
            }
            const CsvTable& t = s.motion.table;
            const auto* time = t.column("time_s");
            if (s.motion.loading) {
                ImGui::TextDisabled("Reading motion.csv...");
            } else if (!time || time->empty()) {
                ImGui::TextDisabled("No time history: %s", s.motion.error.empty() ? "motion.csv is empty" : s.motion.error.c_str());
            } else {
                const double factor = s.transient_x == 1 ? std::abs(model.speed) * 180.0 / kPi : 1.0;
                std::vector<double> x(time->size());
                for (std::size_t i = 0; i < x.size(); ++i) x[i] = (*time)[i] * factor;
                const auto& hu = result_data::unit_of(unit_map, result_data::kUnitFilm, gui::UnitFamily::SmallLength);
                const auto& tu = result_data::unit_of(unit_map, result_data::kUnitTorque, gui::UnitFamily::Torque);
                std::vector<double> h, torque, cavitated, load;
                if (const auto* c = t.column("h_min_m")) h = scaled(*c, hu);
                if (const auto* c = t.column("friction_torque_Nm")) torque = scaled(*c, tu);
                if (const auto* c = t.column("cavitated_fraction")) for (double v : *c) cavitated.push_back(v * 100.0);
                const auto* lx = t.column("load_x_N");
                const auto* ly = t.column("load_y_N");
                if (lx && ly) for (std::size_t i = 0; i < lx->size(); ++i) load.push_back(std::hypot((*lx)[i], (*ly)[i]));
                std::vector<double> ecc;
                if (const auto* c = t.column("eccentricity_ratio")) ecc = *c;
                const std::string h_label = std::string("h_min [") + hu.label + "]";
                const std::string t_label = std::string("Torque [") + tu.label + "]";
                struct Series { const char* label; const std::vector<double>* y; };
                const Series series[] = {{h_label.c_str(), &h}, {"Ecc. ratio", &ecc}, {t_label.c_str(), &torque},
                                         {"Cav. area [%]", &cavitated}, {"Load [N]", &load}};
                if (!steps_.empty())
                    s.cursor = steps_[static_cast<std::size_t>(std::max(shown_step, 0))].time * factor;
                bool dragged = false;
                const float height = std::max(ImGui::GetFontSize() * 30.0f, ImGui::GetContentRegionAvail().y);
                ImGui::BeginGroup();
                if (ImPlot::BeginSubplots("##motion", 5, 1, ImVec2(-1, height), ImPlotSubplotFlags_LinkAllX)) {
                    for (int k = 0; k < 5; ++k) {
                        if (!ImPlot::BeginPlot(series[k].label, ImVec2(), ImPlotFlags_NoLegend | ImPlotFlags_NoTitle)) continue;
                        ImPlot::SetupAxes(k == 4 ? (s.transient_x == 1 ? "Crank angle [deg]" : "Time [s]") : nullptr,
                                          series[k].label,
                                          ImPlotAxisFlags_AutoFit | (k < 4 ? ImPlotAxisFlags_NoTickLabels : 0),
                                          ImPlotAxisFlags_AutoFit);
                        const std::vector<double>& y = *series[k].y;
                        const int n = static_cast<int>(std::min(x.size(), y.size()));
                        ImPlot::SetNextLineStyle(theme::kAccent, 1.5f);
                        if (n > 0) ImPlot::PlotLine("##y", x.data(), y.data(), n);
                        if (ImPlot::DragLineX(0, &s.cursor, theme::kHighlight, 1.5f, ImPlotDragToolFlags_NoFit)) dragged = true;
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
                ImGui::EndGroup();
                remember_plot(s, "time_history");
                if (dragged && !steps_.empty()) {
                    const double t_cursor = s.cursor / factor;
                    int best = 0;
                    for (int i = 1; i <= last_step; ++i)
                        if (std::abs(steps_[static_cast<std::size_t>(i)].time - t_cursor) <
                            std::abs(steps_[static_cast<std::size_t>(best)].time - t_cursor)) best = i;
                    selected_step_ = best == last_step ? -1 : best;
                    s.playing = false;
                }
            }
            ImGui::EndTabItem();
        }
        if ((on >= 0 || off >= 0) && ImGui::BeginTabItem("Comparison")) {
            const RunRecord* runs[2] = {on >= 0 ? &records_[static_cast<std::size_t>(on)] : nullptr,
                                        off >= 0 ? &records_[static_cast<std::size_t>(off)] : nullptr};
            for (int k = 0; k < 2; ++k)
                if (runs[k]) result_data::ensure_compare(s.compare[k], runs[k]->job.output_dir, runs[k]->job.config);
            const auto result_of = [&](int k) -> const gui::KeyResults* {
                if (!runs[k]) return nullptr;
                if (runs[k]->results_ready) return &runs[k]->results;
                return s.compare[k].result.grid ? &s.compare[k].result.results : nullptr;
            };
            const gui::KeyResults* a = result_of(0);
            const gui::KeyResults* b = result_of(1);
            if (ImGui::BeginTable("comparison", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                                       ImGuiTableFlags_SizingStretchProp)) {
                ImGui::TableSetupColumn("Quantity");
                ImGui::TableSetupColumn("Gas on (GC+VC)");
                ImGui::TableSetupColumn("Gas off (VC only)");
                ImGui::TableSetupColumn("Difference [%]");
                ImGui::TableHeadersRow();
                const auto row = [&](const char* name, double gui::KeyResults::*field, gui::UnitFamily family,
                                     const char* id, double factor, const char* format) {
                    const double va = a ? a->*field : kNaN, vb = b ? b->*field : kNaN;
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(name);
                    for (double v : {va, vb}) {
                        ImGui::TableNextColumn();
                        ImGui::TextUnformatted(id ? value_text(v, family, unit_map, id).c_str()
                                                  : plain_text(v * factor, format).c_str());
                    }
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(plain_text(gui::percent_difference(vb, va), "%+.1f").c_str());
                };
                using K = gui::KeyResults;
                row("Minimum film thickness", &K::h_min, gui::UnitFamily::SmallLength, result_data::kUnitFilm, 1, nullptr);
                row("Eccentricity ratio", &K::eccentricity_ratio, gui::UnitFamily::Fixed, nullptr, 1, "%.3f");
                row("Attitude angle [deg]", &K::attitude_deg, gui::UnitFamily::Fixed, nullptr, 1, "%.1f");
                row("Friction torque", &K::friction_torque, gui::UnitFamily::Torque, result_data::kUnitTorque, 1, nullptr);
                row("Friction power loss", &K::friction_power, gui::UnitFamily::Power, result_data::kUnitPower, 1, nullptr);
                row("Peak pressure", &K::p_max, gui::UnitFamily::Pressure, result_data::kUnitPressure, 1, nullptr);
                row("Lowest film pressure", &K::p_min, gui::UnitFamily::Pressure, result_data::kUnitPressure, 1, nullptr);
                row("Cavitated area [%]", &K::cavitated_fraction, gui::UnitFamily::Fixed, nullptr, 100, "%.1f");
                row("Shear share of film force [%]", &K::shear_share, gui::UnitFamily::Fixed, nullptr, 100, "%.1f");
                ImGui::EndTable();
            }
            ImGui::TextDisabled("Difference: (gas on - gas off) / gas off.");
            const auto& pu = result_data::unit_of(unit_map, result_data::kUnitPressure, gui::UnitFamily::Pressure);
            const std::string p_axis = std::string("Mid-plane pressure [") + pu.label + "]";
            const float height = std::max(240.0f, ImGui::GetContentRegionAvail().y);
            ImGui::BeginGroup();
            if (ImPlot::BeginPlot("##comparison_profiles", ImVec2(-1, height))) {
                ImPlot::SetupAxes(x_label.c_str(), p_axis.c_str(), 0, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisLimits(ImAxis_X1, start, start + 360.0, ImPlotCond_Always);
                for (int k = 0; k < 2; ++k) {
                    const CompareResult& c = s.compare[k].result;
                    if (!runs[k] || !c.grid) continue;
                    const double run_zero = result_data::zero_position(angle_zero_, runs[k]->job.model.feed.angle_deg, c.results);
                    const Profile profile = make_profile(*c.grid, run_zero, start, 0.5);
                    const auto p = scaled(profile.p, pu);
                    ImPlot::SetNextLineStyle(k == 0 ? theme::kModelGcVc : theme::kModelVc, 2.0f);
                    if (!p.empty()) ImPlot::PlotLine(k == 0 ? kGasOnLabel : kGasOffLabel, profile.angle.data(), p.data(),
                                                     static_cast<int>(p.size()));
                }
                const double sat = runs[0] ? runs[0]->job.saturation_pressure : p_sat;
                if (sat > 0.0) {
                    const double v = gui::to_display(pu, sat);
                    ImPlot::SetNextLineStyle(theme::kHighlight, 1.0f);
                    ImPlot::PlotInfLines("Saturation pressure", &v, 1, ImPlotInfLinesFlags_Horizontal);
                }
                ImPlot::EndPlot();
            }
            ImGui::EndGroup();
            remember_plot(s, "comparison_profiles");
            for (int k = 0; k < 2; ++k)
                if (runs[k] && !s.compare[k].loading && !s.compare[k].result.error.empty())
                    ImGui::TextDisabled("%s: %s", k == 0 ? kGasOnLabel : kGasOffLabel, s.compare[k].result.error.c_str());
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}
