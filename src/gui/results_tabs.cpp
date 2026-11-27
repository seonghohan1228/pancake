// Tabs of the Results window: Field maps | Polar | Profiles | Transient |
// Comparison. The window (run picker, key-results card, export, step selector)
// is results_window.cpp. Field maps follow the run's model (plan 0.5): FBNS
// shows the void fraction 1 - theta and never a gas field.

#include "gui_app.hpp"

#include "imgui.h"
#include "imgui_internal.h"
#include "implot.h"
#include "result_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace result_data {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr const char* kGasOnLabel = "gas on (GC+VC)";
constexpr const char* kGasOffLabel = "gas off (VC only)";

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

/// Wrapped disabled text.
void note(const char* text) {
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", text);
    ImGui::PopTextWrapPos();
}

/// Display unit of a view: the Results unit choice for p, h and T; as stored otherwise.
const gui::UnitOption* view_unit(const gui::FieldView& view, ui::UnitChoices& units) {
    if (view.id == "p") return &unit_of(units, kUnitPressure, gui::UnitFamily::Pressure);
    if (view.id == "h") return &unit_of(units, kUnitFilm, gui::UnitFamily::SmallLength);
    if (view.id == "T") return &unit_of(units, kUnitTemperature, gui::UnitFamily::Temperature);
    return nullptr;
}

/// Text placed at `pos` unless it overlaps a rectangle already taken; tries a few
/// positions below. Returns false when no free place was found.
bool place_text(ImDrawList* draw, std::vector<ImRect>& taken, const ImRect& area, ImVec2 pos, ImU32 colour,
                const char* text) {
    const ImVec2 size = ImGui::CalcTextSize(text);
    for (int k = 0; k < 4; ++k) {
        ImVec2 at(std::clamp(pos.x, area.Min.x, std::max(area.Min.x, area.Max.x - size.x)),
                  std::clamp(pos.y + k * size.y, area.Min.y, std::max(area.Min.y, area.Max.y - size.y)));
        const ImRect rect(at, ImVec2(at.x + size.x, at.y + size.y));
        bool free = true;
        for (const ImRect& other : taken)
            if (rect.Overlaps(other)) free = false;
        if (!free) continue;
        taken.push_back(rect);
        // Dark backing so white labels stay readable on light parts of the colormap.
        draw->AddRectFilled(ImVec2(rect.Min.x - 2.0f, rect.Min.y), ImVec2(rect.Max.x + 2.0f, rect.Max.y),
                            IM_COL32(0, 0, 0, 140), 3.0f);
        draw->AddText(at, colour, text);
        return true;
    }
    return false;
}

void draw_field_map(PlotCapture& capture, const MapCache& m, const std::string& bar_label,
                    const std::string& x_label, double start, double aspect, double feed_x, double feed_y,
                    double hmin_x) {
    const float font = ImGui::GetFontSize();
    const float bar_w = font * 6.5f;
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float pad_x = font * 4.0f, pad_y = font * 3.6f;
    float inner_w = std::max(avail.x - bar_w - ImGui::GetStyle().ItemSpacing.x - pad_x, 120.0f);
    float inner_h = inner_w * static_cast<float>(aspect);
    const float max_h = std::max(avail.y - ImGui::GetTextLineHeightWithSpacing() * 4.0f - pad_y, 80.0f);
    if (inner_h > max_h) {
        inner_h = max_h;
        inner_w = inner_h / static_cast<float>(aspect);
    }
    const ImVec2 size(inner_w + pad_x, inner_h + pad_y);
    ImGui::BeginGroup();
    ImPlot::PushColormap(ImPlotColormap_Viridis);
    if (ImPlot::BeginPlot("##field_map", size, ImPlotFlags_NoLegend | ImPlotFlags_NoMenus)) {
        ImPlot::SetupAxes(x_label.c_str(), "z / L", ImPlotAxisFlags_NoGridLines, ImPlotAxisFlags_NoGridLines);
        ImPlot::SetupAxisLimits(ImAxis_X1, start, start + 360.0, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
        // A flat map (pi D >> L) leaves little height: fewer z/L ticks so they never collide.
        static const double few_ticks[] = {0.0, 1.0};
        static const double some_ticks[] = {0.0, 0.5, 1.0};
        if (inner_h < 3.0f * font) ImPlot::SetupAxisTicks(ImAxis_Y1, few_ticks, 2);
        else if (inner_h < 7.0f * font) ImPlot::SetupAxisTicks(ImAxis_Y1, some_ticks, 3);
        ImPlot::PlotHeatmap("##values", m.values.data(), m.rows, m.cols, m.vmin, m.vmax, nullptr,
                            ImPlotPoint(m.x_left, 0.0), ImPlotPoint(m.x_left + 360.0, 1.0));
        ImDrawList* draw = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        const ImU32 outline = ImGui::GetColorU32(theme::kHighlight);
        for (std::size_t i = 0; i + 3 < m.outline.size(); i += 4)
            draw->AddLine(ImPlot::PlotToPixels(m.outline[i], m.outline[i + 1]),
                          ImPlot::PlotToPixels(m.outline[i + 2], m.outline[i + 3]), outline, 2.0f);
        if (std::isfinite(hmin_x)) {
            ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 0.9f), 1.5f);
            ImPlot::PlotInfLines("##hmin", &hmin_x, 1);
        }
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Down, 7.0f, ImVec4(1, 1, 1, 1), 1.0f, ImVec4(0, 0, 0, 1));
        ImPlot::PlotScatter("##feed", &feed_x, &feed_y, 1);
        // Labels in pixel space; the second one moves down when they would overlap.
        const ImVec2 plot_pos = ImPlot::GetPlotPos(), plot_size = ImPlot::GetPlotSize();
        const ImRect area(plot_pos.x + 2.0f, plot_pos.y + 2.0f, plot_pos.x + plot_size.x - 2.0f,
                          plot_pos.y + plot_size.y - 2.0f);
        std::vector<ImRect> taken;
        const ImU32 white = IM_COL32(255, 255, 255, 255);
        if (std::isfinite(hmin_x)) {
            const ImVec2 at = ImPlot::PlotToPixels(hmin_x, 0.85);
            place_text(draw, taken, area, ImVec2(at.x + 4.0f, at.y - 0.5f * font), white, "h_min");
        }
        const ImVec2 feed_px = ImPlot::PlotToPixels(feed_x, feed_y);
        const float feed_w = ImGui::CalcTextSize("feed").x;
        place_text(draw, taken, area, ImVec2(feed_px.x - 0.5f * feed_w, feed_px.y - 2.0f * font), white, "feed");
        ImPlot::PopPlotClipRect();
        ImPlot::EndPlot();
    }
    ImGui::SameLine();
    // The colour-bar title runs along the bar; drop the quantity when it does not fit.
    std::string bar = bar_label;
    const std::size_t bracket = bar.find(" [");
    if (ImGui::CalcTextSize(bar.c_str()).x > size.y - font && bracket != std::string::npos) bar = bar.substr(bracket + 1);
    // The colour bar keeps a readable height even beside a flat map.
    const float bar_h = std::max(size.y, 9.0f * font);
    if (ImGui::CalcTextSize(bar.c_str()).x > bar_h - font && bracket != std::string::npos) bar = bar.substr(bracket + 1);
    ImPlot::ColormapScale(bar.c_str(), m.vmin, m.vmax, ImVec2(bar_w, bar_h), "%.3g");
    ImPlot::PopColormap();
    ImGui::EndGroup();
    remember_plot(capture, "field_map");
}

}  // namespace

void remember_plot(PlotCapture& capture, const char* name) {
    const ImVec2 min = ImGui::GetItemRectMin(), max = ImGui::GetItemRectMax();
    capture.min_x = min.x;
    capture.min_y = min.y;
    capture.max_x = max.x;
    capture.max_y = max.y;
    capture.name = name;
}

void draw_field_maps_tab(FieldMapUi& state, const TabContext& tab, PlotCapture& capture) {
    if (!tab.grid) {
        ImGui::TextDisabled(tab.loading ? "Reading the results..." : "No result step loaded.");
        return;
    }
    const SimulationConfig& cfg = *tab.config;
    const gui::CaseModel& model = *tab.model;
    const std::vector<gui::FieldView> views = gui::result_fields(cfg, *tab.grid);
    if (views.empty()) {
        ImGui::TextDisabled("This output has no field the maps can show.");
        return;
    }
    const auto selected = std::find_if(views.begin(), views.end(),
                                       [&](const gui::FieldView& v) { return v.id == state.view_id; });
    if (selected == views.end()) state.view_id = views.front().id;

    // View selector: radio buttons that wrap to the next line instead of running off the window.
    const ImGuiStyle& style = ImGui::GetStyle();
    for (std::size_t i = 0; i < views.size(); ++i) {
        const gui::FieldView& view = views[i];
        if (i > 0) {
            const float width = ImGui::GetFrameHeight() + style.ItemInnerSpacing.x +
                                ImGui::CalcTextSize(view.label.c_str()).x;
            ImGui::SameLine();
            if (ImGui::GetContentRegionAvail().x < width) ImGui::NewLine();
        }
        if (ImGui::RadioButton(view.label.c_str(), state.view_id == view.id)) state.view_id = view.id;
        ui::item_tooltip(view.tooltip.c_str());
    }
    const gui::FieldView& view = *std::find_if(views.begin(), views.end(),
                                               [&](const gui::FieldView& v) { return v.id == state.view_id; });

    const gui::UnitOption* unit = view_unit(view, *tab.units);
    const std::string bar = view.label + " [" + (unit ? std::string(unit->label) : view.unit) + "]";
    MapCache& m = state.map;
    if (m.grid != tab.grid || m.view_id != view.id || m.form != gui::result_form(cfg) || m.zero != tab.zero ||
        m.start != tab.start || m.scale != (unit ? unit->scale : 1.0) || m.offset != (unit ? unit->offset : 0.0))
        build_map(m, *tab.grid, view, cfg, tab.zero, tab.start, unit);
    if (!m.available) {
        ImGui::TextDisabled("This field is not in the output.");
        return;
    }
    const double length = model.bearing_length > 0.0 ? model.bearing_length : 1.0;
    const double aspect = model.bearing_diameter > 0.0 ? length / (kPi * model.bearing_diameter) : 0.3;
    const double feed_x = gui::display_angle(model.feed.angle_deg, tab.zero, tab.start);
    const double feed_y = model.feed.type == gui::FeedType::Hole && model.feed.hole_axial_position >= 0.0
                              ? std::clamp(model.feed.hole_axial_position / length, 0.0, 1.0)
                              : 0.5;
    const double hmin_x = std::isfinite(tab.grid_results.h_min_theta_deg)
                              ? gui::display_angle(tab.grid_results.h_min_theta_deg, tab.zero, tab.start)
                              : kNaN;
    draw_field_map(capture, m, bar, tab.x_label, tab.start, std::clamp(aspect, 0.05, 2.0), feed_x, feed_y, hmin_x);

    const gui::ResultForm form = gui::result_form(cfg);
    note(form == gui::ResultForm::Coupled
             ? "Purple outline: cavitated region (free-gas volume fraction 0.5). White line: minimum film. "
               "Width and height to scale (pi D : L)."
             : "Purple outline: cavitated region (void fraction 0.5). White line: minimum film. "
               "Width and height to scale (pi D : L).");
    if (form == gui::ResultForm::Fbns)
        note("Pressure-plateau model: the cavity is liquid streaks plus a massless void; gas is not tracked.");
}

void draw_polar_tab(const TabContext& tab, PlotCapture& capture) {
    if (!tab.grid) {
        ImGui::TextDisabled("No result step loaded.");
        return;
    }
    const ResultGrid& grid = *tab.grid;
    const gui::CaseModel& model = *tab.model;
    ui::UnitChoices& units = *tab.units;
    const double zero = tab.zero, start = tab.start, p_sat = tab.p_sat;
    const auto* p = grid.field("p");
    if (!p) {
        ImGui::TextDisabled("This output has no pressure field.");
        return;
    }
    const auto order = angle_order(grid, zero, start);
    const auto values = row_values(grid, *p, grid.ny / 2, order);
    const double p_amb = model.ambient_pressure;
    double hi = 0.0, lo = 0.0;
    for (double v : values)
        if (std::isfinite(v)) hi = std::max(hi, v - p_amb), lo = std::min(lo, v - p_amb);
    if (p_sat > 0.0) hi = std::max(hi, p_sat - p_amb), lo = std::min(lo, p_sat - p_amb);

    const ImVec2 avail = ImGui::GetContentRegionAvail();
    const float side = std::max(220.0f, std::min(avail.x, avail.y - ImGui::GetTextLineHeightWithSpacing() * 4.5f));
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
    std::vector<ImRect> taken;  // angle and marker labels
    draw->AddCircleFilled(centre, r_max * 1.05f, ImGui::GetColorU32(ImGuiCol_FrameBg), 128);
    for (int q = 0; q < 12; ++q) {
        const double angle = start + 30.0 * q;
        draw->AddLine(centre, at(angle, r_max * 1.05), grid_colour);
        if (q % 3 != 0) continue;
        char label[16];
        std::snprintf(label, sizeof(label), "%.0f", angle);
        const ImVec2 size = ImGui::CalcTextSize(label);
        const ImVec2 pos = at(angle, r_max * 1.16);
        const ImVec2 corner(pos.x - size.x * 0.5f, pos.y - size.y * 0.5f);
        taken.emplace_back(corner, ImVec2(corner.x + size.x, corner.y + size.y));
        draw->AddText(corner, text_colour, label);
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
    const ImRect area(origin, ImVec2(origin.x + side, origin.y + side));
    const auto marker = [&](double fixed_deg, const char* label) {
        const double display = gui::display_angle(fixed_deg, zero, start);
        const ImVec2 pos = at(display, r_max * 1.05);
        draw->AddCircleFilled(pos, 4.5f, text_colour);
        // Outside the ring, on the side of the marker; stacked outwards if taken.
        const ImVec2 size = ImGui::CalcTextSize(label);
        const float c = static_cast<float>(std::cos(display * kPi / 180.0));
        const float s = static_cast<float>(std::sin(display * kPi / 180.0));
        for (int k = 0; k < 4; ++k) {
            const ImVec2 anchor = at(display, r_max * 1.12 + k * size.y);
            ImVec2 corner(anchor.x + (c - 1.0f) * 0.5f * size.x, anchor.y - (s + 1.0f) * 0.5f * size.y);
            corner.x = std::clamp(corner.x, area.Min.x, area.Max.x - size.x);
            corner.y = std::clamp(corner.y, area.Min.y, area.Max.y - size.y);
            const ImRect rect(corner, ImVec2(corner.x + size.x, corner.y + size.y));
            bool free = true;
            for (const ImRect& other : taken)
                if (rect.Overlaps(other)) free = false;
            if (!free && k < 3) continue;
            taken.push_back(rect);
            draw->AddText(corner, text_colour, label);
            break;
        }
    };
    marker(model.feed.angle_deg, "feed");
    if (std::isfinite(tab.grid_results.h_min_theta_deg)) marker(tab.grid_results.h_min_theta_deg, "h_min");
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
                          value_text(values[best], gui::UnitFamily::Pressure, units, kUnitPressure).c_str());
    }
    ImGui::EndGroup();
    remember_plot(capture, "polar_pressure");
    const std::string legend =
        "Mid-plane pressure (blue): radius grows with the pressure above ambient. Grey circle: ambient pressure " +
        value_text(p_amb, gui::UnitFamily::Pressure, units, kUnitPressure) + ". Purple circle: saturation pressure " +
        (p_sat > 0.0 ? value_text(p_sat, gui::UnitFamily::Pressure, units, kUnitPressure) : std::string("unknown")) +
        ".";
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(legend.c_str());
    ImGui::PopTextWrapPos();
    note("Angles grow counter-clockwise, in the direction of rotation.");
}

void draw_profiles_tab(const TabContext& tab, float& profile_z, PlotCapture& capture) {
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
    ImGui::SliderFloat("Axial position z / L", &profile_z, 0.0f, 1.0f, "%.2f");
    ui::item_tooltip("Where along the bearing length the profiles are taken (0.5 = mid-plane).");
    if (!tab.grid) {
        ImGui::TextDisabled("No result step loaded.");
        return;
    }
    ui::UnitChoices& units = *tab.units;
    const Profile profile = make_profile(*tab.grid, *tab.config, tab.zero, tab.start, profile_z);
    const auto& pu = unit_of(units, kUnitPressure, gui::UnitFamily::Pressure);
    const auto& hu = unit_of(units, kUnitFilm, gui::UnitFamily::SmallLength);
    const auto p = scaled(profile.p, pu), h = scaled(profile.h, hu);
    const std::string p_axis = std::string("Pressure [") + pu.label + "]";
    const std::string h_axis = std::string("Film thickness [") + hu.label + "]";
    const float height = std::max(260.0f, ImGui::GetContentRegionAvail().y);
    ImGui::BeginGroup();
    if (ImPlot::BeginPlot("##profiles", ImVec2(-1, height))) {
        ImPlot::SetupAxis(ImAxis_X1, tab.x_label.c_str());
        const bool data = ui::has_points(profile.angle, p);
        ImPlot::SetupAxis(ImAxis_Y1, p_axis.c_str(), ui::fit_if(data));
        ImPlot::SetupAxis(ImAxis_Y2, h_axis.c_str(), ImPlotAxisFlags_AuxDefault | ui::fit_if(data));
        if (!data) {
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y2, 0.0, 1.0, ImPlotCond_Always);
        }
        ImPlot::SetupAxisLimits(ImAxis_X1, tab.start, tab.start + 360.0, ImPlotCond_Always);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast);
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
        ImPlot::SetNextLineStyle(theme::kAccent, 2.0f);
        if (!p.empty()) ImPlot::PlotLine("Pressure", profile.angle.data(), p.data(), static_cast<int>(p.size()));
        const double ambient = gui::to_display(pu, tab.model->ambient_pressure);
        ImPlot::SetNextLineStyle(theme::kMuted, 1.0f);
        ImPlot::PlotInfLines("Ambient pressure", &ambient, 1, ImPlotInfLinesFlags_Horizontal);
        if (tab.p_sat > 0.0) {
            const double sat = gui::to_display(pu, tab.p_sat);
            ImPlot::SetNextLineStyle(theme::kHighlight, 1.0f);
            ImPlot::PlotInfLines("Saturation pressure", &sat, 1, ImPlotInfLinesFlags_Horizontal);
        }
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
        ImPlot::SetNextLineStyle(theme::kWarn, 2.0f);
        if (!h.empty()) ImPlot::PlotLine("Film thickness", profile.angle.data(), h.data(), static_cast<int>(h.size()));
        ImPlot::EndPlot();
    }
    ImGui::EndGroup();
    remember_plot(capture, "profiles");
}

int draw_transient_tab(TransientUi& state, const TabContext& tab, const std::filesystem::path& motion_csv,
                       const std::vector<StepInfo>& steps, int shown_step, bool loader_busy, PlotCapture& capture) {
    int result = -2;
    const int last_step = static_cast<int>(steps.size()) - 1;
    ensure_motion(state.motion, motion_csv);
    if (ImGui::SmallButton("Reload")) state.motion.path.clear();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 9.0f);
    ImGui::Combo("x axis", &state.x_axis, "Time\0Crank angle\0");
    ui::item_tooltip("Plot against time or against the shaft angle turned since the start.");
    ImGui::SameLine();
    if (ImGui::SmallButton(state.playing ? "Pause" : "Play")) {
        state.playing = !state.playing;
        state.last_tick = ImGui::GetTime();
    }
    ui::item_tooltip("Step the field maps through the saved outputs (about 5 per second).");
    if (state.playing && steps.size() > 1 && !loader_busy && ImGui::GetTime() - state.last_tick >= 0.2) {
        state.last_tick = ImGui::GetTime();
        result = (shown_step + 1) % static_cast<int>(steps.size());
        if (result == last_step) result = -1;
    }
    const CsvTable& t = state.motion.table;
    const auto* time = t.column("time_s");
    if (state.motion.loading) {
        ImGui::TextDisabled("Reading motion.csv...");
        return result;
    }
    if (!time || time->empty()) {
        const std::string text = "No time history: " +
                                 (state.motion.error.empty() ? std::string("motion.csv is empty") : state.motion.error);
        note(text.c_str());
        return result;
    }
    ui::UnitChoices& units = *tab.units;
    const double factor = state.x_axis == 1 ? std::abs(tab.model->speed) * 180.0 / kPi : 1.0;
    std::vector<double> x(time->size());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = (*time)[i] * factor;
    const auto& hu = unit_of(units, kUnitFilm, gui::UnitFamily::SmallLength);
    const auto& tu = unit_of(units, kUnitTorque, gui::UnitFamily::Torque);
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
    if (!steps.empty()) state.cursor = steps[static_cast<std::size_t>(std::max(shown_step, 0))].time * factor;
    bool dragged = false;
    const float height = std::max(ImGui::GetFontSize() * 30.0f, ImGui::GetContentRegionAvail().y);
    ImGui::BeginGroup();
    if (ImPlot::BeginSubplots("##motion", 5, 1, ImVec2(-1, height), ImPlotSubplotFlags_LinkAllX)) {
        for (int k = 0; k < 5; ++k) {
            if (!ImPlot::BeginPlot(series[k].label, ImVec2(), ImPlotFlags_NoLegend | ImPlotFlags_NoTitle)) continue;
            const std::vector<double>& y = *series[k].y;
            const bool data = ui::has_points(x, y);
            ImPlot::SetupAxes(k == 4 ? (state.x_axis == 1 ? "Crank angle [deg]" : "Time [s]") : nullptr,
                              series[k].label, ui::fit_if(data) | (k < 4 ? ImPlotAxisFlags_NoTickLabels : 0),
                              ui::fit_if(data));
            if (!data) {
                ui::empty_plot();
                ImPlot::EndPlot();
                continue;
            }
            const int n = static_cast<int>(std::min(x.size(), y.size()));
            ImPlot::SetNextLineStyle(theme::kAccent, 1.5f);
            if (n > 0) ImPlot::PlotLine("##y", x.data(), y.data(), n);
            if (ImPlot::DragLineX(0, &state.cursor, theme::kHighlight, 1.5f, ImPlotDragToolFlags_NoFit)) dragged = true;
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }
    ImGui::EndGroup();
    remember_plot(capture, "time_history");
    if (dragged && !steps.empty()) {
        const double t_cursor = state.cursor / factor;
        int best = 0;
        for (int i = 1; i <= last_step; ++i)
            if (std::abs(steps[static_cast<std::size_t>(i)].time - t_cursor) <
                std::abs(steps[static_cast<std::size_t>(best)].time - t_cursor))
                best = i;
        result = best == last_step ? -1 : best;
        state.playing = false;
    }
    return result;
}

void draw_comparison_tab(ComparisonSide on, ComparisonSide off, const TabContext& tab, gui::AngleZero angle_zero,
                         PlotCapture& capture) {
    ui::UnitChoices& units = *tab.units;
    const ComparisonSide sides[2] = {on, off};
    for (const ComparisonSide& side : sides)
        if (side.record && side.run) ensure_compare(*side.run, side.record->job.output_dir, side.record->job.config);
    const auto result_of = [&](int k) -> const gui::KeyResults* {
        const ComparisonSide& side = sides[k];
        if (!side.record || !side.run) return nullptr;
        if (side.record->results_ready) return &side.record->results;
        return side.run->result.grid ? &side.run->result.results : nullptr;
    };
    const gui::KeyResults* a = result_of(0);
    const gui::KeyResults* b = result_of(1);
    if (ImGui::BeginTable("comparison", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                               ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Quantity", ImGuiTableColumnFlags_WidthStretch, 0.4f);
        ImGui::TableSetupColumn("Gas on (GC+VC)", ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableSetupColumn("Gas off (VC only)", ImGuiTableColumnFlags_WidthStretch, 0.22f);
        ImGui::TableSetupColumn("Difference [%]", ImGuiTableColumnFlags_WidthStretch, 0.16f);
        ImGui::TableHeadersRow();
        const auto row = [&](const char* name, double gui::KeyResults::*field, gui::UnitFamily family,
                             const char* id, double factor, const char* format) {
            const double va = a ? a->*field : kNaN, vb = b ? b->*field : kNaN;
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextUnformatted(name);
            ImGui::PopTextWrapPos();
            for (double v : {va, vb}) {
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(id ? value_text(v, family, units, id).c_str() : plain_text(v * factor, format).c_str());
            }
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(plain_text(gui::percent_difference(vb, va), "%+.1f").c_str());
        };
        using K = gui::KeyResults;
        row("Minimum film thickness", &K::h_min, gui::UnitFamily::SmallLength, kUnitFilm, 1, nullptr);
        row("Eccentricity ratio", &K::eccentricity_ratio, gui::UnitFamily::Fixed, nullptr, 1, "%.3f");
        row("Attitude angle [deg]", &K::attitude_deg, gui::UnitFamily::Fixed, nullptr, 1, "%.1f");
        row("Friction torque", &K::friction_torque, gui::UnitFamily::Torque, kUnitTorque, 1, nullptr);
        row("Friction power loss", &K::friction_power, gui::UnitFamily::Power, kUnitPower, 1, nullptr);
        row("Peak pressure", &K::p_max, gui::UnitFamily::Pressure, kUnitPressure, 1, nullptr);
        row("Lowest film pressure", &K::p_min, gui::UnitFamily::Pressure, kUnitPressure, 1, nullptr);
        row("Cavitated area (void fraction > 0) [%]", &K::cavitated_fraction, gui::UnitFamily::Fixed, nullptr, 100,
            "%.1f");
        row("Shear share of film force [%]", &K::shear_share, gui::UnitFamily::Fixed, nullptr, 100, "%.1f");
        ImGui::EndTable();
    }
    note("Difference: (gas on - gas off) / gas off.");
    const bool fbns_on = on.record && gui::result_form(on.record->job.config) == gui::ResultForm::Fbns;
    if (fbns_on)
        note("Pressure-plateau model: the cavity is liquid streaks plus a massless void; gas is not tracked. "
             "Gas on and gas off differ only in the cavity pressure (p_sat or the oil vapour pressure).");
    const auto& pu = unit_of(units, kUnitPressure, gui::UnitFamily::Pressure);
    const std::string p_axis = std::string("Mid-plane pressure [") + pu.label + "]";
    const float height = std::max(240.0f, ImGui::GetContentRegionAvail().y - ImGui::GetTextLineHeightWithSpacing());
    ImGui::BeginGroup();
    if (ImPlot::BeginPlot("##comparison_profiles", ImVec2(-1, height))) {
        bool any_run = false;
        for (int k = 0; k < 2; ++k)
            any_run = any_run || (sides[k].record && sides[k].run && sides[k].run->result.grid);
        ImPlot::SetupAxes(tab.x_label.c_str(), p_axis.c_str(), 0, ui::fit_if(any_run));
        ImPlot::SetupAxisLimits(ImAxis_X1, tab.start, tab.start + 360.0, ImPlotCond_Always);
        if (!any_run) ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.0, ImPlotCond_Always);
        for (int k = 0; k < 2; ++k) {
            const ComparisonSide& side = sides[k];
            if (!side.record || !side.run || !side.run->result.grid) continue;
            const CompareResult& c = side.run->result;
            const double run_zero = zero_position(angle_zero, side.record->job.model.feed.angle_deg, c.results);
            const Profile profile = make_profile(*c.grid, side.record->job.config, run_zero, tab.start, 0.5);
            const auto p = scaled(profile.p, pu);
            ImPlot::SetNextLineStyle(k == 0 ? theme::kModelGcVc : theme::kModelVc, 2.0f);
            if (!p.empty())
                ImPlot::PlotLine(k == 0 ? kGasOnLabel : kGasOffLabel, profile.angle.data(), p.data(),
                                 static_cast<int>(p.size()));
        }
        const double sat = on.record ? on.record->job.saturation_pressure : tab.p_sat;
        if (sat > 0.0) {
            const double v = gui::to_display(pu, sat);
            ImPlot::SetNextLineStyle(theme::kHighlight, 1.0f);
            ImPlot::PlotInfLines("Saturation pressure", &v, 1, ImPlotInfLinesFlags_Horizontal);
        }
        ImPlot::EndPlot();
    }
    ImGui::EndGroup();
    remember_plot(capture, "comparison_profiles");
    for (int k = 0; k < 2; ++k) {
        const ComparisonSide& side = sides[k];
        if (side.record && side.run && !side.run->loading && !side.run->result.error.empty()) {
            const std::string text = std::string(k == 0 ? kGasOnLabel : kGasOffLabel) + ": " + side.run->result.error;
            note(text.c_str());
        }
    }
}

}  // namespace result_data
