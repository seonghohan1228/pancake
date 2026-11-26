// Live view window (tabs Geometry | Mesh) and the bearing section drawing.
//
// Fixed frame: x to the right, y up; angles from +x counterclockwise, which is
// the direction of rotation for positive speed. h = c + b_x cos(theta) +
// b_y sin(theta) with b the bearing displacement, so the journal approaches the
// bore at h_min. Everything is to scale apart from the clearance.

#include "gui_app.hpp"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

ImU32 color_u32(const ImVec4& color, float alpha = 1.0f) {
    return ImGui::GetColorU32(ImVec4(color.x, color.y, color.z, color.w * alpha));
}

/// Point at fixed-frame angle `deg` and radius `r` [px] around `centre` (y up).
ImVec2 polar(ImVec2 centre, float r, double deg) {
    return ImVec2(centre.x + r * static_cast<float>(std::cos(deg * kDeg)),
                  centre.y - r * static_cast<float>(std::sin(deg * kDeg)));
}

/// Straight arrow with a filled head at `to`.
void draw_arrow(ImDrawList* dl, ImVec2 from, ImVec2 to, ImU32 color, float thickness, float head) {
    const float dx = to.x - from.x, dy = to.y - from.y;
    const float length = std::sqrt(dx * dx + dy * dy);
    if (length < 1.0f) return;
    const float ux = dx / length, uy = dy / length;
    const ImVec2 base(to.x - ux * head, to.y - uy * head);
    dl->AddLine(from, base, color, thickness);
    dl->AddTriangleFilled(to, ImVec2(base.x - uy * head * 0.5f, base.y + ux * head * 0.5f),
                          ImVec2(base.x + uy * head * 0.5f, base.y - ux * head * 0.5f), color);
}

/// Text placed just outside a circle at fixed-frame angle `deg`.
void label_outside(ImDrawList* dl, ImVec2 anchor, double deg, ImU32 color, const char* text) {
    const ImVec2 size = ImGui::CalcTextSize(text);
    const float c = static_cast<float>(std::cos(deg * kDeg));
    const float s = static_cast<float>(std::sin(deg * kDeg));
    dl->AddText(ImVec2(anchor.x + (c - 1.0f) * 0.5f * size.x, anchor.y - (s + 1.0f) * 0.5f * size.y),
                color, text);
}

/// Largest of 1, 2, 5 x 10^k not above `value` (at least 1).
double nice_factor(double value) {
    if (!(value > 1.0)) return 1.0;
    const double decade = std::pow(10.0, std::floor(std::log10(value)));
    const double mantissa = value / decade;
    const double step = mantissa >= 5.0 ? 5.0 : mantissa >= 2.0 ? 2.0 : 1.0;
    return step * decade;
}

/// Splits `text` into lines no wider than `width` [px] (breaks at spaces).
std::vector<std::string> wrap_lines(const std::string& text, float width) {
    std::vector<std::string> lines;
    std::string line;
    std::size_t start = 0;
    while (start <= text.size()) {
        std::size_t end = text.find(' ', start);
        if (end == std::string::npos) end = text.size();
        const std::string word = text.substr(start, end - start);
        const std::string candidate = line.empty() ? word : line + " " + word;
        if (!line.empty() && ImGui::CalcTextSize(candidate.c_str()).x > width) {
            lines.push_back(line);
            line = word;
        } else {
            line = candidate;
        }
        start = end + 1;
    }
    if (!line.empty()) lines.push_back(line);
    return lines;
}

double wrap_deg(double deg) {
    double value = std::fmod(deg, 360.0);
    if (value < 0.0) value += 360.0;
    return value;
}

}  // namespace

void GuiApp::draw_live_view_window() {
    if (!ImGui::Begin("Live view")) {
        ImGui::End();
        return;
    }

    // live_view_tab_ may be set elsewhere (mesh panel, screenshots): select that tab
    // until the tab bar reports it, then follow the user's clicks again.
    static int synced_tab = -1;
    const bool external = live_view_tab_ != synced_tab;
    int active = -1;
    const float combo_width = ImGui::GetFontSize() * 13.0f;
    const float label_width = ImGui::CalcTextSize("Angle zero").x + ImGui::GetStyle().ItemSpacing.x;
    if (ImGui::BeginTable("##live_header", 2, ImGuiTableFlags_None)) {
        ImGui::TableSetupColumn("tabs", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("zero", ImGuiTableColumnFlags_WidthFixed, combo_width + label_width);
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        if (ImGui::BeginTabBar("##live_tabs")) {
            const char* names[] = {"Geometry", "Mesh"};
            const char* tips[] = {"Bearing section to scale, clearance exaggerated: feed, load, rotation, "
                                  "shaft position and the angle zeros.",
                                  "Unwrapped film (angle x z/L, true proportions) with the cells, the "
                                  "feed and the bearing ends."};
            for (int tab = 0; tab < 2; ++tab) {
                const ImGuiTabItemFlags flags =
                    external && live_view_tab_ == tab ? ImGuiTabItemFlags_SetSelected : 0;
                if (ImGui::BeginTabItem(names[tab], nullptr, flags)) {
                    active = tab;
                    ImGui::EndTabItem();
                }
                ui::item_tooltip(tips[tab]);
            }
            ImGui::EndTabBar();
        }
        ImGui::TableSetColumnIndex(1);
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Angle zero");
        ui::item_tooltip("Where angle 0 lies on every angle axis. Angles grow in the direction of "
                         "rotation. theta = 0 at the maximum film, theta_in = 0 at the feed.");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##angle_zero", gui::angle_zero_name(angle_zero_))) {
            for (int i = 0; i < 5; ++i) {
                const auto zero = static_cast<gui::AngleZero>(i);
                if (ImGui::Selectable(gui::angle_zero_name(zero), zero == angle_zero_)) {
                    angle_zero_ = zero;
                    settings.angle_zero = i;
                }
            }
            ImGui::EndCombo();
        }
        ImGui::EndTable();
    }
    if (active >= 0) {
        if (!external) live_view_tab_ = active;
        if (active == live_view_tab_) synced_tab = active;
    }

    const ImVec2 avail = ImGui::GetContentRegionAvail();
    if (live_view_tab_ == 1) draw_mesh_view(avail);
    else draw_geometry_view(avail);
    ImGui::End();
}

void GuiApp::draw_geometry_view(ImVec2 size) {
    size.x = std::max(size.x, 120.0f);
    size.y = std::max(size.y, 120.0f);
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##geometry_canvas", size);
    const ImVec2 p1(p0.x + size.x, p0.y + size.y);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(p0, p1, true);
    dl->AddRectFilled(p0, p1, IM_COL32(252, 252, 253, 255));

    const float font = ImGui::GetFontSize();
    const ImU32 ink = ImGui::GetColorU32(ImGuiCol_Text);
    const ImU32 muted = color_u32(theme::kMuted);
    const ImU32 accent = color_u32(theme::kAccent);
    const ImU32 highlight = color_u32(theme::kHighlight);
    const ImU32 housing_fill = IM_COL32(214, 218, 225, 255);
    const ImU32 film_fill = IM_COL32(250, 234, 190, 255);
    const ImU32 shaft_fill = IM_COL32(160, 167, 178, 255);

    // Geometry in metres.
    const double radius = 0.5 * std::max(case_.bearing_diameter, 1e-6);
    const double clearance = std::max(case_.radial_clearance, 1e-9);
    const double outer_radius =
        case_.thermal ? (case_.housing.outer_diameter > 0.0 ? 0.5 * case_.housing.outer_diameter
                                                            : 1.1 * radius)
                      : 1.12 * radius;
    const double factor = nice_factor(0.15 * radius / clearance);
    const double clearance_shown = std::min(clearance * factor, 0.3 * radius);

    // Shaft position: last result, the fixed position, or a nominal estimate.
    double eps = 0.6;
    double hmin_deg = case_.load_direction_deg;
    double hmax_deg = hmin_deg + 180.0;
    std::string position_note = "Shaft position: estimate (0.6 c toward the load), no result yet";
    double load_deg = case_.load_direction_deg;
    bool load_known = case_.position == gui::ShaftPosition::LoadBalance || case_.analysis ==
                                                                                gui::AnalysisType::Transient;
    bool load_from_table = false;
    if (case_.analysis == gui::AnalysisType::Transient && !case_.load_table.file.empty()) {
        // The table is the load on the moving bearing; the shaft carries the opposite.
        double sum_x = 0.0, sum_y = 0.0;
        for (const LoadRow& row : load_rows_) {
            if (case_.load_table.format == gui::LoadFormat::MagnitudeDirection) {
                sum_x += row.a * std::cos(row.b * kDeg);
                sum_y += row.a * std::sin(row.b * kDeg);
            } else {
                sum_x += row.a;
                sum_y += row.b;
            }
        }
        sum_x *= case_.load_table.scale;
        sum_y *= case_.load_table.scale;
        load_known = std::hypot(sum_x, sum_y) > 0.0;
        if (load_known) load_deg = wrap_deg(std::atan2(sum_y, sum_x) / kDeg + 180.0);
        load_from_table = true;
        hmin_deg = load_deg;
        hmax_deg = load_deg + 180.0;
    }
    if (case_.analysis == gui::AnalysisType::Steady && case_.position == gui::ShaftPosition::Fixed) {
        eps = std::clamp(case_.fixed_eccentricity_ratio, 0.0, 0.99);
        hmin_deg = case_.fixed_attitude_deg;
        hmax_deg = hmin_deg + 180.0;
        position_note = "Shaft position: fixed by the input";
    }
    if (shown_record_ >= 0 && shown_record_ < static_cast<int>(records_.size())) {
        const RunRecord& record = records_[static_cast<std::size_t>(shown_record_)];
        const gui::KeyResults& r = record.results;
        if (record.results_ready && std::isfinite(r.eccentricity_ratio) && std::isfinite(r.h_min_theta_deg)) {
            eps = std::clamp(r.eccentricity_ratio, 0.0, 0.999);
            hmin_deg = r.h_min_theta_deg;
            hmax_deg = std::isfinite(r.h_max_theta_deg) ? r.h_max_theta_deg : hmin_deg + 180.0;
            char text[128];
            std::snprintf(text, sizeof(text), "Shaft position: from run \"%s\", eccentricity ratio %.3f",
                          record.job.label.c_str(), r.eccentricity_ratio);
            position_note = text;
        }
    }
    const double feed_deg = case_.feed.angle_deg;

    // Captions (top left), wrapped to the canvas width; the arrows are named here so
    // that no text is drawn on the shaft.
    struct Caption {
        std::string text;
        ImU32 color;
    };
    std::vector<Caption> captions;
    {
        char text[96];
        std::snprintf(text, sizeof(text), "; clearance x %.0f (exaggerated), the rest to scale", factor);
        captions.push_back({"D = " + ui::format_value(case_.bearing_diameter, gui::UnitFamily::Length, units(),
                                                      "geom.diameter") +
                                ", c = " + ui::format_value(case_.radial_clearance, gui::UnitFamily::SmallLength,
                                                            units(), "geom.clearance") +
                                text,
                            ink});
        captions.push_back({position_note, muted});
        std::string arrows = "Arrows: blue = rotation " +
                             ui::format_value(std::abs(case_.speed), gui::UnitFamily::Speed, units(), "op.speed", "%.0f");
        const bool fixed = case_.analysis == gui::AnalysisType::Steady && case_.position == gui::ShaftPosition::Fixed;
        if (load_known && !fixed)
            arrows += load_from_table ? std::string(", black = load on the shaft (mean of the load table)")
                                      : ", black = load on the shaft, W = " +
                                            ui::format_value(case_.load, gui::UnitFamily::Force, units(), "op.load",
                                                             "%.0f");
        captions.push_back({arrows, muted});
        if (highlighted_zone_ == 0) captions.push_back({"Bearing ends: see Mesh view", highlight});
    }
    using CaptionLines = std::vector<std::pair<std::string, ImU32>>;
    const auto wrap_captions = [&](float width) {
        CaptionLines lines;
        for (const Caption& caption : captions)
            for (std::string& line : wrap_lines(caption.text, std::max(width, 40.0f)))
                lines.emplace_back(std::move(line), caption.color);
        return lines;
    };
    const auto block_height = [&](const CaptionLines& lines) {
        return static_cast<float>(lines.size()) * 1.2f * font + 0.4f * font;
    };

    // Layout: captions above the drawing, or in a left column when that leaves a
    // larger drawing (wide, short canvas). The ring keeps room for its angle labels.
    const float label_w = ImGui::CalcTextSize("theta_in = 0 (feed)").x + 0.5f * font;
    CaptionLines top_lines = wrap_captions(size.x - font);
    const float top_h = block_height(top_lines);
    const float r_top = std::min(0.5f * size.x - label_w, 0.5f * (size.y - top_h) - 2.6f * font);
    const float column_w = 0.4f * size.x;
    CaptionLines side_lines = wrap_captions(column_w - font);
    const float side_h = block_height(side_lines);
    const float region_w = size.x - column_w;
    const float r_side = std::min(0.5f * region_w - label_w, 0.5f * size.y - 2.6f * font);
    const bool side = r_side > r_top && side_h < size.y - 3.5f * font;
    const CaptionLines caption_lines = side ? std::move(side_lines) : std::move(top_lines);
    const ImVec2 centre = side ? ImVec2(p0.x + column_w + 0.5f * region_w, p0.y + 0.5f * size.y)
                               : ImVec2(p0.x + 0.5f * size.x, p0.y + top_h + 0.5f * (size.y - top_h));
    const float avail_r = side ? r_side : r_top;
    const float scale = std::max(avail_r, 20.0f) / static_cast<float>(outer_radius);  // px per m
    const float r_bore = static_cast<float>(radius) * scale;
    const float r_outer = static_cast<float>(outer_radius) * scale;
    const float c_px = static_cast<float>(clearance_shown) * scale;
    const float r_shaft = r_bore - c_px;
    const ImVec2 shaft_centre = polar(centre, static_cast<float>(eps) * c_px, hmin_deg);
    const int segments = 128;

    // Bearing (or housing) ring, film and shaft.
    dl->AddCircleFilled(centre, r_outer, housing_fill, segments);
    if (case_.thermal) dl->AddCircle(centre, r_outer, ink, segments, 1.2f);
    dl->AddCircleFilled(centre, r_bore, film_fill, segments);

    // Oil feed (cut into the bore).
    const float band = std::max(0.04f * r_bore, 4.0f);
    if (case_.feed.type == gui::FeedType::Hole) {
        const float hole_r = std::max(0.5f * static_cast<float>(case_.feed.hole_diameter) * scale, 3.0f);
        const ImVec2 inner = polar(centre, r_bore - 1.0f, feed_deg);
        const ImVec2 outer = polar(centre, r_outer + 1.0f, feed_deg);
        dl->AddLine(inner, outer, film_fill, 2.0f * hole_r);
        dl->AddCircleFilled(inner, hole_r, film_fill, 24);
    } else {
        const double half = 0.5 * std::max(case_.feed.groove_width_deg, 0.5);
        dl->PathArcTo(centre, r_bore + 0.5f * band - 0.5f, static_cast<float>(-(feed_deg + half) * kDeg),
                      static_cast<float>(-(feed_deg - half) * kDeg), 32);
        dl->PathStroke(film_fill, ImDrawFlags_None, band + 1.0f);
    }
    dl->AddCircle(centre, r_bore, ink, segments, 1.2f);

    dl->AddCircleFilled(shaft_centre, r_shaft, shaft_fill, segments);
    dl->AddCircle(shaft_centre, r_shaft, ink, segments, 1.2f);

    // Zone highlights.
    const float thick = std::max(3.0f, 0.15f * font);
    switch (highlighted_zone_) {
        case 0: dl->AddCircle(centre, r_bore, color_u32(theme::kHighlight, 0.45f), segments, thick); break;
        case 1:
            if (case_.feed.type == gui::FeedType::Hole) {
                const float hole_r =
                    std::max(0.5f * static_cast<float>(case_.feed.hole_diameter) * scale, 3.0f);
                dl->AddCircle(polar(centre, r_bore, feed_deg), hole_r + 3.0f, highlight, 24, thick);
            } else {
                const double half = 0.5 * std::max(case_.feed.groove_width_deg, 0.5);
                dl->PathArcTo(centre, r_bore + 0.5f * band, static_cast<float>(-(feed_deg + half) * kDeg),
                              static_cast<float>(-(feed_deg - half) * kDeg), 32);
                dl->PathStroke(color_u32(theme::kHighlight, 0.8f), ImDrawFlags_None, band + 4.0f);
            }
            break;
        case 2: dl->AddCircle(shaft_centre, r_shaft, highlight, segments, thick); break;
        case 3: dl->AddCircle(centre, r_bore, highlight, segments, thick); break;
        case 4: dl->AddCircle(centre, r_outer, highlight, segments, thick); break;
        default: break;
    }

    // Centres and line of centres (towards h_min).
    dl->AddLine(centre, polar(centre, r_bore, hmin_deg), color_u32(theme::kMuted, 0.6f), 1.0f);
    dl->AddLine(ImVec2(centre.x - 4, centre.y), ImVec2(centre.x + 4, centre.y), ink, 1.0f);
    dl->AddLine(ImVec2(centre.x, centre.y - 4), ImVec2(centre.x, centre.y + 4), ink, 1.0f);
    dl->AddCircleFilled(shaft_centre, 2.5f, ink, 12);

    // Rotation arrow inside the shaft, opposite the load, pointing in the direction of rotation.
    const double rotation = case_.speed >= 0.0 ? 1.0 : -1.0;
    {
        const float r_arc = 0.72f * r_shaft;
        const double mid = load_deg + 180.0;
        const double from_deg = mid - rotation * 45.0, to_deg = mid + rotation * 45.0;
        dl->PathArcTo(shaft_centre, r_arc, static_cast<float>(-from_deg * kDeg),
                      static_cast<float>(-(to_deg - rotation * 6.0) * kDeg), 32);
        dl->PathStroke(accent, ImDrawFlags_None, 2.0f);
        const ImVec2 tip = polar(shaft_centre, r_arc, to_deg);
        const ImVec2 back = polar(shaft_centre, r_arc, to_deg - rotation * 10.0);
        draw_arrow(dl, back, tip, accent, 2.0f, 0.6f * font);
    }

    // Load arrow from the shaft centre.
    if (load_known) {
        const ImVec2 tip = polar(shaft_centre, 0.62f * r_shaft, load_deg);
        draw_arrow(dl, shaft_centre, tip, ink, 2.2f, 0.75f * font);
    }

    // Angle zeros and h_min, labelled outside the ring; nearby labels are stacked outwards.
    struct Mark {
        double deg;
        ImU32 color;
        const char* text;
    };
    const char* hmax_text = "theta = 0 (h_max)";
    const char* feed_text = "theta_in = 0 (feed)";
    const char* hmin_text = "h_min";
    const Mark marks[] = {{hmax_deg, accent, hmax_text},
                          {feed_deg, highlighted_zone_ == 1 ? highlight : ink, feed_text},
                          {hmin_deg, ink, hmin_text}};
    std::vector<double> placed;
    for (const Mark& mark : marks) {
        int stack = 0;
        for (double other : placed) {
            double gap = std::abs(wrap_deg(mark.deg - other));
            gap = std::min(gap, 360.0 - gap);
            if (gap < 18.0) ++stack;
        }
        placed.push_back(mark.deg);
        dl->AddLine(polar(centre, r_bore - 0.5f * c_px, mark.deg), polar(centre, r_outer + 6.0f, mark.deg),
                    mark.color, 1.5f);
        label_outside(dl, polar(centre, r_outer + 8.0f + stack * 1.15f * font, mark.deg), mark.deg,
                      mark.color, mark.text);
    }
    if (highlighted_zone_ == 4 && !case_.thermal)
        label_outside(dl, polar(centre, r_outer + 8.0f, 45.0), 45.0, highlight,
                      "housing: thermal analysis off");

    // Fixed-frame axes, bottom left.
    {
        const float arm = 1.8f * font;
        const ImVec2 origin(p0.x + 0.8f * font, p1.y - 0.8f * font);
        draw_arrow(dl, origin, ImVec2(origin.x + arm, origin.y), muted, 1.2f, 0.45f * font);
        draw_arrow(dl, origin, ImVec2(origin.x, origin.y - arm), muted, 1.2f, 0.45f * font);
        dl->AddText(ImVec2(origin.x + arm + 2.0f, origin.y - 0.55f * font), muted, "x");
        dl->AddText(ImVec2(origin.x - 0.3f * font, origin.y - arm - 1.1f * font), muted, "y");
    }

    // Captions, top left.
    {
        ImVec2 at(p0.x + 0.5f * font, p0.y + 0.3f * font);
        for (const auto& [line, color] : caption_lines) {
            dl->AddText(at, color, line.c_str());
            at.y += 1.2f * font;
        }
    }
    dl->PopClipRect();

    ui::item_tooltip("Bearing section seen along the shaft axis. Angles are measured from +x "
                     "counterclockwise; positive speed turns counterclockwise. The clearance is "
                     "exaggerated so the film and the shaft offset are visible.");
}
