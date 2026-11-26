// Mesh view: the unwrapped film (angle on the selected zero x z/L) in true
// proportions, with the cell boundaries, the feed and the bearing ends.

#include "gui_app.hpp"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

ImU32 color_u32(const ImVec4& color, float alpha = 1.0f) {
    return ImGui::GetColorU32(ImVec4(color.x, color.y, color.z, color.w * alpha));
}

/// "2nd", "3rd", "4th", ... for the thinning note.
std::string ordinal(int k) {
    const int tens = k % 100;
    const char* suffix = "th";
    if (tens < 11 || tens > 13) {
        if (k % 10 == 1) suffix = "st";
        else if (k % 10 == 2) suffix = "nd";
        else if (k % 10 == 3) suffix = "rd";
    }
    return std::to_string(k) + suffix;
}

/// Length in um below 1 mm, else in mm.
std::string format_cell(double metres) {
    char text[32];
    if (metres < 1e-3) std::snprintf(text, sizeof(text), "%.0f um", metres * 1e6);
    else std::snprintf(text, sizeof(text), "%.2f mm", metres * 1e3);
    return text;
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

}  // namespace

void GuiApp::draw_mesh_view(ImVec2 size) {
    size.x = std::max(size.x, 160.0f);
    size.y = std::max(size.y, 120.0f);
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##mesh_canvas", size);
    const ImVec2 p1(p0.x + size.x, p0.y + size.y);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->PushClipRect(p0, p1, true);
    dl->AddRectFilled(p0, p1, IM_COL32(252, 252, 253, 255));

    const float font = ImGui::GetFontSize();
    const ImU32 ink = ImGui::GetColorU32(ImGuiCol_Text);
    const ImU32 muted = color_u32(theme::kMuted);
    const ImU32 accent = color_u32(theme::kAccent);
    const ImU32 highlight = color_u32(theme::kHighlight);
    const ImU32 film_fill = IM_COL32(250, 240, 212, 255);
    const ImU32 feed_fill = IM_COL32(236, 196, 110, 255);

    int n_theta = case_.n_theta, n_z = case_.n_z;
    if (case_.mesh_level != gui::MeshLevel::Custom) gui::mesh_level_cells(case_.mesh_level, n_theta, n_z);
    n_theta = std::max(n_theta, 1);
    n_z = std::max(n_z, 1);
    const double diameter = std::max(case_.bearing_diameter, 1e-6);
    const double length = std::max(case_.bearing_length, 1e-6);
    const double circumference = kPi * diameter;

    // Moving zeros: the shown result, else the estimate of the geometry view.
    const double feed_deg = case_.feed.angle_deg;
    double hmin_deg = case_.load_direction_deg;
    if (case_.analysis == gui::AnalysisType::Transient && !case_.load_table.file.empty()) {
        double sum_x = 0.0, sum_y = 0.0;
        for (const LoadRow& row : load_rows_) {
            const bool polar = case_.load_table.format == gui::LoadFormat::MagnitudeDirection;
            sum_x += polar ? row.a * std::cos(row.b * kDeg) : row.a;
            sum_y += polar ? row.a * std::sin(row.b * kDeg) : row.b;
        }
        sum_x *= case_.load_table.scale;
        sum_y *= case_.load_table.scale;
        if (std::hypot(sum_x, sum_y) > 0.0) hmin_deg = std::atan2(sum_y, sum_x) / kDeg + 180.0;
    }
    if (case_.analysis == gui::AnalysisType::Steady && case_.position == gui::ShaftPosition::Fixed)
        hmin_deg = case_.fixed_attitude_deg;
    double hmax_deg = hmin_deg + 180.0;
    bool estimated = true;
    if (shown_record_ >= 0 && shown_record_ < static_cast<int>(records_.size())) {
        const RunRecord& record = records_[static_cast<std::size_t>(shown_record_)];
        if (record.results_ready && std::isfinite(record.results.h_min_theta_deg)) {
            hmin_deg = record.results.h_min_theta_deg;
            hmax_deg = std::isfinite(record.results.h_max_theta_deg) ? record.results.h_max_theta_deg
                                                                     : hmin_deg + 180.0;
            estimated = false;
        }
    }
    const double zero_deg = gui::angle_zero_position(angle_zero_, feed_deg, hmax_deg, hmin_deg);
    const double start = gui::angle_axis_start(angle_zero_);
    const double rotation = case_.speed >= 0.0 ? 1.0 : -1.0;

    // Cell sizes and the thinning are known before the layout; the info lines go on top.
    const double cell_theta = circumference / n_theta;
    const double cell_z = length / n_z;
    std::vector<std::pair<std::string, ImU32>> info;
    {
        char text[256];
        std::snprintf(text, sizeof(text), "%d x %d cells (%s), cell %s around x %s along", n_theta, n_z,
                      gui::mesh_level_name(case_.mesh_level), format_cell(cell_theta).c_str(),
                      format_cell(cell_z).c_str());
        std::string first = text;
        if (case_.feed.type == gui::FeedType::Hole) {
            std::snprintf(text, sizeof(text), "; hole spans %.0f cells",
                          std::max(1.0, std::round(case_.feed.hole_diameter / cell_theta)));
            first += text;
        }
        info.emplace_back(first, ink);
        std::string second = "Unwrapped film, true proportions (pi D : L).";
        if (estimated) second += " h_max, h_min: estimated, no result yet.";
        info.emplace_back(second, muted);
    }
    if (highlighted_zone_ == 0)
        info[1] = {"Bearing ends: the top and bottom edges (z = L, z = 0), highlighted.", highlight};
    // Info lines wrap to the canvas width; room is kept for the grid note (known after the layout).
    const float info_width = std::max(size.x - font, 40.0f);
    int info_line_count =
        static_cast<int>(wrap_lines("Grid: every 10th line around, every 10th line along shown.", info_width).size());
    for (const auto& entry : info) info_line_count += static_cast<int>(wrap_lines(entry.first, info_width).size());

    // Plot rectangle with true proportions (pi D : L).
    const float line = 1.25f * font;
    const float info_h = static_cast<float>(info_line_count) * line + 0.4f * font;
    const float top = p0.y + info_h;
    const float left = p0.x + 2.6f * font + ImGui::CalcTextSize("0.75").x;  // z/L title and ticks
    const float right = p1.x - 1.0f * font;
    const float bottom = p1.y - 2.4f * line;                 // angle ticks and title
    const float marker_row = 0.5f * font + 3.0f * ImGui::GetTextLineHeight();  // h_max / h_min / feed, up to 3 rows
    const float avail_w = std::max(right - left, 40.0f);
    const float avail_h = std::max(bottom - top - marker_row, 20.0f);
    const float aspect = static_cast<float>(circumference / length);
    float w = avail_w, h = avail_w / aspect;
    if (h > avail_h) {
        h = avail_h;
        w = h * aspect;
    }
    const float x0 = left + 0.5f * (avail_w - w);
    const float y0 = top + marker_row + 0.5f * (avail_h - h);
    const float x1 = x0 + w, y1 = y0 + h;
    const float px_per_m = w / static_cast<float>(circumference);
    const auto x_of = [&](double fixed_deg) {
        const double shown = gui::display_angle(fixed_deg, zero_deg, start, rotation);
        return x0 + static_cast<float>((shown - start) / 360.0) * w;
    };
    const auto y_of = [&](double z_over_l) { return y1 - static_cast<float>(z_over_l) * h; };

    dl->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, y1), film_fill);

    // Cell boundaries, thinned to keep lines at least ~4 px apart.
    const int k_theta = std::max(1, static_cast<int>(std::ceil(4.0f * n_theta / w)));
    const int k_z = std::max(1, static_cast<int>(std::ceil(4.0f * n_z / h)));
    const ImU32 grid = color_u32(theme::kMuted, 0.35f);
    for (int i = 0; i < n_theta; i += k_theta) {
        const float x = x_of(360.0 * i / n_theta);
        dl->AddLine(ImVec2(x, y0), ImVec2(x, y1), grid, 1.0f);
    }
    for (int j = 0; j <= n_z; j += k_z) {
        const float y = y_of(static_cast<double>(j) / n_z);
        dl->AddLine(ImVec2(x0, y), ImVec2(x1, y), grid, 1.0f);
    }

    // Oil feed to scale (wraps across the axis ends).
    dl->PushClipRect(ImVec2(x0, y0), ImVec2(x1, y1), true);
    const float feed_x = x_of(feed_deg);
    const bool feed_lit = highlighted_zone_ == 1;
    for (float shift : {-w, 0.0f, w}) {
        if (case_.feed.type == gui::FeedType::Hole) {
            const double z_c = case_.feed.hole_axial_position >= 0.0 ? case_.feed.hole_axial_position
                                                                      : 0.5 * length;
            const float r = std::max(0.5f * static_cast<float>(case_.feed.hole_diameter) * px_per_m, 2.5f);
            const ImVec2 c(feed_x + shift, y_of(z_c / length));
            dl->AddCircleFilled(c, r, feed_fill, 32);
            dl->AddCircle(c, r, feed_lit ? highlight : ink, 32, feed_lit ? 3.0f : 1.2f);
        } else {
            const float half = std::max(static_cast<float>(case_.feed.groove_width_deg / 360.0) * w, 1.5f) * 0.5f;
            const ImVec2 a(feed_x + shift - half, y0), b(feed_x + shift + half, y1);
            dl->AddRectFilled(a, b, feed_fill);
            dl->AddRect(a, b, feed_lit ? highlight : ink, 0.0f, 0, feed_lit ? 3.0f : 1.2f);
        }
    }
    dl->PopClipRect();

    // Outline; the bearing ends are the bottom (z = 0) and top (z = L) edges.
    dl->AddRect(ImVec2(x0, y0), ImVec2(x1, y1), ink, 0.0f, 0, 1.2f);
    if (highlighted_zone_ == 0) {
        dl->AddLine(ImVec2(x0, y0), ImVec2(x1, y0), highlight, 4.0f);
        dl->AddLine(ImVec2(x0, y1), ImVec2(x1, y1), highlight, 4.0f);
    }
    if (highlighted_zone_ == 2 || highlighted_zone_ == 3)
        dl->AddRect(ImVec2(x0, y0), ImVec2(x1, y1), color_u32(theme::kHighlight, 0.6f), 0.0f, 0, 3.0f);

    // h_max (theta = 0), h_min and the feed above the plot.
    struct Mark {
        double deg;
        ImU32 color;
        const char* text;
    };
    const Mark marks[] = {{hmax_deg, accent, "h_max"}, {hmin_deg, ink, "h_min"}, {feed_deg, ink, "feed"}};
    std::vector<std::pair<float, float>> taken[3];  // label x ranges per row
    for (const Mark& mark : marks) {
        const float x = x_of(mark.deg);
        dl->AddLine(ImVec2(x, y0 - 0.25f * font), ImVec2(x, y1), color_u32(theme::kAccent, 0.35f), 1.0f);
        dl->AddTriangleFilled(ImVec2(x - 4.0f, y0 - 0.35f * font), ImVec2(x + 4.0f, y0 - 0.35f * font),
                              ImVec2(x, y0 - 0.05f * font), mark.color);
        const ImVec2 text_size = ImGui::CalcTextSize(mark.text);
        const float tx = std::clamp(x - 0.5f * text_size.x, x0 - 0.5f * text_size.x, x1 - 0.5f * text_size.x);
        int row = 0;
        for (; row < 2; ++row) {
            bool free = true;
            for (const auto& range : taken[row])
                if (tx < range.second + 4.0f && tx + text_size.x > range.first - 4.0f) free = false;
            if (free) break;
        }
        taken[row].emplace_back(tx, tx + text_size.x);
        dl->AddText(ImVec2(tx, y0 - 0.4f * font - text_size.y * (1.0f + row)), mark.color, mark.text);
    }

    // Angle ticks every 60 deg, z/L ticks every 0.25.
    for (int t = 0; t <= 6; ++t) {
        const float x = x0 + w * t / 6.0f;
        dl->AddLine(ImVec2(x, y1), ImVec2(x, y1 + 4.0f), ink, 1.0f);
        char text[16];
        std::snprintf(text, sizeof(text), "%.0f", start + 60.0 * t);
        const ImVec2 text_size = ImGui::CalcTextSize(text);
        dl->AddText(ImVec2(x - 0.5f * text_size.x, y1 + 5.0f), ink, text);
    }
    {
        const std::string title = gui::angle_axis_label(angle_zero_);
        const ImVec2 text_size = ImGui::CalcTextSize(title.c_str());
        dl->AddText(ImVec2(0.5f * (x0 + x1) - 0.5f * text_size.x, y1 + 5.0f + line), ink, title.c_str());
    }
    for (int t = 0; t <= 4; ++t) {
        const float y = y_of(0.25 * t);
        dl->AddLine(ImVec2(x0 - 4.0f, y), ImVec2(x0, y), ink, 1.0f);
        char text[16];
        std::snprintf(text, sizeof(text), "%.2g", 0.25 * t);
        const ImVec2 text_size = ImGui::CalcTextSize(text);
        dl->AddText(ImVec2(x0 - 6.0f - text_size.x, y - 0.5f * text_size.y), ink, text);
    }
    dl->AddText(ImVec2(p0.x + 0.2f * font, 0.5f * (y0 + y1) - 0.5f * font), ink, "z/L");

    // Info lines.
    if (k_theta > 1 || k_z > 1) {
        std::string grid_note = "Grid: every " + (k_theta > 1 ? ordinal(k_theta) + " line" : std::string("line")) +
                                " around, every " + (k_z > 1 ? ordinal(k_z) + " line" : std::string("line")) +
                                " along shown.";
        info.emplace_back(grid_note, muted);
    }
    float info_y = p0.y + 0.3f * font;
    for (const auto& entry : info)
        for (const std::string& text : wrap_lines(entry.first, info_width)) {
            dl->AddText(ImVec2(p0.x + 0.5f * font, info_y), entry.second, text.c_str());
            info_y += line;
        }
    dl->PopClipRect();

    ui::item_tooltip("Film unwrapped around the bearing: angle along x (zero from the Angle zero "
                     "choice, growing in the direction of rotation), axial position z/L along y. "
                     "Bottom edge z = 0 and top edge z = L are the bearing ends.");
}
