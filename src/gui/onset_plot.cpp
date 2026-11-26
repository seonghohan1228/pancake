#include "onset_plot.hpp"

#include "implot.h"
#include "theme.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace ui {

namespace {

ImVec4 verdict_color(gui::OnsetVerdict verdict) {
    switch (verdict) {
        case gui::OnsetVerdict::Expected: return theme::kRiskHigh;
        case gui::OnsetVerdict::Borderline: return theme::kRiskMedium;
        case gui::OnsetVerdict::NotExpected: return theme::kRiskLow;
        case gui::OnsetVerdict::Unknown: break;
    }
    return theme::kMuted;
}

ImVec4 with_alpha(ImVec4 color, float alpha) {
    color.w = alpha;
    return color;
}

}  // namespace

void onset_chart(const char* id, ImVec2 size, const std::vector<OnsetMarker>& markers,
                 bool compact) {
    const gui::OnsetCurve& curve = gui::tester_onset_curve();

    // Axis window: the chart range, extended to keep every case point visible.
    double x_min = 0.01, x_max = 10.0, y_min = 0.05, y_max = 20.0;
    for (const OnsetMarker& marker : markers) {
        if (std::isfinite(marker.sommerfeld) && marker.sommerfeld > 0.0) {
            x_min = std::min(x_min, 0.5 * marker.sommerfeld);
            x_max = std::max(x_max, 2.0 * marker.sommerfeld);
        }
        if (std::isfinite(marker.margin_ratio) && marker.margin_ratio > 0.0) {
            y_min = std::min(y_min, 0.5 * marker.margin_ratio);
            y_max = std::max(y_max, 2.0 * marker.margin_ratio);
        }
    }

    ImPlotFlags plot_flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoMenus;
    if (compact) plot_flags |= ImPlotFlags_NoLegend | ImPlotFlags_NoTitle;
    if (!ImPlot::BeginPlot(id, size, plot_flags)) return;

    const ImPlotAxisFlags axis_flags = compact ? ImPlotAxisFlags_NoLabel : ImPlotAxisFlags_None;
    ImPlot::SetupAxes(compact ? nullptr : "Sommerfeld number S",
                      compact ? nullptr : "Pressure margin / specific load, (p_amb - p_sat) / P",
                      axis_flags, axis_flags);
    ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
    ImPlot::SetupAxesLimits(x_min, x_max, y_min, y_max, ImPlotCond_Always);
    if (!compact) ImPlot::SetupLegend(ImPlotLocation_NorthEast);

    std::vector<double> xs, ys;
    xs.reserve(curve.points.size());
    ys.reserve(curve.points.size());
    for (const gui::OnsetPoint& point : curve.points) {
        xs.push_back(point.sommerfeld);
        ys.push_back(point.pressure_drop_ratio);
    }
    const int count = static_cast<int>(xs.size());

    if (count >= 2) {
        // Below the curve the margin is smaller than the film pressure drop: gas release.
        ImPlot::SetNextFillStyle(theme::kRiskHigh, 0.12f);
        ImPlot::PlotShaded("Gas release region", xs.data(), ys.data(), count, -INFINITY);
        ImPlot::SetNextLineStyle(theme::kMuted, 2.0f);
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, 3.5f, theme::kMuted, 1.0f, theme::kMuted);
        ImPlot::PlotLine("Onset curve Phi(S)", xs.data(), ys.data(), count);

        const double s_low = 0.1, s_high = 1.0;
        const double phi_low = curve.phi(s_low), phi_high = curve.phi(s_high);
        if (std::isfinite(phi_low))
            ImPlot::PlotText("gas release", s_low, phi_low * 0.3);
        if (std::isfinite(phi_high))
            ImPlot::PlotText("no gas release", s_high, phi_high * 4.0);
    }

    // Case points, one scatter per verdict so the legend names the categories.
    const gui::OnsetVerdict verdicts[] = {gui::OnsetVerdict::Expected, gui::OnsetVerdict::Borderline,
                                          gui::OnsetVerdict::NotExpected,
                                          gui::OnsetVerdict::Unknown};
    const float marker_size = compact ? 4.5f : 6.0f;
    for (gui::OnsetVerdict verdict : verdicts) {
        std::vector<double> mx, my;
        std::vector<double> below_x, below_y;  // margin <= 0: drawn at the bottom edge
        for (const OnsetMarker& marker : markers) {
            if (marker.verdict != verdict) continue;
            if (!std::isfinite(marker.sommerfeld) || marker.sommerfeld <= 0.0 ||
                !std::isfinite(marker.margin_ratio))
                continue;
            if (marker.margin_ratio > 0.0) {
                mx.push_back(marker.sommerfeld);
                my.push_back(marker.margin_ratio);
            } else {
                below_x.push_back(marker.sommerfeld);
                below_y.push_back(y_min * 1.15);
            }
        }
        if (mx.empty() && below_x.empty()) continue;
        const ImVec4 color = verdict_color(verdict);
        const char* name = gui::onset_verdict_text(verdict);
        if (!mx.empty()) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, marker_size, color, 1.5f,
                                       with_alpha(color, 1.0f));
            ImPlot::PlotScatter(name, mx.data(), my.data(), static_cast<int>(mx.size()));
        }
        if (!below_x.empty()) {
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Down, marker_size + 1.0f, color, 1.5f, color);
            ImPlot::PlotScatter(name, below_x.data(), below_y.data(),
                                static_cast<int>(below_x.size()));
        }
    }

    // Marker labels, offset to the upper right of each point.
    if (!compact || markers.size() <= 3) {
        for (const OnsetMarker& marker : markers) {
            if (marker.label.empty() || !std::isfinite(marker.sommerfeld) || marker.sommerfeld <= 0.0 ||
                !std::isfinite(marker.margin_ratio))
                continue;
            const bool below = marker.margin_ratio <= 0.0;
            const double y = below ? y_min * 1.15 : marker.margin_ratio;
            std::string text = marker.label;
            if (below) text += " (margin <= 0)";
            ImPlot::PushStyleColor(ImPlotCol_InlayText, verdict_color(marker.verdict));
            ImPlot::PlotText(text.c_str(), marker.sommerfeld, y, ImVec2(6.0f, -10.0f),
                             ImPlotTextFlags_None);
            ImPlot::PopStyleColor();
        }
    }

    if (!compact) {
        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 plot_size = ImPlot::GetPlotSize();
        const char* caption = "approximate for other bearings";
        const ImVec2 text_size = ImGui::CalcTextSize(caption);
        ImPlot::GetPlotDrawList()->AddText(
            ImVec2(pos.x + plot_size.x - text_size.x - 6.0f, pos.y + plot_size.y - text_size.y - 4.0f),
            ImGui::GetColorU32(theme::kMuted), caption);
    }

    if (ImPlot::IsPlotHovered()) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted("Gaseous-cavitation onset chart");
        ImGui::TextDisabled("%s", curve.source.c_str());
        ImGui::TextUnformatted(
            "Gas release is expected when the case point lies below the curve: the pressure "
            "margin is smaller than the film pressure drop of a gas-free run. The curve is "
            "approximate for other bearings; the run gives the direct check (film p_min vs p_sat).");
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }

    ImPlot::EndPlot();
}

}  // namespace ui
