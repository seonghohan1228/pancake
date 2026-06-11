// Results, Convergence, Solver Log and Run History panels. All plots are
// ImPlot instruments: pan/zoom, hover readouts, labeled axes with units,
// CSV and PNG export.

#include "gui_app.hpp"
#include "implot.h"
#include "win_util.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>

namespace {

struct GridStats {
    const void* grid_key = nullptr;
    std::string field;
    double min = 0, max = 0, mean = 0;
    bool valid = false;
};

GridStats compute_stats(const ResultGrid& grid, const std::string& field_name) {
    GridStats stats;
    stats.grid_key = &grid;
    stats.field = field_name;
    const auto* field = grid.field(field_name);
    if (!field || field->empty()) return stats;
    double sum = 0;
    size_t count = 0;
    stats.min = DBL_MAX;
    stats.max = -DBL_MAX;
    for (double value : *field) {
        if (!std::isfinite(value)) continue;
        stats.min = std::min(stats.min, value);
        stats.max = std::max(stats.max, value);
        sum += value;
        ++count;
    }
    if (count > 0) {
        stats.mean = sum / static_cast<double>(count);
        stats.valid = true;
    }
    return stats;
}

const char* field_unit(const std::string& name) {
    if (name == "p") return "Pa";
    if (name == "h") return "m";
    if (name == "T") return "K";
    if (name == "rho" || name == "rho_liquid_solution") return "kg/m^3";
    if (name == "mu" || name == "mu_liquid_solution") return "Pa.s";
    if (name == "film_content" || name == "alpha_gas" || name == "inlet_indicator") return "-";
    if (name == "heat_generation") return "W/m^3";
    if (name == "dissolved_gas") return "kg/kg";
    if (name == "gas_mass_transfer") return "kg/(s.m^3)";
    if (name.rfind("F", 0) == 0) return "N";
    if (name.rfind("xB", 0) == 0) return "m";
    if (name.rfind("U", 0) == 0) return "m/s";
    return "";
}

void copyable_value(const char* label, double value, const char* unit, const char* fmt = "%.6g") {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    ImGui::TableSetColumnIndex(1);
    char text[64];
    if (std::isfinite(value)) {
        std::snprintf(text, sizeof(text), fmt, value);
    } else {
        std::snprintf(text, sizeof(text), "-");
    }
    ImGui::TextUnformatted(text);
    ImGui::SameLine();
    ImGui::TextDisabled("%s", unit);
    ImGui::TableSetColumnIndex(2);
    ImGui::PushID(label);
    if (std::isfinite(value) && ImGui::SmallButton("Copy")) {
        ImGui::SetClipboardText(text);
    }
    ImGui::PopID();
}

}  // namespace

// ---------------------------------------------------------------------------
// Results window
// ---------------------------------------------------------------------------

void GuiApp::draw_results_panel() {
    if (!ImGui::Begin("Results")) {
        ImGui::End();
        return;
    }

    // Toolbar: step selection + refresh.
    if (ImGui::Button("Refresh")) request_scan(resolve_output_dir());
    ImGui::SetItemTooltip("Rescan the output directory for saved steps.");
    ImGui::SameLine();
    ImGui::BeginDisabled(steps_.empty());
    const int current_index = selected_step_ < 0 ? static_cast<int>(steps_.size()) - 1
                                                 : selected_step_;
    ImGui::BeginDisabled(current_index <= 0);
    if (ImGui::ArrowButton("##prev_step", ImGuiDir_Left)) select_step(current_index - 1);
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 13);
    std::string step_preview = "no steps";
    if (!steps_.empty() && current_index >= 0) {
        const StepInfo& info = steps_[static_cast<size_t>(current_index)];
        step_preview = "step " + std::to_string(info.step) + "  (t = " +
                       std::to_string(info.time).substr(0, 8) + " s)";
        if (selected_step_ < 0) step_preview = "latest: " + step_preview;
    }
    if (ImGui::BeginCombo("##step", step_preview.c_str())) {
        if (ImGui::Selectable("latest", selected_step_ < 0)) select_step(-1);
        for (int i = 0; i < static_cast<int>(steps_.size()); ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "step %d  (t = %.6g s)", steps_[i].step,
                          steps_[i].time);
            if (ImGui::Selectable(label, selected_step_ == i)) select_step(i);
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(steps_.empty() || (selected_step_ < 0));
    if (ImGui::ArrowButton("##next_step", ImGuiDir_Right)) select_step(current_index + 1);
    ImGui::EndDisabled();
    ImGui::EndDisabled();
    if (loader_.busy()) {
        ImGui::SameLine();
        ImGui::TextDisabled("loading...");
    }

    if (!results_error_.empty() && !grid_) {
        ImGui::TextColored(theme::kWarn, "%s", results_error_.c_str());
        if (progress_.state == SolverState::Idle) {
            ImGui::TextDisabled("Run a case or use Refresh after pointing the output directory "
                                "at an existing results folder.");
        }
    } else if (!results_error_.empty()) {
        ImGui::TextColored(theme::kWarn, "Result file problem: %s", results_error_.c_str());
    }

    if (ImGui::BeginTabBar("results_tabs")) {
        if (ImGui::BeginTabItem("Field map")) {
            draw_heatmap_tab();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Profiles")) {
            draw_profiles_tab();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Summary")) {
            draw_summary_tab();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

void GuiApp::draw_heatmap_tab() {
    if (!grid_) {
        ImGui::TextDisabled("No result step loaded yet.");
        return;
    }
    const ResultGrid& grid = *grid_;

    // Field picker over the scalar arrays present in the file.
    std::vector<std::string> field_names;
    for (const auto& [name, values] : grid.fields) {
        if (name == "theta_rad" || name == "z_m") continue;
        field_names.push_back(name);
    }
    if (field_names.empty()) {
        ImGui::TextColored(theme::kWarn, "The result file contains no cell fields.");
        return;
    }
    if (!grid.field(selected_field_)) selected_field_ = field_names.front();

    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12);
    if (ImGui::BeginCombo("Field", selected_field_.c_str())) {
        for (const std::string& name : field_names) {
            std::string label = name;
            const char* unit = field_unit(name);
            if (*unit) label += std::string("  [") + unit + "]";
            if (ImGui::Selectable(label.c_str(), name == selected_field_)) {
                selected_field_ = name;
            }
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    static bool show_cavitation = true;
    ImGui::Checkbox("Cavitation zone", &show_cavitation);
    ImGui::SetItemTooltip("Hatch cells where the film content is below 1 (ruptured film).");
    ImGui::SameLine();
    if (ImGui::SmallButton("CSV")) export_grid_csv();
    ImGui::SameLine();
    bool want_png = ImGui::SmallButton("PNG");

    const auto* field = grid.field(selected_field_);
    if (!field) return;
    GridStats stats = compute_stats(grid, selected_field_);
    if (!stats.valid) {
        ImGui::TextColored(theme::kWarn, "Field '%s' has no finite values.", selected_field_.c_str());
        return;
    }
    double scale_min = stats.min;
    double scale_max = stats.max;
    if (scale_max <= scale_min) scale_max = scale_min + 1.0;

    // Row-flip so z increases upward; NaN cells pinned to the scale minimum.
    static std::vector<double> display;
    display.assign(field->size(), scale_min);
    for (int j = 0; j < grid.ny; ++j) {
        const size_t source_row = static_cast<size_t>(j) * grid.nx;
        const size_t target_row = static_cast<size_t>(grid.ny - 1 - j) * grid.nx;
        for (int i = 0; i < grid.nx; ++i) {
            const double value = (*field)[source_row + i];
            display[target_row + i] = std::isfinite(value) ? value : scale_min;
        }
    }

    const char* unit = field_unit(selected_field_);
    char value_label[96];
    std::snprintf(value_label, sizeof(value_label), "%s%s%s", selected_field_.c_str(),
                  *unit ? " " : "", *unit ? (std::string("[") + unit + "]").c_str() : "");

    ImGui::TextDisabled("step %d  t = %.6g s   min %.6g   max %.6g   mean %.6g %s   (%d x %d cells)",
                        grid.step, grid.time, stats.min, stats.max, stats.mean, unit, grid.nx,
                        grid.ny);

    const float scale_width = ImGui::GetFontSize() * 6.5f;
    const ImVec2 plot_size(ImGui::GetContentRegionAvail().x - scale_width,
                           ImGui::GetContentRegionAvail().y);
    ImVec2 capture_min(0, 0), capture_max(0, 0);

    // Outer edge of the last cell row (cell centers + half spacing).
    const double z_max = grid.z_m.empty() ? 1.0 : grid.z_m.back() + grid.z_m.front();

    ImPlot::PushColormap(ImPlotColormap_Viridis);
    if (ImPlot::BeginPlot("##fieldmap", plot_size, ImPlotFlags_NoLegend)) {
        ImPlot::SetupAxes("bearing angle [deg]", "axial position z [m]");
        ImPlot::SetupAxesLimits(0, 360, 0, z_max, ImPlotCond_Once);
        ImPlot::PlotHeatmap("##values", display.data(), grid.ny, grid.nx, scale_min, scale_max,
                            nullptr, ImPlotPoint(0, 0), ImPlotPoint(360, z_max));

        if (show_cavitation) {
            const auto* film_content = grid.field("film_content");
            if (film_content) {
                static std::vector<double> mask;
                mask.assign(film_content->size(), 0.0);
                bool any = false;
                for (int j = 0; j < grid.ny; ++j) {
                    const size_t source_row = static_cast<size_t>(j) * grid.nx;
                    const size_t target_row = static_cast<size_t>(grid.ny - 1 - j) * grid.nx;
                    for (int i = 0; i < grid.nx; ++i) {
                        const double content = (*film_content)[source_row + i];
                        if (std::isfinite(content) && content < 0.999) {
                            mask[target_row + i] = 1.0;
                            any = true;
                        }
                    }
                }
                if (any) {
                    static ImPlotColormap mask_colormap = -1;
                    if (mask_colormap == -1) {
                        const ImU32 colors[2] = {IM_COL32(0, 0, 0, 0),
                                                 IM_COL32(255, 255, 255, 110)};
                        mask_colormap = ImPlot::AddColormap("##cav_mask", colors, 2, false);
                    }
                    ImPlot::PushColormap(mask_colormap);
                    ImPlot::PlotHeatmap("##cavitation", mask.data(), grid.ny, grid.nx, 0.0, 1.0,
                                        nullptr, ImPlotPoint(0, 0), ImPlotPoint(360, z_max));
                    ImPlot::PopColormap();
                }
            } else {
                ImPlot::Annotation(180, z_max * 0.5, ImVec4(1, 1, 1, 0.6), ImVec2(0, 0), true,
                                   "film_content not in output_fields - cavitation overlay "
                                   "unavailable");
            }
        }

        // Hover readout: angle, z, cell value.
        if (ImPlot::IsPlotHovered()) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            if (mouse.x >= 0 && mouse.x <= 360 && mouse.y >= 0 && mouse.y <= z_max) {
                const int i = std::clamp(static_cast<int>(mouse.x / 360.0 * grid.nx), 0,
                                         grid.nx - 1);
                const int j = std::clamp(static_cast<int>(mouse.y / z_max * grid.ny), 0,
                                         grid.ny - 1);
                const double value = (*field)[static_cast<size_t>(j) * grid.nx + i];
                ImGui::BeginTooltip();
                ImGui::Text("theta = %.1f deg, z = %.4g m", grid.theta_deg[i], grid.z_m[j]);
                ImGui::Text("%s = %.6g %s", selected_field_.c_str(), value, unit);
                const auto* film_content = grid.field("film_content");
                if (film_content) {
                    const double content = (*film_content)[static_cast<size_t>(j) * grid.nx + i];
                    if (std::isfinite(content) && content < 0.999) {
                        ImGui::TextColored(theme::kWarn, "cavitated (film content %.3f)", content);
                    }
                }
                ImGui::EndTooltip();
            }
        }
        capture_min = ImPlot::GetPlotPos();
        const ImVec2 size = ImPlot::GetPlotSize();
        capture_max = ImVec2(capture_min.x + size.x, capture_min.y + size.y);
        ImPlot::EndPlot();
    }
    ImGui::SameLine();
    ImPlot::ColormapScale(value_label, scale_min, scale_max,
                          ImVec2(scale_width - ImGui::GetStyle().ItemSpacing.x, plot_size.y));
    ImPlot::PopColormap();

    if (want_png && capture_max.x > capture_min.x) {
        // Include the colorbar in the capture.
        const std::string name = selected_field_ + "_step" + std::to_string(grid.step) + ".png";
        request_plot_capture(name.c_str(), capture_min.x, capture_min.y,
                             capture_max.x + scale_width, capture_max.y);
    }
}

void GuiApp::draw_profiles_tab() {
    if (!grid_) {
        ImGui::TextDisabled("No result step loaded yet.");
        return;
    }
    const ResultGrid& grid = *grid_;
    const auto* pressure = grid.field("p");
    const auto* thickness = grid.field("h");
    if (!pressure && !thickness) {
        ImGui::TextColored(theme::kWarn, "Profiles need 'p' or 'h' in output_fields.");
        return;
    }

    const int mid = grid.ny / 2;
    int j = profile_z_index_ < 0 ? mid : std::clamp(profile_z_index_, 0, grid.ny - 1);
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14);
    if (ImGui::SliderInt("z slice", &j, 0, grid.ny - 1,
                         ("z = " + std::to_string(grid.z_m[static_cast<size_t>(j)])).c_str())) {
        profile_z_index_ = j;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("Mid-plane")) profile_z_index_ = -1;
    ImGui::SameLine();
    if (ImGui::SmallButton("CSV")) export_profiles_csv();
    ImGui::SameLine();
    const bool want_png = ImGui::SmallButton("PNG");
    int overlay_count = 0;
    for (const RunRecord& record : history_) {
        if (record.overlay && !record.profile_p.empty()) ++overlay_count;
    }
    if (overlay_count > 0) {
        ImGui::SameLine();
        ImGui::TextDisabled("overlaying %d run%s (mid-plane curves)", overlay_count,
                            overlay_count == 1 ? "" : "s");
    }

    static std::vector<double> p_row, h_row;
    if (pressure) {
        p_row.resize(grid.nx);
        for (int i = 0; i < grid.nx; ++i) {
            p_row[i] = (*pressure)[static_cast<size_t>(j) * grid.nx + i];
        }
    }
    if (thickness) {
        h_row.resize(grid.nx);
        for (int i = 0; i < grid.nx; ++i) {
            h_row[i] = (*thickness)[static_cast<size_t>(j) * grid.nx + i];
        }
    }

    ImVec2 capture_min(FLT_MAX, FLT_MAX), capture_max(-FLT_MAX, -FLT_MAX);
    const auto merge_capture = [&]() {
        const ImVec2 pos = ImPlot::GetPlotPos();
        const ImVec2 size = ImPlot::GetPlotSize();
        capture_min.x = std::min(capture_min.x, pos.x);
        capture_min.y = std::min(capture_min.y, pos.y);
        capture_max.x = std::max(capture_max.x, pos.x + size.x);
        capture_max.y = std::max(capture_max.y, pos.y + size.y);
    };

    if (ImPlot::BeginSubplots("##profiles", 2, 1, ImGui::GetContentRegionAvail(),
                              ImPlotSubplotFlags_LinkAllX)) {
        if (ImPlot::BeginPlot("Pressure")) {
            ImPlot::SetupAxes("bearing angle [deg]", "p [Pa]");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, 360, ImPlotCond_Once);
            if (pressure) {
                ImPlot::PlotLine("current", grid.theta_deg.data(), p_row.data(), grid.nx);
            }
            for (const RunRecord& record : history_) {
                if (!record.overlay || record.profile_p.empty()) continue;
                const std::string label = "run " + std::to_string(record.id);
                ImPlot::PlotLine(label.c_str(), record.profile_theta_deg.data(),
                                 record.profile_p.data(),
                                 static_cast<int>(record.profile_p.size()));
            }
            merge_capture();
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Film thickness")) {
            ImPlot::SetupAxes("bearing angle [deg]", "h [m]");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, 360, ImPlotCond_Once);
            if (thickness) {
                ImPlot::PlotLine("current", grid.theta_deg.data(), h_row.data(), grid.nx);
            }
            for (const RunRecord& record : history_) {
                if (!record.overlay || record.profile_h.empty()) continue;
                const std::string label = "run " + std::to_string(record.id);
                ImPlot::PlotLine(label.c_str(), record.profile_theta_deg.data(),
                                 record.profile_h.data(),
                                 static_cast<int>(record.profile_h.size()));
            }
            merge_capture();
            ImPlot::EndPlot();
        }
        ImPlot::EndSubplots();
    }

    if (want_png && capture_max.x > capture_min.x) {
        const std::string name = "profiles_step" + std::to_string(grid.step) + ".png";
        request_plot_capture(name.c_str(), capture_min.x, capture_min.y, capture_max.x,
                             capture_max.y);
    }
}

void GuiApp::draw_summary_tab() {
    if (progress_.state == SolverState::Failed) {
        ImGui::TextColored(theme::kError,
                           "The last run failed (exit code %d). Check the Solver Log tab - the "
                           "last lines usually name the cause (config validation, file "
                           "permissions, or a diverged solve).",
                           progress_.exit_code);
        ImGui::Separator();
    }
    if (!grid_) {
        ImGui::TextDisabled("Key results appear here after a run (or Refresh).");
        return;
    }
    const ResultGrid& grid = *grid_;
    const GridStats h_stats = compute_stats(grid, "h");
    const GridStats p_stats = compute_stats(grid, "p");
    const GridStats film_stats = compute_stats(grid, "film_content");
    const GridStats temperature_stats = compute_stats(grid, "T");

    double cavitated_fraction = std::numeric_limits<double>::quiet_NaN();
    if (const auto* film_content = grid.field("film_content")) {
        size_t cavitated = 0, total = 0;
        for (double value : *film_content) {
            if (!std::isfinite(value)) continue;
            ++total;
            if (value < 0.999) ++cavitated;
        }
        if (total > 0) cavitated_fraction = static_cast<double>(cavitated) / total;
    }

    const RunRecord* last_record = history_.empty() ? nullptr : &history_.back();

    // Prefer the solver's own resultant fields for this step; fall back to the
    // last run record (DETAILED-stdout parse) when the fields are not written.
    const auto first_cell = [&](const char* name) {
        const auto* values = grid.field(name);
        return (values && !values->empty()) ? (*values)[0]
                                            : std::numeric_limits<double>::quiet_NaN();
    };
    double torque = first_cell("friction_torque");
    if (!std::isfinite(torque) && last_record) torque = last_record->friction_torque;
    double attitude = first_cell("bearing_attitude_angle");
    if (!std::isfinite(attitude) && last_record) attitude = last_record->attitude_deg;
    const double omega = last_record ? last_record->config.omega : config_.omega;
    const double power =
        std::isfinite(torque) ? std::abs(torque * omega)
                              : std::numeric_limits<double>::quiet_NaN();
    const double force_x = first_cell("F.x");
    const double force_y = first_cell("F.y");

    ImGui::PushFont(fonts.heading);
    ImGui::Text("Key results - step %d (t = %.6g s)", grid.step, grid.time);
    ImGui::PopFont();
    if (ImGui::BeginTable("##summary", 3,
                          ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("quantity", ImGuiTableColumnFlags_WidthStretch, 0.5f);
        ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch, 0.4f);
        ImGui::TableSetupColumn("copy", ImGuiTableColumnFlags_WidthStretch, 0.1f);
        copyable_value("Minimum film thickness", h_stats.valid ? h_stats.min : NAN, "m");
        copyable_value("Peak pressure", p_stats.valid ? p_stats.max : NAN, "Pa");
        copyable_value("Cavitated area fraction", cavitated_fraction, "-", "%.4f");
        copyable_value("Attitude angle", attitude, "deg", "%.2f");
        copyable_value("Friction torque", torque, "N.m");
        copyable_value("Friction power loss", power, "W");
        if (std::isfinite(force_x) && std::isfinite(force_y)) {
            copyable_value("Fluid force magnitude", std::hypot(force_x, force_y), "N");
        }
        if (last_record) {
            copyable_value("Axial boundary mass flux", last_record->boundary_flux, "kg/s");
        }
        if (temperature_stats.valid) {
            copyable_value("Max temperature", temperature_stats.max, "K", "%.2f");
        }
        if (film_stats.valid) {
            copyable_value("Min film content", film_stats.min, "-", "%.4f");
        }
        ImGui::EndTable();
    }
    if (!std::isfinite(torque)) {
        ImGui::TextDisabled("Friction torque needs 'friction_torque' in output_fields or "
                            "Detailed terminal verbosity.");
    }
    ImGui::Separator();
    ImGui::TextDisabled("Full fields: Results > Field map, or open the run in ParaView "
                        "(Results menu).");
}

// ---------------------------------------------------------------------------
// Convergence window
// ---------------------------------------------------------------------------

void GuiApp::draw_convergence_panel() {
    if (!ImGui::Begin("Convergence")) {
        ImGui::End();
        return;
    }
    if (diagnostics_.empty()) {
        if (config_.diagnostics_interval <= 0) {
            ImGui::TextColored(theme::kWarn,
                               "diagnostics_interval is 0 - enable it (Numerics, advanced) to "
                               "get residual traces.");
        } else {
            ImGui::TextDisabled("Residuals stream here once a solve is running.");
        }
        ImGui::End();
        return;
    }

    if (ImGui::SmallButton("CSV")) export_convergence_csv();
    ImGui::SameLine();
    const bool want_png = ImGui::SmallButton("PNG");
    ImGui::SameLine();
    ImGui::TextDisabled("%zu diagnostic rows", diagnostics_.size());

    static std::vector<double> steps, residuals, theta_changes, cavitated;
    steps.clear(); residuals.clear(); theta_changes.clear(); cavitated.clear();
    for (const DiagRow& row : diagnostics_) {
        steps.push_back(row.step);
        residuals.push_back(std::max(std::abs(row.liquid_residual), 1e-300));
        theta_changes.push_back(std::max(row.max_theta_change, 1e-300));
        cavitated.push_back(row.cavitated_fraction);
    }

    ImVec2 capture_min(0, 0), capture_max(0, 0);
    if (ImPlot::BeginPlot("##convergence", ImGui::GetContentRegionAvail())) {
        ImPlot::SetupAxes("step", "residual (log)");
        ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
        ImPlot::SetupAxis(ImAxis_Y2, "cavitated fraction [-]",
                          ImPlotAxisFlags_AuxDefault);
        ImPlot::SetupAxisLimits(ImAxis_Y2, 0, 1, ImPlotCond_Once);
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
        ImPlot::PlotLine("liquid mass residual", steps.data(), residuals.data(),
                         static_cast<int>(steps.size()));
        ImPlot::PlotLine("max film-content change", steps.data(), theta_changes.data(),
                         static_cast<int>(steps.size()));
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
        ImPlot::PlotLine("cavitated fraction", steps.data(), cavitated.data(),
                         static_cast<int>(steps.size()));
        capture_min = ImPlot::GetPlotPos();
        const ImVec2 size = ImPlot::GetPlotSize();
        capture_max = ImVec2(capture_min.x + size.x, capture_min.y + size.y);
        ImPlot::EndPlot();
    }
    if (want_png && capture_max.x > capture_min.x) {
        request_plot_capture("convergence.png", capture_min.x, capture_min.y, capture_max.x,
                             capture_max.y);
    }
    ImGui::End();
}

// ---------------------------------------------------------------------------
// Solver Log window
// ---------------------------------------------------------------------------

void GuiApp::draw_log_panel() {
    if (!ImGui::Begin("Solver Log")) {
        ImGui::End();
        return;
    }
    if (ImGui::SmallButton("Clear")) log_.clear();
    ImGui::SameLine();
    if (ImGui::SmallButton("Copy all")) {
        std::string text;
        for (const LogLine& line : log_) {
            text += line.text;
            text += '\n';
        }
        ImGui::SetClipboardText(text.c_str());
    }
    ImGui::SameLine();
    ImGui::Checkbox("Auto-scroll", &log_autoscroll_);
    ImGui::SameLine();
    ImGui::TextDisabled("%zu lines", log_.size());

    ImGui::BeginChild("##log_scroll", ImVec2(0, 0), ImGuiChildFlags_None,
                      ImGuiWindowFlags_HorizontalScrollbar);
    ImGui::PushFont(fonts.mono);
    ImGuiListClipper clipper;
    clipper.Begin(static_cast<int>(log_.size()));
    while (clipper.Step()) {
        for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; ++i) {
            const LogLine& line = log_[static_cast<size_t>(i)];
            if (line.color != 0) {
                ImGui::PushStyleColor(ImGuiCol_Text, line.color);
                ImGui::TextUnformatted(line.text.c_str());
                ImGui::PopStyleColor();
            } else {
                ImGui::TextUnformatted(line.text.c_str());
            }
        }
    }
    if (log_autoscroll_ && ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 4.0f) {
        ImGui::SetScrollHereY(1.0f);
    }
    ImGui::PopFont();
    ImGui::EndChild();
    ImGui::End();
}

// ---------------------------------------------------------------------------
// Run History window
// ---------------------------------------------------------------------------

void GuiApp::draw_history_panel() {
    if (!ImGui::Begin("Run History")) {
        ImGui::End();
        return;
    }
    if (history_.empty()) {
        ImGui::TextDisabled("Each solve adds a row here for comparison.\n"
                            "Tick 'overlay' to compare runs on the Profiles plots.");
        ImGui::End();
        return;
    }
    if (ImGui::SmallButton("Export CSV")) export_history_csv();
    ImGui::SameLine();
    ImGui::TextDisabled("%zu runs this session", history_.size());

    const ImGuiTableFlags flags = ImGuiTableFlags_Sortable | ImGuiTableFlags_RowBg |
                                  ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY |
                                  ImGuiTableFlags_ScrollX | ImGuiTableFlags_SizingFixedFit;
    if (ImGui::BeginTable("##history", 11, flags)) {
        ImGui::TableSetupScrollFreeze(1, 1);
        ImGui::TableSetupColumn("run", ImGuiTableColumnFlags_DefaultSort);
        ImGui::TableSetupColumn("outcome");
        ImGui::TableSetupColumn("duration [s]");
        ImGui::TableSetupColumn("omega [rad/s]");
        ImGui::TableSetupColumn("e [m]");
        ImGui::TableSetupColumn("load [N]");
        ImGui::TableSetupColumn("mesh");
        ImGui::TableSetupColumn("h_min [m]");
        ImGui::TableSetupColumn("p_max [Pa]");
        ImGui::TableSetupColumn("overlay", ImGuiTableColumnFlags_NoSort);
        ImGui::TableSetupColumn("actions", ImGuiTableColumnFlags_NoSort);
        ImGui::TableHeadersRow();

        // Sort a copy of the indices per the table sort specs.
        std::vector<int> order(history_.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
        if (ImGuiTableSortSpecs* sort_specs = ImGui::TableGetSortSpecs()) {
            if (sort_specs->SpecsCount > 0) {
                const ImGuiTableColumnSortSpecs& spec = sort_specs->Specs[0];
                const auto key = [&](int index) -> double {
                    const RunRecord& record = history_[static_cast<size_t>(index)];
                    switch (spec.ColumnIndex) {
                        case 0: return record.id;
                        case 1: return static_cast<double>(record.outcome);
                        case 2: return record.duration_s;
                        case 3: return record.config.omega;
                        case 4: return record.config.e;
                        case 5: return std::hypot(record.config.external_load_x,
                                                  record.config.external_load_y);
                        case 6: return record.config.n_theta_global * 1e6 +
                                       record.config.n_z_global;
                        case 7: return record.h_min;
                        case 8: return record.p_max;
                        default: return record.id;
                    }
                };
                std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                    const double ka = key(a), kb = key(b);
                    return spec.SortDirection == ImGuiSortDirection_Ascending ? ka < kb : ka > kb;
                });
            }
        }

        for (int index : order) {
            RunRecord& record = history_[static_cast<size_t>(index)];
            ImGui::PushID(record.id);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d  (%s)", record.id, record.started_at.c_str());
            ImGui::TableSetColumnIndex(1);
            const ImVec4 color = record.outcome == SolverState::Finished ? theme::kOk
                                 : record.outcome == SolverState::Cancelled ? theme::kWarn
                                                                            : theme::kError;
            ImGui::TextColored(color, "%s", to_string(record.outcome));
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.1f", record.duration_s);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%.4g", record.config.omega);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%.4g", record.config.e);
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%.4g", std::hypot(record.config.external_load_x,
                                           record.config.external_load_y));
            ImGui::TableSetColumnIndex(6);
            ImGui::Text("%dx%d", record.config.n_theta_global, record.config.n_z_global);
            ImGui::TableSetColumnIndex(7);
            if (std::isfinite(record.h_min)) ImGui::Text("%.4g", record.h_min);
            else ImGui::TextDisabled("-");
            ImGui::TableSetColumnIndex(8);
            if (std::isfinite(record.p_max)) ImGui::Text("%.4g", record.p_max);
            else ImGui::TextDisabled("-");
            ImGui::TableSetColumnIndex(9);
            ImGui::BeginDisabled(record.profile_p.empty());
            ImGui::Checkbox("##overlay", &record.overlay);
            ImGui::EndDisabled();
            if (record.profile_p.empty()) {
                ImGui::SetItemTooltip("No stored profile (run produced no readable results).");
            }
            ImGui::TableSetColumnIndex(10);
            if (ImGui::SmallButton("Load inputs")) {
                config_ = record.config;
                dirty_ = true;
                sync_raw_from_config();
                set_status("Inputs from run " + std::to_string(record.id) +
                           " loaded into the form.");
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }
    ImGui::End();
}
