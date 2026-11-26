// Workflow window (run bar, workflow tree with status lights, estimates box)
// and the "Case setup" window that hosts the panel of the selected node.

#include "gui_app.hpp"
#include "imgui.h"
#include "onset_plot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr gui::Node kNodes[] = {gui::Node::Geometry,   gui::Node::Fluid, gui::Node::Operating,
                                gui::Node::Boundaries, gui::Node::Mesh,  gui::Node::Solution,
                                gui::Node::Results};

/// True when the value can be shown (finite and, for ratios, not a placeholder 0).
bool shown(double value) { return std::isfinite(value); }

std::string number_text(double value, const char* format = "%.3g") {
    if (!shown(value)) return "-";
    char buffer[48];
    std::snprintf(buffer, sizeof(buffer), format, value);
    return buffer;
}

/// Colour of the onset verdict (design risk: gas release expected is red).
ImVec4 verdict_color(gui::OnsetVerdict verdict) {
    switch (verdict) {
        case gui::OnsetVerdict::Expected: return theme::kRiskHigh;
        case gui::OnsetVerdict::Borderline: return theme::kRiskMedium;
        case gui::OnsetVerdict::NotExpected: return theme::kRiskLow;
        case gui::OnsetVerdict::Unknown: break;
    }
    return theme::kMuted;
}

/// One row of the estimates table: label (with tooltip) | value.
void estimate_row(const char* label, const char* help, const std::string& value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::TextUnformatted(label);
    ui::item_tooltip(help);
    ImGui::TableSetColumnIndex(1);
    ImGui::TextUnformatted(value.c_str());
    ui::item_tooltip(help);
}

/// All reasons of a node, one per line, for the tree tooltip.
std::string joined_reasons(const gui::NodeStatus& status) {
    std::string text;
    for (const std::string& reason : status.reasons) {
        if (!text.empty()) text += "\n";
        text += reason;
    }
    return text;
}

}  // namespace

/// Window "Workflow": run bar, the seven workflow nodes with their lights, estimates.
void GuiApp::draw_workflow_window() {
    if (!ImGui::Begin("Workflow")) {
        ImGui::End();
        return;
    }

    // Run bar
    const bool running = runner_.running();
    const bool can_run = ready_to_run() && !running;
    std::string run_reason;
    for (int i = 0; i < static_cast<int>(gui::Node::Results) && run_reason.empty(); ++i) {
        const gui::NodeStatus& node = status_[static_cast<std::size_t>(i)];
        if (node.light == gui::Light::Red || node.light == gui::Light::Yellow)
            run_reason = std::string(gui::node_name(static_cast<gui::Node>(i))) + ": " +
                         (node.reasons.empty() ? std::string("needs attention.") : node.reasons.front());
    }
    if (running) run_reason = "A run is in progress.";
    if (run_reason.empty()) run_reason = "Run the case (F5). The case is saved first.";

    ImGui::BeginDisabled(!can_run);
    ImGui::PushStyleColor(ImGuiCol_Button, theme::kOk);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.14f, 0.62f, 0.34f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.08f, 0.44f, 0.24f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
    if (ImGui::Button("Run (F5)")) start_run();
    ImGui::PopStyleColor(4);
    ImGui::EndDisabled();
    ui::item_tooltip(run_reason.c_str());

    ImGui::SameLine();
    ImGui::BeginDisabled(!running);
    if (ImGui::Button("Cancel")) cancel_all();
    ImGui::EndDisabled();
    ui::item_tooltip("Stop the running case and drop any queued runs.");

    ImGui::SameLine();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("CPU cores");
    const char* cores_help =
        "Number of processes the solver runs on. More cores = faster; above 1 the solver "
        "needs Microsoft MPI (mpiexec) installed. Typical: 1 to 8.";
    ui::item_tooltip(cores_help);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 5.5f);
    int cores = std::max(1, case_.cores);
    if (ImGui::InputInt("##cores", &cores, 1, 4)) {
        cores = std::clamp(cores, 1, 64);
        if (cores != case_.cores) {
            case_.cores = cores;
            mark_changed();
        }
    }
    ui::item_tooltip(cores_help);

    // Workflow tree
    ui::section("Workflow");
    for (gui::Node node : kNodes) {
        const gui::NodeStatus& status = status_[static_cast<std::size_t>(node)];
        ImGui::PushID(static_cast<int>(node));
        ui::status_dot(status.light);
        ImGui::SameLine();
        if (ImGui::Selectable(gui::node_name(node), selected_node_ == node)) {
            selected_node_ = node;
            if (node == gui::Node::Results) ImGui::SetWindowFocus("Results");
        }
        const std::string tip = joined_reasons(status);
        ui::item_tooltip(tip.empty() ? (status.light == gui::Light::Green ? "Complete." : "")
                                     : tip.c_str());
        ImGui::PopID();
    }

    draw_estimates_box();
    ImGui::End();
}

/// Instant estimates before the run: load, Sommerfeld, saturation, margin, onset verdict.
void GuiApp::draw_estimates_box() {
    ui::section("Estimates (before the run)");
    const gui::Estimates& e = estimates_;
    const bool have_load = e.specific_load > 0.0 && shown(e.specific_load);
    const bool have_fluid = pair() != nullptr;

    if (ImGui::BeginTable("##estimates", 2, ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.5f);
        ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch, 0.5f);
        estimate_row("Specific load P", "Load per projected area, P = W / (L D).",
                     have_load ? ui::format_value(e.specific_load, gui::UnitFamily::Pressure, units(),
                                                  "est.specific_load")
                               : "-");
        estimate_row("Sommerfeld number S",
                     "Bearing duty parameter S = mu N L D (R/c)^2 / W, N in rev/s, mu at ambient "
                     "temperature and the content used.",
                     have_load && e.sommerfeld > 0.0 ? number_text(e.sommerfeld) : "-");
        estimate_row("Refrigerant content", "Refrigerant content used by the case (saturated or given).",
                     have_fluid ? ui::format_value(e.content, gui::UnitFamily::Content, units(),
                                                   "fluid.content")
                                : "-");
        estimate_row("Saturation pressure p_sat",
                     "Pressure below which dissolved refrigerant comes out of solution, at the "
                     "content used and the ambient temperature.",
                     have_fluid ? ui::format_value(e.saturation_pressure, gui::UnitFamily::Pressure,
                                                   units(), "est.p_sat")
                                : "-");
        estimate_row("Pressure margin", "Ambient minus saturation pressure, p_amb - p_sat.",
                     have_fluid ? ui::format_value(e.pressure_margin, gui::UnitFamily::Pressure,
                                                   units(), "est.margin")
                                : "-");
        estimate_row("Margin / P", "Pressure margin per specific load, (p_amb - p_sat) / P.",
                     have_fluid && have_load ? number_text(e.margin_ratio) : "-");
        estimate_row("Onset curve Phi(S)",
                     "Pressure drop ratio Phi = (p_amb - p_min) / P of gas-free runs at this S. "
                     "Gas release starts when margin / P falls below Phi(S).",
                     e.onset_phi > 0.0 ? number_text(e.onset_phi) : "-");
        ImGui::EndTable();
    }

    ImGui::TextColored(verdict_color(e.verdict), "%s", gui::onset_verdict_text(e.verdict));
    ui::item_tooltip(
        "From the onset chart of the tester bearing (L/D = 0.93, axial groove). For other "
        "bearings it is an estimate; the run gives the direct check (film p_min vs p_sat). "
        "Gas release is a design risk.");
    ImGui::PushTextWrapPos(0.0f);
    if (!e.note.empty()) ImGui::TextDisabled("%s", e.note.c_str());
    if (!e.fluid_in_range)
        ImGui::TextColored(theme::kRiskMedium,
                           "The ambient state is outside the fluid table; values are clamped.");
    ImGui::PopTextWrapPos();
    if (e.run_time.typical_s > 0.0 && shown(e.run_time.typical_s)) {
        const double minutes = std::max(1.0, std::round(e.run_time.typical_s / 60.0));
        ImGui::TextDisabled("Run time: about %.0f min (rough)", minutes);
        ui::item_tooltip("Cells x typical iterations x time per cell-iteration, divided by the cores.");
    }

    std::vector<ui::OnsetMarker> markers;
    if (have_load && have_fluid && e.sommerfeld > 0.0 && shown(e.margin_ratio)) {
        ui::OnsetMarker marker;
        marker.label = "this case";
        marker.sommerfeld = e.sommerfeld;
        marker.margin_ratio = e.margin_ratio;
        marker.verdict = e.verdict;
        markers.push_back(marker);
    }
    ui::onset_chart("##mini", ImVec2(-1, 160), markers, true);
}

/// Window "Case setup": heading, Advanced toggle, node reasons, then the node's panel.
void GuiApp::draw_panel_window() {
    if (!ImGui::Begin("Case setup")) {
        ImGui::End();
        return;
    }
    const auto index = static_cast<std::size_t>(selected_node_);
    if (fonts.heading) ImGui::PushFont(fonts.heading);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(gui::node_name(selected_node_));
    if (fonts.heading) ImGui::PopFont();
    if (selected_node_ != gui::Node::Results) ui::advanced_toggle("##advanced", advanced_[index]);

    const gui::NodeStatus& status = status_[index];
    if (status.light != gui::Light::Green && status.light != gui::Light::Idle) {
        ImGui::PushTextWrapPos(0.0f);
        for (const std::string& reason : status.reasons)
            ImGui::TextColored(ui::light_color(status.light), "%s", reason.c_str());
        ImGui::PopTextWrapPos();
    }
    ImGui::Separator();

    switch (selected_node_) {
        case gui::Node::Geometry: draw_geometry_panel(); break;
        case gui::Node::Fluid: draw_fluid_panel(); break;
        case gui::Node::Operating: draw_operating_panel(); break;
        case gui::Node::Boundaries: draw_boundaries_panel(); break;
        case gui::Node::Mesh: draw_mesh_panel(); break;
        case gui::Node::Solution: draw_solution_panel(); break;
        case gui::Node::Results: {
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextUnformatted(
                "Results are shown in the Results window: key results, field maps, polar plot, "
                "profiles and export.");
            ImGui::PopTextWrapPos();
            if (ImGui::Button("Open the Results window")) ImGui::SetWindowFocus("Results");
            if (shown_record_ < 0 || shown_record_ >= static_cast<int>(records_.size())) {
                ImGui::TextDisabled("No run yet in this session.");
                break;
            }
            const RunRecord& record = records_[static_cast<std::size_t>(shown_record_)];
            ImGui::TextDisabled("Run \"%s\" started %s, %.0f s", record.job.label.c_str(),
                                record.started_at.c_str(), record.duration_s);
            if (!record.results_ready) {
                ImGui::TextDisabled("%s", record.results_error.empty() ? "Reading results..."
                                                                       : record.results_error.c_str());
                break;
            }
            const gui::KeyResults& r = record.results;
            if (ui::begin_form("##results_summary")) {
                ui::form_value("Minimum film thickness h_min", "Smallest gap.",
                               shown(r.h_min) ? ui::format_value(r.h_min, gui::UnitFamily::SmallLength,
                                                                 units(), "res.h_min")
                                              : "-");
                ui::form_value("Eccentricity ratio", "Shaft centre offset per clearance, e / c.",
                               number_text(r.eccentricity_ratio, "%.3f"));
                ui::form_value("Friction torque", "Viscous torque on the shaft.",
                               shown(r.friction_torque)
                                   ? ui::format_value(r.friction_torque, gui::UnitFamily::Torque,
                                                      units(), "res.torque")
                                   : "-");
                ui::form_value("Friction power loss", "Torque times angular speed.",
                               shown(r.friction_power)
                                   ? ui::format_value(r.friction_power, gui::UnitFamily::Power,
                                                      units(), "res.power")
                                   : "-");
                ui::form_value("Peak pressure", "Highest film pressure.",
                               shown(r.p_max) ? ui::format_value(r.p_max, gui::UnitFamily::Pressure,
                                                                 units(), "res.p_max")
                                              : "-");
                ui::form_value("Cavitated area", "Fraction of the film with film content below 1.",
                               shown(r.cavitated_fraction)
                                   ? number_text(100.0 * r.cavitated_fraction, "%.1f %%")
                                   : "-");
                ui::form_value("Film thickness ratio Lambda", "h_min / combined roughness.",
                               number_text(r.lambda, "%.2f"));
                ui::end_form();
            }
            ImGui::TextColored(record.converged ? theme::kRiskLow : theme::kRiskHigh, "%s",
                               record.converged ? "Converged." : "Not converged: results are not reliable.");
            break;
        }
    }
    ImGui::End();
}
