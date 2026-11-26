// Solution panel: the three model choices (dissolved-gas release, thermal,
// comparison mode); Advanced: gas model, convergence, time stepping, solver
// form, hidden settings and the exact command line.

#include "gui_app.hpp"
#include "imgui.h"

#include <cstdio>
#include <string>

namespace {

/// Convergence table row with a double tolerance; returns true on edit.
bool tolerance_row(const char* label, const char* meaning, const char* id, double& value, double default_value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    ui::item_tooltip(meaning);
    ImGui::TableSetColumnIndex(1);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    double edited = value;
    const bool changed = ImGui::InputDouble("##v", &edited, 0.0, 0.0, "%.3g") && edited > 0.0;
    if (changed) value = edited;
    ImGui::PopID();
    ImGui::TableSetColumnIndex(2);
    ImGui::TextDisabled("%s (default %.0e)", meaning, default_value);
    return changed;
}

/// Convergence table row with an iteration count; returns true on edit.
bool iterations_row(const char* label, const char* meaning, const char* id, int& value, int default_value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    ui::item_tooltip(meaning);
    ImGui::TableSetColumnIndex(1);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    int edited = value;
    const bool changed = ImGui::InputInt("##i", &edited, 0, 0) && edited >= 1;
    if (changed) value = edited;
    ImGui::PopID();
    ImGui::TableSetColumnIndex(2);
    ImGui::TextDisabled("%s (default %d)", meaning, default_value);
    return changed;
}

}  // namespace

/// Model choices; Advanced: gas model, convergence, time stepping, solver form, command line.
void GuiApp::draw_solution_panel() {
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Solution)];
    const bool transient = case_.analysis == gui::AnalysisType::Transient;
    bool changed = false;

    if (ui::begin_form("##solution")) {
        changed |= ui::form_checkbox(
            "Dissolved-gas release (gaseous cavitation)",
            "On (GC+VC): refrigerant comes out of solution, so the film cannot fall below the saturation "
            "pressure. Off (VC only): the refrigerant stays dissolved and the film cavitates only at the oil "
            "vapour pressure.",
            "sol.gas", case_.dissolved_gas);
        ImGui::SameLine();
        ImGui::TextDisabled(case_.dissolved_gas ? "on (GC+VC): film held at p_sat or above"
                                                : "off (VC only): refrigerant stays dissolved");
        changed |= ui::form_checkbox("Thermal analysis (film, housing and shaft)",
                                     "Solve the film, housing and shaft temperatures instead of a uniform ambient "
                                     "temperature. Steady cases only.",
                                     "sol.thermal", case_.thermal);
        ImGui::SameLine();
        ImGui::TextDisabled(case_.thermal ? "on" : "off (isothermal at ambient)");
        ui::form_checkbox("Comparison mode: run with and without gas release",
                          "Run launches two copies of the case, gas release on and off, and shows the "
                          "differences in %.",
                          "sol.comparison", comparison_mode_);

        if (advanced) {
            const bool forced = case_.thermal && gui::effective_content(case_, pair()) > 0.0;
            ImGui::BeginDisabled(forced || !case_.dissolved_gas);
            changed |= ui::form_radio<gui::GasModel>(
                "Gas model", "How the released refrigerant is modelled.", "sol.gas_model", case_.gas_model,
                {{gui::GasModel::Plateau, "Pressure plateau",
                  "FBNS: the cavitated film is held at the local saturation pressure. Fast; deck B0 and "
                  "rev03 runs."},
                 {gui::GasModel::Transport, "Gas transport",
                  "COUPLED: free gas is carried with the film and redissolves. Slower; deck X-runs and onset "
                  "chart. Steady only."}});
            ImGui::EndDisabled();
            if (forced) {
                ui::form_label("");
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextDisabled("Thermal analysis with refrigerant in the oil always uses gas transport "
                                    "(the plateau model's film energy needs pure oil).");
                ImGui::PopTextWrapPos();
            } else if (!case_.dissolved_gas) {
                ui::form_label("");
                ImGui::TextDisabled("Used only with dissolved-gas release on.");
            }
        }
        ui::end_form();
    }

    if (advanced) {
        const bool coupled = gui::uses_coupled(case_, pair());
        const gui::ConvergenceInput defaults;
        gui::ConvergenceInput& c = case_.convergence;
        ui::section("Convergence");
        if (ImGui::BeginTable("##convergence", 3, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("setting", ImGuiTableColumnFlags_WidthStretch, 0.34f);
            ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch, 0.18f);
            ImGui::TableSetupColumn("meaning", ImGuiTableColumnFlags_WidthStretch, 0.48f);
            if (!transient && case_.position == gui::ShaftPosition::LoadBalance)
                changed |= tolerance_row("Load balance tolerance", "|F + W| / |W|", "conv.force", c.force_tolerance,
                                         defaults.force_tolerance);
            changed |= tolerance_row("Outer tolerance", "pressure / temperature change per iteration", "conv.outer",
                                     c.outer_tolerance, defaults.outer_tolerance);
            changed |= iterations_row("Max outer iterations", "load balance / coupling iterations", "conv.outer_max",
                                      c.outer_max_iterations, defaults.outer_max_iterations);
            changed |= tolerance_row("Inner tolerance", "FBNS film equations residual", "conv.inner",
                                     c.inner_tolerance, defaults.inner_tolerance);
            changed |= iterations_row("Max inner iterations", "FBNS Newton iterations", "conv.inner_max",
                                      c.inner_max_iterations, defaults.inner_max_iterations);
            if (coupled)
                changed |= tolerance_row("Gas transport tolerance", "COUPLED nonlinear residual", "conv.transport",
                                         c.transport_tolerance, defaults.transport_tolerance);
            if (transient) {
                changed |= tolerance_row("Motion tolerance", "shaft position change per step", "conv.motion",
                                         c.motion_tolerance, defaults.motion_tolerance);
                changed |= iterations_row("Max motion iterations", "film force / motion per step", "conv.motion_max",
                                          c.motion_max_iterations, defaults.motion_max_iterations);
            }
            ImGui::EndTable();
        }

        if (transient) {
            ui::section("Time stepping");
            if (ui::begin_form("##time")) {
                changed |= ui::form_int("Steps per revolution", "Time steps per shaft revolution. Typical: 720 to 1440.",
                                        "time.steps", case_.steps_per_revolution, 36, 1000000);
                ui::form_label("Revolutions", "Shaft revolutions to simulate. Typical: 3.");
                ImGui::SetNextItemWidth(-FLT_MIN);
                double revolutions = case_.revolutions;
                if (ImGui::InputDouble("##revolutions", &revolutions, 0.0, 0.0, "%.4g") && revolutions > 0.0) {
                    case_.revolutions = revolutions;
                    changed = true;
                }
                changed |= ui::form_int("Field outputs per revolution",
                                        "Field snapshots written per revolution (views and export). Typical: 36.",
                                        "time.outputs", case_.outputs_per_revolution, 1, 100000);
                ui::end_form();
            }
        }

        ui::section("Solver");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Solver form: %s", gui::solver_form_text(case_, pair()).c_str());
        ImGui::PopTextWrapPos();
    }

    char header[96];
    std::snprintf(header, sizeof(header), "Settings in this file that the GUI does not show (kept on save): %zu###hidden",
                  hidden_keys_.size());
    if (ImGui::CollapsingHeader(header)) {
        if (hidden_keys_.empty()) ImGui::TextDisabled("None.");
        for (const std::string& key : hidden_keys_) ImGui::BulletText("%s", key.c_str());
    }

    ui::section("Command line");
    ImGui::PushTextWrapPos(0.0f);
    if (case_path_.empty()) {
        ImGui::TextDisabled("Save the case first to get the command line.");
    } else {
        const std::string command = command_line_for(case_path_, case_.cores);
        if (fonts.mono) ImGui::PushFont(fonts.mono);
        ImGui::TextUnformatted(command.c_str());
        if (fonts.mono) ImGui::PopFont();
        if (ImGui::Button("Copy command")) ImGui::SetClipboardText(command.c_str());
        ui::item_tooltip("Runs the same case from a terminal.");
    }
    ImGui::PopTextWrapPos();

    if (changed) mark_changed();
}
