// Solution panel: the three model choices (dissolved-gas release, thermal,
// comparison mode) and the run-time estimate; Advanced: gas model, cavity
// shear (read-only), convergence, time stepping, solver form, hidden settings
// and the exact command line.

#include "gui_app.hpp"
#include "imgui.h"

#include <cmath>
#include <cstdio>
#include <string>

namespace {

/// "40 s", "25 min" or "3.1 h".
std::string duration_text(double seconds) {
    char buffer[32];
    if (seconds < 90.0) std::snprintf(buffer, sizeof(buffer), "%.0f s", seconds);
    else if (seconds < 5400.0) std::snprintf(buffer, sizeof(buffer), "%.0f min", seconds / 60.0);
    else std::snprintf(buffer, sizeof(buffer), "%.1f h", seconds / 3600.0);
    return buffer;
}

/// Convergence table row with a double tolerance; returns true on edit.
bool tolerance_row(const char* label, const char* meaning, const char* id, double& value, double default_value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(label);
    ImGui::PopTextWrapPos();
    ui::item_tooltip(meaning);
    ImGui::TableSetColumnIndex(1);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    double edited = value;
    const bool changed = ImGui::InputDouble("##v", &edited, 0.0, 0.0, "%.3g") && edited > 0.0;
    if (changed) value = edited;
    ImGui::PopID();
    ImGui::TableSetColumnIndex(2);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s (default %.0e)", meaning, default_value);
    ImGui::PopTextWrapPos();
    return changed;
}

/// Convergence table row with an iteration count; returns true on edit.
bool iterations_row(const char* label, const char* meaning, const char* id, int& value, int default_value) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(label);
    ImGui::PopTextWrapPos();
    ui::item_tooltip(meaning);
    ImGui::TableSetColumnIndex(1);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    int edited = value;
    const bool changed = ImGui::InputInt("##i", &edited, 0, 0) && edited >= 1;
    if (changed) value = edited;
    ImGui::PopID();
    ImGui::TableSetColumnIndex(2);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s (default %d)", meaning, default_value);
    ImGui::PopTextWrapPos();
    return changed;
}

/// Wrapped disabled text in the current cell (no clipping at the column edge).
void wrapped_note(const char* text) {
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("%s", text);
    ImGui::PopTextWrapPos();
}

/// Form row label that wraps inside the label column (long labels are not cut off).
void wrapped_label(const char* label, const char* help) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextUnformatted(label);
    ImGui::PopTextWrapPos();
    ui::item_tooltip(help);
    ImGui::TableSetColumnIndex(1);
}

/// Form row: wrapped label + checkbox. Returns true on edit.
bool checkbox_row(const char* label, const char* help, const char* id, bool& value) {
    wrapped_label(label, help);
    ImGui::PushID(id);
    const bool changed = ImGui::Checkbox("##c", &value);
    ImGui::PopID();
    ui::item_tooltip(help);
    return changed;
}

/// One radio button on its own line; `disabled_reason` disables it and becomes the tooltip.
/// The label wraps inside the column and is clickable like the button.
bool gas_model_radio(const char* label, const char* help, bool active, const std::string& disabled_reason) {
    const bool disabled = !disabled_reason.empty();
    const char* tip = disabled ? disabled_reason.c_str() : help;
    ImGui::PushID(label);
    ImGui::BeginDisabled(disabled);
    bool clicked = ImGui::RadioButton("##radio", active);
    ImGui::EndDisabled();
    ui::item_tooltip(tip);
    ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
    ImGui::AlignTextToFramePadding();
    ImGui::PushTextWrapPos(0.0f);
    if (disabled) ImGui::TextDisabled("%s", label);
    else ImGui::TextUnformatted(label);
    ImGui::PopTextWrapPos();
    if (ImGui::IsItemClicked()) clicked = true;
    ui::item_tooltip(tip);
    ImGui::PopID();
    return clicked && !disabled;
}

const char* cavity_shear_text(CavitatedFilmModel model) {
    return model == CavitatedFilmModel::STRIATED ? "No shear from the bubble; only the liquid shears"
                                                 : "Whole gap shears (FULL_FILM)";
}

}  // namespace

/// Model choices and run time; Advanced: gas model, cavity shear, convergence, time stepping, solver form.
void GuiApp::draw_solution_panel() {
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Solution)];
    const bool transient = case_.analysis == gui::AnalysisType::Transient;
    bool changed = false;

    if (ui::begin_form("##solution")) {
        changed |= checkbox_row(
            "Dissolved-gas release (gaseous cavitation)",
            "On (GC+VC): refrigerant comes out of solution, so the film cannot fall below the saturation "
            "pressure. Off (VC only): the refrigerant stays dissolved and the film cavitates only at the oil "
            "vapour pressure.",
            "sol.gas", case_.dissolved_gas);
        ImGui::SameLine();
        wrapped_note(case_.dissolved_gas ? "on (GC+VC): film held at p_sat or above"
                                         : "off (VC only): refrigerant stays dissolved");
        changed |= checkbox_row("Thermal analysis (film, housing and shaft)",
                                     "Solve the film, housing and shaft temperatures instead of a uniform ambient "
                                     "temperature. Steady cases use the full thermal model, transient cases a "
                                     "simplified one.",
                                     "sol.thermal", case_.thermal);
        ImGui::SameLine();
        wrapped_note(case_.thermal ? "on" : "off (isothermal at ambient)");
        if (case_.thermal && transient) {
            ui::form_label("");
            wrapped_note("Simplified thermal model for transient cases: polynomial film profile, 2-D housing. "
                         "The finite-volume film and 3-D housing run in steady cases only.");
        }
        checkbox_row("Comparison mode: run with and without gas release",
                          "Run launches two copies of the case, gas release on and off, and shows the "
                          "differences in %.",
                          "sol.comparison", comparison_mode_);

        const gui::RunTime& run_time = estimates_.run_time;
        if (run_time.typical_s > 0.0 && std::isfinite(run_time.typical_s)) {
            std::string text = gui::format_run_time(run_time);
            if (transient && run_time.per_revolution_s > 0.0)
                text += "; " + duration_text(run_time.per_revolution_s) + " per revolution";
            char tip[384];
            std::snprintf(tip, sizeof(tip),
                          "Scaled from measured runs (%s). Typical core-time %.2g core-h (wall time x cores).",
                          run_time.basis.c_str(), run_time.core_hours);
            wrapped_label("Estimated run time", tip);
            wrapped_note(text.c_str());
            ui::item_tooltip(tip);
        }

        if (advanced) {
            const char* gas_off_reason = "Used only with dissolved-gas release on.";
            const bool no_gas_data = pair() && !pair()->has_gas_phase;
            const bool forced = gui::uses_coupled(case_, pair()) && case_.gas_model == gui::GasModel::Plateau;
            const bool transport_shown = case_.gas_model == gui::GasModel::Transport || forced;
            std::string transport_reason;
            if (no_gas_data) transport_reason = "This pair has no gas-phase data; gas transport is not available.";
            else if (forced)
                transport_reason = "A steady thermal case with refrigerant in the oil always uses gas transport "
                                   "(the plateau model's film energy needs pure oil).";
            else if (!case_.dissolved_gas) transport_reason = gas_off_reason;
            else if (transient && case_.gas_model != gui::GasModel::Transport)
                transport_reason = "Gas transport runs in steady cases only; transient cases use the pressure plateau.";
            std::string plateau_reason;
            if (forced) plateau_reason = transport_reason;
            else if (!case_.dissolved_gas) plateau_reason = gas_off_reason;
            wrapped_label("Gas model", "How the released refrigerant is modelled.");
            ImGui::PushID("sol.gas_model");
            if (gas_model_radio("Pressure plateau (fast; deck compressor form)",
                                "FBNS: the cavitated film is held at the local saturation pressure. In the "
                                "cavity the bubble carries no shear; only the liquid (fraction theta) does. "
                                "Fast; the form of the deck's compressor results.",
                                !transport_shown, plateau_reason)) {
                case_.gas_model = gui::GasModel::Plateau;
                changed = true;
            }
            if (gas_model_radio("Gas transport (detailed, slower)",
                                "COUPLED: free gas is carried with the film and redissolves. Slower; deck "
                                "X-runs and onset chart. Steady only.",
                                transport_shown, transport_reason)) {
                case_.gas_model = gui::GasModel::Transport;
                changed = true;
            }
            ImGui::PopID();
            if (!transport_reason.empty()) wrapped_note(transport_reason.c_str());

            const bool fbns = gui::result_form(config_) == gui::ResultForm::Fbns;
            std::string shear = cavity_shear_text(config_.cavitated_film_model);
            if (fbns && config_.cavitated_film_model == CavitatedFilmModel::FULL_FILM && loaded_case_)
                shear += " (kept from the file; new cases give the bubble no shear)";
            wrapped_label("Cavity shear",
                           "Which part of the cavitated film carries the viscous shear. Fixed by the model: "
                           "no shear from the bubble for the pressure plateau, the whole gap for gas transport.");
            wrapped_note(shear.c_str());
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
