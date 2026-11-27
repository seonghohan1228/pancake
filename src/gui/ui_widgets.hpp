#pragma once

// Shared ImGui building blocks for the setup panels and result views: a
// two-column form (label | widget), unit-aware number inputs, status dots,
// risk colours and tooltips. Every panel uses these so labels, units and
// tooltips look the same everywhere.

#include "case_status.hpp"
#include "imgui.h"
#include "key_results.hpp"
#include "units.hpp"

#include <initializer_list>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ui {

/// Display-unit choice per field id (persisted in UserSettings::unit_choices).
using UnitChoices = std::map<std::string, int>;

ImVec4 light_color(gui::Light light);
ImVec4 risk_color(gui::Risk risk);

/// Filled circle inline with text (status light).
void status_dot(gui::Light light);
/// Tooltip on the last item, wrapped.
void item_tooltip(const char* text);
/// Small "(?)" marker with a tooltip.
void help_marker(const char* text);
/// Bold section title with a separator.
void section(const char* title);
/// Right-aligned "Advanced" checkbox on the current line; returns true when toggled.
bool advanced_toggle(const char* id, bool& open);

/// Two-column form. Rows add a label (with tooltip) and leave the cursor in the
/// widget column. Must be closed with end_form() when begin_form() returned true.
bool begin_form(const char* id);
void end_form();
void form_label(const char* label, const char* help = nullptr);

/// Number in display units with a unit selector; edits `si`. Returns true on edit.
bool unit_input(const char* id, double& si, gui::UnitFamily family, UnitChoices& choices,
                const char* format = "%.6g");
/// Form row: label + unit_input.
bool form_number(const char* label, const char* help, const char* id, double& si,
                 gui::UnitFamily family, UnitChoices& choices, const char* format = "%.6g");
bool form_int(const char* label, const char* help, const char* id, int& value, int min_value,
              int max_value);
bool form_checkbox(const char* label, const char* help, const char* id, bool& value);
/// Read-only value row (Expanded "Auto" items).
void form_value(const char* label, const char* help, const std::string& text);
/// Value formatted in the user's display unit, e.g. "24.5 bar".
std::string format_value(double si, gui::UnitFamily family, UnitChoices& choices, const char* id,
                         const char* format = "%.4g");

/// True when at least one (x, y) pair is finite.
bool has_points(const std::vector<double>& x, const std::vector<double>& y);
/// Axis flags that auto-fit only when the plot has data (no tick labels when empty). ImPlot's fit on an
/// empty plot changes the limits every frame, so empty plots must not fit.
int fit_if(bool has_data);
/// For an empty plot, right after its axes are set up: fixed limits and a
/// centred "No data yet". Ends the setup phase.
void empty_plot(const char* text = "No data yet");
/// Soft drop shadow around the current popup or floating window, so it stands
/// out from the windows behind it. Call right after a successful Begin*/
/// BeginPopup*. Dear ImGui has no window shadows; this draws fading rings just
/// outside the window on its own draw list (nothing for docked windows).
void window_shadow();
/// Speed in all three units, e.g. "3000 rpm = 50 Hz = 314.2 rad/s".
std::string speed_all_units(double rad_s);

/// Combo over a small enum set. Options: {value, label, tooltip}.
template <typename Enum>
struct Option {
    Enum value;
    const char* label;
    const char* tooltip = nullptr;
};

template <typename Enum>
bool combo(const char* id, Enum& value, std::initializer_list<Option<Enum>> options) {
    const char* preview = "?";
    for (const auto& option : options)
        if (option.value == value) preview = option.label;
    bool changed = false;
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo(id, preview)) {
        for (const auto& option : options) {
            if (ImGui::Selectable(option.label, option.value == value)) {
                value = option.value;
                changed = true;
            }
            if (option.tooltip) item_tooltip(option.tooltip);
        }
        ImGui::EndCombo();
    }
    return changed;
}

/// Radio buttons on one line over a small enum set.
template <typename Enum>
bool radio_row(const char* id, Enum& value, std::initializer_list<Option<Enum>> options) {
    bool changed = false;
    ImGui::PushID(id);
    bool first = true;
    for (const auto& option : options) {
        if (!first) ImGui::SameLine();
        first = false;
        if (ImGui::RadioButton(option.label, value == option.value)) {
            value = option.value;
            changed = true;
        }
        if (option.tooltip) item_tooltip(option.tooltip);
    }
    ImGui::PopID();
    return changed;
}

template <typename Enum>
bool form_combo(const char* label, const char* help, const char* id, Enum& value,
                std::initializer_list<Option<Enum>> options) {
    form_label(label, help);
    return combo(id, value, options);
}

template <typename Enum>
bool form_radio(const char* label, const char* help, const char* id, Enum& value,
                std::initializer_list<Option<Enum>> options) {
    form_label(label, help);
    return radio_row(id, value, options);
}

}  // namespace ui
