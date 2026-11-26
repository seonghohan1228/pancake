#include "ui_widgets.hpp"

#include "theme.hpp"

#include <algorithm>
#include <cstdio>

namespace ui {

ImVec4 light_color(gui::Light light) {
    switch (light) {
        case gui::Light::Green: return theme::kRiskLow;
        case gui::Light::Yellow: return theme::kRiskMedium;
        case gui::Light::Red: return theme::kRiskHigh;
        case gui::Light::Idle: break;
    }
    return theme::kMuted;
}

ImVec4 risk_color(gui::Risk risk) {
    switch (risk) {
        case gui::Risk::Low: return theme::kRiskLow;
        case gui::Risk::Medium: return theme::kRiskMedium;
        case gui::Risk::High: return theme::kRiskHigh;
        case gui::Risk::None: break;
    }
    return theme::kMuted;
}

void status_dot(gui::Light light) {
    const float size = ImGui::GetTextLineHeight();
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 centre(pos.x + size * 0.5f, pos.y + size * 0.55f);
    const ImU32 colour = ImGui::GetColorU32(light_color(light));
    if (light == gui::Light::Idle) draw->AddCircle(centre, size * 0.3f, colour, 0, 1.5f);
    else draw->AddCircleFilled(centre, size * 0.32f, colour);
    ImGui::Dummy(ImVec2(size, size));
}

void item_tooltip(const char* text) {
    if (!text || !*text) return;
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip | ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void help_marker(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    item_tooltip(text);
}

void section(const char* title) {
    ImGui::Spacing();
    ImGui::SeparatorText(title);
}

bool advanced_toggle(const char* id, bool& open) {
    const char* label = "Advanced";
    const float width = ImGui::CalcTextSize(label).x + ImGui::GetFrameHeight() +
                        ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SameLine();
    const float right = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;
    ImGui::SameLine(std::max(ImGui::GetCursorPosX(), right - width));
    ImGui::PushID(id);
    const bool changed = ImGui::Checkbox(label, &open);
    ImGui::PopID();
    item_tooltip("Show tolerances, exact mesh sizes, tilt and the settings chosen automatically.");
    return changed;
}

bool begin_form(const char* id) {
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingStretchProp)) return false;
    ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.45f);
    ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthStretch, 0.55f);
    return true;
}

void end_form() { ImGui::EndTable(); }

void form_label(const char* label, const char* help) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    item_tooltip(help);
    ImGui::TableSetColumnIndex(1);
}

namespace {
int& choice_for(UnitChoices& choices, const char* id, gui::UnitFamily family) {
    int& choice = choices[id];
    const int count = static_cast<int>(gui::unit_options(family).size());
    if (choice < 0 || choice >= count) choice = 0;
    return choice;
}
}  // namespace

bool unit_input(const char* id, double& si, gui::UnitFamily family, UnitChoices& choices,
                const char* format) {
    const auto& options = gui::unit_options(family);
    int& choice = choice_for(choices, id, family);
    const gui::UnitOption& unit = options[static_cast<std::size_t>(choice)];
    ImGui::PushID(id);
    const float unit_width = options.size() > 1 ? ImGui::GetFontSize() * 5.2f
                                                : ImGui::CalcTextSize(unit.label).x + 8.0f;
    ImGui::SetNextItemWidth(-unit_width);
    double display = gui::to_display(unit, si);
    bool changed = false;
    if (ImGui::InputDouble("##v", &display, 0.0, 0.0, format)) {
        si = gui::to_si(unit, display);
        changed = true;
    }
    ImGui::SameLine();
    if (options.size() == 1) {
        ImGui::TextDisabled("%s", unit.label);
    } else {
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##u", unit.label)) {
            for (int i = 0; i < static_cast<int>(options.size()); ++i)
                if (ImGui::Selectable(options[static_cast<std::size_t>(i)].label, i == choice)) choice = i;
            ImGui::EndCombo();
        }
    }
    ImGui::PopID();
    return changed;
}

bool form_number(const char* label, const char* help, const char* id, double& si,
                 gui::UnitFamily family, UnitChoices& choices, const char* format) {
    form_label(label, help);
    return unit_input(id, si, family, choices, format);
}

bool form_int(const char* label, const char* help, const char* id, int& value, int min_value,
              int max_value) {
    form_label(label, help);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    int edited = value;
    const bool changed = ImGui::InputInt("##i", &edited, 0, 0);
    ImGui::PopID();
    if (changed) value = std::clamp(edited, min_value, max_value);
    return changed;
}

bool form_checkbox(const char* label, const char* help, const char* id, bool& value) {
    form_label(label, help);
    ImGui::PushID(id);
    const bool changed = ImGui::Checkbox("##c", &value);
    ImGui::PopID();
    return changed;
}

void form_value(const char* label, const char* help, const std::string& text) {
    form_label(label, help);
    ImGui::TextDisabled("%s", text.c_str());
}

std::string format_value(double si, gui::UnitFamily family, UnitChoices& choices, const char* id,
                         const char* format) {
    const auto& options = gui::unit_options(family);
    const gui::UnitOption& unit = options[static_cast<std::size_t>(choice_for(choices, id, family))];
    char number[48];
    std::snprintf(number, sizeof(number), format, gui::to_display(unit, si));
    std::string text = number;
    if (unit.label[0]) text += std::string(" ") + unit.label;
    return text;
}

}  // namespace ui
