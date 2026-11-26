// Boundary conditions panel: five named zones. A click lights the zone in the
// live views (highlighted_zone_); a double-click opens the zone's values.

#include "gui_app.hpp"
#include "imgui.h"

#include <string>

namespace {

constexpr int kZoneCount = 5;
constexpr const char* kZoneNames[kZoneCount] = {"Bearing ends", "Oil feed", "Shaft surface", "Bearing surface",
                                                "Housing outer surface"};
constexpr const char* kZoneHelp[kZoneCount] = {
    "Both axial ends of the film: open to the ambient pressure.",
    "Feed hole or groove: oil supplied above the ambient pressure.",
    "Shaft surface in the film: heat flows into the shaft and from it to the ambient (thermal analysis).",
    "Bore surface in the film: heat conducts into the housing (thermal analysis).",
    "Outside of the housing: heat transfer to the ambient (thermal analysis)."};

bool is_thermal_zone(int zone) { return zone >= 2; }

}  // namespace

/// Zone list with one-line summaries; double-click opens a pop-up with the zone's values.
void GuiApp::draw_boundaries_panel() {
    static int popup_zone = -1;
    bool changed = false;
    auto pressure = [&](double si) {
        return ui::format_value(si, gui::UnitFamily::Pressure, units(), "op.ambient_pressure");
    };
    auto temperature = [&](double si) {
        return ui::format_value(si, gui::UnitFamily::Temperature, units(), "op.ambient_temperature");
    };
    auto heat_transfer = [&](double si) {
        return ui::format_value(si, gui::UnitFamily::HeatTransfer, units(), "bc.htc");
    };

    const std::string supply_temperature = case_.feed.supply_temperature < 0.0
                                               ? "ambient (" + temperature(case_.ambient_temperature) + ")"
                                               : temperature(case_.feed.supply_temperature);
    const std::string summaries[kZoneCount] = {
        "pressure = ambient (" + pressure(case_.ambient_pressure) + "); oil entering at ambient temperature",
        "supply pressure = ambient + " +
            ui::format_value(case_.feed.supply_overpressure, gui::UnitFamily::Pressure, units(),
                             "bc.supply_overpressure") +
            "; supply temperature " + supply_temperature,
        "heat transfer to ambient, " + heat_transfer(case_.housing.shaft_heat_transfer),
        "housing conduction, k = " + ui::format_value(case_.housing.conductivity, gui::UnitFamily::Conductivity,
                                                      units(), "bc.conductivity"),
        "heat transfer to ambient, " + heat_transfer(case_.housing.outer_heat_transfer) + "; outer diameter " +
            (case_.housing.outer_diameter > 0.0
                 ? ui::format_value(case_.housing.outer_diameter, gui::UnitFamily::Length, units(),
                                    "bc.housing_diameter")
                 : std::string("default")),
    };

    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("Click a zone to light it in the live view; double-click to edit its values.");
    ImGui::PopTextWrapPos();
    int open_zone = -1;
    if (ImGui::BeginTable("##zones", 2, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("zone", ImGuiTableColumnFlags_WidthStretch, 0.32f);
        ImGui::TableSetupColumn("summary", ImGuiTableColumnFlags_WidthStretch, 0.68f);
        for (int zone = 0; zone < kZoneCount; ++zone) {
            const bool muted = is_thermal_zone(zone) && !case_.thermal;
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::PushID(zone);
            if (muted) ImGui::PushStyleColor(ImGuiCol_Text, theme::kMuted);
            if (ImGui::Selectable(kZoneNames[zone], highlighted_zone_ == zone,
                                  ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowDoubleClick)) {
                highlighted_zone_ = zone;
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) open_zone = zone;
            }
            if (muted) ImGui::PopStyleColor();
            ui::item_tooltip(kZoneHelp[zone]);
            ImGui::PopID();
            ImGui::TableSetColumnIndex(1);
            ImGui::PushTextWrapPos(0.0f);
            if (muted) ImGui::TextDisabled("thermal analysis off");
            else ImGui::TextUnformatted(summaries[zone].c_str());
            ImGui::PopTextWrapPos();
        }
        ImGui::EndTable();
    }
    if (open_zone >= 0) {
        popup_zone = open_zone;
        ImGui::OpenPopup("Zone values");
    }
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 30.0f, 0.0f));
    if (ImGui::BeginPopup("Zone values")) {
        const int zone = popup_zone;
        if (zone >= 0 && zone < kZoneCount) {
            ImGui::TextUnformatted(kZoneNames[zone]);
            ImGui::Separator();
            if (is_thermal_zone(zone) && !case_.thermal)
                ImGui::TextColored(theme::kRiskMedium, "Thermal analysis is off; these values are not used.");
            if (ui::begin_form("##zone_form")) {
                switch (zone) {
                    case 0:
                        ui::form_value("Pressure", "Ambient pressure (Operating conditions).",
                                       pressure(case_.ambient_pressure));
                        ui::form_value("Inflow temperature", "Ambient temperature (Operating conditions).",
                                       temperature(case_.ambient_temperature));
                        break;
                    case 1: {
                        changed |= ui::form_number(
                            "Supply overpressure", "Feed pressure above ambient. Typical: 0 to 1 bar.",
                            "bc.supply_overpressure", case_.feed.supply_overpressure, gui::UnitFamily::Pressure,
                            units());
                        bool same_as_ambient = case_.feed.supply_temperature < 0.0;
                        if (ui::form_checkbox("Supply temperature same as ambient",
                                              "Oil is supplied at the ambient temperature.", "bc.supply_t_amb",
                                              same_as_ambient)) {
                            case_.feed.supply_temperature = same_as_ambient ? -1.0 : case_.ambient_temperature;
                            changed = true;
                        }
                        if (!same_as_ambient)
                            changed |= ui::form_number("Supply temperature", "Temperature of the supplied oil.",
                                                       "bc.supply_temperature", case_.feed.supply_temperature,
                                                       gui::UnitFamily::Temperature, units());
                        if (case_.analysis == gui::AnalysisType::Transient)
                            changed |= ui::form_checkbox(
                                "Feed hole turns with the shaft",
                                "The feed is in the shaft and rotates with it (transient only).",
                                "bc.feed_turns", case_.feed.turns_with_shaft);
                        break;
                    }
                    case 2:
                        changed |= ui::form_number(
                            "Heat transfer coefficient to ambient",
                            "From the exposed shaft to the ambient. Typical: 200 to 1000 W/(m2 K).", "bc.htc",
                            case_.housing.shaft_heat_transfer, gui::UnitFamily::HeatTransfer, units());
                        changed |= ui::form_number(
                            "Exposed area ratio", "Exposed shaft area per film area. 1 = equal areas.",
                            "bc.shaft_area_ratio", case_.housing.shaft_area_ratio, gui::UnitFamily::Fixed,
                            units(), "%.4g");
                        break;
                    case 3:
                        changed |= ui::form_number(
                            "Housing conductivity", "Thermal conductivity of the bearing housing. Steel: 15 to 50.",
                            "bc.conductivity", case_.housing.conductivity, gui::UnitFamily::Conductivity, units());
                        break;
                    case 4: {
                        bool use_default = case_.housing.outer_diameter <= 0.0;
                        if (ui::form_checkbox("Default outer diameter", "Housing outer diameter = 1.1 x D.",
                                              "bc.housing_default", use_default)) {
                            case_.housing.outer_diameter = use_default ? 0.0 : 1.1 * case_.bearing_diameter;
                            changed = true;
                        }
                        if (!use_default)
                            changed |= ui::form_number("Housing outer diameter", "Outer diameter of the housing.",
                                                       "bc.housing_diameter", case_.housing.outer_diameter,
                                                       gui::UnitFamily::Length, units());
                        changed |= ui::form_number(
                            "Heat transfer coefficient", "From the housing outside to the ambient. Typical: 20 to 100.",
                            "bc.htc", case_.housing.outer_heat_transfer, gui::UnitFamily::HeatTransfer, units());
                        break;
                    }
                    default: break;
                }
                ui::end_form();
            }
            if (zone == 0) {
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextDisabled("These follow the operating conditions; change them there.");
                ImGui::PopTextWrapPos();
            }
        }
        ImGui::EndPopup();
    }
    if (changed) mark_changed();
}
