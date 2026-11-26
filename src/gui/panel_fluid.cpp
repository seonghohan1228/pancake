// Fluid panel: oil-refrigerant pair, refrigerant content, surface roughness.

#include "gui_app.hpp"
#include "imgui.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>

namespace {

std::string lower(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return text;
}

/// True when the pair header says it is example data only.
bool is_example_only(const gui::FluidPair& pair) {
    const std::string note = lower(pair.note + " " + pair.source);
    return note.find("synthetic") != std::string::npos || note.find("not for design") != std::string::npos;
}

/// Read-only description of the pair's gas-phase block.
const char* gas_phase_text(const gui::FluidPair& pair) {
    if (!pair.has_gas_phase) return "none: pressure plateau only (no gas transport)";
    switch (pair.solver_species) {
        case DissolvedGasSpecies::R410A: return "R410A real-fluid model";
        case DissolvedGasSpecies::PROPANE: return "Propane, ideal gas";
        case DissolvedGasSpecies::AIR: return "Air, ideal gas";
    }
    return "-";
}

}  // namespace

/// Pair choice, refrigerant content, roughness; Advanced: fluid properties at ambient.
void GuiApp::draw_fluid_panel() {
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Fluid)];
    if (!ui::begin_form("##fluid")) return;
    bool changed = false;

    // Pair
    ui::form_label("Oil-refrigerant pair",
                   "Oil and dissolved refrigerant with their property tables (fluid directory).");
    std::string preview = case_.fluid_pair.empty() ? std::string("(choose a pair)") : case_.fluid_pair;
    if (pair_ok_ && !pair_.name.empty()) preview = pair_.name;
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##pair", preview.c_str())) {
        for (const std::string& id : pair_ids_) {
            const bool selected = id == case_.fluid_pair;
            const std::string label = (selected && pair_ok_ && !pair_.name.empty()) ? pair_.name + "##" + id : id;
            if (ImGui::Selectable(label.c_str(), selected) && !selected) {
                select_pair(id);
                case_.fluid_pair = id;
                changed = true;
            }
        }
        if (pair_ids_.empty()) ImGui::TextDisabled("No pairs in the fluid directory.");
        ImGui::EndCombo();
    }
    ui::form_label("");
    if (ImGui::Button("Open fluid library...")) show_fluid_library_ = true;
    ui::item_tooltip("View, edit, copy or add oil-refrigerant pairs and their tables.");

    if (!pair_error_.empty()) {
        ui::form_label("");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::kRiskHigh, "Cannot load the pair: %s", pair_error_.c_str());
        ImGui::PopTextWrapPos();
    } else if (pair_ok_) {
        if (!pair_.note.empty()) {
            ui::form_label("Note");
            ImGui::PushTextWrapPos(0.0f);
            if (is_example_only(pair_)) ImGui::TextColored(theme::kRiskMedium, "%s", pair_.note.c_str());
            else ImGui::TextUnformatted(pair_.note.c_str());
            ImGui::PopTextWrapPos();
        }
        if (!pair_.source.empty()) {
            ui::form_label("Source", "Where the property data come from.");
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("%s", pair_.source.c_str());
            ImGui::PopTextWrapPos();
        }
    }

    if (pair_ok_) {
        ui::form_label("Gas phase",
                       "Data of the released refrigerant gas, needed for gas transport. Edit it in the fluid "
                       "library.");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("%s", gas_phase_text(pair_));
        ImGui::PopTextWrapPos();
    }

    // Refrigerant content
    ui::form_label("Refrigerant content", "Mass fraction of refrigerant dissolved in the oil.");
    bool saturated_edit = false;
    {
        // Radio buttons that wrap in a narrow panel instead of running past its edge.
        const char* given_label = "Given";
        ImGui::PushID("fluid.saturated");
        if (ImGui::RadioButton("Saturated at ambient pressure", case_.saturated) && !case_.saturated) {
            case_.saturated = true;
            saturated_edit = true;
        }
        ui::item_tooltip("Content at which the oil is in equilibrium with the refrigerant at the ambient pressure "
                         "and temperature.");
        ImGui::SameLine();
        const float given_w = ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x +
                              ImGui::CalcTextSize(given_label).x;
        if (ImGui::GetContentRegionAvail().x < given_w) ImGui::NewLine();
        if (ImGui::RadioButton(given_label, !case_.saturated) && case_.saturated) {
            case_.saturated = false;
            saturated_edit = true;
        }
        ui::item_tooltip("Enter the content directly.");
        ImGui::PopID();
    }
    if (saturated_edit) {
        if (!case_.saturated && case_.content <= 0.0) case_.content = estimates_.content;
        changed = true;
    }
    if (!case_.saturated) {
        changed |= ui::form_number("Given content", "Refrigerant mass fraction in the oil. Typical: 5 to 30 wt%.",
                                   "fluid.content", case_.content, gui::UnitFamily::Content, units());
    }
    if (pair_ok_) {
        const std::string used = ui::format_value(estimates_.content, gui::UnitFamily::Content, units(),
                                                  "fluid.content") +
                                 " at " +
                                 ui::format_value(case_.ambient_pressure, gui::UnitFamily::Pressure, units(),
                                                  "op.ambient_pressure") +
                                 ", " +
                                 ui::format_value(case_.ambient_temperature, gui::UnitFamily::Temperature,
                                                  units(), "op.ambient_temperature");
        ui::form_value("Content used", "Content the run uses, at the ambient pressure and temperature.", used);
    }

    // Roughness
    const char* roughness_help =
        "RMS surface roughness Rq. Used only for the film thickness ratio Lambda = h_min / sigma, "
        "sigma = sqrt(Rq_shaft^2 + Rq_bush^2); there is no contact model. Typical: 0.1 to 0.4 um.";
    changed |= ui::form_number("Surface roughness Rq, shaft", roughness_help, "fluid.rq_shaft",
                               case_.roughness_shaft, gui::UnitFamily::SmallLength, units());
    changed |= ui::form_number("Surface roughness Rq, bush", roughness_help, "fluid.rq_bush",
                               case_.roughness_bush, gui::UnitFamily::SmallLength, units());

    // Range
    if (pair_ok_ && !estimates_.fluid_in_range) {
        ui::form_label("Table range");
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::kRiskMedium,
                           "The ambient state is outside the fluid table; properties are clamped to the "
                           "table ends.");
        ImGui::PopTextWrapPos();
        ui::form_label("");
        if (ImGui::Checkbox("Run outside the table range (results marked)", &allow_out_of_range_))
            refresh_derived();
        ui::item_tooltip("Accept the clamped properties. The results are marked as outside the table range.");
    } else if (pair_ok_) {
        ui::form_label("Table range");
        ImGui::TextColored(theme::kRiskLow, "Within the table range.");
    }

    // Advanced: properties at ambient
    if (advanced && pair_ok_) {
        const double t_amb = case_.ambient_temperature;
        ui::form_value("Viscosity at ambient",
                       "Dynamic viscosity at the ambient temperature and the content used.",
                       ui::format_value(estimates_.viscosity, gui::UnitFamily::Viscosity, units(),
                                        "fluid.viscosity"));
        char text[64];
        std::snprintf(text, sizeof(text), "%.1f kg/m3", pair_.density(t_amb));
        ui::form_value("Density at ambient", "Liquid density, linear law of the pair.", text);
        ui::form_value("Saturation pressure",
                       "Pressure below which dissolved refrigerant comes out of solution, at the content "
                       "used and the ambient temperature.",
                       ui::format_value(estimates_.saturation_pressure, gui::UnitFamily::Pressure, units(),
                                        "est.p_sat"));
        ui::form_value("Oil vapour pressure",
                       "Pressure at which the oil itself boils (vapour cavitation, VC).",
                       ui::format_value(pair_.oil_vapour_pressure, gui::UnitFamily::Pressure, units(),
                                        "fluid.p_vap"));
    }
    ui::end_form();
    if (changed) mark_changed();
}
