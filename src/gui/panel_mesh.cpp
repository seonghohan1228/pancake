// Mesh panel: mesh level with cell counts, resulting cell sizes, and
// (Advanced) custom cell counts and cross-film layers.

#include "gui_app.hpp"
#include "imgui.h"

#include <cmath>
#include <cstdio>
#include <string>

/// Mesh level, cell sizes around and along the film; Advanced: custom counts, cross-film layers.
void GuiApp::draw_mesh_panel() {
    constexpr double kPi = 3.14159265358979323846;
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Mesh)];
    if (!ui::begin_form("##mesh")) return;
    bool changed = false;

    ui::form_label("Mesh level",
                   "Number of cells around x along the film. Finer meshes are more accurate and slower; "
                   "Medium is the usual choice.");
    ImGui::PushID("mesh.level");
    for (gui::MeshLevel level : {gui::MeshLevel::Coarse, gui::MeshLevel::Medium, gui::MeshLevel::Fine}) {
        int n_theta = 0;
        int n_z = 0;
        gui::mesh_level_cells(level, n_theta, n_z);
        char label[64];
        const char* name = level == gui::MeshLevel::Coarse   ? "Coarse"
                           : level == gui::MeshLevel::Medium ? "Medium"
                                                             : "Fine";
        std::snprintf(label, sizeof(label), "%-7s %d x %d", name, n_theta, n_z);
        if (ImGui::RadioButton(label, case_.mesh_level == level)) {
            case_.mesh_level = level;
            changed = true;
        }
        ui::item_tooltip(level == gui::MeshLevel::Coarse   ? "Quick check of a trend."
                         : level == gui::MeshLevel::Medium ? "Design runs."
                                                           : "Check that the result no longer changes with the mesh.");
    }
    if (case_.mesh_level == gui::MeshLevel::Custom)
        ImGui::TextColored(theme::kRiskMedium, "Custom %d x %d (see Advanced)", case_.n_theta, case_.n_z);
    ImGui::PopID();

    int n_theta = case_.n_theta;
    int n_z = case_.n_z;
    if (case_.mesh_level != gui::MeshLevel::Custom) gui::mesh_level_cells(case_.mesh_level, n_theta, n_z);
    if (n_theta > 0 && n_z > 0) {
        const double arc = kPi * case_.bearing_diameter / n_theta;
        const double along = case_.bearing_length / n_z;
        char text[128];
        std::snprintf(text, sizeof(text), "%.3g deg, ", 360.0 / n_theta);
        ui::form_value("Cell size around", "Arc length per cell around the film (angle and length).",
                       text + ui::format_value(arc, gui::UnitFamily::SmallLength, units(), "mesh.cell_size"));
        ui::form_value("Cell size along", "Bearing length per cell, L / n_z.",
                       ui::format_value(along, gui::UnitFamily::SmallLength, units(), "mesh.cell_size"));
        if (case_.feed.type == gui::FeedType::Hole && arc > 0.0 && along > 0.0) {
            std::snprintf(text, sizeof(text), "about %.0f around, %.0f along",
                          case_.feed.hole_diameter / arc, case_.feed.hole_diameter / along);
            ui::form_value("Cells across the feed hole",
                           "How many cells the feed hole spans. At least 4 resolves the hole.", text);
        } else if (case_.feed.type == gui::FeedType::Groove) {
            std::snprintf(text, sizeof(text), "about %.0f around", case_.feed.groove_width_deg / (360.0 / n_theta));
            ui::form_value("Cells across the groove", "How many cells the feed groove spans.", text);
        }
    }

    if (advanced) {
        bool custom = case_.mesh_level == gui::MeshLevel::Custom;
        if (ui::form_checkbox("Custom cell counts", "Enter the cell counts directly.", "mesh.custom", custom)) {
            if (custom) {
                case_.n_theta = n_theta;
                case_.n_z = n_z;
                case_.mesh_level = gui::MeshLevel::Custom;
            } else {
                case_.mesh_level = gui::MeshLevel::Medium;
            }
            changed = true;
        }
        ImGui::BeginDisabled(!custom);
        if (ui::form_int("Cells around n_theta", "Cells around the circumference. Typical: 180 to 1440.",
                         "mesh.n_theta", case_.n_theta, 8, 100000)) {
            case_.mesh_level = gui::MeshLevel::Custom;
            changed = true;
        }
        if (ui::form_int("Cells along n_z", "Cells along the bearing length. Typical: 48 to 384.", "mesh.n_z",
                         case_.n_z, 4, 100000)) {
            case_.mesh_level = gui::MeshLevel::Custom;
            changed = true;
        }
        ImGui::EndDisabled();
        ImGui::BeginDisabled(!case_.thermal);
        changed |= ui::form_int("Cross-film layers",
                                "Layers across the film thickness for the temperature (thermal analysis only). "
                                "Typical: 16.",
                                "mesh.cross_film", case_.cross_film_layers, 2, 256);
        ImGui::EndDisabled();
        if (!case_.thermal) {
            ui::form_label("");
            ImGui::TextDisabled("Used only with thermal analysis.");
        }
    }
    ui::end_form();
    if (changed) mark_changed();
}
