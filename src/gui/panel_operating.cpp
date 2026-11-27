// Operating conditions panel: analysis type, speed, ambient state, load
// (steady value or transient load table) and the shaft position mode.

#include "gui_app.hpp"
#include "exports.hpp"
#include "imgui.h"
#include "implot.h"
#include "win_util.hpp"

#include <cmath>
#include <cstring>
#include <initializer_list>
#include <string>
#include <vector>

namespace {

/// Form row with a plain double and a fixed unit text (no unit selector).
bool plain_number(const char* label, const char* help, const char* id, double& value, const char* unit,
                  const char* format = "%.6g") {
    ui::form_label(label, help);
    ImGui::PushID(id);
    const float unit_width = unit[0] ? ImGui::CalcTextSize(unit).x + 8.0f : 0.0f;
    ImGui::SetNextItemWidth(unit[0] ? -unit_width : -FLT_MIN);
    const bool changed = ImGui::InputDouble("##v", &value, 0.0, 0.0, format);
    if (unit[0]) {
        ImGui::SameLine();
        ImGui::TextDisabled("%s", unit);
    }
    ImGui::PopID();
    return changed;
}

/// Form row of radio buttons that continue on the next line when the column is
/// too narrow (no label is cut off at the panel edge). Returns true on edit.
template <typename Enum>
bool radio_row(const char* label, const char* help, const char* id, Enum& value,
               std::initializer_list<ui::Option<Enum>> options) {
    ui::form_label(label, help);
    ImGui::PushID(id);
    bool changed = false, first = true;
    for (const auto& option : options) {
        if (!first) {
            const float width = ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x +
                                ImGui::CalcTextSize(option.label).x;
            ImGui::SameLine();
            if (ImGui::GetContentRegionAvail().x < width) ImGui::NewLine();
        }
        first = false;
        if (ImGui::RadioButton(option.label, value == option.value) && value != option.value) {
            value = option.value;
            changed = true;
        }
        if (option.tooltip) ui::item_tooltip(option.tooltip);
    }
    ImGui::PopID();
    return changed;
}

}  // namespace

/// Analysis type, speed, ambient state, load; transient load table and plot.
void GuiApp::draw_operating_panel() {
    const bool advanced = advanced_[static_cast<std::size_t>(gui::Node::Operating)];
    if (!ui::begin_form("##operating")) return;
    bool changed = false;

    if (radio_row<gui::AnalysisType>(
            "Analysis type", "Steady: fixed operating point. Transient: load varying in time.", "op.analysis",
            case_.analysis,
            {{gui::AnalysisType::Steady, "Steady", "Fixed operating point; the shaft finds its load balance."},
             {gui::AnalysisType::Transient, "Transient",
              "Load from a table over the cycle; the shaft moves with the film force."}})) {
        reload_load_table();
        changed = true;
    }
    changed |= ui::form_number("Speed", "Shaft rotational speed. Typical compressor: 3000 to 9000 rpm.",
                               "op.speed", case_.speed, gui::UnitFamily::Speed, units());
    changed |= ui::form_number("Ambient pressure",
                               "Pressure around the bearing (chamber or shell), absolute. Typical: 5 to 30 bar.",
                               "op.ambient_pressure", case_.ambient_pressure, gui::UnitFamily::Pressure,
                               units());
    changed |= ui::form_number("Ambient temperature",
                               "Temperature of the surroundings and of the supplied oil. Typical: 40 to 90 degC.",
                               "op.ambient_temperature", case_.ambient_temperature,
                               gui::UnitFamily::Temperature, units());

    if (case_.analysis == gui::AnalysisType::Steady) {
        changed |= ui::form_number("Load", "Force on the bearing, magnitude W.", "op.load", case_.load,
                                   gui::UnitFamily::Force, units());
        changed |= ui::form_number("Load direction",
                                   "Direction of the load acting on the shaft, fixed frame, in the direction "
                                   "of rotation from the x axis; 270 deg = downward.",
                                   "op.load_direction", case_.load_direction_deg, gui::UnitFamily::Angle,
                                   units());
        if (advanced) {
            changed |= radio_row<gui::ShaftPosition>(
                "Shaft position", "How the shaft position is found.", "op.position", case_.position,
                {{gui::ShaftPosition::LoadBalance, "From load balance",
                  "The shaft moves until the film force carries the load."},
                 {gui::ShaftPosition::Fixed, "Fixed (no load balance)",
                  "The shaft stays at the given eccentricity and attitude; the film force is a result."}});
            if (case_.position == gui::ShaftPosition::Fixed) {
                changed |= ui::form_number("Eccentricity ratio", "Shaft centre offset per clearance, e / c (0 to 1).",
                                           "op.fixed_eps", case_.fixed_eccentricity_ratio, gui::UnitFamily::Fixed,
                                           units(), "%.4g");
                changed |= ui::form_number("Attitude angle",
                                           "Angle from the load direction to the line of centres (h_min).",
                                           "op.fixed_attitude", case_.fixed_attitude_deg, gui::UnitFamily::Angle,
                                           units());
            }
        }
        ui::end_form();
        if (changed) mark_changed();
        return;
    }

    // Transient: load table
    static char path_buffer[512] = {};
    static std::string synced_path;
    static bool path_active = false;
    const std::string current = winutil::path_to_utf8(case_.load_table.file);
    ui::form_label("Load table file",
                   "CSV with abscissa (crank angle or time) and the load (magnitude and direction, or x and y "
                   "components) per row.");
    const float browse_width = ImGui::CalcTextSize("Browse...").x + ImGui::GetStyle().FramePadding.x * 2.0f +
                               ImGui::GetStyle().ItemSpacing.x;
    ImGui::SetNextItemWidth(-browse_width);
    if (!path_active && current != synced_path) {
        std::strncpy(path_buffer, current.c_str(), sizeof(path_buffer) - 1);
        path_buffer[sizeof(path_buffer) - 1] = '\0';
        synced_path = current;
    }
    ImGui::InputText("##load_file", path_buffer, sizeof(path_buffer));
    path_active = ImGui::IsItemActive();
    if (ImGui::IsItemDeactivatedAfterEdit()) {
        case_.load_table.file = winutil::utf8_to_path(path_buffer);
        synced_path = path_buffer;
        reload_load_table();
        changed = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Browse...")) {
        const std::filesystem::path file =
            exports::open_file_dialog(hwnd, L"CSV files (*.csv)\0*.csv\0All files\0*.*\0");
        if (!file.empty()) {
            case_.load_table.file = file;
            reload_load_table();
            changed = true;
        }
    }
    changed |= ui::form_combo<gui::LoadAbscissa>(
        "Abscissa", "What the first column of the table is.", "op.load_abscissa", case_.load_table.abscissa,
        {{gui::LoadAbscissa::CrankAngle, "Crank angle", "First column in degrees of crank angle."},
         {gui::LoadAbscissa::Time, "Time", "First column in seconds."}});
    changed |= ui::form_combo<gui::LoadFormat>(
        "Load format", "How the load columns are given.", "op.load_format", case_.load_table.format,
        {{gui::LoadFormat::MagnitudeDirection, "Magnitude and direction",
          "Columns: magnitude [N] and direction [deg] of the load on the shaft, fixed frame."},
         {gui::LoadFormat::Components, "x and y components", "Columns: Fx and Fy [N], fixed frame."}});
    ImGui::BeginDisabled(case_.load_table.abscissa != gui::LoadAbscissa::CrankAngle);
    changed |= plain_number("Cycle frequency",
                            "Load cycles per second for a crank-angle table (usually the shaft speed in Hz).",
                            "op.load_frequency", case_.load_table.frequency_hz, "Hz");
    ImGui::EndDisabled();
    changed |= ui::form_number("Phase offset", "Shift of the table abscissa (crank angle at the run start).",
                               "op.load_phase", case_.load_table.phase_offset_deg, gui::UnitFamily::Angle,
                               units());
    changed |= plain_number("Scale", "Factor applied to every load value of the table. 1 = as given.",
                            "op.load_scale", case_.load_table.scale, "");

    ui::form_label("Load over the cycle");
    if (load_table_ok_ && !load_rows_.empty()) {
        std::vector<double> x;
        std::vector<double> y;
        x.reserve(load_rows_.size());
        y.reserve(load_rows_.size());
        for (const LoadRow& row : load_rows_) {
            x.push_back(row.abscissa);
            const double magnitude = case_.load_table.format == gui::LoadFormat::MagnitudeDirection
                                         ? row.a
                                         : std::hypot(row.a, row.b);
            y.push_back(magnitude * case_.load_table.scale);
        }
        if (ImPlot::BeginPlot("##load_plot", ImVec2(-1, ImGui::GetFontSize() * 11.0f))) {
            ImPlot::SetupAxes(case_.load_table.abscissa == gui::LoadAbscissa::CrankAngle ? "Crank angle [deg]"
                                                                                       : "Time [s]",
                              "Load [N]", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::PlotLine("load magnitude", x.data(), y.data(), static_cast<int>(x.size()));
            ImPlot::EndPlot();
        }
    } else if (case_.load_table.file.empty()) {
        ImGui::TextColored(theme::kRiskMedium, "Choose a load table file.");
    } else {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextColored(theme::kRiskHigh, "Load table: %s", load_table_error_.c_str());
        ImGui::PopTextWrapPos();
    }

    changed |= ui::form_number("Moving mass", "Mass of the part that moves with the bearing film.",
                               "op.moving_mass", case_.moving_mass, gui::UnitFamily::Mass, units());
    changed |= plain_number("Revolutions", "Shaft revolutions to simulate. Typical: 3 (the last one is periodic).",
                            "op.revolutions", case_.revolutions, "rev", "%.4g");
    ui::end_form();
    if (changed) mark_changed();
}
