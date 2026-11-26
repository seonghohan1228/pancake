// Fluid library dialog: list the oil-refrigerant pairs of the fluid directory,
// edit a pair's header and its table (in the file units of solution.csv), save
// it through the solver's loader, and plot what the solver will read.

#include "exports.hpp"
#include "gui_app.hpp"
#include "imgui.h"
#include "implot.h"
#include "win_util.hpp"

#include <shellapi.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

/// Editor state of the dialog (purely visual; the case only sees saved pairs).
struct LibraryEditor {
    std::string editing_id;               // pair shown on the right; empty = none
    bool is_new = false;                  // not yet written to the fluid directory
    bool modified = false;
    gui::FluidPair pair;                  // header being edited
    std::vector<gui::SolutionRow> rows;   // table being edited (file units)
    gui::FluidPair plotted;               // last version accepted by the solver's loader
    bool plotted_ok = false;
    std::string message;
    bool message_error = false;
    char new_name[64] = "my_pair";
    char copy_name[64] = "";
    std::vector<std::string> listed;      // pair ids the synthetic flags were read for
    std::map<std::string, bool> synthetic;
    // Cached saturation curve of `plotted` at the case ambient temperature.
    std::string curve_key;
    std::vector<double> curve_content, curve_pressure;
};

LibraryEditor& editor() {
    static LibraryEditor state;
    return state;
}

int string_resize(ImGuiInputTextCallbackData* data) {
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* text = static_cast<std::string*>(data->UserData);
        text->resize(static_cast<std::size_t>(data->BufTextLen));
        data->Buf = text->data();
    }
    return 0;
}

/// InputText bound to a std::string.
bool input_string(const char* id, std::string& text) {
    return ImGui::InputText(id, text.data(), text.capacity() + 1, ImGuiInputTextFlags_CallbackResize,
                            string_resize, &text);
}

bool form_text(const char* label, const char* help, const char* id, std::string& text) {
    ui::form_label(label, help);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    const bool changed = input_string("##t", text);
    ImGui::PopID();
    return changed;
}

bool is_synthetic_note(const std::string& note) {
    std::string lower = note;
    for (char& ch : lower) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return lower.find("synthetic") != std::string::npos || lower.find("not for design") != std::string::npos;
}

/// Reads only the note line of <root>/<id>/pair.txt (cheap, no table parse).
bool pair_is_synthetic(const fs::path& root, const std::string& id) {
    std::ifstream input(root / id / "pair.txt");
    std::string line;
    while (std::getline(input, line)) {
        const auto start = line.find_first_not_of(" \t");
        if (start == std::string::npos || line.compare(start, 4, "note") != 0) continue;
        if (is_synthetic_note(line)) return true;
    }
    return false;
}

/// A pair id is a folder name: letters, digits, '_', '-', '.'.
bool valid_id(const std::string& id) {
    if (id.empty() || id == "." || id == "..") return false;
    return std::all_of(id.begin(), id.end(), [](char ch) {
        return std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '-' || ch == '.';
    });
}

/// Number cell; blank_is_nan shows NaN as an empty cell (pure-oil saturation pressure).
bool number_cell(const char* id, double& value, bool blank_is_nan) {
    char buffer[48] = "";
    if (std::isfinite(value) || !blank_is_nan) std::snprintf(buffer, sizeof(buffer), "%.6g", value);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (!ImGui::InputText(id, buffer, sizeof(buffer), ImGuiInputTextFlags_CharsScientific)) return false;
    if (buffer[0] == '\0') {
        if (!blank_is_nan) return false;
        value = std::numeric_limits<double>::quiet_NaN();
        return true;
    }
    char* end = nullptr;
    const double parsed = std::strtod(buffer, &end);
    if (end == buffer) return false;
    value = parsed;
    return true;
}

void set_message(LibraryEditor& state, std::string text, bool error) {
    state.message = std::move(text);
    state.message_error = error;
}

}  // namespace

void GuiApp::draw_fluid_library_dialog() {
    LibraryEditor& state = editor();
    const fs::path root = fluid_root();

    // Loads a pair into the editor; a table the solver rejects still opens for repair.
    const auto open_pair = [&](const std::string& id) {
        state.editing_id = id;
        state.is_new = false;
        state.modified = false;
        state.plotted_ok = false;
        state.curve_key.clear();
        std::string error;
        gui::FluidPair loaded;
        if (library_.load(id, loaded, error)) {
            state.pair = loaded;
            state.plotted = loaded;
            state.plotted_ok = true;
            std::string row_error;
            if (!gui::read_solution_rows(loaded.table_file, state.rows, row_error))
                set_message(state, "Table: " + row_error, true);
            else
                state.message.clear();
        } else {
            state.pair = loaded;
            state.pair.id = id;
            std::string row_error;
            gui::read_solution_rows(root / id / "solution.csv", state.rows, row_error);
            set_message(state, "The solver's loader rejected this pair: " + error, true);
        }
        std::snprintf(state.copy_name, sizeof(state.copy_name), "%s_copy", id.c_str());
    };
    const auto start_new = [&](const std::string& id, std::vector<gui::SolutionRow> rows) {
        state.editing_id = id;
        state.is_new = true;
        state.modified = true;
        state.pair = gui::FluidPair{};
        state.pair.id = id;
        state.pair.name = id;
        state.rows = std::move(rows);
        state.plotted_ok = false;
        state.curve_key.clear();
    };

    if (state.listed != pair_ids_) {
        state.listed = pair_ids_;
        state.synthetic.clear();
        for (const std::string& id : pair_ids_) state.synthetic[id] = pair_is_synthetic(root, id);
    }
    if (state.editing_id.empty() && !pair_ids_.empty()) {
        const bool has_case_pair =
            std::find(pair_ids_.begin(), pair_ids_.end(), case_.fluid_pair) != pair_ids_.end();
        open_pair(has_case_pair ? case_.fluid_pair : pair_ids_.front());
    }

    ImGui::SetNextWindowSize(ImVec2(1000, 650), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Fluid library", &show_fluid_library_)) {
        ImGui::End();
        return;
    }

    ImGui::AlignTextToFramePadding();
    ImGui::Text("Fluid directory: %s", winutil::path_to_utf8(root).c_str());
    ImGui::SameLine();
    if (ImGui::Button("Open folder")) {
        std::error_code ec;
        fs::create_directories(root, ec);
        ShellExecuteW(hwnd, L"open", root.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
    }
    ui::item_tooltip("Each pair is a folder with pair.txt (header) and solution.csv (table). "
                     "The folder can be changed in Settings.");
    if (!state.message.empty()) {
        if (state.message_error) ImGui::PushStyleColor(ImGuiCol_Text, theme::kError);
        ImGui::TextWrapped("%s", state.message.c_str());
        if (state.message_error) ImGui::PopStyleColor();
    }
    ImGui::Separator();

    // ---- Left: pairs -------------------------------------------------------
    const float font = ImGui::GetFontSize();
    ImGui::BeginChild("##pairs", ImVec2(16.0f * font, 0), ImGuiChildFlags_Borders);
    ImGui::TextUnformatted("Oil-refrigerant pairs");
    if (ImGui::Button("Refresh")) {
        refresh_pairs();
        state.listed.clear();
    }
    ui::item_tooltip("Read the fluid directory again.");
    ImGui::Separator();
    for (const std::string& id : pair_ids_) {
        const bool synthetic = state.synthetic[id];
        std::string label = id;
        if (synthetic) label += "  (synthetic)";
        if (id == case_.fluid_pair) label += "  [case]";
        if (synthetic) ImGui::PushStyleColor(ImGuiCol_Text, theme::kRiskMedium);
        if (ImGui::Selectable(label.c_str(), id == state.editing_id && !state.is_new)) open_pair(id);
        if (synthetic) ImGui::PopStyleColor();
        if (synthetic) ui::item_tooltip("Synthetic example data, not for design.");
    }
    if (state.is_new) {
        ImGui::PushStyleColor(ImGuiCol_Text, theme::kMuted);
        ImGui::Selectable((state.editing_id + "  (new, not saved)").c_str(), true);
        ImGui::PopStyleColor();
    }
    ImGui::Separator();

    ImGui::TextUnformatted("New pair");
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputText("##new_name", state.new_name, sizeof(state.new_name));
    ui::item_tooltip("Folder name of the new pair (letters, digits, _ - .).");
    const auto exists = [&](const std::string& id) {
        return std::find(pair_ids_.begin(), pair_ids_.end(), id) != pair_ids_.end();
    };
    if (ImGui::Button("New")) {
        const std::string id = state.new_name;
        if (!valid_id(id)) set_message(state, "Use letters, digits, _ - . for the pair name.", true);
        else if (exists(id)) set_message(state, "A pair named '" + id + "' already exists.", true);
        else {
            start_new(id, {});
            set_message(state, "New pair: fill in the header, paste the table, then Save.", false);
        }
    }
    ui::item_tooltip("Start an empty pair; nothing is written until Save.");
    ImGui::SameLine();
    if (ImGui::Button("Import CSV...")) {
        const fs::path file = exports::open_file_dialog(
            hwnd, L"Solution tables (*.csv)\0*.csv\0All files\0*.*\0");
        if (!file.empty()) {
            std::vector<gui::SolutionRow> rows;
            std::string error;
            if (!gui::read_solution_rows(file, rows, error)) {
                set_message(state, "Import failed: " + error, true);
            } else {
                std::string id = winutil::path_to_utf8(file.stem());
                if (!valid_id(id) || exists(id)) id = state.new_name;
                start_new(id, std::move(rows));
                set_message(state, "Imported " + std::to_string(state.rows.size()) +
                                       " rows into a new pair; fill in the header and Save.",
                            false);
            }
        }
    }
    ui::item_tooltip("Start a new pair from a CSV in solution.csv format: content [kg/kg], "
                     "temperature [degC], saturation pressure [MPa abs], viscosity [mm2/s].");

    ImGui::Spacing();
    ImGui::TextUnformatted("Copy selected as");
    ImGui::SetNextItemWidth(-FLT_MIN);
    ImGui::InputText("##copy_name", state.copy_name, sizeof(state.copy_name));
    ImGui::BeginDisabled(state.is_new || state.editing_id.empty());
    if (ImGui::Button("Copy")) {
        const std::string id = state.copy_name;
        std::string error;
        if (!valid_id(id)) set_message(state, "Use letters, digits, _ - . for the pair name.", true);
        else if (!library_.copy(state.editing_id, id, error)) set_message(state, "Copy failed: " + error, true);
        else {
            refresh_pairs();
            open_pair(id);
            set_message(state, "Copied to '" + id + "'.", false);
        }
    }
    ImGui::EndDisabled();
    ui::item_tooltip("Copy the selected pair folder under a new name, e.g. to edit a variant.");
    ImGui::EndChild();

    ImGui::SameLine();

    // ---- Right: the pair being edited --------------------------------------
    ImGui::BeginChild("##pair_editor", ImVec2(0, 0));
    if (state.editing_id.empty()) {
        ImGui::TextDisabled("No pair selected. Create one with New or Import CSV.");
        ImGui::EndChild();
        ImGui::End();
        return;
    }
    gui::FluidPair& pair = state.pair;
    ImGui::PushFont(fonts.heading);
    ImGui::Text("%s%s", state.editing_id.c_str(), state.modified ? " *" : "");
    ImGui::PopFont();
    if (is_synthetic_note(pair.note)) {
        ImGui::SameLine();
        ImGui::TextColored(theme::kRiskMedium, "Synthetic example data, not for design.");
    }

    // Header
    if (ImGui::CollapsingHeader("Header", ImGuiTreeNodeFlags_DefaultOpen) && ui::begin_form("##header")) {
        bool changed = false;
        changed |= form_text("Name", "Name shown in the case pair list.", "name", pair.name);
        changed |= form_text("Oil", "Oil, e.g. grade and base stock.", "oil", pair.oil);
        changed |= form_text("Refrigerant", "Refrigerant dissolved in the oil.", "refrigerant",
                             pair.refrigerant);
        changed |= form_text("Source", "Where the data come from (measurement, paper, supplier sheet).",
                             "source", pair.source);
        changed |= form_text("Note", "Any limitation, e.g. 'example only, not for design'.", "note",
                             pair.note);
        changed |= ui::form_combo(
            "Gas model for gas transport",
            "Equation of state of the free gas used when the gas is transported with the film.",
            "fluid_lib.species", pair.solver_species,
            {{DissolvedGasSpecies::PROPANE, "PROPANE", "Propane (R290) ideal-gas properties."},
             {DissolvedGasSpecies::AIR, "AIR", "Air ideal-gas properties."},
             {DissolvedGasSpecies::R410A, "R410A", "R410A ideal-gas properties."}});
        auto& units_map = units();
        changed |= ui::form_number("Density at 0 degC [kg/m3]",
                                   "Liquid density at 0 degC; the density is linear in temperature. "
                                   "Typical oil: 850 to 1000 kg/m3.",
                                   "fluid_lib.rho0", pair.density_intercept, gui::UnitFamily::Fixed,
                                   units_map);
        changed |= ui::form_number("Density slope [kg/(m3 K)]",
                                   "Change of density per kelvin, negative. Typical: -0.6 to -0.75.",
                                   "fluid_lib.rho_slope", pair.density_slope, gui::UnitFamily::Fixed,
                                   units_map);
        changed |= ui::form_number("Oil vapour pressure",
                                   "Pressure at which the oil itself boils (vaporous cavitation). "
                                   "Typical: a few hundred Pa.",
                                   "fluid_lib.p_vap", pair.oil_vapour_pressure, gui::UnitFamily::Pressure,
                                   units_map);
        changed |= ui::form_number("Specific heat [J/(kg K)]", "Liquid specific heat. Typical oil: about 2000.",
                                   "fluid_lib.cp", pair.specific_heat, gui::UnitFamily::Fixed, units_map);
        changed |= ui::form_number("Thermal conductivity", "Liquid thermal conductivity. Typical oil: 0.13.",
                                   "fluid_lib.k", pair.thermal_conductivity, gui::UnitFamily::Conductivity,
                                   units_map);
        ui::end_form();
        const bool optional_open = ImGui::TreeNodeEx("Optional: gas transport and thermal gas release",
                                                     ImGuiTreeNodeFlags_SpanAvailWidth);
        ui::item_tooltip("Free-gas viscosity and caloric data; needed only for gas transport and "
                         "thermal gas release.");
        if (optional_open) {
            if (ui::begin_form("##header_optional")) {
                changed |= ui::form_number("Free-gas viscosity (optional)",
                                           "Dynamic viscosity of the released gas, used with gas transport. "
                                           "0 = not given. Typical: 1e-5 Pa s.",
                                           "fluid_lib.mu_gas", pair.gas_viscosity, gui::UnitFamily::Viscosity,
                                           units_map);
                changed |= ui::form_number("Dissolved heat capacity [J/(kg K)] (optional)",
                                           "Heat capacity of the dissolved refrigerant, for thermal gas release. "
                                           "0 = not given.",
                                           "fluid_lib.cp_dissolved", pair.dissolved_heat_capacity,
                                           gui::UnitFamily::Fixed, units_map);
                double desorption = std::isfinite(pair.desorption_energy) ? pair.desorption_energy : 0.0;
                if (ui::form_number("Desorption energy [J/kg] (optional)",
                                    "Energy taken up when dissolved refrigerant comes out of solution "
                                    "(u_gas - u_dissolved), for thermal gas release. 0 = not given.",
                                    "fluid_lib.desorption", desorption, gui::UnitFamily::Fixed, units_map)) {
                    pair.desorption_energy = desorption != 0.0 ? desorption : std::numeric_limits<double>::quiet_NaN();
                    changed = true;
                }
                changed |= ui::form_number("Free-gas thermal conductivity (optional)",
                                           "Thermal conductivity of the released gas, for thermal gas release. "
                                           "Saved together with the two values above.",
                                           "fluid_lib.k_gas", pair.gas_thermal_conductivity,
                                           gui::UnitFamily::Conductivity, units_map);
                ui::end_form();
            }
            ImGui::TreePop();
        }
        if (changed) state.modified = true;
    }

    // Table
    if (ImGui::CollapsingHeader("Table (solution.csv)", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("Each content curve lists the temperature at which that content is saturated "
                           "at the given pressure, with the viscosity there.");
        ImGui::PopStyleColor();
        if (ImGui::Button("Add row")) {
            state.rows.push_back(state.rows.empty() ? gui::SolutionRow{} : state.rows.back());
            state.modified = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Sort")) {
            std::stable_sort(state.rows.begin(), state.rows.end(),
                             [](const gui::SolutionRow& a, const gui::SolutionRow& b) {
                                 if (a.content != b.content) return a.content < b.content;
                                 return a.temperature_c < b.temperature_c;
                             });
            state.modified = true;
        }
        ui::item_tooltip("Sort by content, then temperature.");
        ImGui::SameLine();
        if (ImGui::Button("Paste from Excel")) {
            const char* clipboard = ImGui::GetClipboardText();
            std::vector<gui::SolutionRow> pasted;
            std::string error;
            if (!clipboard || !*clipboard) {
                set_message(state, "The clipboard is empty.", true);
            } else if (!gui::parse_pasted_rows(clipboard, pasted, error)) {
                set_message(state, "Paste failed: " + error, true);
            } else {
                state.rows = std::move(pasted);
                state.modified = true;
                set_message(state, "Pasted " + std::to_string(state.rows.size()) + " rows.", false);
            }
        }
        ui::item_tooltip("Replace the rows with four columns copied from a spreadsheet, in the units of "
                         "the column titles. A header line is allowed.");
        ImGui::SameLine();
        ImGui::TextDisabled("%d rows", static_cast<int>(state.rows.size()));

        const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                      ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp;
        if (ImGui::BeginTable("##rows", 5, flags, ImVec2(0, 14.0f * font))) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Content [kg/kg]");
            ImGui::TableSetupColumn("Temperature [degC]");
            ImGui::TableSetupColumn("Saturation pressure [MPa, abs; blank for pure oil]");
            ImGui::TableSetupColumn("Viscosity [mm2/s]");
            ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 1.6f * font);
            ImGui::TableHeadersRow();
            int remove = -1;
            ImGuiListClipper clipper;
            clipper.Begin(static_cast<int>(state.rows.size()));
            while (clipper.Step()) {
                for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; ++i) {
                    gui::SolutionRow& row = state.rows[static_cast<std::size_t>(i)];
                    ImGui::TableNextRow();
                    ImGui::PushID(i);
                    ImGui::TableSetColumnIndex(0);
                    if (number_cell("##c", row.content, false)) state.modified = true;
                    ImGui::TableSetColumnIndex(1);
                    if (number_cell("##t", row.temperature_c, false)) state.modified = true;
                    ImGui::TableSetColumnIndex(2);
                    if (number_cell("##p", row.pressure_mpa, true)) state.modified = true;
                    ImGui::TableSetColumnIndex(3);
                    if (number_cell("##nu", row.viscosity_mm2s, false)) state.modified = true;
                    ImGui::TableSetColumnIndex(4);
                    if (ImGui::SmallButton("x")) remove = i;
                    ui::item_tooltip("Delete this row.");
                    ImGui::PopID();
                }
            }
            ImGui::EndTable();
            if (remove >= 0) {
                state.rows.erase(state.rows.begin() + remove);
                state.modified = true;
            }
        }

        if (ImGui::Button("Save", ImVec2(7.0f * font, 0))) {
            pair.id = state.editing_id;
            std::string error;
            if (!library_.save(pair, state.rows, error)) {
                set_message(state, "Save failed or the solver's loader rejected the table: " + error, true);
                state.plotted_ok = false;
            } else {
                state.is_new = false;
                state.modified = false;
                state.curve_key.clear();
                std::string load_error;
                state.plotted_ok = library_.load(state.editing_id, state.plotted, load_error);
                set_message(state, "Saved '" + state.editing_id + "'; the solver's loader accepted the table.",
                            false);
            }
            refresh_pairs();
            if (state.editing_id == case_.fluid_pair) {
                select_pair(state.editing_id);
                refresh_derived();
            }
        }
        ui::item_tooltip("Write pair.txt and solution.csv, then read them back with the solver's loader.");
        ImGui::SameLine();
        if (ImGui::Button("Revert", ImVec2(7.0f * font, 0))) {
            if (state.is_new) {
                start_new(state.editing_id, {});
                set_message(state, "Cleared the new pair.", false);
            } else {
                open_pair(state.editing_id);
                set_message(state, "Reverted to the saved files.", false);
            }
        }
        ui::item_tooltip("Discard the edits and read the saved files again.");
    }

    // Plots of what the solver reads
    if (ImGui::CollapsingHeader("Plots", ImGuiTreeNodeFlags_DefaultOpen)) {
        const bool case_pair = state.editing_id == case_.fluid_pair;
        if (case_pair && !estimates_.fluid_in_range)
            ImGui::TextColored(theme::kRiskMedium,
                               "The case state (ambient temperature, content) lies outside this table; "
                               "properties are taken at the table edge.");
        if (!state.plotted_ok) {
            ImGui::TextDisabled("Plots appear after the pair loads or saves without errors.");
        } else {
            const gui::FluidPair& shown = state.plotted;
            const double t_amb = case_.ambient_temperature;
            const float plot_w = 0.5f * (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x);
            const ImVec2 plot_size(std::max(plot_w, 150.0f), 18.0f * font);

            if (ImPlot::BeginPlot("Viscosity vs temperature", plot_size, ImPlotFlags_NoBoxSelect)) {
                ImPlot::SetupAxes("Temperature [degC]", "Kinematic viscosity [mm2/s]", ImPlotAxisFlags_AutoFit,
                                  ImPlotAxisFlags_AutoFit);
                ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
                const auto& table = shown.viscosity;
                for (std::size_t i = 0; i < table.curves.size() && i < table.axis.size(); ++i) {
                    std::vector<double> xs, ys;
                    for (const auto& point : table.curves[i]) {
                        xs.push_back(point.x - 273.15);
                        ys.push_back(point.value * 1e6);
                    }
                    char label[48];
                    std::snprintf(label, sizeof(label), "%.3g wt%%", table.axis[i] * 100.0);
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 2.5f);
                    ImPlot::PlotLine(label, xs.data(), ys.data(), static_cast<int>(xs.size()));
                }
                if (case_pair) {
                    const double t_c = t_amb - 273.15;
                    ImPlot::SetNextLineStyle(theme::kMuted, 1.0f);
                    ImPlot::PlotInfLines("Ambient temperature", &t_c, 1);
                }
                ImPlot::EndPlot();
            }
            ImGui::SameLine();

            // Saturation pressure vs content at the case ambient temperature (cached).
            char key[160];
            std::snprintf(key, sizeof(key), "%s|%.6g", shown.id.c_str(), t_amb);
            if (state.curve_key != key) {
                state.curve_key = key;
                state.curve_content.clear();
                state.curve_pressure.clear();
                const double c_lo = shown.content_min(), c_hi = shown.content_max();
                const int samples = 60;
                for (int i = 0; i <= samples && c_hi > c_lo; ++i) {
                    const double c = c_lo + (c_hi - c_lo) * i / samples;
                    try {
                        const gui::PropertyValue p = shown.saturation_pressure(c, t_amb);
                        if (!std::isfinite(p.value)) continue;
                        state.curve_content.push_back(c * 100.0);
                        state.curve_pressure.push_back(p.value * 1e-5);
                    } catch (const std::exception&) {
                    }
                }
            }
            char title[96];
            std::snprintf(title, sizeof(title), "Saturation pressure at %.1f degC###p_sat_plot",
                          t_amb - 273.15);
            if (ImPlot::BeginPlot(title, plot_size, ImPlotFlags_NoBoxSelect)) {
                ImPlot::SetupAxes("Refrigerant content [wt%]", "Saturation pressure [bar, abs]",
                                  ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                ImPlot::SetupLegend(ImPlotLocation_NorthWest, ImPlotLegendFlags_None);
                ImPlot::SetNextLineStyle(theme::kAccent, 2.0f);
                ImPlot::PlotLine("p_sat (case ambient temperature)", state.curve_content.data(),
                                 state.curve_pressure.data(), static_cast<int>(state.curve_content.size()));
                if (case_pair) {
                    const double p_amb = case_.ambient_pressure * 1e-5;
                    ImPlot::SetNextLineStyle(theme::kMuted, 1.5f);
                    ImPlot::PlotInfLines("Ambient pressure", &p_amb, 1, ImPlotInfLinesFlags_Horizontal);
                    if (std::isfinite(estimates_.content) && std::isfinite(estimates_.saturation_pressure)) {
                        const double cx = estimates_.content * 100.0;
                        const double py = estimates_.saturation_pressure * 1e-5;
                        const ImVec4 dark(0.10f, 0.11f, 0.13f, 1.0f);
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Diamond, 7.0f, dark, 1.5f, dark);
                        ImPlot::PlotScatter("Case content", &cx, &py, 1);
                    }
                }
                ImPlot::EndPlot();
            }
            if (!case_pair)
                ImGui::TextDisabled("Ambient lines and the case marker appear for the pair the case uses (%s).",
                                    case_.fluid_pair.empty() ? "none" : case_.fluid_pair.c_str());
        }
    }

    ImGui::EndChild();
    ImGui::End();
}
