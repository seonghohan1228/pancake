// Case Setup panel: run controls, validation summary, schema-generated
// parameter form (with units, live validation, common/advanced split),
// enum/model selectors, inlet editor, and the raw config.txt editor.

#include "gui_app.hpp"
#include "imgui.h"
#include "win_util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>

namespace {

// InputTextMultiline over std::string (resize callback, as in imgui_stdlib).
int string_resize_callback(ImGuiInputTextCallbackData* data) {
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* text = static_cast<std::string*>(data->UserData);
        text->resize(static_cast<size_t>(data->BufTextLen));
        data->Buf = text->data();
    }
    return 0;
}

bool input_text_multiline_string(const char* label, std::string& text, const ImVec2& size) {
    return ImGui::InputTextMultiline(label, text.data(), text.capacity() + 1, size,
                                     ImGuiInputTextFlags_CallbackResize | ImGuiInputTextFlags_AllowTabInput,
                                     string_resize_callback, &text);
}

bool input_text_string(const char* label, std::string& text) {
    return ImGui::InputText(label, text.data(), text.capacity() + 1,
                            ImGuiInputTextFlags_CallbackResize, string_resize_callback, &text);
}

bool contains_ci(const char* haystack, const char* needle) {
    if (!needle || !*needle) return true;
    const size_t haystack_len = std::strlen(haystack);
    const size_t needle_len = std::strlen(needle);
    if (needle_len > haystack_len) return false;
    for (size_t i = 0; i + needle_len <= haystack_len; ++i) {
        size_t j = 0;
        while (j < needle_len &&
               std::tolower(static_cast<unsigned char>(haystack[i + j])) ==
                   std::tolower(static_cast<unsigned char>(needle[j]))) {
            ++j;
        }
        if (j == needle_len) return true;
    }
    return false;
}

void help_tooltip(const char* help, const char* unit) {
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled | ImGuiHoveredFlags_ForTooltip)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26);
        ImGui::TextUnformatted(help);
        if (unit && std::strcmp(unit, "-") != 0) {
            ImGui::TextDisabled("Unit: %s", unit);
        }
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

/// Generic enum combo over explicit (value, label) options.
template <typename Enum>
bool enum_combo(const char* label, Enum& value,
                std::initializer_list<std::pair<Enum, const char*>> options,
                const char* help) {
    const char* preview = "?";
    for (const auto& [option, text] : options) {
        if (option == value) preview = text;
    }
    bool changed = false;
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo(label, preview)) {
        for (const auto& [option, text] : options) {
            if (ImGui::Selectable(text, option == value)) {
                value = option;
                changed = true;
            }
        }
        ImGui::EndCombo();
    }
    help_tooltip(help, nullptr);
    return changed;
}

/// Two-column row: label + tooltip on the left, widget on the right.
void begin_param_row(const char* label, const char* help, const char* unit, bool invalid,
                     const std::vector<ValidationIssue>& issues, const char* key) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    if (invalid) {
        ImGui::TextColored(ImVec4(0.93f, 0.36f, 0.36f, 1.0f), "%s", label);
    } else {
        ImGui::TextUnformatted(label);
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 26);
        ImGui::TextUnformatted(help);
        if (unit && std::strcmp(unit, "-") != 0) ImGui::TextDisabled("Unit: %s", unit);
        if (invalid) {
            for (const ValidationIssue& issue : issues) {
                if (issue.field == key && issue.is_error) {
                    ImGui::TextColored(ImVec4(0.93f, 0.36f, 0.36f, 1.0f), "%s",
                                       issue.message.c_str());
                }
            }
        }
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
    ImGui::TableSetColumnIndex(1);
}

}  // namespace

void GuiApp::draw_setup_panel() {
    if (!ImGui::Begin("Case Setup")) {
        setup_panel_focused_ = false;
        ImGui::End();
        return;
    }
    setup_panel_focused_ = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);

    // --- Run bar -----------------------------------------------------------
    const bool running = runner_.running();
    const bool can_run = !has_errors_ && !running;
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.13f, 0.55f, 0.33f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.16f, 0.63f, 0.39f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.11f, 0.47f, 0.29f, 1.0f));
    ImGui::BeginDisabled(!can_run);
    const float run_width = (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x) * 0.62f;
    if (ImGui::Button(running ? "Running..." : "Run  (F5)", ImVec2(run_width, 0))) start_run();
    ImGui::EndDisabled();
    ImGui::PopStyleColor(3);
    if (!can_run && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        ImGui::SetTooltip("%s", running ? "A solve is already running."
                                        : (first_error() ? first_error() : ""));
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(!running);
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.62f, 0.22f, 0.24f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.72f, 0.27f, 0.29f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.55f, 0.18f, 0.20f, 1.0f));
    if (ImGui::Button("Cancel", ImVec2(-FLT_MIN, 0))) runner_.cancel();
    ImGui::PopStyleColor(3);
    ImGui::EndDisabled();

    if (has_errors_) {
        ImGui::TextColored(theme::kError, "Run disabled: %s", first_error());
    } else {
        int warning_count = 0;
        for (const ValidationIssue& issue : issues_) {
            if (!issue.is_error) ++warning_count;
        }
        if (warning_count > 0) {
            ImGui::TextColored(theme::kWarn, "Ready to run - %d warning%s (hover fields for details).",
                               warning_count, warning_count == 1 ? "" : "s");
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                ImGui::BeginTooltip();
                for (const ValidationIssue& issue : issues_) {
                    if (!issue.is_error) ImGui::BulletText("%s", issue.message.c_str());
                }
                ImGui::EndTooltip();
            }
        } else {
            ImGui::TextColored(theme::kOk, "Ready to run.");
        }
    }
    ImGui::Separator();

    if (ImGui::BeginTabBar("setup_tabs")) {
        if (ImGui::BeginTabItem("Parameters")) {
            draw_form_tab();
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(raw_dirty_ ? "Raw config.txt *###rawtab" : "Raw config.txt###rawtab")) {
            draw_raw_tab();
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

void GuiApp::draw_form_tab() {
    // Toolbar: search, advanced toggle, presets.
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x * 0.45f);
    ImGui::InputTextWithHint("##search", "Filter parameters...", search_buffer_,
                             sizeof(search_buffer_));
    ImGui::SameLine();
    ImGui::Checkbox("Advanced", &show_advanced_);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##presets", "Presets...")) {
        if (ImGui::Selectable("Reset to defaults")) {
            const std::string output_dir = config_.output_dir;
            config_ = SimulationConfig{};
            config_.output_dir = output_dir;
            dirty_ = true;
        }
        if (ImGui::Selectable("High eccentricity (e = 0.9 c)")) {
            config_.e = 0.9 * config_.c;
            dirty_ = true;
        }
        if (ImGui::Selectable("Misaligned journal")) {
            config_.tilt_slope_x = 0.5 * config_.c / std::max(config_.L, 1e-12);
            dirty_ = true;
        }
        if (ImGui::Selectable("Add circular oil hole inlet")) {
            InletConfig inlet;
            inlet.type = InletConfig::Type::CIRCULAR;
            inlet.theta = 180.0;
            inlet.z = 0.5 * config_.L;
            inlet.size = 0.1 * config_.L;
            inlet.p_supply = 2e5;
            config_.inlets.push_back(inlet);
            dirty_ = true;
        }
        ImGui::EndCombo();
    }

    ImGui::BeginChild("##form_scroll", ImVec2(0, 0), ImGuiChildFlags_None);
    const char* filter = search_buffer_;
    const bool filtering = filter[0] != '\0';

    static const ParamGroup kGroupOrder[] = {
        ParamGroup::Geometry, ParamGroup::Operating, ParamGroup::Lubricant,
        ParamGroup::MeshTime, ParamGroup::Cavitation, ParamGroup::Boundaries,
        ParamGroup::Thermal, ParamGroup::FluidProperties, ParamGroup::Motion,
        ParamGroup::Output, ParamGroup::Numerics,
    };
    // Advanced-only groups are hidden entirely unless toggled or searched.
    const auto group_is_advanced = [](ParamGroup group) {
        return group == ParamGroup::FluidProperties || group == ParamGroup::Numerics;
    };

    for (ParamGroup group : kGroupOrder) {
        // Does anything in this group survive the filter / advanced toggle?
        const auto spec_visible = [&](const char* label, const char* key, bool advanced) {
            if (filtering) return contains_ci(label, filter) || contains_ci(key, filter);
            return show_advanced_ || !advanced;
        };
        bool any_visible = false;
        for (const auto& spec : double_specs()) {
            if (spec.group == group && spec_visible(spec.label, spec.key, spec.advanced)) {
                any_visible = true;
            }
        }
        for (const auto& spec : int_specs()) {
            if (spec.group == group && spec_visible(spec.label, spec.key, spec.advanced)) {
                any_visible = true;
            }
        }
        for (const auto& spec : bool_specs()) {
            if (spec.group == group && spec_visible(spec.label, spec.key, spec.advanced)) {
                any_visible = true;
            }
        }
        // Groups with bespoke widgets always have content when not filtering.
        const bool has_bespoke = group == ParamGroup::MeshTime || group == ParamGroup::Cavitation ||
                                 group == ParamGroup::Boundaries || group == ParamGroup::Thermal ||
                                 group == ParamGroup::FluidProperties ||
                                 group == ParamGroup::Motion || group == ParamGroup::Output ||
                                 group == ParamGroup::Numerics;
        if (!any_visible && !(has_bespoke && !filtering)) continue;
        if (!filtering && !show_advanced_ && group_is_advanced(group) && !any_visible) continue;

        ImGui::PushFont(fonts.heading);
        const bool open = ImGui::CollapsingHeader(group_label(group),
                                                  filtering ? ImGuiTreeNodeFlags_DefaultOpen
                                                            : ImGuiTreeNodeFlags_DefaultOpen);
        ImGui::PopFont();
        if (!open) continue;

        if (ImGui::BeginTable((std::string("##params_") + group_label(group)).c_str(), 2,
                              ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.48f);
            ImGui::TableSetupColumn("widget", ImGuiTableColumnFlags_WidthStretch, 0.52f);

            for (const auto& spec : double_specs()) {
                if (spec.group != group || !spec_visible(spec.label, spec.key, spec.advanced)) {
                    continue;
                }
                const bool invalid = field_invalid(spec.key);
                begin_param_row(spec.label, spec.help, spec.unit, invalid, issues_, spec.key);
                if (invalid) ImGui::PushStyleColor(ImGuiCol_FrameBg, theme::kInvalidBg);
                const float unit_width = ImGui::CalcTextSize(spec.unit).x + 8.0f;
                ImGui::SetNextItemWidth(-unit_width);
                double value = config_.*(spec.member);
                if (ImGui::InputDouble((std::string("##") + spec.key).c_str(), &value, 0.0, 0.0,
                                       "%g")) {
                    config_.*(spec.member) = value;
                    dirty_ = true;
                }
                help_tooltip(spec.help, spec.unit);
                if (invalid) ImGui::PopStyleColor();
                ImGui::SameLine();
                ImGui::TextDisabled("%s", spec.unit);
            }
            for (const auto& spec : int_specs()) {
                if (spec.group != group || !spec_visible(spec.label, spec.key, spec.advanced)) {
                    continue;
                }
                const bool invalid = field_invalid(spec.key);
                begin_param_row(spec.label, spec.help, spec.unit, invalid, issues_, spec.key);
                if (invalid) ImGui::PushStyleColor(ImGuiCol_FrameBg, theme::kInvalidBg);
                ImGui::SetNextItemWidth(-FLT_MIN);
                int value = config_.*(spec.member);
                if (ImGui::InputInt((std::string("##") + spec.key).c_str(), &value, 0, 0)) {
                    config_.*(spec.member) = value;
                    dirty_ = true;
                }
                help_tooltip(spec.help, spec.unit);
                if (invalid) ImGui::PopStyleColor();
            }
            for (const auto& spec : bool_specs()) {
                if (spec.group != group || !spec_visible(spec.label, spec.key, spec.advanced)) {
                    continue;
                }
                begin_param_row(spec.label, spec.help, "-", false, issues_, spec.key);
                bool value = config_.*(spec.member);
                if (ImGui::Checkbox((std::string("##") + spec.key).c_str(), &value)) {
                    config_.*(spec.member) = value;
                    dirty_ = true;
                }
                help_tooltip(spec.help, nullptr);
            }
            ImGui::EndTable();
        }

        if (!filtering) draw_enum_widgets(group);
        if (group == ParamGroup::Boundaries && !filtering) draw_inlets_editor();
    }
    ImGui::EndChild();
}

void GuiApp::draw_enum_widgets(ParamGroup group) {
    bool changed = false;
    if (ImGui::BeginTable((std::string("##enums_") + group_label(group)).c_str(), 2,
                          ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthStretch, 0.48f);
        ImGui::TableSetupColumn("widget", ImGuiTableColumnFlags_WidthStretch, 0.52f);
        const auto row = [&](const char* label, const char* help) {
            begin_param_row(label, help, nullptr, false, issues_, "");
        };

        switch (group) {
            case ParamGroup::MeshTime:
                row("Solution mode", "TRANSIENT marches end_t/dt; STEADY_STATE solves one step.");
                changed |= enum_combo("##solution_mode", config_.solution_mode,
                                      {{SolutionMode::TRANSIENT, "Transient"},
                                       {SolutionMode::STEADY_STATE, "Steady state"}},
                                      "Solution mode");
                break;
            case ParamGroup::Cavitation:
                row("Cavitation model",
                    "ELROD_ADAMS is mass-conserving (JFO); GUMBEL clamps negative pressures; "
                    "FULL_SOMMERFELD keeps them.");
                changed |= enum_combo("##cavitation_model", config_.cavitation_model,
                                      {{CavitationModel::ELROD_ADAMS, "Elrod-Adams (JFO)"},
                                       {CavitationModel::GUMBEL, "Gumbel (half Sommerfeld)"},
                                       {CavitationModel::FULL_SOMMERFELD, "Full Sommerfeld"}},
                                      "Cavitation treatment");
                break;
            case ParamGroup::Boundaries:
                row("South boundary type (z = 0)", "Pressure condition at the south face.");
                changed |= enum_combo("##bc_s", config_.bc_z_south_type,
                                      {{BCType::DIRICHLET, "Fixed pressure (Dirichlet)"},
                                       {BCType::NEUMANN, "Zero gradient (Neumann)"},
                                       {BCType::INLET_OUTLET, "Inlet/outlet"}},
                                      "South BC");
                row("South thermal inflow", "Temperature carried by oil entering at z = 0.");
                changed |= enum_combo("##bc_s_th", config_.bc_z_south_thermal,
                                      {{ThermalInflowMode::OPEN, "Open (recirculating oil)"},
                                       {ThermalInflowMode::RESERVOIR, "Reservoir (fixed T)"}},
                                      "South thermal inflow");
                row("North boundary type (z = L)", "Pressure condition at the north face.");
                changed |= enum_combo("##bc_n", config_.bc_z_north_type,
                                      {{BCType::DIRICHLET, "Fixed pressure (Dirichlet)"},
                                       {BCType::NEUMANN, "Zero gradient (Neumann)"},
                                       {BCType::INLET_OUTLET, "Inlet/outlet"}},
                                      "North BC");
                row("North thermal inflow", "Temperature carried by oil entering at z = L.");
                changed |= enum_combo("##bc_n_th", config_.bc_z_north_thermal,
                                      {{ThermalInflowMode::OPEN, "Open (recirculating oil)"},
                                       {ThermalInflowMode::RESERVOIR, "Reservoir (fixed T)"}},
                                      "North thermal inflow");
                break;
            case ParamGroup::Thermal:
                row("Temperature model",
                    "ISOTHERMAL keeps properties at the initial temperature; ENERGY_EQUATION "
                    "transports heat in the film.");
                changed |= enum_combo("##temperature_model", config_.temperature_model,
                                      {{TemperatureModel::ISOTHERMAL, "Isothermal"},
                                       {TemperatureModel::ENERGY_EQUATION, "Energy equation (THD)"}},
                                      "Thermal model");
                break;
            case ParamGroup::FluidProperties:
                if (show_advanced_) {
                    row("Property model", "How density/viscosity react to pressure, temperature "
                                          "and dissolved gas.");
                    changed |= enum_combo(
                        "##fluid_property_model", config_.fluid_property_model,
                        {{FluidPropertyModel::CONSTANT, "Constant"},
                         {FluidPropertyModel::OIL_DISSOLVED_GAS, "Oil + dissolved gas"},
                         {FluidPropertyModel::GAS_CAVITATION_MIXTURE, "Gas cavitation mixture"}},
                        "Fluid property model");
                    row("Dissolved gas species", "Gas dissolved in the lubricant.");
                    changed |= enum_combo("##dissolved_gas_species", config_.dissolved_gas_species,
                                          {{DissolvedGasSpecies::PROPANE, "Propane (R-290)"},
                                           {DissolvedGasSpecies::AIR, "Air"}},
                                          "Dissolved gas");
                    row("Solubility model", "Saturation law c_sat(p, T).");
                    changed |= enum_combo("##oil_gas_solution_model", config_.oil_gas_solution_model,
                                          {{OilGasSolutionModel::HENRY, "Henry's law"},
                                           {OilGasSolutionModel::BUNSEN, "Bunsen coefficient"},
                                           {OilGasSolutionModel::TABLE, "Table"}},
                                          "Solubility model");
                    row("Density model", "Liquid-solution density mixing rule.");
                    changed |= enum_combo("##density_model", config_.density_model,
                                          {{DensityModel::PURE_OIL, "Pure oil"},
                                           {DensityModel::MASS_VOLUME_MIXING, "Mass-volume mixing"},
                                           {DensityModel::TABLE, "Table"}},
                                          "Density model");
                    row("Viscosity model", "Liquid-solution viscosity mixing rule.");
                    changed |= enum_combo("##viscosity_model", config_.viscosity_model,
                                          {{ViscosityModel::PURE_OIL, "Pure oil"},
                                           {ViscosityModel::LOG_MIXING, "Log mixing"},
                                           {ViscosityModel::EMPIRICAL_CORRELATION, "Empirical correlation"},
                                           {ViscosityModel::TABLE, "Table"}},
                                          "Viscosity model");
                    row("Two-phase viscosity", "Effective viscosity once free gas appears.");
                    changed |= enum_combo(
                        "##gas_mixture_viscosity_model", config_.gas_mixture_viscosity_model,
                        {{GasMixtureViscosityModel::EINSTEIN_DILUTE, "Einstein (dilute)"},
                         {GasMixtureViscosityModel::DUKLER_VOID, "Dukler (void fraction)"},
                         {GasMixtureViscosityModel::MCADAMS_QUALITY, "McAdams (quality)"},
                         {GasMixtureViscosityModel::KRIEGER_DOUGHERTY, "Krieger-Dougherty"},
                         {GasMixtureViscosityModel::LINEAR_QUALITY, "Linear quality"}},
                        "Two-phase viscosity model");
                }
                break;
            case ParamGroup::Motion:
                row("Motion model",
                    "STATIC keeps the configured eccentricity; MOVING_BEARING integrates the "
                    "bearing equation of motion under the loads.");
                changed |= enum_combo("##motion_model", config_.motion_model,
                                      {{MotionModel::STATIC, "Static"},
                                       {MotionModel::MOVING_BEARING, "Moving bearing"}},
                                      "Motion model");
                break;
            case ParamGroup::Output: {
                row("Terminal verbosity",
                    "DETAILED adds per-step film/force/torque diagnostics to the log (needed for "
                    "the friction-torque result).");
                changed |= enum_combo("##output_verbosity", config_.output_verbosity,
                                      {{OutputVerbosity::SIMPLIFIED, "Simplified"},
                                       {OutputVerbosity::DETAILED, "Detailed"}},
                                      "Verbosity");
                row("Output directory", "Where VTK results land; relative paths resolve next to "
                                        "the config file.");
                if (input_text_string("##output_dir", config_.output_dir)) dirty_ = true;
                row("Filename prefix", "Prefix of all written VTK files.");
                if (input_text_string("##filename_prefix", config_.filename_prefix)) dirty_ = true;
                break;
            }
            case ParamGroup::Numerics:
                if (show_advanced_) {
                    row("Motion time integrator", "Integrator for the bearing equation of motion.");
                    changed |= enum_combo("##motion_time_method", config_.motion_time_method,
                                          {{TimeSteppingMethod::EULER_IMPLICIT, "Implicit Euler"},
                                           {TimeSteppingMethod::EULER_EXPLICIT, "Explicit Euler"},
                                           {TimeSteppingMethod::CRANK_NICOLSON, "Crank-Nicolson"},
                                           {TimeSteppingMethod::RK2, "Runge-Kutta 2"},
                                           {TimeSteppingMethod::RK4, "Runge-Kutta 4"}},
                                          "Motion integrator");
                    const auto scheme_combo = [&](const char* id, ConvectionScheme& scheme,
                                                  bool with_type_differencing) {
                        if (with_type_differencing) {
                            return enum_combo(id, scheme,
                                              {{ConvectionScheme::UPWIND, "Upwind"},
                                               {ConvectionScheme::TVD_VANLEER, "TVD van Leer"},
                                               {ConvectionScheme::TVD_MINMOD, "TVD minmod"},
                                               {ConvectionScheme::TYPE_DIFFERENCING,
                                                "Type differencing (V&K 1989)"}},
                                              "Convection face interpolation");
                        }
                        return enum_combo(id, scheme,
                                          {{ConvectionScheme::UPWIND, "Upwind"},
                                           {ConvectionScheme::TVD_VANLEER, "TVD van Leer"},
                                           {ConvectionScheme::TVD_MINMOD, "TVD minmod"}},
                                          "Convection face interpolation");
                    };
                    row("Film-content convection", "Scheme for the Elrod film-content transport.");
                    changed |= scheme_combo("##theta_scheme", config_.theta_convection_scheme, true);
                    row("Thermal convection", "Scheme for temperature transport.");
                    changed |= scheme_combo("##thermal_scheme", config_.thermal_convection_scheme,
                                            false);
                    row("Gas convection", "Scheme for dissolved-gas transport.");
                    changed |= scheme_combo("##gas_scheme", config_.gas_convection_scheme, false);
                }
                break;
            default: break;
        }
        ImGui::EndTable();
    }
    if (changed) dirty_ = true;
}

void GuiApp::draw_inlets_editor() {
    ImGui::PushFont(fonts.heading);
    ImGui::SeparatorText("Oil inlets");
    ImGui::PopFont();
    const bool inlets_invalid = field_invalid("inlets");
    if (inlets_invalid) {
        for (const ValidationIssue& issue : issues_) {
            if (issue.field == "inlets" && issue.is_error) {
                ImGui::TextColored(theme::kError, "%s", issue.message.c_str());
            }
        }
    }
    int remove_index = -1;
    for (size_t i = 0; i < config_.inlets.size(); ++i) {
        InletConfig& inlet = config_.inlets[i];
        ImGui::PushID(static_cast<int>(i));
        const bool is_groove = inlet.type == InletConfig::Type::GROOVE;
        std::string title = std::string(is_groove ? "Groove" : "Circular hole") + " inlet " +
                            std::to_string(i + 1);
        if (ImGui::TreeNodeEx(title.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
            int type_index = is_groove ? 1 : 0;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
            if (ImGui::Combo("Type", &type_index, "Circular hole\0Axial groove\0")) {
                inlet.type = type_index == 1 ? InletConfig::Type::GROOVE
                                             : InletConfig::Type::CIRCULAR;
                dirty_ = true;
            }
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
            if (ImGui::InputDouble("Position theta [deg]", &inlet.theta, 0, 0, "%g")) dirty_ = true;
            if (inlet.type == InletConfig::Type::CIRCULAR) {
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
                if (ImGui::InputDouble("Position z [m]", &inlet.z, 0, 0, "%g")) dirty_ = true;
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
                if (ImGui::InputDouble("Radius [m]", &inlet.size, 0, 0, "%g")) dirty_ = true;
            } else {
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
                if (ImGui::InputDouble("Angular width [deg]", &inlet.size, 0, 0, "%g")) dirty_ = true;
            }
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
            if (ImGui::InputDouble("Supply pressure [Pa]", &inlet.p_supply, 0, 0, "%g")) dirty_ = true;
            if (ImGui::Checkbox("Feeds fresh oil (fixed supply temperature)",
                                &inlet.feeds_fresh_oil)) {
                dirty_ = true;
            }
            if (inlet.feeds_fresh_oil) {
                ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
                if (ImGui::InputDouble("Supply temperature [K]", &inlet.t_supply, 0, 0, "%g")) {
                    dirty_ = true;
                }
            }
            if (ImGui::SmallButton("Remove inlet")) remove_index = static_cast<int>(i);
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
    if (remove_index >= 0) {
        config_.inlets.erase(config_.inlets.begin() + remove_index);
        dirty_ = true;
    }
    if (ImGui::Button("Add inlet")) {
        InletConfig inlet;
        inlet.type = InletConfig::Type::GROOVE;
        inlet.theta = 90.0;
        inlet.size = 10.0;
        inlet.p_supply = 2e5;
        config_.inlets.push_back(inlet);
        dirty_ = true;
    }
}

void GuiApp::draw_raw_tab() {
    ImGui::TextDisabled("Direct key = value editing; Apply parses it back into the form.");
    ImGui::SameLine();
    if (ImGui::SmallButton("Sync from form")) sync_raw_from_config();
    ImGui::SameLine();
    ImGui::BeginDisabled(!raw_dirty_);
    if (ImGui::SmallButton("Apply to form")) {
        SimulationConfig parsed;
        parsed.load_from_text(raw_text_);
        config_ = std::move(parsed);
        raw_dirty_ = false;
        dirty_ = true;
        for (const std::string& warning : config_.parse_warnings) {
            log_.push_back({"Config warning: " + warning, IM_COL32(230, 190, 80, 255)});
        }
    }
    ImGui::EndDisabled();
    ImGui::PushFont(fonts.mono);
    if (input_text_multiline_string("##rawtext", raw_text_, ImVec2(-FLT_MIN, -FLT_MIN))) {
        raw_dirty_ = true;
    }
    ImGui::PopFont();
}
