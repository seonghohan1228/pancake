// Case Setup panel: run controls, validation summary, schema-generated
// parameter form (with units, live validation, common/advanced split),
// enum/model selectors with hover descriptions, and the inlet editor.

#include "gui_app.hpp"
#include "imgui.h"
#include "win_util.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>

namespace {

// InputText over std::string (resize callback, as in imgui_stdlib).
int string_resize_callback(ImGuiInputTextCallbackData* data) {
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* text = static_cast<std::string*>(data->UserData);
        text->resize(static_cast<size_t>(data->BufTextLen));
        data->Buf = text->data();
    }
    return 0;
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

/// One choice of an enum combo, with a hover description that helps the user
/// decide (accuracy, robustness, when to use).
template <typename Enum>
struct EnumOption {
    Enum value;
    const char* label;
    const char* desc;
};

template <typename Enum>
bool enum_combo(const char* id, Enum& value,
                std::initializer_list<EnumOption<Enum>> options) {
    const char* preview = "?";
    for (const auto& option : options) {
        if (option.value == value) preview = option.label;
    }
    bool changed = false;
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo(id, preview)) {
        for (const auto& option : options) {
            if (ImGui::Selectable(option.label, option.value == value)) {
                value = option.value;
                changed = true;
            }
            if (option.desc && ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                ImGui::BeginTooltip();
                ImGui::PushTextWrapPos(ImGui::GetFontSize() * 24);
                ImGui::TextUnformatted(option.desc);
                ImGui::PopTextWrapPos();
                ImGui::EndTooltip();
            }
        }
        ImGui::EndCombo();
    }
    return changed;
}

/// Two-column row: label + tooltip on the left, widget on the right.
void begin_param_row(const char* label, const char* help, const char* unit, bool invalid,
                     const std::vector<ValidationIssue>& issues, const char* key) {
    ImGui::TableNextRow();
    ImGui::TableSetColumnIndex(0);
    ImGui::AlignTextToFramePadding();
    if (invalid) {
        ImGui::TextColored(theme::kError, "%s", label);
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
                    ImGui::TextColored(theme::kError, "%s",
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

    // MPI ranks live next to Run: more ranks = faster solve.
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("MPI processes");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
        ImGui::SetTooltip("Number of parallel solver processes. More processes solve faster on\n"
                          "multi-core machines (needs the Microsoft MPI runtime when > 1).");
    }
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
    ImGui::BeginDisabled(running);
    if (ImGui::InputInt("##mpi_ranks", &settings.mpi_ranks)) {
        settings.mpi_ranks = std::clamp(settings.mpi_ranks, 1, 256);
    }
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
    draw_form_tab();
    ImGui::End();
}

void GuiApp::draw_form_tab() {
    // Toolbar: search, advanced toggle, presets.
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x * 0.40f);
    ImGui::InputTextWithHint("##search", "Filter parameters...", search_buffer_,
                             sizeof(search_buffer_));
    ImGui::SameLine();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Advanced");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
        ImGui::SetTooltip("Show expert parameters (numerics, fluid-property models, ...).");
    }
    ImGui::SameLine();
    ImGui::Checkbox("##advanced", &show_advanced_);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##presets", "Presets...")) {
        if (ImGui::Selectable("Reset to defaults")) {
            const std::string output_dir = config_.output_dir;
            config_ = SimulationConfig{};
            // Runnable atmospheric defaults (see default_case in gui_app.cpp).
            config_.p_cav = 101325.0;
            config_.bc_z_south_val = 101325.0;
            config_.bc_z_north_val = 101325.0;
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

    // Display order = edit frequency: the parameters of a typical tweak-run
    // loop sit at the top; set-and-forget groups follow, collapsed by default.
    static const ParamGroup kGroupOrder[] = {
        ParamGroup::Operating, ParamGroup::Geometry, ParamGroup::Lubricant,
        ParamGroup::Cavitation, ParamGroup::TimeStepping, ParamGroup::Mesh,
        ParamGroup::Boundaries, ParamGroup::Thermal, ParamGroup::Motion,
        ParamGroup::FluidProperties, ParamGroup::Output, ParamGroup::Numerics,
    };
    const auto group_default_open = [](ParamGroup group) {
        switch (group) {
            case ParamGroup::Operating:
            case ParamGroup::Geometry:
            case ParamGroup::Lubricant:
            case ParamGroup::Cavitation:
            case ParamGroup::TimeStepping:
            case ParamGroup::Mesh:
                return true;
            default:
                return false;
        }
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
        const bool has_bespoke =
            group == ParamGroup::TimeStepping || group == ParamGroup::Cavitation ||
            group == ParamGroup::Boundaries || group == ParamGroup::Thermal ||
            group == ParamGroup::FluidProperties || group == ParamGroup::Motion ||
            group == ParamGroup::Output || group == ParamGroup::Numerics;
        if (!any_visible && !(has_bespoke && !filtering)) continue;
        if (!filtering && !show_advanced_ && group_is_advanced(group) && !any_visible) continue;

        // While filtering, force every surviving group open so matches show.
        if (filtering) ImGui::SetNextItemOpen(true, ImGuiCond_Always);
        ImGui::PushFont(fonts.heading);
        const bool open = ImGui::CollapsingHeader(
            group_label(group),
            group_default_open(group) ? ImGuiTreeNodeFlags_DefaultOpen : ImGuiTreeNodeFlags_None);
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

                const auto& options = unit_options(spec.family);
                int& choice = unit_choice(spec.key, spec.family);
                const double factor = options[static_cast<size_t>(choice)].to_display;

                if (invalid) ImGui::PushStyleColor(ImGuiCol_FrameBg, theme::kInvalidBg);
                const float unit_width =
                    options.size() == 1 ? ImGui::CalcTextSize(spec.unit).x + 8.0f
                                        : ImGui::GetFontSize() * 4.8f;
                ImGui::SetNextItemWidth(-unit_width);
                double display_value = config_.*(spec.member) * factor;
                if (ImGui::InputDouble((std::string("##") + spec.key).c_str(), &display_value,
                                       0.0, 0.0, "%g")) {
                    config_.*(spec.member) = display_value / factor;
                    dirty_ = true;
                }
                help_tooltip(spec.help, options.size() == 1 ? spec.unit : nullptr);
                if (invalid) ImGui::PopStyleColor();
                ImGui::SameLine();
                if (options.size() == 1) {
                    ImGui::TextDisabled("%s", spec.unit);
                } else {
                    // Display-unit selector; the config file stays solver-native.
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::BeginCombo((std::string("##unit_") + spec.key).c_str(),
                                          options[static_cast<size_t>(choice)].label)) {
                        for (int i = 0; i < static_cast<int>(options.size()); ++i) {
                            if (ImGui::Selectable(options[static_cast<size_t>(i)].label,
                                                  i == choice)) {
                                choice = i;
                            }
                        }
                        ImGui::EndCombo();
                    }
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                        ImGui::SetTooltip("Display unit only - the config file keeps %s.",
                                          spec.unit);
                    }
                }
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
            case ParamGroup::TimeStepping:
                row("Solution mode", "How the simulation advances in time.");
                changed |= enum_combo<SolutionMode>(
                    "##solution_mode", config_.solution_mode,
                    {{SolutionMode::TRANSIENT, "Transient",
                      "Marches from 0 to the end time in dt steps. Use for startup, dynamic "
                      "loads, moving bearings, or to let cavitation settle."},
                     {SolutionMode::STEADY_STATE, "Steady state",
                      "Solves a single equilibrium step. Fastest for fixed operating points "
                      "with a static journal."}});
                break;
            case ParamGroup::Cavitation:
                row("Cavitation model", "How film rupture is treated.");
                changed |= enum_combo<CavitationModel>(
                    "##cavitation_model", config_.cavitation_model,
                    {{CavitationModel::FULL_SOMMERFELD, "Full Sommerfeld (reference)",
                      "No cavitation: sub-cavity (tensile) pressures are kept. The simplest "
                      "model - a reference baseline for verification / overlay plots, not for "
                      "real operating points."},
                     {CavitationModel::GUMBEL, "Gumbel (half Sommerfeld)",
                      "Clamps sub-cavity pressures to p_cav after the solve. Fast and robust "
                      "but not mass-conserving - cavitation extent is approximate. OK for "
                      "quick load-capacity estimates."},
                     {CavitationModel::JFO, "JFO (Elrod-Adams)",
                      "Mass-conserving Jakobsson-Floberg-Olsson cavitation via the Elrod-Adams "
                      "film-content algorithm: rupture AND reformation are physical. The "
                      "recommended production model for accurate extent and oil transport."}});
                break;
            case ParamGroup::Boundaries:
                row("South boundary type (z = 0)", "Pressure condition at the south face.");
                changed |= enum_combo<BCType>(
                    "##bc_s", config_.bc_z_south_type,
                    {{BCType::DIRICHLET, "Fixed pressure (Dirichlet)",
                      "Holds the boundary at the given absolute pressure - a submerged or "
                      "vented bearing end. The usual choice."},
                     {BCType::NEUMANN, "Zero gradient (Neumann)",
                      "No axial pressure gradient - a sealed/symmetry end with no side flow."},
                     {BCType::INLET_OUTLET, "Inlet/outlet",
                      "Fixed pressure where flow enters, zero-gradient where it leaves - an "
                      "open end that should not pull unphysical inflow."}});
                row("South thermal inflow", "Temperature carried by oil entering at z = 0.");
                changed |= enum_combo<ThermalInflowMode>(
                    "##bc_s_th", config_.bc_z_south_thermal,
                    {{ThermalInflowMode::ZERO_GRADIENT, "Zero gradient (recirculating oil)",
                      "Entering oil is the same recirculating film oil (zero-gradient "
                      "temperature). Use for submerged/open ends."},
                     {ThermalInflowMode::CONSTANT, "Constant (fixed supply T)",
                      "Entering oil carries a fixed supply temperature. Use when the end is "
                      "fed from an external reservoir."}});
                row("North boundary type (z = L)", "Pressure condition at the north face.");
                changed |= enum_combo<BCType>(
                    "##bc_n", config_.bc_z_north_type,
                    {{BCType::DIRICHLET, "Fixed pressure (Dirichlet)",
                      "Holds the boundary at the given absolute pressure - a submerged or "
                      "vented bearing end. The usual choice."},
                     {BCType::NEUMANN, "Zero gradient (Neumann)",
                      "No axial pressure gradient - a sealed/symmetry end with no side flow."},
                     {BCType::INLET_OUTLET, "Inlet/outlet",
                      "Fixed pressure where flow enters, zero-gradient where it leaves - an "
                      "open end that should not pull unphysical inflow."}});
                row("North thermal inflow", "Temperature carried by oil entering at z = L.");
                changed |= enum_combo<ThermalInflowMode>(
                    "##bc_n_th", config_.bc_z_north_thermal,
                    {{ThermalInflowMode::ZERO_GRADIENT, "Zero gradient (recirculating oil)",
                      "Entering oil is the same recirculating film oil (zero-gradient "
                      "temperature). Use for submerged/open ends."},
                     {ThermalInflowMode::CONSTANT, "Constant (fixed supply T)",
                      "Entering oil carries a fixed supply temperature. Use when the end is "
                      "fed from an external reservoir."}});
                break;
            case ParamGroup::Thermal:
                row("Temperature model", "Whether heat transport in the film is solved.");
                changed |= enum_combo<TemperatureModel>(
                    "##temperature_model", config_.temperature_model,
                    {{TemperatureModel::ISOTHERMAL, "Isothermal",
                      "Constant film temperature; viscosity stays at its reference value. "
                      "Fast; fine when thermal rise is small or unknown inputs dominate."},
                     {TemperatureModel::ENERGY_EQUATION, "Energy equation (THD)",
                      "Transports viscous heating through the film (film-averaged energy "
                      "equation). Use when temperature rise matters - high speed, low "
                      "clearance, viscosity-temperature coupling."}});
                break;
            case ParamGroup::FluidProperties:
                if (show_advanced_) {
                    row("Property model", "How density/viscosity react to pressure, temperature "
                                          "and dissolved gas.");
                    changed |= enum_combo<FluidPropertyModel>(
                        "##fluid_property_model", config_.fluid_property_model,
                        {{FluidPropertyModel::CONSTANT, "Constant",
                          "Fixed density and viscosity. The default for plain oil bearings."},
                         {FluidPropertyModel::SINGLE_PHASE, "Single phase (oil + dissolved gas)",
                          "Liquid density/viscosity vary with the dissolved-gas fraction (and "
                          "T, p), but gas stays in solution - no free phase. Use for "
                          "refrigerant/oil mixtures below the bubble point."},
                         {FluidPropertyModel::TWO_PHASE, "Two phase (free gas)",
                          "Adds finite-rate gas release/resorption and a free-gas phase in "
                          "cavitated zones. The full gaseous-cavitation model; needs JFO "
                          "cavitation and gas parameters."}});
                    row("Dissolved gas species", "Gas dissolved in the lubricant.");
                    changed |= enum_combo<DissolvedGasSpecies>(
                        "##dissolved_gas_species", config_.dissolved_gas_species,
                        {{DissolvedGasSpecies::PROPANE, "Propane (R-290)",
                          "Refrigerant propane dissolved in oil (compressor bearings)."},
                         {DissolvedGasSpecies::AIR, "Air",
                          "Air dissolved in oil (conventional aerated lubricant)."}});
                    row("Solubility model", "Saturation law c_sat(p, T).");
                    changed |= enum_combo<SolubilityModel>(
                        "##solubility_model", config_.solubility_model,
                        {{SolubilityModel::HENRY, "Henry's law",
                          "c_sat = H(T) p, linear in pressure with optional van't Hoff "
                          "temperature dependence. Good for dilute solutions."},
                         {SolubilityModel::TABLE, "Table",
                          "Measured c_sat(p,T) surface (1-D or 2-D PTSV data). Use supplier "
                          "data for the real saturation curve."}});
                    row("Density model", "Liquid-solution density law.");
                    changed |= enum_combo<DensityModel>(
                        "##density_model", config_.density_model,
                        {{DensityModel::CONSTANT, "Constant",
                          "Dissolved gas does not change liquid density (rho = rho_oil)."},
                         {DensityModel::MIXTURE, "Mixture",
                          "Mass/volume mixing of the two phase densities: oil (rho) and the "
                          "dissolved-gas liquid phase (dissolved_gas_liquid_density), both "
                          "explicit inputs."},
                         {DensityModel::TABLE, "Table",
                          "Measured rho(p,T) surface."}});
                    row("Liquid viscosity model", "Liquid-solution viscosity law.");
                    changed |= enum_combo<LiquidViscosityModel>(
                        "##liquid_viscosity_model", config_.liquid_viscosity_model,
                        {{LiquidViscosityModel::CONSTANT, "Constant",
                          "Dissolved gas does not change liquid viscosity (mu = mu_oil)."},
                         {LiquidViscosityModel::EMPIRICAL, "Empirical",
                          "mu_oil * Andrade(T) * exp(a_c c_d) * Barus(p): exponential thinning "
                          "with dissolved gas plus temperature and pressure terms. Set the T/p "
                          "coefficients to 0 for pure gas-thinning (log mixing)."},
                         {LiquidViscosityModel::TABLE, "Table",
                          "Measured kinematic nu(p,T), converted to dynamic per cell."}});
                    row("Mixture viscosity", "Effective viscosity once free gas appears.");
                    changed |= enum_combo<MixtureViscosityModel>(
                        "##mixture_viscosity_model", config_.mixture_viscosity_model,
                        {{MixtureViscosityModel::DUKLER, "Dukler (legacy)",
                          "Void-fraction-weighted linear blend alpha*mu_g + (1-alpha)*mu_l; "
                          "needs mu_gas. Legacy/deprecated - kept for comparison with prior "
                          "results."},
                         {MixtureViscosityModel::LINEAR, "Linear (Grando)",
                          "Mass-quality-weighted linear blend x*mu_g + (1-x)*mu_l (Grando "
                          "2006); needs mu_gas. Refrigerant-oil validated."},
                         {MixtureViscosityModel::MCADAMS, "McAdams",
                          "Homogeneous two-phase standard: reciprocal blend by mass quality "
                          "1/mu = x/mu_g + (1-x)/mu_l; needs mu_gas. Correct asymptotes; the "
                          "default."}});
                }
                break;
            case ParamGroup::Motion:
                row("Motion model", "Whether the bearing can move.");
                changed |= enum_combo<MotionModel>(
                    "##motion_model", config_.motion_model,
                    {{MotionModel::STATIC, "Static",
                      "Bearing fixed at the configured eccentricity/attitude. Use to evaluate "
                      "an operating point."},
                     {MotionModel::MOVING_BEARING, "Moving bearing",
                      "Integrates the bearing equation of motion under fluid + external "
                      "loads (mass/stiffness/damping). Use for equilibrium search or "
                      "dynamics."}});
                break;
            case ParamGroup::Output: {
                row("Terminal verbosity", "How much the solver prints per step.");
                changed |= enum_combo<OutputVerbosity>(
                    "##output_verbosity", config_.output_verbosity,
                    {{OutputVerbosity::SIMPLIFIED, "Simplified",
                      "One short line per step. Keeps long runs readable."},
                     {OutputVerbosity::DETAILED, "Detailed",
                      "Adds per-step film thickness, temperature range, forces and friction "
                      "torque to the log."}});
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
                    changed |= enum_combo<TimeSteppingMethod>(
                        "##motion_time_method", config_.motion_time_method,
                        {{TimeSteppingMethod::EULER_EXPLICIT, "Explicit Euler (reference)",
                          "1st order, cheapest, but stability-limited - needs small dt. "
                          "Reference/comparison only."},
                         {TimeSteppingMethod::EULER_IMPLICIT, "Implicit Euler",
                          "1st order, unconditionally stable. The robust default."},
                         {TimeSteppingMethod::CRANK_NICOLSON, "Crank-Nicolson (reference)",
                          "2nd order, semi-implicit. More accurate trajectories but can ring "
                          "on stiff supports. Reference/comparison only."},
                         {TimeSteppingMethod::RK2, "Runge-Kutta 2 (reference)",
                          "2nd order explicit midpoint; stability-limited. Reference only."},
                         {TimeSteppingMethod::RK4, "Runge-Kutta 4",
                          "4th order explicit; most accurate per step (4 force evaluations). "
                          "Use for accurate whirl/orbit trajectories."}});

                    // Film-content face interpolation: Upwind or Linear.
                    // (TVD limiters and the thermal/gas scheme keys remain
                    // solver options reachable by editing the config file;
                    // a value loaded from file is preserved and shown.)
                    constexpr const char* kUpwindDesc =
                        "1st-order, unconditionally bounded and most robust, but numerically "
                        "diffusive - smears sharp cavitation fronts. The safe default.";
                    constexpr const char* kLinearDesc =
                        "2nd-order central interpolation. Sharpest in smooth regions but "
                        "UNBOUNDED - can oscillate at the rupture front. A textbook reference "
                        "scheme; use for comparison or smooth fields.";
                    constexpr const char* kVanLeerDesc =
                        "2nd-order van Leer TVD limiter: sharp AND bounded (no oscillations). "
                        "The recommended high-accuracy scheme for cavitation fronts.";
                    constexpr const char* kTypeDesc =
                        "Vijayaraghavan & Keith (1989): central across full-film faces, upwind "
                        "across cavitated faces. The classic Elrod scheme; JFO film-content only.";
                    row("Film-content convection", "Scheme for the JFO film-content transport.");
                    changed |= enum_combo<ConvectionScheme>(
                        "##theta_scheme", config_.theta_convection_scheme,
                        {{ConvectionScheme::UPWIND, "Upwind", kUpwindDesc},
                         {ConvectionScheme::LINEAR, "Linear (central)", kLinearDesc},
                         {ConvectionScheme::VANLEER, "TVD van Leer", kVanLeerDesc},
                         {ConvectionScheme::TYPE_DIFFERENCING, "Type differencing", kTypeDesc}});
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
    // Label-left input with a display-unit selector (shared per choice_key, so
    // all inlets switch units together).
    const auto unit_input = [&](const char* label, const char* choice_key, UnitFamily family,
                                double& stored, const char* fixed_unit) {
        const float label_column = ImGui::GetFontSize() * 11.0f;
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(label);
        ImGui::SameLine(label_column);
        const auto& options = unit_options(family);
        int& choice = unit_choice(choice_key, family);
        const double factor = options[static_cast<size_t>(choice)].to_display;
        double display = stored * factor;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 7.5f);
        if (ImGui::InputDouble((std::string("##v_") + label).c_str(), &display, 0, 0, "%g")) {
            stored = display / factor;
            dirty_ = true;
        }
        ImGui::SameLine();
        if (options.size() == 1) {
            ImGui::TextDisabled("%s", fixed_unit);
        } else {
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 4.4f);
            if (ImGui::BeginCombo((std::string("##u_") + label).c_str(),
                                  options[static_cast<size_t>(choice)].label)) {
                for (int option = 0; option < static_cast<int>(options.size()); ++option) {
                    if (ImGui::Selectable(options[static_cast<size_t>(option)].label,
                                          option == choice)) {
                        choice = option;
                    }
                }
                ImGui::EndCombo();
            }
        }
    };

    int remove_index = -1;
    for (size_t i = 0; i < config_.inlets.size(); ++i) {
        InletConfig& inlet = config_.inlets[i];
        ImGui::PushID(static_cast<int>(i));
        const bool is_groove = inlet.type == InletConfig::Type::GROOVE;
        std::string title = std::string(is_groove ? "Groove" : "Circular hole") + " inlet " +
                            std::to_string(i + 1);
        if (ImGui::TreeNodeEx(title.c_str(), ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted("Type");
            ImGui::SameLine(ImGui::GetFontSize() * 11.0f);
            int type_index = is_groove ? 1 : 0;
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 11);
            if (ImGui::Combo("##type", &type_index, "Circular hole\0Axial groove\0")) {
                inlet.type = type_index == 1 ? InletConfig::Type::GROOVE
                                             : InletConfig::Type::CIRCULAR;
                dirty_ = true;
            }
            unit_input("Position theta", "inlet.theta", UnitFamily::Angle, inlet.theta, "deg");
            if (inlet.type == InletConfig::Type::CIRCULAR) {
                unit_input("Position z", "inlet.z", UnitFamily::Length, inlet.z, "m");
                unit_input("Radius", "inlet.size_r", UnitFamily::Length, inlet.size, "m");
            } else {
                unit_input("Angular width", "inlet.size_w", UnitFamily::Angle, inlet.size, "deg");
            }
            unit_input("Supply pressure", "inlet.p_supply", UnitFamily::Pressure, inlet.p_supply,
                       "Pa");
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted("Feeds fresh oil");
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                ImGui::SetTooltip("Checked: an actual fed inlet pumping fresh oil in at the "
                                  "supply temperature.\nUnchecked: an open pocket of the same "
                                  "recirculating oil (pressure only).");
            }
            ImGui::SameLine(ImGui::GetFontSize() * 11.0f);
            if (ImGui::Checkbox("##feeds", &inlet.feeds_fresh_oil)) dirty_ = true;
            if (inlet.feeds_fresh_oil) {
                unit_input("Supply temperature", "inlet.t_supply", UnitFamily::Fixed,
                           inlet.t_supply, "K");
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
