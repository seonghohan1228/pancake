// Results window: run picker, key-results card with the design checks, PNG,
// CSV and case-report export, the step selector and the tab bar. The tabs
// (Field maps | Polar | Profiles | Transient | Comparison) are drawn by
// results_tabs.cpp. Grids are parsed off the UI thread (ResultsLoader, and
// std::async for motion.csv and the two comparison runs).

#include "gui_app.hpp"

#include "exports.hpp"
#include "imgui.h"
#include "implot.h"
#include "result_data.hpp"
#include "win_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace fs = std::filesystem;
using result_data::Profile;
using result_data::make_profile;

namespace {

constexpr const char* kGasOnLabel = "gas on (GC+VC)";
constexpr const char* kGasOffLabel = "gas off (VC only)";

struct ResultsUi {
    float profile_z = 0.5f;
    result_data::PlotCapture capture;
    fs::path pending_scan;
    result_data::ReloadGuard reload;
    const ResultGrid* results_grid = nullptr;
    int results_record = -2;
    gui::KeyResults grid_results;
    int estimate_record = -2;
    bool estimate_in_range = true;
    result_data::FieldMapUi field_map;
    result_data::TransientUi transient;
    result_data::CompareRun compare[2];
};

ResultsUi& state() {
    static ResultsUi s;
    return s;
}

std::string value_text(double si, gui::UnitFamily family, ui::UnitChoices& units, const char* id,
                       const char* format = "%.4g") {
    return std::isfinite(si) ? ui::format_value(si, family, units, id, format) : std::string("-");
}

std::string plain_text(double value, const char* format) {
    if (!std::isfinite(value)) return "-";
    char text[64];
    std::snprintf(text, sizeof(text), format, value);
    return text;
}

std::string record_label(const RunRecord& r, int index) {
    char text[320];
    std::snprintf(text, sizeof(text), "#%d  %s  [%s]  %s  %s", index + 1, r.job.label.c_str(),
                  r.job.group.c_str(), r.started_at.c_str(), to_string(r.outcome));
    return text;
}

gui::Light risk_light(gui::Risk risk) {
    switch (risk) {
        case gui::Risk::Low: return gui::Light::Green;
        case gui::Risk::Medium: return gui::Light::Yellow;
        case gui::Risk::High: return gui::Light::Red;
        case gui::Risk::None: break;
    }
    return gui::Light::Idle;
}

/// Coloured status mark with one sentence.
void mark(gui::Light light, const std::string& text) {
    ui::status_dot(light);
    ImGui::SameLine();
    ImGui::PushStyleColor(ImGuiCol_Text, ui::light_color(light));
    ImGui::TextWrapped("%s", text.c_str());
    ImGui::PopStyleColor();
}

/// SameLine when an item of `width` still fits on the line, else a new line.
void same_line_if_fits(float width) {
    ImGui::SameLine();
    if (ImGui::GetContentRegionAvail().x < width) ImGui::NewLine();
}

/// Width of a button with this label.
float button_width(const char* label) {
    return ImGui::CalcTextSize(label).x + ImGui::GetStyle().FramePadding.x * 2.0f;
}

void unit_choice_row(const char* label, const char* id, gui::UnitFamily family, ui::UnitChoices& choices) {
    const auto& options = gui::unit_options(family);
    result_data::unit_of(choices, id, family);
    int& choice = choices[id];
    ui::form_label(label, nullptr);
    ImGui::PushID(id);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##unit", options[static_cast<std::size_t>(choice)].label)) {
        for (int i = 0; i < static_cast<int>(options.size()); ++i)
            if (ImGui::Selectable(options[static_cast<std::size_t>(i)].label, i == choice)) choice = i;
        ImGui::EndCombo();
    }
    ImGui::PopID();
}

}  // namespace

void GuiApp::draw_results_window() {
    ResultsUi& s = state();
    s.reload.poll(loader_.busy(), loader_.revision(), steps_, selected_step_, grid_.get(),
                  [&](const StepInfo& info) { loader_.load_step(info); });
    if (!s.pending_scan.empty() && !loader_.busy()) {
        s.reload.arm(loader_.revision(), grid_.get());
        request_scan(s.pending_scan);
        s.pending_scan.clear();
    }
    if (!ImGui::Begin("Results")) {
        ImGui::End();
        return;
    }
    ui::UnitChoices& unit_map = units();
    const int record_count = static_cast<int>(records_.size());
    if (shown_record_ >= record_count) shown_record_ = -1;

    // --- Run picker, angle zero, units ---------------------------------------
    if (record_count > 0) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted("Run");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x * 0.5f, ImGui::GetFontSize() * 30.0f));
        const std::string preview = shown_record_ >= 0
            ? record_label(records_[static_cast<std::size_t>(shown_record_)], shown_record_)
            : std::string("Choose a run");
        if (ImGui::BeginCombo("##run", preview.c_str())) {
            for (int i = record_count - 1; i >= 0; --i) {
                ImGui::PushID(i);
                const RunRecord& r = records_[static_cast<std::size_t>(i)];
                if (ImGui::Selectable(record_label(r, i).c_str(), i == shown_record_) && i != shown_record_) {
                    shown_record_ = i;
                    selected_step_ = -1;
                    s.transient.playing = false;
                    s.pending_scan = r.job.output_dir;
                    refresh_derived();
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        ui::item_tooltip("Runs of this session, newest first: label, group, start time and outcome.");
        same_line_if_fits(button_width("Open results folder..."));
    }
    if (ImGui::Button("Open results folder...")) {
        const fs::path pvd = exports::open_file_dialog(
            hwnd, L"Result collection (results.pvd)\0*.pvd\0All files\0*.*\0");
        if (!pvd.empty()) {
            shown_record_ = -1;
            selected_step_ = -1;
            s.pending_scan = pvd.parent_path();
        }
    }
    ui::item_tooltip("Show the results of an earlier run: pick the results.pvd file in its results folder.");
    same_line_if_fits(ImGui::CalcTextSize("Angle zero").x + ImGui::GetStyle().ItemSpacing.x + ImGui::GetFontSize() * 12.0f);
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted("Angle zero");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
    if (ImGui::BeginCombo("##angle_zero", gui::angle_zero_name(angle_zero_))) {
        for (int z = 0; z <= static_cast<int>(gui::AngleZero::OppositeFeed); ++z) {
            const auto zero = static_cast<gui::AngleZero>(z);
            if (ImGui::Selectable(gui::angle_zero_name(zero), zero == angle_zero_)) {
                angle_zero_ = zero;
                settings.angle_zero = z;
            }
        }
        ImGui::EndCombo();
    }
    ui::item_tooltip("Zero of every angle axis and of the polar plot. The maximum and minimum film zeros "
                     "follow the field shown.");
    same_line_if_fits(button_width("Units"));
    if (ImGui::Button("Units")) ImGui::OpenPopup("result_units");
    ui::item_tooltip("Display units of the results.");
    if (ImGui::BeginPopup("result_units")) {
        if (ui::begin_form("units_form")) {
            unit_choice_row("Pressure", result_data::kUnitPressure, gui::UnitFamily::Pressure, unit_map);
            unit_choice_row("Film thickness", result_data::kUnitFilm, gui::UnitFamily::SmallLength, unit_map);
            unit_choice_row("Friction torque", result_data::kUnitTorque, gui::UnitFamily::Torque, unit_map);
            unit_choice_row("Friction power loss", result_data::kUnitPower, gui::UnitFamily::Power, unit_map);
            unit_choice_row("Temperature", result_data::kUnitTemperature, gui::UnitFamily::Temperature, unit_map);
            ui::end_form();
        }
        ImGui::EndPopup();
    }

    if (!grid_ && steps_.empty() && record_count == 0) {
        ImGui::Spacing();
        ImGui::TextDisabled("%s", results_error_.empty()
                                      ? "No results yet. Run the case (F5) or open a results folder."
                                      : results_error_.c_str());
        ImGui::End();
        return;
    }

    // --- Context of the shown run ---------------------------------------------
    const RunRecord* record = shown_record_ >= 0 ? &records_[static_cast<std::size_t>(shown_record_)] : nullptr;
    const SimulationConfig& cfg = record ? record->job.config : config_;
    const gui::CaseModel& model = record ? record->job.model : case_;
    if (grid_ && (s.results_grid != grid_.get() || s.results_record != shown_record_)) {
        s.results_grid = grid_.get();
        s.results_record = shown_record_;
        s.grid_results = gui::compute_key_results(*grid_, cfg);
    }
    if (s.estimate_record != shown_record_) {
        s.estimate_record = shown_record_;
        s.estimate_in_range = gui::compute_estimates(model, pair()).fluid_in_range;
    }
    const int last_step = static_cast<int>(steps_.size()) - 1;
    const int shown_step = selected_step_ < 0 ? last_step : std::min(selected_step_, last_step);
    const bool latest = shown_step == last_step;
    gui::KeyResults shown;
    bool have_results = false;
    if (record && record->results_ready && latest) {
        shown = record->results;
        have_results = true;
    } else if (grid_) {
        shown = s.grid_results;
        have_results = true;
    }
    const gui::KeyResults& grid_r = grid_ ? s.grid_results : shown;
    const double zero = result_data::zero_position(angle_zero_, model.feed.angle_deg, grid_r);
    const double start = gui::angle_axis_start(angle_zero_);
    const std::string x_label = gui::angle_axis_label(angle_zero_);
    const double p_sat = record ? record->job.saturation_pressure : estimates_.saturation_pressure;

    // --- Key-results card -------------------------------------------------------
    ui::section("Key results");
    if (record && record->outcome == SolverState::Failed)
        mark(gui::Light::Red, "The run failed (exit " + std::to_string(record->exit_code) +
                                  "): the solver log names the cause.");
    else if (record && record->outcome == SolverState::Cancelled)
        mark(gui::Light::Yellow, "The run was cancelled: the last saved step is shown.");
    else if (record && !record->converged)
        mark(gui::Light::Red, "The run did not converge: the results are not reliable.");
    if (!s.estimate_in_range)
        mark(gui::Light::Yellow, "The case state is outside the fluid table: properties were taken at the table edge.");
    if (record && !record->results_ready && !record->results_error.empty())
        ImGui::TextDisabled("Key results not read: %s", record->results_error.c_str());
    if (!have_results) {
        ImGui::TextDisabled(loader_.busy() ? "Reading the results..." : "No result step loaded.");
    } else {
        if (!steps_.empty())
            ImGui::TextDisabled("Output %d of %d, t = %.6g s%s", shown_step + 1, last_step + 1,
                                shown.time, record ? (record->converged ? ", converged" : ", not converged") : "");
        if (ImGui::BeginTable("key_results", 4, ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_RowBg |
                                                     ImGuiTableFlags_SizingStretchProp)) {
            const auto cell = [](const char* name, const char* help, const std::string& value) {
                ImGui::TableNextColumn();
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextUnformatted(name);
                ImGui::PopTextWrapPos();
                ui::item_tooltip(help);
                ImGui::TableNextColumn();
                ImGui::PushTextWrapPos(0.0f);
                ImGui::TextUnformatted(value.c_str());
                ImGui::PopTextWrapPos();
            };
            cell("Minimum film thickness h_min", "Smallest gap in the film.",
                 value_text(shown.h_min, gui::UnitFamily::SmallLength, unit_map, result_data::kUnitFilm, "%.3g"));
            cell("Eccentricity ratio", "Shaft centre offset per radial clearance, e / c.",
                 plain_text(shown.eccentricity_ratio, "%.3f"));
            cell("Attitude angle", "Angle from the load direction to the line of centres.",
                 plain_text(shown.attitude_deg, "%.1f deg"));
            cell("Friction torque", "Viscous torque on the shaft.",
                 value_text(shown.friction_torque, gui::UnitFamily::Torque, unit_map, result_data::kUnitTorque));
            cell("Friction power loss", "Friction torque times angular speed.",
                 value_text(shown.friction_power, gui::UnitFamily::Power, unit_map, result_data::kUnitPower));
            cell("Peak pressure", "Highest film pressure (absolute).",
                 value_text(shown.p_max, gui::UnitFamily::Pressure, unit_map, result_data::kUnitPressure));
            cell("Cavitated area (void fraction > 0)",
                 "Share of the film area where the gap is not completely filled by liquid (film content "
                 "below 1).",
                 plain_text(shown.cavitated_fraction * 100.0, "%.1f %%"));
            if (std::isfinite(shown.lambda))
                cell("Film thickness ratio Lambda", "Minimum film per combined surface roughness.",
                     plain_text(shown.lambda, "%.2f"));
            else
                cell("Lowest film pressure", "Lowest film pressure (absolute).",
                     value_text(shown.p_min, gui::UnitFamily::Pressure, unit_map, result_data::kUnitPressure));
            ImGui::EndTable();
        }
        ImGui::TextWrapped("Shear share of the film force: %s", plain_text(shown.shear_share * 100.0, "%.1f %%").c_str());
        ui::item_tooltip("Viscous (shear) part of the total film force; the rest is pressure.");
        const gui::DesignChecks checks = gui::design_checks(shown, p_sat, model.dissolved_gas);
        ImGui::PushTextWrapPos(0.0f);
        if (!checks.gas_release_text.empty()) mark(risk_light(checks.gas_release_risk), checks.gas_release_text);
        else ImGui::TextDisabled("Gas release check: the saturation pressure is unknown (choose a fluid pair).");
        if (!checks.lambda_text.empty()) mark(risk_light(checks.lambda_risk), checks.lambda_text);
        else ImGui::TextDisabled("Film thickness ratio: enter the surface roughness to check it.");
        ImGui::PopTextWrapPos();
    }
    if (record) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Command");
        ImGui::SameLine();
        if (fonts.mono) ImGui::PushFont(fonts.mono);
        ImGui::TextWrapped("%s", record->command.c_str());
        if (fonts.mono) ImGui::PopFont();
        if (ImGui::SmallButton("Copy command")) {
            ImGui::SetClipboardText(record->command.c_str());
            set_status("Command copied.");
        }
        ui::item_tooltip("Runs the same case from a terminal.");
        same_line_if_fits(button_width("PNG"));
    }

    // --- Export ----------------------------------------------------------------
    if (ImGui::SmallButton("PNG")) {
        const result_data::PlotCapture& c = s.capture;
        if (c.max_x > c.min_x) {
            const fs::path file = exports::save_file_dialog(hwnd, L"PNG image (*.png)\0*.png\0", L"png",
                                                            winutil::utf8_to_wide(c.name + ".png"));
            if (!file.empty()) pending_capture = CaptureRequest{c.min_x, c.min_y, c.max_x, c.max_y, file};
        } else {
            set_status("Show a plot first.");
        }
    }
    ui::item_tooltip("Save the plot shown in the tab below as a PNG image.");
    same_line_if_fits(button_width("CSV: key results"));
    if (ImGui::SmallButton("CSV: key results") && have_results) {
        const fs::path file = exports::save_file_dialog(hwnd, L"CSV file (*.csv)\0*.csv\0", L"csv", L"key_results.csv");
        if (!file.empty()) {
            const auto row = [](const char* name, double value, const char* unit) {
                return std::vector<std::string>{name, plain_text(value, "%.8g"), unit};
            };
            const std::vector<std::vector<std::string>> rows = {
                row("minimum film thickness", shown.h_min, "m"),
                row("eccentricity ratio", shown.eccentricity_ratio, "-"),
                row("attitude angle", shown.attitude_deg, "deg"),
                row("friction torque", shown.friction_torque, "N m"),
                row("friction power loss", shown.friction_power, "W"),
                row("peak pressure", shown.p_max, "Pa"),
                row("lowest film pressure", shown.p_min, "Pa"),
                row("cavitated area", shown.cavitated_fraction, "-"),
                row("shear share of film force", shown.shear_share, "-"),
                row("film thickness ratio", shown.lambda, "-"),
                row("maximum temperature", shown.t_max, "K"),
                row("time", shown.time, "s")};
            std::string error;
            if (exports::write_csv_rows(file, {"quantity", "value", "unit"}, rows, error))
                set_status("Saved " + winutil::path_to_utf8(file));
            else set_status("CSV export failed: " + error);
        }
    }
    ui::item_tooltip("Key results of the step shown, in SI units.");
    same_line_if_fits(button_width("CSV: profiles"));
    if (ImGui::SmallButton("CSV: profiles") && grid_) {
        const fs::path file = exports::save_file_dialog(hwnd, L"CSV file (*.csv)\0*.csv\0", L"csv", L"profiles.csv");
        if (!file.empty()) {
            const Profile profile = make_profile(*grid_, cfg, zero, start, s.profile_z);
            std::string error;
            if (exports::write_csv_columns(file, {"angle_deg (" + x_label + ")", "p_Pa", "h_m",
                                                  result_data::cavity_column(cfg), "T_K"},
                                           {profile.angle, profile.p, profile.h, profile.cavity, profile.t}, error))
                set_status("Saved " + winutil::path_to_utf8(file));
            else set_status("CSV export failed: " + error);
        }
    }
    ui::item_tooltip("Pressure, film thickness, void fraction (free-gas fraction for gas transport) and "
                     "temperature along the angle at the axial position of the Profiles tab, in SI units.");
    same_line_if_fits(button_width("Case report..."));
    if (ImGui::SmallButton("Case report...")) {
        const std::wstring name = (case_path_.empty() ? std::wstring(L"case") : case_path_.stem().wstring()) + L"_report.html";
        const fs::path file = exports::save_file_dialog(hwnd, L"HTML page (*.html)\0*.html\0", L"html", name);
        if (!file.empty()) export_case_report(file);
    }
    ui::item_tooltip("One printable page: inputs, estimates, key results, checks and figures.");

    // --- Step selector -------------------------------------------------------------
    if (steps_.size() > 1) {
        int index = shown_step;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 16.0f);
        if (ImGui::SliderInt("Output step", &index, 0, last_step, "%d")) {
            selected_step_ = index == last_step ? -1 : index;
            s.transient.playing = false;
        }
        ui::item_tooltip("Saved output shown in the field maps, polar plot and profiles.");
        char time_text[48];
        std::snprintf(time_text, sizeof(time_text), "t = %.6g s",
                      steps_[static_cast<std::size_t>(std::max(shown_step, 0))].time);
        same_line_if_fits(ImGui::CalcTextSize(time_text).x);
        ImGui::TextDisabled("%s", time_text);
    }

    // --- Tabs ------------------------------------------------------------------------
    const bool transient = model.analysis == gui::AnalysisType::Transient ||
                           cfg.solution_mode == SolutionMode::TRANSIENT;
    int on = -1, off = -1;
    for (int i = record_count - 1; i >= 0; --i) {
        const RunRecord& r = records_[static_cast<std::size_t>(i)];
        if (r.job.group != "comparison") continue;
        if (on < 0 && r.job.label == kGasOnLabel) on = i;
        if (off < 0 && r.job.label == kGasOffLabel) off = i;
    }
    result_data::TabContext tab;
    tab.grid = grid_.get();
    tab.config = &cfg;
    tab.model = &model;
    tab.grid_results = grid_r;
    tab.zero = zero;
    tab.start = start;
    tab.x_label = x_label;
    tab.p_sat = p_sat;
    tab.units = &unit_map;
    tab.loading = loader_.busy();
    if (ImGui::BeginTabBar("results_tabs")) {
        if (ImGui::BeginTabItem("Field maps")) {
            result_data::draw_field_maps_tab(s.field_map, tab, s.capture);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Polar")) {
            result_data::draw_polar_tab(tab, s.capture);
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Profiles")) {
            result_data::draw_profiles_tab(tab, s.profile_z, s.capture);
            ImGui::EndTabItem();
        }
        const fs::path base = record ? record->job.output_dir
                                     : (results_dir_.filename() == "flat" ? results_dir_.parent_path() : results_dir_);
        if (transient && ImGui::BeginTabItem("Transient")) {
            const int step = result_data::draw_transient_tab(s.transient, tab, base / "motion.csv", steps_, shown_step,
                                                             loader_.busy(), s.capture);
            if (step != -2) selected_step_ = step;
            ImGui::EndTabItem();
        }
        if ((on >= 0 || off >= 0) && ImGui::BeginTabItem("Comparison")) {
            result_data::ComparisonSide side_on, side_off;
            if (on >= 0) side_on = {&records_[static_cast<std::size_t>(on)], &s.compare[0]};
            if (off >= 0) side_off = {&records_[static_cast<std::size_t>(off)], &s.compare[1]};
            result_data::draw_comparison_tab(side_on, side_off, tab, angle_zero_, s.capture);
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}
