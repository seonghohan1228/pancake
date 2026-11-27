// Sweep / batch window: a sweep of one input around the current case (or a list
// of existing case files), case files and run scripts written to disk, the batch
// run through the GUI queue, a live table, result curves that grow as cases
// finish, and the onset chart with every case point.

#include "gui_app.hpp"

#include "exports.hpp"
#include "imgui.h"
#include "implot.h"
#include "onset_plot.hpp"
#include "result_data.hpp"
#include "win_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>

namespace fs = std::filesystem;

namespace {

enum class SweepInput { Speed, Load, Content, AmbientPressure, AmbientTemperature };

struct InputInfo {
    const char* name;
    const char* slug;
    gui::UnitFamily family;
    const char* unit_id;
    const char* help;
};

const InputInfo& info(SweepInput input) {
    static const InputInfo table[] = {
        {"Speed", "speed", gui::UnitFamily::Speed, "batch.speed", "Shaft rotational speed."},
        {"Load", "load", gui::UnitFamily::Force, "batch.load", "Load magnitude; the direction stays as in the case."},
        {"Refrigerant content", "content", gui::UnitFamily::Content, "batch.content",
         "Mass fraction of refrigerant dissolved in the oil (the case is no longer 'saturated')."},
        {"Ambient pressure", "ambient_pressure", gui::UnitFamily::Pressure, "batch.ambient_pressure",
         "Pressure around the bearing."},
        {"Ambient temperature", "ambient_temperature", gui::UnitFamily::Temperature, "batch.ambient_temperature",
         "Temperature of the surroundings and the supplied oil."},
    };
    return table[static_cast<int>(input)];
}

double input_value(const gui::CaseModel& m, SweepInput input, const gui::FluidPair* pair) {
    switch (input) {
        case SweepInput::Speed: return m.speed;
        case SweepInput::Load: return m.load;
        case SweepInput::Content: return gui::effective_content(m, pair);
        case SweepInput::AmbientPressure: return m.ambient_pressure;
        case SweepInput::AmbientTemperature: return m.ambient_temperature;
    }
    return 0.0;
}

void set_input(gui::CaseModel& m, SweepInput input, double value) {
    switch (input) {
        case SweepInput::Speed: m.speed = value; break;
        case SweepInput::Load: m.load = value; break;
        case SweepInput::Content:
            m.saturated = false;
            m.content = value;
            break;
        case SweepInput::AmbientPressure: m.ambient_pressure = value; break;
        case SweepInput::AmbientTemperature: m.ambient_temperature = value; break;
    }
}

struct BatchCase {
    gui::CaseModel model;
    SimulationConfig config;
    double value = 0.0;           // SI (file cases: index)
    std::string label;            // value with unit, or the file name
    fs::path config_file, output_dir;
    gui::Estimates estimates;
};

struct BatchUi {
    SweepInput input = SweepInput::Load;
    int seeded_input = -1;
    bool use_list = false;
    bool files_mode = false;
    double from = 0.0, to = 0.0;
    int steps = 5;
    char list[512] = "";
    int cores = 1;
    std::vector<fs::path> files;
    std::string batch_id;          // group of the jobs of the last "Run batch"
    std::vector<BatchCase> cases;  // of the last written or run batch
    int y_quantity = 0;
};

BatchUi& state() {
    static BatchUi s;
    return s;
}

std::string local_time() {
    SYSTEMTIME now{};
    GetLocalTime(&now);
    char text[16];
    std::snprintf(text, sizeof(text), "%02u:%02u:%02u", now.wHour, now.wMinute, now.wSecond);
    return text;
}

std::string text_or_dash(double value, const char* format) {
    if (!std::isfinite(value)) return "-";
    char text[64];
    std::snprintf(text, sizeof(text), format, value);
    return text;
}

ImVec4 verdict_colour(gui::OnsetVerdict verdict) {
    switch (verdict) {
        case gui::OnsetVerdict::Expected: return theme::kRiskHigh;
        case gui::OnsetVerdict::Borderline: return theme::kRiskMedium;
        case gui::OnsetVerdict::NotExpected: return theme::kRiskLow;
        case gui::OnsetVerdict::Unknown: break;
    }
    return theme::kMuted;
}

/// Swept values in SI from the range or the list (list entries in the display unit).
std::vector<double> sweep_values(const BatchUi& s, const gui::UnitOption& unit) {
    std::vector<double> values;
    if (s.use_list) {
        std::string text = s.list;
        for (char& ch : text)
            if (ch == ';' || ch == '\t') ch = ',';
        std::stringstream cells(text);
        std::string cell;
        while (std::getline(cells, cell, ',')) {
            char* end = nullptr;
            const double v = std::strtod(cell.c_str(), &end);
            if (end != cell.c_str() && std::isfinite(v)) values.push_back(gui::to_si(unit, v));
        }
        return values;
    }
    const int n = std::clamp(s.steps, 1, 200);
    for (int i = 0; i < n; ++i)
        values.push_back(n == 1 ? s.from : s.from + (s.to - s.from) * i / (n - 1));
    return values;
}

bool write_scripts(const fs::path& folder, const std::vector<std::string>& commands, std::string& error) {
    std::ofstream bat(folder / "run_all.bat", std::ios::binary | std::ios::trunc);
    std::ofstream ps1(folder / "run_all.ps1", std::ios::binary | std::ios::trunc);
    if (!bat || !ps1) {
        error = "cannot write the run scripts in " + winutil::path_to_utf8(folder);
        return false;
    }
    bat << "@echo off\r\nrem Batch written by the pancake GUI: the cases run one after another.\r\n"
        << "cd /d \"%~dp0\"\r\n";
    ps1 << "# Batch written by the pancake GUI: the cases run one after another.\r\n"
        << "Set-Location -Path $PSScriptRoot\r\n";
    for (std::size_t i = 0; i < commands.size(); ++i) {
        bat << "echo Case " << i + 1 << " of " << commands.size() << "\r\n" << commands[i] << "\r\n"
            << "if errorlevel 1 echo Case " << i + 1 << " failed.\r\n";
        ps1 << "Write-Host 'Case " << i + 1 << " of " << commands.size() << "'\r\n"
            << (commands[i].front() == '"' ? "& " : "") << commands[i] << "\r\n"
            << "if ($LASTEXITCODE -ne 0) { Write-Host 'Case " << i + 1 << " failed.' }\r\n";
    }
    return bat.good() && ps1.good();
}

struct Quantity {
    const char* name;
    double gui::KeyResults::*field;
    gui::UnitFamily family;
    const char* unit_id;  // null: plain number times factor
    double factor;
};

const Quantity kQuantities[] = {
    {"Minimum film thickness", &gui::KeyResults::h_min, gui::UnitFamily::SmallLength, result_data::kUnitFilm, 1.0},
    {"Eccentricity ratio", &gui::KeyResults::eccentricity_ratio, gui::UnitFamily::Fixed, nullptr, 1.0},
    {"Attitude angle [deg]", &gui::KeyResults::attitude_deg, gui::UnitFamily::Fixed, nullptr, 1.0},
    {"Friction torque", &gui::KeyResults::friction_torque, gui::UnitFamily::Torque, result_data::kUnitTorque, 1.0},
    {"Friction power loss", &gui::KeyResults::friction_power, gui::UnitFamily::Power, result_data::kUnitPower, 1.0},
    {"Cavitated area [%]", &gui::KeyResults::cavitated_fraction, gui::UnitFamily::Fixed, nullptr, 100.0},
    {"Peak pressure", &gui::KeyResults::p_max, gui::UnitFamily::Pressure, result_data::kUnitPressure, 1.0},
};

}  // namespace

void GuiApp::draw_batch_window() {
    BatchUi& s = state();
    ImGui::SetNextWindowSize(ImVec2(980, 720), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Sweep / batch", &show_batch_)) {
        ImGui::End();
        return;
    }
    ui::window_shadow();
    ui::UnitChoices& unit_map = units();
    const InputInfo& in = info(s.input);
    const gui::UnitOption& unit = result_data::unit_of(unit_map, in.unit_id, in.family);
    if (s.seeded_input != static_cast<int>(s.input)) {
        s.seeded_input = static_cast<int>(s.input);
        const double base = input_value(case_, s.input, pair());
        s.from = in.family == gui::UnitFamily::Temperature ? base - 10.0 : base * 0.8;
        s.to = in.family == gui::UnitFamily::Temperature ? base + 10.0 : base * 1.2;
    }

    // --- Definition ---------------------------------------------------------------
    ui::section("Sweep definition");
    if (ImGui::RadioButton("Sweep one input of the current case", !s.files_mode)) s.files_mode = false;
    ui::item_tooltip("Copies of the current case that differ in one input.");
    ImGui::SameLine();
    if (ImGui::GetContentRegionAvail().x < ImGui::GetFrameHeight() + ImGui::GetStyle().ItemInnerSpacing.x +
                                             ImGui::CalcTextSize("Run existing case files").x)
        ImGui::NewLine();
    if (ImGui::RadioButton("Run existing case files", s.files_mode)) s.files_mode = true;
    ui::item_tooltip("Run config files as they are, one after another.");
    if (!s.files_mode && ui::begin_form("sweep_form")) {
        ui::form_combo("Input", "Operating condition varied from case to case; everything else stays as in the "
                       "current case.", "##input", s.input,
                       {{SweepInput::Speed, "Speed", info(SweepInput::Speed).help},
                        {SweepInput::Load, "Load", info(SweepInput::Load).help},
                        {SweepInput::Content, "Refrigerant content", info(SweepInput::Content).help},
                        {SweepInput::AmbientPressure, "Ambient pressure", info(SweepInput::AmbientPressure).help},
                        {SweepInput::AmbientTemperature, "Ambient temperature", info(SweepInput::AmbientTemperature).help}});
        ui::form_label("Values", "A linear range, or a list of values in the unit shown.");
        ImGui::Checkbox("List", &s.use_list);
        ui::item_tooltip("Enter the values as a comma-separated list instead of a range.");
        if (!s.use_list) {
            ui::form_label("From", "First value of the range.");
            ImGui::PushID("from");
            ui::unit_input(in.unit_id, s.from, in.family, unit_map);
            ImGui::PopID();
            ui::form_label("To", "Last value of the range.");
            ImGui::PushID("to");
            ui::unit_input(in.unit_id, s.to, in.family, unit_map);
            ImGui::PopID();
            ui::form_int("Number of cases", "Cases spread evenly from the first to the last value (typical 5).",
                         "##steps", s.steps, 1, 200);
        } else {
            ui::form_label("List", "Comma-separated values, e.g. 2000, 2500, 3000.");
            ImGui::SetNextItemWidth(-ImGui::CalcTextSize(unit.label).x - ImGui::GetFontSize());
            ImGui::InputText("##list", s.list, sizeof(s.list));
            ImGui::SameLine();
            ImGui::TextDisabled("%s", unit.label);
        }
        ui::form_int("Cores per case", "CPU cores used by each solver run (cases run one after another).",
                     "##cores", s.cores, 1, 64);
        ui::end_form();
    }
    if (s.files_mode) {
        if (ImGui::Button("Add case files...")) {
            const fs::path file = exports::open_file_dialog(hwnd, L"Case files (*.txt)\0*.txt\0All files\0*.*\0");
            if (!file.empty() && std::find(s.files.begin(), s.files.end(), file) == s.files.end())
                s.files.push_back(file);
        }
        ui::item_tooltip("Add one config file to the batch; it runs exactly as written.");
        ImGui::SameLine();
        if (ImGui::Button("Clear list")) s.files.clear();
        ImGui::SameLine();
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6.0f);
        ImGui::InputInt("Cores per case", &s.cores);
        s.cores = std::clamp(s.cores, 1, 64);
        for (const fs::path& file : s.files) ImGui::BulletText("%s", winutil::path_to_utf8(file).c_str());
        if (s.files.empty()) ImGui::TextDisabled("No case files added.");
    }

    const std::vector<double> values = sweep_values(s, unit);
    if (!s.files_mode) {
        std::string preview;
        for (std::size_t i = 0; i < values.size() && i < 12; ++i)
            preview += (i ? ", " : "") + text_or_dash(gui::to_display(unit, values[i]), "%.4g");
        if (values.size() > 12) preview += ", ...";
        ImGui::TextDisabled("%d case(s): %s %s", static_cast<int>(values.size()), preview.c_str(), unit.label);
    }

    // Builds the case list (sweep variants written to disk, or the given files).
    const auto build_cases = [&](std::vector<BatchCase>& out, fs::path& folder, std::string& error) -> bool {
        out.clear();
        if (s.files_mode) {
            if (s.files.empty()) return error = "add case files first", false;
            folder = s.files.front().parent_path();
            for (std::size_t i = 0; i < s.files.size(); ++i) {
                BatchCase c;
                std::ifstream file(s.files[i], std::ios::binary);
                std::ostringstream text;
                text << file.rdbuf();
                c.config.load_from_text(text.str());
                c.config.config_dir = s.files[i].parent_path();
                c.model = gui::from_config(c.config, gui::read_gui_meta(text.str()), pair()).model;
                c.config_file = s.files[i];
                fs::path output = c.config.output_dir;
                c.output_dir = output.is_relative() ? s.files[i].parent_path() / output : output;
                c.value = static_cast<double>(i + 1);
                c.label = winutil::path_to_utf8(s.files[i].filename());
                c.estimates = gui::compute_estimates(c.model, pair());
                out.push_back(std::move(c));
            }
            return true;
        }
        if (case_path_.empty()) return error = "save the case first (File > Save)", false;
        if (values.empty()) return error = "enter at least one value", false;
        folder = case_path_.parent_path() / (case_path_.stem().wstring() + L"_sweep_" + winutil::utf8_to_wide(in.slug));
        std::error_code ec;
        fs::create_directories(folder, ec);
        if (ec) return error = "cannot create " + winutil::path_to_utf8(folder), false;
        for (std::size_t i = 0; i < values.size(); ++i) {
            BatchCase c;
            c.model = case_;
            set_input(c.model, s.input, values[i]);
            c.model.cores = s.cores;
            char name[32];
            std::snprintf(name, sizeof(name), "case_%03d.txt", static_cast<int>(i + 1));
            c.config_file = folder / name;
            c.output_dir = output_dir_for(c.config_file);
            if (!write_case_file(c.model, c.config_file, c.output_dir, error)) return false;
            c.config = gui::to_config(c.model, pair(), base_config_, loaded_case_);
            c.value = values[i];
            c.label = ui::format_value(values[i], in.family, unit_map, in.unit_id, "%.4g");
            c.estimates = gui::compute_estimates(c.model, pair());
            out.push_back(std::move(c));
        }
        return true;
    };
    const auto write_all = [&](bool announce) -> bool {
        fs::path folder;
        std::string error;
        std::vector<BatchCase> built;
        if (!build_cases(built, folder, error)) {
            set_status("Batch not written: " + error);
            return false;
        }
        std::vector<std::string> commands;
        for (const BatchCase& c : built) commands.push_back(command_line_for(c.config_file, s.cores));
        if (!write_scripts(folder, commands, error)) {
            set_status("Batch scripts not written: " + error);
            return false;
        }
        s.cases = std::move(built);
        s.batch_id.clear();
        if (announce) set_status("Wrote " + std::to_string(s.cases.size()) + " cases and run_all.ps1 / run_all.bat in " +
                                 winutil::path_to_utf8(folder));
        return true;
    };

    const bool can_sweep = s.files_mode || ready_to_run();
    ImGui::BeginDisabled(!can_sweep);
    if (ImGui::Button("Write cases and scripts")) write_all(true);
    ui::item_tooltip("Write one config file per case and run_all.ps1 / run_all.bat that run them with the exact "
                     "solver command, so the batch can also run from a terminal.");
    ImGui::SameLine();
    if (ui::primary_button("Run batch") && write_all(false)) {
        s.batch_id = "batch " + local_time();
        for (const BatchCase& c : s.cases) {
            RunJob job;
            job.label = c.label;
            job.group = s.batch_id;
            job.config_file = c.config_file;
            job.output_dir = c.output_dir;
            job.cores = s.cores;
            job.model = c.model;
            job.config = c.config;
            job.sweep_value = c.value;
            if (pair())
                job.saturation_pressure =
                    pair_.saturation_pressure(gui::effective_content(c.model, pair()), c.model.ambient_temperature).value;
            enqueue(std::move(job));
        }
        set_status("Queued " + std::to_string(s.cases.size()) + " cases (" + s.batch_id + ").");
    }
    ui::item_tooltip("Write the cases and run them one after another in this window's queue.");
    ImGui::EndDisabled();
    if (!can_sweep) {
        ImGui::SameLine();
        ImGui::TextColored(theme::kRiskMedium, "Complete the case setup first (every workflow light green).");
    }

    // --- Live table ----------------------------------------------------------------
    ui::section("Cases");
    if (s.cases.empty()) {
        ImGui::TextDisabled("No batch yet.");
        ImGui::End();
        return;
    }
    const auto find_record = [&](const BatchCase& c) -> const RunRecord* {
        if (s.batch_id.empty()) return nullptr;
        for (auto it = records_.rbegin(); it != records_.rend(); ++it)
            if (it->job.group == s.batch_id && it->job.config_file == c.config_file) return &*it;
        return nullptr;
    };
    const auto status_of = [&](const BatchCase& c, const RunRecord* r) -> std::string {
        if (s.batch_id.empty()) return "not started";
        if (active_job_ && active_job_->group == s.batch_id && active_job_->config_file == c.config_file) return "running";
        for (const RunJob& job : queue_)
            if (job.group == s.batch_id && job.config_file == c.config_file) return "queued";
        if (!r) return "not queued";
        if (r->outcome == SolverState::Finished) return "done";
        if (r->outcome == SolverState::Failed) return "failed";
        return "cancelled";
    };
    const char* headers[] = {"#", "Value", "Status", "h_min", "Eccentricity", "Torque", "Power", "Cavitated %",
                             "p_max", "Onset (estimate)", "Converged"};
    std::vector<std::vector<std::string>> csv_rows;
    const float table_h = std::min(ImGui::GetTextLineHeightWithSpacing() * (s.cases.size() + 1.5f),
                                   ImGui::GetContentRegionAvail().y * 0.45f);
    if (ImGui::BeginTable("batch_table", 11, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                                 ImGuiTableFlags_ScrollY | ImGuiTableFlags_ScrollX |
                                                 ImGuiTableFlags_SizingFixedFit,
                          ImVec2(0, std::max(table_h, 80.0f)))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        for (const char* header : headers) ImGui::TableSetupColumn(header);
        ImGui::TableHeadersRow();
        for (std::size_t i = 0; i < s.cases.size(); ++i) {
            const BatchCase& c = s.cases[i];
            const RunRecord* r = find_record(c);
            const bool ready = r && r->results_ready;
            const gui::KeyResults* k = ready ? &r->results : nullptr;
            const auto q = [&](double gui::KeyResults::*field, gui::UnitFamily family, const char* id, double factor,
                               const char* format) {
                if (!k || !std::isfinite(k->*field)) return std::string("-");
                return id ? ui::format_value(k->*field, family, unit_map, id, "%.4g") : text_or_dash(k->*field * factor, format);
            };
            const std::string status = status_of(c, r);
            const std::vector<std::string> row = {
                std::to_string(i + 1), c.label, status,
                q(&gui::KeyResults::h_min, gui::UnitFamily::SmallLength, result_data::kUnitFilm, 1, nullptr),
                q(&gui::KeyResults::eccentricity_ratio, gui::UnitFamily::Fixed, nullptr, 1, "%.3f"),
                q(&gui::KeyResults::friction_torque, gui::UnitFamily::Torque, result_data::kUnitTorque, 1, nullptr),
                q(&gui::KeyResults::friction_power, gui::UnitFamily::Power, result_data::kUnitPower, 1, nullptr),
                q(&gui::KeyResults::cavitated_fraction, gui::UnitFamily::Fixed, nullptr, 100, "%.1f"),
                q(&gui::KeyResults::p_max, gui::UnitFamily::Pressure, result_data::kUnitPressure, 1, nullptr),
                gui::onset_verdict_text(c.estimates.verdict),
                r && r->outcome != SolverState::Running ? (r->converged ? "yes" : "no") : "-"};
            csv_rows.push_back(row);
            ImGui::TableNextRow();
            for (int col = 0; col < 11; ++col) {
                ImGui::TableNextColumn();
                const std::string& cell = row[static_cast<std::size_t>(col)];
                if (col == 2 && status == "failed") ImGui::TextColored(theme::kRiskHigh, "%s", cell.c_str());
                else if (col == 2 && status == "running") ImGui::TextColored(theme::kAccent, "%s", cell.c_str());
                else if (col == 9) ImGui::TextColored(verdict_colour(c.estimates.verdict), "%s", cell.c_str());
                else if (col == 10 && cell == "no") ImGui::TextColored(theme::kRiskHigh, "%s", cell.c_str());
                else ImGui::TextUnformatted(cell.c_str());
                if (col == 9) ui::item_tooltip("Gaseous-cavitation onset from the onset chart (before the run).");
            }
        }
        ImGui::EndTable();
    }
    if (ImGui::SmallButton("Export table (CSV)")) {
        const fs::path file = exports::save_file_dialog(hwnd, L"CSV file (*.csv)\0*.csv\0", L"csv", L"batch.csv");
        std::string error;
        if (!file.empty()) {
            if (exports::write_csv_rows(file, std::vector<std::string>(std::begin(headers), std::end(headers)), csv_rows, error))
                set_status("Saved " + winutil::path_to_utf8(file));
            else set_status("CSV export failed: " + error);
        }
    }
    ui::item_tooltip("The table as shown, in the display units.");

    // --- Result curve and onset chart -------------------------------------------------
    const float plot_h = std::max(220.0f, ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing());
    if (ImGui::BeginTable("batch_plots", 2, ImGuiTableFlags_SizingStretchSame)) {
        ImGui::TableNextColumn();
        const char* preview = kQuantities[s.y_quantity].name;
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 14.0f);
        if (ImGui::BeginCombo("Result", preview)) {
            for (int i = 0; i < static_cast<int>(std::size(kQuantities)); ++i)
                if (ImGui::Selectable(kQuantities[i].name, i == s.y_quantity)) s.y_quantity = i;
            ImGui::EndCombo();
        }
        ui::item_tooltip("Result plotted against the swept input; points appear as cases finish.");
        const Quantity& quantity = kQuantities[s.y_quantity];
        std::vector<double> xs, ys;
        for (const BatchCase& c : s.cases) {
            const RunRecord* r = find_record(c);
            if (!r || !r->results_ready || !std::isfinite(r->results.*quantity.field)) continue;
            xs.push_back(s.files_mode ? c.value : gui::to_display(unit, c.value));
            const double v = r->results.*quantity.field;
            ys.push_back(quantity.unit_id
                             ? gui::to_display(result_data::unit_of(unit_map, quantity.unit_id, quantity.family), v)
                             : v * quantity.factor);
        }
        std::string y_label = quantity.name;
        if (quantity.unit_id)
            y_label += std::string(" [") + result_data::unit_of(unit_map, quantity.unit_id, quantity.family).label + "]";
        const std::string x_label = s.files_mode ? std::string("Case #") : std::string(in.name) + " [" + unit.label + "]";
        if (ImPlot::BeginPlot("##batch_curve", ImVec2(-1, plot_h - ImGui::GetFrameHeightWithSpacing()),
                              ImPlotFlags_NoLegend)) {
            const bool data = ui::has_points(xs, ys);
            ImPlot::SetupAxes(x_label.c_str(), y_label.c_str(), ui::fit_if(data), ui::fit_if(data));
            if (!data) ui::empty_plot("Results appear as cases finish");
            ImPlot::SetNextLineStyle(theme::kAccent, 2.0f);
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 4.0f);
            if (data) ImPlot::PlotLine("##curve", xs.data(), ys.data(), static_cast<int>(xs.size()));
            ImPlot::EndPlot();
        }
        ImGui::TableNextColumn();
        std::vector<ui::OnsetMarker> markers;
        for (const BatchCase& c : s.cases)
            markers.push_back({c.label, c.estimates.sommerfeld, c.estimates.margin_ratio, c.estimates.verdict});
        ImGui::TextDisabled("Onset chart (estimate before the run)");
        ui::onset_chart("##batch_onset", ImVec2(-1, plot_h - ImGui::GetFrameHeightWithSpacing()), markers);
        ImGui::EndTable();
    }
    ImGui::End();
}
