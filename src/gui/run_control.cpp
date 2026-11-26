// Run control: the run queue (single, comparison, batch, two-step thermal),
// solver launch and progress, key results parsed off the UI thread, results
// loading, and the scripted --run mode.

#include "gui_app.hpp"

#include "config_text_patch.hpp"
#include "imgui.h"
#include "win_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

constexpr double kPi = 3.14159265358979323846;

fs::path find_mpiexec() {
    wchar_t buffer[MAX_PATH] = {};
    const DWORD env_len = GetEnvironmentVariableW(L"MSMPI_BIN", buffer, MAX_PATH);
    if (env_len > 0 && env_len < MAX_PATH) {
        const fs::path candidate = fs::path(buffer) / L"mpiexec.exe";
        if (fs::exists(candidate)) return candidate;
    }
    for (const wchar_t* prefix : {L"C:/Program Files/Microsoft MPI/Bin/mpiexec.exe",
                                  L"C:/Program Files (x86)/Microsoft SDKs/MPI/Bin/mpiexec.exe"}) {
        if (fs::exists(prefix)) return prefix;
    }
    wchar_t found[MAX_PATH] = {};
    if (SearchPathW(nullptr, L"mpiexec.exe", nullptr, MAX_PATH, found, nullptr) > 0) return found;
    return {};
}

std::string local_time_string() {
    SYSTEMTIME now{};
    GetLocalTime(&now);
    char text[16];
    std::snprintf(text, sizeof(text), "%02u:%02u:%02u", now.wHour, now.wMinute, now.wSecond);
    return text;
}

std::string quoted(const fs::path& path) { return "\"" + winutil::path_to_utf8(path) + "\""; }

}  // namespace

void GuiApp::advance_scripted_run() {
    if (scripted_run_state_ == 1) {
        if (!ready_to_run()) {
            log_.push_back({"--run: the case is not ready to run.", IM_COL32(185, 35, 35, 255)});
            scripted_run_state_ = 3;
            wants_exit = startup_screenshots.empty();
            return;
        }
        start_run();
        scripted_run_state_ = 2;
        return;
    }
    if (runner_.running() || !queue_.empty() || active_job_ || !result_tasks_.empty()) return;
    for (const RunRecord& record : records_)
        write_key_results_json(record, record.job.output_dir / "gui_key_results.json");
    scripted_run_state_ = 3;
    if (!startup_screenshots.empty()) {
        screenshot_dir_ = startup_screenshots;
        exit_after_screenshots_ = true;
        start_screenshot_series();
    } else {
        wants_exit = true;
    }
}

bool GuiApp::write_key_results_json(const RunRecord& record, const fs::path& file) const {
    std::ofstream out(file, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    const gui::KeyResults& r = record.results;
    const char* nl = "\n";
    const auto number = [&](const char* key, double value, bool last = false) {
        out << "  \"" << key << "\": ";
        if (std::isfinite(value)) out << value;
        else out << "null";
        out << (last ? "" : ",") << nl;
    };
    out.precision(10);
    out << "{" << nl << "  \"label\": \"" << record.job.label << "\"," << nl
        << "  \"outcome\": \"" << to_string(record.outcome) << "\"," << nl
        << "  \"converged\": " << (record.converged ? "true" : "false") << "," << nl
        << "  \"results_ready\": " << (record.results_ready ? "true" : "false") << "," << nl
        << "  \"duration_s\": " << record.duration_s << "," << nl;
    number("h_min_m", r.h_min);
    number("eccentricity_ratio", r.eccentricity_ratio);
    number("attitude_deg", r.attitude_deg);
    number("friction_torque_Nm", r.friction_torque);
    number("friction_power_W", r.friction_power);
    number("p_max_Pa", r.p_max);
    number("p_min_Pa", r.p_min);
    number("cavitated_fraction", r.cavitated_fraction);
    number("shear_share", r.shear_share);
    number("lambda", r.lambda, true);
    out << "}" << nl;
    return static_cast<bool>(out);
}


// ===========================================================================
// Run control
// ===========================================================================

std::string GuiApp::command_line_for(const fs::path& config, int cores) const {
    std::string source;
    const fs::path solver = discover_solver(source);
    const std::string exe = solver.empty() ? std::string("pancake.exe") : quoted(solver);
    if (cores > 1) return "mpiexec -n " + std::to_string(cores) + " " + exe + " -c " + quoted(config);
    return exe + " -c " + quoted(config);
}

fs::path GuiApp::discover_solver(std::string& source) const {
    if (!settings.solver_path.empty()) {
        fs::path override_path = winutil::utf8_to_path(settings.solver_path);
        if (override_path.is_relative()) override_path = winutil::executable_dir() / override_path;
        if (fs::exists(override_path)) {
            source = "settings";
            return override_path;
        }
    }
    const fs::path candidate = winutil::executable_dir() / "pancake.exe";
    if (fs::exists(candidate)) {
        source = "next to the GUI";
        return candidate;
    }
    source = "not found";
    return {};
}

RunJob GuiApp::make_job(const gui::CaseModel& model, const std::string& label, const fs::path& file,
                        const std::string& group) const {
    RunJob job;
    job.label = label;
    job.group = group;
    job.config_file = file;
    job.output_dir = output_dir_for(file);
    job.cores = std::max(1, model.cores);
    job.model = model;
    job.config = gui::to_config(model, pair(), base_config_, loaded_case_);
    if (pair()) {
        job.saturation_pressure =
            pair_.saturation_pressure(gui::effective_content(model, pair()), model.ambient_temperature).value;
    }
    return job;
}

void GuiApp::start_run() {
    if (runner_.running() || !ready_to_run()) return;
    // A thermal transient runs for hours; ask first (plan 0.3).
    if (case_.thermal && case_.analysis == gui::AnalysisType::Transient && !thermal_transient_confirmed_) {
        show_thermal_transient_dialog_ = true;
        return;
    }
    thermal_transient_confirmed_ = false;
    // The solver runs the file on disk: save first so the file equals the GUI
    // state and writes its results where the GUI reads them (<case>_results).
    if (!save_case(case_path_.empty())) return;
    if (!comparison_mode_) {
        enqueue(make_job(case_, "case", case_path_, "single"));
        return;
    }
    // Comparison: two copies of the case, gas on and gas off, run one after another.
    for (bool gas : {true, false}) {
        gui::CaseModel model = case_;
        model.dissolved_gas = gas;
        const fs::path file = case_path_.parent_path() /
                              (case_path_.stem().wstring() + (gas ? L"_gas_on.txt" : L"_gas_off.txt"));
        std::string error;
        if (!write_case_file(model, file, output_dir_for(file), error)) {
            set_status("Cannot write the comparison case: " + error);
            return;
        }
        enqueue(make_job(model, gas ? "gas on (GC+VC)" : "gas off (VC only)", file, "comparison"));
    }
}

bool GuiApp::two_step_cases(gui::CaseModel& step1, gui::CaseModel& step2, std::string& reason) const {
    step1 = case_;
    step1.analysis = gui::AnalysisType::Steady;
    step1.position = gui::ShaftPosition::LoadBalance;
    if (!case_.load_table.file.empty()) {
        if (!load_table_ok_ || load_rows_.empty()) {
            reason = "the load table cannot be read";
            return false;
        }
        // Mean load vector of the table. The table holds the load on the moving
        // bearing (as the solver reads it); the shaft load is the opposite.
        double sx = 0.0, sy = 0.0;
        for (const LoadRow& row : load_rows_) {
            if (case_.load_table.format == gui::LoadFormat::MagnitudeDirection) {
                sx += row.a * std::cos(row.b * kPi / 180.0);
                sy += row.a * std::sin(row.b * kPi / 180.0);
            } else {
                sx += row.a;
                sy += row.b;
            }
        }
        sx *= case_.load_table.scale / static_cast<double>(load_rows_.size());
        sy *= case_.load_table.scale / static_cast<double>(load_rows_.size());
        step1.load = std::hypot(sx, sy);
        step1.load_direction_deg = std::fmod(std::atan2(-sy, -sx) * 180.0 / kPi + 360.0, 360.0);
    }
    if (!(step1.load > 0.0)) {
        reason = "the mean load is zero";
        return false;
    }
    const std::string unsupported = gui::unsupported_reason(step1, pair());
    if (!unsupported.empty()) {
        reason = "step 1: " + unsupported;
        return false;
    }
    step2 = case_;
    step2.thermal = false;  // isothermal transient; film temperature set after step 1
    return true;
}

void GuiApp::start_two_step_route() {
    gui::CaseModel step1, step2;
    std::string reason;
    if (!two_step_cases(step1, step2, reason)) {
        set_status("Two-step route not possible: " + reason);
        return;
    }
    if (!save_case(case_path_.empty())) return;
    const fs::path file = case_path_.parent_path() /
                          (case_path_.stem().wstring() + L"_step1_thermal_steady.txt");
    std::string error;
    if (!write_case_file(step1, file, output_dir_for(file), error)) {
        set_status("Cannot write step 1: " + error);
        return;
    }
    RunJob job = make_job(step1, "step 1: thermal steady at the mean load", file, "two-step");
    job.followup = step2;
    enqueue(std::move(job));
}

void GuiApp::draw_thermal_transient_dialog() {
    if (show_thermal_transient_dialog_) {
        ImGui::OpenPopup("Thermal transient");
        show_thermal_transient_dialog_ = false;
    }
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 34.0f, 0.0f), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal("Thermal transient", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return;
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 33.0f);
    const gui::RunTime& direct = estimates_.run_time;
    const auto duration = [](double seconds) {
        char text[32];
        if (seconds < 5400.0) std::snprintf(text, sizeof(text), "%.0f min", seconds / 60.0);
        else std::snprintf(text, sizeof(text), "%.1f h", seconds / 3600.0);
        return std::string(text);
    };
    ImGui::TextUnformatted("This case solves the film, housing and shaft temperatures in time.");
    ImGui::Spacing();
    ImGui::Text("Per revolution: about %s on %d %s.", duration(direct.per_revolution_s).c_str(),
                direct.cores, direct.cores == 1 ? "core" : "cores");
    ImGui::Text("Total for %.0f revolutions: %s.", case_.revolutions, gui::format_run_time(direct).c_str());
    ImGui::Spacing();
    ImGui::TextColored(theme::kRiskMedium,
                       "The housing warms over seconds to minutes; a few revolutions (milliseconds "
                       "each) do not reach the thermal periodic state.");
    ImGui::TextDisabled("Simplified thermal model: polynomial film profile, 2-D housing.");
    ImGui::Spacing();
    ImGui::SeparatorText("Faster route when the thermal state matters");
    gui::CaseModel step1, step2;
    std::string reason;
    const bool route = two_step_cases(step1, step2, reason);
    if (route) {
        const gui::RunTime t1 = gui::estimate_run_time(step1, pair());
        const gui::RunTime t2 = gui::estimate_run_time(step2, pair());
        ImGui::TextUnformatted("1. Thermal steady case at the mean load (finite-volume film, 3-D housing).");
        ImGui::TextDisabled("   %s", gui::format_run_time(t1).c_str());
        ImGui::TextUnformatted("2. Isothermal transient at the mean film temperature of step 1.");
        ImGui::TextDisabled("   %s", gui::format_run_time(t2).c_str());
    } else {
        ImGui::TextColored(theme::kRiskHigh, "Not available: %s.", reason.c_str());
    }
    ImGui::PopTextWrapPos();
    ImGui::Separator();
    if (ImGui::Button("Run thermal transient")) {
        thermal_transient_confirmed_ = true;
        ImGui::CloseCurrentPopup();
        start_run();
    }
    ImGui::SameLine();
    ImGui::BeginDisabled(!route);
    if (ImGui::Button("Run the two-step route")) {
        ImGui::CloseCurrentPopup();
        start_two_step_route();
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
}

void GuiApp::enqueue(RunJob job) {
    queue_.push_back(std::move(job));
    if (!runner_.running() && !active_job_) launch_next();
}

void GuiApp::cancel_all() {
    queue_.clear();
    runner_.cancel();
}

void GuiApp::launch_next() {
    if (queue_.empty() || runner_.running()) return;
    RunJob job = std::move(queue_.front());
    queue_.pop_front();

    std::string source;
    SolverRunner::LaunchSpec spec;
    spec.solver_exe = discover_solver(source);
    if (spec.solver_exe.empty()) {
        set_status("pancake.exe not found. Place it next to this program or set it in Tools > Settings.");
        log_.push_back({"Run aborted: pancake.exe not found.", IM_COL32(185, 35, 35, 255)});
        queue_.clear();
        return;
    }
    if (job.cores > 1) {
        spec.mpiexec = find_mpiexec();
        if (spec.mpiexec.empty()) {
            set_status("mpiexec.exe not found. Install Microsoft MPI or use 1 core.");
            log_.push_back({"Run aborted: mpiexec.exe not found (needed for more than 1 core).",
                            IM_COL32(185, 35, 35, 255)});
            queue_.clear();
            return;
        }
        spec.ranks = job.cores;
    }
    spec.config_file = job.config_file;
    spec.working_dir = job.config_file.parent_path();
    spec.diagnostics_csv = job.output_dir / "diagnostics.csv";
    spec.end_t = job.config.solution_mode == SolutionMode::TRANSIENT ? job.config.end_t : 0.0;
    spec.steady_state = job.config.solution_mode == SolutionMode::STEADY_STATE;
    std::error_code ec;
    fs::remove(spec.diagnostics_csv, ec);

    diagnostics_.clear();
    monitors_.clear();
    stage_text_ = "Starting";
    log_.push_back({"> " + command_line_for(job.config_file, job.cores), IM_COL32(40, 90, 170, 255)});
    std::string error;
    if (!runner_.start(spec, error)) {
        set_status("Launch failed: " + error);
        log_.push_back({"Launch failed: " + error, IM_COL32(185, 35, 35, 255)});
        return;
    }
    active_started_at_ = local_time_string();
    active_job_ = std::move(job);
    results_dir_ = active_job_->output_dir;
    refresh_derived();
}

void GuiApp::process_runner_events() {
    const std::size_t before = log_.size();
    runner_.drain_log(log_);
    runner_.drain_diagnostics(diagnostics_);
    runner_.drain_monitors(monitors_);
    const SolverState previous = progress_.state;
    progress_ = runner_.snapshot();

    // Plain-words stage from the solver log.
    for (std::size_t i = before; i < log_.size(); ++i) {
        const std::string& line = log_[i].text;
        char text[160];
        if (const auto pos = line.find("Equilibrium iter="); pos != std::string::npos) {
            const int iteration = std::atoi(line.c_str() + pos + 17);
            const auto rf = line.find("R_F=");
            const double residual = rf != std::string::npos ? std::atof(line.c_str() + rf + 4) : NAN;
            std::snprintf(text, sizeof(text), "Load balance: iteration %d, force residual %.1e",
                          iteration, residual);
            stage_text_ = text;
        } else if (line.find("Coupled FBNS initialization") != std::string::npos) {
            stage_text_ = "Initial solution (gas-free start)";
        } else if (const auto step = line.find("Step "); step != std::string::npos && active_job_ &&
                                                          active_job_->config.solution_mode ==
                                                              SolutionMode::TRANSIENT) {
            const auto tpos = line.find("t=", step);
            if (tpos != std::string::npos && active_job_->model.speed != 0.0) {
                const double t = std::atof(line.c_str() + tpos + 2);
                const double period = 2.0 * kPi / std::abs(active_job_->model.speed);
                const double turns = t / period;
                std::snprintf(text, sizeof(text), "Revolution %d of %.0f, angle %.0f deg",
                              static_cast<int>(turns) + 1, active_job_->model.revolutions,
                              std::fmod(turns, 1.0) * 360.0);
                stage_text_ = text;
            }
        } else if (line.find("finished successfully") != std::string::npos) {
            stage_text_ = "Finished";
        }
    }
    if (log_.size() > 60000) log_.erase(log_.begin(), log_.begin() + (log_.size() - 50000));

    if (previous == SolverState::Running && progress_.state != SolverState::Running)
        on_job_finished(progress_.state);
}

void GuiApp::on_job_finished(SolverState outcome) {
    if (!active_job_) return;
    RunRecord record;
    record.job = std::move(*active_job_);
    active_job_.reset();
    record.outcome = outcome;
    record.exit_code = progress_.exit_code;
    record.duration_s = progress_.elapsed_s;
    record.started_at = active_started_at_;
    record.command = command_line_for(record.job.config_file, record.job.cores);
    record.converged = outcome == SolverState::Finished &&
                       (diagnostics_.empty() || diagnostics_.back().converged);
    records_.push_back(record);
    const int index = static_cast<int>(records_.size()) - 1;
    shown_record_ = index;

    if (outcome == SolverState::Finished || outcome == SolverState::Cancelled) {
        // Key results of the last step, parsed off the UI thread.
        const fs::path dir = record.job.output_dir;
        const std::string prefix = record.job.config.filename_prefix;
        const SimulationConfig cfg = record.job.config;
        result_tasks_.emplace_back(index, std::async(std::launch::async, [dir, prefix, cfg]() {
            ParsedResult parsed;
            std::string error;
            fs::path folder = dir / "flat";
            std::error_code ec;
            if (!fs::exists(folder / "results.pvd", ec)) folder = dir;
            const auto steps = parse_pvd(folder, prefix, error);
            if (steps.empty()) {
                parsed.error = error.empty() ? "no result steps found" : error;
                return parsed;
            }
            const auto grid = parse_step_grid(steps.back(), error);
            if (!grid) {
                parsed.error = error;
                return parsed;
            }
            parsed.results = gui::compute_key_results(*grid, cfg);
            parsed.ok = true;
            return parsed;
        }));
        selected_step_ = -1;
        request_scan(record.job.output_dir);
    }
    if (outcome == SolverState::Failed)
        set_status("The run failed (exit " + std::to_string(progress_.exit_code) +
                   "); the solver log names the cause.");
    else if (outcome == SolverState::Finished)
        set_status("Run finished in " + std::to_string(static_cast<int>(progress_.elapsed_s)) + " s.");
    refresh_derived();
    if (outcome != SolverState::Cancelled) launch_next();
    else queue_.clear();
}

void GuiApp::poll_result_tasks() {
    for (auto it = result_tasks_.begin(); it != result_tasks_.end();) {
        if (it->second.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
            ++it;
            continue;
        }
        const ParsedResult parsed = it->second.get();
        RunRecord& record = records_[static_cast<std::size_t>(it->first)];
        record.results_ready = parsed.ok;
        record.results = parsed.results;
        record.results_error = parsed.error;
        if (record.job.followup) {
            // Two-step route: step 2 at the mean film temperature of step 1.
            gui::CaseModel step2 = *record.job.followup;
            if (parsed.ok && std::isfinite(parsed.results.t_mean) && record.converged) {
                step2.film_temperature = parsed.results.t_mean;
                const fs::path file = record.job.config_file.parent_path() /
                                      (record.job.config_file.stem().wstring() + L"_step2_isothermal_transient.txt");
                std::string error;
                if (write_case_file(step2, file, output_dir_for(file), error)) {
                    char label[96];
                    std::snprintf(label, sizeof(label), "step 2: isothermal transient at %.1f degC",
                                  parsed.results.t_mean - 273.15);
                    enqueue(make_job(step2, label, file, "two-step"));
                } else {
                    set_status("Cannot write step 2: " + error);
                }
            } else {
                log_.push_back({"Two-step route stopped: step 1 gave no converged mean film temperature.",
                                IM_COL32(185, 35, 35, 255)});
            }
        }
        it = result_tasks_.erase(it);
    }
}

int GuiApp::latest_record(const std::string& group) const {
    for (int i = static_cast<int>(records_.size()) - 1; i >= 0; --i)
        if (records_[static_cast<std::size_t>(i)].job.group == group) return i;
    return -1;
}

// ===========================================================================
// Results loading
// ===========================================================================

void GuiApp::request_scan(const fs::path& output_dir) {
    // The GUI reads the unwrapped (flat) output when present.
    std::error_code ec;
    const fs::path flat = output_dir / "flat";
    results_dir_ = fs::exists(flat / "results.pvd", ec) ? flat : output_dir;
    loader_.scan(results_dir_, config_.filename_prefix.empty() ? "solution" : config_.filename_prefix);
}

void GuiApp::poll_loader() {
    const unsigned revision = loader_.revision();
    if (revision == loader_revision_seen_) {
        if (!loader_.busy() && !steps_.empty()) {
            const int index = selected_step_ < 0 ? static_cast<int>(steps_.size()) - 1 : selected_step_;
            const StepInfo& desired = steps_[static_cast<std::size_t>(std::min(index, static_cast<int>(steps_.size()) - 1))];
            // Compare the step file, not the number: steady snapshots are all step 0.
            if (!grid_ || grid_->step != desired.step || loaded_step_file_ != desired.pvts) {
                loaded_step_file_ = desired.pvts;
                loader_.load_step(desired);
            }
        }
        return;
    }
    loader_revision_seen_ = revision;
    steps_ = loader_.steps();
    grid_ = loader_.grid();
    results_error_ = loader_.error();
    if (selected_step_ >= static_cast<int>(steps_.size())) selected_step_ = -1;
}

