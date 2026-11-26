// Application shell: lifecycle, dock layout, menus, status bar, case files,
// fluid library, run queue and the screenshot series. Window bodies live in
// workflow_tree.cpp, panel_*.cpp, view_*.cpp, dialog_fluid_library.cpp,
// results_window.cpp, run_monitor.cpp, batch.cpp and case_report.cpp.

#include "gui_app.hpp"

#include "config_text_patch.hpp"
#include "exports.hpp"
#include "imgui.h"
#include "imgui_internal.h"  // DockBuilder
#include "win_util.hpp"

#include <shellapi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

const wchar_t* kConfigFilter = L"Case files (*.txt)\0*.txt\0All files\0*.*\0";
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

std::string read_text(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    std::ostringstream text;
    text << file.rdbuf();
    return text.str();
}

bool write_text(const fs::path& path, const std::string& text, std::string& error) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        error = "cannot write " + winutil::path_to_utf8(path);
        return false;
    }
    file << text;
    return static_cast<bool>(file);
}

std::string quoted(const fs::path& path) { return "\"" + winutil::path_to_utf8(path) + "\""; }

}  // namespace

// ===========================================================================
// Lifecycle
// ===========================================================================

void GuiApp::init() {
    angle_zero_ = static_cast<gui::AngleZero>(std::clamp(settings.angle_zero, 0, 4));
    library_.set_root(fluid_root());
    refresh_pairs();
    std::snprintf(solver_path_buffer_, sizeof(solver_path_buffer_), "%s", settings.solver_path.c_str());
    std::snprintf(fluid_dir_buffer_, sizeof(fluid_dir_buffer_), "%s", settings.fluid_dir.c_str());

    fs::path start = startup_case;
    if (start.empty() && !settings.last_config.empty()) start = winutil::utf8_to_path(settings.last_config);
    std::error_code ec;
    if (!start.empty() && fs::exists(start, ec)) load_case(start);
    else new_case();
    if (startup_run) {
        scripted_run_state_ = 1;  // screenshots, if requested, follow the run
    } else if (!startup_screenshots.empty()) {
        screenshot_dir_ = startup_screenshots;
        exit_after_screenshots_ = true;
        start_screenshot_series();
    }
}

void GuiApp::shutdown() {
    runner_.cancel();
    settings.angle_zero = static_cast<int>(angle_zero_);
    settings.save(paths.settings_file);
}

void GuiApp::capture_finished(bool ok, const std::string& message) {
    set_status(ok ? "Saved " + message : "PNG export failed: " + message);
}

void GuiApp::draw() {
    process_runner_events();
    poll_loader();
    poll_result_tasks();

    draw_dockspace_and_menu();
    draw_workflow_window();
    draw_panel_window();
    draw_live_view_window();
    draw_results_window();
    draw_monitor_window();
    draw_log_window();
    if (show_batch_) draw_batch_window();
    if (show_fluid_library_) draw_fluid_library_dialog();
    draw_status_bar();
    draw_modals();

    ImGuiIO& io = ImGui::GetIO();
    if (ImGui::IsKeyPressed(ImGuiKey_F5, false) && ready_to_run() && !runner_.running()) start_run();
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_S, false)) save_case(false);
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_O, false)) open_case_dialog();
    if (scripted_run_state_ > 0 && scripted_run_state_ < 3) advance_scripted_run();
    if (screenshot_step_ >= 0) advance_screenshot_series();
}

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
// Case files
// ===========================================================================

void GuiApp::mark_changed() {
    dirty_ = true;
    refresh_derived();
}

void GuiApp::refresh_derived() {
    if (case_.fluid_pair != pair_.id || (!pair_ok_ && pair_error_.empty())) select_pair(case_.fluid_pair);
    config_ = gui::to_config(case_, pair(), base_config_, loaded_case_);
    if (pair()) {
        config_.solubility_curve_table = pair_.solubility;
        config_.kinematic_composition_table = pair_.viscosity;
    }
    if (!case_path_.empty()) config_.config_dir = case_path_.parent_path();
    estimates_ = gui::compute_estimates(case_, pair());

    gui::StatusInputs inputs;
    inputs.model = &case_;
    inputs.pair = pair();
    inputs.pair_error = pair_error_;
    inputs.estimates = &estimates_;
    inputs.config = &config_;
    inputs.allow_out_of_range = allow_out_of_range_;
    inputs.load_table_found = load_table_ok_;
    inputs.load_table_error = load_table_error_;
    const int record = shown_record_;
    if (runner_.running()) {
        inputs.results_light = gui::Light::Idle;
        inputs.results_reason = "Running.";
    } else if (record >= 0 && record < static_cast<int>(records_.size())) {
        const RunRecord& r = records_[static_cast<std::size_t>(record)];
        if (r.outcome == SolverState::Failed) {
            inputs.results_light = gui::Light::Red;
            inputs.results_reason = "The last run failed; see the solver log.";
        } else if (!r.converged) {
            inputs.results_light = gui::Light::Yellow;
            inputs.results_reason = "The last run did not converge; results are not reliable.";
        } else {
            inputs.results_light = gui::Light::Green;
            inputs.results_reason = "Results of the last run.";
        }
    }
    status_ = gui::compute_status(inputs);
}

void GuiApp::new_case() {
    case_ = gui::CaseModel{};
    if (!pair_ids_.empty()) case_.fluid_pair = pair_ids_.front();
    loaded_case_.reset();
    base_config_ = gui::new_case_template();
    base_text_.clear();
    case_path_.clear();
    hidden_keys_.clear();
    load_notes_.clear();
    allow_out_of_range_ = false;
    reload_load_table();
    refresh_derived();
    dirty_ = false;
    set_status("New case.");
}

void GuiApp::open_case_dialog() {
    const fs::path picked = exports::open_file_dialog(hwnd, kConfigFilter);
    if (!picked.empty()) load_case(picked);
}

void GuiApp::load_case(const fs::path& path) {
    const std::string text = read_text(path);
    SimulationConfig cfg;
    cfg.load_from_text(text);
    cfg.config_dir = path.parent_path();
    for (const std::string& error : cfg.parse_errors)
        log_.push_back({"Case file: " + error, IM_COL32(200, 60, 60, 255)});

    const gui::GuiMeta meta = gui::read_gui_meta(text);
    // Pair: from the GUI metadata, else the library folder whose table this config uses.
    std::string pair_id;
    if (const auto it = meta.find("fluid_pair"); it != meta.end()) pair_id = it->second;
    if (pair_id.empty()) pair_id = library_.find_by_table(cfg.solution_table_file);
    select_pair(pair_id);

    const gui::CaseRead read = gui::from_config(cfg, meta, pair());
    case_ = read.model;
    if (case_.fluid_pair.empty()) case_.fluid_pair = pair_id;
    loaded_case_ = case_;
    base_config_ = cfg;
    base_text_ = text;
    const bool same_case = path == case_path_;
    case_path_ = path;
    hidden_keys_ = read.hidden_keys;
    load_notes_ = read.notes;
    if (!same_case) allow_out_of_range_ = false;  // keep it across a re-save
    for (const std::string& note : load_notes_)
        log_.push_back({"Case note: " + note, IM_COL32(160, 110, 10, 255)});
    if (!hidden_keys_.empty())
        log_.push_back({"This case uses " + std::to_string(hidden_keys_.size()) +
                            " settings the GUI does not show; they are kept on save.",
                        IM_COL32(160, 110, 10, 255)});
    if (pair_id.empty() && !cfg.solution_table_file.empty())
        log_.push_back({"The fluid table of this case (" + cfg.solution_table_file.generic_string() +
                            ") is not in the fluid directory. Choose a pair to run it here.",
                        IM_COL32(160, 110, 10, 255)});
    settings.last_config = winutil::path_to_utf8(path);
    reload_load_table();
    refresh_derived();
    dirty_ = false;
    set_status("Opened " + winutil::path_to_utf8(path));
}

bool GuiApp::write_case_file(const gui::CaseModel& model, const fs::path& path,
                             const fs::path& output_dir, std::string& error) {
    SimulationConfig cfg = gui::to_config(model, pair(), base_config_, loaded_case_);
    std::error_code ec;
    fs::path relative = fs::relative(output_dir, path.parent_path(), ec);
    cfg.output_dir = (ec || relative.empty()) ? winutil::path_to_utf8(output_dir)
                                              : relative.generic_string();
    const gui::GuiMeta meta = gui::case_meta(model);
    const std::string text = base_text_.empty()
                                 ? gui::new_config_text(cfg, meta)
                                 : gui::patch_config_text(base_text_, base_config_, cfg, meta);
    return write_text(path, text, error);
}

bool GuiApp::save_case(bool save_as) {
    fs::path target = case_path_;
    if (save_as || target.empty()) {
        target = exports::save_file_dialog(hwnd, kConfigFilter, L"txt", L"case.txt");
        if (target.empty()) return false;
    }
    std::string error;
    if (!write_case_file(case_, target, output_dir_for(target), error)) {
        set_status("Save failed: " + error);
        return false;
    }
    // The saved file becomes the base: later saves patch it.
    load_case(target);
    set_status("Saved " + winutil::path_to_utf8(target));
    return true;
}

void GuiApp::duplicate_case() {
    fs::path suggestion = case_path_.empty() ? fs::path("case_copy.txt")
                                             : case_path_.parent_path() /
                                                   (case_path_.stem().wstring() + L"_copy.txt");
    const fs::path target = exports::save_file_dialog(hwnd, kConfigFilter, L"txt",
                                                      suggestion.filename().wstring());
    if (target.empty()) return;
    std::string error;
    if (!write_case_file(case_, target, output_dir_for(target), error)) {
        set_status("Duplicate failed: " + error);
        return;
    }
    load_case(target);
    set_status("Duplicated as " + winutil::path_to_utf8(target));
}

fs::path GuiApp::output_dir_for(const fs::path& config_file) const {
    return config_file.parent_path() / (config_file.stem().wstring() + L"_results");
}

// ===========================================================================
// Fluid library
// ===========================================================================

fs::path GuiApp::fluid_root() const {
    if (!startup_fluids.empty()) return startup_fluids;
    if (!settings.fluid_dir.empty()) return winutil::utf8_to_path(settings.fluid_dir);
    return winutil::executable_dir() / "fluids";
}

void GuiApp::refresh_pairs() {
    library_.set_root(fluid_root());
    pair_ids_ = library_.list();
}

void GuiApp::select_pair(const std::string& id) {
    pair_ok_ = false;
    pair_error_.clear();
    pair_ = gui::FluidPair{};
    pair_.id = id;
    if (id.empty()) return;
    gui::FluidPair loaded;
    if (library_.load(id, loaded, pair_error_)) {
        pair_ = std::move(loaded);
        pair_ok_ = true;
    }
}

void GuiApp::reload_load_table() {
    load_rows_.clear();
    load_table_ok_ = true;
    load_table_error_.clear();
    if (case_.analysis != gui::AnalysisType::Transient || case_.load_table.file.empty()) return;
    fs::path file = case_.load_table.file;
    if (file.is_relative() && !case_path_.empty()) file = case_path_.parent_path() / file;
    std::ifstream input(file);
    if (!input) {
        load_table_ok_ = false;
        load_table_error_ = "file not found (" + winutil::path_to_utf8(file) + ")";
        return;
    }
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        for (char& ch : line)
            if (ch == ',' || ch == ';' || ch == '\t') ch = ' ';
        std::istringstream cells(line);
        LoadRow row{};
        if (!(cells >> row.abscissa >> row.a >> row.b)) continue;  // header or text line
        load_rows_.push_back(row);
    }
    if (load_rows_.size() < 2) {
        load_table_ok_ = false;
        load_table_error_ = "fewer than two numeric rows";
    }
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

void GuiApp::start_run() {
    if (runner_.running() || !ready_to_run()) return;
    // The solver runs the file on disk: save first so the file equals the GUI
    // state and writes its results where the GUI reads them (<case>_results).
    if (!save_case(case_path_.empty())) return;
    const auto make_job = [&](const gui::CaseModel& model, const std::string& label,
                              const fs::path& file) {
        RunJob job;
        job.label = label;
        job.group = comparison_mode_ ? "comparison" : "single";
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
    };
    if (!comparison_mode_) {
        enqueue(make_job(case_, "case", case_path_));
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
        enqueue(make_job(model, gas ? "gas on (GC+VC)" : "gas off (VC only)", file));
    }
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
            if (!grid_ || grid_->step != desired.step) loader_.load_step(desired);
        }
        return;
    }
    loader_revision_seen_ = revision;
    steps_ = loader_.steps();
    grid_ = loader_.grid();
    results_error_ = loader_.error();
    if (selected_step_ >= static_cast<int>(steps_.size())) selected_step_ = -1;
}

// ===========================================================================
// Shell: menus, layout, status bar, log, modals
// ===========================================================================

void GuiApp::draw_dockspace_and_menu() {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float status_height = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y;
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x, viewport->WorkSize.y - status_height));
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("##host", nullptr,
                 ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus |
                     ImGuiWindowFlags_MenuBar);
    ImGui::PopStyleVar(3);

    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New case")) new_case();
            if (ImGui::MenuItem("Open case...", "Ctrl+O")) open_case_dialog();
            if (ImGui::MenuItem("Save", "Ctrl+S")) save_case(false);
            if (ImGui::MenuItem("Save as...")) save_case(true);
            if (ImGui::MenuItem("Duplicate case...")) duplicate_case();
            if (ImGui::MenuItem("Reload from disk", nullptr, false, !case_path_.empty())) load_case(case_path_);
            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) wants_exit = true;
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Run")) {
            if (ImGui::MenuItem("Run", "F5", false, ready_to_run() && !runner_.running())) start_run();
            if (ImGui::MenuItem("Cancel", nullptr, false, runner_.running())) cancel_all();
            ImGui::Separator();
            if (ImGui::MenuItem("Sweep / batch...")) show_batch_ = true;
            if (ImGui::MenuItem("Copy command line", nullptr, false, !case_path_.empty()))
                ImGui::SetClipboardText(command_line_for(case_path_, case_.cores).c_str());
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Results")) {
            const bool have = !results_dir_.empty() && fs::exists(results_dir_);
            if (ImGui::MenuItem("Open results folder", nullptr, false, have))
                ShellExecuteW(hwnd, L"open", results_dir_.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
            const fs::path pvd = results_dir_ / "results.pvd";
            if (ImGui::MenuItem("Open in ParaView", nullptr, false, have && fs::exists(pvd)))
                ShellExecuteW(hwnd, L"open", pvd.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("View")) {
            ImGui::MenuItem("Fluid library", nullptr, &show_fluid_library_);
            ImGui::MenuItem("Sweep / batch", nullptr, &show_batch_);
            if (ImGui::BeginMenu("Angle zero")) {
                for (gui::AngleZero zero : {gui::AngleZero::Fixed, gui::AngleZero::MaxFilm,
                                            gui::AngleZero::MinFilm, gui::AngleZero::Feed,
                                            gui::AngleZero::OppositeFeed}) {
                    if (ImGui::MenuItem(gui::angle_zero_name(zero), nullptr, angle_zero_ == zero)) {
                        angle_zero_ = zero;
                        settings.angle_zero = static_cast<int>(zero);
                    }
                }
                ImGui::EndMenu();
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Reset window layout")) reset_layout_ = true;
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Tools")) {
            if (ImGui::MenuItem("Settings...")) show_settings_ = true;
            if (ImGui::MenuItem("Save screenshots of every screen...")) {
                const fs::path target = exports::save_file_dialog(hwnd, L"PNG images (*.png)\0*.png\0",
                                                                  L"png", L"screens.png");
                if (!target.empty()) {
                    screenshot_dir_ = target.parent_path();
                    start_screenshot_series();
                }
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Help")) {
            if (ImGui::MenuItem("About")) show_about_ = true;
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    const ImGuiID dockspace_id = ImGui::GetID("pancake_dockspace_v2");
    if (reset_layout_ || ImGui::DockBuilderGetNode(dockspace_id) == nullptr) {
        build_default_layout(dockspace_id);
        reset_layout_ = false;
    }
    ImGui::DockSpace(dockspace_id, ImVec2(0, 0), ImGuiDockNodeFlags_None);
    ImGui::End();
}

void GuiApp::build_default_layout(unsigned dockspace_id) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);
    ImGuiID rest = dockspace_id;
    const ImGuiID left = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Left, 0.19f, nullptr, &rest);
    ImGuiID right = 0;
    ImGuiID middle = ImGui::DockBuilderSplitNode(rest, ImGuiDir_Left, 0.42f, nullptr, &right);
    ImGuiID middle_bottom = 0;
    const ImGuiID middle_top = ImGui::DockBuilderSplitNode(middle, ImGuiDir_Up, 0.56f, nullptr, &middle_bottom);
    ImGuiID right_bottom = 0;
    const ImGuiID right_top = ImGui::DockBuilderSplitNode(right, ImGuiDir_Up, 0.74f, nullptr, &right_bottom);
    ImGui::DockBuilderDockWindow("Workflow", left);
    ImGui::DockBuilderDockWindow("Case setup", middle_top);
    ImGui::DockBuilderDockWindow("Live view", middle_bottom);
    ImGui::DockBuilderDockWindow("Results", right_top);
    ImGui::DockBuilderDockWindow("Run monitor", right_top);
    ImGui::DockBuilderDockWindow("Sweep / batch", right_top);
    ImGui::DockBuilderDockWindow("Solver log", right_bottom);
    ImGui::DockBuilderFinish(dockspace_id);
}

void GuiApp::draw_status_bar() {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float height = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y;
    ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x, viewport->WorkPos.y + viewport->WorkSize.y - height));
    ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x, height));
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(10, 3));
    ImGui::Begin("##statusbar", nullptr,
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoBringToFrontOnFocus);

    // First yellow/red reason of the setup, or the run state.
    std::string message;
    gui::Light light = gui::Light::Green;
    for (int i = 0; i < static_cast<int>(gui::Node::Results); ++i) {
        const gui::NodeStatus& node = status_[static_cast<std::size_t>(i)];
        if (node.light == gui::Light::Red || (node.light == gui::Light::Yellow && light != gui::Light::Red)) {
            if (node.light == gui::Light::Red || message.empty()) {
                light = node.light;
                message = std::string(gui::node_name(static_cast<gui::Node>(i))) + ": " +
                          (node.reasons.empty() ? "" : node.reasons.front());
            }
            if (node.light == gui::Light::Red) break;
        }
    }
    if (runner_.running()) {
        ImGui::TextColored(theme::kAccent, "Running: %s", stage_text_.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("%.0f s", progress_.elapsed_s);
        if (!queue_.empty()) {
            ImGui::SameLine();
            ImGui::TextDisabled("(%zu queued)", queue_.size());
        }
    } else if (!message.empty()) {
        ui::status_dot(light);
        ImGui::SameLine();
        ImGui::TextColored(ui::light_color(light), "%s", message.c_str());
    } else {
        ui::status_dot(gui::Light::Green);
        ImGui::SameLine();
        ImGui::TextColored(theme::kOk, "Ready to run.");
    }

    std::string right = (!status_message_.empty() && ImGui::GetTime() < status_message_until_)
                            ? status_message_
                            : (case_path_.empty() ? std::string("unsaved case")
                                                  : winutil::path_to_utf8(case_path_)) +
                                  (dirty_ ? " *" : "");
    const float width = ImGui::CalcTextSize(right.c_str()).x;
    ImGui::SameLine(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() - width - 12.0f));
    ImGui::TextDisabled("%s", right.c_str());
    ImGui::End();
    ImGui::PopStyleVar(2);
}

void GuiApp::draw_log_window() {
    if (!ImGui::Begin("Solver log")) {
        ImGui::End();
        return;
    }
    if (ImGui::SmallButton("Clear")) log_.clear();
    ImGui::SameLine();
    if (ImGui::SmallButton("Copy all")) {
        std::string all;
        for (const LogLine& line : log_) all += line.text + "\n";
        ImGui::SetClipboardText(all.c_str());
    }
    ImGui::BeginChild("##log", ImVec2(0, 0), ImGuiChildFlags_None, ImGuiWindowFlags_HorizontalScrollbar);
    ImGui::PushFont(fonts.mono);
    ImGuiListClipper clipper;
    clipper.Begin(static_cast<int>(log_.size()));
    while (clipper.Step()) {
        for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; ++i) {
            const LogLine& line = log_[static_cast<std::size_t>(i)];
            if (line.color) ImGui::PushStyleColor(ImGuiCol_Text, line.color);
            ImGui::TextUnformatted(line.text.c_str());
            if (line.color) ImGui::PopStyleColor();
        }
    }
    ImGui::PopFont();
    if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY() - 20.0f) ImGui::SetScrollHereY(1.0f);
    ImGui::EndChild();
    ImGui::End();
}

void GuiApp::draw_modals() {
    if (show_settings_) {
        ImGui::OpenPopup("Settings");
        show_settings_ = false;
    }
    if (ImGui::BeginPopupModal("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("Solver program (pancake.exe)");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 30);
        ImGui::InputTextWithHint("##solver", "pancake.exe next to this program", solver_path_buffer_,
                                 sizeof(solver_path_buffer_));
        ImGui::SameLine();
        if (ImGui::Button("Browse...##solver")) {
            const fs::path picked = exports::open_file_dialog(hwnd, L"Programs (*.exe)\0*.exe\0");
            if (!picked.empty())
                std::snprintf(solver_path_buffer_, sizeof(solver_path_buffer_), "%s",
                              winutil::path_to_utf8(picked).c_str());
        }
        std::string source;
        const fs::path solver = discover_solver(source);
        if (solver.empty()) ImGui::TextColored(theme::kError, "Solver not found.");
        else ImGui::TextColored(theme::kOk, "Using %s (%s)", winutil::path_to_utf8(solver).c_str(), source.c_str());

        ImGui::Spacing();
        ImGui::TextUnformatted("Fluid directory");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 30);
        ImGui::InputTextWithHint("##fluids", "fluids next to this program", fluid_dir_buffer_,
                                 sizeof(fluid_dir_buffer_));
        ImGui::TextDisabled("Now: %s", winutil::path_to_utf8(fluid_root()).c_str());
        ImGui::Separator();
        if (ImGui::Button("OK", ImVec2(120, 0))) {
            settings.solver_path = solver_path_buffer_;
            settings.fluid_dir = fluid_dir_buffer_;
            settings.save(paths.settings_file);
            refresh_pairs();
            select_pair(case_.fluid_pair);
            refresh_derived();
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (show_about_) {
        ImGui::OpenPopup("About pancake");
        show_about_ = false;
    }
    if (ImGui::BeginPopupModal("About pancake", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("pancake GUI 3.0 - journal bearings with refrigerant-oil mixtures");
        ImGui::TextDisabled("Every case is a config file: pancake.exe -c <case.txt> reproduces it.");
        ImGui::TextDisabled("Dear ImGui + ImPlot, Win32 + Direct3D 11. Fonts: Roboto, JetBrains Mono.");
        ImGui::TextDisabled("Settings: %s", winutil::path_to_utf8(paths.data_dir).c_str());
        if (ImGui::Button("Close", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
}

void GuiApp::set_status(const std::string& message) {
    status_message_ = message;
    status_message_until_ = ImGui::GetTime() + 6.0;
}

// ===========================================================================
// Screenshot series (user guide figures)
// ===========================================================================
//
// Steps: each setup node, then the fluid library, the batch window and the
// results. Each step selects a screen, waits a few frames for layout, then
// captures the whole window into <dir>/NN_<name>.png.

namespace {
struct ScreenStep {
    const char* name;
    int node;          // gui::Node to select, or -1
    int live_tab;      // live view tab
    bool fluid_library;
    bool batch;
};
constexpr ScreenStep kScreens[] = {
    {"geometry", 0, 0, false, false},      {"fluid", 1, 0, false, false},
    {"operating", 2, 0, false, false},     {"boundaries", 3, 0, false, false},
    {"mesh", 4, 1, false, false},          {"solution", 5, 0, false, false},
    {"results", 6, 0, false, false},       {"fluid_library", 1, 0, true, false},
    {"batch", 0, 0, false, true},
};
constexpr int kFramesPerScreen = 12;
}  // namespace

void GuiApp::start_screenshot_series() {
    std::error_code ec;
    fs::create_directories(screenshot_dir_, ec);
    screenshot_step_ = 0;
}

void GuiApp::advance_screenshot_series() {
    static int frame = 0;
    const int count = static_cast<int>(sizeof(kScreens) / sizeof(kScreens[0]));
    if (screenshot_step_ >= count) {
        screenshot_step_ = -1;
        show_fluid_library_ = false;
        show_batch_ = false;
        set_status("Screenshots saved to " + winutil::path_to_utf8(screenshot_dir_));
        if (exit_after_screenshots_) wants_exit = true;
        return;
    }
    const ScreenStep& screen = kScreens[screenshot_step_];
    if (frame == 0) {
        selected_node_ = static_cast<gui::Node>(screen.node);
        live_view_tab_ = screen.live_tab;
        show_fluid_library_ = screen.fluid_library;
        show_batch_ = screen.batch;
        if (screen.node == static_cast<int>(gui::Node::Results)) ImGui::SetWindowFocus("Results");
        if (screen.batch) ImGui::SetWindowFocus("Sweep / batch");
    }
    if (++frame == kFramesPerScreen) {
        char name[64];
        std::snprintf(name, sizeof(name), "%02d_%s.png", screenshot_step_ + 1, screen.name);
        CaptureRequest request;  // whole window
        request.path = screenshot_dir_ / name;
        pending_capture = request;
    } else if (frame > kFramesPerScreen + 1) {
        frame = 0;
        ++screenshot_step_;
    }
}
