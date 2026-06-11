// Application shell: lifecycle, dock layout, menu, status bar, run control,
// config file handling, exports. Panel bodies live in panels_setup.cpp and
// panels_results.cpp.

#include "gui_app.hpp"

#include "exports.hpp"
#include "imgui.h"
#include "imgui_internal.h"  // DockBuilder
#include "win_util.hpp"

#include <shellapi.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace fs = std::filesystem;

namespace {

const wchar_t* kConfigFilter = L"Config files (*.txt)\0*.txt\0All files\0*.*\0";

std::filesystem::path find_mpiexec() {
    wchar_t buffer[MAX_PATH] = {};
    const DWORD env_len = GetEnvironmentVariableW(L"MSMPI_BIN", buffer, MAX_PATH);
    if (env_len > 0 && env_len < MAX_PATH) {
        const fs::path candidate = fs::path(buffer) / L"mpiexec.exe";
        if (fs::exists(candidate)) return candidate;
    }
    for (const wchar_t* prefix :
         {L"C:/Program Files/Microsoft MPI/Bin/mpiexec.exe",
          L"C:/Program Files (x86)/Microsoft SDKs/MPI/Bin/mpiexec.exe"}) {
        if (fs::exists(prefix)) return prefix;
    }
    wchar_t found[MAX_PATH] = {};
    if (SearchPathW(nullptr, L"mpiexec.exe", nullptr, MAX_PATH, found, nullptr) > 0) {
        return found;
    }
    return {};
}

std::string local_time_string() {
    SYSTEMTIME now{};
    GetLocalTime(&now);
    char text[16];
    std::snprintf(text, sizeof(text), "%02u:%02u:%02u", now.wHour, now.wMinute, now.wSecond);
    return text;
}

/// Runnable defaults for a typical journal bearing. The raw SimulationConfig
/// defaults fail the solver's validation (p_cav and the axial boundary
/// pressures default to 0; the solver requires absolute pressures > 0).
SimulationConfig default_case() {
    SimulationConfig config;
    config.p_cav = 101325.0;
    config.bc_z_south_val = 101325.0;
    config.bc_z_north_val = 101325.0;
    return config;
}

}  // namespace

void GuiApp::init() {
    config_ = default_case();
    settings.load(paths.settings_file);
    if (!settings.last_config.empty()) {
        const fs::path last = winutil::utf8_to_path(settings.last_config);
        if (fs::exists(last)) load_config_from(last);
    }
    if (config_path_.empty()) {
        // Fall back to a config.txt next to the exe if one is shipped there.
        const fs::path sibling = winutil::executable_dir() / "config.txt";
        if (fs::exists(sibling)) load_config_from(sibling);
    }
    sync_raw_from_config();
    std::snprintf(solver_path_buffer_, sizeof(solver_path_buffer_), "%s",
                  settings.solver_path.c_str());

    // Pick up results from a previous session right away.
    const fs::path output_dir = resolve_output_dir();
    std::error_code ec;
    if (fs::exists(output_dir / "results.pvd", ec)) request_scan(output_dir);
}

void GuiApp::shutdown() {
    runner_.cancel();
    settings.save(paths.settings_file);
}

void GuiApp::capture_finished(bool ok, const std::string& message) {
    set_status(ok ? "Saved " + message : "PNG export failed: " + message);
}

void GuiApp::draw() {
    process_runner_events();
    poll_loader();

    issues_ = validate_config(config_);
    has_errors_ = std::any_of(issues_.begin(), issues_.end(),
                              [](const ValidationIssue& issue) { return issue.is_error; });

    draw_dockspace_and_menu();
    draw_setup_panel();
    draw_results_panel();
    draw_convergence_panel();
    draw_log_panel();
    draw_history_panel();
    draw_status_bar();
    draw_modals();

    // Keyboard flow: F5 / Enter (outside text inputs) runs, Ctrl+S saves.
    ImGuiIO& io = ImGui::GetIO();
    const bool can_run = !has_errors_ && !runner_.running();
    if (ImGui::IsKeyPressed(ImGuiKey_F5, false) && can_run) start_run();
    if (setup_panel_focused_ && !io.WantTextInput &&
        ImGui::IsKeyPressed(ImGuiKey_Enter, false) && can_run) {
        start_run();
    }
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_S, false)) save_config(false);
}

void GuiApp::draw_dockspace_and_menu() {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const float status_height = ImGui::GetFrameHeight() + ImGui::GetStyle().WindowPadding.y;
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x, viewport->WorkSize.y - status_height));
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    const ImGuiWindowFlags host_flags =
        ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse |
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus |
        ImGuiWindowFlags_MenuBar;
    ImGui::Begin("##host", nullptr, host_flags);
    ImGui::PopStyleVar(3);

    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New case", nullptr)) new_case();
            if (ImGui::MenuItem("Open config...", "Ctrl+O")) open_config_dialog();
            if (ImGui::MenuItem("Save", "Ctrl+S")) save_config(false);
            if (ImGui::MenuItem("Save as...")) save_config(true);
            if (ImGui::MenuItem("Reload from disk", nullptr, false, !config_path_.empty())) {
                load_config_from(config_path_);
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) wants_exit = true;
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Run")) {
            if (ImGui::MenuItem("Run", "F5", false, !has_errors_ && !runner_.running())) {
                start_run();
            }
            if (ImGui::MenuItem("Cancel run", nullptr, false, runner_.running())) {
                runner_.cancel();
            }
            ImGui::Separator();
            ImGui::SetNextItemWidth(ImGui::GetFontSize() * 6);
            if (ImGui::InputInt("MPI ranks", &settings.mpi_ranks)) {
                settings.mpi_ranks = std::clamp(settings.mpi_ranks, 1, 256);
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Results")) {
            if (ImGui::MenuItem("Refresh results")) {
                request_scan(resolve_output_dir());
            }
            const bool have_dir = !results_dir_.empty() && fs::exists(results_dir_);
            if (ImGui::MenuItem("Open results folder", nullptr, false, have_dir)) {
                ShellExecuteW(hwnd, L"open", results_dir_.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
            }
            const fs::path pvd = results_dir_ / "results.pvd";
            if (ImGui::MenuItem("Open in ParaView", nullptr, false, have_dir && fs::exists(pvd))) {
                ShellExecuteW(hwnd, L"open", pvd.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Tools")) {
            if (ImGui::MenuItem("Settings...")) show_settings_ = true;
            if (ImGui::MenuItem("Open data folder (settings, layout)")) {
                ShellExecuteW(hwnd, L"open", paths.data_dir.c_str(), nullptr, nullptr,
                              SW_SHOWNORMAL);
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Reset window layout")) reset_layout_ = true;
            ImGui::MenuItem("Show advanced parameters", nullptr, &show_advanced_);
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Help")) {
            if (ImGui::MenuItem("About")) show_about_ = true;
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    const ImGuiID dockspace_id = ImGui::GetID("pancake_dockspace");
    if (reset_layout_ || ImGui::DockBuilderGetNode(dockspace_id) == nullptr) {
        build_default_layout(dockspace_id);
        reset_layout_ = false;
    }
    ImGui::DockSpace(dockspace_id, ImVec2(0, 0), ImGuiDockNodeFlags_None);
    ImGui::End();

    if (ImGui::IsKeyPressed(ImGuiKey_O, false) && ImGui::GetIO().KeyCtrl) open_config_dialog();
}

void GuiApp::build_default_layout(unsigned dockspace_id) {
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);

    ImGuiID center = dockspace_id;
    const ImGuiID left = ImGui::DockBuilderSplitNode(center, ImGuiDir_Left, 0.30f, nullptr, &center);
    const ImGuiID right = ImGui::DockBuilderSplitNode(center, ImGuiDir_Right, 0.26f, nullptr, &center);
    const ImGuiID bottom = ImGui::DockBuilderSplitNode(center, ImGuiDir_Down, 0.30f, nullptr, &center);

    ImGui::DockBuilderDockWindow("Case Setup", left);
    ImGui::DockBuilderDockWindow("Results", center);
    ImGui::DockBuilderDockWindow("Convergence", center);
    ImGui::DockBuilderDockWindow("Run History", right);
    ImGui::DockBuilderDockWindow("Solver Log", bottom);
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
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking |
                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                     ImGuiWindowFlags_NoBringToFrontOnFocus);

    // Solver state chip: idle / running step N / finished / failed / cancelled.
    ImVec4 chip_color = theme::kMuted;
    std::string state_text = to_string(progress_.state);
    switch (progress_.state) {
        case SolverState::Running:
            chip_color = theme::kAccent;
            state_text = "running";
            if (progress_.step >= 0) state_text += "  step " + std::to_string(progress_.step);
            break;
        case SolverState::Finished: chip_color = theme::kOk; state_text = "converged / finished"; break;
        case SolverState::Failed:
            chip_color = theme::kError;
            state_text = "failed (exit " + std::to_string(progress_.exit_code) + ")";
            break;
        case SolverState::Cancelled: chip_color = theme::kWarn; break;
        case SolverState::Idle: break;
    }
    ImGui::TextColored(chip_color, "[%s]", state_text.c_str());

    ImGui::SameLine();
    {
        std::string solver_source;
        if (discover_solver(solver_source).empty()) {
            ImGui::TextColored(theme::kWarn, "solver not found");
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                ImGui::SetTooltip("pancake.exe was not found.\nPlace it next to this program or "
                                  "set the path in Tools > Settings.");
            }
            ImGui::SameLine();
        }
    }
    if (progress_.state == SolverState::Running || progress_.elapsed_s > 0.0) {
        ImGui::TextDisabled("elapsed %.1f s", progress_.elapsed_s);
        ImGui::SameLine();
    }
    if (std::isfinite(progress_.last_residual)) {
        ImGui::TextDisabled("residual %.2e", progress_.last_residual);
        ImGui::SameLine();
    }
    if (progress_.state == SolverState::Running && progress_.end_t > 0.0) {
        const float fraction =
            static_cast<float>(std::clamp(progress_.sim_t / progress_.end_t, 0.0, 1.0));
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12);
        char overlay[48];
        std::snprintf(overlay, sizeof(overlay), "t = %.4g / %.4g s", progress_.sim_t,
                      progress_.end_t);
        ImGui::ProgressBar(fraction, ImVec2(ImGui::GetFontSize() * 12, 0), overlay);
        ImGui::SameLine();
    }

    // Right side: transient status message, then config path.
    std::string right_text;
    if (!status_message_.empty() && ImGui::GetTime() < status_message_until_) {
        right_text = status_message_;
    } else {
        right_text = config_path_.empty() ? std::string("unsaved case")
                                          : winutil::path_to_utf8(config_path_);
        if (dirty_) right_text += " *";
    }
    const float text_width = ImGui::CalcTextSize(right_text.c_str()).x;
    ImGui::SameLine(std::max(ImGui::GetCursorPosX(),
                             ImGui::GetWindowWidth() - text_width - 12.0f));
    ImGui::TextDisabled("%s", right_text.c_str());

    ImGui::End();
    ImGui::PopStyleVar(2);
}

void GuiApp::draw_modals() {
    if (show_settings_) {
        ImGui::OpenPopup("Settings");
        show_settings_ = false;
    }
    if (ImGui::BeginPopupModal("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("Solver executable");
        ImGui::SetNextItemWidth(ImGui::GetFontSize() * 28);
        ImGui::InputTextWithHint("##solverpath", "auto-discover (pancake.exe next to this app)",
                                 solver_path_buffer_, sizeof(solver_path_buffer_));
        ImGui::SameLine();
        if (ImGui::Button("Browse...")) {
            const fs::path picked = exports::open_file_dialog(
                hwnd, L"Executables (*.exe)\0*.exe\0All files\0*.*\0");
            if (!picked.empty()) {
                std::snprintf(solver_path_buffer_, sizeof(solver_path_buffer_), "%s",
                              winutil::path_to_utf8(picked).c_str());
            }
        }
        std::string source;
        const fs::path discovered = discover_solver(source);
        if (discovered.empty()) {
            ImGui::TextColored(theme::kError,
                               "Solver not found. Place pancake.exe next to this program or set "
                               "a path above.");
        } else {
            ImGui::TextColored(theme::kOk, "Using: %s (%s)",
                               winutil::path_to_utf8(discovered).c_str(), source.c_str());
        }
        ImGui::Separator();
        if (ImGui::Button("OK", ImVec2(120, 0))) {
            settings.solver_path = solver_path_buffer_;
            settings.save(paths.settings_file);
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            std::snprintf(solver_path_buffer_, sizeof(solver_path_buffer_), "%s",
                          settings.solver_path.c_str());
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    if (show_about_) {
        ImGui::OpenPopup("About Pancake GUI");
        show_about_ = false;
    }
    if (ImGui::BeginPopupModal("About Pancake GUI", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted("Pancake GUI 2.0");
        ImGui::TextDisabled("Journal-bearing thin-film solver front end.");
        ImGui::Separator();
        ImGui::TextDisabled("Dear ImGui (docking) + ImPlot, Win32 + Direct3D 11.");
        ImGui::TextDisabled("Fonts: Roboto (Apache-2.0), JetBrains Mono (OFL-1.1), embedded.");
        ImGui::TextDisabled("PNG export: stb_image_write (public domain).");
        ImGui::TextDisabled("Settings & layout: %s", winutil::path_to_utf8(paths.data_dir).c_str());
        if (ImGui::Button("Close", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
}

// ---------------------------------------------------------------------------
// Run lifecycle
// ---------------------------------------------------------------------------

void GuiApp::process_runner_events() {
    const size_t log_before = log_.size();
    runner_.drain_log(log_);
    runner_.drain_diagnostics(diagnostics_);
    progress_ = runner_.snapshot();

    // Friction torque only appears in DETAILED per-step lines ("| M=<x> Nm").
    for (size_t i = log_before; i < log_.size(); ++i) {
        const std::string& line = log_[i].text;
        const size_t torque_pos = line.find("| M=");
        if (torque_pos != std::string::npos) {
            last_detailed_torque_ = std::atof(line.c_str() + torque_pos + 4);
        }
    }
    if (log_.size() > 60000) log_.erase(log_.begin(), log_.begin() + (log_.size() - 50000));

    if (progress_.state != previous_state_) {
        if (previous_state_ == SolverState::Running) on_run_finished(progress_.state);
        previous_state_ = progress_.state;
    }
}

void GuiApp::start_run() {
    if (runner_.running() || has_errors_) return;

    // The solver reads a config file from disk; persist the case first.
    if (config_path_.empty()) {
        const fs::path case_dir = paths.data_dir / "case";
        std::error_code ec;
        fs::create_directories(case_dir, ec);
        config_path_ = case_dir / "config.txt";
    }
    if (!save_config(false)) return;

    std::string source;
    SolverRunner::LaunchSpec spec;
    spec.solver_exe = discover_solver(source);
    if (spec.solver_exe.empty()) {
        set_status("Solver not found - place pancake.exe next to this program or set the path "
                   "in Tools > Settings.");
        log_.push_back({"Run aborted: pancake.exe not found. Place it next to this program or "
                        "set the path in Tools > Settings.",
                        IM_COL32(235, 100, 100, 255)});
        return;
    }
    if (settings.mpi_ranks > 1) {
        spec.mpiexec = find_mpiexec();
        if (spec.mpiexec.empty()) {
            set_status("mpiexec.exe not found - install Microsoft MPI or set ranks to 1.");
            log_.push_back({"Run aborted: mpiexec.exe not found (needed for MPI ranks > 1). "
                            "Install the Microsoft MPI runtime or set ranks to 1.",
                            IM_COL32(235, 100, 100, 255)});
            return;
        }
        spec.ranks = settings.mpi_ranks;
    }

    spec.config_file = config_path_;
    spec.working_dir = config_path_.parent_path();
    active_output_dir_ = resolve_output_dir();
    spec.diagnostics_csv = active_output_dir_ / "diagnostics.csv";
    spec.end_t = config_.end_t;
    spec.steady_state = config_.solution_mode == SolutionMode::STEADY_STATE;

    // Remove a stale diagnostics file so the tail only sees this run's rows.
    std::error_code ec;
    fs::remove(spec.diagnostics_csv, ec);

    diagnostics_.clear();
    last_detailed_torque_ = std::numeric_limits<double>::quiet_NaN();
    running_config_ = config_;
    run_started_at_ = local_time_string();

    std::string error;
    if (!runner_.start(spec, error)) {
        set_status("Launch failed: " + error);
        log_.push_back({"Launch failed: " + error, IM_COL32(235, 100, 100, 255)});
        return;
    }
    set_status("Solver started.");
}

void GuiApp::on_run_finished(SolverState outcome) {
    RunRecord record;
    record.id = ++run_counter_;
    record.started_at = run_started_at_;
    record.duration_s = progress_.elapsed_s;
    record.outcome = outcome;
    record.exit_code = progress_.exit_code;
    record.config = running_config_;
    record.friction_torque = last_detailed_torque_;
    if (std::isfinite(last_detailed_torque_)) {
        record.friction_power = std::abs(last_detailed_torque_ * running_config_.omega);
    }
    if (!diagnostics_.empty()) {
        record.cavitated_fraction = diagnostics_.back().cavitated_fraction;
        record.boundary_flux = diagnostics_.back().liquid_boundary_flux;
        record.conv_step.reserve(diagnostics_.size());
        record.conv_residual.reserve(diagnostics_.size());
        for (const DiagRow& row : diagnostics_) {
            record.conv_step.push_back(row.step);
            record.conv_residual.push_back(std::abs(row.liquid_residual));
        }
    }
    history_.push_back(std::move(record));

    if (outcome == SolverState::Finished || outcome == SolverState::Cancelled) {
        fill_history_when_loaded_ = true;
        selected_step_ = -1;  // jump to the latest step of the fresh run
        request_scan(active_output_dir_);
    }
    if (outcome == SolverState::Failed) {
        set_status("Solver failed (exit " + std::to_string(progress_.exit_code) +
                   ") - see the Solver Log for the error.");
    } else if (outcome == SolverState::Finished) {
        set_status("Run finished in " + std::to_string(progress_.elapsed_s).substr(0, 5) + " s.");
    }
}

void GuiApp::fill_history_results(RunRecord& record) {
    if (!grid_) return;
    const auto* h = grid_->field("h");
    const auto* p = grid_->field("p");
    if (h && !h->empty()) {
        double h_min = std::numeric_limits<double>::infinity();
        for (double value : *h) {
            if (std::isfinite(value)) h_min = std::min(h_min, value);
        }
        if (std::isfinite(h_min)) record.h_min = h_min;
    }
    if (p && !p->empty()) {
        double p_max = -std::numeric_limits<double>::infinity();
        for (double value : *p) {
            if (std::isfinite(value)) p_max = std::max(p_max, value);
        }
        if (std::isfinite(p_max)) record.p_max = p_max;
    }

    // The solver writes these resultants as constant cell fields when enabled.
    const auto first_cell = [&](const char* name) {
        const auto* values = grid_->field(name);
        return (values && !values->empty()) ? (*values)[0]
                                            : std::numeric_limits<double>::quiet_NaN();
    };
    const double attitude = first_cell("bearing_attitude_angle");
    record.attitude_deg = std::isfinite(attitude) ? attitude : record.config.attitude_angle_deg;
    const double torque = first_cell("friction_torque");
    if (std::isfinite(torque)) {
        record.friction_torque = torque;
        record.friction_power = std::abs(torque * record.config.omega);
    }

    // Mid-plane profiles for overlay comparison.
    if (p && grid_->nx > 0 && grid_->ny > 0) {
        const int j = grid_->ny / 2;
        record.profile_theta_deg = grid_->theta_deg;
        record.profile_p.resize(grid_->nx);
        for (int i = 0; i < grid_->nx; ++i) {
            record.profile_p[i] = (*p)[static_cast<size_t>(j) * grid_->nx + i];
        }
        if (h) {
            record.profile_h.resize(grid_->nx);
            for (int i = 0; i < grid_->nx; ++i) {
                record.profile_h[i] = (*h)[static_cast<size_t>(j) * grid_->nx + i];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Config file handling
// ---------------------------------------------------------------------------

void GuiApp::new_case() {
    config_ = default_case();
    config_path_.clear();
    dirty_ = false;
    sync_raw_from_config();
    set_status("New case with default parameters.");
}

void GuiApp::open_config_dialog() {
    const fs::path picked = exports::open_file_dialog(hwnd, kConfigFilter);
    if (!picked.empty()) load_config_from(picked);
}

void GuiApp::load_config_from(const fs::path& path) {
    SimulationConfig fresh;
    fresh.load_from_file(path);
    config_ = std::move(fresh);
    config_path_ = path;
    dirty_ = false;
    sync_raw_from_config();
    settings.last_config = winutil::path_to_utf8(path);
    for (const std::string& warning : config_.parse_warnings) {
        log_.push_back({"Config warning: " + warning, IM_COL32(230, 190, 80, 255)});
    }
    set_status("Loaded " + winutil::path_to_utf8(path));
}

bool GuiApp::save_config(bool save_as) {
    fs::path target = config_path_;
    if (save_as || target.empty()) {
        target = exports::save_file_dialog(hwnd, kConfigFilter, L"txt", L"config.txt");
        if (target.empty()) return false;
    }
    try {
        config_.save_to_file(target);
    } catch (const std::exception& ex) {
        set_status(std::string("Save failed: ") + ex.what());
        return false;
    }
    config_path_ = target;
    dirty_ = false;
    settings.last_config = winutil::path_to_utf8(target);
    set_status("Saved " + winutil::path_to_utf8(target));
    return true;
}

void GuiApp::sync_raw_from_config() {
    raw_text_ = config_.to_config_text();
    raw_dirty_ = false;
}

// ---------------------------------------------------------------------------
// Results helpers
// ---------------------------------------------------------------------------

void GuiApp::request_scan(const fs::path& output_dir) {
    results_dir_ = output_dir;
    loader_.scan(output_dir, config_.filename_prefix);
}

void GuiApp::poll_loader() {
    const unsigned revision = loader_.revision();
    if (revision == loader_revision_seen_) {
        // Still kick a pending step load if the worker is free.
        if (!loader_.busy() && !steps_.empty()) {
            const int desired_index =
                selected_step_ < 0 ? static_cast<int>(steps_.size()) - 1 : selected_step_;
            const StepInfo& desired = steps_[static_cast<size_t>(desired_index)];
            if (!grid_ || grid_->step != desired.step) loader_.load_step(desired);
        }
        return;
    }
    loader_revision_seen_ = revision;
    steps_ = loader_.steps();
    grid_ = loader_.grid();
    results_error_ = loader_.error();
    if (selected_step_ >= static_cast<int>(steps_.size())) selected_step_ = -1;

    if (fill_history_when_loaded_ && grid_ && !history_.empty() && !steps_.empty() &&
        grid_->step == steps_.back().step) {
        fill_history_results(history_.back());
        fill_history_when_loaded_ = false;
    }
}

void GuiApp::select_step(int index) {
    selected_step_ = index;
    if (steps_.empty()) return;
    const int clamped =
        index < 0 ? static_cast<int>(steps_.size()) - 1
                  : std::clamp(index, 0, static_cast<int>(steps_.size()) - 1);
    loader_.load_step(steps_[static_cast<size_t>(clamped)]);
}

// ---------------------------------------------------------------------------
// Misc helpers
// ---------------------------------------------------------------------------

void GuiApp::set_status(const std::string& message) {
    status_message_ = message;
    status_message_until_ = ImGui::GetTime() + 6.0;
}

bool GuiApp::field_invalid(const char* key) const {
    for (const ValidationIssue& issue : issues_) {
        if (issue.is_error && issue.field == key) return true;
    }
    return false;
}

const char* GuiApp::first_error() const {
    for (const ValidationIssue& issue : issues_) {
        if (issue.is_error) return issue.message.c_str();
    }
    return nullptr;
}

std::filesystem::path GuiApp::discover_solver(std::string& source) const {
    if (!settings.solver_path.empty()) {
        fs::path override_path = winutil::utf8_to_path(settings.solver_path);
        if (override_path.is_relative()) {
            override_path = winutil::executable_dir() / override_path;
        }
        if (fs::exists(override_path)) {
            source = "settings override";
            return override_path;
        }
    }
    const fs::path exe_dir = winutil::executable_dir();
    for (const wchar_t* name : {L"pancake.exe", L"solver.exe"}) {
        const fs::path candidate = exe_dir / name;
        if (fs::exists(candidate)) {
            source = "next to the GUI";
            return candidate;
        }
    }
    source = "not found";
    return {};
}

std::filesystem::path GuiApp::resolve_output_dir() const {
    fs::path output_dir = winutil::utf8_to_path(config_.output_dir);
    if (output_dir.is_relative()) {
        const fs::path base =
            config_path_.empty() ? winutil::executable_dir() : config_path_.parent_path();
        output_dir = base / output_dir;
    }
    return output_dir;
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------

void GuiApp::export_history_csv() {
    const fs::path target =
        exports::save_file_dialog(hwnd, L"CSV files (*.csv)\0*.csv\0", L"csv", L"run_history.csv");
    if (target.empty()) return;
    std::vector<std::vector<std::string>> rows;
    for (const RunRecord& record : history_) {
        char cell[64];
        std::vector<std::string> row;
        row.push_back(std::to_string(record.id));
        row.push_back(record.started_at);
        row.push_back(to_string(record.outcome));
        const auto push_double = [&](double value) {
            std::snprintf(cell, sizeof(cell), "%.9g", value);
            row.push_back(std::isfinite(value) ? cell : "");
        };
        push_double(record.duration_s);
        push_double(record.config.omega);
        push_double(record.config.e);
        push_double(record.config.mu);
        push_double(std::hypot(record.config.external_load_x, record.config.external_load_y));
        row.push_back(std::to_string(record.config.n_theta_global) + "x" +
                      std::to_string(record.config.n_z_global));
        push_double(record.h_min);
        push_double(record.p_max);
        push_double(record.attitude_deg);
        push_double(record.friction_torque);
        push_double(record.friction_power);
        push_double(record.cavitated_fraction);
        rows.push_back(std::move(row));
    }
    std::string error;
    if (exports::write_csv_rows(target,
                                {"run", "started", "outcome", "duration_s", "omega_rad_s", "e_m",
                                 "mu_Pa_s", "load_N", "mesh", "h_min_m", "p_max_Pa",
                                 "attitude_deg", "friction_torque_Nm", "friction_power_W",
                                 "cavitated_fraction"},
                                rows, error)) {
        set_status("Saved " + winutil::path_to_utf8(target));
    } else {
        set_status(error);
    }
}

void GuiApp::export_grid_csv() {
    if (!grid_) return;
    const auto* field = grid_->field(selected_field_);
    if (!field) return;
    const fs::path target = exports::save_file_dialog(
        hwnd, L"CSV files (*.csv)\0*.csv\0", L"csv",
        winutil::utf8_to_wide(selected_field_ + "_step" + std::to_string(grid_->step) + ".csv"));
    if (target.empty()) return;
    std::vector<double> theta_column, z_column, value_column;
    theta_column.reserve(field->size());
    for (int j = 0; j < grid_->ny; ++j) {
        for (int i = 0; i < grid_->nx; ++i) {
            theta_column.push_back(grid_->theta_deg[i]);
            z_column.push_back(grid_->z_m[j]);
            value_column.push_back((*field)[static_cast<size_t>(j) * grid_->nx + i]);
        }
    }
    std::string error;
    if (exports::write_csv_columns(target, {"theta_deg", "z_m", selected_field_},
                                   {theta_column, z_column, value_column}, error)) {
        set_status("Saved " + winutil::path_to_utf8(target));
    } else {
        set_status(error);
    }
}

void GuiApp::export_profiles_csv() {
    if (!grid_) return;
    const fs::path target = exports::save_file_dialog(
        hwnd, L"CSV files (*.csv)\0*.csv\0", L"csv",
        winutil::utf8_to_wide("profiles_step" + std::to_string(grid_->step) + ".csv"));
    if (target.empty()) return;
    const int j = profile_z_index_ < 0 ? grid_->ny / 2
                                       : std::clamp(profile_z_index_, 0, grid_->ny - 1);
    std::vector<std::string> headers = {"theta_deg"};
    std::vector<std::vector<double>> columns = {grid_->theta_deg};
    for (const char* name : {"p", "h", "film_content", "T"}) {
        const auto* field = grid_->field(name);
        if (!field) continue;
        std::vector<double> column(grid_->nx);
        for (int i = 0; i < grid_->nx; ++i) {
            column[i] = (*field)[static_cast<size_t>(j) * grid_->nx + i];
        }
        headers.push_back(name);
        columns.push_back(std::move(column));
    }
    for (const RunRecord& record : history_) {
        if (!record.overlay || record.profile_p.empty()) continue;
        headers.push_back("p_run" + std::to_string(record.id));
        columns.push_back(record.profile_p);
    }
    std::string error;
    if (exports::write_csv_columns(target, headers, columns, error)) {
        set_status("Saved " + winutil::path_to_utf8(target));
    } else {
        set_status(error);
    }
}

void GuiApp::export_convergence_csv() {
    if (diagnostics_.empty()) return;
    const fs::path target = exports::save_file_dialog(hwnd, L"CSV files (*.csv)\0*.csv\0", L"csv",
                                                      L"convergence.csv");
    if (target.empty()) return;
    std::vector<double> steps, times, outer, residuals, cavitated;
    for (const DiagRow& row : diagnostics_) {
        steps.push_back(row.step);
        times.push_back(row.time);
        outer.push_back(row.outer_iters);
        residuals.push_back(row.liquid_residual);
        cavitated.push_back(row.cavitated_fraction);
    }
    std::string error;
    if (exports::write_csv_columns(
            target, {"step", "time_s", "outer_iters", "liquid_residual", "cavitated_fraction"},
            {steps, times, outer, residuals, cavitated}, error)) {
        set_status("Saved " + winutil::path_to_utf8(target));
    } else {
        set_status(error);
    }
}

void GuiApp::request_plot_capture(const char* suggested_name, float min_x, float min_y,
                                  float max_x, float max_y) {
    const fs::path target =
        exports::save_file_dialog(hwnd, L"PNG images (*.png)\0*.png\0", L"png",
                                  winutil::utf8_to_wide(suggested_name));
    if (target.empty()) return;
    CaptureRequest request;
    request.min_x = min_x;
    request.min_y = min_y;
    request.max_x = max_x;
    request.max_y = max_y;
    request.path = target;
    pending_capture = request;
}
