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


}  // namespace

// ===========================================================================
// Lifecycle
// ===========================================================================

void GuiApp::init() {
    {
        // Physical cores and threads, for the run-time estimate.
        DWORD length = 0;
        GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);
        std::vector<unsigned char> buffer(length);
        int physical = 0, logical = 0;
        if (length && GetLogicalProcessorInformationEx(
                          RelationProcessorCore,
                          reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data()),
                          &length)) {
            for (DWORD offset = 0; offset < length;) {
                const auto* info = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(
                    buffer.data() + offset);
                ++physical;
                for (WORD g = 0; g < info->Processor.GroupCount; ++g) {
                    KAFFINITY mask = info->Processor.GroupMask[g].Mask;
                    while (mask) {
                        logical += static_cast<int>(mask & 1);
                        mask >>= 1;
                    }
                }
                offset += info->Size;
            }
        }
        gui::set_machine_cores(physical, logical);
    }
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
        // Leave room for the right-aligned text; cut the message with "..." and
        // show it whole in a tooltip.
        const float room = ImGui::GetContentRegionAvail().x * 0.62f;
        std::string shown = message;
        if (ImGui::CalcTextSize(shown.c_str()).x > room) {
            while (shown.size() > 4 && ImGui::CalcTextSize((shown + "...").c_str()).x > room)
                shown.pop_back();
            shown += "...";
        }
        ImGui::TextColored(ui::light_color(light), "%s", shown.c_str());
        if (shown != message) ui::item_tooltip(message.c_str());
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
    // Right-aligned text fits the space left after the message (cut from the left).
    ImGui::SameLine();
    const float available = ImGui::GetContentRegionAvail().x - 12.0f;
    std::string right_shown = right;
    if (ImGui::CalcTextSize(right_shown.c_str()).x > available) {
        while (right_shown.size() > 4 && ImGui::CalcTextSize(("..." + right_shown).c_str()).x > available)
            right_shown.erase(right_shown.begin());
        right_shown = "..." + right_shown;
    }
    const float width = ImGui::CalcTextSize(right_shown.c_str()).x;
    ImGui::SameLine(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() - width - 12.0f));
    ImGui::TextDisabled("%s", right_shown.c_str());
    if (right_shown != right) ui::item_tooltip(right.c_str());
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
    draw_thermal_transient_dialog();
    if (show_settings_) {
        ImGui::OpenPopup("Settings");
        show_settings_ = false;
    }
    if (ImGui::BeginPopupModal("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ui::window_shadow();
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
        ui::window_shadow();
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
    {"batch", 0, 0, false, true},          {"run_monitor", 5, 0, false, false},
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
        if (std::string(screen.name) == "run_monitor") ImGui::SetWindowFocus("Run monitor");
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
