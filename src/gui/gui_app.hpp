#pragma once

// GuiApp owns all application state (no globals) and draws every window.
// Platform plumbing (window, D3D11, DPI, frame loop, PNG capture) lives in
// main_gui.cpp. Window bodies are split by screen; each member function below
// names the file that implements it.
//
// Data flow: panels edit `case_` (engineer inputs) and call mark_changed().
// refresh_derived() then maps the case to `config_` (gui::to_config, onto the
// loaded file's config so values the GUI does not show are kept), validates
// it with the solver's own SimulationConfig::validate(), and recomputes the
// estimates and the workflow status lights. Run is enabled only when every
// setup node is green. The config file is the only case format: saving a
// loaded case patches its text (gui::patch_config_text).

#include "case_model.hpp"
#include "case_status.hpp"
#include "config.hpp"
#include "angle_reference.hpp"
#include "estimates.hpp"
#include "fluid_library.hpp"
#include "key_results.hpp"
#include "results.hpp"
#include "settings.hpp"
#include "solver_runner.hpp"
#include "theme.hpp"
#include "ui_widgets.hpp"

#include <array>
#include <deque>
#include <filesystem>
#include <future>
#include <memory>
#include <optional>
#include <string>
#include <vector>

/// One solver launch: a config file on disk plus what the GUI needs to label it.
struct RunJob {
    std::string label;                    // e.g. "gas on", "load 2500 N"
    std::string group;                    // "single", "comparison", or a batch id
    std::filesystem::path config_file;
    std::filesystem::path output_dir;     // absolute
    int cores = 1;
    gui::CaseModel model;                 // case as run
    SimulationConfig config;              // config as run (tables attached)
    double saturation_pressure = 0.0;     // Pa, for the gas-release check
    double sweep_value = 0.0;             // batch: the swept input (SI)
    /// Two-step thermal route: the isothermal transient queued when this thermal
    /// steady job finishes, at its mean film temperature.
    std::optional<gui::CaseModel> followup;
};

/// A finished (or failed) run and its key results.
struct RunRecord {
    RunJob job;
    SolverState outcome = SolverState::Idle;
    int exit_code = 0;
    double duration_s = 0.0;
    std::string started_at;               // HH:MM:SS
    std::string command;                  // exact command line
    bool converged = false;               // solver reported convergence
    bool results_ready = false;           // key results parsed from the last step
    std::string results_error;
    gui::KeyResults results;
};

class GuiApp {
public:
    /// Wired by main_gui.cpp before the first frame.
    HWND hwnd = nullptr;
    AppPaths paths;
    UserSettings settings;
    Fonts fonts;
    bool wants_exit = false;
    /// Command line: case to open at start, and the screenshot-series folder.
    std::filesystem::path startup_case;
    std::filesystem::path startup_screenshots;
    bool startup_run = false;               // --run: run the opened case, write key results
    /// Scripted session (--run or --screenshots): reads the user's settings and
    /// layout but never writes them, so a test run leaves the user's GUI as it was.
    bool scripted() const { return startup_run || !startup_screenshots.empty(); }
    std::filesystem::path startup_fluids;   // --fluids: fluid directory for this session only

    /// Screen-region PNG capture handshake with the frame loop: the app requests,
    /// the main loop captures the next rendered frame and reports back.
    struct CaptureRequest {
        float min_x = 0, min_y = 0, max_x = 0, max_y = 0;  // max_x <= 0: whole window
        std::filesystem::path path;
    };
    std::optional<CaptureRequest> pending_capture;
    void capture_finished(bool ok, const std::string& message);

    void init();      // load settings, fluid library, last case
    void shutdown();  // persist settings, stop solver
    void draw();      // entire UI, called once per frame

    bool solver_running() const { return runner_.running(); }

private:
    // ===================================================================
    // Case (gui_app.cpp)
    // ===================================================================
    gui::CaseModel case_;                       // engineer inputs (edited by panels)
    std::optional<gui::CaseModel> loaded_case_;  // as read from the file; nullopt = new case
    SimulationConfig base_config_;               // config read from the file, or the template
    std::string base_text_;                      // file text as read; empty = new case
    std::filesystem::path case_path_;            // empty = unsaved
    bool dirty_ = false;
    std::vector<std::string> hidden_keys_;       // settings in the file the GUI does not show
    std::vector<std::string> load_notes_;        // what could not be represented

    /// Call after every edit of case_: marks dirty and refreshes derived state.
    void mark_changed();
    void new_case();
    void open_case_dialog();
    void load_case(const std::filesystem::path& path);
    bool save_case(bool save_as);
    void duplicate_case();
    /// Writes `model` as a config file (new text, or patched when it is the
    /// loaded case). Used by Run, comparison and batch.
    bool write_case_file(const gui::CaseModel& model, const std::filesystem::path& path,
                         const std::filesystem::path& output_dir, std::string& error);

    // Derived state, refreshed by refresh_derived().
    SimulationConfig config_;                    // to_config(case_) with tables attached
    gui::Estimates estimates_;
    gui::StatusArray status_;
    bool allow_out_of_range_ = false;            // user accepted an out-of-table state
    void refresh_derived();
    bool ready_to_run() const { return gui::ready_to_run(status_); }

    // ===================================================================
    // Fluid library (gui_app.cpp)
    // ===================================================================
    gui::FluidLibrary library_;
    std::vector<std::string> pair_ids_;
    gui::FluidPair pair_;
    bool pair_ok_ = false;
    std::string pair_error_;
    const gui::FluidPair* pair() const { return pair_ok_ ? &pair_ : nullptr; }
    std::filesystem::path fluid_root() const;
    void refresh_pairs();
    void select_pair(const std::string& id);    // sets case_.fluid_pair and loads it

    // Transient load table (gui_app.cpp): parsed rows for the plot and checks.
    struct LoadRow { double abscissa, a, b; };
    std::vector<LoadRow> load_rows_;
    bool load_table_ok_ = true;
    std::string load_table_error_;
    void reload_load_table();

    // ===================================================================
    // UI state shared between windows
    // ===================================================================
    gui::Node selected_node_ = gui::Node::Geometry;
    std::array<bool, gui::kNodeCount> advanced_{};  // per-panel Advanced toggle
    gui::AngleZero angle_zero_ = gui::AngleZero::MaxFilm;
    int live_view_tab_ = 0;                      // 0 geometry, 1 mesh
    int highlighted_zone_ = -1;                  // boundary zone lit in the live views
    bool comparison_mode_ = false;               // Run launches gas on + gas off
    bool show_fluid_library_ = false;
    bool show_batch_ = false;
    bool show_settings_ = false;
    bool show_about_ = false;
    bool reset_layout_ = false;
    std::string status_message_;
    double status_message_until_ = 0.0;
    char solver_path_buffer_[512] = {};
    char fluid_dir_buffer_[512] = {};

    /// Display-unit choices (persisted).
    ui::UnitChoices& units() { return settings.unit_choices; }

    // ===================================================================
    // Run control (gui_app.cpp)
    // ===================================================================
    SolverRunner runner_;
    SolverProgress progress_;
    std::vector<LogLine> log_;
    std::vector<DiagRow> diagnostics_;
    std::vector<MonitorRow> monitors_;
    std::deque<RunJob> queue_;
    std::optional<RunJob> active_job_;
    std::string active_started_at_;
    std::vector<RunRecord> records_;             // every finished run this session
    /// Key results parsed off the UI thread after a run: (record index, task).
    struct ParsedResult {
        bool ok = false;
        gui::KeyResults results;
        std::string error;
    };
    std::vector<std::pair<int, std::future<ParsedResult>>> result_tasks_;
    void poll_result_tasks();
    std::string stage_text_;                     // plain-words stage of the active run

    /// Run the current case (or the gas on / gas off pair in comparison mode).
    /// A thermal transient case asks for confirmation first.
    void start_run();
    RunJob make_job(const gui::CaseModel& model, const std::string& label,
                    const std::filesystem::path& file, const std::string& group) const;
    // Thermal transient confirmation (run_control.cpp, plan 0.3).
    bool show_thermal_transient_dialog_ = false;
    bool thermal_transient_confirmed_ = false;
    void draw_thermal_transient_dialog();
    /// Step 1 thermal steady at the mean load, step 2 isothermal transient at
    /// the mean film temperature of step 1. Returns false when not possible.
    bool two_step_cases(gui::CaseModel& step1, gui::CaseModel& step2, std::string& reason) const;
    void start_two_step_route();
    void enqueue(RunJob job);
    void cancel_all();
    void process_runner_events();
    void launch_next();
    void on_job_finished(SolverState outcome);
    /// Exact command line for a config and core count (Copy command).
    std::string command_line_for(const std::filesystem::path& config, int cores) const;
    std::filesystem::path discover_solver(std::string& source) const;
    std::filesystem::path output_dir_for(const std::filesystem::path& config_file) const;
    /// Index of the latest record of a group, or -1.
    int latest_record(const std::string& group) const;

    // ===================================================================
    // Shell (gui_app.cpp)
    // ===================================================================
    void draw_dockspace_and_menu();
    void build_default_layout(unsigned dockspace_id);
    void draw_status_bar();
    void draw_modals();
    void draw_log_window();
    void set_status(const std::string& message);

    // Screenshots of every screen for the user guide (gui_app.cpp).
    int screenshot_step_ = -1;
    int screenshot_results_tab_ = -1;     // results tab a screenshot step selects (-1: leave)
    bool screenshot_units_popup_ = false;  // screenshot step with the Units popup open
    int screenshot_zone_popup_ = -1;       // screenshot step with this boundary zone's values open
    std::filesystem::path screenshot_dir_;
    bool exit_after_screenshots_ = false;    // --screenshots: quit when done
    void start_screenshot_series();
    // --run (gui_app.cpp): 0 idle, 1 waiting to start, 2 running, 3 done.
    int scripted_run_state_ = 0;
    void advance_scripted_run();
    bool write_key_results_json(const RunRecord& record, const std::filesystem::path& file) const;
    void advance_screenshot_series();

    // ===================================================================
    // Setup: workflow tree and panels (workflow_tree.cpp, panel_*.cpp)
    // ===================================================================
    void draw_workflow_window();     // workflow_tree.cpp: tree with lights + estimates box
    void draw_estimates_box();       // workflow_tree.cpp
    void draw_panel_window();        // workflow_tree.cpp: dispatches on selected_node_
    void draw_geometry_panel();      // panel_geometry.cpp
    void draw_fluid_panel();         // panel_fluid.cpp
    void draw_operating_panel();     // panel_operating.cpp
    void draw_boundaries_panel();    // panel_boundaries.cpp (zones + pop-ups)
    void draw_mesh_panel();          // panel_mesh.cpp
    void draw_solution_panel();      // panel_solution.cpp

    // ===================================================================
    // Live views and fluid library dialog (view_geometry.cpp, view_mesh.cpp,
    // dialog_fluid_library.cpp)
    // ===================================================================
    void draw_live_view_window();    // view_geometry.cpp: tabs Geometry | Mesh
    void draw_geometry_view(ImVec2 size);  // view_geometry.cpp
    void draw_mesh_view(ImVec2 size);      // view_mesh.cpp
    void draw_fluid_library_dialog();      // dialog_fluid_library.cpp

    // ===================================================================
    // Results, monitors, batch, report (results_*.cpp, run_monitor.cpp,
    // batch.cpp, case_report.cpp)
    // ===================================================================
    // Results browsing state
    ResultsLoader loader_;
    unsigned loader_revision_seen_ = 0;
    std::vector<StepInfo> steps_;
    std::shared_ptr<const ResultGrid> grid_;
    std::string results_error_;
    std::filesystem::path results_dir_;
    int selected_step_ = -1;                     // -1 = latest
    std::filesystem::path loaded_step_file_;     // step file last requested from the loader
    int shown_record_ = -1;                      // index into records_ shown in Results
    void request_scan(const std::filesystem::path& output_dir);
    void poll_loader();                          // gui_app.cpp

    void draw_results_window();      // results_window.cpp: card, maps, polar, profiles, export
    void draw_monitor_window();      // run_monitor.cpp: live plots, stage, command line
    void draw_batch_window();        // batch.cpp: sweep definition, live table, curves
    bool export_case_report(const std::filesystem::path& html);  // case_report.cpp
};
