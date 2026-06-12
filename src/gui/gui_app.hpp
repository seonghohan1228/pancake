#pragma once

// GuiApp owns all application state (single AppState aggregate, no globals)
// and draws every panel. Platform plumbing (window, D3D11, DPI, frame loop)
// lives in main_gui.cpp and talks to GuiApp through this interface.

#include "config.hpp"
#include "results.hpp"
#include "schema.hpp"
#include "settings.hpp"
#include "solver_runner.hpp"
#include "theme.hpp"

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

/// Snapshot of the most recent solve, feeding the key-results summary.
struct RunRecord {
    int id = 0;
    std::string started_at;  // wall clock HH:MM:SS
    double duration_s = 0.0;
    SolverState outcome = SolverState::Finished;
    int exit_code = 0;
    SimulationConfig config;  // input snapshot

    // Key scalar results (NaN = unavailable).
    double h_min = std::numeric_limits<double>::quiet_NaN();          // m
    double p_max = std::numeric_limits<double>::quiet_NaN();          // Pa
    double attitude_deg = std::numeric_limits<double>::quiet_NaN();   // deg
    double friction_torque = std::numeric_limits<double>::quiet_NaN();  // N.m
    double friction_power = std::numeric_limits<double>::quiet_NaN();   // W
    double cavitated_fraction = std::numeric_limits<double>::quiet_NaN();
    double boundary_flux = std::numeric_limits<double>::quiet_NaN();  // kg/s
};

class GuiApp {
public:
    /// Wired by main_gui.cpp before the first frame.
    HWND hwnd = nullptr;
    AppPaths paths;
    UserSettings settings;
    Fonts fonts;
    bool wants_exit = false;

    /// Plot-region PNG capture handshake with the frame loop: panels request,
    /// main loop captures the next rendered frame and reports back.
    struct CaptureRequest {
        float min_x = 0, min_y = 0, max_x = 0, max_y = 0;
        std::filesystem::path path;
    };
    std::optional<CaptureRequest> pending_capture;
    void capture_finished(bool ok, const std::string& message);

    void init();      // load settings, last config
    void shutdown();  // persist settings, stop solver
    void draw();      // entire UI, called once per frame

    bool solver_running() const { return runner_.running(); }

private:
    // --- case being edited ---
    SimulationConfig config_;
    std::filesystem::path config_path_;  // empty = unsaved case
    bool dirty_ = false;
    std::vector<ValidationIssue> issues_;
    bool has_errors_ = false;
    bool show_advanced_ = false;
    char search_buffer_[96] = {};

    // --- solver run ---
    SolverRunner runner_;
    SolverProgress progress_;
    std::vector<LogLine> log_;
    bool log_autoscroll_ = true;
    std::vector<DiagRow> diagnostics_;
    std::vector<MonitorRow> monitors_;  // physical per-step quantities (DETAILED runs)
    SolverState previous_state_ = SolverState::Idle;
    std::filesystem::path active_output_dir_;
    int run_counter_ = 0;
    std::string run_started_at_;
    double last_detailed_torque_ = std::numeric_limits<double>::quiet_NaN();

    // --- results browsing ---
    ResultsLoader loader_;
    unsigned loader_revision_seen_ = 0;
    std::vector<StepInfo> steps_;
    std::shared_ptr<const ResultGrid> grid_;
    std::string results_error_;
    int selected_step_ = -1;       // index into steps_; -1 = latest
    std::string selected_field_ = "p";
    std::filesystem::path results_dir_;
    int profile_z_index_ = -1;     // -1 = mid-plane
    bool fill_history_when_loaded_ = false;

    // --- last-run records (for the key-results summary) ---
    std::vector<RunRecord> history_;
    SimulationConfig running_config_;  // snapshot taken at launch
    bool setup_panel_focused_ = false;
    bool heatmap_fit_requested_ = false;

    // --- transient UI state ---
    std::string status_message_;
    double status_message_until_ = 0.0;
    bool show_about_ = false;
    bool show_settings_ = false;
    bool reset_layout_ = false;
    char solver_path_buffer_[512] = {};

    // panels
    void draw_dockspace_and_menu();
    void build_default_layout(unsigned dockspace_id);
    void draw_setup_panel();
    void draw_form_tab();
    void draw_inlets_editor();
    void draw_enum_widgets(ParamGroup group);
    void draw_results_panel();
    void draw_heatmap_tab();
    void draw_profiles_tab();
    void draw_summary_tab();
    void draw_convergence_panel();
    void draw_log_panel();
    void draw_status_bar();
    void draw_modals();

    // run lifecycle
    void process_runner_events();
    void start_run();
    void on_run_finished(SolverState outcome);
    void fill_history_results(RunRecord& record);

    // config file handling
    void new_case();
    void open_config_dialog();
    void load_config_from(const std::filesystem::path& path);
    bool save_config(bool save_as);

    // results helpers
    void request_scan(const std::filesystem::path& output_dir);
    void poll_loader();
    void select_step(int index);

    // misc helpers
    void set_status(const std::string& message);
    bool field_invalid(const char* key) const;
    const char* first_error() const;
    std::filesystem::path discover_solver(std::string& source) const;
    std::filesystem::path resolve_output_dir() const;
    /// Selected display-unit index for a field (persisted in settings).
    int& unit_choice(const char* key, UnitFamily family);
    void export_grid_csv();
    void export_profiles_csv();
    void export_monitors_csv();
    /// Opens a save dialog and schedules a capture of the given screen rect
    /// (the plot just drawn) for the next rendered frame.
    void request_plot_capture(const char* suggested_name, float min_x, float min_y,
                              float max_x, float max_y);
};
