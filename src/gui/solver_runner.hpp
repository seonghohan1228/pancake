#pragma once

// Launches the solver as a subprocess with redirected stdout/stderr and
// tails its diagnostics.csv. All blocking work happens on worker threads;
// the UI thread only takes cheap mutex-protected snapshots.

#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <atomic>
#include <chrono>
#include <deque>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

enum class SolverState { Idle, Running, Finished, Failed, Cancelled };

const char* to_string(SolverState state);

/// One log line with an optional ANSI-derived color (IM_COL32 value; 0 = default).
struct LogLine {
    std::string text;
    unsigned color = 0;
};

/// One row of the solver's diagnostics.csv.
struct DiagRow {
    int step = 0;
    double time = 0.0;
    int outer_iters = 0;
    int flag_flips = 0;
    double max_theta_change = 0.0;
    bool converged = true;
    double liquid_mass = 0.0;
    double liquid_boundary_flux = 0.0;
    double inlet_mass_source = 0.0;
    double liquid_residual = 0.0;
    double gas_residual = 0.0;
    double cavitated_fraction = 0.0;
};

/// UI-facing snapshot of the run state. Copied under mutex once per frame.
struct SolverProgress {
    SolverState state = SolverState::Idle;
    double sim_t = 0.0;     // latest "t=" parsed from stdout
    double end_t = 0.0;     // from config at launch; 0 for steady state
    int step = -1;          // latest "Step N" parsed
    double elapsed_s = 0.0;
    double last_residual = std::numeric_limits<double>::quiet_NaN();
    int exit_code = 0;
    std::string command_line;  // what was launched (for the log/report)
};

class SolverRunner {
public:
    struct LaunchSpec {
        std::filesystem::path solver_exe;
        std::filesystem::path mpiexec;  // empty => direct launch
        int ranks = 1;
        std::filesystem::path config_file;
        std::filesystem::path working_dir;
        std::filesystem::path diagnostics_csv;  // tailed while running; may not exist
        double end_t = 0.0;
        bool steady_state = false;
    };

    ~SolverRunner();

    /// Starts the subprocess and worker threads. Returns false with a message
    /// if CreateProcess fails. Must not be called while running().
    bool start(const LaunchSpec& spec, std::string& error);

    /// Terminates the subprocess (TerminateProcess; the solver has no
    /// cooperative shutdown channel). Safe to call when not running.
    void cancel();

    bool running() const { return state_.load() == SolverState::Running; }

    SolverProgress snapshot() const;

    /// Moves any new log lines / diagnostic rows into the caller's buffers.
    void drain_log(std::vector<LogLine>& out);
    void drain_diagnostics(std::vector<DiagRow>& out);

private:
    void reader_thread_main();
    void diagnostics_thread_main();
    void append_output(const char* data, size_t size);
    void finish_line(std::string&& line);
    void join_threads();
    void close_handles();

    LaunchSpec spec_;
    HANDLE process_ = nullptr;
    HANDLE pipe_read_ = nullptr;
    std::thread reader_thread_;
    std::thread diag_thread_;
    std::atomic<SolverState> state_{SolverState::Idle};
    std::atomic<bool> stop_requested_{false};
    std::chrono::steady_clock::time_point start_time_{};

    mutable std::mutex mutex_;
    std::deque<LogLine> pending_log_;
    std::deque<DiagRow> pending_diag_;
    std::string partial_line_;
    unsigned current_color_ = 0;  // sticky ANSI color across chunks
    SolverProgress progress_;
};
