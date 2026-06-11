#include "solver_runner.hpp"

#include "imgui.h"
#include "win_util.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

const char* to_string(SolverState state) {
    switch (state) {
        case SolverState::Idle: return "idle";
        case SolverState::Running: return "running";
        case SolverState::Finished: return "finished";
        case SolverState::Failed: return "failed";
        case SolverState::Cancelled: return "cancelled";
    }
    return "?";
}

namespace {

unsigned ansi_code_to_color(int code) {
    switch (code) {
        case 31: case 91: return IM_COL32(235, 100, 100, 255);  // red
        case 32: case 92: return IM_COL32(120, 200, 120, 255);  // green
        case 33: case 93: return IM_COL32(230, 190, 80, 255);   // yellow
        case 34: case 94: return IM_COL32(110, 160, 240, 255);  // blue
        case 35: case 95: return IM_COL32(200, 120, 220, 255);  // magenta
        case 36: case 96: return IM_COL32(90, 200, 210, 255);   // cyan
        case 0: default: return 0;
    }
}

std::string quote_path(const std::filesystem::path& path) {
    std::string text = winutil::path_to_utf8(path);
    std::string quoted = "\"";
    for (char ch : text) {
        if (ch == '"') quoted += '\\';
        quoted += ch;
    }
    quoted += '"';
    return quoted;
}

}  // namespace

SolverRunner::~SolverRunner() {
    cancel();
    join_threads();
    close_handles();
}

bool SolverRunner::start(const LaunchSpec& spec, std::string& error) {
    if (running()) {
        error = "A solve is already running.";
        return false;
    }
    join_threads();
    close_handles();

    spec_ = spec;
    stop_requested_ = false;

    std::string command;
    if (!spec.mpiexec.empty() && spec.ranks > 1) {
        command = quote_path(spec.mpiexec) + " -n " + std::to_string(spec.ranks) + " " +
                  quote_path(spec.solver_exe) + " -c " + quote_path(spec.config_file);
    } else {
        command = quote_path(spec.solver_exe) + " -c " + quote_path(spec.config_file);
    }

    SECURITY_ATTRIBUTES security{};
    security.nLength = sizeof(security);
    security.bInheritHandle = TRUE;

    HANDLE pipe_write = nullptr;
    if (!CreatePipe(&pipe_read_, &pipe_write, &security, 0)) {
        error = "CreatePipe failed.";
        return false;
    }
    SetHandleInformation(pipe_read_, HANDLE_FLAG_INHERIT, 0);

    STARTUPINFOW startup{};
    startup.cb = sizeof(startup);
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = pipe_write;
    startup.hStdError = pipe_write;
    startup.hStdInput = INVALID_HANDLE_VALUE;

    std::wstring command_wide = winutil::utf8_to_wide(command);
    std::vector<wchar_t> command_buffer(command_wide.begin(), command_wide.end());
    command_buffer.push_back(L'\0');

    PROCESS_INFORMATION process_info{};
    const BOOL created = CreateProcessW(
        nullptr, command_buffer.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr,
        spec.working_dir.empty() ? nullptr : spec.working_dir.c_str(), &startup,
        &process_info);
    CloseHandle(pipe_write);  // child owns the write end now

    if (!created) {
        const DWORD last_error = GetLastError();
        CloseHandle(pipe_read_);
        pipe_read_ = nullptr;
        error = "CreateProcess failed (error " + std::to_string(last_error) + ").";
        return false;
    }
    CloseHandle(process_info.hThread);
    process_ = process_info.hProcess;
    start_time_ = std::chrono::steady_clock::now();

    {
        std::lock_guard lock(mutex_);
        progress_ = SolverProgress{};
        progress_.state = SolverState::Running;
        progress_.end_t = spec.steady_state ? 0.0 : spec.end_t;
        progress_.command_line = command;
        partial_line_.clear();
        current_color_ = 0;
        pending_log_.push_back({"> " + command, IM_COL32(140, 150, 165, 255)});
    }
    state_ = SolverState::Running;

    reader_thread_ = std::thread(&SolverRunner::reader_thread_main, this);
    diag_thread_ = std::thread(&SolverRunner::diagnostics_thread_main, this);
    return true;
}

void SolverRunner::cancel() {
    if (!running()) return;
    stop_requested_ = true;
    if (process_) TerminateProcess(process_, 130);
}

SolverProgress SolverRunner::snapshot() const {
    std::lock_guard lock(mutex_);
    SolverProgress copy = progress_;
    if (copy.state == SolverState::Running) {
        copy.elapsed_s = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start_time_)
                             .count();
    }
    return copy;
}

void SolverRunner::drain_log(std::vector<LogLine>& out) {
    std::lock_guard lock(mutex_);
    while (!pending_log_.empty()) {
        out.push_back(std::move(pending_log_.front()));
        pending_log_.pop_front();
    }
}

void SolverRunner::drain_diagnostics(std::vector<DiagRow>& out) {
    std::lock_guard lock(mutex_);
    while (!pending_diag_.empty()) {
        out.push_back(pending_diag_.front());
        pending_diag_.pop_front();
    }
}

void SolverRunner::reader_thread_main() {
    char buffer[4096];
    DWORD bytes_read = 0;
    while (ReadFile(pipe_read_, buffer, sizeof(buffer), &bytes_read, nullptr) &&
           bytes_read > 0) {
        append_output(buffer, bytes_read);
    }

    WaitForSingleObject(process_, INFINITE);
    DWORD exit_code = 0;
    GetExitCodeProcess(process_, &exit_code);

    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time_)
            .count();

    SolverState final_state;
    if (stop_requested_) {
        final_state = SolverState::Cancelled;
    } else if (exit_code == 0) {
        final_state = SolverState::Finished;
    } else {
        final_state = SolverState::Failed;
    }

    {
        std::lock_guard lock(mutex_);
        if (!partial_line_.empty()) finish_line(std::move(partial_line_));
        progress_.state = final_state;
        progress_.exit_code = static_cast<int>(exit_code);
        progress_.elapsed_s = elapsed;
        if (final_state == SolverState::Finished) {
            if (progress_.end_t > 0.0) progress_.sim_t = progress_.end_t;
        }
        char summary[160];
        if (final_state == SolverState::Cancelled) {
            std::snprintf(summary, sizeof(summary),
                          "-- cancelled by user after %.1f s (solver terminated; the last "
                          "output file may be incomplete)",
                          elapsed);
        } else {
            std::snprintf(summary, sizeof(summary), "-- solver exited with code %lu after %.1f s",
                          static_cast<unsigned long>(exit_code), elapsed);
        }
        pending_log_.push_back(
            {summary, final_state == SolverState::Finished
                          ? IM_COL32(120, 200, 120, 255)
                          : IM_COL32(235, 150, 90, 255)});
    }
    state_ = final_state;
}

void SolverRunner::append_output(const char* data, size_t size) {
    std::lock_guard lock(mutex_);
    for (size_t i = 0; i < size; ++i) {
        const char ch = data[i];
        if (ch == '\r') continue;
        if (ch == '\n') {
            finish_line(std::move(partial_line_));
            partial_line_.clear();
        } else {
            partial_line_ += ch;
        }
    }
}

// Strips ANSI escapes (keeping the color), updates progress, queues the line.
// Caller holds mutex_.
void SolverRunner::finish_line(std::string&& line) {
    std::string clean;
    clean.reserve(line.size());
    unsigned color = current_color_;
    for (size_t i = 0; i < line.size();) {
        if (line[i] == '\033' && i + 1 < line.size() && line[i + 1] == '[') {
            size_t j = i + 2;
            int code = 0;
            while (j < line.size() && isdigit(static_cast<unsigned char>(line[j]))) {
                code = code * 10 + (line[j] - '0');
                ++j;
            }
            if (j < line.size() && line[j] == 'm') {
                current_color_ = ansi_code_to_color(code);
                if (clean.empty()) color = current_color_;
                i = j + 1;
                continue;
            }
        }
        clean += line[i];
        ++i;
    }

    // Progress: "Step <n> t=<t>" lines (indented by the solver logger), plus
    // the final success marker.
    const size_t first_char = clean.find_first_not_of(' ');
    if (first_char != std::string::npos && clean.compare(first_char, 5, "Step ") == 0) {
        progress_.step = std::atoi(clean.c_str() + first_char + 5);
        const size_t t_pos = clean.find("t=", first_char);
        if (t_pos != std::string::npos) {
            progress_.sim_t = std::atof(clean.c_str() + t_pos + 2);
        }
    }
    if (clean.find("finished successfully") != std::string::npos && progress_.end_t > 0.0) {
        progress_.sim_t = progress_.end_t;
    }

    pending_log_.push_back({std::move(clean), color});
    if (pending_log_.size() > 50000) pending_log_.pop_front();
}

// Tails diagnostics.csv while the solver runs: re-reads new bytes every 300 ms,
// parses complete rows, and feeds residuals to the progress snapshot.
void SolverRunner::diagnostics_thread_main() {
    const std::filesystem::path csv_path = spec_.diagnostics_csv;
    if (csv_path.empty()) return;

    std::streamoff offset = 0;
    std::string carry;
    bool header_skipped = false;

    const auto parse_pending = [&](const std::string& chunk) {
        carry += chunk;
        size_t start = 0;
        while (true) {
            const size_t newline = carry.find('\n', start);
            if (newline == std::string::npos) break;
            std::string row = carry.substr(start, newline - start);
            start = newline + 1;
            if (!row.empty() && row.back() == '\r') row.pop_back();
            if (row.empty()) continue;
            if (!header_skipped) {
                header_skipped = true;
                if (row.rfind("step,", 0) == 0) continue;
            }
            DiagRow diag;
            // Columns: step,time,outer_iters,flag_flips,max_theta_change,converged,
            // liquid_mass,liquid_boundary_flux,inlet_mass_source,clamp_mass_liquid,
            // liquid_residual,gas_mass,gas_boundary_flux,clamp_mass_gas,gas_residual,
            // cavitated_fraction
            std::stringstream ss(row);
            std::string cell;
            std::vector<double> values;
            while (std::getline(ss, cell, ',')) {
                values.push_back(cell.empty() ? 0.0 : std::atof(cell.c_str()));
            }
            if (values.size() < 16) continue;
            diag.step = static_cast<int>(values[0]);
            diag.time = values[1];
            diag.outer_iters = static_cast<int>(values[2]);
            diag.flag_flips = static_cast<int>(values[3]);
            diag.max_theta_change = values[4];
            diag.converged = values[5] != 0.0;
            diag.liquid_mass = values[6];
            diag.liquid_boundary_flux = values[7];
            diag.inlet_mass_source = values[8];
            diag.liquid_residual = values[10];
            diag.gas_residual = values[14];
            diag.cavitated_fraction = values[15];

            std::lock_guard lock(mutex_);
            progress_.last_residual = diag.liquid_residual;
            pending_diag_.push_back(diag);
            if (pending_diag_.size() > 100000) pending_diag_.pop_front();
        }
        carry.erase(0, start);
    };

    while (state_.load() == SolverState::Running) {
        std::ifstream file(csv_path, std::ios::binary);
        if (file.is_open()) {
            file.seekg(0, std::ios::end);
            const std::streamoff size = file.tellg();
            if (size < offset) {  // truncated: a fresh run reuses the file
                offset = 0;
                carry.clear();
                header_skipped = false;
            }
            if (size > offset) {
                file.seekg(offset);
                std::string chunk(static_cast<size_t>(size - offset), '\0');
                file.read(chunk.data(), static_cast<std::streamsize>(chunk.size()));
                offset = size;
                parse_pending(chunk);
            }
        }
        for (int i = 0; i < 6 && state_.load() == SolverState::Running; ++i) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }
}

void SolverRunner::join_threads() {
    if (reader_thread_.joinable()) reader_thread_.join();
    if (diag_thread_.joinable()) diag_thread_.join();
}

void SolverRunner::close_handles() {
    if (pipe_read_) {
        CloseHandle(pipe_read_);
        pipe_read_ = nullptr;
    }
    if (process_) {
        CloseHandle(process_);
        process_ = nullptr;
    }
}
