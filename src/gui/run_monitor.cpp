// Run monitor: the stage in plain words, progress, elapsed time, queue, the exact
// command, live plots (minimum film, force residual, friction torque, cavitated
// area; motion.csv for transient runs) and the live refresh of the Results
// field views from the solver's snapshots. motion.csv is tailed on a worker
// thread; the UI thread only drains parsed lines.

#include "gui_app.hpp"

#include "imgui.h"
#include "implot.h"
#include "result_data.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <thread>

namespace fs = std::filesystem;
using result_data::CsvTable;

namespace {

/// Follows a growing CSV file: polls every 0.5 s and hands over complete new lines.
/// A file whose modification time equals `stale_time` (left by an earlier run) is ignored.
class CsvTailer {
public:
    ~CsvTailer() { stop(); }

    void start(const fs::path& path, fs::file_time_type stale_time) {
        stop();
        {
            std::lock_guard lock(mutex_);
            header_.clear();
            lines_.clear();
            reset_ = true;
        }
        stop_ = false;
        worker_ = std::thread([this, path, stale_time]() { run(path, stale_time); });
    }

    void stop() {
        stop_ = true;
        if (worker_.joinable()) worker_.join();
    }

    /// Moves the lines read so far into `table` (cleared when the file started over).
    void drain(CsvTable& table) {
        std::lock_guard lock(mutex_);
        if (reset_) {
            table = {};
            reset_ = false;
        }
        if (table.headers.empty() && !header_.empty()) table.set_header(header_);
        if (table.headers.empty()) return;
        for (const std::string& line : lines_) table.add_row(line);
        lines_.clear();
    }

private:
    void run(fs::path path, fs::file_time_type stale_time) {
        std::streamoff offset = 0;
        std::string partial;
        bool have_header = false, fresh = false;
        while (!stop_) {
            std::error_code ec;
            if (!fresh && fs::exists(path, ec)) {
                const auto time = fs::last_write_time(path, ec);
                fresh = !ec && time != stale_time;
            }
            std::ifstream file;
            if (fresh) file.open(path, std::ios::binary);
            if (file) {
                file.seekg(0, std::ios::end);
                const std::streamoff size = file.tellg();
                if (size < offset) {  // rewritten from the start
                    offset = 0;
                    partial.clear();
                    have_header = false;
                    std::lock_guard lock(mutex_);
                    header_.clear();
                    lines_.clear();
                    reset_ = true;
                }
                if (size > offset) {
                    std::string chunk(static_cast<std::size_t>(size - offset), '\0');
                    file.seekg(offset);
                    file.read(chunk.data(), static_cast<std::streamsize>(chunk.size()));
                    chunk.resize(static_cast<std::size_t>(file.gcount()));
                    offset += static_cast<std::streamoff>(chunk.size());
                    partial += chunk;
                    std::vector<std::string> complete;
                    std::string header;
                    std::size_t newline;
                    while ((newline = partial.find('\n')) != std::string::npos) {
                        std::string line = partial.substr(0, newline);
                        partial.erase(0, newline + 1);
                        if (!line.empty() && line.back() == '\r') line.pop_back();
                        if (line.empty()) continue;
                        if (!have_header) {
                            header = line;
                            have_header = true;
                        } else {
                            complete.push_back(std::move(line));
                        }
                    }
                    std::lock_guard lock(mutex_);
                    if (!header.empty()) header_ = header;
                    for (auto& line : complete) lines_.push_back(std::move(line));
                }
            }
            for (int i = 0; i < 10 && !stop_; ++i) std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }

    std::thread worker_;
    std::atomic<bool> stop_{false};
    std::mutex mutex_;
    std::string header_;
    std::vector<std::string> lines_;
    bool reset_ = false;
};

struct MonitorUi {
    std::string job_key;
    bool transient = false;
    std::size_t log_cursor = 0;
    std::vector<double> eq_iteration, eq_residual;
    CsvTailer motion_tailer;
    CsvTable motion;
    bool follow_live = true;
    double last_live_check = -10.0;
    fs::path last_snapshot;
    fs::file_time_type last_pvd_time{};
    result_data::ReloadGuard reload;
};

MonitorUi& state() {
    static MonitorUi s;
    return s;
}

std::string duration_text(double seconds) {
    const int total = static_cast<int>(std::max(seconds, 0.0));
    char text[32];
    if (total >= 3600) std::snprintf(text, sizeof(text), "%d h %02d min", total / 3600, (total / 60) % 60);
    else std::snprintf(text, sizeof(text), "%d min %02d s", total / 60, total % 60);
    return text;
}

/// One subplot row: a line, NaN points skipped by ImPlot.
void plot_row(const char* title, const char* x_label, const char* y_label, const std::vector<double>& x,
              const std::vector<double>& y, bool log_y) {
    if (!ImPlot::BeginPlot(title, ImVec2(), ImPlotFlags_NoLegend | ImPlotFlags_NoTitle)) return;
    const bool data = ui::has_points(x, y);
    ImPlot::SetupAxes(x_label, y_label, ui::fit_if(data) | (x_label ? 0 : ImPlotAxisFlags_NoTickLabels),
                      ui::fit_if(data));
    if (!data) {
        ui::empty_plot("Waiting for data");
        ImPlot::EndPlot();
        return;
    }
    if (log_y) ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
    const int n = static_cast<int>(std::min(x.size(), y.size()));
    ImPlot::SetNextLineStyle(theme::kAccent, 1.5f);
    ImPlot::PlotLine("##y", x.data(), y.data(), n);
    ImPlot::EndPlot();
}

}  // namespace

void GuiApp::draw_monitor_window() {
    MonitorUi& s = state();
    s.reload.poll(loader_.busy(), loader_.revision(), steps_, selected_step_, grid_.get(),
                  [&](const StepInfo& info) { loader_.load_step(info); });

    // New job: reset the parsed series and start tailing its motion.csv.
    const std::string key = active_job_ ? active_job_->config_file.string() + "|" + active_started_at_ : std::string();
    if (key != s.job_key) {
        s.job_key = key;
        if (active_job_) {
            s.eq_iteration.clear();
            s.eq_residual.clear();
            s.log_cursor = 0;
            for (std::size_t i = log_.size(); i-- > 0;)
                if (log_[i].text.rfind("> ", 0) == 0) {
                    s.log_cursor = i;
                    break;
                }
            s.last_snapshot.clear();
            s.last_pvd_time = {};
            s.transient = active_job_->config.solution_mode == SolutionMode::TRANSIENT;
            if (s.transient) {
                const fs::path motion = active_job_->output_dir / "motion.csv";
                std::error_code ec;
                const auto stale = fs::exists(motion, ec) ? fs::last_write_time(motion, ec) : fs::file_time_type::min();
                s.motion = {};
                s.motion_tailer.start(motion, stale);
            } else {
                s.motion_tailer.stop();
                s.motion = {};
            }
        } else {
            s.motion_tailer.drain(s.motion);
            s.motion_tailer.stop();
        }
    }
    if (active_job_ && s.transient) s.motion_tailer.drain(s.motion);

    // Force residual of the load balance from "Equilibrium iter=K R_F=..." lines.
    if (s.log_cursor > log_.size()) s.log_cursor = log_.size();
    if (active_job_) {
        for (; s.log_cursor < log_.size(); ++s.log_cursor) {
            const std::string& line = log_[s.log_cursor].text;
            const auto pos = line.find("Equilibrium iter=");
            if (pos == std::string::npos) continue;
            const auto rf = line.find("R_F=", pos);
            if (rf == std::string::npos) continue;
            s.eq_iteration.push_back(std::atof(line.c_str() + pos + 17));
            s.eq_residual.push_back(std::abs(std::atof(line.c_str() + rf + 4)));
        }
    }

    // Live field views: newest steady snapshot, or the growing transient output.
    const double now = ImGui::GetTime();
    if (active_job_ && s.follow_live && now - s.last_live_check >= 1.0 && !loader_.busy()) {
        s.last_live_check = now;
        std::error_code ec;
        if (!s.transient) {
            fs::path newest;
            for (fs::directory_iterator it(active_job_->output_dir / "iterations", ec), end; !ec && it != end;
                 it.increment(ec))
                if (it->is_directory(ec) && (newest.empty() || it->path().filename() > newest.filename()))
                    newest = it->path();
            if (!newest.empty() && newest != s.last_snapshot &&
                (fs::exists(newest / "results.pvd", ec) || fs::exists(newest / "flat" / "results.pvd", ec))) {
                s.last_snapshot = newest;
                shown_record_ = -1;
                selected_step_ = -1;
                s.reload.arm(loader_.revision(), grid_.get());
                request_scan(newest);
            }
        } else {
            fs::path pvd = active_job_->output_dir / "flat" / "results.pvd";
            if (!fs::exists(pvd, ec)) pvd = active_job_->output_dir / "results.pvd";
            const auto time = fs::last_write_time(pvd, ec);
            if (!ec && time != s.last_pvd_time) {
                s.last_pvd_time = time;
                shown_record_ = -1;
                selected_step_ = -1;
                request_scan(active_job_->output_dir);
            }
        }
    }

    if (!ImGui::Begin("Run monitor")) {
        ImGui::End();
        return;
    }
    if (!active_job_) {
        ImGui::TextDisabled("No run active.");
        if (!records_.empty()) {
            const RunRecord& last = records_.back();
            ImGui::TextDisabled("Last run: %s, %s after %s.", last.job.label.c_str(), to_string(last.outcome),
                                duration_text(last.duration_s).c_str());
        }
        if (!queue_.empty()) ImGui::TextDisabled("Queued: %d", static_cast<int>(queue_.size()));
    } else {
        const RunJob& job = *active_job_;
        const bool load_balance = !s.transient && job.model.position == gui::ShaftPosition::LoadBalance;
        ImGui::Text("Stage: %s", stage_text_.empty() ? "Starting" : stage_text_.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("(%s)", job.label.c_str());

        float fraction = 0.0f;
        char overlay[128];
        const char* help = nullptr;
        if (s.transient && progress_.end_t > 0.0) {
            fraction = static_cast<float>(std::clamp(progress_.sim_t / progress_.end_t, 0.0, 1.0));
            std::snprintf(overlay, sizeof(overlay), "t = %.4g of %.4g s (%.0f %%)", progress_.sim_t, progress_.end_t,
                          fraction * 100.0f);
            help = "Simulated time of the transient run.";
        } else if (load_balance) {
            const int max_iterations = std::max(1, job.model.convergence.outer_max_iterations);
            const int iteration = s.eq_iteration.empty() ? 0 : static_cast<int>(s.eq_iteration.back());
            fraction = std::clamp(static_cast<float>(iteration) / max_iterations, 0.0f, 1.0f);
            std::snprintf(overlay, sizeof(overlay), "Load balance iteration %d of at most %d", iteration, max_iterations);
            help = "The bar counts against the maximum number of iterations, an upper bound: a run usually "
                   "converges well before it.";
        } else {
            std::snprintf(overlay, sizeof(overlay), "Iterating until the residuals fall below the tolerances");
            help = "Fixed shaft position: no fixed number of iterations is known in advance.";
        }
        ImGui::ProgressBar(fraction, ImVec2(-FLT_MIN, 0), overlay);
        ui::item_tooltip(help);
        ImGui::Text("Elapsed %s", duration_text(progress_.elapsed_s).c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("|  %d job(s) queued", static_cast<int>(queue_.size()));
        ImGui::SameLine();
        if (ImGui::SmallButton("Cancel")) cancel_all();
        ui::item_tooltip("Stop the running solver and clear the queue.");

        const std::string command = command_line_for(job.config_file, job.cores);
        ImGui::AlignTextToFramePadding();
        ImGui::TextDisabled("Command");
        ImGui::SameLine();
        if (fonts.mono) ImGui::PushFont(fonts.mono);
        ImGui::TextWrapped("%s", command.c_str());
        if (fonts.mono) ImGui::PopFont();
        if (ImGui::SmallButton("Copy command")) {
            ImGui::SetClipboardText(command.c_str());
            set_status("Command copied.");
        }
        ui::item_tooltip("Runs the same case from a terminal.");
        ImGui::SameLine();
        ImGui::Checkbox("Show live fields in Results", &s.follow_live);
        ui::item_tooltip("Load the newest field output into the Results window while the run goes on.");
        if (s.transient) {
            ImGui::TextDisabled("Field views refresh at every field output.");
        } else if (job.config.steady_output_interval > 0) {
            ImGui::TextDisabled("Field views refresh every %d outer iterations.", job.config.steady_output_interval);
        } else {
            ImGui::TextDisabled("Field views appear when the run ends.");
        }
    }

    // Live plots (the last run stays visible after it ends).
    const float height = std::max(320.0f, ImGui::GetContentRegionAvail().y);
    const auto& hu = result_data::unit_of(units(), result_data::kUnitFilm, gui::UnitFamily::SmallLength);
    const auto& tu = result_data::unit_of(units(), result_data::kUnitTorque, gui::UnitFamily::Torque);
    const std::string h_label = std::string("h_min [") + hu.label + "]";
    const std::string t_label = std::string("Torque [") + tu.label + "]";
    const auto* motion_time = s.motion.column("time_s");
    if (s.transient && motion_time && !motion_time->empty()) {
        std::vector<double> h, torque, cavitated, ecc;
        if (const auto* c = s.motion.column("h_min_m")) h = result_data::scaled(*c, hu);
        if (const auto* c = s.motion.column("friction_torque_Nm")) torque = result_data::scaled(*c, tu);
        if (const auto* c = s.motion.column("cavitated_fraction"))
            for (double v : *c) cavitated.push_back(v * 100.0);
        if (const auto* c = s.motion.column("eccentricity_ratio")) ecc = *c;
        if (ImPlot::BeginSubplots("##live_transient", 4, 1, ImVec2(-1, height), ImPlotSubplotFlags_LinkAllX)) {
            plot_row("h_min", nullptr, h_label.c_str(), *motion_time, h, false);
            plot_row("ecc", nullptr, "Ecc. ratio", *motion_time, ecc, false);
            plot_row("torque", nullptr, t_label.c_str(), *motion_time, torque, false);
            plot_row("cav", "Time [s]", "Cav. area [%]", *motion_time, cavitated, false);
            ImPlot::EndSubplots();
        }
        ImGui::End();
        return;
    }
    std::vector<double> mx, mh, mt, dx, dc;
    for (const MonitorRow& row : monitors_) {
        mx.push_back(s.transient ? row.t : row.step);
        mh.push_back(gui::to_display(hu, row.h_min));
        mt.push_back(gui::to_display(tu, row.torque));
    }
    for (const DiagRow& row : diagnostics_) {
        dx.push_back(s.transient ? row.time : row.step);
        dc.push_back(row.cavitated_fraction * 100.0);
    }
    const char* x_label = s.transient ? "Time [s]" : "Outer iteration / step";
    const int rows = s.transient ? 3 : 4;
    if (ImPlot::BeginSubplots("##live", rows, 1, ImVec2(-1, height), ImPlotSubplotFlags_LinkAllX)) {
        plot_row("h_min", nullptr, h_label.c_str(), mx, mh, false);
        if (!s.transient) plot_row("force", nullptr, "Force resid.", s.eq_iteration, s.eq_residual, true);
        plot_row("torque", nullptr, t_label.c_str(), mx, mt, false);
        plot_row("cav", x_label, "Cav. area [%]", dx, dc, false);
        ImPlot::EndSubplots();
    }
    ImGui::End();
}
