#pragma once

// Asynchronous reader for the solver's VTK XML output. Result files are
// treated as untrusted input: any malformed content becomes a readable error
// string, never a crash. All parsing happens on a worker thread; the UI takes
// shared_ptr snapshots.

#include <atomic>
#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

/// One saved timestep discovered from results.pvd / *.pvts files.
struct StepInfo {
    int step = 0;
    double time = 0.0;
    std::filesystem::path pvts;  // parallel collection file for this step
};

/// Assembled global cell-centered grid for one timestep.
/// Fields are row-major over cells: value[j * nx + i], i in theta, j in z.
struct ResultGrid {
    int nx = 0;  // theta cells
    int ny = 0;  // z cells
    int step = 0;
    double time = 0.0;
    std::vector<double> theta_deg;  // cell-center angle, size nx
    std::vector<double> z_m;        // cell-center axial position, size ny
    std::map<std::string, std::vector<double>> fields;

    const std::vector<double>* field(const std::string& name) const {
        const auto it = fields.find(name);
        return it == fields.end() ? nullptr : &it->second;
    }
};

/// Discovers steps and loads grids on a worker thread.
class ResultsLoader {
public:
    ~ResultsLoader();

    /// Discover available steps under output_dir (worker thread).
    void scan(const std::filesystem::path& output_dir, const std::string& prefix);

    /// Parse one step's grid (worker thread). No-op if already loading.
    void load_step(const StepInfo& info);

    bool busy() const { return busy_.load(); }

    /// Snapshot accessors (cheap copies / shared_ptr).
    std::vector<StepInfo> steps() const;
    std::shared_ptr<const ResultGrid> grid() const;
    std::string error() const;
    /// Increments whenever steps/grid/error change, so the UI can react.
    unsigned revision() const { return revision_.load(); }

private:
    void join();

    std::thread worker_;
    std::atomic<bool> busy_{false};
    std::atomic<unsigned> revision_{0};
    mutable std::mutex mutex_;
    std::vector<StepInfo> steps_;
    std::shared_ptr<const ResultGrid> grid_;
    std::string error_;
};

/// Parses one .pvts + its piece .vts files into a grid. Returns nullptr and
/// fills `error` on malformed input. Synchronous; used by ResultsLoader and
/// by run-completion summaries.
std::shared_ptr<ResultGrid> parse_step_grid(const StepInfo& info, std::string& error);

/// Parses a results.pvd collection into step infos (sorted by step).
std::vector<StepInfo> parse_pvd(const std::filesystem::path& output_dir,
                                const std::string& prefix, std::string& error);
