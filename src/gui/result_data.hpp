#pragma once

// Data helpers shared by the results window, run monitor, batch window and case
// report: numeric CSV tables (motion.csv, equilibrium_history.csv), the column
// order of a result grid by display angle, mid-plane profiles, the gas-fraction
// field, display units of result quantities, field-map and profile extraction,
// async readers for motion.csv and comparison runs, and a reload guard for the
// results loader. Header-only (no CMake entry needed); no solver calls.

#include "angle_reference.hpp"
#include "config.hpp"
#include "key_results.hpp"
#include "results.hpp"
#include "units.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace result_data {

/// Display-unit ids of result quantities (the Results "Units" pop-up edits them).
inline constexpr const char* kUnitPressure = "result.pressure";
inline constexpr const char* kUnitFilm = "result.film";
inline constexpr const char* kUnitTorque = "result.torque";
inline constexpr const char* kUnitPower = "result.power";
inline constexpr const char* kUnitTemperature = "result.temperature";

/// The display unit chosen for `id` (first option when unset).
inline const gui::UnitOption& unit_of(std::map<std::string, int>& choices, const char* id,
                                      gui::UnitFamily family) {
    const auto& options = gui::unit_options(family);
    int& choice = choices[id];
    if (choice < 0 || choice >= static_cast<int>(options.size())) choice = 0;
    return options[static_cast<std::size_t>(choice)];
}

/// Numeric CSV table stored by column; the first line names the columns.
struct CsvTable {
    std::vector<std::string> headers;
    std::vector<std::vector<double>> columns;

    std::size_t rows() const { return columns.empty() ? 0 : columns.front().size(); }

    const std::vector<double>* column(const std::string& name) const {
        for (std::size_t i = 0; i < headers.size(); ++i)
            if (headers[i] == name) return &columns[i];
        return nullptr;
    }

    void set_header(std::string line) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        headers.clear();
        std::size_t start = 0;
        while (start <= line.size()) {
            const std::size_t comma = line.find(',', start);
            std::string cell = line.substr(start, comma == std::string::npos ? std::string::npos
                                                                            : comma - start);
            while (!cell.empty() && (cell.back() == ' ' || cell.back() == '\t')) cell.pop_back();
            while (!cell.empty() && (cell.front() == ' ' || cell.front() == '\t')) cell.erase(0, 1);
            headers.push_back(cell);
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        columns.assign(headers.size(), {});
    }

    /// Appends one data line; false (and nothing added) when it has too few cells.
    bool add_row(const std::string& line) {
        if (headers.empty()) return false;
        std::vector<double> values;
        values.reserve(headers.size());
        const char* cursor = line.c_str();
        while (*cursor && values.size() < headers.size()) {
            char* end = nullptr;
            const double value = std::strtod(cursor, &end);
            if (end == cursor) {
                // Non-numeric cell (e.g. "true"/"nan" text): read as NaN up to the next comma.
                const char* comma = std::strchr(cursor, ',');
                values.push_back(std::numeric_limits<double>::quiet_NaN());
                if (!comma) break;
                cursor = comma + 1;
                continue;
            }
            values.push_back(value);
            cursor = end;
            while (*cursor == ' ' || *cursor == '\t') ++cursor;
            if (*cursor == ',') ++cursor;
            else break;
        }
        if (values.size() < headers.size()) return false;
        for (std::size_t i = 0; i < headers.size(); ++i) columns[i].push_back(values[i]);
        return true;
    }
};

/// Reads a whole CSV file (header + numeric rows). Call off the UI thread for large files.
inline bool read_csv(const std::filesystem::path& path, CsvTable& table, std::string& error) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        error = "cannot open " + path.filename().string();
        return false;
    }
    std::string line;
    if (!std::getline(file, line)) {
        error = path.filename().string() + " is empty";
        return false;
    }
    table.set_header(line);
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) table.add_row(line);
    }
    return true;
}

/// Gas fraction of the film: alpha_gas when written and nonzero somewhere (gas transport
/// runs), else 1 - film content (FBNS runs write alpha_gas = 0 and carry the void in it).
inline std::vector<double> gas_fraction(const ResultGrid& grid) {
    if (const auto* alpha = grid.field("alpha_gas"))
        if (std::any_of(alpha->begin(), alpha->end(), [](double v) { return std::isfinite(v) && v > 0.0; }))
            return *alpha;
    std::vector<double> out;
    if (const auto* film = grid.field("film_content")) {
        out.resize(film->size());
        for (std::size_t i = 0; i < film->size(); ++i) out[i] = 1.0 - (*film)[i];
    }
    return out;
}

/// Columns of a grid ordered by display angle (ascending).
struct AngleOrder {
    std::vector<int> columns;    // grid column index per display position
    std::vector<double> angles;  // display angle per position [deg]
    double cell_deg = 1.0;       // angular cell size
};

inline AngleOrder angle_order(const ResultGrid& grid, double zero_deg, double start_deg) {
    AngleOrder order;
    const int nx = grid.nx;
    order.columns.resize(static_cast<std::size_t>(std::max(nx, 0)));
    std::iota(order.columns.begin(), order.columns.end(), 0);
    std::vector<double> display(order.columns.size());
    for (int i = 0; i < nx; ++i)
        display[static_cast<std::size_t>(i)] = gui::display_angle(
            grid.theta_deg[static_cast<std::size_t>(i)], zero_deg, start_deg);
    std::sort(order.columns.begin(), order.columns.end(),
              [&](int a, int b) { return display[static_cast<std::size_t>(a)] <
                                         display[static_cast<std::size_t>(b)]; });
    order.angles.resize(order.columns.size());
    for (std::size_t k = 0; k < order.columns.size(); ++k)
        order.angles[k] = display[static_cast<std::size_t>(order.columns[k])];
    order.cell_deg = nx > 0 ? 360.0 / nx : 1.0;
    return order;
}

/// Grid row nearest to an axial fraction z / L in [0, 1].
inline int row_at(const ResultGrid& grid, double z_fraction) {
    if (grid.ny <= 0) return 0;
    return std::clamp(static_cast<int>(std::floor(z_fraction * grid.ny)), 0, grid.ny - 1);
}

/// Values of `field` along grid row j in display order.
inline std::vector<double> row_values(const ResultGrid& grid, const std::vector<double>& field,
                                      int j, const AngleOrder& order) {
    std::vector<double> out(order.columns.size(), std::numeric_limits<double>::quiet_NaN());
    const std::size_t base = static_cast<std::size_t>(j) * static_cast<std::size_t>(grid.nx);
    for (std::size_t k = 0; k < order.columns.size(); ++k) {
        const std::size_t index = base + static_cast<std::size_t>(order.columns[k]);
        if (index < field.size()) out[k] = field[index];
    }
    return out;
}

/// Fixed-frame angle of the chosen zero; moving zeros come from the results of the grid shown.
inline double zero_position(gui::AngleZero zero, double feed_deg, const gui::KeyResults& results) {
    const double h_max = std::isfinite(results.h_max_theta_deg) ? results.h_max_theta_deg : 0.0;
    const double h_min = std::isfinite(results.h_min_theta_deg) ? results.h_min_theta_deg : 0.0;
    return gui::angle_zero_position(zero, feed_deg, h_max, h_min);
}

/// The results loader skips a reload when the new step has the same step number as the grid
/// shown (steady snapshots are all step 0). After a scan this guard loads the step itself.
struct ReloadGuard {
    bool active = false;
    unsigned revision = 0;
    const ResultGrid* grid = nullptr;

    void arm(unsigned loader_revision, const ResultGrid* shown) {
        active = true;
        revision = loader_revision;
        grid = shown;
    }
    /// Call every frame; `load` is called with the step to show once the scan is done.
    template <typename LoadFn>
    void poll(bool busy, unsigned loader_revision, const std::vector<StepInfo>& steps,
              int selected_step, const ResultGrid* shown, LoadFn load) {
        if (!active || busy || loader_revision == revision) return;
        active = false;
        if (shown != grid || steps.empty()) return;  // the loader already brought a new grid
        const int last = static_cast<int>(steps.size()) - 1;
        const int index = selected_step < 0 ? last : std::min(selected_step, last);
        load(steps[static_cast<std::size_t>(index)]);
    }
};

/// Field shown in a map.
enum class MapField { Pressure, Film, Gas, Temperature };

/// A field in display order, ready for a heatmap, with the cavitated-region outline.
struct MapCache {
    const ResultGrid* grid = nullptr;
    MapField field = MapField::Pressure;
    double zero = std::numeric_limits<double>::quiet_NaN(), start = std::numeric_limits<double>::quiet_NaN(), scale = 1.0, offset = 0.0;
    bool available = false;
    std::vector<double> values;  // display order, first row at z = L (heatmap top)
    int rows = 0, cols = 0;
    double vmin = 0.0, vmax = 1.0, x_left = 0.0;
    std::vector<double> outline;  // x0, y0, x1, y1 per cell edge (plot coordinates)
};

/// Values along the angle (display order) at one axial position [SI].
struct Profile {
    std::vector<double> angle, p, h, gas, t;
};

/// The grid values of a map field (gas fraction computed into `gas`).
inline const std::vector<double>* map_source(const ResultGrid& grid, MapField field,
                                      std::vector<double>& gas) {
    switch (field) {
        case MapField::Pressure: return grid.field("p");
        case MapField::Film: return grid.field("h");
        case MapField::Temperature: return grid.field("T");
        case MapField::Gas: break;
    }
    gas = gas_fraction(grid);
    return gas.empty() ? nullptr : &gas;
}

/// Fills `m` for one grid, field, angle zero and display unit (null unit = as stored).
inline void build_map(MapCache& m, const ResultGrid& grid, MapField field, double zero, double start,
               const gui::UnitOption* unit) {
    m.grid = &grid;
    m.field = field;
    m.zero = zero;
    m.start = start;
    m.scale = unit ? unit->scale : 1.0;
    m.offset = unit ? unit->offset : 0.0;
    m.available = false;
    m.values.clear();
    m.outline.clear();
    std::vector<double> gas_store;
    const std::vector<double>* source = map_source(grid, field, gas_store);
    const int nx = grid.nx, ny = grid.ny;
    const std::size_t cells = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    if (!source || nx <= 0 || ny <= 0 || source->size() < cells) return;
    const auto order = angle_order(grid, zero, start);
    m.rows = ny;
    m.cols = nx;
    m.values.resize(cells);
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (int r = 0; r < ny; ++r) {
        const std::size_t j = static_cast<std::size_t>(ny - 1 - r);
        for (int k = 0; k < nx; ++k) {
            double v = (*source)[j * nx + static_cast<std::size_t>(order.columns[k])];
            if (unit) v = gui::to_display(*unit, v);
            m.values[static_cast<std::size_t>(r) * nx + k] = v;
            if (std::isfinite(v)) {
                lo = std::min(lo, v);
                hi = std::max(hi, v);
            }
        }
    }
    if (!std::isfinite(lo)) lo = 0.0, hi = 1.0;
    if (hi <= lo) hi = lo + std::max(std::abs(lo) * 1e-6, 1e-12);
    m.vmin = lo;
    m.vmax = hi;
    m.x_left = order.angles.front() - 0.5 * order.cell_deg;

    // Cavitated-region outline: cell edges where the gas fraction crosses 0.5.
    const std::vector<double> gas = field == MapField::Gas ? gas_store : gas_fraction(grid);
    if (gas.size() >= cells) {
        const double dx = order.cell_deg, dy = 1.0 / ny;
        const auto cavitated = [&](int k, int j) {
            return gas[static_cast<std::size_t>(j) * nx + static_cast<std::size_t>(order.columns[k])] >= 0.5;
        };
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nx; ++k) {
                const bool here = cavitated(k, j);
                if (k + 1 < nx && here != cavitated(k + 1, j)) {
                    const double x = m.x_left + (k + 1) * dx;
                    m.outline.insert(m.outline.end(), {x, j * dy, x, (j + 1) * dy});
                }
                if (j + 1 < ny && here != cavitated(k, j + 1)) {
                    const double y = (j + 1) * dy;
                    m.outline.insert(m.outline.end(), {m.x_left + k * dx, y, m.x_left + (k + 1) * dx, y});
                }
            }
        }
    }
    m.available = true;
}

/// Profiles at axial fraction z / L, in display-angle order.
inline Profile make_profile(const ResultGrid& grid, double zero, double start, double z_fraction) {
    Profile out;
    const auto order = angle_order(grid, zero, start);
    const int j = row_at(grid, z_fraction);
    out.angle = order.angles;
    if (const auto* p = grid.field("p")) out.p = row_values(grid, *p, j, order);
    if (const auto* h = grid.field("h")) out.h = row_values(grid, *h, j, order);
    if (const auto* t = grid.field("T")) out.t = row_values(grid, *t, j, order);
    const auto gas = gas_fraction(grid);
    if (!gas.empty()) out.gas = row_values(grid, gas, j, order);
    return out;
}

/// Values converted to a display unit.
inline std::vector<double> scaled(const std::vector<double>& values, const gui::UnitOption& unit) {
    std::vector<double> out(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) out[i] = gui::to_display(unit, values[i]);
    return out;
}

/// motion.csv of one run, read off the UI thread.
struct MotionCache {
    std::filesystem::path path;
    std::future<std::pair<CsvTable, std::string>> task;
    bool loading = false;
    CsvTable table;
    std::string error;
};

/// Last step of a run and its key results.
struct CompareResult {
    std::shared_ptr<const ResultGrid> grid;
    gui::KeyResults results;
    std::string error;
};

/// Async read of a run for the comparison view.
struct CompareRun {
    std::filesystem::path dir;
    std::future<CompareResult> task;
    bool loading = false;
    CompareResult result;
};

/// Starts (or collects) the async read of a motion.csv; re-reads when `path` changes.
inline void ensure_motion(MotionCache& m, const std::filesystem::path& path) {
    if (m.loading && m.task.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        auto result = m.task.get();
        m.table = std::move(result.first);
        m.error = std::move(result.second);
        m.loading = false;
    }
    if (m.path == path || m.loading) return;
    m.path = path;
    m.table = {};
    m.error.clear();
    m.loading = true;
    m.task = std::async(std::launch::async, [path]() {
        CsvTable table;
        std::string error;
        read_csv(path, table, error);
        return std::make_pair(std::move(table), error);
    });
}

/// Starts (or collects) the async read of the last step of a run folder and its key results.
inline void ensure_compare(CompareRun& run, const std::filesystem::path& output_dir,
                           const SimulationConfig& config) {
    if (run.loading && run.task.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        run.result = run.task.get();
        run.loading = false;
    }
    if (run.dir == output_dir || run.loading) return;
    run.dir = output_dir;
    run.result = {};
    run.loading = true;
    const std::filesystem::path dir = run.dir;
    const std::string prefix =
        config.filename_prefix.empty() ? "solution" : config.filename_prefix;
    const SimulationConfig cfg = config;
    run.task = std::async(std::launch::async, [dir, prefix, cfg]() {
        CompareResult out;
        std::error_code ec;
        std::filesystem::path folder = dir / "flat";
        if (!std::filesystem::exists(folder / "results.pvd", ec)) folder = dir;
        const auto steps = parse_pvd(folder, prefix, out.error);
        if (steps.empty()) {
            if (out.error.empty()) out.error = "no result steps found";
            return out;
        }
        auto grid = parse_step_grid(steps.back(), out.error);
        if (!grid) return out;
        out.results = gui::compute_key_results(*grid, cfg);
        out.grid = std::move(grid);
        return out;
    });
}


}  // namespace result_data
