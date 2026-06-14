#include "property_table_io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>

namespace FluidProperties {
namespace {

std::string strip_bom_and_cr(std::string line) {
    if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF &&
        static_cast<unsigned char>(line[1]) == 0xBB &&
        static_cast<unsigned char>(line[2]) == 0xBF) {
        line.erase(0, 3);  // UTF-8 BOM
    }
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return line;
}

std::string trim_copy(std::string s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return {};
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) out.push_back(trim_copy(cell));
    return out;
}

// Index of value v within a small axis vector (tolerant compare); -1 if absent.
long axis_index(const std::vector<double>& axis, double v) {
    for (std::size_t i = 0; i < axis.size(); ++i) {
        const double tol = 1.0e-6 * std::max(1.0, std::fabs(axis[i]));
        if (std::fabs(axis[i] - v) <= tol) return static_cast<long>(i);
    }
    return -1;
}

}  // namespace

PropertyTable2D load_table_2d_csv(const std::filesystem::path& path,
                                  const std::string& value_header,
                                  double value_scale,
                                  std::string& error_out) {
    error_out.clear();
    PropertyTable2D table;
    const std::string where = "'" + path.string() + "'";

    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        error_out = "cannot open " + where;
        return table;
    }

    std::string raw;
    if (!std::getline(file, raw)) {
        error_out = where + " is empty";
        return table;
    }
    const std::vector<std::string> header = split_csv(strip_bom_and_cr(raw));
    auto find_col = [&](const std::string& name) -> long {
        const std::string want = to_lower(name);
        for (std::size_t i = 0; i < header.size(); ++i)
            if (to_lower(header[i]) == want) return static_cast<long>(i);
        return -1;
    };
    const long col_p = find_col("P_MPa");
    const long col_t = find_col("T_C");
    const long col_v = find_col(value_header);
    if (col_p < 0 || col_t < 0 || col_v < 0) {
        error_out = where + " is missing a required header (need P_MPa, T_C, " + value_header + ")";
        return table;
    }
    const std::size_t need = static_cast<std::size_t>(std::max({col_p, col_t, col_v})) + 1;

    struct Row { double p, t, v; };
    std::vector<Row> rows;
    std::vector<double> p_vals, t_vals;
    long line_no = 1;
    while (std::getline(file, raw)) {
        ++line_no;
        const std::string line = strip_bom_and_cr(raw);
        const std::string trimmed = trim_copy(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        const std::vector<std::string> cells = split_csv(line);
        if (cells.size() < need) {
            error_out = where + " line " + std::to_string(line_no) + ": too few columns";
            return table;
        }
        Row r{};
        try {
            r.p = std::stod(cells[static_cast<std::size_t>(col_p)]) * 1.0e6;     // MPa -> Pa
            r.t = std::stod(cells[static_cast<std::size_t>(col_t)]) + 273.15;    // C  -> K
            r.v = std::stod(cells[static_cast<std::size_t>(col_v)]) * value_scale;  // -> SI
        } catch (const std::exception&) {
            error_out = where + " line " + std::to_string(line_no) + ": non-numeric value";
            return table;
        }
        if (!std::isfinite(r.p) || !std::isfinite(r.t) || !std::isfinite(r.v)) {
            error_out = where + " line " + std::to_string(line_no) + ": non-finite value";
            return table;
        }
        rows.push_back(r);
        if (axis_index(p_vals, r.p) < 0) p_vals.push_back(r.p);
        if (axis_index(t_vals, r.t) < 0) t_vals.push_back(r.t);
    }

    if (rows.empty()) {
        error_out = where + " has no data rows";
        return table;
    }
    std::sort(p_vals.begin(), p_vals.end());
    std::sort(t_vals.begin(), t_vals.end());
    if (p_vals.size() < 2 || t_vals.size() < 2) {
        error_out = where + " needs >= 2 distinct pressures and temperatures";
        return table;
    }
    const std::size_t np = p_vals.size();
    const std::size_t nt = t_vals.size();
    if (rows.size() != np * nt) {
        error_out = where + " is not a complete rectangular grid (" + std::to_string(rows.size()) +
                    " rows, expected " + std::to_string(np * nt) + " = " + std::to_string(np) +
                    "x" + std::to_string(nt) + ")";
        return table;
    }

    std::vector<double> grid(np * nt, 0.0);
    std::vector<char> filled(np * nt, 0);
    for (const Row& r : rows) {
        const long ip = axis_index(p_vals, r.p);
        const long it = axis_index(t_vals, r.t);
        if (ip < 0 || it < 0) {
            error_out = where + ": internal grid placement error";
            return table;
        }
        const std::size_t idx = static_cast<std::size_t>(ip) * nt + static_cast<std::size_t>(it);
        if (filled[idx]) {
            error_out = where + ": duplicate (P,T) entry";
            return table;
        }
        grid[idx] = r.v;
        filled[idx] = 1;
    }
    if (std::find(filled.begin(), filled.end(), 0) != filled.end()) {
        error_out = where + ": incomplete grid (a (P,T) node has no value)";
        return table;
    }

    table.p = std::move(p_vals);
    table.t = std::move(t_vals);
    table.values = std::move(grid);
    return table;
}

void load_config_property_tables(SimulationConfig& cfg, std::vector<std::string>& errors) {
    auto resolve = [&](const std::filesystem::path& f) -> std::filesystem::path {
        if (f.empty() || f.is_absolute()) return f;
        if (!cfg.config_dir.empty()) {
            std::error_code ec;
            const std::filesystem::path cand = cfg.config_dir / f;
            if (std::filesystem::exists(cand, ec)) return cand;
        }
        return f;  // CWD-relative fallback
    };
    auto load_one = [&](const std::filesystem::path& f, const char* header, double scale,
                        PropertyTable2D& dst, const char* label) {
        if (f.empty()) return;
        std::string err;
        PropertyTable2D loaded = load_table_2d_csv(resolve(f), header, scale, err);
        if (!err.empty()) {
            errors.push_back(std::string(label) + "_table_file: " + err);
            return;
        }
        dst = std::move(loaded);
    };
    load_one(cfg.solubility_table_file, "solubility_pct", 0.01, cfg.solubility_table_2d, "solubility");
    load_one(cfg.viscosity_table_file, "viscosity_mm2s", 1.0e-6, cfg.viscosity_table_2d, "viscosity");
    load_one(cfg.density_table_file, "density_kgm3", 1.0, cfg.density_table_2d, "density");
}

}  // namespace FluidProperties
