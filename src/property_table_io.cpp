#include "property_table_io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <fstream>
#include <map>
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

static PropertyTable2D load_rectangular_csv(const std::filesystem::path& path,
                                  const std::string& value_header,
                                  double value_scale,
                                  std::string& error_out,
                                  const std::string& axis_header, double axis_scale) {
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
    const long col_p = find_col(axis_header);
    const long col_t = find_col("T_C");
    const long col_v = find_col(value_header);
    if (col_p < 0 || col_t < 0 || col_v < 0) {
        error_out = where + " is missing a required header (need " + axis_header + ", T_C, " + value_header + ")";
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
            r.p = std::stod(cells[static_cast<std::size_t>(col_p)]) * axis_scale;
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
        error_out = where + " needs >= 2 distinct first-axis values and temperatures";
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
            error_out = where + ": duplicate (" + axis_header + ",T) entry";
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

PropertyTable2D load_table_2d_csv(const std::filesystem::path& path,const std::string& value_header,
                               double value_scale,std::string& error_out) {
    return load_rectangular_csv(path,value_header,value_scale,error_out,"P_MPa",1e6);
}
CompositionViscosityTable load_viscosity_composition_csv(const std::filesystem::path& path,std::string& error_out) {
    auto grid=load_rectangular_csv(path,"mu_Pa_s",1,error_out,"c_mass_fraction",1);
    if (!error_out.empty()) return {};
    if (grid.p.front()<0 || grid.p.back()>=1 || grid.t.front()<=0 ||
        std::any_of(grid.values.begin(),grid.values.end(),[](double mu) { return mu<=0; })) {
        error_out="Composition viscosity requires 0 <= c < 1, T > 0 K and mu > 0 Pa s"; return {};
    }
    for (double& value : grid.values) value=std::log(value);
    return {std::move(grid.p),std::move(grid.t),std::move(grid.values)};
}

void load_solution_table_csv(const std::filesystem::path& path, PropertyCurveTable& solubility,
                             PropertyCurveTable& viscosity, std::string& error_out) {
    solubility = {};
    viscosity = {};
    error_out.clear();
    std::map<double, std::vector<PropertyTablePoint>> pressure_curves, concentration_curves;
    std::ifstream file(path);
    std::string line;
    std::size_t line_number = 1;
    try {
        if (!file || !std::getline(file, line)) throw std::runtime_error("cannot read table");
        auto header = split_csv(strip_bom_and_cr(line));
        const auto column = [&](const char* name) {
            for (std::size_t i = 0; i < header.size(); ++i)
                if (to_lower(header[i]) == to_lower(name)) return i;
            throw std::runtime_error(std::string("missing column ") + name);
        };
        const auto c_column = column("c_mass_fraction"), t_column = column("T_C");
        const auto p_column = column("P_MPa"), nu_column = column("viscosity_mm2s");
        const auto needed = std::max({c_column, t_column, p_column, nu_column}) + 1;
        while (std::getline(file, line)) {
            ++line_number;
            line = trim_copy(strip_bom_and_cr(line));
            if (line.empty() || line[0] == '#') continue;
            auto cells = split_csv(line);
            // A reordered blank oil-pressure column may be the last cell.
            if (line.back() == ',') cells.emplace_back();
            if (cells.size() < needed) throw std::runtime_error("too few columns");
            const auto number = [&](std::size_t index) {
                std::size_t end;
                const double value = std::stod(cells[index], &end);
                if (end != cells[index].size() || !std::isfinite(value))
                    throw std::runtime_error("invalid numeric value");
                return value;
            };
            const double c = number(c_column), t = number(t_column) + 273.15;
            const double nu = number(nu_column) * 1e-6;
            if (c < 0 || c > 1 || t <= 0 || nu <= 0 || !std::isfinite(t) || !std::isfinite(nu))
                throw std::runtime_error("require 0 <= c <= 1, T > 0 K, nu > 0");
            concentration_curves[c].push_back({t, nu});
            if (c == 0) {
                if (!cells[p_column].empty()) throw std::runtime_error("pure-oil rows need blank pressure");
            } else {
                const double pressure = number(p_column) * 1e6;
                if (pressure <= 0 || !std::isfinite(pressure)) throw std::runtime_error("require positive absolute pressure");
                pressure_curves[pressure].push_back({t, c});
            }
        }
        const auto finish = [](auto& groups, PropertyCurveTable& table, bool saturation) {
            if (groups.size() < 2) throw std::runtime_error("need at least two pressure/concentration curves");
            for (auto& [axis, points] : groups) {
                std::sort(points.begin(), points.end(), [](auto a, auto b) { return a.x < b.x; });
                if (points.size() < 2) throw std::runtime_error("each curve needs at least two temperatures");
                for (std::size_t i = 1; i < points.size(); ++i) {
                    if (points[i].x <= points[i - 1].x) throw std::runtime_error("duplicate temperature in curve");
                    if (saturation && points[i].value >= points[i - 1].value)
                        throw std::runtime_error("solubility must decrease with T along each isobar");
                }
                table.axis.push_back(axis);
                table.curves.push_back(std::move(points));
            }
        };
        finish(pressure_curves, solubility, true);
        finish(concentration_curves, viscosity, false);
    } catch (const std::exception& error) {
        solubility = {};
        viscosity = {};
        error_out = path.string() + " line " + std::to_string(line_number) + ": " + error.what();
    }
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
    cfg.solubility_curve_table = {};
    cfg.kinematic_composition_table = {};
    if (!cfg.solution_table_file.empty()) {
        std::string error;
        load_solution_table_csv(resolve(cfg.solution_table_file), cfg.solubility_curve_table,
                                cfg.kinematic_composition_table, error);
        if (!error.empty()) errors.push_back("solution_table_file: " + error);
    }
    if (!cfg.viscosity_composition_file.empty()) {
        std::string error;
        cfg.viscosity_composition_table=load_viscosity_composition_csv(resolve(cfg.viscosity_composition_file),error);
        if (!error.empty()) errors.push_back("viscosity_composition_file: "+error);
    }
}

}  // namespace FluidProperties
