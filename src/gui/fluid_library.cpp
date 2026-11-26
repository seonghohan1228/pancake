#include "fluid_library.hpp"

#include "property_table_io.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>

namespace fs = std::filesystem;

namespace gui {

namespace {

std::string trim(const std::string& text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::map<std::string, std::string> read_key_values(const fs::path& path, std::string& error) {
    std::map<std::string, std::string> values;
    std::ifstream file(path);
    if (!file) {
        error = "cannot read " + path.string();
        return values;
    }
    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.rfind("\xEF\xBB\xBF", 0) == 0) line = line.substr(3);
        if (line.empty() || line[0] == '#') continue;
        const auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        values[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }
    return values;
}

bool parse_species(const std::string& text, DissolvedGasSpecies& species) {
    if (text == "PROPANE") species = DissolvedGasSpecies::PROPANE;
    else if (text == "AIR") species = DissolvedGasSpecies::AIR;
    else if (text == "R410A") species = DissolvedGasSpecies::R410A;
    else return false;
    return true;
}

/// Solver-side configuration holding only the pair's curve tables, so that
/// saturation_pressure() runs the solver's own inversion.
SimulationConfig table_config(const FluidPair& pair) {
    SimulationConfig cfg;
    cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
    cfg.solubility_model = SolubilityModel::CURVE_TABLE;
    cfg.solubility_curve_table = pair.solubility;
    cfg.kinematic_composition_table = pair.viscosity;
    cfg.solution_table_clamp = true;
    return cfg;
}

}  // namespace

const char* species_name(DissolvedGasSpecies species) {
    switch (species) {
        case DissolvedGasSpecies::PROPANE: return "PROPANE";
        case DissolvedGasSpecies::AIR: return "AIR";
        case DissolvedGasSpecies::R410A: return "R410A";
    }
    return "PROPANE";
}

double FluidPair::density(double temperature_k) const {
    return density_intercept + density_slope * (temperature_k - 273.15);
}

PropertyValue FluidPair::kinematic_viscosity(double content, double temperature_k) const {
    PropertyValue out;
    if (viscosity.empty()) {
        out.value = std::numeric_limits<double>::quiet_NaN();
        out.in_range = false;
        return out;
    }
    out.in_range = content >= viscosity.axis.front() && content <= viscosity.axis.back();
    double t_min = 0.0, t_max = 0.0;
    temperature_range(content, t_min, t_max);
    out.in_range = out.in_range && temperature_k >= t_min && temperature_k <= t_max;
    out.value = viscosity.interpolate(content, temperature_k, true, true);
    return out;
}

PropertyValue FluidPair::dynamic_viscosity(double content, double temperature_k) const {
    PropertyValue out = kinematic_viscosity(content, temperature_k);
    out.value *= density(temperature_k);
    return out;
}

PropertyValue FluidPair::saturation_pressure(double content, double temperature_k) const {
    PropertyValue out;
    if (solubility.empty() || content <= 0.0) {
        out.value = 0.0;
        out.in_range = !solubility.empty();
        return out;
    }
    const SimulationConfig cfg = table_config(*this);
    out.value = cfg.saturation_pressure(content, temperature_k);
    // In range when the content is reached on the isotherm, i.e. the inverted
    // pressure lies strictly inside the table pressure axis or hits a node.
    std::vector<double> c_iso;
    for (std::size_t i = 0; i < solubility.axis.size(); ++i)
        c_iso.push_back(solubility.isotherm_at(i, temperature_k, true));
    const double c_lo = *std::min_element(c_iso.begin(), c_iso.end());
    const double c_hi = *std::max_element(c_iso.begin(), c_iso.end());
    double t_lo = std::numeric_limits<double>::infinity(), t_hi = 0.0;
    for (const auto& curve : solubility.curves) {
        t_lo = std::min(t_lo, curve.front().x);
        t_hi = std::max(t_hi, curve.back().x);
    }
    out.in_range = content >= c_lo && content <= c_hi && temperature_k >= t_lo &&
                   temperature_k <= t_hi;
    return out;
}

PropertyValue FluidPair::saturated_content(double pressure_pa, double temperature_k) const {
    PropertyValue out;
    if (solubility.empty()) {
        out.value = std::numeric_limits<double>::quiet_NaN();
        out.in_range = false;
        return out;
    }
    out.value = solubility.interpolate(pressure_pa, temperature_k, true);
    out.in_range = pressure_pa >= solubility.axis.front() && pressure_pa <= solubility.axis.back();
    if (out.in_range) {
        try {
            (void)solubility.interpolate(pressure_pa, temperature_k, false);
        } catch (const std::exception&) {
            out.in_range = false;
        }
    }
    return out;
}

double FluidPair::content_min() const {
    return viscosity.empty() ? 0.0 : viscosity.axis.front();
}

double FluidPair::content_max() const {
    return viscosity.empty() ? 0.0 : viscosity.axis.back();
}

void FluidPair::temperature_range(double content, double& t_min, double& t_max) const {
    t_min = t_max = std::numeric_limits<double>::quiet_NaN();
    if (viscosity.empty()) return;
    // The curves bracketing `content` both have to cover the query temperature.
    std::size_t lo = 0, hi = 0;
    for (std::size_t i = 0; i < viscosity.axis.size(); ++i) {
        if (viscosity.axis[i] <= content) lo = i;
    }
    hi = std::min(lo + (viscosity.axis[lo] < content ? 1 : 0), viscosity.axis.size() - 1);
    t_min = std::max(viscosity.curves[lo].front().x, viscosity.curves[hi].front().x);
    t_max = std::min(viscosity.curves[lo].back().x, viscosity.curves[hi].back().x);
}

bool load_fluid_pair(const fs::path& folder, FluidPair& out, std::string& error) {
    error.clear();
    FluidPair pair;
    pair.id = folder.filename().string();
    const auto values = read_key_values(folder / "pair.txt", error);
    if (!error.empty()) return false;
    const auto text = [&](const char* key, const std::string& fallback = {}) {
        const auto it = values.find(key);
        return it == values.end() ? fallback : it->second;
    };
    const auto number = [&](const char* key, double fallback) {
        const auto it = values.find(key);
        if (it == values.end()) return fallback;
        try {
            return std::stod(it->second);
        } catch (const std::exception&) {
            error = std::string("pair.txt: ") + key + " is not a number";
            return fallback;
        }
    };
    pair.name = text("name", pair.id);
    pair.oil = text("oil");
    pair.refrigerant = text("refrigerant");
    pair.source = text("source");
    pair.note = text("note");
    // gas_model (rev. 3); solver_species is the earlier name of the same key.
    const std::string gas_model =
        !text("gas_model").empty() ? text("gas_model") : text("solver_species");
    pair.has_gas_phase = !gas_model.empty();
    if (pair.has_gas_phase && !parse_species(gas_model, pair.solver_species)) {
        error = "pair.txt: gas_model must be R410A, PROPANE or AIR";
        return false;
    }
    pair.density_intercept = number("density_intercept_kg_m3", pair.density_intercept);
    pair.density_slope = number("density_slope_kg_m3_K", pair.density_slope);
    pair.oil_vapour_pressure = number("oil_vapour_pressure_Pa", pair.oil_vapour_pressure);
    pair.specific_heat = number("specific_heat_J_kgK", pair.specific_heat);
    pair.thermal_conductivity = number("thermal_conductivity_W_mK", pair.thermal_conductivity);
    pair.dissolved_heat_capacity = number("dissolved_heat_capacity_J_kgK", pair.dissolved_heat_capacity);
    pair.desorption_energy = number("desorption_energy_J_kg", pair.desorption_energy);
    pair.gas_heat_capacity_cv = number("gas_heat_capacity_cv_J_kgK", pair.gas_heat_capacity_cv);
    pair.gas_thermal_conductivity = number("gas_thermal_conductivity_W_mK", pair.gas_thermal_conductivity);
    pair.gas_viscosity = number("gas_viscosity_Pa_s", pair.gas_viscosity);
    if (!error.empty()) return false;
    pair.table_file = fs::absolute(folder / text("table", "solution.csv"));
    FluidProperties::load_solution_table_csv(pair.table_file, pair.solubility, pair.viscosity,
                                             error);
    if (!error.empty()) return false;
    if (pair.solubility.empty()) {
        error = "solution table has no saturation curves (content > 0 rows)";
        return false;
    }
    out = std::move(pair);
    return true;
}

std::vector<std::string> FluidLibrary::list() const {
    std::vector<std::string> ids;
    std::error_code ec;
    if (root_.empty() || !fs::is_directory(root_, ec)) return ids;
    for (const auto& entry : fs::directory_iterator(root_, ec)) {
        if (entry.is_directory() && fs::exists(entry.path() / "pair.txt"))
            ids.push_back(entry.path().filename().string());
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

bool FluidLibrary::load(const std::string& id, FluidPair& out, std::string& error) const {
    return load_fluid_pair(root_ / id, out, error);
}

std::string FluidLibrary::find_by_table(const fs::path& table) const {
    if (table.empty()) return "";
    std::error_code ec;
    for (const std::string& id : list()) {
        const fs::path candidate = root_ / id / "solution.csv";
        if (fs::equivalent(candidate, table, ec)) return id;
    }
    return "";
}

bool read_solution_rows(const fs::path& csv, std::vector<SolutionRow>& rows, std::string& error) {
    rows.clear();
    error.clear();
    std::ifstream file(csv);
    if (!file) {
        error = "cannot read " + csv.string();
        return false;
    }
    std::string line;
    std::getline(file, line);  // header (column order fixed by write_solution_rows)
    std::size_t number = 1;
    while (std::getline(file, line)) {
        ++number;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cells;
        std::stringstream stream(line);
        std::string cell;
        while (std::getline(stream, cell, ',')) cells.push_back(trim(cell));
        if (line.back() == ',') cells.emplace_back();
        if (cells.size() < 4) {
            error = "line " + std::to_string(number) + ": expected 4 columns";
            return false;
        }
        try {
            SolutionRow row;
            row.content = std::stod(cells[0]);
            row.temperature_c = std::stod(cells[1]);
            row.pressure_mpa = cells[2].empty() ? std::numeric_limits<double>::quiet_NaN()
                                                : std::stod(cells[2]);
            row.viscosity_mm2s = std::stod(cells[3]);
            rows.push_back(row);
        } catch (const std::exception&) {
            error = "line " + std::to_string(number) + ": not a number";
            return false;
        }
    }
    return true;
}

bool write_solution_rows(const fs::path& csv, const std::vector<SolutionRow>& rows,
                         std::string& error) {
    std::ofstream file(csv, std::ios::binary);
    if (!file) {
        error = "cannot write " + csv.string();
        return false;
    }
    file << "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n" << std::setprecision(10);
    for (const SolutionRow& row : rows) {
        file << row.content << ',' << row.temperature_c << ',';
        if (std::isfinite(row.pressure_mpa) && row.content > 0.0) file << row.pressure_mpa;
        file << ',' << row.viscosity_mm2s << '\n';
    }
    return static_cast<bool>(file);
}

bool parse_pasted_rows(const std::string& text, std::vector<SolutionRow>& rows,
                       std::string& error) {
    rows.clear();
    error.clear();
    std::stringstream stream(text);
    std::string line;
    std::size_t number = 0;
    while (std::getline(stream, line)) {
        ++number;
        line = trim(line);
        if (line.empty()) continue;
        std::vector<std::string> cells;
        std::string cell;
        for (char ch : line) {
            if (ch == '\t' || ch == ',' || ch == ';') {
                cells.push_back(trim(cell));
                cell.clear();
            } else {
                cell.push_back(ch);
            }
        }
        cells.push_back(trim(cell));
        const auto is_number = [](const std::string& s) {
            if (s.empty()) return false;
            char* end = nullptr;
            std::strtod(s.c_str(), &end);
            return end && *end == '\0';
        };
        if (cells.size() < 4 || !is_number(cells[0])) {
            if (number == 1) continue;  // header line
            error = "row " + std::to_string(number) +
                    ": expected content, temperature, saturation pressure, viscosity";
            return false;
        }
        SolutionRow row;
        row.content = std::stod(cells[0]);
        row.temperature_c = std::stod(cells[1]);
        row.pressure_mpa = is_number(cells[2]) ? std::stod(cells[2])
                                               : std::numeric_limits<double>::quiet_NaN();
        if (!is_number(cells[3])) {
            error = "row " + std::to_string(number) + ": viscosity is not a number";
            return false;
        }
        row.viscosity_mm2s = std::stod(cells[3]);
        rows.push_back(row);
    }
    if (rows.empty()) error = "no data rows found";
    return error.empty();
}

bool FluidLibrary::save(const FluidPair& pair, const std::vector<SolutionRow>& rows,
                        std::string& error) const {
    error.clear();
    const fs::path folder = root_ / pair.id;
    std::error_code ec;
    fs::create_directories(folder, ec);
    if (ec) {
        error = "cannot create " + folder.string();
        return false;
    }
    std::ofstream header(folder / "pair.txt", std::ios::binary);
    if (!header) {
        error = "cannot write pair.txt";
        return false;
    }
    header << std::setprecision(10)
           << "# Oil-refrigerant pair for the pancake GUI. Format: fluids/README.md\n"
           << "name = " << pair.name << "\n"
           << "oil = " << pair.oil << "\n"
           << "refrigerant = " << pair.refrigerant << "\n"
           << "source = " << pair.source << "\n"
           << "note = " << pair.note << "\n"
           << (pair.has_gas_phase
                   ? std::string("gas_model = ") + species_name(pair.solver_species) + "\n"
                   : std::string("# no gas_model: FBNS only (no gas transport)\n"))
           << "density_intercept_kg_m3 = " << pair.density_intercept << "\n"
           << "density_slope_kg_m3_K = " << pair.density_slope << "\n"
           << "oil_vapour_pressure_Pa = " << pair.oil_vapour_pressure << "\n"
           << "specific_heat_J_kgK = " << pair.specific_heat << "\n"
           << "thermal_conductivity_W_mK = " << pair.thermal_conductivity << "\n"
           << "table = solution.csv\n";
    if (pair.gas_viscosity > 0.0) header << "gas_viscosity_Pa_s = " << pair.gas_viscosity << "\n";
    if (pair.has_gas_phase) {
        if (pair.gas_heat_capacity_cv > 0.0)
            header << "gas_heat_capacity_cv_J_kgK = " << pair.gas_heat_capacity_cv << "\n";
        if (pair.gas_thermal_conductivity > 0.0)
            header << "gas_thermal_conductivity_W_mK = " << pair.gas_thermal_conductivity << "\n";
    }
    if (pair.has_caloric_data()) {
        header << "dissolved_heat_capacity_J_kgK = " << pair.dissolved_heat_capacity << "\n"
               << "desorption_energy_J_kg = " << pair.desorption_energy << "\n";
    }
    if (!header) {
        error = "cannot write pair.txt";
        return false;
    }
    header.close();
    if (!write_solution_rows(folder / "solution.csv", rows, error)) return false;
    // Load back with the solver's loader: a table it rejects is reported now.
    FluidPair check;
    return load_fluid_pair(folder, check, error);
}

bool FluidLibrary::copy(const std::string& from_id, const std::string& to_id,
                        std::string& error) const {
    std::error_code ec;
    const fs::path target = root_ / to_id;
    if (fs::exists(target, ec)) {
        error = "a pair named '" + to_id + "' already exists";
        return false;
    }
    fs::copy(root_ / from_id, target, fs::copy_options::recursive, ec);
    if (ec) {
        error = "copy failed: " + ec.message();
        return false;
    }
    return true;
}

}  // namespace gui
