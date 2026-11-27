#pragma once

// Oil-refrigerant pairs stored as plain files in the fluid directory:
//
//   <fluids>/<pair>/pair.txt      key = value header (names, source, density law, ...)
//   <fluids>/<pair>/solution.csv  c_mass_fraction,T_C,P_MPa,viscosity_mm2s
//
// solution.csv is the solver's own `solution_table_file` format: each content
// curve lists the temperature at which that content is saturated at P_MPa,
// with the kinematic viscosity there. Pure-oil rows (content 0) leave P_MPa
// blank. One table therefore gives both views the GUI shows: viscosity vs
// temperature per content, and saturation pressure vs content. The tables are
// read with the solver's loader, and saturation pressure is inverted with the
// solver's own SimulationConfig::saturation_pressure, so GUI estimates and
// solver runs cannot disagree. No toolkit dependency: unit-tested on its own.

#include "config.hpp"

#include <cmath>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace gui {

/// Range check result of one property query. `in_range` is false when the
/// query left the digitized table; the value is then the clamped end value.
struct PropertyValue {
    double value = 0.0;
    bool in_range = true;
};

/// One oil-refrigerant pair (header + loaded tables).
struct FluidPair {
    // Header (pair.txt)
    std::string id;             // folder name
    std::string name;           // display name
    std::string oil;
    std::string refrigerant;
    std::string source;         // citation or "synthetic"
    std::string note;           // e.g. "example only, not for design"
    // Optional gas-phase block: the solver's gas data for released refrigerant
    // (R410A real-fluid EOS; PROPANE, AIR ideal gas). Without it the pair runs
    // FBNS only (the pressure plateau needs only the saturation pressure).
    bool has_gas_phase = false;
    DissolvedGasSpecies solver_species = DissolvedGasSpecies::PROPANE;  // valid if has_gas_phase
    double density_intercept = 954.9;     // kg/m^3 at 0 degC
    double density_slope = -0.7284;       // kg/(m^3 K)
    double oil_vapour_pressure = 500.0;   // Pa (vaporous cavitation pressure)
    double specific_heat = 2010.0;        // J/(kg K)
    double thermal_conductivity = 0.13;   // W/(m K)
    std::filesystem::path table_file;     // absolute path of solution.csv
    // Caloric data for thermal gas release (COUPLED); optional.
    double dissolved_heat_capacity = -1.0;  // J/(kg K) of dissolved refrigerant; <= 0 unknown
    double desorption_energy = std::numeric_limits<double>::quiet_NaN();  // J/kg, u_gas - u_dissolved
    double gas_heat_capacity_cv = -1.0;     // J/(kg K) free gas; <= 0 from the species model
    double gas_thermal_conductivity = 0.0;  // W/(m K) free gas
    double gas_viscosity = 0.0;             // Pa s free gas (gas release with transport)
    bool has_caloric_data() const {
        return dissolved_heat_capacity > 0.0 && std::isfinite(desorption_energy);
    }

    // Loaded tables (SI: c in kg/kg, T in K, p in Pa, nu in m^2/s)
    PropertyCurveTable solubility;  // c_sat(p, T)
    PropertyCurveTable viscosity;   // nu(c, T)

    bool loaded() const { return !solubility.empty() && !viscosity.empty(); }

    /// Liquid density, linear in temperature [kg/m^3].
    double density(double temperature_k) const;
    /// Kinematic viscosity [m^2/s] at content c [kg/kg] and temperature T [K],
    /// interpolated in log(nu) as the solver does by default.
    PropertyValue kinematic_viscosity(double content, double temperature_k) const;
    /// Dynamic viscosity [Pa s] = density x kinematic viscosity.
    PropertyValue dynamic_viscosity(double content, double temperature_k) const;
    /// Saturation (bubble-point) pressure [Pa] of content c at T.
    PropertyValue saturation_pressure(double content, double temperature_k) const;
    /// Content [kg/kg] in equilibrium with refrigerant at pressure p and T.
    PropertyValue saturated_content(double pressure_pa, double temperature_k) const;

    /// Content range of the table [kg/kg] (curves with a saturation pressure).
    double content_min() const;
    double content_max() const;
    /// Temperature range [K] of the curve nearest to `content`.
    void temperature_range(double content, double& t_min, double& t_max) const;
};

/// Rows of solution.csv in display units, for the table editor.
struct SolutionRow {
    double content = 0.0;        // kg/kg
    double temperature_c = 0.0;  // degC
    double pressure_mpa = 0.0;   // MPa abs; NaN for pure oil
    double viscosity_mm2s = 0.0; // mm^2/s
};

/// Fluid directory: list, load, save, copy pairs.
class FluidLibrary {
public:
    explicit FluidLibrary(std::filesystem::path root = {}) : root_(std::move(root)) {}

    const std::filesystem::path& root() const { return root_; }
    void set_root(std::filesystem::path root) { root_ = std::move(root); }

    /// Pair folder names (those containing pair.txt), sorted.
    std::vector<std::string> list() const;
    /// Loads one pair; on failure returns false and fills `error`.
    bool load(const std::string& id, FluidPair& out, std::string& error) const;
    /// Finds the pair whose table file has the same file content name as `table`
    /// (matched by folder/file name); empty if none.
    std::string find_by_table(const std::filesystem::path& table) const;

    /// Writes pair.txt and solution.csv into <root>/<pair.id>/.
    bool save(const FluidPair& pair, const std::vector<SolutionRow>& rows,
              std::string& error) const;
    /// Copies a pair folder under a new id.
    bool copy(const std::string& from_id, const std::string& to_id, std::string& error) const;

private:
    std::filesystem::path root_;
};

/// Loads a pair header (pair.txt) and its table; standalone for tests and the dialog.
bool load_fluid_pair(const std::filesystem::path& folder, FluidPair& out, std::string& error);
/// Reads solution.csv rows for the table editor.
bool read_solution_rows(const std::filesystem::path& csv, std::vector<SolutionRow>& rows,
                        std::string& error);
/// Writes solution.csv rows (header included).
bool write_solution_rows(const std::filesystem::path& csv, const std::vector<SolutionRow>& rows,
                         std::string& error);
/// Parses tab- or comma-separated text pasted from a spreadsheet into rows.
/// Accepts an optional header line; returns false with a message on a bad row.
bool parse_pasted_rows(const std::string& text, std::vector<SolutionRow>& rows,
                       std::string& error);

const char* species_name(DissolvedGasSpecies species);

}  // namespace gui
