#pragma once

// Engineer-level case model and its two mappings to SimulationConfig.
//
// Panels edit a CaseModel only. `to_config` writes the solver keys derived from
// the case onto a base config; `from_config` reads a case back from any config.
// A loaded config keeps every value the GUI does not show: `to_config` rewrites
// only the key groups whose case inputs changed relative to the case that was
// read from the file (`original`). For a new case (no original) every group is
// written. Model choices follow docs/GUI_LGE_PLAN.md section 4 and the deck
// runs (FBNS plateau for dissolved gas, COUPLED for gas transport).
//
// GUI-only settings (pair id, angle zero, mesh level label, cores, ...) live in
// `#@gui key = value` comment lines; the solver ignores comments, so the config
// file remains the only case format. No toolkit dependency: unit-tested.

#include "config.hpp"
#include "fluid_library.hpp"

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace gui {

enum class AnalysisType { Steady, Transient };
enum class ShaftPosition { LoadBalance, Fixed };
enum class GasModel { Plateau, Transport };  // FBNS + LOCAL_PSAT | COUPLED
enum class FeedType { Hole, Groove };
enum class MeshLevel { Coarse, Medium, Fine, Custom };
enum class LoadAbscissa { CrankAngle, Time };
enum class LoadFormat { MagnitudeDirection, Components };

/// Cells around x along the film for a mesh level (Custom returns {0, 0}).
void mesh_level_cells(MeshLevel level, int& n_theta, int& n_z);
const char* mesh_level_name(MeshLevel level);

struct FeedInput {
    FeedType type = FeedType::Hole;
    double angle_deg = 90.0;                 // centre, fixed frame, direction of rotation
    double hole_diameter = 2.0e-3;           // m
    double hole_axial_position = -1.0;       // m from z = 0; < 0 means mid-length
    double groove_width_deg = 10.0;
    double supply_overpressure = 0.3e5;      // Pa above ambient
    double supply_temperature = -1.0;        // K; < 0 means ambient
    bool turns_with_shaft = false;           // transient FBNS only
    bool operator==(const FeedInput&) const = default;
};

struct LoadTableInput {
    std::filesystem::path file;
    LoadAbscissa abscissa = LoadAbscissa::CrankAngle;
    LoadFormat format = LoadFormat::MagnitudeDirection;
    double frequency_hz = 0.0;     // crank-angle tables: cycles per second
    double phase_offset_deg = 0.0;
    double scale = 1.0;
    bool operator==(const LoadTableInput&) const = default;
};

struct HousingInput {
    double outer_diameter = 0.0;           // m; <= 0 means 2.2 R (solver default 1.1 R)
    double conductivity = 15.0;            // W/(m K)
    double outer_heat_transfer = 40.0;     // W/(m^2 K)
    double shaft_heat_transfer = 500.0;    // W/(m^2 K)
    double shaft_area_ratio = 1.0;
    bool operator==(const HousingInput&) const = default;
};

struct ConvergenceInput {
    double force_tolerance = 1e-3;     // load balance |F + W| / |W|
    double outer_tolerance = 1e-6;     // pressure / temperature change per outer iteration
    int outer_max_iterations = 200;
    double inner_tolerance = 1e-10;    // film equations (FBNS residual)
    int inner_max_iterations = 200;
    double transport_tolerance = 1e-10;  // COUPLED nonlinear tolerance
    double motion_tolerance = 1e-6;      // transient load-position coupling
    int motion_max_iterations = 30;
    bool operator==(const ConvergenceInput&) const = default;
};

/// Engineer-level inputs (SI units).
struct CaseModel {
    // Geometry
    double bearing_diameter = 0.021;   // D = 2R [m]
    double bearing_length = 0.020;     // L [m]
    double radial_clearance = 30e-6;   // c [m]
    FeedInput feed;
    double tilt_angle = 0.0;           // gamma [rad] (slope dh/dz)
    double tilt_direction_deg = 0.0;

    // Operating conditions
    AnalysisType analysis = AnalysisType::Steady;
    ShaftPosition position = ShaftPosition::LoadBalance;
    double speed = 816.8140899333462;  // rad/s
    double load = 0.0;                 // N (steady)
    double load_direction_deg = 90.0;  // fixed frame, as the feed angle
    LoadTableInput load_table;         // transient
    double moving_mass = 1.0;          // kg (transient)
    double fixed_eccentricity_ratio = 0.5;  // fixed shaft position mode
    double fixed_attitude_deg = 0.0;
    double ambient_pressure = 101325.0;     // Pa abs
    double ambient_temperature = 313.15;    // K

    // Fluid
    std::string fluid_pair;           // id in the fluid directory
    bool saturated = true;            // content = saturated at ambient p, T
    double content = 0.0;             // kg/kg when not saturated
    double roughness_shaft = 0.0;     // RMS [m]
    double roughness_bush = 0.0;

    // Model
    bool dissolved_gas = true;
    GasModel gas_model = GasModel::Plateau;
    bool thermal = false;
    HousingInput housing;

    // Mesh and time
    MeshLevel mesh_level = MeshLevel::Medium;
    int n_theta = 360;
    int n_z = 96;
    int cross_film_layers = 16;
    int steps_per_revolution = 1440;
    double revolutions = 3.0;
    int outputs_per_revolution = 36;

    ConvergenceInput convergence;
    int cores = 1;
};

/// Result of reading a config into a case.
struct CaseRead {
    CaseModel model;
    /// Config keys whose value differs from what the GUI would write for this case
    /// on a new file: settings the GUI does not show. They are kept on save.
    std::vector<std::string> hidden_keys;
    /// One-sentence notes about what could not be represented (shown in the log).
    std::vector<std::string> notes;
};

/// GUI metadata stored in `#@gui key = value` comment lines.
using GuiMeta = std::map<std::string, std::string>;
GuiMeta read_gui_meta(const std::string& config_text);

/// Template config for a new GUI case (deck-run conventions: DETAILED log,
/// diagnostics every step, flat VTK, default output fields).
SimulationConfig new_case_template();

/// Writes the solver keys of `model` onto `base`. Only groups that differ from
/// `original` are written; with no original, all groups are written. `pair`
/// must be loaded when the fluid state is written (else the keys stay as in base).
SimulationConfig to_config(const CaseModel& model, const FluidPair* pair,
                           const SimulationConfig& base,
                           const std::optional<CaseModel>& original = std::nullopt);

/// Reads a case from a config. `pair` (optional) resolves "saturated" when no
/// metadata says so. Fills hidden_keys by comparing with a fresh GUI case.
CaseRead from_config(const SimulationConfig& cfg, const GuiMeta& meta, const FluidPair* pair);

/// Metadata lines for a case (pair id, saturated flag, mesh level, cores).
GuiMeta case_meta(const CaseModel& model);

/// True when the case runs the COUPLED solver: gas transport chosen, or a thermal
/// case with refrigerant in the oil (FBNS film energy needs pure oil).
bool uses_coupled(const CaseModel& model, const FluidPair* pair);
/// One sentence when the model choices form a combination the solver does not
/// support (the Solution node turns red); empty when supported.
std::string unsupported_reason(const CaseModel& model, const FluidPair* pair);

/// Content actually used by the case [kg/kg] (saturated or given).
double effective_content(const CaseModel& model, const FluidPair* pair);

/// Splits config text into key -> value (repeated inlet lines joined by '\n').
std::map<std::string, std::string> config_key_values(const std::string& text);

}  // namespace gui
