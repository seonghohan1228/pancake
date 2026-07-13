#pragma once

// Data-driven description of SimulationConfig for the case-setup UI:
// labels, units, help text, bounds, grouping, and common-vs-advanced split.
// Enums and inlets have bespoke widgets in the setup panel; everything
// numeric/bool/string is generated from these specs.

#include "config.hpp"

#include <string>
#include <vector>

// Listed in display order: frequently tweaked parameters first so the common
// edit -> run loop needs no scrolling.
enum class ParamGroup {
    Operating,
    Geometry,
    Lubricant,
    Cavitation,
    TimeStepping,
    Mesh,
    Boundaries,
    Thermal,
    Motion,
    FluidProperties,
    Output,
    Numerics,
};

const char* group_label(ParamGroup group);

/// Display-unit family of a numeric field. The config file always stores the
/// solver-native value (first option of each family); the form converts for
/// display only.
enum class UnitFamily {
    Fixed,         // single unit, shown as plain text
    Pressure,      // Pa, kPa, MPa, bar       (stored: Pa)
    BulkModulus,   // Pa, kPa, MPa, GPa       (stored: Pa)
    Angle,         // deg, rad                (stored: deg)
    AngularSpeed,  // rad/s, rpm              (stored: rad/s)
    Length,        // m, mm, um               (stored: m)
    Time,          // s, ms                   (stored: s)
    TiltSlope,     // m/m, mrad               (stored: m/m, small-angle slope)
    Temperature,   // K, degC                 (stored: K)
    MassFraction,  // kg/kg, %                (stored: kg/kg)
    Stiffness,     // N/m, kN/m, MN/m         (stored: N/m)
    Damping,       // N.s/m, kN.s/m           (stored: N.s/m)
};

/// One display unit. display = stored * to_display + offset (offset is non-zero
/// only for affine scales like Celsius); stored = (display - offset) / to_display.
struct UnitOption {
    const char* label;
    double to_display;
    double offset = 0.0;
};

/// Affine display<->stored conversions shared by every unit-aware widget.
inline double unit_to_display(const UnitOption& unit, double stored) {
    return stored * unit.to_display + unit.offset;
}
inline double unit_to_stored(const UnitOption& unit, double display) {
    return (display - unit.offset) / unit.to_display;
}

/// Options for a family; size 1 for Fixed (label unused).
const std::vector<UnitOption>& unit_options(UnitFamily family);

struct DoubleSpec {
    const char* key;    // config file key (= validation highlight id)
    const char* label;
    const char* unit;   // solver-native unit, SI; "-" for dimensionless
    const char* help;   // tooltip; constraints and engineering hints
    double SimulationConfig::* member;
    ParamGroup group;
    bool advanced;
    double min;
    double max;
    bool min_exclusive;  // true: value must be > min, false: >= min
    UnitFamily family = UnitFamily::Fixed;
};

struct IntSpec {
    const char* key;
    const char* label;
    const char* unit;
    const char* help;
    int SimulationConfig::* member;
    ParamGroup group;
    bool advanced;
    int min;
    int max;
};

struct BoolSpec {
    const char* key;
    const char* label;
    const char* help;
    bool SimulationConfig::* member;
    ParamGroup group;
    bool advanced;
};

const std::vector<DoubleSpec>& double_specs();
const std::vector<IntSpec>& int_specs();
const std::vector<BoolSpec>& bool_specs();

struct ValidationIssue {
    bool is_error = true;  // false = warning (does not block Run)
    std::string field;     // spec key, for red highlighting; may be empty
    std::string message;
};

/// Range checks from the specs plus cross-field physical consistency rules.
std::vector<ValidationIssue> validate_config(const SimulationConfig& config);
