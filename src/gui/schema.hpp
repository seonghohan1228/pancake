#pragma once

// Data-driven description of SimulationConfig for the case-setup UI:
// labels, units, help text, bounds, grouping, and common-vs-advanced split.
// Enums and inlets have bespoke widgets in the setup panel; everything
// numeric/bool/string is generated from these specs.

#include "config.hpp"

#include <string>
#include <vector>

enum class ParamGroup {
    Geometry,
    Operating,
    Lubricant,
    MeshTime,
    Cavitation,
    Boundaries,
    Thermal,
    FluidProperties,
    Motion,
    Output,
    Numerics,
};

const char* group_label(ParamGroup group);

struct DoubleSpec {
    const char* key;    // config file key (= validation highlight id)
    const char* label;
    const char* unit;   // display unit, SI; "-" for dimensionless
    const char* help;   // tooltip; constraints and engineering hints
    double SimulationConfig::* member;
    ParamGroup group;
    bool advanced;
    double min;
    double max;
    bool min_exclusive;  // true: value must be > min, false: >= min
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
