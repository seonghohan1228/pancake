#pragma once

// Display units per field. The case model and config stay SI; each unit-aware
// widget converts for display only (affine: display = si * scale + offset).
// The first option of each family is the engineering default shown to a new user.

#include <vector>

namespace gui {

enum class UnitFamily {
    Fixed,          // dimensionless or single unit
    Pressure,       // bar, MPa, kPa, Pa
    Length,         // mm, m
    SmallLength,    // um, mm, m (clearance, film thickness, roughness)
    Angle,          // deg, rad
    TiltAngle,      // mrad, deg (stored: rad = small-angle slope)
    Speed,          // rpm, Hz, rad/s (stored: rad/s)
    Temperature,    // degC, K
    Content,        // wt%, kg/kg
    Force,          // N, kN
    Mass,           // kg, g
    Time,           // s, ms
    Viscosity,      // mPa s, Pa s
    HeatTransfer,   // W/(m^2 K)
    Conductivity,   // W/(m K)
    Torque,         // N m, N mm
    Power,          // W, kW
};

struct UnitOption {
    const char* label;
    double scale;        // display = si * scale + offset
    double offset = 0.0;
};

const std::vector<UnitOption>& unit_options(UnitFamily family);

inline double to_display(const UnitOption& unit, double si) { return si * unit.scale + unit.offset; }
inline double to_si(const UnitOption& unit, double display) {
    return (display - unit.offset) / unit.scale;
}

}  // namespace gui
