#include "units.hpp"

namespace gui {

const std::vector<UnitOption>& unit_options(UnitFamily family) {
    constexpr double kPi = 3.14159265358979323846;
    static const std::vector<UnitOption> fixed = {{"", 1.0}};
    static const std::vector<UnitOption> pressure = {
        {"bar", 1e-5}, {"MPa", 1e-6}, {"kPa", 1e-3}, {"Pa", 1.0}};
    static const std::vector<UnitOption> length = {{"mm", 1e3}, {"m", 1.0}};
    static const std::vector<UnitOption> small_length = {{"um", 1e6}, {"mm", 1e3}, {"m", 1.0}};
    static const std::vector<UnitOption> angle = {{"deg", 1.0}, {"rad", kPi / 180.0}};
    static const std::vector<UnitOption> tilt = {{"mrad", 1e3}, {"deg", 180.0 / kPi}};
    static const std::vector<UnitOption> speed = {
        {"rpm", 60.0 / (2.0 * kPi)}, {"Hz", 1.0 / (2.0 * kPi)}, {"rad/s", 1.0}};
    static const std::vector<UnitOption> temperature = {{"degC", 1.0, -273.15}, {"K", 1.0}};
    static const std::vector<UnitOption> content = {{"wt%", 100.0}, {"kg/kg", 1.0}};
    static const std::vector<UnitOption> force = {{"N", 1.0}, {"kN", 1e-3}};
    static const std::vector<UnitOption> mass = {{"kg", 1.0}, {"g", 1e3}};
    static const std::vector<UnitOption> time = {{"s", 1.0}, {"ms", 1e3}};
    static const std::vector<UnitOption> viscosity = {{"mPa s", 1e3}, {"Pa s", 1.0}};
    static const std::vector<UnitOption> heat_transfer = {{"W/(m2 K)", 1.0}};
    static const std::vector<UnitOption> conductivity = {{"W/(m K)", 1.0}};
    static const std::vector<UnitOption> torque = {{"N m", 1.0}, {"N mm", 1e3}};
    static const std::vector<UnitOption> power = {{"W", 1.0}, {"kW", 1e-3}};
    switch (family) {
        case UnitFamily::Pressure: return pressure;
        case UnitFamily::Length: return length;
        case UnitFamily::SmallLength: return small_length;
        case UnitFamily::Angle: return angle;
        case UnitFamily::TiltAngle: return tilt;
        case UnitFamily::Speed: return speed;
        case UnitFamily::Temperature: return temperature;
        case UnitFamily::Content: return content;
        case UnitFamily::Force: return force;
        case UnitFamily::Mass: return mass;
        case UnitFamily::Time: return time;
        case UnitFamily::Viscosity: return viscosity;
        case UnitFamily::HeatTransfer: return heat_transfer;
        case UnitFamily::Conductivity: return conductivity;
        case UnitFamily::Torque: return torque;
        case UnitFamily::Power: return power;
        case UnitFamily::Fixed: break;
    }
    return fixed;
}

}  // namespace gui
