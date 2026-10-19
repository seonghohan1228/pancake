#include "schema.hpp"

#include <algorithm>
#include <cmath>
#include <set>

const char* group_label(ParamGroup group) {
    switch (group) {
        case ParamGroup::Operating: return "Operating conditions";
        case ParamGroup::Geometry: return "Geometry";
        case ParamGroup::Lubricant: return "Lubricant";
        case ParamGroup::Cavitation: return "Cavitation";
        case ParamGroup::TimeStepping: return "Time stepping";
        case ParamGroup::Mesh: return "Mesh";
        case ParamGroup::Boundaries: return "Axial boundaries & inlets";
        case ParamGroup::Thermal: return "Thermal model";
        case ParamGroup::Motion: return "Bearing motion";
        case ParamGroup::FluidProperties: return "Fluid property models";
        case ParamGroup::Output: return "Output";
        case ParamGroup::Numerics: return "Numerics";
    }
    return "?";
}

const std::vector<UnitOption>& unit_options(UnitFamily family) {
    static const std::vector<UnitOption> fixed = {{"", 1.0}};
    static const std::vector<UnitOption> pressure = {
        {"Pa", 1.0}, {"kPa", 1e-3}, {"MPa", 1e-6}, {"bar", 1e-5}};
    static const std::vector<UnitOption> bulk_modulus = {
        {"Pa", 1.0}, {"kPa", 1e-3}, {"MPa", 1e-6}, {"GPa", 1e-9}};
    static const std::vector<UnitOption> angle = {
        {"deg", 1.0}, {"rad", 3.14159265358979323846 / 180.0}};
    static const std::vector<UnitOption> angular_speed = {
        {"rad/s", 1.0}, {"rpm", 60.0 / (2.0 * 3.14159265358979323846)}};
    static const std::vector<UnitOption> length = {{"m", 1.0}, {"mm", 1e3}, {"um", 1e6}};
    static const std::vector<UnitOption> time = {{"s", 1.0}, {"ms", 1e3}};
    // Misalignment slope (dh/dz). For the small tilts of a real bearing the slope
    // and the tilt angle coincide (tan x ~ x), so 1 m/m ~ 1 rad = 1000 mrad.
    static const std::vector<UnitOption> tilt_slope = {{"m/m", 1.0}, {"mrad", 1e3}};
    static const std::vector<UnitOption> temperature = {{"K", 1.0, 0.0}, {"degC", 1.0, -273.15}};
    static const std::vector<UnitOption> mass_fraction = {{"kg/kg", 1.0}, {"%", 100.0}};
    static const std::vector<UnitOption> stiffness = {
        {"N/m", 1.0}, {"kN/m", 1e-3}, {"MN/m", 1e-6}};
    static const std::vector<UnitOption> damping = {{"N.s/m", 1.0}, {"kN.s/m", 1e-3}};
    switch (family) {
        case UnitFamily::Pressure: return pressure;
        case UnitFamily::BulkModulus: return bulk_modulus;
        case UnitFamily::Angle: return angle;
        case UnitFamily::AngularSpeed: return angular_speed;
        case UnitFamily::Length: return length;
        case UnitFamily::Time: return time;
        case UnitFamily::TiltSlope: return tilt_slope;
        case UnitFamily::Temperature: return temperature;
        case UnitFamily::MassFraction: return mass_fraction;
        case UnitFamily::Stiffness: return stiffness;
        case UnitFamily::Damping: return damping;
        case UnitFamily::Fixed: break;
    }
    return fixed;
}

namespace {
constexpr double kInf = 1e300;
using SC = SimulationConfig;
constexpr UnitFamily kPressure = UnitFamily::Pressure;
constexpr UnitFamily kBulkModulus = UnitFamily::BulkModulus;
constexpr UnitFamily kAngle = UnitFamily::Angle;
constexpr UnitFamily kAngularSpeed = UnitFamily::AngularSpeed;
constexpr UnitFamily kLength = UnitFamily::Length;
constexpr UnitFamily kTime = UnitFamily::Time;
constexpr UnitFamily kTiltSlope = UnitFamily::TiltSlope;
constexpr UnitFamily kTemperature = UnitFamily::Temperature;
constexpr UnitFamily kMassFraction = UnitFamily::MassFraction;
constexpr UnitFamily kStiffness = UnitFamily::Stiffness;
constexpr UnitFamily kDamping = UnitFamily::Damping;
}  // namespace

const std::vector<DoubleSpec>& double_specs() {
    static const std::vector<DoubleSpec> specs = {
        // Geometry
        {"R", "Journal radius", "m", "Shaft radius R. Typical 5-100 mm.", &SC::R,
         ParamGroup::Geometry, false, 0.0, kInf, true, kLength},
        {"c", "Radial clearance", "m", "Nominal radial clearance c. Typical 0.5-2 per mille of R (e.g. 1e-4 m = 100 um).",
         &SC::c, ParamGroup::Geometry, false, 0.0, kInf, true, kLength},
        {"e", "Eccentricity", "m", "Journal-bearing center offset. Must satisfy 0 <= e < c. Ignored as initial condition when the motion model derives it from the attitude angle.",
         &SC::e, ParamGroup::Geometry, false, 0.0, kInf, false, kLength},
        {"L", "Bearing length", "m", "Axial length L. L/D ratio typically 0.25-1.", &SC::L,
         ParamGroup::Geometry, false, 0.0, kInf, true, kLength},
        {"attitude_angle_deg", "Attitude angle", "deg",
         "Initial angle between load line and line of centers.", &SC::attitude_angle_deg,
         ParamGroup::Geometry, false, -360.0, 360.0, false, kAngle},
        {"load_angle_deg", "Load angle", "deg", "Direction of the external load in the x-y plane.",
         &SC::load_angle_deg, ParamGroup::Geometry, false, -360.0, 360.0, false, kAngle},
        {"tilt_slope_x", "Tilt slope x", "m/m", "Journal misalignment slope along x "
         "(displayable as mrad of tilt; equal for the small slopes of a real bearing).",
         &SC::tilt_slope_x, ParamGroup::Geometry, true, -kInf, kInf, false, kTiltSlope},
        {"tilt_slope_y", "Tilt slope y", "m/m", "Journal misalignment slope along y "
         "(displayable as mrad of tilt; equal for the small slopes of a real bearing).",
         &SC::tilt_slope_y, ParamGroup::Geometry, true, -kInf, kInf, false, kTiltSlope},

        // Operating conditions
        {"omega", "Rotational speed", "rad/s",
         "Journal angular velocity. 100 rad/s = 955 rpm.", &SC::omega, ParamGroup::Operating,
         false, -kInf, kInf, false, kAngularSpeed},
        {"omega_ramp_time", "Speed ramp time", "s",
         "Linear spin-up duration from rest; 0 starts at full speed.", &SC::omega_ramp_time,
         ParamGroup::Operating, true, 0.0, kInf, false, kTime},

        // Lubricant
        {"mu", "Dynamic viscosity", "Pa.s",
         "Lubricant dynamic viscosity at reference conditions. ISO VG32 at 40 C is about 0.027 Pa.s.",
         &SC::mu, ParamGroup::Lubricant, false, 0.0, kInf, true},
        {"rho", "Density", "kg/m^3", "Lubricant density. Mineral oils are about 850-900 kg/m^3.",
         &SC::rho, ParamGroup::Lubricant, false, 0.0, kInf, true},
        {"p_cav", "Cavitation pressure", "Pa",
         "Absolute pressure at which the film ruptures. About 101325 Pa for atmospheric "
         "cavitation; inputs are absolute, not gauge.", &SC::p_cav,
         ParamGroup::Lubricant, false, 0.0, kInf, true, kPressure},
        {"bulk_modulus", "Bulk modulus", "Pa",
         "Liquid compressibility for the Elrod-Adams model. Oils are about 1-2 GPa.",
         &SC::bulk_modulus, ParamGroup::Lubricant, true, 0.0, kInf, true, kBulkModulus},

        // Time stepping
        {"end_t", "End time", "s", "Simulated duration (transient mode).", &SC::end_t,
         ParamGroup::TimeStepping, false, 0.0, kInf, true, kTime},
        {"dt", "Time step", "s", "Time step size. Must not exceed the end time.", &SC::dt,
         ParamGroup::TimeStepping, false, 0.0, kInf, true, kTime},
        {"write_interval", "Write interval", "s",
         "Simulated time between saved output steps.", &SC::write_interval,
         ParamGroup::TimeStepping, false, 0.0, kInf, true, kTime},

        // Cavitation
        {"outer_tol", "Outer tolerance", "-",
         "Elrod-Adams outer-iteration convergence tolerance on the film content.", &SC::outer_tol,
         ParamGroup::Cavitation, true, 0.0, kInf, true},
        {"theta_min", "Minimum film content", "-",
         "Lower clamp for the liquid fraction in cavitated cells.", &SC::theta_min,
         ParamGroup::Cavitation, true, 0.0, 1.0, true},

        // Boundaries
        {"bc_z_south_val", "South boundary pressure", "Pa",
         "Pressure value at z = 0 (DIRICHLET / INLET_OUTLET).", &SC::bc_z_south_val,
         ParamGroup::Boundaries, false, -kInf, kInf, false, kPressure},
        {"bc_z_north_val", "North boundary pressure", "Pa",
         "Pressure value at z = L (DIRICHLET / INLET_OUTLET).", &SC::bc_z_north_val,
         ParamGroup::Boundaries, false, -kInf, kInf, false, kPressure},

        // Thermal
        {"temperature_initial", "Initial temperature", "K", "Initial film temperature.",
         &SC::temperature_initial, ParamGroup::Thermal, false, 0.0, kInf, true, kTemperature},
        {"temperature_reference", "Reference temperature", "K",
         "Reference temperature for property models.", &SC::temperature_reference,
         ParamGroup::Thermal, true, 0.0, kInf, true, kTemperature},
        {"journal_wall_temperature", "Journal wall temperature", "K",
         "Journal surface temperature for wall heat exchange.", &SC::journal_wall_temperature,
         ParamGroup::Thermal, true, 0.0, kInf, true, kTemperature},
        {"bearing_wall_temperature", "Bearing wall temperature", "K",
         "Bearing surface temperature for wall heat exchange.", &SC::bearing_wall_temperature,
         ParamGroup::Thermal, true, 0.0, kInf, true, kTemperature},
        {"rho_cp", "Volumetric heat capacity", "J/(m^3.K)",
         "rho*cp of the lubricant. Mineral oil is about 1.9e6.", &SC::rho_cp, ParamGroup::Thermal,
         true, 0.0, kInf, true},
        {"thermal_conductivity", "Thermal conductivity", "W/(m.K)",
         "Lubricant conductivity. Mineral oil is about 0.13-0.15.", &SC::thermal_conductivity,
         ParamGroup::Thermal, true, 0.0, kInf, true},
        {"gas_heat_capacity_cv", "Free gas heat capacity (cv)", "J/(kg.K)",
         "COUPLED ideal-gas heat capacity; R410A obtains this from its real-fluid EOS.",
         &SC::gas_heat_capacity_cv, ParamGroup::Thermal, true, -1.0, kInf, false},
        {"dissolved_heat_capacity", "Dissolved species heat capacity", "J/(kg.K)",
         "Partial caloric property of dissolved species. Requires mixture data.",
         &SC::dissolved_heat_capacity, ParamGroup::Thermal, true, -1.0, kInf, false},
        {"desorption_internal_energy", "Desorption internal energy", "J/kg",
         "Free minus dissolved specific internal energy at reference p,T. Positive values cool the film on release. No FW68D value is supplied.",
         &SC::desorption_internal_energy, ParamGroup::Thermal, true, -kInf, kInf, false},
        {"gas_thermal_conductivity", "Free phase thermal conductivity", "W/(m.K)",
         "COUPLED volume-weighted free-phase conduction; zero neglects this contribution.",
         &SC::gas_thermal_conductivity, ParamGroup::Thermal, true, 0.0, kInf, false},
        {"oil_vapor_temperature_coeff", "Oil vapor pressure temperature slope", "K",
         "COUPLED Clausius-Clapeyron slope. Zero uses constant p_cav; fit from the actual oil vapor-pressure data.",
         &SC::oil_vapor_temperature_coeff, ParamGroup::Thermal, true, 0.0, kInf, false},
        {"journal_heat_transfer", "Journal heat transfer coeff.", "W/(m^2.K)",
         "Wall heat-transfer coefficient to the journal; 0 = adiabatic.",
         &SC::journal_heat_transfer, ParamGroup::Thermal, true, 0.0, kInf, false},
        {"bearing_heat_transfer", "Bearing heat transfer coeff.", "W/(m^2.K)",
         "Wall heat-transfer coefficient to the bearing; 0 = adiabatic.",
         &SC::bearing_heat_transfer, ParamGroup::Thermal, true, 0.0, kInf, false},

        // Fluid property models (all advanced)
        {"dissolved_gas_initial", "Initial dissolved gas", "kg/kg",
         "Initial dissolved gas mass fraction (displayable as %).", &SC::dissolved_gas_initial,
         ParamGroup::FluidProperties, true, 0.0, 1.0, false, kMassFraction},
        {"dissolved_gas_max", "Max dissolved gas", "kg/kg",
         "Upper bound of the dissolved gas mass fraction (displayable as %).",
         &SC::dissolved_gas_max, ParamGroup::FluidProperties, true, 0.0, 1.0, false, kMassFraction},
        {"dissolved_gas_henry_coeff", "Henry coefficient", "1/Pa",
         "Henry's-law slope at the reference temperature: c_sat = H p.",
         &SC::dissolved_gas_henry_coeff, ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"dissolved_gas_henry_temp_coeff", "Henry temperature coeff.", "K",
         "van't Hoff exponent for H(T).", &SC::dissolved_gas_henry_temp_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"dissolved_gas_liquid_density", "Dissolved gas liquid density", "kg/m^3",
         "Liquid-phase density of the dissolved gas, the second phase density in the MIXTURE "
         "density model.", &SC::dissolved_gas_liquid_density, ParamGroup::FluidProperties, true,
         0.0, kInf, false},
        {"gas_mass_transfer_rate", "Gas mass-transfer rate", "1/s",
         "Finite-rate gas release/resorption coefficient.", &SC::gas_mass_transfer_rate,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"gas_absorption_rate", "Gas absorption rate", "1/s",
         "-1 inherits release rate; zero preserves bubbles during compression without dissolution.",
         &SC::gas_absorption_rate, ParamGroup::FluidProperties, true, -1.0, kInf, false},
        {"dissolved_gas_diffusivity", "Dissolved gas diffusivity", "m^2/s",
         "Fickian diffusivity of dissolved gas.", &SC::dissolved_gas_diffusivity,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"solution_density_gas_coeff", "Solution density gas coeff.", "-",
         "Density correction by dissolved-gas concentration.", &SC::solution_density_gas_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"density_temperature_intercept", "Density at 0 C", "kg/m^3",
         "LINEAR_TEMPERATURE: rho = intercept + slope * T[degC].",
         &SC::density_temperature_intercept, ParamGroup::FluidProperties, false, -kInf, kInf, false},
        {"density_temperature_slope", "Density temperature slope", "kg/(m^3 K)",
         "LINEAR_TEMPERATURE slope; FW68D default is -0.7284.",
         &SC::density_temperature_slope, ParamGroup::FluidProperties, false, -kInf, kInf, false},
        {"solution_viscosity_gas_coeff", "Solution viscosity gas coeff.", "-",
         "Log-mixing viscosity exponent a_c.", &SC::solution_viscosity_gas_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_temperature_coeff", "Viscosity temperature coeff.", "K",
         "Andrade slope E_mu: mu(T) = mu exp(E_mu (1/T - 1/T_ref)).",
         &SC::viscosity_temperature_coeff, ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_pressure_coeff", "Piezo-viscosity coeff.", "1/Pa",
         "Barus coefficient alpha_p: mu *= exp(alpha_p (p - p_ref)).",
         &SC::viscosity_pressure_coeff, ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"viscosity_poly_k0", "Ferron k0", "-",
         "POLYNOMIAL constant term: mu(T) = mu (k0 - k1 Tbar + k2 Tbar^2), Tbar = T[C]/T_ref[C]. "
         "Ferron 1983 VG32: 3.287. Defaults to 1 (identity).",
         &SC::viscosity_poly_k0, ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_poly_k1", "Ferron k1", "-",
         "POLYNOMIAL linear term, SUBTRACTED. Ferron 1983 VG32: 3.064.",
         &SC::viscosity_poly_k1, ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_poly_k2", "Ferron k2", "-",
         "POLYNOMIAL quadratic term. Ferron 1983 VG32: 0.777. k0 - k1 + k2 should be 1, "
         "since mu is defined at the reference temperature.",
         &SC::viscosity_poly_k2, ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"gas_alpha_max", "Max free-gas fraction", "-",
         "Maximum resolved free-gas volume fraction.", &SC::gas_alpha_max,
         ParamGroup::FluidProperties, true, 0.0, 1.0, true},
        {"gas_pressure_floor", "Gas pressure floor", "Pa",
         "Lower pressure bound for the ideal-gas density.", &SC::gas_pressure_floor,
         ParamGroup::FluidProperties, true, 0.0, kInf, true, kPressure},
        {"mu_gas", "Free-gas viscosity", "Pa.s",
         "Gas dynamic viscosity, used by the two-phase viscosity models.", &SC::mu_gas,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"property_reference_pressure", "Reference pressure", "Pa",
         "Reference pressure for the property models.", &SC::property_reference_pressure,
         ParamGroup::FluidProperties, true, 0.0, kInf, true, kPressure},
        {"property_reference_temperature", "Property reference temperature", "K",
         "Reference temperature for the property models.", &SC::property_reference_temperature,
         ParamGroup::FluidProperties, true, 0.0, kInf, true, kTemperature},

        {"coupled_nonlinear_tolerance", "Coupled nonlinear tolerance", "",
         "Scaled inner pressure/species residual limit. Separate from outer thermal convergence.",
         &SC::coupled_nonlinear_tolerance, ParamGroup::Numerics, true, 0.0, 1.0, true},

        // Bearing motion
        {"external_load_x", "External load x", "N", "Applied load on the bearing, x component.",
         &SC::external_load_x, ParamGroup::Motion, false, -kInf, kInf, false},
        {"external_load_y", "External load y", "N", "Applied load on the bearing, y component.",
         &SC::external_load_y, ParamGroup::Motion, false, -kInf, kInf, false},
        {"external_load_z", "External load z", "N", "Applied load on the bearing, z component.",
         &SC::external_load_z, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_mass", "Bearing mass", "kg", "Moving bearing mass (MOVING_BEARING model).",
         &SC::bearing_mass, ParamGroup::Motion, true, 0.0, kInf, true},
        // Isotropic support: one value applied to every direction. The GUI mirrors
        // it onto the y/z config members when this row is *edited*, never on load
        // or save, so an externally authored anisotropic config survives a
        // load/save cycle intact (it is flagged in the log instead). The initial
        // position comes from the bearing geometry (eccentricity + attitude angle)
        // and the initial velocity is zero, so neither is exposed here; a config
        // that sets them explicitly is likewise left alone.
        {"bearing_stiffness_x", "Support stiffness", "N/m",
         "Isotropic support spring stiffness (same in all directions).",
         &SC::bearing_stiffness_x, ParamGroup::Motion, true, 0.0, kInf, false, kStiffness},
        {"bearing_damping_x", "Support damping", "N.s/m",
         "Isotropic support damping (same in all directions).",
         &SC::bearing_damping_x, ParamGroup::Motion, true, 0.0, kInf, false, kDamping},
        {"min_film_thickness", "Minimum film thickness", "m",
         "Simulation guard: stop/clamp when the film drops below this.", &SC::min_film_thickness,
         ParamGroup::Motion, true, 0.0, kInf, true, kLength},

        // Numerics
        {"linear_rtol", "Linear solver tolerance", "-",
         "Relative tolerance of the inner Krylov solves; bounds the reported mass residual.",
         &SC::linear_rtol, ParamGroup::Numerics, true, 0.0, 1.0, true},
    };
    return specs;
}

const std::vector<IntSpec>& int_specs() {
    static const std::vector<IntSpec> specs = {
        {"n_theta_global", "Cells around circumference", "-",
         "Circumferential mesh resolution. 120-360 is typical.", &SC::n_theta_global,
         ParamGroup::Mesh, false, 3, 100000},
        {"n_z_global", "Cells along axis", "-", "Axial mesh resolution. 20-100 is typical.",
         &SC::n_z_global, ParamGroup::Mesh, false, 1, 100000},
        {"max_outer_iters", "Max outer iterations", "-",
         "Elrod-Adams cavitation flag-update iteration limit per step.", &SC::max_outer_iters,
         ParamGroup::Cavitation, true, 1, 100000},
        {"diagnostics_interval", "Diagnostics interval", "steps",
         "Write a diagnostics.csv row every N steps; 0 disables. The convergence plot needs this on.",
         &SC::diagnostics_interval, ParamGroup::Numerics, true, 0, 1000000},
    };
    return specs;
}

const std::vector<BoolSpec>& bool_specs() {
    static const std::vector<BoolSpec> specs = {
        {"stop_on_nonpositive_film", "Stop on film contact",
         "Abort the run when the film thickness reaches the minimum.",
         &SC::stop_on_nonpositive_film, ParamGroup::Motion, true},
        {"log_outer_iters", "Log outer iterations",
         "Print Elrod-Adams outer-iteration progress to the solver log.", &SC::log_outer_iters,
         ParamGroup::Numerics, true},
        {"output_write_3d", "Write 3D (cylindrical) VTK",
         "Write the wrapped cylindrical surface output for ParaView.", &SC::output_write_3d,
         ParamGroup::Output, true},
        {"output_write_flat", "Write flat (unwrapped) VTK",
         "Write the unwrapped bearing surface output. The GUI heatmap reads either form.",
         &SC::output_write_flat, ParamGroup::Output, true},
    };
    return specs;
}

namespace {

/// Maps a solver validation message to a config key for field highlighting.
/// Solver messages start with the key they complain about ("p_cav must be...",
/// "inlet[0].p_supply must be..."); unknown shapes map to "" (list-only).
std::string field_for_message(const std::string& message) {
    static std::vector<std::string> keys = [] {
        std::vector<std::string> all;
        for (const auto& spec : double_specs()) all.push_back(spec.key);
        for (const auto& spec : int_specs()) all.push_back(spec.key);
        for (const auto& spec : bool_specs()) all.push_back(spec.key);
        // Longest first so "dissolved_gas_max" wins over a hypothetical prefix.
        std::sort(all.begin(), all.end(),
                  [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
        return all;
    }();
    if (message.rfind("inlet[", 0) == 0) return "inlets";
    if (message.rfind("n_theta_global", 0) == 0) return "n_theta_global";
    for (const std::string& key : keys) {
        if (message.size() > key.size() && message.compare(0, key.size(), key) == 0 &&
            message[key.size()] == ' ') {
            return key;
        }
    }
    return "";
}

void check_range(const DoubleSpec& spec, const SimulationConfig& config,
                 std::vector<ValidationIssue>& issues) {
    const double value = config.*(spec.member);
    if (!std::isfinite(value)) {
        issues.push_back({true, spec.key, std::string(spec.label) + " must be a finite number."});
        return;
    }
    if (spec.min_exclusive ? value <= spec.min : value < spec.min) {
        issues.push_back({true, spec.key,
                          std::string(spec.label) + " must be " +
                              (spec.min_exclusive ? "greater than " : "at least ") +
                              (spec.min == 0.0 ? std::string("0") : std::to_string(spec.min)) +
                              " " + spec.unit + "."});
    }
    if (value > spec.max) {
        issues.push_back({true, spec.key,
                          std::string(spec.label) + " must be at most " +
                              std::to_string(spec.max) + " " + spec.unit + "."});
    }
}

}  // namespace

std::vector<ValidationIssue> validate_config(const SimulationConfig& config) {
    std::vector<ValidationIssue> issues;

    for (const auto& spec : double_specs()) check_range(spec, config, issues);
    for (const auto& spec : int_specs()) {
        const int value = config.*(spec.member);
        if (value < spec.min || value > spec.max) {
            issues.push_back({true, spec.key,
                              std::string(spec.label) + " must be between " +
                                  std::to_string(spec.min) + " and " + std::to_string(spec.max) +
                                  "."});
        }
    }

    // GUI-only ergonomics the solver does not check.
    if (config.c > 0.0 && config.R > 0.0 && config.c >= config.R) {
        issues.push_back({true, "c", "Radial clearance must be much smaller than the radius."});
    }
    if (config.solution_mode == SolutionMode::TRANSIENT && config.dt > 0.0 &&
        config.end_t > 0.0) {
        if (config.dt > config.end_t) {
            issues.push_back({false, "dt", "Time step exceeds the end time; only one step runs."});
        }
        if (config.end_t / config.dt > 2e6) {
            issues.push_back({false, "dt",
                              "More than 2 million time steps configured; this run will take "
                              "very long."});
        }
    }
    if (config.output_dir.empty()) {
        issues.push_back({true, "output_dir", "Output directory must not be empty."});
    }
    if (config.filename_prefix.empty()) {
        issues.push_back({true, "filename_prefix", "Filename prefix must not be empty."});
    }
    if (!config.output_write_3d && !config.output_write_flat) {
        issues.push_back({false, "output_write_3d",
                          "Both VTK outputs are disabled; the results view will stay empty."});
    }
    if (config.diagnostics_interval <= 0) {
        issues.push_back({false, "diagnostics_interval",
                          "Diagnostics are disabled; the Monitors plots will stay empty."});
    }
    for (size_t i = 0; i < config.inlets.size(); ++i) {
        const InletConfig& inlet = config.inlets[i];
        if (inlet.type == InletConfig::Type::CIRCULAR &&
            (inlet.z < 0.0 || inlet.z > config.L)) {
            issues.push_back({true, "inlets",
                              "Inlet " + std::to_string(i + 1) +
                                  ": axial position is outside [0, L]."});
        }
    }

    // The solver's own validation is the source of truth: anything it rejects
    // must disable Run, with the field highlighted when the message names it.
    std::set<std::string> flagged_fields;
    for (const ValidationIssue& issue : issues) {
        if (issue.is_error && !issue.field.empty()) flagged_fields.insert(issue.field);
    }
    std::vector<std::string> solver_errors, solver_warnings;
    config.validate(solver_errors, solver_warnings);
    const auto merge = [&](const std::vector<std::string>& messages, bool is_error) {
        for (const std::string& message : messages) {
            std::string field = field_for_message(message);
            if (is_error && !field.empty() && flagged_fields.count(field) > 0) {
                continue;  // a range check already reported this field
            }
            if (is_error && !field.empty()) flagged_fields.insert(field);
            issues.push_back({is_error, field, message});
        }
    };
    merge(solver_errors, true);
    merge(solver_warnings, false);
    return issues;
}
