#include "schema.hpp"

#include <cmath>

const char* group_label(ParamGroup group) {
    switch (group) {
        case ParamGroup::Geometry: return "Geometry";
        case ParamGroup::Operating: return "Operating conditions";
        case ParamGroup::Lubricant: return "Lubricant";
        case ParamGroup::MeshTime: return "Mesh & time";
        case ParamGroup::Cavitation: return "Cavitation";
        case ParamGroup::Boundaries: return "Axial boundaries & inlets";
        case ParamGroup::Thermal: return "Thermal model";
        case ParamGroup::FluidProperties: return "Fluid property models";
        case ParamGroup::Motion: return "Bearing motion";
        case ParamGroup::Output: return "Output";
        case ParamGroup::Numerics: return "Numerics";
    }
    return "?";
}

namespace {
constexpr double kInf = 1e300;
using SC = SimulationConfig;
}  // namespace

const std::vector<DoubleSpec>& double_specs() {
    static const std::vector<DoubleSpec> specs = {
        // Geometry
        {"R", "Journal radius", "m", "Shaft radius R. Typical 5-100 mm.", &SC::R,
         ParamGroup::Geometry, false, 0.0, kInf, true},
        {"c", "Radial clearance", "m", "Nominal radial clearance c. Typical 0.5-2 per mille of R (e.g. 1e-4 m = 100 um).",
         &SC::c, ParamGroup::Geometry, false, 0.0, kInf, true},
        {"e", "Eccentricity", "m", "Journal-bearing center offset. Must satisfy 0 <= e < c. Ignored as initial condition when the motion model derives it from the attitude angle.",
         &SC::e, ParamGroup::Geometry, false, 0.0, kInf, false},
        {"L", "Bearing length", "m", "Axial length L. L/D ratio typically 0.25-1.", &SC::L,
         ParamGroup::Geometry, false, 0.0, kInf, true},
        {"attitude_angle_deg", "Attitude angle", "deg",
         "Initial angle between load line and line of centers.", &SC::attitude_angle_deg,
         ParamGroup::Geometry, false, -360.0, 360.0, false},
        {"load_angle_deg", "Load angle", "deg", "Direction of the external load in the x-y plane.",
         &SC::load_angle_deg, ParamGroup::Geometry, false, -360.0, 360.0, false},
        {"tilt_slope_x", "Tilt slope x", "m/m", "Journal misalignment slope along x.",
         &SC::tilt_slope_x, ParamGroup::Geometry, true, -kInf, kInf, false},
        {"tilt_slope_y", "Tilt slope y", "m/m", "Journal misalignment slope along y.",
         &SC::tilt_slope_y, ParamGroup::Geometry, true, -kInf, kInf, false},

        // Operating conditions
        {"omega", "Rotational speed", "rad/s",
         "Journal angular velocity. 100 rad/s = 955 rpm.", &SC::omega, ParamGroup::Operating,
         false, -kInf, kInf, false},
        {"omega_ramp_time", "Speed ramp time", "s",
         "Linear spin-up duration from rest; 0 starts at full speed.", &SC::omega_ramp_time,
         ParamGroup::Operating, true, 0.0, kInf, false},

        // Lubricant
        {"mu", "Dynamic viscosity", "Pa.s",
         "Lubricant dynamic viscosity at reference conditions. ISO VG32 at 40 C is about 0.027 Pa.s.",
         &SC::mu, ParamGroup::Lubricant, false, 0.0, kInf, true},
        {"rho", "Density", "kg/m^3", "Lubricant density. Mineral oils are about 850-900 kg/m^3.",
         &SC::rho, ParamGroup::Lubricant, false, 0.0, kInf, true},
        {"p_cav", "Cavitation pressure", "Pa",
         "Absolute pressure at which the film ruptures. 0 for gauge-zero rupture.", &SC::p_cav,
         ParamGroup::Lubricant, false, -kInf, kInf, false},
        {"bulk_modulus", "Bulk modulus", "Pa",
         "Liquid compressibility for the Elrod-Adams model. Oils are about 1-2 GPa.",
         &SC::bulk_modulus, ParamGroup::Lubricant, true, 0.0, kInf, true},

        // Mesh & time
        {"end_t", "End time", "s", "Simulated duration (transient mode).", &SC::end_t,
         ParamGroup::MeshTime, false, 0.0, kInf, true},
        {"dt", "Time step", "s", "Time step size. Must not exceed the end time.", &SC::dt,
         ParamGroup::MeshTime, false, 0.0, kInf, true},
        {"write_interval", "Write interval", "s",
         "Simulated time between saved output steps.", &SC::write_interval, ParamGroup::MeshTime,
         false, 0.0, kInf, true},

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
         ParamGroup::Boundaries, false, -kInf, kInf, false},
        {"bc_z_north_val", "North boundary pressure", "Pa",
         "Pressure value at z = L (DIRICHLET / INLET_OUTLET).", &SC::bc_z_north_val,
         ParamGroup::Boundaries, false, -kInf, kInf, false},

        // Thermal
        {"temperature_initial", "Initial temperature", "K", "Initial film temperature.",
         &SC::temperature_initial, ParamGroup::Thermal, false, 0.0, kInf, true},
        {"temperature_reference", "Reference temperature", "K",
         "Reference temperature for property models.", &SC::temperature_reference,
         ParamGroup::Thermal, true, 0.0, kInf, true},
        {"journal_wall_temperature", "Journal wall temperature", "K",
         "Journal surface temperature for wall heat exchange.", &SC::journal_wall_temperature,
         ParamGroup::Thermal, true, 0.0, kInf, true},
        {"bearing_wall_temperature", "Bearing wall temperature", "K",
         "Bearing surface temperature for wall heat exchange.", &SC::bearing_wall_temperature,
         ParamGroup::Thermal, true, 0.0, kInf, true},
        {"rho_cp", "Volumetric heat capacity", "J/(m^3.K)",
         "rho*cp of the lubricant. Mineral oil is about 1.9e6.", &SC::rho_cp, ParamGroup::Thermal,
         true, 0.0, kInf, true},
        {"thermal_conductivity", "Thermal conductivity", "W/(m.K)",
         "Lubricant conductivity. Mineral oil is about 0.13-0.15.", &SC::thermal_conductivity,
         ParamGroup::Thermal, true, 0.0, kInf, true},
        {"journal_heat_transfer", "Journal heat transfer coeff.", "W/(m^2.K)",
         "Wall heat-transfer coefficient to the journal; 0 = adiabatic.",
         &SC::journal_heat_transfer, ParamGroup::Thermal, true, 0.0, kInf, false},
        {"bearing_heat_transfer", "Bearing heat transfer coeff.", "W/(m^2.K)",
         "Wall heat-transfer coefficient to the bearing; 0 = adiabatic.",
         &SC::bearing_heat_transfer, ParamGroup::Thermal, true, 0.0, kInf, false},

        // Fluid property models (all advanced)
        {"dissolved_gas_initial", "Initial dissolved gas", "kg/kg",
         "Initial dissolved gas mass fraction.", &SC::dissolved_gas_initial,
         ParamGroup::FluidProperties, true, 0.0, 1.0, false},
        {"dissolved_gas_max", "Max dissolved gas", "kg/kg",
         "Upper bound of the dissolved gas mass fraction.", &SC::dissolved_gas_max,
         ParamGroup::FluidProperties, true, 0.0, 1.0, false},
        {"dissolved_gas_henry_coeff", "Henry coefficient", "1/Pa",
         "Henry's-law slope at the reference temperature: c_sat = H p.",
         &SC::dissolved_gas_henry_coeff, ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"dissolved_gas_henry_temp_coeff", "Henry temperature coeff.", "K",
         "van't Hoff exponent for H(T).", &SC::dissolved_gas_henry_temp_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"dissolved_gas_bunsen_coeff", "Bunsen coefficient", "-",
         "Bunsen volume-based solubility coefficient.", &SC::dissolved_gas_bunsen_coeff,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"dissolved_gas_liquid_density", "Dissolved gas liquid density", "kg/m^3",
         "Liquid-phase density of the dissolved gas for volume mixing; 0 uses ideal gas.",
         &SC::dissolved_gas_liquid_density, ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"gas_mass_transfer_rate", "Gas mass-transfer rate", "1/s",
         "Finite-rate gas release/resorption coefficient.", &SC::gas_mass_transfer_rate,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"dissolved_gas_diffusivity", "Dissolved gas diffusivity", "m^2/s",
         "Fickian diffusivity of dissolved gas.", &SC::dissolved_gas_diffusivity,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"solution_density_gas_coeff", "Solution density gas coeff.", "-",
         "Density correction by dissolved-gas concentration.", &SC::solution_density_gas_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"solution_viscosity_gas_coeff", "Solution viscosity gas coeff.", "-",
         "Log-mixing viscosity exponent a_c.", &SC::solution_viscosity_gas_coeff,
         ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_temperature_coeff", "Viscosity temperature coeff.", "K",
         "Andrade slope E_mu: mu(T) = mu exp(E_mu (1/T - 1/T_ref)).",
         &SC::viscosity_temperature_coeff, ParamGroup::FluidProperties, true, -kInf, kInf, false},
        {"viscosity_pressure_coeff", "Piezo-viscosity coeff.", "1/Pa",
         "Barus coefficient alpha_p: mu *= exp(alpha_p (p - p_ref)).",
         &SC::viscosity_pressure_coeff, ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"gas_alpha_max", "Max free-gas fraction", "-",
         "Maximum resolved free-gas volume fraction.", &SC::gas_alpha_max,
         ParamGroup::FluidProperties, true, 0.0, 1.0, true},
        {"gas_pressure_floor", "Gas pressure floor", "Pa",
         "Lower pressure bound for the ideal-gas density.", &SC::gas_pressure_floor,
         ParamGroup::FluidProperties, true, 0.0, kInf, true},
        {"mu_gas", "Free-gas viscosity", "Pa.s",
         "Gas dynamic viscosity, used by the two-phase viscosity models.", &SC::mu_gas,
         ParamGroup::FluidProperties, true, 0.0, kInf, false},
        {"property_reference_pressure", "Reference pressure", "Pa",
         "Reference pressure for the property models.", &SC::property_reference_pressure,
         ParamGroup::FluidProperties, true, 0.0, kInf, true},
        {"property_reference_temperature", "Property reference temperature", "K",
         "Reference temperature for the property models.", &SC::property_reference_temperature,
         ParamGroup::FluidProperties, true, 0.0, kInf, true},

        // Bearing motion
        {"external_load_x", "External load x", "N", "Applied load on the bearing, x component.",
         &SC::external_load_x, ParamGroup::Motion, false, -kInf, kInf, false},
        {"external_load_y", "External load y", "N", "Applied load on the bearing, y component.",
         &SC::external_load_y, ParamGroup::Motion, false, -kInf, kInf, false},
        {"external_load_z", "External load z", "N", "Applied load on the bearing, z component.",
         &SC::external_load_z, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_mass", "Bearing mass", "kg", "Moving bearing mass (MOVING_BEARING model).",
         &SC::bearing_mass, ParamGroup::Motion, true, 0.0, kInf, true},
        {"bearing_initial_x", "Initial position x", "m", "Initial bearing offset, x.",
         &SC::bearing_initial_x, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_initial_y", "Initial position y", "m", "Initial bearing offset, y.",
         &SC::bearing_initial_y, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_initial_z", "Initial position z", "m", "Initial bearing offset, z.",
         &SC::bearing_initial_z, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_initial_vx", "Initial velocity x", "m/s", "Initial bearing velocity, x.",
         &SC::bearing_initial_vx, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_initial_vy", "Initial velocity y", "m/s", "Initial bearing velocity, y.",
         &SC::bearing_initial_vy, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_initial_vz", "Initial velocity z", "m/s", "Initial bearing velocity, z.",
         &SC::bearing_initial_vz, ParamGroup::Motion, true, -kInf, kInf, false},
        {"bearing_stiffness_x", "Support stiffness x", "N/m", "Support spring stiffness, x.",
         &SC::bearing_stiffness_x, ParamGroup::Motion, true, 0.0, kInf, false},
        {"bearing_stiffness_y", "Support stiffness y", "N/m", "Support spring stiffness, y.",
         &SC::bearing_stiffness_y, ParamGroup::Motion, true, 0.0, kInf, false},
        {"bearing_stiffness_z", "Support stiffness z", "N/m", "Support spring stiffness, z.",
         &SC::bearing_stiffness_z, ParamGroup::Motion, true, 0.0, kInf, false},
        {"bearing_damping_x", "Support damping x", "N.s/m", "Support damping, x.",
         &SC::bearing_damping_x, ParamGroup::Motion, true, 0.0, kInf, false},
        {"bearing_damping_y", "Support damping y", "N.s/m", "Support damping, y.",
         &SC::bearing_damping_y, ParamGroup::Motion, true, 0.0, kInf, false},
        {"bearing_damping_z", "Support damping z", "N.s/m", "Support damping, z.",
         &SC::bearing_damping_z, ParamGroup::Motion, true, 0.0, kInf, false},
        {"min_film_thickness", "Minimum film thickness", "m",
         "Simulation guard: stop/clamp when the film drops below this.", &SC::min_film_thickness,
         ParamGroup::Motion, true, 0.0, kInf, true},

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
         ParamGroup::MeshTime, false, 3, 100000},
        {"n_z_global", "Cells along axis", "-", "Axial mesh resolution. 20-100 is typical.",
         &SC::n_z_global, ParamGroup::MeshTime, false, 1, 100000},
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
        {"bearing_initial_from_attitude", "Initial position from attitude angle",
         "Derive the initial bearing position from eccentricity and attitude angle instead of explicit coordinates.",
         &SC::bearing_initial_from_attitude, ParamGroup::Motion, true},
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

    // Cross-field physical consistency.
    if (config.e >= config.c) {
        issues.push_back({true, "e",
                          "Eccentricity must be smaller than the radial clearance (e < c), "
                          "otherwise the surfaces touch."});
    }
    if (config.c >= config.R) {
        issues.push_back({true, "c", "Radial clearance must be much smaller than the radius."});
    }
    if (config.solution_mode == SolutionMode::TRANSIENT) {
        if (config.dt > config.end_t) {
            issues.push_back({true, "dt", "Time step exceeds the end time."});
        }
        if (config.write_interval < config.dt) {
            issues.push_back({false, "write_interval",
                              "Write interval is below the time step; every step will be saved."});
        }
        const double steps = config.end_t / config.dt;
        if (steps > 2e6) {
            issues.push_back({false, "dt",
                              "More than 2 million time steps configured; this run will take very long."});
        }
    }
    if (config.cavitation_model == CavitationModel::ELROD_ADAMS) {
        if (config.bulk_modulus <= 0.0) {
            issues.push_back({true, "bulk_modulus",
                              "Elrod-Adams cavitation requires a positive bulk modulus."});
        }
    }
    const auto check_bc = [&](BCType type, double value, const char* key, const char* side) {
        if (type == BCType::DIRICHLET && value < config.p_cav) {
            issues.push_back({false, key,
                              std::string(side) +
                                  " boundary pressure is below the cavitation pressure; the "
                                  "boundary will cavitate."});
        }
    };
    check_bc(config.bc_z_south_type, config.bc_z_south_val, "bc_z_south_val", "South");
    check_bc(config.bc_z_north_type, config.bc_z_north_val, "bc_z_north_val", "North");

    if (config.temperature_model == TemperatureModel::ENERGY_EQUATION) {
        if (config.rho_cp <= 0.0) {
            issues.push_back({true, "rho_cp",
                              "The energy equation requires a positive volumetric heat capacity."});
        }
        if (config.thermal_conductivity <= 0.0) {
            issues.push_back({true, "thermal_conductivity",
                              "The energy equation requires a positive thermal conductivity."});
        }
    }
    if (config.fluid_property_model != FluidPropertyModel::CONSTANT) {
        if (config.dissolved_gas_initial > config.dissolved_gas_max) {
            issues.push_back({true, "dissolved_gas_initial",
                              "Initial dissolved gas exceeds the configured maximum."});
        }
    }
    if (config.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE) {
        const bool needs_mu_gas =
            config.gas_mixture_viscosity_model == GasMixtureViscosityModel::DUKLER_VOID ||
            config.gas_mixture_viscosity_model == GasMixtureViscosityModel::MCADAMS_QUALITY ||
            config.gas_mixture_viscosity_model == GasMixtureViscosityModel::LINEAR_QUALITY;
        if (needs_mu_gas && config.mu_gas <= 0.0) {
            issues.push_back({true, "mu_gas",
                              "The selected two-phase viscosity model needs a positive free-gas "
                              "viscosity."});
        }
    }
    if (config.motion_model == MotionModel::MOVING_BEARING && config.bearing_mass <= 0.0) {
        issues.push_back({true, "bearing_mass", "The moving bearing needs a positive mass."});
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
                          "Diagnostics are disabled; the convergence plot will stay empty."});
    }
    for (size_t i = 0; i < config.inlets.size(); ++i) {
        const InletConfig& inlet = config.inlets[i];
        const std::string label = "Inlet " + std::to_string(i + 1);
        if (inlet.size <= 0.0) {
            issues.push_back({true, "inlets", label + ": size must be positive."});
        }
        if (inlet.p_supply < config.p_cav) {
            issues.push_back({false, "inlets",
                              label + ": supply pressure is below the cavitation pressure."});
        }
        if (inlet.type == InletConfig::Type::CIRCULAR &&
            (inlet.z < 0.0 || inlet.z > config.L)) {
            issues.push_back({true, "inlets", label + ": axial position is outside [0, L]."});
        }
    }
    return issues;
}
