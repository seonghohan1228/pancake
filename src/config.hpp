#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

enum class CavitationModel { GUMBEL, ELROD_ADAMS, FULL_SOMMERFELD };
enum class BCType { DIRICHLET, NEUMANN, INLET_OUTLET };
enum class MotionModel { STATIC, MOVING_BEARING };
enum class TimeSteppingMethod { EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICOLSON, RK2, RK4 };
enum class SolutionMode { TRANSIENT, STEADY_STATE };
enum class TemperatureModel { ISOTHERMAL, ENERGY_EQUATION };

inline std::string normalise_config_token(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    text = text.substr(first, last - first + 1);
    std::replace(text.begin(), text.end(), '-', '_');
    std::replace(text.begin(), text.end(), ' ', '_');
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return text;
}

inline const char* to_config_value(BCType type) {
    switch (type) {
        case BCType::DIRICHLET: return "DIRICHLET";
        case BCType::NEUMANN: return "NEUMANN";
        case BCType::INLET_OUTLET: return "INLET_OUTLET";
    }
    return "DIRICHLET";
}

inline const char* to_config_value(CavitationModel model) {
    switch (model) {
        case CavitationModel::GUMBEL: return "GUMBEL";
        case CavitationModel::ELROD_ADAMS: return "ELROD_ADAMS";
        case CavitationModel::FULL_SOMMERFELD: return "FULL_SOMMERFELD";
    }
    return "ELROD_ADAMS";
}

inline const char* to_config_value(MotionModel model) {
    switch (model) {
        case MotionModel::STATIC: return "STATIC";
        case MotionModel::MOVING_BEARING: return "MOVING_BEARING";
    }
    return "STATIC";
}

inline const char* to_config_value(TimeSteppingMethod method) {
    switch (method) {
        case TimeSteppingMethod::EULER_EXPLICIT: return "EULER_EXPLICIT";
        case TimeSteppingMethod::EULER_IMPLICIT: return "EULER_IMPLICIT";
        case TimeSteppingMethod::CRANK_NICOLSON: return "CRANK_NICOLSON";
        case TimeSteppingMethod::RK2: return "RK2";
        case TimeSteppingMethod::RK4: return "RK4";
    }
    return "EULER_IMPLICIT";
}

inline const char* to_config_value(SolutionMode mode) {
    switch (mode) {
        case SolutionMode::TRANSIENT: return "TRANSIENT";
        case SolutionMode::STEADY_STATE: return "STEADY_STATE";
    }
    return "TRANSIENT";
}

inline const char* to_config_value(TemperatureModel model) {
    switch (model) {
        case TemperatureModel::ISOTHERMAL: return "ISOTHERMAL";
        case TemperatureModel::ENERGY_EQUATION: return "ENERGY_EQUATION";
    }
    return "ISOTHERMAL";
}

inline MotionModel motion_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "MOVING_BEARING" || value == "DYNAMIC" || value == "BEARING") {
        return MotionModel::MOVING_BEARING;
    }
    return MotionModel::STATIC;
}

inline TimeSteppingMethod time_stepping_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "EULER_EXPLICIT" || value == "EXPLICIT_EULER" || value == "EXPLICIT") {
        return TimeSteppingMethod::EULER_EXPLICIT;
    }
    if (value == "EULER_IMPLICIT" || value == "IMPLICIT_EULER" || value == "IMPLICIT") {
        return TimeSteppingMethod::EULER_IMPLICIT;
    }
    if (value == "CRANK_NICOLSON" || value == "CRANK_NICHOLSON" ||
        value == "SEMI_IMPLICIT" || value == "SEMIIMPLICIT" || value == "IMEX") {
        return TimeSteppingMethod::CRANK_NICOLSON;
    }
    if (value == "RK2" || value == "RUNGE_KUTTA_2" || value == "MIDPOINT") {
        return TimeSteppingMethod::RK2;
    }
    if (value == "RK4" || value == "RUNGE_KUTTA_4") {
        return TimeSteppingMethod::RK4;
    }
    return TimeSteppingMethod::EULER_IMPLICIT;
}

inline SolutionMode solution_mode_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "STEADY_STATE" || value == "STEADY" || value == "SINGLE_STEP") {
        return SolutionMode::STEADY_STATE;
    }
    return SolutionMode::TRANSIENT;
}

inline TemperatureModel temperature_model_from_config(std::string value) {
    value = normalise_config_token(value);
    if (value == "ENERGY_EQUATION" || value == "ENERGY" || value == "THERMAL" || value == "THD") {
        return TemperatureModel::ENERGY_EQUATION;
    }
    return TemperatureModel::ISOTHERMAL;
}

struct InletConfig {
    enum class Type { CIRCULAR, GROOVE } type = Type::CIRCULAR;
    double theta = 0.0;     // Center angular position (deg)
    double z = 0.0;         // Center axial position (m)
    double size = 0.0;      // Radius for CIRCULAR (m), angular width for GROOVE (deg)
    double p_supply = 0.0;  // Supply pressure (Pa)
};

struct SimulationConfig {
    // Geometry
    double R = 0.01;
    double c = 0.0001;
    double e = 0.00008;
    double L = 0.02;
    double attitude_angle_deg = -90.0;
    double load_angle_deg = 0.0;

    // Grid
    int n_theta_global = 120;
    int n_z_global = 40;

    // Time
    double end_t = 1.0;
    double dt = 0.01;
    double write_interval = 0.1;
    SolutionMode solution_mode = SolutionMode::TRANSIENT;
    TimeSteppingMethod pressure_time_method = TimeSteppingMethod::EULER_IMPLICIT;
    TimeSteppingMethod motion_time_method = TimeSteppingMethod::EULER_IMPLICIT;
    TimeSteppingMethod temperature_time_method = TimeSteppingMethod::EULER_IMPLICIT;

    // Physics
    double omega = 100.0;
    double mu = 0.01;
    double rho = 900.0;
    double p_cav = 0.0;
    double bulk_modulus = 1e9;

    // Misalignment (tilt)
    double tilt_slope_x = 0.0; // Slope in x-direction (m/m)
    double tilt_slope_y = 0.0; // Slope in y-direction (m/m)

    // Energy equation / thermal model
    TemperatureModel temperature_model = TemperatureModel::ISOTHERMAL;
    double temperature_initial = 300.0;
    double temperature_reference = 300.0;
    double journal_wall_temperature = 300.0;
    double bearing_wall_temperature = 300.0;
    double rho_cp = 1.92e6;              // Volumetric heat capacity [J/(m^3 K)]
    double thermal_conductivity = 0.15;  // Lubricant thermal conductivity [W/(m K)]
    double journal_heat_transfer = 0.0;  // Wall heat-transfer coefficient [W/(m^2 K)]
    double bearing_heat_transfer = 0.0;  // Wall heat-transfer coefficient [W/(m^2 K)]

    // Bearing motion: outer bearing moves while the shaft remains fixed.
    MotionModel motion_model = MotionModel::STATIC;
    bool bearing_initial_from_attitude = true;
    double bearing_initial_x = 0.0;
    double bearing_initial_y = 0.0;
    double bearing_initial_z = 0.0;
    double bearing_initial_vx = 0.0;
    double bearing_initial_vy = 0.0;
    double bearing_initial_vz = 0.0;
    double bearing_mass = 1.0;
    double bearing_stiffness_x = 0.0;
    double bearing_stiffness_y = 0.0;
    double bearing_stiffness_z = 0.0;
    double bearing_damping_x = 0.0;
    double bearing_damping_y = 0.0;
    double bearing_damping_z = 0.0;
    double external_load_x = 0.0;
    double external_load_y = 0.0;
    double external_load_z = 0.0;
    double min_film_thickness = 1e-9;
    bool stop_on_nonpositive_film = true;

    // Axial Boundary Conditions
    BCType bc_z_south_type = BCType::DIRICHLET;
    double bc_z_south_val  = 0.0;
    double bc_z_south_theta = 1.0; // Phase selection for INLET_OUTLET inflow
    BCType bc_z_north_type = BCType::DIRICHLET;
    double bc_z_north_val  = 0.0;
    double bc_z_north_theta = 1.0; // Phase selection for INLET_OUTLET inflow

    // Inlets
    std::vector<InletConfig> inlets;

    // Cavitation
    CavitationModel cavitation_model = CavitationModel::ELROD_ADAMS;
    int    max_outer_iters = 50;    // Elrod-Adams outer flag-update iteration limit
    double outer_tol       = 1e-6;  // Convergence tolerance on max change between outer iters
    double theta_min       = 1e-6;  // Minimum film content
    bool   log_outer_iters = true;  // Whether to print outer iteration progress

    // Output
    std::string output_dir = "results";
    std::string filename_prefix = "solution";
    bool output_write_3d = true;
    bool output_write_flat = true;
    std::vector<std::string> output_fields = default_output_fields();

    void load_from_file(const std::filesystem::path& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open config file '" << path
                      << "'. Using defaults.\n";
            return;
        }

        load_from_stream(file);
    }

    void load_from_text(const std::string& text) {
        std::istringstream input(text);
        load_from_stream(input);
    }

    void load_from_stream(std::istream& input) {
        inlets.clear();

        std::string line;
        while (std::getline(input, line)) {
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);

            trim(line);
            if (line.empty()) continue;

            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            std::string key = line.substr(0, eq_pos);
            std::string value = line.substr(eq_pos + 1);

            trim(key);
            trim(value);

            if (key.empty()) continue;
            if (value.empty() && key != "output_fields") continue;

            if (key == "R") R = std::stod(value);
            else if (key == "c") c = std::stod(value);
            else if (key == "e") e = std::stod(value);
            else if (key == "L") L = std::stod(value);
            else if (key == "attitude_angle_deg") attitude_angle_deg = std::stod(value);
            else if (key == "load_angle_deg") load_angle_deg = std::stod(value);
            else if (key == "n_theta_global") n_theta_global = std::stoi(value);
            else if (key == "n_z_global") n_z_global = std::stoi(value);
            else if (key == "end_t") end_t = std::stod(value);
            else if (key == "dt") dt = std::stod(value);
            else if (key == "write_interval") write_interval = std::stod(value);
            else if (key == "solution_mode") solution_mode = solution_mode_from_config(value);
            else if (key == "pressure_time_method") pressure_time_method = time_stepping_from_config(value);
            else if (key == "motion_time_method") motion_time_method = time_stepping_from_config(value);
            else if (key == "temperature_time_method") temperature_time_method = time_stepping_from_config(value);
            else if (key == "omega") omega = std::stod(value);
            else if (key == "mu") mu = std::stod(value);
            else if (key == "rho") rho = std::stod(value);
            else if (key == "p_cav") p_cav = std::stod(value);
            else if (key == "bulk_modulus") bulk_modulus = std::stod(value);
            else if (key == "tilt_slope_x") tilt_slope_x = std::stod(value);
            else if (key == "tilt_slope_y") tilt_slope_y = std::stod(value);
            else if (key == "temperature_model") temperature_model = temperature_model_from_config(value);
            else if (key == "temperature_initial") temperature_initial = std::stod(value);
            else if (key == "temperature_reference") temperature_reference = std::stod(value);
            else if (key == "journal_wall_temperature") journal_wall_temperature = std::stod(value);
            else if (key == "bearing_wall_temperature") bearing_wall_temperature = std::stod(value);
            else if (key == "rho_cp") rho_cp = std::stod(value);
            else if (key == "thermal_conductivity") thermal_conductivity = std::stod(value);
            else if (key == "journal_heat_transfer") journal_heat_transfer = std::stod(value);
            else if (key == "bearing_heat_transfer") bearing_heat_transfer = std::stod(value);
            else if (key == "motion_model") motion_model = motion_model_from_config(value);
            else if (key == "bearing_initial_from_attitude") {
                if (is_true(value)) bearing_initial_from_attitude = true;
                else if (is_false(value)) bearing_initial_from_attitude = false;
            }
            else if (key == "bearing_initial_x") bearing_initial_x = std::stod(value);
            else if (key == "bearing_initial_y") bearing_initial_y = std::stod(value);
            else if (key == "bearing_initial_z") bearing_initial_z = std::stod(value);
            else if (key == "bearing_initial_vx") bearing_initial_vx = std::stod(value);
            else if (key == "bearing_initial_vy") bearing_initial_vy = std::stod(value);
            else if (key == "bearing_initial_vz") bearing_initial_vz = std::stod(value);
            else if (key == "bearing_mass") bearing_mass = std::stod(value);
            else if (key == "bearing_stiffness_x") bearing_stiffness_x = std::stod(value);
            else if (key == "bearing_stiffness_y") bearing_stiffness_y = std::stod(value);
            else if (key == "bearing_stiffness_z") bearing_stiffness_z = std::stod(value);
            else if (key == "bearing_damping_x") bearing_damping_x = std::stod(value);
            else if (key == "bearing_damping_y") bearing_damping_y = std::stod(value);
            else if (key == "bearing_damping_z") bearing_damping_z = std::stod(value);
            else if (key == "external_load_x") external_load_x = std::stod(value);
            else if (key == "external_load_y") external_load_y = std::stod(value);
            else if (key == "external_load_z") external_load_z = std::stod(value);
            else if (key == "min_film_thickness") min_film_thickness = std::stod(value);
            else if (key == "stop_on_nonpositive_film") {
                if (is_true(value)) stop_on_nonpositive_film = true;
                else if (is_false(value)) stop_on_nonpositive_film = false;
            }
            else if (key == "bc_z_south_type") {
                if (value == "DIRICHLET") bc_z_south_type = BCType::DIRICHLET;
                else if (value == "NEUMANN") bc_z_south_type = BCType::NEUMANN;
                else if (value == "INLET_OUTLET") bc_z_south_type = BCType::INLET_OUTLET;
            }
            else if (key == "bc_z_south_val") bc_z_south_val = std::stod(value);
            else if (key == "bc_z_south_theta") bc_z_south_theta = std::stod(value);
            else if (key == "bc_z_north_type") {
                if (value == "DIRICHLET") bc_z_north_type = BCType::DIRICHLET;
                else if (value == "NEUMANN") bc_z_north_type = BCType::NEUMANN;
                else if (value == "INLET_OUTLET") bc_z_north_type = BCType::INLET_OUTLET;
            }
            else if (key == "bc_z_north_val") bc_z_north_val = std::stod(value);
            else if (key == "bc_z_north_theta") bc_z_north_theta = std::stod(value);
            else if (key == "cavitation_model") {
                if (value == "GUMBEL") cavitation_model = CavitationModel::GUMBEL;
                else if (value == "ELROD_ADAMS") cavitation_model = CavitationModel::ELROD_ADAMS;
                else if (value == "FULL_SOMMERFELD") cavitation_model = CavitationModel::FULL_SOMMERFELD;
            }
            else if (key == "max_outer_iters") max_outer_iters = std::stoi(value);
            else if (key == "outer_tol") outer_tol = std::stod(value);
            else if (key == "theta_min") theta_min = std::stod(value);
            else if (key == "log_outer_iters") {
                if (is_true(value)) log_outer_iters = true;
                else if (is_false(value)) log_outer_iters = false;
            }
            else if (key == "output_dir") output_dir = value;
            else if (key == "filename_prefix") filename_prefix = value;
            else if (key == "output_write_3d") {
                if (is_true(value)) output_write_3d = true;
                else if (is_false(value)) output_write_3d = false;
            }
            else if (key == "output_write_flat") {
                if (is_true(value)) output_write_flat = true;
                else if (is_false(value)) output_write_flat = false;
            }
            else if (key == "output_fields") output_fields = parse_output_fields(value);
            else if (key == "inlet_circular") {
                InletConfig inlet;
                inlet.type = InletConfig::Type::CIRCULAR;
                std::stringstream ss(value);
                ss >> inlet.theta >> inlet.z >> inlet.size >> inlet.p_supply;
                inlets.push_back(inlet);
            }
            else if (key == "inlet_groove") {
                InletConfig inlet;
                inlet.type = InletConfig::Type::GROOVE;
                std::stringstream ss(value);
                ss >> inlet.theta >> inlet.size >> inlet.p_supply;
                inlets.push_back(inlet);
            }
        }
    }

    void save_to_file(const std::filesystem::path& path) const {
        std::ofstream file(path, std::ios::trunc);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open config file for writing: " + path.string());
        }
        file << to_config_text();
    }

    bool output_field_enabled(const std::string& name) const {
        return std::find(output_fields.begin(), output_fields.end(), normalise_key(name)) != output_fields.end();
    }

    double initial_pressure() const {
        return inlets.empty() ? p_cav : inlets.front().p_supply;
    }

    std::string output_fields_text() const {
        std::ostringstream out;
        for (size_t i = 0; i < output_fields.size(); ++i) {
            if (i > 0) out << ", ";
            out << output_fields[i];
        }
        return out.str();
    }

    static std::vector<std::string> default_output_fields() {
        return {
            "pressure",
            "film_content",
            "h",
            "rho",
            "temperature",
            "heat_generation",
            "inlet_indicator",
            "velocity",
            "pressure_force_x",
            "pressure_force_y",
            "pressure_force_z",
            "viscous_force_x",
            "viscous_force_y",
            "viscous_force_z",
            "fluid_force_x",
            "fluid_force_y",
            "fluid_force_z",
            "external_load_x",
            "external_load_y",
            "external_load_z",
            "bearing_x",
            "bearing_y",
            "bearing_z",
            "bearing_attitude_angle",
            "friction_torque"};
    }

    std::string to_config_text() const {
        std::ostringstream out;
        out << std::setprecision(17);

        out << "# Geometry\n"
            << "R = " << R << "\n"
            << "c = " << c << "\n"
            << "e = " << e << "\n"
            << "L = " << L << "\n"
            << "attitude_angle_deg = " << attitude_angle_deg << "\n"
            << "load_angle_deg = " << load_angle_deg << "\n"
            << "tilt_slope_x = " << tilt_slope_x << "\n"
            << "tilt_slope_y = " << tilt_slope_y << "\n\n";

        out << "# Grid\n"
            << "n_theta_global = " << n_theta_global << "\n"
            << "n_z_global = " << n_z_global << "\n\n";

        out << "# Time\n"
            << "end_t = " << end_t << "\n"
            << "dt = " << dt << "\n"
            << "write_interval = " << write_interval << "\n"
            << "# solution_mode options: TRANSIENT, STEADY_STATE\n"
            << "solution_mode = " << to_config_value(solution_mode) << "\n"
            << "# Options: EULER_EXPLICIT, EULER_IMPLICIT, CRANK_NICOLSON, RK2, RK4\n"
            << "pressure_time_method = " << to_config_value(pressure_time_method) << "\n"
            << "motion_time_method = " << to_config_value(motion_time_method) << "\n"
            << "temperature_time_method = " << to_config_value(temperature_time_method) << "\n\n";

        out << "# Physics\n"
            << "omega = " << omega << "\n"
            << "mu = " << mu << "\n"
            << "rho = " << rho << "\n"
            << "p_cav = " << p_cav << "\n"
            << "bulk_modulus = " << bulk_modulus << "\n\n";

        out << "# Energy Equation\n"
            << "# temperature_model options: ISOTHERMAL, ENERGY_EQUATION\n"
            << "temperature_model = " << to_config_value(temperature_model) << "\n"
            << "temperature_initial = " << temperature_initial << "\n"
            << "temperature_reference = " << temperature_reference << "\n"
            << "journal_wall_temperature = " << journal_wall_temperature << "\n"
            << "bearing_wall_temperature = " << bearing_wall_temperature << "\n"
            << "rho_cp = " << rho_cp << "\n"
            << "thermal_conductivity = " << thermal_conductivity << "\n"
            << "journal_heat_transfer = " << journal_heat_transfer << "\n"
            << "bearing_heat_transfer = " << bearing_heat_transfer << "\n\n";

        out << "# Bearing Motion\n"
            << "# motion_model options: STATIC, MOVING_BEARING\n"
            << "motion_model = " << to_config_value(motion_model) << "\n"
            << "bearing_initial_from_attitude = " << (bearing_initial_from_attitude ? "true" : "false") << "\n"
            << "bearing_initial_x = " << bearing_initial_x << "\n"
            << "bearing_initial_y = " << bearing_initial_y << "\n"
            << "bearing_initial_z = " << bearing_initial_z << "\n"
            << "bearing_initial_vx = " << bearing_initial_vx << "\n"
            << "bearing_initial_vy = " << bearing_initial_vy << "\n"
            << "bearing_initial_vz = " << bearing_initial_vz << "\n"
            << "bearing_mass = " << bearing_mass << "\n"
            << "bearing_stiffness_x = " << bearing_stiffness_x << "\n"
            << "bearing_stiffness_y = " << bearing_stiffness_y << "\n"
            << "bearing_stiffness_z = " << bearing_stiffness_z << "\n"
            << "bearing_damping_x = " << bearing_damping_x << "\n"
            << "bearing_damping_y = " << bearing_damping_y << "\n"
            << "bearing_damping_z = " << bearing_damping_z << "\n"
            << "external_load_x = " << external_load_x << "\n"
            << "external_load_y = " << external_load_y << "\n"
            << "external_load_z = " << external_load_z << "\n"
            << "min_film_thickness = " << min_film_thickness << "\n"
            << "stop_on_nonpositive_film = " << (stop_on_nonpositive_film ? "true" : "false") << "\n\n";

        out << "# Axial Boundary Conditions (Outlets)\n"
            << "# Options: DIRICHLET, NEUMANN, INLET_OUTLET\n"
            << "bc_z_south_type = " << to_config_value(bc_z_south_type) << "\n"
            << "bc_z_south_val = " << bc_z_south_val << "\n"
            << "bc_z_south_theta = " << bc_z_south_theta << "\n"
            << "bc_z_north_type = " << to_config_value(bc_z_north_type) << "\n"
            << "bc_z_north_val = " << bc_z_north_val << "\n"
            << "bc_z_north_theta = " << bc_z_north_theta << "\n\n";

        out << "# Cavitation\n"
            << "# Options: FULL_SOMMERFELD, GUMBEL, ELROD_ADAMS\n"
            << "cavitation_model = " << to_config_value(cavitation_model) << "\n"
            << "max_outer_iters = " << max_outer_iters << "\n"
            << "outer_tol = " << outer_tol << "\n"
            << "theta_min = " << theta_min << "\n"
            << "log_outer_iters = " << (log_outer_iters ? "true" : "false") << "\n\n";

        out << "# Inlets\n"
            << "# inlet_circular = theta(deg) z(m) radius(m) p_supply(Pa)\n"
            << "# inlet_groove   = theta(deg) width(deg) p_supply(Pa)\n";
        if (inlets.empty()) {
            out << "# inlet_groove = 90.0 10.0 2e5\n";
        } else {
            for (const auto& inlet : inlets) {
                if (inlet.type == InletConfig::Type::CIRCULAR) {
                    out << "inlet_circular = " << inlet.theta << " " << inlet.z << " "
                        << inlet.size << " " << inlet.p_supply << "\n";
                } else {
                    out << "inlet_groove = " << inlet.theta << " "
                        << inlet.size << " " << inlet.p_supply << "\n";
                }
            }
        }
        out << "\n";

        out << "# Output\n"
            << "# output_fields accepts comma-separated field names.\n"
            << "# Options: pressure, film_content, h, rho, temperature, heat_generation,\n"
            << "#          inlet_indicator, velocity,\n"
            << "#          pressure_force_x, pressure_force_y, pressure_force_z,\n"
            << "#          load_x, load_y, load_z (legacy pressure-force aliases),\n"
            << "#          viscous_force_x, viscous_force_y, viscous_force_z,\n"
            << "#          fluid_force_x, fluid_force_y, fluid_force_z,\n"
            << "#          external_load_x, external_load_y, external_load_z,\n"
            << "#          bearing_x, bearing_y, bearing_z, bearing_attitude_angle,\n"
            << "#          friction_torque\n"
            << "output_dir = " << output_dir << "\n"
            << "filename_prefix = " << filename_prefix << "\n"
            << "output_write_3d = " << (output_write_3d ? "true" : "false") << "\n"
            << "output_write_flat = " << (output_write_flat ? "true" : "false") << "\n"
            << "output_fields = " << output_fields_text() << "\n";

        return out.str();
    }

private:
    static std::string normalise_key(std::string text) {
        trim(text);
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });
        return text;
    }

    static std::vector<std::string> parse_output_fields(std::string value) {
        std::replace(value.begin(), value.end(), ',', ' ');
        std::replace(value.begin(), value.end(), ';', ' ');
        std::replace(value.begin(), value.end(), '|', ' ');

        std::vector<std::string> fields;
        std::stringstream input(value);
        std::string token;
        while (input >> token) {
            const std::string key = normalise_key(token);
            if (!key.empty() && std::find(fields.begin(), fields.end(), key) == fields.end()) {
                fields.push_back(key);
            }
        }
        return fields;
    }

    static void trim(std::string& text) {
        const auto first = text.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            text.clear();
            return;
        }
        const auto last = text.find_last_not_of(" \t\r\n");
        text = text.substr(first, last - first + 1);
    }

    static bool is_true(const std::string& value) {
        return value == "true" || value == "TRUE" || value == "True" || value == "1";
    }

    static bool is_false(const std::string& value) {
        return value == "false" || value == "FALSE" || value == "False" || value == "0";
    }
};
