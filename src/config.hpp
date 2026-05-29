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

    // Physics
    double omega = 100.0;
    double mu = 0.01;
    double rho = 900.0;
    double p_cav = 0.0;
    double bulk_modulus = 1e9;

    // Misalignment (tilt)
    double tilt_slope_x = 0.0; // Slope in x-direction (m/m)
    double tilt_slope_y = 0.0; // Slope in y-direction (m/m)

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
            else if (key == "omega") omega = std::stod(value);
            else if (key == "mu") mu = std::stod(value);
            else if (key == "rho") rho = std::stod(value);
            else if (key == "p_cav") p_cav = std::stod(value);
            else if (key == "bulk_modulus") bulk_modulus = std::stod(value);
            else if (key == "tilt_slope_x") tilt_slope_x = std::stod(value);
            else if (key == "tilt_slope_y") tilt_slope_y = std::stod(value);
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
            "inlet_indicator",
            "velocity",
            "load_x",
            "load_y",
            "load_z",
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
            << "write_interval = " << write_interval << "\n\n";

        out << "# Physics\n"
            << "omega = " << omega << "\n"
            << "mu = " << mu << "\n"
            << "rho = " << rho << "\n"
            << "p_cav = " << p_cav << "\n"
            << "bulk_modulus = " << bulk_modulus << "\n\n";

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
            << "# Options: pressure, film_content, h, rho, inlet_indicator, velocity,\n"
            << "#          load_x, load_y, load_z, friction_torque\n"
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
