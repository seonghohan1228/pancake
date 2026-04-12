#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

enum class CavitationModel { GUMBEL, ELROD_ADAMS, FULL_SOMMERFELD };
enum class BCType { DIRICHLET, NEUMANN, INLET_OUTLET };

struct InletConfig {
    enum class Type { CIRCULAR, GROOVE } type;
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
    double bulk_modulus = 1e5;

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
    int    max_outer_iters = 50;    // Elrod-Adams outer (flag-update) iteration limit
    double outer_tol       = 1e-6;  // Convergence tolerance on max |Δθ| between outer iters
    double theta_min       = 1e-6;  // Minimum film content (prevents unphysical θ ≤ 0)
    bool   log_outer_iters = true;  // Whether to print outer iteration progress

    // Output
    std::string output_dir = "results";
    std::string filename_prefix = "solution";

    void load_from_file(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Warning: Could not open config file '" << path << "'. Using defaults.\n";
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            // Remove comments
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);

            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            if (line.empty()) continue;

            // Simple split by '='
            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            std::string key = line.substr(0, eq_pos);
            std::string value = line.substr(eq_pos + 1);

            // Trim key and value
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            if (key.empty() || value.empty()) continue;

            if (key == "R") R = std::stod(value);
            else if (key == "c") c = std::stod(value);
            else if (key == "e") e = std::stod(value);
            else if (key == "L") L = std::stod(value);
            else if (key == "attitude_angle_deg") attitude_angle_deg = std::stod(value);
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
            else if (key == "max_outer_iters") max_outer_iters = std::stoi(value);            else if (key == "outer_tol") outer_tol = std::stod(value);
            else if (key == "theta_min") theta_min = std::stod(value);
            else if (key == "log_outer_iters") {
                if (value == "true" || value == "1") log_outer_iters = true;
                else if (value == "false" || value == "0") log_outer_iters = false;
            }
            else if (key == "output_dir") output_dir = value;
            else if (key == "filename_prefix") filename_prefix = value;
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
};
