#include "bubble_config.hpp"

void BubbleConfig::load_from_file(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open bubble config file '" << path
                  << "'. Using defaults.\n";
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Strip comments
        size_t comment = line.find('#');
        if (comment != std::string::npos) line = line.substr(0, comment);

        // Trim whitespace
        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(" \t\r\n"));
            s.erase(s.find_last_not_of(" \t\r\n") + 1);
        };
        trim(line);
        if (line.empty()) continue;

        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key   = line.substr(0, eq);
        std::string value = line.substr(eq + 1);
        trim(key);
        trim(value);
        if (key.empty() || value.empty()) continue;

        // Geometry
        if      (key == "R_bubble")  R_bubble  = std::stod(value);
        else if (key == "L_box")     L_box     = std::stod(value);
        else if (key == "n_u")       n_u       = std::stoi(value);
        else if (key == "n_v")       n_v       = std::stoi(value);
        // Film thickness
        else if (key == "h_initial") h_initial = std::stod(value);
        else if (key == "h_rim")     h_rim     = std::stod(value);
        else if (key == "h_min")     h_min     = std::stod(value);
        else if (key == "h_rupture") h_rupture = std::stod(value);
        // Liquid properties
        else if (key == "rho_l")   rho_l   = std::stod(value);
        else if (key == "mu")      mu      = std::stod(value);
        else if (key == "sigma_0") sigma_0 = std::stod(value);
        else if (key == "gamma_T") gamma_T = std::stod(value);
        else if (key == "k_l")     k_l     = std::stod(value);
        else if (key == "c_p")     c_p     = std::stod(value);
        // Gas and evaporation
        else if (key == "rho_gas")   rho_gas   = std::stod(value);
        else if (key == "D_AB")      D_AB      = std::stod(value);
        else if (key == "omega_inf") omega_inf = std::stod(value);
        else if (key == "T_ref")     T_ref     = std::stod(value);
        else if (key == "T_ambient") T_ambient = std::stod(value);
        else if (key == "T_rim")     T_rim     = std::stod(value);
        else if (key == "delta_H_v") delta_H_v = std::stod(value);
        else if (key == "h_conv")    h_conv    = std::stod(value);
        // Disjoining pressure
        else if (key == "A_hamaker") A_hamaker = std::stod(value);
        // Gravity
        else if (key == "g_x") g_x = std::stod(value);
        else if (key == "g_y") g_y = std::stod(value);
        else if (key == "g_z") g_z = std::stod(value);
        // Time stepping
        else if (key == "dt")     dt     = std::stod(value);
        else if (key == "dt_max") dt_max = std::stod(value);
        else if (key == "end_t")  end_t  = std::stod(value);
        else if (key == "cfl")    cfl    = std::stod(value);
        // SIMPLE solver
        else if (key == "alpha_p")          alpha_p          = std::stod(value);
        else if (key == "alpha_u")          alpha_u          = std::stod(value);
        else if (key == "max_simple_iters") max_simple_iters = std::stoi(value);
        else if (key == "simple_tol")       simple_tol       = std::stod(value);
        else if (key == "advection_scheme") {
            if      (value == "UPWIND")       advection_scheme = ConvectionScheme::UPWIND;
            else if (value == "TVD_VANLEER")  advection_scheme = ConvectionScheme::TVD_VANLEER;
            else if (value == "TVD_MINMOD")   advection_scheme = ConvectionScheme::TVD_MINMOD;
        }
        // Physics toggles
        else if (key == "enable_evaporation") {
            enable_evaporation = (value == "true" || value == "1");
        }
        else if (key == "enable_marangoni") {
            enable_marangoni = (value == "true" || value == "1");
        }
        else if (key == "enable_disjoining") {
            enable_disjoining = (value == "true" || value == "1");
        }
        else if (key == "enable_gravity") {
            enable_gravity = (value == "true" || value == "1");
        }
        // Output
        else if (key == "output_dir")     output_dir     = value;
        else if (key == "write_interval") write_interval = std::stod(value);
    }
}
