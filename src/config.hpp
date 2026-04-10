#pragma once

#include <string>

struct SimulationConfig {
    // Geometry
    double R = 0.01;
    double c = 0.001;
    double e = 0.0008;
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

    // Output
    std::string output_dir = "results";
    std::string filename_prefix = "solution";
};
