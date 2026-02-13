#pragma once
#include <string>

struct SimulationConfig {
    // TODO: Add a parser to read these from a file (e.g., .json or .ini)

    // Geometry
    double R = 0.01;
    double c = 0.001;
    double e = 0.0001;
    double L = 0.02;

    // Grid Dimensions (Global)
    int n_theta_global = 30;
    int n_z_global = 10;

    // Time Integration
    double end_t = 1.0;
    double dt = 0.1;
    double write_interval = 0.1;

    // Output Settings
    std::string output_dir = "results";
    std::string filename_prefix = "solution";
};
