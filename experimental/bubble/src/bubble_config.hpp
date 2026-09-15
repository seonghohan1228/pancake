#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

enum class ConvectionScheme { UPWIND, TVD_VANLEER, TVD_MINMOD };

/// Configuration for the hemispherical soap bubble thin-film solver.
/// All parameters have physically reasonable defaults for a 1 cm radius bubble
/// in ambient air at 20 °C. Load from a text file with load_from_file().
struct BubbleConfig {
    // --- Geometry ---
    double R_bubble = 0.01;    // Bubble radius [m]
    double L_box    = 0.012;   // Half-width of 2D square stereographic domain [m]
    int    n_u = 64;           // Grid cells in u direction
    int    n_v = 64;           // Grid cells in v direction

    // --- Film thickness ---
    double h_initial = 1.0e-6; // Initial uniform film thickness [m]
    double h_rim     = 1.0e-6; // Rim Dirichlet boundary condition on h [m]
    double h_min     = 1.0e-8; // Minimum allowed film thickness (clamping floor) [m]
    double h_rupture = 1.0e-8; // Film thickness below which rupture is declared [m]

    // --- Liquid properties ---
    double rho_l   = 1000.0;   // Liquid density [kg/m^3]
    double mu      = 1.0e-3;   // Dynamic viscosity [Pa·s]
    double sigma_0 = 0.025;    // Surface tension at T_ref [N/m]
    double gamma_T = -1.7e-4;  // dσ/dT [N/(m·K)]
    double k_l     = 0.6;      // Thermal conductivity [W/(m·K)]
    double c_p     = 4186.0;   // Specific heat capacity [J/(kg·K)]

    // --- Gas and evaporation properties ---
    double rho_gas   = 1.2;    // Ambient air density [kg/m^3]
    double D_AB      = 2.5e-5; // Water vapour diffusivity in air [m^2/s]
    double omega_inf = 0.01;   // Ambient humidity mass fraction [-]
    double T_ref     = 293.15; // Reference temperature [K]
    double T_ambient = 293.15; // Ambient air temperature [K]
    double T_rim     = 293.15; // Rim temperature Dirichlet BC [K]
    double delta_H_v = 2.26e6; // Latent heat of vaporisation [J/kg]
    double h_conv    = 10.0;   // Convective heat transfer coefficient [W/(m^2·K)]

    // --- Disjoining pressure ---
    double A_hamaker = 1.0e-20; // Hamaker constant [J]

    // --- Gravity ---
    double g_x = 0.0;
    double g_y = 0.0;
    double g_z = -9.81;        // Gravitational acceleration components [m/s^2]

    // --- Time stepping ---
    double dt          = 1.0e-6; // Initial timestep [s]
    double dt_max      = 1.0e-3; // Maximum allowed timestep [s]
    double end_t       = 10.0;   // Simulation end time [s]
    double cfl         = 0.5;    // CFL number for adaptive timestep

    // --- SIMPLE solver ---
    double alpha_p         = 0.3;  // Pressure under-relaxation factor
    double alpha_u         = 0.7;  // Velocity under-relaxation factor
    int    max_simple_iters = 50;  // Maximum SIMPLE iterations per timestep
    double simple_tol      = 1.0e-4; // SIMPLE convergence tolerance

    ConvectionScheme advection_scheme = ConvectionScheme::TVD_VANLEER;

    // --- Physics toggles ---
    bool enable_evaporation = true;
    bool enable_marangoni   = true;
    bool enable_disjoining  = true;
    bool enable_gravity     = true;

    // --- Output ---
    std::string output_dir    = "results";
    double      write_interval = 0.01; // Time between VTK snapshots [s]

    void load_from_file(const std::string& path);
};
