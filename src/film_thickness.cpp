#include "film_thickness.hpp"

#include <cmath>

namespace FilmThickness {

void compute_static(Field& h, const Mesh& mesh, const SimulationConfig& cfg) {
    const double psi_rad = cfg.attitude_angle_deg * (M_PI / 180.0);
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();

    // Eccentricity components at the center (z = L/2)
    const double ex_c = cfg.e * std::cos(psi_rad);
    const double ey_c = cfg.e * std::sin(psi_rad);

    for (int j = 0; j < mesh.n_z_local; ++j) {
        const double z = (j + 0.5) * d_z;
        const double dz = z - 0.5 * cfg.L;

        // Misaligned eccentricity components
        const double ex = ex_c + cfg.tilt_slope_x * dz;
        const double ey = ey_c + cfg.tilt_slope_y * dz;

        for (int i = 0; i < mesh.n_theta_local; ++i) {
            double theta = (mesh.offset_theta + i + 0.5) * d_theta;
            // h = c - ex*cos(theta) - ey*sin(theta)
            h(i, j) = cfg.c - ex * std::cos(theta) - ey * std::sin(theta);
        }
    }
}

void compute_inlet_indicator(Field& indicator, const Mesh& mesh, const SimulationConfig& cfg) {
    indicator.fill(0.0);
    if (cfg.inlets.empty()) return;

    const int n_theta = mesh.n_theta_local;
    const int n_z = mesh.n_z_local;
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double R = mesh.R;

    for (const auto& inlet : cfg.inlets) {
        const double theta_in_rad = inlet.theta * M_PI / 180.0;
        const double size_rad = inlet.size * M_PI / 180.0;
        
        for (int i = 0; i < n_theta; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - theta_in_rad);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            for (int j = 0; j < n_z; ++j) {
                const double z_c = (j + 0.5) * d_z;
                bool inside = false;

                if (inlet.type == InletConfig::Type::CIRCULAR) {
                    const double dist_sq = (R * dth) * (R * dth) + (z_c - inlet.z) * (z_c - inlet.z);
                    if (dist_sq <= inlet.size * inlet.size) inside = true;
                } else if (inlet.type == InletConfig::Type::GROOVE) {
                    if (dth <= 0.5 * size_rad) inside = true;
                }

                if (inside) indicator(i, j) = 1.0;
            }
        }
    }
}

}  // namespace FilmThickness
