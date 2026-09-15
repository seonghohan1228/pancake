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

void compute_moving_bearing(Field& h, Field& dh_dt, const Mesh& mesh,
                            const SimulationConfig& cfg,
                            const JournalMotion::BearingState& state) {
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();

    for (int j = 0; j < mesh.n_z_local; ++j) {
        const double z = (j + 0.5) * d_z;
        const double dz = z - 0.5 * cfg.L;

        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta = (mesh.offset_theta + i + 0.5) * d_theta;
            const double cos_t = std::cos(theta);
            const double sin_t = std::sin(theta);

            h(i, j) = cfg.c
                + state.x * cos_t + state.y * sin_t
                - cfg.tilt_slope_x * dz * cos_t
                - cfg.tilt_slope_y * dz * sin_t;
            dh_dt(i, j) = state.vx * cos_t + state.vy * sin_t;
        }
    }
}

void compute_inlet_indicator(Field& indicator, const Mesh& mesh, const SimulationConfig& cfg) {
    indicator.fill(0.0);
    if (cfg.inlets.empty()) return;

    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            indicator(i,j) = mesh.inlet_at(cfg,i,j) ? 1.0 : 0.0;
}

}  // namespace FilmThickness
