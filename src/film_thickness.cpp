#include "film_thickness.hpp"

#include <cmath>

namespace FilmThickness {

void compute_static(Field& h, const Mesh& mesh, const SimulationConfig& cfg) {
    const double psi_rad = cfg.attitude_angle_deg * (M_PI / 180.0);
    const double d_theta = mesh.get_d_theta();

    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            double theta_coord = (mesh.offset_theta + i + 0.5) * d_theta;
            h(i, j) = cfg.c - cfg.e * std::cos(theta_coord - psi_rad);
        }
    }
}

}  // namespace FilmThickness
