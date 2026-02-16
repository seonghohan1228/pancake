#include "field.hpp"

Field::Field(std::string name, const Mesh& mesh, int ghost_layers, GridLocation loc)
    : name(name), location(loc), ng(ghost_layers) {

    n_theta_phys = mesh.n_theta_local;
    n_z_phys = mesh.n_z_local;

    // Calculate total size including ghost layers
    int n_theta_total = n_theta_phys + 2 * ng;
    int n_z_total = n_z_phys + 2 * ng;

    stride_theta = n_theta_total;

    // Allocate and initialize to 0.0
    data.resize(n_theta_total * n_z_total, 0.0);
}

void Field::fill(double value) {
    std::fill(data.begin(), data.end(), value);
}
