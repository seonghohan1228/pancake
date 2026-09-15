#include "mesh.hpp"

#include <cmath>
#include <numbers>

Mesh::Mesh(const SimulationConfig& cfg, MPI_Comm comm)
    : R(cfg.R), L(cfg.L) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    n_theta_global = cfg.n_theta_global;
    n_z_global = cfg.n_z_global;
    n_z_local = n_z_global;
    offset_z = 0;

    // Domain decomposition (1D in theta direction)
    n_theta_local = n_theta_global / size;
    int remainder = n_theta_global % size;
    offset_theta = rank * n_theta_local;

    if (rank < remainder) {
        n_theta_local++;
        offset_theta += rank;
    } else {
        offset_theta += remainder;
    }
}

bool Mesh::inside_inlet(const InletConfig& inlet, int i, int j) const {
    const double angle = (offset_theta + i + 0.5) * get_d_theta();
    const double delta = std::remainder(angle - inlet.theta * std::numbers::pi / 180.0,
                                        2.0 * std::numbers::pi);
    if (inlet.type == InletConfig::Type::GROOVE)
        return std::abs(delta) <= inlet.size * std::numbers::pi / 360.0;
    const double dz = (j + 0.5) * get_d_z() - inlet.z;
    return R * R * delta * delta + dz * dz <= inlet.size * inlet.size;
}

const InletConfig* Mesh::inlet_at(const SimulationConfig& cfg, int i, int j) const {
    for (auto it = cfg.inlets.rbegin(); it != cfg.inlets.rend(); ++it)
        if (inside_inlet(*it, i, j)) return &*it;
    return nullptr;
}
