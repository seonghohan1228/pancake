#include "mesh.hpp"

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
