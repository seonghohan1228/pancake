#include "mesh.hpp"

Mesh::Mesh(const SimulationConfig& cfg, MPI_Comm comm)
    : R(cfg.R), L(cfg.L) {

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    n_theta_global = cfg.n_theta_global;
    n_z_global = cfg.n_z_global;

    // --- Domain Decomposition (1D Split in Theta) ---
    // We split the CELLS. Points will be derived from cell offsets.
    n_theta_local = n_theta_global / size;
    int remainder = n_theta_global % size;
    offset_theta = rank * n_theta_local;

    // Distribute remainder to first few ranks to balance load
    if (rank < remainder) {
        n_theta_local++;
        offset_theta += rank;
    } else {
        offset_theta += remainder;
    }

    // Z-direction is not split in this version (each rank has full height)
    n_z_local = n_z_global;
    offset_z = 0;
}
