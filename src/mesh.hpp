#pragma once
#include <mpi.h>
#include <cmath>
#include "config.hpp"

class Mesh {
public:
    int rank, size;

    // Global Dimensions
    int n_theta_global;
    int n_z_global;

    // Local Dimensions (Owned by this rank)
    int n_theta_local;
    int n_z_local;

    // Offsets (Where this rank starts in global index)
    int offset_theta;
    int offset_z;

    // Geometry parameters
    double R, L;

    Mesh(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);

    // Returns the number of CELLS owned by this rank
    size_t local_cell_count() const { return n_theta_local * n_z_local; }

    // Coordinate Helpers
    double get_d_theta() const { return (2.0 * M_PI) / n_theta_global; }
    double get_d_z() const { return L / n_z_global; }
};
