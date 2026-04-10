#pragma once

#include <math.h>
#include <mpi.h>

#include "config.hpp"

class Mesh {
public:
    int rank, size;
    int n_theta_global, n_z_global;
    int n_theta_local, n_z_local;
    int offset_theta, offset_z;
    double R, L;

    Mesh(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);

    double get_d_theta() const { return (2.0 * M_PI) / n_theta_global; }
    double get_d_z() const { return L / n_z_global; }
    double cell_volume() const { return (R * get_d_theta()) * get_d_z(); }
};
