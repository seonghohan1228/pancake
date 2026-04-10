#pragma once

#include <mpi.h>
#include <vector>

#include "field.hpp"
#include "mesh.hpp"

class Communicator {
private:
    const Mesh& mesh;

    // Buffers for packing/unpacking non-contiguous ghost data (size: n_z_local * n_ghost)
    std::vector<double> send_buf_left, recv_buf_left;
    std::vector<double> send_buf_right, recv_buf_right;

public:
    Communicator(const Mesh& mesh_ref);

    // Update ghost cells for ALL fields in the container
    void update_ghosts(Fields& fields);

    // Update ghost cells for a SINGLE field
    void update_ghosts(Field& field);
};
