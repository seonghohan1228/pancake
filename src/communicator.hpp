#pragma once

#include <mpi.h>
#include <vector>

#include "field.hpp"
#include "stereo_mesh.hpp"

/// MPI ghost-cell exchange for the bubble solver on a StereoMesh.
///
/// Domain decomposition is 1D in u (same as bearing solver in theta).
/// u boundaries are OPEN (non-periodic):
///   - Left ghost of rank 0 receives zero-gradient fill (copied from physical cell 0).
///   - Right ghost of rank size-1 receives zero-gradient fill from the last cell.
///   - Interior rank boundaries exchange via MPI_Sendrecv.
/// v (j) direction is not decomposed; top/bottom ghost cells receive zero-gradient fill.
class Communicator {
public:
    Communicator(const StereoMesh& mesh);

    /// Exchange ghost cells for all fields in the container.
    void update_ghosts(Fields& fields);

    /// Exchange ghost cells for a single field.
    void update_ghosts(Field& field);

private:
    const StereoMesh& mesh;
    std::vector<double> send_buf_left, recv_buf_left;
    std::vector<double> send_buf_right, recv_buf_right;

    void fill_v_ghosts(Field& field);
};
