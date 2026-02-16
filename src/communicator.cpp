#include "communicator.hpp"

Communicator::Communicator(const Mesh& mesh_ref) : mesh(mesh_ref) {
    // Max buffer size needed: Z-height * max expected ghost layers (e.g. 2)
    // We assume 2 layers for safety, but will resize dynamically if needed in update
}

void Communicator::update_ghosts(Fields& fields) {
    for (auto& [name, field] : fields) {
        update_ghosts(*field);
    }
}

void Communicator::update_ghosts(Field& field) {
    int ng = field.ng;
    int nz = field.n_z_phys;
    int n_theta = field.n_theta_phys;
    int buffer_size = nz * ng;

    // Resize buffers if needed
    if (send_buf_left.size() < buffer_size) {
        send_buf_left.resize(buffer_size);
        recv_buf_left.resize(buffer_size);
        send_buf_right.resize(buffer_size);
        recv_buf_right.resize(buffer_size);
    }

    // --- 1. Pack Data ---
    // Send the "Left Edge" to the Left Neighbor
    // and the "Right Edge" to the Right Neighbor.

    int idx = 0;
    for (int j = 0; j < nz; ++j) {
        // Pack Left Edge (Physical indices 0 to ng-1)
        for (int k = 0; k < ng; ++k) {
            send_buf_left[idx + k] = field(k, j);
        }
        // Pack Right Edge (Physical indices N-ng to N-1)
        for (int k = 0; k < ng; ++k) {
            send_buf_right[idx + k] = field(n_theta - ng + k, j);
        }
        idx += ng;
    }

    // --- 2. MPI Exchange ---
    // Periodic boundaries:
    // Left neighbor of Rank 0 is Rank Size-1
    // Right neighbor of Rank Size-1 is Rank 0
    int left_rank = (mesh.rank - 1 + mesh.size) % mesh.size;
    int right_rank = (mesh.rank + 1) % mesh.size;

    MPI_Status status;

    // Send Left, Receive from Right (filling Right Ghost)
    MPI_Sendrecv(send_buf_left.data(), buffer_size, MPI_DOUBLE, left_rank, 0,
                 recv_buf_right.data(), buffer_size, MPI_DOUBLE, right_rank, 0,
                 MPI_COMM_WORLD, &status);

    // Send Right, Receive from Left (filling Left Ghost)
    MPI_Sendrecv(send_buf_right.data(), buffer_size, MPI_DOUBLE, right_rank, 1,
                 recv_buf_left.data(), buffer_size, MPI_DOUBLE, left_rank, 1,
                 MPI_COMM_WORLD, &status);

    // --- 3. Unpack Data ---
    idx = 0;
    for (int j = 0; j < nz; ++j) {
        // Unpack into Left Ghost (indices -ng to -1)
        // Note: recv_buf_left contains data sent by Left Neighbor's Right Edge
        for (int k = 0; k < ng; ++k) {
            // field(-ng + k, j)
            field(-ng + k, j) = recv_buf_left[idx + k];
        }

        // Unpack into Right Ghost (indices N to N+ng-1)
        for (int k = 0; k < ng; ++k) {
            // field(n_theta + k, j)
            field(n_theta + k, j) = recv_buf_right[idx + k];
        }
        idx += ng;
    }
}
