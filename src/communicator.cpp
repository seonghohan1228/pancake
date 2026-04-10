#include "communicator.hpp"

Communicator::Communicator(const Mesh& mesh_ref) : mesh(mesh_ref) {}

void Communicator::update_ghosts(Fields& fields) {
    for (auto& [name, field] : fields) {
        update_ghosts(*field);
    }
}

void Communicator::update_ghosts(Field& field) {
    int n_ghost = field.n_ghost;
    int n_z = field.n_z_phys;
    int n_theta = field.n_theta_phys;
    int buffer_size = n_z * n_ghost;

    // Resize buffers if needed
    if (send_buf_left.size() < buffer_size) {
        send_buf_left.resize(buffer_size);
        recv_buf_left.resize(buffer_size);
        send_buf_right.resize(buffer_size);
        recv_buf_right.resize(buffer_size);
    }

    // Pack data
    // Send the left edge to the left neighborand the right edge to the right neighbor.
    int idx = 0;
    for (int j = 0; j < n_z; ++j) {
        // Pack left edge (physical indices 0 to n_ghost-1)
        for (int k = 0; k < n_ghost; ++k) send_buf_left[idx + k] = field(k, j);
        // Pack right edge (physical indices N-n_ghost to N-1)
        for (int k = 0; k < n_ghost; ++k) send_buf_right[idx + k] = field(n_theta - n_ghost + k, j);
        idx += n_ghost;
    }

    // MPI exchange
    // Periodic boundaries:
    //   - Left neighbor of Rank 0 is Rank Size-1
    //   - Right neighbor of Rank Size-1 is Rank 0
    int left_rank = (mesh.rank - 1 + mesh.size) % mesh.size;
    int right_rank = (mesh.rank + 1) % mesh.size;

    MPI_Status status;

    // Send left, receive from right (filling right ghost)
    MPI_Sendrecv(send_buf_left.data(), buffer_size, MPI_DOUBLE, left_rank, 0,
                 recv_buf_right.data(), buffer_size, MPI_DOUBLE, right_rank, 0,
                 MPI_COMM_WORLD, &status);

    // Send right, receive from left (filling left ghost)
    MPI_Sendrecv(send_buf_right.data(), buffer_size, MPI_DOUBLE, right_rank, 1,
                 recv_buf_left.data(), buffer_size, MPI_DOUBLE, left_rank, 1,
                 MPI_COMM_WORLD, &status);

    // Unpack data
    idx = 0;
    for (int j = 0; j < n_z; ++j) {
        // Unpack into left ghost (indices -n_ghost to -1)
        // Note: recv_buf_left contains data sent by left neighbor's right edge
        for (int k = 0; k < n_ghost; ++k) {
            // field(-n_ghost + k, j)
            field(-n_ghost + k, j) = recv_buf_left[idx + k];
        }

        // Unpack into right ghost (indices N to N+n_ghost-1)
        for (int k = 0; k < n_ghost; ++k) {
            field(n_theta + k, j) = recv_buf_right[idx + k];
        }
        idx += n_ghost;
    }
}
