#include "communicator.hpp"

Communicator::Communicator(const StereoMesh& mesh_ref) : mesh(mesh_ref) {}

void Communicator::update_ghosts(Fields& fields) {
    for (auto& [name, field_ptr] : fields)
        update_ghosts(*field_ptr);
}

void Communicator::update_ghosts(Field& field) {
    const int n_ghost  = field.n_ghost;
    const int n_u      = field.n_theta_phys;
    const int n_v      = field.n_z_phys;
    const int buf_size = n_v * n_ghost;

    // Resize send/recv buffers if needed
    if ((int)send_buf_left.size() < buf_size) {
        send_buf_left .resize(buf_size);
        recv_buf_left .resize(buf_size);
        send_buf_right.resize(buf_size);
        recv_buf_right.resize(buf_size);
    }

    // --- u-direction exchange ---

    // Pack: left physical strip (i=0..n_ghost-1) and right strip (i=n_u-n_ghost..n_u-1)
    int idx = 0;
    for (int j = 0; j < n_v; ++j) {
        for (int k = 0; k < n_ghost; ++k) send_buf_left [idx + k] = field(k,           j);
        for (int k = 0; k < n_ghost; ++k) send_buf_right[idx + k] = field(n_u - n_ghost + k, j);
        idx += n_ghost;
    }

    const int left_rank  = mesh.rank - 1;
    const int right_rank = mesh.rank + 1;

    MPI_Status status;

    // Exchange with right neighbour: send our right strip, receive their left strip.
    if (right_rank < mesh.size) {
        MPI_Sendrecv(send_buf_right.data(), buf_size, MPI_DOUBLE, right_rank, 0,
                     recv_buf_right.data(), buf_size, MPI_DOUBLE, right_rank, 0,
                     MPI_COMM_WORLD, &status);
    } else {
        // Rightmost rank: copy last physical cell into right ghost (zero-gradient)
        int k = 0;
        for (int j = 0; j < n_v; ++j)
            for (int g = 0; g < n_ghost; ++g)
                recv_buf_right[k++] = field(n_u - 1, j);
    }

    // Exchange with left neighbour: send our left strip, receive their right strip.
    if (left_rank >= 0) {
        MPI_Sendrecv(send_buf_left.data(),  buf_size, MPI_DOUBLE, left_rank,  0,
                     recv_buf_left.data(),  buf_size, MPI_DOUBLE, left_rank,  0,
                     MPI_COMM_WORLD, &status);
    } else {
        // Leftmost rank: copy first physical cell into left ghost (zero-gradient)
        int k = 0;
        for (int j = 0; j < n_v; ++j)
            for (int g = 0; g < n_ghost; ++g)
                recv_buf_left[k++] = field(0, j);
    }

    // Unpack received data into ghost cells
    idx = 0;
    for (int j = 0; j < n_v; ++j) {
        for (int k = 0; k < n_ghost; ++k) {
            field(-n_ghost + k, j) = recv_buf_left [idx + k];
            field(n_u + k,      j) = recv_buf_right[idx + k];
        }
        idx += n_ghost;
    }

    // --- v-direction: zero-gradient fill (not decomposed) ---
    fill_v_ghosts(field);
}

void Communicator::fill_v_ghosts(Field& field) {
    const int n_ghost = field.n_ghost;
    const int n_u     = field.n_theta_phys;
    const int n_v     = field.n_z_phys;

    for (int i = -n_ghost; i < n_u + n_ghost; ++i) {
        for (int g = 0; g < n_ghost; ++g) {
            field(i, -g - 1) = field(i, 0);       // South ghosts
            field(i, n_v + g) = field(i, n_v - 1); // North ghosts
        }
    }
}
