#include "stereo_mesh.hpp"

#include <algorithm>
#include <stdexcept>

StereoMesh::StereoMesh(const BubbleConfig& cfg, MPI_Comm comm) {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    R          = cfg.R_bubble;
    L_box      = cfg.L_box;
    n_u_global = cfg.n_u;
    n_v_global = cfg.n_v;
    n_v_local  = n_v_global;
    du         = 2.0 * L_box / n_u_global;
    dv         = 2.0 * L_box / n_v_global;

    // 1D decomposition in u
    int base  = n_u_global / size;
    int rem   = n_u_global % size;
    n_u_local = base + (rank < rem ? 1 : 0);
    offset_u  = rank * base + std::min(rank, rem);

    if (n_u_local < 1)
        throw std::runtime_error("StereoMesh: too many MPI ranks for grid size");

    compute_geometry(cfg.g_x, cfg.g_y, cfg.g_z);
}

double StereoMesh::east_face_len(int gu, int j) const {
    int li = gu - offset_u;
    if (li >= 0 && li < n_u_local) return face_len_e[idx(li, j)];
    double u_face = -L_box + (gu + 1) * du;
    double v0     = vertex_v(j);
    return (project(u_face, v0 + dv) - project(u_face, v0)).norm();
}

double StereoMesh::d_u_east(int gu, int j) const {
    int li = gu - offset_u;
    if (li >= 0 && li < n_u_local) return d_u[idx(li, j)];
    double uc = -L_box + (gu + 0.5) * du;
    double vc = cell_v(j);
    return (project(uc + du, vc) - project(uc, vc)).norm();
}

Vec3 StereoMesh::project(double u, double v) const {
    double r2 = u*u + v*v;
    double D  = r2 + R*R;
    return {2.0*R*R*u/D, 2.0*R*R*v/D, R*(R*R - r2)/D};
}

void StereoMesh::compute_geometry(double gx, double gy, double gz) {
    const int n = n_u_local * n_v_global;
    cell_area.resize(n);
    face_len_e.resize(n);
    face_len_n.resize(n);
    d_u.resize(n);
    d_v.resize(n);
    cell_nx.resize(n);
    cell_ny.resize(n);
    cell_nz.resize(n);
    g_u.resize(n);
    g_v.resize(n);

    const double eps = du * 1.0e-5;  // Small step for tangent vector estimation

    const Vec3 gvec{gx, gy, gz};

    for (int i = 0; i < n_u_local; ++i) {
        for (int j = 0; j < n_v_global; ++j) {
            const double uc = cell_u(i);
            const double vc = cell_v(j);
            const double u0 = vertex_u(i);
            const double u1 = u0 + du;
            const double v0 = vertex_v(j);
            const double v1 = v0 + dv;

            // 4 corners of the cell in 3D (SW, SE, NE, NW)
            const Vec3 p00 = project(u0, v0);
            const Vec3 p10 = project(u1, v0);
            const Vec3 p11 = project(u1, v1);
            const Vec3 p01 = project(u0, v1);

            // Cell area = 0.5 * |diag1 x diag2|
            cell_area[idx(i,j)] = 0.5 * (p11 - p00).cross(p01 - p10).norm();

            // Face lengths (chord in 3D)
            face_len_e[idx(i,j)] = (p11 - p10).norm();  // East face (v-aligned)
            face_len_n[idx(i,j)] = (p11 - p01).norm();  // North face (u-aligned)

            // Centre-to-centre 3D distances
            const Vec3 c_cur   = project(uc, vc);
            const Vec3 c_east  = project(uc + du, vc);
            const Vec3 c_north = project(uc, vc + dv);
            d_u[idx(i,j)] = (c_east  - c_cur).norm();
            d_v[idx(i,j)] = (c_north - c_cur).norm();

            // Surface normal: outward unit normal on sphere = position / radius
            const Vec3 n_hat = c_cur.normalized();
            cell_nx[idx(i,j)] = n_hat.x;
            cell_ny[idx(i,j)] = n_hat.y;
            cell_nz[idx(i,j)] = n_hat.z;

            // Surface gravity: project out normal component
            const double g_dot_n = gvec.dot(n_hat);
            const Vec3 gs = gvec - n_hat * g_dot_n;

            // Local tangent vectors via central difference of stereographic map
            const Vec3 dxdu = (project(uc + eps, vc) - project(uc - eps, vc)) * (0.5 / eps);
            const Vec3 dxdv = (project(uc, vc + eps) - project(uc, vc - eps)) * (0.5 / eps);
            const Vec3 eu   = dxdu.normalized();
            const Vec3 ev   = dxdv.normalized();

            g_u[idx(i,j)] = gs.dot(eu);
            g_v[idx(i,j)] = gs.dot(ev);
        }
    }
}
