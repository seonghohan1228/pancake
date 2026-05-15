/// Phase B2 — Stereographic mesh geometry tests.
///
/// Checks: area sum ≈ 2πR², unit normals, gravity at pole/equator, positive geometry.
#include <petsc.h>
#include <cmath>
#include <iostream>

#include "bubble_config.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"

static int g_pass = 0, g_fail = 0;

static void check(const std::string& name, bool cond, int rank) {
    if (rank == 0)
        std::cout << (cond ? "[PASS] " : "[FAIL] ") << name << "\n";
    if (cond) ++g_pass; else ++g_fail;
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    BubbleConfig cfg;
    cfg.n_u = 64; cfg.n_v = 64;
    cfg.R_bubble = 0.01; cfg.L_box = 0.012;
    cfg.g_x = 0.0; cfg.g_y = 0.0; cfg.g_z = -9.81;

    StereoMesh mesh(cfg);
    Mask mask(mesh, cfg.R_bubble);

    // --- North pole projection ---
    {
        Vec3 p = mesh.project(0.0, 0.0);
        check("North pole: x=0", std::abs(p.x) < 1e-15, rank);
        check("North pole: y=0", std::abs(p.y) < 1e-15, rank);
        check("North pole: z=R", std::abs(p.z - cfg.R_bubble) < 1e-14, rank);
    }

    // --- Equator projection (r_2D = R maps to z=0) ---
    {
        Vec3 p = mesh.project(cfg.R_bubble, 0.0);
        check("Equator: z≈0",   std::abs(p.z) < 1e-14, rank);
        check("Equator: |p|=R", std::abs(p.norm() - cfg.R_bubble) < 1e-13, rank);
    }

    // --- Total hemisphere area ≈ 2πR² ---
    {
        double local_area = 0.0;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j)
                if (mask.is_active(i, j))
                    local_area += mesh.cell_area[mesh.idx(i, j)];
        double total_area;
        MPI_Allreduce(&local_area, &total_area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double analytic = 2.0 * M_PI * cfg.R_bubble * cfg.R_bubble;
        check("Area sum within 1% of 2πR²",
              std::abs(total_area - analytic) / analytic < 0.01, rank);
    }

    // --- Normals are unit vectors ---
    {
        bool ok = true;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j) {
                if (!mask.is_active(i, j)) continue;
                double nx = mesh.cell_nx[mesh.idx(i, j)];
                double ny = mesh.cell_ny[mesh.idx(i, j)];
                double nz = mesh.cell_nz[mesh.idx(i, j)];
                if (std::abs(nx*nx + ny*ny + nz*nz - 1.0) > 1e-12) { ok = false; break; }
            }
        check("All normals are unit vectors", ok, rank);
    }

    // --- Gravity at north pole: |g_s| = 0 ---
    {
        Vec3 n = mesh.project(0.0, 0.0).normalized();
        Vec3 g{cfg.g_x, cfg.g_y, cfg.g_z};
        Vec3 gs = g - n * g.dot(n);
        check("|g_s| = 0 at north pole", gs.norm() < 1e-12, rank);
    }

    // --- Gravity at equator: |g_s| = |g| ---
    {
        Vec3 n = mesh.project(cfg.R_bubble, 0.0).normalized();
        Vec3 g{cfg.g_x, cfg.g_y, cfg.g_z};
        Vec3 gs = g - n * g.dot(n);
        check("|g_s| = |g| at equator",
              std::abs(gs.norm() - g.norm()) < 1e-12, rank);
    }

    // --- All geometry values positive for active cells ---
    {
        bool ok = true;
        for (int i = 0; i < mesh.n_u_local && ok; ++i)
            for (int j = 0; j < mesh.n_v_global && ok; ++j) {
                if (!mask.is_active(i, j)) continue;
                int idx = mesh.idx(i, j);
                if (mesh.cell_area[idx] <= 0 ||
                    mesh.face_len_e[idx] <= 0 || mesh.face_len_n[idx] <= 0 ||
                    mesh.d_u[idx] <= 0       || mesh.d_v[idx] <= 0)
                    ok = false;
            }
        check("All areas, face lengths, distances > 0", ok, rank);
    }

    // --- east_face_len helper: local and analytic agree ---
    {
        // For a local cell, east_face_len(gu, j) should match face_len_e[idx(i,j)]
        int i = 0, j = mesh.n_v_global / 2;
        int gu = mesh.global_u(i);
        double stored  = mesh.face_len_e[mesh.idx(i, j)];
        double computed = mesh.east_face_len(gu, j);
        check("east_face_len helper matches precomputed for local cell",
              std::abs(stored - computed) < 1e-14, rank);
    }

    if (rank == 0)
        std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n";

    PetscFinalize();
    return (g_fail > 0) ? 1 : 0;
}
