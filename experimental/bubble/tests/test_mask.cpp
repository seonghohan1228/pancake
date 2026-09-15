/// Phase B2 — Circular mask tests.
///
/// Checks: active cell count, rim detection, enforce_rim_bc, blank_masked.
#include <petsc.h>
#include <cmath>
#include <iostream>

#include "bubble_config.hpp"
#include "field.hpp"
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
    cfg.n_u = 32; cfg.n_v = 32;
    cfg.R_bubble = 0.01; cfg.L_box = 0.012;

    StereoMesh mesh(cfg);
    Mask mask(mesh, cfg.R_bubble);

    // --- Active cells: all with r_2D <= R ---
    {
        int local_active = 0, local_inactive = 0;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j) {
                double uc = mesh.cell_u(i);
                double vc = mesh.cell_v(j);
                double r2 = uc*uc + vc*vc;
                bool expected = (r2 <= cfg.R_bubble * cfg.R_bubble);
                bool actual   = mask.is_active(i, j);
                if (actual == expected) ++local_active; else ++local_inactive;
            }
        int total_wrong;
        MPI_Allreduce(&local_inactive, &total_wrong, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        check("Active flags match r_2D <= R condition", total_wrong == 0, rank);
    }

    // --- Cells at origin are active ---
    {
        // Only rank that owns the origin cell can check
        int global_u_origin = (mesh.n_u_global) / 2;  // roughly at u=0
        int j_origin = mesh.n_v_global / 2;
        bool ok = true;
        for (int i = 0; i < mesh.n_u_local; ++i) {
            if (mesh.global_u(i) == global_u_origin) {
                // Cell near center should be active (u≈0, v≈0 << R)
                ok = mask.is_active(i, j_origin);
                break;
            }
        }
        int ok_int = ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &ok_int, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        check("Center cell is active", ok_int == 1, rank);
    }

    // --- Rim cells are active and have blanked neighbour ---
    {
        int bad = 0;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j) {
                if (!mask.is_rim(i, j)) continue;
                if (!mask.is_active(i, j)) { ++bad; continue; }
                // Must have at least one inactive neighbour
                int gu = mesh.global_u(i);
                bool has_inactive =
                    (gu - 1 < 0   || !mask.is_active(i-1, j)) ||
                    (gu + 1 >= mesh.n_u_global || !mask.is_active(i+1, j)) ||
                    (j - 1 < 0   || !mask.is_active(i, j-1)) ||
                    (j + 1 >= mesh.n_v_global  || !mask.is_active(i, j+1));
                // Note: for cross-rank neighbors we can't check directly,
                // so only count cells we can verify
                if (gu > 0 && gu < mesh.n_u_global - 1 && j > 0 && j < mesh.n_v_global - 1)
                    if (!has_inactive) ++bad;
            }
        int total_bad;
        MPI_Allreduce(&bad, &total_bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        check("Rim cells are active and adjacent to blank", total_bad == 0, rank);
    }

    // --- enforce_rim_bc sets rim cells to given value ---
    {
        Field f("test", mesh.n_u_local, mesh.n_v_global);
        f.fill(0.0);
        mask.enforce_rim_bc(f, 42.0);

        bool ok = true;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j) {
                if (mask.is_rim(i, j) && std::abs(f(i,j) - 42.0) > 1e-14) { ok = false; break; }
                if (!mask.is_rim(i,j) && !mask.is_active(i,j) && f(i,j) != 0.0) { ok = false; break; }
            }
        check("enforce_rim_bc sets rim cells correctly", ok, rank);
    }

    // --- blank_masked zeroes blanked cells, leaves active unchanged ---
    {
        Field f("test2", mesh.n_u_local, mesh.n_v_global);
        f.fill(1.0);
        mask.blank_masked(f, -1.0);

        bool ok = true;
        for (int i = 0; i < mesh.n_u_local; ++i)
            for (int j = 0; j < mesh.n_v_global; ++j) {
                if (!mask.is_active(i,j) && f(i,j) != -1.0) { ok = false; break; }
                if ( mask.is_active(i,j) && f(i,j) != 1.0)  { ok = false; break; }
            }
        check("blank_masked fills only masked cells", ok, rank);
    }

    if (rank == 0)
        std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n";

    PetscFinalize();
    return (g_fail > 0) ? 1 : 0;
}
