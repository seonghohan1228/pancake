/// Phase B3 — SurfaceFVM operator unit tests.
///
/// Checks:
///   ddt:       a_P = A/dt, source = A/dt * phi_old for each active cell.
///   laplacian: uniform field → net LHS = 0 for interior active cells.
///   divergence:constant flux field → source contribution matches analytic sum.
///   add_source:source *= A_cell.
#include <petsc.h>
#include <cmath>
#include <iostream>

#include "bubble_config.hpp"
#include "bubble_linear_system.hpp"
#include "communicator.hpp"
#include "field.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"
#include "surface_operators.hpp"

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
    cfg.g_z = -9.81;

    StereoMesh   mesh(cfg);
    Mask         mask(mesh, cfg.R_bubble);
    Communicator comm(mesh);

    // Scope sys so its PETSc objects are destroyed before PetscFinalize.
    {
    BubbleLinearSystem sys(mesh);

    const int nu = mesh.n_u_local, nv = mesh.n_v_global;
    const double dt = 1.0e-3;

    // --- ddt: coefficients match A_cell / dt, source = A_cell/dt * phi_old ---
    {
        Field phi("phi", nu, nv);
        phi.fill(3.14);
        phi.store_old_time();

        sys.reset();
        SurfaceFVM::ddt(sys, phi, dt, mesh, mask);

        bool ap_ok = true, src_ok = true;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                if (!mask.is_active(i, j)) continue;
                const double A   = mesh.cell_area[mesh.idx(i, j)];
                const double expected_ap  = A / dt;
                const double expected_src = (A / dt) * phi.old(i, j);
                if (std::abs(sys.a_p(i,j) - expected_ap) / expected_ap > 1e-12) ap_ok = false;
                if (std::abs(sys.source(i,j) - expected_src) / std::abs(expected_src) > 1e-12) src_ok = false;
            }
        }
        check("ddt: a_P = A_cell / dt", ap_ok, rank);
        check("ddt: source = A_cell / dt * phi_old", src_ok, rank);
    }

    // --- ddt_weighted: coefficients match w(i,j)*A_cell/dt ---
    {
        Field phi("phi", nu, nv);
        phi.fill(1.0);
        phi.store_old_time();

        Field weight("w", nu, nv);
        // Vary weight to make the test non-trivial
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j < nv; ++j)
                weight(i,j) = 2.0 + 0.1 * (i + j);

        sys.reset();
        SurfaceFVM::ddt_weighted(sys, phi, weight, dt, mesh, mask);

        bool ok = true;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                if (!mask.is_active(i,j)) continue;
                const double expected = weight(i,j) * mesh.cell_area[mesh.idx(i,j)] / dt;
                if (std::abs(sys.a_p(i,j) - expected) / expected > 1e-12) { ok = false; break; }
            }
        }
        check("ddt_weighted: a_P = w * A_cell / dt", ok, rank);
    }

    // --- laplacian: uniform field → net LHS sum ≈ 0 for interior active cells ---
    {
        // For a uniform gamma and phi=const, the laplacian should give zero net flux.
        // Equivalently: for each interior cell, a_P * 1 - a_E * 1 - a_W * 1 - a_N * 1 - a_S * 1 = 0.
        Field gamma("g", nu, nv);
        gamma.fill(1.0);
        comm.update_ghosts(gamma);

        sys.reset();
        SurfaceFVM::laplacian(sys, gamma, mesh, mask);

        double max_residual = 0.0;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                if (!mask.is_active(i,j) || mask.is_rim(i,j)) continue;
                int gu = mesh.global_u(i);
                // Interior non-rim active cell: all 4 neighbours are active
                bool all_nbrs_active =
                    (gu+1 < mesh.n_u_global) && mask.is_active(i+1, j) && !mask.is_rim(i+1, j) &&
                    (gu-1 >= 0)              && mask.is_active(i-1, j) && !mask.is_rim(i-1, j) &&
                    (j+1 < nv)               && mask.is_active(i, j+1) && !mask.is_rim(i, j+1) &&
                    (j-1 >= 0)               && mask.is_active(i, j-1) && !mask.is_rim(i, j-1);
                if (!all_nbrs_active) continue;
                // For uniform gamma and phi: a_P = a_E + a_W + a_N + a_S
                double residual = sys.a_p(i,j) - sys.a_e(i,j) - sys.a_w(i,j)
                                               - sys.a_n(i,j) - sys.a_s(i,j);
                max_residual = std::max(max_residual, std::abs(residual));
            }
        }
        double global_max;
        MPI_Allreduce(&max_residual, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        check("laplacian: a_P = Σ(neighbours) for interior cells", global_max < 1e-12, rank);
    }

    // --- add_source: contribution = field * A_cell ---
    {
        Field src_field("src", nu, nv);
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j < nv; ++j)
                src_field(i,j) = 2.5;

        sys.reset();
        SurfaceFVM::add_source(sys, src_field, mesh, mask);

        bool ok = true;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                if (!mask.is_active(i,j)) continue;
                const double expected = 2.5 * mesh.cell_area[mesh.idx(i,j)];
                if (std::abs(sys.source(i,j) - expected) / expected > 1e-12)
                    { ok = false; break; }
            }
        }
        check("add_source: source += field * A_cell", ok, rank);
    }

    // --- divergence: zero flux field → all coefficients zero ---
    {
        Field fu("fu", nu, nv);  fu.fill(0.0);
        Field fv("fv", nu, nv);  fv.fill(0.0);
        Field phi("ph", nu, nv); phi.fill(1.0); phi.store_old_time();
        comm.update_ghosts(fu);  comm.update_ghosts(fv);  comm.update_ghosts(phi);

        sys.reset();
        SurfaceFVM::divergence(sys, fu, fv, phi, ConvectionScheme::UPWIND, mesh, mask);

        bool ok = true;
        for (int i = 0; i < nu; ++i)
            for (int j = 0; j < nv; ++j)
                if (sys.a_p(i,j) != 0 || sys.a_e(i,j) != 0 || sys.a_w(i,j) != 0 ||
                    sys.a_n(i,j) != 0 || sys.a_s(i,j) != 0)
                    { ok = false; break; }
        check("divergence: zero flux → all coefficients zero", ok, rank);
    }

    } // destroy BubbleLinearSystem before PetscFinalize

    if (rank == 0)
        std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n";

    PetscFinalize();
    return (g_fail > 0) ? 1 : 0;
}
