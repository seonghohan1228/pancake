/// Phase B3 — ThicknessTransport unit tests.
///
/// Checks:
///   compute_gravity_velocity: zero gravity → zero velocity; non-zero g_z produces
///                             nonzero velocity only where g_s != 0 (not at pole).
///   compute_face_fluxes:      uniform velocity → face flux = u * face_len.
///   solve (zero velocity):    h unchanged after one step when there are no fluxes
///                             and no evaporation (transient + zero divergence).
///   solve (h_min clamp):      negative initial h is clamped to h_min after solve.
#include <petsc.h>
#include <cmath>
#include <iostream>

#include "bubble_config.hpp"
#include "bubble_linear_system.hpp"
#include "communicator.hpp"
#include "field.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"
#include "thickness_transport.hpp"

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
    cfg.g_x = 0.0; cfg.g_y = 0.0; cfg.g_z = -9.81;
    cfg.rho_l = 1000.0;
    cfg.mu    = 1.0e-3;
    cfg.h_rim = 1.0e-6;
    cfg.h_min = 1.0e-9;
    cfg.h_initial = 1.0e-5;
    cfg.advection_scheme = ConvectionScheme::UPWIND;

    StereoMesh   mesh(cfg);
    Mask         mask(mesh, cfg.R_bubble);
    Communicator comm(mesh);

    // Scope sys so its PETSc objects are destroyed before PetscFinalize.
    {
    BubbleLinearSystem sys(mesh);

    const int nu = mesh.n_u_local;
    const int nv = mesh.n_v_global;
    const double dt = 1.0e-3;

    // --- compute_gravity_velocity: zero gravity → zero velocity ---
    {
        BubbleConfig cfg0 = cfg;
        cfg0.g_x = 0.0; cfg0.g_y = 0.0; cfg0.g_z = 0.0;
        StereoMesh mesh0(cfg0);
        Mask       mask0(mesh0, cfg0.R_bubble);

        Field h0("h", nu, nv);   h0.fill(cfg0.h_initial);
        Field uu("uu", nu, nv),  uv("uv", nu, nv);

        ThicknessTransport::compute_gravity_velocity(uu, uv, h0, mesh0, mask0, cfg0);

        bool ok = true;
        for (int i = 0; i < nu && ok; ++i)
            for (int j = 0; j < nv && ok; ++j)
                if (uu(i,j) != 0.0 || uv(i,j) != 0.0)
                    ok = false;
        check("compute_gravity_velocity: zero g → zero u_s", ok, rank);
    }

    // --- compute_gravity_velocity: nonzero g_z → nonzero surface velocity except at pole ---
    {
        Field h("h", nu, nv);  h.fill(cfg.h_initial);
        Field uu("uu", nu, nv), uv("uv", nu, nv);
        ThicknessTransport::compute_gravity_velocity(uu, uv, h, mesh, mask, cfg);

        // Pole cell (r≈0): g_s = 0 → velocity should be zero.
        // Off-pole cells: |g_s| > 0 → velocity should be nonzero.
        bool pole_ok = true, offpole_ok = true;
        for (int i = 0; i < nu; ++i) {
            for (int j = 0; j < nv; ++j) {
                if (!mask.is_active(i,j) || mask.is_rim(i,j)) continue;
                const double r2 = mesh.cell_u(i)*mesh.cell_u(i)
                                + mesh.cell_v(j)*mesh.cell_v(j);
                const bool near_pole = r2 < (0.01 * cfg.R_bubble * cfg.R_bubble);
                const double speed2  = uu(i,j)*uu(i,j) + uv(i,j)*uv(i,j);
                if (near_pole && speed2 > 1e-30) pole_ok    = false;
                if (!near_pole && speed2 == 0.0) offpole_ok = false;
            }
        }
        int pole_int = pole_ok ? 1 : 0, offpole_int = offpole_ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &pole_int,    1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &offpole_int, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        check("compute_gravity_velocity: pole velocity is zero", pole_int == 1, rank);
        check("compute_gravity_velocity: off-pole velocity nonzero", offpole_int == 1, rank);
    }

    // --- compute_face_fluxes: uniform velocity → flux = u * face_len_e ---
    {
        const double u_const = 0.5;
        Field uu("uu", nu, nv, 2), uv("uv", nu, nv, 2);
        uu.fill(u_const);  uv.fill(0.0);
        comm.update_ghosts(uu);  comm.update_ghosts(uv);

        Field fu("fu", nu, nv), fv("fv", nu, nv);
        ThicknessTransport::compute_face_fluxes(fu, fv, uu, uv, mesh, mask);

        bool ok = true;
        for (int i = 0; i < nu && ok; ++i) {
            for (int j = 0; j < nv && ok; ++j) {
                const double expected = u_const * mesh.face_len_e[mesh.idx(i,j)];
                if (std::abs(fu(i,j) - expected) > 1e-12 * std::abs(expected) + 1e-30)
                    ok = false;
            }
        }
        check("compute_face_fluxes: uniform u → flux = u*face_len_e", ok, rank);
    }

    // --- solve: zero velocity + zero evap → h unchanged (transient only) ---
    {
        const double h_val = cfg.h_initial;
        Field h("h", nu, nv);   h.fill(h_val);
        mask.enforce_rim_bc(h, cfg.h_rim);
        mask.blank_masked(h, 0.0);
        h.store_old_time();

        // Zero face fluxes
        Field fu("fu", nu, nv), fv("fv", nu, nv);
        fu.fill(0.0);  fv.fill(0.0);
        comm.update_ghosts(fu);  comm.update_ghosts(fv);

        ThicknessTransport::solve(h, fu, fv, nullptr, sys, mesh, mask, cfg, dt);

        // Interior active non-rim cells: h should remain h_val (within solver tolerance)
        bool ok = true;
        for (int i = 0; i < nu && ok; ++i) {
            for (int j = 0; j < nv && ok; ++j) {
                if (!mask.is_active(i,j) || mask.is_rim(i,j)) continue;
                if (std::abs(h(i,j) - h_val) / h_val > 1e-6)
                    ok = false;
            }
        }
        int ok_int = ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &ok_int, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        check("solve: zero flux → h unchanged", ok_int == 1, rank);
    }

    // --- solve: h below h_min is clamped upward ---
    {
        const double h_neg = -1.0e-5;   // deliberately below h_min
        Field h("h", nu, nv);  h.fill(h_neg);
        mask.enforce_rim_bc(h, cfg.h_rim);
        mask.blank_masked(h, 0.0);
        h.store_old_time();

        Field fu("fu", nu, nv), fv("fv", nu, nv);
        fu.fill(0.0);  fv.fill(0.0);

        ThicknessTransport::solve(h, fu, fv, nullptr, sys, mesh, mask, cfg, dt);

        bool ok = true;
        for (int i = 0; i < nu && ok; ++i)
            for (int j = 0; j < nv && ok; ++j)
                if (mask.is_active(i,j) && h(i,j) < cfg.h_min - 1e-30)
                    ok = false;
        int ok_int = ok ? 1 : 0;
        MPI_Allreduce(MPI_IN_PLACE, &ok_int, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        check("solve: h clamped to h_min after solve", ok_int == 1, rank);
    }

    } // destroy BubbleLinearSystem before PetscFinalize

    if (rank == 0)
        std::cout << "\n" << g_pass << " passed, " << g_fail << " failed.\n";

    PetscFinalize();
    return (g_fail > 0) ? 1 : 0;
}
