#include <petsc.h>
#include <cmath>
#include <string>

#include "bubble_config.hpp"
#include "bubble_io.hpp"
#include "bubble_linear_system.hpp"
#include "communicator.hpp"
#include "field.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"
#include "thickness_transport.hpp"
#include "utils.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    // --- Configuration ---
    BubbleConfig cfg;
    std::string config_path = "bubble_config.txt";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-c" || arg == "--config") && i + 1 < argc)
            config_path = argv[++i];
    }
    cfg.load_from_file(config_path);

    Utils::log("Loaded bubble configuration from: " + config_path);

    // --- Mesh & geometry ---
    StereoMesh mesh(cfg);
    Mask       mask(mesh, cfg.R_bubble);

    Utils::log("Grid: " + std::to_string(cfg.n_u) + " x " + std::to_string(cfg.n_v)
               + "  R = " + std::to_string(cfg.R_bubble) + " m");

    // --- Fields ---
    Fields fields;
    const int nu = mesh.n_u_local;
    const int nv = mesh.n_v_global;

    // Film thickness
    Field& h = fields.add("h", nu, nv);
    h.fill(cfg.h_initial);

    // Geometry diagnostic fields (for VTK output and verification)
    Field& cell_area_fld = fields.add("cell_area", nu, nv);
    Field& mask_fld      = fields.add("mask",      nu, nv);
    Field& g_u_fld       = fields.add("g_u",       nu, nv);
    Field& g_v_fld       = fields.add("g_v",       nu, nv);
    Field& g_mag_fld     = fields.add("g_s_mag",   nu, nv);

    for (int i = 0; i < nu; ++i) {
        for (int j = 0; j < nv; ++j) {
            cell_area_fld(i, j) = mesh.cell_area[mesh.idx(i, j)];
            mask_fld     (i, j) = mask.is_active(i, j) ? 1.0 : 0.0;
            g_u_fld      (i, j) = mesh.g_u[mesh.idx(i, j)];
            g_v_fld      (i, j) = mesh.g_v[mesh.idx(i, j)];

            const double gu = mesh.g_u[mesh.idx(i, j)];
            const double gv = mesh.g_v[mesh.idx(i, j)];
            g_mag_fld(i, j) = std::sqrt(gu*gu + gv*gv);
        }
    }

    // Apply rim BC on h and blank outside cells
    mask.enforce_rim_bc(h, cfg.h_rim);
    mask.blank_masked(h, 0.0);

    // Velocity and flux fields (Phase 3)
    Field& u_u   = fields.add("u_u",    nu, nv);
    Field& u_v   = fields.add("u_v",    nu, nv);
    Field& flux_u = fields.add("flux_u", nu, nv);
    Field& flux_v = fields.add("flux_v", nu, nv);

    // --- Communicator ---
    Communicator comm(mesh);
    comm.update_ghosts(fields);

    // --- Verification: total cell area vs analytic 2*pi*R^2 ---
    double local_area = 0.0;
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
            if (mask.is_active(i, j))
                local_area += mesh.cell_area[mesh.idx(i, j)];

    double total_area = 0.0;
    MPI_Allreduce(&local_area, &total_area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const double analytic_area = 2.0 * M_PI * cfg.R_bubble * cfg.R_bubble;
    const double area_error    = std::abs(total_area - analytic_area) / analytic_area;

    Utils::log("Hemispherical area: computed = " + std::to_string(total_area)
               + "  analytic = " + std::to_string(analytic_area)
               + "  error = " + std::to_string(area_error * 100.0) + " %");

    if (area_error > 0.01)
        Utils::log("WARNING: area error > 1%", Utils::Color::YELLOW);

    // --- Output ---
    BubbleIO::prepare_output_directory(cfg);
    BubbleIO::write_timestep(0.0, 0, mesh, mask, fields, cfg);
    Utils::log("Wrote initial VTK output.", Utils::Color::GREEN);

    // --- Linear system (Phase 3) ---
    BubbleLinearSystem sys(mesh);

    // --- Time loop ---
    double t        = 0.0;
    double out_time = cfg.write_interval;
    int    step     = 1;

    Utils::log("Starting bubble simulation...");

    while (t < cfg.end_t) {
        t += cfg.dt;

        // Phase 3: thickness transport
        h.store_old_time();
        ThicknessTransport::compute_gravity_velocity(u_u, u_v, h, mesh, mask, cfg);
        comm.update_ghosts(u_u);
        comm.update_ghosts(u_v);
        ThicknessTransport::compute_face_fluxes(flux_u, flux_v, u_u, u_v, mesh, mask);
        ThicknessTransport::solve(h, flux_u, flux_v, nullptr, sys, mesh, mask, cfg, cfg.dt);
        comm.update_ghosts(h);

        // TODO Phase 4: SIMPLE momentum solve
        // TODO Phase 5: energy, evaporation, Marangoni
        // TODO Phase 6: disjoining pressure, rupture detection

        if (t >= out_time) {
            BubbleIO::write_timestep(t, step, mesh, mask, fields, cfg);
            Utils::log("Step " + std::to_string(step) + "  t = " + std::to_string(t),
                       Utils::Color::GREEN, 1);
            out_time += cfg.write_interval;
            ++step;
        }
    }

    Utils::log("Bubble simulation finished successfully.");
    PetscFinalize();
    return 0;
}
