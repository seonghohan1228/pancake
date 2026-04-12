#include <fstream>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "io.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"
#include "utils.hpp"

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    SimulationConfig cfg;
    std::string config_path = "config.txt";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" || arg == "--config") {
            if (i + 1 < argc) {
                config_path = argv[i + 1];
                i++;
            }
        }
    }

    cfg.load_from_file(config_path);
    Utils::log("Loaded configuration from: " + config_path);

    Mesh mesh(cfg);
    Communicator comm(mesh);
    IO::prepare_output_directory(cfg, MPI_COMM_WORLD);

    {
        // Initialize fields
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("theta",    mesh).fill(1.0);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("rho",      mesh).fill(cfg.rho);
        fields.add("inlet_indicator", mesh).fill(0.0);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z",     mesh, 2, GridLocation::FACE_Z);

        Utils::log("Starting Simulation...");

        // Geometry initialization
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        FilmThickness::compute_inlet_indicator(fields["inlet_indicator"], mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);

        // Initialize old_time to current state
        fields["pressure"].store_old_time();
        fields["theta"].store_old_time();

        double t = 0.0;
        double output_time = 0;
        int step = 0;

        // Time loop
        while (t <= cfg.end_t || std::abs(t - cfg.end_t) < 1e-8) {
            // Re-initialise to full film at each step (flooded guess)
            // This prevents the "starvation" trap where theta=0 persists.
            fields["theta"].fill(1.0);
            fields["pressure"].fill(cfg.p_cav);
            comm.update_ghosts(fields);

            if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS)
                Reynolds::solve_elrod(fields, sys, mesh, cfg);
            else
                Reynolds::solve(fields, sys, mesh, cfg);

            // Calculate velocities
            comm.update_ghosts(fields["pressure"]);
            comm.update_ghosts(fields["theta"]);
            Reynolds::calculate_velocities(fields, mesh, cfg);

            // Output
            if (t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0) {
                comm.update_ghosts(fields["velocity_theta"]);
                comm.update_ghosts(fields["velocity_z"]);
                IO::write_timestep(t, step, mesh, fields, cfg);
                Utils::log("Step " + std::to_string(step) + " t=" + std::to_string(t),
                           Utils::Color::GREEN, 1);
                output_time += cfg.write_interval;
            }

            // Advance time
            fields["pressure"].store_old_time();
            fields["theta"].store_old_time();
            t += cfg.dt;
            step++;
        }

        Utils::log("Simulation finished successfully.");
    }  // LinearSystem and Fields destroyed before PetscFinalize

    PetscFinalize();
    return 0;
}
