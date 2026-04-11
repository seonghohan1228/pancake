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
    Mesh mesh(cfg);
    Communicator comm(mesh);
    IO::prepare_output_directory(cfg, MPI_COMM_WORLD);

    {
        // Initialize fields
        Fields fields;
        fields.add("pressure", mesh).fill(101325.0);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho);

        Utils::log("Starting Simulation...");

        // Geometry initialization
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);

        // Initialize old_time to current state
        fields["pressure"].store_old_time();

        double t = 0.0;
        double output_time = 0;
        int step = 0;

        // Time loop
        while (t <= cfg.end_t || std::abs(t - cfg.end_t) < 1e-8) {
            // TODO: update h for dynamic bearing (Phase 4+)
            comm.update_ghosts(fields);

            Reynolds::solve(fields, sys, mesh, cfg);

            // Output
            if (t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0) {
                IO::write_timestep(t, step, mesh, fields, cfg);
                Utils::log("Step " + std::to_string(step) + " t=" + std::to_string(t),
                           Utils::Color::GREEN, 1);
                output_time += cfg.write_interval;
            }

            // Advance time
            fields["pressure"].store_old_time();
            t += cfg.dt;
            step++;
        }

        // Close the PVD file
        if (mesh.rank == 0) {
            std::ofstream pvd(cfg.output_dir + "/results.pvd", std::ios::app);
            pvd << "  </Collection>\n"
                << "</VTKFile>\n";
        }

        Utils::log("Simulation finished successfully.");
    }  // LinearSystem and Fields destroyed before PetscFinalize

    PetscFinalize();
    return 0;
}
