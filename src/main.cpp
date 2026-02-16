#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>

#include "config.hpp"
#include "field.hpp"
#include "logger.hpp"
#include "mesh.hpp"
#include "fields.hpp"
#include "io.hpp"
#include "communicator.hpp"

#define NEAR(a, b) (std::abs((a) - (b)) < 1e-8)

void print_info() {
    Logger::print("pancake", Logger::Color::CYAN); // Removed BOLD_CYAN (not in your logger enum)
    Logger::print("Thin film fluid flow solver.");
    Logger::print("");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    print_info();

    // Configuration
    SimulationConfig config;
    // TODO: Load config from file

    // I/O setup
    IO::prepare_output_directory(config);

    // Mesh initialization
    Logger::print("Initializing Mesh...");
    Mesh mesh(config);

    // Communicator initialization
    Communicator comm(mesh);

    // Solution fields
    Fields fields;
    fields.add("pressure", mesh, 2, GridLocation::CENTER);
    fields.add("film_thickness", mesh, 2, GridLocation::CENTER);
    fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA); // Staggered
    fields.add("velocity_z", mesh, 2, GridLocation::FACE_Z);

    fields["pressure"].fill(101325.0);
    fields["film_thickness"].fill(1.0);
    fields["velocity_theta"].fill(0.0);
    fields["velocity_z"].fill(0.0);

    // Decomposition fields
    fields.add("processor", mesh);
    fields["processor"].fill(static_cast<double>(rank));

    // Apply initial ghost cell exchange
    comm.update_ghosts(fields);

    Logger::print("Starting time loop.");

    // Time stepping loop
    double t = 0;
    double next_output_time = 0;
    int step_index = 0;

    while (t <= config.end_t || NEAR(t, config.end_t)) {

        // --- SOLVER STEP ---
        // TODO: Implement solver step
        // Update ghosts after every change
        // comm.update_ghosts(fields);

        // --- I/O STEP ---
        if (t >= next_output_time || NEAR(t, next_output_time)) {
            std::string msg = "Writing output t = " + std::to_string(t);
            Logger::print(msg, Logger::Color::GREEN, 1);

            // Ensure ghosts are valid before interpolating for IO (optional but good practice)
            comm.update_ghosts(fields);

            IO::write_timestep(t, step_index, mesh, fields, config);

            next_output_time += config.write_interval;
            step_index++;
        }

        t += config.dt;
    }

    // Finalization
    if (rank == 0) {
        // Close the PVD file tags
        std::ofstream pvd(config.output_dir + "/results.pvd", std::ios::app);
        pvd << "  </Collection>\n</VTKFile>\n";
        Logger::print("\nSimulation finished successfully.\n");
    }

    MPI_Finalize();
    return 0;
}
