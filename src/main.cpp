#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>
#include <vector>

#include "config.hpp"
#include "mesh.hpp"
#include "io.hpp"

#define NEAR(a, b) (std::abs((a) - (b)) < 1e-8)

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) std::cout << "[pancake] Simulation started.\n";

    // 1. Configuration
    SimulationConfig config;
    // TODO: Load config from file here if needed

    // 2. Mesh Initialization
    if(rank == 0) std::cout << "[pancake] Initializing Mesh...\n";
    Mesh mesh(config);

    // 3. I/O Setup (Clean start)
    IO::prepare_output_directory(config);

    // 4. State Initialization
    // TODO: Create a proper 'State' or 'Fields' class if variables increase (u, v, p, T...)
    std::vector<double> solution_field(mesh.local_cell_count());

    // TODO: Apply actual initial conditions here
    for(auto& val : solution_field) {
        val = 0.0;
    }

    if(rank == 0) std::cout << "[pancake] Initialization complete. Starting time loop.\n";

    // 5. Time Stepping Loop
    double t = 0;
    double next_output_time = 0;
    int step_index = 0;

    while (t <= config.end_t || NEAR(t, config.end_t)) {

        // --- SOLVER STEP ---
        // TODO: Implement your Finite Volume / Finite Difference update here.
        // Example:
        // calculate_fluxes(...);
        // update_variables(...);
        // exchange_ghost_cells(...);

        // Dummy update for visualization testing
        for(auto& val : solution_field) {
            val += 0.01 * (mesh.rank + 1);
        }

        // --- I/O STEP ---
        if (t >= next_output_time || NEAR(t, next_output_time)) {
            if(rank == 0) std::cout << "  -> Writing output t=" << t << std::endl;

            IO::write_timestep(t, step_index, mesh, solution_field, config);

            next_output_time += config.write_interval;
            step_index++;
        }

        t += config.dt;
    }

    // 6. Finalization
    if (rank == 0) {
        // Close the PVD file tags
        std::ofstream pvd(config.output_dir + "/results.pvd", std::ios::app);
        pvd << "  </Collection>\n</VTKFile>\n";
        std::cout << "[pancake] Simulation finished successfully.\n";
    }

    MPI_Finalize();
    return 0;
}
