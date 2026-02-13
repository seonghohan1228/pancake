#pragma once
#include <vector>
#include "config.hpp"
#include "mesh.hpp"

namespace IO {
    // Prepares output directory structure. Rank 0 wipes old data.
    void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);

    // Main output function: writes local .vts and updates global .pvts/.pvd
    void write_timestep(double time, int step_index, const Mesh& mesh,
                       const std::vector<double>& data,
                       const SimulationConfig& cfg);
}
