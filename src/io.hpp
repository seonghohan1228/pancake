#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

namespace IO {
    void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);
    void write_timestep(double time, int step_index, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg);
}
