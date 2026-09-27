#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

namespace IO {
    /// Create output directories and replace named solver products, preserving
    /// unrelated files. Setup failures are reported to all ranks before throwing.
    void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);
    void write_timestep(double time, int step_index, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg);
    /// Portable, rank-independent initial guesses; geometry and finite values are checked collectively.
    /// This stores no physical time history and does not assert convergence of the new solve.
    void read_steady_initial_state(const std::filesystem::path& path, Fields& fields,
                                   const Mesh& mesh, const SimulationConfig& cfg);
    void write_steady_initial_state(const std::filesystem::path& path, const Fields& fields,
                                    const Mesh& mesh, const SimulationConfig& cfg);
}
