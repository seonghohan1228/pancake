#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

namespace IO {
    /// A directory holding this file is a run archive root. prepare_output_directory
    /// refuses to use it as output_dir; runs must write to a campaign folder inside it.
    inline constexpr const char* archive_root_marker = ".pancake_archive_root";

    /// Create output directories and replace named solver products, preserving
    /// unrelated files. Only files directly in output_dir, processorN/, flat/ and
    /// flat/processorN/ are cleared; other subfolders are never entered.
    /// Refuses an output_dir that contains archive_root_marker. Setup failures
    /// are reported to all ranks before throwing.
    void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm = MPI_COMM_WORLD);
    void write_timestep(double time, int step_index, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg);
    /// Portable, rank-independent initial guesses; geometry and finite values are checked collectively.
    /// This stores no physical time history and does not assert convergence of the new solve.
    void read_steady_initial_state(const std::filesystem::path& path, Fields& fields,
                                   const Mesh& mesh, const SimulationConfig& cfg);
    void write_steady_initial_state(const std::filesystem::path& path, const Fields& fields,
                                    const Mesh& mesh, const SimulationConfig& cfg);
}
