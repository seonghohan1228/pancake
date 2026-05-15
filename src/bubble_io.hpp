#pragma once

#include <mpi.h>
#include <string>

#include "bubble_config.hpp"
#include "field.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"

/// VTK unstructured-grid (vtu/pvtu/pvd) output for the hemispherical bubble solver.
///
/// Output layout:
///   <output_dir>/results.pvd           — time collection (rank 0)
///   <output_dir>/step_N.pvtu           — parallel descriptor (rank 0)
///   <output_dir>/rank_R/step_N.vtu     — per-rank piece
///
/// Each active cell is written as a VTK_QUAD (type 9) with 4 projected 3D vertices.
/// Blanked cells (mask.is_active == false) are omitted.
namespace BubbleIO {

    /// Create output directories and initialise the PVD file. Must be called once before
    /// the time loop. Collective (all ranks participate; rank 0 creates files).
    void prepare_output_directory(const BubbleConfig& cfg,
                                  MPI_Comm comm = MPI_COMM_WORLD);

    /// Write one timestep. Rank 0 also updates the PVD collection file.
    void write_timestep(double time, int step_index,
                        const StereoMesh& mesh, const Mask& mask,
                        Fields& fields, const BubbleConfig& cfg);

}  // namespace BubbleIO
