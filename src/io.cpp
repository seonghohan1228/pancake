#include "io.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <mpi.h>

namespace fs = std::filesystem;

namespace IO {

// --- Private Helpers ---

void write_local_vts(int step_index, const Mesh& mesh, const std::vector<double>& data, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/processor" << mesh.rank << "/"
       << cfg.filename_prefix << "_" << step_index << "_" << mesh.rank << ".vts";

    std::ofstream file(ss.str());

    // VTK Structured Grid defines EXTENTS based on POINTS, not cells.
    // If we own N cells, we define N+1 points in the visualization.
    int t_start = mesh.offset_theta;
    int t_end = mesh.offset_theta + mesh.n_theta_local;
    int z_start = mesh.offset_z;
    int z_end = mesh.offset_z + mesh.n_z_local;

    // Global Extents
    int t_global_max = cfg.n_theta_global;
    int z_global_max = cfg.n_z_global;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <StructuredGrid WholeExtent=\"0 " << t_global_max << " 0 " << z_global_max << " 0 0\">\n";
    file << "    <Piece Extent=\"" << t_start << " " << t_end << " " << z_start << " " << z_end << " 0 0\">\n";

    // 1. Points
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    double d_theta = mesh.get_d_theta();
    double d_z = mesh.get_d_z();

    for (int j = z_start; j <= z_end; ++j) {
        double z = j * d_z;
        for (int i = t_start; i <= t_end; ++i) {
            double theta = i * d_theta;
            file << (mesh.R * std::cos(theta)) << " "
                 << (mesh.R * std::sin(theta)) << " "
                 << z << " ";
        }
    }
    file << "\n        </DataArray>\n      </Points>\n";

    // 2. Cell Data
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Solution\" format=\"ascii\">\n";
    for (const auto& val : data) file << val << " ";
    file << "\n        </DataArray>\n      </CellData>\n";

    file << "    </Piece>\n  </StructuredGrid>\n</VTKFile>\n";
}

void write_pvts(int step_index, int total_ranks, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/" << cfg.filename_prefix << "_" << step_index << ".pvts";
    std::ofstream file(ss.str());

    int t_global_max = cfg.n_theta_global;
    int z_global_max = cfg.n_z_global;

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << t_global_max << " 0 " << z_global_max << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData>\n"
         << "      <PDataArray type=\"Float64\" Name=\"Solution\"/>\n"
         << "    </PCellData>\n";

    // Reconstruct extents for all ranks to write valid Piece entries
    for (int r = 0; r < total_ranks; ++r) {
        int n_local = t_global_max / total_ranks;
        int remainder = t_global_max % total_ranks;
        int offset = r * n_local;

        if (r < remainder) {
            n_local++;
            offset += r;
        } else {
            offset += remainder;
        }

        file << "    <Piece Extent=\"" << offset << " " << (offset + n_local) << " 0 " << z_global_max << " 0 0\" "
             << "Source=\"processor" << r << "/" << cfg.filename_prefix << "_" << step_index << "_" << r << ".vts\"/>\n";
    }

    file << "  </PStructuredGrid>\n</VTKFile>\n";
}

void update_pvd(double time, int step_index, const SimulationConfig& cfg) {
    std::ofstream file(cfg.output_dir + "/results.pvd", std::ios::app);
    file << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" "
         << "file=\"" << cfg.filename_prefix << "_" << step_index << ".pvts\"/>\n";
}

// --- Public Interface ---

void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        if (fs::exists(cfg.output_dir)) fs::remove_all(cfg.output_dir);
        fs::create_directory(cfg.output_dir);

        // Initialize PVD header
        std::ofstream pvd(cfg.output_dir + "/results.pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n";
    }
    MPI_Barrier(comm);

    // Each rank creates its own subdirectory
    std::string my_subfolder = cfg.output_dir + "/processor" + std::to_string(rank);
    if (!fs::exists(my_subfolder)) fs::create_directory(my_subfolder);
    MPI_Barrier(comm);
}

void write_timestep(double time, int step_index, const Mesh& mesh,
                   const std::vector<double>& data, const SimulationConfig& cfg) {
    write_local_vts(step_index, mesh, data, cfg);

    if (mesh.rank == 0) {
        write_pvts(step_index, mesh.size, cfg);
        update_pvd(time, step_index, cfg);
    }
}

}
