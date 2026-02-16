#include "io.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <mpi.h>
#include <iostream>

namespace fs = std::filesystem;

namespace IO {

// --- Helper Functions ---

// Interpolate face-based data to cell centers for visualization
std::vector<double> get_centered_data(const Field& f, const Mesh& mesh) {
    std::vector<double> result;
    result.reserve(mesh.n_theta_local * mesh.n_z_local);

    if (f.location == GridLocation::CENTER) {
        // Copy directly
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                result.push_back(f(i, j));
            }
        }
    }
    else if (f.location == GridLocation::FACE_THETA) {
        // Interpolate u(i) and u(i+1) -> center(i)
        // Assumption: u(i) is left face, u(i+1) is right face
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                double val = 0.5 * (f(i, j) + f(i + 1, j));
                result.push_back(val);
            }
        }
    }
    else {
        // Fallback (e.g. FACE_Z): just take current index for now
        // A proper implementation would average (j) and (j+1)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                result.push_back(f(i, j));
            }
        }
    }
    return result;
}

void write_local_vts(int step_index, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/processor" << mesh.rank << "/"
       << cfg.filename_prefix << "_" << step_index << "_" << mesh.rank << ".vts";

    std::ofstream file(ss.str());

    int t_start = mesh.offset_theta;
    int t_end = mesh.offset_theta + mesh.n_theta_local;
    int z_start = mesh.offset_z;
    int z_end = mesh.offset_z + mesh.n_z_local;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <StructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\">\n";
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
            file << (mesh.R * std::cos(theta)) << " " << (mesh.R * std::sin(theta)) << " " << z << " ";
        }
    }
    file << "\n        </DataArray>\n      </Points>\n";

    // 2. Cell Data
    file << "      <CellData>\n";
    for(auto& [name, field_ptr] : fields) {
         // Interpolate/Copy data
         std::vector<double> centered = get_centered_data(*field_ptr, mesh);

         file << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
         for (const auto& val : centered) file << val << " ";
         file << "\n        </DataArray>\n";
    }
    file << "      </CellData>\n";
    file << "    </Piece>\n  </StructuredGrid>\n</VTKFile>\n";
}

void write_pvts(int step_index, int total_ranks, const SimulationConfig& cfg, Fields& fields) {
    std::stringstream ss;
    ss << cfg.output_dir << "/" << cfg.filename_prefix << "_" << step_index << ".pvts";
    std::ofstream file(ss.str());

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData>\n";

    // List all fields so ParaView knows what to expect
    for(auto& [name, field] : fields) {
        file << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"/>\n";
    }

    file << "    </PCellData>\n";

    // Write Piece entries for each rank
    for (int r = 0; r < total_ranks; ++r) {
        int n_local = cfg.n_theta_global / total_ranks;
        int remainder = cfg.n_theta_global % total_ranks;
        int offset = r * n_local;
        if (r < remainder) { n_local++; offset += r; } else { offset += remainder; }

        file << "    <Piece Extent=\"" << offset << " " << (offset + n_local) << " 0 " << cfg.n_z_global << " 0 0\" "
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

    // 1. Rank 0 cleans up and creates the root folder
    if (rank == 0) {
        if (fs::exists(cfg.output_dir)) {
            fs::remove_all(cfg.output_dir);
        }
        fs::create_directory(cfg.output_dir);

        // Initialize PVD file
        std::ofstream pvd(cfg.output_dir + "/results.pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n";
    }

    // 2. Wait for Rank 0 to finish filesystem ops
    MPI_Barrier(comm);

    // 3. Each rank creates its own subfolder
    std::string my_subfolder = cfg.output_dir + "/processor" + std::to_string(rank);
    if (!fs::exists(my_subfolder)) {
        fs::create_directory(my_subfolder);
    }

    // 4. Final sync to ensure all folders exist before any writing happens
    MPI_Barrier(comm);
}

void write_timestep(double time, int step_index, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    write_local_vts(step_index, mesh, fields, cfg);
    if (mesh.rank == 0) {
        write_pvts(step_index, mesh.size, cfg, fields);
        update_pvd(time, step_index, cfg);
    }
}

} // namespace IO
