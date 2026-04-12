#include "io.hpp"

#include <filesystem>
#include <fstream>
#include <mpi.h>
#include <sstream>

namespace fs = std::filesystem;

namespace IO {
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
        // u(i) is west face, u(i+1) is east face
        // We use field(i+1, j) which is a ghost cell if i == n_theta-1
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                double val = 0.5 * (f(i, j) + f(i + 1, j));
                result.push_back(val);
            }
        }
    }
    else if (f.location == GridLocation::FACE_Z) {
        // Interpolate v(j) and v(j+1) -> center(j)
        // v(j) is south face, v(j+1) is north face
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                double val = 0.5 * (f(i, j) + f(i, j + 1));
                result.push_back(val);
            }
        }
    }
    return result;
}

void update_pvd(double time, int step, const SimulationConfig& cfg) {
    // We rewrite the footer every time to keep the XML valid
    // This is slightly less efficient than keeping it open but allows live viewing.
    std::ofstream file(cfg.output_dir + "/results.pvd", std::ios::in | std::ios::out | std::ios::ate);
    
    // Move back before the closing tags
    // <Collection> (13) + </Collection> (14) + </VTKFile> (11) + newlines
    // Easier approach: Just keep the file valid by closing it every time.
    
    // Read the whole file except the last two lines
    std::ifstream infile(cfg.output_dir + "/results.pvd");
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.find("</Collection>") == std::string::npos && 
            line.find("</VTKFile>") == std::string::npos) {
            lines.push_back(line);
        }
    }
    infile.close();

    std::ofstream outfile(cfg.output_dir + "/results.pvd");
    for (const auto& l : lines) outfile << l << "\n";
    
    outfile << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" "
            << "file=\"" << cfg.filename_prefix << "_" << step << ".pvts\"/>\n"
            << "  </Collection>\n"
            << "</VTKFile>\n";
}

void update_flat_pvd(double time, int step, const SimulationConfig& cfg) {
    std::string path = cfg.output_dir + "/flat/results.pvd";
    std::ifstream infile(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(infile, line)) {
        if (line.find("</Collection>") == std::string::npos && 
            line.find("</VTKFile>") == std::string::npos) {
            lines.push_back(line);
        }
    }
    infile.close();

    std::ofstream outfile(path);
    for (const auto& l : lines) outfile << l << "\n";
    outfile << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" "
            << "file=\"" << cfg.filename_prefix << "_" << step << ".pvts\"/>\n"
            << "  </Collection>\n"
            << "</VTKFile>\n";
}

void write_flat_vts(int step, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/flat/processor" << mesh.rank << "/" << cfg.filename_prefix << "_" << step << "_" << mesh.rank << ".vts";

    std::ofstream file(ss.str());
    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <StructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\">\n"
         << "    <Piece Extent=\"" << mesh.offset_theta << " " << mesh.offset_theta+mesh.n_theta_local << " 0 " << mesh.n_z_local << " 0 0\">\n";

    // Points (Flat: X = theta*R, Y = z, Z = 0)
    file << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int j = 0; j<=mesh.n_z_local; ++j) {
        double z = j * mesh.get_d_z();
        for(int i = 0; i<=mesh.n_theta_local; ++i) {
            double theta = (mesh.offset_theta + i) * mesh.get_d_theta();
            file << (mesh.R * theta) << " " << z << " 0 ";
        }
    }
    file << "\n        </DataArray>\n      </Points>\n";

    // Cell data (Same as 3D mesh)
    file << "      <CellData>\n";
    for(auto& [name, field_ptr] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        std::vector<double> centered = get_centered_data(*field_ptr, mesh);
        std::string output_name = (name == "theta") ? "film_content" : name;
        file << "        <DataArray type=\"Float64\" Name=\"" << output_name << "\" format=\"ascii\">\n";
        for (const auto& val : centered) file << val << " ";
        file << "\n        </DataArray>\n";
    }
    
    // Coordinates
    file << "        <DataArray type=\"Float64\" Name=\"theta_rad\" format=\"ascii\">\n";
    for (int j = 0; j < mesh.n_z_local; ++j)
        for (int i = 0; i < mesh.n_theta_local; ++i)
            file << (mesh.offset_theta + i + 0.5) * mesh.get_d_theta() << " ";
    file << "\n        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"z_m\" format=\"ascii\">\n";
    for (int j = 0; j < mesh.n_z_local; ++j)
        for (int i = 0; i < mesh.n_theta_local; ++i)
            file << (j + 0.5) * mesh.get_d_z() << " ";
    file << "\n        </DataArray>\n";

    // Velocity vectors (In flat mesh, we can export local components as UVW)
    if (fields.has("velocity_theta") && fields.has("velocity_z")) {
        std::vector<double> u_th = get_centered_data(fields["velocity_theta"], mesh);
        std::vector<double> u_z  = get_centered_data(fields["velocity_z"], mesh);
        file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t idx = 0; idx < u_th.size(); ++idx) {
            file << u_th[idx] << " " << u_z[idx] << " 0 ";
        }
        file << "\n        </DataArray>\n";
    }

    file << "      </CellData>\n    </Piece>\n  </StructuredGrid>\n</VTKFile>\n";
}

void write_flat_pvts(int step, int total_ranks, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/flat/" << cfg.filename_prefix << "_" << step << ".pvts";
    std::string path = ss.str();
    std::ofstream file(path);
    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData>\n";

    for(auto& [name, field] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        std::string output_name = (name == "theta") ? "film_content" : name;
        file << "      <PDataArray type=\"Float64\" Name=\"" << output_name << "\"/>\n";
    }
    file << "      <PDataArray type=\"Float64\" Name=\"theta_rad\"/>\n"
         << "      <PDataArray type=\"Float64\" Name=\"z_m\"/>\n"
         << "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n"
         << "    </PCellData>\n";

    for (int r = 0; r < total_ranks; ++r) {
        int n_local = cfg.n_theta_global / total_ranks;
        int remainder = cfg.n_theta_global % total_ranks;
        int offset = r * n_local;
        if (r < remainder) { n_local++; offset += r; } else { offset += remainder; }

        file << "    <Piece Extent=\"" << offset << " " << (offset + n_local) << " 0 " << cfg.n_z_global << " 0 0\" "
             << "Source=\"processor" << r << "/" << cfg.filename_prefix << "_" << step << "_" << r << ".vts\"/>\n";
    }
    file << "  </PStructuredGrid>\n</VTKFile>\n";
}

void write_vts(int step, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/processor" << mesh.rank << "/" << cfg.filename_prefix << "_" << step << "_" << mesh.rank << ".vts";

    // Headers
    std::ofstream file(ss.str());
    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <StructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\">\n"
         << "    <Piece Extent=\"" << mesh.offset_theta << " " << mesh.offset_theta+mesh.n_theta_local << " 0 " << mesh.n_z_local << " 0 0\">\n";

    // Points
    file << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(int j = 0; j<=mesh.n_z_local; ++j) {
        double z = j * mesh.get_d_z();
        for(int i = 0; i<=mesh.n_theta_local; ++i) {
            double theta = (mesh.offset_theta + i) * mesh.get_d_theta();
            file << (mesh.R * cos(theta)) << " " << (mesh.R * sin(theta)) << " " << z << " ";
        }
    }
    file << "\n"
         << "        </DataArray>\n"
         << "      </Points>\n";

    // Cell data
    file << "      <CellData>\n";
    
    // Write scalar fields (excluding individual velocity components which we'll group)
    for(auto& [name, field_ptr] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;

        std::vector<double> centered = get_centered_data(*field_ptr, mesh);
        
        // Use 'film_content' for Elrod-Adams 'theta' to avoid coordinate confusion
        std::string output_name = (name == "theta") ? "film_content" : name;

        file << "        <DataArray type=\"Float64\" Name=\"" << output_name << "\" format=\"ascii\">\n";
        for (const auto& val : centered) file << val << " ";
        file << "\n        </DataArray>\n";
    }

    // Write explicit coordinates as scalar fields for easy unwrapping/plotting
    file << "        <DataArray type=\"Float64\" Name=\"theta_rad\" format=\"ascii\">\n";
    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            file << (mesh.offset_theta + i + 0.5) * mesh.get_d_theta() << " ";
        }
    }
    file << "\n        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"z_m\" format=\"ascii\">\n";
    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            file << (j + 0.5) * mesh.get_d_z() << " ";
        }
    }
    file << "\n        </DataArray>\n";

    // Write Velocity as a Vector field
    if (fields.has("velocity_theta") && fields.has("velocity_z")) {
        std::vector<double> u_th = get_centered_data(fields["velocity_theta"], mesh);
        std::vector<double> u_z  = get_centered_data(fields["velocity_z"], mesh);

        file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        int idx = 0;
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                double theta = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
                double vx = -u_th[idx] * sin(theta);
                double vy =  u_th[idx] * cos(theta);
                double vz =  u_z[idx];
                file << vx << " " << vy << " " << vz << " ";
                idx++;
            }
        }
        file << "\n        </DataArray>\n";
    }

    file << "      </CellData>\n"
         << "    </Piece>\n"
         << "  </StructuredGrid>\n"
         << "</VTKFile>\n";
}

void write_pvts(int step, int total_ranks, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/" << cfg.filename_prefix << "_" << step << ".pvts";

    std::ofstream file(ss.str());
    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData>\n";

    for(auto& [name, field] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        std::string output_name = (name == "theta") ? "film_content" : name;
        file << "      <PDataArray type=\"Float64\" Name=\"" << output_name << "\"/>\n";
    }
    // Explicit coordinates
    file << "      <PDataArray type=\"Float64\" Name=\"theta_rad\"/>\n"
         << "      <PDataArray type=\"Float64\" Name=\"z_m\"/>\n";

    // Vector field header
    file << "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";

    file << "    </PCellData>\n";

    // Write Piece entries for each rank
    for (int r = 0; r < total_ranks; ++r) {
        int n_local = cfg.n_theta_global / total_ranks;
        int remainder = cfg.n_theta_global % total_ranks;
        int offset = r * n_local;
        if (r < remainder) { n_local++; offset += r; } else { offset += remainder; }

        file << "    <Piece Extent=\"" << offset << " " << (offset + n_local) << " 0 " << cfg.n_z_global << " 0 0\" "
             << "Source=\"processor" << r << "/" << cfg.filename_prefix << "_" << step << "_" << r << ".vts\"/>\n";
    }
    file << "  </PStructuredGrid>\n"
         << "</VTKFile>\n";
}

void prepare_output_directory(const SimulationConfig& cfg, MPI_Comm comm) {
    int rank; MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        // Directory and file initialization
        if (fs::exists(cfg.output_dir)) fs::remove_all(cfg.output_dir);
        fs::create_directory(cfg.output_dir);
        fs::create_directory(cfg.output_dir + "/flat");

        // 3D PVD
        std::ofstream pvd(cfg.output_dir + "/results.pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n"
            << "  </Collection>\n"
            << "</VTKFile>\n";

        // Flat PVD
        std::ofstream flat_pvd(cfg.output_dir + "/flat/results.pvd");
        flat_pvd << "<?xml version=\"1.0\"?>\n"
                 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                 << "  <Collection>\n"
                 << "  </Collection>\n"
                 << "</VTKFile>\n";
    }
    MPI_Barrier(comm);

    // Subdirectory for each processor
    fs::create_directory(cfg.output_dir + "/processor" + std::to_string(rank));
    fs::create_directory(cfg.output_dir + "/flat/processor" + std::to_string(rank));
    MPI_Barrier(comm);
}

void write_timestep(double time, int step, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    // 3D output
    write_vts(step, mesh, fields, cfg);
    if (mesh.rank == 0) {
        write_pvts(step, mesh.size, fields, cfg);
        update_pvd(time, step, cfg);
    }

    // Flat output
    write_flat_vts(step, mesh, fields, cfg);
    if (mesh.rank == 0) {
        write_flat_pvts(step, mesh.size, fields, cfg);
        update_flat_pvd(time, step, cfg);
    }
}
}
