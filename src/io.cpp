#include "io.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "output_naming.hpp"

namespace fs = std::filesystem;

namespace IO {
namespace {
// Selection key (what config output_fields / output_field_enabled use). theta is
// selected as film_content for backward compatibility.
std::string selection_key_for_field(const std::string& name) {
    return (name == "theta") ? "film_content" : name;
}

bool vector_group_enabled(const SimulationConfig& cfg, const OutputNaming::VectorGroup& g) {
    return cfg.output_field_enabled(g.x) || cfg.output_field_enabled(g.y) ||
           cfg.output_field_enabled(g.z);
}
}

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

void write_vector_data_array(std::ofstream& file, const std::string& name,
                             const std::vector<double>& x, const std::vector<double>& y,
                             const std::vector<double>& z) {
    file << "        <DataArray type=\"Float64\" Name=\"" << name
         << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t idx = 0; idx < x.size(); ++idx) {
        file << x[idx] << " " << y[idx] << " " << z[idx] << " ";
    }
    file << "\n        </DataArray>\n";
}

// Flat (unwrapped theta-z plane): the natural Cartesian vector is (u_theta, u_z, 0).
void write_flat_velocity_arrays(std::ofstream& file, const Mesh& mesh, Fields& fields) {
    std::vector<double> u_theta = get_centered_data(fields["velocity_theta"], mesh);
    std::vector<double> u_z = get_centered_data(fields["velocity_z"], mesh);
    std::vector<double> zero(u_theta.size(), 0.0);
    write_vector_data_array(file, OutputNaming::velocity_name(), u_theta, u_z, zero);
}

// Curved (real cylinder): Cartesian (vx, vy, vz) from the curvilinear components.
void write_curved_velocity_arrays(std::ofstream& file, const Mesh& mesh, Fields& fields) {
    std::vector<double> u_theta = get_centered_data(fields["velocity_theta"], mesh);
    std::vector<double> u_z = get_centered_data(fields["velocity_z"], mesh);
    std::vector<double> vx(u_theta.size()), vy(u_theta.size());
    size_t idx = 0;
    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
            vx[idx] = -u_theta[idx] * std::sin(theta);
            vy[idx] =  u_theta[idx] * std::cos(theta);
            ++idx;
        }
    }
    write_vector_data_array(file, OutputNaming::velocity_name(), vx, vy, u_z);
}

// Cartesian force / load / bearing resultant vectors (constant fields).
void write_vector_groups(std::ofstream& file, const Mesh& mesh, Fields& fields,
                         const SimulationConfig& cfg) {
    for (const auto& g : OutputNaming::vector_groups()) {
        if (!vector_group_enabled(cfg, g)) continue;
        if (!fields.has(g.x) || !fields.has(g.y) || !fields.has(g.z)) continue;
        write_vector_data_array(file, g.out,
                                get_centered_data(fields[g.x], mesh),
                                get_centered_data(fields[g.y], mesh),
                                get_centered_data(fields[g.z], mesh));
    }
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

    const bool write_velocity = cfg.output_field_enabled("velocity") &&
        fields.has("velocity_theta") && fields.has("velocity_z");

    // Cell data (Same as 3D mesh)
    file << "      <CellData";
    if (write_velocity) file << " Vectors=\"" << OutputNaming::velocity_name() << "\"";
    file << ">\n";
    for(auto& [name, field_ptr] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        if (OutputNaming::is_consumed_component(name)) continue;
        if (!cfg.output_field_enabled(selection_key_for_field(name))) continue;
        std::vector<double> centered = get_centered_data(*field_ptr, mesh);
        file << "        <DataArray type=\"Float64\" Name=\"" << OutputNaming::scalar(name) << "\" format=\"ascii\">\n";
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

    write_vector_groups(file, mesh, fields, cfg);
    if (write_velocity) write_flat_velocity_arrays(file, mesh, fields);

    file << "      </CellData>\n    </Piece>\n  </StructuredGrid>\n</VTKFile>\n";
}

void write_flat_pvts(int step, int total_ranks, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/flat/" << cfg.filename_prefix << "_" << step << ".pvts";
    std::string path = ss.str();
    std::ofstream file(path);
    const bool write_velocity = cfg.output_field_enabled("velocity") &&
        fields.has("velocity_theta") && fields.has("velocity_z");

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData";
    if (write_velocity) file << " Vectors=\"" << OutputNaming::velocity_name() << "\"";
    file << ">\n";

    for(auto& [name, field] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        if (OutputNaming::is_consumed_component(name)) continue;
        if (!cfg.output_field_enabled(selection_key_for_field(name))) continue;
        file << "      <PDataArray type=\"Float64\" Name=\"" << OutputNaming::scalar(name) << "\"/>\n";
    }
    file << "      <PDataArray type=\"Float64\" Name=\"theta_rad\"/>\n"
         << "      <PDataArray type=\"Float64\" Name=\"z_m\"/>\n";
    for (const auto& g : OutputNaming::vector_groups()) {
        if (vector_group_enabled(cfg, g) && fields.has(g.x) && fields.has(g.y) && fields.has(g.z))
            file << "      <PDataArray type=\"Float64\" Name=\"" << g.out << "\" NumberOfComponents=\"3\"/>\n";
    }
    if (write_velocity) {
        file << "      <PDataArray type=\"Float64\" Name=\"" << OutputNaming::velocity_name()
             << "\" NumberOfComponents=\"3\"/>\n";
    }
    file << "    </PCellData>\n";

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

    const bool write_velocity = cfg.output_field_enabled("velocity") &&
        fields.has("velocity_theta") && fields.has("velocity_z");

    // Cell data
    file << "      <CellData";
    if (write_velocity) file << " Vectors=\"" << OutputNaming::velocity_name() << "\"";
    file << ">\n";

    // Write scalar fields (vector components are grouped into Cartesian vectors below)
    for(auto& [name, field_ptr] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        if (OutputNaming::is_consumed_component(name)) continue;
        if (!cfg.output_field_enabled(selection_key_for_field(name))) continue;

        std::vector<double> centered = get_centered_data(*field_ptr, mesh);

        file << "        <DataArray type=\"Float64\" Name=\"" << OutputNaming::scalar(name) << "\" format=\"ascii\">\n";
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

    write_vector_groups(file, mesh, fields, cfg);
    if (write_velocity) write_curved_velocity_arrays(file, mesh, fields);

    file << "      </CellData>\n"
         << "    </Piece>\n"
         << "  </StructuredGrid>\n"
         << "</VTKFile>\n";
}

void write_pvts(int step, int total_ranks, Fields& fields, const SimulationConfig& cfg) {
    std::stringstream ss;
    ss << cfg.output_dir << "/" << cfg.filename_prefix << "_" << step << ".pvts";

    std::ofstream file(ss.str());
    const bool write_velocity = cfg.output_field_enabled("velocity") &&
        fields.has("velocity_theta") && fields.has("velocity_z");

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PStructuredGrid WholeExtent=\"0 " << cfg.n_theta_global << " 0 " << cfg.n_z_global << " 0 0\" GhostLevel=\"0\">\n"
         << "    <PPoints>\n"
         << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
         << "    </PPoints>\n"
         << "    <PCellData";
    if (write_velocity) file << " Vectors=\"" << OutputNaming::velocity_name() << "\"";
    file << ">\n";

    for(auto& [name, field] : fields) {
        if (name == "velocity_theta" || name == "velocity_z") continue;
        if (OutputNaming::is_consumed_component(name)) continue;
        if (!cfg.output_field_enabled(selection_key_for_field(name))) continue;
        file << "      <PDataArray type=\"Float64\" Name=\"" << OutputNaming::scalar(name) << "\"/>\n";
    }
    // Explicit coordinates
    file << "      <PDataArray type=\"Float64\" Name=\"theta_rad\"/>\n"
         << "      <PDataArray type=\"Float64\" Name=\"z_m\"/>\n";

    // Cartesian vector field headers
    for (const auto& g : OutputNaming::vector_groups()) {
        if (vector_group_enabled(cfg, g) && fields.has(g.x) && fields.has(g.y) && fields.has(g.z))
            file << "      <PDataArray type=\"Float64\" Name=\"" << g.out << "\" NumberOfComponents=\"3\"/>\n";
    }
    if (write_velocity) {
        file << "      <PDataArray type=\"Float64\" Name=\"" << OutputNaming::velocity_name()
             << "\" NumberOfComponents=\"3\"/>\n";
    }

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

    std::string setup_error;
    if (rank == 0) {
        std::error_code ec;
        const fs::path root(cfg.output_dir);

        if (fs::exists(root, ec)) {
            fs::remove_all(root, ec);
            if (ec) {
                setup_error = "Cannot clear output_dir '" + root.string() + "': " +
                    ec.message() + ". Close any viewer using the results, or choose a new output_dir.";
            }
        } else if (ec) {
            setup_error = "Cannot inspect output_dir '" + root.string() + "': " + ec.message();
        }

        if (setup_error.empty()) {
            fs::create_directories(root, ec);
            if (ec) setup_error = "Cannot create output_dir '" + root.string() + "': " + ec.message();
        }

        if (setup_error.empty() && cfg.output_write_3d) {
            std::ofstream pvd(root / "results.pvd");
            if (!pvd) {
                setup_error = "Cannot write '" + (root / "results.pvd").string() + "'.";
            } else {
                pvd << "<?xml version=\"1.0\"?>\n"
                    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                    << "  <Collection>\n"
                    << "  </Collection>\n"
                    << "</VTKFile>\n";
            }
        }

        if (setup_error.empty() && cfg.output_write_flat) {
            const fs::path flat_dir = root / "flat";
            fs::create_directories(flat_dir, ec);
            if (ec) {
                setup_error = "Cannot create output_dir '" + flat_dir.string() + "': " + ec.message();
            } else {
                std::ofstream flat_pvd(flat_dir / "results.pvd");
                if (!flat_pvd) {
                    setup_error = "Cannot write '" + (flat_dir / "results.pvd").string() + "'.";
                } else {
                    flat_pvd << "<?xml version=\"1.0\"?>\n"
                             << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                             << "  <Collection>\n"
                             << "  </Collection>\n"
                             << "</VTKFile>\n";
                }
            }
        }
    }

    int root_error_len = static_cast<int>(setup_error.size());
    MPI_Bcast(&root_error_len, 1, MPI_INT, 0, comm);
    if (root_error_len > 0) {
        if (rank != 0) setup_error.assign(static_cast<size_t>(root_error_len), '\0');
        MPI_Bcast(setup_error.data(), root_error_len, MPI_CHAR, 0, comm);
        throw std::runtime_error(setup_error);
    }
    MPI_Barrier(comm);

    // Subdirectory for each processor
    std::string rank_error;
    std::error_code ec;
    if (cfg.output_write_3d) {
        const fs::path dir = fs::path(cfg.output_dir) / ("processor" + std::to_string(rank));
        fs::create_directories(dir, ec);
        if (ec) rank_error = "Cannot create output_dir '" + dir.string() + "': " + ec.message();
    }
    if (rank_error.empty() && cfg.output_write_flat) {
        const fs::path dir = fs::path(cfg.output_dir) / "flat" / ("processor" + std::to_string(rank));
        fs::create_directories(dir, ec);
        if (ec) rank_error = "Cannot create output_dir '" + dir.string() + "': " + ec.message();
    }

    const int local_failed = rank_error.empty() ? 0 : 1;
    int any_failed = 0;
    MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, comm);
    if (any_failed) {
        if (!rank_error.empty()) throw std::runtime_error(rank_error);
        throw std::runtime_error("Cannot create one or more processor output directories.");
    }
    MPI_Barrier(comm);
}

void write_timestep(double time, int step, const Mesh& mesh, Fields& fields, const SimulationConfig& cfg) {
    // 3D output
    if (cfg.output_write_3d) {
        write_vts(step, mesh, fields, cfg);
        if (mesh.rank == 0) {
            write_pvts(step, mesh.size, fields, cfg);
            update_pvd(time, step, cfg);
        }
    }

    // Flat output
    if (cfg.output_write_flat) {
        write_flat_vts(step, mesh, fields, cfg);
        if (mesh.rank == 0) {
            write_flat_pvts(step, mesh.size, fields, cfg);
            update_flat_pvd(time, step, cfg);
        }
    }
}
}
