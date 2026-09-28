#include "io.hpp"

#include <cmath>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <array>
#include <iomanip>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <string>

#include "output_naming.hpp"

namespace fs = std::filesystem;

namespace IO {
namespace {
constexpr std::array<const char*, 4> initial_state_fields = {
    "pressure", "theta", "temperature", "cross_film_collocation"};

std::vector<std::string> initial_state_names(const Fields& fields) {
    std::vector<std::string> names(initial_state_fields.begin(), initial_state_fields.end());
    for (int layer = 0; ; ++layer) {
        const std::string name = "cross_film_layer_" + std::to_string(layer);
        if (!fields.has(name)) break;
        names.push_back(name);
    }
    return names;
}

void initial_state_error(std::string error, const Mesh& mesh) {
    int length = mesh.rank == 0 ? static_cast<int>(error.size()) : 0;
    MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    error.resize(length);
    if (length) {
        MPI_Bcast(error.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
        throw std::runtime_error("Steady initial state: " + error);
    }
}

void initial_state_layout(const Mesh& mesh, std::vector<int>& counts, std::vector<int>& offsets) {
    const int count = mesh.n_theta_local * mesh.n_z_local;
    const int offset = mesh.offset_theta * mesh.n_z_local;
    counts.resize(mesh.size);
    offsets.resize(mesh.size);
    MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, MPI_COMM_WORLD);
}

bool admissible_initial_value(int field, double value, const SimulationConfig& cfg) {
    return std::isfinite(value) &&
        (field != 0 || value >= cfg.effective_p_cav() - 1e-8) &&
        (field != 1 || (value >= 0.0 && value <= 1.0 + 1e-12)) &&
        ((field != 2 && field < 4) || value > 0.0);
}

bool processor_directory(const std::string& name) {
    return name.starts_with("processor") && name.size() > 9 &&
        std::all_of(name.begin() + 9, name.end(), [](unsigned char c) { return std::isdigit(c); });
}

// Remove only products with the solver's naming convention. A rerun must not
// recursively erase reports, input files or unrelated content in output_dir.
void clear_solver_products(const fs::path& root, const std::string& prefix, std::error_code& ec) {
    fs::recursive_directory_iterator it(root, fs::directory_options::none, ec), end;
    for (; !ec && it != end; it.increment(ec)) {
        if (it->is_symlink(ec)) { it.disable_recursion_pending(); continue; }
        if (it->is_directory(ec)) {
            const std::string name = it->path().filename().string();
            if (!processor_directory(name) && !(it.depth() == 0 && name == "flat"))
                it.disable_recursion_pending();
            continue;
        }
        const std::string name = it->path().filename().string();
        const std::string stem = it->path().stem().string();
        const std::string ext = it->path().extension().string();
        bool generated = name == "results.pvd" || (it.depth() == 0 && name == "diagnostics.csv");
        if ((ext == ".vts" || ext == ".pvts") && stem.starts_with(prefix + "_")) {
            const std::string suffix = stem.substr(prefix.size() + 1);
            generated = !suffix.empty() && std::all_of(suffix.begin(), suffix.end(),
                [](unsigned char c) { return std::isdigit(c) || c == '_'; });
        }
        if (generated) fs::remove(it->path(), ec);
    }
}

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

void write_steady_initial_state(const fs::path& path, const Fields& fields,
                                const Mesh& mesh, const SimulationConfig& cfg) {
    std::vector<int> counts, offsets;
    initial_state_layout(mesh, counts, offsets);
    const size_t total = static_cast<size_t>(mesh.n_theta_global) * mesh.n_z_global;
    const auto names = initial_state_names(fields);
    std::vector<std::vector<double>> gathered(names.size());
    int invalid = 0;
    for (size_t k = 0; k < names.size(); ++k) {
        std::vector<double> local(counts[mesh.rank]);
        const std::string& name = names[k];
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j) {
            const double value = fields.has(name) ? fields[name](i, j) : 0.0;
            local[static_cast<size_t>(i) * mesh.n_z_local + j] = value;
            if (!admissible_initial_value(k, value, cfg)) invalid = 1;
        }
        if (mesh.rank == 0) gathered[k].resize(total);
        MPI_Gatherv(local.data(), counts[mesh.rank], MPI_DOUBLE, gathered[k].data(),
                    counts.data(), offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    int any_invalid = 0;
    MPI_Allreduce(&invalid, &any_invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    std::string error;
    if (mesh.rank == 0) {
        if (any_invalid) error = "refusing to write inadmissible fields";
        else {
            std::ofstream stream(path);
            stream << std::setprecision(17) << "PANCAKE_STEADY_GUESS " << (names.size() > 4 ? 2 : 1) << '\n'
                   << mesh.n_theta_global << ' ' << mesh.n_z_global << ' '
                   << cfg.R << ' ' << cfg.c << ' ' << cfg.L << '\n';
            if (names.size() > 4) stream << names.size() << '\n';
            for (size_t k = 0; k < names.size(); ++k) {
                stream << names[k] << '\n';
                for (double value : gathered[k]) stream << value << '\n';
            }
            stream.close();
            if (!stream) error = "could not write " + path.string();
        }
    }
    initial_state_error(error, mesh);
}

void read_steady_initial_state(const fs::path& path, Fields& fields,
                               const Mesh& mesh, const SimulationConfig& cfg) {
    const size_t total = static_cast<size_t>(mesh.n_theta_global) * mesh.n_z_global;
    const auto names = initial_state_names(fields);
    std::vector<std::vector<double>> gathered(names.size());
    int stored_count = 4;
    std::string error;
    if (mesh.rank == 0) {
        std::ifstream stream(path);
        std::string magic, name;
        int version = 0, nt = 0, nz = 0;
        double radius = 0.0, clearance = 0.0, length = 0.0;
        stream >> magic >> version >> nt >> nz >> radius >> clearance >> length;
        if (!stream || magic != "PANCAKE_STEADY_GUESS" || (version != 1 && version != 2) ||
            nt != mesh.n_theta_global || nz != mesh.n_z_global ||
            radius != cfg.R || clearance != cfg.c || length != cfg.L) {
            error = "missing, malformed, or incompatible file " + path.string();
        }
        if (error.empty() && version == 2) {
            stream >> stored_count;
            if (!stream || stored_count != static_cast<int>(names.size())) error = "incompatible cross-film layer count";
        }
        for (int k = 0; k < stored_count && error.empty(); ++k) {
            stream >> name;
            if (name != names[k]) { error = "unexpected field " + name; break; }
            gathered[k].resize(total);
            for (double& value : gathered[k]) {
                if (!(stream >> value) || !admissible_initial_value(k, value, cfg)) {
                    error = "invalid or missing value in " + name;
                    break;
                }
            }
        }
        if (error.empty() && (stream >> name)) error = "unexpected trailing data";
    }
    initial_state_error(error, mesh);
    MPI_Bcast(&stored_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> counts, offsets;
    initial_state_layout(mesh, counts, offsets);
    for (int k = 0; k < stored_count; ++k) {
        std::vector<double> local(counts[mesh.rank]);
        MPI_Scatterv(gathered[k].data(), counts.data(), offsets.data(), MPI_DOUBLE,
                     local.data(), counts[mesh.rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!fields.has(names[k])) continue;
        Field& field = fields[names[k]];
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j)
            field(i, j) = local[static_cast<size_t>(i) * mesh.n_z_local + j];
    }
    // Legacy mean-only states seed every layer uniformly; they are initial guesses.
    for (size_t k = stored_count; k < names.size(); ++k)
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j)
            fields[names[k]](i, j) = fields["temperature"](i, j);
    if (fields.has("cross_film_layer_initialized")) {
        fields["cross_film_layer_initialized"].fill(1.0);
        for (const char* name : {"cross_film_t_shaft", "cross_film_t_bush"})
            if (fields.has(name))
                for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j)
                    fields[name](i,j) = fields["temperature"](i,j);
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
    file.precision(17);
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
    file.precision(17);
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
    file.precision(17);
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
    file.precision(17);
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
            clear_solver_products(root, cfg.filename_prefix, ec);
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
