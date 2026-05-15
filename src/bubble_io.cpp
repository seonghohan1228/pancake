#include "bubble_io.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace BubbleIO {

// --- Internal helpers ---

static std::string rank_dir(const BubbleConfig& cfg, int rank) {
    return cfg.output_dir + "/rank_" + std::to_string(rank);
}

static std::string vtu_path(const BubbleConfig& cfg, int rank, int step) {
    return rank_dir(cfg, rank) + "/step_" + std::to_string(step) + ".vtu";
}

static std::string pvtu_path(const BubbleConfig& cfg, int step) {
    return cfg.output_dir + "/step_" + std::to_string(step) + ".pvtu";
}

static std::string pvtu_relative(const BubbleConfig& cfg, int rank, int step) {
    return "rank_" + std::to_string(rank) + "/step_" + std::to_string(step) + ".vtu";
}

static void update_pvd(double time, int step, const BubbleConfig& cfg) {
    const std::string path = cfg.output_dir + "/results.pvd";

    // Re-read, strip closing tags, append new entry
    std::ifstream in(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("</Collection>") == std::string::npos &&
            line.find("</VTKFile>")    == std::string::npos)
            lines.push_back(line);
    }
    in.close();

    std::ofstream out(path);
    for (const auto& l : lines) out << l << "\n";
    out << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" "
        << "file=\"step_" << step << ".pvtu\"/>\n"
        << "  </Collection>\n"
        << "</VTKFile>\n";
}

static void write_vtu(int step, const StereoMesh& mesh, const Mask& mask,
                      Fields& fields, const BubbleConfig& cfg)
{
    // Collect active cell data for this rank
    struct Cell {
        Vec3 p[4];  // SW, SE, NE, NW corners in 3D
        int  i, j;
    };

    std::vector<Cell> cells;
    cells.reserve(mesh.n_u_local * mesh.n_v_global);

    for (int i = 0; i < mesh.n_u_local; ++i) {
        for (int j = 0; j < mesh.n_v_global; ++j) {
            if (!mask.is_active(i, j)) continue;

            const double u0 = mesh.vertex_u(i);
            const double u1 = u0 + mesh.du;
            const double v0 = mesh.vertex_v(j);
            const double v1 = v0 + mesh.dv;

            Cell c;
            c.p[0] = mesh.project(u0, v0);  // SW
            c.p[1] = mesh.project(u1, v0);  // SE
            c.p[2] = mesh.project(u1, v1);  // NE
            c.p[3] = mesh.project(u0, v1);  // NW
            c.i    = i;
            c.j    = j;
            cells.push_back(c);
        }
    }

    const int n_cells  = (int)cells.size();
    const int n_points = n_cells * 4;

    std::ofstream f(vtu_path(cfg, mesh.rank, step));
    if (!f.is_open())
        throw std::runtime_error("BubbleIO: cannot open " + vtu_path(cfg, mesh.rank, step));

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << n_points
      << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" Name=\"Points\""
      << " NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& c : cells)
        for (int k = 0; k < 4; ++k)
            f << c.p[k].x << " " << c.p[k].y << " " << c.p[k].z << " ";
    f << "\n        </DataArray>\n"
      << "      </Points>\n";

    // Cells
    f << "      <Cells>\n";

    // Connectivity: 0 1 2 3,  4 5 6 7, ...
    f << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int c = 0; c < n_cells; ++c)
        f << c*4 << " " << c*4+1 << " " << c*4+2 << " " << c*4+3 << " ";
    f << "\n        </DataArray>\n";

    // Offsets: 4, 8, 12, ...
    f << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int c = 1; c <= n_cells; ++c) f << c*4 << " ";
    f << "\n        </DataArray>\n";

    // Types: VTK_QUAD = 9
    f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int c = 0; c < n_cells; ++c) f << "9 ";
    f << "\n        </DataArray>\n";

    f << "      </Cells>\n";

    // Cell data
    f << "      <CellData>\n";
    for (const auto& [name, field_ptr] : fields) {
        const Field& fld = *field_ptr;
        if (fld.location != GridLocation::CENTER) continue;

        f << "        <DataArray type=\"Float64\" Name=\"" << name
          << "\" format=\"ascii\">\n";
        for (const auto& c : cells)
            f << fld(c.i, c.j) << " ";
        f << "\n        </DataArray>\n";
    }
    f << "      </CellData>\n";

    f << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";
}

static void write_pvtu(int step, int n_ranks, Fields& fields,
                       const BubbleConfig& cfg)
{
    std::ofstream f(pvtu_path(cfg, step));
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <PUnstructuredGrid GhostLevel=\"0\">\n"
      << "    <PPoints>\n"
      << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n"
      << "    </PPoints>\n"
      << "    <PCellData>\n";

    for (const auto& [name, field_ptr] : fields) {
        if (field_ptr->location != GridLocation::CENTER) continue;
        f << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"/>\n";
    }

    f << "    </PCellData>\n";

    for (int r = 0; r < n_ranks; ++r)
        f << "    <Piece Source=\"" << pvtu_relative(cfg, r, step) << "\"/>\n";

    f << "  </PUnstructuredGrid>\n"
      << "</VTKFile>\n";
}

// --- Public API ---

void prepare_output_directory(const BubbleConfig& cfg, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == 0) {
        if (fs::exists(cfg.output_dir)) fs::remove_all(cfg.output_dir);
        fs::create_directories(cfg.output_dir);

        std::ofstream pvd(cfg.output_dir + "/results.pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n"
            << "  </Collection>\n"
            << "</VTKFile>\n";
    }
    MPI_Barrier(comm);

    fs::create_directories(rank_dir(cfg, rank));
    MPI_Barrier(comm);
}

void write_timestep(double time, int step_index,
                    const StereoMesh& mesh, const Mask& mask,
                    Fields& fields, const BubbleConfig& cfg)
{
    write_vtu(step_index, mesh, mask, fields, cfg);

    if (mesh.rank == 0) {
        write_pvtu(step_index, mesh.size, fields, cfg);
        update_pvd(time, step_index, cfg);
    }
}

}  // namespace BubbleIO
