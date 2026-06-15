#include "results.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

std::string read_file(const fs::path& path, std::string& error) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        error = "Cannot open " + path.string();
        return {};
    }
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

/// Extracts the value of attribute `name` from the tag starting at `tag_pos`.
std::string tag_attribute(const std::string& xml, size_t tag_pos, const std::string& name) {
    const size_t tag_end = xml.find('>', tag_pos);
    if (tag_end == std::string::npos) return {};
    const std::string needle = name + "=\"";
    const size_t attr = xml.find(needle, tag_pos);
    if (attr == std::string::npos || attr > tag_end) return {};
    const size_t value_start = attr + needle.size();
    const size_t value_end = xml.find('"', value_start);
    if (value_end == std::string::npos) return {};
    return xml.substr(value_start, value_end - value_start);
}

/// Parses whitespace-separated doubles between `from` and the next closing tag.
bool parse_numbers(const std::string& xml, size_t from, size_t to,
                   std::vector<double>& out) {
    const char* p = xml.data() + from;
    const char* limit = xml.data() + to;
    while (p < limit) {
        while (p < limit && std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (p >= limit) break;
        char* end = nullptr;
        const double value = std::strtod(p, &end);
        if (end == p) return false;  // junk that is not a number
        out.push_back(value);
        p = end;
    }
    return true;
}

struct PieceData {
    int extent[6] = {0, 0, 0, 0, 0, 0};  // point extents: theta lo/hi, z lo/hi
    std::map<std::string, std::vector<double>> arrays;  // scalars + vector components
};

bool parse_extent(const std::string& text, int extent[6]) {
    std::stringstream ss(text);
    for (int i = 0; i < 6; ++i) {
        if (!(ss >> extent[i])) return false;
    }
    return true;
}

/// Parses the CellData arrays of one serial .vts piece.
bool parse_vts_piece(const fs::path& path, PieceData& piece, std::string& error) {
    std::string xml = read_file(path, error);
    if (xml.empty() && !error.empty()) return false;

    const size_t piece_tag = xml.find("<Piece");
    if (piece_tag == std::string::npos) {
        error = path.filename().string() + ": no <Piece> element";
        return false;
    }
    if (!parse_extent(tag_attribute(xml, piece_tag, "Extent"), piece.extent)) {
        error = path.filename().string() + ": bad Piece Extent";
        return false;
    }

    const size_t cell_data = xml.find("<CellData");
    if (cell_data == std::string::npos) {
        error = path.filename().string() + ": no <CellData> section";
        return false;
    }

    const int cells = (piece.extent[1] - piece.extent[0]) * (piece.extent[3] - piece.extent[2]);
    if (cells <= 0 || cells > 50'000'000) {
        error = path.filename().string() + ": implausible cell count";
        return false;
    }

    size_t cursor = cell_data;
    while (true) {
        const size_t array_tag = xml.find("<DataArray", cursor);
        if (array_tag == std::string::npos) break;
        const std::string name = tag_attribute(xml, array_tag, "Name");
        const std::string components = tag_attribute(xml, array_tag, "NumberOfComponents");
        const int n_components = components.empty() ? 1 : std::atoi(components.c_str());
        const size_t content_start = xml.find('>', array_tag);
        const size_t content_end = xml.find("</DataArray>", array_tag);
        if (content_start == std::string::npos || content_end == std::string::npos) {
            error = path.filename().string() + ": truncated DataArray '" + name + "'";
            return false;
        }
        cursor = content_end + 1;
        if (name.empty()) continue;

        std::vector<double> values;
        values.reserve(static_cast<size_t>(cells) * static_cast<size_t>(std::max(1, n_components)));
        if (!parse_numbers(xml, content_start + 1, content_end, values)) {
            error = path.filename().string() + ": non-numeric data in '" + name + "'";
            return false;
        }
        if (n_components == 1) {
            if (static_cast<int>(values.size()) != cells) {
                error = path.filename().string() + ": '" + name + "' has " +
                        std::to_string(values.size()) + " values, expected " +
                        std::to_string(cells);
                return false;
            }
            piece.arrays[name] = std::move(values);
        } else if (n_components == 3) {
            if (static_cast<int>(values.size()) != cells * 3) {
                error = path.filename().string() + ": vector '" + name + "' has wrong size";
                return false;
            }
            std::vector<double> x(cells), y(cells), z(cells);
            for (int cell = 0; cell < cells; ++cell) {
                x[cell] = values[cell * 3 + 0];
                y[cell] = values[cell * 3 + 1];
                z[cell] = values[cell * 3 + 2];
            }
            piece.arrays[name + ".x"] = std::move(x);
            piece.arrays[name + ".y"] = std::move(y);
            piece.arrays[name + ".z"] = std::move(z);
        }
    }
    return true;
}

}  // namespace

std::vector<StepInfo> parse_pvd(const fs::path& output_dir, const std::string& prefix,
                                std::string& error) {
    std::vector<StepInfo> steps;
    const fs::path pvd_path = output_dir / "results.pvd";
    std::string ignored;
    const std::string xml = read_file(pvd_path, ignored);

    if (!xml.empty()) {
        size_t cursor = 0;
        while (true) {
            const size_t dataset = xml.find("<DataSet", cursor);
            if (dataset == std::string::npos) break;
            cursor = dataset + 8;
            StepInfo info;
            info.time = std::atof(tag_attribute(xml, dataset, "timestep").c_str());
            const std::string file = tag_attribute(xml, dataset, "file");
            if (file.empty()) continue;
            // file name is "<prefix>_<step>.pvts"
            const size_t underscore = file.rfind('_');
            const size_t dot = file.rfind(".pvts");
            if (underscore == std::string::npos || dot == std::string::npos) continue;
            info.step = std::atoi(file.substr(underscore + 1, dot - underscore - 1).c_str());
            info.pvts = output_dir / file;
            steps.push_back(info);
        }
    }

    if (steps.empty()) {
        // No/empty pvd: fall back to globbing pvts files (e.g. run cancelled
        // before the collection footer was rewritten).
        std::error_code ec;
        for (const auto& entry : fs::directory_iterator(output_dir, ec)) {
            const std::string name = entry.path().filename().string();
            if (name.size() < prefix.size() + 6) continue;
            if (name.rfind(prefix + "_", 0) != 0 || entry.path().extension() != ".pvts") continue;
            StepInfo info;
            info.step = std::atoi(name.substr(prefix.size() + 1).c_str());
            info.time = static_cast<double>(info.step);
            info.pvts = entry.path();
            steps.push_back(info);
        }
        if (ec) error = "Cannot list " + output_dir.string();
    }

    std::sort(steps.begin(), steps.end(),
              [](const StepInfo& a, const StepInfo& b) { return a.step < b.step; });
    if (steps.empty() && error.empty()) {
        error = "No results found in " + output_dir.string();
    }
    return steps;
}

std::shared_ptr<ResultGrid> parse_step_grid(const StepInfo& info, std::string& error) {
    std::string xml = read_file(info.pvts, error);
    if (xml.empty()) {
        if (error.empty()) error = info.pvts.string() + " is empty";
        return nullptr;
    }

    const size_t grid_tag = xml.find("<PStructuredGrid");
    int whole[6] = {};
    if (grid_tag == std::string::npos ||
        !parse_extent(tag_attribute(xml, grid_tag, "WholeExtent"), whole)) {
        error = info.pvts.filename().string() + ": missing/invalid WholeExtent";
        return nullptr;
    }

    auto grid = std::make_shared<ResultGrid>();
    grid->nx = whole[1] - whole[0];
    grid->ny = whole[3] - whole[2];
    grid->step = info.step;
    grid->time = info.time;
    if (grid->nx <= 0 || grid->ny <= 0 || grid->nx > 100000 || grid->ny > 100000) {
        error = info.pvts.filename().string() + ": implausible grid extent";
        return nullptr;
    }

    // Collect piece sources.
    std::vector<fs::path> sources;
    size_t cursor = 0;
    while (true) {
        const size_t piece = xml.find("<Piece", cursor);
        if (piece == std::string::npos) break;
        cursor = piece + 6;
        const std::string source = tag_attribute(xml, piece, "Source");
        if (!source.empty()) sources.push_back(info.pvts.parent_path() / source);
    }
    if (sources.empty()) {
        error = info.pvts.filename().string() + ": no piece sources";
        return nullptr;
    }

    const size_t total_cells = static_cast<size_t>(grid->nx) * static_cast<size_t>(grid->ny);
    for (const auto& source : sources) {
        PieceData piece;
        if (!parse_vts_piece(source, piece, error)) return nullptr;
        const int i0 = piece.extent[0];
        const int local_nx = piece.extent[1] - piece.extent[0];
        const int local_ny = piece.extent[3] - piece.extent[2];
        if (i0 < 0 || i0 + local_nx > grid->nx || local_ny != grid->ny) {
            error = source.filename().string() + ": piece extent outside whole extent";
            return nullptr;
        }
        for (auto& [name, values] : piece.arrays) {
            auto& field = grid->fields[name];
            if (field.empty()) field.assign(total_cells, std::nan(""));
            for (int j = 0; j < local_ny; ++j) {
                for (int i = 0; i < local_nx; ++i) {
                    field[static_cast<size_t>(j) * grid->nx + (i0 + i)] =
                        values[static_cast<size_t>(j) * local_nx + i];
                }
            }
        }
    }

    // Cell-center coordinates come from the always-written helper fields.
    const auto* theta_rad = grid->field("theta_rad");
    const auto* z_m = grid->field("z_m");
    grid->theta_deg.resize(grid->nx);
    grid->z_m.resize(grid->ny);
    for (int i = 0; i < grid->nx; ++i) {
        grid->theta_deg[i] = theta_rad ? (*theta_rad)[i] * 180.0 / 3.14159265358979323846
                                       : (i + 0.5) * 360.0 / grid->nx;
    }
    for (int j = 0; j < grid->ny; ++j) {
        grid->z_m[j] = z_m ? (*z_m)[static_cast<size_t>(j) * grid->nx] : j + 0.5;
    }
    return grid;
}

ResultsLoader::~ResultsLoader() { join(); }

void ResultsLoader::join() {
    if (worker_.joinable()) worker_.join();
}

void ResultsLoader::scan(const fs::path& output_dir, const std::string& prefix) {
    if (busy_.exchange(true)) return;
    join();
    worker_ = std::thread([this, output_dir, prefix]() {
        std::string error;
        auto steps = parse_pvd(output_dir, prefix, error);
        {
            std::lock_guard lock(mutex_);
            steps_ = std::move(steps);
            error_ = error;
        }
        ++revision_;
        busy_ = false;
    });
}

void ResultsLoader::load_step(const StepInfo& info) {
    if (busy_.exchange(true)) return;
    join();
    worker_ = std::thread([this, info]() {
        std::string error;
        auto grid = parse_step_grid(info, error);
        {
            std::lock_guard lock(mutex_);
            if (grid) grid_ = grid;
            error_ = error;
        }
        ++revision_;
        busy_ = false;
    });
}

std::vector<StepInfo> ResultsLoader::steps() const {
    std::lock_guard lock(mutex_);
    return steps_;
}

std::shared_ptr<const ResultGrid> ResultsLoader::grid() const {
    std::lock_guard lock(mutex_);
    return grid_;
}

std::string ResultsLoader::error() const {
    std::lock_guard lock(mutex_);
    return error_;
}
