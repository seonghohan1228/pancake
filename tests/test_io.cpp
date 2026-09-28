#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include <mpi.h>

#include "config.hpp"
#include "field.hpp"
#include "io.hpp"
#include "mesh.hpp"

namespace fs = std::filesystem;

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static std::string read_text(const fs::path& path) {
    std::ifstream file(path);
    check(file.is_open(), "expected VTK output file was not written");
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
}

static bool has_array(const std::string& text, const std::string& name) {
    return text.find("Name=\"" + name + "\"") != std::string::npos;
}

static bool has_vector_array(const std::string& text, const std::string& name) {
    const std::string token = "Name=\"" + name + "\"";
    const size_t pos = text.find(token);
    if (pos == std::string::npos) return false;
    const size_t tag_end = text.find('>', pos);
    if (tag_end == std::string::npos) return false;
    const std::string tag = text.substr(pos, tag_end - pos);
    return tag.find("NumberOfComponents=\"3\"") != std::string::npos;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    SimulationConfig cfg;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 3;
    cfg.output_dir = "test_io_results";
    cfg.filename_prefix = "io";
    cfg.output_write_3d = true;
    cfg.output_write_flat = true;
    cfg.output_fields = {"velocity"};

    Mesh mesh(cfg);
    Fields fields;
    fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA).fill(2.0);
    fields.add("velocity_z", mesh, 2, GridLocation::FACE_Z).fill(3.0);

    IO::prepare_output_directory(cfg);
    if (mesh.rank == 0) {
        std::ofstream(fs::path(cfg.output_dir) / "experiment_notes.txt") << "preserve me";
        std::ofstream(fs::path(cfg.output_dir) / "io_999.pvts") << "old solver output";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    IO::prepare_output_directory(cfg);
    if (mesh.rank == 0) {
        check(fs::exists(fs::path(cfg.output_dir) / "experiment_notes.txt"),
              "rerunning must preserve unrelated output-directory files");
        check(!fs::exists(fs::path(cfg.output_dir) / "io_999.pvts"),
              "rerunning must remove stale solver products");
    }
    IO::write_timestep(0.0, 0, mesh, fields, cfg);
    MPI_Barrier(MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        const std::string curved_vts = read_text(fs::path(cfg.output_dir) / "processor0" / "io_0_0.vts");
        const std::string curved_pvts = read_text(fs::path(cfg.output_dir) / "io_0.pvts");
        const std::string flat_vts = read_text(fs::path(cfg.output_dir) / "flat" / "processor0" / "io_0_0.vts");
        const std::string flat_pvts = read_text(fs::path(cfg.output_dir) / "flat" / "io_0.pvts");

        check(curved_vts.find("<CellData Vectors=\"U\">") != std::string::npos,
              "curved VTS should mark U as active cell vector");
        check(curved_pvts.find("<PCellData Vectors=\"U\">") != std::string::npos,
              "curved PVTS should mark U as active cell vector");
        check(flat_vts.find("<CellData Vectors=\"U\">") != std::string::npos,
              "flat VTS should mark U as active cell vector");
        check(flat_pvts.find("<PCellData Vectors=\"U\">") != std::string::npos,
              "flat PVTS should mark U as active cell vector");

        check(has_vector_array(curved_vts, "U"), "curved VTS U should be a 3-component vector");
        check(has_vector_array(curved_pvts, "U"), "curved PVTS U should be a 3-component vector");
        check(has_vector_array(flat_vts, "U"), "flat VTS U should be a 3-component vector");
        check(has_vector_array(flat_pvts, "U"), "flat PVTS U should be a 3-component vector");

        // Components are unified into the Cartesian vector; no per-component scalar arrays.
        check(!has_array(curved_vts, "velocity_x") && !has_array(curved_vts, "velocity_y") &&
              !has_array(curved_vts, "velocity_z") && !has_array(curved_vts, "velocity_theta") &&
              !has_array(curved_vts, "velocity"),
              "curved VTS should not emit velocity component scalars or legacy 'velocity'");
        check(!has_array(flat_vts, "velocity_theta") && !has_array(flat_vts, "velocity_z") &&
              !has_array(flat_vts, "velocity"),
              "flat VTS should not emit velocity component scalars or legacy 'velocity'");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // A warm start is a global field ordering, independent of the MPI theta split.
    Fields initial;
    const char* names[] = {"pressure", "theta", "temperature", "cross_film_collocation"};
    for (const char* name : names) initial.add(name, mesh);
    auto expected = [](int k, int index) {
        const double starts[] = {200000.0, 0.5, 300.0, -20.0};
        const double steps[] = {1.0, 0.001, 0.01, 0.1};
        return starts[k] + steps[k] * index;
    };
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j)
            initial[names[k]](i, j) = expected(k, (mesh.offset_theta + i) * mesh.n_z_global + j);
    const fs::path guess = fs::path(cfg.output_dir) / "steady_state.dat";
    IO::write_steady_initial_state(guess, initial, mesh, cfg);
    for (const char* name : names) initial[name].fill(0.0);
    IO::read_steady_initial_state(guess, initial, mesh, cfg);
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j)
            check(initial[names[k]](i, j) == expected(k, (mesh.offset_theta + i) * mesh.n_z_global + j),
                  "steady initial guess must round-trip every physical cell exactly");
    initial.add("cross_film_layer_0", mesh).fill(0.0);
    initial.add("cross_film_layer_1", mesh).fill(0.0);
    IO::read_steady_initial_state(guess, initial, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j) {
        check(initial["cross_film_layer_0"](i, j) == initial["temperature"](i, j),
              "legacy mean-only guesses must initialize every resolved layer");
        initial["cross_film_layer_1"](i, j) = initial["temperature"](i, j) + 2.0;
    }
    IO::write_steady_initial_state(guess, initial, mesh, cfg);
    initial["cross_film_layer_0"].fill(0.0);
    initial["cross_film_layer_1"].fill(0.0);
    IO::read_steady_initial_state(guess, initial, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j) {
        check(initial["cross_film_layer_0"](i, j) == initial["temperature"](i, j) &&
              initial["cross_film_layer_1"](i, j) == initial["temperature"](i, j) + 2.0,
              "resolved-layer guesses must round-trip every layer without overwriting the mean");
    }
    SimulationConfig mismatch = cfg;
    mismatch.c *= 2.0;
    bool rejected = false;
    try { IO::read_steady_initial_state(guess, initial, mesh, mismatch); }
    catch (const std::exception&) { rejected = true; }
    check(rejected, "steady initial guess must reject incompatible geometry on every rank");
    if (mesh.rank == 0) std::ofstream(guess) << "PANCAKE_STEADY_GUESS 1\n8 3\n";
    MPI_Barrier(MPI_COMM_WORLD);
    rejected = false;
    try { IO::read_steady_initial_state(guess, initial, mesh, cfg); }
    catch (const std::exception&) { rejected = true; }
    check(rejected, "steady initial guess must reject truncated input collectively");

    SimulationConfig legacy_cfg;
    legacy_cfg.n_theta_global = 8;
    legacy_cfg.n_z_global = 3;
    legacy_cfg.output_dir = "test_io_legacy_results";
    legacy_cfg.filename_prefix = "legacy";
    legacy_cfg.output_write_3d = true;
    legacy_cfg.output_write_flat = false;
    legacy_cfg.output_fields = {"load_x"};

    Mesh legacy_mesh(legacy_cfg);
    Fields legacy_fields;
    legacy_fields.add("pressure_force_x", legacy_mesh).fill(42.0);
    legacy_fields.add("pressure_force_y", legacy_mesh).fill(7.0);
    legacy_fields.add("pressure_force_z", legacy_mesh).fill(0.0);

    IO::prepare_output_directory(legacy_cfg);
    IO::write_timestep(0.0, 0, legacy_mesh, legacy_fields, legacy_cfg);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        const std::string legacy_vts =
            read_text(fs::path(legacy_cfg.output_dir) / "processor0" / "legacy_0_0.vts");
        // load_x selects the pressure force, now written as the Cartesian vector Fp.
        check(has_vector_array(legacy_vts, "Fp"),
              "legacy load_x request should write the pressure-force vector Fp");
        check(!has_array(legacy_vts, "pressure_force_x") && !has_array(legacy_vts, "load_x"),
              "pressure-force components should be unified into Fp, not written as scalars");

        fs::remove_all(cfg.output_dir);
        fs::remove_all(legacy_cfg.output_dir);
        std::cout << "PASS: test_io\n";
    }

    MPI_Finalize();
    return 0;
}
