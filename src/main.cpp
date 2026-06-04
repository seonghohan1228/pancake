#include <petsc.h>

#include <algorithm>
#include <filesystem>

#include "communicator.hpp"
#include "config.hpp"
#include "energy.hpp"
#include "equation_of_state.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "io.hpp"
#include "journal_motion.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"
#include "utils.hpp"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {
namespace fs = std::filesystem;

fs::path executable_directory(const char* argv0) {
#ifdef _WIN32
    wchar_t path[MAX_PATH] = {};
    const DWORD length = GetModuleFileNameW(nullptr, path, MAX_PATH);
    if (length > 0 && length < MAX_PATH) {
        return fs::path(path).parent_path();
    }
#endif

    if (argv0 != nullptr && argv0[0] != '\0') {
        fs::path exe_path(argv0);
        if (exe_path.has_parent_path()) {
            return fs::absolute(exe_path).parent_path();
        }
    }
    return fs::current_path();
}

double global_min_field(const Field& field, const Mesh& mesh) {
    double local_min = 1e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local_min = std::min(local_min, field(i, j));
        }
    }
    double global_min = local_min;
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global_min;
}

void update_film_geometry(Fields& fields, const Mesh& mesh,
                          const SimulationConfig& cfg,
                          const JournalMotion::BearingState& bearing_state) {
    if (cfg.motion_model == MotionModel::MOVING_BEARING) {
        FilmThickness::compute_moving_bearing(fields["h"], fields["dh_dt"], mesh, cfg, bearing_state);
    } else {
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        fields["dh_dt"].fill(0.0);
    }
    JournalMotion::write_state_fields(fields, bearing_state, cfg);
}
}

int main(int argc, char** argv) {
    fs::path config_path = executable_directory(argc > 0 ? argv[0] : nullptr) / "config.txt";
    int petsc_argc = argc > 0 ? 1 : 0;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" || arg == "--config") {
            if (i + 1 < argc) {
                config_path = argv[i + 1];
                i++;
            }
        } else {
            argv[petsc_argc++] = argv[i];
        }
    }
    argc = petsc_argc;

    PetscInitialize(&argc, &argv, nullptr, nullptr);

    SimulationConfig cfg;

    cfg.load_from_file(config_path);
    Utils::log("Loaded configuration from: " + config_path.string());

    Mesh mesh(cfg);
    Communicator comm(mesh);
    IO::prepare_output_directory(cfg, MPI_COMM_WORLD);

    {
        // Initialize fields
        Fields fields;
        const double pressure_initial = cfg.initial_pressure();
        const double theta_initial =
            std::max(EOS::theta_from_pressure(pressure_initial, cfg.p_cav, cfg.bulk_modulus), cfg.theta_min);
        fields.add("pressure", mesh).fill(pressure_initial);
        fields.add("theta",    mesh).fill(theta_initial);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("dh_dt",    mesh).fill(0.0);
        fields.add("rho",      mesh).fill(cfg.rho);
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("heat_generation", mesh).fill(0.0);
        fields.add("inlet_indicator", mesh).fill(0.0);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z",     mesh, 2, GridLocation::FACE_Z);
        fields.add("pressure_force_x", mesh).fill(0.0);
        fields.add("pressure_force_y", mesh).fill(0.0);
        fields.add("pressure_force_z", mesh).fill(0.0);
        fields.add("load_x", mesh).fill(0.0);
        fields.add("load_y", mesh).fill(0.0);
        fields.add("load_z", mesh).fill(0.0);
        fields.add("viscous_force_x", mesh).fill(0.0);
        fields.add("viscous_force_y", mesh).fill(0.0);
        fields.add("viscous_force_z", mesh).fill(0.0);
        fields.add("fluid_force_x", mesh).fill(0.0);
        fields.add("fluid_force_y", mesh).fill(0.0);
        fields.add("fluid_force_z", mesh).fill(0.0);
        fields.add("external_load_x", mesh).fill(cfg.external_load_x);
        fields.add("external_load_y", mesh).fill(cfg.external_load_y);
        fields.add("external_load_z", mesh).fill(cfg.external_load_z);
        fields.add("bearing_x", mesh).fill(0.0);
        fields.add("bearing_y", mesh).fill(0.0);
        fields.add("bearing_z", mesh).fill(0.0);
        fields.add("bearing_attitude_angle", mesh).fill(cfg.attitude_angle_deg);
        fields.add("friction_torque", mesh).fill(0.0);

        Utils::log("Starting Simulation...");

        // Geometry initialization
        JournalMotion::BearingState bearing_state = JournalMotion::initial_state(cfg);
        update_film_geometry(fields, mesh, cfg, bearing_state);
        FilmThickness::compute_inlet_indicator(fields["inlet_indicator"], mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);

        // Initialize old_time to current state
        fields["pressure"].store_old_time();
        fields["theta"].store_old_time();
        fields["h"].store_old_time();
        fields["temperature"].store_old_time();

        double t = 0.0;
        double output_time = 0;
        int step = 0;
        const bool steady_state = cfg.solution_mode == SolutionMode::STEADY_STATE;

        // Time loop
        while (steady_state ? (step == 0) : (t <= cfg.end_t || std::abs(t - cfg.end_t) < 1e-8)) {
            update_film_geometry(fields, mesh, cfg, bearing_state);
            const double h_min = global_min_field(fields["h"], mesh);
            if (cfg.stop_on_nonpositive_film && h_min <= cfg.min_film_thickness) {
                Utils::log("Stopping: minimum film thickness reached " + std::to_string(h_min),
                           Utils::Color::RED, 1);
                break;
            }

            // Re-initialise to inlet pressure at each step (flooded guess).
            // This prevents the starvation trap where theta=0 persists.
            fields["theta"].fill(theta_initial);
            fields["pressure"].fill(pressure_initial);
            comm.update_ghosts(fields);

            if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS)
                Reynolds::solve_elrod(fields, sys, mesh, cfg);
            else
                Reynolds::solve(fields, sys, mesh, cfg);

            // Calculate velocities and macroscopic properties
            comm.update_ghosts(fields["pressure"]);
            comm.update_ghosts(fields["theta"]);
            Reynolds::calculate_velocities(fields, mesh, cfg);
            const Reynolds::ForceComponents forces =
                Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
            JournalMotion::write_state_fields(fields, bearing_state, cfg);

            comm.update_ghosts(fields["h"]);
            comm.update_ghosts(fields["temperature"]);
            Energy::solve(fields, sys, mesh, cfg);

            // Output
            if (t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0) {
                comm.update_ghosts(fields["velocity_theta"]);
                comm.update_ghosts(fields["velocity_z"]);
                IO::write_timestep(t, step, mesh, fields, cfg);
                Utils::log("Step " + std::to_string(step) + " t=" + std::to_string(t),
                           Utils::Color::GREEN, 1);
                output_time += cfg.write_interval;
            }

            // Advance time
            fields["pressure"].store_old_time();
            fields["theta"].store_old_time();
            fields["h"].store_old_time();
            fields["temperature"].store_old_time();
            if (!steady_state) {
                JournalMotion::advance(
                    bearing_state,
                    {forces.fluid_x, forces.fluid_y, forces.fluid_z},
                    cfg,
                    cfg.dt);
                t += cfg.dt;
            }
            step++;
        }

        Utils::log("Simulation finished successfully.");
    }  // LinearSystem and Fields destroyed before PetscFinalize

    PetscFinalize();
    return 0;
}
