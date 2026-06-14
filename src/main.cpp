#include <petsc.h>

#include <algorithm>
#include <exception>
#include <filesystem>
#include <sstream>
#include <vector>

#include "communicator.hpp"
#include "config.hpp"
#include "diagnostics.hpp"
#include "energy.hpp"
#include "equation_of_state.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "fluid_properties.hpp"
#include "gas_transport.hpp"
#include "io.hpp"
#include "journal_motion.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "property_table_io.hpp"
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

void global_min_max_field(const Field& field, const Mesh& mesh, double& out_min, double& out_max) {
    double local_min = 1e300, local_max = -1e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            local_min = std::min(local_min, field(i, j));
            local_max = std::max(local_max, field(i, j));
        }
    }
    MPI_Allreduce(&local_min, &out_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &out_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

// Copy physical (non-ghost) cells of src into dst, for outer-iteration snapshots.
void copy_physical(const Field& src, Field& dst, const Mesh& mesh) {
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            dst(i, j) = src(i, j);
}

// Global max |cur - prev| over physical cells, normalized by ref (pressure space).
double normalized_max_change(const Field& cur, const Field& prev, const Mesh& mesh, double ref) {
    double local_max = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            local_max = std::max(local_max, std::abs(cur(i, j) - prev(i, j)));
    double global_max = local_max;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_max / std::max(std::abs(ref), 1.0);
}

// Under-relax field toward its pre-sweep snapshot: field = omega*field + (1-omega)*prev.
void relax_field(Field& field, const Field& prev, double omega, const Mesh& mesh) {
    if (omega >= 1.0) return;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            field(i, j) = omega * field(i, j) + (1.0 - omega) * prev(i, j);
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
    for (const auto& warning : cfg.parse_warnings) {
        Utils::log("Config warning: " + warning, Utils::Color::YELLOW, 1);
    }

    std::vector<std::string> config_errors;
    std::vector<std::string> config_warnings;
    // Load external 2-D PTSV property tables before validation so the saturation
    // context sees them; validate() clears config_errors, so merge afterwards.
    std::vector<std::string> table_errors;
    FluidProperties::load_config_property_tables(cfg, table_errors);
    cfg.validate(config_errors, config_warnings);
    config_errors.insert(config_errors.end(), table_errors.begin(), table_errors.end());
    for (const auto& warning : config_warnings) {
        Utils::log("Config warning: " + warning, Utils::Color::YELLOW, 1);
    }
    if (!config_errors.empty()) {
        Utils::log("Configuration failed validation:", Utils::Color::RED, 1);
        for (const auto& error : config_errors) {
            Utils::log("- " + error, Utils::Color::RED, 2);
        }
        PetscFinalize();
        return 1;
    }

    Mesh mesh(cfg);
    Communicator comm(mesh);
    try {
        IO::prepare_output_directory(cfg, MPI_COMM_WORLD);
    } catch (const std::exception& ex) {
        Utils::log(std::string("Output directory setup failed: ") + ex.what(),
                   Utils::Color::RED, 1);
        PetscFinalize();
        return 3;
    }

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
        fields.add("mu",       mesh).fill(cfg.mu);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("mu_liquid_solution", mesh).fill(cfg.mu);
        fields.add("cp_liquid_solution", mesh).fill(cfg.rho_cp / cfg.rho);
        fields.add("rho_cp_liquid_solution", mesh).fill(cfg.rho_cp);
        fields.add("k_liquid_solution", mesh).fill(cfg.thermal_conductivity);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("gas_mass_transfer", mesh).fill(0.0);
        fields.add("rho_gas", mesh).fill(0.0);
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("heat_generation", mesh).fill(0.0);
        fields.add("inlet_indicator", mesh).fill(0.0);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA);
        fields.add("velocity_z",     mesh, 2, GridLocation::FACE_Z);
        fields.add("pressure_force_x", mesh).fill(0.0);
        fields.add("pressure_force_y", mesh).fill(0.0);
        fields.add("pressure_force_z", mesh).fill(0.0);
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
        fields.add("cavitation_clamp_mass", mesh).fill(0.0);
        fields.add("gas_clamp_mass", mesh).fill(0.0);

        Utils::log("Starting Simulation...");

        // Geometry initialization
        JournalMotion::BearingState bearing_state = JournalMotion::initial_state(cfg);
        update_film_geometry(fields, mesh, cfg, bearing_state);
        FilmThickness::compute_inlet_indicator(fields["inlet_indicator"], mesh, cfg);
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);

        LinearSystem sys(mesh);
        Field p_outer_prev("p_outer_prev", mesh);  // pre-sweep snapshots for the
        Field T_outer_prev("T_outer_prev", mesh);  // outer (Picard) coupling residual

        // Initialize old_time to current state
        fields["pressure"].store_old_time();
        fields["theta"].store_old_time();
        fields["h"].store_old_time();
        fields["temperature"].store_old_time();
        fields["dissolved_gas"].store_old_time();
        fields["free_gas_mass"].store_old_time();

        double t = 0.0;
        double output_time = 0;
        int step = 0;
        const bool steady_state = cfg.solution_mode == SolutionMode::STEADY_STATE;

        // Time loop
        while (steady_state ? (step == 0) : (t <= cfg.end_t || std::abs(t - cfg.end_t) < 1e-8)) {
            SimulationConfig step_cfg = cfg;
            step_cfg.omega = cfg.omega_at_time(t);

            update_film_geometry(fields, mesh, cfg, bearing_state);
            const double h_min = global_min_field(fields["h"], mesh);
            if (cfg.stop_on_nonpositive_film && h_min <= cfg.min_film_thickness) {
                Utils::log("Stopping: minimum film thickness reached " + std::to_string(h_min),
                           Utils::Color::RED, 1);
                break;
            }

            // Outer (Picard) coupling: in STEADY_STATE repeat the segregated sweep
            // until the pressure/temperature residuals converge (a real fixed
            // point, reported in pressure space). TRANSIENT runs exactly one
            // sweep, identical to before.
            Reynolds::ElrodStats elrod_stats;
            Reynolds::ForceComponents forces;
            const double gas_dt = steady_state ? 0.0 : cfg.dt;
            const int max_coupling = steady_state ? std::max(1, cfg.coupling_max_iters) : 1;
            const double p_ref = std::max(cfg.property_reference_pressure, 1.0);
            const double T_ref = std::max(cfg.temperature_reference, 1.0);
            int coupling_iters = 0;
            double coupling_R_p = 0.0;
            double coupling_R_T = 0.0;

            for (int outer = 0; outer < max_coupling; ++outer) {
                if (steady_state) {
                    copy_physical(fields["pressure"], p_outer_prev, mesh);
                    copy_physical(fields["temperature"], T_outer_prev, mesh);
                }

                FluidProperties::update_solution_fields(fields, mesh, cfg);
                comm.update_ghosts(fields);

                if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS)
                    elrod_stats = Reynolds::solve_elrod(fields, sys, mesh, step_cfg);
                else
                    Reynolds::solve(fields, sys, mesh, step_cfg);

                FluidProperties::update_solution_fields(fields, mesh, cfg);

                // Sync the flow + property fields the gas transport and velocity solve read.
                comm.update_ghosts(fields["pressure"]);
                comm.update_ghosts(fields["theta"]);
                comm.update_ghosts(fields["rho"]);
                comm.update_ghosts(fields["mu"]);
                comm.update_ghosts(fields["rho_liquid_solution"]);
                comm.update_ghosts(fields["mu_liquid_solution"]);
                comm.update_ghosts(fields["rho_cp_liquid_solution"]);
                comm.update_ghosts(fields["k_liquid_solution"]);

                // Gaseous cavitation (operator splitting): transport dissolved gas with the
                // liquid film mass flux, then exchange dissolved <-> free gas locally, then
                // refresh mixture properties (rho, mu) for the velocity/energy solves. In
                // STEADY_STATE with steady_gas_model = EQUILIBRIUM, flash to the local
                // saturation state instead (transport is meaningless without time advance).
                if (steady_state && cfg.steady_gas_model == SteadyGasModel::EQUILIBRIUM) {
                    FluidProperties::equilibrate_gas(fields, mesh, cfg);
                } else {
                    comm.update_ghosts(fields["dissolved_gas"]);
                    comm.update_ghosts(fields["free_gas_mass"]);
                    GasTransport::solve(fields, sys, mesh, step_cfg, gas_dt);
                    FluidProperties::update_gas_state(fields, mesh, cfg, gas_dt);
                }
                FluidProperties::update_solution_fields(fields, mesh, cfg);
                comm.update_ghosts(fields["rho"]);
                comm.update_ghosts(fields["mu"]);

                // Calculate velocities and macroscopic properties
                Reynolds::calculate_velocities(fields, mesh, step_cfg);
                forces = Reynolds::calculate_macroscopic_properties(fields, mesh, step_cfg);
                JournalMotion::write_state_fields(fields, bearing_state, cfg);

                comm.update_ghosts(fields["h"]);
                comm.update_ghosts(fields["temperature"]);
                Energy::solve(fields, sys, mesh, step_cfg);

                coupling_iters = outer + 1;
                if (steady_state) {
                    coupling_R_p = normalized_max_change(fields["pressure"], p_outer_prev, mesh, p_ref);
                    coupling_R_T = normalized_max_change(fields["temperature"], T_outer_prev, mesh, T_ref);
                    relax_field(fields["pressure"], p_outer_prev, cfg.coupling_relaxation, mesh);
                    relax_field(fields["temperature"], T_outer_prev, cfg.coupling_relaxation, mesh);
                    if (coupling_R_p < cfg.coupling_tolerance && coupling_R_T < cfg.coupling_tolerance)
                        break;
                }
            }
            if (steady_state && cfg.log_outer_iters) {
                std::ostringstream cl;
                cl << "Coupling: " << coupling_iters << " outer iters, R_p=" << coupling_R_p
                   << " R_T=" << coupling_R_T;
                Utils::log(cl.str(), Utils::Color::CYAN, 1);
            }

            // Conservation/convergence diagnostics. Must run before
            // store_old_time so the step's storage terms are still formed
            // against the previous timestep.
            const bool diagnostics_now = cfg.diagnostics_interval > 0 &&
                step % cfg.diagnostics_interval == 0;
            Diagnostics::MassBalance balance;
            if (diagnostics_now) {
                balance = Diagnostics::mass_balance(fields, mesh, step_cfg,
                                                    steady_state ? 0.0 : cfg.dt);
                Diagnostics::append_csv(
                    std::filesystem::path(cfg.output_dir) / "diagnostics.csv",
                    t, step, balance, elrod_stats);
            }

            // Output: write VTK at the configured save cadence ...
            const bool save_now =
                (t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0);
            if (save_now) {
                comm.update_ghosts(fields["velocity_theta"]);
                comm.update_ghosts(fields["velocity_z"]);
                IO::write_timestep(t, step, mesh, fields, cfg);
                output_time += cfg.write_interval;
            }

            // ... but log every timestep. SIMPLIFIED is a compact line; DETAILED adds
            // per-step diagnostics. "[saved]" marks the steps written to disk.
            {
                std::ostringstream line;
                line << "Step " << step << " t=" << std::to_string(t);
                if (cfg.output_verbosity == OutputVerbosity::DETAILED) {
                    double t_min = 0.0, t_max = 0.0;
                    global_min_max_field(fields["temperature"], mesh, t_min, t_max);
                    line << " | h_min=" << h_min
                         << " | T=[" << t_min << ", " << t_max << "] K"
                         << " | F=(" << forces.fluid_x << ", " << forces.fluid_y << ", "
                         << forces.fluid_z << ") N"
                         << " | M=" << forces.friction_torque << " Nm";
                    if (cfg.cavitation_model == CavitationModel::ELROD_ADAMS) {
                        line << " | outer=" << elrod_stats.outer_iters
                             << " flips=" << elrod_stats.flag_flips;
                    }
                    if (diagnostics_now) {
                        line << " | r_l=" << balance.liquid_residual;
                        if (cfg.fluid_property_model == FluidPropertyModel::GAS_CAVITATION_MIXTURE) {
                            line << " r_gas=" << balance.gas_residual;
                        }
                    }
                }
                if (save_now) line << "  [saved]";
                Utils::log(line.str(), Utils::Color::GREEN, 1);
            }

            // Advance time
            fields["pressure"].store_old_time();
            fields["theta"].store_old_time();
            fields["h"].store_old_time();
            fields["temperature"].store_old_time();
            fields["dissolved_gas"].store_old_time();
            fields["free_gas_mass"].store_old_time();
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
