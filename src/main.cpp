#include <petsc.h>

#include <algorithm>
#include <exception>
#include <filesystem>
#include <limits>
#include <sstream>
#include <vector>

#include "communicator.hpp"
#include "config.hpp"
#include "coupling.hpp"
#include "coupled_cavitation.hpp"
#include "diagnostics.hpp"
#include "cross_film.hpp"
#include "energy.hpp"
#include "equation_of_state.hpp"
#include "fbns.hpp"
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
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double change = std::abs(cur(i, j) - prev(i, j));
            local_max = !std::isfinite(change) ? std::numeric_limits<double>::infinity()
                                              : std::max(local_max, change);
        }
    }
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

    const PetscErrorCode initialization_error=PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (initialization_error) return initialization_error;

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

    bool nonlinear_failed = false;
    {
        // Initialize fields
        Fields fields;
        const double pressure_initial = cfg.initial_pressure();
        // FBNS carries no bulk modulus, so its film content starts at the
        // full-film value; the solver detects that as a cold start and builds its
        // own uncavitated initial guess.
        const bool coupled = cfg.cavitation_model == CavitationModel::COUPLED;
        const double theta_initial = (FBNS::active(cfg) || coupled)
            ? 1.0
            : std::max(EOS::theta_from_pressure(pressure_initial, cfg.effective_p_cav(), cfg.bulk_modulus), cfg.theta_min);
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
        // Lagged, under-relaxed cavitated-zone plateau. Only VOID_COUPLED needs it;
        // its presence is what switches Reynolds from recomputing p_void every
        // assembly (which limit-cycles) to relaxing it once per outer iteration.
        if (cfg.gas_pressure_coupling == GasPressureCoupling::VOID_COUPLED) {
            fields.add("p_void", mesh).fill(cfg.effective_p_cav());
        }
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("heat_generation", mesh).fill(0.0);
        CrossFilm::register_fields(fields, mesh, cfg);
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
        // Per-cell continuity residual [kg/s]. Registered only when selected, so
        // an unselected run allocates nothing and takes the original code path.
        if (cfg.output_field_enabled("mass_residual")) {
            fields.add("mass_residual", mesh).fill(0.0);
        }

        Utils::log("Starting Simulation...");

        // Geometry initialization
        JournalMotion::BearingState bearing_state = JournalMotion::initial_state(cfg);
        update_film_geometry(fields, mesh, cfg, bearing_state);
        FilmThickness::compute_inlet_indicator(fields["inlet_indicator"], mesh, cfg);
        if (!cfg.steady_initial_state_file.empty()) {
            try {
                IO::read_steady_initial_state(cfg.steady_initial_state_file, fields, mesh, cfg);
            } catch (const std::exception& ex) {
                Utils::log(ex.what(), Utils::Color::RED, 1);
                MPI_Abort(MPI_COMM_WORLD, 3);
            }
        }
        if (!coupled) FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);

        std::unique_ptr<CoupledCavitation::Solver> coupled_solver;
        if (coupled) {
            coupled_solver = std::make_unique<CoupledCavitation::Solver>(mesh,comm,cfg.temperature_model==TemperatureModel::ENERGY_EQUATION);
            coupled_solver->initialize(fields,cfg);
        }

        LinearSystem sys(mesh, linear_solver_options_from_config(cfg));
        Field p_outer_prev("p_outer_prev", mesh);  // pre-sweep snapshots for the
        Field T_outer_prev("T_outer_prev", mesh);  // outer (Picard) coupling residual
        Field c_outer_prev("c_outer_prev", mesh);
        Field alpha_outer_prev("alpha_outer_prev", mesh);

        // Initialize old_time to current state
        fields["pressure"].store_old_time();
        fields["theta"].store_old_time();
        fields["h"].store_old_time();
        fields["rho_liquid_solution"].store_old_time();
        fields["temperature"].store_old_time();
        fields["dissolved_gas"].store_old_time();
        fields["free_gas_mass"].store_old_time();

        double t = 0.0;
        double output_time = 0;
        int step = 0;
        const bool steady_state = cfg.solution_mode == SolutionMode::STEADY_STATE;
        Reynolds::ForceComponents previous_forces;
        if (!steady_state) {
            SimulationConfig initial_cfg = cfg;
            initial_cfg.omega = cfg.omega_at_time(0);
            Reynolds::calculate_velocities(fields, mesh, initial_cfg);
            previous_forces = Reynolds::calculate_macroscopic_properties(fields, mesh, initial_cfg);
            JournalMotion::write_state_fields(fields, bearing_state, cfg);
            comm.update_ghosts(fields["velocity_theta"]);
            comm.update_ghosts(fields["velocity_z"]);
            IO::write_timestep(0, 0, mesh, fields, cfg);
            output_time = cfg.write_interval;
        }

        // Time loop
        while (steady_state ? (step == 0) : (t < cfg.end_t)) {
            SimulationConfig step_cfg = cfg;
            if (!steady_state) {
                const double remaining = cfg.end_t - t;
                const double time_roundoff = 1e-12 * std::max(cfg.dt, cfg.end_t);
                step_cfg.dt = remaining <= cfg.dt + time_roundoff ? remaining : cfg.dt;
                t = step_cfg.dt == remaining ? cfg.end_t : t + step_cfg.dt;
                ++step;
                // Explicit fluid-force predictor: geometry and the field solution
                // now share the target time. The force stays lagged by one step.
                JournalMotion::advance(bearing_state,
                    {previous_forces.fluid_x, previous_forces.fluid_y, previous_forces.fluid_z},
                    cfg, step_cfg.dt);
            }
            step_cfg.omega = cfg.omega_at_time(t);

            update_film_geometry(fields, mesh, cfg, bearing_state);
            const double h_min = global_min_field(fields["h"], mesh);
            if (cfg.stop_on_nonpositive_film && h_min <= cfg.min_film_thickness) {
                nonlinear_failed = true;
                Utils::log("Stopping: minimum film thickness reached " + std::to_string(h_min),
                           Utils::Color::RED, 1);
                break;
            }

            // Outer (Picard) coupling: in STEADY_STATE repeat the segregated sweep
            // until the pressure/temperature residuals converge (a real fixed
            // point, reported in pressure space). TRANSIENT runs exactly one
            // sweep, identical to before.
            Reynolds::ElrodStats elrod_stats;
            CoupledCavitation::Stats coupled_stats;
            Reynolds::ForceComponents forces;
            const double gas_dt = steady_state ? 0.0 : step_cfg.dt;
            const int max_coupling = steady_state ? std::max(1, cfg.coupling_max_iters) : 1;
            const double p_ref = std::max(cfg.property_reference_pressure, 1.0);
            const double T_ref = std::max(cfg.temperature_reference, 1.0);
            int coupling_iters = 0;
            double coupling_R_p = 0.0;
            double coupling_R_T = 0.0;
            CrossFilm::WallState cross_film_state;

            // Outer-loop control: convergence test, relaxation schedule and
            // best-iterate bookkeeping. The snapshot below is what makes a
            // cap-exhausted run deliver its lowest-residual iterate instead of
            // whatever the limit cycle happened to land on.
            Coupling::Controller coupling_control(cfg);
            Coupling::Snapshot coupling_best;
            Reynolds::ElrodStats best_elrod_stats;
            Reynolds::ForceComponents best_forces;
            CrossFilm::WallState best_cross_film_state;
            const bool track_best =
                steady_state && cfg.coupling_keep_best_iterate && max_coupling > 1;

            for (int outer = 0; outer < max_coupling; ++outer) {
                if (coupled) {
                    try {
                        coupled_stats = coupled_solver->solve(fields,step_cfg,step_cfg.dt,cfg.coupled_phase_equilibrium);
                    } catch (const std::exception& ex) {
                        nonlinear_failed = true;
                        Utils::log(ex.what(),Utils::Color::RED,1);
                        break;
                    }
                    elrod_stats.converged = coupled_stats.converged;
                    elrod_stats.outer_iters = coupled_stats.iterations;
                    if (!coupled_stats.converged) {
                        nonlinear_failed = true;
                        Utils::log("Coupled conservation/volume solve failed; residual = " + std::to_string(coupled_stats.residual),Utils::Color::RED,1);
                        break;
                    }
                    forces = Reynolds::calculate_macroscopic_properties(fields,mesh,step_cfg);
                    JournalMotion::write_state_fields(fields,bearing_state,cfg);
                    coupling_iters = 1;
                    break;
                }
                if (steady_state) {
                    copy_physical(fields["pressure"], p_outer_prev, mesh);
                    copy_physical(fields["temperature"], T_outer_prev, mesh);
                    copy_physical(fields["dissolved_gas"], c_outer_prev, mesh);
                    copy_physical(fields["alpha_gas"], alpha_outer_prev, mesh);
                }

                FluidProperties::update_solution_fields(fields, mesh, cfg);
                comm.update_ghosts(fields);

                // Cross-film profile, conjugate walls and the Dowson F-integrals.
                // Runs before Reynolds so the generalized coefficients it produces
                // are the ones the pressure solve uses this iteration, and the wall
                // temperatures it solves are the ones the energy solve closes on.
                // No-op unless cross_film_thermal_model is enabled.
                if (CrossFilm::active(cfg)) {
                    try {
                        cross_film_state = CrossFilm::update(fields, mesh, cfg);
                    } catch (const std::exception& ex) {
                        nonlinear_failed = true;
                        Utils::log(ex.what(), Utils::Color::RED, 1);
                        break;
                    }
                    // The Reynolds assembly reads these on both theta ghost layers.
                    comm.update_ghosts(fields["cross_film_f0"]);
                    comm.update_ghosts(fields["cross_film_f1"]);
                    comm.update_ghosts(fields["cross_film_f2"]);
                    comm.update_ghosts(fields["cross_film_couette_h"]);
                }

                if (FBNS::active(cfg)) {
                    // The FBNS Newton stats are reported through the same
                    // diagnostics columns: outer_iters is the Newton count and
                    // flag_flips stays 0 (there is no switch function to flip).
                    const FBNS::Stats fbns_stats =
                        FBNS::solve(fields, sys, comm, mesh, step_cfg);
                    elrod_stats.outer_iters = fbns_stats.newton_iters;
                    elrod_stats.flag_flips = 0;
                    elrod_stats.max_theta_change = fbns_stats.max_theta_change;
                    elrod_stats.converged = fbns_stats.converged;
                    if (!fbns_stats.converged) {
                        std::ostringstream msg;
                        msg << "FBNS did not converge in " << fbns_stats.newton_iters
                            << " Newton iterations: |G|=" << fbns_stats.residual_mass
                            << " |FB|=" << fbns_stats.residual_complementarity;
                        Utils::log(msg.str(), Utils::Color::YELLOW, 1);
                    }
                } else if (cfg.cavitation_model == CavitationModel::JFO) {
                    elrod_stats = Reynolds::solve_elrod(fields, sys, mesh, step_cfg);
                } else {
                    Reynolds::solve(fields, sys, mesh, step_cfg);
                }
                if (!steady_state && cfg.mass_conserving_cavitation() && !elrod_stats.converged)
                    nonlinear_failed = true;

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
                // feed flash instead; conservative steady species transport is still absent.
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
                // The free-gas fraction sets the full-film ceiling (FullFilm::bound)
                // under TWO_PHASE, and the energy, velocity and diagnostics passes
                // below evaluate that ceiling on theta ghost cells. Without this the
                // ceiling would be one gas update stale across every rank boundary.
                if (cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE) {
                    comm.update_ghosts(fields["alpha_gas"]);
                    comm.update_ghosts(fields["theta"]);
                    comm.update_ghosts(fields["rho_liquid_solution"]);
                    comm.update_ghosts(fields["mu_liquid_solution"]);
                    comm.update_ghosts(fields["rho_cp_liquid_solution"]);
                    comm.update_ghosts(fields["k_liquid_solution"]);
                }

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
                    const double coupling_R_g = cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE
                        ? std::max(normalized_max_change(fields["dissolved_gas"], c_outer_prev, mesh, 1.0) /
                                       std::max(cfg.dissolved_gas_initial, 1e-6),
                                   normalized_max_change(fields["alpha_gas"], alpha_outer_prev, mesh, 1.0)) : 0.0;
                    bool energy_balanced = !CrossFilm::active(cfg) ||
                        (std::isfinite(cross_film_state.profile_residual) &&
                         cross_film_state.profile_residual < cfg.coupling_tolerance);
                    double coupling_R_E = std::numeric_limits<double>::quiet_NaN();
                    if (cfg.groove_mixing_model == GrooveMixingModel::FLUX_BALANCED) {
                        if (energy_balanced &&
                            std::max(coupling_R_p, coupling_R_T) < cfg.coupling_tolerance) {
                            comm.update_ghosts(fields["temperature"]);
                            const auto ledger = Energy::balance(fields, mesh, step_cfg);
                            coupling_R_E = std::abs(ledger.residual_W) / std::max(ledger.scale_W, 1.0);
                            energy_balanced = std::isfinite(ledger.residual_W) &&
                                std::isfinite(ledger.scale_W) && ledger.scale_W >= 0.0 &&
                                coupling_R_E < cfg.coupling_tolerance;
                        } else {
                            energy_balanced = false;
                        }
                    }
                    const bool tolerance_met = coupling_control.record(
                        outer, coupling_R_p, coupling_R_T, coupling_R_g, energy_balanced);
                    if (cfg.log_coupling_history) {
                        std::ostringstream h;
                        h << "Coupling outer=" << (outer + 1) << " R_p=" << coupling_R_p
                          << " R_T=" << coupling_R_T
                          << " R_g=" << coupling_R_g
                          << " R_X=" << cross_film_state.profile_residual
                          << " R_E=" << coupling_R_E
                          << " omega=" << coupling_control.relaxation()
                          << (coupling_control.is_best() ? " *best" : "");
                        Utils::log(h.str(), Utils::Color::BLUE, 2);
                    }
                    // Snapshot BEFORE relaxation: the unrelaxed end-of-sweep state is
                    // the one the solvers produced and the one the residual just
                    // recorded describes.
                    if (track_best && coupling_control.is_best()) {
                        coupling_best.capture(fields);
                        best_elrod_stats = elrod_stats;
                        best_forces = forces;
                        best_cross_film_state = cross_film_state;
                    }
                    if (tolerance_met) break;
                    relax_field(fields["pressure"], p_outer_prev, coupling_control.relaxation(), mesh);
                    relax_field(fields["temperature"], T_outer_prev, coupling_control.relaxation(), mesh);
                }
            }
            if (steady_state && !coupled) {
                coupling_control.finish(coupling_iters);
                if (coupling_control.status().restored_best && coupling_best.valid()) {
                    // Ghost layers are part of the snapshot, so the restored state
                    // is bit-identical to the captured sweep and needs no exchange.
                    coupling_best.restore(fields);
                    elrod_stats = best_elrod_stats;
                    forces = best_forces;
                    cross_film_state = best_cross_film_state;
                }
            }
            // The non-convergence verdict is not optional output: a run that ends
            // on the cap without meeting coupling_tolerance must say so whether or
            // not log_outer_iters is on.
            if (steady_state && !coupled && !coupling_control.status().converged) {
                nonlinear_failed = true;
                Utils::log(coupling_control.summary_line(), Utils::Color::RED, 1);
                Utils::log(coupling_control.rollback_line(), Utils::Color::RED, 1);
            }
            if (cfg.mass_conserving_cavitation() && !elrod_stats.converged) nonlinear_failed = true;
            // A rejected coupled solve leaves inventories at the last accepted
            // time. Do not label that retained state with the failed target time.
            if (coupled && nonlinear_failed) break;
            if (steady_state && cfg.log_outer_iters) {
                if (coupling_control.status().converged) {
                    Utils::log(coupling_control.summary_line(), Utils::Color::CYAN, 1);
                }
                if (CrossFilm::active(cfg)) {
                    std::ostringstream xf;
                    xf << "CrossFilm: T_shaft=" << cross_film_state.t_shaft
                       << " K T_bush_mean=" << cross_film_state.t_bush_mean
                       << " K Q_shaft=" << cross_film_state.q_shaft_net
                       << " W Q_bush=" << cross_film_state.q_bush_net
                       << " W core_rise_max=" << cross_film_state.core_rise_max << " K";
                    Utils::log(xf.str(), Utils::Color::CYAN, 1);
                }
            }

            // Conservation/convergence diagnostics. Must run before
            // store_old_time so the step's storage terms are still formed
            // against the previous timestep.
            // No-op unless mass_residual was selected as an output field.
            if (!coupled) Diagnostics::mass_residual(fields, mesh, step_cfg,
                                                     steady_state ? 0.0 : step_cfg.dt);
            const bool diagnostics_now = cfg.diagnostics_interval > 0 &&
                step % cfg.diagnostics_interval == 0;
            Diagnostics::MassBalance balance;
            if (diagnostics_now) {
                if (coupled) {
                    balance.liquid_mass = coupled_stats.liquid_mass;
                    balance.gas_mass = coupled_stats.gas_mass;
                    balance.liquid_boundary_flux = coupled_stats.liquid_boundary_flux;
                    balance.gas_boundary_flux = coupled_stats.gas_boundary_flux;
                    balance.inlet_mass_source = coupled_stats.liquid_inlet_source;
                    balance.liquid_residual = coupled_stats.liquid_balance*step_cfg.dt/std::max(coupled_stats.liquid_mass,1e-30);
                    balance.gas_residual = coupled_stats.species_balance*step_cfg.dt/std::max(coupled_stats.gas_mass,1e-30);
                    balance.gas_volume_excess = coupled_stats.volume_error;
                    if (cfg.temperature_model==TemperatureModel::ENERGY_EQUATION) {
                        balance.energy_residual_w=coupled_stats.energy_balance;
                        balance.internal_energy_j=coupled_stats.internal_energy;
                    }
                    double local_cavity = 0, total_cavity = 0;
                    for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
                        if (fields["theta"](i,j)<1-1e-8) local_cavity += 1;
                    MPI_Allreduce(&local_cavity,&total_cavity,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                    balance.cavitated_fraction = total_cavity/(mesh.n_theta_global*mesh.n_z_global);
                } else balance = Diagnostics::mass_balance(fields, mesh, step_cfg,
                                                            steady_state ? 0.0 : step_cfg.dt);
                if (cfg.groove_mixing_model == GrooveMixingModel::FLUX_BALANCED) {
                    comm.update_ghosts(fields["temperature"]);
                    balance.energy = Energy::balance(fields, mesh, step_cfg);
                }
                if (CrossFilm::active(cfg)) {
                    balance.cross_film_profile_residual = cross_film_state.profile_residual;
                    balance.cross_film_profile_t_min_K = cross_film_state.profile_t_min;
                    balance.cross_film_profile_t_max_K = cross_film_state.profile_t_max;
                }
                Diagnostics::append_csv(
                    std::filesystem::path(cfg.output_dir) / "diagnostics.csv",
                    t, step, balance, elrod_stats,
                    steady_state && !coupled ? &coupling_control.status() : nullptr);
            }

            // Output: write VTK at the configured save cadence ...
            const bool save_now =
                (t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0 ||
                 (!steady_state && t == cfg.end_t));
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
                    if (cfg.cavitation_model == CavitationModel::JFO) {
                        line << " | outer=" << elrod_stats.outer_iters
                             << " flips=" << elrod_stats.flag_flips;
                    } else if (FBNS::active(cfg)) {
                        line << " | newton=" << elrod_stats.outer_iters;
                    }
                    if (diagnostics_now) {
                        line << " | r_l=" << balance.liquid_residual;
                        if (cfg.fluid_property_model == FluidPropertyModel::TWO_PHASE) {
                            line << " r_gas=" << balance.gas_residual;
                            line << " gas_volume_excess=" << balance.gas_volume_excess
                                 << " gas_eos_defect=" << balance.gas_eos_defect;
                        }
                    }
                }
                if (save_now) line << "  [saved]";
                Utils::log(line.str(), Utils::Color::GREEN, 1);
            }

            if (nonlinear_failed) break;

            if (steady_state && cfg.cavitation_model == CavitationModel::FBNS &&
                cfg.fluid_property_model == FluidPropertyModel::SINGLE_PHASE &&
                cfg.gas_pressure_coupling == GasPressureCoupling::NONE &&
                cfg.cavitation_threshold == CavitationThreshold::SCALAR_PCAV &&
                cfg.dissolved_gas_initial == 0.0) {
                try {
                    IO::write_steady_initial_state(fs::path(cfg.output_dir) / "steady_state.dat", fields, mesh, cfg);
                } catch (const std::exception& ex) {
                    Utils::log(ex.what(), Utils::Color::RED, 1);
                    nonlinear_failed = true;
                    break;
                }
            }

            // Store the completed state for the next physical interval.
            fields["pressure"].store_old_time();
            fields["theta"].store_old_time();
            fields["h"].store_old_time();
            fields["rho_liquid_solution"].store_old_time();
            fields["temperature"].store_old_time();
            fields["dissolved_gas"].store_old_time();
            fields["free_gas_mass"].store_old_time();
            if (coupled) coupled_solver->store_old_time(fields);
            previous_forces = forces;
            if (steady_state) ++step;
        }

        Utils::log(nonlinear_failed ? "Simulation stopped with an unconverged or inadmissible state."
                                    : "Simulation finished successfully.",
                   nonlinear_failed ? Utils::Color::RED : Utils::Color::GREEN);
    }  // LinearSystem and Fields destroyed before PetscFinalize

    PetscFinalize();
    return nonlinear_failed ? 4 : 0;
}
