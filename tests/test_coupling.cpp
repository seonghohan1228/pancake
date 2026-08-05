#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "config.hpp"
#include "coupling.hpp"
#include "field.hpp"
#include "fluid_properties.hpp"
#include "mesh.hpp"

namespace {

void check(bool condition, const char* message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

bool has_substr(const std::vector<std::string>& items, const std::string& needle) {
    return std::any_of(items.begin(), items.end(),
                       [&](const std::string& s) { return s.find(needle) != std::string::npos; });
}

SimulationConfig gas_steady_config() {
    SimulationConfig cfg;
    cfg.p_cav = 1.0e5;
    cfg.bc_z_south_val = 1.0e6;
    cfg.bc_z_north_val = 1.0e6;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    cfg.solubility_model = SolubilityModel::HENRY;
    cfg.density_model = DensityModel::CONSTANT;
    cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
    cfg.dissolved_gas_initial = 0.2175;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.175e-7;
    cfg.gas_mass_transfer_rate = 500.0;
    cfg.gas_alpha_max = 0.6;
    cfg.mu_gas = 8.7e-6;
    cfg.mixture_viscosity_model = MixtureViscosityModel::MCADAMS;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;
    cfg.temperature_initial = 313.15;
    return cfg;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // --- steady + gas + FROZEN honesty guard --------------------------------
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig frozen = gas_steady_config();
        frozen.steady_gas_model = SteadyGasModel::FROZEN;
        frozen.allow_frozen_gas_steady = false;
        frozen.validate(errors, warnings);
        check(has_substr(errors, "freezes the gas model"),
              "guard: STEADY_STATE + gas + FROZEN must error");

        SimulationConfig equil = gas_steady_config();
        equil.steady_gas_model = SteadyGasModel::EQUILIBRIUM;
        errors.clear(); warnings.clear();
        equil.validate(errors, warnings);
        check(!has_substr(errors, "freezes the gas model"), "guard: EQUILIBRIUM is allowed");
        check(errors.empty(), "guard: EQUILIBRIUM steady gas config is otherwise valid");

        SimulationConfig allowed = gas_steady_config();
        allowed.steady_gas_model = SteadyGasModel::FROZEN;
        allowed.allow_frozen_gas_steady = true;
        errors.clear(); warnings.clear();
        allowed.validate(errors, warnings);
        check(!has_substr(errors, "freezes the gas model"),
              "guard: allow_frozen_gas_steady opts out of the error");
    }

    // --- WP-6 cavitated_film_model = STRIATED requires Elrod-Adams ----------
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig striated = gas_steady_config();
        striated.steady_gas_model = SteadyGasModel::EQUILIBRIUM;  // avoid the unrelated steady guard
        striated.cavitated_film_model = CavitatedFilmModel::STRIATED;
        striated.cavitation_model = CavitationModel::GUMBEL;
        striated.validate(errors, warnings);
        check(has_substr(errors, "STRIATED requires cavitation_model = JFO"),
              "WP-6: STRIATED + GUMBEL must error");

        striated.cavitation_model = CavitationModel::JFO;
        errors.clear(); warnings.clear();
        striated.validate(errors, warnings);
        check(!has_substr(errors, "STRIATED requires"), "WP-6: STRIATED + JFO is valid");

        // FBNS resolves the same free boundary, so it carries a film content
        // STRIATED can weight by. (Its own phase-1 errors are checked in
        // test_fbns; only the STRIATED gate matters here.)
        striated.cavitation_model = CavitationModel::FBNS;
        errors.clear(); warnings.clear();
        striated.validate(errors, warnings);
        check(!has_substr(errors, "STRIATED requires"), "WP-6: STRIATED + FBNS is valid");
    }

    // --- coupling config round-trip -----------------------------------------
    {
        SimulationConfig cfg;
        cfg.load_from_text(
            "coupling_max_iters = 12\n"
            "coupling_tolerance = 1e-5\n"
            "coupling_relaxation = 0.7\n"
            "steady_gas_model = EQUILIBRIUM\n"
            "allow_frozen_gas_steady = true\n");
        check(cfg.coupling_max_iters == 12, "parse coupling_max_iters");
        check(close(cfg.coupling_tolerance, 1e-5, 1e-12), "parse coupling_tolerance");
        check(close(cfg.coupling_relaxation, 0.7, 1e-12), "parse coupling_relaxation");
        check(cfg.steady_gas_model == SteadyGasModel::EQUILIBRIUM, "parse steady_gas_model");
        check(cfg.allow_frozen_gas_steady, "parse allow_frozen_gas_steady");
    }

    // --- best-iterate config defaults and round-trip ------------------------
    {
        SimulationConfig defaults;
        check(defaults.coupling_keep_best_iterate,
              "best-iterate rollback is on by default (it only fires on cap exhaustion)");
        check(defaults.coupling_relaxation_mode == CouplingRelaxationMode::FIXED,
              "relaxation schedule defaults to FIXED, i.e. the historical behaviour");

        SimulationConfig cfg;
        cfg.load_from_text(
            "coupling_keep_best_iterate = false\n"
            "coupling_relaxation_mode = ADAPTIVE\n"
            "coupling_relaxation_min = 0.05\n"
            "coupling_relaxation_backoff = 0.4\n"
            "coupling_relaxation_growth = 1.25\n"
            "log_coupling_history = true\n");
        check(cfg.log_coupling_history, "parse log_coupling_history");
        check(!defaults.log_coupling_history, "per-sweep coupling history is off by default");
        check(!cfg.coupling_keep_best_iterate, "parse coupling_keep_best_iterate");
        check(cfg.coupling_relaxation_mode == CouplingRelaxationMode::ADAPTIVE,
              "parse coupling_relaxation_mode");
        check(close(cfg.coupling_relaxation_min, 0.05, 1e-12), "parse coupling_relaxation_min");
        check(close(cfg.coupling_relaxation_backoff, 0.4, 1e-12), "parse coupling_relaxation_backoff");
        check(close(cfg.coupling_relaxation_growth, 1.25, 1e-12), "parse coupling_relaxation_growth");

        // ADAPTIVE bounds are validated; a floor above the ceiling is an error.
        std::vector<std::string> errors, warnings;
        SimulationConfig bad;
        bad.solution_mode = SolutionMode::STEADY_STATE;
        bad.coupling_relaxation_mode = CouplingRelaxationMode::ADAPTIVE;
        bad.coupling_relaxation = 0.5;
        bad.coupling_relaxation_min = 0.8;
        bad.validate(errors, warnings);
        check(has_substr(errors, "coupling_relaxation_min"),
              "ADAPTIVE: relaxation floor above the ceiling must error");
    }

    // --- Controller: a converging loop returns its LAST iterate --------------
    // This is the bit-identity guard. A monotonically decreasing history that
    // meets the tolerance must never trigger a rollback, so a converged run is
    // byte-for-byte what it was before best-iterate tracking existed.
    {
        SimulationConfig cfg;
        cfg.coupling_tolerance = 1e-6;
        cfg.coupling_max_iters = 50;
        Coupling::Controller c(cfg);

        const double history[] = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
        int outer = 0;
        bool stopped = false;
        for (double r : history) {
            if (c.record(outer, r, 0.1 * r)) { stopped = true; ++outer; break; }
            check(c.is_best(), "monotone: every sweep is the new best");
            check(close(c.relaxation(), 1.0, 1e-15), "FIXED: relaxation stays at 1");
            ++outer;
        }
        c.finish(outer);
        check(stopped, "monotone: the loop stops when the tolerance is met");
        check(c.status().converged, "monotone: reported converged");
        check(!c.status().cap_exhausted, "monotone: did not exhaust the cap");
        check(c.status().best_iter == c.status().iters,
              "monotone: best iterate IS the last iterate");
        check(!c.should_restore_best(), "monotone: no rollback");
        check(!c.status().restored_best, "monotone: status records no rollback");
    }

    // --- Controller: an oscillating loop returns its BEST iterate ------------
    // The measured signature of the conjugate Ferron cases: the residual falls to
    // a minimum part-way through and then rises and oscillates until the cap.
    {
        SimulationConfig cfg;
        cfg.coupling_tolerance = 1e-9;   // unreachable by this history -> cap exhaustion
        cfg.coupling_max_iters = 10;
        Coupling::Controller c(cfg);

        // Minimum at index 3 (1-based outer 4), then a rising oscillation.
        const double history[] = {1.1e-4, 8.0e-6, 9.0e-7, 1.2e-7,
                                  3.0e-7, 4.5e-7, 3.2e-7, 4.4e-7, 3.1e-7, 4.5e-7};
        int outer = 0;
        for (double r : history) {
            const bool stop = c.record(outer, r, 0.0);
            check(!stop, "oscillation: never meets the tolerance");
            ++outer;
        }
        c.finish(outer);
        check(c.status().iters == 10, "oscillation: ran to the cap");
        check(!c.status().converged, "oscillation: reported NOT converged");
        check(c.status().cap_exhausted, "oscillation: cap exhaustion flagged");
        check(c.status().best_iter == 4, "oscillation: best iterate is outer 4");
        check(close(c.status().best_residual, 1.2e-7, 1e-20), "oscillation: best residual recorded");
        check(close(c.status().R_p, 4.5e-7, 1e-20), "oscillation: last residual recorded");
        check(c.should_restore_best(), "oscillation: rollback requested");
        check(c.status().restored_best, "oscillation: status records the rollback");
        check(c.status().best_to_last_ratio() > 3.0,
              "oscillation: the final iterate is measurably worse than the best");
        check(c.rollback_line().find("NON-CONVERGED") != std::string::npos,
              "oscillation: the rollback line says NON-CONVERGED out loud");
        check(c.summary_line().find("NOT CONVERGED") != std::string::npos,
              "oscillation: the summary line says NOT CONVERGED out loud");

        // Opting out restores the historical (last-iterate) behaviour exactly.
        SimulationConfig off = cfg;
        off.coupling_keep_best_iterate = false;
        Coupling::Controller c2(off);
        outer = 0;
        for (double r : history) { c2.record(outer, r, 0.0); ++outer; }
        c2.finish(outer);
        check(!c2.should_restore_best(), "coupling_keep_best_iterate = false: no rollback");
        check(!c2.status().converged, "coupling_keep_best_iterate = false: still flagged non-converged");
    }

    // --- Controller: a loop whose LAST iterate is its best needs no rollback --
    {
        SimulationConfig cfg;
        cfg.coupling_tolerance = 1e-12;  // unreachable
        cfg.coupling_max_iters = 4;
        Coupling::Controller c(cfg);
        const double history[] = {1e-3, 1e-4, 1e-5, 1e-6};
        int outer = 0;
        for (double r : history) { c.record(outer, r, 0.0); ++outer; }
        c.finish(outer);
        check(!c.status().converged, "slow: not converged");
        check(c.status().cap_exhausted, "slow: cap exhausted");
        check(!c.should_restore_best(), "slow but monotone: the last iterate is the best");
        check(c.rollback_line().find("NON-CONVERGED") != std::string::npos,
              "slow: still reported NON-CONVERGED");
    }

    // --- Controller: ADAPTIVE relaxation backs off on a rising residual -------
    {
        SimulationConfig cfg;
        cfg.coupling_tolerance = 1e-12;
        cfg.coupling_max_iters = 8;
        cfg.coupling_relaxation = 1.0;
        cfg.coupling_relaxation_mode = CouplingRelaxationMode::ADAPTIVE;
        cfg.coupling_relaxation_min = 0.2;
        cfg.coupling_relaxation_backoff = 0.5;
        cfg.coupling_relaxation_growth = 1.1;
        Coupling::Controller c(cfg);

        check(close(c.relaxation(), 1.0, 1e-15), "adaptive: starts at the ceiling");
        c.record(0, 1e-3, 0.0);
        check(close(c.relaxation(), 1.0, 1e-15), "adaptive: first sweep does not move omega");
        c.record(1, 1e-4, 0.0);
        check(close(c.relaxation(), 1.0, 1e-15), "adaptive: falling residual cannot exceed the ceiling");
        c.record(2, 5e-4, 0.0);
        check(close(c.relaxation(), 0.5, 1e-15), "adaptive: rising residual halves omega");
        c.record(3, 6e-4, 0.0);
        check(close(c.relaxation(), 0.25, 1e-15), "adaptive: a second rise halves it again");
        c.record(4, 7e-4, 0.0);
        check(close(c.relaxation(), 0.2, 1e-15), "adaptive: omega is clamped at the floor");
        c.record(5, 1e-5, 0.0);
        check(close(c.relaxation(), 0.22, 1e-12), "adaptive: a falling residual lets omega recover");
        c.finish(6);
        check(close(c.status().relaxation_min_used, 0.2, 1e-15),
              "adaptive: the minimum omega used is reported");
    }

    // --- Snapshot: restoring an iterate is exact, ghosts included ------------
    {
        SimulationConfig cfg;
        cfg.n_theta_global = 12;
        cfg.n_z_global = 6;
        cfg.R = 0.02;
        cfg.c = 4.25e-5;
        cfg.L = 0.04;
        Mesh mesh(cfg);

        Fields fields;
        fields.add("pressure", mesh).fill(1.0);
        fields.add("temperature", mesh).fill(313.15);
        fields.add("velocity_theta", mesh, 2, GridLocation::FACE_THETA).fill(2.0);

        // Write a distinguishable pattern into physical cells AND a ghost cell,
        // then store an old-time state, so the snapshot has something to get wrong.
        for (int i = 0; i < mesh.n_theta_local; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                fields["pressure"](i, j) = 100.0 + i + 0.5 * j;
        fields["pressure"](-1, 0) = -7.25;
        fields["pressure"].store_old_time();
        fields["temperature"](0, 0) = 400.0;

        Coupling::Snapshot best;
        check(!best.valid(), "snapshot: empty until captured");
        best.capture(fields);
        check(best.valid(), "snapshot: valid after capture");
        check(best.bytes() > 0, "snapshot: reports its size");

        // Trash every field the way a further 150 outer iterations would.
        fields["pressure"].fill(-1.0);
        fields["temperature"].fill(-1.0);
        fields["velocity_theta"].fill(-1.0);

        best.restore(fields);
        bool exact = true;
        for (int i = 0; i < mesh.n_theta_local; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                exact = exact && fields["pressure"](i, j) == 100.0 + i + 0.5 * j;
        check(exact, "snapshot: physical cells restored bit-exactly");
        check(fields["pressure"](-1, 0) == -7.25, "snapshot: ghost cells restored too");
        check(fields["pressure"].old(0, 0) == 100.0, "snapshot: old_data restored too");
        check(fields["temperature"](0, 0) == 400.0, "snapshot: every field is restored");
        check(fields["velocity_theta"](0, 0) == 2.0, "snapshot: face fields are restored");
    }

    // --- equilibrate_gas local flash ----------------------------------------
    {
        SimulationConfig cfg = gas_steady_config();
        cfg.steady_gas_model = SteadyGasModel::EQUILIBRIUM;
        cfg.n_theta_global = 8;
        cfg.n_z_global = 4;
        cfg.R = 0.02;
        cfg.c = 4.25e-5;
        cfg.L = 0.04;
        Mesh mesh(cfg);

        Fields fields;
        fields.add("pressure", mesh).fill(5.0e5);  // below p_sat (1.0 MPa) -> outgassing
        fields.add("temperature", mesh).fill(313.15);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("theta", mesh).fill(1.0);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("gas_mass_transfer", mesh).fill(0.0);
        fields.add("rho_gas", mesh).fill(0.0);
        fields.add("inlet_indicator", mesh).fill(0.0);

        FluidProperties::equilibrate_gas(fields, mesh, cfg);
        // c_eq(5e5) = 2.175e-7 * 5e5 = 0.108750 < feed -> some R290 comes out of solution.
        check(close(fields["dissolved_gas"](0, 0), 0.10875, 1e-5),
              "equilibrate_gas: dissolved = local saturation below p_sat");
        check(fields["free_gas_mass"](0, 0) > 1e-9, "equilibrate_gas: free gas released");
        check(fields["alpha_gas"](0, 0) > 0.0 && fields["alpha_gas"](0, 0) <= cfg.gas_alpha_max + 1e-12,
              "equilibrate_gas: void fraction within [0, alpha_max]");

        // At/above p_sat the feed stays fully dissolved (no free gas).
        fields["pressure"].fill(1.0e6);
        fields["dissolved_gas"].fill(cfg.dissolved_gas_initial);
        fields["free_gas_mass"].fill(0.0);
        fields["alpha_gas"].fill(0.0);
        FluidProperties::equilibrate_gas(fields, mesh, cfg);
        check(close(fields["dissolved_gas"](0, 0), cfg.dissolved_gas_initial, 1e-6),
              "equilibrate_gas: fully dissolved at p_sat");
        check(close(fields["free_gas_mass"](0, 0), 0.0, 1e-12), "equilibrate_gas: no free gas at p_sat");
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_coupling: OK\n";
    MPI_Finalize();
    return 0;
}
