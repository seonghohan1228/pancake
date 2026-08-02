// Cross-film temperature profile, conjugate walls, generalized Reynolds, and
// hot-oil carry-over at grooves.
//
// The tests are ordered by what they protect:
//   1. every new config key parses AND survives a save/reload round trip -- a
//      silently-ignored key produced a fake "no effect" result once already;
//   2. the F-integrals reduce to h^3/(12 mu) and h exactly when viscosity does
//      not vary across the gap, so the generalized branch cannot quietly change
//      an isothermal answer;
//   3. an adiabatic shaft carries zero NET heat while carrying nonzero LOCAL
//      flux -- the upstream heat path that the constant-T_wall closure lacks;
//   4. groove mixing at lambda = 0 delivers the upstream film temperature, and
//      at lambda = 1 delivers the supply temperature (the historical pin).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "cross_film.hpp"
#include "energy.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

static double global_max_abs(double local) {
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}

static const char* kKeys =
    "cross_film_thermal_model = POLYNOMIAL\n"
    "cross_film_profile_order = 2\n"
    "cross_film_quadrature_points = 5\n"
    "cross_film_relaxation = 0.45\n"
    "shaft_ambient_heat_transfer = 37.5\n"
    "shaft_ambient_area_ratio = 2.25\n"
    "thermal_ambient_temperature = 305.5\n"
    "bush_thermal_conductivity = 46.5\n"
    "bush_outer_radius = 0.075\n"
    "bush_outer_heat_transfer = 82.5\n"
    "groove_mixing_model = FLOW_WEIGHTED\n"
    "inlet_groove = 12.0 18.0 2e5 313.15 0.4 1.5e-4\n";

static void assert_keys_present(const SimulationConfig& cfg, const char* where) {
    const std::string w = where;
    check(cfg.cross_film_thermal_model == CrossFilmModel::POLYNOMIAL,
          (w + ": cross_film_thermal_model").c_str());
    check(cfg.cross_film_profile_order == 2, (w + ": cross_film_profile_order").c_str());
    check(cfg.cross_film_quadrature_points == 5, (w + ": cross_film_quadrature_points").c_str());
    check(std::abs(cfg.cross_film_relaxation - 0.45) < 1e-12,
          (w + ": cross_film_relaxation").c_str());
    check(std::abs(cfg.shaft_ambient_heat_transfer - 37.5) < 1e-12,
          (w + ": shaft_ambient_heat_transfer").c_str());
    check(std::abs(cfg.shaft_ambient_area_ratio - 2.25) < 1e-12,
          (w + ": shaft_ambient_area_ratio").c_str());
    check(std::abs(cfg.thermal_ambient_temperature - 305.5) < 1e-12,
          (w + ": thermal_ambient_temperature").c_str());
    check(std::abs(cfg.bush_thermal_conductivity - 46.5) < 1e-12,
          (w + ": bush_thermal_conductivity").c_str());
    check(std::abs(cfg.bush_outer_radius - 0.075) < 1e-12,
          (w + ": bush_outer_radius").c_str());
    check(std::abs(cfg.bush_outer_heat_transfer - 82.5) < 1e-12,
          (w + ": bush_outer_heat_transfer").c_str());
    check(cfg.groove_mixing_model == GrooveMixingModel::FLOW_WEIGHTED,
          (w + ": groove_mixing_model").c_str());
    check(cfg.inlets.size() == 1, (w + ": inlet count").c_str());
    check(cfg.inlets[0].feeds_fresh_oil, (w + ": inlet feeds_fresh_oil").c_str());
    check(std::abs(cfg.inlets[0].t_supply - 313.15) < 1e-12, (w + ": inlet t_supply").c_str());
    check(std::abs(cfg.inlets[0].mixing_lambda - 0.4) < 1e-12, (w + ": inlet mixing_lambda").c_str());
    check(std::abs(cfg.inlets[0].supply_flow_rate - 1.5e-4) < 1e-12,
          (w + ": inlet supply_flow_rate").c_str());
}

// Every key must parse, and must still be there after the config is written out
// and read back: an emit-side omission is exactly as damaging as a parse-side one,
// because the emitted file is what a reproduction run is built from.
static void test_config_round_trip() {
    SimulationConfig cfg;
    cfg.load_from_text(kKeys);
    assert_keys_present(cfg, "parse");

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) return;
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "pancake_cross_film_round_trip.txt";
    cfg.save_to_file(path);

    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    in.close();
    SimulationConfig reloaded;
    reloaded.load_from_text(buffer.str());
    assert_keys_present(reloaded, "round trip");
    std::filesystem::remove(path);

    // Defaults must be the historical behaviour, or "off" is not off.
    SimulationConfig defaults;
    check(defaults.cross_film_thermal_model == CrossFilmModel::NONE,
          "cross_film_thermal_model must default to NONE");
    check(defaults.groove_mixing_model == GrooveMixingModel::PINNED,
          "groove_mixing_model must default to PINNED");

    // An inlet given only a supply temperature must keep lambda = 1 (a full pin).
    // std::istream leaves a failed extraction target zeroed, which would silently
    // turn every unannotated groove into a pure carry-over groove.
    SimulationConfig short_form;
    short_form.load_from_text("inlet_groove = 0.0 18.0 2e5 320.0\n");
    check(short_form.inlets.size() == 1, "short-form inlet parse");
    check(std::abs(short_form.inlets[0].mixing_lambda - 1.0) < 1e-12,
          "omitted mixing_lambda must stay 1.0");
    check(short_form.inlets[0].supply_flow_rate == 0.0,
          "omitted supply_flow_rate must stay 0.0");
}

static void build_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    fields.add("pressure", mesh).fill(cfg.p_cav);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("temperature", mesh).fill(cfg.temperature_initial);
    fields.add("heat_generation", mesh).fill(0.0);
    fields.add("mu", mesh).fill(cfg.mu);
    fields.add("rho_cp_liquid_solution", mesh).fill(cfg.rho_cp);
    fields.add("k_liquid_solution", mesh).fill(cfg.thermal_conductivity);
    fields.add("dissolved_gas", mesh).fill(0.0);
    CrossFilm::register_fields(fields, mesh, cfg);
    fields["temperature"].store_old_time();
}

// With liquid_viscosity_model = CONSTANT the reconstructed T(y) cannot change mu,
// so the Dowson integrals must collapse onto the classical coefficients. This is
// the guard that the generalized branch is a generalisation and not a change.
static void test_uniform_viscosity_reduction() {
    SimulationConfig cfg;
    cfg.n_theta_global = 24;
    cfg.n_z_global = 6;
    cfg.omega = 800.0;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    cfg.temperature_initial = 320.0;
    cfg.temperature_reference = 320.0;
    cfg.thermal_ambient_temperature = 320.0;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    build_fields(fields, mesh, cfg);
    // A non-uniform temperature so the profile coefficients are genuinely nonzero.
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            fields["temperature"](i, j) = 320.0 + 5.0 * (mesh.offset_theta + i) / cfg.n_theta_global;
        }
    }
    comm.update_ghosts(fields);
    CrossFilm::update(fields, mesh, cfg);

    const double f2_exact = cfg.c * cfg.c * cfg.c / (12.0 * cfg.mu);
    double err_f2 = 0.0, err_hc = 0.0, max_b = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            err_f2 = std::max(err_f2, std::abs(fields["cross_film_f2"](i, j) / f2_exact - 1.0));
            err_hc = std::max(err_hc, std::abs(fields["cross_film_couette_h"](i, j) / cfg.c - 1.0));
            max_b = std::max(max_b, std::abs(fields["cross_film_b"](i, j)));
        }
    }
    check(global_max_abs(err_f2) < 1e-12, "F2 must reduce to h^3/(12 mu) at uniform mu");
    check(global_max_abs(err_hc) < 1e-12, "2 F1/F0 must reduce to h at uniform mu");
    check(global_max_abs(max_b) > 1e-6, "profile coefficients should not be trivially zero");
}

// The whole point of the conjugate shaft: with no path to ambient, the shaft
// exchanges heat with the film locally but takes none away in total. That net-zero
// with nonzero local flux is what moves heat upstream against the entrainment.
static void test_adiabatic_shaft_balance() {
    SimulationConfig cfg;
    cfg.n_theta_global = 32;
    cfg.n_z_global = 4;
    cfg.omega = 1200.0;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    cfg.shaft_ambient_heat_transfer = 0.0;   // adiabatic shaft
    cfg.bush_outer_heat_transfer = 0.0;      // adiabatic bush
    cfg.temperature_initial = 330.0;
    cfg.temperature_reference = 330.0;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    build_fields(fields, mesh, cfg);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double x = 2.0 * M_PI * (mesh.offset_theta + i) / cfg.n_theta_global;
        for (int j = 0; j < mesh.n_z_local; ++j) {
            fields["temperature"](i, j) = 330.0 + 8.0 * std::sin(x);
        }
    }
    comm.update_ghosts(fields);
    const CrossFilm::WallState state = CrossFilm::update(fields, mesh, cfg);

    const double scale = std::abs(state.q_shaft_net) + 1.0;
    check(std::abs(state.q_shaft_net) / scale < 1e-9,
          "adiabatic shaft must carry zero net heat");
    check(std::abs(state.q_bush_net) < 1e-9,
          "adiabatic bush must carry zero net heat");

    double max_local_flux = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            max_local_flux = std::max(max_local_flux, std::abs(fields["cross_film_wall_flux"](i, j)));
        }
    }
    check(global_max_abs(max_local_flux) > 1.0,
          "zero net shaft heat must come from cancelling local fluxes, not from no flux");
    check(state.t_shaft > 320.0 && state.t_shaft < 340.0,
          "shaft temperature should land inside the film temperature range");
}

// lambda = 0 is pure carry-over (groove delivers the upstream film), lambda = 1 is
// the historical hard pin. Both are checked on the solved field, not on the closure,
// so a mixing value that never reaches the linear system would fail here.
static void test_groove_mixing_limits(double lambda, bool expect_supply) {
    SimulationConfig cfg;
    cfg.n_theta_global = 40;
    cfg.n_z_global = 4;
    cfg.omega = 600.0;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.groove_mixing_model = GrooveMixingModel::FIXED_LAMBDA;
    cfg.temperature_initial = 310.0;
    cfg.temperature_reference = 310.0;
    InletConfig inlet;
    inlet.type = InletConfig::Type::GROOVE;
    inlet.theta = 180.0;
    inlet.size = 18.0;
    inlet.p_supply = 2.0e5;
    inlet.t_supply = 300.0;
    inlet.feeds_fresh_oil = true;
    inlet.mixing_lambda = lambda;
    cfg.inlets.push_back(inlet);

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    build_fields(fields, mesh, cfg);
    // Ramp in theta so the upstream cell has a distinctive temperature.
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            fields["temperature"](i, j) = 310.0 + 0.5 * (mesh.offset_theta + i);
        }
    }
    comm.update_ghosts(fields);

    // The groove spans 18 deg around 180 deg on 40 cells (9 deg each): cells 19 and
    // 20 by centre angle, so the upstream sample is cell 18 at 310 + 9 = 319 K.
    const double t_upstream = 310.0 + 0.5 * 18.0;
    const double expected = expect_supply ? inlet.t_supply : t_upstream;

    LinearSystem sys(mesh);
    Energy::solve(fields, sys, mesh, cfg);

    double err = 0.0;
    int groove_cells = 0;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const int gi = mesh.offset_theta + i;
        if (gi != 19 && gi != 20) continue;
        for (int j = 0; j < mesh.n_z_local; ++j) {
            err = std::max(err, std::abs(fields["temperature"](i, j) - expected));
            ++groove_cells;
        }
    }
    int total_cells = 0;
    MPI_Allreduce(&groove_cells, &total_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    check(total_cells == 2 * cfg.n_z_global, "groove cell count");
    check(global_max_abs(err) < 1e-6,
          expect_supply ? "lambda = 1 must reproduce the hard supply pin"
                        : "lambda = 0 must deliver the upstream film temperature");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_config_round_trip();
    test_uniform_viscosity_reduction();
    test_adiabatic_shaft_balance();
    test_groove_mixing_limits(0.0, false);
    test_groove_mixing_limits(1.0, true);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_cross_film: all checks passed\n";

    PetscFinalize();
    MPI_Finalize();
    return 0;
}
