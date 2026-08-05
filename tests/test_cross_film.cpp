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
//      at lambda = 1 delivers the supply temperature (the historical pin);
//   5. the optional film-surface outputs are the closed-form profile evaluated at
//      the two walls (temperatures) and the conduction slope there (fluxes,
//      positive INTO the fluid), they cost nothing and change nothing when they
//      are not selected, and selecting them without the model is a config error
//      rather than a field of zeros.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "cross_film.hpp"
#include "energy.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "output_naming.hpp"

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
    cfg.output_fields.push_back("heat_flux_journal");
    cfg.output_fields.push_back("heat_flux_bearing");

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
    // The output fluxes see the same limit from the other side: an adiabatic bush
    // conducts nothing anywhere, cell by cell, while the adiabatic shaft still
    // exchanges locally and only cancels in the area integral.
    double max_bearing_flux = 0.0, max_journal_flux = 0.0, journal_integral = 0.0;
    const double cell_area = mesh.R * mesh.get_d_theta() * mesh.get_d_z();
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            max_local_flux = std::max(max_local_flux, std::abs(fields["cross_film_wall_flux"](i, j)));
            const double q_j = fields["heat_flux_journal"](i, j);
            max_bearing_flux = std::max(max_bearing_flux, std::abs(fields["heat_flux_bearing"](i, j)));
            max_journal_flux = std::max(max_journal_flux, std::abs(q_j));
            journal_integral += cell_area * q_j;
        }
    }
    check(global_max_abs(max_local_flux) > 1.0,
          "zero net shaft heat must come from cancelling local fluxes, not from no flux");
    // "Identically zero" up to round-off of the quartic collocation term, which is
    // O(10^3 K) here and cancels between b and a; the journal flux on the same
    // wall sets the scale that residual is measured against.
    const double flux_scale = global_max_abs(max_journal_flux);
    check(global_max_abs(max_bearing_flux) < 1e-10 * flux_scale,
          "an adiabatic bush must give an identically zero bearing heat-flux field");
    check(flux_scale > 1.0,
          "the journal heat-flux field must stay locally nonzero on an adiabatic shaft");
    double total_journal = 0.0;
    MPI_Allreduce(&journal_integral, &total_journal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    check(std::abs(total_journal) < 1e-9 * (1.0 + std::abs(state.q_shaft_net)),
          "the journal heat-flux field must integrate to the zero net adiabatic shaft heat");
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

// -------------------------------------------------------------------------
// Film-surface outputs: temperature_journal / temperature_bearing [K] and
// heat_flux_journal / heat_flux_bearing [W/m^2, positive INTO the fluid].
// -------------------------------------------------------------------------

// The profile of PHYSICS.md 13.1.2, written out independently of CrossFilm so
// the check is against the documented closed form and not against the code.
static double profile_at(double t_bar, double a, double b, double c, double e) {
    return t_bar + a * (e - 0.5) + b * (e - e * e - 1.0 / 6.0) +
           c * (e * e * (1.0 - e) * (1.0 - e) - 1.0 / 30.0);
}

// d T / d eta by the 5-point central stencil, which is exact for a quartic up to
// round-off, so the flux check never has to reuse the analytic wall slopes.
static double profile_slope_at(double t_bar, double a, double b, double c, double e) {
    const double d = 0.01;
    return (-profile_at(t_bar, a, b, c, e + 2 * d) + 8.0 * profile_at(t_bar, a, b, c, e + d) -
            8.0 * profile_at(t_bar, a, b, c, e - d) + profile_at(t_bar, a, b, c, e - 2 * d)) /
           (12.0 * d);
}

static SimulationConfig surface_config() {
    SimulationConfig cfg;
    cfg.n_theta_global = 24;
    cfg.n_z_global = 6;
    cfg.omega = 900.0;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    cfg.cross_film_profile_order = 4;        // exercise the quartic term too
    cfg.temperature_initial = 330.0;
    cfg.temperature_reference = 330.0;
    cfg.thermal_ambient_temperature = 310.0; // cooled bush: the fluxes are nonzero
    cfg.bush_outer_heat_transfer = 80.0;
    return cfg;
}

static void ramp_temperature(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg) {
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        const double x = 2.0 * M_PI * (mesh.offset_theta + i) / cfg.n_theta_global;
        for (int j = 0; j < mesh.n_z_local; ++j) {
            fields["temperature"](i, j) = 330.0 + 6.0 * std::sin(x) + 0.5 * j;
        }
    }
}

// The two surface temperatures must be the profile evaluated at eta = 0 and 1,
// and the two fluxes must be the conduction flux -k dT/dy there with the sign
// flipped to "positive into the fluid". A wrong sign is the failure mode this
// guards: the module's own q_shaft / q_bush are positive OUT of the film.
static void test_surface_fields_match_the_profile() {
    SimulationConfig cfg = surface_config();
    for (const char* name : OutputNaming::surface_fields()) cfg.output_fields.push_back(name);

    Mesh mesh(cfg);
    Communicator comm(mesh);
    Fields fields;
    build_fields(fields, mesh, cfg);
    ramp_temperature(fields, mesh, cfg);
    comm.update_ghosts(fields);
    CrossFilm::update(fields, mesh, cfg);

    for (const char* name : OutputNaming::surface_fields()) {
        check(fields.has(name), (std::string("selected field ") + name + " must be allocated").c_str());
    }

    const double kh = cfg.thermal_conductivity / cfg.c;
    double err_t = 0.0, err_q = 0.0, err_walls = 0.0;
    double max_split = 0.0, max_flux = 0.0, worst_bearing_sign = -1.0e300;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            const double t_bar = fields["temperature"](i, j);
            const double a = fields["cross_film_a"](i, j);
            const double b = fields["cross_film_b"](i, j);
            const double c = fields["cross_film_c"](i, j);

            const double t_j = fields["temperature_journal"](i, j);
            const double t_b = fields["temperature_bearing"](i, j);
            err_t = std::max(err_t, std::abs(t_j - profile_at(t_bar, a, b, c, 0.0)));
            err_t = std::max(err_t, std::abs(t_b - profile_at(t_bar, a, b, c, 1.0)));

            // The wall values are what closes the profile, so the surface fields
            // must also reproduce the conjugate wall temperatures.
            err_walls = std::max(err_walls, std::abs(t_j - fields["cross_film_t_shaft"](i, j)));
            err_walls = std::max(err_walls, std::abs(t_b - fields["cross_film_t_bush"](i, j)));

            const double q_j = fields["heat_flux_journal"](i, j);
            const double q_b = fields["heat_flux_bearing"](i, j);
            // Into the fluid: at eta = 0 the fluid lies at larger eta, at eta = 1
            // at smaller eta, hence the opposite signs on the same slope.
            const double q_j_expected = -kh * profile_slope_at(t_bar, a, b, c, 0.0);
            const double q_b_expected = kh * profile_slope_at(t_bar, a, b, c, 1.0);
            err_q = std::max(err_q, std::abs(q_j - q_j_expected));
            err_q = std::max(err_q, std::abs(q_b - q_b_expected));

            max_split = std::max(max_split, std::abs(t_j - t_b));
            max_flux = std::max(max_flux, std::max(std::abs(q_j), std::abs(q_b)));
            worst_bearing_sign = std::max(worst_bearing_sign, q_b);
        }
    }
    check(global_max_abs(err_t) < 1e-9,
          "surface temperatures must equal the closed-form profile at eta = 0 and 1");
    check(global_max_abs(err_walls) < 1e-9,
          "surface temperatures must reproduce the conjugate wall temperatures");
    check(global_max_abs(err_q) < 1e-6 * (1.0 + global_max_abs(max_flux)),
          "surface heat fluxes must equal -k dT/dy at the wall, positive into the fluid");
    check(global_max_abs(max_split) > 1e-3 && global_max_abs(max_flux) > 1.0,
          "the surface fields must not be trivially equal or zero in this setup");
    // A film hotter than ambient loses heat through the cooled bush, so the
    // into-the-fluid flux there is negative everywhere. This is the assertion a
    // sign-flip regression trips first.
    check(worst_bearing_sign < 0.0,
          "a cooled bush must draw heat out of the film (negative into-fluid flux)");
}

// The fields are opt-in: not selecting them must leave every other field
// bit-for-bit unchanged, so an existing run's output cannot move because the
// feature exists. Compared on the raw storage, ghosts included.
static void test_surface_fields_are_opt_in() {
    SimulationConfig with = surface_config();
    for (const char* name : OutputNaming::surface_fields()) with.output_fields.push_back(name);
    SimulationConfig without = surface_config();

    Mesh mesh(with);
    Communicator comm(mesh);
    Fields fields_with, fields_without;
    build_fields(fields_with, mesh, with);
    build_fields(fields_without, mesh, without);
    ramp_temperature(fields_with, mesh, with);
    ramp_temperature(fields_without, mesh, without);
    comm.update_ghosts(fields_with);
    comm.update_ghosts(fields_without);
    CrossFilm::update(fields_with, mesh, with);
    CrossFilm::update(fields_without, mesh, without);

    for (const char* name : OutputNaming::surface_fields()) {
        check(!fields_without.has(name),
              (std::string("unselected field ") + name + " must not be allocated").c_str());
    }
    int compared = 0;
    for (auto& [name, field] : fields_without) {
        check(fields_with.has(name), (name + " missing from the selected run").c_str());
        const Field& other = fields_with[name];
        check(field->data == other.data, (name + " must be byte-identical").c_str());
        check(field->old_data == other.old_data, (name + " old_data must be byte-identical").c_str());
        ++compared;
    }
    check(compared > 5, "the byte-identity check must actually compare fields");
}

// Selecting a surface field without the model that defines it is a config error,
// not a silently written field of zeros.
static void test_surface_field_selection_is_validated() {
    auto errors_for = [](const SimulationConfig& cfg) {
        std::vector<std::string> errors, warnings;
        SimulationConfig copy = cfg;
        copy.validate(errors, warnings);
        return errors;
    };

    SimulationConfig none;  // defaults: cross_film_thermal_model = NONE
    const size_t baseline = errors_for(none).size();
    for (const char* name : OutputNaming::surface_fields()) {
        SimulationConfig cfg = none;
        cfg.output_fields.push_back(name);
        const std::vector<std::string> errors = errors_for(cfg);
        check(errors.size() > baseline,
              (std::string(name) + " under cross_film NONE must be a config error").c_str());
        bool named = false;
        for (const auto& message : errors) named = named || message.find(name) != std::string::npos;
        check(named, (std::string("the error message must name ") + name).c_str());
    }

    SimulationConfig active = surface_config();
    const size_t active_baseline = errors_for(active).size();
    for (const char* name : OutputNaming::surface_fields()) active.output_fields.push_back(name);
    check(errors_for(active).size() == active_baseline,
          "the surface fields must validate cleanly with the cross-film model on");
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    test_config_round_trip();
    test_uniform_viscosity_reduction();
    test_adiabatic_shaft_balance();
    test_groove_mixing_limits(0.0, false);
    test_groove_mixing_limits(1.0, true);
    test_surface_fields_match_the_profile();
    test_surface_fields_are_opt_in();
    test_surface_field_selection_is_validated();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "test_cross_film: all checks passed\n";

    PetscFinalize();
    MPI_Finalize();
    return 0;
}
