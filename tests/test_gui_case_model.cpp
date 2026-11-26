// GUI case model: fluid tables, config mappings, text-preserving save,
// estimates and angle references. Plain C++, no MPI or toolkit.
//
// Test 1 - Fluid pair: the FW68D table gives the deck-run state (compressor:
//   saturated content and viscosity at 24.5 bar, 70 degC; tester: p_sat at 27.4 %).
// Test 2 - A new compressor case maps to the deck B0 config's physics keys.
// Test 3 - Round trip (acceptance test 1): every deck config opens and saves
//   back equal in every solver key, with untouched lines byte-identical.
// Test 4 - An edit changes only the keys it owns.
// Test 5 - Estimates: Sommerfeld number and saturation pressure of the deck.
// Test 6 - Angle references.
// Test 7 - New case defaults (1 core) and the case -> config -> case mapping.

#include "angle_reference.hpp"
#include "case_model.hpp"
#include "config_text_patch.hpp"
#include "estimates.hpp"
#include "fluid_library.hpp"
#include "key_results.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace gui;

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        ++failures;
    }
}

bool close(double a, double b, double rel) {
    return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b));
}

std::string read_file(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    std::ostringstream text;
    text << file.rdbuf();
    return text.str();
}

const fs::path kData = fs::path(PANCAKE_SOURCE_DIR) / "tests" / "data" / "lge";

FluidPair fw68d() {
    FluidPair pair;
    std::string error;
    const bool ok = load_fluid_pair(kData / "fw68d_pair", pair, error);
    check(ok, "FW68D test pair loads: " + error);
    return pair;
}

std::vector<std::string> lines_of(const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) lines.push_back(line);
    return lines;
}

// Deck compressor B0 mean state (parameters.json of the campaign).
constexpr double kShellPressure = 2451662.5;
constexpr double kShellTemperature = 343.15;
constexpr double kDeckContent = 0.22535420248737287;
constexpr double kDeckViscosity = 0.0050849678203342035;

CaseModel compressor_case() {
    CaseModel m;
    m.bearing_diameter = 0.021;
    m.bearing_length = 0.020;
    m.radial_clearance = 30e-6;
    m.feed.type = FeedType::Hole;
    m.feed.angle_deg = 90.0;
    m.feed.hole_diameter = 2.7e-3;
    m.feed.supply_overpressure = 0.3e5;
    m.speed = 2.0 * 3.14159265358979323846 * 130.0;
    m.load = 3450.192492003545;
    m.load_direction_deg = 270.0;  // on the shaft, downward (bearing load +y in the solver)
    m.ambient_pressure = kShellPressure;
    m.ambient_temperature = kShellTemperature;
    m.fluid_pair = "fw68d_pair";
    m.saturated = true;
    m.dissolved_gas = true;
    m.mesh_level = MeshLevel::Fine;
    return m;
}

void test_fluid_pair() {
    const FluidPair pair = fw68d();
    const PropertyValue c = pair.saturated_content(kShellPressure, kShellTemperature);
    check(close(c.value, kDeckContent, 1e-9) && c.in_range,
          "saturated content at 24.5 bar, 70 degC equals the deck value 0.225354");
    const PropertyValue mu = pair.dynamic_viscosity(kDeckContent, kShellTemperature);
    check(close(mu.value, kDeckViscosity, 1e-9) && mu.in_range,
          "viscosity at the deck state equals the driver value 5.085 mPa s");
    const PropertyValue p_sat = pair.saturation_pressure(kDeckContent, kShellTemperature);
    check(close(p_sat.value, kShellPressure, 1e-6), "p_sat inverts c_sat at the deck state");
    const PropertyValue out = pair.saturated_content(4.5e6, kShellTemperature);
    check(!out.in_range, "a pressure above the table is flagged out of range");
    // The shipped synthetic example loads with the solver's loader and is usable
    // at a compressor-like state (20 bar, 70 degC).
    FluidPair example;
    std::string error;
    check(load_fluid_pair(fs::path(PANCAKE_SOURCE_DIR) / "fluids" / "example_vg68_rx", example, error),
          "shipped example pair loads: " + error);
    const PropertyValue example_c = example.saturated_content(20e5, 343.15);
    check(example_c.in_range && example_c.value > 0.05 && example_c.value < 0.5,
          "example pair: saturated content at 20 bar, 70 degC is inside the table");
    check(example.note.find("NOT FOR DESIGN") != std::string::npos, "example pair is labelled");

    const PropertyValue hot = pair.saturation_pressure(0.3, 473.15);
    check(!hot.in_range, "200 degC is above the table and flagged out of range");
    const PropertyValue rich = pair.saturation_pressure(0.95, 293.15);
    check(rich.in_range, "95 wt % at 20 degC lies inside the FW68D table");
}

void test_new_compressor_case_matches_deck() {
    const FluidPair pair = fw68d();
    const CaseModel m = compressor_case();
    const SimulationConfig cfg = to_config(m, &pair, new_case_template());
    SimulationConfig deck;
    deck.load_from_text(read_file(kData / "configs" / "compressor_gcvc_steady.txt"));
    check(cfg.R == deck.R && cfg.L == deck.L && cfg.c == deck.c, "geometry keys");
    check(close(cfg.omega, deck.omega, 1e-12), "speed");
    check(close(cfg.mu, deck.mu, 1e-9) && close(cfg.rho, deck.rho, 1e-12), "mu and rho");
    check(cfg.p_cav == deck.p_cav, "oil vapour pressure");
    check(close(cfg.dissolved_gas_initial, deck.dissolved_gas_initial, 1e-9), "content");
    check(cfg.bc_z_south_val == deck.bc_z_south_val && cfg.bc_z_north_val == deck.bc_z_north_val,
          "end pressures");
    check(cfg.cavitation_model == CavitationModel::FBNS &&
              cfg.cavitation_threshold == CavitationThreshold::LOCAL_PSAT &&
              cfg.fluid_property_model == FluidPropertyModel::SINGLE_PHASE,
          "GC+VC plateau model form");
    check(cfg.dissolved_gas_species == DissolvedGasSpecies::PROPANE,
          "FBNS keeps the propane label of the deck runs");
    check(cfg.liquid_viscosity_model == LiquidViscosityModel::CONSTANT &&
              cfg.density_model == DensityModel::CONSTANT,
          "isothermal properties fixed at the ambient state");
    check(cfg.inlets.size() == 1 && cfg.inlets[0].type == InletConfig::Type::CIRCULAR &&
              cfg.inlets[0].theta == 90.0 && close(cfg.inlets[0].size, 1.35e-3, 1e-12) &&
              close(cfg.inlets[0].z, 0.01, 1e-12) && cfg.inlets[0].p_supply == 2481662.5 &&
              cfg.inlets[0].t_supply == kShellTemperature,
          "feed hole matches the deck inlet line");
    check(cfg.n_theta_global == 720 && cfg.n_z_global == 192, "fine mesh is the deck grid");
    check(cfg.motion_model == MotionModel::STEADY_EQUILIBRIUM &&
              close(cfg.external_load_y, 3450.192492003545, 1e-12) &&
              std::abs(cfg.external_load_x) < 1e-9,
          "steady load balance with the load along +y");
    std::vector<std::string> errors, warnings;
    cfg.validate(errors, warnings);
    check(errors.empty(), "the new compressor case validates: " +
                              (errors.empty() ? std::string() : errors.front()));
}

void test_round_trip_of_deck_configs() {
    const FluidPair pair = fw68d();
    for (const auto& entry : fs::directory_iterator(kData / "configs")) {
        const std::string name = entry.path().filename().string();
        const std::string text = read_file(entry.path());
        SimulationConfig cfg;
        cfg.load_from_text(text);
        check(cfg.parse_errors.empty(), name + ": parses");
        const GuiMeta meta = read_gui_meta(text);
        const CaseRead read = from_config(cfg, meta, &pair);
        const SimulationConfig updated = to_config(read.model, &pair, cfg, read.model);
        const std::string saved = patch_config_text(text, cfg, updated, case_meta(read.model));
        SimulationConfig reloaded;
        reloaded.load_from_text(saved);
        const auto before = config_key_values(cfg.to_config_text());
        const auto after = config_key_values(reloaded.to_config_text());
        std::size_t differing = 0;
        for (const auto& [key, value] : before) {
            const auto it = after.find(key);
            if (it == after.end() || it->second != value) {
                ++differing;
                std::cerr << "  " << name << ": key " << key << " changed\n";
            }
        }
        check(differing == 0, name + ": every solver key equal after open and save");
        // Untouched lines stay byte-identical (CRLF and afml-r6 paths included).
        const auto original_lines = lines_of(text);
        const auto saved_lines = lines_of(saved);
        check(saved_lines.size() >= original_lines.size() &&
                  std::equal(original_lines.begin(), original_lines.end(), saved_lines.begin()),
              name + ": original lines kept byte-identical, GUI lines appended");
        std::cout << "  round trip " << name << ": " << before.size() << " keys equal, "
                  << read.hidden_keys.size() << " not shown in the GUI\n";
    }
}

void test_edit_changes_only_its_keys() {
    const FluidPair pair = fw68d();
    const std::string text = read_file(kData / "configs" / "compressor_gcvc_steady.txt");
    SimulationConfig cfg;
    cfg.load_from_text(text);
    const CaseRead read = from_config(cfg, read_gui_meta(text), &pair);
    CaseModel edited = read.model;
    edited.speed *= 0.5;
    const SimulationConfig updated = to_config(edited, &pair, cfg, read.model);
    const auto before = config_key_values(cfg.to_config_text());
    const auto after = config_key_values(updated.to_config_text());
    std::vector<std::string> changed;
    for (const auto& [key, value] : after)
        if (before.at(key) != value) changed.push_back(key);
    check(changed == std::vector<std::string>{"omega"},
          "changing the speed of a steady case rewrites only omega");

    CaseModel gas_off = read.model;
    gas_off.dissolved_gas = false;
    const SimulationConfig vc = to_config(gas_off, &pair, cfg, read.model);
    SimulationConfig deck_vc;
    deck_vc.load_from_text(read_file(kData / "configs" / "compressor_vc_steady.txt"));
    check(vc.cavitation_threshold == CavitationThreshold::SCALAR_PCAV &&
              vc.p_cav == deck_vc.p_cav && vc.dissolved_gas_initial == deck_vc.dissolved_gas_initial,
          "dissolved gas off gives the deck VC-only form (SCALAR_PCAV, content kept)");
}

void test_estimates() {
    const FluidPair pair = fw68d();
    // Tester, 500 rpm, 2,610 N, c = 0.274 at the P5g film temperature.
    CaseModel tester;
    tester.bearing_diameter = 0.043;
    tester.bearing_length = 0.040;
    tester.radial_clearance = 42.5e-6;
    tester.speed = 52.35987755982988;
    tester.load = 2610.276612252356;
    tester.ambient_pressure = 1657756.124499086;
    tester.ambient_temperature = 297.5390033315531;
    tester.saturated = false;
    tester.content = 0.274;
    const Estimates e = compute_estimates(tester, &pair);
    std::cout << "  tester estimates: S = " << e.sommerfeld << ", p_sat = "
              << e.saturation_pressure / 1e5 << " bar, margin/P = " << e.margin_ratio << "\n";
    check(e.sommerfeld > 0.0150 && e.sommerfeld < 0.0165, "tester S rounds to the deck's 0.016");
    check(close(e.saturation_pressure, 9.653e5, 2e-3), "tester p_sat ~ 9.65 bar at 27.4 wt %");
    check(e.verdict == OnsetVerdict::Expected, "tester 2,610 N: gas release expected (deck p43)");

    const Estimates compressor = compute_estimates(compressor_case(), &pair);
    check(std::abs(compressor.pressure_margin) < 1.0, "saturated feed has zero pressure margin");
    check(compressor.verdict == OnsetVerdict::Expected,
          "compressor example (zero margin, S below the chart): gas release expected");

    check(onset_verdict(2.0, 1.0) == OnsetVerdict::NotExpected, "margin above curve: not expected");
    check(onset_verdict(1.1, 1.0) == OnsetVerdict::Borderline, "within 20 %: borderline");
    check(close(tester_onset_curve().phi(0.2), 1.2137, 1e-9), "curve node value");
}

void test_angle_references() {
    // Compressor-like view: feed at 90 deg, h_min at 109.6 + 180 deg.
    const double feed = 90.0, hmax = 109.6, hmin = 289.6;
    const double zero_feed = angle_zero_position(AngleZero::Feed, feed, hmax, hmin);
    check(display_angle(feed, zero_feed, angle_axis_start(AngleZero::Feed)) == 0.0,
          "feed reference puts the feed at 0");
    check(close(display_angle(hmin, zero_feed, -240.0), -160.4, 1e-9),
          "h_min line at theta_in = -160.4 deg in the deck window [-240, 120)");
    const double zero_hmax = angle_zero_position(AngleZero::MaxFilm, feed, hmax, hmin);
    check(close(display_angle(hmin, zero_hmax, 0.0), 180.0, 1e-9),
          "h_min is at theta = 180 from the maximum film");
    check(close(display_angle(10.0, 0.0, 0.0, -1.0), 350.0, 1e-9),
          "reverse rotation flips the positive direction");
}

void test_new_case_defaults_and_mapping() {
    const CaseModel defaults;
    check(defaults.cores == 1, "a new case uses 1 core");
    const FluidPair pair = fw68d();
    CaseModel m = compressor_case();
    m.thermal = true;
    m.tilt_angle = 1e-4;
    m.tilt_direction_deg = 30.0;
    m.analysis = AnalysisType::Transient;
    m.load_table.file = kData / "os_load_130hz.csv";
    m.load_table.frequency_hz = 130.0;
    m.moving_mass = 0.686;
    m.revolutions = 2.0;
    m.steps_per_revolution = 1440;
    const SimulationConfig cfg = to_config(m, &pair, new_case_template());
    const CaseRead back = from_config(cfg, case_meta(m), &pair);
    const CaseModel& r = back.model;
    check(close(r.bearing_diameter, m.bearing_diameter, 1e-12) && r.feed == m.feed,
          "geometry and feed read back");
    check(close(r.tilt_angle, m.tilt_angle, 1e-9) && close(r.tilt_direction_deg, 30.0, 1e-9),
          "tilt reads back");
    check(r.analysis == AnalysisType::Transient && r.load_table == m.load_table &&
              r.moving_mass == 0.686 && close(r.revolutions, 2.0, 1e-9) &&
              r.steps_per_revolution == 1440,
          "transient inputs read back");
    check(r.thermal && r.saturated && r.dissolved_gas && r.mesh_level == MeshLevel::Fine,
          "model choices read back");
    check(close(cfg.dt, 1.0 / 130.0 / 1440.0, 1e-12), "1440 steps per revolution at 130 Hz");
    // Thermal transient (rev. 3): simplified thermal model, FBNS, validates.
    check(unsupported_reason(m, &pair).empty(), "thermal transient is supported");
    check(cfg.cross_film_thermal_model == CrossFilmModel::POLYNOMIAL &&
              cfg.bush_thermal_model == BushThermalModel::CONDUCTION_2D &&
              cfg.cavitation_model == CavitationModel::FBNS,
          "thermal transient uses the simplified thermal model with FBNS");
    std::vector<std::string> errors, warnings;
    cfg.validate(errors, warnings);
    check(errors.empty(), "thermal transient validates: " +
                              (errors.empty() ? std::string() : errors.front()));
    check(cfg.cavitated_film_model == CavitatedFilmModel::STRIATED,
          "every FBNS case uses the STRIATED cavity shear (plan R3)");

    // Supported thermal forms validate: gas-free FBNS (O-runs) and COUPLED (X-runs).
    for (bool gas : {false, true}) {
        CaseModel steady = compressor_case();
        steady.thermal = true;
        steady.dissolved_gas = gas;
        check(uses_coupled(steady, &pair), "thermal with refrigerant in the oil runs COUPLED");
        const SimulationConfig thermal_cfg = to_config(steady, &pair, new_case_template());
        std::vector<std::string> e2, w2;
        thermal_cfg.validate(e2, w2);
        check(e2.empty(), std::string("steady thermal case validates (gas ") + (gas ? "on" : "off") +
                              "): " + (e2.empty() ? std::string() : e2.front()));
    }
    // Transient isothermal plateau case (compressor rev03 form) validates.
    CaseModel transient = compressor_case();
    transient.analysis = AnalysisType::Transient;
    transient.load_table.file = kData / "os_load_130hz.csv";
    transient.load_table.frequency_hz = 130.0;
    transient.moving_mass = 0.686;
    const SimulationConfig transient_cfg = to_config(transient, &pair, new_case_template());
    std::vector<std::string> e3, w3;
    transient_cfg.validate(e3, w3);
    check(e3.empty(), "transient plateau case validates: " + (e3.empty() ? std::string() : e3.front()));

    // A pair without a gas block runs FBNS only: gas transport is red with the reason.
    FluidPair example;
    std::string error;
    load_fluid_pair(fs::path(PANCAKE_SOURCE_DIR) / "fluids" / "example_vg68_rx", example, error);
    check(!example.has_gas_phase, "the shipped example pair has no gas block (plan R5)");
    CaseModel transport = compressor_case();
    transport.fluid_pair = example.id;
    transport.gas_model = GasModel::Transport;
    check(unsupported_reason(transport, &example).find("no gas-phase data") != std::string::npos,
          "gas transport with a pair without gas data is red");
    CaseModel thermal_example = transport;
    thermal_example.gas_model = GasModel::Plateau;
    thermal_example.thermal = true;
    check(unsupported_reason(thermal_example, &example).find("gas-phase data") != std::string::npos,
          "steady thermal with refrigerant needs gas data: red for the example pair");
    check(pair.has_gas_phase && pair.solver_species == DissolvedGasSpecies::R410A,
          "the FW68D test pair names the R410A gas model");
}

// Test 11: an FBNS result never shows a nonzero gas field.
void test_fbns_result_fields() {
    ResultGrid grid;
    grid.nx = 4;
    grid.ny = 1;
    grid.theta_deg = {45, 135, 225, 315};
    grid.z_m = {0.01};
    grid.fields["p"] = {2e6, 3e6, 2.4e6, 2.4e6};
    grid.fields["h"] = {3e-5, 2e-5, 1e-5, 2e-5};
    grid.fields["film_content"] = {1.0, 1.0, 0.4, 0.7};
    grid.fields["alpha_gas"] = {0, 0, 0, 0};             // as written by FBNS
    grid.fields["dissolved_gas"] = {0.2, 0.2, 0.2, 0.2};  // constant in FBNS
    grid.fields["rho"] = {900, 900, 900, 900};
    SimulationConfig fbns;
    fbns.cavitation_model = CavitationModel::FBNS;
    fbns.rho = 900.0;
    bool void_found = false;
    for (const FieldView& view : result_fields(fbns, grid)) {
        check(view.id != "alpha_gas" && view.id != "dissolved_gas" && view.id != "rho" &&
                  view.id != "gas_mass_transfer",
              "FBNS result offers no gas field (" + view.id + ")");
        if (view.id == "void") {
            void_found = true;
            const auto values = field_view_values(view, grid, fbns);
            check(values.size() == 4 && close(values[2], 0.6, 1e-12),
                  "void fraction is 1 - theta");
        }
    }
    check(void_found, "FBNS result offers the void fraction");
    check(field_view_values({"alpha_gas", "", "", ""}, grid, fbns).empty(),
          "asking an FBNS result for alpha_gas returns nothing");
    SimulationConfig coupled;
    coupled.cavitation_model = CavitationModel::COUPLED;
    bool gas_found = false;
    for (const FieldView& view : result_fields(coupled, grid)) gas_found |= view.id == "alpha_gas";
    check(gas_found, "COUPLED result offers the free-gas fraction");
}

// Run-time model (plan 0.2): scaling with cells, cores and steps per revolution.
void test_run_time_model() {
    const FluidPair pair = fw68d();
    CaseModel m = compressor_case();
    m.cores = 1;
    const RunTime one = estimate_run_time(m, &pair);
    check(close(one.typical_s, 600.0, 1e-9), "isothermal FBNS steady at 720 x 192, 1 core: 10 min");
    m.cores = 8;
    const RunTime eight = estimate_run_time(m, &pair);
    check(close(eight.typical_s, 600.0 * 1.14 / 8.0, 1e-9), "8 cores with 2 % efficiency loss each");
    m.mesh_level = MeshLevel::Medium;
    check(close(estimate_run_time(m, &pair).typical_s, eight.typical_s / 4.0, 1e-9),
          "a quarter of the cells, a quarter of the time");
    CaseModel t = compressor_case();
    t.analysis = AnalysisType::Transient;
    t.revolutions = 1.0;
    t.cores = 1;
    t.steps_per_revolution = 1440;
    check(close(estimate_run_time(t, &pair).typical_s / 3600.0, 3.7, 0.01),
          "transient at 0.25 deg: 3.7 core-h per revolution");
    t.steps_per_revolution = 360;
    check(close(estimate_run_time(t, &pair).typical_s / 3600.0, 0.5, 1e-9),
          "transient at 1 deg: 0.5 core-h per revolution");
}

}  // namespace

// Writes the shipped example cases (synthetic pair; no LGE data):
//   test_gui_case_model --write-examples <dir>
int write_examples(const fs::path& dir) {
    FluidPair pair;
    std::string error;
    if (!load_fluid_pair(fs::path(PANCAKE_SOURCE_DIR) / "fluids" / "example_vg68_rx", pair, error)) {
        std::cerr << error << "\n";
        return 1;
    }
    CaseModel steady;  // generic defaults: D 40 mm, L 20 mm, c 20 um, 3000 rpm
    steady.fluid_pair = pair.id;
    steady.load = 2000.0;
    steady.ambient_pressure = 15e5;
    steady.ambient_temperature = 333.15;
    steady.feed.type = FeedType::Hole;
    steady.feed.hole_diameter = 3e-3;
    steady.mesh_level = MeshLevel::Medium;
    CaseModel transient = steady;
    transient.analysis = AnalysisType::Transient;
    transient.load_table.file = "load_table_example.csv";
    transient.load_table.frequency_hz = 50.0;
    transient.moving_mass = 0.5;
    transient.revolutions = 2.0;
    transient.steps_per_revolution = 720;
    transient.mesh_level = MeshLevel::Coarse;
    const std::pair<const char*, CaseModel> cases[] = {{"example_steady.txt", steady},
                                                       {"example_transient.txt", transient}};
    for (const auto& [name, model] : cases) {
        SimulationConfig cfg = to_config(model, &pair, new_case_template());
        cfg.solution_table_file = fs::path("../fluids/example_vg68_rx/solution.csv");
        cfg.output_dir = fs::path(name).stem().string() + "_results";
        std::ofstream out(dir / name, std::ios::binary);
        out << new_config_text(cfg, case_meta(model));
        std::cout << "wrote " << (dir / name).string() << "\n";
    }
    return 0;
}

// Writes the acceptance-test cases (test 2) through the GUI case model with the
// FW68D test pair. Inputs are the deck-run values (tests/data/lge/README.md):
//   test_gui_case_model --write-acceptance <dir>
int write_acceptance(const fs::path& dir) {
    const FluidPair pair = fw68d();
    if (!pair.loaded()) return 1;
    CaseModel gcvc = compressor_case();
    gcvc.convergence.force_tolerance = 1e-4;  // the driver's tolerance (part_b.py FORCE_TOL)
    CaseModel vc = gcvc;
    vc.dissolved_gas = false;
    CaseModel transient = gcvc;
    transient.analysis = AnalysisType::Transient;
    transient.load_table.file = kData / "os_load_130hz.csv";
    transient.load_table.frequency_hz = 130.0;
    transient.moving_mass = 0.686;
    transient.revolutions = 3.0;           // the deck accepted the third revolution (rev03)
    transient.steps_per_revolution = 1440;
    transient.outputs_per_revolution = 36;
    // Tester, 500 rpm, 2,610 N, c = 0.274 (P5g inputs), in the GUI's plateau form (D5).
    CaseModel tester = compressor_case();
    tester.bearing_diameter = 0.043;
    tester.bearing_length = 0.040;
    tester.radial_clearance = 42.5e-6;
    tester.feed.type = FeedType::Groove;
    tester.feed.angle_deg = 90.0;
    tester.feed.groove_width_deg = 13.298312524795712;
    tester.feed.supply_overpressure = 0.0;
    tester.speed = 52.35987755982988;
    tester.load = 2610.276612252356;
    tester.ambient_pressure = 1657756.124499086;
    tester.ambient_temperature = 297.5390033315531;
    tester.saturated = false;
    tester.content = 0.274;
    tester.convergence.force_tolerance = 1e-3;  // run_fbns.py FORCE_TOL
    // Test 10: a thermal steady case on this machine. Pure oil (content 0) with
    // the synthetic example pair, deck O-run form (FBNS, finite-volume film,
    // 3-D housing), medium mesh.
    FluidPair example;
    std::string error;
    load_fluid_pair(fs::path(PANCAKE_SOURCE_DIR) / "fluids" / "example_vg68_rx", example, error);
    CaseModel thermal;
    thermal.fluid_pair = example.id;
    thermal.saturated = false;
    thermal.content = 0.0;
    thermal.dissolved_gas = false;
    thermal.thermal = true;
    thermal.load = 2000.0;
    thermal.ambient_pressure = 15e5;
    thermal.feed.hole_diameter = 3e-3;
    thermal.housing.outer_diameter = 0.060;
    thermal.mesh_level = MeshLevel::Medium;
    thermal.cores = 8;
    const std::pair<const char*, CaseModel> cases[] = {
        {"compressor_gcvc_steady.txt", gcvc}, {"compressor_vc_steady.txt", vc},
        {"compressor_transient.txt", transient}, {"tester_500rpm_2610N.txt", tester},
        {"example_thermal_steady_pure_oil.txt", thermal}};
    for (const auto& [name, model] : cases) {
        const FluidPair& used = model.fluid_pair == example.id ? example : pair;
        SimulationConfig cfg = to_config(model, &used, new_case_template());
        cfg.output_dir = fs::path(name).stem().string() + "_results";
        std::ofstream out(dir / name, std::ios::binary);
        out << new_config_text(cfg, case_meta(model));
        std::cout << "wrote " << (dir / name).string() << "\n";
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "--write-examples") return write_examples(argv[2]);
    if (argc == 3 && std::string(argv[1]) == "--write-acceptance") return write_acceptance(argv[2]);
    std::cout << "Running test_gui_case_model...\n";
    test_fluid_pair();
    test_new_compressor_case_matches_deck();
    test_round_trip_of_deck_configs();
    test_edit_changes_only_its_keys();
    test_estimates();
    test_angle_references();
    test_new_case_defaults_and_mapping();
    test_fbns_result_fields();
    test_run_time_model();
    if (failures) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "PASS: test_gui_case_model\n";
    return 0;
}
