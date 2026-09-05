// Config-key parsing contract.
//
// Test 1 — The shipped silent drop: `cross_film_model = NONE` (the key the
//   Ferron THD case templates carried through a whole validation campaign) is
//   rejected by name and points at `cross_film_thermal_model`.
// Test 2 — Deliberate legacy branches keep working: value aliases, key aliases,
//   and the warn-and-drop keys for removed models stay warnings, not errors.
// Test 3 — Manifest vs parser: every key in known_config_keys() is accepted.
// Test 4 — Round trip: to_config_text() writes only keys the parser knows.
// Test 5 — Every config file shipped in the tree parses without an unknown key.
//
// Run as: mpirun -n <N> ./test_config_keys

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"

namespace {

namespace fs = std::filesystem;

void check(bool condition, const std::string& message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

bool mentions(const std::vector<std::string>& messages, const std::string& text) {
    return std::any_of(messages.begin(), messages.end(), [&](const std::string& m) {
        return m.find(text) != std::string::npos;
    });
}

// The exact line that shipped in cases/01_ferron_thd/config_template.txt: the
// parser only ever knew `cross_film_thermal_model`, so this was dropped without
// a word on every run and on ~100 generated fb_iter_*/config.txt copies.
void test_shipped_silent_drop_is_now_fatal() {
    SimulationConfig cfg;
    cfg.load_from_text("cross_film_model = NONE\n");

    check(cfg.parse_errors.size() == 1,
          "cross_film_model must produce exactly one parse error");
    check(mentions(cfg.parse_errors, "cross_film_model"),
          "the unknown-key error must name the offending key");
    check(mentions(cfg.parse_errors, "cross_film_thermal_model"),
          "the unknown-key error must suggest the closest known key");

    std::vector<std::string> errors, warnings;
    cfg.validate(errors, warnings);
    check(mentions(errors, "cross_film_model"),
          "validate() must promote the unknown key to a hard error");

    // A close misspelling of a long key still resolves; an unrelated token does
    // not attract a bogus suggestion.
    SimulationConfig typo;
    typo.load_from_text("bush_outer_radius_ = 0.04\ntotally_made_up_key = 3\n");
    check(typo.parse_errors.size() == 2, "both unknown keys must be reported");
    check(mentions(typo.parse_errors, "did you mean 'bush_outer_radius'"),
          "a one-character typo must suggest the intended key");
    check(!mentions({typo.parse_errors[1]}, "did you mean"),
          "an unrelated key must not be given a spurious suggestion");
}

void test_legacy_branches_still_parse() {
    SimulationConfig cfg;
    cfg.load_from_text(
        "viscosity_model = CONSTANT\n"              // alias of liquid_viscosity_model
        "oil_gas_solution_model = HENRY\n"          // alias of solubility_model
        "gas_mixture_viscosity_model = DUKLER\n"    // alias of mixture_viscosity_model
        "cavitation_model = ELROD_ADAMS\n"          // value alias of JFO
        "pressure_time_method = CRANK_NICOLSON\n"   // removed: warn and drop
        "temperature_time_method = RK4\n"           // removed: warn and drop
        "dissolved_gas_bunsen_coeff = 0\n");        // removed: warn and drop

    check(cfg.parse_errors.empty(),
          "deliberate legacy branches must not be reported as unknown keys");
    check(cfg.parse_warnings.size() == 3,
          "the three warn-and-drop keys must still warn");
    check(cfg.cavitation_model == CavitationModel::JFO,
          "ELROD_ADAMS must still map to JFO");
    check(cfg.liquid_viscosity_model == LiquidViscosityModel::CONSTANT,
          "the viscosity_model alias must still apply");
}

void test_manifest_matches_parser() {
    for (const std::string& key : SimulationConfig::known_config_keys()) {
        SimulationConfig cfg;
        cfg.load_from_text(key + (key == "dissolved_gas_species" ? " = AIR\n" : " = 1\n"));
        check(cfg.parse_errors.empty(),
              "known_config_keys() lists '" + key + "' but the parser rejects it");
    }
}

void test_written_config_round_trips() {
    SimulationConfig cfg;
    SimulationConfig reloaded;
    reloaded.load_from_text(cfg.to_config_text());
    check(reloaded.parse_errors.empty(),
          "to_config_text() must write only keys the parser accepts: " +
              (reloaded.parse_errors.empty() ? std::string() : reloaded.parse_errors[0]));
}

void test_shipped_configs_parse() {
#ifdef PANCAKE_SOURCE_DIR
    const fs::path root(PANCAKE_SOURCE_DIR);
    int checked = 0;
    const auto check_file = [&](const fs::path& path) {
        if (path.extension() != ".txt") return;
        if (path.filename().string().find("config") == std::string::npos) return;
        std::ifstream file(path);
        if (!file.is_open()) return;
        std::string text((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
        // *_template.txt files carry Python "{placeholder}" values the drivers
        // substitute. Replace them with a parseable number so the key names are
        // still policed -- the templates are exactly where the dropped
        // `cross_film_model` line lived.
        for (size_t open = text.find('{'); open != std::string::npos;
             open = text.find('{', open)) {
            const size_t close = text.find('}', open);
            if (close == std::string::npos) break;
            text.replace(open, close - open + 1, "1");
        }
        SimulationConfig cfg;
        cfg.load_from_text(text);
        check(cfg.parse_errors.empty(),
              path.string() + ": " +
                  (cfg.parse_errors.empty() ? std::string() : cfg.parse_errors[0]));
        ++checked;
    };
    // Shipped starting points at the root, then every validation-case config.
    for (const auto& entry : fs::directory_iterator(root)) {
        if (entry.is_regular_file()) check_file(entry.path());
    }
    if (fs::exists(root / "cases")) {
        for (const auto& entry : fs::recursive_directory_iterator(root / "cases")) {
            if (entry.is_regular_file()) check_file(entry.path());
        }
    }
    check(checked > 0, "no shipped config files were found to check");

    // config.txt is the new user's starting point: it must not merely parse, it
    // must be a case the solver agrees to run. This is what pins the shipped
    // cavitation default against the FBNS scope rules in validate().
    SimulationConfig shipped;
    shipped.load_from_file(root / "config.txt");
    std::vector<std::string> errors, warnings;
    shipped.validate(errors, warnings);
    check(errors.empty(), "config.txt must pass validate(): " +
                              (errors.empty() ? std::string() : errors[0]));
#endif
}

}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_config_keys...\n";

    test_shipped_silent_drop_is_now_fatal();
    test_legacy_branches_still_parse();
    test_manifest_matches_parser();
    test_written_config_round_trips();
    test_shipped_configs_parse();

    if (rank == 0) std::cout << "PASS: test_config_keys\n";
    PetscFinalize();
    return 0;
}
