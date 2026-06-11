// Console harness that exercises the GUI's result-file parser against a real
// solver output directory. Not part of any shipped binary; used during
// development/acceptance checks.
//   results_parser_check.exe <output_dir> [prefix]

#include "../src/gui/results.hpp"

#include <cstdio>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: results_parser_check <output_dir> [prefix]\n");
        return 2;
    }
    const std::filesystem::path output_dir = argv[1];
    const std::string prefix = argc > 2 ? argv[2] : "solution";

    std::string error;
    const auto steps = parse_pvd(output_dir, prefix, error);
    std::printf("steps: %zu  (error: '%s')\n", steps.size(), error.c_str());
    for (const StepInfo& info : steps) {
        std::printf("  step %d t=%g  %s\n", info.step, info.time, info.pvts.string().c_str());
    }
    if (steps.empty()) return 1;

    error.clear();
    const auto grid = parse_step_grid(steps.back(), error);
    if (!grid) {
        std::printf("parse_step_grid FAILED: %s\n", error.c_str());
        return 1;
    }
    std::printf("grid %dx%d step %d t=%g, %zu fields\n", grid->nx, grid->ny, grid->step,
                grid->time, grid->fields.size());
    for (const auto& [name, values] : grid->fields) {
        double lo = 1e300, hi = -1e300;
        for (double v : values) {
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        std::printf("  %-22s [%zu] min %.6g max %.6g\n", name.c_str(), values.size(), lo, hi);
    }
    std::printf("theta_deg[0]=%.3f theta_deg[last]=%.3f z[0]=%.6g z[last]=%.6g\n",
                grid->theta_deg.front(), grid->theta_deg.back(), grid->z_m.front(),
                grid->z_m.back());

    // Malformed-input check: truncated piece file must yield an error, not a crash.
    const std::filesystem::path bad_dir = output_dir / "malformed_check";
    std::filesystem::create_directories(bad_dir / "processor0");
    std::filesystem::copy_file(steps.back().pvts, bad_dir / steps.back().pvts.filename(),
                               std::filesystem::copy_options::overwrite_existing);
    {
        FILE* f = _wfopen((bad_dir / "processor0" /
                           (prefix + "_" + std::to_string(steps.back().step) + "_0.vts"))
                              .c_str(),
                          L"wb");
        std::fputs("<VTKFile><StructuredGrid WholeExtent=\"0 80 0 24 0 0\"><Piece Extent=\"0 80 "
                   "0 24 0 0\"><CellData><DataArray Name=\"p\">1 2 garbage",
                   f);
        std::fclose(f);
    }
    StepInfo bad = steps.back();
    bad.pvts = bad_dir / steps.back().pvts.filename();
    error.clear();
    const auto bad_grid = parse_step_grid(bad, error);
    std::printf("malformed test: grid=%s error='%s'\n", bad_grid ? "NOT-NULL (BAD)" : "null (ok)",
                error.c_str());
    return bad_grid ? 1 : 0;
}
