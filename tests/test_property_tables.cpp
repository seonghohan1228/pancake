#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "config.hpp"
#include "fluid_properties.hpp"
#include "property_table_io.hpp"

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

void write_file(const std::filesystem::path& p, const std::string& text) {
    std::ofstream out(p, std::ios::binary);
    out << text;
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // --- PropertyTable2D bilinear interpolation (non-uniform axis, clamp) ----
    {
        PropertyTable2D tbl;
        tbl.p = {1.0, 2.0, 4.0};  // non-uniform pressure axis
        tbl.t = {300.0, 400.0};
        tbl.values = {10.0, 20.0,    // p = 1
                      30.0, 40.0,    // p = 2
                      70.0, 80.0};   // p = 4

        check(close(tbl.interpolate(2.0, 400.0, -1.0), 40.0, 1e-12), "2D node-exact");
        check(close(tbl.interpolate(1.0, 300.0, -1.0), 10.0, 1e-12), "2D corner node");
        check(close(tbl.interpolate(1.0, 350.0, -1.0), 15.0, 1e-12), "2D linear in T");
        check(close(tbl.interpolate(3.0, 300.0, -1.0), 50.0, 1e-12), "2D linear in P (non-uniform)");
        check(close(tbl.interpolate(3.0, 350.0, -1.0), 55.0, 1e-12), "2D bilinear interior");
        check(close(tbl.interpolate(0.5, 350.0, -1.0), 15.0, 1e-12), "2D clamp low P");
        check(close(tbl.interpolate(9.0, 300.0, -1.0), 70.0, 1e-12), "2D clamp high P");
        check(close(tbl.interpolate(2.0, 250.0, -1.0), 30.0, 1e-12), "2D clamp low T");

        PropertyTable2D empty;
        check(close(empty.interpolate(1.0, 1.0, 42.0), 42.0, 1e-12), "2D empty -> fallback");

        // isotherm extraction used by saturation_pressure inversion
        check(close(tbl.isotherm_at(2, 300.0), 70.0, 1e-12), "2D isotherm node");
        check(close(tbl.isotherm_at(2, 350.0), 75.0, 1e-12), "2D isotherm interp in T");
    }

    // --- invert_csat: first-crossing, flat tail, below-range, noise ---------
    {
        const std::vector<double> p = {1.0, 2.0, 3.0, 4.0, 5.0};
        const std::vector<double> c = {0.0, 0.10, 0.20, 0.20, 0.20};  // rise then flat
        check(close(invert_csat(p, c, 0.05), 1.5, 1e-12), "invert mid-rise");
        check(close(invert_csat(p, c, 0.10), 2.0, 1e-12), "invert at node");
        check(close(invert_csat(p, c, 0.20), 3.0, 1e-12), "invert onset of flat tail");
        check(close(invert_csat(p, c, 0.30), 5.0, 1e-12), "invert above range -> p_max");
        check(close(invert_csat(p, c, 0.0), 1.0, 1e-12), "invert at/below front -> p_min");

        const std::vector<double> cn = {0.0, 0.10, 0.099, 0.20, 0.20};  // 1-LSB dip
        check(close(invert_csat(p, cn, 0.05), 1.5, 1e-12), "invert robust to noise blip");

        check(invert_csat({1.0}, {0.1}, 0.05) == 0.0, "invert degenerate table -> 0");
    }

    // --- 2-D tables consumed by the property layer (config integration) ------
    {
        SimulationConfig cfg;
        cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
        cfg.solubility_model = SolubilityModel::TABLE;
        cfg.dissolved_gas_max = 1.0;
        cfg.property_reference_temperature = 300.0;
        cfg.solubility_table_2d.p = {1.0e5, 2.0e5, 3.0e5};
        cfg.solubility_table_2d.t = {300.0, 350.0};
        cfg.solubility_table_2d.values = {0.05, 0.02,    // p = 1e5 : (T=300, T=350)
                                          0.10, 0.04,    // p = 2e5
                                          0.10, 0.06};   // p = 3e5

        // equilibrium_dissolved_gas reads the 2-D saturation surface.
        check(close(FluidProperties::equilibrium_dissolved_gas(2.0e5, 300.0, cfg), 0.10, 1e-12),
              "2D solubility at a node");
        check(close(FluidProperties::equilibrium_dissolved_gas(1.5e5, 300.0, cfg), 0.075, 1e-12),
              "2D solubility interp in P");

        // saturation_pressure inverts the local isotherm to the bubble-point onset.
        check(close(cfg.saturation_pressure(0.10, 300.0), 2.0e5, 1e-6), "2D p_sat at plateau onset");
        check(close(cfg.saturation_pressure(0.075, 300.0), 1.5e5, 1e-6), "2D p_sat mid-rise");
        check(close(cfg.saturation_pressure(0.05, 350.0), 2.5e5, 1e-6), "2D p_sat on a hotter isotherm");
    }

    // --- effective_p_cav threshold (default SCALAR_PCAV vs LOCAL_PSAT) -------
    {
        SimulationConfig cfg;
        cfg.p_cav = 3.0e4;  // oil cavitation pressure
        cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
        cfg.solubility_model = SolubilityModel::HENRY;
        cfg.dissolved_gas_henry_coeff = 2.175e-7;
        cfg.dissolved_gas_initial = 0.2175;
        cfg.dissolved_gas_max = 0.5;
        cfg.temperature_initial = 313.15;
        cfg.property_reference_temperature = 313.15;

        // Default reproduces the oil value exactly (decoupled-limit regression).
        check(cfg.effective_p_cav() == cfg.p_cav, "effective_p_cav: SCALAR_PCAV == p_cav");
        // LOCAL_PSAT takes the propane bubble point (0.2175 / 2.175e-7 = 1.0 MPa),
        // i.e. the higher-release species governs onset, not the oil p_cav.
        cfg.cavitation_threshold = CavitationThreshold::LOCAL_PSAT;
        check(close(cfg.effective_p_cav(), 1.0e6, 1.0), "effective_p_cav: LOCAL_PSAT = p_sat (1.0 MPa)");
        // No dissolved gas -> nothing outgasses first -> falls back to p_cav.
        cfg.dissolved_gas_initial = 0.0;
        check(cfg.effective_p_cav() == cfg.p_cav, "effective_p_cav: LOCAL_PSAT no-gas -> p_cav");
    }

    // --- CSV loader (rank 0 only: temp-file I/O) -----------------------------
    if (rank == 0) {
        const auto dir = std::filesystem::temp_directory_path();
        const auto good = dir / "pancake_prop_good.csv";
        write_file(good,
                   "P_MPa,T_C,solubility_pct,viscosity_mm2s\n"
                   "0.1,-60,100.0,0.3472\n"
                   "0.2,-60,100.0,0.3472\n"
                   "0.1,40,5.0,12.0\n"
                   "0.2,40,8.0,10.0\n");

        std::string err;
        PropertyTable2D sol = FluidProperties::load_table_2d_csv(good, "solubility_pct", 0.01, err);
        check(err.empty(), "loader: clean load reports no error");
        check(sol.p.size() == 2 && sol.t.size() == 2, "loader: 2x2 axes");
        check(close(sol.p.front(), 0.1e6, 1.0), "loader: P MPa -> Pa");
        check(close(sol.t.front(), 213.15, 1e-6), "loader: T C -> K");
        check(close(sol.interpolate(0.1e6, 313.15, -1.0), 0.05, 1e-9), "loader: solubility %% -> fraction");

        std::string verr;
        PropertyTable2D vis = FluidProperties::load_table_2d_csv(good, "viscosity_mm2s", 1.0e-6, verr);
        check(verr.empty() && close(vis.interpolate(0.2e6, 313.15, -1.0), 10.0e-6, 1e-12),
              "loader: viscosity mm2/s -> m2/s, column-mapped by header");

        // CRLF + BOM + comment + reordered columns must all be tolerated.
        const auto quirky = dir / "pancake_prop_quirky.csv";
        write_file(quirky,
                   "\xEF\xBB\xBF" "T_C,viscosity_mm2s,P_MPa,solubility_pct\r\n"
                   "# a comment line\r\n"
                   "-60,0.3472,0.1,100.0\r\n"
                   "\r\n"
                   "40,12.0,0.1,5.0\r\n"
                   "-60,0.3472,0.2,100.0\r\n"
                   "40,10.0,0.2,8.0\r\n");
        std::string qerr;
        PropertyTable2D q = FluidProperties::load_table_2d_csv(quirky, "solubility_pct", 0.01, qerr);
        check(qerr.empty(), "loader: BOM/CRLF/comment/reordered tolerated");
        check(close(q.interpolate(0.2e6, 313.15, -1.0), 0.08, 1e-9), "loader: reordered columns value");

        // Malformed: missing required header.
        const auto bad_hdr = dir / "pancake_prop_badhdr.csv";
        write_file(bad_hdr, "P_MPa,T_C,foo\n0.1,-60,1.0\n0.2,-60,1.0\n");
        std::string e_hdr;
        PropertyTable2D th = FluidProperties::load_table_2d_csv(bad_hdr, "solubility_pct", 0.01, e_hdr);
        check(!e_hdr.empty() && th.empty(), "loader: missing header -> error");

        // Malformed: incomplete (ragged) grid.
        const auto ragged = dir / "pancake_prop_ragged.csv";
        write_file(ragged,
                   "P_MPa,T_C,solubility_pct\n0.1,-60,1.0\n0.2,-60,1.0\n0.1,40,5.0\n");
        std::string e_rag;
        PropertyTable2D tr = FluidProperties::load_table_2d_csv(ragged, "solubility_pct", 0.01, e_rag);
        check(!e_rag.empty() && tr.empty(), "loader: ragged grid -> error");

        // Malformed: non-numeric cell.
        const auto nan_cell = dir / "pancake_prop_nan.csv";
        write_file(nan_cell,
                   "P_MPa,T_C,solubility_pct\n0.1,-60,abc\n0.2,-60,1.0\n0.1,40,5.0\n0.2,40,8.0\n");
        std::string e_nan;
        PropertyTable2D tn = FluidProperties::load_table_2d_csv(nan_cell, "solubility_pct", 0.01, e_nan);
        check(!e_nan.empty() && tn.empty(), "loader: non-numeric cell -> error");

        // Missing file.
        std::string e_miss;
        PropertyTable2D tm = FluidProperties::load_table_2d_csv(dir / "pancake_does_not_exist.csv",
                                                                "solubility_pct", 0.01, e_miss);
        check(!e_miss.empty() && tm.empty(), "loader: missing file -> error");

        for (const auto& f : {good, quirky, bad_hdr, ragged, nan_cell}) {
            std::error_code ec;
            std::filesystem::remove(f, ec);
        }

        // --- Real PTSV data file, when available (spot-check vs README) ------
#ifdef PANCAKE_DATA_DIR
        const std::filesystem::path data =
            std::filesystem::path(PANCAKE_DATA_DIR) / "R290_PZ68S" / "r290_pz68s.csv";
        std::error_code ec;
        if (std::filesystem::exists(data, ec)) {
            std::string e_real;
            PropertyTable2D s = FluidProperties::load_table_2d_csv(data, "solubility_pct", 0.01, e_real);
            check(e_real.empty(), "loader: real PTSV data loads");
            check(s.p.size() >= 2 && s.t.size() >= 2, "loader: real PTSV axes present");
            // (1.0 MPa, 40 C) is a grid node = 21.746 % (README spot-check).
            check(close(s.interpolate(1.0e6, 313.15, -1.0), 0.21746, 1e-6),
                  "loader: real PTSV (1MPa,40C) solubility = 21.746%");
            std::string e_v;
            PropertyTable2D v = FluidProperties::load_table_2d_csv(data, "viscosity_mm2s", 1.0e-6, e_v);
            check(e_v.empty() && close(v.interpolate(1.0e6, 313.15, -1.0), 6.843e-6, 1e-9),
                  "loader: real PTSV (1MPa,40C) viscosity = 6.843 mm2/s");
        }
#endif
        std::cout << "test_property_tables: OK\n";
    }

    MPI_Finalize();
    return 0;
}
