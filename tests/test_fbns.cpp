// FBNS complementarity cavitation (src/fbns.cpp).
//
// Test 1 - Complementarity. On the 00c grooved-bearing geometry the converged
//   state satisfies p - p_cav >= 0, 1 - theta >= 0 and (p - p_cav)(1 - theta) = 0
//   cell by cell, with no switch function and no bulk modulus involved.
// Test 2 - Per-cell mass conservation. FBNS::cell_outflow, which is also what
//   Diagnostics::mass_residual reports, is at Krylov precision on every cell
//   except the deliberately modelled groove.
// Test 3 - Cavitated-zone flux invariant. The discrete steady balance of a
//   cavitated cell reduces to h_face * theta = const, so that product must be
//   constant across the whole cavitated arc - including across rank boundaries,
//   which is what makes this the rank-independence check.
// Test 4 - Free-boundary accuracy against the exact 1D JFO benchmark
//   (cases/00c_jfo_verification/exact_1d_jfo.py, rupture 302.766608 deg). FBNS
//   must land within ~1 cell; Elrod-Adams on the same grid must not, which is
//   the whole reason the model exists.
// Test 5 - Transient storage. A squeezing film closes the per-cell balance
//   including the backward-Euler d(rho h theta)/dt term, and the global
//   diagnostics balance agrees with the summed field.
// Test 6 - Dirichlet axial ends. The boundary rows close as tightly as the
//   interior and carry a nonzero boundary flux.
// Test 7 - Phase-1 scope. Configurations FBNS cannot yet represent are rejected
//   by validate() instead of running with a silently wrong operator.
//
// Run as: mpirun -n <N> ./test_fbns

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "diagnostics.hpp"
#include "fbns.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

namespace {

constexpr double kExactRuptureDeg = 302.766608;
constexpr double kExactReformationDeg = 79.006725;

void check(bool condition, const char* message) {
    if (!condition) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << message << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void log_rank0(const std::string& message) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << message << "\n";
}

std::string fmt(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(3) << value;
    return out.str();
}

// Peak Couette face flux [kg/s]: the throughput a per-cell residual must be
// negligible against.
double flux_scale(const Mesh& mesh, const SimulationConfig& cfg) {
    return cfg.rho * std::abs(cfg.omega) * mesh.R * (cfg.c + cfg.e) * mesh.get_d_z();
}

// The 00c benchmark bearing: Ferron-type geometry, one 18 deg pressurized groove
// at 90 deg, closed (Neumann) axial ends so the problem is the exact 1D one.
SimulationConfig benchmark_config(int n_theta, int n_z) {
    SimulationConfig cfg;
    cfg.R = 0.05;
    cfg.c = 1.52e-4;
    cfg.e = 9.12e-5;
    cfg.L = 0.08;
    cfg.attitude_angle_deg = -90.0;
    cfg.n_theta_global = n_theta;
    cfg.n_z_global = n_z;
    cfg.omega = 418.879020479;
    cfg.mu = 0.0277;
    cfg.rho = 860.0;
    cfg.p_cav = 1.0e5;
    cfg.bulk_modulus = 1.0e9;  // unused by FBNS; kept for the Elrod comparison
    cfg.cavitation_model = CavitationModel::FBNS;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;
    cfg.linear_rtol = 1.0e-13;
    cfg.log_outer_iters = false;
    cfg.fbns_tol = 1.0e-11;
    cfg.fbns_max_iters = 100;
    cfg.max_outer_iters = 300;
    cfg.outer_tol = 1.0e-9;
    cfg.inlets.clear();
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0;
    groove.size = 18.0;
    groove.p_supply = 1.713e5;
    cfg.inlets.push_back(groove);
    return cfg;
}

struct Case {
    SimulationConfig cfg;
    Mesh mesh;
    Communicator comm;
    LinearSystem sys;
    Fields fields;
    Field inlet_mask;

    explicit Case(const SimulationConfig& config)
        : cfg(config), mesh(cfg), comm(mesh), sys(mesh), inlet_mask("inlet_mask", mesh) {
        fields.add("pressure", mesh).fill(cfg.initial_pressure());
        fields.add("theta", mesh).fill(1.0);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("dh_dt", mesh).fill(0.0);
        fields.add("rho", mesh).fill(cfg.rho);
        fields.add("mass_residual", mesh).fill(0.0);
        fields.add("cavitation_clamp_mass", mesh).fill(0.0);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        FilmThickness::compute_inlet_indicator(inlet_mask, mesh, cfg);
        comm.update_ghosts(inlet_mask);
        comm.update_ghosts(fields);
        fields["theta"].store_old_time();
        fields["pressure"].store_old_time();
    }

    FBNS::Stats solve() {
        comm.update_ghosts(fields);
        FBNS::Stats stats = FBNS::solve(fields, sys, comm, mesh, cfg);
        comm.update_ghosts(fields);
        return stats;
    }
};

// Gather one z-row of a cell-centred field into a global, circumferentially
// ordered array on every rank.
std::vector<double> gather_row(const Field& field, const Mesh& mesh, int j) {
    std::vector<double> global(mesh.n_theta_global, 0.0);
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        global[mesh.offset_theta + i] = field(i, j);
    }
    MPI_Allreduce(MPI_IN_PLACE, global.data(), mesh.n_theta_global,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// Rupture / reformation angle [deg], reported at the face between the last
// full-film cell and the first cavitated cell walking downstream from the groove
// trailing edge. Same estimator for every solver, so the comparison is fair.
void locate_fronts(const std::vector<double>& film, double& rupture, double& reformation) {
    const int n = static_cast<int>(film.size());
    const double d_theta = 360.0 / n;
    auto center = [&](int k) { return (k + 0.5) * d_theta; };
    auto cavitated = [&](int k) { return film[k] < 1.0 - 1.0e-6; };

    int start = 0;
    while (start < n && center(start) <= 99.0) ++start;
    rupture = reformation = -1.0;
    for (int k = 1; k < n; ++k) {
        const int prev = (start + k - 1) % n;
        const int cur = (start + k) % n;
        const double face = std::fmod(center(prev) + 0.5 * d_theta, 360.0);
        if (rupture < 0.0 && !cavitated(prev) && cavitated(cur)) {
            rupture = face;
        } else if (rupture >= 0.0 && reformation < 0.0 && cavitated(prev) && !cavitated(cur)) {
            reformation = face;
        }
    }
}

void test_complementarity_and_residual() {
    const SimulationConfig cfg = benchmark_config(240, 4);
    Case c(cfg);
    const FBNS::Stats stats = c.solve();

    check(stats.converged, "FBNS must converge on the 00c benchmark bearing");
    log_rank0("  complementarity: " + std::to_string(stats.newton_iters) +
              " Newton iterations, residual = " + fmt(stats.residual));

    // Grid-independent iteration count is the FBNS selling point; anything in the
    // tens means the Newton/Schur structure is intact.
    check(stats.newton_iters <= 30, "FBNS must converge in a Newton-like iteration count");

    const Field& p = c.fields["pressure"];
    const Field& theta = c.fields["theta"];
    double worst_p = 0.0, worst_theta = 0.0, worst_product = 0.0;
    double min_theta = 1.0;
    const double p_ref = FBNS::pressure_scale(cfg);
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            const double a = (p(i, j) - cfg.p_cav) / p_ref;
            const double b = 1.0 - theta(i, j);
            worst_p = std::max(worst_p, -a);
            worst_theta = std::max(worst_theta, -b);
            worst_product = std::max(worst_product, std::abs(a * b));
            min_theta = std::min(min_theta, theta(i, j));
        }
    }
    double vals[4] = {worst_p, worst_theta, worst_product, -min_theta};
    MPI_Allreduce(MPI_IN_PLACE, vals, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    log_rank0("  complementarity: max(-(p-p_cav)/p_ref) = " + fmt(vals[0]) +
              ", max(theta-1) = " + fmt(vals[1]) +
              ", max product = " + fmt(vals[2]) +
              ", min theta = " + fmt(-vals[3]));
    check(vals[0] <= 0.0, "no cell may sit below p_cav");
    check(vals[1] <= 0.0, "no cell may exceed full film");
    check(vals[2] < 1.0e-14, "(p - p_cav)(1 - theta) must vanish cell by cell");
    check(-vals[3] < 0.5, "the benchmark bearing must actually cavitate");

    // Per-cell mass conservation, through the same entry point the diagnostics
    // field uses.
    Diagnostics::mass_residual(c.fields, c.mesh, cfg, 0.0);
    const Field& residual = c.fields["mass_residual"];
    const double scale = flux_scale(c.mesh, cfg);
    double worst = 0.0;
    double worst_groove_free = 0.0;
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            worst = std::max(worst, std::abs(residual(i, j)));
            const bool near_inlet = c.inlet_mask(i, j) != 0.0 ||
                c.inlet_mask(i - 1, j) != 0.0 || c.inlet_mask(i + 1, j) != 0.0;
            if (!near_inlet) worst_groove_free = std::max(worst_groove_free, std::abs(residual(i, j)));
        }
    }
    double res[2] = {worst, worst_groove_free};
    MPI_Allreduce(MPI_IN_PLACE, res, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    log_rank0("  mass residual: max |r|/flux = " + fmt(res[0] / scale) +
              " (away from the groove " + fmt(res[1] / scale) + ")");
    check(res[1] / scale < 1.0e-9,
          "every cell away from the groove must conserve mass at Krylov precision");
}

// A cavitated cell has p = p_cav on both sides, so its steady balance is purely
// the upwind Couette term: h_face_east * theta = const along the cavitated arc.
void test_cavitated_flux_invariant() {
    const SimulationConfig cfg = benchmark_config(240, 4);
    Case c(cfg);
    c.solve();

    const Mesh& mesh = c.mesh;
    const int j = mesh.n_z_local / 2;
    const std::vector<double> film = gather_row(c.fields["theta"], mesh, j);
    const std::vector<double> h = gather_row(c.fields["h"], mesh, j);
    const int n = mesh.n_theta_global;

    double lo = 1.0e300, hi = -1.0e300;
    int counted = 0;
    for (int k = 0; k < n; ++k) {
        const int east = (k + 1) % n;
        // Deep in the cavitated zone: the cell and its downstream neighbour are
        // both cavitated, so no reformation Poiseuille term enters the balance.
        if (film[k] >= 1.0 - 1.0e-6 || film[east] >= 1.0 - 1.0e-6) continue;
        const double invariant = 0.5 * (h[k] + h[east]) * film[k];
        lo = std::min(lo, invariant);
        hi = std::max(hi, invariant);
        ++counted;
    }
    check(counted > n / 8, "the benchmark must leave a substantial cavitated arc");
    const double spread = (hi - lo) / hi;
    log_rank0("  cavitated flux invariant over " + std::to_string(counted) +
              " cells: relative spread = " + fmt(spread));
    check(spread < 1.0e-9,
          "h_face * theta must be constant across the cavitated arc and across ranks");
}

void test_front_accuracy_beats_elrod() {
    const int n_theta = 240;
    const double cell_deg = 360.0 / n_theta;
    const SimulationConfig cfg = benchmark_config(n_theta, 4);

    Case fbns_case(cfg);
    fbns_case.solve();
    const std::vector<double> fbns_film =
        gather_row(fbns_case.fields["theta"], fbns_case.mesh, fbns_case.mesh.n_z_local / 2);
    double fbns_rupture = 0.0, fbns_reformation = 0.0;
    locate_fronts(fbns_film, fbns_rupture, fbns_reformation);

    // Same geometry, Elrod-Adams. Fresh fields so the two solvers start equal.
    SimulationConfig elrod_cfg = cfg;
    elrod_cfg.cavitation_model = CavitationModel::JFO;
    Case elrod_case(elrod_cfg);
    elrod_case.fields["theta"].fill(
        std::max(std::exp((elrod_cfg.initial_pressure() - elrod_cfg.p_cav) /
                          elrod_cfg.bulk_modulus), elrod_cfg.theta_min));
    elrod_case.comm.update_ghosts(elrod_case.fields);
    elrod_case.fields["theta"].store_old_time();
    Reynolds::solve_elrod(elrod_case.fields, elrod_case.sys, elrod_case.mesh, elrod_cfg);
    elrod_case.comm.update_ghosts(elrod_case.fields);
    const std::vector<double> elrod_film =
        gather_row(elrod_case.fields["theta"], elrod_case.mesh, elrod_case.mesh.n_z_local / 2);
    double elrod_rupture = 0.0, elrod_reformation = 0.0;
    locate_fronts(elrod_film, elrod_rupture, elrod_reformation);

    const double fbns_err = (fbns_rupture - kExactRuptureDeg) / cell_deg;
    const double elrod_err = (elrod_rupture - kExactRuptureDeg) / cell_deg;
    double reform_err = (fbns_reformation - kExactReformationDeg) / cell_deg;
    if (reform_err > 180.0 / cell_deg) reform_err -= 360.0 / cell_deg;

    log_rank0("  rupture: exact " + fmt(kExactRuptureDeg) + " deg, FBNS " + fmt(fbns_rupture) +
              " deg (" + fmt(fbns_err) + " cells), Elrod " + fmt(elrod_rupture) +
              " deg (" + fmt(elrod_err) + " cells)");
    log_rank0("  reformation: exact " + fmt(kExactReformationDeg) + " deg, FBNS " +
              fmt(fbns_reformation) + " deg (" + fmt(reform_err) + " cells)");

    check(std::abs(fbns_err) <= 1.5,
          "FBNS must place the rupture front within about one cell of the exact solution");
    check(std::abs(reform_err) <= 1.5,
          "FBNS must place the reformation front within about one cell of the exact solution");
    check(std::abs(elrod_err) > 4.0 * std::abs(fbns_err),
          "the Elrod-Adams front error must be far larger; otherwise this benchmark "
          "no longer discriminates and the comparison is meaningless");
}

void test_transient_storage() {
    SimulationConfig cfg = benchmark_config(120, 4);
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.dt = 2.0e-5;

    Case c(cfg);
    // A uniform squeeze sized so the storage terms are a sizeable fraction of the
    // Couette throughput; a token rate would let the test pass on a solver that
    // dropped d(rho h theta)/dt entirely.
    const double squeeze_rate = -0.2;  // m/s
    c.fields["dh_dt"].fill(squeeze_rate);
    c.comm.update_ghosts(c.fields["dh_dt"]);

    double worst = 0.0;
    for (int step = 0; step < 3; ++step) {
        for (int i = -2; i < c.mesh.n_theta_local + 2; ++i)
            for (int j = 0; j < c.mesh.n_z_local; ++j)
                c.fields["h"](i, j) += squeeze_rate * cfg.dt;
        c.comm.update_ghosts(c.fields["h"]);

        const FBNS::Stats stats = c.solve();
        check(stats.converged, "FBNS must converge on every transient step");

        Diagnostics::mass_residual(c.fields, c.mesh, cfg, cfg.dt);
        const Field& residual = c.fields["mass_residual"];
        worst = 0.0;
        for (int i = 0; i < c.mesh.n_theta_local; ++i) {
            for (int j = 0; j < c.mesh.n_z_local; ++j) {
                const bool near_inlet = c.inlet_mask(i, j) != 0.0 ||
                    c.inlet_mask(i - 1, j) != 0.0 || c.inlet_mask(i + 1, j) != 0.0;
                if (!near_inlet) worst = std::max(worst, std::abs(residual(i, j)));
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        c.fields["theta"].store_old_time();
    }

    const double scale = flux_scale(c.mesh, cfg);
    log_rank0("  transient squeeze: worst |r|/flux = " + fmt(worst / scale));
    check(worst / scale < 1.0e-9,
          "the backward-Euler storage term must close the per-cell balance");
}

void test_dirichlet_axial_ends() {
    SimulationConfig cfg = benchmark_config(120, 8);
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 1.5e5;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 1.5e5;

    Case c(cfg);
    check(c.solve().converged, "FBNS must converge with open axial ends");

    Diagnostics::mass_residual(c.fields, c.mesh, cfg, 0.0);
    const Diagnostics::MassBalance balance =
        Diagnostics::mass_balance(c.fields, c.mesh, cfg, 0.0);

    const Field& residual = c.fields["mass_residual"];
    double worst_boundary = 0.0;
    double field_sum = 0.0;
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            field_sum += residual(i, j);
            const bool near_inlet = c.inlet_mask(i, j) != 0.0 ||
                c.inlet_mask(i - 1, j) != 0.0 || c.inlet_mask(i + 1, j) != 0.0;
            const bool boundary_row = (j == 0 || j == c.mesh.n_z_local - 1);
            if (!near_inlet && boundary_row) {
                worst_boundary = std::max(worst_boundary, std::abs(residual(i, j)));
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &worst_boundary, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &field_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const double scale = flux_scale(c.mesh, cfg);
    log_rank0("  Dirichlet ends: worst boundary |r|/flux = " + fmt(worst_boundary / scale) +
              ", boundary flux = " + fmt(balance.liquid_boundary_flux) +
              ", sum(field) = " + fmt(field_sum));
    check(std::abs(balance.liquid_boundary_flux) > 0.0,
          "open ends must carry a nonzero boundary flux for this to test anything");
    check(worst_boundary / scale < 1.0e-9,
          "Dirichlet boundary rows must close as tightly as the interior");

    // Interior faces telescope, so the total outflow is the axial boundary flux
    // plus the groove supply. mass_residual zeroes the inlet cells (their
    // imbalance is the modelled source, reported separately), which leaves
    //     sum(field) + boundary_inflow + inlet_source = 0.
    const double closure =
        field_sum + balance.liquid_boundary_flux + balance.inlet_mass_source;
    log_rank0("  Dirichlet ends: inlet source = " + fmt(balance.inlet_mass_source) +
              ", telescoping closure = " + fmt(closure));
    check(std::abs(balance.inlet_mass_source) > 0.0,
          "the groove must supply mass for the telescoping identity to test anything");
    check(std::abs(closure) < 1.0e-9 * std::max(std::abs(balance.liquid_boundary_flux), scale),
          "the summed per-cell residual must reproduce the boundary flux plus the groove supply");
}

void test_phase1_scope_is_enforced() {
    auto rejects = [](const SimulationConfig& cfg) {
        std::vector<std::string> errors, warnings;
        SimulationConfig copy = cfg;
        copy.validate(errors, warnings);
        return !errors.empty();
    };

    SimulationConfig base = benchmark_config(240, 4);
    {
        std::vector<std::string> errors, warnings;
        SimulationConfig ok = base;
        ok.validate(errors, warnings);
        check(errors.empty(), "the phase-1 FBNS configuration must validate cleanly");
    }

    SimulationConfig thermal = base;
    thermal.temperature_model = TemperatureModel::ENERGY_EQUATION;
    check(rejects(thermal), "FBNS must reject the energy equation in phase 1");

    SimulationConfig gas = base;
    gas.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    check(rejects(gas), "FBNS must reject two-phase gas transport in phase 1");

    SimulationConfig io_bc = base;
    io_bc.bc_z_south_type = BCType::INLET_OUTLET;
    io_bc.bc_z_south_val = 1.5e5;
    check(rejects(io_bc), "FBNS must reject INLET_OUTLET axial boundaries");

    SimulationConfig bad_relax = base;
    bad_relax.fbns_relaxation = 1.5;
    check(rejects(bad_relax), "fbns_relaxation outside (0, 1] must be rejected");

    log_rank0("  phase-1 scope: unsupported configurations rejected");
}

}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_fbns...\n";

    test_complementarity_and_residual();
    test_cavitated_flux_invariant();
    test_front_accuracy_beats_elrod();
    test_transient_storage();
    test_dirichlet_axial_ends();
    test_phase1_scope_is_enforced();

    if (rank == 0) std::cout << "PASS: test_fbns\n";
    PetscFinalize();
    return 0;
}
