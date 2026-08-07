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
// Test 7 - Scope. Configurations FBNS cannot yet represent are rejected by
//   validate() instead of running with a silently wrong operator, and the ones
//   phase 2 added (energy equation, SINGLE_PHASE properties) are accepted --
//   including with the dissolved-gas coefficients left at zero, which is the
//   correct input for a pure-oil thermal run.
// Test 8 - Property coupling is inert when switched off. SINGLE_PHASE with every
//   property coefficient at zero must reproduce the CONSTANT solve bit for bit,
//   so enabling the coupling machinery cannot perturb a run that does not use it.
// Test 9 - Viscosity feedback. Replaying the converged coupled mu field through an
//   isothermal FBNS solve as a prescribed coefficient must reproduce the coupled
//   pressure and film content bit for bit, which is the exact statement that the
//   face mobilities read mu per cell.
// Test 10 - Thermohydrodynamic coupling. On a fine effectively-1D grid, FBNS run
//   inside the segregated Picard loop (flow -> properties -> energy) must land on
//   the same temperature field as Elrod-Adams run through the identical loop,
//   within 0.5 K wherever both hold a full film, and must not shift the front gap
//   the two models already have isothermally. The fronts themselves cannot be
//   asked to agree; see the comment on the test.
// Test 11 - Gaseous cavitation (phase 3). With a resolved free-gas phase the
//   complementarity bound moves to theta <= 1 - alpha_gas. FBNS must satisfy it
//   cell by cell, keep per-cell mass conservation at Krylov precision, and
//   release gas in the same place and order of magnitude as the JFO path. Fronts
//   are logged for the record and asserted against nothing.
// Test 12 - Zero-gas identity. TWO_PHASE with alpha_gas = 0 everywhere must
//   reproduce the phase-2 SINGLE_PHASE solve bit for bit, so the moved bound
//   cannot perturb a run that has no free gas in it.
// Test 13 - The moved ceiling is FBNS-only. FullFilm::bound returns exactly 1 for
//   every other cavitation model whatever alpha_gas holds, which is what makes the
//   Elrod and Gumbel paths byte-identical by construction rather than by luck.
// Test 14 - Cross-film generalized (Dowson) coefficients. With
//   cross_film_thermal_model = POLYNOMIAL the FBNS face mobilities become rho_l F2
//   and 2 F1/F0, i.e. mu(T(y)) inside the film integrals. The coefficients must
//   actually differ from h^3/(12 mu) and h; FBNS and Elrod driven through the
//   identical Picard loop with the identical conjugate-wall closure must land on
//   the same temperature over the shared full-film arc; and per-cell mass
//   conservation must hold including on the fronts, which is what pins the
//   diagnostics operator to the assembled one under BOTH cavitation models.
// Test 15 - Inlet/axial-boundary precedence. A full-length pressure groove
//   overwrites the axial boundary rows at their overlap. Diagnostics must not
//   restore those discarded faces or its inlet source grows as 1/dz.
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
#include "cross_film.hpp"
#include "diagnostics.hpp"
#include "energy.hpp"
#include "equation_of_state.hpp"
#include "fbns.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "fluid_properties.hpp"
#include "full_film.hpp"
#include "gas_transport.hpp"
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

double manufactured_full_length_groove_source(int n_z) {
    SimulationConfig cfg = benchmark_config(120, n_z);
    cfg.omega = 0.0;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = cfg.p_cav;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = cfg.p_cav;

    Case c(cfg);
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            c.fields["pressure"](i, j) = c.inlet_mask(i, j) != 0.0
                ? cfg.inlets.front().p_supply : cfg.p_cav;
        }
    }
    c.comm.update_ghosts(c.fields);
    return Diagnostics::mass_balance(c.fields, c.mesh, cfg, 0.0).inlet_mass_source;
}

void test_inlet_overwrites_axial_boundary() {
    const double source_8 = manufactured_full_length_groove_source(8);
    const double source_16 = manufactured_full_length_groove_source(16);
    const double rel = std::abs(source_16 - source_8) /
        std::max(std::abs(source_8), 1.0e-30);
    log_rank0("  full-length groove: source(nz=8) = " + fmt(source_8) +
              ", source(nz=16) = " + fmt(source_16) + ", rel = " + fmt(rel));
    check(std::abs(source_8) > 0.0,
          "the manufactured groove must carry a nonzero circumferential source");
    check(rel < 1.0e-12,
          "an overwritten axial boundary face must not contaminate the inlet source");
}

void test_scope_is_enforced() {
    auto rejects = [](const SimulationConfig& cfg) {
        std::vector<std::string> errors, warnings;
        SimulationConfig copy = cfg;
        copy.validate(errors, warnings);
        return !errors.empty();
    };
    auto accepts = [](const SimulationConfig& cfg) {
        std::vector<std::string> errors, warnings;
        SimulationConfig copy = cfg;
        copy.validate(errors, warnings);
        if (!errors.empty()) {
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) for (const auto& e : errors) std::cerr << "    unexpected error: " << e << "\n";
        }
        return errors.empty();
    };

    SimulationConfig base = benchmark_config(240, 4);
    check(accepts(base), "the isothermal FBNS configuration must validate cleanly");

    // Phase 2: the energy equation and the single-phase liquid property model are
    // in scope. They reach the FBNS operator only through its per-cell
    // coefficients, which the outer Picard loop refreshes between Newton solves.
    SimulationConfig thermal = base;
    thermal.temperature_model = TemperatureModel::ENERGY_EQUATION;
    check(accepts(thermal), "FBNS must accept the energy equation in phase 2");

    SimulationConfig thd = thermal;
    thd.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
    thd.liquid_viscosity_model = LiquidViscosityModel::EMPIRICAL;
    thd.viscosity_temperature_coeff = 3700.0;
    thd.property_reference_temperature = 313.15;
    // The pure-oil THD case: no gas species participates, so both dissolved-gas
    // coefficients are legitimately zero. Demanding a non-zero Henry slope or a
    // non-zero log-mixing exponent here would be demanding an invented number.
    thd.dissolved_gas_initial = 0.0;
    thd.dissolved_gas_henry_coeff = 0.0;
    thd.solution_viscosity_gas_coeff = 0.0;
    check(accepts(thd),
          "a SINGLE_PHASE pure-oil THD run must validate with both dissolved-gas "
          "coefficients at zero");

    // ... but once a gas species really does participate, a zero Henry slope is a
    // genuine error again: c_sat = 0 would flash the whole charge at once.
    SimulationConfig thd_with_gas = thd;
    thd_with_gas.dissolved_gas_initial = 1.0e-3;
    check(rejects(thd_with_gas),
          "dissolved_gas_henry_coeff = 0 must still be rejected when dissolved gas is present");

    // Phase 3: free gas is in scope. It moves the full-film bound to
    // theta <= 1 - alpha_gas, which is affine in theta with the same derivative,
    // so the Fischer-Burmeister composition and both Jacobian blocks are intact.
    SimulationConfig gas = base;
    gas.solution_mode = SolutionMode::TRANSIENT;
    gas.dt = 1.0e-6;
    gas.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    gas.gas_mass_transfer_rate = 1.0e3;
    gas.dissolved_gas_initial = 1.0e-3;
    gas.dissolved_gas_henry_coeff = 1.2e-9;
    gas.mu_gas = 8.7e-6;
    check(accepts(gas), "FBNS must accept TWO_PHASE gaseous cavitation in phase 3");

    // ... but the void-coupled plateau closure is not. It makes the cavitation
    // datum a per-cell field that only solve_elrod's outer loop advances, so FBNS
    // would run it frozen at its initial value.
    SimulationConfig void_coupled = gas;
    void_coupled.gas_pressure_coupling = GasPressureCoupling::VOID_COUPLED;
    check(rejects(void_coupled),
          "FBNS must reject gas_pressure_coupling = VOID_COUPLED: the per-cell "
          "plateau has no update site outside solve_elrod");

    // The generalized (Dowson) coefficients ARE in the FBNS operator now: the
    // face mobilities read rho_l F2 and the Couette thickness 2 F1/F0 from the
    // same ReynoldsCoeff definition the Elrod path uses, and CrossFilm::update
    // rebuilds them once per outer iteration, which keeps them frozen inside one
    // Newton solve exactly as mu is.
    SimulationConfig cross = thermal;
    cross.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    check(accepts(cross), "FBNS must accept the cross-film generalized coefficients");

    SimulationConfig io_bc = base;
    io_bc.bc_z_south_type = BCType::INLET_OUTLET;
    io_bc.bc_z_south_val = 1.5e5;
    check(rejects(io_bc), "FBNS must reject INLET_OUTLET axial boundaries");

    SimulationConfig bad_relax = base;
    bad_relax.fbns_relaxation = 1.5;
    check(rejects(bad_relax), "fbns_relaxation outside (0, 1] must be rejected");

    log_rank0("  scope: phase-2 configurations accepted, out-of-scope ones rejected");
}

// ---------------------------------------------------------------------------
// Thermohydrodynamic coupling
// ---------------------------------------------------------------------------

// Every field FluidProperties::update_solution_fields and Energy::solve read or
// write, so the coupled driver below sees the same field set main.cpp builds.
struct ThermalCase {
    SimulationConfig cfg;
    Mesh mesh;
    Communicator comm;
    LinearSystem sys;
    Fields fields;

    explicit ThermalCase(const SimulationConfig& config)
        : cfg(config), mesh(cfg), comm(mesh), sys(mesh) {
        const double p0 = cfg.initial_pressure();
        fields.add("pressure", mesh).fill(p0);
        fields.add("theta", mesh).fill(
            cfg.cavitation_model == CavitationModel::FBNS
                ? 1.0
                : std::max(std::exp((p0 - cfg.effective_p_cav()) / cfg.bulk_modulus), cfg.theta_min));
        fields.add("h", mesh).fill(cfg.c);
        fields.add("dh_dt", mesh).fill(0.0);
        fields.add("rho", mesh).fill(cfg.rho);
        fields.add("mu", mesh).fill(cfg.mu);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("mu_liquid_solution", mesh).fill(cfg.mu);
        fields.add("cp_liquid_solution", mesh).fill(cfg.rho_cp / cfg.rho);
        fields.add("rho_cp_liquid_solution", mesh).fill(cfg.rho_cp);
        fields.add("k_liquid_solution", mesh).fill(cfg.thermal_conductivity);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("rho_gas", mesh).fill(0.0);
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("heat_generation", mesh).fill(0.0);
        fields.add("inlet_indicator", mesh).fill(0.0);
        fields.add("mass_residual", mesh).fill(0.0);
        fields.add("cavitation_clamp_mass", mesh).fill(0.0);
        // Allocation is unconditional and inert; only CrossFilm::update is gated,
        // exactly as main.cpp does it.
        CrossFilm::register_fields(fields, mesh, cfg);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        FilmThickness::compute_inlet_indicator(fields["inlet_indicator"], mesh, cfg);
        comm.update_ghosts(fields);
        fields["theta"].store_old_time();
        fields["pressure"].store_old_time();
        fields["temperature"].store_old_time();
    }
};

// Segregated Picard sweep, deliberately mirroring the STEADY_STATE outer loop of
// main.cpp: properties -> cavitation solve -> properties -> energy, repeated
// until the normalized pressure and temperature changes both fall below
// coupling_tolerance. Both cavitation models are driven through this one
// function, so any difference the comparison finds is the models', not the
// driver's. Gas transport, velocities and forces are omitted; they feed back into
// nothing this test measures.
int solve_coupled(ThermalCase& c) {
    Field p_prev("p_prev", c.mesh);
    Field t_prev("t_prev", c.mesh);
    const double p_ref = std::max(c.cfg.property_reference_pressure, 1.0);
    const double t_ref = std::max(c.cfg.temperature_reference, 1.0);

    auto max_change = [&](const Field& cur, const Field& prev, double ref) {
        double local = 0.0;
        for (int i = 0; i < c.mesh.n_theta_local; ++i)
            for (int j = 0; j < c.mesh.n_z_local; ++j)
                local = std::max(local, std::abs(cur(i, j) - prev(i, j)));
        MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return local / ref;
    };

    for (int outer = 0; outer < c.cfg.coupling_max_iters; ++outer) {
        for (int i = 0; i < c.mesh.n_theta_local; ++i) {
            for (int j = 0; j < c.mesh.n_z_local; ++j) {
                p_prev(i, j) = c.fields["pressure"](i, j);
                t_prev(i, j) = c.fields["temperature"](i, j);
            }
        }

        FluidProperties::update_solution_fields(c.fields, c.mesh, c.cfg);
        c.comm.update_ghosts(c.fields);

        // Cross-film profile, conjugate walls and the Dowson F-integrals, before
        // the flow solve and with the same ghost exchange main.cpp performs, so
        // the generalized coefficients this iteration uses are the current ones.
        // No-op unless cross_film_thermal_model is enabled.
        if (CrossFilm::active(c.cfg)) {
            CrossFilm::update(c.fields, c.mesh, c.cfg);
            c.comm.update_ghosts(c.fields["cross_film_f0"]);
            c.comm.update_ghosts(c.fields["cross_film_f1"]);
            c.comm.update_ghosts(c.fields["cross_film_f2"]);
            c.comm.update_ghosts(c.fields["cross_film_couette_h"]);
        }

        if (FBNS::active(c.cfg)) {
            const FBNS::Stats stats = FBNS::solve(c.fields, c.sys, c.comm, c.mesh, c.cfg);
            check(stats.converged, "FBNS must converge on every coupled outer iteration");
        } else {
            Reynolds::solve_elrod(c.fields, c.sys, c.mesh, c.cfg);
        }

        FluidProperties::update_solution_fields(c.fields, c.mesh, c.cfg);
        c.comm.update_ghosts(c.fields);
        Energy::solve(c.fields, c.sys, c.mesh, c.cfg);
        c.comm.update_ghosts(c.fields["temperature"]);

        const double r_p = max_change(c.fields["pressure"], p_prev, p_ref);
        const double r_t = max_change(c.fields["temperature"], t_prev, t_ref);
        if (r_p < c.cfg.coupling_tolerance && r_t < c.cfg.coupling_tolerance) {
            return outer + 1;
        }
    }
    return c.cfg.coupling_max_iters;
}

// Load capacity [N] from the pressure field, for the isothermal/thermal sanity
// direction. Only the magnitude of the gauge pressure integral is needed.
double load_capacity(const ThermalCase& c) {
    const Field& p = c.fields["pressure"];
    const double area = c.mesh.R * c.mesh.get_d_theta() * c.mesh.get_d_z();
    double fx = 0.0, fy = 0.0;
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        const double angle = (c.mesh.offset_theta + i + 0.5) * c.mesh.get_d_theta();
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            const double gauge = (p(i, j) - c.cfg.effective_p_cav()) * area;
            fx -= gauge * std::cos(angle);
            fy -= gauge * std::sin(angle);
        }
    }
    double f[2] = {fx, fy};
    MPI_Allreduce(MPI_IN_PLACE, f, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return std::sqrt(f[0] * f[0] + f[1] * f[1]);
}

// Effectively-1D thermal bearing: the film and both axial boundaries are
// z-independent, so the solution is a pure function of theta and the comparison
// is not contaminated by side-leakage differences.
SimulationConfig thermal_config(int n_theta, CavitationModel model, double ecc_ratio = 0.6) {
    SimulationConfig cfg = benchmark_config(n_theta, 4);
    cfg.e = ecc_ratio * cfg.c;
    cfg.cavitation_model = model;
    cfg.temperature_model = TemperatureModel::ENERGY_EQUATION;
    cfg.temperature_initial = 313.15;
    cfg.temperature_reference = 313.15;
    cfg.journal_wall_temperature = 313.15;
    cfg.bearing_wall_temperature = 313.15;
    cfg.journal_heat_transfer = 1000.0;
    cfg.bearing_heat_transfer = 1000.0;
    cfg.rho_cp = 1.72e6;
    cfg.thermal_conductivity = 0.13;
    // Andrade viscosity in Kelvin: mu = mu_ref exp[E_mu (1/T - 1/T_ref)].
    cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
    cfg.liquid_viscosity_model = LiquidViscosityModel::EMPIRICAL;
    cfg.viscosity_temperature_coeff = 3700.0;
    cfg.property_reference_temperature = 313.15;
    cfg.property_reference_pressure = 1.013e5;
    // Cavitated cells hold only theta*h of liquid: the D14 weighting, which FBNS
    // must honour through exactly the same "theta" field the Elrod path fills.
    cfg.cavitated_film_model = CavitatedFilmModel::STRIATED;
    cfg.coupling_max_iters = 200;
    cfg.coupling_tolerance = 1e-8;
    cfg.inlets[0].feeds_fresh_oil = true;
    cfg.inlets[0].t_supply = 313.15;
    return cfg;
}

// FBNS and Elrod-Adams do NOT place the free boundary in the same cell, and they
// are not supposed to: test_front_accuracy_beats_elrod above measures Elrod
// missing the exact 1-D rupture angle by an order of magnitude more cells than
// FBNS, which is the entire reason this solver exists. That gap is a property of
// Elrod's barotropic smearing, not of the grid or of the coupling: on this
// bearing it stays near 13 deg from n_theta = 240 to 720, and raising
// bulk_modulus to sharpen the switch makes it worse and stops the Elrod outer
// loop converging at all. So a coupled comparison cannot ask the two fronts to
// coincide.
//
// What it can ask -- and what tests the coupling rather than the cavitation
// discretization -- is that turning the energy equation on introduces no NEW
// disagreement:
// wherever both models agree the film is full, their temperatures must agree to
// within 0.5 K. The front positions are reported for the record only.
// The eccentricity is held at 0.3 for this comparison because the residual
// temperature difference is proportional to the thermal rise (measured 0.24 K at
// e/c = 0.2, 0.44 K at 0.3, 0.70 K at 0.4, 1.42 K at 0.6): it is Elrod's
// misplaced front leaking into its pressure field and hence into the Poiseuille
// share of the enthalpy transport, which is a quarter of the total here.
void test_energy_coupling_matches_elrod() {
    const int n_theta = 720;
    const double cell_deg = 360.0 / n_theta;
    const double ecc = 0.3;

    auto isothermal = [&](CavitationModel model) {
        SimulationConfig cfg = thermal_config(n_theta, model, ecc);
        cfg.temperature_model = TemperatureModel::ISOTHERMAL;
        return cfg;
    };

    ThermalCase fbns(thermal_config(n_theta, CavitationModel::FBNS, ecc));
    const int fbns_iters = solve_coupled(fbns);
    ThermalCase elrod(thermal_config(n_theta, CavitationModel::JFO, ecc));
    const int elrod_iters = solve_coupled(elrod);
    ThermalCase fbns_iso(isothermal(CavitationModel::FBNS));
    solve_coupled(fbns_iso);
    ThermalCase elrod_iso(isothermal(CavitationModel::JFO));
    solve_coupled(elrod_iso);

    check(fbns_iters < fbns.cfg.coupling_max_iters,
          "the FBNS thermal Picard loop must reach its fixed point");
    check(elrod_iters < elrod.cfg.coupling_max_iters,
          "the Elrod thermal Picard loop must reach its fixed point");

    const int j = fbns.mesh.n_z_local / 2;
    auto rupture_of = [&](ThermalCase& c) {
        double rup = 0.0, ref = 0.0;
        locate_fronts(gather_row(c.fields["theta"], c.mesh, j), rup, ref);
        return rup;
    };
    const double fbns_shift = rupture_of(fbns) - rupture_of(fbns_iso);
    const double elrod_shift = rupture_of(elrod) - rupture_of(elrod_iso);
    const double gap_thd = rupture_of(fbns) - rupture_of(elrod);
    const double gap_iso = rupture_of(fbns_iso) - rupture_of(elrod_iso);

    // Temperature, compared only where both models hold a full film: in the
    // cavitated arc the two are transporting enthalpy over arcs of different
    // length, which measures Elrod's front error, not the coupling.
    const std::vector<double> t_fbns = gather_row(fbns.fields["temperature"], fbns.mesh, j);
    const std::vector<double> t_elrod = gather_row(elrod.fields["temperature"], elrod.mesh, j);
    const std::vector<double> th_fbns = gather_row(fbns.fields["theta"], fbns.mesh, j);
    const std::vector<double> th_elrod = gather_row(elrod.fields["theta"], elrod.mesh, j);
    double worst_full_film = 0.0, worst_any = 0.0;
    double t_max_fbns = 0.0, t_max_elrod = 0.0;
    int counted = 0;
    for (int k = 0; k < n_theta; ++k) {
        const double diff = std::abs(t_fbns[k] - t_elrod[k]);
        worst_any = std::max(worst_any, diff);
        if (th_fbns[k] >= 1.0 - 1.0e-6 && th_elrod[k] >= 1.0 - 1.0e-6) {
            worst_full_film = std::max(worst_full_film, diff);
            ++counted;
        }
        t_max_fbns = std::max(t_max_fbns, t_fbns[k]);
        t_max_elrod = std::max(t_max_elrod, t_elrod[k]);
    }

    log_rank0("  THD coupling: FBNS " + std::to_string(fbns_iters) + " Picard iters, Elrod " +
              std::to_string(elrod_iters) + " Picard iters");
    log_rank0("  THD coupling: heating moves the rupture front by " + fmt(fbns_shift / cell_deg) +
              " cells (FBNS) and " + fmt(elrod_shift / cell_deg) + " cells (Elrod); the "
              "FBNS-Elrod gap goes " + fmt(gap_iso / cell_deg) + " -> " + fmt(gap_thd / cell_deg) +
              " cells");
    log_rank0("  THD coupling: T_max " + fmt(t_max_fbns) + " K vs " + fmt(t_max_elrod) +
              " K; max |dT| = " + fmt(worst_full_film) + " K over " + std::to_string(counted) +
              " shared full-film cells (" + fmt(worst_any) + " K everywhere)");

    check(t_max_fbns > fbns.cfg.temperature_initial + 1.0,
          "the coupled run must actually heat the film, or the comparison is vacuous");
    check(counted > n_theta / 4,
          "the shared full-film arc must be substantial for the comparison to mean anything");
    // The front positions above are logged, not asserted. Elrod cannot serve as a
    // front reference in either regime -- it sits ~12 deg from the exact solution
    // isothermally, and under heating its front quantizes onto a cell and does not
    // move at all while the FBNS front shifts 2 deg. The exact 1-D JFO benchmark in
    // test_front_accuracy_beats_elrod is the reference that means something, and it
    // is unaffected by anything this test changes.
    check(worst_full_film < 0.5,
          "FBNS and Elrod temperatures must agree within 0.5 K wherever both hold a full film");

    // Direction check: the Andrade law thins the oil, which must cost load capacity.
    const std::vector<double> mu_hot = gather_row(fbns.fields["mu"], fbns.mesh, j);
    double mu_lo = 1.0e300, mu_hi = -1.0e300;
    for (int k = 0; k < n_theta; ++k) {
        mu_lo = std::min(mu_lo, mu_hot[k]);
        mu_hi = std::max(mu_hi, mu_hot[k]);
    }
    const double w_hot = load_capacity(fbns);
    const double w_iso = load_capacity(fbns_iso);
    log_rank0("  THD coupling: mu in [" + fmt(mu_lo) + ", " + fmt(mu_hi) + "] Pa s, W = " +
              fmt(w_hot) + " N thermal vs " + fmt(w_iso) + " N isothermal");
    check(mu_hi / mu_lo > 1.1,
          "the Andrade law must produce a genuinely varying viscosity field");
    check(w_hot < 0.95 * w_iso,
          "thinning the oil must reduce the load capacity relative to the isothermal run");
}

// Same bearing with the cross-film model on. The conjugate shaft/bush closure
// replaces the Newton-cooling walls entirely (SimulationConfig::validate warns
// about this), so journal_heat_transfer / bearing_heat_transfer are zeroed here
// to say so explicitly rather than leave inert values in the config.
SimulationConfig cross_film_config(int n_theta, CavitationModel model, double ecc_ratio) {
    SimulationConfig cfg = thermal_config(n_theta, model, ecc_ratio);
    cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    cfg.cross_film_profile_order = 4;
    cfg.cross_film_quadrature_points = 6;
    cfg.journal_heat_transfer = 0.0;
    cfg.bearing_heat_transfer = 0.0;
    cfg.shaft_ambient_heat_transfer = 0.0;  // globally adiabatic shaft
    cfg.bush_thermal_conductivity = 50.0;
    cfg.bush_outer_radius = 0.10;
    cfg.bush_outer_heat_transfer = 1000.0;
    cfg.thermal_ambient_temperature = 313.15;
    return cfg;
}

// Cross-film generalized (Dowson) coefficients under FBNS.
//
// The point of the model is that mu(T(y)) acts INSIDE the film integrals: the
// Poiseuille mobility becomes rho_l F2 and the Couette thickness 2 F1/F0, not
// h^3/(12 mu_avg) and h. Three things have to hold at once for that to be real
// rather than decorative.
//
//  1. The coefficients must actually differ from their classical values, or
//     every assertion below is vacuous.
//  2. FBNS and Elrod-Adams, driven through the identical Picard loop with the
//     identical cross-film closure, must land on the same temperature wherever
//     both hold a full film -- the same statement
//     test_energy_coupling_matches_elrod makes without the cross-film model, so
//     a regression here is the generalized coefficients disagreeing between the
//     two operators rather than anything about cavitation.
//  3. Per-cell mass conservation must survive, INCLUDING on the rupture and
//     reformation fronts. This is the sharp test of the diagnostics fix: the
//     residual is re-evaluated by Diagnostics from ReynoldsCoeff, so if the
//     diagnostics still used h^3/(12 mu) while the operator used F2 the front
//     cells would light up by the size of the difference measured in (1).
void test_cross_film_generalized_coefficients() {
    const int n_theta = 720;
    const double ecc = 0.3;

    ThermalCase fbns(cross_film_config(n_theta, CavitationModel::FBNS, ecc));
    const int fbns_iters = solve_coupled(fbns);
    ThermalCase elrod(cross_film_config(n_theta, CavitationModel::JFO, ecc));
    const int elrod_iters = solve_coupled(elrod);

    check(fbns_iters < fbns.cfg.coupling_max_iters,
          "the FBNS cross-film Picard loop must reach its fixed point");
    check(elrod_iters < elrod.cfg.coupling_max_iters,
          "the Elrod cross-film Picard loop must reach its fixed point");

    // (1) The generalized coefficients must be genuinely generalized.
    double worst_f2 = 0.0, worst_couette = 0.0, core_rise = 0.0;
    for (int i = 0; i < fbns.mesh.n_theta_local; ++i) {
        for (int j = 0; j < fbns.mesh.n_z_local; ++j) {
            const double h_val = fbns.fields["h"](i, j);
            const double classical = h_val * h_val * h_val / (12.0 * fbns.fields["mu"](i, j));
            worst_f2 = std::max(worst_f2,
                std::abs(fbns.fields["cross_film_f2"](i, j) / classical - 1.0));
            worst_couette = std::max(worst_couette,
                std::abs(fbns.fields["cross_film_couette_h"](i, j) / h_val - 1.0));
            core_rise = std::max(core_rise, std::abs(fbns.fields["cross_film_core_rise"](i, j)));
        }
    }
    double gen[3] = {worst_f2, worst_couette, core_rise};
    MPI_Allreduce(MPI_IN_PLACE, gen, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    log_rank0("  cross-film: max |F2/(h^3/12mu) - 1| = " + fmt(gen[0]) +
              ", max |2F1/F0/h - 1| = " + fmt(gen[1]) +
              ", max core rise = " + fmt(gen[2]) + " K");
    check(gen[0] > 1.0e-3,
          "the cross-film F2 must differ measurably from h^3/(12 mu), or mu(T(y)) is "
          "not reaching the film integrals");
    check(gen[2] > 0.1,
          "the reconstructed profile must carry a real core-to-mean temperature "
          "difference for the comparison to test anything");

    // (2) Same temperature as Elrod on the shared full-film arc.
    const int j = fbns.mesh.n_z_local / 2;
    const std::vector<double> t_fbns = gather_row(fbns.fields["temperature"], fbns.mesh, j);
    const std::vector<double> t_elrod = gather_row(elrod.fields["temperature"], elrod.mesh, j);
    const std::vector<double> th_fbns = gather_row(fbns.fields["theta"], fbns.mesh, j);
    const std::vector<double> th_elrod = gather_row(elrod.fields["theta"], elrod.mesh, j);
    const std::vector<double> tb_fbns = gather_row(fbns.fields["cross_film_t_bush"], fbns.mesh, j);
    const std::vector<double> tb_elrod = gather_row(elrod.fields["cross_film_t_bush"], elrod.mesh, j);
    double worst_full_film = 0.0, worst_bush = 0.0;
    double t_max_fbns = 0.0, t_max_elrod = 0.0;
    int counted = 0;
    for (int k = 0; k < n_theta; ++k) {
        if (th_fbns[k] >= 1.0 - 1.0e-6 && th_elrod[k] >= 1.0 - 1.0e-6) {
            worst_full_film = std::max(worst_full_film, std::abs(t_fbns[k] - t_elrod[k]));
            worst_bush = std::max(worst_bush, std::abs(tb_fbns[k] - tb_elrod[k]));
            ++counted;
        }
        t_max_fbns = std::max(t_max_fbns, t_fbns[k]);
        t_max_elrod = std::max(t_max_elrod, t_elrod[k]);
    }
    log_rank0("  cross-film: " + std::to_string(fbns_iters) + " Picard iters (FBNS) vs " +
              std::to_string(elrod_iters) + " (Elrod); T_max " + fmt(t_max_fbns) + " K vs " +
              fmt(t_max_elrod) + " K; max |dT| = " + fmt(worst_full_film) + " K, max |dT_bush| = " +
              fmt(worst_bush) + " K over " + std::to_string(counted) + " shared full-film cells");
    check(t_max_fbns > fbns.cfg.temperature_initial + 1.0,
          "the cross-film coupled run must actually heat the film");
    check(counted > n_theta / 4,
          "the shared full-film arc must be substantial for the comparison to mean anything");
    check(worst_full_film < 0.5,
          "FBNS and Elrod cross-film temperatures must agree within 0.5 K wherever both "
          "hold a full film");
    // The bush surface temperature is affine in the depth-averaged temperature
    // through the conjugate closure, so its residual disagreement must be
    // EXPLAINED by the depth-averaged one rather than be an independent error of
    // the wall solve. Measured ratio is 1.05; the bound is loose enough not to be
    // a tolerance to tune and tight enough to catch a wall closure that reacted to
    // the cavitation model on its own.
    check(worst_bush < 2.0 * worst_full_film,
          "the conjugate bush surface temperature must not disagree by more than the "
          "depth-averaged temperature it is derived from");

    // (3) Per-cell mass conservation under the generalized operator, fronts
    // included. Only the groove and its immediate neighbours are excluded: their
    // imbalance is the modelled supply, which the global balance reports
    // separately.
    auto residual_report = [&](ThermalCase& c, const char* label, double tol) {
        Diagnostics::mass_residual(c.fields, c.mesh, c.cfg, 0.0);
        const Field& r = c.fields["mass_residual"];
        const Field& mask = c.fields["inlet_indicator"];
        const Field& film = c.fields["theta"];
        double worst = 0.0, worst_front = 0.0;
        double fronts = 0.0;
        for (int ii = 0; ii < c.mesh.n_theta_local; ++ii) {
            for (int jj = 0; jj < c.mesh.n_z_local; ++jj) {
                if (mask(ii, jj) != 0.0 || mask(ii - 1, jj) != 0.0 || mask(ii + 1, jj) != 0.0)
                    continue;
                worst = std::max(worst, std::abs(r(ii, jj)));
                const bool full = film(ii, jj) >= 1.0 - 1.0e-6;
                const bool front = (film(ii - 1, jj) >= 1.0 - 1.0e-6) != full ||
                                   (film(ii + 1, jj) >= 1.0 - 1.0e-6) != full;
                if (front) {
                    worst_front = std::max(worst_front, std::abs(r(ii, jj)));
                    fronts += 1.0;
                }
            }
        }
        double vals[3] = {worst, worst_front, fronts};
        MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, vals + 2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const double scale = flux_scale(c.mesh, c.cfg);
        log_rank0(std::string("  cross-film ") + label + ": worst |r|/flux = " +
                  fmt(vals[0] / scale) + " over all non-groove cells, " +
                  fmt(vals[1] / scale) + " on the " + std::to_string((int)vals[2]) +
                  " front cells");
        check(vals[2] > 0.0, "the comparison needs actual rupture/reformation fronts");
        check(vals[0] / scale < tol,
              "the generalized operator must conserve mass cell by cell, fronts included");
    };
    residual_report(fbns, "FBNS", 1.0e-9);
    // Elrod's barotropic residual is measured against the same Couette throughput
    // for comparability; it is looser because the theta-formulation's Poiseuille
    // conductance exceeds the throughput by O(beta/p), so a round-off-level row
    // residual is a larger fraction of the flux.
    residual_report(elrod, "Elrod", 1.0e-6);
}

// The viscosity feedback has to reach the FBNS *face mobilities*, not merely the
// energy equation, and it has to be the same mu the property update produced --
// the whole frozen-coefficient Picard argument rests on that. Take the converged
// coupled mu field, hand it straight to an isothermal FBNS solve as a prescribed
// coefficient, and require the pressure and film content to come back. If any
// cell silently fell back to the scalar cfg.mu the two differ by the ~40% the
// Andrade law moved mu; agreement is instead at the level of the outer coupling
// tolerance, because the mu read out after the last energy solve is one Picard
// step ahead of the one the last FBNS solve actually consumed.
void test_viscosity_reaches_the_operator() {
    const int n_theta = 240;
    ThermalCase coupled(thermal_config(n_theta, CavitationModel::FBNS));
    const int iters = solve_coupled(coupled);
    check(iters < coupled.cfg.coupling_max_iters, "the coupled solve must converge");

    // Same geometry, no property model at all: mu is whatever we put in it.
    SimulationConfig frozen_cfg = benchmark_config(n_theta, 4);
    frozen_cfg.cavitated_film_model = CavitatedFilmModel::STRIATED;
    Case frozen(frozen_cfg);
    frozen.fields.add("mu", frozen.mesh).fill(frozen_cfg.mu);
    for (int i = 0; i < frozen.mesh.n_theta_local; ++i)
        for (int j = 0; j < frozen.mesh.n_z_local; ++j)
            frozen.fields["mu"](i, j) = coupled.fields["mu"](i, j);
    frozen.solve();

    double worst_p = 0.0, worst_theta = 0.0, mu_spread = 0.0;
    for (int i = 0; i < frozen.mesh.n_theta_local; ++i) {
        for (int j = 0; j < frozen.mesh.n_z_local; ++j) {
            worst_p = std::max(worst_p,
                std::abs(frozen.fields["pressure"](i, j) - coupled.fields["pressure"](i, j)));
            worst_theta = std::max(worst_theta,
                std::abs(frozen.fields["theta"](i, j) - coupled.fields["theta"](i, j)));
            mu_spread = std::max(mu_spread, std::abs(coupled.fields["mu"](i, j) - frozen_cfg.mu));
        }
    }
    double vals[3] = {worst_p, worst_theta, mu_spread};
    MPI_Allreduce(MPI_IN_PLACE, vals, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    const double p_ref = FBNS::pressure_scale(frozen_cfg);
    log_rank0("  frozen-mu replay: max |dp|/p_ref = " + fmt(vals[0] / p_ref) +
              ", max |dtheta| = " + fmt(vals[1]) + ", max |mu - mu_ref|/mu_ref = " +
              fmt(vals[2] / frozen_cfg.mu));
    check(vals[2] > 0.1 * frozen_cfg.mu,
          "the coupled viscosity field must differ substantially from the reference value, "
          "or replaying it proves nothing");
    check(vals[0] / p_ref < 1.0e-6 && vals[1] < 1.0e-6,
          "replaying the converged mu field through FBNS must reproduce the coupled solution: "
          "the face mobilities read that field per cell, not the scalar cfg.mu");
}

// Enabling the property machinery must not perturb a run that does not use it:
// SINGLE_PHASE with every coefficient at zero is the CONSTANT model, and the FBNS
// solve must return the same bits.
void test_property_coupling_is_inert_when_off() {
    SimulationConfig constant_cfg = benchmark_config(240, 4);
    SimulationConfig single_cfg = constant_cfg;
    single_cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
    single_cfg.liquid_viscosity_model = LiquidViscosityModel::EMPIRICAL;
    single_cfg.density_model = DensityModel::CONSTANT;
    single_cfg.viscosity_temperature_coeff = 0.0;
    single_cfg.viscosity_pressure_coeff = 0.0;
    single_cfg.solution_viscosity_gas_coeff = 0.0;
    single_cfg.dissolved_gas_initial = 0.0;
    single_cfg.dissolved_gas_henry_coeff = 0.0;

    ThermalCase a(constant_cfg);
    ThermalCase b(single_cfg);
    solve_coupled(a);
    solve_coupled(b);

    double worst_p = 0.0, worst_theta = 0.0;
    for (int i = 0; i < a.mesh.n_theta_local; ++i) {
        for (int j = 0; j < a.mesh.n_z_local; ++j) {
            worst_p = std::max(worst_p,
                std::abs(a.fields["pressure"](i, j) - b.fields["pressure"](i, j)));
            worst_theta = std::max(worst_theta,
                std::abs(a.fields["theta"](i, j) - b.fields["theta"](i, j)));
        }
    }
    double vals[2] = {worst_p, worst_theta};
    MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    log_rank0("  inert coupling: max |dp| = " + fmt(vals[0]) + " Pa, max |dtheta| = " + fmt(vals[1]));
    check(vals[0] == 0.0 && vals[1] == 0.0,
          "SINGLE_PHASE with zero coefficients must reproduce CONSTANT bit for bit");
}

}  // namespace

// --- Gaseous cavitation (phase 3) -------------------------------------------

// The 00c bearing with a saturated propane charge: the groove feed sits exactly
// at its bubble point, so every cell below the supply pressure is supersaturated
// and outgasses. Transient, because gas transport and release need a time step.
SimulationConfig outgassing_config(int n_theta, int n_z, CavitationModel model) {
    SimulationConfig cfg = benchmark_config(n_theta, n_z);
    cfg.cavitation_model = model;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.dt = 1.0e-5;
    cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    cfg.dissolved_gas_species = DissolvedGasSpecies::PROPANE;
    cfg.solubility_model = SolubilityModel::HENRY;
    cfg.dissolved_gas_initial = 1.0e-3;
    cfg.dissolved_gas_max = 5.0e-3;
    // Henry slope calibrated so the groove feed is exactly saturated: any pressure
    // drop below p_supply releases gas, nothing releases above it.
    cfg.dissolved_gas_henry_coeff = cfg.dissolved_gas_initial / 1.713e5;
    cfg.gas_mass_transfer_rate = 1.0e4;
    cfg.gas_alpha_max = 0.6;
    cfg.mixture_viscosity_model = MixtureViscosityModel::MCADAMS;
    cfg.mu_gas = 8.7e-6;
    cfg.gas_pressure_floor = 1.0e5;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;
    cfg.temperature_initial = 313.15;
    return cfg;
}

// One bearing carrying the full gas field set, advanced through the same
// segregated sequence main.cpp runs: properties -> flow -> properties -> gas
// transport -> gas state. Either cavitation model can drive it.
struct GasCase {
    SimulationConfig cfg;
    Mesh mesh;
    Communicator comm;
    LinearSystem sys;
    Fields fields;
    Field inlet_mask;
    int newton_iters = 0;
    bool converged = true;

    explicit GasCase(const SimulationConfig& config)
        : cfg(config), mesh(cfg), comm(mesh), sys(mesh), inlet_mask("inlet_mask", mesh) {
        fields.add("pressure", mesh).fill(cfg.initial_pressure());
        fields.add("theta", mesh).fill(FBNS::active(cfg)
            ? 1.0
            : std::max(EOS::theta_from_pressure(cfg.initial_pressure(), cfg.p_cav, cfg.bulk_modulus),
                       cfg.theta_min));
        fields.add("h", mesh).fill(cfg.c);
        fields.add("dh_dt", mesh).fill(0.0);
        fields.add("rho", mesh).fill(cfg.rho);
        fields.add("mu", mesh).fill(cfg.mu);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("mu_liquid_solution", mesh).fill(cfg.mu);
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("gas_mass_transfer", mesh).fill(0.0);
        fields.add("rho_gas", mesh).fill(0.0);
        fields.add("gas_clamp_mass", mesh).fill(0.0);
        fields.add("cavitation_clamp_mass", mesh).fill(0.0);
        fields.add("mass_residual", mesh).fill(0.0);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        FilmThickness::compute_inlet_indicator(inlet_mask, mesh, cfg);
        fields.add("inlet_indicator", mesh);
        for (int i = -inlet_mask.n_ghost; i < mesh.n_theta_local + inlet_mask.n_ghost; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                fields["inlet_indicator"](i, j) = inlet_mask(i, j);
        comm.update_ghosts(inlet_mask);
        comm.update_ghosts(fields);
        for (const char* name : {"pressure", "theta", "h", "dissolved_gas", "free_gas_mass"})
            fields[name].store_old_time();
    }

    // Establish the developed, cavitated film before any gas is released, so the
    // transient steps below exercise the moved bound on a bearing that actually
    // has a rupture front rather than on the first few microseconds of spin-up.
    void prime_steady() {
        SimulationConfig steady = cfg;
        steady.solution_mode = SolutionMode::STEADY_STATE;
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);
        if (FBNS::active(cfg)) FBNS::solve(fields, sys, comm, mesh, steady);
        else Reynolds::solve_elrod(fields, sys, mesh, steady);
        comm.update_ghosts(fields);
        for (const char* name : {"pressure", "theta", "h", "dissolved_gas", "free_gas_mass"})
            fields[name].store_old_time();
    }

    // Flow solve only, leaving the gas state exactly as the previous step left it.
    // The complementarity assertions have to be made here: after the gas update
    // alpha_gas has moved and theta has not, so the pair is one lag out of phase
    // by construction (that lag is what makes alpha_gas a coefficient).
    void solve_flow() {
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);
        if (FBNS::active(cfg)) {
            const FBNS::Stats stats = FBNS::solve(fields, sys, comm, mesh, cfg);
            newton_iters = stats.newton_iters;
            converged = stats.converged;
        } else {
            const Reynolds::ElrodStats stats = Reynolds::solve_elrod(fields, sys, mesh, cfg);
            newton_iters = stats.outer_iters;
            converged = stats.converged;
        }
        comm.update_ghosts(fields);
    }

    void advance_gas() {
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);
        GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);
        FluidProperties::update_gas_state(fields, mesh, cfg, cfg.dt);
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        comm.update_ghosts(fields);
        for (const char* name : {"pressure", "theta", "h", "dissolved_gas", "free_gas_mass"})
            fields[name].store_old_time();
    }

    double global_max_of(const char* name) const {
        double local = -1.0e300;
        for (int i = 0; i < mesh.n_theta_local; ++i)
            for (int j = 0; j < mesh.n_z_local; ++j)
                local = std::max(local, fields[name](i, j));
        MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return local;
    }
};

// Test 11 - Gaseous cavitation under FBNS. The invariants asserted are physical
// ones the solver owns: complementarity against the *moved* bound
// theta <= 1 - alpha_gas, the film content staying inside [theta_min, bound], and
// per-cell mass conservation at Krylov precision. The JFO run beside it is a
// qualitative reference only -- both must release gas in the same place, and the
// fronts are logged for the record, but nothing is asserted about their equality.
// Elrod's front on this bearing sits 12-18 deg from the exact 1-D solution (7.3.5),
// so agreeing with it would be evidence of a defect, not of correctness.
void test_gaseous_cavitation() {
    const int n_theta = 240, n_z = 4;
    const int n_steps = 5;
    GasCase fbns(outgassing_config(n_theta, n_z, CavitationModel::FBNS));
    GasCase jfo(outgassing_config(n_theta, n_z, CavitationModel::JFO));
    fbns.prime_steady();
    jfo.prime_steady();

    const double p_ref = FBNS::pressure_scale(fbns.cfg);
    const double scale = flux_scale(fbns.mesh, fbns.cfg);
    double worst_below_pcav = 0.0, worst_above_bound = 0.0;
    double worst_product = 0.0, worst_residual = 0.0;
    int newton_max = 0;

    for (int step = 0; step < n_steps; ++step) {
        fbns.solve_flow();
        check(fbns.converged, "FBNS must converge on every gaseous-cavitation step");
        newton_max = std::max(newton_max, fbns.newton_iters);

        // Complementarity against the gas-adjusted bound, in the same state the
        // solver left, and per-cell conservation through the same coefficients.
        Diagnostics::mass_residual(fbns.fields, fbns.mesh, fbns.cfg, fbns.cfg.dt);
        const Field& p = fbns.fields["pressure"];
        const Field& theta = fbns.fields["theta"];
        const Field& residual = fbns.fields["mass_residual"];
        for (int i = 0; i < fbns.mesh.n_theta_local; ++i) {
            for (int j = 0; j < fbns.mesh.n_z_local; ++j) {
                const double bound = FullFilm::bound(fbns.fields, fbns.cfg, i, j);
                const double a = (p(i, j) - fbns.cfg.effective_p_cav()) / p_ref;
                const double b = bound - theta(i, j);
                worst_below_pcav = std::max(worst_below_pcav, -a);
                worst_above_bound = std::max(worst_above_bound, -b);
                worst_product = std::max(worst_product, std::abs(a * b));
                const bool near_inlet = fbns.inlet_mask(i, j) != 0.0 ||
                    fbns.inlet_mask(i - 1, j) != 0.0 || fbns.inlet_mask(i + 1, j) != 0.0;
                if (!near_inlet) worst_residual = std::max(worst_residual, std::abs(residual(i, j)));
            }
        }
        fbns.advance_gas();

        jfo.solve_flow();
        jfo.advance_gas();
    }

    double vals[4] = {worst_below_pcav, worst_above_bound, worst_product, worst_residual};
    MPI_Allreduce(MPI_IN_PLACE, vals, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double alpha_fbns = fbns.global_max_of("alpha_gas");
    const double alpha_jfo = jfo.global_max_of("alpha_gas");
    const double gas_fbns = fbns.global_max_of("free_gas_mass");
    const double gas_jfo = jfo.global_max_of("free_gas_mass");

    log_rank0("  gaseous FBNS: max Newton = " + std::to_string(newton_max) +
              ", max alpha_gas = " + fmt(alpha_fbns) + " (JFO " + fmt(alpha_jfo) + ")");
    log_rank0("  gaseous FBNS: max free_gas_mass = " + fmt(gas_fbns) +
              " kg/m^2 (JFO " + fmt(gas_jfo) + ")");
    log_rank0("  gaseous FBNS: max(-(p-p_cav)/p_ref) = " + fmt(vals[0]) +
              ", max(theta-(1-alpha_gas)) = " + fmt(vals[1]) +
              ", max product = " + fmt(vals[2]));
    log_rank0("  gaseous FBNS: max |r|/flux away from the groove = " + fmt(vals[3] / scale));

    // Both models must actually release gas, in the same order of magnitude. This
    // is the qualitative agreement the two paths owe each other; the fields
    // themselves are not comparable because Elrod's gas is volumetrically inert.
    check(alpha_fbns > 1.0e-3, "FBNS must release free gas on a supersaturated feed");
    check(alpha_jfo > 1.0e-3, "JFO must release free gas on the same feed");
    check(alpha_fbns < 10.0 * alpha_jfo && alpha_jfo < 10.0 * alpha_fbns,
          "the two models must release gas within an order of magnitude of each other");

    // The invariants that are actually asserted.
    check(vals[0] <= 0.0, "no cell may sit below p_cav");
    check(vals[1] <= 1.0e-12, "no cell may hold more liquid than 1 - alpha_gas");
    check(vals[2] < 1.0e-12, "(p - p_cav)((1 - alpha_gas) - theta) must vanish cell by cell");
    check(vals[3] / scale < 1.0e-6, "per-cell mass conservation must hold at Krylov precision");

    // Fronts: logged, never asserted against each other.
    std::vector<double> film_fbns = gather_row(fbns.fields["theta"], fbns.mesh, n_z / 2);
    std::vector<double> film_jfo = gather_row(jfo.fields["theta"], jfo.mesh, n_z / 2);
    std::vector<double> alpha_row = gather_row(fbns.fields["alpha_gas"], fbns.mesh, n_z / 2);
    // locate_fronts asks "is this cell below full film", so the FBNS row has to be
    // normalized by its own moved ceiling first. Feeding it raw theta would call a
    // gassy but unbroken cell cavitated -- the exact confusion the shared bound
    // exists to prevent.
    int cavitated_cells = 0;
    for (size_t k = 0; k < film_fbns.size(); ++k) {
        film_fbns[k] /= std::max(1.0 - alpha_row[k], 1.0e-30);
        if (film_fbns[k] < 1.0 - 1.0e-6) ++cavitated_cells;
    }
    double rupture_fbns = 0.0, reform_fbns = 0.0, rupture_jfo = 0.0, reform_jfo = 0.0;
    locate_fronts(film_fbns, rupture_fbns, reform_fbns);
    locate_fronts(film_jfo, rupture_jfo, reform_jfo);
    // locate_fronts reports the first transition pair walking downstream from the
    // groove, so a single re-flooded cell just past rupture shows up as an early
    // reformation; the cavitated-cell count is the reliable extent measure.
    log_rank0("  gaseous fronts (logged, not asserted): FBNS rupture " +
              std::to_string(rupture_fbns) + " deg reformation " + std::to_string(reform_fbns) +
              " deg over " + std::to_string(cavitated_cells) + " cavitated cells | JFO rupture " +
              std::to_string(rupture_jfo) + " deg reformation " + std::to_string(reform_jfo) + " deg");
}

// Test 12 - The gas path is inert when there is no gas. TWO_PHASE with
// alpha_gas = 0 everywhere must reproduce the phase-2 (SINGLE_PHASE) FBNS solve
// bit for bit, over a cold start and two warm ones, so admitting the moved bound
// cannot perturb a run that has no free gas in it.
void test_zero_gas_reproduces_phase_two() {
    SimulationConfig single = benchmark_config(240, 4);
    single.solution_mode = SolutionMode::TRANSIENT;
    single.dt = 1.0e-5;
    single.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;

    SimulationConfig two_phase = single;
    two_phase.fluid_property_model = FluidPropertyModel::TWO_PHASE;

    Case a(single);
    Case b(two_phase);
    // Both carry an explicit zero free-gas field, so the TWO_PHASE run really does
    // walk the clamp in FullFilm::gas_volume_fraction rather than short-circuit on
    // a missing field.
    a.fields.add("alpha_gas", a.mesh).fill(0.0);
    b.fields.add("alpha_gas", b.mesh).fill(0.0);

    for (int step = 0; step < 3; ++step) {
        a.solve();
        b.solve();
        int mismatches = 0;
        for (int i = 0; i < a.mesh.n_theta_local; ++i) {
            for (int j = 0; j < a.mesh.n_z_local; ++j) {
                if (a.fields["pressure"](i, j) != b.fields["pressure"](i, j)) ++mismatches;
                if (a.fields["theta"](i, j) != b.fields["theta"](i, j)) ++mismatches;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &mismatches, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        check(mismatches == 0,
              "TWO_PHASE with alpha_gas = 0 must reproduce the SINGLE_PHASE FBNS solve bit for bit");
        a.fields["theta"].store_old_time();
        a.fields["pressure"].store_old_time();
        b.fields["theta"].store_old_time();
        b.fields["pressure"].store_old_time();
    }
    log_rank0("  zero-gas identity: 3 steps, pressure and theta bit-identical to phase 2");
}

// Test 13 - The Elrod and Gumbel paths cannot see the moved ceiling. This is the
// by-construction half of the byte-identity acceptance: FullFilm::bound is the
// only place the ceiling lives, and it returns exactly 1 for every model except
// FBNS + TWO_PHASE, whatever the free-gas field holds.
void test_ceiling_is_fbns_only() {
    SimulationConfig cfg = benchmark_config(16, 4);
    Mesh mesh(cfg);
    Fields fields;
    fields.add("alpha_gas", mesh).fill(0.25);

    const CavitationModel models[] = {CavitationModel::FULL_SOMMERFELD, CavitationModel::GUMBEL,
                                      CavitationModel::JFO, CavitationModel::FBNS};
    const FluidPropertyModel properties[] = {FluidPropertyModel::CONSTANT,
                                             FluidPropertyModel::SINGLE_PHASE,
                                             FluidPropertyModel::TWO_PHASE};
    const GasPressureCoupling couplings[] = {GasPressureCoupling::NONE,
                                             GasPressureCoupling::VOID_COUPLED};
    for (CavitationModel model : models) {
        for (FluidPropertyModel property : properties) {
            for (GasPressureCoupling coupling : couplings) {
                SimulationConfig probe = cfg;
                probe.cavitation_model = model;
                probe.fluid_property_model = property;
                probe.gas_pressure_coupling = coupling;
                const double bound = FullFilm::bound(fields, probe, 0, 0);
                const bool moved = model == CavitationModel::FBNS &&
                    property == FluidPropertyModel::TWO_PHASE;
                check(bound == (moved ? 0.75 : 1.0),
                      "FullFilm::bound may only move under FBNS + TWO_PHASE");
            }
        }
    }
    log_rank0("  ceiling: 1.0 exactly for every model but FBNS + TWO_PHASE (24 combinations)");
}

// The Ferron polynomial viscosity law under both energy couplings.
//
// benchmark_config is already the Ferron, Frene & Boncompain (1983) bearing
// (R = 50 mm, L = 80 mm, c = 152 um, mu = 0.0277 Pa s, rho = 860 kg/m3), so the
// paper's own coefficients belong on it. The law is
//   mu(T) = mu (k0 - k1 Tbar + k2 Tbar^2),  Tbar = T[degC]/T_ref[degC],
// with Tbar a CELSIUS ratio: k0 - k1 + k2 = 1.000 at the 40 degC reference.
//
// Three things have to hold, and each is a separate assertion below:
//   1. the polynomial actually moves mu (otherwise the rest proves nothing);
//   2. that mu field reaches the Reynolds operator under FBNS AND under Elrod,
//      shown by the pressure moving away from the frozen-mu solution in both;
//   3. the polynomial and the Andrade exponential are genuinely different laws,
//      so their converged fields differ.
SimulationConfig ferron_polynomial_config(int n_theta, CavitationModel model) {
    SimulationConfig cfg = thermal_config(n_theta, model);
    cfg.liquid_viscosity_model = LiquidViscosityModel::POLYNOMIAL;
    cfg.viscosity_temperature_coeff = 0.0;  // Andrade off; the polynomial replaces it
    cfg.viscosity_poly_k0 = 3.287;
    cfg.viscosity_poly_k1 = 3.064;
    cfg.viscosity_poly_k2 = 0.777;
    cfg.property_reference_temperature = 313.15;  // 40 degC
    return cfg;
}

void test_polynomial_viscosity_couples() {
    const int n_theta = 240;
    const double mu_ref = 0.0277;
    const double p_ref =
        FBNS::pressure_scale(ferron_polynomial_config(n_theta, CavitationModel::FBNS));

    // Per cavitation model: how far the polynomial moves mu, and how far the
    // resulting pressure moves away from the same run with mu frozen at mu_ref.
    auto spread_and_gap = [&](CavitationModel model, double& mu_spread, double& dp) {
        ThermalCase poly(ferron_polynomial_config(n_theta, model));
        const int poly_iters = solve_coupled(poly);
        check(poly_iters < poly.cfg.coupling_max_iters,
              "the POLYNOMIAL coupled solve must reach its fixed point");

        SimulationConfig frozen_cfg = ferron_polynomial_config(n_theta, model);
        frozen_cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
        ThermalCase frozen(frozen_cfg);
        solve_coupled(frozen);

        mu_spread = 0.0;
        dp = 0.0;
        for (int i = 0; i < poly.mesh.n_theta_local; ++i) {
            for (int j = 0; j < poly.mesh.n_z_local; ++j) {
                mu_spread = std::max(mu_spread, std::abs(poly.fields["mu"](i, j) - poly.cfg.mu));
                dp = std::max(dp,
                    std::abs(poly.fields["pressure"](i, j) - frozen.fields["pressure"](i, j)));
            }
        }
        double vals[2] = {mu_spread, dp};
        MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        mu_spread = vals[0];
        dp = vals[1];
    };

    double fbns_mu = 0.0, fbns_dp = 0.0, elrod_mu = 0.0, elrod_dp = 0.0;
    spread_and_gap(CavitationModel::FBNS, fbns_mu, fbns_dp);
    spread_and_gap(CavitationModel::JFO, elrod_mu, elrod_dp);
    log_rank0("  Ferron polynomial (FBNS):  max |mu - mu_ref|/mu_ref = " + fmt(fbns_mu / mu_ref) +
              ", max |dp|/p_ref vs frozen mu = " + fmt(fbns_dp / p_ref));
    log_rank0("  Ferron polynomial (Elrod): max |mu - mu_ref|/mu_ref = " + fmt(elrod_mu / mu_ref) +
              ", max |dp|/p_ref vs frozen mu = " + fmt(elrod_dp / p_ref));

    check(fbns_mu > 0.05 * mu_ref && elrod_mu > 0.05 * mu_ref,
          "the polynomial must move mu appreciably under both cavitation models");
    check(fbns_dp > 1.0e-3 * p_ref,
          "the polynomial mu field must reach the FBNS Reynolds operator");
    check(elrod_dp > 1.0e-3 * p_ref,
          "the polynomial mu field must reach the Elrod Reynolds operator");

    // Same reference viscosity, same reference temperature, different law: the
    // Andrade exponential and the Ferron quadratic must not land on each other.
    ThermalCase poly(ferron_polynomial_config(n_theta, CavitationModel::FBNS));
    ThermalCase andrade(thermal_config(n_theta, CavitationModel::FBNS));
    solve_coupled(poly);
    solve_coupled(andrade);
    double dmu = 0.0, dp = 0.0;
    for (int i = 0; i < poly.mesh.n_theta_local; ++i) {
        for (int j = 0; j < poly.mesh.n_z_local; ++j) {
            dmu = std::max(dmu, std::abs(poly.fields["mu"](i, j) - andrade.fields["mu"](i, j)));
            dp = std::max(dp,
                std::abs(poly.fields["pressure"](i, j) - andrade.fields["pressure"](i, j)));
        }
    }
    double gap[2] = {dmu, dp};
    MPI_Allreduce(MPI_IN_PLACE, gap, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    log_rank0("  POLYNOMIAL vs EMPIRICAL: max |dmu|/mu_ref = " + fmt(gap[0] / mu_ref) +
              ", max |dp|/p_ref = " + fmt(gap[1] / p_ref));
    check(gap[0] > 0.01 * mu_ref && gap[1] > 1.0e-4 * p_ref,
          "POLYNOMIAL and EMPIRICAL are different laws and must give different fields");
}

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
    test_inlet_overwrites_axial_boundary();
    test_scope_is_enforced();
    test_property_coupling_is_inert_when_off();
    test_viscosity_reaches_the_operator();
    test_energy_coupling_matches_elrod();
    test_cross_film_generalized_coefficients();
    test_polynomial_viscosity_couples();
    test_gaseous_cavitation();
    test_zero_gas_reproduces_phase_two();
    test_ceiling_is_fbns_only();

    if (rank == 0) std::cout << "PASS: test_fbns\n";
    PetscFinalize();
    return 0;
}
