// Test #3 of cases/00c_jfo_verification: the per-cell continuity residual field.
//
// The global mass balance closes to round-off, which by itself cannot tell a
// cell-by-cell conservative scheme apart from large local errors that cancel.
// These tests pin down the per-cell field instead.
//
// Test 1 - Steady grooved bearing (the 00c geometry: JFO, Neumann axial ends,
//   one pressurized groove). Every cell that is not an inlet, not clamped and
//   not on the rupture front sits at the linear solver tolerance, and the field
//   sums to the global closure mass_balance() reports.
// Test 2 - Open axial ends. The Dirichlet boundary rows carry their boundary
//   link, so they are at tolerance too; the sum identity holds with a nonzero
//   boundary flux and transient storage.
// Test 3 - Localization. Perturbing one converged cell lights up exactly that
//   cell and its four neighbours; every other cell stays at tolerance and the
//   sum identity still holds. Guards against a field that is merely zeroed.
// Test 4 - Absent field. mass_residual() is a no-op and mass_balance() is
//   unchanged, so a run that does not select the field pays nothing.
//
// Run as: mpirun -n <N> ./test_mass_residual

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "diagnostics.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

namespace {

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

// std::to_string is fixed-point and prints every residual here as 0.000000.
std::string fmt(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(3) << value;
    return out.str();
}

// Upper bound on a single Couette face flux [kg/s]: the mass throughput the
// residual has to be small against to be useful as a diagnostic.
double flux_scale(const Mesh& mesh, const SimulationConfig& cfg) {
    return cfg.rho * std::abs(cfg.omega) * mesh.R * (2.0 * cfg.c) * mesh.get_d_z();
}

// Size of the individual terms whose near-cancellation forms the residual. The
// barotropic Poiseuille conductance rho_0 beta h^3/(12 mu) * face/distance beats
// the Couette flux by ~beta/p, so it, not the throughput, sets what "round-off"
// means for one row.
double row_scale(const Mesh& mesh, const SimulationConfig& cfg) {
    const double h = cfg.c + cfg.e;  // thickest film in the gap
    const double gamma = cfg.rho * cfg.bulk_modulus * h * h * h / (12.0 * cfg.mu);
    const double d_theta = mesh.get_d_theta();
    const double d_z = mesh.get_d_z();
    const double conductance =
        gamma * std::max(d_z / (mesh.R * d_theta), mesh.R * d_theta / d_z);
    return std::max(conductance, flux_scale(mesh, cfg));
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
        const double theta0 = std::max(
            std::exp((cfg.initial_pressure() - cfg.p_cav) / cfg.bulk_modulus), cfg.theta_min);
        fields.add("pressure", mesh).fill(cfg.initial_pressure());
        fields.add("theta", mesh).fill(theta0);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho * theta0);
        fields.add("cavitation_clamp_mass", mesh).fill(0.0);
        fields.add("mass_residual", mesh).fill(0.0);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        FilmThickness::compute_inlet_indicator(inlet_mask, mesh, cfg);
        comm.update_ghosts(inlet_mask);
        comm.update_ghosts(fields);
        fields["theta"].store_old_time();
    }

    void solve() {
        comm.update_ghosts(fields);
        Reynolds::solve_elrod(fields, sys, mesh, cfg);
        comm.update_ghosts(fields);
    }
};

// Signed sum, magnitude sum, and the worst value over "clean" cells: cells that
// are not an inlet or next to one, carry no clamp mass, and are full film with
// full-film neighbours (so they are away from the rupture/reformation fronts).
// Those are exactly the cells whose residual must be at the solver tolerance.
struct ResidualStats {
    double sum = 0.0;
    double sum_abs = 0.0;
    double worst_clean = 0.0;
    double clean_cells = 0.0;
};

ResidualStats collect(const Case& c) {
    const Field& r = c.fields["mass_residual"];
    const Field& theta = c.fields["theta"];
    const Field& clamp = c.fields["cavitation_clamp_mass"];
    const Mesh& mesh = c.mesh;

    ResidualStats s;
    for (int i = 0; i < mesh.n_theta_local; ++i) {
        for (int j = 0; j < mesh.n_z_local; ++j) {
            s.sum += r(i, j);
            s.sum_abs += std::abs(r(i, j));

            const bool no_inlet = c.inlet_mask(i, j) == 0.0 &&
                c.inlet_mask(i - 1, j) == 0.0 && c.inlet_mask(i + 1, j) == 0.0;
            const bool full_film = theta(i, j) >= 1.0 && theta(i - 1, j) >= 1.0 &&
                theta(i + 1, j) >= 1.0 &&
                (j == 0 || theta(i, j - 1) >= 1.0) &&
                (j == mesh.n_z_local - 1 || theta(i, j + 1) >= 1.0);
            if (no_inlet && full_film && clamp(i, j) == 0.0) {
                s.worst_clean = std::max(s.worst_clean, std::abs(r(i, j)));
                s.clean_cells += 1.0;
            }
        }
    }

    double vals[4] = {s.sum, s.sum_abs, s.worst_clean, s.clean_cells};
    MPI_Allreduce(MPI_IN_PLACE, vals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, vals + 2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, vals + 3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    s.sum = vals[0];
    s.sum_abs = vals[1];
    s.worst_clean = vals[2];
    s.clean_cells = vals[3];
    return s;
}

// Interior faces telescope, so summing the field must reproduce the global
// closure. Compared in the normalized units mass_balance() reports (r_l), which
// is the number the diagnostics CSV and the per-step log actually carry.
void check_sum_identity(const Case& c, const ResidualStats& s,
                        const Diagnostics::MassBalance& balance, double dt,
                        const char* label) {
    const bool transient = c.cfg.solution_mode != SolutionMode::STEADY_STATE && dt > 0.0;
    const double step = transient ? dt : std::max(c.cfg.dt, 1.0e-30);
    const double from_field = s.sum * step / std::max(balance.liquid_mass, 1.0e-30);
    const double mismatch = std::abs(from_field - balance.liquid_residual);
    log_rank0(std::string("  ") + label + " r_l from field = " + fmt(from_field) +
              ", from diagnostics = " + fmt(balance.liquid_residual) +
              ", |difference| = " + fmt(mismatch));
    check(mismatch < 1.0e-13,
          "mass_residual must sum to the global closure the diagnostics report");
}

SimulationConfig base_config() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.n_theta_global = 64;
    cfg.n_z_global = 8;
    cfg.p_cav = 1.0e5;
    cfg.bulk_modulus = 1.0e9;
    cfg.omega = 100.0;
    cfg.e = 0.6 * cfg.c;
    cfg.outer_tol = 1.0e-12;
    cfg.max_outer_iters = 80;
    cfg.linear_rtol = 1.0e-13;
    cfg.log_outer_iters = false;
    cfg.inlets.clear();
    return cfg;
}

// 00c geometry: steady JFO, closed axial ends, one pressurized groove. Without a
// supply the closed domain runs cavitated almost everywhere and leaves no
// full-film interior to check.
SimulationConfig groove_config() {
    SimulationConfig cfg = base_config();
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0;
    groove.size = 18.0;
    groove.p_supply = 1.713e5;
    cfg.inlets.push_back(groove);
    return cfg;
}

void test_steady_groove_residual() {
    const SimulationConfig cfg = groove_config();

    Case c(cfg);
    c.solve();
    Diagnostics::mass_residual(c.fields, c.mesh, cfg, 0.0);
    const Diagnostics::MassBalance balance =
        Diagnostics::mass_balance(c.fields, c.mesh, cfg, 0.0);

    const ResidualStats s = collect(c);
    const double row = row_scale(c.mesh, cfg);
    log_rank0("  steady groove: clean cells = " + std::to_string((int)s.clean_cells) +
              ", worst clean |r|/row = " + fmt(s.worst_clean / row) +
              ", worst clean |r|/flux = " + fmt(s.worst_clean / flux_scale(c.mesh, cfg)));
    check(s.clean_cells > 0.5 * cfg.n_theta_global * cfg.n_z_global,
          "most of a moderately loaded bearing must be clean full-film cells");
    check(s.worst_clean / row < 1.0e-12,
          "clean full-film cells must be at round-off of the assembled row");
    check_sum_identity(c, s, balance, 0.0, "steady groove");
}

// Open axial ends: the Dirichlet link belongs to the boundary row's balance, so
// those rows must close as well as the interior.
void test_open_boundary_residual() {
    SimulationConfig cfg = base_config();
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.dt = 1.0e-4;
    cfg.e = 0.4 * cfg.c;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 2.0e5;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 2.0e5;

    Case c(cfg);
    for (int step = 0; step < 3; ++step) {
        c.solve();
        Diagnostics::mass_residual(c.fields, c.mesh, cfg, cfg.dt);
        const Diagnostics::MassBalance balance =
            Diagnostics::mass_balance(c.fields, c.mesh, cfg, cfg.dt);
        const ResidualStats s = collect(c);
        const double row = row_scale(c.mesh, cfg);
        if (step == 2) {
            log_rank0("  open ends: worst clean |r|/row = " + fmt(s.worst_clean / row) +
                      ", boundary flux = " + fmt(balance.liquid_boundary_flux));
            check(std::abs(balance.liquid_boundary_flux) > 0.0,
                  "Dirichlet ends must carry a nonzero boundary flux for this to test anything");
            check(s.worst_clean / row < 1.0e-12,
                  "Dirichlet boundary rows must close as tightly as the interior");
            check_sum_identity(c, s, balance, cfg.dt, "open ends");
        }
        c.fields["theta"].store_old_time();
    }
}

// A converged state perturbed in one cell: the residual must localize on that
// cell and its four neighbours and stay at tolerance everywhere else.
void test_perturbation_is_localized() {
    const SimulationConfig cfg = groove_config();

    Case c(cfg);
    c.solve();
    const Diagnostics::MassBalance before =
        Diagnostics::mass_balance(c.fields, c.mesh, cfg, 0.0);

    // Pick a full-film cell with full-film neighbours and no groove cell in the
    // stencil (a groove cell's residual is zeroed by construction), on the lowest
    // rank that has one, so exactly one cell is perturbed across the communicator.
    const Field& theta_ro = c.fields["theta"];
    int rank = 0, n_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
    int local_i = -1, local_j = -1;
    for (int i = 0; i < c.mesh.n_theta_local && local_i < 0; ++i) {
        for (int j = 1; j < c.mesh.n_z_local - 1; ++j) {
            if (c.inlet_mask(i, j) != 0.0 || c.inlet_mask(i - 1, j) != 0.0 ||
                c.inlet_mask(i + 1, j) != 0.0) {
                continue;
            }
            if (theta_ro(i, j) >= 1.0 && theta_ro(i - 1, j) >= 1.0 &&
                theta_ro(i + 1, j) >= 1.0 && theta_ro(i, j - 1) >= 1.0 &&
                theta_ro(i, j + 1) >= 1.0) {
                local_i = i;
                local_j = j;
                break;
            }
        }
    }
    int candidate = (local_i >= 0) ? rank : n_ranks;
    int owner = candidate;
    MPI_Allreduce(&candidate, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    check(owner < n_ranks, "the perturbation test needs one full-film interior cell");

    if (rank == owner) c.fields["theta"](local_i, local_j) *= 1.02;
    c.comm.update_ghosts(c.fields["theta"]);

    Diagnostics::mass_residual(c.fields, c.mesh, cfg, 0.0);
    const Diagnostics::MassBalance balance =
        Diagnostics::mass_balance(c.fields, c.mesh, cfg, 0.0);

    const Field& r = c.fields["mass_residual"];
    const double scale = flux_scale(c.mesh, cfg);
    double hit = 0.0;      // largest |r| among the perturbed cell + neighbours
    double clean = 0.0;    // largest |r| elsewhere
    double signed_sum = 0.0;
    double abs_sum = 0.0;
    for (int i = 0; i < c.mesh.n_theta_local; ++i) {
        for (int j = 0; j < c.mesh.n_z_local; ++j) {
            signed_sum += r(i, j);
            abs_sum += std::abs(r(i, j));
            const bool touched = (rank == owner) &&
                ((i == local_i && std::abs(j - local_j) <= 1) ||
                 (j == local_j && std::abs(i - local_i) <= 1));
            // A perturbed cell owned by another rank shows up here as a theta
            // ghost, so also treat the two theta-edge columns as possibly touched.
            const bool edge_column = (i == 0 || i == c.mesh.n_theta_local - 1);
            if (touched) hit = std::max(hit, std::abs(r(i, j)));
            else if (!edge_column) clean = std::max(clean, std::abs(r(i, j)));
        }
    }
    double reduced[4] = {hit, clean, signed_sum, abs_sum};
    MPI_Allreduce(MPI_IN_PLACE, reduced, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, reduced + 2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    log_rank0("  perturbed cell: |r|_hit/flux = " + fmt(reduced[0] / scale) +
              ", |r|_elsewhere/flux = " + fmt(reduced[1] / scale));
    check(reduced[0] / scale > 1.0e-6,
          "the perturbed cell must show a large residual (field is not a stub)");
    check(reduced[1] < 1.0e-6 * reduced[0],
          "the perturbation must stay local: untouched cells keep their tolerance");

    // The reason the field exists: the five disturbed cells' errors telescope, so
    // the global balance is unchanged while the field localizes the defect.
    log_rank0("  global r_l before = " + fmt(before.liquid_residual) +
              ", after = " + fmt(balance.liquid_residual));
    check(std::abs(balance.liquid_residual - before.liquid_residual) < 1.0e-13,
          "the global balance must stay blind to the local defect the field reveals");

    ResidualStats s;
    s.sum = reduced[2];
    s.sum_abs = reduced[3];
    check_sum_identity(c, s, balance, 0.0, "perturbed");
}

// Not selecting the field must cost nothing and change nothing.
void test_absent_field_is_noop() {
    SimulationConfig cfg = base_config();
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);
    Fields fields;
    const double theta0 = std::max(
        std::exp((cfg.initial_pressure() - cfg.p_cav) / cfg.bulk_modulus), cfg.theta_min);
    fields.add("pressure", mesh).fill(cfg.initial_pressure());
    fields.add("theta", mesh).fill(theta0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho * theta0);
    fields.add("cavitation_clamp_mass", mesh).fill(0.0);
    FilmThickness::compute_static(fields["h"], mesh, cfg);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();

    comm.update_ghosts(fields);
    Reynolds::solve_elrod(fields, sys, mesh, cfg);
    const Diagnostics::MassBalance before =
        Diagnostics::mass_balance(fields, mesh, cfg, 0.0);
    Diagnostics::mass_residual(fields, mesh, cfg, 0.0);  // must not throw or add a field
    const Diagnostics::MassBalance after =
        Diagnostics::mass_balance(fields, mesh, cfg, 0.0);

    check(!fields.has("mass_residual"), "mass_residual() must not register the field itself");
    check(before.liquid_mass == after.liquid_mass &&
          before.liquid_residual == after.liquid_residual,
          "mass_residual() must not perturb the solution or the global balance");
    log_rank0("  absent field: no-op confirmed");
}

}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_mass_residual...\n";

    test_steady_groove_residual();
    test_open_boundary_residual();
    test_perturbation_is_localized();
    test_absent_field_is_noop();

    if (rank == 0) std::cout << "PASS: test_mass_residual\n";
    PetscFinalize();
    return 0;
}
