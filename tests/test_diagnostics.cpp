// WP-4 validation: global mass balances, clamp accounting, and regime guards.
//
// Test 1 — Closed-domain conservation: periodic theta, zero-gradient z, no
//   inlets, fixed-bearing transient. Per-step liquid residual at the linear
//   solver tolerance.
// Test 2 — Open-domain accounting: Dirichlet z ends; balance closes after
//   boundary-flux accounting.
// Test 3 — Inlet accounting: penalty-pinned supply cells appear as an explicit
//   inlet mass source and the balance still closes.
// Test 4 — Clamp visibility: a strongly cavitating case with a high theta_min
//   reports nonzero clamp mass and the balance closes including it.
// Test 5 — Gas exchange closure: a no-flow supersaturated cell releases gas;
//   total propane is conserved and r_gas vanishes.
// Test 6 — Regime guards: Taylor-number warning and saturation-pressure
//   context warnings from validate().
//
// Run as: mpirun -n <N> ./test_diagnostics

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "diagnostics.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "fluid_properties.hpp"
#include "gas_transport.hpp"
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

bool contains_text(const std::vector<std::string>& values, const std::string& needle) {
    return std::any_of(values.begin(), values.end(), [&](const std::string& value) {
        return value.find(needle) != std::string::npos;
    });
}

void log_rank0(const std::string& message) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << message << "\n";
}

Fields make_elrod_fields(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields fields;
    const double theta0 = std::max(
        std::exp((cfg.initial_pressure() - cfg.p_cav) / cfg.bulk_modulus), cfg.theta_min);
    fields.add("pressure", mesh).fill(cfg.initial_pressure());
    fields.add("theta", mesh).fill(theta0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho * theta0);
    fields.add("cavitation_clamp_mass", mesh).fill(0.0);
    FilmThickness::compute_static(fields["h"], mesh, cfg);
    return fields;
}

// Run n_steps of the Elrod solve and return the worst per-step |r_l| plus the
// last balance for inspection.
double run_and_track_liquid_residual(SimulationConfig& cfg, int n_steps,
                                     Diagnostics::MassBalance& last_balance) {
    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);
    Fields fields = make_elrod_fields(mesh, cfg);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();

    double worst = 0.0;
    for (int step = 0; step < n_steps; ++step) {
        comm.update_ghosts(fields);
        Reynolds::solve_elrod(fields, sys, mesh, cfg);
        last_balance = Diagnostics::mass_balance(fields, mesh, cfg, cfg.dt);
        worst = std::max(worst, std::abs(last_balance.liquid_residual));
        fields["theta"].store_old_time();
    }
    return worst;
}

SimulationConfig base_transient_config() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.n_theta_global = 48;
    cfg.n_z_global = 12;
    cfg.p_cav = 1.0e5;
    cfg.bulk_modulus = 1.0e5;
    cfg.dt = 1.0e-3;
    cfg.omega = 100.0;
    cfg.outer_tol = 1.0e-12;
    cfg.max_outer_iters = 60;
    cfg.linear_rtol = 1.0e-13;
    cfg.log_outer_iters = false;
    cfg.inlets.clear();
    return cfg;
}

void test_closed_domain_conservation() {
    SimulationConfig cfg = base_transient_config();
    cfg.e = 0.6 * cfg.c;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;

    Diagnostics::MassBalance balance;
    const double worst = run_and_track_liquid_residual(cfg, 5, balance);
    log_rank0("  closed-domain worst |r_l| = " + std::to_string(worst));
    check(worst < 1.0e-10, "closed-domain liquid residual must sit at solver tolerance");
    check(std::abs(balance.liquid_boundary_flux) < 1.0e-30,
          "Neumann ends must report zero boundary flux");
}

void test_open_domain_accounting() {
    SimulationConfig cfg = base_transient_config();
    cfg.e = 0.3 * cfg.c;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 2.0e5;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 2.0e5;

    Diagnostics::MassBalance balance;
    const double worst = run_and_track_liquid_residual(cfg, 5, balance);
    log_rank0("  open-domain worst |r_l| = " + std::to_string(worst));
    check(worst < 1.0e-10, "open-domain balance must close after boundary-flux accounting");
}

void test_inlet_accounting() {
    SimulationConfig cfg = base_transient_config();
    cfg.e = 0.2 * cfg.c;
    cfg.bc_z_south_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = 2.0e5;
    cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_north_val = 2.0e5;
    InletConfig groove;
    groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90.0;
    groove.size = 20.0;
    groove.p_supply = 3.0e5;
    cfg.inlets.push_back(groove);

    Diagnostics::MassBalance balance;
    const double worst = run_and_track_liquid_residual(cfg, 5, balance);
    log_rank0("  inlet-case worst |r_l| = " + std::to_string(worst) +
              ", inlet source = " + std::to_string(balance.inlet_mass_source));
    check(balance.inlet_mass_source != 0.0,
          "penalty-pinned supply cells must report an inlet mass source");
    check(worst < 1.0e-8, "balance must close with inlet sources accounted");
}

void test_clamp_visibility() {
    SimulationConfig cfg = base_transient_config();
    cfg.e = 0.9 * cfg.c;
    cfg.theta_min = 0.5;
    cfg.dt = 5.0e-3;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;

    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);
    Fields fields = make_elrod_fields(mesh, cfg);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();

    bool clamp_seen = false;
    double worst = 0.0;
    for (int step = 0; step < 40 && !(clamp_seen && step > 10); ++step) {
        comm.update_ghosts(fields);
        Reynolds::solve_elrod(fields, sys, mesh, cfg);
        const Diagnostics::MassBalance balance =
            Diagnostics::mass_balance(fields, mesh, cfg, cfg.dt);
        if (balance.clamp_mass_liquid != 0.0) clamp_seen = true;
        worst = std::max(worst, std::abs(balance.liquid_residual));
        fields["theta"].store_old_time();
    }
    log_rank0("  clamp case worst |r_l| = " + std::to_string(worst) +
              ", clamp seen = " + std::to_string(clamp_seen));
    check(clamp_seen, "theta_min clamp mass must be reported, not silent");
    check(worst < 1.0e-10, "balance must close including clamp mass");
}

void test_gas_exchange_closure() {
    SimulationConfig cfg;
    cfg.cavitation_model = CavitationModel::JFO;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    cfg.solubility_model = SolubilityModel::HENRY;
    cfg.n_theta_global = 8;
    cfg.n_z_global = 2;
    cfg.p_cav = 1.0e5;
    cfg.omega = 0.0;
    cfg.dt = 1.0e-3;
    cfg.rho = 960.0;
    cfg.c = 1.0e-4;
    cfg.dissolved_gas_initial = 0.3;
    cfg.dissolved_gas_max = 0.5;
    cfg.dissolved_gas_henry_coeff = 2.0e-7;
    cfg.gas_mass_transfer_rate = 50.0;
    cfg.gas_alpha_max = 1.0;
    cfg.gas_pressure_floor = 1.0e5;
    cfg.temperature_initial = 313.15;
    cfg.bc_z_south_type = BCType::NEUMANN;
    cfg.bc_z_north_type = BCType::NEUMANN;
    cfg.linear_rtol = 1.0e-13;
    cfg.inlets.clear();

    Mesh mesh(cfg);
    Communicator comm(mesh);
    LinearSystem sys(mesh);

    Fields fields;
    fields.add("pressure", mesh).fill(1.0e6);
    fields.add("theta", mesh).fill(1.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);
    fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
    fields.add("mu", mesh).fill(cfg.mu);
    fields.add("temperature", mesh).fill(cfg.temperature_initial);
    fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
    fields.add("free_gas_mass", mesh).fill(0.0);
    fields.add("alpha_gas", mesh).fill(0.0);
    fields.add("gas_mass_transfer", mesh).fill(0.0);
    fields.add("gas_clamp_mass", mesh).fill(0.0);
    fields.add("cavitation_clamp_mass", mesh).fill(0.0);
    comm.update_ghosts(fields);
    fields["theta"].store_old_time();
    fields["dissolved_gas"].store_old_time();
    fields["free_gas_mass"].store_old_time();

    Diagnostics::MassBalance initial =
        Diagnostics::mass_balance(fields, mesh, cfg, 0.0);
    const double gas_mass_start = initial.gas_mass;

    double worst_gas_residual = 0.0;
    for (int step = 0; step < 5; ++step) {
        comm.update_ghosts(fields["dissolved_gas"]);
        comm.update_ghosts(fields["free_gas_mass"]);
        GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);
        FluidProperties::update_gas_state(fields, mesh, cfg, cfg.dt);
        const Diagnostics::MassBalance balance =
            Diagnostics::mass_balance(fields, mesh, cfg, cfg.dt);
        worst_gas_residual = std::max(worst_gas_residual, std::abs(balance.gas_residual));
        fields["dissolved_gas"].store_old_time();
        fields["free_gas_mass"].store_old_time();
    }

    const Diagnostics::MassBalance final_balance =
        Diagnostics::mass_balance(fields, mesh, cfg, 0.0);
    const double drift =
        std::abs(final_balance.gas_mass - gas_mass_start) / std::max(gas_mass_start, 1e-30);

    double released = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j)
            released = std::max(released, fields["free_gas_mass"](i, j));
    MPI_Allreduce(MPI_IN_PLACE, &released, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    log_rank0("  gas closure: worst |r_gas| = " + std::to_string(worst_gas_residual) +
              ", total drift = " + std::to_string(drift) +
              ", released m_g = " + std::to_string(released));
    check(released > 0.0, "supersaturated cell must release free gas");
    check(worst_gas_residual < 1.0e-10,
          "release/resorption exchange must cancel in the total gas balance");
    check(drift < 1.0e-12, "total propane mass must be conserved in a closed cell");
}

// WP-11 follow-up: under consistent_boundary_flux a cavitated axial-boundary cell
// re-floods at the reformation rate. The reflood must convect the reservoir
// dissolved-gas composition into the cell (it carries nothing without the flag),
// and the diagnostics gas-boundary flux must mirror the same reflood (nonzero with
// the flag, zero without -- a cavitated full-film face carries none). omega = 0,
// uniform cavitated film, no kinetics/diffusion, so the axial boundary inflow is
// the only thing that can change the composition. The film content theta is held
// fixed to isolate the boundary convection, so the *global* gas residual is not a
// closure check here (frozen theta against a net inflow breaks liquid continuity);
// the gas balance closure under solved continuity is covered by the fixed-theta
// closed cell in test_gas_exchange_closure.
void test_reformation_boundary_reservoir_state() {
    auto run = [](bool reflood, double& boundary_c, double& interior_c,
                  double& gas_boundary_flux) {
        SimulationConfig cfg;
        cfg.cavitation_model = CavitationModel::JFO;
        cfg.solution_mode = SolutionMode::TRANSIENT;
        cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
        cfg.solubility_model = SolubilityModel::HENRY;
        cfg.density_model = DensityModel::CONSTANT;
        cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
        cfg.n_theta_global = 8;
        cfg.n_z_global = 4;
        cfg.p_cav = 1.0e5;
        cfg.bulk_modulus = 1.0e9;
        cfg.omega = 0.0;                       // no Couette: isolate the axial inflow
        cfg.dt = 1.0e-3;
        cfg.rho = 960.0;
        cfg.c = 1.0e-4;
        cfg.dissolved_gas_initial = 0.2175;    // reservoir composition c_in
        cfg.dissolved_gas_max = 0.5;
        cfg.dissolved_gas_henry_coeff = 2.175e-7;
        cfg.gas_mass_transfer_rate = 0.0;      // disable 0-D kinetics
        cfg.dissolved_gas_diffusivity = 0.0;   // disable diffusion
        cfg.gas_alpha_max = 0.6;
        cfg.gas_pressure_floor = 1.0e5;
        cfg.temperature_initial = 313.15;
        cfg.bc_z_south_type = BCType::DIRICHLET;
        cfg.bc_z_south_val = 1.0e6;            // supply > p_cav: reflood drives inflow
        cfg.bc_z_north_type = BCType::DIRICHLET;
        cfg.bc_z_north_val = 1.0e6;
        cfg.consistent_boundary_flux = reflood;
        cfg.linear_rtol = 1.0e-13;
        cfg.inlets.clear();

        Mesh mesh(cfg);
        Communicator comm(mesh);
        LinearSystem sys(mesh);

        const double theta_cav = 0.5;          // cavitated film at the boundary
        const double c_start = 0.05;           // depleted, below c_in
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);            // cavitated plateau
        fields.add("theta", mesh).fill(theta_cav);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho * theta_cav);       // Elrod density rho_l*theta
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);   // full liquid density
        fields.add("mu", mesh).fill(cfg.mu);
        fields.add("temperature", mesh).fill(cfg.temperature_initial);
        fields.add("dissolved_gas", mesh).fill(c_start);
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("alpha_gas", mesh).fill(0.0);
        fields.add("gas_mass_transfer", mesh).fill(0.0);
        fields.add("gas_clamp_mass", mesh).fill(0.0);
        comm.update_ghosts(fields);
        fields["dissolved_gas"].store_old_time();
        fields["free_gas_mass"].store_old_time();
        fields["theta"].store_old_time();

        Diagnostics::MassBalance balance;
        for (int step = 0; step < 20; ++step) {
            comm.update_ghosts(fields);
            GasTransport::solve(fields, sys, mesh, cfg, cfg.dt);
            balance = Diagnostics::mass_balance(fields, mesh, cfg, cfg.dt);
            fields["dissolved_gas"].store_old_time();
            fields["free_gas_mass"].store_old_time();
        }
        gas_boundary_flux = balance.gas_boundary_flux;  // already MPI-reduced

        double bc = 0.0, ic = 0.0;             // south boundary (j=0) vs interior (j=1)
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            bc = std::max(bc, fields["dissolved_gas"](i, 0));
            ic = std::max(ic, fields["dissolved_gas"](i, 1));
        }
        MPI_Allreduce(&bc, &boundary_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&ic, &interior_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    };

    double c_bound_off = 0, c_int_off = 0, flux_off = 0;
    double c_bound_on = 0, c_int_on = 0, flux_on = 0;
    run(false, c_bound_off, c_int_off, flux_off);
    run(true, c_bound_on, c_int_on, flux_on);

    log_rank0("  reformation BC: c_boundary OFF=" + std::to_string(c_bound_off) +
              " ON=" + std::to_string(c_bound_on) + " (c_in=0.2175); gas_boundary_flux OFF=" +
              std::to_string(flux_off) + " ON=" + std::to_string(flux_on));

    const double c_in = 0.2175, c_start = 0.05;
    check(std::abs(c_bound_off - c_start) < 1e-6,
          "reformation BC OFF: cavitated boundary cell receives no reflood (stays depleted)");
    check(c_bound_on > c_start + 1e-4,
          "reformation BC ON: reflood re-saturates the boundary cell toward c_in");
    check(c_bound_on <= c_in + 1e-9, "reformation BC ON: boundary cell stays bounded by c_in");
    check(c_bound_on > c_int_on + 1e-6,
          "reformation BC ON: boundary cell gains more gas than the interior");
    // Diagnostics mirror: a cavitated boundary carries no full-film flux (OFF), and
    // the diagnostic picks up the reflood inflow (ON) -- same formula as the solver.
    check(std::abs(flux_off) < 1e-30,
          "reformation BC OFF: diagnostics report zero gas boundary flux at a cavitated end");
    check(flux_on > 0.0,
          "reformation BC ON: diagnostics mirror the reflood as a positive gas inflow");
}

void test_regime_and_saturation_guards() {
    std::vector<std::string> errors;
    std::vector<std::string> warnings;

    // Taylor-number guard: low viscosity, fast shaft.
    SimulationConfig turbulent;
    turbulent.p_cav = 1.0e5;
    turbulent.bc_z_south_val = 2.0e5;
    turbulent.bc_z_north_val = 2.0e5;
    turbulent.mu = 1.0e-6;
    turbulent.omega = 1000.0;
    turbulent.validate(errors, warnings);
    check(errors.empty(), "low-viscosity case should warn, not error");
    check(contains_text(warnings, "Taylor"),
          "validate() must warn when the Taylor-vortex limit is exceeded");

    // Saturation context for the shipped R290 constants: p_sat = 1 MPa and the
    // 1 MPa boundary feed enters exactly saturated.
    SimulationConfig r290;
    r290.p_cav = 1.0e5;
    r290.bc_z_south_val = 1.0e6;
    r290.bc_z_north_val = 1.0e6;
    r290.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    r290.solubility_model = SolubilityModel::HENRY;
    r290.dissolved_gas_initial = 0.2175;
    r290.dissolved_gas_max = 0.5;
    r290.dissolved_gas_henry_coeff = 2.175e-7;
    r290.temperature_initial = 313.15;
    r290.property_reference_temperature = 313.15;
    r290.gas_mass_transfer_rate = 500.0;
    r290.gas_alpha_max = 0.6;
    r290.mixture_viscosity_model = MixtureViscosityModel::MCADAMS;
    r290.mu_gas = 8.7e-6;
    r290.validate(errors, warnings);
    check(errors.empty(), "R290 saturation-context case should validate");
    check(std::abs(r290.saturation_pressure(r290.dissolved_gas_initial,
                                            r290.temperature_initial) - 1.0e6) < 1.0,
          "p_sat for the shipped R290 constants must be 1.0 MPa");
    check(contains_text(warnings, "p_sat"),
          "validate() must log the saturation-pressure context");
    check(contains_text(warnings, "saturated"),
          "validate() must warn when boundary oil enters saturated");
}

}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Running test_diagnostics...\n";

    test_regime_and_saturation_guards();
    test_closed_domain_conservation();
    test_open_domain_accounting();
    test_inlet_accounting();
    test_clamp_visibility();
    test_gas_exchange_closure();
    test_reformation_boundary_reservoir_state();

    if (rank == 0) std::cout << "PASS: test_diagnostics\n";
    PetscFinalize();
    return 0;
}
