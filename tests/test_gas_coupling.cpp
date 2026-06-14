// WP-1 acceptance gates for the two-phase (VOID_COUPLED) gas-pressure coupling.
//   - 0-D closed cell: dissolved->free gas release conserves total propane and
//     reaches a steady state (the finite-rate kinetics that drive the coupling).
//   - 0-D resorption: under-saturated liquid reabsorbs free gas.
//   - Vaporous limit: with no gas, VOID_COUPLED reproduces the single-phase
//     Elrod solve exactly (the decoupled-limit gate the AUDIT_PLAN requires).
//
// Run as: mpiexec -n <N> ./test_gas_coupling

#include <algorithm>
#include <cmath>
#include <iostream>

#include <petsc.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "fluid_properties.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
static bool close_rel(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (std::fabs(b) + 1e-30);
}

static SimulationConfig gas_config() {
    SimulationConfig cfg;
    cfg.p_cav = 1.0e5;
    cfg.fluid_property_model = FluidPropertyModel::GAS_CAVITATION_MIXTURE;
    cfg.dissolved_gas_species = DissolvedGasSpecies::PROPANE;
    cfg.oil_gas_solution_model = OilGasSolutionModel::HENRY;
    cfg.density_model = DensityModel::PURE_OIL;   // rho_l = cfg.rho deterministically
    cfg.viscosity_model = ViscosityModel::PURE_OIL;
    cfg.dissolved_gas_henry_coeff = 2.175e-7;     // p_sat(0.2175) = 1.0 MPa
    cfg.dissolved_gas_initial = 0.2175;
    cfg.dissolved_gas_max = 0.5;
    cfg.gas_mass_transfer_rate = 500.0;
    cfg.gas_alpha_max = 0.6;
    cfg.mu_gas = 8.7e-6;
    cfg.gas_pressure_floor = 1.0e5;
    cfg.property_reference_pressure = 1.0e6;
    cfg.property_reference_temperature = 313.15;
    cfg.temperature_initial = 313.15;
    return cfg;
}

static Fields make_cell_fields(const Mesh& mesh, const SimulationConfig& cfg,
                               double pressure, double c_d, double m_g) {
    Fields fields;
    fields.add("pressure", mesh).fill(pressure);
    fields.add("temperature", mesh).fill(cfg.temperature_initial);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("theta", mesh).fill(1.0);
    fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
    fields.add("dissolved_gas", mesh).fill(c_d);
    fields.add("free_gas_mass", mesh).fill(m_g);
    fields.add("alpha_gas", mesh).fill(0.0);
    fields.add("gas_mass_transfer", mesh).fill(0.0);
    fields.add("rho_gas", mesh).fill(0.0);
    fields.add("inlet_indicator", mesh).fill(0.0);
    return fields;
}

// Total dissolved + free propane per unit area in cell (0,0) [kg/m^2].
static double propane_per_area(const Fields& fields, const SimulationConfig& cfg) {
    const double c_d = fields["dissolved_gas"](0, 0);
    const double m_g = fields["free_gas_mass"](0, 0);
    return c_d * cfg.rho * 1.0 * cfg.c + m_g;  // theta = 1, h = cfg.c
}

// Gate 1: sub-saturated cell releases gas; total propane is conserved and steady.
static void test_release_conservation(const Mesh& mesh) {
    SimulationConfig cfg = gas_config();
    // p = 0.5 MPa < p_sat (1.0 MPa): the 21.75% feed is supersaturated and outgasses.
    Fields fields = make_cell_fields(mesh, cfg, 5.0e5, cfg.dissolved_gas_initial, 0.0);

    const double total0 = propane_per_area(fields, cfg);
    const double c_d0 = fields["dissolved_gas"](0, 0);
    const double dt = 1.0e-5;
    for (int step = 0; step < 4000; ++step)
        FluidProperties::update_gas_state(fields, mesh, cfg, dt);

    const double c_d1 = fields["dissolved_gas"](0, 0);
    const double m_g1 = fields["free_gas_mass"](0, 0);
    const double total1 = propane_per_area(fields, cfg);

    check(close_rel(total1, total0, 1e-9), "release: total propane conserved");
    check(c_d1 < c_d0 - 1e-6, "release: dissolved gas drops below the feed");
    check(m_g1 > 0.0, "release: free gas appears");
    check(fields["alpha_gas"](0, 0) <= cfg.gas_alpha_max + 1e-12, "release: void fraction <= alpha_max");

    // Steady: another block of steps does not move the state.
    for (int step = 0; step < 500; ++step)
        FluidProperties::update_gas_state(fields, mesh, cfg, dt);
    check(close_rel(fields["dissolved_gas"](0, 0), c_d1, 1e-6), "release: reaches a steady dissolved fraction");
    check(close_rel(propane_per_area(fields, cfg), total0, 1e-9), "release: mass still conserved at steady");
}

// Gate 2: under-saturated liquid with free gas reabsorbs it (mass conserved).
static void test_resorption(const Mesh& mesh) {
    SimulationConfig cfg = gas_config();
    // p = 1.0 MPa -> c_eq = 0.2175; start under-saturated with free gas present.
    Fields fields = make_cell_fields(mesh, cfg, 1.0e6, 0.05, 1.0e-5);
    const double total0 = propane_per_area(fields, cfg);
    const double c_d0 = fields["dissolved_gas"](0, 0);
    const double m_g0 = fields["free_gas_mass"](0, 0);

    const double dt = 1.0e-5;
    for (int step = 0; step < 4000; ++step)
        FluidProperties::update_gas_state(fields, mesh, cfg, dt);

    check(close_rel(propane_per_area(fields, cfg), total0, 1e-9), "resorption: total propane conserved");
    check(fields["dissolved_gas"](0, 0) > c_d0 + 1e-9, "resorption: dissolved gas rises");
    check(fields["free_gas_mass"](0, 0) < m_g0 - 1e-12, "resorption: free gas reabsorbs");
}

// Gate 3 (vaporous limit): with no gas, VOID_COUPLED == single-phase Elrod (NONE).
static void test_vaporous_limit(const Mesh& mesh, Communicator& comm) {
    auto solve_theta = [&](GasPressureCoupling coupling) {
        SimulationConfig cfg;                       // default journal bearing (cavitates)
        cfg.solution_mode = SolutionMode::STEADY_STATE;
        cfg.cavitation_model = CavitationModel::ELROD_ADAMS;
        cfg.gas_pressure_coupling = coupling;
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("theta", mesh).fill(1.0);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("rho", mesh).fill(cfg.rho);
        fields.add("alpha_gas", mesh).fill(0.0);    // no free gas
        fields.add("free_gas_mass", mesh).fill(0.0);
        fields.add("dissolved_gas", mesh).fill(0.0);
        FilmThickness::compute_static(fields["h"], mesh, cfg);
        comm.update_ghosts(fields);
        LinearSystem sys(mesh);
        Reynolds::solve_elrod(fields, sys, mesh, cfg);
        return fields;  // copy out theta
    };

    Fields none = solve_theta(GasPressureCoupling::NONE);
    Fields coupled = solve_theta(GasPressureCoupling::VOID_COUPLED);

    double max_dtheta = 0.0, max_dp = 0.0;
    for (int i = 0; i < mesh.n_theta_local; ++i)
        for (int j = 0; j < mesh.n_z_local; ++j) {
            max_dtheta = std::max(max_dtheta, std::fabs(none["theta"](i, j) - coupled["theta"](i, j)));
            max_dp = std::max(max_dp, std::fabs(none["pressure"](i, j) - coupled["pressure"](i, j)));
        }
    MPI_Allreduce(MPI_IN_PLACE, &max_dtheta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &max_dp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    check(max_dtheta < 1e-12, "vaporous limit: VOID_COUPLED theta == NONE with no gas");
    check(max_dp < 1e-6, "vaporous limit: VOID_COUPLED pressure == NONE with no gas");
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    {
        SimulationConfig cfg;
        Mesh mesh(cfg);
        Communicator comm(mesh);

        test_release_conservation(mesh);
        test_resorption(mesh);
        test_vaporous_limit(mesh, comm);

        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << "PASS: test_gas_coupling\n";
    }
    PetscFinalize();
    return 0;
}
