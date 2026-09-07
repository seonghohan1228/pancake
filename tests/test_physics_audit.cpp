#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <petsc.h>

#include "communicator.hpp"
#include "coupling.hpp"
#include "diagnostics.hpp"
#include "fbns.hpp"
#include "fluid_properties.hpp"
#include "film_thickness.hpp"
#include "gas_transport.hpp"
#include "reynolds.hpp"

namespace {
int failures = 0;

void check(bool condition, const char* message) {
    int bad = condition ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (bad && rank == 0) std::cerr << "FAIL: " << message << '\n';
    failures += bad;
}

bool close(double a, double b, double tolerance = 1e-10) {
    return std::abs(a - b) <= tolerance * std::max(std::abs(b), 1e-20);
}

SimulationConfig config() {
    SimulationConfig cfg;
    cfg.n_theta_global = 24;
    cfg.n_z_global = 4;
    cfg.cavitation_model = CavitationModel::FBNS;
    cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
    cfg.solution_mode = SolutionMode::TRANSIENT;
    cfg.density_model = DensityModel::CONSTANT;
    cfg.p_cav = 1e5;
    cfg.dissolved_gas_henry_coeff = 1e-10;
    cfg.dissolved_gas_initial = 0.002;
    cfg.dissolved_gas_max = 0.5;
    cfg.gas_mass_transfer_rate = 20;
    cfg.gas_alpha_max = 0.6;
    cfg.omega = 0;
    cfg.bc_z_south_type = cfg.bc_z_north_type = BCType::NEUMANN;
    cfg.linear_rtol = 1e-12;
    return cfg;
}

Fields fields_for(const Mesh& mesh, const SimulationConfig& cfg) {
    Fields f;
    f.add("pressure", mesh).fill(1e6);
    f.add("temperature", mesh).fill(300);
    f.add("theta", mesh).fill(0.5);
    f.add("h", mesh).fill(cfg.c);
    f.add("rho_liquid_solution", mesh).fill(cfg.rho);
    f.add("rho", mesh).fill(0.5 * cfg.rho);
    f.add("mu", mesh).fill(cfg.mu);
    f.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
    f.add("free_gas_mass", mesh).fill(0);
    f.add("alpha_gas", mesh).fill(0);
    f.add("gas_mass_transfer", mesh).fill(0);
    f.add("rho_gas", mesh).fill(0);
    f.add("inlet_indicator", mesh).fill(0);
    f.add("gas_clamp_mass", mesh).fill(0);
    for (auto& [name, field] : f) field->store_old_time();
    return f;
}

double liquid_mass(const Fields& f) {
    return f["rho_liquid_solution"](0, 0) * f["theta"](0, 0) * f["h"](0, 0);
}

void test_phase_exchange(const Mesh& mesh) {
    auto cfg = config();
    for (auto density : {DensityModel::CONSTANT, DensityModel::MIXTURE}) {
        cfg.density_model = density;
        cfg.dissolved_gas_liquid_density = 500;
        auto f = fields_for(mesh, cfg);
        FluidProperties::update_solution_fields(f, mesh, cfg);
        f.add("mass_residual", mesh).fill(0);
        for (auto& [name, field] : f) field->store_old_time();
        const double liquid0 = liquid_mass(f);
        const double gas0 = liquid0 * f["dissolved_gas"](0, 0);
        const double oil0 = liquid0 - gas0;
        FluidProperties::update_gas_state(f, mesh, cfg, 0.001);
        FluidProperties::update_solution_fields(f, mesh, cfg);
        const double liquid1 = liquid_mass(f);
        Diagnostics::mass_residual(f, mesh, cfg, 0.001);
        check(std::abs(f["mass_residual"](0, 0)) < 1e-12,
              "liquid residual must include the conservative phase-transfer source");
        const double gas1 = liquid1 * f["dissolved_gas"](0, 0) + f["free_gas_mass"](0, 0);
        check(close(gas1, gas0), "phase exchange must conserve dissolved plus free gas under FBNS");
        check(close(liquid1 * (1 - f["dissolved_gas"](0, 0)), oil0),
              "phase exchange must conserve oil as solution mass changes");
        check(close(liquid1 + f["free_gas_mass"](0, 0), liquid0),
              "phase exchange must conserve total mixture mass");
    }

    cfg.density_model = DensityModel::CONSTANT;
    auto f = fields_for(mesh, cfg);
    FluidProperties::update_gas_state(f, mesh, cfg, 1.0);
    const double equilibrium = FluidProperties::equilibrium_dissolved_gas(1e6, 300, cfg);
    check(f["dissolved_gas"](0, 0) >= equilibrium - 1e-14,
          "stiff release must approach equilibrium without crossing it");
    const double gas_density = FluidProperties::gas_density(1e6, 300, cfg);
    check(close(f["free_gas_mass"](0, 0), gas_density * cfg.c * f["alpha_gas"](0, 0)),
          "void fraction must be derived from the current gas inventory and EOS");

    f["pressure"].fill(2e6);
    cfg.gas_mass_transfer_rate = 0;
    FluidProperties::update_gas_state(f, mesh, cfg, 0.01);
    check(close(f["free_gas_mass"](0, 0),
                FluidProperties::gas_density(2e6, 300, cfg) * cfg.c * f["alpha_gas"](0, 0)),
          "pressure changes must update gas volume without a fictitious relaxation time");
}

void test_separate_absorption(const Mesh& mesh) {
    auto cfg = config();
    cfg.dissolved_gas_species = DissolvedGasSpecies::AIR;
    cfg.dissolved_gas_henry_coeff = 1e-8;
    cfg.gas_absorption_rate = 0;
    auto f = fields_for(mesh, cfg);
    f["pressure"].fill(1e5);
    FluidProperties::update_gas_state(f, mesh, cfg, 0.01);
    const double free = f["free_gas_mass"](0, 0);
    const double alpha = f["alpha_gas"](0, 0);
    const double oil = liquid_mass(f) * (1 - f["dissolved_gas"](0, 0));
    const double gas = liquid_mass(f) * f["dissolved_gas"](0, 0) + free;
    check(free > 0, "zero absorption must still allow supersaturated oil to release air");
    f["pressure"].fill(4e5);
    FluidProperties::update_gas_state(f, mesh, cfg, 0.01);
    check(close(f["free_gas_mass"](0, 0), free), "compression at zero absorption must conserve free gas mass");
    check(close(f["alpha_gas"](0, 0), alpha / 4), "quadrupling ideal-gas pressure must quarter its volume");
    cfg.gas_mass_transfer_rate = 0;
    cfg.gas_absorption_rate = 3;
    const double ml = liquid_mass(f);
    const double md = f["dissolved_gas"](0, 0) * ml;
    const double ceq = 0.004;
    const double transfer = (md - ceq * ml) * (-std::expm1(-(1-ceq) * 3 * 0.01)) / (1-ceq);
    FluidProperties::update_gas_state(f, mesh, cfg, 0.01);
    check(close(f["free_gas_mass"](0, 0), free + transfer),
          "absorption-only kinetics must match the exact conservative reaction");
    check(close(liquid_mass(f) * (1 - f["dissolved_gas"](0, 0)), oil), "absorption must conserve oil");
    check(close(liquid_mass(f) * f["dissolved_gas"](0, 0) + f["free_gas_mass"](0, 0), gas),
          "absorption must conserve dissolved plus free gas");
}

void test_conservative_transport(const Mesh& mesh) {
    auto cfg = config();
    auto f = fields_for(mesh, cfg);
    // A prescribed carrier compression halves concentration. With no face flux,
    // changing capacity cannot create dissolved gas.
    f["h"].fill(2 * cfg.c);
    Communicator comm(mesh);
    comm.update_ghosts(f);
    LinearSystem sys(mesh);
    GasTransport::solve(f, sys, mesh, cfg, 1e-3);
    check(close(f["dissolved_gas"](0, 0), 0.5 * cfg.dissolved_gas_initial, 1e-9),
          "gas transport must conserve old concentration times old liquid mass");
    const auto balance = Diagnostics::mass_balance(f, mesh, cfg, 1e-3);
    check(std::abs(balance.gas_residual) < 1e-9,
          "gas diagnostics must measure the complete stored mass product");
}

void test_carrier_storage(const Mesh& mesh) {
    auto cfg = config();
    cfg.dt = 1e-4;
    for (auto model : {CavitationModel::JFO, CavitationModel::FBNS}) {
        cfg.cavitation_model = model;
        auto f = fields_for(mesh, cfg);
        f["pressure"].fill(cfg.p_cav);
        f["h"].fill(2 * cfg.c);
        f["rho_liquid_solution"].fill(0.8 * cfg.rho);
        Communicator comm(mesh);
        comm.update_ghosts(f);
        LinearSystem sys(mesh);
        if (model == CavitationModel::FBNS) FBNS::solve(f, sys, comm, mesh, cfg);
        else Reynolds::solve_elrod(f, sys, mesh, cfg);
        check(close(f["theta"](0, 0), 0.5 / 1.6, 1e-8),
              "flow must conserve liquid mass when density and film thickness change");
    }
}

void test_closed_bearing_gas_inventory(const Mesh& mesh) {
    auto cfg = config();
    cfg.omega = 200;
    cfg.dt = 1e-5;
    auto f = fields_for(mesh, cfg);
    for (int j = 0; j < mesh.n_z_local; ++j) {
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double angle = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
            f["h"](i, j) *= 1 + 0.3 * std::cos(angle);
            f["dissolved_gas"](i, j) *= 1 + 0.2 * std::sin(angle);
        }
    }
    auto inventories = [&]() {
        std::vector<double> masses(2, 0);
        for (int j = 0; j < mesh.n_z_local; ++j) {
            for (int i = 0; i < mesh.n_theta_local; ++i) {
                const double liquid = f["rho_liquid_solution"](i, j) * f["theta"](i, j) * f["h"](i, j);
                const double dissolved = liquid * f["dissolved_gas"](i, j);
                masses[0] += liquid - dissolved;
                masses[1] += dissolved + f["free_gas_mass"](i, j);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, masses.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return masses;
    };
    Communicator comm(mesh);
    LinearSystem sys(mesh);
    const auto initial = inventories();
    double worst = 0;
    for (int step = 0; step < 20; ++step) {
        FluidProperties::update_solution_fields(f, mesh, cfg);
        comm.update_ghosts(f);
        for (auto& [name, field] : f) field->store_old_time();
        FBNS::solve(f, sys, comm, mesh, cfg);
        comm.update_ghosts(f);
        GasTransport::solve(f, sys, mesh, cfg, cfg.dt);
        FluidProperties::update_gas_state(f, mesh, cfg, cfg.dt);
        const auto current = inventories();
        worst = std::max({worst, std::abs(current[0] / initial[0] - 1),
                                 std::abs(current[1] / initial[1] - 1)});
    }
    check(worst < 1e-8, "closed moving film must conserve oil and gas through flow, transport and reaction");
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Closed bearing worst relative inventory drift: " << worst << '\n';
}

void test_fbns_reservoir_transport(const Mesh& mesh) {
    auto cfg = config();
    cfg.omega = 200;
    cfg.dt = 1e-5;
    cfg.inlets.push_back({InletConfig::Type::GROOVE, 90, 0, 30, 2e5});
    cfg.bc_z_south_type = cfg.bc_z_north_type = BCType::DIRICHLET;
    cfg.bc_z_south_val = cfg.bc_z_north_val = cfg.p_cav;
    auto f = fields_for(mesh, cfg);
    f["theta"].fill(1);
    f["rho"].fill(cfg.rho);
    f["pressure"].fill(cfg.p_cav);
    FilmThickness::compute_inlet_indicator(f["inlet_indicator"], mesh, cfg);
    for (auto& [name, field] : f) field->store_old_time();
    Communicator comm(mesh);
    comm.update_ghosts(f);
    LinearSystem sys(mesh);
    FBNS::solve(f, sys, comm, mesh, cfg);
    comm.update_ghosts(f);
    GasTransport::solve(f, sys, mesh, cfg, cfg.dt);
    double error = 0;
    for (int j = 0; j < mesh.n_z_local; ++j)
        for (int i = 0; i < mesh.n_theta_local; ++i)
            error = std::max(error, std::abs(f["dissolved_gas"](i, j) - cfg.dissolved_gas_initial));
    check(error < 1e-9, "FBNS reservoir flow must preserve a uniform feed concentration throughout the domain");
    check(std::abs(Diagnostics::mass_balance(f, mesh, cfg, cfg.dt).gas_residual) < 1e-8,
          "open gas transport must account for the imposed reservoir source and axial fluxes");
}

void test_gas_volume_diagnostics(const Mesh& mesh) {
    const auto cfg = config();
    auto f = fields_for(mesh, cfg);
    f["theta"].fill(0.9);
    f["rho_gas"].fill(2);
    f["free_gas_mass"].fill(0.3 * 2 * cfg.c);
    f["alpha_gas"].fill(0.1);
    const auto balance = Diagnostics::mass_balance(f, mesh, cfg, cfg.dt);
    check(close(balance.gas_volume_excess, 0.2), "volume diagnostics must see gas hidden by alpha clipping");
    check(close(balance.gas_eos_defect, 0.2), "EOS diagnostics must compare stored alpha with gas inventory");
}

void test_species_and_mass_fraction() {
    for (const char* species : {"NITROGEN", "AR", "PROAPNE"}) {
        SimulationConfig cfg;
        cfg.load_from_text(std::string("dissolved_gas_species = ") + species + "\n");
        std::vector<std::string> errors, warnings;
        cfg.validate(errors, warnings);
        check(!cfg.parse_errors.empty(), "unsupported gas species must not silently become propane");
    }
    auto cfg = config();
    cfg.dissolved_gas_max = 1.1;
    std::vector<std::string> errors, warnings;
    cfg.validate(errors, warnings);
    check(!errors.empty(), "dissolved gas mass fractions above one must be rejected");
    cfg.coupling_tolerance = 1e-6;
    Coupling::Controller controller(cfg);
    check(!controller.record(0, 0, 0, 1e-3), "steady gas changes must prevent a false coupled convergence");
    check(controller.record(1, 0, 0, 0), "converged gas and thermal fields should pass the coupling gate");
}

void test_tabulated_viscosity(const Mesh& mesh) {
    auto cfg = config();
    cfg.liquid_viscosity_model = LiquidViscosityModel::TABLE;
    cfg.viscosity_table_2d.p = {1e5, 1e6};
    cfg.viscosity_table_2d.t = {300, 350};
    cfg.viscosity_table_2d.values = {10e-6, 20e-6, 5e-6, 10e-6};
    auto f = fields_for(mesh, cfg);
    FluidProperties::update_solution_fields(f, mesh, cfg);
    const double reference = cfg.rho * cfg.viscosity_table_2d.interpolate(1e6, 300, 0);
    check(close(FluidProperties::liquid_solution_viscosity(0.002, 1e6, 300, cfg), reference),
          "cross-film property calls must use the same tabulated viscosity as the mean film");
    check(close(f["mu"](0, 0), reference), "tabulated kinematic viscosity must convert to dynamic viscosity");
}
}  // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    {
        const auto cfg = config();
        Mesh mesh(cfg);
        test_phase_exchange(mesh);
        test_separate_absorption(mesh);
        test_conservative_transport(mesh);
        test_carrier_storage(mesh);
        test_closed_bearing_gas_inventory(mesh);
        test_fbns_reservoir_transport(mesh);
        test_gas_volume_diagnostics(mesh);
        test_species_and_mass_fraction();
        test_tabulated_viscosity(mesh);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::cout << "Physics audit failures: " << failures << '\n';
    PetscFinalize();
    return failures ? 1 : 0;
}
