#include <petsc.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "steady_equilibrium.hpp"
#include "film_thickness.hpp"

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void recover(bool gas, bool chart = false) {
    SimulationConfig cfg;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.motion_model = MotionModel::STEADY_EQUILIBRIUM;
    cfg.cavitation_model = CavitationModel::COUPLED;
    cfg.fluid_property_model = gas ? FluidPropertyModel::TWO_PHASE : FluidPropertyModel::SINGLE_PHASE;
    cfg.mixture_viscosity_model = MixtureViscosityModel::DUKLER;
    cfg.coupled_phase_equilibrium = gas;
    cfg.n_theta_global = 24; cfg.n_z_global = 4;
    cfg.R = .05; cfg.c = 1.52e-4; cfg.L = .08;
    cfg.mu = .0277; cfg.rho = 860; cfg.omega = 10; cfg.p_cav = 0;
    cfg.bc_z_south_val = cfg.bc_z_north_val = 1e6;
    cfg.dissolved_gas_initial = gas ? .03 : 0;
    cfg.dissolved_gas_henry_coeff = 1e-6;
    if (!gas) cfg.solubility_model = SolubilityModel::CURVE_TABLE; // Disabled exchange must not query an absent table.
    if (chart) {
        cfg.solubility_model = SolubilityModel::CURVE_TABLE;
        cfg.solution_table_clamp = false;
        cfg.solubility_curve_table = {{1e5,2e5}, {{{200,.02},{500,.02}},{{200,.04},{500,.04}}}};
        require(std::abs(cfg.saturation_pressure(.03,cfg.temperature_initial)-1.5e5)<1e-8,
                "strict chart bubble pressure remains defined for high-pressure undersaturated liquid");
    }
    cfg.steady_force_tolerance = 1e-7;
    InletConfig groove; groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90; groove.size = 30; groove.p_supply = 1.2e6;
    cfg.inlets = {groove};
    Mesh mesh(cfg); Communicator comm(mesh);
    LinearSystem sys(mesh, linear_solver_options_from_config(cfg));
    Fields fields;
    for (auto name : {"pressure","h","dh_dt","theta","temperature","dissolved_gas","free_gas_mass",
                      "alpha_gas","rho","mu","rho_liquid_solution","mu_liquid_solution"}) fields.add(name,mesh);
    fields["pressure"].fill(1e6); fields["h"].fill(cfg.c); fields["theta"].fill(1);
    fields["temperature"].fill(cfg.temperature_initial);
    fields["dissolved_gas"].fill(cfg.dissolved_gas_initial);
    CoupledCavitation::Solver solver(mesh,comm);
    JournalMotion::BearingState target; target.x = -.25*cfg.c; target.y = .4*cfg.c;
    auto geometry = [&](const JournalMotion::BearingState& state) {
        FilmThickness::compute_moving_bearing(fields["h"], fields["dh_dt"], mesh, cfg, state);
        fields["dh_dt"].fill(0); comm.update_ghosts(fields);
    };
    geometry(target); solver.initialize(fields,cfg);
    require(solver.solve(fields,cfg,cfg.dt,gas).converged,"known coupled equilibrium pressure");
    const auto target_forces = Reynolds::calculate_macroscopic_properties(fields,mesh,cfg);
    cfg.external_load_x = -target_forces.fluid_x; cfg.external_load_y = -target_forces.fluid_y;
    auto state = target; state.x += .03*cfg.c; state.y -= .02*cfg.c;
    geometry(state);
    require(solver.solve(fields,cfg,cfg.dt,gas).converged,"displaced coupled base pressure");
    CoupledCavitation::Solver::store_old_time(fields);
    const auto base_h = fields["h"].data;
    const auto old_oil = fields["oil_mass"].old_data;
    const auto old_dissolved = fields["dissolved_mass"].old_data;
    const auto old_free = fields["free_gas_mass"].old_data;
    const auto base_state = state;
    cfg.steady_position_fd_step = 2;
    const auto rejected = SteadyEquilibrium::update(fields,sys,comm,mesh,cfg,state,&solver);
    require(!rejected.accepted && state.x == base_state.x && state.y == base_state.y,"invalid coupled probes retain position");
    require(fields["h"].data == base_h,"rejected coupled probes restore geometry including ghosts");
    require(fields["oil_mass"].old_data == old_oil && fields["dissolved_mass"].old_data == old_dissolved &&
            fields["free_gas_mass"].old_data == old_free,"rejected coupled probes preserve every conserved history");
    cfg.steady_position_fd_step = .001;
    SteadyEquilibrium::Stats stats;
    for (int iteration=0; iteration<30; ++iteration) {
        stats = SteadyEquilibrium::update(fields,sys,comm,mesh,cfg,state,&solver);
        if (stats.force_residual < cfg.steady_force_tolerance) break;
    }
    require(stats.coupled.converged && stats.force_residual < cfg.steady_force_tolerance,"coupled load closure");
    require(std::hypot(state.x-target.x,state.y-target.y)/cfg.c < 1e-6,"recover known coupled Cartesian position");
    require(std::abs(stats.coupled.oil_balance) < 1e-10 && std::abs(stats.coupled.species_balance) < 1e-10,
            "accepted coupled equilibrium conserves both inventories");
    for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
        require(fields["temperature"](i,j) == cfg.temperature_initial && fields["dh_dt"](i,j) == 0,
                "coupled pressure probes freeze thermal state and squeeze");
}
}

int main(int argc, char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int failed = 0;
    try { recover(false); recover(true); recover(true,true); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; failed = 1; }
    int any_failed = 0; MPI_Allreduce(&failed,&any_failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (!any_failed) std::cout << "PASS: coupled load closure, gas-off/on, frozen energy, conservative rollback\n";
    PetscFinalize(); return any_failed;
}
