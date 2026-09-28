#include <cmath>
#include <iostream>
#include <petsc.h>
#include "steady_equilibrium.hpp"
#include "film_thickness.hpp"
#include "coupling.hpp"

void check(bool condition, const char* message) {
    if (!condition) { std::cerr << message << '\n'; MPI_Abort(MPI_COMM_WORLD, 1); }
}

JournalMotion::BearingState recover(double dt) {
    SimulationConfig cfg;
    cfg.motion_model = MotionModel::STEADY_EQUILIBRIUM;
    cfg.solution_mode = SolutionMode::STEADY_STATE;
    cfg.cavitation_model = CavitationModel::FBNS;
    cfg.fluid_property_model = FluidPropertyModel::SINGLE_PHASE;
    cfg.n_theta_global = 64; cfg.n_z_global = 4;
    cfg.R = .05; cfg.c = 1.52e-4; cfg.L = .08; cfg.mu = .0277; cfg.rho = 860;
    cfg.omega = 418.879; cfg.p_cav = 1e5; cfg.dt = dt;
    cfg.fbns_tol = 1e-10; cfg.fbns_max_iters = 100; cfg.linear_rtol = 1e-12;
    cfg.bc_z_south_val = cfg.bc_z_north_val = cfg.p_cav;
    cfg.steady_force_tolerance = 1e-7;
    InletConfig groove; groove.type = InletConfig::Type::GROOVE;
    groove.theta = 90; groove.size = 18; groove.p_supply = 1.713e5;
    cfg.inlets = {groove};
    Mesh mesh(cfg); Communicator comm(mesh); LinearSystem sys(mesh, linear_solver_options_from_config(cfg));
    Fields fields;
    fields.add("pressure",mesh).fill(cfg.p_cav); fields.add("theta",mesh).fill(1);
    fields.add("h",mesh).fill(cfg.c); fields.add("dh_dt",mesh).fill(0);
    fields.add("rho",mesh).fill(cfg.rho); fields.add("rho_liquid_solution",mesh).fill(cfg.rho);
    fields.add("mu",mesh).fill(cfg.mu); fields.add("temperature",mesh).fill(300);
    for (auto name : {"cross_film_f0","cross_film_f1","cross_film_f2","cross_film_couette_h"}) fields.add(name,mesh);
    cfg.cross_film_thermal_model = CrossFilmModel::POLYNOMIAL;
    JournalMotion::BearingState target; target.x = -.25*cfg.c; target.y = .4*cfg.c;
    auto geometry = [&](const JournalMotion::BearingState& state) {
        FilmThickness::compute_moving_bearing(fields["h"],fields["dh_dt"],mesh,cfg,state);
        for (int i=0; i<mesh.n_theta_local; ++i) for (int j=0; j<mesh.n_z_local; ++j) {
            const double h=fields["h"](i,j);
            fields["cross_film_f0"](i,j)=h/cfg.mu;
            fields["cross_film_f1"](i,j)=h*h/(2*cfg.mu);
            fields["cross_film_f2"](i,j)=h*h*h/(12*cfg.mu);
            fields["cross_film_couette_h"](i,j)=h;
        }
        comm.update_ghosts(fields);
    };
    geometry(target);
    check(FBNS::solve(fields,sys,comm,mesh,cfg).converged,"target pressure converged");
    const auto force = Reynolds::calculate_macroscopic_properties(fields,mesh,cfg);
    cfg.external_load_x = -force.fluid_x; cfg.external_load_y = -force.fluid_y;
    auto state=target; state.x += .09*cfg.c; state.y -= .07*cfg.c;
    state.vx=7; state.vy=-5;
    geometry(state);
    fields["h"].store_old_time(); fields["pressure"].store_old_time(); fields["theta"].store_old_time();
    const auto displaced = state;
    const auto displaced_h = fields["h"].data;
    cfg.steady_position_fd_step = 2.0;
    const auto rejected = SteadyEquilibrium::update(fields,sys,comm,mesh,cfg,state);
    check(!rejected.accepted && state.x==displaced.x && state.y==displaced.y,"invalid Jacobian probes retain base position");
    check(fields["h"].data==displaced_h,"invalid probes restore physical and ghost geometry");
    cfg.steady_position_fd_step = .001;
    auto first=SteadyEquilibrium::update(fields,sys,comm,mesh,cfg,state);
    check(first.force_residual > cfg.steady_force_tolerance,"one iteration must not falsely close displaced bearing");
    cfg.coupling_max_iters=1;
    Coupling::Controller control(cfg);
    check(!control.record(0,0,0,0,first.force_residual < cfg.steady_force_tolerance),"thermal stationarity cannot bypass force closure");
    control.finish(1); check(!control.status().converged && control.status().cap_exhausted,"small cap fails honestly");
    SteadyEquilibrium::Stats stats;
    int accepted_updates = first.accepted ? 1 : 0;
    for (int iteration=0; iteration<60; ++iteration) {
        stats=SteadyEquilibrium::update(fields,sys,comm,mesh,cfg,state);
        if (stats.accepted) ++accepted_updates;
        if (stats.force_residual < cfg.steady_force_tolerance) break;
    }
    check(stats.force_residual < cfg.steady_force_tolerance,"recover displaced known equilibrium");
    check(accepted_updates > 1,"position updates repeatedly inside steady iteration");
    check(std::hypot(state.x-target.x,state.y-target.y)/cfg.c < 1e-6,"recover target Cartesian position");
    check(state.vx==0 && state.vy==0 && state.vz==0,"no physical velocities in algebraic equilibrium");
    check(FBNS::solve(fields,sys,comm,mesh,cfg).converged,"final pressure re-solve converges at final geometry");
    const auto final_force=Reynolds::calculate_macroscopic_properties(fields,mesh,cfg);
    check(SteadyEquilibrium::force_residual(final_force,cfg)<cfg.steady_force_tolerance,"independent force at final geometry");
    for (int i=0; i<mesh.n_theta_local; ++i) for(int j=0; j<mesh.n_z_local; ++j) {
        const double h=fields["h"](i,j);
        check(fields["dh_dt"](i,j)==0,"no squeeze term");
        check(std::abs(fields["cross_film_f2"](i,j)/(h*h*h/(12*cfg.mu))-1)<1e-12,"probe and accepted geometry retain Dowson powers");
        check(fields["temperature"](i,j)==300,"pressure probes do not solve energy");
    }
    return state;
}
int main(int argc,char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    { const auto a=recover(1e-8),b=recover(1e3);
      check(std::hypot(a.x-b.x,a.y-b.y)<1e-13,"equilibrium independent of physical timestep"); }
    int rank=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0) std::cout << "PASS: steady equilibrium recovery, frozen profile, force closure, no dt, honest cap\n";
    PetscFinalize();
}

