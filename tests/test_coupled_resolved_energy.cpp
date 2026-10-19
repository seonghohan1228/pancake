#include <petscsys.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "coupled_cavitation.hpp"
#include "cross_film.hpp"
#include "resolved_cross_film.hpp"
#include "fluid_properties.hpp"

namespace {
void require(bool value,const char* message) {if(!value)throw std::runtime_error(message);}
void near(double a,double b,double tol,const char* message) {
    if(std::abs(a-b)>tol) {std::cerr<<message<<": "<<a<<" vs "<<b<<'\n';throw std::runtime_error(message);}
}
SimulationConfig config() {
    SimulationConfig c;
    c.n_theta_global=24;c.n_z_global=4;c.e=0;c.R=.0215;c.L=.04;c.c=42.5e-6;
    c.solution_mode=SolutionMode::STEADY_STATE;c.cavitation_model=CavitationModel::COUPLED;
    c.fluid_property_model=FluidPropertyModel::SINGLE_PHASE;
    c.temperature_model=TemperatureModel::ENERGY_EQUATION;
    c.cross_film_thermal_model=CrossFilmModel::FINITE_VOLUME;c.cross_film_layers=16;
    c.mixture_viscosity_model=MixtureViscosityModel::DUKLER;
    c.groove_mixing_model=GrooveMixingModel::PINNED;
    c.mu=.04;c.mu_gas=1.8e-5;c.rho=920;c.liquid_specific_heat=2010;c.rho_cp=1; // explicit cp must win
    c.dissolved_heat_capacity=1656.929131763835;c.gas_heat_capacity_cv=717.5;
    c.desorption_internal_energy=192486.2474970321;
    c.thermal_conductivity=.13;c.gas_thermal_conductivity=.013;
    c.temperature_initial=c.temperature_reference=c.property_reference_temperature=300;
    c.thermal_ambient_temperature=c.bush_outer_temperature=c.shaft_ambient_temperature=300;
    c.bush_thermal_model=BushThermalModel::CONDUCTION_3D;c.bush_radial_cells=4;c.bush_axial_cells=4;
    c.bush_outer_radius=.0425;c.bush_thermal_conductivity=15;c.bush_outer_heat_transfer=65;
    c.bush_outer_boundary=BushOuterBoundary::CONVECTION;
    c.shaft_ambient_heat_transfer=1300;c.shaft_ambient_area_ratio=.5375;
    c.omega=100;c.p_cav=1;c.bc_z_south_val=c.bc_z_north_val=1e5;
    c.linear_rtol=1e-12;c.dissolved_gas_species=DissolvedGasSpecies::AIR;
    return c;
}
void fields(Fields& f,const Mesh& mesh,const SimulationConfig& c,double concentration) {
    for(auto name:{"pressure","h","theta","dissolved_gas","free_gas_mass","temperature","alpha_gas","rho","mu","rho_liquid_solution","mu_liquid_solution"})f.add(name,mesh);
    f["pressure"].fill(.5*(c.bc_z_south_val+c.bc_z_north_val));f["h"].fill(c.c);f["theta"].fill(1);
    f["dissolved_gas"].fill(concentration);f["free_gas_mass"].fill(0);f["temperature"].fill(c.temperature_initial);
    f["rho"].fill(c.rho);f["mu"].fill(c.mu);CrossFilm::register_fields(f,mesh,c);
}
void uniform_limit() {
    auto c=config();Mesh mesh(c);Communicator comm(mesh);Fields f;fields(f,mesh,c,0);
    CoupledCavitation::Solver solver(mesh,comm,false);solver.initialize(f,c);
    CrossFilm::update(f,mesh,c);
    auto s=solver.solve(f,c,.01,false);require(s.converged,"zero-gas isothermal hydraulic limit");
    ResolvedCrossFilm::solve(f,mesh,c);
    const auto wall=ResolvedCrossFilm::diagnostics(f,mesh,c);
    const auto ledger=ResolvedCrossFilm::balance(f,mesh,c);
    const double heat=c.mu*std::pow(c.omega*c.R,2)/c.c;
    const double gs=c.shaft_ambient_heat_transfer*c.shaft_ambient_area_ratio;
    const double gb=1/(c.R*std::log(c.bush_outer_radius/c.R)/c.bush_thermal_conductivity+c.R/(c.bush_outer_radius*c.bush_outer_heat_transfer));
    const double g=c.thermal_conductivity/c.c,det=(gs+g)*(gb+g)-g*g;
    near(wall.t_shaft,300+heat*.5*(gb+2*g)/det,2e-6,"zero-gas uniform shaft analytic balance");
    near(wall.t_bush_mean,300+heat*.5*(gs+2*g)/det,2e-6,"zero-gas 3D housing analytic balance");
    require(std::abs(ledger.residual_W)/ledger.scale_W<1e-8,"zero-gas conservative energy balance");
    // Each film-solid interface enters the two equations with opposite signs.
    const double surface=2*std::acos(-1.0)*c.R*c.L;
    near(wall.q_shaft_net,gs*surface*(wall.t_shaft-300),1e-7,"shaft interface cancellation");
    near(wall.q_bush_net,gb*surface*(wall.t_bush_mean-300),1e-7,"housing interface cancellation");
    const auto u=CoupledCavitation::internal_energies(310,1e5,c);
    near(u[0],20100,1e-8,"liquid_specific_heat overrides legacy rho_cp/rho");
    near(u[1],10*c.dissolved_heat_capacity-c.desorption_internal_energy,1e-8,"explicit dissolved caloric binding analogue");
    if(mesh.rank==0)std::cout<<"zero_gas_energy_residual_W="<<ledger.residual_W<<" shaft_K="<<wall.t_shaft<<" housing_K="<<wall.t_bush_mean<<'\n';
}
double phase_transport(double binding) {
    auto c=config();c.fluid_property_model=FluidPropertyModel::TWO_PHASE;c.omega=0;
    c.bc_z_south_val=2e5;c.bc_z_north_val=1e5;
    c.bc_z_south_thermal=c.bc_z_north_thermal=ThermalInflowMode::CONSTANT;
    c.bc_z_south_temperature=c.bc_z_north_temperature=300;
    c.bc_z_south_dissolved_gas=c.bc_z_north_dissolved_gas=.03;
    c.dissolved_gas_initial=.03;c.dissolved_gas_henry_coeff=1e-7;
    c.desorption_internal_energy=binding;c.gas_mass_transfer_rate=c.gas_absorption_rate=0;
    Mesh mesh(c);Communicator comm(mesh);Fields f;fields(f,mesh,c,.03);
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<mesh.n_z_local;++j)
        f["pressure"](i,j)=c.bc_z_south_val+(c.bc_z_north_val-c.bc_z_south_val)*(j+.5)/mesh.n_z_local;
    CoupledCavitation::Solver solver(mesh,comm,false);solver.initialize(f,c);
    for(int iteration=0;iteration<40;++iteration) {
        CrossFilm::update(f,mesh,c);
        const auto s=solver.solve(f,c,.01,true);require(s.converged,"steady phase-equilibrium hydraulic convergence");
        require(std::abs(s.oil_balance)<1e-8 && std::abs(s.species_balance)<1e-8,"oil and total-species conservation");
        ResolvedCrossFilm::solve(f,mesh,c);
        const auto wall=ResolvedCrossFilm::diagnostics(f,mesh,c);
        if(wall.profile_residual<1e-10)break;
        require(iteration<39,"phase-energy Picard convergence");
    }
    const auto ledger=ResolvedCrossFilm::balance(f,mesh,c);
    require(std::abs(ledger.residual_W)/std::max(ledger.scale_W,1e-12)<1e-8,"phase transport and pressure-work energy balance");
    double free=0,mean=0;
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<mesh.n_z_local;++j) {
        free+=f["free_gas_mass"](i,j);mean+=f["temperature"](i,j);
        for(const char* axis:{"theta","z"}) {
            double sum=0;
            for(int phase=0;phase<3;++phase)for(const char* part:{"pressure","moving","diffusion"})
                sum+=f["coupled_"+std::string(part)+"_"+std::to_string(phase)+"_"+axis](i,j);
            near(sum,f["mixture_mass_flux_"+std::string(axis)](i,j),1e-14,"resolved phase fluxes exactly match hydraulic transport");
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&free,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&mean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    require(free>0,"instantaneous equilibrium must release gas despite zero finite transfer rates");
    mean/=mesh.n_theta_global*mesh.n_z_global;
    if(mesh.rank==0)std::cout<<"binding_J_kg="<<binding<<" mean_K="<<mean<<" phase_energy_residual_W="<<ledger.residual_W<<'\n';
    return mean;
}
}
int main(int argc,char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int failed=0;
    try {uniform_limit();const double unbound=phase_transport(0),bound=phase_transport(192486.2474970321);require(bound<unbound,"desorption binding energy must cool the conserved phase-energy solution");}
    catch(const std::exception& e) {std::cerr<<e.what()<<'\n';failed=1;}
    MPI_Allreduce(MPI_IN_PLACE,&failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    PetscFinalize();return failed;
}
