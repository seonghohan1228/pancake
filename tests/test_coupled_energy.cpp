#include <petscsys.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "coupled_cavitation.hpp"
#include "fluid_properties.hpp"
#include "r410a.hpp"

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double value,double expected,double tolerance,const char* message) {
    if (std::abs(value-expected)>tolerance) {
        std::cerr << message << ": " << value << " vs " << expected << '\n';
        throw std::runtime_error(message);
    }
}
}
int main(int argc,char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int failed=0;
    try {
        SimulationConfig cfg;
        cfg.n_theta_global=16; cfg.n_z_global=4; cfg.e=0;
        cfg.cavitation_model=CavitationModel::COUPLED;
        cfg.temperature_model=TemperatureModel::ENERGY_EQUATION;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        cfg.rho=850; cfg.rho_cp=1.7e6; cfg.mu=.03; cfg.mu_gas=1.8e-5;
        cfg.mixture_viscosity_model=MixtureViscosityModel::DUKLER;
        cfg.omega=200; cfg.p_cav=1; cfg.dissolved_gas_species=DissolvedGasSpecies::AIR;
        cfg.temperature_initial=cfg.temperature_reference=cfg.property_reference_temperature=313.15;
        cfg.dissolved_heat_capacity=2000; cfg.gas_heat_capacity_cv=717.5;
        cfg.desorption_internal_energy=1e5; cfg.property_reference_pressure=1e5;
        Mesh mesh(cfg); Communicator comm(mesh); Fields f;
        for (auto name : {"pressure","h","theta","dissolved_gas","free_gas_mass","temperature","alpha_gas"}) f.add(name,mesh);
        CoupledCavitation::Solver solver(mesh,comm,true);
        const auto reset=[&](double theta,double c,double pg) {
            f["pressure"].fill(pg+cfg.p_cav); f["h"].fill(cfg.c); f["theta"].fill(theta);
            f["temperature"].fill(cfg.temperature_initial); f["dissolved_gas"].fill(c);
            const double rg=cfg.dissolved_gas_species==DissolvedGasSpecies::R410A ?
                R410A::gas_density(pg,cfg.temperature_initial) : pg/(287.05*cfg.temperature_initial);
            f["free_gas_mass"].fill(rg*cfg.c*(1-theta)); solver.initialize(f,cfg);
        };
        const double dt=.01;
        // Sealed constant-gap pure oil has a pressure nullspace. Supply one
        // full-width pressure patch; its net mass source is exactly zero here.
        InletConfig groove; groove.type=InletConfig::Type::GROOVE; groove.theta=90;
        groove.size=30; groove.p_supply=1e5+cfg.p_cav; cfg.inlets.push_back(groove);
        reset(1,0,1e5);
        auto result=solver.solve(f,cfg,dt);
        require(result.converged,"Adiabatic Couette solve");
        const double heat=cfg.mu*std::pow(cfg.omega*cfg.R,2)/cfg.c;
        const double exact=cfg.temperature_initial+dt*heat/(cfg.rho_cp*cfg.c);
        near(f["temperature"](0,0),exact,1e-7,"Exact Couette temperature rise");
        near(result.energy_balance,0,1e-7,"Global heat balance");
        const double once=f["temperature"](0,0);
        result=solver.solve(f,cfg,dt);
        require(result.converged,"Repeated thermal nonlinear solve");
        near(f["temperature"](0,0),once,1e-8,"No repeated physical heat advance");
        cfg.journal_heat_transfer=200; cfg.bearing_heat_transfer=300;
        cfg.journal_wall_temperature=cfg.bearing_wall_temperature=300;
        reset(1,0,1e5); result=solver.solve(f,cfg,dt);
        require(result.converged,"Couette wall cooling solve");
        const double capacity=cfg.rho_cp*cfg.c;
        near(f["temperature"](0,0),(capacity*cfg.temperature_initial+dt*(heat+500*300))/(capacity+dt*500),1e-7,
            "Independent implicit wall-cooling temperature");
        cfg.journal_heat_transfer=cfg.bearing_heat_transfer=0; cfg.omega=0; cfg.inlets.clear();
        cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE;
        cfg.gas_mass_transfer_rate=0; cfg.gas_absorption_rate=0;
        reset(.5,0,1e5);
        const double oil=f["oil_mass"](0,0), gas=f["free_gas_mass"](0,0);
        const double new_h=.75*cfg.c, vg=new_h-oil/cfg.rho;
        const double cv=oil*cfg.rho_cp/cfg.rho+gas*cfg.gas_heat_capacity_cv;
        const double dh=new_h-cfg.c;
        const double compressed=(cv*cfg.temperature_initial-cfg.p_cav*dh)/(cv+gas*287.05*dh/vg);
        f["h"].fill(new_h); result=solver.solve(f,cfg,dt);
        require(result.converged,"Adiabatic bubble compression");
        near(f["temperature"](0,0),compressed,1e-7,"Independent implicit compression temperature");
        near(f["pressure"](0,0),cfg.p_cav+gas*287.05*compressed/vg,.001,"Hot compressed bubble pressure");
        require(compressed>cfg.temperature_initial,"Compression must heat without absorption");
        cfg.fluid_property_model=FluidPropertyModel::CONSTANT;
        cfg.omega=200; cfg.mu_gas=0; cfg.oil_vapor_temperature_coeff=3000;
        reset(.8,0,0);
        result=solver.solve(f,cfg,dt);
        require(result.converged,"Heated pure-vapor cavity");
        near(f["temperature"](0,0),exact,1e-7,"Liquid capacity and shear scale consistently in a vapor cavity");
        near(f["pressure"](0,0),cfg.oil_vapor_pressure(exact),1e-6,"Oil vapor pressure follows film temperature");
        cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE;
        cfg.omega=0; cfg.mu_gas=1.8e-5; cfg.oil_vapor_temperature_coeff=0;
        cfg.dissolved_gas_henry_coeff=1e-8; cfg.gas_mass_transfer_rate=100;
        reset(.9,.002,1e5);
        const double old_energy=f["internal_energy"](0,0);
        result=solver.solve(f,cfg,dt);
        require(result.converged,"Adiabatic desorption");
        near(f["internal_energy"](0,0),old_energy,1e-7,"Closed desorption conserves internal energy");
        require(f["temperature"](0,0)<cfg.temperature_initial,"Positive desorption energy cools mixture");
        require(f["gas_mass_transfer"](0,0)>0,"Published release rate must be nonzero");
        cfg.fluid_property_model=FluidPropertyModel::CONSTANT;
        cfg.gas_mass_transfer_rate=0; cfg.thermal_conductivity=0;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::DIRICHLET;
        cfg.bc_z_south_val=2e5; cfg.bc_z_north_val=1e5;
        cfg.bc_z_south_thermal=cfg.bc_z_north_thermal=ThermalInflowMode::CONSTANT;
        cfg.bc_z_south_temperature=300; cfg.bc_z_north_temperature=340;
        for (bool reverse : {false,true}) {
            if (reverse) std::swap(cfg.bc_z_south_val,cfg.bc_z_north_val);
            reset(1,0,1.5e5);
            result=solver.solve(f,cfg,dt);
            require(result.converged,"Bidirectional thermal reservoir solve");
            const double velocity=cfg.c*cfg.c*1e5/(12*cfg.mu*cfg.L);
            const double courant=dt*velocity/mesh.get_d_z();
            const double rise=dt*velocity*1e5/(cfg.L*cfg.rho_cp);
            double expected=reverse ? 340 : 300;
            for (int index=0;index<mesh.n_z_local;++index) {
                expected=(cfg.temperature_initial+courant*expected+rise)/(1+courant);
                const int j=reverse ? mesh.n_z_local-index-1 : index;
                near(f["temperature"](0,j),expected,1e-7,"Independent inlet advection and pressure dissipation");
            }
            near(result.energy_balance,0,1e-7,"Axial energy conservation");
        }
        // Feed temperature enters an enthalpy source; it must never pin drainage.
        cfg.bc_z_south_val=cfg.bc_z_north_val=2e5;
        cfg.bc_z_south_temperature=cfg.bc_z_north_temperature=cfg.temperature_initial;
        for (auto geometry : {InletConfig::Type::GROOVE,InletConfig::Type::CIRCULAR}) {
            InletConfig inlet; inlet.type=geometry; inlet.theta=90; inlet.z=cfg.L/2;
            inlet.size=geometry==InletConfig::Type::GROOVE ? 30 : cfg.L*.3;
            inlet.feeds_fresh_oil=true; inlet.t_supply=280;
            for (bool feeding : {false,true}) {
                inlet.p_supply=feeding ? 3e5 : 1e5; cfg.inlets={inlet}; reset(1,0,2e5);
                result=solver.solve(f,cfg,dt);
                require(result.converged,"Groove/hole thermal transport");
                near(result.energy_balance,0,1e-7,"Pressure patch enthalpy conservation");
                if (!feeding) for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
                    require(f["temperature"](i,j)>=cfg.temperature_initial-1e-7,"Drain must not impose its cold supply temperature");
            }
        }
        cfg.inlets.clear(); cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::DIRICHLET;
        cfg.bc_z_south_val=2e5; cfg.bc_z_north_val=1e5;
        cfg.fluid_property_model=FluidPropertyModel::SINGLE_PHASE;
        cfg.liquid_viscosity_model=LiquidViscosityModel::EMPIRICAL;
        cfg.viscosity_temperature_coeff=2500;
        cfg.bc_z_north_temperature=350; reset(1,0,1.5e5);
        result=solver.solve(f,cfg,dt); require(result.converged,"Variable-viscosity outflow");
        const auto reference_pressure=f["pressure"].data, reference_temperature=f["temperature"].data;
        cfg.bc_z_north_temperature=450; reset(1,0,1.5e5);
        result=solver.solve(f,cfg,dt); require(result.converged,"Outflow reservoir sensitivity");
        for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j) {
            const int index=(j+f["temperature"].n_ghost)*f["temperature"].stride_theta+i+f["temperature"].n_ghost;
            near(f["temperature"](i,j),reference_temperature[index],1e-8,"External outflow temperature must not change hydraulic resistance");
            near(f["pressure"](i,j),reference_pressure[index],1e-6,"External outflow temperature must not change pressure");
        }
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        cfg.liquid_viscosity_model=LiquidViscosityModel::CONSTANT;
        cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE; cfg.gas_mass_transfer_rate=100;
        cfg.thermal_conductivity=.13;
        // Refrigerant thermal closure uses real caloric properties, not air cp.
        cfg.dissolved_gas_species=DissolvedGasSpecies::R410A;
        cfg.dissolved_gas_henry_coeff=1e-7;
        for (double t : {293.15,313.15,333.15,353.15}) {
            cfg.temperature_initial=cfg.temperature_reference=cfg.property_reference_temperature=t;
            reset(.9,.1,6e5); const double before=f["internal_energy"](0,0);
            result=solver.solve(f,cfg,.001);
            require(result.converged,"R410A thermal release at 20--80 C");
            near(f["internal_energy"](0,0),before,1e-7,"R410A closed thermal energy balance");
            near(f["theta"](0,0)+f["alpha_gas"](0,0),1,1e-8,"R410A thermal volume closure");
        }
        // Gas-free oil may flow above the pure-refrigerant dew pressure. The
        // Jacobian's tiny trial gas perturbation is not an accepted free phase.
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::DIRICHLET;
        cfg.bc_z_south_val=2e6; cfg.bc_z_north_val=1.8e6;
        cfg.temperature_initial=cfg.temperature_reference=cfg.property_reference_temperature=293.15;
        cfg.bc_z_south_temperature=cfg.bc_z_north_temperature=293.15;
        cfg.dissolved_gas_initial=.1;
        reset(1,.1,6e5); f["pressure"].fill(1.9e6);
        result=solver.solve(f,cfg,.0001,true);
        require(result.converged,"Gas-free R410A solution above the dew pressure must permit Newton perturbations");
        for (int j=0;j<mesh.n_z_local;++j) {
            near(f["pressure"](0,j),2e6-2e5*(j+.5)/mesh.n_z_local,.01,"High-pressure dissolved solution has the linear pressure profile");
            near(f["free_gas_mass"](0,j),0,1e-10,"High-pressure undersaturated solution remains gas-free");
        }
        // Cooling a uniform, pressure-fed expanding liquid draws in cold oil.
        // Supply flow work must use the reservoir density at its own temperature.
        cfg.dissolved_gas_species=DissolvedGasSpecies::AIR;
        cfg.fluid_property_model=FluidPropertyModel::SINGLE_PHASE;
        cfg.density_model=DensityModel::LINEAR_TEMPERATURE;
        cfg.dissolved_gas_initial=cfg.dissolved_gas_henry_coeff=0;
        cfg.gas_mass_transfer_rate=cfg.omega=0;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        cfg.temperature_initial=333.15;
        cfg.journal_heat_transfer=10000; cfg.journal_wall_temperature=280;
        InletConfig cold_feed; cold_feed.type=InletConfig::Type::GROOVE;
        cold_feed.size=360; cold_feed.p_supply=1e6; cold_feed.feeds_fresh_oil=true;
        cold_feed.t_supply=293.15; cfg.inlets={cold_feed};
        reset(1,0,cold_feed.p_supply-cfg.p_cav);
        const double mass_before=f["oil_mass"](0,0), energy_before=f["internal_energy"](0,0);
        result=solver.solve(f,cfg,dt);
        require(result.converged,"Thermally expanding liquid reservoir solve");
        const double temperature=f["temperature"](0,0), mass_after=f["oil_mass"](0,0);
        const double rho_feed=FluidProperties::liquid_solution_density(0,cold_feed.p_supply,cold_feed.t_supply,cfg);
        const double supply_enthalpy=cfg.rho_cp/cfg.rho*(cold_feed.t_supply-cfg.property_reference_temperature)
            +cold_feed.p_supply/rho_feed;
        require(mass_after>mass_before,"Cooling must draw liquid from the pressure reservoir");
        near(f["internal_energy"](0,0)-energy_before+dt*cfg.journal_heat_transfer*(temperature-280)
                -(mass_after-mass_before)*supply_enthalpy,0,1e-6,
             "Reservoir flow work uses supply temperature, not receiving cell temperature");
        if (mesh.rank==0) std::cout << "PASS: coupled heat, compression, dissolution and R410A energy\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; failed=1; }
    int any=0; MPI_Allreduce(&failed,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    PetscFinalize(); return any;
}
