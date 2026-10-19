#include <petscsys.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "coupled_cavitation.hpp"
#include "r410a.hpp"

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main(int argc, char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int failed=0;
    try {
        SimulationConfig cfg;
        cfg.n_theta_global=16; cfg.n_z_global=4; cfg.e=0; cfg.omega=0;
        cfg.rho=850; cfg.mu=.03; cfg.mu_gas=1.8e-5; cfg.p_cav=100;
        cfg.temperature_reference=313.15; cfg.dissolved_gas_species=DissolvedGasSpecies::AIR;
        cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE;
        cfg.dissolved_gas_initial=.002; cfg.dissolved_gas_henry_coeff=1e-8;
        cfg.gas_mass_transfer_rate=20; cfg.gas_absorption_rate=0;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        Mesh mesh(cfg); Communicator comm(mesh); Fields fields;
        for (auto name : {"pressure","h","theta","dissolved_gas","free_gas_mass","temperature","alpha_gas"}) fields.add(name,mesh);
        fields["h"].fill(cfg.c); fields["theta"].fill(.9);
        fields["pressure"].fill(1e5+cfg.p_cav); fields["temperature"].fill(cfg.temperature_reference);
        fields["dissolved_gas"].fill(cfg.dissolved_gas_initial);
        const double mg=1e5/(287.05*cfg.temperature_reference)*cfg.c*.1;
        fields["free_gas_mass"].fill(mg);
        CoupledCavitation::Solver solver(mesh,comm); solver.initialize(fields,cfg);
        const double mo=fields["oil_mass"](0,0), total=fields["dissolved_mass"](0,0)+mg;
        // Runtime iteration limits must survive per-step default setup. A failed
        // attempt must preserve the last accepted pressure and inventories.
        for (const auto* limit : {"0", "1"}) {
            PetscOptionsSetValue(nullptr,"-coupled_snes_max_it",limit);
            const auto stopped=solver.solve(fields,cfg,.001);
            require(!stopped.converged && stopped.iterations==std::stoi(limit),
                    "Coupled runtime iteration limit and failure count");
            require(fields["pressure"](0,0)==1e5+cfg.p_cav && fields["oil_mass"](0,0)==mo,
                    "Rejected solve must not overwrite accepted state");
        }
        PetscOptionsClearValue(nullptr,"-coupled_snes_max_it");
        for (int step=0;step<10;++step) {
            const auto stats=solver.solve(fields,cfg,.001);
            if (mesh.rank==0) std::cout << "step " << step << " converged=" << stats.converged << " residual=" << stats.residual << " iterations=" << stats.iterations << '\n';
            require(stats.converged,"Closed-cell coupled solve failed");
            if (step==0) {
                const double once=fields["pressure"](0,0);
                const auto repeat=solver.solve(fields,cfg,.001);
                require(repeat.converged && std::abs(fields["pressure"](0,0)-once)<1e-6,
                        "Nonlinear repeats must retain the same previous physical time");
            }
            for(int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j) {
                require(std::abs(fields["oil_mass"](i,j)/mo-1)<1e-9,"Oil conservation");
                require(std::abs((fields["dissolved_mass"](i,j)+fields["free_gas_mass"](i,j))/total-1)<1e-9,"Species conservation");
                require(std::abs(fields["theta"](i,j)+fields["alpha_gas"](i,j)-1)<1e-8,"Gap volume closure");
            }
            solver.store_old_time(fields);
        }
        require(fields["pressure"](0,0)>100100,"Release must raise pressure in a closed fixed gap");
        // With exchange disabled, gas compression must obey Boyle's law while
        // preserving oil and gas separately. The new gap is ten percent smaller.
        cfg.gas_mass_transfer_rate=0;
        const double p0=fields["pressure"](0,0)-cfg.p_cav;
        const double ml=fields["oil_mass"](0,0)+fields["dissolved_mass"](0,0);
        const double v0=cfg.c-ml/cfg.rho;
        const double new_h=ml/cfg.rho+v0/2;
        fields["h"].fill(new_h);
        const auto compression=solver.solve(fields,cfg,.001);
        require(compression.converged,"Closed compression solve");
        require(std::abs((fields["pressure"](0,0)-cfg.p_cav)/(2*p0)-1)<1e-8,"Boyle pressure-volume relation");
        require(std::abs(R410A::gas_density(2e6,313.15)/75.71095-1)<1e-6,"R410A density reference");
        bool rejected=false;
        try { R410A::gas_density(2e6,293.15); } catch (const std::domain_error&) { rejected=true; }
        require(rejected,"R410A condensed state rejection");
        // Flooded axial pressure-driven flow must be the exact linear profile.
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::DIRICHLET;
        cfg.bc_z_south_val=2e5; cfg.bc_z_north_val=1e5;
        cfg.dissolved_gas_initial=0; cfg.dissolved_gas_henry_coeff=0;
        fields["h"].fill(cfg.c); fields["theta"].fill(1);
        fields["dissolved_gas"].fill(0); fields["free_gas_mass"].fill(0);
        fields["pressure"].fill(1.5e5); solver.initialize(fields,cfg);
        auto axial=solver.solve(fields,cfg,.001);
        require(axial.converged,"Axial full-film solve");
        for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
            require(std::abs(fields["pressure"](i,j)-(2e5-1e5*(j+.5)/mesh.n_z_local))<.01,"Linear reservoir profile");
        // A rotating eccentric film must admit a vapor cavity with no gas mass.
        cfg.omega=1000; cfg.bc_z_south_val=cfg.bc_z_north_val=1e5;
        cfg.fluid_property_model=FluidPropertyModel::CONSTANT;
        fields["pressure"].fill(1e5); fields["theta"].fill(1);
        for (int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j)
            fields["h"](i,j)=cfg.c*(1+.6*std::cos((mesh.offset_theta+i+.5)*mesh.get_d_theta()));
        solver.initialize(fields,cfg);
        for (int step=0;step<20;++step) {
            const auto bearing=solver.solve(fields,cfg,.001);
            if (mesh.rank==0) std::cout << "bearing " << step << " converged=" << bearing.converged << " residual=" << bearing.residual << '\n';
            require(bearing.converged,"Rotating vaporous bearing solve");
            for (int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j) {
                require(fields["pressure"](i,j)>=cfg.p_cav-.01,"Vapor pressure floor");
                require(fields["theta"](i,j)<=1+1e-8,"Liquid volume bound");
                require(std::abs((fields["pressure"](i,j)-cfg.p_cav)*(1-fields["theta"](i,j)))<.01,"JFO complementarity");
            }
            solver.store_old_time(fields);
        }
        double local_min=1, global_min=1;
        for (int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j)
            local_min=std::min(local_min,fields["theta"](i,j));
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        require(global_min<.99,"Vapor benchmark must actually cavitate");
        cfg.solution_mode=SolutionMode::STEADY_STATE;
        const auto steady_bearing=solver.solve(fields,cfg,.001);
        require(steady_bearing.converged,"Steady conservative vapor bearing");
        // Dissolved transport cannot create a free phase when its equation has no source or incoming mass.
        const double original_p_cav=cfg.p_cav;
        cfg.fluid_property_model=FluidPropertyModel::SINGLE_PHASE;
        cfg.dissolved_gas_initial=.3; cfg.p_cav=0;
        fields["dissolved_gas"].fill(.3); fields["free_gas_mass"].fill(0);
        for (double eccentricity : {.6,.62,.58,.6}) {
            for (int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j)
                fields["h"](i,j)=cfg.c*(1+eccentricity*std::cos((mesh.offset_theta+i+.5)*mesh.get_d_theta()));
            solver.initialize(fields,cfg);
            const auto gas_off=solver.solve(fields,cfg,.001);
            require(gas_off.converged,"Gas-off dissolved mixture pressure solve");
            require(std::abs(gas_off.oil_balance)<1e-10 && std::abs(gas_off.species_balance)<1e-10,
                    "Gas-off mixture preserves oil and total species");
            double bounds[3]={0,0,INFINITY};
            for (int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j) {
                bounds[0]=std::min(bounds[0],fields["free_gas_mass"](i,j));
                bounds[1]=std::min(bounds[1],-fields["free_gas_mass"](i,j));
                bounds[2]=std::min(bounds[2],fields["pressure"](i,j));
            }
            MPI_Allreduce(MPI_IN_PLACE,bounds,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
            if(mesh.rank==0)std::cout<<std::setprecision(17)<<"gas_off_invariant eccentricity="<<eccentricity
                <<" free_min="<<bounds[0]<<" free_max="<<-bounds[1]<<" p_min="<<bounds[2]
                <<" volume_error="<<gas_off.volume_error<<" oil_balance="<<gas_off.oil_balance
                <<" species_balance="<<gas_off.species_balance<<" residual="<<gas_off.residual<<'\n';
            require(bounds[0]==0 && bounds[1]==0,"Source-free free-gas inventory remains exactly zero");
            require(bounds[2]>=-1e-6,"Gas-off mixture preserves the strict absolute pressure floor");
        }
        cfg.p_cav=original_p_cav;
        cfg.solution_mode=SolutionMode::TRANSIENT;

        // Independent equilibrium root from Henry's law, gas EOS and conserved
        // oil/total species. This is a physical-volume flash, not a fixed-p flash.
        cfg.omega=0; cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE;
        cfg.bc_z_south_type=cfg.bc_z_north_type=BCType::NEUMANN;
        cfg.dissolved_gas_initial=.002; cfg.dissolved_gas_henry_coeff=1e-8;
        fields["h"].fill(cfg.c); fields["theta"].fill(.9);
        fields["dissolved_gas"].fill(.002); fields["free_gas_mass"].fill(mg);
        fields["pressure"].fill(1e5+cfg.p_cav); solver.initialize(fields,cfg);
        const double oil=fields["oil_mass"](0,0), gas=fields["dissolved_mass"](0,0)+mg;
        const double rt=287.05*cfg.temperature_reference, henry=cfg.dissolved_gas_henry_coeff;
        const double qa=cfg.c*henry, qb=cfg.c-oil/cfg.rho+rt*henry*(gas+oil), qc=rt*gas;
        const double exact_p=2*qc/(qb+std::sqrt(qb*qb-4*qa*qc));
        const auto flash=solver.solve(fields,cfg,.001,true);
        require(flash.converged,"Equilibrium coupled flash");
        require(std::abs((fields["pressure"](0,0)-cfg.p_cav)/exact_p-1)<1e-8,"Independent equilibrium pressure root");
        require(std::abs(fields["dissolved_gas"](0,0)-henry*exact_p)<1e-10,"Equilibrium solubility");

        // The real R410A EOS participates inside the nonlinear volume solve.
        // Constitutive oil constants below are synthetic verification inputs.
        cfg.dissolved_gas_species=DissolvedGasSpecies::R410A;
        cfg.dissolved_gas_initial=.10; cfg.dissolved_gas_henry_coeff=1e-7;
        cfg.gas_mass_transfer_rate=20;
        for (double temperature : {293.15,313.15,333.15,353.15}) {
            cfg.temperature_reference=temperature;
            const double gas_mass=R410A::gas_density(6e5,temperature)*cfg.c*.1;
            fields["temperature"].fill(temperature); fields["theta"].fill(.9);
            fields["dissolved_gas"].fill(.1); fields["free_gas_mass"].fill(gas_mass);
            fields["pressure"].fill(6e5+cfg.p_cav); solver.initialize(fields,cfg);
            const double total_gas=fields["dissolved_mass"](0,0)+gas_mass;
            const auto refrigerant=solver.solve(fields,cfg,.001);
            require(refrigerant.converged,"R410A coupled finite-rate solve");
            const double pg=fields["pressure"](0,0)-cfg.p_cav;
            const double alpha=fields["free_gas_mass"](0,0)/(R410A::gas_density(pg,temperature)*cfg.c);
            require(std::abs(fields["theta"](0,0)+alpha-1)<1e-8,"R410A real-fluid volume closure");
            require(std::abs((fields["dissolved_mass"](0,0)+fields["free_gas_mass"](0,0))/total_gas-1)<1e-9,"R410A inventory conservation");
        }
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; failed=1; }
    int global_failed=0; MPI_Allreduce(&failed,&global_failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    PetscFinalize(); return global_failed;
}
