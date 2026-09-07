#include <petscsys.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "coupled_cavitation.hpp"

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main(int argc, char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int failed=0;
    try {
        SimulationConfig cfg;
        cfg.n_theta_global=24; cfg.n_z_global=8; cfg.e=0; cfg.omega=0;
        cfg.rho=850; cfg.mu=.03; cfg.mu_gas=1.8e-5; cfg.p_cav=100;
        cfg.temperature_reference=313.15; cfg.dissolved_gas_species=DissolvedGasSpecies::AIR;
        cfg.fluid_property_model=FluidPropertyModel::TWO_PHASE;
        cfg.dissolved_gas_initial=.001;
        cfg.bc_z_south_val=2e5; cfg.bc_z_north_val=1e5;
        cfg.bc_z_south_dissolved_gas=.003; cfg.bc_z_north_dissolved_gas=.005;
        Mesh mesh(cfg); Communicator comm(mesh); Fields fields;
        for (auto name : {"pressure","h","theta","dissolved_gas","free_gas_mass","temperature","alpha_gas"}) fields.add(name,mesh);
        CoupledCavitation::Solver solver(mesh,comm);
        const auto reset=[&] {
            fields["h"].fill(cfg.c); fields["theta"].fill(1);
            fields["pressure"].fill((cfg.bc_z_south_val+cfg.bc_z_north_val)/2);
            fields["temperature"].fill(cfg.temperature_reference);
            fields["dissolved_gas"].fill(cfg.dissolved_gas_initial); fields["free_gas_mass"].fill(0);
            solver.initialize(fields,cfg);
        };
        constexpr double dt=.003;
        for (bool reverse : {false,true}) {
            if (reverse) std::swap(cfg.bc_z_south_val,cfg.bc_z_north_val);
            reset();
            const auto stats=solver.solve(fields,cfg,dt);
            require(stats.converged,"Reverse-capable axial solve");
            if (mesh.rank==0) std::cout << "axial residual " << stats.residual << " oil " << stats.oil_balance << " species " << stats.species_balance << '\n';
            const double courant=dt*cfg.c*cfg.c*1e5/(12*cfg.mu*cfg.L*mesh.get_d_z());
            double expected=reverse ? cfg.bc_z_north_dissolved_gas : cfg.bc_z_south_dissolved_gas;
            for (int index=0;index<mesh.n_z_local;++index) {
                expected=(cfg.dissolved_gas_initial+courant*expected)/(1+courant);
                const int j=reverse ? mesh.n_z_local-1-index : index;
                for(int i=0;i<mesh.n_theta_local;++i)
                    require(std::abs(fields["dissolved_gas"](i,j)-expected)<1e-10,"Independent implicit upwind concentration profile");
            }
            require(std::abs(stats.oil_balance)<1e-12 && std::abs(stats.species_balance)<1e-12,"Boundary species conservation");
        }
        // An internal supply below end pressure drains local composition. Its
        // nominal feed concentration is relevant only after its flow reverses.
        cfg.bc_z_south_val=cfg.bc_z_north_val=2e5;
        cfg.bc_z_south_dissolved_gas=cfg.bc_z_north_dissolved_gas=.001;
        for (auto geometry : {InletConfig::Type::GROOVE,InletConfig::Type::CIRCULAR}) {
            InletConfig inlet; inlet.type=geometry; inlet.theta=90; inlet.z=cfg.L/2;
            inlet.size=geometry==InletConfig::Type::GROOVE ? 30 : cfg.L*.22;
            inlet.dissolved_gas=.009;
            for (bool feeding : {false,true}) {
                inlet.p_supply=feeding ? 3e5 : 1e5; cfg.inlets={inlet}; reset();
                const auto stats=solver.solve(fields,cfg,dt);
                if (mesh.rank==0) std::cout << "inlet geometry=" << static_cast<int>(geometry) << " feeding=" << feeding << " residual=" << stats.residual << " iterations=" << stats.iterations << '\n';
                require(stats.converged,"Groove/hole reservoir solve");
                require(feeding ? stats.liquid_inlet_source>0 : stats.liquid_inlet_source<0,"Internal reservoir flow sign");
                require(std::abs(stats.oil_balance)<1e-12 && std::abs(stats.species_balance)<1e-12,"Internal reservoir conservation");
                if (!feeding) for(int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j)
                    require(std::abs(fields["dissolved_gas"](i,j)-.001)<1e-10,"Outflow must retain local composition");
            }
        }
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; failed=1; }
    int global_failed=0; MPI_Allreduce(&failed,&global_failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    PetscFinalize(); return global_failed;
}
