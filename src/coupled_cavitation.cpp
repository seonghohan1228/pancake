#include "coupled_cavitation.hpp"

#include <petscsnes.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "fluid_properties.hpp"
#include "r410a.hpp"

namespace CoupledCavitation {
namespace {
constexpr std::array<const char*,5> names{{"pressure", "oil_mass", "dissolved_mass", "free_gas_mass", "temperature"}};
void check(PetscErrorCode code) {
    if (code) throw std::runtime_error("Coupled cavitation PETSc error " + std::to_string(code));
}
double harmonic(double a, double b) { return 2*a*b/(a+b); }
double fb(double a, double b) {
    const double length=std::hypot(a,b), sum=a+b;
    return sum>0 ? 2*a*b/(sum+length) : sum-length;
}
struct State {
    double p, h, t, rho_l, mu_l, rho_g, theta, mobility, energy=0, conductivity=0, couette_h=0;
    std::array<double,3> u{}, volume{};
    std::array<double,3> mass;
    double total() const { return mass[0]+mass[1]+mass[2]; }
    double concentration() const { return mass[1]/(mass[0]+mass[1]); }
};
}
struct Face {
    std::array<double,3> mass{}, pressure{}, moving{}, diffusion{};
    double volume=0, energy=0, dissipation=0, conduction=0, diffusion_work=0;
    double operator[](int k) const { return mass[k]; }
};
struct Solver::Impl {
    const Mesh& mesh;
    Communicator& comm;
    const bool thermal;
    const int dof;
    Fields trial;
    Fields* fields = nullptr;
    const SimulationConfig* cfg = nullptr;
    SNES snes = nullptr;
    Vec x = nullptr, residual = nullptr;
    Mat jacobian = nullptr;
    MatFDColoring coloring = nullptr;
    double dt = 0, p_scale = 1e5, m_scale = 1, c_scale = .001, g_scale = 1;
    bool equilibrium = false;
    bool steady_initializing = false;
    Stats stats;
    std::vector<State> states;
    std::vector<double> oil_residual, species_residual, liquid_residual, energy_residual, heat, reaction_field;

    Impl(const Mesh& m, Communicator& c, bool energy) : mesh(m), comm(c), thermal(energy), dof(energy ? 5 : 4) {
        states.resize((mesh.n_theta_local+2)*mesh.n_z_local);
        for (auto name : names) trial.add(name,mesh);
        const PetscInt local = dof*mesh.n_theta_local*mesh.n_z_local;
        const PetscInt global = dof*mesh.n_theta_global*mesh.n_z_global;
        check(VecCreateMPI(MPI_COMM_WORLD,local,global,&x));
        check(VecDuplicate(x,&residual));
        check(MatCreateAIJ(MPI_COMM_WORLD,local,local,global,global,5*dof,nullptr,5*dof,nullptr,&jacobian));
        // Every residual depends on all active states in the five-cell stencil.
        for (int i=0; i<mesh.n_theta_local; ++i) for (int j=0; j<mesh.n_z_local; ++j) {
            const int gi=mesh.offset_theta+i;
            std::vector<PetscInt> columns;
            for (auto offset : {std::array<int,2>{0,0}, {1,0}, {-1,0}, {0,1}, {0,-1}}) {
                const int jj=j+offset[1];
                if (jj<0 || jj>=mesh.n_z_global) continue;
                const int ii=(gi+offset[0]+mesh.n_theta_global)%mesh.n_theta_global;
                for (int k=0;k<dof;++k) columns.push_back(dof*(ii*mesh.n_z_global+jj)+k);
            }
            std::vector<PetscScalar> values(columns.size(),1);
            for (int k=0;k<dof;++k) {
                const PetscInt row=dof*(gi*mesh.n_z_global+j)+k;
                check(MatSetValues(jacobian,1,&row,static_cast<PetscInt>(columns.size()),columns.data(),values.data(),INSERT_VALUES));
            }
        }
        check(MatAssemblyBegin(jacobian,MAT_FINAL_ASSEMBLY));
        check(MatAssemblyEnd(jacobian,MAT_FINAL_ASSEMBLY));
        MatColoring mc; ISColoring colors;
        check(MatColoringCreate(jacobian,&mc));
        check(MatColoringSetType(mc,MATCOLORINGSL));
        check(MatColoringSetDistance(mc,2));
        check(MatColoringApply(mc,&colors));
        check(MatColoringDestroy(&mc));
        check(MatFDColoringCreate(jacobian,colors,&coloring));
#if PETSC_VERSION_GE(3,24,0)
        check(MatFDColoringSetFunction(coloring,fd_function,this));
#else
        // Older PETSc headers expose this four-argument callback untyped.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif
        check(MatFDColoringSetFunction(coloring,reinterpret_cast<PetscErrorCode(*)(void)>(fd_function),this));
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif
        check(MatFDColoringSetUp(jacobian,colors,coloring));
        check(ISColoringDestroy(&colors));
        check(SNESCreate(MPI_COMM_WORLD,&snes));
        check(SNESSetOptionsPrefix(snes,"coupled_"));
        check(SNESSetFunction(snes,residual,function,this));
        check(SNESSetJacobian(snes,jacobian,jacobian,SNESComputeJacobianDefaultColor,coloring));
        check(SNESSetType(snes,SNESNEWTONLS));
        KSP ksp; PC pc;
        check(SNESGetKSP(snes,&ksp));
        check(KSPSetType(ksp,KSPGMRES));
        check(KSPGMRESSetRestart(ksp,100));
        check(KSPSetTolerances(ksp,1e-9,PETSC_DEFAULT,PETSC_DEFAULT,1000));
        check(KSPGetPC(ksp,&pc));
        check(PCSetType(pc,PCASM));
        check(PCASMSetOverlap(pc,4));
        check(SNESSetFromOptions(snes));
    }
    ~Impl() {
        SNESDestroy(&snes); MatFDColoringDestroy(&coloring);
        MatDestroy(&jacobian); VecDestroy(&x); VecDestroy(&residual);
    }
    double scale(int k) const {
        return k==0 ? p_scale : k==1 ? m_scale : k==2 ? m_scale*c_scale : k==3 ? g_scale : 10;
    }
    double gas_density(double p, double t) const {
        if (cfg->dissolved_gas_species != DissolvedGasSpecies::R410A || p <= 0)
            return p/(cfg->gas_specific_constant_for_species()*t);
        // A bounded continuation lets Newton visit gas-free high-pressure cells.
        // Accepted cells containing free refrigerant undergo the strict EOS check.
        const double limit=std::min(R410A::maximum_pressure,R410A::dew_pressure(t));
        if (p<=limit) return R410A::gas_density(p,t);
        const double density=R410A::gas_density(limit,t);
        return density+(p-limit)/R410A::pressure_derivative(density,t);
    }
    void unpack(Vec input) {
        const PetscScalar* values;
        check(VecGetArrayRead(input,&values));
        int n=0;
        for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
            for (int k=0;k<dof;++k) trial[names[k]](i,j)=values[n++]*scale(k)+(k==0 ? cfg->p_cav : 0);
        check(VecRestoreArrayRead(input,&values));
        comm.update_ghosts(trial);
    }
    State state(int i, int j) const {
        State s{};
        s.p=trial["pressure"](i,j); s.h=(*fields)["h"](i,j);
        s.t=thermal ? trial["temperature"](i,j) : (*fields)["temperature"](i,j);
        for (int k=0;k<3;++k) s.mass[k]=trial[names[k+1]](i,j);
        properties(s);
        if (cfg->cross_film_thermal_model==CrossFilmModel::FINITE_VOLUME) {
            s.mobility=s.total()/s.h*(*fields)["cross_film_f2"](i,j);
            s.couette_h=(*fields)["cross_film_couette_h"](i,j);
        }
        return s;
    }
    const State& at(int i, int j) const { return states[(i+1)*mesh.n_z_local+j]; }
    void properties(State& s) const {
        if (!std::isfinite(s.t) || s.t<200 || s.t>500) throw std::domain_error("Coupled model requires 200 <= T <= 500 K");
        const double ml=s.mass[0]+s.mass[1];
        if (!(s.h>0) || !(ml>0) || !(s.total()>0)) throw std::domain_error("Nonpositive liquid inventory or gap");
        s.rho_l=FluidProperties::liquid_solution_density(s.concentration(),s.p,s.t,*cfg);
        s.mu_l=FluidProperties::liquid_solution_viscosity(s.concentration(),s.p,s.t,*cfg);
        s.theta=ml/(s.rho_l*s.h);
        // During Newton iteration the dilute EOS is continued through zero.
        // The volume complementarity equation enforces nonnegative partial p.
        s.rho_g=gas_density(s.p-cfg->oil_vapor_pressure(s.t),s.t);
        const double mu_mix=s.theta*s.mu_l+std::max(0.0,1-s.theta)*cfg->mu_gas;
        s.mobility=s.total()*s.h*s.h/(12*mu_mix);
        s.couette_h=s.h;
        const double liquid_volume = cfg->density_model == DensityModel::LINEAR_TEMPERATURE ? 1/s.rho_l : 1/cfg->rho;
        s.volume={liquid_volume, cfg->density_model==DensityModel::MIXTURE && cfg->dissolved_gas_liquid_density>0 ?
            1/cfg->dissolved_gas_liquid_density : liquid_volume, s.rho_g>0 ? 1/s.rho_g : 0};
        if (thermal || cfg->cross_film_thermal_model==CrossFilmModel::FINITE_VOLUME) {
            s.energy=0;
            for(int k=0;k<3;++k)s.u[k]=s.mass[k]!=0?partial_internal_energy(k,s.t,s.p,*cfg):0;
            for (int k=0;k<3;++k) s.energy+=s.mass[k]*s.u[k];
            s.conductivity=s.h*(s.theta*cfg->thermal_conductivity+
                std::max(0.0,1-s.theta)*cfg->gas_thermal_conductivity);
        }
    }
    State reservoir(const State& cell, bool south) const {
        State s=cell;
        s.p=south ? cfg->bc_z_south_val : cfg->bc_z_north_val;
        const bool incoming=s.p>cell.p;
        const auto thermal=south ? cfg->bc_z_south_thermal : cfg->bc_z_north_thermal;
        if (incoming && thermal==ThermalInflowMode::CONSTANT) s.t=cfg->axial_inflow_temperature(south);
        const double c=incoming ? cfg->axial_inflow_dissolved_gas(south) : cell.concentration();
        const double q=incoming ? (south ? cfg->bc_z_south_free_gas_mass_fraction : cfg->bc_z_north_free_gas_mass_fraction) : cell.mass[2]/cell.total();
        s.rho_l=FluidProperties::liquid_solution_density(c,s.p,s.t,*cfg);
        if (s.p<cfg->oil_vapor_pressure(s.t)) throw std::domain_error("Reservoir pressure is below oil vapor pressure at its temperature");
        const double rg=gas_density(s.p-cfg->oil_vapor_pressure(s.t),s.t);
        // Only prescribed incoming gas is a fixed physical reservoir state.
        // Outgoing trial gas must permit the same continuation as interior cells.
        if (incoming && q>0 && cfg->dissolved_gas_species==DissolvedGasSpecies::R410A)
            R410A::gas_density(s.p-cfg->oil_vapor_pressure(s.t),s.t);
        if (q>0 && rg<=0) throw std::domain_error("Reservoir with free gas requires p > p_cav");
        const double mt=s.h/((1-q)/s.rho_l+(q>0 ? q/rg : 0));
        s.mass={mt*(1-q)*(1-c),mt*(1-q)*c,mt*q};
        properties(s);
        return s;
    }
    const InletConfig* inlet(int i, int j) const { return mesh.inlet_at(*cfg,i,j); }
    std::array<double,3> internal_energies(double t, double p, double) const {
        return CoupledCavitation::internal_energies(t,p,*cfg);
    }
    Face flux(const State& a, const State& b, double distance, double width, double speed) const {
        const double fp=-harmonic(a.mobility,b.mobility)*(b.p-a.p)/distance*width;
        const State& donor=fp>=0 ? a : b;
        const State& moving=speed>=0 ? a : b;
        Face f;
        const double qp=fp*donor.h/donor.total();
        const double h_face=.5*(a.couette_h+b.couette_h);
        f.volume=qp+speed*h_face*width;
        f.dissipation=-qp*(b.p-a.p);
        for (int k=0;k<3;++k) {
            const double pressure_mass=fp*donor.mass[k]/donor.total();
            const double moving_mass=speed*h_face*moving.mass[k]/moving.h*width;
            f.pressure[k]=pressure_mass; f.moving[k]=moving_mass;
            f.mass[k]=pressure_mass+moving_mass;
            f.energy+=pressure_mass*donor.u[k]+moving_mass*moving.u[k];
        }
        const double fd=-cfg->dissolved_gas_diffusivity*harmonic(a.mass[0]+a.mass[1],b.mass[0]+b.mass[1])*
            (b.concentration()-a.concentration())/distance*width;
        f.mass[1]+=fd; f.mass[0]-=fd;
        f.diffusion={-fd,fd,0};
        f.diffusion_work=fd*.5*(a.p*(a.volume[1]-a.volume[0])+b.p*(b.volume[1]-b.volume[0]));
        // Fickian diffusion transports partial enthalpy, including displacement work.
        f.energy+=fd*.5*(a.u[1]-a.u[0]+a.p*(a.volume[1]-a.volume[0])+
            b.u[1]-b.u[0]+b.p*(b.volume[1]-b.volume[0]));
        if (thermal && a.conductivity+b.conductivity>0)
            f.conduction=-harmonic(a.conductivity,b.conductivity)*(b.t-a.t)/distance*width;
        return f;
    }
    static PetscErrorCode function(SNES solver, Vec input, Vec output, void* context) {
        auto& self=*static_cast<Impl*>(context);
        // All ranks finish ghost exchange before any local property-domain check.
        self.unpack(input);
        int failed=0;
        try { self.evaluate(output); } catch (const std::domain_error&) { failed=1; }
        int any_failed=0;
        MPI_Allreduce(&failed,&any_failed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        if (any_failed) return SNESSetFunctionDomainError(solver);
        return 0;
    }
    static PetscErrorCode fd_function(void* solver, Vec input, Vec output, void* context) {
        return function(static_cast<SNES>(solver),input,output,context);
    }
    void evaluate(Vec output) {
        std::vector<double> result(dof*mesh.n_theta_local*mesh.n_z_local);
        oil_residual.resize(result.size()/dof); species_residual.resize(result.size()/dof); liquid_residual.resize(result.size()/dof);
        energy_residual.resize(result.size()/dof); heat.resize(result.size()/dof); reaction_field.resize(result.size()/dof);
        stats={};
        const double dx=mesh.R*mesh.get_d_theta(), dz=mesh.get_d_z(), area=dx*dz;
        const double speed=cfg->omega*mesh.R/2;
        for (int i=-1;i<=mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
            states[(i+1)*mesh.n_z_local+j]=state(i,j);
        int n=0;
        for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j) {
            const State& s=at(i,j);
            const auto east=flux(s,at(i+1,j),dx,dz,speed);
            const auto west=flux(at(i-1,j),s,dx,dz,speed);
            Face north{},south{};
            if (j+1<mesh.n_z_local) north=flux(s,at(i,j+1),dz,dx,0);
            else if (cfg->bc_z_north_type!=BCType::NEUMANN) north=flux(s,reservoir(s,false),dz/2,dx,0);
            if (j>0) south=flux(at(i,j-1),s,dz,dx,0);
            else if (cfg->bc_z_south_type!=BCType::NEUMANN) south=flux(reservoir(s,true),s,dz/2,dx,0);
            if (j==0) south.conduction=0;
            if (j+1==mesh.n_z_local) north.conduction=0;
            std::array<double,3> balance;
            for (int k=0;k<3;++k) balance[k]=
                (cfg->solution_mode==SolutionMode::TRANSIENT ? s.mass[k]-(*fields)[names[k+1]].old(i,j) : 0)+
                dt*(east[k]-west[k]+north[k]-south[k])/area;
            // The inverse chart needs only the saturation boundary at the actual
            // composition, even in undersaturated high-pressure liquid cells.
            const bool bubble_constraint=equilibrium && cfg->solubility_model==SolubilityModel::CURVE_TABLE;
            const bool exchange=cfg->fluid_property_model==FluidPropertyModel::TWO_PHASE &&
                (cfg->gas_mass_transfer_rate>0 || cfg->effective_gas_absorption_rate()>0);
            const double ceq=!bubble_constraint && (equilibrium || exchange) ?
                std::min(1-1e-12,FluidProperties::equilibrium_dissolved_gas(s.p-cfg->oil_vapor_pressure(s.t),s.t,*cfg)) : 0;
            const double md_eq=std::min(s.mass[1]+s.mass[2],ceq/(1-ceq)*s.mass[0]);
            const double rate=cfg->fluid_property_model!=FluidPropertyModel::TWO_PHASE ? 0 :
                (s.mass[1]>=md_eq ? cfg->gas_mass_transfer_rate : cfg->effective_gas_absorption_rate());
            const double reaction=rate*(1-ceq)*(s.mass[1]-md_eq);
            double source_o=0, source_d=0, source_g=0;
            const auto* supply=inlet(i,j);
            if (supply) {
                const double q=balance[0]+balance[1]+balance[2];
                const double c=supply->dissolved_gas>=0 ? supply->dissolved_gas : cfg->dissolved_gas_initial;
                const double quality=supply->free_gas_mass_fraction;
                source_o=q*(q>=0 ? (1-quality)*(1-c) : s.mass[0]/s.total());
                source_d=q*(q>=0 ? (1-quality)*c : s.mass[1]/s.total());
                source_g=q*(q>=0 ? quality : s.mass[2]/s.total());
            }
            {
                const std::array<double,3> source{{source_o,source_d,source_g}};
                double source_work=0;
                for (int k=0;k<3;++k) {
                    (*fields)["coupled_source_"+std::to_string(k)](i,j)=source[k]*area/dt;
                    double volume=s.volume[k];
                    if (supply && source_o+source_d+source_g>=0 && supply->feeds_fresh_oil && k<2 && cfg->density_model==DensityModel::LINEAR_TEMPERATURE)
                        volume=1/FluidProperties::liquid_solution_density(0,s.p,supply->t_supply,*cfg);
                    source_work+=source[k]*area/dt*s.p*volume;
                }
                (*fields)["coupled_source_work"](i,j)=source_work;
            }
            const double a=s.rho_g/(g_scale/cfg->c), b=1-s.theta, g=s.mass[2]/(g_scale*s.h/cfg->c);
            const double radical=a*a+b*b+2*g;
            if (radical<0) throw std::domain_error("Negative gas-volume complementarity radical");
            const int cell=n/dof;
            oil_residual[cell]=(balance[0]-source_o)*area/dt;
            species_residual[cell]=(balance[1]+balance[2]-source_d-source_g)*area/dt;
            const double accepted_reaction=equilibrium ? -(balance[1]-source_d)/dt : reaction;
            reaction_field[cell]=accepted_reaction/s.h;
            const double mu_mix=s.theta*s.mu_l+std::max(0.0,1-s.theta)*cfg->mu_gas;
            const double couette_heat=cfg->cross_film_thermal_model==CrossFilmModel::FINITE_VOLUME ?
                4*speed*speed/(*fields)["cross_film_f0"](i,j) : mu_mix*4*speed*speed/s.h;
            heat[cell]=couette_heat+(.5*(east.dissipation+west.dissipation)+
                (j==0 ? 1 : .5)*south.dissipation+(j+1==mesh.n_z_local ? 1 : .5)*north.dissipation)/area;
            double energy_balance=0;
            if (thermal) {
                double source_energy=0;
                if (supply) {
                    State feed=s;
                    if (source_o+source_d+source_g>=0 && supply->feeds_fresh_oil) feed.t=supply->t_supply;
                    const double pg=feed.p-cfg->oil_vapor_pressure(feed.t);
                    if (pg<0) throw std::domain_error("Feed pressure is below oil vapor pressure at its temperature");
                    if (source_o+source_d+source_g>=0 && supply->free_gas_mass_fraction>0 &&
                        cfg->dissolved_gas_species==DissolvedGasSpecies::R410A) R410A::gas_density(pg,feed.t);
                    const double rg=gas_density(pg,feed.t);
                    const auto u=internal_energies(feed.t,feed.p,rg);
                    const double gas_enthalpy=u[2]+(rg>0 ? feed.p/rg : 0);
                    auto volume=s.volume;
                    if (cfg->density_model==DensityModel::LINEAR_TEMPERATURE)
                        volume[0]=volume[1]=1/FluidProperties::liquid_solution_density(0,feed.p,feed.t,*cfg);
                    source_energy=source_o*(u[0]+feed.p*volume[0])+
                        source_d*(u[1]+feed.p*volume[1])+source_g*gas_enthalpy;
                }
                const double storage=cfg->solution_mode==SolutionMode::TRANSIENT ?
                    s.energy-(*fields)["internal_energy"].old(i,j) : 0;
                const double dh=cfg->solution_mode==SolutionMode::TRANSIENT ? s.h-(*fields)["h"].old(i,j) : 0;
                const double work=s.p*(dh+dt*(east.volume-west.volume+north.volume-south.volume)/area);
                const double wall=cfg->journal_heat_transfer*(s.t-cfg->journal_wall_temperature)+
                    cfg->bearing_heat_transfer*(s.t-cfg->bearing_wall_temperature);
                energy_balance=storage+dt*(east.energy-west.energy+north.energy-south.energy+
                    east.conduction-west.conduction+north.conduction-south.conduction)/area+
                    work-dt*(heat[cell]-wall)-source_energy;
                stats.energy_balance+=energy_balance*area/dt;
                stats.internal_energy+=s.energy*area;
            }
            energy_residual[cell]=energy_balance*area/dt;
            liquid_residual[cell]=(balance[0]+balance[1]+dt*accepted_reaction-source_o-source_d)*area/dt;
            stats.oil_balance+=oil_residual[cell]; stats.species_balance+=species_residual[cell];
            stats.liquid_balance+=liquid_residual[cell];
            stats.liquid_mass+=(s.mass[0]+s.mass[1])*area;
            stats.gas_mass+=(s.mass[1]+s.mass[2])*area;
            stats.liquid_inlet_source+=(source_o+source_d)*area/dt;
            stats.gas_inlet_source+=(source_d+source_g)*area/dt;
            if (j==0) { stats.liquid_boundary_flux+=south[0]+south[1]; stats.gas_boundary_flux+=south[1]+south[2]; }
            if (j+1==mesh.n_z_local) { stats.liquid_boundary_flux-=north[0]+north[1]; stats.gas_boundary_flux-=north[1]+north[2]; }
            result[n++]=supply ? (s.p-supply->p_supply)/p_scale : balance[0]/m_scale;
            const double sum=a+b, root=std::sqrt(radical);
            result[n++]=sum>0 ? 2*(a*b-g)/(sum+root) : sum-root;
            const double total_species=(balance[1]+balance[2]-source_d-source_g)/(m_scale*c_scale);
            if (equilibrium) {
                result[n++]=total_species;
                const double undersaturation=bubble_constraint ?
                    (s.p-cfg->oil_vapor_pressure(s.t)-cfg->saturation_pressure(s.concentration(),s.t))/p_scale :
                    (ceq-s.concentration())/c_scale;
                result[n++]=fb(undersaturation,s.mass[2]/g_scale);
            } else {
                result[n++]=(balance[1]+dt*reaction-source_d)/(m_scale*c_scale);
                result[n++]=total_species;
            }
            if (thermal) result[n++]=energy_balance/(m_scale*(cfg->liquid_specific_heat>0 ? cfg->liquid_specific_heat : cfg->rho_cp/cfg->rho)*10);
        }
        PetscScalar* values;
        check(VecGetArray(output,&values));
        std::copy(result.begin(),result.end(),values);
        check(VecRestoreArray(output,&values));
    }
};

double partial_internal_energy(int phase,double t,double p,const SimulationConfig& cfg) {
    const double cp=cfg.liquid_specific_heat>0 ? cfg.liquid_specific_heat : cfg.rho_cp/cfg.rho;
    const double delta=t-cfg.property_reference_temperature;
    if(phase==0)return cp*delta;
    if(cfg.fluid_property_model==FluidPropertyModel::CONSTANT)return phase==1?cp*delta:0;
    if(phase==1)return cfg.dissolved_heat_capacity*delta-cfg.desorption_internal_energy;
    if(cfg.dissolved_gas_species!=DissolvedGasSpecies::R410A)return cfg.gas_heat_capacity_cv*delta;
    const double pg=p-cfg.oil_vapor_pressure(t);
    const double limit=std::min(R410A::maximum_pressure,R410A::dew_pressure(t));
    double rg=0;
    if(pg<=limit)rg=R410A::gas_density(std::max(0.0,pg),t);
    else {
        const double rlim=R410A::gas_density(limit,t);
        rg=rlim+(pg-limit)/R410A::pressure_derivative(rlim,t);
    }
    thread_local double last_p=-1,last_t=-1,reference_energy=0;
    if(last_p!=cfg.property_reference_pressure||last_t!=cfg.property_reference_temperature) {
        const double reference_rho=R410A::gas_density(cfg.property_reference_pressure,cfg.property_reference_temperature);
        reference_energy=R410A::internal_energy(reference_rho,cfg.property_reference_temperature);
        last_p=cfg.property_reference_pressure;last_t=cfg.property_reference_temperature;
    }
    return R410A::internal_energy(rg,t)-reference_energy;
}
std::array<double,3> internal_energies(double t,double p,const SimulationConfig& cfg) {
    return {partial_internal_energy(0,t,p,cfg),partial_internal_energy(1,t,p,cfg),partial_internal_energy(2,t,p,cfg)};
}

Solver::Solver(const Mesh& mesh, Communicator& comm, bool thermal) : impl(std::make_unique<Impl>(mesh,comm,thermal)) {}
Solver::~Solver() = default;
void Solver::initialize(Fields& fields, const SimulationConfig& cfg) {
    const auto& m=impl->mesh;
    if (impl->thermal != (cfg.temperature_model==TemperatureModel::ENERGY_EQUATION && cfg.cross_film_thermal_model!=CrossFilmModel::FINITE_VOLUME))
        throw std::invalid_argument("Coupled solver thermal mode differs from configuration");
    impl->cfg=&cfg; impl->fields=&fields;
    for (auto name : {"internal_energy","energy_residual","heat_generation","gas_mass_transfer"})
        if (!fields.has(name)) fields.add(name,m);
    for (auto name : {"mixture_mass_flux_theta","mixture_volume_flux_theta"})
        if (!fields.has(name)) fields.add(name,m,2,GridLocation::FACE_THETA);
    for (auto name : {"mixture_mass_flux_z","mixture_volume_flux_z"})
        if (!fields.has(name)) fields.add(name,m,2,GridLocation::FACE_Z);
    {
        if (cfg.cross_film_thermal_model==CrossFilmModel::FINITE_VOLUME && cfg.solution_mode!=SolutionMode::STEADY_STATE) throw std::invalid_argument("Resolved coupled energy requires steady state");
        for (const char* axis : {"theta","z"}) {
            const auto location=std::string(axis)=="theta" ? GridLocation::FACE_THETA : GridLocation::FACE_Z;
            for (int k=0;k<3;++k) for (const char* part : {"pressure","moving","diffusion"})
                if(!fields.has("coupled_"+std::string(part)+"_"+std::to_string(k)+"_"+axis)) fields.add("coupled_"+std::string(part)+"_"+std::to_string(k)+"_"+axis,m,2,location);
            if(!fields.has("coupled_diffusion_work_"+std::string(axis)))fields.add("coupled_diffusion_work_"+std::string(axis),m,2,location);
        }
        for (int k=0;k<3;++k) if(!fields.has("coupled_source_"+std::to_string(k)))fields.add("coupled_source_"+std::to_string(k),m);
        if(!fields.has("coupled_source_work"))fields.add("coupled_source_work",m);
    }
    for (auto name : {"oil_mass","dissolved_mass"}) if (!fields.has(name)) fields.add(name,m);
    for (int i=0;i<m.n_theta_local;++i) for (int j=0;j<m.n_z_local;++j) {
        const double c=fields["dissolved_gas"](i,j);
        const double rho=FluidProperties::liquid_solution_density(c,fields["pressure"](i,j),fields["temperature"](i,j),cfg);
        const double ml=rho*fields["h"](i,j)*fields["theta"](i,j);
        fields["oil_mass"](i,j)=ml*(1-c); fields["dissolved_mass"](i,j)=ml*c;
        const double mu=FluidProperties::liquid_solution_viscosity(c,fields["pressure"](i,j),fields["temperature"](i,j),cfg);
        if (fields.has("rho_liquid_solution")) fields["rho_liquid_solution"](i,j)=rho;
        if (fields.has("mu_liquid_solution")) fields["mu_liquid_solution"](i,j)=mu;
        if (fields.has("rho")) fields["rho"](i,j)=(ml+fields["free_gas_mass"](i,j))/fields["h"](i,j);
        if (fields.has("mu")) fields["mu"](i,j)=fields["theta"](i,j)*mu+std::max(0.0,1-fields["theta"](i,j))*cfg.mu_gas;
    }
    if (impl->thermal || cfg.cross_film_thermal_model==CrossFilmModel::FINITE_VOLUME) {
        for (int i=0;i<m.n_theta_local;++i) for (int j=0;j<m.n_z_local;++j) {
            for (int k=0;k<impl->dof;++k) impl->trial[names[k]](i,j)=fields[names[k]](i,j);
            impl->trial["temperature"](i,j)=fields["temperature"](i,j);
            fields["internal_energy"](i,j)=impl->state(i,j).energy;
        }
    }
    impl->comm.update_ghosts(fields);
    store_old_time(fields);
}
void Solver::store_old_time(Fields& fields) {
    for (int k=1;k<4;++k) fields[names[k]].store_old_time();
    for (auto name : {"h","temperature","internal_energy"}) fields[name].store_old_time();
}
Stats Solver::solve(Fields& fields, const SimulationConfig& cfg, double dt, bool equilibrium) {
    if (!(dt>0)) throw std::invalid_argument("Coupled cavitation requires positive time step");
    if (cfg.solution_mode==SolutionMode::STEADY_STATE && !impl->steady_initializing) {
        struct Reset { bool& flag; ~Reset() { flag=false; } } reset{impl->steady_initializing};
        impl->steady_initializing=true;
        Stats result=solve(fields,cfg,dt,equilibrium);
        if (result.converged) return result;
        // Conservative pseudo-time continuation supplies a wet/dry initial
        // guess. Only the final zero-storage residual establishes steady state.
        auto transient=cfg;
        transient.solution_mode=SolutionMode::TRANSIENT;
        double pseudo_dt=std::min(dt,1e-4);
        int iterations=result.iterations;
        store_old_time(fields);
        for (int warm=0;warm<300;++warm) {
            const Stats advance=solve(fields,transient,pseudo_dt,equilibrium);
            iterations+=advance.iterations;
            if (!advance.converged) {
                pseudo_dt*=.5;
                if (pseudo_dt<1e-12) break;
                continue;
            }
            store_old_time(fields);
            if (warm%5==4) {
                result=solve(fields,cfg,dt,equilibrium);
                iterations+=result.iterations;
                if (result.converged) { result.iterations=iterations; return result; }
            }
            pseudo_dt=std::min(pseudo_dt*1.3,.05);
        }
        result.converged=false; result.iterations=iterations;
        return result;
    }
    auto& s=*impl;
    s.fields=&fields; s.cfg=&cfg; s.dt=dt; s.equilibrium=equilibrium;
    s.p_scale=std::max({1e5,cfg.bc_z_south_val,cfg.bc_z_north_val,cfg.mu*std::abs(cfg.omega)*cfg.R*cfg.R/(cfg.c*cfg.c)});
    s.m_scale=cfg.rho*cfg.c;
    s.c_scale=std::max(cfg.dissolved_gas_initial,.001);
    s.g_scale=s.p_scale/(cfg.gas_specific_constant_for_species()*cfg.temperature_reference)*cfg.c;
    s.comm.update_ghosts(fields["h"]); s.comm.update_ghosts(fields["temperature"]);
    PetscScalar* values;
    check(VecGetArray(s.x,&values));
    int n=0;
    for (int i=0;i<s.mesh.n_theta_local;++i) for (int j=0;j<s.mesh.n_z_local;++j) {
        const auto* supply=s.inlet(i,j);
        for (int k=0;k<s.dof;++k) {
            const double value=k==0 && supply ? supply->p_supply : fields[names[k]](i,j);
            values[n++]=(value-(k==0 ? cfg.p_cav : 0))/s.scale(k);
        }
    }
    check(VecRestoreArray(s.x,&values));
    // Early steady attempts only test whether continuation reached Newton's
    // basin. Long, slowly improving attempts delay the next useful warm step.
    PetscReal atol=1e-10, rtol=0, stol=1e-14;
    PetscInt max_iterations=cfg.solution_mode==SolutionMode::STEADY_STATE ? 30 : 100, max_functions=10000;
    // Read only the stopping controls here; repeating SNESSetFromOptions would
    // reinstall monitors and reset solver/preconditioner options each timestep.
    check(PetscOptionsGetReal(nullptr,"coupled_","-snes_atol",&atol,nullptr));
    check(PetscOptionsGetReal(nullptr,"coupled_","-snes_rtol",&rtol,nullptr));
    check(PetscOptionsGetReal(nullptr,"coupled_","-snes_stol",&stol,nullptr));
    check(PetscOptionsGetInt(nullptr,"coupled_","-snes_max_it",&max_iterations,nullptr));
    check(PetscOptionsGetInt(nullptr,"coupled_","-snes_max_funcs",&max_functions,nullptr));
    check(SNESSetTolerances(s.snes,atol,rtol,stol,max_iterations,max_functions));
    check(SNESSolve(s.snes,nullptr,s.x));
    // Query failure before reevaluating: a rejected EOS trial may be outside its domain.
    SNESConvergedReason solve_reason;
    check(SNESGetConvergedReason(s.snes,&solve_reason));
    if (solve_reason<=0) {
        Stats failure; PetscInt iterations;
        check(SNESGetIterationNumber(s.snes,&iterations));
        failure.iterations=iterations; failure.residual=INFINITY; return failure;
    }
    s.unpack(s.x); s.evaluate(s.residual);
    PetscReal norm; PetscInt iterations; SNESConvergedReason reason;
    check(VecNorm(s.residual,NORM_INFINITY,&norm));
    check(SNESGetIterationNumber(s.snes,&iterations));
    check(SNESGetConvergedReason(s.snes,&reason));
    Stats result; result.converged=reason>0 && norm<1e-10; result.iterations=iterations; result.residual=norm;
    if (!result.converged) return result;
    std::array<double,9> local{{s.stats.oil_balance,s.stats.species_balance,s.stats.liquid_balance,
        s.stats.liquid_mass,s.stats.gas_mass,s.stats.liquid_boundary_flux,s.stats.gas_boundary_flux,
        s.stats.liquid_inlet_source,s.stats.gas_inlet_source}}, global{};
    MPI_Allreduce(local.data(),global.data(),9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    result.oil_balance=global[0]; result.species_balance=global[1]; result.liquid_balance=global[2];
    result.liquid_mass=global[3]; result.gas_mass=global[4];
    result.liquid_boundary_flux=global[5]; result.gas_boundary_flux=global[6];
    result.liquid_inlet_source=global[7]; result.gas_inlet_source=global[8];
    double energy_local[2]={s.stats.energy_balance,s.stats.internal_energy}, energy_global[2];
    MPI_Allreduce(energy_local,energy_global,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    result.energy_balance=energy_global[0]; result.internal_energy=energy_global[1];
    int invalid=0;
    std::array<double,5> local_bounds{{1e300,1e300,1e300,1e300,1e300}}, global_bounds{};
    double local_volume_error=0;
    for (int i=0;i<s.mesh.n_theta_local;++i) for (int j=0;j<s.mesh.n_z_local;++j) {
        const State& cell=s.at(i,j);
        const double pv=cfg.oil_vapor_pressure(cell.t);
        local_bounds[0]=std::min(local_bounds[0],cell.p-pv);
        local_bounds[1]=std::min(local_bounds[1],cell.mass[0]/s.m_scale);
        local_bounds[2]=std::min(local_bounds[2],cell.mass[1]/(s.m_scale*s.c_scale));
        local_bounds[3]=std::min(local_bounds[3],cell.mass[2]/s.g_scale);
        local_bounds[4]=std::min(local_bounds[4],1-cell.theta);
        if (cell.p<pv-1e-6 || cell.mass[0]<0 || cell.mass[1]<-1e-12*s.m_scale || cell.mass[2]<-1e-9*s.g_scale || cell.theta>1+1e-8) invalid=1;
        if (cfg.dissolved_gas_species==DissolvedGasSpecies::R410A && cell.mass[2]>1e-10*s.g_scale) {
            try { R410A::gas_density(cell.p-pv,cell.t); } catch (const std::domain_error&) { invalid=1; }
        }
        if (cell.mass[2]>1e-10*s.g_scale && cell.rho_g>0)
            local_volume_error=std::max(local_volume_error,std::abs(cell.theta+cell.mass[2]/(cell.rho_g*cell.h)-1));
        if (cfg.solubility_model==SolubilityModel::TABLE) {
            const double pg=cell.p-pv;
            const auto& table=cfg.solubility_table_2d;
            if (!table.empty()) {
                if (pg<table.p.front() || pg>table.p.back() || cell.t<table.t.front() || cell.t>table.t.back()) invalid=1;
            } else if (!cfg.solubility_table.empty() &&
                (pg<cfg.solubility_table.front().x || pg>cfg.solubility_table.back().x)) invalid=1;
        }
    }
    int any_invalid=0;
    MPI_Allreduce(&invalid,&any_invalid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(local_bounds.data(),global_bounds.data(),5,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&local_volume_error,&result.volume_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (any_invalid) {
        std::ostringstream message;
        message << "Coupled state violates positivity, the solubility table bounds, or the stable free-gas EOS domain; state was not accepted. Minima: p-p_v="
            << global_bounds[0] << " Pa, scaled oil=" << global_bounds[1] << ", dissolved=" << global_bounds[2]
            << ", free=" << global_bounds[3] << ", void=" << global_bounds[4];
        throw std::domain_error(message.str());
    }
    for (int i=0;i<s.mesh.n_theta_local;++i) for (int j=0;j<s.mesh.n_z_local;++j) {
        const State& cell=s.at(i,j);
        for (int k=0;k<s.dof;++k) fields[names[k]](i,j)=s.trial[names[k]](i,j);
        fields["theta"](i,j)=cell.theta;
        fields["dissolved_gas"](i,j)=cell.concentration();
        if (fields.has("alpha_gas")) fields["alpha_gas"](i,j)=cell.rho_g>0 ? cell.mass[2]/(cell.rho_g*cell.h) : 0;
        if (fields.has("rho")) fields["rho"](i,j)=cell.total()/cell.h;
        if (fields.has("rho_gas")) fields["rho_gas"](i,j)=cell.rho_g;
        if (fields.has("rho_liquid_solution")) fields["rho_liquid_solution"](i,j)=cell.rho_l;
        if (fields.has("mu_liquid_solution")) fields["mu_liquid_solution"](i,j)=cell.mu_l;
        if (fields.has("mu")) fields["mu"](i,j)=cell.theta*cell.mu_l+std::max(0.0,1-cell.theta)*cfg.mu_gas;
        fields["internal_energy"](i,j)=cell.energy;
        const int index=i*s.mesh.n_z_local+j;
        fields["energy_residual"](i,j)=s.energy_residual[index];
        fields["heat_generation"](i,j)=s.heat[index];
        fields["gas_mass_transfer"](i,j)=s.reaction_field[index];
        if (fields.has("mass_residual")) fields["mass_residual"](i,j)=s.liquid_residual[i*s.mesh.n_z_local+j];
    }
    const double dx=s.mesh.R*s.mesh.get_d_theta(), dz=s.mesh.get_d_z();
    auto export_flux=[&](const Face& face,const char* axis,int i,int j) {
        for (int k=0;k<3;++k) {
            fields["coupled_pressure_"+std::to_string(k)+"_"+axis](i,j)=face.pressure[k];
            fields["coupled_moving_"+std::to_string(k)+"_"+axis](i,j)=face.moving[k];
            fields["coupled_diffusion_"+std::to_string(k)+"_"+axis](i,j)=face.diffusion[k];
        }
        fields["coupled_diffusion_work_"+std::string(axis)](i,j)=face.diffusion_work;
    };
    for (int i=0;i<s.mesh.n_theta_local;++i) {
        for (int j=0;j<s.mesh.n_z_local;++j) {
            const auto& a=s.at(i-1,j);
            const auto& b=s.at(i,j);
            const auto f=s.flux(a,b,dx,dz,cfg.omega*cfg.R/2);
            export_flux(f,"theta",i,j);
            fields["mixture_mass_flux_theta"](i,j)=f[0]+f[1]+f[2];
            fields["mixture_volume_flux_theta"](i,j)=f.volume;
            if (fields.has("velocity_theta")) fields["velocity_theta"](i,j)=f.volume/(dz*.5*(a.h+b.h));
        }
        for (int j=0;j<=s.mesh.n_z_local;++j) {
            const auto& cell=s.at(i,std::min(j,s.mesh.n_z_local-1));
            Face f; double h_face=cell.h;
            if (j==0 && cfg.bc_z_south_type!=BCType::NEUMANN)
                f=s.flux(s.reservoir(cell,true),cell,dz/2,dx,0);
            else if (j==s.mesh.n_z_local && cfg.bc_z_north_type!=BCType::NEUMANN)
                f=s.flux(cell,s.reservoir(cell,false),dz/2,dx,0);
            else if (j>0 && j<s.mesh.n_z_local) {
                const auto& a=s.at(i,j-1);
                f=s.flux(a,cell,dz,dx,0); h_face=.5*(a.h+cell.h);
            }
            export_flux(f,"z",i,j);
            fields["mixture_mass_flux_z"](i,j)=f[0]+f[1]+f[2];
            fields["mixture_volume_flux_z"](i,j)=f.volume;
            if (fields.has("velocity_z")) fields["velocity_z"](i,j)=f.volume/(dx*h_face);
        }
    }
    s.comm.update_ghosts(fields);
    return result;
}
}
