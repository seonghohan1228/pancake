#include "resolved_cross_film.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <petscksp.h>

#include "communicator.hpp"
#include "coupled_cavitation.hpp"
#include "fbns.hpp"
#include "fluid_properties.hpp"
#include "full_film.hpp"
#include "reynolds_coefficients.hpp"

namespace ResolvedCrossFilm {
namespace {
std::string layer_name(int k) { return "cross_film_layer_" + std::to_string(k); }
std::string viscosity_name(int k) { return "cross_film_inv_mu_" + std::to_string(k); }
double value(const Fields& f, const char* name, double fallback, int i, int j) {
    return f.has(name) ? f[name](i,j) : fallback;
}
void require_valid(bool invalid, const char* message) {
    int local = invalid ? 1 : 0, global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global) throw std::runtime_error(message);
}
const InletConfig* inlet_at(const Mesh& mesh, const SimulationConfig& cfg, int i, int j) {
    const auto* inlet = mesh.inlet_at(cfg,i,j);
    return inlet && inlet->feeds_fresh_oil ? inlet : nullptr;
}
double liquid(const Fields& f, const SimulationConfig& cfg, int i, int j) {
    if (cfg.cavitation_model==CavitationModel::COUPLED) return 1;
    return cfg.cavitated_film_model == CavitatedFilmModel::STRIATED
        && FullFilm::gate(f,cfg,i,j) <= 0.0 ? std::clamp(f["theta"](i,j),0.0,1.0) : 1.0;
}
double conductivity(const Fields& f,const SimulationConfig& cfg,int i,int j) {
    if (cfg.cavitation_model!=CavitationModel::COUPLED)
        return value(f,"k_liquid_solution",cfg.thermal_conductivity,i,j)*liquid(f,cfg,i,j);
    const double fill=f["theta"](i,j);
    return fill*cfg.thermal_conductivity+(1-fill)*cfg.gas_thermal_conductivity;
}

double film_conductance(const Fields& f,const SimulationConfig& cfg,int i,int j) {
    return 2.0*cfg.cross_film_layers*conductivity(f,cfg,i,j)/f["h"](i,j);
}
double outer_temperature(const SimulationConfig& cfg) {
    return cfg.bush_outer_temperature > 0 ? cfg.bush_outer_temperature : cfg.thermal_ambient_temperature;
}
double shaft_temperature(const SimulationConfig& cfg) {
    return cfg.shaft_ambient_temperature > 0 ? cfg.shaft_ambient_temperature : cfg.thermal_ambient_temperature;
}
double bush_conductance(const SimulationConfig& cfg) {
    const double ro=cfg.bush_outer_radius > cfg.R ? cfg.bush_outer_radius : 1.1*cfg.R;
    double resistance=cfg.R*std::log(ro/cfg.R)/cfg.bush_thermal_conductivity;
    if (cfg.bush_outer_boundary == BushOuterBoundary::CONVECTION) {
        if (cfg.bush_outer_heat_transfer <= 0) return 0;
        resistance += cfg.R/(ro*cfg.bush_outer_heat_transfer);
    }
    return 1.0/resistance;
}

// Churchill & Chu (1975), vertical-plate analogue with height 2*R_outer.
// The annular geometry is an approximation; properties are at the film temperature.
double end_heat_transfer(double node_temperature,double dz,const SimulationConfig& cfg,bool south) {
    if (!cfg.bush_end_natural_convection) return cfg.bush_end_heat_transfer;
    const double ambient=outer_temperature(cfg),length=2*cfg.bush_outer_radius;
    double surface=node_temperature,h=0;
    for(int iteration=0;iteration<20;++iteration) {
        const double film=0.5*(surface+ambient),p=cfg.property_reference_pressure;
        const double rho=FluidProperties::liquid_solution_density(0,p,film,cfg);
        const double mu=FluidProperties::liquid_solution_viscosity(cfg.axial_inflow_dissolved_gas(south),p,film,cfg);
        const double nu=mu/rho,alpha=cfg.thermal_conductivity/(rho*cfg.liquid_specific_heat);
        const double pr=nu/alpha,beta=-cfg.density_temperature_slope/rho;
        const double ra=9.80665*std::max(beta,0.0)*std::abs(surface-ambient)*length*length*length/(nu*alpha);
        h=cfg.thermal_conductivity/length*std::pow(.825+.387*std::pow(ra,1.0/6.0)/
          std::pow(1+std::pow(.492/pr,9.0/16.0),8.0/27.0),2);
        const double next=ambient+(node_temperature-ambient)/(1+h*dz/(2*cfg.bush_thermal_conductivity));
        if(std::abs(next-surface)<1e-8)break;
        surface=next;
    }
    return h;
}

// Exact layer integrals for the piecewise-constant viscosity: Couette velocity
// U(1-I0/F0) and Poiseuille velocity dp/ds*(I1-F1*I0/F0).
struct Shape { std::vector<double> couette, pressure; };
Shape velocity_shape(const Fields& f,const SimulationConfig& cfg,int i,int j) {
    const int n=cfg.cross_film_layers;
    Shape shape{std::vector<double>(n),std::vector<double>(n)};
    const double h=f["h"](i,j), dy=h/n;
    const double f0=f["cross_film_f0"](i,j), f1=f["cross_film_f1"](i,j), f2=f["cross_film_f2"](i,j);
    double i0=0,i1=0;
    for(int k=0;k<n;++k) {
        const double inv=f[viscosity_name(k)](i,j), y=k*dy;
        shape.couette[k]=(dy*(1.0-i0/f0)-inv*dy*dy/(2*f0))/(f1/f0);
        shape.pressure[k]=-(dy*(i1-f1*i0/f0)+inv*((y-f1/f0)*dy*dy/2+dy*dy*dy/6))/f2;
        i0+=inv*dy; i1+=inv*dy*(y+0.5*dy);
    }
    return shape;
}
struct Fluxes {
    int nt,nz,n;
    std::vector<double> east,north;
    Fluxes(const Fields& f,const Mesh& mesh,const SimulationConfig& cfg)
        :nt(mesh.n_theta_local),nz(mesh.n_z_local),n(cfg.cross_film_layers),
         east(static_cast<size_t>(nt+1)*nz*n),north(static_cast<size_t>(nt)*(nz+1)*n) {
        if (cfg.cavitation_model==CavitationModel::COUPLED) return;
        std::vector<Shape> shapes;
        shapes.reserve(static_cast<size_t>(nt+2)*nz);
        for(int i=-1;i<=nt;++i) for(int j=0;j<nz;++j) shapes.push_back(velocity_shape(f,cfg,i,j));
        auto shape=[&](int i,int j)->const Shape& {return shapes[(i+1)*nz+j];};
        for(int i=-1;i<nt;++i) for(int j=0;j<nz;++j) {
            const double c=ReynoldsCoeff::couette_face_flux(f,mesh,cfg,i,j);
            const double mc=c*f["theta"](c>=0?i:i+1,j);
            const double mp=FBNS::mass_flux_theta(f,mesh,cfg,i,j)-mc;
            for(int k=0;k<n;++k) e(i,j,k)=mc*0.5*(shape(i,j).couette[k]+shape(i+1,j).couette[k])
                +mp*0.5*(shape(i,j).pressure[k]+shape(i+1,j).pressure[k]);
        }
        for(int i=0;i<nt;++i) for(int j=-1;j<nz;++j) {
            const int js=std::clamp(j,0,nz-1),jn=std::clamp(j+1,0,nz-1);
            const double mass=FBNS::mass_flux_z(f,mesh,cfg,i,j);
            for(int k=0;k<n;++k) z(i,j,k)=mass*0.5*(shape(i,js).pressure[k]+shape(i,jn).pressure[k]);
        }
    }
    double& e(int i,int j,int k) {return east[((i+1)*nz+j)*n+k];}
    double e(int i,int j,int k) const {return east[((i+1)*nz+j)*n+k];}
    double& z(int i,int j,int k) {return north[(i*(nz+1)+j+1)*n+k];}
    double z(int i,int j,int k) const {return north[(i*(nz+1)+j+1)*n+k];}
    double divergence(int i,int j,int k,bool feed) const {
        return e(i,j,k)-e(i-1,j,k)+(j+1<nz||!feed?z(i,j,k):0)-(j>0||!feed?z(i,j-1,k):0);
    }
};

// The three conservative flux components can counterflow, so each retains its
// own upwind caloric state. Their layer sums equal the accepted hydraulic flux.
struct CoupledFluxes {
    const Fields& f; const Mesh& mesh; const SimulationConfig& cfg;
    int n,nz,nt;
    std::vector<Shape> shapes;
    CoupledFluxes(const Fields& fields,const Mesh& m,const SimulationConfig& c):f(fields),mesh(m),cfg(c),n(c.cross_film_layers),nz(m.n_z_local),nt(m.n_theta_local) {
        if(c.cavitation_model==CavitationModel::COUPLED)
            for(int i=-1;i<=nt;++i)for(int j=0;j<nz;++j)shapes.push_back(velocity_shape(f,cfg,i,j));
    }
    const Shape& shape(int i,int j) const {return shapes[(i+1)*nz+j];}
    double face(bool theta,int i,int j,int k,int phase,int part) const {
        const char* component=part==0?"pressure":part==1?"moving":"diffusion";
        const char* axis=theta?"theta":"z";
        const double mass=f["coupled_"+std::string(component)+"_"+std::to_string(phase)+"_"+axis](i,j);
        if(part==2)return mass/n;
        const auto& a=shape(theta?i-1:i,theta?j:std::max(j-1,0));
        const auto& b=shape(i,theta?j:std::min(j,nz-1));
        return mass*.5*(part==0?a.pressure[k]+b.pressure[k]:a.couette[k]+b.couette[k]);
    }
    double divergence(int i,int j,int k,int phase) const {
        double net=0;
        for(int part=0;part<3;++part) net+=face(true,i+1,j,k,phase,part)-face(true,i,j,k,phase,part)+
            face(false,i,j+1,k,phase,part)-face(false,i,j,k,phase,part);
        return net;
    }
    double work(int i,int j) const {
        return f["pressure"](i,j)*(f["mixture_volume_flux_theta"](i+1,j)-f["mixture_volume_flux_theta"](i,j)+
            f["mixture_volume_flux_z"](i,j+1)-f["mixture_volume_flux_z"](i,j))+
            f["coupled_diffusion_work_theta"](i+1,j)-f["coupled_diffusion_work_theta"](i,j)+
            f["coupled_diffusion_work_z"](i,j+1)-f["coupled_diffusion_work_z"](i,j)-f["coupled_source_work"](i,j);
    }
};
struct Caloric { double slope,offset; };
Caloric caloric(int phase,double t,double p,const SimulationConfig& cfg) {
    const double u=CoupledCavitation::partial_internal_energy(phase,t,p,cfg);
    double slope=0;
    if(phase==0||cfg.fluid_property_model==FluidPropertyModel::CONSTANT)
        slope=phase<2?(cfg.liquid_specific_heat>0?cfg.liquid_specific_heat:cfg.rho_cp/cfg.rho):0;
    else if(phase==1)slope=cfg.dissolved_heat_capacity;
    else if(cfg.dissolved_gas_species!=DissolvedGasSpecies::R410A)slope=cfg.gas_heat_capacity_cv;
    else slope=(CoupledCavitation::partial_internal_energy(phase,t+.01,p,cfg)-
        CoupledCavitation::partial_internal_energy(phase,t-.01,p,cfg))/.02;
    return {slope,u-slope*(t-cfg.temperature_reference)};
}

// Owning PETSc objects ensures exceptions do not leak a production-size matrix.
struct System {
    Mat a=nullptr; Vec rhs=nullptr,x=nullptr; KSP ksp=nullptr;
    ~System() {if(ksp)KSPDestroy(&ksp);if(x)VecDestroy(&x);if(rhs)VecDestroy(&rhs);if(a)MatDestroy(&a);}
};
}

void register_fields(Fields& f,const Mesh& mesh,const SimulationConfig& cfg) {
    for(int k=0;k<cfg.cross_film_layers;++k) {
        f.add(layer_name(k),mesh).fill(cfg.temperature_initial);
        f.add(viscosity_name(k),mesh).fill(1.0/cfg.mu);
    }
    f.add("cross_film_layer_initialized",mesh).fill(0);
    f.add("cross_film_layer_residual",mesh).fill(0);
    if(cfg.bush_thermal_model!=BushThermalModel::LOCAL) {
        for(int k=0;k<cfg.bush_radial_cells;++k)
            f.add("cross_film_solid_"+std::to_string(k),mesh).fill(outer_temperature(cfg));
        f.add("temperature_housing_outer",mesh).fill(outer_temperature(cfg));
    }
}

CrossFilm::WallState update(Fields& f,const Mesh& mesh,const SimulationConfig& cfg,double relaxation) {
    const int n=cfg.cross_film_layers;
    const double area=mesh.cell_volume();
    bool invalid=false;
    const bool damp=cfg.cavitation_model==CavitationModel::COUPLED &&
        cfg.solution_mode==SolutionMode::STEADY_STATE;
    require_valid(damp && !(std::isfinite(relaxation)&&relaxation>0&&relaxation<=1),
                  "Resolved viscosity relaxation must lie in (0, 1]");
    double reductions[5]={0,0,0,0,0};
    double tmin=std::numeric_limits<double>::infinity(),tmax=-tmin;
    for(int i=0;i<mesh.n_theta_local;++i) for(int j=0;j<mesh.n_z_local;++j) {
        const double h=f["h"](i,j),bar=f["temperature"](i,j);
        if(f["cross_film_layer_initialized"](i,j)==0) {
            for(int k=0;k<n;++k) f[layer_name(k)](i,j)=bar;
            f["cross_film_t_shaft"](i,j)=bar;f["cross_film_t_bush"](i,j)=bar;
        }
        double f0=0,f1=0,g2=0;
        const double mu=value(f,"mu",cfg.mu,i,j),p=f["pressure"](i,j);
        const double dissolved=value(f,"dissolved_gas",0,i,j);
        double ref=0;
        try {ref=FluidProperties::liquid_solution_viscosity(dissolved,p,bar,cfg);} catch(...) {invalid=true;}
        for(int k=0;k<n;++k) {
            const double t=f[layer_name(k)](i,j),lo=double(k)/n,hi=double(k+1)/n;
            double inv=0;
            try {
                const double mul=FluidProperties::liquid_solution_viscosity(dissolved,p,t,cfg);
                inv=cfg.cavitation_model==CavitationModel::COUPLED ?
                    1/(f["theta"](i,j)*mul+(1-f["theta"](i,j))*cfg.mu_gas) : ref/(mu*mul);
            } catch(...) {invalid=true;}
            invalid|=!(std::isfinite(t)&&t>0&&std::isfinite(inv)&&inv>0&&std::isfinite(h)&&h>0);
            if(damp && f["cross_film_layer_initialized"](i,j)!=0 && relaxation<1) {
                const double previous=f[viscosity_name(k)](i,j);
                invalid|=!(std::isfinite(previous)&&previous>0);
                inv=relaxation*inv+(1-relaxation)*previous;
            }
            f[viscosity_name(k)](i,j)=inv;
            f0+=h*(hi-lo)*inv;f1+=h*h*(hi*hi-lo*lo)*inv/2;
            g2+=h*h*h*(hi*hi*hi-lo*lo*lo)*inv/3;
            tmin=std::min(tmin,t);tmax=std::max(tmax,t);
        }
        const double f2=g2-f1*f1/f0;
        invalid|=!(std::isfinite(f2)&&f2>0);
        f["cross_film_f0"](i,j)=f0;f["cross_film_f1"](i,j)=f1;
        f["cross_film_f2"](i,j)=f2;f["cross_film_couette_h"](i,j)=2*f1/f0;
        const double ts=f["cross_film_t_shaft"](i,j),tb=f["cross_film_t_bush"](i,j);
        const double g=film_conductance(f,cfg,i,j);
        const double qs=g*(f[layer_name(0)](i,j)-ts),qb=g*(f[layer_name(n-1)](i,j)-tb);
        f["cross_film_wall_flux"](i,j)=qs+qb;
        if(f.has("temperature_journal"))f["temperature_journal"](i,j)=ts;
        if(f.has("temperature_bearing"))f["temperature_bearing"](i,j)=tb;
        if(f.has("heat_flux_journal"))f["heat_flux_journal"](i,j)=-qs;
        if(f.has("heat_flux_bearing"))f["heat_flux_bearing"](i,j)=-qb;
        tmin=std::min({tmin,ts,tb});tmax=std::max({tmax,ts,tb});
        reductions[0]+=ts*area;reductions[1]+=tb*area;reductions[2]+=qs*area;reductions[3]+=qb*area;
        reductions[4]=std::max(reductions[4],f["cross_film_layer_residual"](i,j));
    }
    require_valid(invalid,"Resolved cross-film temperature or viscosity is nonphysical");
    MPI_Allreduce(MPI_IN_PLACE,reductions,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,reductions+4,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    Communicator comm(mesh);
    for(int k=0;k<n;++k) {comm.update_ghosts(f[layer_name(k)]);comm.update_ghosts(f[viscosity_name(k)]);}
    for(const char* name:{"cross_film_f0","cross_film_f1","cross_film_f2","cross_film_couette_h"})comm.update_ghosts(f[name]);
    CrossFilm::WallState state;
    const double surface=2*M_PI*mesh.R*mesh.L;
    state.t_shaft=reductions[0]/surface;state.t_bush_mean=reductions[1]/surface;
    state.q_shaft_net=reductions[2];state.q_bush_net=reductions[3];state.profile_residual=reductions[4];
    state.profile_t_min=tmin;state.profile_t_max=tmax;
    state.bush_2d_solved=cfg.bush_thermal_model!=BushThermalModel::LOCAL;
    return state;
}

void solve(Fields& f,const Mesh& mesh,const SimulationConfig& cfg) {
    const int n=cfg.cross_film_layers,nz=mesh.n_z_local,nt=mesh.n_theta_local;
    const bool solid=cfg.bush_thermal_model!=BushThermalModel::LOCAL;
    const int ns=cfg.bush_thermal_model==BushThermalModel::CONDUCTION_3D?cfg.bush_axial_cells:1;
    const int nr=solid?cfg.bush_radial_cells:0,block=nz*n+(solid?ns*(nr+1):0);
    const PetscInt shaft=mesh.n_theta_global*block,total=shaft+1;
    const PetscInt local=nt*block+(mesh.rank==mesh.size-1?1:0);
    const double area=mesh.cell_volume(),cp=cfg.liquid_specific_heat,ref=cfg.temperature_reference;
    const double ds=mesh.R*mesh.get_d_theta(),dz=mesh.get_d_z(),u=cfg.omega*mesh.R;
    const double ro=cfg.bush_outer_radius>mesh.R?cfg.bush_outer_radius:1.1*mesh.R;
    const double dr=solid?(ro-mesh.R)/nr:0;
    const double solid_dz=mesh.L/ns;
    Fluxes flux(f,mesh,cfg);
    CoupledFluxes phase_flux(f,mesh,cfg);
    const bool coupled=cfg.cavitation_model==CavitationModel::COUPLED;
    System system;
    std::vector<PetscInt> diagonal(local,9),off_diagonal(local,3);
    if(solid)for(int i=0;i<nt;++i)for(int slab=0;slab<ns;++slab)
        diagonal[i*block+nz*n+slab*(nr+1)]=nz/ns+2;
    if(mesh.rank==mesh.size-1) {
        diagonal.back()=nt*nz+1;
        off_diagonal.back()=(mesh.n_theta_global-nt)*nz;
    }
    MatCreateAIJ(PETSC_COMM_WORLD,local,local,total,total,0,diagonal.data(),0,off_diagonal.data(),&system.a);
    MatSetOption(system.a,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
    VecCreateMPI(PETSC_COMM_WORLD,local,total,&system.rhs);VecDuplicate(system.rhs,&system.x);
    VecSet(system.rhs,0.0);VecSet(system.x,0.0);
    auto index=[&](int i,int j,int k)->PetscInt {
        const int gi=(mesh.offset_theta+i+mesh.n_theta_global)%mesh.n_theta_global;
        return gi*block+j*n+k;
    };
    auto bore=[&](int i,int slab)->PetscInt{return index(i,0,0)+nz*n+slab*(nr+1);};
    auto radial=[&](int i,int slab,int k)->PetscInt{return bore(i,slab)+1+k;};
    auto add=[&](PetscInt row,PetscInt col,double x){MatSetValue(system.a,row,col,x,ADD_VALUES);};
    double shaft_diag=0;
    std::vector<double> previous(static_cast<size_t>(nt)*nz*n);
    std::vector<double> phase_divergence(3*n),phase_before(3*n),phase_total(3);
    for(int i=0;i<nt;++i) for(int j=0;j<nz;++j) {
        const bool feed=inlet_at(mesh,cfg,i,j)!=nullptr;
        double net=0;for(int k=0;k<n;++k)net+=flux.divergence(i,j,k,feed);
        double transverse=0,heat_total=0;
        const double h=f["h"](i,j),dy=h/n,kth=conductivity(f,cfg,i,j);
        const double fill=liquid(f,cfg,i,j),gfilm=film_conductance(f,cfg,i,j)*area;
        const double ghalf=kth/dy*area;
        const double gp=(f["pressure"](i+1,j)-f["pressure"](i-1,j))/(2*ds);
        const double gz=nz>1?(f["pressure"](i,std::min(j+1,nz-1))-f["pressure"](i,std::max(j-1,0)))/((std::min(j+1,nz-1)-std::max(j-1,0))*dz):0;
        const double f0=f["cross_film_f0"](i,j),ybar=f["cross_film_f1"](i,j)/f0;
        double raw_heat=0;
        if(coupled)for(int k=0;k<n;++k) {
            const double y=(k+.5)*dy,tau=-u/f0+gp*(y-ybar),tauz=gz*(y-ybar);
            raw_heat+=f[viscosity_name(k)](i,j)*dy*(tau*tau+tauz*tauz+(gp*gp+gz*gz)*dy*dy/12);
        }
        const double accepted_heat=f["heat_generation"](i,j);
        if(coupled)for(int phase=0;phase<3;++phase) {
            double phase_net=0;
            for(int k=0;k<n;++k) {
                const double divergence=phase_flux.divergence(i,j,k,phase);
                phase_divergence[phase*n+k]=divergence;
                phase_net+=divergence;
            }
            phase_total[phase]=phase_net;
            double before=0;
            for(int k=0;k<n;++k) {
                phase_before[phase*n+k]=before;
                before+=phase_net/n-phase_divergence[phase*n+k];
            }
        }
        for(int k=0;k<n;++k) {
            const PetscInt row=index(i,j,k);
            const double old=f[layer_name(k)](i,j);
            previous[(i*nz+j)*n+k]=old;
            double diag=0,source=0;
            auto advect=[&](double mass,PetscInt neighbour,double inflow) {
                const double rate=cp*mass;
                if(rate>=0)diag+=rate;
                else if(neighbour>=0)add(row,neighbour,rate);
                else source-=rate*(inflow-ref);
            };
            if(coupled) {
                auto phase_advect=[&](int phase,double mass,PetscInt neighbour,int ii,int jj,double boundary_t=0,double boundary_p=0,int donor_layer=-1) {
                    if(mass==0)return;
                    const bool own=mass>=0;
                    const double t=own?old:boundary_t>0?boundary_t:f[layer_name(donor_layer>=0?donor_layer:k)](ii,jj);
                    const double p=own?f["pressure"](i,j):boundary_t>0?boundary_p:f["pressure"](ii,jj);
                    const auto u=caloric(phase,t,p,cfg);
                    if(own)diag+=mass*u.slope;
                    else if(neighbour>=0)add(row,neighbour,mass*u.slope);
                    else source-=mass*u.slope*(t-ref);
                    source-=mass*u.offset;
                };
                for(int phase=0;phase<3;++phase) {
                    const double phase_net=phase_total[phase];
                    for(int part=0;part<3;++part) {
                        phase_advect(phase,phase_flux.face(true,i+1,j,k,phase,part),index(i+1,j,k),i+1,j);
                        phase_advect(phase,-phase_flux.face(true,i,j,k,phase,part),index(i-1,j,k),i-1,j);
                        phase_advect(phase,phase_flux.face(false,i,j+1,k,phase,part),j+1<nz?index(i,j+1,k):cfg.bc_z_north_thermal==ThermalInflowMode::CONSTANT?-1:row,
                            i,std::min(j+1,nz-1),j+1==nz&&cfg.bc_z_north_thermal==ThermalInflowMode::CONSTANT?cfg.axial_inflow_temperature(false):0,cfg.bc_z_north_val);
                        phase_advect(phase,-phase_flux.face(false,i,j,k,phase,part),j>0?index(i,j-1,k):cfg.bc_z_south_thermal==ThermalInflowMode::CONSTANT?-1:row,
                            i,std::max(j-1,0),j==0&&cfg.bc_z_south_thermal==ThermalInflowMode::CONSTANT?cfg.axial_inflow_temperature(true):0,cfg.bc_z_south_val);
                    }
                    const double before=phase_before[phase*n+k];
                    const double after=before+phase_net/n-phase_divergence[phase*n+k];
                    if(k>0)phase_advect(phase,-before,index(i,j,k-1),i,j,0,0,k-1);
                    if(k+1<n)phase_advect(phase,after,index(i,j,k+1),i,j,0,0,k+1);
                    const double mass=f["coupled_source_"+std::to_string(phase)](i,j)/n;
                    const auto* inlet=mesh.inlet_at(cfg,i,j);
                    if(mass>0) {
                        const double t=inlet&&inlet->feeds_fresh_oil?inlet->t_supply:old;
                        source+=mass*CoupledCavitation::partial_internal_energy(phase,t,f["pressure"](i,j),cfg);
                    } else if(mass<0) {
                        const auto u=caloric(phase,old,f["pressure"](i,j),cfg);
                        diag-=mass*u.slope;source+=mass*u.offset;
                    }
                }
                source-=phase_flux.work(i,j)/n;
            } else {
            advect(flux.e(i,j,k),index(i+1,j,k),0);
            advect(-flux.e(i-1,j,k),index(i-1,j,k),0);
            if(j+1<nz)advect(flux.z(i,j,k),index(i,j+1,k),0);
            else if(!feed)advect(flux.z(i,j,k),cfg.bc_z_north_thermal==ThermalInflowMode::CONSTANT?-1:row,cfg.axial_inflow_temperature(false));
            if(j>0)advect(-flux.z(i,j-1,k),index(i,j-1,k),0);
            else if(!feed)advect(-flux.z(i,j-1,k),cfg.bc_z_south_thermal==ThermalInflowMode::CONSTANT?-1:row,cfg.axial_inflow_temperature(true));
            if(k>0)advect(-transverse,index(i,j,k-1),0);
            transverse+=net/n-flux.divergence(i,j,k,feed);
            if(k+1<n)advect(transverse,index(i,j,k+1),0);
            if(feed) {
                const auto* inlet=inlet_at(mesh,cfg,i,j);
                if(net>=0)source+=cp*net/n*(inlet->t_supply-ref);
                else diag-=cp*net/n;
            } else diag-=cp*net/n;
            }
            auto conduct=[&](PetscInt neighbour,double g) {diag+=g;add(row,neighbour,-g);};
            if(k>0)conduct(index(i,j,k-1),ghalf);
            if(k+1<n)conduct(index(i,j,k+1),ghalf);
            if(k==0) {
                conduct(shaft,gfilm);add(shaft,row,-gfilm);shaft_diag+=gfilm;
            }
            if(k==n-1) {
                if(solid) {
                    const auto wall=bore(i,j*ns/nz);
                    conduct(wall,gfilm);add(wall,row,-gfilm);add(wall,wall,gfilm);
                } else {
                    const double gb=bush_conductance(cfg)*area;
                    const double g=gb>0?gfilm*gb/(gfilm+gb):0;
                    diag+=g;source+=g*(outer_temperature(cfg)-ref);
                }
            }
            // Thin-film surface conduction, one conservative link per face.
            const double own=kth*h/n;
            auto inplane=[&](int ii,int jj,double length,double distance) {
                const double other=conductivity(f,cfg,ii,jj)*f["h"](ii,jj)/n;
                conduct(index(ii,jj,k),0.5*(own+other)*length/distance);
            };
            inplane(i+1,j,dz,ds);inplane(i-1,j,dz,ds);
            if(j+1<nz)inplane(i,j+1,ds,dz);
            if(j>0)inplane(i,j-1,ds,dz);
            // Integral of the full squared Couette/Poiseuille stress. Its cross
            // term cancels on summation, giving precisely U^2/F0 + F2|grad p|^2.
            const double y=(k+0.5)*dy;
            const double tau=-u*std::sqrt(fill)/f0+gp*(y-ybar),tauz=gz*(y-ybar);
            const double raw=f[viscosity_name(k)](i,j)*dy*(tau*tau+tauz*tauz+(gp*gp+gz*gz)*dy*dy/12);
            const double heat=coupled?(raw_heat>0?accepted_heat*raw/raw_heat:accepted_heat/n):raw;
            heat_total+=heat;source+=heat*area;
            add(row,row,diag);VecSetValue(system.rhs,row,source,INSERT_VALUES);
        }
        f["heat_generation"](i,j)=heat_total;
    }
    add(shaft,shaft,shaft_diag);
    const double sink=cfg.shaft_ambient_heat_transfer*cfg.shaft_ambient_area_ratio*2*M_PI*mesh.R*mesh.L;
    if(mesh.rank==mesh.size-1) {
        add(shaft,shaft,sink);VecSetValue(system.rhs,shaft,sink*(shaft_temperature(cfg)-ref),INSERT_VALUES);
    }
    if(solid)for(int i=0;i<nt;++i)for(int slab=0;slab<ns;++slab) {
        const double kb=cfg.bush_thermal_conductivity,dtheta=mesh.get_d_theta();
        const double inner=kb*dtheta*solid_dz/std::log((mesh.R+0.5*dr)/mesh.R);
        add(bore(i,slab),bore(i,slab),inner);add(bore(i,slab),radial(i,slab,0),-inner);
        add(radial(i,slab,0),radial(i,slab,0),inner);add(radial(i,slab,0),bore(i,slab),-inner);
        for(int k=0;k<nr;++k) {
            const PetscInt row=radial(i,slab,k);const double r=mesh.R+(k+0.5)*dr;
            const double angle=kb*dr*solid_dz/(r*dtheta);
            add(row,row,2*angle);add(row,radial(i+1,slab,k),-angle);add(row,radial(i-1,slab,k),-angle);
            const double end_area=r*dr*dtheta;
            if(slab+1<ns) {
                const double g=kb*end_area/solid_dz;
                add(row,row,g);add(row,radial(i,slab+1,k),-g);
                add(radial(i,slab+1,k),radial(i,slab+1,k),g);add(radial(i,slab+1,k),row,-g);
            }
            double source=0;
            if(cfg.bush_thermal_model==BushThermalModel::CONDUCTION_3D && (slab==0||slab==ns-1)) {
                const double old=f["cross_film_solid_"+std::to_string(k)](i,slab*nz/ns);
                const double he=end_heat_transfer(old,solid_dz,cfg,slab==0);
                if(he>0) {
                    const double g=end_area/(solid_dz/(2*kb)+1/he)*(ns==1?2:1);
                    add(row,row,g);source+=g*(outer_temperature(cfg)-ref);
                }
            }
            if(k+1<nr) {
                const double g=kb*dtheta*solid_dz/std::log((r+dr)/r);
                add(row,row,g);add(row,radial(i,slab,k+1),-g);
                add(radial(i,slab,k+1),radial(i,slab,k+1),g);add(radial(i,slab,k+1),row,-g);
            } else {
                double resistance=std::log(ro/r)/(kb*dtheta*solid_dz);
                if(cfg.bush_outer_boundary==BushOuterBoundary::CONVECTION) {
                    resistance=cfg.bush_outer_heat_transfer>0 ?
                        resistance+1.0/(cfg.bush_outer_heat_transfer*ro*dtheta*solid_dz) :
                        std::numeric_limits<double>::infinity();
                }
                const double g=1/resistance;
                add(row,row,g);source+=g*(outer_temperature(cfg)-ref);
            }
            VecSetValue(system.rhs,row,source,INSERT_VALUES);
        }
    }
    MatAssemblyBegin(system.a,MAT_FINAL_ASSEMBLY);MatAssemblyEnd(system.a,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(system.rhs);VecAssemblyEnd(system.rhs);
    KSPCreate(PETSC_COMM_WORLD,&system.ksp);KSPSetOperators(system.ksp,system.a,system.a);
    KSPSetOptionsPrefix(system.ksp,"cross_film_");KSPSetType(system.ksp,KSPBCGSL);
    PC pc;KSPGetPC(system.ksp,&pc);PCSetType(pc,PCBJACOBI);
    KSPSetTolerances(system.ksp,cfg.linear_rtol,1e-14,PETSC_DEFAULT,2000);KSPSetFromOptions(system.ksp);
    KSPSolve(system.ksp,system.rhs,system.x);
    KSPConvergedReason reason;KSPGetConvergedReason(system.ksp,&reason);
    require_valid(reason<0,"Resolved cross-film energy linear solve did not converge");
    PetscInt first,last;VecGetOwnershipRange(system.x,&first,&last);
    const PetscScalar* values=nullptr;VecGetArrayRead(system.x,&values);
    double tshaft=mesh.rank==mesh.size-1?PetscRealPart(values[shaft-first])+ref:0;
    MPI_Bcast(&tshaft,1,MPI_DOUBLE,mesh.size-1,MPI_COMM_WORLD);
    bool invalid=false;
    for(int i=0;i<nt;++i)for(int j=0;j<nz;++j) {
        double mean=0,residual=0;
        for(int k=0;k<n;++k) {
            const double t=PetscRealPart(values[index(i,j,k)-first])+ref;
            invalid|=!(std::isfinite(t)&&t>0);
            residual=std::max(residual,std::abs(t-previous[(i*nz+j)*n+k])/std::max(std::abs(ref),1.0));
            f[layer_name(k)](i,j)=t;mean+=t/n;
        }
        f["temperature"](i,j)=mean;f["cross_film_layer_residual"](i,j)=residual;
        f["cross_film_layer_initialized"](i,j)=1;f["cross_film_t_shaft"](i,j)=tshaft;
        const double gfilm=film_conductance(f,cfg,i,j);
        if(solid) {
            const int slab=j*ns/nz;
            f["cross_film_t_bush"](i,j)=PetscRealPart(values[bore(i,slab)-first])+ref;
            for(int k=0;k<nr;++k) {
                auto& field=f["cross_film_solid_"+std::to_string(k)];
                const double t=PetscRealPart(values[radial(i,slab,k)-first])+ref;
                invalid|=!(std::isfinite(t)&&t>0);
                residual=std::max(residual,std::abs(t-field(i,j))/std::max(std::abs(ref),1.0));
                field(i,j)=t;
            }
            f["cross_film_layer_residual"](i,j)=residual;
            const double last=f["cross_film_solid_"+std::to_string(nr-1)](i,j);
            const double factor=cfg.bush_outer_boundary==BushOuterBoundary::FIXED_TEMPERATURE ? 0 :
                1/(1+cfg.bush_outer_heat_transfer*ro/cfg.bush_thermal_conductivity*std::log(ro/(ro-.5*dr)));
            f["temperature_housing_outer"](i,j)=outer_temperature(cfg)+(last-outer_temperature(cfg))*factor;
        } else {
            const double gb=bush_conductance(cfg);
            f["cross_film_t_bush"](i,j)=(gfilm*f[layer_name(n-1)](i,j)+gb*outer_temperature(cfg))/(gfilm+gb);
        }
    }
    VecRestoreArrayRead(system.x,&values);
    require_valid(invalid,"Resolved cross-film energy returned a nonphysical temperature");
    // Wall fluxes/outputs use this solution immediately. Preserve the frozen
    // viscosity integrals until the next pressure solve, as in the Picard model.
    for(int i=0;i<nt;++i)for(int j=0;j<nz;++j) {
        const double g=film_conductance(f,cfg,i,j);
        const double qs=g*(f[layer_name(0)](i,j)-f["cross_film_t_shaft"](i,j));
        const double qb=g*(f[layer_name(n-1)](i,j)-f["cross_film_t_bush"](i,j));
        f["cross_film_wall_flux"](i,j)=qs+qb;
        if(f.has("temperature_journal"))f["temperature_journal"](i,j)=f["cross_film_t_shaft"](i,j);
        if(f.has("temperature_bearing"))f["temperature_bearing"](i,j)=f["cross_film_t_bush"](i,j);
        if(f.has("heat_flux_journal"))f["heat_flux_journal"](i,j)=-qs;
        if(f.has("heat_flux_bearing"))f["heat_flux_bearing"](i,j)=-qb;
    }
    Communicator comm(mesh);comm.update_ghosts(f["temperature"]);
    for(int k=0;k<n;++k)comm.update_ghosts(f[layer_name(k)]);
}

Energy::Balance coupled_balance(const Fields& f,const Mesh& mesh,const SimulationConfig& cfg) {
    CoupledFluxes flux(f,mesh,cfg);
    const int n=cfg.cross_film_layers,nz=mesh.n_z_local;
    const double area=mesh.cell_volume();
    double sums[7]={}; // heat, wall, transport, feed, work, feed mass, local residual
    auto energy_face=[&](bool theta,int i,int j,int k,int phase,int part) {
        const double mass=flux.face(theta,i,j,k,phase,part);
        if(mass==0)return 0.0;
        int ii=theta?(mass>=0?i-1:i):i,jj=theta?j:(mass>=0?j-1:j);
        double t=0,p=0;
        if(jj<0||jj>=nz) {
            const bool south=jj<0;
            const auto mode=south?cfg.bc_z_south_thermal:cfg.bc_z_north_thermal;
            t=mode==ThermalInflowMode::CONSTANT?cfg.axial_inflow_temperature(south):f[layer_name(k)](ii,std::clamp(jj,0,nz-1));
            p=mode==ThermalInflowMode::CONSTANT?(south?cfg.bc_z_south_val:cfg.bc_z_north_val):f["pressure"](ii,std::clamp(jj,0,nz-1));
        } else {t=f[layer_name(k)](ii,jj);p=f["pressure"](ii,jj);}
        return mass*CoupledCavitation::partial_internal_energy(phase,t,p,cfg);
    };
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<nz;++j) {
        double transport=0,feed=0;
        for(int k=0;k<n;++k)for(int phase=0;phase<3;++phase) {
            for(int part=0;part<3;++part)transport+=energy_face(true,i+1,j,k,phase,part)-energy_face(true,i,j,k,phase,part)+
                energy_face(false,i,j+1,k,phase,part)-energy_face(false,i,j,k,phase,part);
            const double mass=f["coupled_source_"+std::to_string(phase)](i,j)/n;
            const auto* inlet=mesh.inlet_at(cfg,i,j);
            const double t=mass>0&&inlet&&inlet->feeds_fresh_oil?inlet->t_supply:f[layer_name(k)](i,j);
            if(mass!=0)feed+=mass*CoupledCavitation::partial_internal_energy(phase,t,f["pressure"](i,j),cfg);
        }
        const double heat=f["heat_generation"](i,j)*area,wall=f["cross_film_wall_flux"](i,j)*area,work=flux.work(i,j);
        sums[0]+=heat;sums[1]+=wall;sums[2]+=transport;sums[3]+=feed;sums[4]+=work;
        for(int phase=0;phase<3;++phase)sums[5]+=f["coupled_source_"+std::to_string(phase)](i,j);
        sums[6]+=transport+wall-heat-feed+work;
    }
    MPI_Allreduce(MPI_IN_PLACE,sums,7,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    Energy::Balance result;
    result.viscous_heat_W=sums[0];result.wall_heat_W=sums[1];result.axial_enthalpy_out_W=sums[2];
    result.groove_enthalpy_in_W=sums[3];result.pressure_work_W=sums[4];result.feed_mass_kg_s=sums[5];
    result.residual_W=sums[6];result.scale_W=std::abs(sums[0])+std::abs(sums[1])+std::abs(sums[2])+std::abs(sums[3])+std::abs(sums[4]);
    return result;
}

Energy::Balance balance(const Fields& f,const Mesh& mesh,const SimulationConfig& cfg) {
    if(cfg.cavitation_model==CavitationModel::COUPLED)return coupled_balance(f,mesh,cfg);
    Fluxes flux(f,mesh,cfg);
    const int n=cfg.cross_film_layers,nz=mesh.n_z_local;
    const double cp=cfg.liquid_specific_heat,ref=cfg.temperature_reference,area=mesh.cell_volume();
    double values[6]={0,0,0,0,0,0}; // heat, wall, axial, feed, continuity defect, mass
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<nz;++j) {
        const auto* inlet=inlet_at(mesh,cfg,i,j);double net=0;
        for(int k=0;k<n;++k)net+=flux.divergence(i,j,k,inlet!=nullptr);
        values[0]+=f["heat_generation"](i,j)*area;
        values[1]+=f["cross_film_wall_flux"](i,j)*area;
        for(int k=0;k<n;++k) {
            const double t=f[layer_name(k)](i,j)-ref;
            if(!inlet) {
                values[4]+=cp*net/n*t;
                if(j==0) {
                    const double mass=-flux.z(i,-1,k);
                    const double tf=mass>=0||cfg.bc_z_south_thermal!=ThermalInflowMode::CONSTANT?t:cfg.axial_inflow_temperature(true)-ref;
                    values[2]+=cp*mass*tf;
                }
                if(j==nz-1) {
                    const double mass=flux.z(i,j,k);
                    const double tf=mass>=0||cfg.bc_z_north_thermal!=ThermalInflowMode::CONSTANT?t:cfg.axial_inflow_temperature(false)-ref;
                    values[2]+=cp*mass*tf;
                }
            } else values[3]+=cp*net/n*(net>=0?inlet->t_supply-ref:t);
        }
        if(inlet)values[5]+=net;
    }
    MPI_Allreduce(MPI_IN_PLACE,values,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    Energy::Balance result;
    result.viscous_heat_W=values[0];result.wall_heat_W=values[1];result.axial_enthalpy_out_W=values[2];
    result.groove_enthalpy_in_W=values[3];result.mass_imbalance_correction_W=values[4];result.feed_mass_kg_s=values[5];
    result.residual_W=values[2]+values[1]-values[0]-values[3]-values[4];
    result.scale_W=std::abs(values[0])+std::abs(values[1])+std::abs(values[2])+std::abs(values[3])+std::abs(values[4]);
    return result;
}

CrossFilm::WallState diagnostics(const Fields& f,const Mesh& mesh,const SimulationConfig& cfg) {
    const double area=mesh.cell_volume();
    double sums[4]={0,0,0,0},change=0;
    double viscosity_residual=0;
    bool invalid=false;
    double lo=std::numeric_limits<double>::infinity(),hi=-lo;
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<mesh.n_z_local;++j) {
        const double ts=f["cross_film_t_shaft"](i,j),tb=f["cross_film_t_bush"](i,j);
        lo=std::min({lo,ts,tb});hi=std::max({hi,ts,tb});
        for(int k=0;k<cfg.cross_film_layers;++k) {
            const double t=f[layer_name(k)](i,j);lo=std::min(lo,t);hi=std::max(hi,t);
            if(cfg.cavitation_model==CavitationModel::COUPLED && cfg.solution_mode==SolutionMode::STEADY_STATE) {
                try {
                    const double mul=FluidProperties::liquid_solution_viscosity(f["dissolved_gas"](i,j),f["pressure"](i,j),t,cfg);
                    const double theta=f["theta"](i,j),target=1/(theta*mul+(1-theta)*cfg.mu_gas);
                    const double used=f[viscosity_name(k)](i,j);
                    invalid|=!(std::isfinite(target)&&target>0&&std::isfinite(used)&&used>0);
                    viscosity_residual=std::max(viscosity_residual,std::abs(used-target)/target);
                } catch(...) {invalid=true;}
            }
        }
        if(cfg.bush_thermal_model!=BushThermalModel::LOCAL)
            for(int k=0;k<cfg.bush_radial_cells;++k) {
                const double t=f["cross_film_solid_"+std::to_string(k)](i,j);
                lo=std::min(lo,t);hi=std::max(hi,t);
            }
        const double g=film_conductance(f,cfg,i,j);
        sums[0]+=ts*area;sums[1]+=tb*area;
        sums[2]+=g*(f[layer_name(0)](i,j)-ts)*area;
        sums[3]+=g*(f[layer_name(cfg.cross_film_layers-1)](i,j)-tb)*area;
        change=std::max(change,f["cross_film_layer_residual"](i,j));
    }
    require_valid(invalid,"Resolved layer constitutive viscosity is nonphysical");
    MPI_Allreduce(MPI_IN_PLACE,sums,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&lo,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&hi,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&change,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&viscosity_residual,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    CrossFilm::WallState state;
    const double surface=2*M_PI*mesh.R*mesh.L;
    state.t_shaft=sums[0]/surface;state.t_bush_mean=sums[1]/surface;
    state.q_shaft_net=sums[2];state.q_bush_net=sums[3];
    state.profile_t_min=lo;state.profile_t_max=hi;state.profile_residual=change;
    state.viscosity_residual=viscosity_residual;
    state.bush_2d_solved=cfg.bush_thermal_model!=BushThermalModel::LOCAL;
    return state;
}
}
