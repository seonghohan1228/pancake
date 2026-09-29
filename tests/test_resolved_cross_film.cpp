#include <cmath>
#include <iostream>
#include <petsc.h>
#include "communicator.hpp"
#include "cross_film.hpp"
#include "resolved_cross_film.hpp"
#include "fbns.hpp"
#include "film_thickness.hpp"

void check(bool good, const char* message) {
    if (!good) { std::cerr << message << '\n'; MPI_Abort(MPI_COMM_WORLD,1); }
}

SimulationConfig config(int layers) {
    SimulationConfig cfg;
    cfg.n_theta_global=24;cfg.n_z_global=4;cfg.R=.0215;cfg.L=.04;cfg.c=42.5e-6;
    cfg.omega=157.07963267948966;cfg.mu=.04;cfg.rho=920;
    cfg.solution_mode=SolutionMode::STEADY_STATE;cfg.cavitation_model=CavitationModel::FBNS;
    cfg.cross_film_thermal_model=CrossFilmModel::FINITE_VOLUME;cfg.cross_film_layers=layers;
    cfg.temperature_model=TemperatureModel::ENERGY_EQUATION;
    cfg.groove_mixing_model=GrooveMixingModel::FLUX_BALANCED;
    cfg.temperature_initial=cfg.temperature_reference=cfg.thermal_ambient_temperature=300;
    cfg.bush_outer_temperature=cfg.shaft_ambient_temperature=300;
    cfg.thermal_conductivity=.13;cfg.liquid_specific_heat=2010;
    cfg.bush_outer_radius=.0425;cfg.bush_thermal_conductivity=15;cfg.bush_outer_heat_transfer=65;
    cfg.bush_outer_boundary=BushOuterBoundary::CONVECTION;
    cfg.bush_thermal_model=BushThermalModel::CONDUCTION_2D;cfg.bush_radial_cells=8;
    cfg.shaft_ambient_heat_transfer=1300;cfg.shaft_ambient_area_ratio=.5375;
    cfg.p_cav=0;cfg.bc_z_south_val=cfg.bc_z_north_val=101325;
    cfg.linear_rtol=1e-11;cfg.fbns_tol=1e-10;cfg.fbns_max_iters=150;
    return cfg;
}

void register_fields(Fields& f, const Mesh& mesh, const SimulationConfig& cfg) {
    f.add("h",mesh).fill(cfg.c);f.add("dh_dt",mesh).fill(0);
    f.add("pressure",mesh).fill(101325);f.add("theta",mesh).fill(1);
    f.add("rho",mesh).fill(cfg.rho);f.add("rho_liquid_solution",mesh).fill(cfg.rho);
    f.add("mu",mesh).fill(cfg.mu);f.add("temperature",mesh).fill(300);
    f.add("k_liquid_solution",mesh).fill(cfg.thermal_conductivity);
    f.add("heat_generation",mesh).fill(0);f.add("inlet_indicator",mesh).fill(0);
    CrossFilm::register_fields(f,mesh,cfg);
}

double analytic_test(int layers) {
    auto cfg=config(layers);Mesh mesh(cfg);Communicator comm(mesh);Fields fields;
    register_fields(fields,mesh,cfg);comm.update_ghosts(fields);
    CrossFilm::update(fields,mesh,cfg);ResolvedCrossFilm::solve(fields,mesh,cfg);
    const auto state=ResolvedCrossFilm::diagnostics(fields,mesh,cfg);
    const auto ledger=ResolvedCrossFilm::balance(fields,mesh,cfg);
    const double speed=cfg.omega*cfg.R,heat=cfg.mu*speed*speed/cfg.c;
    const double gs=cfg.shaft_ambient_heat_transfer*cfg.shaft_ambient_area_ratio;
    const double gb=1/(cfg.R*std::log(cfg.bush_outer_radius/cfg.R)/cfg.bush_thermal_conductivity+
                      cfg.R/(cfg.bush_outer_radius*cfg.bush_outer_heat_transfer));
    const double g=cfg.thermal_conductivity/cfg.c,det=(gs+g)*(gb+g)-g*g;
    const double ts=300+heat*.5*(gb+2*g)/det,tb=300+heat*.5*(gs+2*g)/det;
    check(std::abs(state.t_shaft-ts)<2e-6 && std::abs(state.t_bush_mean-tb)<2e-6,
          "uniform Couette conjugate wall temperatures must match the analytic resistance balance");
    check(std::abs(ledger.residual_W)/ledger.scale_W<1e-8,"resolved energy must conserve heat");
    check(state.profile_t_min>=300-1e-8,"positive heating with ambient sinks must not undershoot ambient");
    double error=0;
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<mesh.n_z_local;++j)for(int k=0;k<layers;++k) {
        const double eta=(k+.5)/layers;
        const double exact=ts+(tb-ts)*eta+cfg.mu*speed*speed/(2*cfg.thermal_conductivity)*eta*(1-eta);
        error=std::max(error,std::abs(fields["cross_film_layer_"+std::to_string(k)](i,j)-exact));
    }
    MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    CrossFilm::update(fields,mesh,cfg);ResolvedCrossFilm::solve(fields,mesh,cfg);
    check(ResolvedCrossFilm::diagnostics(fields,mesh,cfg).profile_residual<1e-8,
          "constant-property resolved thermal fixed point must remain unchanged");
    if(mesh.rank==0)std::cout<<"layers="<<layers<<" analytic_error_K="<<error<<" energy_residual_W="<<ledger.residual_W<<'\n';
    return error;
}

void groove_test(double omega) {
    auto cfg=config(8);cfg.n_theta_global=48;cfg.n_z_global=8;cfg.omega=omega;
    cfg.cavitated_film_model=CavitatedFilmModel::STRIATED;
    InletConfig groove;groove.type=InletConfig::Type::GROOVE;
    groove.theta=90;groove.size=15;groove.p_supply=101325;groove.feeds_fresh_oil=true;groove.t_supply=300;
    cfg.inlets={groove};Mesh mesh(cfg);Communicator comm(mesh);Fields fields;
    register_fields(fields,mesh,cfg);
    for(int i=0;i<mesh.n_theta_local;++i)for(int j=0;j<mesh.n_z_local;++j) {
        const double theta=(mesh.offset_theta+i+.5)*mesh.get_d_theta();
        fields["h"](i,j)=cfg.c*(1+.65*std::sin(theta)-.15*std::cos(theta));
    }
    FilmThickness::compute_inlet_indicator(fields["inlet_indicator"],mesh,cfg);
    comm.update_ghosts(fields);CrossFilm::update(fields,mesh,cfg);
    LinearSystem sys(mesh,linear_solver_options_from_config(cfg));
    check(FBNS::solve(fields,sys,comm,mesh,cfg).converged,"groove pressure regression must converge");
    comm.update_ghosts(fields);ResolvedCrossFilm::solve(fields,mesh,cfg);
    const auto state=ResolvedCrossFilm::diagnostics(fields,mesh,cfg);
    const auto ledger=ResolvedCrossFilm::balance(fields,mesh,cfg);
    check(state.profile_t_min>=300-1e-7,"groove/reverse-flow thermal solution must not create subambient cells");
    check(std::isfinite(state.profile_t_max)&&state.profile_t_max<400,"groove thermal field must remain finite");
    check(std::abs(ledger.residual_W)/ledger.scale_W<1e-8,"groove/layer/boundary enthalpy must close globally");
    if(mesh.rank==0)std::cout<<"omega="<<omega<<" min_K="<<state.profile_t_min<<" max_K="<<state.profile_t_max
                           <<" energy_residual_W="<<ledger.residual_W<<'\n';
}

void housing_ends_test() {
    auto cfg=config(8);cfg.n_z_global=16;cfg.bush_axial_cells=16;
    cfg.bush_thermal_model=BushThermalModel::CONDUCTION_3D;
    Mesh mesh(cfg);Communicator comm(mesh);Fields fields;
    register_fields(fields,mesh,cfg);comm.update_ghosts(fields);
    CrossFilm::update(fields,mesh,cfg);ResolvedCrossFilm::solve(fields,mesh,cfg);
    const double insulated=ResolvedCrossFilm::diagnostics(fields,mesh,cfg).t_bush_mean;
    cfg.bush_end_heat_transfer=80;
    CrossFilm::update(fields,mesh,cfg);ResolvedCrossFilm::solve(fields,mesh,cfg);
    const auto state=ResolvedCrossFilm::diagnostics(fields,mesh,cfg);
    check(state.t_bush_mean<insulated,"exposed housing ends must cool the housing");
    check(fields["cross_film_solid_0"](0,8)>fields["cross_film_solid_0"](0,0),
          "housing midplane must be warmer than cooled ends under uniform heating");
    const double dr=(cfg.bush_outer_radius-cfg.R)/cfg.bush_radial_cells;
    double loss=0;
    for(int i=0;i<mesh.n_theta_local;++i) {
        for(int j=0;j<mesh.n_z_local;++j)
            loss+=cfg.bush_outer_heat_transfer*(fields["temperature_housing_outer"](i,j)-300)*
                  cfg.bush_outer_radius*mesh.get_d_theta()*mesh.get_d_z();
        for(int k=0;k<cfg.bush_radial_cells;++k)for(int j:{0,mesh.n_z_local-1})
            loss+=(fields["cross_film_solid_"+std::to_string(k)](i,j)-300)*
                (cfg.R+(k+.5)*dr)*dr*mesh.get_d_theta()/
                (mesh.get_d_z()/(2*cfg.bush_thermal_conductivity)+1/cfg.bush_end_heat_transfer);
    }
    MPI_Allreduce(MPI_IN_PLACE,&loss,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    check(std::abs(loss-state.q_bush_net)/loss<1e-8,
          "housing film heat input must equal cylindrical plus annular-end heat loss");
    cfg.bush_end_natural_convection=true;cfg.density_model=DensityModel::LINEAR_TEMPERATURE;
    cfg.density_temperature_intercept=954.9;cfg.density_temperature_slope=-.7284;
    for(int iteration=0;iteration<12;++iteration) {
        CrossFilm::update(fields,mesh,cfg);ResolvedCrossFilm::solve(fields,mesh,cfg);
    }
    const auto natural=ResolvedCrossFilm::diagnostics(fields,mesh,cfg);
    check(natural.profile_residual<1e-8 && natural.profile_t_min>=300,
          "natural convection at both housing ends must converge without subambient temperatures");
    if(mesh.rank==0)std::cout<<"housing_heat_W="<<loss<<" natural_profile_residual="<<natural.profile_residual<<'\n';
}

int main(int argc,char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    {
        const double e8=analytic_test(8),e16=analytic_test(16),e32=analytic_test(32);
        check(e16<.27*e8&&e32<.27*e16,"through-film temperature must converge at second order");
        groove_test(157.07963267948966);groove_test(-157.07963267948966);
        housing_ends_test();
    }
    PetscFinalize();return 0;
}
