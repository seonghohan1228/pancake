#include "fluid_properties.hpp"
#include "property_table_io.hpp"
#include <array>
#include <cmath>
#include <iostream>

int main() {
    SimulationConfig cfg;
    cfg.fluid_property_model=FluidPropertyModel::SINGLE_PHASE;
    cfg.liquid_viscosity_model=LiquidViscosityModel::COMPOSITION_TABLE;
    cfg.viscosity_composition_file=std::filesystem::path(PANCAKE_SOURCE_DIR)/"cases/ref_data/fvc68b_r410a_1999_viscosity.csv";
    std::vector<std::string> errors;
    FluidProperties::load_config_property_tables(cfg,errors);
    if (!errors.empty()) { std::cerr << errors.front(); return 1; }
    // Published FVC68B data verify loading and interpolation, not FW68D identity.
    const double measured[3][4]={{.16132,.027,.018,.0078},{.05159,.012,.008,.0044},{.02183,.0065,.0048,.0029}};
    const double concentrations[]={0,.15,.2,.3}, temperatures[]={296.15,318.15,338.15};
    for (int j=0;j<3;++j) for (int i=0;i<4;++i) {
        const double mu=FluidProperties::liquid_solution_viscosity(concentrations[i],1e5,temperatures[j],cfg);
        if (std::abs(mu/measured[j][i]-1)>1e-13) return 2;
    }
    const double expected=std::pow(.027*.018*.012*.008,.25);
    if (std::abs(FluidProperties::liquid_solution_viscosity(.175,1e5,307.15,cfg)/expected-1)>1e-13) return 3;
    for (auto query : {std::array<double,2>{.1,293.15}, {.1,353.15}, {.31,313.15}}) {
        bool rejected=false;
        try { FluidProperties::liquid_solution_viscosity(query[0],1e5,query[1],cfg); }
        catch (const std::domain_error&) { rejected=true; }
        if (!rejected) return 4;
    }
    SimulationConfig copy; copy.load_from_text(cfg.to_config_text());
    if (copy.liquid_viscosity_model!=cfg.liquid_viscosity_model || copy.viscosity_composition_file!=cfg.viscosity_composition_file) return 5;
    std::cout << "PASS: actual-composition viscosity, exact source nodes, log interpolation and domain rejection\n";
}
