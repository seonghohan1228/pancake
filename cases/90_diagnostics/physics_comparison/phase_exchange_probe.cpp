#include <fstream>
#include <iomanip>
#include <petsc.h>

#include "config.hpp"
#include "field.hpp"
#include "fluid_properties.hpp"
#include "mesh.hpp"

// Isolate local exchange at prescribed p,T, without boundary transport. The same
// source is compiled against both revisions; it measures inventories directly.
int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (argc < 2) { PetscFinalize(); return 1; }
    {
        SimulationConfig cfg;
        cfg.R = 0.02; cfg.L = 0.04; cfg.c = 1e-4;
        cfg.n_theta_global = 24; cfg.n_z_global = 4;
        cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
        cfg.cavitation_model = CavitationModel::FBNS;
        cfg.dissolved_gas_species = DissolvedGasSpecies::AIR;
        cfg.density_model = DensityModel::CONSTANT;
        cfg.liquid_viscosity_model = LiquidViscosityModel::CONSTANT;
        cfg.rho = 950; cfg.mu = 0.05; cfg.mu_gas = 1.8e-5;
        cfg.dissolved_gas_initial = 0.002;
        cfg.dissolved_gas_henry_coeff = 1e-8;
        cfg.dissolved_gas_henry_temp_coeff = 0;
        cfg.dissolved_gas_max = 0.5;
        cfg.gas_mass_transfer_rate = 20;
        cfg.gas_alpha_max = 0.6;
        cfg.dt = 1e-4;
        Mesh mesh(cfg);
        Fields fields;
        fields.add("pressure", mesh).fill(1e5);
        fields.add("temperature", mesh).fill(313.15);
        fields.add("h", mesh).fill(cfg.c);
        fields.add("theta", mesh).fill(0.5);
        fields.add("rho_liquid_solution", mesh).fill(cfg.rho);
        fields.add("rho", mesh).fill(0.5 * cfg.rho);
        fields.add("mu", mesh).fill(cfg.mu);
        fields.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
        fields.add("free_gas_mass", mesh).fill(0);
        fields.add("alpha_gas", mesh).fill(0);
        fields.add("rho_gas", mesh).fill(0);
        fields.add("gas_mass_transfer", mesh).fill(0);
        fields.add("inlet_indicator", mesh).fill(0);
        FluidProperties::update_solution_fields(fields, mesh, cfg);
        std::ofstream output;
        if (rank == 0) {
            output.open(argv[1]);
            output << std::setprecision(17)
                   << "time,liquid_mass_area,dissolved_mass_area,free_mass_area,oil_mass_area,"
                      "total_mass_area,total_gas_area,alpha_gas,theta,concentration\n";
        }
        for (int step = 0; step <= 200; ++step) {
            if (rank == 0) {
                const double liquid = fields["rho_liquid_solution"](0, 0) * cfg.c * fields["theta"](0, 0);
                const double dissolved = fields["dissolved_gas"](0, 0) * liquid;
                const double free_mass = fields["free_gas_mass"](0, 0);
                output << step * cfg.dt << ',' << liquid << ',' << dissolved << ',' << free_mass
                       << ',' << liquid - dissolved << ',' << liquid + free_mass
                       << ',' << dissolved + free_mass << ',' << fields["alpha_gas"](0, 0)
                       << ',' << fields["theta"](0, 0) << ',' << fields["dissolved_gas"](0, 0) << '\n';
            }
            if (step == 200) break;
            FluidProperties::update_gas_state(fields, mesh, cfg, cfg.dt);
            FluidProperties::update_solution_fields(fields, mesh, cfg);
        }
    }
    PetscFinalize();
    return 0;
}
