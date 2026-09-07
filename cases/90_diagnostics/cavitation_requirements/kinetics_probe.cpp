#include <fstream>
#include <iomanip>
#include <petsc.h>
#include "fluid_properties.hpp"

// Prescribed-pressure parcel: isolates kinetic choice from the unresolved
// bearing pressure/volume coupling. Rates are illustrative, not FW68D data.
int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (argc < 2) { PetscFinalize(); return 1; }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {
        std::ofstream out;
        if (!rank) {
            out.open(argv[1]);
            out << std::setprecision(17) << "absorption_rate,time,pressure,free_mass,alpha_gas,oil_mass,total_gas_mass\n";
        }
        for (double absorption : {20., 0.}) {
            SimulationConfig cfg;
            cfg.n_theta_global = 24; cfg.n_z_global = 4;
            cfg.c = 1e-4; cfg.rho = 950;
            cfg.fluid_property_model = FluidPropertyModel::TWO_PHASE;
            cfg.cavitation_model = CavitationModel::FBNS;
            cfg.dissolved_gas_species = DissolvedGasSpecies::AIR;
            cfg.density_model = DensityModel::CONSTANT;
            cfg.dissolved_gas_initial = 0.002; cfg.dissolved_gas_max = 0.5;
            cfg.dissolved_gas_henry_coeff = 1e-8;
            cfg.gas_mass_transfer_rate = 20;
            cfg.gas_absorption_rate = absorption;
            cfg.gas_alpha_max = 0.9;
            Mesh mesh(cfg);
            Fields f;
            f.add("pressure", mesh).fill(1e5);
            f.add("temperature", mesh).fill(313.15);
            f.add("h", mesh).fill(cfg.c);
            f.add("theta", mesh).fill(0.5);
            f.add("rho_liquid_solution", mesh).fill(cfg.rho);
            f.add("rho", mesh).fill(cfg.rho * 0.5);
            f.add("dissolved_gas", mesh).fill(cfg.dissolved_gas_initial);
            f.add("free_gas_mass", mesh).fill(0);
            f.add("alpha_gas", mesh).fill(0);
            for (int step = 0; step <= 400; ++step) {
                const double time = step * 1e-4;
                // Resolve a 10-ms pressure recovery after 20 ms of release.
                const double pressure = 1e5 + 3e5 * std::clamp((time - 0.02) / 0.01, 0., 1.);
                f["pressure"].fill(pressure);
                FluidProperties::update_gas_state(f, mesh, cfg, step ? 1e-4 : 0);
                if (!rank) {
                    const double ml = f["rho_liquid_solution"](0,0) * f["theta"](0,0) * cfg.c;
                    const double md = ml * f["dissolved_gas"](0,0);
                    out << absorption << ',' << time << ',' << pressure << ',' << f["free_gas_mass"](0,0)
                        << ',' << f["alpha_gas"](0,0) << ',' << ml - md << ','
                        << md + f["free_gas_mass"](0,0) << '\n';
                }
            }
        }
    }
    PetscFinalize();
}
