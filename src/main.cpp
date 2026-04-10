#include <cmath>
#include <fstream>
#include <mpi.h>

#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include "utils.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    SimulationConfig cfg;
    Mesh mesh(cfg);
    Communicator comm(mesh);
    IO::prepare_output_directory(cfg, MPI_COMM_WORLD);

    // Initialize fields
    Fields fields;
    fields.add("pressure", mesh).fill(101325.0);
    fields.add("h", mesh).fill(cfg.c);
    fields.add("rho", mesh).fill(cfg.rho);

    Utils::log("Starting Simulation...");

    // Geometry initialization
    double psi_rad = cfg.attitude_angle_deg * (M_PI / 180.0);
    Field& h = fields["h"];
    for(int j = 0; j < mesh.n_z_local; ++j) {
        for(int i = 0; i < mesh.n_theta_local; ++i) {
            double theta_coord = (mesh.offset_theta + i + 0.5) * mesh.get_d_theta();
            h(i,j) = cfg.c - cfg.e * std::cos(theta_coord - psi_rad);
        }
    }
    comm.update_ghosts(fields);

    // Initialize old_time to current state
    fields["pressure"].store_old_time();

    double t = 0.0;
    double output_time = 0;
    int step = 0;

    // Time loop
    while(t <= cfg.end_t || std::abs(t - cfg.end_t) < 1e-8) {
        // Geometry & physics update
        Field& p = fields["pressure"];

        // Update h and sync
        // TODO: Solve position and update h
        comm.update_ghosts(fields);

        // Assemble & solve equations
        // Reynolds equation
        // d(rho*theta)/dt + div(flux_u * theta) - div(gamma * grad(theta)) = 0

        // Solve

        // Output
        if(t >= output_time || std::abs(t - output_time) < 1e-8 || step == 0) {
             IO::write_timestep(t, step, mesh, fields, cfg);
             Utils::log("Step " + std::to_string(step) + " t=" + std::to_string(t), Utils::Color::GREEN, 1);
             output_time += cfg.write_interval;
        }

        // Advance time
        p.store_old_time(); // Save current solution for next ddt step
        t += cfg.dt;
        step++;
    }

    // Finalization
    // Close the PVD file tags
    if (mesh.rank == 0) { // Consistency with IO
        std::ofstream pvd(cfg.output_dir + "/results.pvd", std::ios::app);
        pvd << "  </Collection>\n"
            << "</VTKFile>\n";
    }

    Utils::log("Simulation finished successfully.");
    MPI_Finalize();
    return 0;
}
