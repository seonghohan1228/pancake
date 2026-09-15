#include <cmath>
#include <iostream>
#include <limits>

#include <mpi.h>
#include <petsc.h>

#include "config.hpp"
#include "communicator.hpp"
#include "field.hpp"
#include "film_thickness.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

static void check(bool cond, const char* msg) {
    if (!cond) {
        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cerr << "FAIL: " << msg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);

    {
        SimulationConfig cfg;
        cfg.cavitation_model = CavitationModel::GUMBEL;
        cfg.omega = 0.0; // Pure diffusion test
        cfg.bc_z_south_type = BCType::DIRICHLET;
        cfg.bc_z_south_val = 1.1e5;
        cfg.bc_z_north_type = BCType::DIRICHLET;
        cfg.bc_z_north_val = 1.2e5;
        
        // Add an axial groove at theta = 180 deg
        InletConfig groove;
        groove.type = InletConfig::Type::GROOVE;
        groove.theta = 180.0;
        groove.size = 10.0; // 10 degrees wide
        groove.p_supply = 2.0e5;
        cfg.inlets.push_back(groove);

        // Add a circular inlet at theta = 90 deg, z = L/2
        InletConfig hole;
        hole.type = InletConfig::Type::CIRCULAR;
        hole.theta = 90.0;
        hole.z = 0.5 * cfg.L;
        hole.size = 0.002; // 2mm radius
        hole.p_supply = 3.0e5;
        cfg.inlets.push_back(hole);

        check(std::abs(cfg.initial_pressure() - cfg.bc_z_south_val) < 1e-12,
              "Initial pressure should use the outlet boundary pressure, not the first inlet supply pressure");

        cfg.omega = 120.0;
        cfg.omega_ramp_time = 0.4;
        cfg.solution_mode = SolutionMode::TRANSIENT;
        check(std::abs(cfg.omega_at_time(0.0)) < 1e-14,
              "Transient omega ramp should start from zero speed");
        check(std::abs(cfg.omega_at_time(0.2) - 60.0) < 1e-12,
              "Transient omega ramp should interpolate linearly");
        check(std::abs(cfg.omega_at_time(0.4) - 120.0) < 1e-12,
              "Transient omega ramp should reach target speed at ramp time");
        cfg.solution_mode = SolutionMode::STEADY_STATE;
        check(std::abs(cfg.omega_at_time(0.0) - 120.0) < 1e-12,
              "Steady-state mode should ignore startup ramp");
        cfg.omega = 0.0;
        cfg.omega_ramp_time = 0.0;
        cfg.solution_mode = SolutionMode::TRANSIENT;

        Mesh mesh(cfg);
        Fields fields;
        fields.add("pressure", mesh).fill(cfg.p_cav);
        fields.add("h",        mesh).fill(cfg.c);
        fields.add("rho",      mesh).fill(cfg.rho);

        FilmThickness::compute_static(fields["h"], mesh, cfg);
        LinearSystem sys(mesh);

        Reynolds::solve(fields, sys, mesh, cfg);

        const Field& pressure = fields["pressure"];
        const double d_theta = mesh.get_d_theta();
        const double d_z = mesh.get_d_z();
        const double R = mesh.R;

        int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Directly prescribed rows are still solved to the Krylov tolerance.
        const double pin_rtol = 1e-6;

        // Check cells inside the groove
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - M_PI);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            if (dth <= 0.5 * (groove.size * M_PI / 180.0)) {
                for (int j = 0; j < mesh.n_z_local; ++j) {
                    check(std::abs(pressure(i, j) - groove.p_supply) < pin_rtol * groove.p_supply,
                          "Groove pressure mismatch");
                }
            }
        }

        // Check cells inside the hole
        for (int i = 0; i < mesh.n_theta_local; ++i) {
            const double theta_c = (mesh.offset_theta + i + 0.5) * d_theta;
            double dth = std::abs(theta_c - 0.5 * M_PI);
            if (dth > M_PI) dth = 2.0 * M_PI - dth;

            for (int j = 0; j < mesh.n_z_local; ++j) {
                const double z_c = (j + 0.5) * d_z;
                const double dist_sq = (R * dth) * (R * dth) + (z_c - hole.z) * (z_c - hole.z);
                if (dist_sq <= hole.size * hole.size) {
                    check(std::abs(pressure(i, j) - hole.p_supply) < pin_rtol * hole.p_supply,
                          "Circular hole pressure mismatch");
                }
            }
        }

        const auto reference = pressure.data;
        Field indicator("indicator",mesh), shifted("shifted",mesh);
        FilmThickness::compute_inlet_indicator(indicator,mesh,cfg);
        for (double turns : {-720.0,1440.0}) {
            auto wrapped=cfg;
            for (auto& inlet : wrapped.inlets) inlet.theta+=turns;
            FilmThickness::compute_inlet_indicator(shifted,mesh,wrapped);
            check(indicator.data==shifted.data,"Inlet mask must be invariant under full turns");
            fields["pressure"].fill(cfg.p_cav);
            fields["rho"].fill(cfg.rho);
            Reynolds::solve(fields,sys,mesh,wrapped);
            for (std::size_t k=0;k<reference.size();++k)
                check(std::abs(pressure.data[k]-reference[k])<0.01,"Wrapped inlet pressure must be unchanged");
        }
        auto overlapping=cfg;
        overlapping.inlets.push_back(groove);
        overlapping.inlets.back().p_supply=4e5;
        for (int i=0;i<mesh.n_theta_local;++i) for (int j=0;j<mesh.n_z_local;++j)
            if (mesh.inside_inlet(groove,i,j))
                check(mesh.inlet_at(overlapping,i,j)==&overlapping.inlets.back(),"Last matching reservoir wins");
        // An undersized partition must fail on every rank before communication.
        auto small=cfg; small.n_theta_global=mesh.size;
        Mesh small_mesh(small); Field small_field("small",small_mesh); Communicator comm(small_mesh);
        bool rejected=false;
        try { comm.update_ghosts(small_field); } catch (const std::invalid_argument&) { rejected=true; }
        check(rejected,"Ghost exchange must reject fewer physical cells than ghost layers");
        auto invalid=cfg;
        invalid.inlets.front().theta=std::numeric_limits<double>::quiet_NaN();
        std::vector<std::string> errors,warnings;
        invalid.validate(errors,warnings);
        bool bad_angle=false;
        for (const auto& error : errors) bad_angle |= error.find(".theta must be a finite")!=std::string::npos;
        check(bad_angle,"Nonfinite inlet angles must fail validation");
        if (rank == 0) std::cout << "PASS: test_inlets\n";
    }

    PetscFinalize();
    return 0;
}
