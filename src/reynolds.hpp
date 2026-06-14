#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Reynolds equation solvers for journal bearings.
///
/// GUMBEL is pressure-primary and can clamp pressure after the solve.
/// Required fields: "pressure", "h", "rho"; optional "mu".
///
/// ELROD_ADAMS is theta-primary, where theta is the JFO liquid-content
/// variable. The Elrod-Adams switch disables Poiseuille diffusion in cavitated
/// cells. Required fields: "theta", "pressure", "h", "rho"; optional
/// field-valued "rho_liquid_solution" and "mu".
///
/// Ghost cells for geometry, flow, and optional property fields must be synced
/// by Communicator before calling the solvers.
namespace Reynolds {

    struct ForceComponents {
        double pressure_x = 0.0;
        double pressure_y = 0.0;
        double pressure_z = 0.0;
        double viscous_x = 0.0;
        double viscous_y = 0.0;
        double viscous_z = 0.0;
        double fluid_x = 0.0;
        double fluid_y = 0.0;
        double fluid_z = 0.0;
        double friction_torque = 0.0;
    };

    /// Convergence summary of the Elrod-Adams outer (flag-update) iteration.
    /// flag_flips counts switch-function changes summed over all outer iterations.
    struct ElrodStats {
        int outer_iters = 0;
        int flag_flips = 0;
        double max_theta_change = 0.0;
        bool converged = false;
    };

    /// Gumbel pressure-primary solver. See PHYSICS.md for the discretization.
    void solve(Fields& fields, LinearSystem& sys,
               const Mesh& mesh, const SimulationConfig& cfg);

    /// Elrod-Adams theta-primary solver. Sets theta, pressure, and rho.
    /// If a "cavitation_clamp_mass" field is present, it receives the liquid mass
    /// [kg per cell] added by the final theta_min clamp and full-film snap so the
    /// diagnostics mass balance can account for them.
    ElrodStats solve_elrod(Fields& fields, LinearSystem& sys,
                           const Mesh& mesh, const SimulationConfig& cfg);

    /// Calculate velocities on a staggered grid.
    /// u_theta is on FACE_THETA, u_z is on FACE_Z.
    void calculate_velocities(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

    /// Calculate integrated bearing-surface forces and shaft friction torque.
    /// pressure_* is pressure force on the moving bearing, viscous_* is shear
    /// force on the moving bearing, and fluid_* is their sum.
    ForceComponents calculate_macroscopic_properties(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace Reynolds
