#pragma once
#include "fbns.hpp"
#include "coupled_cavitation.hpp"
#include "journal_motion.hpp"
#include "reynolds.hpp"

/// Algebraic load closure at frozen thermal/property coefficients; no physical time integration.
namespace SteadyEquilibrium {
struct Stats {
    double force_residual = 0.0;
    double position_change = 0.0; // correction / clearance
    int pressure_solves = 0;
    bool accepted = false;
    FBNS::Stats pressure;
    CoupledCavitation::Stats coupled;
    Reynolds::ForceComponents forces;
};
/// Solve base pressure and accept at most one bounded Cartesian Newton correction.
/// Rejected probes restore all fields and the position. Temperature is never solved here.
Stats update(Fields& fields, LinearSystem& sys, Communicator& comm, const Mesh& mesh,
             const SimulationConfig& cfg, JournalMotion::BearingState& state,
             CoupledCavitation::Solver* coupled_solver = nullptr);
/// Relative residual of fluid force plus independent external load.
double force_residual(const Reynolds::ForceComponents& forces, const SimulationConfig& cfg);
}
