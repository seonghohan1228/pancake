#pragma once

#include "config.hpp"
#include "field.hpp"

namespace JournalMotion {

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct BearingState {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
};

/// Build the initial moving-bearing state. When configured from attitude, the
/// state is the bearing-center displacement that reproduces the static film gap.
BearingState initial_state(const SimulationConfig& cfg);

/// External load is constant in the global frame and independent of bearing motion.
Vector3 external_load(const SimulationConfig& cfg);

/// Read the integrated fluid force field used to move the bearing.
Vector3 fluid_force_from_fields(const Fields& fields);

/// Advance the bearing state by one timestep using cfg.motion_time_method.
void advance(BearingState& state, Vector3 fluid_force,
             const SimulationConfig& cfg, double dt);

/// Equivalent shaft eccentricity magnitude for the current moving-bearing state.
double equivalent_eccentricity(const BearingState& state);

/// Equivalent attitude angle of the minimum gap, in degrees.
double equivalent_attitude_angle_deg(const BearingState& state,
                                     const SimulationConfig& cfg);

/// Fill uniform output fields for bearing state and independent external load.
void write_state_fields(Fields& fields, const BearingState& state,
                        const SimulationConfig& cfg);

}  // namespace JournalMotion
