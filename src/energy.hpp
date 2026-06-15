#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Film-averaged energy equation for the lubricant.
///
/// The temperature field is cell-centered. `heat_generation` stores the local
/// viscous dissipation per bearing surface area [W/m^2]. When the thermal model
/// is ISOTHERMAL, the solver only refreshes `heat_generation` and leaves
/// `temperature` at the configured initial/reference value.
namespace Energy {

    void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace Energy
