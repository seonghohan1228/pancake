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

    struct Balance {
        double residual_W = 0.0;
        double scale_W = 0.0;
        double viscous_heat_W = 0.0;
        double wall_heat_W = 0.0;
        double pressure_work_W = 0.0; // p div(Q) + diffusion displacement - feed pV [W]
        double axial_enthalpy_out_W = 0.0;
        double groove_enthalpy_in_W = 0.0;
        double mass_imbalance_correction_W = 0.0;
        double feed_mass_kg_s = 0.0;
    };

    void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh, const SimulationConfig& cfg);

    /// Discrete steady FBNS energy ledger, with enthalpy referenced to
    /// temperature_reference. Requires current temperature ghosts after solve.
    Balance balance(const Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace Energy
