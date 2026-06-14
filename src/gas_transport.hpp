#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Transport of dissolved and released gas in the liquid film.
///
/// This is the advective-diffusive transport half of the gaseous-cavitation model
/// (JFO Phase B.2). It moves the dissolved-gas mass fraction `dissolved_gas` (c_d,
/// kg gas / kg liquid solution) with the SAME film mass flux that the Elrod-Adams
/// solve conserves; `free_gas_mass` is transported with the film-averaged
/// velocity when present. Pressure inlets impose the configured incoming mixture:
/// `dissolved_gas_initial` and zero free gas. Axial pressure-boundary outflow
/// convects gas out, while inflow carries the same initial dissolved composition.
/// The local release/resorption exchange between dissolved and free gas
/// (Henry's-law-driven) is handled separately by FluidProperties::update_gas_state
/// (operator splitting: transport, then react).
namespace GasTransport {

    /// Advect-diffuse gas over one step. No-op unless the fluid model is
    /// GAS_CAVITATION_MIXTURE. Reads `pressure`, `h`, `theta`, `rho`, `mu`; updates
    /// `dissolved_gas` and optional `free_gas_mass` in place.
    void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh,
               const SimulationConfig& cfg, double dt);

}  // namespace GasTransport
