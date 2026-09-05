#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Transport of dissolved and released gas in the liquid film.
///
/// This is the advective-diffusive transport half of the gaseous-cavitation model
/// (JFO Phase B.2). It moves the dissolved-gas mass fraction `dissolved_gas` (c_d,
/// kg gas / kg liquid solution) conservatively. FBNS uses its exact discrete
/// liquid mass flux; JFO uses the reconstructed liquid velocity. Free gas uses a
/// homogeneous no-slip approximation to the film-averaged
/// velocity when present. Pressure inlets impose the configured incoming mixture:
/// `dissolved_gas_initial` and zero free gas. Axial pressure-boundary outflow
/// convects gas out, while inflow carries the same initial dissolved composition.
/// The local release/resorption exchange between dissolved and free gas
/// (Henry's-law-driven) is handled separately by FluidProperties::update_gas_state
/// (operator splitting: transport, then react).
namespace GasTransport {

    /// Local-rank net gas inflow [kg/s] through axial boundaries, using the same
    /// dissolved and free-gas fluxes as transport. A completed step retains its
    /// boundary flux before the subsequent phase/property update.
    double axial_gas_inflow(const Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

    /// Advect-diffuse gas over one step. No-op unless the fluid model is
    /// TWO_PHASE. Reads `pressure`, `h`, `theta`, `rho`, `mu`; updates
    /// `dissolved_gas` and optional `free_gas_mass` in place.
    void solve(Fields& fields, LinearSystem& sys, const Mesh& mesh,
               const SimulationConfig& cfg, double dt);

}  // namespace GasTransport
