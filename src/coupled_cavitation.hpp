#pragma once

#include <memory>
#include <array>
#include "communicator.hpp"
#include "config.hpp"
#include "field.hpp"

/// Conservative homogeneous oil/dissolved/free-species film model.
/// Oil vapor is in equilibrium at p_v(T) with negligible mass. Gas occupies the
/// remaining volume at partial pressure p-p_v(T); no prescribed gas plateau.
/// Optional film-averaged temperature follows conservative internal energy.
namespace CoupledCavitation {
/// Partial internal energies [J/kg], relative to the configured caloric datum.
/// The dissolved offset is an explicit binding-energy analogue, not latent heat data.
double partial_internal_energy(int phase, double temperature, double pressure, const SimulationConfig&);
std::array<double,3> internal_energies(double temperature, double pressure, const SimulationConfig&);
/// Scaled round-off slack for admissibility checks on a state that the coupled
/// Newton solve accepted at coupled_nonlinear_tolerance, clamped to [1e-9, 1e-8].
double admissibility_roundoff(double coupled_nonlinear_tolerance);
/// True when p-p_v >= -admissibility_roundoff*pressure_scale [Pa]. The coupled
/// residual measures pressure in units of pressure_scale, so a converged
/// cavitated cell may sit that far below p_v without being physically invalid.
bool pressure_admissible(double pressure, double vapor_pressure, double pressure_scale,
                         double coupled_nonlinear_tolerance);
struct Stats {
    bool converged = false;
    int iterations = 0;
    double residual = 0;
    double volume_error = 0;
    double oil_balance = 0;       // global storage + boundary flux [kg/s]
    double species_balance = 0;   // dissolved + free [kg/s]
    double liquid_mass = 0, gas_mass = 0;
    double liquid_boundary_flux = 0, gas_boundary_flux = 0;
    double liquid_inlet_source = 0, gas_inlet_source = 0;
    double liquid_balance = 0;
    double energy_balance = 0;    // storage + transport + work - heat [W]
    double internal_energy = 0;   // global mixture internal energy [J]
};

class Solver {
public:
    Solver(const Mesh& mesh, Communicator& comm, bool thermal = false);
    ~Solver();
    Solver(const Solver&) = delete;
    Solver& operator=(const Solver&) = delete;
    /// Build conserved inventories from theta, c and free_gas_mass. Store the
    /// initial history. Calling solve never advances old-time storage itself.
    void initialize(Fields& fields, const SimulationConfig& cfg);
    /// Backward Euler conservation and volume closure, solved simultaneously.
    /// equilibrium=true imposes local phase equilibrium instead of finite rates.
    /// Disable steady continuation for disposable position trials; rejected
    /// trials then return without entering the conservative pseudo-time march.
    Stats solve(Fields& fields, const SimulationConfig& cfg, double dt, bool equilibrium = false,
                bool allow_steady_continuation = true);
    /// Commit accepted inventories as the next step's fixed old state.
    static void store_old_time(Fields& fields);
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
}
