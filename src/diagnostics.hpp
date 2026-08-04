#pragma once

#include <filesystem>

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"
#include "reynolds.hpp"

/// Per-step global conservation and convergence diagnostics.
///
/// The liquid balance mirrors the discrete operators the Elrod solve assembles
/// (units kg/s per cell): storage from theta vs theta.old, boundary fluxes from
/// the same Dirichlet/INLET_OUTLET formulas, inlet penalty cells re-evaluated as
/// physical operator residuals. In a closed domain the residual sits at the
/// linear-solver tolerance, not O(1). Clamp/snap operations are accounted
/// through the "cavitation_clamp_mass" / "gas_clamp_mass" accumulator fields
/// written by the Reynolds and gas-transport solves.
///
/// Call mass_balance() after the step's solves but BEFORE store_old_time(), so
/// theta.old / dissolved_gas.old / free_gas_mass.old still hold the previous
/// timestep.
///
/// Known accounting limits (deliberately visible, they motivate WP-1/WP-2):
///  - variable-property runs recompute rho_l after the solve, so the balance
///    reflects the property lag of the segregated sweep;
///  - the dissolved-gas transport is solved in material (non-conservative)
///    form, so r_gas in advecting cases includes that discretization choice;
///  - inlet sources are evaluated for the JFO path only (GUMBEL is not
///    mass-conserving by construction, which the residual then shows).
namespace Diagnostics {

struct MassBalance {
    double liquid_mass = 0.0;          // ∫ rho h dA [kg], rho = rho_l theta
    double gas_mass = 0.0;             // ∫ (c_d rho h + m_g) dA [kg]
    double liquid_boundary_flux = 0.0; // net liquid inflow through z boundaries [kg/s]
    double gas_boundary_flux = 0.0;    // net total-gas inflow through z boundaries [kg/s]
    double inlet_mass_source = 0.0;    // net liquid injected by inlet penalty cells [kg/s]
    double clamp_mass_liquid = 0.0;    // liquid mass added by theta clamp/snap this step [kg]
    double clamp_mass_gas = 0.0;       // gas mass added by clamps/supply overwrites this step [kg]
    double liquid_residual = 0.0;      // normalized step-closure residual r_l
    double gas_residual = 0.0;         // normalized step-closure residual r_gas
    double cavitated_fraction = 0.0;   // fraction of cells with theta < 1 (cavitation extent)
};

/// Global mass balance for the step just completed. dt is the step size used by
/// the solves (steady mode passes cfg.dt; storage terms are then zero).
MassBalance mass_balance(const Fields& fields, const Mesh& mesh,
                         const SimulationConfig& cfg, double dt);

/// Per-cell liquid-continuity residual of the converged JFO state, written into
/// the optional "mass_residual" field. No-op when that field is absent, so a run
/// that does not select it keeps the original code path exactly.
///
/// Units kg/s (cell-integrated, NOT area-normalized), the same units as every
/// MassBalance rate. For one cell,
///
///   r = (Poiseuille + Couette flux out of the four faces)
///     + d(rho_l h theta)/dt + rho_l theta_old dh/dt        [transient only]
///     - m_clamp / step,
///
/// formed from the same expressions solve_elrod assembles: harmonically averaged
/// Gamma_g = Gamma_base g(theta) on interior faces, cfg.elrod_boundary_outflow on
/// active axial Dirichlet / INLET_OUTLET faces, and the scheme's Couette face
/// value. Inlet penalty cells carry only the clamp term: their flux imbalance is
/// the deliberately modelled supply source, which mass_balance() reports
/// separately as inlet_mass_source.
///
/// Interior faces telescope, so the residual field reproduces the global closure:
///   sum over ranks and cells of r == liquid_residual * liquid_mass / step,
/// with step = dt in transient runs and cfg.dt in steady runs. A cell that is not
/// an inlet, not clamped and not at the rupture front must sit at the linear
/// solver tolerance; a larger value is a defect in the assembled operator rather
/// than a physical source.
///
/// JFO only. Under GUMBEL the solved variable is pressure, not the theta balance
/// this residual describes, so the field is left at zero (validate() warns).
///
/// Same placement rule as mass_balance(): call after the step's solves but BEFORE
/// store_old_time(), so the storage term still sees the previous timestep.
void mass_residual(Fields& fields, const Mesh& mesh,
                   const SimulationConfig& cfg, double dt);

/// Rank-0 CSV appender; writes a header row when creating the file.
void append_csv(const std::filesystem::path& path, double time, int step,
                const MassBalance& balance, const Reynolds::ElrodStats& stats);

}  // namespace Diagnostics
