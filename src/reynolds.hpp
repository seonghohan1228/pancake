#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Steady-state Reynolds equation solver for journal bearing.
///
/// Solves: ∂/∂θ[ρh³/(12μR²) ∂p/∂θ] + ∂/∂z[ρh³/(12μ) ∂p/∂z] = ω/2 · ∂(ρh)/∂θ
/// Boundary conditions: Dirichlet p = p_cav at z = 0 and z = L; periodic in θ.
///
/// Uses lagged density ρ (from fields["rho"]) as the diffusion coefficient.
/// After solving, applies the Gumbel cavitation condition (p ≥ p_cav) and
/// updates fields["rho"] via the barotropic EOS.
///
/// Prerequisite: ghost cells of "h" and "rho" must be synced via Communicator.
/// Required fields: "pressure", "h", "rho".
namespace Reynolds {

    void solve(Fields& fields, LinearSystem& sys,
               const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace Reynolds
