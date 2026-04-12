#pragma once

#include "config.hpp"
#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

/// Reynolds equation solvers for journal bearing.
///
/// Two formulations are available, selected via cfg.cavitation_model:
///
/// GUMBEL — pressure primary variable, post-solve pressure clamp.
///   Solves: ∂/∂θ[ρh³/(12μR²) ∂p/∂θ] + ∂/∂z[ρh³/(12μ) ∂p/∂z] = ω/2 · ∂(ρh)/∂θ
///   Required fields: "pressure", "h", "rho".
///
/// ELROD_ADAMS — universal variable θ (film content) primary variable.
///   Solves: div(Γ·∇θ) = div(F_couette·θ)  where Γ = ρ₀βh³/(12μ), F = ωRh/2·Δz
///   Full-film formulation (Phase A.1): switch function g(θ)=1 everywhere.
///   After solve: recovers pressure and density from EOS.
///   Required fields: "theta", "pressure", "h", "rho".
///
/// Prerequisite: ghost cells of "h", "rho" (GUMBEL) or "h", "theta" (ELROD_ADAMS)
/// must be synced via Communicator before calling.
namespace Reynolds {

    /// Gumbel (pressure-primary) solver. See Phase 3 notes in PHYSICS.md.
    void solve(Fields& fields, LinearSystem& sys,
               const Mesh& mesh, const SimulationConfig& cfg);

    /// Elrod-Adams (θ-primary) solver. Phase A.1: full-film only (g=1 everywhere).
    /// Sets fields["theta"], fields["pressure"], and fields["rho"].
    void solve_elrod(Fields& fields, LinearSystem& sys,
                     const Mesh& mesh, const SimulationConfig& cfg);

    /// Calculate velocities on a staggered grid.
    /// u_theta is on FACE_THETA, u_z is on FACE_Z.
    void calculate_velocities(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace Reynolds
