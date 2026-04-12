#pragma once

#include "field.hpp"
#include "linear_system.hpp"
#include "mesh.hpp"

enum class ConvectionScheme { UPWIND, TVD_VANLEER, TVD_MINMOD };

/// Finite-volume operators for the cylindrical (theta x z) mesh.
/// Each operator adds its contribution to a LinearSystem's stencil coefficients
/// and source terms. Operators are additive: call multiple to build an equation.
///
/// Sign convention: operators represent terms that appear on the LHS of the PDE
/// with a POSITIVE sign. The assembled equation is:
///   aP * phi_P = aE * phi_E + aW * phi_W + aN * phi_N + aS * phi_S + source
///
/// Ghost cells for phi and flux fields must be synced via Communicator before calling.
/// Z-boundary faces (j=0 south, j=n_z-1 north) are omitted; callers handle those BCs.
namespace FVM {

    /// Time derivative: d(phi)/dt.
    /// Adds vol/dt to aP and vol/dt * phi.old(i,j) to source.
    void ddt(LinearSystem& sys, const Field& phi, double dt, const Mesh& mesh);

    /// Weighted time derivative: d(weight * phi)/dt with static weight field.
    /// Adds weight(i,j)*vol/dt to aP and weight(i,j)*vol/dt * phi.old(i,j) to source.
    /// Use for transient terms of the form ∂(h·θ)/∂t where h is the weight.
    void ddt_weighted(LinearSystem& sys, const Field& phi, const Field& weight,
                      double dt, const Mesh& mesh);

    /// Laplacian: div(gamma * grad(phi)).
    /// gamma is cell-centered; face values use harmonic averaging.
    /// Adds to aP, aE, aW, aN, aS.
    void laplacian(LinearSystem& sys, const Field& gamma, const Mesh& mesh);

    /// Divergence of convective flux: div(F * phi) where F is a known face flux.
    /// flux_theta(i,j) = outward volumetric flux through the east face of cell (i,j).
    /// flux_z(i,j)     = outward volumetric flux through the north face of cell (i,j).
    /// TVD schemes use phi.old for deferred correction (requires outer iteration).
    void divergence(LinearSystem& sys,
                    const Field& flux_theta, const Field& flux_z, const Field& phi,
                    ConvectionScheme scheme, const Mesh& mesh);

    /// Explicit source: adds source_field(i,j) * cell_volume to source(i,j).
    void add_source(LinearSystem& sys, const Field& source_field, const Mesh& mesh);

}  // namespace FVM
