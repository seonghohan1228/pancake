#pragma once

#include "bubble_linear_system.hpp"
#include "field.hpp"
#include "bubble_config.hpp"  // ConvectionScheme
#include "mask.hpp"
#include "stereo_mesh.hpp"

/// Metric-aware finite-volume operators for the hemispherical surface mesh.
///
/// Sign convention: operators contribute to the LHS of the discrete PDE.
/// Assembled equation: a_P*phi_P = a_E*phi_E + a_W*phi_W + a_N*phi_N + a_S*phi_S + source.
///
/// Each operator skips cells whose neighbour is masked (outside the active domain)
/// so that no unphysical diffusion/convection crosses the domain boundary.
/// Rim cells are NOT skipped; BubbleLinearSystem::apply_mask() enforces their BCs.
///
/// Ghost cells in phi and flux fields must be synced via Communicator before calling.
namespace SurfaceFVM {

    /// d(phi)/dt: adds A_cell/dt to a_P and A_cell/dt * phi.old(i,j) to source.
    void ddt(BubbleLinearSystem& sys, const Field& phi, double dt,
             const StereoMesh& mesh, const Mask& mask);

    /// d(weight * phi)/dt: adds weight(i,j)*A_cell/dt to a_P.
    /// Use for terms of the form d(rho*h*phi)/dt where rho*h is the weight.
    void ddt_weighted(BubbleLinearSystem& sys, const Field& phi, const Field& weight,
                      double dt, const StereoMesh& mesh, const Mask& mask);

    /// div(gamma * grad(phi)):
    ///   face coefficient = gamma_f * L_face / d_PN  (harmonic-averaged gamma).
    ///   Face geometry taken from precomputed StereoMesh arrays; boundary faces use
    ///   analytic computation via StereoMesh::east_face_len / d_u_east.
    void laplacian(BubbleLinearSystem& sys, const Field& gamma,
                   const StereoMesh& mesh, const Mask& mask);

    /// div(flux * phi) where flux_{u,v} are volumetric face fluxes [m^2/s]:
    ///   flux_u(i,j) = outward east-face flux of cell (i,j)   (positive = eastward).
    ///   flux_v(i,j) = outward north-face flux of cell (i,j)  (positive = northward).
    /// UPWIND:     fully implicit.
    /// TVD_*:      upwind implicit + deferred TVD correction added to source (uses phi.old).
    void divergence(BubbleLinearSystem& sys,
                    const Field& flux_u, const Field& flux_v, const Field& phi,
                    ConvectionScheme scheme,
                    const StereoMesh& mesh, const Mask& mask);

    /// Explicit source: adds source_field(i,j) * A_cell(i,j) to source(i,j).
    void add_source(BubbleLinearSystem& sys, const Field& source_field,
                    const StereoMesh& mesh, const Mask& mask);

    /// Cell-centred gradient via central difference of the 3D chord distances.
    ///   grad_u(i,j) = (phi(i+1,j) - phi(i-1,j)) / (d_u(i,j) + d_u(i-1,j))
    /// Masked cells get zero gradient.
    void gradient(const Field& phi, Field& grad_u, Field& grad_v,
                  const StereoMesh& mesh, const Mask& mask);

}  // namespace SurfaceFVM
