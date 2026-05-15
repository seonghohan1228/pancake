#pragma once

#include "bubble_config.hpp"
#include "bubble_linear_system.hpp"
#include "field.hpp"
#include "bubble_config.hpp"  // ConvectionScheme
#include "mask.hpp"
#include "stereo_mesh.hpp"

/// Thickness transport and gravity-driven velocity for the bubble solver.
namespace ThicknessTransport {

    /// Compute cell-centred velocity driven by thin-film Stokes drainage:
    ///   u_s = (rho_l * h^2 / (3 * mu)) * g_s
    /// Rim and masked cells are zeroed. Results in u_u and u_v (overwritten).
    void compute_gravity_velocity(Field& u_u, Field& u_v,
                                  const Field& h,
                                  const StereoMesh& mesh, const Mask& mask,
                                  const BubbleConfig& cfg);

    /// Interpolate cell-centred velocities to face fluxes:
    ///   flux_u(i,j) = 0.5*(u_u(i,j) + u_u(i+1,j)) * face_len_e(i,j)
    ///   flux_v(i,j) = 0.5*(u_v(i,j) + u_v(i,j+1)) * face_len_n(i,j)
    /// Ghost layers of u_u and u_v must be exchanged before calling.
    void compute_face_fluxes(Field& flux_u, Field& flux_v,
                              const Field& u_u, const Field& u_v,
                              const StereoMesh& mesh, const Mask& mask);

    /// Assemble and solve the thickness transport equation:
    ///   d(h)/dt + div(h * u_s) = evap_rate (optional, zero if nullptr)
    /// Applies rim Dirichlet BC (h = h_rim) and masks blanked cells.
    /// Clamps result to h >= h_min after solve.
    void solve(Field& h,
               const Field& flux_u, const Field& flux_v,
               const Field* evap_rate,            // nullptr to disable evaporation source
               BubbleLinearSystem& sys,
               const StereoMesh& mesh, const Mask& mask,
               const BubbleConfig& cfg, double dt);

}  // namespace ThicknessTransport
