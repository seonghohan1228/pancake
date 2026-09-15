#pragma once

#include <vector>

#include "config.hpp"

/// Steady two-dimensional (r, theta) heat conduction in the bush annulus.
///
/// This is the optional replacement for the LOCAL 1-D radial bush closure of
/// `CrossFilm` (PHYSICS.md 13.1.2). The local closure sends every watt that
/// enters the bore straight out through the wall underneath it, so the bush
/// surface temperature can only peak where the film is hottest. Ferron (1983)
/// instead solves conduction in the bush itself, which lets heat picked up in
/// the loaded arc leave through a different angular position, and it is that
/// circumferential redistribution that turns the measured bush temperature down
/// before the groove.
///
/// The module is deliberately free of MPI and of `Fields`: it takes the
/// z-averaged film-side flux law on the global theta grid and returns the bore
/// temperature on the same grid, so it can be unit-tested on its own and every
/// rank can run it redundantly on identical input and obtain identical output.
///
/// Governing equation of this reduced, z-collapsed bush model,
///
///     (1/r) d/dr (r k_b dT/dr) + (1/r^2) d/dtheta (k_b dT/dtheta) = 0,
///     R <= r <= R2, theta periodic.
///
/// Ferron (1983) instead retains the axial term and applies free convection on
/// the two axial end faces. This module cannot represent that local r-theta-z
/// boundary-value problem; its input and output are explicitly z-averaged.
///
/// Inner boundary: flux continuity with the film. The cross-film profile makes
/// the film-side flux linear in the bore temperature, so the coupling is a
/// Robin condition and the whole bush problem stays linear (see `Coupling`).
/// Outer boundary: convection to `bush_outer_temperature` when positive
/// (otherwise `thermal_ambient_temperature`) with `bush_outer_heat_transfer`,
/// or a fixed temperature, per
/// `cfg.bush_outer_boundary`.
namespace BushConduction {

/// Per-theta, z-averaged film-side flux law at the bore, positive for heat
/// LEAVING the film and entering the bush:
///
///     q(theta) = a0(theta) + a1(theta) * T_shaft - slope(theta) * T_bush(theta)
///
/// The split into a shaft-independent part and a shaft slope is what keeps the
/// conjugate system linear: the bush is solved once for each part and the shaft
/// balance then closes in one step, exactly as in the LOCAL closure.
///
/// All arrays are global in theta (length `n_theta`), and `slope` must be
/// non-negative (it is (k/h)(1 + b_slope/2) > 0 for the polynomial profile).
struct Coupling {
    std::vector<double> a0;     ///< [W/m^2]
    std::vector<double> a1;     ///< [W/(m^2 K)]
    std::vector<double> slope;  ///< [W/(m^2 K)]
};

/// Bore temperature as an affine function of the shaft temperature.
struct Solution {
    std::vector<double> alpha;  ///< T_bush(theta) at T_shaft = 0 [K]
    std::vector<double> beta;   ///< d T_bush / d T_shaft [-]
    double residual = 0.0;      ///< relative residual of the factorisation (round-off)
    bool converged = true;      ///< false only if the direct solve broke down
};

/// Solve the annulus. `r_inner` is the bore radius (the mesh radius R); the
/// outer radius, conductivity, outer boundary condition, ambient temperature and
/// radial cell count come from `cfg`. `n_theta` must equal the length of the
/// coupling arrays.
///
/// Discretisation: finite volume, uniform in theta, uniform in r, with the exact
/// logarithmic conductance k_b d(theta)/ln(r_2/r_1) on every radial face
/// (including the two half-cells at the boundaries). A purely radial problem is
/// therefore reproduced exactly at any radial resolution.
///
/// Linear solve: DIRECT, not iterative. The system is block tridiagonal in theta
/// with n_r x n_r blocks and one periodic wrap, so it is solved by block Thomas
/// on the first n_theta-1 blocks carried with an extra n_r right-hand-side
/// columns for the wrap, followed by one dense n_r x n_r Schur solve for the last
/// block. Cost is O(n_theta * n_r^3), about 10^7 flops at 720 x 16 - negligible
/// beside one film pressure solve, and with no convergence tolerance to tune. An
/// iterative solver was tried first and rejected: the bush grid is close to
/// isotropic (radial and circumferential face spacings are within a factor of a
/// few), which is exactly the regime where line relaxation stalls.
Solution solve(const Coupling& coupling, int n_theta, double r_inner,
               const SimulationConfig& cfg);

/// Relative residual above which the direct solve is reported as broken down.
constexpr double kResidualTolerance = 1.0e-8;

}  // namespace BushConduction
