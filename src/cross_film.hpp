#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

/// Cross-film (through-thickness) temperature profile, conjugate wall coupling,
/// and the Dowson viscosity integrals of the generalized Reynolds equation.
///
/// The whole module is gated on `cfg.cross_film_thermal_model`. With the default
/// `NONE` nothing here is called and the solver keeps its depth-averaged,
/// uniform-across-the-gap behaviour exactly.
///
/// With `POLYNOMIAL`, `update()` must be called once per outer (Picard) iteration,
/// before the Reynolds solve, with the temperature/pressure/property ghosts
/// current. It writes the fields listed in `register_fields()`; the caller must
/// exchange ghosts for `cross_film_f2` and `cross_film_couette_h` afterwards
/// because the Reynolds assembly reads them on the theta ghost layers.
namespace CrossFilm {

/// Names of every field this module writes, so the caller can allocate them.
/// Allocation is unconditional and harmless; only `update()` is gated.
void register_fields(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

/// Result of one conjugate-wall solve, for logging and diagnostics.
struct WallState {
    double t_shaft = 0.0;      ///< Circumferentially uniform shaft temperature [K]
    double t_bush_mean = 0.0;  ///< Area-mean bush surface temperature [K]
    double q_shaft_net = 0.0;  ///< Net film -> shaft heat flow [W] (0 when adiabatic)
    double q_bush_net = 0.0;   ///< Net film -> bush heat flow [W]
    double core_rise_max = 0.0;///< Largest midplane-minus-mean temperature [K]
};

/// Rebuild the cross-film profile, the conjugate wall temperatures, the wall
/// heat-flux closure, and the F-integrals. No-op unless the model is POLYNOMIAL.
WallState update(Fields& fields, const Mesh& mesh, const SimulationConfig& cfg);

/// True when the generalized-Reynolds coefficients are available and should be used.
inline bool active(const SimulationConfig& cfg) {
    return cfg.cross_film_thermal_model != CrossFilmModel::NONE;
}

}  // namespace CrossFilm
