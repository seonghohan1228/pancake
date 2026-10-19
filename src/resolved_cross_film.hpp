#pragma once
#include "cross_film.hpp"
#include "energy.hpp"

/// Steady conservative energy in equal through-film finite-volume layers,
/// coupled to an isothermal shaft and local, r-theta or r-theta-z housing conduction.
namespace ResolvedCrossFilm {
void register_fields(Fields&, const Mesh&, const SimulationConfig&);
/// Blend steady COUPLED layer inverse viscosities before integrating all flow coefficients.
CrossFilm::WallState update(Fields&, const Mesh&, const SimulationConfig&, double relaxation = 1.0);
CrossFilm::WallState diagnostics(const Fields&, const Mesh&, const SimulationConfig&);
void solve(Fields&, const Mesh&, const SimulationConfig&);
Energy::Balance balance(const Fields&, const Mesh&, const SimulationConfig&);
}
