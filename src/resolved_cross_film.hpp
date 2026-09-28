#pragma once
#include "cross_film.hpp"
#include "energy.hpp"

/// Steady conservative energy in equal through-film finite-volume layers,
/// coupled to an isothermal shaft and local, r-theta or r-theta-z housing conduction.
namespace ResolvedCrossFilm {
void register_fields(Fields&, const Mesh&, const SimulationConfig&);
CrossFilm::WallState update(Fields&, const Mesh&, const SimulationConfig&);
CrossFilm::WallState diagnostics(const Fields&, const Mesh&, const SimulationConfig&);
void solve(Fields&, const Mesh&, const SimulationConfig&);
Energy::Balance balance(const Fields&, const Mesh&, const SimulationConfig&);
}
