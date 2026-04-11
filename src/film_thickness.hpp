#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

/// Compute film thickness fields for the journal bearing geometry.
namespace FilmThickness {

    /// Static eccentric geometry: h(theta) = c - e * cos(theta - psi)
    /// where theta is measured from some reference, psi is the attitude angle.
    void compute_static(Field& h, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace FilmThickness
