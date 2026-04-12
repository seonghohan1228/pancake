#pragma once

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

/// Compute film thickness fields for the journal bearing geometry.
namespace FilmThickness {

    /// Static eccentric geometry: h(theta) = c - e * cos(theta - psi)
    /// where theta is measured from some reference, psi is the attitude angle.
    void compute_static(Field& h, const Mesh& mesh, const SimulationConfig& cfg);

    /// Populate the inlet indicator field (1 if cell is inside an inlet, 0 otherwise).
    void compute_inlet_indicator(Field& indicator, const Mesh& mesh, const SimulationConfig& cfg);

}  // namespace FilmThickness
