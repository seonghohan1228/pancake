#pragma once

#include <vector>

#include "field.hpp"
#include "stereo_mesh.hpp"

/// Circular domain mask for the hemispherical bubble solver.
///
/// Cells with 2D stereographic radius r > R_mask are "blanked" (inactive).
/// The mask is computed from the StereoMesh cell centres at construction.
///
/// Rim cells are active cells with at least one blanked neighbour; they carry
/// Dirichlet boundary conditions (pinned h, T, no-slip velocity).
class Mask {
public:
    Mask(const StereoMesh& mesh, double R_mask);

    /// True if cell (local i, j) is inside the circular domain.
    bool is_active(int i, int j) const { return active[mesh.idx(i, j)]; }

    /// True if cell is active and borders at least one blanked cell.
    bool is_rim(int i, int j)    const { return rim[mesh.idx(i, j)]; }

    /// Set all masked cells in field to fill value.
    void blank_masked(Field& field, double fill = 0.0) const;

    /// Enforce Dirichlet BC on rim cells: set field(i,j) = value.
    void enforce_rim_bc(Field& field, double value) const;

    /// Enforce no-slip on rim: set both velocity components to zero.
    void enforce_rim_noslip(Field& u_field, Field& v_field) const;

private:
    const StereoMesh& mesh;
    int n_u, n_v;
    std::vector<bool> active;
    std::vector<bool> rim;
};
