#include "mask.hpp"

Mask::Mask(const StereoMesh& m, double R_mask)
    : mesh(m), n_u(m.n_u_local), n_v(m.n_v_global), radius_squared(R_mask * R_mask) {
    active.assign(n_u * n_v, false);
    rim.assign(n_u * n_v, false);

    // Mark active cells: those whose 2D centre radius is within R_mask
    for (int i = 0; i < n_u; ++i) {
        for (int j = 0; j < n_v; ++j) {
            double uc = mesh.cell_u(i);
            double vc = mesh.cell_v(j);
            double r2 = uc*uc + vc*vc;
            active[mesh.idx(i, j)] = (r2 <= R_mask * R_mask);
        }
    }

    // Mark rim cells: active cells with at least one blanked (or out-of-domain) neighbour.
    // Neighbours are the 4 cardinal directions; cells at the domain edge are treated as blanked.
    for (int i = 0; i < n_u; ++i) {
        for (int j = 0; j < n_v; ++j) {
            if (!active[mesh.idx(i, j)]) continue;

            auto neighbour_active = [&](int ni, int nj) -> bool {
                // Domain edges count as blanked
                int global_ni = mesh.global_u(ni);
                if (global_ni < 0 || global_ni >= mesh.n_u_global) return false;
                if (nj < 0 || nj >= n_v) return false;
                // A neighbour across an MPI boundary is not in the local array.
                const double u = mesh.cell_u(ni);
                const double v = mesh.cell_v(nj);
                return u*u + v*v <= R_mask * R_mask;
            };

            bool blanked_nbr =
                !neighbour_active(i - 1, j) ||
                !neighbour_active(i + 1, j) ||
                !neighbour_active(i, j - 1) ||
                !neighbour_active(i, j + 1);

            rim[mesh.idx(i, j)] = blanked_nbr;
        }
    }
}

bool Mask::is_active(int i, int j) const {
    if (mesh.global_u(i) < 0 || mesh.global_u(i) >= mesh.n_u_global ||
        j < 0 || j >= n_v) return false;
    if (i >= 0 && i < n_u) return active[mesh.idx(i, j)];
    const double u = mesh.cell_u(i), v = mesh.cell_v(j);
    return u*u + v*v <= radius_squared;
}

bool Mask::is_rim(int i, int j) const {
    if (!is_active(i, j)) return false;
    if (i >= 0 && i < n_u) return rim[mesh.idx(i, j)];
    return !is_active(i-1, j) || !is_active(i+1, j) ||
           !is_active(i, j-1) || !is_active(i, j+1);
}

void Mask::blank_masked(Field& field, double fill) const {
    for (int i = 0; i < n_u; ++i)
        for (int j = 0; j < n_v; ++j)
            if (!active[mesh.idx(i, j)])
                field(i, j) = fill;
}

void Mask::enforce_rim_bc(Field& field, double value) const {
    for (int i = 0; i < n_u; ++i)
        for (int j = 0; j < n_v; ++j)
            if (rim[mesh.idx(i, j)])
                field(i, j) = value;
}

void Mask::enforce_rim_noslip(Field& u_field, Field& v_field) const {
    for (int i = 0; i < n_u; ++i)
        for (int j = 0; j < n_v; ++j)
            if (rim[mesh.idx(i, j)]) {
                u_field(i, j) = 0.0;
                v_field(i, j) = 0.0;
            }
}
