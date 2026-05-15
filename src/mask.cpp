#include "mask.hpp"

Mask::Mask(const StereoMesh& m, double R_mask) : mesh(m), n_u(m.n_u_local), n_v(m.n_v_global) {
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
                return active[mesh.idx(ni, nj)];
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
