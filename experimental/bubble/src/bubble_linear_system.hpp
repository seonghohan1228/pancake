#pragma once

#include <petsc.h>
#include <vector>

#include "field.hpp"
#include "mask.hpp"
#include "stereo_mesh.hpp"

/// Sparse linear system for 2D structured FVM on a StereoMesh.
///
/// Global row index: (offset_u + i) * n_v_global + j  (no periodic wrap).
///
/// Stencil equation (same sign convention as LinearSystem):
///   a_P * phi_P = a_E * phi_E + a_W * phi_W + a_N * phi_N + a_S * phi_S + source
///
/// Before calling solve(), call apply_mask() to inject identity (Dirichlet = 0) for
/// blanked cells and Dirichlet BC for rim cells.
class BubbleLinearSystem {
public:
    BubbleLinearSystem(const StereoMesh& mesh);
    ~BubbleLinearSystem();

    /// Zero all stencil coefficients and source terms.
    void reset();

    /// Accessors for stencil coefficients at local physical cell (i, j).
    double& a_p(int i, int j);
    double& a_e(int i, int j);
    double& a_w(int i, int j);
    double& a_n(int i, int j);
    double& a_s(int i, int j);
    double& source(int i, int j);

    /// Overwrite masked and rim cell equations.
    ///   - Masked (blanked): identity, source = 0  → phi = 0.
    ///   - Rim:              identity, source = bc_field(i,j)  → phi = bc_value.
    void apply_mask(const Mask& mask, const Field& bc_field);

    /// Assemble PETSc matrix and RHS, solve with KSP, write result into field.
    void solve(Field& result, double rtol = 1e-8, int max_iters = 1000);

private:
    const StereoMesh& mesh_ref;
    int n_u, n_v;

    std::vector<double> coeff_ap, coeff_ae, coeff_aw, coeff_an, coeff_as, coeff_src;

    Mat A;
    Vec b, x;
    KSP ksp;

    int global_row(int gu, int j) const { return gu * mesh_ref.n_v_global + j; }
    int local_index(int i, int j) const { return i * n_v + j; }
};
