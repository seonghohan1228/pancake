#pragma once

#include <petsc.h>
#include <vector>

#include "field.hpp"
#include "mesh.hpp"

/// Sparse linear system for 2D structured FVM on a cylindrical (theta x z) mesh.
/// Each physical cell (i,j) maps to a global row index: global_theta * n_z_global + j.
/// Coefficients aP, aE, aW, aN, aS represent the standard FVM stencil:
///   aP * phi_P = aE * phi_E + aW * phi_W + aN * phi_N + aS * phi_S + source
/// Assembled into a PETSc MPIAIJ matrix and solved in parallel with KSP.
class LinearSystem {
public:
    LinearSystem(const Mesh& mesh);
    ~LinearSystem();

    /// Zero all stencil coefficients and source terms.
    void reset();

    /// Accessors for stencil coefficients at local physical cell (i,j).
    double& a_p(int i, int j);
    double& a_e(int i, int j);
    double& a_w(int i, int j);
    double& a_n(int i, int j);
    double& a_s(int i, int j);
    double& source(int i, int j);

    /// Assemble the PETSc matrix and RHS, then solve. Writes result into field.
    void solve(Field& result, double rtol = 1e-8, int max_iters = 1000);

private:
    const Mesh& mesh_ref;
    int n_theta, n_z;

    // Flat arrays indexed by local_index(i,j) = i * n_z + j
    std::vector<double> coeff_ap, coeff_ae, coeff_aw, coeff_an, coeff_as, coeff_src;

    Mat A;
    Vec b, x;
    KSP ksp;

    int global_row(int global_theta_idx, int j) const;
    int local_index(int i, int j) const;
};
