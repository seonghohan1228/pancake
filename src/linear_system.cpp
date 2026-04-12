#include "linear_system.hpp"

LinearSystem::LinearSystem(const Mesh& mesh)
    : mesh_ref(mesh), n_theta(mesh.n_theta_local), n_z(mesh.n_z_local)
{
    int n_local = n_theta * n_z;
    coeff_ap.resize(n_local, 0.0);
    coeff_ae.resize(n_local, 0.0);
    coeff_aw.resize(n_local, 0.0);
    coeff_an.resize(n_local, 0.0);
    coeff_as.resize(n_local, 0.0);
    coeff_src.resize(n_local, 0.0);

    int n_global = mesh.n_theta_global * mesh.n_z_global;

    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, n_local, n_local, n_global, n_global);
    MatSetType(A, MATMPIAIJ);
    // d_nz=5: diagonal block non-zeros (center + up to 4 neighbors on same rank)
    // o_nz=2: off-diagonal block non-zeros (east/west neighbors across rank boundary)
    MatMPIAIJSetPreallocation(A, 5, nullptr, 2, nullptr);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetUp(A);

    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, n_local, n_global);
    VecSetFromOptions(b);
    VecDuplicate(b, &x);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetType(ksp, KSPBCGS);
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCBJACOBI);
    KSPSetFromOptions(ksp);
}

LinearSystem::~LinearSystem() {
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
    KSPDestroy(&ksp);
}

void LinearSystem::reset() {
    std::fill(coeff_ap.begin(),  coeff_ap.end(),  0.0);
    std::fill(coeff_ae.begin(),  coeff_ae.end(),  0.0);
    std::fill(coeff_aw.begin(),  coeff_aw.end(),  0.0);
    std::fill(coeff_an.begin(),  coeff_an.end(),  0.0);
    std::fill(coeff_as.begin(),  coeff_as.end(),  0.0);
    std::fill(coeff_src.begin(), coeff_src.end(), 0.0);
}

int LinearSystem::global_row(int global_theta_idx, int j) const {
    return global_theta_idx * mesh_ref.n_z_global + j;
}

int LinearSystem::local_index(int i, int j) const {
    return i * n_z + j;
}

double& LinearSystem::a_p(int i, int j)    { return coeff_ap [local_index(i, j)]; }
double& LinearSystem::a_e(int i, int j)    { return coeff_ae [local_index(i, j)]; }
double& LinearSystem::a_w(int i, int j)    { return coeff_aw [local_index(i, j)]; }
double& LinearSystem::a_n(int i, int j)    { return coeff_an [local_index(i, j)]; }
double& LinearSystem::a_s(int i, int j)    { return coeff_as [local_index(i, j)]; }
double& LinearSystem::source(int i, int j) { return coeff_src[local_index(i, j)]; }

void LinearSystem::solve(Field& result, double rtol, int max_iters) {
    MatZeroEntries(A);
    VecZeroEntries(b);

    const int n_theta_g = mesh_ref.n_theta_global;
    const int off       = mesh_ref.offset_theta;

    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            const int idx = local_index(i, j);
            const PetscInt row = global_row(off + i, j);

            // Diagonal: a_P
            PetscScalar ap = coeff_ap[idx];
            MatSetValue(A, row, row, ap, INSERT_VALUES);

            // East: column wraps periodically in theta
            if (coeff_ae[idx] != 0.0) {
                PetscInt col_e = global_row((off + i + 1) % n_theta_g, j);
                MatSetValue(A, row, col_e, -coeff_ae[idx], INSERT_VALUES);
            }

            // West: column wraps periodically in theta
            if (coeff_aw[idx] != 0.0) {
                PetscInt col_w = global_row((off + i - 1 + n_theta_g) % n_theta_g, j);
                MatSetValue(A, row, col_w, -coeff_aw[idx], INSERT_VALUES);
            }

            // North: no periodicity in z; skip at top boundary (a_n should be 0 there)
            if (coeff_an[idx] != 0.0 && j + 1 < n_z) {
                PetscInt col_n = global_row(off + i, j + 1);
                MatSetValue(A, row, col_n, -coeff_an[idx], INSERT_VALUES);
            }

            // South: no periodicity in z; skip at bottom boundary
            if (coeff_as[idx] != 0.0 && j - 1 >= 0) {
                PetscInt col_s = global_row(off + i, j - 1);
                MatSetValue(A, row, col_s, -coeff_as[idx], INSERT_VALUES);
            }

            VecSetValue(b, row, coeff_src[idx], INSERT_VALUES);
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    KSPSetTolerances(ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters);
    KSPSetOperators(ksp, A, A);
    VecZeroEntries(x);
    KSPSolve(ksp, b, x);

    // Copy local solution back to field
    const PetscScalar* arr;
    VecGetArrayRead(x, &arr);
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_z; ++j)
            result(i, j) = arr[local_index(i, j)];
    VecRestoreArrayRead(x, &arr);
}
