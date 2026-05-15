#include "bubble_linear_system.hpp"

BubbleLinearSystem::BubbleLinearSystem(const StereoMesh& mesh)
    : mesh_ref(mesh), n_u(mesh.n_u_local), n_v(mesh.n_v_global)
{
    const int n_local  = n_u * n_v;
    const int n_global = mesh.n_u_global * mesh.n_v_global;

    coeff_ap.resize(n_local, 0.0);
    coeff_ae.resize(n_local, 0.0);
    coeff_aw.resize(n_local, 0.0);
    coeff_an.resize(n_local, 0.0);
    coeff_as.resize(n_local, 0.0);
    coeff_src.resize(n_local, 0.0);

    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, n_local, n_local, n_global, n_global);
    MatSetType(A, MATMPIAIJ);
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

BubbleLinearSystem::~BubbleLinearSystem() {
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
    KSPDestroy(&ksp);
}

void BubbleLinearSystem::reset() {
    std::fill(coeff_ap.begin(),  coeff_ap.end(),  0.0);
    std::fill(coeff_ae.begin(),  coeff_ae.end(),  0.0);
    std::fill(coeff_aw.begin(),  coeff_aw.end(),  0.0);
    std::fill(coeff_an.begin(),  coeff_an.end(),  0.0);
    std::fill(coeff_as.begin(),  coeff_as.end(),  0.0);
    std::fill(coeff_src.begin(), coeff_src.end(), 0.0);
}

double& BubbleLinearSystem::a_p(int i, int j)    { return coeff_ap [local_index(i, j)]; }
double& BubbleLinearSystem::a_e(int i, int j)    { return coeff_ae [local_index(i, j)]; }
double& BubbleLinearSystem::a_w(int i, int j)    { return coeff_aw [local_index(i, j)]; }
double& BubbleLinearSystem::a_n(int i, int j)    { return coeff_an [local_index(i, j)]; }
double& BubbleLinearSystem::a_s(int i, int j)    { return coeff_as [local_index(i, j)]; }
double& BubbleLinearSystem::source(int i, int j) { return coeff_src[local_index(i, j)]; }

void BubbleLinearSystem::apply_mask(const Mask& mask, const Field& bc_field) {
    for (int i = 0; i < n_u; ++i) {
        for (int j = 0; j < n_v; ++j) {
            if (mask.is_active(i, j) && !mask.is_rim(i, j)) continue;
            const int idx = local_index(i, j);
            coeff_ae [idx] = 0.0;
            coeff_aw [idx] = 0.0;
            coeff_an [idx] = 0.0;
            coeff_as [idx] = 0.0;
            coeff_ap [idx] = 1.0;
            coeff_src[idx] = mask.is_rim(i, j) ? bc_field(i, j) : 0.0;
        }
    }
}

void BubbleLinearSystem::solve(Field& result, double rtol, int max_iters) {
    MatZeroEntries(A);
    VecZeroEntries(b);

    for (int i = 0; i < n_u; ++i) {
        for (int j = 0; j < n_v; ++j) {
            const int     idx = local_index(i, j);
            const int     gu  = mesh_ref.global_u(i);
            const PetscInt row = global_row(gu, j);

            MatSetValue(A, row, row, coeff_ap[idx], INSERT_VALUES);

            // East (no periodic wrap)
            if (gu + 1 < mesh_ref.n_u_global && coeff_ae[idx] != 0.0)
                MatSetValue(A, row, global_row(gu + 1, j), -coeff_ae[idx], INSERT_VALUES);

            // West
            if (gu - 1 >= 0 && coeff_aw[idx] != 0.0)
                MatSetValue(A, row, global_row(gu - 1, j), -coeff_aw[idx], INSERT_VALUES);

            // North
            if (j + 1 < n_v && coeff_an[idx] != 0.0)
                MatSetValue(A, row, global_row(gu, j + 1), -coeff_an[idx], INSERT_VALUES);

            // South
            if (j - 1 >= 0 && coeff_as[idx] != 0.0)
                MatSetValue(A, row, global_row(gu, j - 1), -coeff_as[idx], INSERT_VALUES);

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

    const PetscScalar* arr;
    VecGetArrayRead(x, &arr);
    for (int i = 0; i < n_u; ++i)
        for (int j = 0; j < n_v; ++j)
            result(i, j) = arr[local_index(i, j)];
    VecRestoreArrayRead(x, &arr);
}
