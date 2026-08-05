#include "linear_system.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

/// Escalation ladder for a solve the configured solver could not finish.
///
/// The rungs are deliberately different paths rather than progressively tighter
/// versions of one path. The systems that defeat a given path here are not
/// uniformly hard: they are marginal, and they appear when a complementarity
/// pivot (FBNS swap, Elrod film-content switch) sits close enough to its
/// threshold that a last-bit difference flips it and leaves one near-singular
/// row behind. Empirically such a system is solved without trouble by a
/// different preconditioner, so rung 1 keeps additive Schwarz and strengthens
/// it, while rung 2 changes the decomposition itself.
struct FallbackSpec {
    const char* prefix;             ///< PETSc options prefix owned by this rung.
    const char* ksp_type;           ///< BiCGStab(ell) throughout: it survives the erratic
                                    ///  residual histories that break plain BiCGStab.
    const char* pc_type;
    int         asm_overlap;        ///< Ignored unless pc_type is PCASM.
    const char* sub_factor_levels;  ///< ILU fill level for the subdomain solves.
};

constexpr FallbackSpec kFallbackLadder[] = {
    {"pancake_ksp_fallback1_", KSPBCGSL, PCASM,     2, "1"},
    {"pancake_ksp_fallback2_", KSPBCGSL, PCBJACOBI, 0, "1"},
};
constexpr int kFallbackLadderSize =
    static_cast<int>(sizeof(kFallbackLadder) / sizeof(kFallbackLadder[0]));

}  // namespace

LinearSystem::LinearSystem(const Mesh& mesh, const LinearSolverOptions& opts)
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
    KSPSetType(ksp, opts.ksp_type.c_str());
    PC pc;
    KSPGetPC(ksp, &pc);
    // On one rank there is nothing to decompose: additive Schwarz degenerates to a
    // single block covering the whole matrix, which is block Jacobi. The two are
    // the same method there but not the same code path -- PETSc reaches the
    // subdomain ILU through MatCreateSubMatrices for ASM and through the diagonal
    // block for block Jacobi, and the resulting ILU(0) factors are not identical.
    // On this operator the ASM path is measurably the weaker of the two (a serial
    // Ferron THD force-balance iterate takes 2.8x longer and needs 25 fallback
    // solves), so use the block-Jacobi path in serial. This also keeps serial runs
    // bit-for-bit on the pre-fix Krylov path.
    const std::string pc_name =
        (mesh.size == 1 && opts.pc_type == PCASM) ? std::string(PCBJACOBI) : opts.pc_type;
    PCSetType(pc, pc_name.c_str());
    // Overlap is meaningless for the other preconditioners and PETSc errors on
    // the call, so gate it on the type actually selected.
    if (pc_name == PCASM) PCASMSetOverlap(pc, opts.asm_overlap);
    // Runtime -ksp_*/-pc_* options still win, so a diverging production run can
    // be re-run with a different Krylov path without rebuilding.
    KSPSetFromOptions(ksp);

    n_fallback_levels = std::clamp(opts.fallback_levels, 0, kFallbackLadderSize);
    fallback_ksps.assign(static_cast<size_t>(n_fallback_levels), nullptr);
}

LinearSystem::~LinearSystem() {
    for (KSP& k : fallback_ksps) {
        if (k) KSPDestroy(&k);
    }
    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
    KSPDestroy(&ksp);
}

KSP LinearSystem::fallback_solver(int level) {
    KSP& slot = fallback_ksps[static_cast<size_t>(level)];
    if (slot) return slot;

    const FallbackSpec& spec = kFallbackLadder[level];
    KSPCreate(PETSC_COMM_WORLD, &slot);
    KSPSetOptionsPrefix(slot, spec.prefix);
    KSPSetOperators(slot, A, A);
    KSPSetType(slot, spec.ksp_type);
    PC pc;
    KSPGetPC(slot, &pc);
    // Same serial degeneracy as in the constructor: rung 1 differs from rung 2
    // only by the decomposition, so on one rank it collapses onto rung 2.
    const std::string pc_name = (mesh_ref.size > 1) ? std::string(spec.pc_type)
                                                    : std::string(PCBJACOBI);
    PCSetType(pc, pc_name.c_str());
    if (pc_name == PCASM) PCASMSetOverlap(pc, spec.asm_overlap);
    // The subdomain KSPs do not exist until PCSetUp runs inside the first solve,
    // so the subdomain fill level has to be requested through the options
    // database under this rung's own prefix rather than set on an object.
    const std::string key = std::string("-") + spec.prefix + "sub_pc_factor_levels";
    PetscOptionsSetValue(nullptr, key.c_str(), spec.sub_factor_levels);
    KSPSetFromOptions(slot);
    return slot;
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

void LinearSystem::solve(Field& result, double rtol, int max_iters, const char* label) {
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

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // A diverged KSP or a NaN must not be mistaken for a converged solve: without
    // this check the bad values propagate into the fields and the run still exits
    // successfully.
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    PetscInt n_iters = 0;

    auto attempt = [&](KSP solver) {
        KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_iters);
        KSPSetOperators(solver, A, A);
        VecZeroEntries(x);
        KSPSolve(solver, b, x);
        KSPGetConvergedReason(solver, &reason);
        KSPGetIterationNumber(solver, &n_iters);

        int local_nonfinite = 0;
        const PetscScalar* arr;
        VecGetArrayRead(x, &arr);
        for (int k = 0; k < n_theta * n_z; ++k) {
            if (!std::isfinite(arr[k])) { local_nonfinite = 1; break; }
        }
        VecRestoreArrayRead(x, &arr);
        int any_nonfinite = 0;
        MPI_Allreduce(&local_nonfinite, &any_nonfinite, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        return reason >= 0 && any_nonfinite == 0;
    };

    auto describe = [&](std::ostream& out) {
        if (reason < 0)
            out << "KSP diverged (" << KSPConvergedReasons[reason]
                << ", reason " << static_cast<int>(reason) << ")";
        else
            out << "solution contains non-finite values";
    };

    bool ok = attempt(ksp);
    for (int level = 0; !ok && level < n_fallback_levels; ++level) {
        if (rank == 0) {
            std::cerr << "\033[33m[WARN] linear solve '" << label << "' failed after "
                      << static_cast<long long>(n_iters) << " iterations: ";
            describe(std::cerr);
            std::cerr << "; retrying on fallback solver " << (level + 1) << "/"
                      << n_fallback_levels << "\033[0m" << std::endl;
        }
        ok = attempt(fallback_solver(level));
        if (ok && rank == 0) {
            std::cerr << "\033[33m[WARN] linear solve '" << label
                      << "' recovered on fallback solver " << (level + 1) << " after "
                      << static_cast<long long>(n_iters) << " iterations\033[0m" << std::endl;
        }
    }

    if (!ok) {
        if (rank == 0) {
            std::cerr << "\033[31m[FATAL] linear solve '" << label << "' failed after "
                      << static_cast<long long>(n_iters) << " iterations: ";
            describe(std::cerr);
            if (n_fallback_levels > 0)
                std::cerr << "; all " << n_fallback_levels << " fallback solvers also failed";
            std::cerr << "\033[0m" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Copy local solution back to field
    const PetscScalar* arr;
    VecGetArrayRead(x, &arr);
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_z; ++j) {
            result(i, j) = arr[local_index(i, j)];
        }
    }
    VecRestoreArrayRead(x, &arr);
}
