#pragma once

#include <petsc.h>
#include <string>
#include <vector>

#include "field.hpp"
#include "mesh.hpp"

/// Krylov solver and preconditioner selection for the inner linear solves.
///
/// The defaults are chosen for rank-independent robustness rather than for the
/// lowest cost per solve. The mesh is decomposed in theta only, so every rank
/// boundary cuts the flow direction of the convection-dominated operators (the
/// steady energy equation and the Elrod film-content equation both have cell
/// Peclet numbers of order 1e4-1e6). Block Jacobi discards the coupling across
/// exactly those cuts, so its quality degrades as the rank count rises and the
/// solve can diverge at one rank count while converging at another. Additive
/// Schwarz with one cell of overlap keeps that coupling and removes the
/// rank dependence; subdomain solves use PETSc's default ILU(0).
///
/// Set `pc_type = PCBJACOBI` and `asm_overlap = 0` to recover the pre-fix
/// behaviour exactly.
/// Removing the rank dependence is necessary but not sufficient: the FBNS Schur
/// complement and the Elrod film-content operator both go near-singular whenever
/// a complementarity pivot sits on the fence, and whether a given Krylov path
/// survives such a system is not predictable from the rank count. `fallback_levels`
/// therefore lets a failed solve be retried on a genuinely different path (see
/// the ladder in linear_system.cpp) before the run is aborted. Retries cost
/// nothing on a solve that converges. Set it to 0 to abort on first failure.
struct LinearSolverOptions {
    std::string ksp_type = KSPBCGS;  ///< PETSc KSP name, e.g. KSPBCGS, KSPBCGSL, KSPGMRES.
    std::string pc_type  = PCASM;    ///< PETSc PC name, e.g. PCASM, PCBJACOBI, PCJACOBI.
    int asm_overlap      = 1;        ///< Subdomain overlap in cells; used only when pc_type == PCASM.
    int fallback_levels  = 2;        ///< Escalation rungs tried after the configured solver fails.
};

/// Translate the `linear_ksp_type` / `linear_pc_type` / `linear_pc_overlap`
/// configuration keys into the PETSc names LinearSystem hands to KSP and PC.
inline LinearSolverOptions linear_solver_options_from_config(const SimulationConfig& cfg) {
    LinearSolverOptions opts;
    switch (cfg.linear_ksp_type) {
        case KrylovMethod::BCGS:   opts.ksp_type = KSPBCGS;   break;
        case KrylovMethod::BCGSL:  opts.ksp_type = KSPBCGSL;  break;
        case KrylovMethod::GMRES:  opts.ksp_type = KSPGMRES;  break;
        case KrylovMethod::FGMRES: opts.ksp_type = KSPFGMRES; break;
    }
    switch (cfg.linear_pc_type) {
        case LinearPreconditioner::ASM:     opts.pc_type = PCASM;     break;
        case LinearPreconditioner::BJACOBI: opts.pc_type = PCBJACOBI; break;
        case LinearPreconditioner::JACOBI:  opts.pc_type = PCJACOBI;  break;
        case LinearPreconditioner::NONE:    opts.pc_type = PCNONE;    break;
    }
    opts.asm_overlap = cfg.linear_pc_overlap;
    opts.fallback_levels = cfg.linear_solver_fallbacks;
    return opts;
}

/// Sparse linear system for 2D structured FVM on a cylindrical (theta x z) mesh.
/// Each physical cell (i,j) maps to a global row index: global_theta * n_z_global + j.
/// Coefficients aP, aE, aW, aN, aS represent the standard FVM stencil:
///   aP * phi_P = aE * phi_E + aW * phi_W + aN * phi_N + aS * phi_S + source
/// Assembled into a PETSc MPIAIJ matrix and solved in parallel with KSP.
class LinearSystem {
public:
    /// `opts` selects the Krylov method and preconditioner; the default is the
    /// rank-robust ASM(1)/ILU(0) + BiCGStab combination documented above.
    explicit LinearSystem(const Mesh& mesh, const LinearSolverOptions& opts = {});
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
    /// A diverged or non-finite solve is retried on the escalation ladder (see
    /// LinearSolverOptions::fallback_levels) and only aborts the run once every
    /// rung has failed; `label` names the equation in that message (e.g.
    /// "pressure", "temperature"). Each retry that succeeds is logged on rank 0.
    void solve(Field& result, double rtol = 1e-8, int max_iters = 1000,
               const char* label = "unnamed");

private:
    const Mesh& mesh_ref;
    int n_theta, n_z;

    // Flat arrays indexed by local_index(i,j) = i * n_z + j
    std::vector<double> coeff_ap, coeff_ae, coeff_aw, coeff_an, coeff_as, coeff_src;

    Mat A;
    Vec b, x;
    KSP ksp;

    // Escalation solvers, built on first use so a run that never fails never pays
    // for them. Entry `level` is null until fallback_solver(level) creates it.
    std::vector<KSP> fallback_ksps;
    int n_fallback_levels = 0;

    int global_row(int global_theta_idx, int j) const;
    int local_index(int i, int j) const;

    /// KSP for escalation rung `level`, created on demand.
    KSP fallback_solver(int level);
};
