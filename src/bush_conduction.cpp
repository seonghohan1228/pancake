#include "bush_conduction.hpp"

#include <algorithm>
#include <cmath>

namespace BushConduction {
namespace {

// -------------------------------------------------------------------------
// Small dense kernels. The blocks are n_r x n_r with n_r = 16 by default, so
// these are deliberately plain: row-major, partial pivoting, no blocking.
// -------------------------------------------------------------------------

/// LU with partial pivoting, in place. Returns false on an exactly singular
/// pivot, which for this operator means the problem was posed with neither a
/// bore Robin term nor an outer sink.
bool lu_factor(std::vector<double>& a, int n, std::vector<int>& piv) {
    piv.resize(n);
    for (int k = 0; k < n; ++k) {
        int best = k;
        double best_abs = std::abs(a[static_cast<size_t>(k) * n + k]);
        for (int r = k + 1; r < n; ++r) {
            const double v = std::abs(a[static_cast<size_t>(r) * n + k]);
            if (v > best_abs) { best_abs = v; best = r; }
        }
        piv[k] = best;
        if (!(best_abs > 0.0)) return false;
        if (best != k) {
            for (int c = 0; c < n; ++c) {
                std::swap(a[static_cast<size_t>(k) * n + c], a[static_cast<size_t>(best) * n + c]);
            }
        }
        const double pivot = a[static_cast<size_t>(k) * n + k];
        for (int r = k + 1; r < n; ++r) {
            const double f = a[static_cast<size_t>(r) * n + k] / pivot;
            a[static_cast<size_t>(r) * n + k] = f;
            if (f == 0.0) continue;
            for (int c = k + 1; c < n; ++c) {
                a[static_cast<size_t>(r) * n + c] -= f * a[static_cast<size_t>(k) * n + c];
            }
        }
    }
    return true;
}

/// Solve L U X = B in place for `nrhs` right-hand sides stored row-major
/// (B[row * nrhs + col]).
void lu_solve(const std::vector<double>& lu, int n, const std::vector<int>& piv,
              std::vector<double>& b, int nrhs) {
    for (int k = 0; k < n; ++k) {
        if (piv[k] != k) {
            for (int c = 0; c < nrhs; ++c) {
                std::swap(b[static_cast<size_t>(k) * nrhs + c],
                          b[static_cast<size_t>(piv[k]) * nrhs + c]);
            }
        }
    }
    for (int r = 1; r < n; ++r) {
        for (int k = 0; k < r; ++k) {
            const double f = lu[static_cast<size_t>(r) * n + k];
            if (f == 0.0) continue;
            for (int c = 0; c < nrhs; ++c) {
                b[static_cast<size_t>(r) * nrhs + c] -= f * b[static_cast<size_t>(k) * nrhs + c];
            }
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        for (int k = r + 1; k < n; ++k) {
            const double f = lu[static_cast<size_t>(r) * n + k];
            if (f == 0.0) continue;
            for (int c = 0; c < nrhs; ++c) {
                b[static_cast<size_t>(r) * nrhs + c] -= f * b[static_cast<size_t>(k) * nrhs + c];
            }
        }
        const double d = lu[static_cast<size_t>(r) * n + r];
        for (int c = 0; c < nrhs; ++c) b[static_cast<size_t>(r) * nrhs + c] /= d;
    }
}

/// Assembled coefficients of the annulus, shared by both right-hand sides.
/// The theta coupling is the same diagonal matrix diag(gt) on every face, which
/// is what makes the block sweep below cheap: the off-diagonal blocks never need
/// a matrix product, only a row scaling.
struct Operator {
    int n_theta = 0;
    int n_r = 0;
    std::vector<double> gr;            ///< radial face conductance, size n_r-1 [W/(m K)]
    std::vector<double> gt;            ///< theta face conductance per radial level [W/(m K)]
    std::vector<double> inner_gain;    ///< lambda*S*slope added to the m = 0 diagonal [W/(m K)]
    std::vector<double> inner_lambda;  ///< half-cell attenuation of the bore flux
    double outer_u = 0.0;              ///< outer-face conductance to ambient [W/(m K)]
    double c_half_inner = 0.0;         ///< bore-to-first-node conductance [W/(m K)]
    double s_inner = 0.0;              ///< bore face length per unit axial length [m]
    std::vector<double> diag;          ///< size n_theta*n_r

    size_t at(int i, int m) const { return static_cast<size_t>(i) * n_r + m; }

    /// Dense block B_i: radial tridiagonal plus the full diagonal.
    void block(int i, std::vector<double>& out) const {
        out.assign(static_cast<size_t>(n_r) * n_r, 0.0);
        for (int m = 0; m < n_r; ++m) {
            out[static_cast<size_t>(m) * n_r + m] = diag[at(i, m)];
            if (m + 1 < n_r) {
                out[static_cast<size_t>(m) * n_r + m + 1] = -gr[m];
                out[static_cast<size_t>(m + 1) * n_r + m] = -gr[m];
            }
        }
    }
};

/// Residual of the assembled system, relative to the source scale. A direct
/// solve should return round-off; anything larger means the factorisation was
/// not the operator the physics assembled.
double relative_residual(const Operator& op, const std::vector<double>& src,
                         const std::vector<double>& x) {
    double scale = 0.0;
    for (double v : src) scale = std::max(scale, std::abs(v));
    if (!(scale > 0.0)) return 0.0;
    double worst = 0.0;
    for (int i = 0; i < op.n_theta; ++i) {
        const int im = (i - 1 + op.n_theta) % op.n_theta;
        const int ip = (i + 1) % op.n_theta;
        for (int m = 0; m < op.n_r; ++m) {
            double lhs = op.diag[op.at(i, m)] * x[op.at(i, m)];
            lhs -= op.gt[m] * (x[op.at(im, m)] + x[op.at(ip, m)]);
            if (m > 0) lhs -= op.gr[m - 1] * x[op.at(i, m - 1)];
            if (m + 1 < op.n_r) lhs -= op.gr[m] * x[op.at(i, m + 1)];
            worst = std::max(worst, std::abs(src[op.at(i, m)] - lhs));
        }
    }
    return worst / scale;
}

/// Dense fallback for the degenerate theta resolutions (n_theta <= 2), where the
/// two theta neighbours of a cell are the same cell and the block-tridiagonal
/// structure does not exist. Never taken on a real mesh.
bool solve_dense(const Operator& op, std::vector<std::vector<double>>& rhs) {
    const int n = op.n_theta * op.n_r;
    std::vector<double> a(static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < op.n_theta; ++i) {
        const int im = (i - 1 + op.n_theta) % op.n_theta;
        const int ip = (i + 1) % op.n_theta;
        for (int m = 0; m < op.n_r; ++m) {
            const size_t row = op.at(i, m);
            a[row * n + row] += op.diag[row];
            a[row * n + op.at(im, m)] -= op.gt[m];
            a[row * n + op.at(ip, m)] -= op.gt[m];
            if (m > 0) a[row * n + op.at(i, m - 1)] -= op.gr[m - 1];
            if (m + 1 < op.n_r) a[row * n + op.at(i, m + 1)] -= op.gr[m];
        }
    }
    std::vector<int> piv;
    if (!lu_factor(a, n, piv)) return false;
    const int nrhs = static_cast<int>(rhs.size());
    std::vector<double> b(static_cast<size_t>(n) * nrhs, 0.0);
    for (int c = 0; c < nrhs; ++c) {
        for (int row = 0; row < n; ++row) b[static_cast<size_t>(row) * nrhs + c] = rhs[c][row];
    }
    lu_solve(a, n, piv, b, nrhs);
    for (int c = 0; c < nrhs; ++c) {
        for (int row = 0; row < n; ++row) rhs[c][row] = b[static_cast<size_t>(row) * nrhs + c];
    }
    return true;
}

/// Direct solve of the periodic block-tridiagonal system for every supplied
/// right-hand side at once (they share one factorisation).
///
/// Blocks 0 .. N-2 are swept with block Thomas, carrying n_r extra right-hand-side
/// columns that express their dependence on the last block:
///     X_i = Y_i + Z_i X_{N-1}.
/// Substituting into the last equation leaves one dense n_r x n_r Schur system
///     (B_{N-1} - diag(gt)(Z_{N-2} + Z_0)) X_{N-1} = D_{N-1} + diag(gt)(Y_{N-2} + Y_0),
/// which closes the periodic wrap exactly.
bool solve_blocks(const Operator& op, std::vector<std::vector<double>>& rhs) {
    const int N = op.n_theta;
    const int nr = op.n_r;
    if (N <= 2) return solve_dense(op, rhs);

    const int n_src = static_cast<int>(rhs.size());
    const int nrhs = n_src + nr;
    const int sub = N - 1;  // blocks 0 .. N-2

    std::vector<double> e_store(static_cast<size_t>(sub) * nr * nr, 0.0);
    std::vector<double> f_store(static_cast<size_t>(sub) * nr * nrhs, 0.0);
    std::vector<double> w(static_cast<size_t>(nr) * nr), block(static_cast<size_t>(nr) * nr);
    std::vector<double> work(static_cast<size_t>(nr) * nrhs);
    std::vector<double> lu(static_cast<size_t>(nr) * nr), c_block(static_cast<size_t>(nr) * nr);
    std::vector<int> piv;

    for (int i = 0; i < sub; ++i) {
        op.block(i, block);
        double* e_prev = (i > 0) ? &e_store[static_cast<size_t>(i - 1) * nr * nr] : nullptr;
        double* f_prev = (i > 0) ? &f_store[static_cast<size_t>(i - 1) * nr * nrhs] : nullptr;
        // W_i = B_i - A E_{i-1}, with A = -diag(gt): a row scaling, no matrix product.
        w = block;
        if (e_prev) {
            for (int m = 0; m < nr; ++m) {
                for (int c = 0; c < nr; ++c) {
                    w[static_cast<size_t>(m) * nr + c] += op.gt[m] * e_prev[static_cast<size_t>(m) * nr + c];
                }
            }
        }
        // Right-hand side columns: the sources, then the coupling to X_{N-1},
        // which touches only the first (through A_0) and last (through C_{N-2}) blocks.
        std::fill(work.begin(), work.end(), 0.0);
        for (int m = 0; m < nr; ++m) {
            for (int c = 0; c < n_src; ++c) {
                work[static_cast<size_t>(m) * nrhs + c] = rhs[c][op.at(i, m)];
            }
            if (i == 0 || i == sub - 1) {
                work[static_cast<size_t>(m) * nrhs + n_src + m] += op.gt[m];
            }
            if (f_prev) {
                for (int c = 0; c < nrhs; ++c) {
                    work[static_cast<size_t>(m) * nrhs + c] +=
                        op.gt[m] * f_prev[static_cast<size_t>(m) * nrhs + c];
                }
            }
        }
        lu = w;
        if (!lu_factor(lu, nr, piv)) return false;
        lu_solve(lu, nr, piv, work, nrhs);
        std::copy(work.begin(), work.end(), f_store.begin() + static_cast<size_t>(i) * nr * nrhs);
        if (i + 1 < sub) {
            // E_i = W_i^{-1} C_i with C_i = -diag(gt).
            std::fill(c_block.begin(), c_block.end(), 0.0);
            for (int m = 0; m < nr; ++m) c_block[static_cast<size_t>(m) * nr + m] = -op.gt[m];
            lu_solve(lu, nr, piv, c_block, nr);
            std::copy(c_block.begin(), c_block.end(),
                      e_store.begin() + static_cast<size_t>(i) * nr * nr);
        }
    }

    // Back-substitute the sub-system: X_i = F_i - E_i X_{i+1}.
    for (int i = sub - 2; i >= 0; --i) {
        double* f_i = &f_store[static_cast<size_t>(i) * nr * nrhs];
        const double* f_next = &f_store[static_cast<size_t>(i + 1) * nr * nrhs];
        const double* e_i = &e_store[static_cast<size_t>(i) * nr * nr];
        for (int m = 0; m < nr; ++m) {
            for (int k = 0; k < nr; ++k) {
                const double g = e_i[static_cast<size_t>(m) * nr + k];
                if (g == 0.0) continue;
                for (int c = 0; c < nrhs; ++c) {
                    f_i[static_cast<size_t>(m) * nrhs + c] -=
                        g * f_next[static_cast<size_t>(k) * nrhs + c];
                }
            }
        }
    }

    // Schur complement for the last block.
    op.block(N - 1, block);
    const double* f0 = &f_store[0];
    const double* fl = &f_store[static_cast<size_t>(sub - 1) * nr * nrhs];
    std::vector<double> schur = block;
    std::vector<double> schur_rhs(static_cast<size_t>(nr) * n_src, 0.0);
    for (int m = 0; m < nr; ++m) {
        for (int c = 0; c < nr; ++c) {
            schur[static_cast<size_t>(m) * nr + c] -=
                op.gt[m] * (f0[static_cast<size_t>(m) * nrhs + n_src + c] +
                            fl[static_cast<size_t>(m) * nrhs + n_src + c]);
        }
        for (int c = 0; c < n_src; ++c) {
            schur_rhs[static_cast<size_t>(m) * n_src + c] =
                rhs[c][op.at(N - 1, m)] +
                op.gt[m] * (f0[static_cast<size_t>(m) * nrhs + c] + fl[static_cast<size_t>(m) * nrhs + c]);
        }
    }
    if (!lu_factor(schur, nr, piv)) return false;
    lu_solve(schur, nr, piv, schur_rhs, n_src);

    for (int c = 0; c < n_src; ++c) {
        for (int m = 0; m < nr; ++m) {
            rhs[c][op.at(N - 1, m)] = schur_rhs[static_cast<size_t>(m) * n_src + c];
        }
    }
    for (int i = 0; i < sub; ++i) {
        const double* f_i = &f_store[static_cast<size_t>(i) * nr * nrhs];
        for (int m = 0; m < nr; ++m) {
            for (int c = 0; c < n_src; ++c) {
                double v = f_i[static_cast<size_t>(m) * nrhs + c];
                for (int k = 0; k < nr; ++k) {
                    v += f_i[static_cast<size_t>(m) * nrhs + n_src + k] *
                         schur_rhs[static_cast<size_t>(k) * n_src + c];
                }
                rhs[c][op.at(i, m)] = v;
            }
        }
    }
    return true;
}

}  // namespace

Solution solve(const Coupling& coupling, int n_theta, double r_inner,
               const SimulationConfig& cfg) {
    Solution out;
    const int n_r = std::max(1, cfg.bush_radial_cells);
    out.alpha.assign(std::max(0, n_theta), 0.0);
    out.beta.assign(std::max(0, n_theta), 0.0);
    if (n_theta <= 0) return out;

    const double k_b = std::max(cfg.bush_thermal_conductivity, 1.0e-30);
    const double r_outer = cfg.bush_outer_radius > r_inner ? cfg.bush_outer_radius
                                                           : 1.1 * r_inner;
    const double d_theta = 2.0 * M_PI / n_theta;
    const double d_r = (r_outer - r_inner) / n_r;
    const double t_amb = cfg.thermal_ambient_temperature;

    Operator op;
    op.n_theta = n_theta;
    op.n_r = n_r;
    std::vector<double> r_c(n_r);
    for (int m = 0; m < n_r; ++m) r_c[m] = r_inner + (m + 0.5) * d_r;

    // Radial faces use the exact logarithmic shell conductance, so a purely
    // radial solution is reproduced to round-off at any radial resolution.
    op.gr.assign(std::max(0, n_r - 1), 0.0);
    for (int m = 0; m + 1 < n_r; ++m) op.gr[m] = k_b * d_theta / std::log(r_c[m + 1] / r_c[m]);
    op.gt.assign(n_r, 0.0);
    for (int m = 0; m < n_r; ++m) op.gt[m] = k_b * d_r / (r_c[m] * d_theta);

    op.s_inner = r_inner * d_theta;
    op.c_half_inner = k_b * d_theta / std::log(r_c[0] / r_inner);
    const double c_half_outer = k_b * d_theta / std::log(r_outer / r_c[n_r - 1]);
    if (cfg.bush_outer_boundary == BushOuterBoundary::FIXED_TEMPERATURE) {
        op.outer_u = c_half_outer;  // the h_out -> infinity limit of the Robin form
    } else if (cfg.bush_outer_heat_transfer > 0.0) {
        const double film = cfg.bush_outer_heat_transfer * r_outer * d_theta;
        op.outer_u = 1.0 / (1.0 / c_half_outer + 1.0 / film);
    } else {
        op.outer_u = 0.0;  // adiabatic outer surface
    }

    // Inner Robin condition. The bore flux q = A - B T_bore is put in series with
    // the bore-to-first-node half cell, which eliminates the bore temperature
    // from the linear system without losing accuracy:
    //   Q_into_cell0 = lambda * S * (A - B T_0),   lambda = C_half/(C_half + S B).
    op.inner_gain.assign(n_theta, 0.0);
    op.inner_lambda.assign(n_theta, 0.0);
    for (int i = 0; i < n_theta; ++i) {
        const double slope = std::max(coupling.slope[i], 0.0);
        const double denom = op.c_half_inner + op.s_inner * slope;
        const double lambda = denom > 0.0 ? op.c_half_inner / denom : 0.0;
        op.inner_lambda[i] = lambda;
        op.inner_gain[i] = lambda * op.s_inner * slope;
    }

    op.diag.assign(static_cast<size_t>(n_theta) * n_r, 0.0);
    for (int i = 0; i < n_theta; ++i) {
        for (int m = 0; m < n_r; ++m) {
            double d = 2.0 * op.gt[m];
            if (m > 0) d += op.gr[m - 1];
            if (m + 1 < n_r) d += op.gr[m];
            if (m == 0) d += op.inner_gain[i];
            if (m == n_r - 1) d += op.outer_u;
            op.diag[op.at(i, m)] = d;
        }
    }

    // Two right-hand sides: the shaft-independent part and the shaft slope. One
    // factorisation serves both, because the film-side slope B does not depend on
    // the shaft temperature.
    const size_t n = static_cast<size_t>(n_theta) * n_r;
    std::vector<std::vector<double>> rhs(2, std::vector<double>(n, 0.0));
    for (int i = 0; i < n_theta; ++i) {
        rhs[0][op.at(i, 0)] += op.inner_lambda[i] * op.s_inner * coupling.a0[i];
        rhs[1][op.at(i, 0)] += op.inner_lambda[i] * op.s_inner * coupling.a1[i];
        rhs[0][op.at(i, n_r - 1)] += op.outer_u * t_amb;  // ambient enters the constant part only
    }
    const std::vector<double> src0 = rhs[0];
    const std::vector<double> src1 = rhs[1];

    out.converged = solve_blocks(op, rhs);
    if (!out.converged) {
        out.residual = 1.0;
        return out;
    }
    out.residual = std::max(relative_residual(op, src0, rhs[0]),
                            relative_residual(op, src1, rhs[1]));
    out.converged = out.residual <= kResidualTolerance;

    // Recover the bore temperature from the first node and the half-cell flux,
    // T_bore = T_0 + Q/C_half with Q = lambda S (A - B T_0).
    for (int i = 0; i < n_theta; ++i) {
        const double slope = std::max(coupling.slope[i], 0.0);
        const double lambda = op.inner_lambda[i];
        const double t0_alpha = rhs[0][op.at(i, 0)];
        const double t0_beta = rhs[1][op.at(i, 0)];
        const double q_alpha = lambda * op.s_inner * (coupling.a0[i] - slope * t0_alpha);
        const double q_beta = lambda * op.s_inner * (coupling.a1[i] - slope * t0_beta);
        out.alpha[i] = t0_alpha + q_alpha / op.c_half_inner;
        out.beta[i] = t0_beta + q_beta / op.c_half_inner;
    }
    return out;
}

}  // namespace BushConduction
