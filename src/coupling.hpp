#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "config.hpp"
#include "field.hpp"
#include "mesh.hpp"

/// Control of the STEADY_STATE outer (Picard) coupling loop: convergence test,
/// under-relaxation schedule, best-iterate bookkeeping and the honesty flag that
/// says whether the delivered state actually met coupling_tolerance.
///
/// Why this is a module and not inline in main(): the loop has three separable
/// failure modes (limit cycle, slow convergence, divergence) and only the first
/// one is what the conjugate Ferron cases hit. Isolating the decision logic makes
/// it testable against a synthetic residual history without running a solver.
///
/// Usage, once per outer iteration, after the segregated sweep has finished and
/// the raw (unrelaxed) pressure/temperature residuals have been formed:
///
///     const bool stop = control.record(outer, R_p, R_T);
///     if (control.is_best()) snapshot.capture(fields, elrod_stats, wall_state);
///     relax_field(p, p_prev, control.relaxation(), mesh);
///     relax_field(T, T_prev, control.relaxation(), mesh);
///     if (stop) break;
///
/// and after the loop:
///
///     if (control.should_restore_best()) snapshot.restore(fields, ...);
///
/// The snapshot is taken BEFORE relaxation on purpose: the unrelaxed end-of-sweep
/// state is the one the solvers actually produced and the one the recorded
/// residual describes. A relaxed blend of two sweeps is neither.
namespace Coupling {

/// Scalar used to rank iterates. Both residuals are normalized and both are
/// tested against the same coupling_tolerance, so the max is the natural metric:
/// it is exactly the quantity that has to fall below the tolerance.
inline double residual_metric(double R_p, double R_T, double R_g = 0.0) {
    if (!std::isfinite(R_p) || !std::isfinite(R_T) || !std::isfinite(R_g))
        return std::numeric_limits<double>::infinity();
    return std::max({R_p, R_T, R_g});
}

/// Outcome of one outer loop, reported to the log, the diagnostics CSV and the
/// case drivers.
struct Status {
    int    iters = 0;              ///< outer sweeps performed
    int    max_iters = 0;          ///< the cap that was in force
    double R_p = 0.0;              ///< residuals of the LAST sweep
    double R_T = 0.0;
    double R_g = 0.0;              ///< dissolved concentration / gas-volume change
    bool   converged = false;      ///< the last sweep met coupling_tolerance
    bool   cap_exhausted = false;  ///< the loop stopped because it ran out of iterations
    int    best_iter = 0;          ///< 1-based index of the lowest-residual sweep
    double best_R_p = 0.0;
    double best_R_T = 0.0;
    double best_R_g = 0.0;
    double best_residual = 0.0;    ///< residual_metric at best_iter
    bool   restored_best = false;  ///< the delivered state was rolled back to best_iter
    double relaxation = 1.0;       ///< relaxation factor in force at the end
    double relaxation_min_used = 1.0;

    /// Residual of the last sweep under the ranking metric.
    double last_residual() const { return residual_metric(R_p, R_T, R_g); }
    /// How much worse the final iterate is than the best one found.
    double best_to_last_ratio() const {
        return best_residual > 0.0 ? last_residual() / best_residual : 1.0;
    }
};

/// Convergence test + relaxation schedule + best-iterate bookkeeping.
///
/// FIXED mode reproduces the historical behaviour exactly: relaxation() always
/// returns cfg.coupling_relaxation and the only new behaviour is the bookkeeping.
/// ADAPTIVE mode backs the factor off whenever the residual rises (the limit-cycle
/// signature) and lets it recover while the residual falls.
class Controller {
public:
    explicit Controller(const SimulationConfig& cfg);

    /// Record the residuals of the sweep numbered `outer` (0-based). Returns true
    /// when coupling_tolerance and any additional conservation check are met.
    bool record(int outer, double R_p, double R_T, double R_g = 0.0,
                bool additional_balance_converged = true);

    /// True when the sweep just recorded is the best seen so far, i.e. the caller
    /// should snapshot it. Always true for the first sweep.
    bool is_best() const { return is_best_; }

    /// Relaxation factor to apply to this sweep's pressure/temperature update.
    double relaxation() const { return omega_; }

    /// True when the loop ended on the iteration cap without meeting tolerance
    /// AND a strictly better iterate than the last one was seen. Only then is a
    /// rollback both necessary and a change of behaviour; a converged run returns
    /// its last iterate bit-identically to before this module existed.
    bool should_restore_best() const;

    /// Call once after the loop; finalizes iters/cap_exhausted/restored_best.
    void finish(int iters_done);

    const Status& status() const { return status_; }

    /// One-line summary for the run log. Loud (and marked) when non-converged.
    std::string summary_line() const;

    /// Extra diagnostic naming best vs final, emitted only when they differ.
    std::string rollback_line() const;

private:
    double tolerance_ = 1e-6;
    int    max_iters_ = 1;
    bool   keep_best_ = true;
    bool   adaptive_ = false;
    double omega_ = 1.0;
    double omega_max_ = 1.0;
    double omega_min_ = 1.0;
    double backoff_ = 0.5;
    double growth_ = 1.1;
    double prev_residual_ = -1.0;
    bool   is_best_ = false;
    bool   any_recorded_ = false;
    Status status_;
};

/// Deep copy of every registered field (data + old_data), plus the per-sweep
/// solver stats that travel with it, so a whole outer iterate can be restored.
///
/// Copies the raw storage including ghost layers, so a restored state is
/// bit-identical to the captured one and needs no ghost exchange to be usable.
class Snapshot {
public:
    /// Capture the physical + ghost storage of every field in `fields`.
    void capture(const Fields& fields);
    /// Restore into `fields`. Fields absent from the snapshot are left alone;
    /// the field set is fixed after registration, so in practice it matches.
    void restore(Fields& fields) const;

    bool valid() const { return valid_; }
    /// Bytes held, for the log line (this is per rank).
    std::size_t bytes() const;

private:
    struct Entry {
        std::string name;
        std::vector<double> data;
        std::vector<double> old_data;
        bool history_initialized = false;
    };
    std::vector<Entry> entries_;
    bool valid_ = false;
};

}  // namespace Coupling
