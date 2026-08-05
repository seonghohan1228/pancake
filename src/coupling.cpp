#include "coupling.hpp"

#include <algorithm>
#include <sstream>

namespace Coupling {

Controller::Controller(const SimulationConfig& cfg) {
    tolerance_ = cfg.coupling_tolerance;
    max_iters_ = std::max(1, cfg.coupling_max_iters);
    keep_best_ = cfg.coupling_keep_best_iterate;
    adaptive_ = cfg.coupling_relaxation_mode == CouplingRelaxationMode::ADAPTIVE;
    omega_max_ = std::clamp(cfg.coupling_relaxation, 0.0, 1.0);
    omega_min_ = adaptive_ ? std::clamp(cfg.coupling_relaxation_min, 0.0, omega_max_) : omega_max_;
    backoff_ = cfg.coupling_relaxation_backoff;
    growth_ = cfg.coupling_relaxation_growth;
    omega_ = omega_max_;
    status_.max_iters = max_iters_;
    status_.relaxation = omega_;
    status_.relaxation_min_used = omega_;
}

bool Controller::record(int outer, double R_p, double R_T) {
    const double r = residual_metric(R_p, R_T);

    // Relaxation schedule. FIXED leaves omega_ at its constructed value, so the
    // sequence of updates is byte-for-byte the historical one.
    if (adaptive_ && any_recorded_) {
        omega_ = (r > prev_residual_) ? std::max(omega_min_, omega_ * backoff_)
                                      : std::min(omega_max_, omega_ * growth_);
    }
    status_.relaxation = omega_;
    status_.relaxation_min_used = std::min(status_.relaxation_min_used, omega_);

    // Best-iterate bookkeeping. Strictly-less keeps the FIRST iterate that
    // attained the minimum, which is the cheapest one on a limit cycle.
    is_best_ = !any_recorded_ || r < status_.best_residual;
    if (is_best_) {
        status_.best_iter = outer + 1;
        status_.best_R_p = R_p;
        status_.best_R_T = R_T;
        status_.best_residual = r;
    }

    status_.iters = outer + 1;
    status_.R_p = R_p;
    status_.R_T = R_T;
    status_.converged = (R_p < tolerance_) && (R_T < tolerance_);

    prev_residual_ = r;
    any_recorded_ = true;
    return status_.converged;
}

bool Controller::should_restore_best() const {
    if (!keep_best_ || !any_recorded_) return false;
    if (status_.converged) return false;          // the last iterate is acceptable
    return status_.best_iter != status_.iters;    // a strictly better one exists
}

void Controller::finish(int iters_done) {
    status_.iters = iters_done;
    status_.cap_exhausted = !status_.converged && iters_done >= max_iters_;
    status_.restored_best = should_restore_best();
}

std::string Controller::summary_line() const {
    std::ostringstream s;
    s << "Coupling: " << status_.iters << " outer iters, R_p=" << status_.R_p
      << " R_T=" << status_.R_T
      << " | " << (status_.converged ? "CONVERGED" : "NOT CONVERGED")
      << " (tol=" << tolerance_ << ")";
    if (!status_.converged && status_.cap_exhausted) {
        s << " -- hit coupling_max_iters=" << max_iters_;
    }
    if (adaptive_) {
        s << " | omega=" << status_.relaxation << " (min " << status_.relaxation_min_used << ")";
    }
    return s.str();
}

std::string Controller::rollback_line() const {
    std::ostringstream s;
    if (status_.restored_best) {
        s << "Coupling rollback: delivering outer " << status_.best_iter
          << " (residual " << status_.best_residual << ", R_p=" << status_.best_R_p
          << " R_T=" << status_.best_R_T << ") instead of outer " << status_.iters
          << " (residual " << status_.last_residual() << "); the final iterate is "
          << status_.best_to_last_ratio() << "x worse. The loop is not converging - "
          << "treat this result as NON-CONVERGED.";
    } else if (!status_.converged) {
        s << "Coupling did not converge and no better iterate than the last was seen "
          << "(best outer " << status_.best_iter << ", residual " << status_.best_residual
          << "). Treat this result as NON-CONVERGED.";
    }
    return s.str();
}

void Snapshot::capture(const Fields& fields) {
    entries_.clear();
    for (const auto& kv : fields) {
        Entry e;
        e.name = kv.first;
        e.data = kv.second->data;
        e.old_data = kv.second->old_data;
        entries_.push_back(std::move(e));
    }
    valid_ = true;
}

void Snapshot::restore(Fields& fields) const {
    if (!valid_) return;
    for (const auto& e : entries_) {
        if (!fields.has(e.name)) continue;
        Field& f = fields[e.name];
        if (f.data.size() == e.data.size()) f.data = e.data;
        if (f.old_data.size() == e.old_data.size()) f.old_data = e.old_data;
    }
}

std::size_t Snapshot::bytes() const {
    std::size_t n = 0;
    for (const auto& e : entries_) n += (e.data.size() + e.old_data.size()) * sizeof(double);
    return n;
}

}  // namespace Coupling
