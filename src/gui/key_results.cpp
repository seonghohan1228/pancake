#include "key_results.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace gui {

namespace {
constexpr double kPi = 3.14159265358979323846;

double wrap180(double deg) {
    double out = std::fmod(deg + 180.0, 360.0);
    if (out < 0.0) out += 360.0;
    return out - 180.0;
}

double first_cell(const ResultGrid& grid, const char* name) {
    const auto* values = grid.field(name);
    return values && !values->empty() ? (*values)[0] : KeyResults::kNaN;
}
}  // namespace

KeyResults compute_key_results(const ResultGrid& grid, const SimulationConfig& cfg) {
    KeyResults r;
    r.time = grid.time;
    r.step = grid.step;
    if (const auto* h = grid.field("h")) {
        double h_min = std::numeric_limits<double>::infinity();
        for (double v : *h)
            if (std::isfinite(v)) h_min = std::min(h_min, v);
        if (std::isfinite(h_min)) r.h_min = h_min;
        // Angular positions of the thinnest and thickest film along the mid-plane.
        const int j = grid.ny / 2;
        int i_min = -1, i_max = -1;
        for (int i = 0; i < grid.nx; ++i) {
            const double v = (*h)[static_cast<std::size_t>(j) * grid.nx + i];
            if (!std::isfinite(v)) continue;
            if (i_min < 0 || v < (*h)[static_cast<std::size_t>(j) * grid.nx + i_min]) i_min = i;
            if (i_max < 0 || v > (*h)[static_cast<std::size_t>(j) * grid.nx + i_max]) i_max = i;
        }
        if (i_min >= 0) r.h_min_theta_deg = grid.theta_deg[i_min];
        if (i_max >= 0) r.h_max_theta_deg = grid.theta_deg[i_max];
    }
    if (const auto* p = grid.field("p")) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (double v : *p) {
            if (!std::isfinite(v)) continue;
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        if (std::isfinite(lo)) r.p_min = lo;
        if (std::isfinite(hi)) r.p_max = hi;
    }
    if (const auto* film = grid.field("film_content")) {
        std::size_t cavitated = 0, total = 0;
        for (double v : *film) {
            if (!std::isfinite(v)) continue;
            ++total;
            if (v < 1.0 - 1e-8) ++cavitated;
        }
        if (total) r.cavitated_fraction = static_cast<double>(cavitated) / total;
    }
    if (const auto* t = grid.field("T")) {
        double hi = -std::numeric_limits<double>::infinity();
        for (double v : *t)
            if (std::isfinite(v)) hi = std::max(hi, v);
        if (std::isfinite(hi)) r.t_max = hi;
    }

    r.friction_torque = first_cell(grid, "friction_torque");
    if (std::isfinite(r.friction_torque)) r.friction_power = std::abs(r.friction_torque * cfg.omega);

    const double fx = first_cell(grid, "F.x"), fy = first_cell(grid, "F.y");
    if (std::isfinite(fx) && std::isfinite(fy)) {
        r.film_force = std::hypot(fx, fy);
        const double vx = first_cell(grid, "Fv.x"), vy = first_cell(grid, "Fv.y");
        if (std::isfinite(vx) && std::isfinite(vy) && r.film_force > 0.0)
            r.shear_share = std::hypot(vx, vy) / r.film_force;
    }

    const double bx = first_cell(grid, "xB.x"), by = first_cell(grid, "xB.y");
    if (std::isfinite(bx) && std::isfinite(by) && cfg.c > 0.0) {
        r.eccentricity_ratio = std::hypot(bx, by) / cfg.c;
        // Journal displacement relative to the bearing is -b; the load on the
        // bearing (+W) acts on the journal as -W.
        const double load_x = std::isfinite(first_cell(grid, "Fext.x")) ? first_cell(grid, "Fext.x")
                                                                         : cfg.external_load_x;
        const double load_y = std::isfinite(first_cell(grid, "Fext.y")) ? first_cell(grid, "Fext.y")
                                                                         : cfg.external_load_y;
        if (std::hypot(load_x, load_y) > 0.0 && r.eccentricity_ratio > 0.0) {
            const double journal_dir = std::atan2(-by, -bx) * 180.0 / kPi;
            const double load_on_journal = std::atan2(-load_y, -load_x) * 180.0 / kPi;
            r.attitude_deg = std::abs(wrap180(journal_dir - load_on_journal));
        }
    }

    const double sigma = std::hypot(cfg.roughness_rq_shaft, cfg.roughness_rq_bush);
    if (sigma > 0.0 && std::isfinite(r.h_min)) r.lambda = r.h_min / sigma;
    return r;
}

DesignChecks design_checks(const KeyResults& r, double p_sat, bool gas_model_on) {
    DesignChecks d;
    char text[160];
    if (std::isfinite(r.lambda)) {
        d.lambda_risk = r.lambda < 3.0 ? Risk::High : (r.lambda < 5.0 ? Risk::Medium : Risk::Low);
        std::snprintf(text, sizeof(text), "%s (Lambda = %.1f). Contact is not modelled.",
                      r.lambda < 3.0 ? "Mixed lubrication likely"
                                     : (r.lambda < 5.0 ? "Thin film" : "Full film"),
                      r.lambda);
        d.lambda_text = text;
    }
    if (p_sat > 0.0 && std::isfinite(r.p_min)) {
        // Gas-off run: the film falls below p_sat -> gas would be released.
        // Gas-on run: the cavitated film sits at p_sat; any cavitated area means release.
        const bool release = gas_model_on ? (std::isfinite(r.cavitated_fraction) &&
                                             r.cavitated_fraction > 0.0 && r.p_min <= p_sat * 1.0001)
                                          : r.p_min < p_sat;
        d.gas_release_risk = release ? Risk::High : Risk::Low;
        std::snprintf(text, sizeof(text), "%s: film p_min %.2f bar, p_sat %.2f bar.",
                      release ? "Gas release in the film" : "No gas release",
                      r.p_min * 1e-5, p_sat * 1e-5);
        d.gas_release_text = text;
    }
    return d;
}

double percent_difference(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b) || a == 0.0) return KeyResults::kNaN;
    return (b - a) / std::abs(a) * 100.0;
}

}  // namespace gui
