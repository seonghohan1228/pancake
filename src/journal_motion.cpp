#include "journal_motion.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace JournalMotion {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct AxisParams {
    double mass = 1.0;
    double damping = 0.0;
    double stiffness = 0.0;
    double force = 0.0;
};

struct AxisState {
    double x = 0.0;
    double v = 0.0;
};

double acceleration(AxisState state, AxisParams params) {
    return (params.force - params.damping * state.v - params.stiffness * state.x) / params.mass;
}

AxisState explicit_euler(AxisState state, AxisParams params, double dt) {
    return {state.x + dt * state.v,
            state.v + dt * acceleration(state, params)};
}

AxisState implicit_euler(AxisState state, AxisParams params, double dt) {
    const double denom = 1.0
        + dt * params.damping / params.mass
        + dt * dt * params.stiffness / params.mass;
    const double numer = state.v
        + dt * (params.force - params.stiffness * state.x) / params.mass;
    const double v_new = numer / denom;
    return {state.x + dt * v_new, v_new};
}

AxisState crank_nicolson(AxisState state, AxisParams params, double dt) {
    const double denom = 1.0
        + dt * (0.5 * params.damping + 0.25 * params.stiffness * dt) / params.mass;
    const double numer = state.v
        + dt * (params.force
                - 0.5 * params.damping * state.v
                - params.stiffness * state.x
                - 0.25 * params.stiffness * dt * state.v) / params.mass;
    const double v_new = numer / denom;
    return {state.x + 0.5 * dt * (state.v + v_new), v_new};
}

AxisState rk2(AxisState state, AxisParams params, double dt) {
    const AxisState k1{state.v, acceleration(state, params)};
    const AxisState mid{state.x + 0.5 * dt * k1.x, state.v + 0.5 * dt * k1.v};
    const AxisState k2{mid.v, acceleration(mid, params)};
    return {state.x + dt * k2.x, state.v + dt * k2.v};
}

AxisState rk4(AxisState state, AxisParams params, double dt) {
    const AxisState k1{state.v, acceleration(state, params)};
    const AxisState s2{state.x + 0.5 * dt * k1.x, state.v + 0.5 * dt * k1.v};
    const AxisState k2{s2.v, acceleration(s2, params)};
    const AxisState s3{state.x + 0.5 * dt * k2.x, state.v + 0.5 * dt * k2.v};
    const AxisState k3{s3.v, acceleration(s3, params)};
    const AxisState s4{state.x + dt * k3.x, state.v + dt * k3.v};
    const AxisState k4{s4.v, acceleration(s4, params)};

    return {
        state.x + dt * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) / 6.0,
        state.v + dt * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v) / 6.0};
}

AxisState advance_axis(AxisState state, AxisParams params,
                       TimeSteppingMethod method, double dt) {
    if (params.mass <= 0.0) {
        throw std::runtime_error("bearing_mass must be positive when moving-bearing motion is enabled");
    }

    switch (method) {
        case TimeSteppingMethod::EULER_EXPLICIT:
            return explicit_euler(state, params, dt);
        case TimeSteppingMethod::EULER_IMPLICIT:
            return implicit_euler(state, params, dt);
        case TimeSteppingMethod::CRANK_NICOLSON:
            return crank_nicolson(state, params, dt);
        case TimeSteppingMethod::RK2:
            return rk2(state, params, dt);
        case TimeSteppingMethod::RK4:
            return rk4(state, params, dt);
    }
    return implicit_euler(state, params, dt);
}

void fill_if_present(Fields& fields, const std::string& name, double value) {
    if (fields.has(name)) fields[name].fill(value);
}

}  // namespace

BearingState initial_state(const SimulationConfig& cfg) {
    BearingState state;
    if (cfg.bearing_initial_from_attitude) {
        const double psi = cfg.attitude_angle_deg * kPi / 180.0;
        state.x = -cfg.e * std::cos(psi);
        state.y = -cfg.e * std::sin(psi);
    } else {
        state.x = cfg.bearing_initial_x;
        state.y = cfg.bearing_initial_y;
    }
    state.z = cfg.bearing_initial_z;
    state.vx = cfg.bearing_initial_vx;
    state.vy = cfg.bearing_initial_vy;
    state.vz = cfg.bearing_initial_vz;
    return state;
}

namespace {

bool parse_load_row(const std::string& raw, double values[3]) {
    std::string line = raw.substr(0, raw.find('#'));
    std::replace(line.begin(), line.end(), ',', ' ');
    std::replace(line.begin(), line.end(), ';', ' ');
    std::replace(line.begin(), line.end(), '\t', ' ');
    std::istringstream in(line);
    std::string token;
    int count = 0;
    while (in >> token) {
        if (count == 3) return false;
        try {
            std::size_t used = 0;
            values[count] = std::stod(token, &used);
            if (used != token.size()) return false;
        } catch (const std::exception&) {
            return false;
        }
        ++count;
    }
    return count == 3;
}

bool blank_or_comment(const std::string& raw) {
    const std::string line = raw.substr(0, raw.find('#'));
    return line.find_first_not_of(" \t\r\n,;") == std::string::npos;
}

}  // namespace

std::vector<LoadTableRow> read_load_table(const std::filesystem::path& path,
                                          const SimulationConfig& cfg, std::string& error) {
    std::vector<LoadTableRow> rows;
    std::ifstream file(path);
    if (!file) {
        error = "cannot open " + path.string();
        return {};
    }
    const double offset = cfg.load_table_angle_offset_deg * kPi / 180.0;
    std::string line;
    int line_number = 0;
    bool data_started = false;
    while (std::getline(file, line)) {
        ++line_number;
        if (blank_or_comment(line)) continue;
        double v[3];
        if (!parse_load_row(line, v)) {
            if (!data_started) continue;  // header lines before the first data row
            error = path.string() + " line " + std::to_string(line_number) +
                    ": expected three numbers (abscissa and two load values)";
            return {};
        }
        data_started = true;
        double x = 0.0, y = 0.0;
        if (cfg.load_table_format == LoadTableFormat::POLAR) {
            const double angle = v[2] * kPi / 180.0;
            x = v[1] * std::cos(angle);
            y = v[1] * std::sin(angle);
        } else {
            x = v[1];
            y = v[2];
        }
        x *= cfg.load_table_scale;
        y *= cfg.load_table_scale;
        rows.push_back({v[0], x * std::cos(offset) - y * std::sin(offset),
                        x * std::sin(offset) + y * std::cos(offset)});
        if (!std::isfinite(v[0]) || !std::isfinite(rows.back().load_x) || !std::isfinite(rows.back().load_y)) {
            error = path.string() + " line " + std::to_string(line_number) + ": non-finite value";
            return {};
        }
    }
    if (rows.empty()) {
        error = path.string() + ": no data rows";
        return {};
    }
    for (std::size_t k = 1; k < rows.size(); ++k) {
        if (!(rows[k].abscissa > rows[k - 1].abscissa)) {
            error = path.string() + ": abscissa must increase strictly (row " + std::to_string(k + 1) + ")";
            return {};
        }
    }
    if (cfg.load_table_abscissa == LoadTableAbscissa::PHASE_DEG &&
        (rows.front().abscissa < 0.0 || rows.back().abscissa >= 360.0)) {
        error = path.string() + ": PHASE_DEG rows must lie in [0, 360); the 360 deg row is the 0 deg row of the next cycle";
        return {};
    }
    return rows;
}

void load_config_load_table(SimulationConfig& cfg, std::vector<std::string>& errors) {
    cfg.load_table.clear();
    if (cfg.load_model != LoadModel::TABLE || cfg.load_table_file.empty()) return;
    std::filesystem::path path = cfg.load_table_file;
    if (!path.is_absolute() && !cfg.config_dir.empty()) {
        std::error_code ec;
        const std::filesystem::path candidate = cfg.config_dir / path;
        if (std::filesystem::exists(candidate, ec)) path = candidate;
    }
    std::string error;
    cfg.load_table = read_load_table(path, cfg, error);
    if (!error.empty()) errors.push_back("load_table_file: " + error);
}

double load_phase_deg(const SimulationConfig& cfg, double time) {
    if (cfg.load_model != LoadModel::TABLE || cfg.load_table_abscissa != LoadTableAbscissa::PHASE_DEG)
        return std::numeric_limits<double>::quiet_NaN();
    double phase = std::fmod(cfg.load_phase_offset_deg + 360.0 * cfg.load_frequency * time, 360.0);
    if (phase < 0.0) phase += 360.0;
    return phase >= 360.0 ? 0.0 : phase;
}

Vector3 external_load(const SimulationConfig& cfg, double time) {
    Vector3 load{cfg.external_load_x, cfg.external_load_y, cfg.external_load_z};
    if (cfg.load_model == LoadModel::ROTATING) {
        const double angle = cfg.load_rotation_start_deg * kPi / 180.0 + cfg.load_rotation_speed * time;
        load.x = cfg.load_magnitude * std::cos(angle);
        load.y = cfg.load_magnitude * std::sin(angle);
    } else if (cfg.load_model == LoadModel::TABLE) {
        const auto& rows = cfg.load_table;
        if (rows.empty()) throw std::runtime_error("TABLE load selected but load_table is empty");
        if (rows.size() == 1) {
            load.x = rows.front().load_x;
            load.y = rows.front().load_y;
            return load;
        }
        // Interval [a, b] containing s and the two rows it joins. Interpolating
        // the Cartesian parts keeps a direction that crosses 360 deg continuous.
        double s = 0.0, a = 0.0, b = 0.0;
        const LoadTableRow* lo = nullptr;
        const LoadTableRow* hi = nullptr;
        if (cfg.load_table_abscissa == LoadTableAbscissa::PHASE_DEG) {
            s = load_phase_deg(cfg, time);
            if (s < rows.front().abscissa || s >= rows.back().abscissa) {
                lo = &rows.back();
                hi = &rows.front();
                a = rows.back().abscissa;
                b = rows.front().abscissa + 360.0;
                if (s < rows.front().abscissa) s += 360.0;
            }
        } else {
            s = time;
            const double tol = 1e-12 * std::max(1.0, std::abs(rows.back().abscissa));
            if (s < rows.front().abscissa - tol || s > rows.back().abscissa + tol)
                throw std::runtime_error("time " + std::to_string(time) + " s lies outside the TIME_S load table");
            s = std::clamp(s, rows.front().abscissa, rows.back().abscissa);
        }
        if (!lo) {
            auto upper = std::upper_bound(rows.begin(), rows.end(), s,
                [](double value, const LoadTableRow& row) { return value < row.abscissa; });
            if (upper == rows.end()) --upper;
            if (upper == rows.begin()) ++upper;
            hi = &*upper;
            lo = &*(upper - 1);
            a = lo->abscissa;
            b = hi->abscissa;
        }
        const double w = (s - a) / (b - a);
        load.x = (1.0 - w) * lo->load_x + w * hi->load_x;
        load.y = (1.0 - w) * lo->load_y + w * hi->load_y;
    }
    return load;
}

Vector3 fluid_force_from_fields(const Fields& fields) {
    Vector3 result;
    if (fields.has("fluid_force_x")) result.x = fields["fluid_force_x"](0, 0);
    if (fields.has("fluid_force_y")) result.y = fields["fluid_force_y"](0, 0);
    if (fields.has("fluid_force_z")) result.z = fields["fluid_force_z"](0, 0);
    return result;
}

void advance(BearingState& state, Vector3 fluid_force,
             const SimulationConfig& cfg, double dt, double time) {
    if (cfg.motion_model != MotionModel::MOVING_BEARING) return;

    const Vector3 applied = external_load(cfg, time);
    const AxisParams px{cfg.bearing_mass, cfg.bearing_damping_x, cfg.bearing_stiffness_x,
                        fluid_force.x + applied.x};
    const AxisParams py{cfg.bearing_mass, cfg.bearing_damping_y, cfg.bearing_stiffness_y,
                        fluid_force.y + applied.y};
    const AxisParams pz{cfg.bearing_mass, cfg.bearing_damping_z, cfg.bearing_stiffness_z,
                        fluid_force.z + applied.z};

    const AxisState sx = advance_axis({state.x, state.vx}, px, cfg.motion_time_method, dt);
    const AxisState sy = advance_axis({state.y, state.vy}, py, cfg.motion_time_method, dt);
    const AxisState sz = advance_axis({state.z, state.vz}, pz, cfg.motion_time_method, dt);

    state.x = sx.x;
    state.vx = sx.v;
    state.y = sy.x;
    state.vy = sy.v;
    state.z = sz.x;
    state.vz = sz.v;
}

double equivalent_eccentricity(const BearingState& state) {
    return std::hypot(state.x, state.y);
}

double equivalent_attitude_angle_deg(const BearingState& state,
                                     const SimulationConfig& cfg) {
    const double e = equivalent_eccentricity(state);
    if (e <= 1e-30) return cfg.attitude_angle_deg;
    return std::atan2(-state.y, -state.x) * 180.0 / kPi;
}

void write_state_fields(Fields& fields, const BearingState& state,
                        const SimulationConfig& cfg, double time) {
    fill_if_present(fields, "bearing_x", state.x);
    fill_if_present(fields, "bearing_y", state.y);
    fill_if_present(fields, "bearing_z", state.z);
    fill_if_present(fields, "bearing_attitude_angle",
                    equivalent_attitude_angle_deg(state, cfg));
    const Vector3 load = external_load(cfg, time);
    fill_if_present(fields, "external_load_x", load.x);
    fill_if_present(fields, "external_load_y", load.y);
    fill_if_present(fields, "external_load_z", load.z);
}

ImplicitCoupling::ImplicitCoupling(const SimulationConfig& cfg) : cfg_(cfg) {}

BearingState ImplicitCoupling::begin_step(const BearingState& start, double dt, double time_new) {
    if (!(dt > 0.0)) throw std::runtime_error("implicit motion coupling needs dt > 0");
    start_ = start;
    dt_ = dt;
    load_ = external_load(cfg_, time_new);
    evaluations_ = 0;
    backtracks_ = 0;
    residual_ = std::numeric_limits<double>::infinity();
    correction_ = std::numeric_limits<double>::infinity();
    phase_ = Phase::BASE;
    if (jacobian_valid_ && std::abs(dt - jacobian_dt_) > 1e-12 * dt) jacobian_valid_ = false;
    BearingState trial;
    set_trial(trial, start.x + dt * start.vx, start.y + dt * start.vy, {});
    return trial;
}

void ImplicitCoupling::residual_at(const BearingState& trial, Vector3 force, double r[2]) const {
    const double m = cfg_.bearing_mass;
    const double vx = (trial.x - start_.x) / dt_;
    const double vy = (trial.y - start_.y) / dt_;
    r[0] = m * (vx - start_.vx) / dt_ + cfg_.bearing_damping_x * vx
         + cfg_.bearing_stiffness_x * trial.x - force.x - load_.x;
    r[1] = m * (vy - start_.vy) / dt_ + cfg_.bearing_damping_y * vy
         + cfg_.bearing_stiffness_y * trial.y - force.y - load_.y;
}

void ImplicitCoupling::set_trial(BearingState& trial, double x, double y, Vector3 force) const {
    trial.x = x;
    trial.y = y;
    trial.vx = (x - start_.x) / dt_;
    trial.vy = (y - start_.y) / dt_;
    const AxisState z = implicit_euler({start_.z, start_.vz},
        {cfg_.bearing_mass, cfg_.bearing_damping_z, cfg_.bearing_stiffness_z, force.z + load_.z}, dt_);
    trial.z = z.x;
    trial.vz = z.v;
}

void ImplicitCoupling::start_probe(BearingState& trial, Vector3 force) {
    phase_ = Phase::PROBE_X;
    backtracks_ = 0;
    set_trial(trial, base_x_ + cfg_.motion_coupling_fd_step * cfg_.c, base_y_, force);
}

bool ImplicitCoupling::newton_step(BearingState& trial, Vector3 force) {
    const double det = jacobian_[0][0] * jacobian_[1][1] - jacobian_[0][1] * jacobian_[1][0];
    const double norm = std::max({std::abs(jacobian_[0][0]), std::abs(jacobian_[0][1]),
                                  std::abs(jacobian_[1][0]), std::abs(jacobian_[1][1])});
    if (!std::isfinite(det) || std::abs(det) <= 1e-14 * norm * norm) return false;
    step_[0] = -(jacobian_[1][1] * base_r_[0] - jacobian_[0][1] * base_r_[1]) / det;
    step_[1] = -(-jacobian_[1][0] * base_r_[0] + jacobian_[0][0] * base_r_[1]) / det;
    const double length = std::hypot(step_[0], step_[1]);
    const double limit = cfg_.motion_coupling_step_limit * cfg_.c;
    if (length > limit) {
        step_[0] *= limit / length;
        step_[1] *= limit / length;
    }
    phase_ = Phase::NEWTON;
    set_trial(trial, base_x_ + step_[0], base_y_ + step_[1], force);
    return true;
}

bool ImplicitCoupling::update(BearingState& trial, Vector3 force) {
    ++evaluations_;
    double r[2];
    residual_at(trial, force, r);
    const double scale = std::max({std::hypot(load_.x, load_.y), std::hypot(force.x, force.y), 1e-12});
    residual_ = std::hypot(r[0], r[1]) / scale;
    if (!std::isfinite(residual_)) throw std::runtime_error("implicit motion coupling: non-finite force residual");
    if (residual_ < cfg_.motion_coupling_tolerance) return true;
    // The film force carries the noise of the field-solve tolerance; once the
    // Newton correction it implies is below the position tolerance, the step has
    // converged as far as the force can tell.
    if (jacobian_valid_) {
        const double det = jacobian_[0][0] * jacobian_[1][1] - jacobian_[0][1] * jacobian_[1][0];
        if (std::isfinite(det) && det != 0.0) {
            const double sx = (jacobian_[1][1] * r[0] - jacobian_[0][1] * r[1]) / det;
            const double sy = (-jacobian_[1][0] * r[0] + jacobian_[0][0] * r[1]) / det;
            correction_ = std::hypot(sx, sy) / cfg_.c;
            if (correction_ < cfg_.motion_coupling_position_tolerance) return true;
        }
    }
    if (evaluations_ >= cfg_.motion_coupling_max_iters)
        throw std::runtime_error("implicit motion coupling did not converge in " +
            std::to_string(evaluations_) + " film solves; relative force residual " + std::to_string(residual_));

    const double h = cfg_.motion_coupling_fd_step * cfg_.c;
    switch (phase_) {
        case Phase::BASE:
            base_x_ = trial.x;
            base_y_ = trial.y;
            base_r_[0] = r[0];
            base_r_[1] = r[1];
            if (!jacobian_valid_ || !newton_step(trial, force)) start_probe(trial, force);
            break;
        case Phase::PROBE_X:
            jacobian_[0][0] = (r[0] - base_r_[0]) / h;
            jacobian_[1][0] = (r[1] - base_r_[1]) / h;
            phase_ = Phase::PROBE_Y;
            set_trial(trial, base_x_, base_y_ + h, force);
            break;
        case Phase::PROBE_Y:
            jacobian_[0][1] = (r[0] - base_r_[0]) / h;
            jacobian_[1][1] = (r[1] - base_r_[1]) / h;
            jacobian_valid_ = true;
            jacobian_dt_ = dt_;
            ++jacobian_builds_;
            if (!newton_step(trial, force))
                throw std::runtime_error("implicit motion coupling: singular finite-difference Jacobian");
            break;
        case Phase::NEWTON:
            if (std::hypot(r[0], r[1]) < std::hypot(base_r_[0], base_r_[1])) {
                const double s[2] = {trial.x - base_x_, trial.y - base_y_};
                const double ss = s[0] * s[0] + s[1] * s[1];
                if (ss > 0.0) {
                    for (int row = 0; row < 2; ++row) {
                        const double y = r[row] - base_r_[row];
                        const double js = jacobian_[row][0] * s[0] + jacobian_[row][1] * s[1];
                        jacobian_[row][0] += (y - js) * s[0] / ss;
                        jacobian_[row][1] += (y - js) * s[1] / ss;
                    }
                }
                base_x_ = trial.x;
                base_y_ = trial.y;
                base_r_[0] = r[0];
                base_r_[1] = r[1];
                backtracks_ = 0;
                if (!newton_step(trial, force)) start_probe(trial, force);
            } else if (++backtracks_ <= 3) {
                step_[0] *= 0.5;
                step_[1] *= 0.5;
                set_trial(trial, base_x_ + step_[0], base_y_ + step_[1], force);
            } else {
                jacobian_valid_ = false;
                start_probe(trial, force);
            }
            break;
    }
    return false;
}

}  // namespace JournalMotion
