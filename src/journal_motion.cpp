#include "journal_motion.hpp"

#include <cmath>
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

Vector3 external_load(const SimulationConfig& cfg) {
    return {cfg.external_load_x, cfg.external_load_y, cfg.external_load_z};
}

Vector3 fluid_force_from_fields(const Fields& fields) {
    Vector3 result;
    if (fields.has("fluid_force_x")) result.x = fields["fluid_force_x"](0, 0);
    if (fields.has("fluid_force_y")) result.y = fields["fluid_force_y"](0, 0);
    if (fields.has("fluid_force_z")) result.z = fields["fluid_force_z"](0, 0);
    return result;
}

void advance(BearingState& state, Vector3 fluid_force,
             const SimulationConfig& cfg, double dt) {
    if (cfg.motion_model != MotionModel::MOVING_BEARING) return;

    const Vector3 applied = external_load(cfg);
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
                        const SimulationConfig& cfg) {
    fill_if_present(fields, "bearing_x", state.x);
    fill_if_present(fields, "bearing_y", state.y);
    fill_if_present(fields, "bearing_z", state.z);
    fill_if_present(fields, "bearing_attitude_angle",
                    equivalent_attitude_angle_deg(state, cfg));
    fill_if_present(fields, "external_load_x", cfg.external_load_x);
    fill_if_present(fields, "external_load_y", cfg.external_load_y);
    fill_if_present(fields, "external_load_z", cfg.external_load_z);
}

}  // namespace JournalMotion
