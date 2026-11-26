#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "config.hpp"
#include "field.hpp"

namespace JournalMotion {

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct BearingState {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
};

/// Build the initial moving-bearing state. When configured from attitude, the
/// state is the bearing-center displacement that reproduces the static film gap.
BearingState initial_state(const SimulationConfig& cfg);

/// Read a load table (cfg.load_table_file) into Cartesian rows, applying
/// load_table_format, load_table_scale and load_table_angle_offset_deg.
/// Format: one row per line, "abscissa, a, b" separated by commas or blanks;
/// '#' starts a comment and non-numeric lines (headers) are skipped. The
/// abscissa must increase strictly. PHASE_DEG rows must lie in [0, 360) and are
/// wrapped periodically (last row -> first row + 360). One row is a constant load.
/// Returns an empty vector and sets `error` on failure.
std::vector<LoadTableRow> read_load_table(const std::filesystem::path& path,
                                          const SimulationConfig& cfg, std::string& error);

/// Resolve cfg.load_table_file (config directory first, then CWD) and fill
/// cfg.load_table. No-op unless load_model = TABLE. Appends to `errors`.
void load_config_load_table(SimulationConfig& cfg, std::vector<std::string>& errors);

/// Phase of the load cycle at time t [deg, in [0, 360)] for a PHASE_DEG table.
/// NaN for every other load option.
double load_phase_deg(const SimulationConfig& cfg, double time);

/// External load on the moving bearing at time t (global frame). CONSTANT
/// ignores t; ROTATING and TABLE follow the load_* keys. external_load_z is
/// always added.
Vector3 external_load(const SimulationConfig& cfg, double time = 0.0);

/// Read the integrated fluid force field used to move the bearing.
Vector3 fluid_force_from_fields(const Fields& fields);

/// Advance the bearing state by one timestep using cfg.motion_time_method with a
/// given (lagged) fluid force; the external load is evaluated at `time`, the end
/// of the step. This is the EXPLICIT coupling.
void advance(BearingState& state, Vector3 fluid_force,
             const SimulationConfig& cfg, double dt, double time = 0.0);

/// Equivalent shaft eccentricity magnitude for the current moving-bearing state.
double equivalent_eccentricity(const BearingState& state);

/// Equivalent attitude angle of the minimum gap, in degrees.
double equivalent_attitude_angle_deg(const BearingState& state,
                                     const SimulationConfig& cfg);

/// Fill uniform output fields for bearing state and the external load at `time`.
void write_state_fields(Fields& fields, const BearingState& state,
                        const SimulationConfig& cfg, double time = 0.0);

/// Implicit (backward-Euler) coupling of the bearing motion to the film force.
///
/// Within a step from t_n to t_{n+1} = t_n + dt the in-plane position x solves
///   R(x) = m [(x - x_n)/dt - v_n]/dt + C (x - x_n)/dt + K x - F(x) - W(t_{n+1}) = 0,
/// with v = (x - x_n)/dt and F(x) the film force of a field solve at the gap of
/// x (the film storage term uses h(x) - h_n, so the squeeze force is implicit too).
/// R is driven to zero by a quasi-Newton iteration: a finite-difference Jacobian
/// (two probe solves), Broyden rank-one updates, step limiting and backtracking.
/// The Jacobian is kept across steps while dt is unchanged. The axial axis is
/// advanced by backward Euler with the latest axial film force (no Jacobian).
///
/// Usage per step: trial = begin_step(state_n, dt, t_new); then repeat
/// { solve film at trial; if (update(trial, F)) break; restore the step-start
/// fields } -- update() writes the next trial into its argument.
class ImplicitCoupling {
public:
    explicit ImplicitCoupling(const SimulationConfig& cfg);

    /// Begin a step from the accepted state at t_n. Returns the first trial
    /// state (constant-velocity predictor).
    BearingState begin_step(const BearingState& start, double dt, double time_new);

    /// Record the film force at `trial`. Returns true when |R| / scale is below
    /// motion_coupling_tolerance, or when the Newton correction |J^-1 R| is below
    /// motion_coupling_position_tolerance * c (trial is then the accepted state). Otherwise
    /// writes the next trial into `trial` and returns false. Throws
    /// std::runtime_error when motion_coupling_max_iters solves are spent.
    bool update(BearingState& trial, Vector3 fluid_force);

    /// Film solves spent in the current step.
    int evaluations() const { return evaluations_; }
    /// Last relative force residual |R| / max(|W|, |F|).
    double residual() const { return residual_; }
    /// Last Newton correction |J^-1 R| / clearance (infinite before a Jacobian exists).
    double correction() const { return correction_; }
    /// Jacobian builds (finite-difference probes) since construction.
    int jacobian_builds() const { return jacobian_builds_; }

private:
    enum class Phase { BASE, PROBE_X, PROBE_Y, NEWTON };

    void residual_at(const BearingState& trial, Vector3 force, double r[2]) const;
    void set_trial(BearingState& trial, double x, double y, Vector3 force) const;
    bool newton_step(BearingState& trial, Vector3 force);
    void start_probe(BearingState& trial, Vector3 force);

    const SimulationConfig& cfg_;
    BearingState start_;
    Vector3 load_;
    double dt_ = 0.0;
    double jacobian_dt_ = 0.0;
    bool jacobian_valid_ = false;
    double jacobian_[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double base_x_ = 0.0, base_y_ = 0.0;
    double base_r_[2] = {0.0, 0.0};
    double step_[2] = {0.0, 0.0};
    int backtracks_ = 0;
    Phase phase_ = Phase::BASE;
    int evaluations_ = 0;
    int jacobian_builds_ = 0;
    double residual_ = 0.0;
    double correction_ = 0.0;
};

}  // namespace JournalMotion
