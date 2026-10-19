#include "steady_equilibrium.hpp"
#include "coupling.hpp"
#include "film_thickness.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SteadyEquilibrium {
double force_residual(const Reynolds::ForceComponents& f, const SimulationConfig& cfg) {
    return std::hypot(f.fluid_x + cfg.external_load_x, f.fluid_y + cfg.external_load_y)
        / std::hypot(cfg.external_load_x, cfg.external_load_y);
}

Stats update(Fields& fields, LinearSystem& sys, Communicator& comm, const Mesh& mesh,
             const SimulationConfig& cfg, JournalMotion::BearingState& state,
             CoupledCavitation::Solver* coupled_solver) {
    Stats result;
    if ((cfg.cavitation_model == CavitationModel::COUPLED) != (coupled_solver != nullptr))
        throw std::invalid_argument("Equilibrium pressure solver does not match cavitation_model");
    auto solve_pressure = [&] {
        ++result.pressure_solves;
        if (coupled_solver) {
            result.coupled = coupled_solver->solve(fields, cfg, cfg.dt, cfg.coupled_phase_equilibrium);
            return result.coupled.converged;
        }
        result.pressure = FBNS::solve(fields, sys, comm, mesh, cfg);
        return result.pressure.converged;
    };
    state.vx = state.vy = state.vz = 0.0;
    fields["dh_dt"].fill(0.0);
    if (!solve_pressure()) throw std::runtime_error("Equilibrium base pressure solve did not converge");
    result.forces = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
    result.force_residual = force_residual(result.forces, cfg);
    if (!std::isfinite(result.force_residual)) throw std::runtime_error("Equilibrium force residual is nonfinite");
    if (result.force_residual < cfg.steady_force_tolerance) return result;

    const auto base_state = state;
    const auto base_force = result.forces;
    const auto base_pressure = result.pressure;
    const auto base_coupled = result.coupled;
    Coupling::Snapshot base;
    base.capture(fields);
    const std::vector<double> base_h = fields["h"].data;
    auto restore = [&] {
        base.restore(fields); state = base_state;
        result.pressure = base_pressure; result.coupled = base_coupled; result.forces = base_force;
    };
    auto solve_position = [&](const JournalMotion::BearingState& trial) {
        restore();
        const double dz = 0.5*cfg.L;
        const double max_offset = std::max(std::hypot(trial.x-cfg.tilt_slope_x*dz,trial.y-cfg.tilt_slope_y*dz),
            std::hypot(trial.x+cfg.tilt_slope_x*dz,trial.y+cfg.tilt_slope_y*dz));
        if (!std::isfinite(max_offset) || cfg.c-max_offset <= cfg.min_film_thickness) return false;
        state = trial;
        FilmThickness::compute_moving_bearing(fields["h"], fields["dh_dt"], mesh, cfg, state);
        fields["dh_dt"].fill(0.0);
        double local_min = std::numeric_limits<double>::infinity();
        int invalid = 0;
        for (int i = 0; i < mesh.n_theta_local; ++i) for (int j = 0; j < mesh.n_z_local; ++j) {
            const double h = fields["h"](i,j);
            if (!std::isfinite(h)) invalid = 1;
            local_min = std::min(local_min, h);
        }
        double global_min = 0.0;
        int any_invalid = 0;
        MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&invalid, &any_invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_invalid || global_min <= cfg.min_film_thickness) return false;
        // Frozen dimensionless viscosity profile: Dowson integrals retain their h powers.
        if (fields.has("cross_film_f0")) {
            for (std::size_t k = 0; k < fields["h"].data.size(); ++k) {
                if (base_h[k] <= 0.0) continue; // ghost geometry is exchanged below
                const double ratio = fields["h"].data[k] / base_h[k];
                fields["cross_film_f0"].data[k] *= ratio;
                fields["cross_film_f1"].data[k] *= ratio * ratio;
                fields["cross_film_f2"].data[k] *= ratio * ratio * ratio;
                fields["cross_film_couette_h"].data[k] *= ratio;
            }
        }
        comm.update_ghosts(fields);
        try {
            if (!solve_pressure()) return false;
        } catch (const std::exception&) {
            return false;
        }
        result.forces = Reynolds::calculate_macroscopic_properties(fields, mesh, cfg);
        return std::isfinite(force_residual(result.forces, cfg));
    };
    double jacobian[2][2]{};
    for (int axis = 0; axis < 2; ++axis) {
        double delta = cfg.steady_position_fd_step * cfg.c;
        auto trial = base_state;
        (axis == 0 ? trial.x : trial.y) += delta;
        if (!solve_position(trial)) {
            delta = -delta;
            trial = base_state;
            (axis == 0 ? trial.x : trial.y) += delta;
            if (!solve_position(trial)) {
                restore(); result.pressure = base_pressure; result.forces = base_force;
                return result;
            }
        }
        jacobian[0][axis] = (result.forces.fluid_x - base_force.fluid_x) / delta;
        jacobian[1][axis] = (result.forces.fluid_y - base_force.fluid_y) / delta;
    }
    const double a = jacobian[0][0], b = jacobian[0][1];
    const double c = jacobian[1][0], d = jacobian[1][1];
    const double det = a*d-b*c;
    const double jac_scale = std::max({std::abs(a),std::abs(b),std::abs(c),std::abs(d)});
    restore(); result.pressure = base_pressure; result.forces = base_force;
    if (!std::isfinite(det) || std::abs(det) <= 1e-12*jac_scale*jac_scale) return result;
    const double rx = base_force.fluid_x + cfg.external_load_x;
    const double ry = base_force.fluid_y + cfg.external_load_y;
    double dx = (-d*rx+b*ry)/det, dy = (c*rx-a*ry)/det;
    if (!std::isfinite(dx) || !std::isfinite(dy)) return result;
    const double scale = std::min(1.0, cfg.steady_position_step_limit*cfg.c / std::max(std::hypot(dx,dy),1e-30));
    dx *= scale; dy *= scale;
    for (double fraction = 1.0; fraction >= 1.0/128.0; fraction *= 0.5) {
        auto trial = base_state;
        trial.x += fraction*dx; trial.y += fraction*dy;
        if (!solve_position(trial)) continue;
        const double residual = force_residual(result.forces, cfg);
        if (residual < result.force_residual) {
            result.force_residual = residual;
            result.position_change = fraction*std::hypot(dx,dy)/cfg.c;
            result.accepted = true;
            JournalMotion::write_state_fields(fields, state, cfg);
            return result;
        }
    }
    restore(); result.pressure = base_pressure; result.forces = base_force;
    return result;
}
}
