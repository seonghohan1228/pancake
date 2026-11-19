"""Part A solver checks of the load input and the implicit motion coupling.

Plain bearing (no inlet, both ends at ambient pressure, FBNS SCALAR_PCAV, constant properties),
OS-bearing geometry and speed on a coarse 180 x 48 grid. Stages:
  static     STEADY_EQUILIBRIUM at the static load W (reference eps, attitude, journal torque)
  rotating   TRANSIENT, IMPLICIT, ROTATING load of size W turning at +omega (shaft direction)
             Expect: eps equal to the static case; attitude mirrored; journal torque equal to
             the static BEARING torque T_j - W e sin(phi) (the rotating surface swaps, see README).
  agree      constant load from an off-equilibrium start, mass 10 kg, dt = 2 us:
             EXPLICIT and IMPLICIT trajectories agree within 0.5 %
  stiff      real mass 0.686 kg, dt = 20 us (above the explicit limit 2m/c):
             IMPLICIT settles at the static equilibrium; EXPLICIT is run and expected to fail

Usage: python part_a_checks.py --work DIR --solver EXE [--ranks 4 --mpi mpiexec] [--stages ...]
Writes DIR/part_a_summary.json.
"""
import argparse
import math
from pathlib import Path

import numpy as np

import common

R, L, C = 0.0105, 0.020, 30e-6
OMEGA = 2 * math.pi * 130.0
MU, RHO = 5.0e-3, 900.0
P_AMB, P_CAV = 1.0e5, 500.0
W = 1000.0
GRID = dict(n_theta_global=180, n_z_global=48)


def plain(**overrides):
    values = common.base_values(
        R=R, c=C, L=L, omega=OMEGA, mu=MU, rho=RHO, **GRID,
        cavitation_threshold="SCALAR_PCAV", p_cav=P_CAV,
        bc_z_south_val=P_AMB, bc_z_north_val=P_AMB,
        bc_z_south_temperature=343.15, bc_z_north_temperature=343.15,
        temperature_initial=343.15, temperature_reference=343.15,
        external_load_x=0.0, external_load_y=-W)
    values.update(overrides)
    return values


def static_stage(args):
    run_dir = args.work / "static"
    values = plain(e=0.7 * C, attitude_angle_deg=0.0, solution_mode="STEADY_STATE",
                   dt=1e-5, end_t=1e-5, write_interval=1e-5, motion_model="STEADY_EQUILIBRIUM",
                   bearing_initial_from_attitude="true", steady_force_tolerance=1e-5,
                   coupling_max_iters=100, coupling_tolerance=1e-8)
    record = common.run(run_dir, values, args.solver, args.ranks, args.mpi)
    if record["returncode"]:
        raise RuntimeError(f"static equilibrium failed: {common.last_log_line(run_dir, 'Equilibrium')}")
    history = list((run_dir / "results/equilibrium_history.csv").read_text().splitlines())
    last = dict(zip(history[0].split(","), map(float, history[-1].split(","))))
    fields = common.read_fields(run_dir)
    x, y = last["x_m"], last["y_m"]
    torque = float(fields["friction_torque"][0])
    e_perp = abs(x * (-W) - y * 0.0) / W          # |e x W_hat|, W along -y
    return dict(x_m=x, y_m=y, eps=math.hypot(x, y) / C,
                attitude_deg=math.degrees(math.atan2(-y, -x)),
                journal_torque_Nm=torque, bearing_torque_Nm=torque - W * e_perp,
                force_residual=last["force_residual"], h_min_m=float(np.min(fields["h"])))


def rotating_stage(args, static):
    run_dir = args.work / "rotating"
    revolutions = args.revolutions
    dt = math.radians(0.5) / OMEGA
    values = plain(e=0.0, attitude_angle_deg=0.0, solution_mode="TRANSIENT", dt=dt,
                   end_t=revolutions * 2 * math.pi / OMEGA, write_interval=2 * math.pi / OMEGA,
                   motion_model="MOVING_BEARING", bearing_mass=0.686,
                   bearing_initial_from_attitude="false",
                   bearing_initial_x=static["x_m"], bearing_initial_y=static["y_m"],
                   motion_time_method="EULER_IMPLICIT", motion_coupling="IMPLICIT",
                   motion_coupling_tolerance=1e-6, diagnostics_interval=0,
                   load_model="ROTATING", load_magnitude=W, load_rotation_start_deg=-90.0,
                   load_rotation_speed=OMEGA)
    record = common.run(run_dir, values, args.solver, args.ranks, args.mpi)
    motion = common.read_motion(run_dir)
    period = 2 * math.pi / OMEGA
    last = motion["time_s"] > motion["time_s"][-1] - period * (1 - 1e-9)
    prev = (motion["time_s"] > motion["time_s"][-1] - 2 * period) & ~last
    eps = motion["eccentricity_ratio"]
    # Angle from the load direction to the minimum-gap direction (-bearing position).
    load_angle = np.degrees(np.arctan2(motion["load_y_N"], motion["load_x_N"]))
    phi = common.wrap180(motion["attitude_angle_deg"] - load_angle)
    out = dict(returncode=record["returncode"], wall_s=record["wall_s"], steps=len(eps),
               eps_last_mean=float(eps[last].mean()), eps_last_min=float(eps[last].min()),
               eps_last_max=float(eps[last].max()), eps_prev_mean=float(eps[prev].mean()),
               attitude_from_load_last_mean_deg=float(phi[last].mean()),
               journal_torque_last_mean_Nm=float(motion["friction_torque_Nm"][last].mean()),
               coupling_solves_mean=float(motion["coupling_solves"].mean()),
               coupling_solves_max=int(motion["coupling_solves"].max()))
    out["eps_rel_diff"] = out["eps_last_mean"] / static["eps"] - 1
    out["torque_rel_diff_vs_static_bearing"] = (abs(out["journal_torque_last_mean_Nm"]) /
                                                abs(static["bearing_torque_Nm"]) - 1)
    out["torque_rel_diff_vs_static_journal"] = (abs(out["journal_torque_last_mean_Nm"]) /
                                                abs(static["journal_torque_Nm"]) - 1)
    out["static_attitude_from_load_deg"] = common.wrap180(static["attitude_deg"] - (-90.0))
    out["pass_eps"] = abs(out["eps_rel_diff"]) < 0.01
    out["pass_torque"] = abs(out["torque_rel_diff_vs_static_bearing"]) < 0.01
    out["pass_mirrored_attitude"] = abs(out["attitude_from_load_last_mean_deg"] +
                                        out["static_attitude_from_load_deg"]) < 1.0
    return out


def transient(run_dir, args, coupling, mass, dt, end_t, eps0, psi0):
    values = plain(e=eps0 * C, attitude_angle_deg=psi0, solution_mode="TRANSIENT", dt=dt, end_t=end_t,
                   write_interval=end_t, motion_model="MOVING_BEARING", bearing_mass=mass,
                   bearing_initial_from_attitude="true", motion_time_method="EULER_IMPLICIT",
                   motion_coupling=coupling, motion_coupling_tolerance=1e-7, diagnostics_interval=0)
    record = common.run(run_dir, values, args.solver, args.ranks, args.mpi)
    motion = common.read_motion(run_dir) if (run_dir / "results/motion.csv").exists() else None
    return record, motion


def agree_stage(args, static):
    dt, end_t = 2e-6, 2e-3
    psi0 = static["attitude_deg"] - 40.0
    rec_e, mot_e = transient(args.work / "agree_explicit", args, "EXPLICIT", 10.0, dt, end_t, 0.3, psi0)
    rec_i, mot_i = transient(args.work / "agree_implicit", args, "IMPLICIT", 10.0, dt, end_t, 0.3, psi0)
    n = min(len(mot_e["time_s"]), len(mot_i["time_s"]))
    dx = np.hypot(mot_e["bearing_x_m"][:n] - mot_i["bearing_x_m"][:n], mot_e["bearing_y_m"][:n] - mot_i["bearing_y_m"][:n])
    scale = max(np.hypot(static["x_m"], static["y_m"]), 1e-30)
    out = dict(explicit_returncode=rec_e["returncode"], implicit_returncode=rec_i["returncode"],
               steps_compared=int(n), steps_explicit=len(mot_e["time_s"]), steps_implicit=len(mot_i["time_s"]),
               max_position_diff_rel=float(dx.max() / scale),
               final_eps_explicit=float(mot_e["eccentricity_ratio"][n - 1]),
               final_eps_implicit=float(mot_i["eccentricity_ratio"][n - 1]),
               max_torque_diff_rel=float(np.max(np.abs(mot_e["friction_torque_Nm"][:n] - mot_i["friction_torque_Nm"][:n])) /
                                         abs(static["journal_torque_Nm"])),
               implicit_solves_mean=float(mot_i["coupling_solves"].mean()),
               wall_s_explicit=rec_e["wall_s"], wall_s_implicit=rec_i["wall_s"])
    out["pass"] = out["max_position_diff_rel"] < 0.005 and rec_e["returncode"] == 0 and rec_i["returncode"] == 0
    return out


def stiff_stage(args, static):
    dt, end_t = 2e-5, 1e-2
    psi0 = static["attitude_deg"] - 40.0
    rec_i, mot_i = transient(args.work / "stiff_implicit", args, "IMPLICIT", 0.686, dt, end_t, 0.3, psi0)
    rec_e, mot_e = transient(args.work / "stiff_explicit", args, "EXPLICIT", 0.686, dt, end_t, 0.3, psi0)
    out = dict(implicit_returncode=rec_i["returncode"], explicit_returncode=rec_e["returncode"],
               implicit_steps=len(mot_i["time_s"]) if mot_i else 0,
               explicit_steps=len(mot_e["time_s"]) if mot_e else 0,
               explicit_last_line=common.last_log_line(args.work / "stiff_explicit", "Step "),
               explicit_stop_line=common.last_log_line(args.work / "stiff_explicit", "Stopping"))
    if mot_i:
        out.update(final_eps=float(mot_i["eccentricity_ratio"][-1]),
                   final_attitude_deg=float(mot_i["attitude_angle_deg"][-1]),
                   final_eps_rel_diff=float(mot_i["eccentricity_ratio"][-1] / static["eps"] - 1),
                   implicit_solves_mean=float(mot_i["coupling_solves"].mean()),
                   implicit_solves_max=int(mot_i["coupling_solves"].max()))
        out["pass"] = rec_i["returncode"] == 0 and abs(out["final_eps_rel_diff"]) < 1e-3
    if mot_e:
        out["explicit_eps_max"] = float(mot_e["eccentricity_ratio"].max())
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--mpi", default="mpiexec")
    parser.add_argument("--revolutions", type=int, default=6)
    parser.add_argument("--stages", nargs="+", default=["static", "rotating", "agree", "stiff"])
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    summary_path = args.work / "part_a_summary.json"
    summary = {}
    if summary_path.exists():
        import json
        summary = json.loads(summary_path.read_text())
    summary["static"] = static_stage(args)
    common.save_json(summary_path, summary)
    print("static", summary["static"], flush=True)
    for stage, function in (("rotating", rotating_stage), ("agree", agree_stage), ("stiff", stiff_stage)):
        if stage in args.stages:
            summary[stage] = function(args, summary["static"])
            common.save_json(summary_path, summary)
            print(stage, summary[stage], flush=True)


if __name__ == "__main__":
    main()
