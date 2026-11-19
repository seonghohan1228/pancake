"""Part C check: the Part B case in the frame of the orbiting scroll (non-rotating, translating
with the crank-pin centre). The OS bore is at rest (solver bearing, moving body); the pin turns at
+omega (solver shaft); the feed hole turns with the pin (inlet_rotation_speed = omega); the load
on the OS turns with the crank (TABLE with direction -90 deg + crank angle).

Map from the Part B (mirrored crank) frame theta' to this frame phi: phi = -theta' + omega t.
At t = 0 this is the cell reversal i -> N-1-i; bearing state (x, y) -> (x, -y); the velocity gains
omega e_z x delta because delta is fixed in the crank frame. The start state is the Part B state
at the start of its accepted revolution, so both runs describe the same cycle.

Comparison (prompt: within 2 %): cycle minimum h_min, and cycle mean friction torque. The solver
torque is the shear torque on its rotating surface: the OS bore in Part B, the pin here. The
like-for-like quantity is the dissipation, D = T omega - F . v (F film force on the OS, v its
velocity relative to the pin), reported as D / omega next to the raw torque.

Usage: python part_c.py --work CAMPAIGN_DIR --solver EXE --table fw68d.csv [--ranks 24 --mpi mpirun]
       [--revs 2] [--end-steps N (smoke)]
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import common
import part_b

OMEGA = part_b.OMEGA


def mirror_state_file(source, target, n_theta, n_z):
    """Reverse the theta index of every field in a PANCAKE_STEADY_GUESS file."""
    lines = Path(source).read_text().split("\n")
    header = lines[:2]
    body = lines[2:]
    out = list(header)
    index = 0
    if lines[0].split()[1] == "2":
        out.append(body[0])
        index = 1
    total = n_theta * n_z
    while index < len(body) and body[index].strip():
        name = body[index]
        values = np.array(body[index + 1:index + 1 + total], dtype=object).reshape(n_theta, n_z)
        out.append(name)
        out.extend(values[::-1].reshape(-1).tolist())
        index += 1 + total
    Path(target).write_text("\n".join(out) + "\n")


def build_table(params, path):
    rows = [line for line in Path(params["load_table_file"]).read_text().splitlines() if line and line[0].isdigit()]
    with Path(path).open("w") as stream:
        stream.write("# Part C: load on the OS in the OS frame, direction -90 deg + crank angle\n")
        stream.write("crank_angle_deg,F_OS_N,direction_deg\n")
        for row in rows:
            angle, load, _ = row.split(",")
            stream.write(f"{angle},{load},{(-90.0 + float(angle)) % 360!r}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--ranks", type=int, default=24)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--revs", type=int, default=2)
    parser.add_argument("--end-steps", type=int, default=0, help="smoke: stop after N steps")
    args = parser.parse_args()
    work = args.work.resolve()
    params = json.loads((work / "parameters.json").read_text())
    b2 = json.loads((work / "b2" / "summary.json").read_text())
    accepted = b2["revolutions"][-1]
    start_b = accepted["start"]
    out = work / "part_c"
    out.mkdir(exist_ok=True)
    table = out / "os_load_130hz_os_frame.csv"
    build_table(params, table)
    state = out / "start_state_mirrored.dat"
    mirror_state_file(start_b["state_file"], state, part_b.N_THETA, part_b.N_Z)
    x, y = start_b["x"], -start_b["y"]
    vx = start_b["vx"] - OMEGA * y
    vy = -start_b["vy"] + OMEGA * x
    start = dict(state_file=str(state), x=x, y=y, vx=vx, vy=vy)
    dt = part_b.deg_to_time(part_b.DT_BASE_DEG)
    hole_phi = (-part_b.HOLE_THETA) % 360
    summary = dict(start_b=start_b, start_c=start, revolutions=[])
    for k in range(1, args.revs + 1):
        run_dir = out / f"rev{k:02d}"
        end_t = args.end_steps * dt if args.end_steps else part_b.PERIOD
        values = part_b.transient_values(params, start, dt, end_t, end_t, "TABLE")
        values.update(load_table_file=str(table), inlet_rotation_speed=OMEGA,
                      inlet_circular=f"{hole_phi:g} {part_b.HOLE_Z!r} {part_b.HOLE_RADIUS!r} {part_b.P_FEED!r} {part_b.T_K!r}")
        record = common.run(run_dir, values, args.solver, args.ranks, args.mpi)
        if record["returncode"]:
            raise SystemExit(f"Part C revolution {k} failed; see {run_dir / 'run.log'}")
        m = common.read_motion(run_dir)
        diag = common.steady_driver.last_diagnostic(run_dir)
        import csv
        rows = list(csv.DictReader((run_dir / "results/diagnostics.csv").open(newline="")))
        mass = max(abs(float(r["liquid_residual_rel"])) for r in rows)
        dissipation = m["friction_torque_Nm"] - (m["fluid_force_x_N"] * m["bearing_vx_m_s"] +
                                                  m["fluid_force_y_N"] * m["bearing_vy_m_s"]) / OMEGA
        entry = dict(rev=k, steps=len(m["time_s"]), wall_s=record["wall_s"],
                     h_min_cycle_min_m=float(m["h_min_m"].min()), h_min_cycle_mean_m=float(m["h_min_m"].mean()),
                     torque_raw_mean_Nm=float(m["friction_torque_Nm"].mean()),
                     torque_dissipation_mean_Nm=float(dissipation.mean()),
                     eps_mean=float(m["eccentricity_ratio"].mean()),
                     coupling_solves_mean=float(m["coupling_solves"].mean()),
                     liquid_residual_rel_max=mass, last_diagnostic_converged=diag.get("converged"))
        summary["revolutions"].append(entry)
        print(f"[c] rev {k}: {entry}", flush=True)
        common.save_json(out / "summary.json", summary)
        if args.end_steps:
            return
        start = part_b.end_state(run_dir)
    ref = common.read_motion(Path(b2["last_revolution_dir"]))
    last = summary["revolutions"][-1]
    compare = dict(
        h_min_cycle_min_B_m=float(ref["h_min_m"].min()), h_min_cycle_min_C_m=last["h_min_cycle_min_m"],
        torque_mean_B_Nm=float(ref["friction_torque_Nm"].mean()),
        torque_dissipation_mean_C_Nm=last["torque_dissipation_mean_Nm"], torque_raw_mean_C_Nm=last["torque_raw_mean_Nm"],
        eps_mean_B=float(ref["eccentricity_ratio"].mean()), eps_mean_C=last["eps_mean"])
    compare["h_min_rel_diff"] = compare["h_min_cycle_min_C_m"] / compare["h_min_cycle_min_B_m"] - 1
    compare["torque_rel_diff"] = compare["torque_dissipation_mean_C_Nm"] / compare["torque_mean_B_Nm"] - 1
    compare["torque_raw_rel_diff"] = compare["torque_raw_mean_C_Nm"] / compare["torque_mean_B_Nm"] - 1
    compare["pass"] = abs(compare["h_min_rel_diff"]) < 0.02 and abs(compare["torque_rel_diff"]) < 0.02
    summary["compare_with_B"] = compare
    common.save_json(out / "summary.json", summary)
    print(f"[c] compare: {compare}", flush=True)


if __name__ == "__main__":
    main()
