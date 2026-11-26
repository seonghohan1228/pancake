"""Set T: transient OS motion over the 130 Hz load cycle at a fixed refrigerant-rich feed.

Same procedure as the 11/19 B2 run (part_b.stage_b2): moving OS (0.686 kg), load table against
crank angle, Part B crank frame, dt = 0.25 deg; chained one-revolution runs from the steady state
at the mean load until two revolutions repeat (h_min curve within 1 %), then the last revolution
again at half the time step. Contact model off (variant N, full film only).

Contact criterion: h_min < H_CONTACT (0.5 um, a fixed criterion set by Seongho, not a roughness).
The runs are not stopped at the criterion: F2 already starts below it (steady mean load
0.489 um), so a stop would give no answer. Each revolution reports the crank angles and the
share of the cycle with h_min < H_CONTACT; those parts of the cycle count as contact, and the
full-film values there are not physical.

Fluids, steady starts and model settings come from the flooded-start campaign (--ref, read only).

Usage: python flooded_transient.py --work NEW_CAMPAIGN --ref FLOODED_CAMPAIGN --solver EXE --fluid F1
       [--ranks 12 --mpi mpirun --rankfile FILE --max-revs 6]
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import common
import flooded_start as fs
import part_b

H_CONTACT = 0.5e-6
SET, VARIANT = "A", "N"


def contact_report(motion, rev_start_time):
    """Cycle share and crank-angle intervals with h_min below H_CONTACT."""
    phase = ((motion["time_s"] - rev_start_time) / part_b.PERIOD * 360.0) % 360.0
    below = motion["h_min_m"] < H_CONTACT
    intervals, start = [], None
    for k, flag in enumerate(below):
        if flag and start is None:
            start = phase[k]
        if not flag and start is not None:
            intervals.append([float(start), float(phase[k - 1])])
            start = None
    if start is not None:
        intervals.append([float(start), float(phase[-1])])
    return dict(contact_share_of_cycle=float(below.mean()), contact_intervals_deg=intervals,
                h_min_cycle_min_m=float(motion["h_min_m"].min()), h_min_cycle_max_m=float(motion["h_min_m"].max()),
                h_min_min_at_deg=float(phase[int(np.argmin(motion["h_min_m"]))]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--fluid", required=True, choices=list(fs.FLUIDS))
    parser.add_argument("--ranks", type=int, default=12)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--rankfile")
    parser.add_argument("--max-revs", type=int, default=6)
    parser.add_argument("--step-deg", type=float, default=part_b.DT_BASE_DEG)
    args = parser.parse_args()
    args.work = args.work.resolve()
    out_dir = args.work / f"T_{args.fluid}"
    out_dir.mkdir(parents=True, exist_ok=True)
    params = json.loads((args.ref / "parameters.json").read_text())
    if not (args.work / "parameters.json").exists():
        common.save_json(args.work / "parameters.json", dict(params, source=str(args.ref / "parameters.json"),
                                                             h_contact_m=H_CONTACT))
    ref = argparse.Namespace(work=args.ref)
    steady = {k: fs.accepted_row(ref, SET, VARIANT, args.fluid, k) for k in ("min", "mean", "max")}
    row = steady["mean"]
    state = Path(row["run_dir"]) / "results" / "steady_state.dat"
    if not state.exists():
        raise SystemExit(f"steady_state.dat missing in {row['run_dir']}")
    eps, psi = row["eps"], math.radians(row["attitude_deg"])
    start = dict(state_file=str(state), x=-eps * part_b.C * math.cos(psi), y=-eps * part_b.C * math.sin(psi),
                 vx=0.0, vy=0.0)
    f = params["fluids"][args.fluid]
    fluid_params = dict(params, mu_Pa_s=f["mu_Pa_s"], rho_kg_m3=f["rho_kg_m3"])
    probe = fs.values(params, args.fluid, SET, VARIANT, 0.5, 0.0)
    keep = {k: probe[k] for k in probe if k.startswith(("cavitation_threshold", "p_cav", "dissolved_gas_initial", "bc_z_",
                                                         "inlet_", "temperature_", "omega"))}
    dt = part_b.deg_to_time(args.step_deg)

    def transient(start_state, step, write_interval=part_b.PERIOD):
        v = part_b.transient_values(fluid_params, start_state, step, part_b.PERIOD, write_interval, "TABLE", phase0=0.0)
        v.update(keep)
        return v

    summary = dict(fluid=args.fluid, mu_Pa_s=f["mu_Pa_s"], h_contact_m=H_CONTACT, step_deg=args.step_deg,
                   steady={k: dict(load_N=r["load_N"], h_min_um=r["h_min_um"], eps=r["eps"]) for k, r in steady.items()},
                   revolutions=[])
    prev, runs = None, []
    for k in range(1, args.max_revs + 1):
        run_dir = out_dir / f"rev{k:02d}"
        record = common.run(run_dir, transient(start, dt), args.solver, args.ranks, args.mpi, rankfile=args.rankfile)
        if record["returncode"]:
            raise SystemExit(f"{args.fluid} revolution {k} failed (exit {record['returncode']}); see {run_dir / 'run.log'}")
        m = common.read_motion(run_dir)
        curve = m["h_min_m"]
        entry = dict(rev=k, steps=len(curve), wall_s=record["wall_s"], start=start,
                     torque_cycle_mean_Nm=float(m["friction_torque_Nm"].mean()),
                     eps_cycle_max=float(m["eccentricity_ratio"].max()),
                     coupling_solves_mean=float(m["coupling_solves"].mean()),
                     coupling_solves_max=int(m["coupling_solves"].max()),
                     force_imbalance_max_N=float(np.max(np.hypot(m["fluid_force_x_N"] + m["load_x_N"],
                                                                 m["fluid_force_y_N"] + m["load_y_N"]))),
                     **contact_report(m, float(m["time_s"][0]) - dt))
        if prev is not None and len(prev) == len(curve):
            entry["h_min_change_max_rel"] = float(np.max(np.abs(curve / prev - 1)))
        summary["revolutions"].append(entry)
        common.save_json(out_dir / "summary.json", summary)
        print(f"[{args.fluid}] rev {k}: " + json.dumps({q: v for q, v in entry.items() if q != "start"}), flush=True)
        runs.append((run_dir, start))
        prev = curve
        start = part_b.end_state(run_dir)
        if entry.get("h_min_change_max_rel", 1.0) < 0.01:
            summary["revolutions_to_repeat"] = k
            break
    else:
        summary["revolutions_to_repeat"] = None
        common.save_json(out_dir / "summary.json", summary)
        raise SystemExit(f"{args.fluid}: cycle did not repeat within {args.max_revs} revolutions")
    last_dir, last_start = runs[-1]
    half_dir = out_dir / "last_half_dt"
    record = common.run(half_dir, transient(last_start, dt / 2), args.solver, args.ranks, args.mpi, rankfile=args.rankfile)
    m_full, m_half = common.read_motion(last_dir), common.read_motion(half_dir)
    half = dict(returncode=record["returncode"], wall_s=record["wall_s"],
                h_min_cycle_min_full_m=float(m_full["h_min_m"].min()), h_min_cycle_min_half_m=float(m_half["h_min_m"].min()),
                torque_mean_full_Nm=float(m_full["friction_torque_Nm"].mean()),
                torque_mean_half_Nm=float(m_half["friction_torque_Nm"].mean()),
                contact_share_half=contact_report(m_half, float(m_half["time_s"][0]) - dt / 2)["contact_share_of_cycle"])
    half["h_min_rel_change"] = half["h_min_cycle_min_half_m"] / half["h_min_cycle_min_full_m"] - 1
    half["torque_rel_change"] = half["torque_mean_half_Nm"] / half["torque_mean_full_Nm"] - 1
    half["pass"] = record["returncode"] == 0 and abs(half["h_min_rel_change"]) < 0.02 and abs(half["torque_rel_change"]) < 0.02
    summary["half_dt"] = half
    summary["last_revolution_dir"] = str(last_dir)
    last = summary["revolutions"][-1]
    summary["result"] = dict(
        h_min_cycle_min_um=last["h_min_cycle_min_m"] * 1e6, h_min_cycle_max_um=last["h_min_cycle_max_m"] * 1e6,
        h_min_min_at_deg=last["h_min_min_at_deg"],
        steady_h_min_um={k: v["h_min_um"] for k, v in summary["steady"].items()},
        cycle_min_vs_steady_max_load_rel=last["h_min_cycle_min_m"] * 1e6 / summary["steady"]["max"]["h_min_um"] - 1,
        contact_share_of_cycle=last["contact_share_of_cycle"], contact_intervals_deg=last["contact_intervals_deg"])
    common.save_json(out_dir / "summary.json", summary)
    print(f"[{args.fluid}] half dt: {json.dumps(half)}", flush=True)
    print(f"[{args.fluid}] result: {json.dumps(summary['result'])}", flush=True)


if __name__ == "__main__":
    main()
