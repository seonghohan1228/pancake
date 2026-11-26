"""Part B: compressor orbiting-scroll (OS) bearing at 130 Hz, GC+VC (FBNS, LOCAL_PSAT, SINGLE_PHASE).

Frame (planning/decisions.md, README): the solver runs in the frame that turns with the crank,
mirrored so the OS bore surface moves at +omega. In that frame the crank pin (with the feed
hole) is at rest and is the mesh surface; the OS bore is the solver's moving body. The load on
the OS points to theta' = 90 deg (fluid force on the OS = -W e_y'), the feed hole sits at
theta' = 90 deg (unloaded side), and only the load size changes with crank angle.

Stages (each reuses finished runs; a failed check stops the script):
  prepare  load table from the LGE sheet export, oil state, parameters.json
  b0       steady load balance (Newton on eps, psi over STATIC solves) at min, mean, max load
  b1       constant-load transient from the steady mean state: drift and time-step check
  b2       table transient, chained one-revolution runs until the cycle repeats, then the
           last revolution at half dt and once more with fields every 5 deg
  report   per-crank-angle post-processing and summary.json

Usage: python part_b.py --work DIR --solver EXE --table fw68d.csv --load-csv load_vs_crank_angle.csv
       [--ranks 24 --mpi mpirun --rankfile FILE] [--stages prepare b0 b1 b2 report]
"""
import argparse
import bisect
import csv
import json
import math
import shutil
import sys
import time
from pathlib import Path

import numpy as np

import common

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02f_lge_mixture_gas"))
import run_gas as gas  # noqa: E402

R, L, C = 0.0105, 0.020, 30e-6                  # bore 21 mm, length 20 mm, loosest radial clearance
FREQ = 130.0
OMEGA = 2 * math.pi * FREQ
PERIOD = 1.0 / FREQ
P_D = 25.0 * 98066.5                            # 25 kgf/cm2 A = shell (discharge) pressure, both ends
P_FEED = P_D + 0.3e5                            # forced-vortex estimate
T_K, T_C = 343.15, 70.0
RHO = 954.9 - 0.7284 * T_C                      # oil density law of the FW68D campaigns at 70 C
MASS = 0.686                                    # OS mass [kg]
HOLE_THETA, HOLE_Z, HOLE_RADIUS = 90.0, L / 2, 1.35e-3
LOAD_THETA = 90.0                               # load on the OS, solver frame
N_THETA, N_Z = 720, 192
P_CAV_OIL = 500.0
SPECIES = "PROPANE"                             # stand-in: SINGLE_PHASE never forms free gas
FORCE_TOL, ANGLE_TOL_DEG = 1e-4, 0.01
DT_BASE_DEG = 0.25


def deg_to_time(deg):
    return deg / 360.0 * PERIOD


# ---------------------------------------------------------------- prepare
def prepare(args):
    out = args.work / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    source = Path(args.load_csv)
    copy = out / "load_vs_crank_angle.csv"
    if not copy.exists():
        shutil.copy2(source, copy)
    rows = [r for r in csv.DictReader(copy.open(encoding="utf-8-sig")) if r["condition"] == "Sweep, 130 Hz"]
    angles = [float(r["crank_angle_deg"]) for r in rows]
    loads = [float(r["F_OS_N"]) for r in rows]
    if angles[:360] != [float(k) for k in range(360)]:
        raise ValueError("130 Hz rows must start with crank angles 0..359 deg")
    radial = [float(r["real_FR_N"]) for r in rows[:360]]
    if max(abs(v) for v in radial) != 0.0:
        raise ValueError("130 Hz rows carry a radial bearing load; the fixed-direction frame does not hold")
    table = out / "os_load_130hz.csv"
    with table.open("w", encoding="utf-8", newline="") as stream:
        stream.write("# F_OS_N of 'Sweep, 130 Hz', crank angle 0..359 deg (row 360 dropped: next cycle's 0 deg)\n")
        stream.write("# direction: load on the OS, solver frame theta' = 90 deg\n")
        stream.write("crank_angle_deg,F_OS_N,direction_deg\n")
        for a, f in zip(angles[:360], loads[:360]):
            stream.write(f"{a:g},{f!r},{LOAD_THETA:g}\n")
    cycle = np.array(loads[:360])
    viscosity, solubility = gas.property_curves(Path(args.table))
    c_feed = gas.curve_value(solubility, P_D / (1 + args.psat_margin), T_C)
    p_sat = gas.bubble_pressure(solubility, c_feed, T_C)
    nu = gas.curve_value(viscosity, c_feed, T_C, logarithmic=True) * 1e-6
    params = dict(
        R_m=R, L_m=L, clearance_m=C, frequency_Hz=FREQ, omega_rad_s=OMEGA, p_shell_Pa=P_D, p_feed_Pa=P_FEED,
        temperature_K=T_K, rho_kg_m3=RHO, nu_m2_s=nu, mu_Pa_s=nu * RHO, c_feed=c_feed,
        c_saturated_at_p_shell=gas.curve_value(solubility, P_D, T_C), p_sat_Pa=p_sat,
        psat_margin=args.psat_margin, p_shell_over_p_sat_minus_1=P_D / p_sat - 1,
        mass_kg=MASS, hole=dict(theta_deg=HOLE_THETA, z_m=HOLE_Z, radius_m=HOLE_RADIUS),
        load_direction_deg=LOAD_THETA, grid=[N_THETA, N_Z],
        cells_across_hole=[2 * HOLE_RADIUS / (2 * math.pi * R / N_THETA), 2 * HOLE_RADIUS / (L / N_Z)],
        load_min_N=float(cycle.min()), load_mean_N=float(cycle.mean()), load_max_N=float(cycle.max()),
        load_min_deg=int(cycle.argmin()), load_max_deg=int(cycle.argmax()),
        jumps=[dict(from_deg=int(k), to_deg=int((k + 1) % 360), from_N=float(cycle[k]), to_N=float(cycle[(k + 1) % 360]))
               for k in range(360) if abs(cycle[(k + 1) % 360] - cycle[k]) > 100.0],
        load_csv_sha256=common.steady_driver.sha256(copy), table_sha256=common.steady_driver.sha256(Path(args.table)),
        load_table_file=str(table), solution_table_file=str(Path(args.table).resolve()), species_stand_in=SPECIES)
    common.save_json(args.work / "parameters.json", params)
    return params


def params_of(args):
    return json.loads((args.work / "parameters.json").read_text())


# ---------------------------------------------------------------- configs
def gc_values(params, **overrides):
    values = common.base_values(
        R=R, c=C, L=L, n_theta_global=N_THETA, n_z_global=N_Z, omega=OMEGA,
        mu=params["mu_Pa_s"], rho=params["rho_kg_m3"],
        cavitation_threshold="LOCAL_PSAT", p_cav=P_CAV_OIL,
        bc_z_south_val=P_D, bc_z_north_val=P_D,
        bc_z_south_temperature=T_K, bc_z_north_temperature=T_K,
        temperature_initial=T_K, temperature_reference=T_K,
        dissolved_gas_species=SPECIES, dissolved_gas_initial=params["c_feed"], dissolved_gas_max=1.0,
        solubility_model="CURVE_TABLE", solution_table_file=params["solution_table_file"],
        dissolved_gas_henry_coeff=0, dissolved_gas_diffusivity=0,
        bc_z_south_dissolved_gas=params["c_feed"], bc_z_north_dissolved_gas=params["c_feed"],
        inlet_circular=f"{HOLE_THETA:g} {HOLE_Z!r} {HOLE_RADIUS!r} {P_FEED!r} {T_K!r}",
        inlet_species=f"0 {params['c_feed']!r} 0")
    values.update(overrides)
    return values


def static_values(params, eps, psi, overrides=None):
    """STATIC GC+VC solve at (eps, psi); `overrides` replaces keys (e.g. the VC-only model)."""
    values = gc_values(params, e=eps * C, attitude_angle_deg=psi % 360, solution_mode="STEADY_STATE",
                       dt=1e-5, end_t=1e-5, write_interval=1e-5, motion_model="STATIC",
                       bearing_initial_from_attitude="true")
    values.update(overrides or {})
    return values


def transient_values(params, start, dt, end_t, write_interval, load, phase0=0.0):
    """start: dict(state_file, x, y, vx, vy). load: float (constant) or 'TABLE'."""
    values = gc_values(params, e=0.0, attitude_angle_deg=0.0, solution_mode="TRANSIENT", dt=dt, end_t=end_t,
                       write_interval=write_interval, motion_model="MOVING_BEARING", bearing_mass=MASS,
                       bearing_initial_from_attitude="false", bearing_initial_x=start["x"],
                       bearing_initial_y=start["y"], bearing_initial_vx=start["vx"], bearing_initial_vy=start["vy"],
                       steady_initial_state_file=start["state_file"], motion_time_method="EULER_IMPLICIT",
                       motion_coupling="IMPLICIT", motion_coupling_tolerance=1e-6,
                       motion_coupling_position_tolerance=1e-6, diagnostics_interval=1,
                       output_fields="pressure,film_content,h")
    if load == "TABLE":
        values.update(load_model="TABLE", load_table_file=params["load_table_file"], load_table_abscissa="PHASE_DEG",
                      load_table_format="POLAR", load_frequency=FREQ, load_phase_offset_deg=phase0 % 360)
    else:
        angle = math.radians(LOAD_THETA)
        values.update(external_load_x=load * math.cos(angle), external_load_y=load * math.sin(angle))
    return values


def run(args, run_dir, values):
    record = common.run(run_dir, values, args.solver, args.ranks, args.mpi, rankfile=args.rankfile)
    return record


# ---------------------------------------------------------------- b0 steady
def check_static(run_dir, record):
    if record["returncode"]:
        raise ValueError(f"solver exit code {record['returncode']}")
    diag = common.steady_driver.last_diagnostic(run_dir)
    if int(diag["converged"]) != 1:
        raise ValueError("FBNS did not converge")
    mass_rel = float(diag["liquid_residual_rel"])
    if not math.isfinite(mass_rel) or abs(mass_rel) > 1e-8:
        raise ValueError(f"liquid balance residual {mass_rel}")
    f = common.read_fields(run_dir)
    force = f["F"].reshape(-1, 3)[0]
    return dict(Fx_N=float(force[0]), Fy_N=float(force[1]), torque_Nm=float(f["friction_torque"][0]),
                h_min_m=float(np.min(f["h"])), p_min_Pa=float(np.min(f["p"])), p_max_Pa=float(np.max(f["p"])),
                film_content_min=float(np.min(f["film_content"])),
                cavitated_fraction=float(np.mean(f["film_content"] < 1 - 1e-8)), mass_residual_rel=mass_rel)


def force_errors(result, load):
    """Fluid force on the OS must equal -W e_load."""
    fx, fy = result["Fx_N"], result["Fy_N"]
    magnitude = math.hypot(fx, fy) / load - 1
    angle = math.degrees(math.atan2(fy, fx)) - (LOAD_THETA + 180.0)
    return magnitude, common.wrap180(angle)


def balance(args, params, label, load, eps0, psi0, overrides=None, root=None):
    """Newton on (eps, psi) over STATIC solves until the film force balances the load.
    Runs go to (root or work/b0)/label."""
    work = (root or args.work / "b0") / label
    count = [0]
    history = []

    def solve(eps, psi, tag):
        if not 0 < eps < 1:
            raise ValueError(f"inadmissible eps {eps}")
        count[0] += 1
        run_dir = work / f"{count[0]:03d}_{tag}"
        record = run(args, run_dir, static_values(params, eps, psi, overrides))
        out = dict(run_dir=str(run_dir), eps=eps, psi_deg=psi % 360, wall_s=record["wall_s"])
        try:
            out.update(check_static(run_dir, record), valid=True)
            out["mag_err"], out["ang_err_deg"] = force_errors(out, load)
        except Exception as error:
            out.update(valid=False, error=str(error))
        history.append(out)
        print(f"[b0 {label}] {run_dir.name} eps={eps:.6f} psi={psi % 360:.4f} " +
              (f"dF={out['mag_err']:+.2e} dang={out['ang_err_deg']:+.4f} hmin={out['h_min_m'] * 1e6:.4f}um "
               f"{out['wall_s']:.0f}s" if out["valid"] else f"FAILED {out['error']}"), flush=True)
        return out

    def residual(rec):
        return np.array([rec["mag_err"], math.radians(rec["ang_err_deg"])])

    x = np.array([eps0, psi0])
    current = solve(*x, "newton0")
    if not current["valid"]:
        raise ValueError(f"{label}: initial solve failed: {current['error']}")
    jacobian = None
    for iteration in range(1, 31):
        r = residual(current)
        if abs(current["mag_err"]) < FORCE_TOL and abs(current["ang_err_deg"]) < ANGLE_TOL_DEG:
            return current, history
        if jacobian is None:
            cols = []
            for index, step in ((0, 2e-4), (1, 0.05)):
                probe = x.copy()
                probe[index] += step
                rec = solve(*probe, f"fd{iteration}_{index}")
                if not rec["valid"]:
                    raise ValueError(f"{label}: FD solve failed: {rec['error']}")
                cols.append((residual(rec) - r) / step)
            jacobian = np.column_stack(cols)
        delta = -np.linalg.solve(jacobian, r)
        delta *= min(1.0, 0.01 / abs(delta[0]) if delta[0] else 1.0, 5.0 / abs(delta[1]) if delta[1] else 1.0)
        for _ in range(5):
            candidate = solve(*(x + delta), f"newton{iteration}")
            if candidate["valid"] and np.linalg.norm(residual(candidate)) < np.linalg.norm(r):
                break
            delta *= 0.5
        else:
            jacobian = None
            continue
        dr = residual(candidate) - r
        jacobian = jacobian + np.outer(dr - jacobian @ delta, delta) / float(delta @ delta)
        x, current = x + delta, candidate
    raise ValueError(f"{label}: load balance did not converge")


def stage_b0(args):
    params = params_of(args)
    loads = dict(min=params["load_min_N"], mean=params["load_mean_N"], max=params["load_max_N"])
    results = {}
    seed = (args.eps0, args.psi0)
    for label in ("mean", "max", "min"):
        accepted, history = balance(args, params, label, loads[label], *seed)
        run_dir = Path(accepted["run_dir"])
        state = run_dir / "results" / "steady_state.dat"
        if not state.exists():
            raise ValueError(f"{label}: steady_state.dat missing in {run_dir}")
        fields = common.read_fields(run_dir)
        results[label] = dict(load_N=loads[label], accepted=accepted, solves=len(history),
                              state_file=str(state), cavity=cavity_report(fields),
                              h_min_um=accepted["h_min_m"] * 1e6,
                              effective_p_cav_Pa=accepted["p_min_Pa"],
                              p_shell_over_p_cav_minus_1=P_D / accepted["p_min_Pa"] - 1)
        common.save_json(args.work / "b0" / f"{label}.json", dict(result=results[label], history=history))
        seed = (accepted["eps"], accepted["psi_deg"])
        print(f"[b0 {label}] accepted eps={accepted['eps']:.6f} psi={accepted['psi_deg']:.4f} "
              f"hmin={accepted['h_min_m'] * 1e6:.4f}um torque={accepted['torque_Nm']:.5f}Nm", flush=True)
    common.save_json(args.work / "b0" / "summary.json", results)
    return results


# ---------------------------------------------------------------- shared transient helpers
def steady_start(args, label="mean"):
    res = json.loads((args.work / "b0" / "summary.json").read_text())[label]
    eps, psi = res["accepted"]["eps"], math.radians(res["accepted"]["psi_deg"])
    return dict(state_file=res["state_file"], x=-eps * C * math.cos(psi), y=-eps * C * math.sin(psi), vx=0.0, vy=0.0)


def end_state(run_dir):
    motion = common.read_motion(run_dir)
    state_file = Path(run_dir) / "results" / "final_state.dat"
    if not state_file.exists():
        raise ValueError(f"{run_dir}: final_state.dat missing (run failed?)")
    return dict(state_file=str(state_file), x=float(motion["bearing_x_m"][-1]), y=float(motion["bearing_y_m"][-1]),
                vx=float(motion["bearing_vx_m_s"][-1]), vy=float(motion["bearing_vy_m_s"][-1]),
                time=float(motion["time_s"][-1]))


def attitude_from_load(motion):
    return common.wrap180(motion["attitude_angle_deg"] - (LOAD_THETA + 180.0))


# ---------------------------------------------------------------- b1 constant-load transient
def stage_b1(args):
    params = params_of(args)
    start = steady_start(args)
    load = params["load_mean_N"]
    out = {}
    for step_deg in args.b1_steps:
        run_dir = args.work / "b1" / f"dt_{step_deg:g}deg"
        record = run(args, run_dir, transient_values(params, start, deg_to_time(step_deg), PERIOD, PERIOD, load))
        entry = dict(step_deg=step_deg, dt_s=deg_to_time(step_deg), returncode=record["returncode"], wall_s=record["wall_s"])
        if (run_dir / "results/motion.csv").exists():
            m = common.read_motion(run_dir)
            eps, phi = m["eccentricity_ratio"], attitude_from_load(m)
            eps0 = math.hypot(start["x"], start["y"]) / C
            phi0 = common.wrap180(math.degrees(math.atan2(-start["y"], -start["x"])) - (LOAD_THETA + 180.0))
            entry.update(steps=len(eps), eps_start=eps0, eps_end=float(eps[-1]),
                         eps_drift_max_rel=float(np.max(np.abs(eps / eps0 - 1))),
                         attitude_drift_max_rel=float(np.max(np.abs(phi - phi0)) / abs(phi0)),
                         attitude_drift_max_deg=float(np.max(np.abs(phi - phi0))),
                         coupling_solves_mean=float(m["coupling_solves"].mean()),
                         coupling_solves_max=int(m["coupling_solves"].max()),
                         force_imbalance_max_rel=float(np.max(np.hypot(m["fluid_force_x_N"] + m["load_x_N"],
                                                                       m["fluid_force_y_N"] + m["load_y_N"])) / load))
            entry["stable"] = record["returncode"] == 0 and len(eps) == round(PERIOD / entry["dt_s"])
            entry["pass_drift"] = entry["stable"] and entry["eps_drift_max_rel"] < 0.01 and entry["attitude_drift_max_rel"] < 0.01
        else:
            entry.update(stable=False, pass_drift=False)
        out[f"{step_deg:g}"] = entry
        common.save_json(args.work / "b1" / "summary.json", out)
        print(f"[b1] {step_deg:g} deg: {entry}", flush=True)
    stable = [e["step_deg"] for e in out.values() if e["stable"] and e["pass_drift"]]
    out["largest_stable_step_deg"] = max(stable) if stable else None
    common.save_json(args.work / "b1" / "summary.json", out)
    if f"{DT_BASE_DEG:g}" not in out or not out[f"{DT_BASE_DEG:g}"]["pass_drift"]:
        raise SystemExit("B1 check failed at the base time step; stopping before B2")
    return out


# ---------------------------------------------------------------- b2 table transient
def h_min_by_angle(motion, rev_start_time):
    phase = (motion["time_s"] - rev_start_time) / PERIOD * 360.0
    return phase, motion["h_min_m"]


def stage_b2(args):
    params = params_of(args)
    step_deg = args.step_deg
    dt = deg_to_time(step_deg)
    start = steady_start(args)
    revs = []
    prev_curve = None
    out_dir = args.work / "b2"
    summary = dict(step_deg=step_deg, revolutions=[])
    for k in range(1, args.max_revs + 1):
        run_dir = out_dir / f"rev{k:02d}"
        record = run(args, run_dir, transient_values(params, start, dt, PERIOD, PERIOD, "TABLE", phase0=0.0))
        if record["returncode"]:
            raise SystemExit(f"B2 revolution {k} failed (exit {record['returncode']}); see {run_dir / 'run.log'}")
        m = common.read_motion(run_dir)
        curve = m["h_min_m"]
        entry = dict(rev=k, steps=len(curve), h_min_cycle_min_m=float(curve.min()),
                     torque_cycle_mean_Nm=float(m["friction_torque_Nm"].mean()),
                     coupling_solves_mean=float(m["coupling_solves"].mean()), wall_s=record["wall_s"],
                     start=start)
        if prev_curve is not None and len(prev_curve) == len(curve):
            entry["h_min_change_max_rel"] = float(np.max(np.abs(curve / prev_curve - 1)))
        summary["revolutions"].append(entry)
        common.save_json(out_dir / "summary.json", summary)
        print(f"[b2] rev {k}: {entry}", flush=True)
        revs.append((run_dir, start))
        prev_curve = curve
        start = end_state(run_dir)
        if entry.get("h_min_change_max_rel", 1.0) < 0.01:
            summary["revolutions_to_repeat"] = k
            break
    else:
        raise SystemExit("B2: cycle did not repeat within max revolutions")
    last_dir, last_start = revs[-1]
    common.save_json(out_dir / "summary.json", summary)
    # B2b: last revolution at half dt from the same start state.
    half_dir = out_dir / "last_half_dt"
    record = run(args, half_dir, transient_values(params, last_start, dt / 2, PERIOD, PERIOD, "TABLE"))
    m_full, m_half = common.read_motion(last_dir), common.read_motion(half_dir)
    b2b = dict(returncode=record["returncode"],
               h_min_cycle_min_full_m=float(m_full["h_min_m"].min()), h_min_cycle_min_half_m=float(m_half["h_min_m"].min()),
               torque_mean_full_Nm=float(m_full["friction_torque_Nm"].mean()),
               torque_mean_half_Nm=float(m_half["friction_torque_Nm"].mean()))
    b2b["h_min_rel_change"] = b2b["h_min_cycle_min_half_m"] / b2b["h_min_cycle_min_full_m"] - 1
    b2b["torque_rel_change"] = b2b["torque_mean_half_Nm"] / b2b["torque_mean_full_Nm"] - 1
    b2b["pass"] = record["returncode"] == 0 and abs(b2b["h_min_rel_change"]) < 0.02 and abs(b2b["torque_rel_change"]) < 0.02
    summary["b2b_half_dt"] = b2b
    # Last revolution again with fields every 5 deg (deterministic repeat of last_dir).
    frames_dir = out_dir / "last_fields"
    values = transient_values(params, last_start, dt, PERIOD, deg_to_time(5.0), "TABLE")
    values["output_fields"] = "pressure,film_content,h,inlet_indicator"
    record = run(args, frames_dir, values)
    m_frames = common.read_motion(frames_dir)
    summary["last_fields_repeat_max_h_min_diff_rel"] = float(np.max(np.abs(m_frames["h_min_m"] / m_full["h_min_m"] - 1)))
    summary["last_revolution_dir"] = str(last_dir)
    summary["last_fields_dir"] = str(frames_dir)
    common.save_json(out_dir / "summary.json", summary)
    print(f"[b2] half dt: {b2b}", flush=True)
    if not b2b["pass"]:
        raise SystemExit("B2b half-dt check failed")
    return summary


# ---------------------------------------------------------------- report
def grid_of(fields):
    nt, nz = len(np.unique(fields["theta_rad"])), len(np.unique(fields["z_m"]))
    order = np.lexsort((fields["theta_rad"], fields["z_m"]))
    shape = lambda key: fields[key][order].reshape(nz, nt)  # noqa: E731
    theta = np.degrees(shape("theta_rad")[0])
    z = shape("z_m")[:, 0]
    return nt, nz, theta, z, shape


def cavity_report(fields):
    """Bubble region (void = 1 - film content >= 0.5): angular extent measured from the hole and
    from the thinnest film, area fraction, and the void around the hole."""
    nt, nz, theta, z, shape = grid_of(fields)
    void = 1 - shape("film_content")
    h = shape("h")
    theta_hmin = float(theta[np.argmin(h.min(axis=0))])
    columns = (void >= 0.5).any(axis=0)
    dth, dz = np.radians(theta[1] - theta[0]) * R, z[1] - z[0]
    tt, zz = np.meshgrid(np.radians(theta) * R, z)
    s_hole = np.radians(HOLE_THETA) * R
    ds = np.minimum(np.abs(tt - s_hole), 2 * math.pi * R - np.abs(tt - s_hole))
    hole = np.hypot(ds, zz - HOLE_Z) <= HOLE_RADIUS
    ring = (np.hypot(ds, zz - HOLE_Z) <= HOLE_RADIUS + 2 * max(dth, dz)) & ~hole
    bubble = void >= 0.5
    gap = (float(np.min(np.hypot(ds, zz - HOLE_Z)[bubble])) - HOLE_RADIUS) * 1e3 if bubble.any() else None
    out = dict(theta_hmin_deg=theta_hmin, h_min_m=float(h.min()), bubble_area_fraction=float(np.mean(bubble)),
               void_max=float(void.max()), void_hole_mean=float(void[hole].mean()),
               void_ring_mean=float(void[ring].mean()), void_ring_max=float(void[ring].max()),
               bubble_touches_hole_ring=bool((void[ring] >= 0.5).any()),
               bubble_gap_to_hole_edge_mm=gap, p_max_Pa=float(shape("p").max()))
    if not columns.any():
        return out
    if columns.all():
        out.update(bubble_start_deg=None, bubble_end_deg=None, note="bubble region spans the full circumference")
        return out
    shift = int(np.argmin(columns))
    rolled = np.roll(columns, -shift)
    runs, start = [], None
    for index, flag in enumerate(np.append(rolled, False)):
        if flag and start is None:
            start = index
        elif not flag and start is not None:
            runs.append((start, index - 1))
            start = None
    first, last = max(runs, key=lambda item: item[1] - item[0])
    step = 360.0 / nt
    start_deg = float((theta[(first + shift) % nt] - step / 2) % 360)
    end_deg = float((theta[(last + shift) % nt] + step / 2) % 360)
    out.update(bubble_start_deg=start_deg, bubble_end_deg=end_deg, bubble_regions=len(runs),
               bubble_start_from_hole_deg=common.wrap180(start_deg - HOLE_THETA),
               bubble_end_from_hole_deg=common.wrap180(end_deg - HOLE_THETA),
               bubble_start_from_hmin_deg=(start_deg - theta_hmin) % 360,
               bubble_end_from_hmin_deg=(end_deg - theta_hmin) % 360)
    return out


def frame_list(run_dir):
    import xml.etree.ElementTree as ET
    pvd = Path(run_dir) / "results/flat/results.pvd"
    return [(float(d.get("timestep")), pvd.parent / d.get("file")) for d in ET.parse(pvd).getroot().iter("DataSet")]


def read_frame(pvts):
    import xml.etree.ElementTree as ET
    arrays = {}
    for piece in ET.parse(pvts).getroot().iter("Piece"):
        source = pvts.parent / piece.get("Source")
        cell = next(ET.parse(source).getroot().iter("CellData"))
        for array in cell.findall("DataArray"):
            arrays.setdefault(array.get("Name"), []).append(np.fromstring(array.text or "", sep=" "))
    return {key: np.concatenate(values) for key, values in arrays.items()}


def stage_report(args):
    params = params_of(args)
    b0 = json.loads((args.work / "b0" / "summary.json").read_text())
    b1 = json.loads((args.work / "b1" / "summary.json").read_text())
    b2 = json.loads((args.work / "b2" / "summary.json").read_text())
    frames_dir = Path(b2["last_fields_dir"])
    motion = common.read_motion(frames_dir)
    diag = list(csv.DictReader((frames_dir / "results/diagnostics.csv").open(newline="")))
    hole_flow = np.array([float(r["inlet_mass_source"]) for r in diag])
    per_angle = []
    for t, pvts in frame_list(frames_dir):
        if t <= 0:
            continue
        fields = read_frame(pvts)
        rep = cavity_report(fields)
        index = int(np.argmin(np.abs(motion["time_s"] - t)))
        rep.update(time_s=t, crank_angle_deg=round(t / PERIOD * 360.0, 6), load_N=float(np.hypot(
            motion["load_x_N"][index], motion["load_y_N"][index])), hole_mass_flow_kg_s=float(hole_flow[index]),
            eps=float(motion["eccentricity_ratio"][index]), torque_Nm=float(motion["friction_torque_Nm"][index]))
        per_angle.append(rep)
    with (args.work / "per_crank_angle.csv").open("w", newline="") as stream:
        keys = sorted({k for row in per_angle for k in row})
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(per_angle)
    m = common.read_motion(Path(b2["last_revolution_dir"]))
    imbalance = np.hypot(m["fluid_force_x_N"] + m["load_x_N"], m["fluid_force_y_N"] + m["load_y_N"])
    load = np.hypot(m["load_x_N"], m["load_y_N"])
    reaches = [r["crank_angle_deg"] for r in per_angle if r["bubble_touches_hole_ring"]]

    def stats(values):
        return dict(min=float(np.min(values)), mean=float(np.mean(values)), max=float(np.max(values)))

    summary = dict(
        parameters=params, steady=b0, b1=b1, b2=b2,
        cycle=dict(h_min_m=stats(m["h_min_m"]), torque_Nm=stats(m["friction_torque_Nm"]), p_max_Pa=stats(m["p_max_Pa"]),
                   cavitated_fraction=stats(m["cavitated_fraction"]), eps=stats(m["eccentricity_ratio"]),
                   friction_power_W=stats(m["friction_torque_Nm"] * OMEGA),
                   hole_mass_flow_kg_s=stats(hole_flow)),
        b2c_force_follows_load=dict(max_abs_N=float(imbalance.max()), max_rel=float(np.max(imbalance / load)),
                                    at_crank_deg=float(((m["time_s"][np.argmax(imbalance)]) / PERIOD * 360.0) % 360)),
        bubble_reaches_hole=bool(reaches), bubble_reaches_hole_at_crank_deg=reaches,
        bubble_gap_to_hole_edge_mm=stats([r["bubble_gap_to_hole_edge_mm"] for r in per_angle
                                          if r["bubble_gap_to_hole_edge_mm"] is not None]),
        h_min_below_1um=bool(np.min(m["h_min_m"]) < 1e-6))
    common.save_json(args.work / "summary.json", summary)
    print(json.dumps(summary["cycle"], indent=1))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--table", required=True, help="fw68d.csv")
    parser.add_argument("--load-csv", help="load_vs_crank_angle.csv (prepare stage)")
    parser.add_argument("--ranks", type=int, default=24)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--rankfile")
    parser.add_argument("--psat-margin", type=float, default=0.0, help="p_shell / p_sat - 1 of the feed oil")
    parser.add_argument("--eps0", type=float, default=0.93)
    parser.add_argument("--psi0", type=float, default=300.0)
    parser.add_argument("--b1-steps", type=float, nargs="+", default=[1.0, 0.5, 0.25])
    parser.add_argument("--step-deg", type=float, default=DT_BASE_DEG)
    parser.add_argument("--max-revs", type=int, default=15)
    parser.add_argument("--stages", nargs="+", default=["prepare", "b0", "b1", "b2", "report"])
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    stages = dict(prepare=prepare, b0=stage_b0, b1=stage_b1, b2=stage_b2, report=stage_report)
    for name in args.stages:
        started = time.monotonic()
        stages[name](args)
        print(f"[{name}] done in {time.monotonic() - started:.0f} s", flush=True)


if __name__ == "__main__":
    main()
