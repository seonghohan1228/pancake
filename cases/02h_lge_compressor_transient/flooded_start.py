"""Flooded start: the OS bearing with a refrigerant-rich feed, steady sweep in concentration.

Prompt and plan: LG final talk folder, 2026-11-24 (pancake prompt "flooded start sweep", run plan).
Steady STATIC solves with a Newton load balance, as B0 of part_b.py, for five fluids (F0-F4),
two speeds (Set A 130 Hz, Set B 30 Hz), three loads and three contact variants:

  N   no contact model (as B0)
  C2  GREENWOOD_TRIPP + PATIR_CHENG, composite roughness 0.2 um (Rq 0.141 um per surface)
  C4  the same, composite roughness 0.4 um

Model per run: FBNS, SINGLE_PHASE, isothermal at the feed temperature, constant mu and rho from
the chart; SCALAR_PCAV with p_cav = shell pressure and no dissolved gas (saturated liquid; the
chart has no p_sat for rich mixtures at 20 C). Feed hole at shell + 0.3 bar x (f / 130 Hz)^2
(forced-vortex head scales with omega^2).

Load balance: the solver's fluid_force already includes the asperity contact force (reynolds.cpp),
so the Newton iteration balances film + contact against the load without any driver addition.
Newton runs on (ln h_geo, psi), h_geo = c (1 - eps), so it stays well scaled at thin films; the
film is not allowed below H_FLOOR. A case whose film force at H_FLOOR is still below the load
has no solution: it is recorded with its last valid state.

Stages:
  repro      F0, LOCAL_PSAT, N, 130 Hz mean load: must reproduce B0 (stop otherwise)
  chain      one variant chain (--variant): Set A then Set B, F0..F4, loads mean, max, min
  checks     F0 SCALAR vs LOCAL_PSAT, shell pressure, flow factors vs contact, grid
  transient  F2, C4, 30 deg of the 130 Hz load table from the steady mean state, dt 0.25/0.125 deg
  report     sweep.csv and summary.json

Usage: python flooded_start.py --work CAMPAIGN --ref CAMPAIGN_1119 --solver EXE --stages ...
       [--variant N|C2|C4] [--ranks 8 --mpi mpirun]
"""
import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

import common
import part_b
from part_b import C, HOLE_THETA, L, LOAD_THETA, P_D, R

FLUIDS = {  # name: (refrigerant mass fraction, temperature C, expected mu mPa s)
    "F0": (None, 70.0, 5.085),   # None = B0 feed (saturated at the shell pressure)
    "F1": (0.50, 20.0, 2.62),
    "F2": (0.60, 20.0, 1.28),
    "F3": (0.80, 20.0, 0.417),
    "F4": (0.95, 20.0, 0.159),
}
VARIANTS = {"N": None, "C2": 1.41421e-07, "C4": 2.82843e-07}   # Rq per surface [m]
CONTACT = dict(asperity_contact_model="GREENWOOD_TRIPP", flow_factor_model="PATIR_CHENG",
               asperity_eta_beta_sigma=0.04, asperity_sigma_over_beta=0.001,
               contact_modulus=107237503566.48, boundary_friction_coefficient=0.1)
LOAD_ORDER = ("mean", "max", "min")
DP_FEED_130 = 0.3e5
P_SHELL_CHECK = 1.44e6
H_FLOOR = 2e-9
FIELDS = "pressure,film_content,h,inlet_indicator,fluid_force_x,fluid_force_y,fluid_force_z,friction_torque,bearing_x,bearing_y,bearing_z"
CONTACT_FIELDS = ",contact_pressure,contact_force_x,contact_force_y,contact_force_z,contact_friction_torque"


class NoSolution(Exception):
    pass


# ---------------------------------------------------------------- inputs
def prepare(args):
    ref = json.loads((args.ref / "parameters.json").read_text())
    viscosity, _ = part_b.gas.property_curves(Path(ref["solution_table_file"]) if Path(ref["solution_table_file"]).exists()
                                              else Path(args.table))
    fluids = {}
    for name, (c, t_c, expected) in FLUIDS.items():
        c = ref["c_feed"] if c is None else c
        nu = part_b.gas.curve_value(viscosity, c, t_c, logarithmic=True) * 1e-6
        rho = 954.9 - 0.7284 * t_c
        fluids[name] = dict(c=c, T_C=t_c, T_K=t_c + 273.15, nu_m2_s=nu, rho_kg_m3=rho, mu_Pa_s=nu * rho,
                            mu_expected_mPa_s=expected, mu_check_rel=nu * rho * 1e3 / expected - 1)
        if abs(fluids[name]["mu_check_rel"]) > 0.01:
            raise SystemExit(f"{name}: chart viscosity {nu * rho * 1e3:.4f} mPa s differs from {expected} by > 1 %")
    rows = [r for r in csv.DictReader((args.ref / "inputs" / "load_vs_crank_angle.csv").open(encoding="utf-8-sig"))
            if r["condition"] == "P2, 30 Hz"][:360]
    f30 = np.array([float(r["F_OS_N"]) for r in rows])
    swing = np.array([float(r["OS_load_angle_from_tangential_deg"]) for r in rows])
    sets = dict(
        A=dict(frequency_Hz=130.0, loads=dict(min=ref["load_min_N"], mean=ref["load_mean_N"], max=ref["load_max_N"])),
        B=dict(frequency_Hz=30.0, loads=dict(min=float(f30.min()), mean=float(f30.mean()), max=float(f30.max())),
               load_min_deg=int(f30.argmin()), load_max_deg=int(f30.argmax()),
               load_direction_swing_max_deg=float(np.abs(swing).max())))
    for s in sets.values():
        s["omega_rad_s"] = 2 * math.pi * s["frequency_Hz"]
        s["dp_feed_Pa"] = DP_FEED_130 * (s["frequency_Hz"] / 130.0) ** 2
    params = dict(ref, source=str(args.ref / "parameters.json"), fluids=fluids, sets=sets,
                  variants={k: dict(rq_per_surface_m=v, sigma_composite_m=None if v is None else v * math.sqrt(2))
                            for k, v in VARIANTS.items()},
                  contact_keys=CONTACT, h_floor_m=H_FLOOR, p_shell_check_Pa=P_SHELL_CHECK)
    common.save_json(args.work / "parameters.json", params)
    return params


def load_params(args):
    path = args.work / "parameters.json"
    return json.loads(path.read_text()) if path.exists() else prepare(args)


def values(params, fluid, set_key, variant, eps, psi, model="SCALAR", p_shell=P_D, grid=None, extra=None):
    f = params["fluids"][fluid]
    s = params["sets"][set_key]
    fluid_params = dict(params, mu_Pa_s=f["mu_Pa_s"], rho_kg_m3=f["rho_kg_m3"])
    t_k = f["T_K"]
    over = dict(omega=s["omega_rad_s"], bc_z_south_val=p_shell, bc_z_north_val=p_shell,
                bc_z_south_temperature=t_k, bc_z_north_temperature=t_k, temperature_initial=t_k,
                temperature_reference=t_k,
                inlet_circular=f"{HOLE_THETA:g} {part_b.HOLE_Z!r} {part_b.HOLE_RADIUS!r} {p_shell + s['dp_feed_Pa']!r} {t_k!r}",
                output_fields=FIELDS)
    if model == "SCALAR":
        over.update(cavitation_threshold="SCALAR_PCAV", p_cav=float(p_shell), dissolved_gas_initial=0.0,
                    bc_z_south_dissolved_gas=0.0, bc_z_north_dissolved_gas=0.0, inlet_species="0 0 0")
    elif fluid != "F0" or p_shell != P_D or set_key != "A":
        raise ValueError("LOCAL_PSAT is the B0 model: F0, 130 Hz, shell pressure only")
    rq = VARIANTS[variant]
    if rq is not None:
        over.update(CONTACT, roughness_rq_shaft=rq, roughness_rq_bush=rq)
        over["output_fields"] += CONTACT_FIELDS
    if grid:
        over.update(n_theta_global=grid[0], n_z_global=grid[1])
    over.update(extra or {})
    return part_b.static_values(fluid_params, eps, psi, over)


# ---------------------------------------------------------------- load balance
def balance(args, label, load, eps0, psi0, make_values, max_iters=40):
    """Newton on (u = ln h_geo, psi) over STATIC solves. Returns (accepted, history); raises
    NoSolution (with .history) when the film cannot carry the load above H_FLOOR or the
    iteration fails."""
    work = args.work / "runs" / label
    history = []
    count = [0]
    u_floor = math.log(H_FLOOR)

    def solve(u, psi, tag):
        eps = 1 - math.exp(u) / C
        count[0] += 1
        run_dir = work / f"{count[0]:03d}_{tag}"
        out = dict(run_dir=str(run_dir), eps=eps, psi_deg=psi % 360, h_geo_m=math.exp(u))
        try:
            record = common.run(run_dir, make_values(eps, psi), args.solver, args.ranks, args.mpi)
            out["wall_s"] = record["wall_s"]
            out.update(part_b.check_static(run_dir, record), valid=True)
            out["mag_err"], out["ang_err_deg"] = part_b.force_errors(out, load)
        except Exception as error:   # noqa: BLE001 - any failed solve is an invalid point
            out.update(valid=False, error=str(error)[:300])
        history.append(out)
        print(f"[{label}] {run_dir.name} eps={eps:.7f} h_geo={math.exp(u) * 1e9:.2f}nm psi={psi % 360:.4f} " +
              (f"dF={out['mag_err']:+.2e} dang={out['ang_err_deg']:+.4f}" if out["valid"] else f"FAILED {out['error'][:120]}"),
              flush=True)
        return out

    def residual(rec):
        return np.array([rec["mag_err"], math.radians(rec["ang_err_deg"])])

    def fail(reason):
        error = NoSolution(reason)
        error.history = history
        raise error

    x = np.array([math.log(C * (1 - eps0)), psi0])
    current = solve(*x, "newton0")
    if not current["valid"]:
        fail(f"initial solve failed: {current['error']}")
    jacobian = None
    for iteration in range(1, max_iters + 1):
        r = residual(current)
        if abs(current["mag_err"]) < part_b.FORCE_TOL and abs(current["ang_err_deg"]) < part_b.ANGLE_TOL_DEG:
            return current, history
        if jacobian is None:
            cols = []
            for index, step in ((0, -0.02), (1, 0.05)):
                probe = x.copy()
                probe[index] += step
                rec = solve(*probe, f"fd{iteration}_{index}")
                if not rec["valid"]:
                    fail(f"FD solve failed: {rec['error']}")
                cols.append((residual(rec) - r) / step)
            jacobian = np.column_stack(cols)
        delta = -np.linalg.solve(jacobian, r)
        delta *= min(1.0, 0.4 / abs(delta[0]) if delta[0] else 1.0, 5.0 / abs(delta[1]) if delta[1] else 1.0)
        if x[0] + delta[0] < u_floor:
            if x[0] - u_floor < 1e-6 and current["mag_err"] < 0:
                fail(f"no film solution above h_geo = {H_FLOOR * 1e9:g} nm: film + contact force is "
                     f"{(1 + current['mag_err']) * 100:.1f} % of the load there")
            delta *= (u_floor - x[0]) / delta[0]
        for _ in range(6):
            candidate = solve(*(x + delta), f"newton{iteration}")
            if candidate["valid"] and np.linalg.norm(residual(candidate)) < np.linalg.norm(r):
                break
            delta *= 0.5
        else:
            if jacobian is not None and iteration > 1:
                jacobian = None
                continue
            fail("line search failed")
        dr = residual(candidate) - r
        jacobian = jacobian + np.outer(dr - jacobian @ delta, delta) / float(delta @ delta)
        x, current = x + delta, candidate
    fail(f"load balance did not converge in {max_iters} iterations")


# ---------------------------------------------------------------- results
def details(params, accepted, fluid, set_key, variant, load):
    run_dir = Path(accepted["run_dir"])
    fields = common.read_fields(run_dir)
    diag = common.steady_driver.last_diagnostic(run_dir)
    cav = part_b.cavity_report(fields)
    f = params["fluids"][fluid]
    s = params["sets"][set_key]
    force = fields["F"].reshape(-1, 3)[0][:2]
    contact = np.array([fields["contact_force_x"][0], fields["contact_force_y"][0]]) if "contact_force_x" in fields else np.zeros(2)
    t_total = float(fields["friction_torque"][0])
    t_contact = float(fields["contact_friction_torque"][0]) if "contact_friction_torque" in fields else 0.0
    h_min = float(np.min(fields["h"]))
    sommerfeld = f["mu_Pa_s"] * s["frequency_Hz"] / (load / (2 * R * L)) * (R / C) ** 2
    row = dict(set=set_key, frequency_Hz=s["frequency_Hz"], fluid=fluid, variant=variant, load_N=load,
               c=f["c"], T_C=f["T_C"], mu_mPa_s=f["mu_Pa_s"] * 1e3, rho_kg_m3=f["rho_kg_m3"], sommerfeld=sommerfeld,
               converged=True, reason="", eps=accepted["eps"], attitude_deg=accepted["psi_deg"], h_min_um=h_min * 1e6,
               h_min_angle_from_hole_deg=common.wrap180(cav["theta_hmin_deg"] - HOLE_THETA),
               lambda_0p2=h_min / (VARIANTS["C2"] * math.sqrt(2)), lambda_0p4=h_min / (VARIANTS["C4"] * math.sqrt(2)),
               p_max_bar=float(np.max(fields["p"])) / 1e5,
               contact_p_max_bar=float(np.max(fields["contact_pressure"])) / 1e5 if "contact_pressure" in fields else 0.0,
               film_load_N=float(np.hypot(*(force - contact))), contact_load_N=float(np.hypot(*contact)),
               contact_share_pct=float(contact @ force / (force @ force)) * 100,
               torque_film_Nm=t_total - t_contact, torque_contact_Nm=t_contact, torque_total_Nm=t_total,
               friction_power_W=t_total * s["omega_rad_s"],
               cavitated_area_pct=accepted["cavitated_fraction"] * 100,
               cavity_start_from_hole_deg=cav.get("bubble_start_from_hole_deg"),
               cavity_end_from_hole_deg=cav.get("bubble_end_from_hole_deg"),
               bubble_gap_to_hole_edge_mm=cav.get("bubble_gap_to_hole_edge_mm"),
               hole_mass_flow_g_s=float(diag["inlet_mass_source"]) * 1e3,
               liquid_residual_rel=float(diag["liquid_residual_rel"]), valid=accepted["valid"],
               force_mag_err=accepted["mag_err"], force_angle_err_deg=accepted["ang_err_deg"], run_dir=str(run_dir))
    return row


def failed_row(params, history, reason, fluid, set_key, variant, load):
    valid = [h for h in history if h.get("valid")]
    row = dict(set=set_key, frequency_Hz=params["sets"][set_key]["frequency_Hz"], fluid=fluid, variant=variant,
               load_N=load, converged=False, reason=reason, solves=len(history))
    if valid:
        last = min(valid, key=lambda h: abs(h["mag_err"]) + abs(math.radians(h["ang_err_deg"])))
        row.update(last_valid_eps=last["eps"], last_valid_h_min_um=last["h_min_m"] * 1e6,
                   last_valid_force_ratio=1 + last["mag_err"], last_valid_run_dir=last["run_dir"],
                   min_h_reached_um=min(h["h_min_m"] for h in valid) * 1e6)
    return row


def result_path(args, set_key, variant, fluid, load_key):
    return args.work / "sweep" / set_key / variant / f"{fluid}_{load_key}.json"


def run_case(args, params, set_key, variant, fluid, load_key, seeds):
    path = result_path(args, set_key, variant, fluid, load_key)
    if path.exists():
        return json.loads(path.read_text())
    load = params["sets"][set_key]["loads"][load_key]
    label = f"{set_key}_{variant}_{fluid}_{load_key}"
    make = lambda eps, psi: values(params, fluid, set_key, variant, eps, psi)  # noqa: E731
    start = time.monotonic()
    result = None
    for seed in seeds:
        if seed is None:
            continue
        try:
            accepted, history = balance(args, label, load, *seed, make)
            result = dict(row=details(params, accepted, fluid, set_key, variant, load), seed=seed, history=history)
            result["row"]["solves"] = len(history)
            break
        except NoSolution as error:
            result = dict(row=failed_row(params, error.history, str(error), fluid, set_key, variant, load),
                          seed=seed, history=error.history)
            if not str(error).startswith("initial solve failed"):
                break
    result["row"]["wall_s"] = time.monotonic() - start
    path.parent.mkdir(parents=True, exist_ok=True)
    common.save_json(path, result)
    row = result["row"]
    print(f"[{label}] " + (f"eps={row['eps']:.6f} h_min={row['h_min_um']:.4f}um share={row['contact_share_pct']:.2f}% "
                           f"T={row['torque_total_Nm']:.5f}Nm" if row["converged"] else f"NOT CONVERGED: {row['reason']}"),
          flush=True)
    return result


def seed_of(result):
    if result is None:
        return None
    row = result["row"]
    if row["converged"]:
        return (row["eps"], row["attitude_deg"])
    history = [h for h in result["history"] if h.get("valid")]
    if not history:
        return None
    last = min(history, key=lambda h: abs(h["mag_err"]))
    return (last["eps"], last["psi_deg"])


def nearest_load(params, set_from, set_to, load_key):
    load = params["sets"][set_to]["loads"][load_key]
    loads = params["sets"][set_from]["loads"]
    return min(loads, key=lambda k: abs(loads[k] - load))


# ---------------------------------------------------------------- stages
def stage_repro(args, params):
    b0 = json.loads((args.ref / "b0" / "summary.json").read_text())["mean"]["accepted"]
    make = lambda eps, psi: values(params, "F0", "A", "N", eps, psi, model="LOCAL")  # noqa: E731
    accepted, history = balance(args, "R_F0_mean_localpsat", params["load_mean_N"], b0["eps"], b0["psi_deg"], make)
    row = details(params, accepted, "F0", "A", "N", params["load_mean_N"])
    target = dict(eps=0.943750, h_min_um=1.68766, torque_total_Nm=0.077335)
    exact = dict(eps=b0["eps"], h_min_um=b0["h_min_m"] * 1e6, torque_total_Nm=b0["torque_Nm"])
    rel = {k: row[k] / exact[k] - 1 for k in exact}
    rel_prompt = {k: row[k] / target[k] - 1 for k in target}
    out = dict(row=row, b0=exact, rel_diff=rel, rel_diff_prompt_rounded=rel_prompt,
               passed=all(abs(v) < 1e-4 for v in rel_prompt.values()), solves=len(history))
    common.save_json(args.work / "checks" / "reproduction.json", out)
    print("[repro]", json.dumps(dict(rel=rel, rel_prompt=rel_prompt, passed=out["passed"])), flush=True)
    if not out["passed"]:
        raise SystemExit("reproduction failed")


def stage_chain(args, params):
    variant = args.variant
    b0 = json.loads((args.ref / "b0" / "summary.json").read_text())
    done = {}
    for set_key in ("A", "B"):
        for index, fluid in enumerate(FLUIDS):
            for load_key in LOAD_ORDER:
                previous = list(FLUIDS)[index - 1] if index else None
                if set_key == "A" and fluid == "F0":
                    seeds = [(b0[load_key]["accepted"]["eps"], b0[load_key]["accepted"]["psi_deg"])]
                elif set_key == "A":
                    seeds = [seed_of(done.get(("A", previous, "mean"))) if load_key == "mean"
                             else seed_of(done.get(("A", fluid, "mean")))]
                else:   # Set B: Set A of the same fluid at the nearest load first (prompt), then Set B neighbours
                    seeds = [seed_of(done.get(("A", fluid, nearest_load(params, "A", "B", load_key)))),
                             seed_of(done.get(("B", previous, load_key))) if previous else None,
                             seed_of(done.get(("B", fluid, "mean"))) if load_key != "mean" else None]
                seeds = [s for s in seeds if s is not None] or [(0.94, 290.0)]
                done[(set_key, fluid, load_key)] = run_case(args, params, set_key, variant, fluid, load_key, seeds)


def accepted_row(args, set_key, variant, fluid, load_key):
    path = result_path(args, set_key, variant, fluid, load_key)
    return json.loads(path.read_text())["row"] if path.exists() else None


def check_case(args, params, label, fluid, set_key, variant, seed, load, **kw):
    path = args.work / "checks" / f"{label}.json"
    if path.exists():
        return json.loads(path.read_text())["row"]
    make = lambda eps, psi: values(params, fluid, set_key, variant, eps, psi, **kw)  # noqa: E731
    try:
        accepted, history = balance(args, label, load, *seed, make)
        row = details(params, accepted, fluid, set_key, variant, load)
        row["solves"] = len(history)
    except NoSolution as error:
        row = failed_row(params, error.history, str(error), fluid, set_key, variant, load)
    common.save_json(path, dict(row=row, overrides={k: str(v) for k, v in kw.items()}))
    return row


def rel(a, b, key):
    return a[key] / b[key] - 1


def stage_checks(args, params):
    out = {}
    mean_a = params["sets"]["A"]["loads"]["mean"]
    repro = json.loads((args.work / "checks" / "reproduction.json").read_text())["row"]
    f0n = accepted_row(args, "A", "N", "F0", "mean")
    out["scalar_vs_localpsat"] = {k: rel(f0n, repro, k) for k in ("eps", "h_min_um", "torque_total_Nm", "hole_mass_flow_g_s")}
    f2n = accepted_row(args, "A", "N", "F2", "mean")
    if f2n and f2n["converged"]:
        low = check_case(args, params, "shell_1p44MPa_A_N_F2_mean", "F2", "A", "N", (f2n["eps"], f2n["attitude_deg"]),
                         mean_a, p_shell=P_SHELL_CHECK)
        out["shell_pressure"] = dict(p_shell_Pa=P_SHELL_CHECK, **{k: rel(low, f2n, k) for k in ("h_min_um", "torque_total_Nm", "eps")})
        out["shell_pressure"]["passed"] = all(abs(out["shell_pressure"][k]) < 1e-4 for k in ("h_min_um", "torque_total_Nm"))
    split = {}
    for variant in ("C2", "C4"):
        full = accepted_row(args, "A", variant, "F0", "mean")
        contact_only = check_case(args, params, f"ff_none_A_{variant}_F0_mean", "F0", "A", variant,
                                  (full["eps"], full["attitude_deg"]), mean_a, extra=dict(flow_factor_model="NONE"))
        split[variant] = dict(
            total_vs_N={k: rel(full, f0n, k) for k in ("h_min_um", "torque_total_Nm")},
            contact_only_vs_N={k: rel(contact_only, f0n, k) for k in ("h_min_um", "torque_total_Nm")},
            flow_factors_part={k: rel(full, contact_only, k) for k in ("h_min_um", "torque_total_Nm")},
            contact_share_pct=full["contact_share_pct"])
    out["contact_at_thick_film"] = split
    rows = [r for p in sorted((args.work / "sweep").rglob("*.json")) for r in [json.loads(p.read_text())["row"]] if r["converged"]]
    thin = min(rows, key=lambda r: r["h_min_um"])
    fine = check_case(args, params, f"grid_1440x384_{thin['set']}_{thin['variant']}_{thin['fluid']}_{thin['load_N']:.0f}N",
                      thin["fluid"], thin["set"], thin["variant"], (thin["eps"], thin["attitude_deg"]), thin["load_N"],
                      grid=(1440, 384))
    out["grid"] = dict(case=dict((k, thin[k]) for k in ("set", "variant", "fluid", "load_N")),
                       coarse=dict((k, thin[k]) for k in ("h_min_um", "torque_total_Nm", "contact_share_pct", "eps")),
                       fine=dict((k, fine.get(k)) for k in ("h_min_um", "torque_total_Nm", "contact_share_pct", "eps", "converged")),
                       rel={k: rel(fine, thin, k) for k in ("h_min_um", "torque_total_Nm", "eps")} if fine["converged"] else None)
    common.save_json(args.work / "checks" / "summary.json", out)
    print(json.dumps(out, indent=1, default=str), flush=True)


def stage_transient(args, params):
    fluid, variant = "F2", "C4"
    row = accepted_row(args, "A", variant, fluid, "mean")
    state = Path(row["run_dir"]) / "results" / "steady_state.dat"
    if not state.exists():
        raise SystemExit(f"steady_state.dat missing in {row['run_dir']}")
    eps, psi = row["eps"], math.radians(row["attitude_deg"])
    start = dict(state_file=str(state), x=-eps * C * math.cos(psi), y=-eps * C * math.sin(psi), vx=0.0, vy=0.0)
    up = [j for j in params["jumps"] if j["to_N"] > j["from_N"]]
    phase0 = (up[0]["from_deg"] - 10) % 360 if up else 0.0
    f = params["fluids"][fluid]
    fluid_params = dict(params, mu_Pa_s=f["mu_Pa_s"], rho_kg_m3=f["rho_kg_m3"])
    probe = values(params, fluid, "A", variant, 0.5, 0.0)
    keep = {k: probe[k] for k in probe if k.startswith(("asperity", "flow_factor", "roughness", "contact", "boundary_friction",
                                                          "cavitation_threshold", "p_cav", "dissolved_gas_initial",
                                                          "bc_z_", "inlet_", "temperature_"))}
    out = dict(start_phase_deg=phase0, window_deg=30.0, steady=dict((k, row[k]) for k in ("eps", "h_min_um", "contact_share_pct",
                                                                                        "contact_load_N", "torque_total_Nm")))
    for step_deg in (0.25, 0.125):
        run_dir = args.work / "transient" / f"F2_C4_dt_{step_deg:g}deg"
        v = part_b.transient_values(fluid_params, start, part_b.deg_to_time(step_deg), part_b.deg_to_time(30.0),
                                    part_b.deg_to_time(1.0), "TABLE", phase0)
        v.update(keep, output_fields="pressure,film_content,h,fluid_force_x,fluid_force_y,friction_torque" + CONTACT_FIELDS)
        record = common.run(run_dir, v, args.solver, args.ranks, args.mpi)
        entry = dict(returncode=record["returncode"], wall_s=record["wall_s"])
        if (run_dir / "results/motion.csv").exists():
            m = common.read_motion(run_dir)
            frames = part_b.frame_list(run_dir)
            history = []
            for t, pvts in frames:
                fr = part_b.read_frame(pvts)
                history.append(dict(time_s=t, phase_deg=phase0 + t / part_b.PERIOD * 360.0,
                                    h_min_um=float(np.min(fr["h"])) * 1e6,
                                    contact_x_N=float(fr["contact_force_x"][0]), contact_y_N=float(fr["contact_force_y"][0]),
                                    contact_N=float(np.hypot(fr["contact_force_x"][0], fr["contact_force_y"][0])),
                                    torque_Nm=float(fr["friction_torque"][0]),
                                    contact_torque_Nm=float(fr["contact_friction_torque"][0])))
            entry.update(steps=len(m["time_s"]), expected_steps=round(30.0 / step_deg),
                         h_min_um_min=float(m["h_min_m"].min()) * 1e6, h_min_um_end=float(m["h_min_m"][-1]) * 1e6,
                         eps_end=float(m["eccentricity_ratio"][-1]),
                         coupling_solves_mean=float(m["coupling_solves"].mean()),
                         coupling_solves_max=int(m["coupling_solves"].max()),
                         force_imbalance_max_N=float(np.max(np.hypot(m["fluid_force_x_N"] + m["load_x_N"],
                                                                     m["fluid_force_y_N"] + m["load_y_N"]))),
                         contact_N_max=max(h["contact_N"] for h in history) if history else None,
                         frames=history)
        out[f"{step_deg:g}"] = entry
        common.save_json(args.work / "transient" / "summary.json", out)
        print(f"[transient {step_deg:g}] " + json.dumps({k: v for k, v in entry.items() if k != "frames"}), flush=True)
    a, b = out.get("0.25", {}), out.get("0.125", {})
    if "h_min_um_min" in a and "h_min_um_min" in b:
        out["dt_difference"] = dict(h_min_min_rel=a["h_min_um_min"] / b["h_min_um_min"] - 1,
                                    contact_max_rel=(a["contact_N_max"] / b["contact_N_max"] - 1) if b["contact_N_max"] else None)
    common.save_json(args.work / "transient" / "summary.json", out)


def stage_report(args, params):
    rows = [json.loads(p.read_text())["row"] for p in sorted((args.work / "sweep").rglob("*.json"))]
    order = {k: i for i, k in enumerate(FLUIDS)}
    rows.sort(key=lambda r: (r["set"], list(VARIANTS).index(r["variant"]), order[r["fluid"]], r["load_N"]))
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with (args.work / "sweep.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    checks = {}
    for name in ("reproduction", "summary"):
        p = args.work / "checks" / f"{name}.json"
        if p.exists():
            checks[name] = json.loads(p.read_text())
    tr = args.work / "transient" / "summary.json"
    summary = dict(solver=str(args.solver), solver_sha256=common.steady_driver.sha256(args.solver),
                   cases=len(rows), converged=sum(r["converged"] for r in rows),
                   not_converged=[dict((k, r.get(k)) for k in ("set", "variant", "fluid", "load_N", "reason",
                                                                "last_valid_h_min_um", "last_valid_force_ratio"))
                                  for r in rows if not r["converged"]],
                   reproduction=checks.get("reproduction"), checks=checks.get("summary"),
                   transient=json.loads(tr.read_text()) if tr.exists() else None, rows=rows)
    common.save_json(args.work / "summary.json", summary)
    print(f"report: {len(rows)} rows, {summary['converged']} converged", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--table", default=str(Path(__file__).resolve().parents[2] / "data/R410A_FW68D/fw68d.csv"))
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--variant", choices=list(VARIANTS))
    parser.add_argument("--stages", nargs="+", default=["repro"])
    args = parser.parse_args()
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    params = load_params(args)
    for stage in args.stages:
        dict(repro=stage_repro, chain=stage_chain, checks=stage_checks, transient=stage_transient,
             report=stage_report)[stage](args, params)


if __name__ == "__main__":
    main()
