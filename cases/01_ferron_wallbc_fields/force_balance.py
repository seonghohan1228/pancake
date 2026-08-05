"""Case 01-WALLBC-FIELDS: four wall treatments at the Ferron 1983 Fig. 4 point,
with the wall heat-flux fields written out.

Method, sign convention, grid, tolerances and the two-pass groove closure are
frozen from cases/01_ferron_thd_crossfilm/force_balance.py. The ONLY thing that
varies between the four groups is the thermal wall block and the output_fields
line; both are supplied here so the rest of the config cannot drift.

Groups
------
rerun_xf_adiabatic   cross-film T(y), LOCAL 1-D bush, both walls adiabatic.
rerun_xf_bush40_J40  cross-film T(y), LOCAL 1-D bush, both surfaces driven onto
                     40 C by the large-conductance limit (the case-01 mapping).
favg_walls40         NO cross-film: depth-averaged T with Newton-cooling walls at
                     40 C (journal h = 10000, bearing h = 7200, the case-01-THD
                     values). The four film-surface fields do NOT exist under
                     cross_film_thermal_model = NONE and requesting them is a
                     config error, so they are omitted here and the wall fluxes
                     are post-derived as q = h (T_wall - T_bar), positive INTO
                     the fluid, in analyse_wallbc.py.
rerun_cht_robin10    conjugate 2-D bush conduction (k = 250 W/(m K), R2 = 0.10 m,
                     16 radial cells) with a Robin outer surface at
                     h = 10 W/(m^2 K) to 313.15 K, shaft globally adiabatic.
                     NOTE: the CHT *verdict* campaign (cases/01_ferron_conjugate)
                     used the paper's h = 80; this production group uses the
                     user-selected h = 10. Both are traceable; the measured
                     difference is +0.2 K and -0.5 deg (see that case's README).

Warm starts: each group starts from the converged equilibrium of its closest
archived twin instead of the cold (0.60, -40) guess, which removes one or two
Newton steps per pass. The starting point cannot change the converged answer --
it is the same fixed point, found from closer in.

Usage: python force_balance.py [group ...]   (default: all four)
"""

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solver_io as S

CASE_DIR = Path(__file__).resolve().parent
PROGRESS = CASE_DIR / "progress.log"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text()

# Outer (Picard) coupling schedule. The defaults ARE the frozen protocol
# (200 sweeps, no under-relaxation); the environment overrides exist so the
# same driver can reproduce the relaxation study in PICARD_FIX_REPORT.md
# without forking it.
COUPLING_MAX_ITERS = int(os.environ.get("PANCAKE_COUPLING_MAX_ITERS", "200"))
COUPLING_RELAXATION = float(os.environ.get("PANCAKE_COUPLING_RELAXATION", "1.0"))
# Extra config lines appended verbatim to every generated config, ";"-separated
# (e.g. PANCAKE_COUPLING_EXTRA="coupling_relaxation_mode = ADAPTIVE;coupling_relaxation_min = 0.1").
# Empty by default, so an unset environment writes a byte-identical config.
_extra = [t.strip() for t in os.environ.get("PANCAKE_COUPLING_EXTRA", "").split(";") if t.strip()]
COUPLING_EXTRA = ("\n" + "\n".join(_extra) + "\n") if _extra else ""

CLEARANCE = 1.52e-4
N_THETA, N_Z = 720, 192
W_TARGET = 4000.0
ANGLE_TARGET = -90.0
W_TOL_REL = 1.0e-3
ANGLE_TOL_DEG = 0.05
D_EPS = 5.0e-3
D_PSI = 0.5

FIXED_WALL_HTC = 1.0e7   # case-01 mapping of a "fixed 40 C" surface
FIXED_WALL_K = 1.0e7

# Wall-flux fields are only defined when the cross-film model is on.
XF_FIELDS = ("temperature_bearing, temperature_journal, heat_flux_bearing, "
             "heat_flux_journal, cross_film_t_bush, cross_film_t_shaft, "
             "cross_film_core_rise, ")
BASE_FIELDS = ("pressure, temperature, film_content, h, mu, mass_residual, "
               "pressure_force_x, pressure_force_y, fluid_force_x, "
               "fluid_force_y, friction_torque")

COMMON_THERMAL = """# Energy Equation
temperature_model = ENERGY_EQUATION
temperature_initial = 313.15
temperature_reference = 313.15
rho_cp = 1.72e6
thermal_conductivity = 0.13
"""

XF_COMMON = """cross_film_thermal_model = POLYNOMIAL
cross_film_profile_order = 4
cross_film_quadrature_points = 6
cross_film_relaxation = 1.0
thermal_ambient_temperature = 313.15
shaft_ambient_area_ratio = 1.0
# Ignored under the conjugate closure; pinned to inert values so nothing in the
# file looks effective when it is not.
journal_wall_temperature = 313.15
bearing_wall_temperature = 313.15
journal_heat_transfer = 0
bearing_heat_transfer = 0
"""

GROUPS = {
    "rerun_xf_adiabatic": dict(
        warm=(0.608595, -37.5462),
        fields=XF_FIELDS + BASE_FIELDS,
        thermal=COMMON_THERMAL + XF_COMMON + f"""shaft_ambient_heat_transfer = 0
bush_thermal_model = LOCAL
bush_thermal_conductivity = 50
bush_outer_radius = 0.10
bush_outer_boundary = CONVECTION
bush_outer_heat_transfer = 0
""",
    ),
    "rerun_xf_bush40_J40": dict(
        warm=(0.564609, -36.1411),
        fields=XF_FIELDS + BASE_FIELDS,
        thermal=COMMON_THERMAL + XF_COMMON + f"""shaft_ambient_heat_transfer = {FIXED_WALL_HTC!r}
bush_thermal_model = LOCAL
bush_thermal_conductivity = {FIXED_WALL_K!r}
bush_outer_radius = 0.10
bush_outer_boundary = CONVECTION
bush_outer_heat_transfer = {FIXED_WALL_HTC!r}
""",
    ),
    "favg_walls40": dict(
        warm=(0.567497, -36.1187),
        fields=BASE_FIELDS,
        thermal=COMMON_THERMAL + """journal_wall_temperature = 313.15
bearing_wall_temperature = 313.15
journal_heat_transfer = 10000
bearing_heat_transfer = 7200
cross_film_thermal_model = NONE
""",
    ),
    "rerun_cht_robin10": dict(
        warm=(0.593400, -37.5036),
        fields=XF_FIELDS + BASE_FIELDS,
        thermal=COMMON_THERMAL + XF_COMMON + """shaft_ambient_heat_transfer = 0
bush_thermal_model = CONDUCTION_2D
bush_thermal_conductivity = 250
bush_outer_radius = 0.10
bush_radial_cells = 16
bush_outer_boundary = CONVECTION
bush_outer_heat_transfer = 10
""",
    ),
}

Q_SUPPLY_PLACEHOLDER = 1.0e-4

METRIC_KEYS = ("Fx", "Fy", "W", "force_angle_deg", "p_max", "h_min", "cav_frac",
               "friction_torque", "T_max", "T_min", "T_max_theta_deg", "mu_min",
               "mu_max", "q_supply", "liquid_residual", "coupling")


def log_path(group):
    return CASE_DIR / group / "fb_log.json"


def load_log(group):
    p = log_path(group)
    return json.loads(p.read_text()) if p.exists() else {"iters": []}


def save_log(group, log):
    p = log_path(group)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(log, indent=2))


def evaluate(group, log, cfgkeys, eps, psi, tag):
    key = json.dumps(cfgkeys, sort_keys=True)
    for rec in log["iters"]:
        if (rec.get("cfgkey") == key and abs(rec["eps"] - eps) < 1e-12
                and abs(rec["psi"] - psi) < 1e-12 and rec["ok"]):
            return rec
    n = len(log["iters"])
    run_dir = CASE_DIR / group / f"fb_iter_{n:02d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    spec = GROUPS[group]
    (run_dir / "config.txt").write_text(TEMPLATE.format(
        e=repr(eps * CLEARANCE), attitude=repr(psi),
        n_theta=N_THETA, n_z=N_Z,
        thermal_block=spec["thermal"], output_fields=spec["fields"],
        q_supply=repr(cfgkeys["q_supply"]),
        groove_mixing_model=cfgkeys["groove_mixing_model"],
        coupling_max_iters=COUPLING_MAX_ITERS,
        coupling_relaxation=repr(COUPLING_RELAXATION)) + COUPLING_EXTRA)
    rc, wall = S.run_case(run_dir)
    m = S.frame_metrics(run_dir) if rc == 0 else None
    rec = {"dir": f"{group}/{run_dir.name}", "cfgkey": key, "eps": eps, "psi": psi,
           "tag": tag, "rc": rc, "wall_s": round(wall, 1),
           "ok": bool(m and "W" in m)}
    if rec["ok"]:
        rec.update({k: m[k] for k in METRIC_KEYS if k in m})
    log["iters"].append(rec)
    save_log(group, log)
    if rec["ok"]:
        cpl = rec.get("coupling") or {}
        S.progress(PROGRESS,
                   f"{rec['dir']} [{tag}] eps={eps:.6f} psi={psi:.4f} -> "
                   f"W={rec['W']:.2f} N angle={rec['force_angle_deg']:.4f} deg "
                   f"Tmax={rec.get('T_max', float('nan'))-273.15:.2f} C "
                   f"picard={cpl.get('iters')} ({wall:.0f} s)")
    else:
        S.progress(PROGRESS, f"{rec['dir']} [{tag}] FAILED rc={rc}")
    return rec


def converged(rec, w_target):
    d_ang = (rec["force_angle_deg"] - ANGLE_TARGET + 180.0) % 360.0 - 180.0
    return (abs(rec["W"] - w_target) / w_target < W_TOL_REL
            and abs(d_ang) < ANGLE_TOL_DEG), d_ang


def newton(group, label, cfgkeys, w_target, eps0, psi0, max_iter=12):
    log = load_log(group)
    target = (0.0, -w_target)
    eps, psi = eps0, psi0
    rec = evaluate(group, log, cfgkeys, eps, psi, f"{label}/x0")
    for it in range(max_iter):
        if not rec["ok"]:
            raise RuntimeError(f"{group}/{label}: solver failed at x{it}")
        ok, d_ang = converged(rec, w_target)
        S.progress(PROGRESS, f"{group}/{label} iter {it}: eps={eps:.6f} psi={psi:.4f} "
                             f"W={rec['W']:.2f} (dW={rec['W']-w_target:+.2f} N) "
                             f"dAngle={d_ang:+.4f} deg")
        if ok:
            S.progress(PROGRESS, f"{group}/{label} CONVERGED eps={eps:.6f} psi={psi:.4f}")
            return {"group": group, "label": label, "cfg": cfgkeys, "W_target": w_target,
                    "eps": eps, "psi": psi, "dir": rec["dir"], "metrics": rec,
                    "newton_steps": it}
        r = (rec["Fx"] - target[0], rec["Fy"] - target[1])
        re_ = evaluate(group, log, cfgkeys, eps + D_EPS, psi, f"{label}/J_eps")
        rp_ = evaluate(group, log, cfgkeys, eps, psi + D_PSI, f"{label}/J_psi")
        if not (re_["ok"] and rp_["ok"]):
            raise RuntimeError(f"{group}/{label}: Jacobian probe failed")
        j = [[(re_["Fx"] - rec["Fx"]) / D_EPS, (rp_["Fx"] - rec["Fx"]) / D_PSI],
             [(re_["Fy"] - rec["Fy"]) / D_EPS, (rp_["Fy"] - rec["Fy"]) / D_PSI]]
        det = j[0][0] * j[1][1] - j[0][1] * j[1][0]
        if abs(det) < 1e-12:
            raise RuntimeError(f"{group}/{label}: singular Jacobian")
        d_e = (-r[0] * j[1][1] + r[1] * j[0][1]) / det
        d_p = (-j[0][0] * r[1] + j[1][0] * r[0]) / det
        scale = min(1.0, 0.15 / max(abs(d_e), 1e-12), 20.0 / max(abs(d_p), 1e-12))
        eps = min(0.95, max(0.02, eps + scale * d_e))
        psi = psi + scale * d_p
        rec = evaluate(group, log, cfgkeys, eps, psi, f"{label}/x{it+1}")
    raise RuntimeError(f"{group}/{label}: not converged in {max_iter} Newton steps")


def solve_point(group, w_target, eps0, psi0):
    pinned = dict(groove_mixing_model="PINNED", q_supply=Q_SUPPLY_PLACEHOLDER)
    a = newton(group, "pinned", pinned, w_target, eps0, psi0)
    q_s = a["metrics"]["q_supply"]
    S.progress(PROGRESS, f"{group}: derived Q_supply = {q_s:.6e} m^3/s "
                         f"({q_s*6e4:.2f} L/min) from the PINNED steady mass balance")
    flow = dict(groove_mixing_model="FLOW_WEIGHTED", q_supply=q_s)
    b = newton(group, "flow", flow, w_target, a["eps"], a["psi"])
    return {"pinned": a, "flow": b, "q_supply": q_s}


def main():
    wanted = sys.argv[1:] or list(GROUPS)
    out_path = CASE_DIR / "equilibrium.json"
    results = json.loads(out_path.read_text()) if out_path.exists() else {}
    failures = {}

    for group in wanted:
        if group not in GROUPS:
            raise SystemExit(f"unknown group {group!r}; choose from {list(GROUPS)}")
        if group in results:
            continue
        eps0, psi0 = GROUPS[group]["warm"]
        S.progress(PROGRESS, f"{group}: warm start from eps={eps0:.6f} psi={psi0:.4f}")
        # One variant failing must not stall the batch.
        try:
            results[group] = solve_point(group, W_TARGET, eps0, psi0)
        except Exception as exc:
            failures[group] = repr(exc)
            S.progress(PROGRESS, f"{group}: ABANDONED -- {exc!r}")
            continue
        out_path.write_text(json.dumps(results, indent=2))

    if failures:
        (CASE_DIR / "failures.json").write_text(json.dumps(failures, indent=2))
    print(json.dumps({k: {"eps": v["flow"]["eps"], "psi": v["flow"]["psi"]}
                      for k, v in results.items()}, indent=2))


if __name__ == "__main__":
    main()
