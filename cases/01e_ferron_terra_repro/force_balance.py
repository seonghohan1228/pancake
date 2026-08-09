"""Isolated adaptive two-pass Ferron driver for case 01e."""

import json
import math
import os
import sys
import time
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
RUN_ROOT = CASE_DIR / "runs" / "production"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text(encoding="utf-8")
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S

CLEARANCE = 1.52e-4
Q_PLACEHOLDER = 1.0e-4
W_TOL_REL, ANGLE_TOL = 1.0e-3, 0.05
D_EPS, D_PSI = 5.0e-3, 0.5
COMMON_XF = """cross_film_thermal_model = POLYNOMIAL
cross_film_profile_order = 4
cross_film_quadrature_points = 6
cross_film_relaxation = 1.0
journal_wall_temperature = 313.15
bearing_wall_temperature = 313.15
journal_heat_transfer = 0
bearing_heat_transfer = 0
thermal_ambient_temperature = 313.15
shaft_ambient_area_ratio = 1.0
"""
BASE_FIELDS = "pressure, temperature, film_content, h, mu, mass_residual, pressure_force_x, pressure_force_y, fluid_force_x, fluid_force_y, friction_torque"
XF_FIELDS = BASE_FIELDS + ", temperature_bearing, temperature_journal"

def conjugate(boundary, h=0.0, shaft_h=0.0, shaft_t=313.15, bush_t=313.15):
    return COMMON_XF + f"""shaft_ambient_heat_transfer = {shaft_h}
shaft_ambient_temperature = {shaft_t}
bush_thermal_model = CONDUCTION_2D
bush_thermal_conductivity = 250
bush_outer_radius = 0.10
bush_radial_cells = 16
bush_outer_boundary = {boundary}
bush_outer_heat_transfer = {h}
bush_outer_temperature = {bush_t}
"""

def depth_average(journal_h):
    return f"""journal_wall_temperature = 313.15
bearing_wall_temperature = 313.15
journal_heat_transfer = {journal_h}
bearing_heat_transfer = 7200
cross_film_thermal_model = NONE
"""

def crossfilm_local(shaft_h, bush_h, k=50):
    return COMMON_XF + f"""shaft_ambient_heat_transfer = {shaft_h}
bush_thermal_model = LOCAL
bush_thermal_conductivity = {k}
bush_outer_radius = 0.10
bush_outer_boundary = CONVECTION
bush_outer_heat_transfer = {bush_h}
"""

def spec(label, grid, thermal, target=4000.0, omega=209.43951, warm=(0.60, -40.0), kind="MODEL-CLOSURE", fields=XF_FIELDS):
    return dict(label=label, grid=grid, thermal=thermal, target=target, omega=omega, warm=warm, kind=kind, fields=fields)

SPECS = {
    "r0_g6": spec("R0 G6 legacy h=10", (720,192), conjugate("CONVECTION", 10), warm=(0.6373,-40.11)),
    "r0_g7": spec("R0 G7 legacy h=10", (1080,288), conjugate("CONVECTION", 10), warm=(0.6364,-40.04)),
    "r1_g6": spec("R1 G6 outer fixed 40 C", (720,192), conjugate("FIXED_TEMPERATURE"), warm=(0.60,-37.5)),
    "r1_g7": spec("R1 G7 outer fixed 40 C", (1080,288), conjugate("FIXED_TEMPERATURE"), warm=(0.60,-37.5)),
    "r2_g6": spec("R2 G6 fixed 40 C + scalar shaft 43 C", (720,192), conjugate("FIXED_TEMPERATURE", shaft_h=1e7, shaft_t=316.15), warm=(0.60,-37.5), kind="SURROGATE + MODEL-CLOSURE"),
    "r2_g7": spec("R2 G7 fixed 40 C + scalar shaft 43 C", (1080,288), conjugate("FIXED_TEMPERATURE", shaft_h=1e7, shaft_t=316.15), warm=(0.60,-37.5), kind="SURROGATE + MODEL-CLOSURE"),
    "d_grid_g1": spec("D-grid G1 current R1 closure", (60,16), conjugate("FIXED_TEMPERATURE"), warm=(0.61,-38.7), kind="MODEL-CLOSURE"),
    "d_grid_g2": spec("D-grid G2 current R1 closure", (120,32), conjugate("FIXED_TEMPERATURE"), warm=(0.61,-38.7), kind="MODEL-CLOSURE"),
    "d_grid_g3": spec("D-grid G3 current R1 closure", (240,64), conjugate("FIXED_TEMPERATURE"), warm=(0.61,-38.7), kind="MODEL-CLOSURE"),
    "d_grid_g4": spec("D-grid G4 current R1 closure", (360,96), conjugate("FIXED_TEMPERATURE"), warm=(0.61,-38.7), kind="MODEL-CLOSURE"),
    "d_grid_g5": spec("D-grid G5 current R1 closure", (480,128), conjugate("FIXED_TEMPERATURE"), warm=(0.60,-37.5), kind="MODEL-CLOSURE"),
    "d_grid_g6": spec("D-grid G6 current R1 closure", (720,192), conjugate("FIXED_TEMPERATURE"), warm=(0.60,-37.5), kind="MODEL-CLOSURE"),
    "d_grid_g7": spec("D-grid G7 current R1 closure", (1080,288), conjugate("FIXED_TEMPERATURE"), warm=(0.60,-37.5), kind="MODEL-CLOSURE"),
    "d_grid_g8": spec("D-grid G8 current R1 closure", (1620,432), conjugate("FIXED_TEMPERATURE"), warm=(0.609,-38.68), kind="MODEL-CLOSURE"),
    "c_fig4_jadiabatic": spec("C Fig.4 J adiabatic", (720,192), depth_average(0), warm=(0.5738,-36.32), fields=BASE_FIELDS),
    "c_fig4_jiso40": spec("C Fig.4 J iso40", (720,192), depth_average(10000), warm=(0.5675,-36.12), fields=BASE_FIELDS),
    "c_xf_adiabatic": spec("C cross-film adiabatic", (720,192), crossfilm_local(0,0), warm=(0.6086,-37.55)),
    "c_xf_bush40_j40": spec("C cross-film bush40/J40", (720,192), crossfilm_local(1e7,1e7,1e7), warm=(0.5646,-36.14)),
    "c_xf_bush40_jad": spec("C cross-film bush40/J adiabatic", (720,192), crossfilm_local(0,1e7,1e7), warm=(0.5860,-36.8)),
    "c_conj_fixed": spec("C conjugate outer fixed 40 C", (720,192), conjugate("FIXED_TEMPERATURE"), warm=(0.5934,-37.50)),
    "c_conj_robin80": spec("C conjugate outer Robin h=80", (720,192), conjugate("CONVECTION",80), warm=(0.6028,-37.70)),
}
for w in (2000,3000,4000,6000,8000,10000):
    SPECS[f"c_fig5_w{w}"] = spec(f"C Fig.5 W={w} N", (720,192), depth_average(0), target=float(w), warm=(0.58,-36.0), fields=BASE_FIELDS)

MODES = {
    "R0": ["r0_g6", "r0_g7"], "R1": ["r1_g6", "r1_g7"], "R2": ["r2_g6", "r2_g7"],
    "D_GRID": [f"d_grid_g{i}" for i in range(1, 9)],
    "C": ["c_fig4_jadiabatic", "c_fig4_jiso40", *[f"c_fig5_w{w}" for w in (2000,3000,4000,6000,8000,10000)], "c_xf_adiabatic", "c_xf_bush40_jad", "c_xf_bush40_j40", "c_conj_fixed", "c_conj_robin80"],
}

def progress(message):
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp}  {message}"
    with (CASE_DIR / "progress.log").open("a", encoding="utf-8") as handle: handle.write(line + "\n")
    print(line, flush=True)

def converged_coupling(metrics):
    c = metrics.get("coupling") or {}
    return c.get("R_p", math.inf) < 1e-6 and c.get("R_T", math.inf) < 1e-6

def log_file(name): return RUN_ROOT / name / "fb_log.json"
def load_log(name): return json.loads(log_file(name).read_text()) if log_file(name).exists() else {"iters": []}
def save_log(name, log):
    log_file(name).parent.mkdir(parents=True, exist_ok=True)
    log_file(name).write_text(json.dumps(log, indent=2))

def evaluate(name, cfg, eps, psi, tag):
    log = load_log(name)
    key = json.dumps(cfg, sort_keys=True)
    for row in log["iters"]:
        if row.get("cfg") == key and row["eps"] == eps and row["psi"] == psi and row["ok"]: return row
    n = len(log["iters"]); run_dir = RUN_ROOT / name / f"fb_iter_{n:02d}"; run_dir.mkdir(parents=True, exist_ok=True)
    s = SPECS[name]
    (run_dir / "config.txt").write_text(TEMPLATE.format(e=eps*CLEARANCE, attitude=psi, n_theta=s["grid"][0], n_z=s["grid"][1], omega=s["omega"], thermal_block=s["thermal"], q_supply=cfg["q_supply"], mixing=cfg["mixing"], output_fields=s["fields"]))
    rc, wall = S.run_case(run_dir); metrics = S.frame_metrics(run_dir) if rc == 0 else None
    row = {"dir": run_dir.relative_to(CASE_DIR).as_posix(), "cfg": key, "eps": eps, "psi": psi, "tag":tag, "rc":rc, "wall_s":round(wall,1), "ok":bool(metrics and "W" in metrics)}
    if metrics: row.update(metrics); row["coupling_ok"] = converged_coupling(metrics); row["ok"] &= row["coupling_ok"]
    log["iters"].append(row); save_log(name,log); progress(f"{row['dir']} {tag}: rc={rc}, coupling_ok={row.get('coupling_ok',False)}")
    return row

def newton(name, cfg, eps, psi, target):
    target_vector=(0.0,-target)
    for step in range(12):
        row=evaluate(name,cfg,eps,psi,f"x{step}")
        if not row["ok"]: raise RuntimeError(f"{name}: failed/adaptive coupling not converged")
        angle=(row["force_angle_deg"]+90+180)%360-180
        if abs(row["W"]-target)/target < W_TOL_REL and abs(angle)<ANGLE_TOL: return {"eps":eps,"psi":psi,"dir":row["dir"],"metrics":row,"steps":step}
        eprobe=evaluate(name,cfg,eps+D_EPS,psi,"J_eps"); pprobe=evaluate(name,cfg,eps,psi+D_PSI,"J_psi")
        if not eprobe["ok"] or not pprobe["ok"]: raise RuntimeError(f"{name}: Jacobian probe failed")
        j00=(eprobe["Fx"]-row["Fx"])/D_EPS; j01=(pprobe["Fx"]-row["Fx"])/D_PSI; j10=(eprobe["Fy"]-row["Fy"])/D_EPS; j11=(pprobe["Fy"]-row["Fy"])/D_PSI
        det=j00*j11-j01*j10
        if abs(det)<1e-12: raise RuntimeError(f"{name}: singular force Jacobian")
        rx,ry=row["Fx"]-target_vector[0],row["Fy"]-target_vector[1]
        de=(-rx*j11+ry*j01)/det; dp=(-j00*ry+j10*rx)/det; scale=min(1.0,0.15/max(abs(de),1e-12),20/max(abs(dp),1e-12))
        eps=min(.95,max(.02,eps+scale*de)); psi+=scale*dp
    raise RuntimeError(f"{name}: force balance not closed")

def solve(name):
    s=SPECS[name]; eps,psi=s["warm"]
    pinned=newton(name,{"mixing":"PINNED","q_supply":Q_PLACEHOLDER},eps,psi,s["target"])
    q=pinned["metrics"]["q_supply"]; progress(f"{name}: source-derived Q={q*6e4:.4f} L/min")
    flow=newton(name,{"mixing":"FLOW_WEIGHTED","q_supply":q},pinned["eps"],pinned["psi"],s["target"])
    return {"label":s["label"],"kind":s["kind"],"grid":s["grid"],"target_N":s["target"],"q_supply_m3s":q,"pinned":pinned,"flow":flow}

def main():
    global RUN_ROOT
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument("mode", choices=MODES); parser.add_argument("--smoke",action="store_true"); args=parser.parse_args()
    result_path=CASE_DIR/("smoke_equilibrium.json" if args.smoke else "campaign_equilibrium.json")
    RUN_ROOT = CASE_DIR / "runs" / ("smoke" if args.smoke else "production")
    data=json.loads(result_path.read_text()) if result_path.exists() else {}
    for name in MODES[args.mode]:
        if name in data and "failed" not in data[name]: progress(f"{name}: already archived, skipped"); continue
        data.pop(name, None)
        if args.smoke: SPECS[name]={**SPECS[name],"grid":(180,48)}
        try: data[name]=solve(name); result_path.write_text(json.dumps(data,indent=2))
        except Exception as exc: data[name]={"failed":str(exc),"label":SPECS[name]["label"]}; result_path.write_text(json.dumps(data,indent=2)); raise

if __name__ == "__main__": main()
