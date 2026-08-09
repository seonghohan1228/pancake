"""Run the D-dt physical-time ladder sequentially after D-grid G6."""
import csv
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import force_balance as F
S = F.S
LADDER = {"d_dt_tn2": 3.2e-4, "d_dt_tn1": 1.6e-4,
          "d_dt_t0": 8e-5, "d_dt_t1": 4e-5, "d_dt_t2": 2e-5,
          "d_dt_t3": 1e-5, "d_dt_t4": 5e-6}

def dynamic_config(dt, q_supply, grid=(720, 192), end_t=0.06, write_interval=5e-4):
    s = F.SPECS["d_grid_g6"]
    text = F.TEMPLATE.format(e=.1 * F.CLEARANCE, attitude=-90., n_theta=grid[0], n_z=grid[1],
        omega=209.43951, thermal_block=s["thermal"], q_supply=q_supply,
        mixing="FLOW_WEIGHTED", output_fields=F.XF_FIELDS + ", bearing_x, bearing_y, bearing_attitude_angle")
    text = text.replace("end_t = 1.0", f"end_t = {end_t}").replace("dt = 1e-5", f"dt = {dt}")
    text = text.replace("write_interval = 1.0", f"write_interval = {write_interval}")
    text = text.replace("solution_mode = STEADY_STATE", "solution_mode = TRANSIENT")
    text = text.replace("motion_model = STATIC", "motion_model = MOVING_BEARING")
    text = text.replace("bearing_mass = 1.64", "bearing_mass = 10")
    return text.replace("external_load_y = 0", "external_load_y = 4000")

def dynamic_metrics(run_dir):
    metrics = S.frame_metrics(run_dir)
    if metrics is None: return None
    final_path = S.final_pvts(run_dir)
    arrays = S._read_pvts(final_path)
    position = arrays.get("xB", ([float("nan")] * 3, 3))[0]
    x, y = position[0], position[1]
    metrics["eps_final"] = math.hypot(x, y) / F.CLEARANCE
    metrics["attitude_final_deg"] = math.degrees(math.atan2(-y, -x))
    pvd = Path(run_dir) / "results" / "flat" / "results.pvd"
    metrics["final_time_s"] = max(
        float(item.get("timestep")) for item in ET.parse(pvd).getroot().iter("DataSet")
    )
    diagnostics_path = Path(run_dir) / "results" / "diagnostics.csv"
    with diagnostics_path.open(encoding="utf-8") as handle:
        diagnostics = list(csv.DictReader(handle))
    metrics["flow_converged_all_steps"] = bool(diagnostics) and all(
        row.get("converged") == "1" for row in diagnostics
    )
    metrics["diagnostic_steps"] = len(diagnostics)
    return metrics

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        ladder = {"d_dt_smoke": 2e-5}
        q_supply = 8e-5
        grid, end_t, write_interval = (180, 48), 4e-5, 2e-5
        output = CASE_DIR / "smoke_d_timestep_runs.json"
        run_root = CASE_DIR / "runs" / "smoke"
    else:
        eq_path = CASE_DIR / "campaign_equilibrium.json"
        if not eq_path.exists(): raise SystemExit("run D_GRID first: no local equilibrium ledger")
        equilibria = json.loads(eq_path.read_text())
        if "d_grid_g6" not in equilibria: raise SystemExit("run D_GRID first: D-grid G6 source missing")
        q_supply = equilibria["d_grid_g6"]["q_supply_m3s"]
        ladder = LADDER
        grid, end_t, write_interval = (720, 192), 0.06, 5e-4
        output = CASE_DIR / "d_timestep_runs.json"
        run_root = CASE_DIR / "runs" / "production"
    runs = json.loads(output.read_text()) if output.exists() else {}
    for name, dt in ladder.items():
        if name in runs:
            F.progress(f"{name}: already archived, skipped"); continue
        run_dir = run_root / name
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "config.txt").write_text(
            dynamic_config(dt, q_supply, grid, end_t, write_interval), encoding="utf-8"
        )
        F.progress(f"{name}: start dt={dt:.1e} s, source-derived Q={q_supply*6e4:.4f} L/min")
        rc, wall = S.run_case(run_dir); metrics = dynamic_metrics(run_dir) if rc == 0 else None
        runs[name] = {"dt_s": dt, "run_dir": run_dir.relative_to(CASE_DIR).as_posix(), "rc": rc, "wall_s": round(wall, 1),
                      "q_supply_m3s": q_supply, "metrics": metrics}
        output.write_text(json.dumps(runs, indent=2)); F.progress(f"{name}: rc={rc}")
        if rc != 0 or metrics is None:
            F.progress(f"{name}: archived failed stability/completion verdict; continuing ladder")

if __name__ == "__main__": main()
