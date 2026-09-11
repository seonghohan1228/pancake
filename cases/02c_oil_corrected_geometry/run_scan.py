"""Run archived steady FBNS evaluations and close vector load by damped Newton.

Usage: python3 run_scan.py pilot|scan|refine|low. settings.json defines the geometry.
Each rpm is independent; evaluations within a load continuation remain serial.
"""
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import numpy as np

ROOT = Path(__file__).resolve().parent
CFG = json.loads((ROOT / "settings.json").read_text())
TEMPLATE = (ROOT / "config_template.txt").read_text()


def read_fields(run_dir):
    """Read final flat VTK frame using physical coordinates, independent of ranks."""
    pvd = run_dir / "results/flat/results.pvd"
    frames = list(ET.parse(pvd).getroot().iter("DataSet"))
    frame = max(frames, key=lambda x: float(x.get("timestep")))
    pvts = pvd.parent / frame.get("file")
    arrays = {}
    for piece in ET.parse(pvts).getroot().iter("Piece"):
        for cell in ET.parse(pvts.parent / piece.get("Source")).getroot().iter("CellData"):
            for da in cell.findall("DataArray"):
                arrays.setdefault(da.get("Name"), []).append(
                    np.fromstring(da.text or "", sep=" "))
    return {key: np.concatenate(val) for key, val in arrays.items()}


def evaluate(group, eps, psi, rpm, grid=None, p_cav=None):
    """Run one prescribed-position solve; reject unconverged or invalid fields."""
    grid = grid or CFG["grid"]
    floor = CFG["p_cav"] if p_cav is None else p_cav
    args = dict(eps=float(eps), psi=float(psi), rpm=rpm, grid=grid, p_cav=floor)
    identity = hashlib.sha256(json.dumps([CFG, args], sort_keys=True).encode()).hexdigest()[:14]
    run_dir = ROOT / "runs" / group / identity
    result_path = run_dir / "metrics.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    run_dir.mkdir(parents=True, exist_ok=True)
    thermal = "temperature_model = ISOTHERMAL\nfluid_property_model = CONSTANT\nliquid_viscosity_model = CONSTANT\nviscosity_temperature_coeff = 0\ncross_film_thermal_model = NONE\n"
    config = TEMPLATE.format(c=CFG["c"], e=eps*CFG["c"], L=CFG["L"],
        n_theta=grid[0], n_z=grid[1], omega=rpm*math.pi/30, mu=CFG["mu"],
        p_cav=floor, p_end=CFG["p_end"], thermal_block=thermal)
    config = re.sub(r"(?m)^R = .*", f'R = {CFG["R"]}', config)
    config = re.sub(r"(?m)^rho = .*", f'rho = {CFG["rho"]}', config)
    config = re.sub(r"(?m)^attitude_angle_deg = .*", f"attitude_angle_deg = {psi}", config)
    config += f'\ninlet_groove = {CFG["groove_angle"]} {CFG["groove_width"]} {CFG["p_end"]} 293.15\n'
    (run_dir / "config.txt").write_text(config)
    started = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMPI_MCA_hwloc_base_binding_policy="none")
    with (run_dir / "run.log").open("w") as log:
        process = subprocess.run(["mpiexec", "-n", str(CFG["ranks"]), CFG["executable"],
            "--config", "config.txt"], cwd=run_dir, stdout=log, stderr=subprocess.STDOUT,
            env=env, timeout=180)
    log = (run_dir / "run.log").read_text(errors="replace")
    residuals = re.findall(r"Coupling:.*?R_p=([\deE+.\-]+).*?R_T=([\deE+.\-]+)", log)
    if process.returncode or not residuals or float(residuals[-1][0]) > 1e-6:
        raise RuntimeError(f"Invalid steady solve {run_dir}; rc={process.returncode}; residual={residuals[-1:]}")
    fbns = re.findall(r"fbns iter=\d+ \|G\|=([\deE+.\-]+) \|FB\|=([\deE+.\-]+)", log)
    if not fbns or max(map(float,fbns[-1])) > 1e-7:
        raise RuntimeError(f"FBNS residual failed in {run_dir}")
    a = read_fields(run_dir)
    p, fc = a["p"], a["film_content"]
    if not np.isfinite(p).all() or not np.isfinite(fc).all() or p.min() < floor-1 or fc.min() < -1e-6:
        raise RuntimeError(f"Nonphysical fields in {run_dir}")
    theta, z = a["theta_rad"], a["z_m"]
    xs, zs = np.unique(theta), np.unique(z)
    planes = {}
    for field in ("p", "film_content", "h", "mu"):
        plane = np.full((len(zs), len(xs)), np.nan)
        plane[np.searchsorted(zs,z), np.searchsorted(xs,theta)] = a[field]
        if not np.isfinite(plane).all():
            raise RuntimeError("Incomplete VTK plane")
        planes[field] = plane
    np.savez_compressed(run_dir / "fields.npz", theta_deg=xs*180/math.pi,
                       z_over_L=zs/CFG["L"]-.5, **planes)
    mask = fc < .999
    midmask = planes["film_content"][len(zs)//2] < .999
    fx, fy = a["Fp"][:2]
    area = CFG["R"] * 2*math.pi / len(xs) * CFG["L"] / len(zs)
    cavity_volume = np.sum((1-np.minimum(fc,1))*a["h"])*area
    diagnostics = run_dir / "results/diagnostics.csv"
    diag = list(csv.DictReader(diagnostics.open()))[-1] if diagnostics.exists() else {}
    row = dict(**args, dir=run_dir.relative_to(ROOT).as_posix(), Fx=float(fx), Fy=float(fy),
        W=float(math.hypot(fx,fy)), force_angle_deg=math.degrees(math.atan2(fy,fx)),
        p_min_bar=float(p.min()/1e5), p_max_bar=float(p.max()/1e5),
        cav_area_pct=float(mask.mean()*100), mid_arc_deg=float(midmask.sum()*360/len(xs)),
        void_max=float(1-fc.min()), void_volume_mm3=float(cavity_volume*1e9),
        h_min_um=float(a["h"].min()*1e6), h_min_exact_um=(1-eps)*CFG["c"]*1e6,
        torque_Nm=float(a["friction_torque"][0]),
        R_p=float(residuals[-1][0]), R_T=float(residuals[-1][1]),
        wall_s=time.monotonic()-started, diagnostics=diag)
    result_path.write_text(json.dumps(row, indent=2))
    print(f"eval {group}: eps={eps:.5f} psi={psi:.3f} W={row['W']:.1f} arc={row['mid_arc_deg']:.1f} wall={row['wall_s']:.1f}s", flush=True)
    return row


def balance(load, rpm, warm, grid=None, p_cav=None, prefix="scan"):
    """Close both force components to 0.2 percent of the imposed load."""
    group = f"{prefix}_n{rpm}_w{load}"
    result_path = ROOT / f"{group}.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    angle = math.radians(CFG["target_force_angle"])
    target = load*np.array([math.cos(angle), math.sin(angle)])
    eps, psi = warm
    history = []
    for iteration in range(18):
        row = evaluate(group, eps, psi, rpm, grid, p_cav)
        force = np.array([row["Fx"],row["Fy"]])
        residual = force-target
        error = float(np.linalg.norm(residual)/load)
        history.append(dict(eps=eps, psi=psi, force_error=error, dir=row["dir"]))
        if error < .002:
            row.update(load_target=load, force_error_rel=error, history=history)
            result_path.write_text(json.dumps(row, indent=2))
            print(f"DONE {group}: error={100*error:.3f}% eps={eps:.5f} arc={row['mid_arc_deg']:.1f}", flush=True)
            return row
        de = min(.002, (.997-eps)*.2)
        dp = .25
        reval = evaluate(group, eps+de, psi, rpm, grid, p_cav)
        pval = evaluate(group, eps, psi+dp, rpm, grid, p_cav)
        jac = np.column_stack(((np.array([reval["Fx"],reval["Fy"]])-force)/de,
                               (np.array([pval["Fx"],pval["Fy"]])-force)/dp))
        step = np.linalg.solve(jac,-residual)
        scale = min(1., .08/max(abs(step[0]),1e-12), 12/max(abs(step[1]),1e-12))
        eps = float(np.clip(eps+scale*step[0], .02, .995))
        psi = float(psi+scale*step[1])
    raise RuntimeError(f"Force balance failed {group}: {history[-1]}")


def scan_rpm(rpm):
    """Continue increasing loads at a fixed speed."""
    warm = (min(.55,.10*1000/rpm), -10.)
    for load in CFG["loads"]:
        row = balance(load,rpm,warm)
        warm = row["eps"],row["psi"]


def main():
    mode = sys.argv[1] if len(sys.argv)>1 else "pilot"
    if mode == "pilot":
        print(json.dumps(evaluate("pilot",.10,-10,1000),indent=2))
    elif mode == "scan":
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
            for result in pool.map(scan_rpm,CFG["rpms"]):
                pass
    elif mode == "refine":
        points = [(1000,2000),(3000,250),(3000,500),(3000,2000),(5000,500),(5000,250)]
        for load,rpm in points:
            row = json.loads((ROOT/f"scan_n{rpm}_w{load}.json").read_text())
            balance(load,rpm,(row["eps"],row["psi"]),[480,128],prefix="fine")
    elif mode == "low":
        row = balance(1000,100,(.55,-30),prefix="low")
        balance(3000,100,(row["eps"],row["psi"]),prefix="low")
    else:
        raise SystemExit("Choose pilot, scan, refine or low")


if __name__ == "__main__":
    main()
