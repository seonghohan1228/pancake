"""Run one pancake steady case and read its final-frame integral metrics.

Case 01-THD copy of the 01-prelim helper. Two differences: it runs on 4 ranks
(the machine is otherwise idle for this batch) and it also reports the thermal
fields the energy-equation runs produce.
"""

import math
import os
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path

EXE = Path(__file__).resolve().parent / "bin" / "pancake.exe"
MINGW_BIN = r"C:\msys64\mingw64\bin"
NPROC = 4
RHO = 860.0  # liquid density [kg/m^3]; converts inlet_mass_source to a volume rate


def progress(progress_file, message):
    """Append one timestamped line to a case progress.log and echo it."""
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp}  {message}"
    with open(progress_file, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def run_case(work_dir, cfg="config.txt", log_name="run.log"):
    """Run the solver in work_dir. Returns (returncode, wall_seconds)."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = MINGW_BIN + os.pathsep + env["PATH"]
    t0 = time.perf_counter()
    with open(work_dir / log_name, "w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.run(
            ["mpiexec", "-n", str(NPROC), str(EXE), "--config", cfg],
            cwd=str(work_dir), stdout=log, stderr=subprocess.STDOUT,
            env=env, text=True, encoding="utf-8", errors="replace",
        )
    return proc.returncode, time.perf_counter() - t0


def _piece_arrays(vts_path):
    root = ET.parse(vts_path).getroot()
    out = {}
    for cd in root.iter("CellData"):
        for da in cd.findall("DataArray"):
            ncomp = int(da.get("NumberOfComponents", "1"))
            out[da.get("Name")] = ([float(v) for v in (da.text or "").split()], ncomp)
    return out


def _read_pvts(pvts_path):
    """Merge every rank piece of a .pvts, in the piece order written by IO."""
    root = ET.parse(pvts_path).getroot()
    merged = {}
    for piece in root.iter("Piece"):
        src = pvts_path.parent / piece.get("Source")
        for name, (vals, ncomp) in _piece_arrays(src).items():
            merged.setdefault(name, ([], ncomp))[0].extend(vals)
    return merged


def final_pvts(run_dir):
    """Path of the last flat .pvts frame, or None."""
    pvd = Path(run_dir) / "results" / "flat" / "results.pvd"
    if not pvd.exists():
        return None
    root = ET.parse(pvd).getroot()
    items = [(float(ds.get("timestep")), pvd.parent / ds.get("file"))
             for ds in root.iter("DataSet")]
    if not items:
        return None
    items.sort(key=lambda x: x[0])
    return items[-1][1]


def frame_metrics(run_dir):
    """Integral metrics of the final frame: forces, p_max, torque, cavitation."""
    path = final_pvts(run_dir)
    if path is None:
        return None
    a = _read_pvts(path)
    m = {}
    if "Fp" in a:
        fx, fy = a["Fp"][0][0], a["Fp"][0][1]
        m["Fx"], m["Fy"] = fx, fy
        m["W"] = math.hypot(fx, fy)
        m["force_angle_deg"] = math.degrees(math.atan2(fy, fx))
    if "p" in a:
        m["p_max"] = max(a["p"][0])
        m["p_min"] = min(a["p"][0])
    if "film_content" in a:
        fc = a["film_content"][0]
        m["cav_frac"] = sum(1 for v in fc if v < 0.999) / len(fc)
        m["film_content_min"] = min(fc)
    if "friction_torque" in a:
        m["friction_torque"] = a["friction_torque"][0][0]
    if "h" in a:
        m["h_min"] = min(a["h"][0])
    if "T" in a:
        t = a["T"][0]
        m["T_max"] = max(t)
        m["T_min"] = min(t)
        # theta of the hottest cell, degrees measured from the groove centre.
        i_hot = max(range(len(t)), key=lambda k: t[k])
        m["T_max_theta_deg"] = (math.degrees(a["theta_rad"][0][i_hot]) - 90.0) % 360.0
    if "mu" in a:
        m["mu_min"] = min(a["mu"][0])
        m["mu_max"] = max(a["mu"][0])
    d = diagnostics_tail(run_dir)
    if "inlet_mass_source" in d:
        # Steady state: the net liquid injected by the groove penalty cells is
        # the fresh-oil supply rate. Volume rate = mass rate / rho.
        m["q_supply"] = float(d["inlet_mass_source"]) / RHO
    if "liquid_residual" in d:
        m["liquid_residual"] = float(d["liquid_residual"])
    m["coupling"] = coupling_health(run_dir)
    return m


def coupling_health(run_dir, log_name="run.log"):
    """Parse the last `Coupling: N outer iters, R_p=.. R_T=..` line of a run log."""
    path = Path(run_dir) / log_name
    if not path.exists():
        return None
    last = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Coupling:" in line:
            last = line
    if last is None:
        return None
    body = last.split("Coupling:", 1)[1]
    try:
        iters = int(body.split("outer iters")[0].strip())
        r_p = float(body.split("R_p=")[1].split()[0])
        r_t = float(body.split("R_T=")[1].split()[0].rstrip("\x1b[0m"))
    except (IndexError, ValueError):
        return None
    return {"iters": iters, "R_p": r_p, "R_T": r_t}


def midplane_profile(run_dir, name="p"):
    """Final-frame mid-span field(theta): (theta_deg list, value list), sorted.

    Uses the solver's own theta_rad / z_m cell fields rather than assuming a
    piece ordering: pieces are theta-blocked but cells inside a piece run
    theta-fastest, so a naive reshape of the concatenation is wrong.
    """
    path = final_pvts(run_dir)
    if path is None:
        return None
    a = _read_pvts(path)
    if name not in a:
        return None
    theta = a["theta_rad"][0]
    z = a["z_m"][0]
    v = a[name][0]
    z_sorted = sorted(set(z))
    z_mid = z_sorted[len(z_sorted) // 2]
    rows = [(math.degrees(t) % 360.0, vv)
            for t, zv, vv in zip(theta, z, v) if zv == z_mid]
    rows.sort()
    return [r[0] for r in rows], [r[1] for r in rows]


def midplane_pressure(run_dir):
    """Final-frame mid-span p(theta). Kept for the 01-prelim call signature."""
    return midplane_profile(run_dir, "p")


def diagnostics_tail(run_dir):
    """Last row of results/diagnostics.csv as a dict of strings."""
    csv_path = Path(run_dir) / "results" / "diagnostics.csv"
    if not csv_path.exists():
        return {}
    lines = [l for l in csv_path.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        return {}
    return dict(zip(lines[0].split(","), lines[-1].split(",")))
