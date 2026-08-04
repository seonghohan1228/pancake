"""Metric extraction for case 00 (grid / timestep independence).

Reads each run's results/flat/results.pvd series (ASCII .vts pieces) plus
results/diagnostics.csv and writes metrics.json in this directory.

Metrics per run:
  W              load magnitude sqrt(Fpx^2 + Fpy^2)                       [N]
  force_angle    atan2(Fpy, Fpx) of the pressure-force resultant          [deg]
  p_max          peak cell pressure in the flat VTK                       [Pa]
  cav_frac       area fraction of cells with film_content < 0.999         [-]
  torque         friction_torque (uniform field)                          [N m]
  T_max, T_mean  temperature extrema (THD runs only)                      [K]
  ecc_ratio      sqrt(xB^2 + yB^2) / c from the final state               [-]
  bearing_angle  atan2 of the bearing centre offset                       [deg]
  history        (time, xB, yB, W, p_max) over the whole VTK series
Diagnostics flags: max outer_iters reached, non-converged steps, worst
|liquid_residual|, cavitated_fraction drift.
"""

import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
CLEARANCE = 1.52e-4  # c [m], frozen Ferron hot clearance; same for every run


def _cfg_value(cfg_path, key, cast=float):
    txt = cfg_path.read_text()
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", txt, re.M)
    return cast(m.group(1)) if m else None


def _piece_arrays(vts_path):
    """Return {array_name: flat list of floats} for the CellData of one .vts."""
    root = ET.parse(vts_path).getroot()
    out = {}
    for cd in root.iter("CellData"):
        for da in cd.findall("DataArray"):
            name = da.get("Name")
            ncomp = int(da.get("NumberOfComponents", "1"))
            vals = [float(v) for v in (da.text or "").split()]
            out[name] = (vals, ncomp)
    return out


def _read_pvts(pvts_path):
    """Merge every piece of a .pvts. Returns {name: (values, ncomp)}."""
    root = ET.parse(pvts_path).getroot()
    merged = {}
    for piece in root.iter("Piece"):
        src = pvts_path.parent / piece.get("Source")
        for name, (vals, ncomp) in _piece_arrays(src).items():
            if name not in merged:
                merged[name] = ([], ncomp)
            merged[name][0].extend(vals)
    return merged


def _series(run_dir):
    """[(time, pvts_path)] from results/flat/results.pvd, time-ordered."""
    pvd = run_dir / "results" / "flat" / "results.pvd"
    if not pvd.exists():
        return []
    root = ET.parse(pvd).getroot()
    items = []
    for ds in root.iter("DataSet"):
        items.append((float(ds.get("timestep")), pvd.parent / ds.get("file")))
    items.sort(key=lambda x: x[0])
    return items


def _uniform(arrs, name, comp=0):
    """Value of a field the solver fills uniformly (forces, torque, position)."""
    if name not in arrs:
        return None
    vals, ncomp = arrs[name]
    return vals[comp] if ncomp == 1 else vals[comp]


def _frame_metrics(pvts_path):
    a = _read_pvts(pvts_path)
    m = {}
    if "Fp" in a:
        vals, ncomp = a["Fp"]
        fx, fy = vals[0], vals[1]
        m["Fx"], m["Fy"] = fx, fy
        m["W"] = math.hypot(fx, fy)
        m["force_angle_deg"] = math.degrees(math.atan2(fy, fx))
    if "F" in a:
        vals, _ = a["F"]
        m["F_fluid_x"], m["F_fluid_y"] = vals[0], vals[1]
        m["W_fluid"] = math.hypot(vals[0], vals[1])
    if "p" in a:
        p = a["p"][0]
        m["p_max"] = max(p)
        m["p_min"] = min(p)
    if "film_content" in a:
        fc = a["film_content"][0]
        m["cav_frac_0999"] = sum(1 for v in fc if v < 0.999) / len(fc)
        m["cav_frac_1"] = sum(1 for v in fc if v < 1.0) / len(fc)
        m["film_content_min"] = min(fc)
    if "friction_torque" in a:
        m["friction_torque"] = a["friction_torque"][0][0]
    if "h" in a:
        m["h_min"] = min(a["h"][0])
    if "T" in a:
        T = a["T"][0]
        m["T_max"] = max(T)
        m["T_mean"] = sum(T) / len(T)
        m["T_min"] = min(T)
    if "xB" in a:
        vals, _ = a["xB"]
        xb, yb = vals[0], vals[1]
        m["bearing_x"], m["bearing_y"] = xb, yb
        m["ecc_ratio"] = math.hypot(xb, yb) / CLEARANCE
        m["bearing_angle_deg"] = math.degrees(math.atan2(yb, xb))
    if "bearing_attitude_angle" in a:
        m["bearing_attitude_angle"] = a["bearing_attitude_angle"][0][0]
    return m


def _diagnostics(run_dir):
    csv_path = run_dir / "results" / "diagnostics.csv"
    if not csv_path.exists():
        return {}
    lines = [l for l in csv_path.read_text().splitlines() if l.strip()]
    head = lines[0].split(",")
    rows = [dict(zip(head, l.split(","))) for l in lines[1:]]
    if not rows:
        return {}
    outer = [int(r["outer_iters"]) for r in rows]
    conv = [int(r["converged"]) for r in rows]
    res = [abs(float(r["liquid_residual"])) for r in rows]
    mass = [float(r["liquid_mass"]) for r in rows]
    cav = [float(r["cavitated_fraction"]) for r in rows]
    return {
        "n_rows": len(rows),
        "outer_iters_max": max(outer),
        "outer_iters_final": outer[-1],
        "n_not_converged": sum(1 for c in conv if c == 0),
        "liquid_residual_max": max(res),
        "liquid_residual_final": res[-1],
        "liquid_mass_first": mass[0],
        "liquid_mass_final": mass[-1],
        "liquid_mass_rel_drift": (mass[-1] - mass[0]) / mass[0] if mass[0] else None,
        "cavitated_fraction_final": cav[-1],
        "final_time": float(rows[-1]["time"]),
        "final_step": int(rows[-1]["step"]),
    }


RUNS = ["grid_G1", "grid_G2", "grid_G3", "grid_G4", "grid_G5",
        "grid_G3_thd", "grid_G4_thd",
        "dt_T1", "dt_T2", "dt_T3", "dt_T4", "dt_T5"]


def main():
    out = {}
    for name in RUNS:
        run_dir = CASE_DIR / name
        if not run_dir.exists():
            continue
        cfg = run_dir / "config.txt"
        entry = {
            "n_theta": _cfg_value(cfg, "n_theta_global", int),
            "n_z": _cfg_value(cfg, "n_z_global", int),
            "dt": _cfg_value(cfg, "dt"),
            "end_t": _cfg_value(cfg, "end_t"),
            "solution_mode": _cfg_value(cfg, "solution_mode", str),
            "temperature_model": _cfg_value(cfg, "temperature_model", str),
            "max_outer_iters": _cfg_value(cfg, "max_outer_iters", int),
        }
        series = _series(run_dir)
        if series:
            entry["n_frames"] = len(series)
            entry["final_time_vtk"] = series[-1][0]
            entry["final"] = _frame_metrics(series[-1][1])
            hist = []
            for t, path in series:
                fm = _frame_metrics(path)
                hist.append({
                    "t": t,
                    "ecc_ratio": fm.get("ecc_ratio"),
                    "bearing_x": fm.get("bearing_x"),
                    "bearing_y": fm.get("bearing_y"),
                    "bearing_angle_deg": fm.get("bearing_angle_deg"),
                    "W": fm.get("W"),
                    "p_max": fm.get("p_max"),
                    "h_min": fm.get("h_min"),
                })
            entry["history"] = hist
        entry["diagnostics"] = _diagnostics(run_dir)
        out[name] = entry
        print(f"[ok] {name}", flush=True)
    (CASE_DIR / "metrics.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
