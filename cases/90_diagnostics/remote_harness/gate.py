"""Standalone sanity-gate extractor: prints W, p_max and friends for one run dir.

Usage: python gate.py <run_dir>   (run_dir contains results/flat/results.pvd)
Mirrors cases/00_grid_dt_independence/extract.py field conventions.
"""
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

CLEARANCE = 1.52e-4


def piece_arrays(vts_path):
    root = ET.parse(vts_path).getroot()
    out = {}
    for cd in root.iter("CellData"):
        for da in cd.findall("DataArray"):
            ncomp = int(da.get("NumberOfComponents", "1"))
            out[da.get("Name")] = ([float(v) for v in (da.text or "").split()], ncomp)
    return out


def read_pvts(pvts_path):
    root = ET.parse(pvts_path).getroot()
    merged = {}
    for piece in root.iter("Piece"):
        src = pvts_path.parent / piece.get("Source")
        for name, (vals, ncomp) in piece_arrays(src).items():
            merged.setdefault(name, ([], ncomp))[0].extend(vals)
    return merged


def series(run_dir):
    pvd = Path(run_dir) / "results" / "flat" / "results.pvd"
    root = ET.parse(pvd).getroot()
    items = [(float(ds.get("timestep")), pvd.parent / ds.get("file"))
             for ds in root.iter("DataSet")]
    items.sort(key=lambda x: x[0])
    return items


def frame(pvts_path):
    a = read_pvts(pvts_path)
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
        m["cav_frac_0999"] = sum(1 for v in fc if v < 0.999) / len(fc)
        m["film_content_min"] = min(fc)
    if "friction_torque" in a:
        m["friction_torque"] = a["friction_torque"][0][0]
    if "h" in a:
        m["h_min"] = min(a["h"][0])
    if "xB" in a:
        xb, yb = a["xB"][0][0], a["xB"][0][1]
        m["bearing_x"], m["bearing_y"] = xb, yb
        m["ecc_ratio"] = math.hypot(xb, yb) / CLEARANCE
        m["bearing_angle_deg"] = math.degrees(math.atan2(yb, xb))
    return m


if __name__ == "__main__":
    run_dir = sys.argv[1]
    s = series(run_dir)
    t, path = s[-1]
    m = frame(path)
    print(f"run_dir   {run_dir}")
    print(f"n_frames  {len(s)}")
    print(f"time      {t:.6g}")
    for k in sorted(m):
        print(f"{k:20s} {m[k]!r}")
