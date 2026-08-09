"""Extract lightweight bearing-eccentricity histories from remote VTK output."""

import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
CLEARANCE = 1.52e-4
RUNS = {
    "d_dt_tn2": 320.0,
    "d_dt_tn1": 160.0,
    "d_dt_t0": 80.0,
    "d_dt_t1": 40.0,
    "d_dt_t2": 20.0,
    "d_dt_t3": 10.0,
    "d_dt_t4": 5.0,
}
RUNS_100MS = {
    "d_dt100_t1": 40.0,
    "d_dt100_t2": 20.0,
    "d_dt100_t3": 10.0,
    "d_dt100_t4": 5.0,
}


def first_vector(path, name):
    blob = path.read_bytes()
    match = re.search(rb'Name="%s"[^>]*>' % name.encode(), blob)
    if not match:
        return None
    values = blob[match.end():match.end() + 400].split(b"<")[0].split()
    return [float(value) for value in values[:3]]


def series(run_dir):
    pvd = run_dir / "results" / "flat" / "results.pvd"
    root = ET.parse(pvd).getroot()
    rows = []
    for dataset in root.iter("DataSet"):
        pvts = pvd.parent / dataset.get("file")
        pvts_root = ET.parse(pvts).getroot()
        piece = next(pvts_root.iter("Piece"))
        x_b = first_vector(pvts.parent / piece.get("Source"), "xB")
        rows.append({
            "time_s": float(dataset.get("timestep")),
            "eps": math.hypot(x_b[0], x_b[1]) / CLEARANCE,
        })
    return sorted(rows, key=lambda row: row["time_s"])


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--extended-100ms", action="store_true")
    args = parser.parse_args()
    runs = RUNS_100MS if args.extended_100ms else RUNS
    output_path = CASE_DIR / ("d_timestep_100ms_history.json"
                              if args.extended_100ms else "d_timestep_history.json")
    output = {}
    for name, dt_us in runs.items():
        run_dir = CASE_DIR / "runs" / "production" / name
        if not (run_dir / "results" / "flat" / "results.pvd").exists():
            continue
        output[name] = {"dt_us": dt_us, "history": series(run_dir)}
        print(f"{name}: {len(output[name]['history'])} frames")
    output_path.write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
