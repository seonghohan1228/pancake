"""Verify a synthetic oil-air transient across MPI counts; not a calibrated case.

Requires NumPy. Paths may be native Windows or Linux. All outputs stay below
--output. Numerical pass criteria do not certify gas pressure-volume closure.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

import numpy as np


def final_fields(case):
    collection = ET.parse(case / "results/flat/results.pvd").getroot()
    dataset = list(collection.iter("DataSet"))[-1]
    parallel = case / "results/flat" / dataset.get("file")
    parts = []
    for piece in ET.parse(parallel).getroot().iter("Piece"):
        data = ET.parse(parallel.parent / piece.get("Source")).getroot().find(".//Piece")
        extent = list(map(int, data.get("Extent").split()))
        nx, nz = extent[1] - extent[0], extent[3] - extent[2]
        arrays = {d.get("Name"): np.fromstring(d.text, sep=" ").reshape(nz, nx)
                  for d in data.find("CellData").iter("DataArray")
                  if int(d.get("NumberOfComponents", "1")) == 1}
        parts.append((extent[0], arrays))
    parts.sort(key=lambda part: part[0])
    return {key: np.concatenate([part[1][key] for part in parts], axis=1)
            for key in parts[0][1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    solver, output = args.solver.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    template = (Path(__file__).parent / "air_smoke.txt").read_text(encoding="utf-8")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    reference = None
    runs = []
    for ranks in (1, 2, 4):
        case = output / f"mpi{ranks}"
        case.mkdir(exist_ok=True)
        config = re.sub(r"(?m)^output_dir\s*=.*$",
                        f"output_dir = {(case / 'results').as_posix()}", template)
        (case / "config.txt").write_text(config, encoding="utf-8")
        process = subprocess.run([args.mpiexec, "-n", str(ranks), str(solver),
                                  "-c", str(case / "config.txt")], cwd=case,
                                 env=env, capture_output=True, timeout=120)
        (case / "solver.log").write_bytes(process.stdout + process.stderr)
        if process.returncode:
            raise RuntimeError(f"MPI-{ranks} exited {process.returncode}; see {case / 'solver.log'}")
        with (case / "results/diagnostics.csv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        arrays = final_fields(case)
        if reference is None:
            reference = arrays
        difference = {key: float(np.max(abs(value - reference[key])) /
                                max(1e-30, np.max(abs(reference[key]))))
                      for key, value in arrays.items()}
        gas_residual = max(abs(float(row["gas_residual"])) for row in rows)
        run = {"ranks": ranks, "diagnostic_rows": len(rows),
               "max_abs_gas_residual": gas_residual,
               "max_gas_volume_excess": max(float(row["gas_volume_excess"]) for row in rows),
               "max_gas_eos_defect": max(float(row["gas_eos_defect"]) for row in rows),
               "final_ranges": {key: [float(value.min()), float(value.max())]
                                for key, value in arrays.items()},
               "relative_difference_vs_mpi1_at_vtk_precision": difference,
               "numerical_pass": gas_residual < 1e-8 and max(difference.values()) < 1e-5
                                 and all(row["converged"] == "1" for row in rows)}
        runs.append(run)
        print(json.dumps(run), flush=True)
    result = {"solver_sha256": hashlib.sha256(solver.read_bytes()).hexdigest(),
              "scope": "Synthetic transient regression; volume errors are reported, not certified",
              "runs": runs, "numerical_pass": all(run["numerical_pass"] for run in runs)}
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0 if result["numerical_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
