"""Rerun the existing exact 1-D FBNS benchmark in an isolated output directory.

Usage: python run_verification.py --solver /path/to/pancake --output /path/to/run
Requires NumPy, as does the existing benchmark extractor. Runs sequentially.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    benchmark = root / "cases/00c_jfo_verification"
    solver = args.solver.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("front", benchmark / "extract_front.py")
    front = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(front)
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    runs = []
    for cells in (120, 240, 480, 960):
        case = output / f"n{cells}"
        case.mkdir(exist_ok=True)
        config = (benchmark / f"oneD_n{cells}_fbns/config.txt").read_text(encoding="utf-8")
        config = re.sub(r"(?m)^output_dir\s*=.*$", f"output_dir = {(case / 'results').as_posix()}", config)
        path = case / "config.txt"
        path.write_text(config, encoding="utf-8")
        started = time.perf_counter()
        process = subprocess.run([str(solver), "-c", str(path)], cwd=case,
                                 env=environment, capture_output=True, timeout=60)
        elapsed = time.perf_counter() - started
        (case / "solver.log").write_bytes(process.stdout + process.stderr)
        if process.returncode:
            raise RuntimeError(f"n={cells} failed with exit {process.returncode}; see {case / 'solver.log'}")
        result = front.analyse(case)
        result["wall_seconds"] = elapsed
        result["rupture_error_cells"] = abs(result["rupture"] - 302.766608) / result["cell_deg"]
        result["reformation_error_cells"] = abs(result["reformation"] - 79.006725) / result["cell_deg"]
        result["peak_error_percent"] = 100 * (result["p_max"] / 7.597141e6 - 1)
        result["load_error_percent"] = 100 * (result["load"] / 495463.6 - 1)
        result["pass"] = (result["rupture_error_cells"] <= 1 and result["reformation_error_cells"] <= 1
                          and abs(result["peak_error_percent"]) < 1 and abs(result["load_error_percent"]) < 1)
        runs.append(result)
        print(json.dumps(result), flush=True)
    summary = {"solver": str(solver), "solver_sha256": hashlib.sha256(solver.read_bytes()).hexdigest(),
               "scope": "Numerical verification of steady constant-property equilibrium JFO cavitation via FBNS",
               "runs": runs, "pass": all(run["pass"] for run in runs)}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0 if summary["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
