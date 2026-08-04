"""Sequential runner for case 00 (grid / timestep independence).

Runs each case directory with `mpiexec -n 4 pancake.exe --config <cfg>` from
inside that directory, records wall-clock time, and appends the result to
run_times.json. Cases already recorded are skipped unless --force is given.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
EXE = Path(r"C:\Dev\pancake\build-windows-mingw\pancake.exe")
MINGW_BIN = r"C:\msys64\mingw64\bin"
TIMES_FILE = CASE_DIR / "run_times.json"
PROGRESS_FILE = CASE_DIR / "progress.log"


def progress(message):
    """Append one timestamped line to progress.log and echo it."""
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp}  {message}"
    with open(PROGRESS_FILE, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)

DEFAULT_RUNS = [
    "grid_G1", "grid_G2", "grid_G3", "grid_G4", "grid_G5",
    "grid_G3_thd", "grid_G4_thd",
    "dt_T1", "dt_T2", "dt_T3", "dt_T4", "dt_T5",
]


def load_times():
    if TIMES_FILE.exists():
        return json.loads(TIMES_FILE.read_text())
    return {}


def save_times(times):
    TIMES_FILE.write_text(json.dumps(times, indent=2))


def run_one(name, cfg="config.txt", tag=None, log_name=None):
    tag = tag or name
    work = CASE_DIR / name
    log_name = log_name or ("run.log" if cfg == "config.txt" else f"run_{Path(cfg).stem}.log")
    env = dict(os.environ)
    env["PATH"] = MINGW_BIN + os.pathsep + env["PATH"]
    log_path = work / log_name
    progress(f"case {tag} STARTED  (mpiexec -n 4 pancake.exe --config {cfg})")
    t0 = time.perf_counter()
    # Tee: the solver's combined stdout/stderr goes to run.log inside the case
    # directory and to this runner's console at the same time.
    with open(log_path, "w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.Popen(
            ["mpiexec", "-n", "4", str(EXE), "--config", cfg],
            cwd=str(work), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            env=env, text=True, encoding="utf-8", errors="replace", bufsize=1,
        )
        for line in proc.stdout:
            log.write(line)
            log.flush()
        proc.wait()
    wall = time.perf_counter() - t0
    if proc.returncode == 0:
        progress(f"case {tag} FINISHED in {wall/60.0:.2f} min ({wall:.1f} s)")
    else:
        progress(f"case {tag} FAILED rc={proc.returncode} after {wall/60.0:.2f} min "
                 f"(see {work.name}/{log_name})")
    return {"returncode": proc.returncode, "wall_s": wall, "config": cfg, "log": log_name}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", default=None)
    ap.add_argument("--config", default="config.txt")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    runs = args.runs or DEFAULT_RUNS
    times = load_times()
    for name in runs:
        tag = args.tag or (name if args.config == "config.txt" else f"{name}:{Path(args.config).stem}")
        if tag in times and not args.force:
            progress(f"case {tag} SKIPPED (already recorded)")
            continue
        times[tag] = run_one(name, args.config, tag)
        save_times(times)
    progress(f"batch complete: {', '.join(runs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
