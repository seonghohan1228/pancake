"""Run three pilots, the remaining fourteen cases, and three focused checks."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if (root / "execution.json").exists():
        raise FileExistsError("Refusing to replace an existing campaign")
    rankfiles = []
    for slot in range(6):
        socket, first = (0 if slot < 3 else 1), (slot % 3)*4
        path = root / f"slot_{slot}.rankfile"
        with path.open("x", newline="\n") as stream:
            stream.writelines(f"rank {rank}=localhost slot={socket}:{first+rank}\n" for rank in range(4))
        rankfiles.append(path)
    worker = root / "source/cases/02e_lge_mixture/run_mixture.py"
    record = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
              "jobs": {}, "phases": [], "mpi_ranks_per_case": 4, "maximum_physical_cores": 24}

    def save():
        temporary = root / "execution.json.tmp"
        temporary.write_text(json.dumps(record, indent=2) + "\n")
        temporary.replace(root / "execution.json")

    def group(phase, cases):
        record["phases"].append(phase)
        pending = list(cases)
        active = {}
        failed = False
        while pending or active:
            for slot in range(6):
                if slot in active or not pending:
                    continue
                case = pending.pop(0)
                key = case if phase in ("pilots", "production") else f"{phase}_{case}"
                work = root / ("production" if phase in ("pilots", "production") else phase)
                command = [sys.executable, str(worker), "--inputs", str(root / "inputs"),
                           "--case", case, "--work", str(work), "--solver", str(root / "source/build-linux/pancake"),
                           "--rankfile", str(rankfiles[slot])]
                if phase.startswith("refinement"):
                    command += ["--tight"]
                if phase == "refinement_layers32":
                    command += ["--layers", "32"]
                log = (root / f"{key}_driver.log").open("x")
                started = time.time()
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=root)
                active[slot] = (process, log, started, key)
                record["jobs"][key] = {"state": "running", "slot": slot, "command": command,
                    "start_time": started, "pid": process.pid, "returncode": None,
                    "rankfile": str(rankfiles[slot]), "phase": phase}
                save()
            for slot, (process, log, started, key) in list(active.items()):
                result = process.poll()
                if result is None:
                    continue
                log.close()
                record["jobs"][key].update(state="complete" if result == 0 else "failed",
                    returncode=result, elapsed_seconds=time.time()-started)
                failed |= result != 0
                del active[slot]
                save()
            if active:
                time.sleep(10)
        return not failed

    save()
    pilots = ["M02", "M11", "M16"]
    remaining = [f"M{i:02d}" for i in range(1, 18) if f"M{i:02d}" not in pilots]
    success = group("pilots", pilots)
    if success:
        success = group("production", remaining)
    if success:
        success = group("refinement", ["M02", "M16"])
    if success:
        success = group("refinement_layers32", ["M02"])
    record.update(status="complete" if success else "failed_needs_diagnosis",
                  finished_utc=datetime.now(timezone.utc).isoformat())
    save()
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
