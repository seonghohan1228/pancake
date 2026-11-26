"""Run the frozen selection and its focused checks on 24 disjoint cores."""

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
    worker = root / "source/cases/02f_lge_mixture_gas/run_gas.py"
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
                case, branch = pending.pop(0)
                key = f"{phase}_{branch}_{case}"
                work = root / ("production" if phase in ("pilots", "production") else phase)
                command = [sys.executable, str(worker), "--inputs", str(root / "inputs"),
                           "--case", case, "--branch", branch, "--work", str(work), "--solver", str(root / "source/build-linux/pancake"),
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
                    "rankfile": str(rankfiles[slot]), "phase": phase, "case_id": case, "branch": branch}
                save()
            for slot, (process, log, started, key) in list(active.items()):
                result = process.poll()
                if result is None:
                    continue
                log.close()
                record["jobs"][key].update(state="complete" if result == 0 else "failed",
                    returncode=result, elapsed_seconds=time.time()-started)
                failed |= result != 0
                if result != 0:
                    record.setdefault("not_started", []).extend(
                        {"phase": phase, "case_id": case, "branch": branch} for case, branch in pending)
                    pending.clear()
                del active[slot]
                save()
            if active:
                time.sleep(10)
        return not failed

    save()
    experiment = json.loads((root / "inputs/experiment.json").read_text())
    selected = [case["case_id"] for case in experiment["cases"]]
    pilots = experiment.get("pilot_case_ids", ["M02", "M11", "M16"])
    checks = experiment.get("numerical_check_cases", {"refinement": ["M02", "M16"], "refinement_layers32": ["M02"]})
    assert len(set(selected)) == len(selected) and set(pilots) <= set(selected)
    assert all(set(keys) <= set(selected) for keys in checks.values())
    remaining = [key for key in selected if key not in pilots]
    pairs = lambda cases: [(case, branch) for case in cases for branch in ("gas_off", "gas_on")]
    success = group("pilots", pairs(pilots))
    if success:
        success = group("production", pairs(remaining))
    if success:
        success = group("refinement", pairs(checks["refinement"]))
    if success:
        success = group("refinement_layers32", pairs(checks["refinement_layers32"]))
    record.update(status="complete" if success else "failed_needs_diagnosis",
                  finished_utc=datetime.now(timezone.utc).isoformat())
    save()
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
