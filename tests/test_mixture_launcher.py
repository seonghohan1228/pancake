"""Exercise real subprocess lifecycle and pilot gates with a disposable fake worker."""

import json
import os
from pathlib import Path
import runpy
import sys
from tempfile import TemporaryDirectory
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
launcher = runpy.run_path(str(ROOT / "cases/02e_lge_mixture/launch_campaign.py"))

def main():
    sleep = time.sleep
    with TemporaryDirectory() as temporary:
        for failure in ("", "M02"):
            root = Path(temporary) / ("failure" if failure else "success")
            worker = root / "source/cases/02e_lge_mixture/run_mixture.py"
            worker.parent.mkdir(parents=True)
            worker.write_text("import argparse,os,sys\n"
                "p=argparse.ArgumentParser();p.add_argument('--case')\n"
                "a,_=p.parse_known_args()\n"
                "sys.exit(3 if os.environ.get('PANCAKE_TEST_FAIL_CASE')==a.case else 0)\n")
            with patch.object(sys, "argv", ["launch_campaign.py", "--root", str(root)]), \
                 patch.object(time, "sleep", lambda _: sleep(.01)), \
                 patch.dict(os.environ, {"PANCAKE_TEST_FAIL_CASE": failure}):
                assert launcher["main"]() == (1 if failure else 0)
                saved = (root / "execution.json").read_bytes()
                record = json.loads(saved)
                jobs = record["jobs"]
                pilots = {"M02", "M11", "M16"}
                expected = pilots if failure else {f"M{i:02d}" for i in range(1, 18)} | {
                    "refinement_M02", "refinement_M16", "refinement_layers32_M02"}
                assert set(jobs) == expected
                assert record["status"] == ("failed_needs_diagnosis" if failure else "complete")
                assert record["phases"] == (["pilots"] if failure else
                    ["pilots", "production", "refinement", "refinement_layers32"])
                for key, job in jobs.items():
                    failed = bool(failure and key == failure)
                    assert job["returncode"] == (3 if failed else 0)
                    assert job["state"] == ("failed" if failed else "complete")
                    command = job["command"]
                    assert command[command.index("--case") + 1] == key.split("_")[-1]
                    assert ("--tight" in command) == key.startswith("refinement")
                    assert ("--layers" in command) == key.startswith("refinement_layers32")
                    if "--layers" in command:
                        assert command[command.index("--layers") + 1] == "32"
                intervals = sorted(jobs.values(), key=lambda job: job["start_time"])
                for index, job in enumerate(intervals):
                    for other in intervals[index + 1:]:
                        if other["slot"] == job["slot"]:
                            assert other["start_time"] >= job["start_time"] + job["elapsed_seconds"]
                if not failure:
                    for before, after in zip(record["phases"], record["phases"][1:]):
                        assert max(j["start_time"] + j["elapsed_seconds"] for j in jobs.values()
                                   if j["phase"] == before) <= min(j["start_time"] for j in jobs.values()
                                                                 if j["phase"] == after)
                cores = []
                for rankfile in sorted(root.glob("slot_*.rankfile")):
                    raw = rankfile.read_bytes()
                    assert b"\r" not in raw and raw.endswith(b"\n")
                    lines = raw.decode().splitlines()
                    assert len(lines) == 4
                    assert [line.split("=")[0] for line in lines] == [f"rank {i}" for i in range(4)]
                    cores.extend(line.split("slot=")[1] for line in lines)
                assert len(cores) == len(set(cores)) == 24
                try:
                    launcher["main"]()
                except FileExistsError:
                    pass
                else:
                    raise AssertionError("existing campaign execution was accepted")
                assert (root / "execution.json").read_bytes() == saved
    print("Mixture launcher passed: 20 jobs, pilot failure gate, distinct cores, LF rankfiles and preservation")

if __name__ == "__main__":
    main()
