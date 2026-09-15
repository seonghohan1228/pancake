"""Queue the 0.10 s timestep ladder behind the independently running G8 case."""

import json
import subprocess
import sys
import time
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
REPO_DIR = CASE_DIR.parent.parent


def run(script, *arguments):
    command = [sys.executable, str(CASE_DIR / script), *arguments]
    subprocess.run(command, cwd=REPO_DIR, check=True)


def main():
    ledger = CASE_DIR / "campaign_equilibrium.json"
    while True:
        data = json.loads(ledger.read_text()) if ledger.exists() else {}
        result = data.get("d_grid_g8")
        if result and "failed" in result:
            raise RuntimeError(f"G8 failed: {result['failed']}")
        if result and "flow" in result:
            break
        print("waiting for G8", flush=True)
        time.sleep(30)
    run("analyse_campaign.py")
    run("export_visual_data.py")
    run("run_timestep.py", "--extended-100ms")
    run("analyse_discretization.py", "--extended-100ms")
    run("export_timestep_history.py", "--extended-100ms")


if __name__ == "__main__":
    main()
