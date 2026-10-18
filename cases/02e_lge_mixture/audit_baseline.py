"""Compare the accepted source snapshot and configurations without rerunning CFD."""

import argparse
import difflib
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cases/02d_fw68d_steady_six"))
import run_steady as driver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    baseline = args.baseline_root / "resolved_results"
    snapshot = args.baseline_root / "resolved_snapshot"
    manifest = json.loads((snapshot / "source_manifest.json").read_text())
    archive = snapshot / "source.tar.gz"
    assert driver.sha256(archive) == manifest["sha256"]
    differences = []
    with tarfile.open(archive) as packed:
        files = {m.name.removeprefix("source/"): m for m in packed if m.isfile()}
        native = {p.relative_to(ROOT).as_posix() for p in (ROOT / "src").rglob("*")
                  if p.suffix in (".cpp", ".hpp")}
        native.add("CMakeLists.txt")
        for name in sorted(native):
            old = packed.extractfile(files[name]).read().replace(b"\r\n", b"\n") if name in files else b""
            current = (ROOT / name).read_bytes().replace(b"\r\n", b"\n")
            if old != current:
                differences.append({"path": name, "diff": "".join(difflib.unified_diff(
                    old.decode().splitlines(True), current.decode().splitlines(True)))})
        name = "cases/02d_fw68d_steady_six/run_steady.py"
        old_driver = packed.extractfile(files[name]).read()
    old = {"__name__": "archived_driver", "__file__": str(ROOT / name)}
    exec(compile(old_driver, "archived_driver", "exec"), old)
    experiment = json.loads(driver.EXPERIMENT.read_text())
    setup = json.loads(driver.SETUP.read_text())
    cases = {}
    for case in experiment["cases"]:
        key = case["case_id"]
        config_path = baseline / f"{key}_config.json"
        config = json.loads(config_path.read_text())
        arguments = (case, experiment, setup["specified_boundaries"], setup["material_settings"],
                     float(config["e"]) / float(config["c"]), float(config["attitude_angle_deg"]),
                     float(config["bush_outer_heat_transfer"]), float(config["shaft_ambient_heat_transfer"]),
                     Path("audit-output"), driver.TABLE)
        rendered = driver.render_config(*arguments, coupling_max_iters=int(config["coupling_max_iters"]))
        assert rendered == old["render_config"](*arguments, coupling_max_iters=int(config["coupling_max_iters"]))
        values = dict(line.split(" = ", 1) for line in rendered.splitlines() if " = " in line)
        ignored = {"e", "output_dir", "solution_table_file", "steady_initial_state_file"}
        mismatch = {k: [v, values.get(k)] for k, v in config.items() if k not in ignored and v != values.get(k)}
        assert not mismatch, (key, mismatch)
        cases[key] = {"config_sha256": driver.sha256(config_path),
                      "fields_sha256": driver.sha256(baseline / f"{key}_fields.npz"),
                      "render_byte_identical": True, "config_differences": mismatch}
    result = {"status": "accepted_baseline_identified_not_new_predictions",
              "current_base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "accepted_snapshot_base_commit": manifest["base_commit"],
              "accepted_snapshot_state": manifest["state"],
              "accepted_source_archive_sha256": driver.sha256(archive),
              "accepted_executable_sha256": (snapshot / "executable.sha256").read_text().split()[0],
              "native_differences": differences, "cases": cases,
              "note": "The archive includes verified uncommitted changes; its base Git commit alone does not reproduce the accepted baseline."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"source_archive_verified": True, "cases_verified": list(cases),
                      "native_difference_paths": [d["path"] for d in differences]}, indent=2))


if __name__ == "__main__":
    main()
