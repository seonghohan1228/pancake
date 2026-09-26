"""Bounded MPI scaling and concurrent-throughput benchmark for one fixed case."""

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
import platform
import re
import socket
import statistics
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def set_config(template, changes):
    lines = template.splitlines()
    for key, value in changes.items():
        pattern = re.compile(r"^\s*" + re.escape(key) + r"\s*=")
        matches = [i for i, line in enumerate(lines) if pattern.match(line)]
        if len(matches) > 1:
            raise ValueError(f"duplicate config key: {key}")
        if matches:
            lines[matches[0]] = f"{key} = {value}"
        else:
            lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def config_value(template, key):
    matches = re.findall(r"^\s*" + re.escape(key) + r"\s*=\s*([^#\s]+)", template, re.M)
    if len(matches) != 1:
        raise ValueError(f"expected one {key} in template")
    return float(matches[0])


def absolute_table(template, source_dir):
    match = re.search(r"^\s*solution_table_file\s*=\s*(\S+)", template, re.M)
    if not match:
        return template
    table = Path(match.group(1))
    if not table.is_absolute():
        table = source_dir / table
    table = table.resolve()
    if not table.is_file():
        raise FileNotFoundError(table)
    return set_config(template, {"solution_table_file": table.as_posix()})


def final_fields(results):
    pvd = results / "flat" / "results.pvd"
    if not pvd.exists():
        pvd = results / "results.pvd"
    if not pvd.exists():
        return {"time_s": None, "fields": {}}
    frames = ET.parse(pvd).findall(".//DataSet")
    if not frames:
        return {"time_s": None, "fields": {}}
    last = max(frames, key=lambda item: float(item.get("timestep")))
    pvts = pvd.parent / last.get("file")
    values = {}
    for piece in ET.parse(pvts).findall(".//Piece"):
        source = pvts.parent / piece.get("Source")
        for array in ET.parse(source).findall(".//CellData/DataArray"):
            name = array.get("Name")
            if name in {"p", "T", "h", "film_content", "F", "xB",
                        "T_bearing", "T_journal", "friction_torque"}:
                numbers = [float(v) for v in (array.text or "").split()]
                if numbers:
                    components = int(array.get("NumberOfComponents", "1"))
                    if len(numbers) % components:
                        raise ValueError(f"truncated VTK field {name}")
                    for component in range(components):
                        key = name if components == 1 else f"{name}[{component}]"
                        values.setdefault(key, []).extend(numbers[component::components])
    fields = {name: {"n": len(a), "sum": sum(a), "min": min(a), "max": max(a)}
              for name, a in values.items()}
    return {"time_s": float(last.get("timestep")), "fields": fields}


def diagnostics(results):
    path = results / "diagnostics.csv"
    if not path.exists():
        return {"rows": 0, "all_converged": False}
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {"rows": 0, "all_converged": False}
    residuals = {}
    for key in ("liquid_residual", "gas_residual", "energy_residual_W"):
        if key in rows[0]:
            numbers = [abs(float(row[key])) for row in rows if row[key]]
            if numbers:
                residuals[key] = max(numbers)
    return {"rows": len(rows), "last_time": rows[-1].get("time"),
            "all_converged": ("converged" in rows[0]
                              and all(row["converged"] == "1" for row in rows)
                              and all(row.get("coupling_converged", "") in ("", "1") for row in rows)),
            "max_abs_residual": residuals}


def make_run(work, template, dt, steps, name):
    run = work / name
    run.mkdir(parents=True, exist_ok=False)
    duration = dt * steps
    config = set_config(template, {"end_t": f"{duration:.12g}",
                                   "write_interval": f"{duration:.12g}",
                                   "diagnostics_interval": "1",
                                   "output_write_flat": "true",
                                   "output_write_3d": "false",
                                   "output_dir": (run / "results").as_posix()})
    path = run / "config.txt"
    path.write_text(config, encoding="utf-8")
    return run, path


def command(solver, config, ranks, cpus, rankfile):
    if ranks < 1 or ranks > 28 or len(cpus) != ranks or len(set(cpus)) != ranks:
        raise ValueError("ranks need distinct physical CPU IDs within 0-27")
    if min(cpus) < 0 or max(cpus) > 27:
        raise ValueError("CPU set must use physical CPU IDs 0-27")
    return ["mpirun", "--rankfile", str(rankfile), "--report-bindings", "-n", str(ranks),
            str(solver), "--config", str(config)]


def launch(solver, config, run, ranks, cpus):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    rankfile = run / "rankfile"
    rankfile.write_text("".join(f"rank {i}={socket.gethostname()} slot={cpu}\n"
                                for i, cpu in enumerate(cpus)), encoding="utf-8")
    cmd = command(solver, config, ranks, cpus, rankfile)
    log = (run / "run.log").open("w", encoding="utf-8")
    start = time.perf_counter()
    try:
        process = subprocess.Popen(cmd, cwd=run, stdout=log, stderr=subprocess.STDOUT, env=env)
    except Exception:
        log.close()
        raise
    return {"process": process, "log": log, "start": start, "run": run,
            "config": config, "ranks": ranks, "cpus": cpus, "command": cmd}


def finish(job, timeout):
    if "end" not in job:
        wait_job(job, timeout)
    run = job["run"]
    return {"name": run.name, "ranks": job["ranks"], "cpus": job["cpus"],
            "command": job["command"], "config_sha256": sha256(job["config"]),
            "return_code": job["rc"], "wall_s": job["end"] - job["start"],
            "core_s": job["ranks"] * (job["end"] - job["start"]),
            "diagnostics": diagnostics(run / "results"),
            "final_fields": final_fields(run / "results") if job["rc"] == 0 else {"time_s": None, "fields": {}}}


def wait_job(job, timeout):
    try:
        rc = job["process"].wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        job["process"].terminate()
        try:
            job["process"].wait(timeout=10)
        except subprocess.TimeoutExpired:
            job["process"].kill()
            job["process"].wait()
        rc = 124
    job["end"] = time.perf_counter()
    job["rc"] = rc
    job["log"].close()


def valid(rec):
    fields = rec["final_fields"]
    numbers = [v for summary in fields["fields"].values() for key, v in summary.items() if key != "n"]
    required = {"p", "T", "h", "film_content", "friction_torque",
                "T_bearing", "T_journal", "F[0]", "F[1]", "F[2]",
                "xB[0]", "xB[1]", "xB[2]"}
    return (rec["return_code"] == 0 and rec["diagnostics"]["all_converged"]
            and rec["diagnostics"]["rows"] == rec["steps"]
            and fields["time_s"] is not None
            and math.isclose(fields["time_s"], rec["expected_time_s"], rel_tol=1e-7, abs_tol=1e-10)
            and required <= fields["fields"].keys()
            and all(fields["fields"][name]["n"] == rec["expected_cells"] for name in required)
            and fields["fields"]["h"]["min"] > 0
            and fields["fields"]["T"]["min"] > 0
            and fields["fields"]["film_content"]["min"] >= 0
            and all(math.isfinite(x) for x in numbers))


def summarize(records):
    paired = {}
    for rec in records:
        if rec.get("mode") != "scale":
            continue
        paired.setdefault((rec["ranks"], rec["repeat"]), {})[rec["steps"]] = rec
    summaries = []
    for (ranks, repeat), trials in sorted(paired.items()):
        if len(trials) != 2:
            continue
        short, long = sorted(trials)
        a, b = trials[short], trials[long]
        passed = all(valid(x) for x in (a, b))
        rate = (b["wall_s"] - a["wall_s"]) / (long - short)
        summaries.append({"ranks": ranks, "repeat": repeat, "valid": passed and rate > 0,
                          "solver_s_per_step": rate, "core_s_per_step": ranks * rate,
                          "startup_s_estimate": a["wall_s"] - short * rate})
    return summaries


def self_check():
    template = "dt = 0.1\nend_t = 1\noutput_dir = old\n"
    changed = set_config(template, {"end_t": "2", "output_dir": "/tmp/new"})
    assert changed.count("end_t =") == 1 and config_value(changed, "end_t") == 2
    assert "output_dir = /tmp/new" in changed
    try:
        set_config("dt=1\ndt=2", {"dt": "3"})
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate key was accepted")
    rows = [{"mode": "scale", "ranks": 4, "repeat": 0, "steps": n,
             "wall_s": t, "return_code": 0, "expected_time_s": n * 0.1, "expected_cells": 1,
             "diagnostics": {"all_converged": True, "rows": n},
             "final_fields": {"time_s": n * 0.1,
                              "fields": {name: {"n": 1, "sum": 1., "min": 1., "max": 1.}
                                         for name in ("p", "T", "h", "film_content", "friction_torque",
                                                      "T_bearing", "T_journal", "F[0]", "F[1]", "F[2]",
                                                      "xB[0]", "xB[1]", "xB[2]")}}}
            for n, t in ((5, 3.0), (25, 11.0))]
    assert abs(summarize(rows)[0]["solver_s_per_step"] - 0.4) < 1e-12
    assert command("solver", "config", 2, [0, 14], "rankfile")[2] == "rankfile"
    print("self-check passed")


def tool_version(command):
    return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).splitlines()[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--source-sha", default="unknown")
    parser.add_argument("--mode", choices=("scale", "throughput"), default="scale")
    parser.add_argument("--ranks", default="1,2,4,7,14,24,28")
    parser.add_argument("--groups", default="7,7,7,7")
    parser.add_argument("--short-steps", type=int, default=2)
    parser.add_argument("--long-steps", type=int, default=6)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=600)
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return 0
    if not args.solver or not args.template or not args.work or not args.source_dir:
        parser.error("--solver, --template, --source-dir, and --work are required")
    if not (0 < args.short_steps < args.long_steps and args.repeats > 0 and args.timeout > 0):
        parser.error("invalid steps, repeats, or timeout")
    solver, template_path = args.solver.resolve(), args.template.resolve()
    template = absolute_table(template_path.read_text(encoding="utf-8"), args.source_dir.resolve())
    dt = config_value(template, "dt")
    if dt <= 0:
        parser.error("dt must be positive")
    args.work = args.work.resolve()
    args.work.mkdir(parents=True, exist_ok=False)
    manifest = {"host": socket.gethostname(), "platform": platform.platform(),
                "python": sys.version, "source_sha": args.source_sha,
                "solver_sha256": sha256(solver), "template_sha256": sha256(template_path),
                "solution_table_sha256": sha256(args.source_dir.resolve() / "data/R410A_FW68D/fw68d.csv"),
                "gcc": tool_version(["g++", "--version"]),
                "cmake": tool_version(["cmake", "--version"]),
                "open_mpi": tool_version(["mpirun", "--version"]),
                "petsc": tool_version(["pkg-config", "--modversion", "PETSc"]),
                "build_type": next((line.split("=", 1)[1] for line in
                                    (solver.parent / "CMakeCache.txt").read_text().splitlines()
                                    if line.startswith("CMAKE_BUILD_TYPE:STRING=")), "unknown"),
                "mode": args.mode, "dt": dt, "runs": []}
    expected_cells = int(config_value(template, "n_theta_global") * config_value(template, "n_z_global"))
    output = args.work / "manifest.json"
    def save():
        manifest["paired_scaling"] = summarize(manifest["runs"])
        output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    save()
    if args.mode == "scale":
        ranks = [int(x) for x in args.ranks.split(",")]
        if not ranks or len(set(ranks)) != len(ranks) or any(n < 1 or n > 28 for n in ranks):
            parser.error("ranks must be distinct integers in 1..28")
        for repeat in range(args.repeats):
            for n in ranks:
                for steps in (args.short_steps, args.long_steps):
                    name = f"scale_n{n}_s{steps}_r{repeat}"
                    run, cfg = make_run(args.work, template, dt, steps, name)
                    rec = finish(launch(solver, cfg, run, n, list(range(n))), args.timeout)
                    rec.update(mode="scale", repeat=repeat, steps=steps, expected_time_s=dt * steps,
                               expected_cells=expected_cells)
                    manifest["runs"].append(rec)
                    save()
                    print(json.dumps({"name": name, "rc": rec["return_code"],
                                      "wall_s": rec["wall_s"], "valid": valid(rec)}), flush=True)
                    if not valid(rec):
                        return 1
    else:
        groups = [int(x) for x in args.groups.split(",")]
        if not groups or any(n < 1 for n in groups) or sum(groups) > 28:
            parser.error("groups must be positive ranks totaling at most 28")
        jobs = []
        offset = 0
        for i, n in enumerate(groups):
            name = f"throughput_g{i}_n{n}"
            run, cfg = make_run(args.work, template, dt, args.long_steps, name)
            jobs.append(launch(solver, cfg, run, n, list(range(offset, offset + n))))
            offset += n
        batch_start = min(job["start"] for job in jobs)
        with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
            futures = [pool.submit(wait_job, job, args.timeout) for job in jobs]
            for future in futures:
                future.result()
        for job in jobs:
            rec = finish(job, args.timeout)
            rec.update(mode="throughput", steps=args.long_steps, expected_time_s=dt * args.long_steps,
                       expected_cells=expected_cells)
            manifest["runs"].append(rec)
            save()
        batch_wall = max(job["end"] for job in jobs) - batch_start
        manifest["throughput"] = {"groups": groups, "batch_wall_s": batch_wall,
                                  "case_steps_per_s": len(groups) * args.long_steps / batch_wall,
                                  "all_valid": all(valid(r) for r in manifest["runs"])}
        save()
        print(json.dumps(manifest["throughput"]), flush=True)
    return 0 if all(valid(r) for r in manifest["runs"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
