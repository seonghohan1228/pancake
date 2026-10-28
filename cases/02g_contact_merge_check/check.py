"""Merge gate for the October 2026 contact-model integration (codex/integrate-2026-10).

Prepares, runs and compares the three check sets in README.md. Nothing runs
unless asked: every run subcommand accepts --dry-run, which only writes the
configs and the manifest.

  python check.py af-smooth  --solver S --work W/af_smooth            # A-F, contact off, full closure
  python check.py af-contact --solver S --work W/af_contact --closed W/af_smooth --cases A,B
  python check.py mixture    --solver S --work W/mixture --base <density-off production root>
  python check.py compare    --work W

Every run directory gets manifest.json with the git commit, branch, dirty flag,
solver SHA-256, parameter file SHA-256 and the exact command. A dirty source
tree is refused unless --allow-dirty is given (the diff is then saved).
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
STEADY = REPO / "cases/02d_fw68d_steady_six"
EXPERIMENT = STEADY / "experiment_recalibrated_20261023.json"
SETUP_DEFAULT = STEADY / "steady_setup.json"
SETUP_TIGHT = HERE / "steady_setup_fbns1e-8.json"   # 2026-10-23: A-E used this (hash-checked), F the default
TIGHT_CASES = {"A", "B", "C", "D", "E"}
PARAMS = REPO / "cases/contact_parameters_sourced_20261022.json"
MIXTURE_CASES = ["M01", "M03", "M05", "M11", "M16"]
CONTACT_FIELDS = "contact_pressure,contact_force_x,contact_force_y,contact_force_z,contact_friction_torque"

sys.path.insert(0, str(STEADY))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True).stdout.strip()


def provenance(args, work):
    """Record the source revision; refuse a dirty tree unless --allow-dirty."""
    dirty = bool(git("status", "--porcelain", "--untracked-files=no"))
    record = {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
              "dirty": dirty, "created_utc": datetime.now(timezone.utc).isoformat(),
              "params_file": str(PARAMS), "params_sha256": sha256(PARAMS)}
    if getattr(args, "solver", None):
        record["solver"] = str(Path(args.solver).resolve())
        record["solver_sha256"] = sha256(args.solver)
    if dirty:
        if not args.allow_dirty:
            sys.exit("Source tree is dirty: commit first or pass --allow-dirty (the diff is then saved).")
        work.mkdir(parents=True, exist_ok=True)
        (work / "source.diff").write_text(git("diff", "HEAD"))
    return record


def read_config(path):
    return dict(line.split(" = ", 1) for line in Path(path).read_text().splitlines() if " = " in line)


def write_run(run_dir, config, manifest, command):
    run_dir.mkdir(parents=True, exist_ok=True)
    config = dict(config, output_dir=str(run_dir / "results"))
    (run_dir / "config.txt").write_text("".join(f"{key} = {value}\n" for key, value in config.items()))
    (run_dir / "manifest.json").write_text(json.dumps(dict(manifest, command=command), indent=2) + "\n")


def run_native(run_dir, args):
    command = [args.mpiexec, "-n", str(args.ranks), str(Path(args.solver).resolve()), "--config", "config.txt"]
    if args.dry_run:
        return None
    started = time.time()
    with open(run_dir / "run.log", "w") as log:
        code = subprocess.call(command, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    manifest.update(returncode=code, wall_s=round(time.time() - started))
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return code


def contact_variants(params):
    tilt = str(params["tilt_slope_y"])
    return {
        "smooth_aligned": {},
        "smooth_tilt": {"tilt_slope_y": tilt},
        "mixed_tilt_sourced": {"tilt_slope_y": tilt, **params["nominal"]},
        "mixed_tilt_upper": {"tilt_slope_y": tilt, **params["upper"]},
    }


def apply_variant(config, variant):
    config = dict(config)
    variant = {key: str(value) for key, value in variant.items()}
    if "tilt_slope_y" in variant:
        # Same treatment as the 2026-10-23 A/B check: keep the edge clear of the wall.
        c, length, e = float(config["c"]), float(config["L"]), float(config["e"])
        tilt = float(variant["tilt_slope_y"])
        config["e"] = str(min(e, 0.95 * (c - tilt * length / 2)))
        config.pop("steady_initial_state_file", None)
    if "asperity_contact_model" in variant and CONTACT_FIELDS not in config.get("output_fields", ""):
        config["output_fields"] = config.get("output_fields", "") + "," + CONTACT_FIELDS
    config.update(variant)
    return config


def af_smooth(args):
    work = Path(args.work).resolve()
    manifest = provenance(args, work)
    work.mkdir(parents=True, exist_ok=True)
    jobs = []
    for case in args.cases.split(","):
        setup = SETUP_TIGHT if case in TIGHT_CASES else SETUP_DEFAULT
        command = [sys.executable, str(STEADY / "run_steady.py"), "--solver", str(Path(args.solver).resolve()),
                   "--ranks", str(args.ranks), "--work", str(work), "--cases", case,
                   "--experiment", str(EXPERIMENT), "--materials", str(setup)]
        jobs.append({"case": case, "setup": str(setup), "setup_sha256": sha256(setup), "command": command})
    (work / "manifest.json").write_text(json.dumps(dict(manifest, jobs=jobs), indent=2) + "\n")
    for job in jobs:
        print(" ".join(job["command"]))
        if not args.dry_run:
            subprocess.call(job["command"], cwd=STEADY)


def af_contact(args):
    work = Path(args.work).resolve()
    manifest = provenance(args, work)
    params = json.loads(PARAMS.read_text())
    for case in args.cases.split(","):
        result_files = glob.glob(str(Path(args.closed) / case / "result_*.json"))
        if not result_files:
            sys.exit(f"No closed result for {case} under {args.closed}; run af-smooth first.")
        result = json.loads(Path(result_files[0]).read_text())
        if result.get("status") != "load_and_thermal_closed":
            sys.exit(f"{case}: smooth run not closed ({result.get('status')}).")
        base = read_config(Path(result["run_dir"]) / "config.txt")
        for tag, variant in contact_variants(params).items():
            run_dir = work / case / tag
            write_run(run_dir, apply_variant(base, variant),
                      dict(manifest, case=case, variant=tag, source_run=result["run_dir"],
                           smooth_closed_torque_Nm=result.get("torque_Nm")), "pending")
            print(run_dir, "rc =", run_native(run_dir, args))


def mixture(args):
    work = Path(args.work).resolve()
    manifest = provenance(args, work)
    params = json.loads(PARAMS.read_text())
    variants = {"smooth": {}, "sourced": params["nominal"], "upper": params["upper"]}
    for case in args.cases.split(","):
        candidates = sorted((Path(args.base) / case / "runs" / case).glob("0*"))
        if not candidates:
            sys.exit(f"No accepted run for {case} under {args.base}.")
        source = candidates[-1]
        base = read_config(source / "config.txt")
        for tag, variant in variants.items():
            run_dir = work / f"{case}_{tag}"
            write_run(run_dir, apply_variant(base, variant),
                      dict(manifest, case=case, variant=tag, source_run=str(source)), "pending")
            print(run_dir, "rc =", run_native(run_dir, args))


def read_result(run_dir):
    from run_steady import read_fields
    fields = read_fields(run_dir)
    row = {"torque_Nm": float(fields["friction_torque"][0]), "h_min_um": float(min(fields["h"])) * 1e6}
    row["contact_torque_Nm"] = float(fields["contact_friction_torque"][0]) if "contact_friction_torque" in fields else 0.0
    return row


def compare(args):
    work = Path(args.work).resolve()
    rows = []
    reference_af = {r["case"]: float(r["M_cfd_new"]) for r in csv.DictReader(open(HERE / "references/smooth_recalibrated_AF.csv"))}
    for result_file in sorted(glob.glob(str(work / "af_smooth/*/result_*.json"))):
        result = json.loads(Path(result_file).read_text())
        case = Path(result_file).parent.name
        rows.append({"set": "af_smooth", "case": case, "variant": "closed", "torque_Nm": result.get("torque_Nm"),
                     "reference_Nm": reference_af.get(case), "status": result.get("status")})
    reference_ab = {(r["case"], r["tag"]): float(r["M"]) for r in csv.DictReader(open(HERE / "references/mixed_lubrication_AB.csv"))}
    reference_ab_alias = {"mixed_tilt_sourced": "mixed_tilt_nominal"}
    for manifest_file in sorted(glob.glob(str(work / "af_contact/*/*/manifest.json"))):
        run_dir = Path(manifest_file).parent
        manifest = json.loads(Path(manifest_file).read_text())
        row = {"set": "af_contact", "case": manifest["case"], "variant": manifest["variant"],
               "reference_Nm": reference_ab.get((manifest["case"], reference_ab_alias.get(manifest["variant"], manifest["variant"])))}
        try:
            row.update(read_result(run_dir))
        except Exception as error:
            row["status"] = f"no result: {error!r}"[:120]
        rows.append(row)
    reference_mix = {(r["case"], r["variant"]): r["torque_Nm"] for r in
                     json.loads((HERE / "references/mixed_sensitivity_summary_20261023.json").read_text())}
    for manifest_file in sorted(glob.glob(str(work / "mixture/*/manifest.json"))):
        run_dir = Path(manifest_file).parent
        manifest = json.loads(Path(manifest_file).read_text())
        # Only the smooth control has an exact 2026-10-23 counterpart; s0.4 is the nearest rough one (fit F5/2, generic K).
        key = (manifest["case"], {"smooth": "smooth", "sourced": "s0.4"}.get(manifest["variant"], ""))
        row = {"set": "mixture", "case": manifest["case"], "variant": manifest["variant"], "reference_Nm": reference_mix.get(key)}
        try:
            row.update(read_result(run_dir))
        except Exception as error:
            row["status"] = f"no result: {error!r}"[:120]
        rows.append(row)
    for row in rows:
        if row.get("torque_Nm") is not None and row.get("reference_Nm"):
            row["change_pct"] = 100.0 * (row["torque_Nm"] / row["reference_Nm"] - 1.0)
    keys = ["set", "case", "variant", "torque_Nm", "reference_Nm", "change_pct", "h_min_um", "contact_torque_Nm", "status"]
    with open(work / "comparison.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print({key: row.get(key) for key in keys if row.get(key) is not None})
    print(f"wrote {work / 'comparison.csv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, default_cases in (("af-smooth", "A,B,C,D,E,F"), ("af-contact", "A,B"), ("mixture", ",".join(MIXTURE_CASES))):
        p = sub.add_parser(name)
        p.add_argument("--solver", required=True)
        p.add_argument("--work", required=True)
        p.add_argument("--cases", default=default_cases)
        p.add_argument("--ranks", type=int, default=2)
        p.add_argument("--mpiexec", default=shutil.which("mpiexec") or "mpiexec")
        p.add_argument("--timeout", type=int, default=5400, help="seconds per native solve")
        p.add_argument("--dry-run", action="store_true", help="write configs and manifests only")
        p.add_argument("--allow-dirty", action="store_true")
        if name == "af-contact":
            p.add_argument("--closed", required=True, help="af-smooth work root with closed result_*.json")
        if name == "mixture":
            p.add_argument("--base", required=True, help="density-off production root: <base>/<case>/runs/<case>/0*/config.txt")
    p = sub.add_parser("compare")
    p.add_argument("--work", required=True, help="parent of af_smooth/, af_contact/, mixture/")
    args = parser.parse_args()
    {"af-smooth": af_smooth, "af-contact": af_contact, "mixture": mixture, "compare": compare}[args.command](args)


if __name__ == "__main__":
    main()
