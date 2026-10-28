"""Prepare the eight thermal sensitivities; this script never launches CFD."""

import argparse
import copy
import csv
import json
import math
from pathlib import Path
import tarfile
from types import SimpleNamespace

import run_steady as driver

ROOT = Path(__file__).resolve().parents[2]
SUPPLY = driver.THERMAL_BOUNDARY_KEYS[:3]
BATHS = driver.THERMAL_BOUNDARY_KEYS[3:]
SCENARIOS = (("supply_plus_1K", "E", SUPPLY),
             ("cooling_plus_1K", "E", BATHS),
             ("supply_TC2", "ABCDEF", SUPPLY))


def parse_config(text):
    return dict(line.split(" = ", 1) for line in text.splitlines() if " = " in line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-csv", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True, help="resolved_results directory")
    parser.add_argument("--archive-source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="new preparation directory")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "four_cores.rankfile").write_text(
        "".join(f"rank {i}=afml-r6 slot={i}\n" for i in range(4)), encoding="utf-8", newline="\n")
    setup = json.loads(driver.SETUP.read_text())
    with args.experiment_csv.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 6 or {row["case_id"] for row in rows} != set("ABCDEF"):
        raise ValueError("experimental table must contain exactly one row for each A-F")
    exp = {row["case_id"]: row for row in rows}
    campaign = driver.Campaign(SimpleNamespace(materials=driver.SETUP, solver=None,
        work=None, trial=None, rankfile=None))
    geometry = campaign.experiment["geometry"]
    assert (geometry["journal_radius_m"], geometry["bearing_length_m"],
            geometry["radial_clearance_m"]) == (.0215, .04, 42.5e-6)
    with tarfile.open(args.archive_source) as archive:
        source = archive.extractfile("source/cases/02d_fw68d_steady_six/run_steady.py").read()
        native = {m.name.removeprefix("source/"): m for m in archive.getmembers()
                  if m.isfile() and (m.name.startswith("source/src/") or m.name == "source/CMakeLists.txt")}
        current = {p.relative_to(ROOT).as_posix() for p in (ROOT / "src").rglob("*") if p.is_file()}
        if current | {"CMakeLists.txt"} != set(native):
            raise ValueError("native source inventory changed; establish a matched baseline")
        for name, member in native.items():
            old = archive.extractfile(member).read().replace(b"\r\n", b"\n")
            now = (ROOT / name).read_bytes().replace(b"\r\n", b"\n")
            if name == "src/resolved_cross_film.hpp":
                old = old.replace(b"/// coupled to the isothermal shaft and the existing r-theta housing model.",
                                  b"/// coupled to an isothermal shaft and local, r-theta or r-theta-z housing conduction.")
            if old != now:
                raise ValueError(f"native source changed at {name}; establish a matched baseline")
    old_driver = {"__name__": "archived_driver", "__file__": str(driver.HERE / "run_steady.py")}
    exec(compile(source, "archived_run_steady.py", "exec"), old_driver)
    baseline = {row["case"]: row for row in json.loads((args.baseline / "results.json").read_text())}
    baseline_rows = []
    for key in "ABCDEF":
        case = campaign.cases[key]
        row = exp[key]
        for a, b in ((float(row["tc0_c"]), case["tc0_boundary_temperature_c"]),
                     (float(row["load_n"]), case["load_mean_n"]),
                     (float(row["rpm"]), case["rpm_nominal"])):
            assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12), (key, a, b)
        # Compare the complete no-scenario rendering, including numerical references.
        render_args = (case, campaign.experiment, campaign.boundaries, campaign.materials,
                       *setup["force_balance"]["initial_positions"][key], 60., 100.,
                       Path("baseline-output"), driver.TABLE)
        assert driver.render_config(*render_args) == old_driver["render_config"](*render_args), key
        saved = baseline[key]
        m = saved["metrics"]
        assert saved["status"] == "accepted" and m["status"] == "load_and_thermal_closed"
        assert m["returncode"] == 0 and m["force_error_rel"] < driver.FORCE_BALANCE_TOL
        assert m["h_relative_mismatch"] < 1e-3 and m["T_change_K"] < .01
        assert m["energy_balance_residual_rel"] <= driver.ENERGY_BALANCE_TOL
        d = m["diagnostics"]
        assert abs(float(d["liquid_residual_rel"])) <= driver.MASS_BALANCE_TOL
        assert float(d["coupling_R_T"]) < 1e-6 and int(d["coupling_converged"]) == 1
        assert int(d["coupling_restored_best"]) == 0 and int(d["converged"]) == 1
        assert m["cross_film_profile_residual"] < 1e-6
        assert m["experiment_sha256"] == driver.sha256(driver.EXPERIMENT)
        assert m["property_table_sha256"] == driver.sha256(driver.TABLE)
        experimental_torque = float(row["new_exp_torque_nm"])
        assert experimental_torque > 0 and math.isclose(saved["torque_Nm"], float(row["cfd_torque_nm"]), rel_tol=1e-12)
        baseline_rows.append({"case": key, "torque_Nm": saved["torque_Nm"],
            "new_exp_torque_nm": experimental_torque,
            "torque_error_percent": 100 * (saved["torque_Nm"] / experimental_torque - 1),
            **{f"TC{i}_error_K": saved["TC_C"][i-1] - float(row[f"tc{i}_exp_c"]) for i in (1, 2, 3)}})
    manifest = {"equilibria": 8, "baseline_repeats": 0,
        "baseline_render_matches_archive": True,
        "experiment_csv_sha256": driver.sha256(args.experiment_csv),
        "archive_source_sha256": driver.sha256(args.archive_source),
        "baseline_results_sha256": driver.sha256(args.baseline / "results.json"),
        "driver_sha256": driver.sha256(driver.__file__),
        "source_sha256": driver.source_sha256(), "scenarios": []}
    # Preserve the corrected comparison input as evidence, without editing experiment.json.
    (args.out / "experimental_comparison.csv").write_bytes(args.experiment_csv.read_bytes())
    driver.save_json(args.out / "baseline_comparison.json", baseline_rows)
    manifest["baseline_comparison_sha256"] = driver.sha256(args.out / "baseline_comparison.json")
    for name, keys, boundaries in SCENARIOS:
        scenario = copy.deepcopy(setup)
        scenario["case_ids"] = list(keys)
        scenario["sensitivity_scenario"] = {"name": name,
            "interpretation": "Hypothetical sensitivity amplitude, not measured uncertainty or inlet temperature",
            "TC2_independent_validation": name != "supply_TC2"}
        scenario["thermal_boundaries_by_case"] = {}
        changes = {}
        for key in keys:
            case = campaign.cases[key]
            temperature = (float(exp[key]["tc2_exp_c"]) + 273.15 if name == "supply_TC2"
                           else case["tc0_boundary_temperature_k"] + 1.)
            overrides = dict.fromkeys(boundaries, temperature)
            scenario["thermal_boundaries_by_case"][key] = overrides
            altered = {**case, "thermal_boundaries_K": overrides, "sensitivity_diagnostics": True}
            render_args = (campaign.experiment, campaign.boundaries, campaign.materials,
                           *setup["force_balance"]["initial_positions"][key], 60., 100.,
                           Path("preview-output"), driver.TABLE)
            before = parse_config(driver.render_config(case, *render_args))
            text = driver.render_config(altered, *render_args)
            after = parse_config(text)
            delta = {k: {"baseline": before[k], "scenario": after[k]} for k in before if before[k] != after[k]}
            expected = {"inlet_groove" if k == "groove_supply_temperature" else k for k in boundaries}
            assert set(delta) == expected | {"output_fields"}, (name, key, delta)
            assert driver.warm_start_identity(text) != driver.warm_start_identity(driver.render_config(case, *render_args))
            changes[key] = delta
            (args.out / f"{name}_{key}_preview.txt").write_text(text, encoding="utf-8")
        path = args.out / f"{name}.json"
        driver.save_json(path, scenario)
        driver.save_json(args.out / f"{name}_config_diff.json", changes)
        manifest["scenarios"].append({"name": name, "cases": list(keys),
            "setup_sha256": driver.sha256(path), "equilibria": len(keys)})
    # Include every native coefficient-closure call, not only the final solve.
    runtimes = {key: [] for key in "ABCDEF"}
    with tarfile.open(args.baseline / "final_results.tar.gz") as archive:
        for member in archive.getmembers():
            if member.name.startswith("history/") and member.name.endswith("/metrics.json"):
                record = json.load(archive.extractfile(member))
                if record.get("case_id") in runtimes:
                    assert record["returncode"] == 0 and record["mpi_ranks"] == 4
                    runtimes[record["case_id"]].append(record["wall_seconds"])
    seconds = {key: sum(values) for key, values in runtimes.items()}
    assert all(seconds.values())
    total = sum(seconds.values()) + 2 * seconds["E"]
    manifest["cost"] = {"baseline_native_calls": {key: len(v) for key, v in runtimes.items()},
        "baseline_seconds_at_4_ranks": seconds,
        "estimated_native_calls": sum(map(len, runtimes.values())) + 2 * len(runtimes["E"]),
        "sequential_wall_hours": total / 3600, "core_hours": total * 4 / 3600,
        "planning_wall_hours_1_to_2x": [total / 3600, 2 * total / 3600],
        "planning_core_hours_1_to_2x": [total * 4 / 3600, total * 8 / 3600],
        "maximum_native_calls": 8 * 30,
        "timeout_note": "7200 seconds per native solve plus 120 seconds stop grace; not a campaign deadline",
        "price": "No workstation billing rate supplied; monetary cost is unknown"}
    driver.save_json(args.out / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
