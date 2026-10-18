"""Prepare or run the fixed-chart-30% LGE surrogate using the established driver."""

import argparse
import copy
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cases/02d_fw68d_steady_six"))
import run_steady as driver
from thermal_expansion import parse_config, solid_observables


def prepare(preflight, destination):
    """Freeze user-confirmed reservoirs and explicitly conditional surrogate inputs."""
    destination.mkdir(parents=True, exist_ok=False)
    rows = list(csv.DictReader((preflight / "case_matrix.csv").open(encoding="utf-8")))
    assert [r["selected_id"] for r in rows] == [f"M{i:02d}" for i in range(1, 18)]
    points = [{"temperature_c": float(r["temperature_c"]), "nu_mm2_s": float(r["nu_mm2_s"])}
              for r in csv.DictReader((preflight / "viscosity_nodes.csv").open())]
    geometry = dict(journal_radius_m=.0215, radial_clearance_m=42.5e-6,
                    bearing_length_m=.040, bush_outer_radius_m=.0425, sensor_radius_m=.026)
    experiment = {"geometry": geometry, "viscosity": {"points": points}, "cases": []}
    reference = json.loads(driver.SETUP.read_text())
    decisions = {
        "status": "conditional_surrogate_authorized_for_execution",
        "user_confirmation": "Pressure sensor inside chamber; groove and both ends share chamber pressure. Exposed surfaces in liquid at roughly TC0. Fluid is undersaturated, with measured roughly30% before run; concentration basis remains unverified.",
        "concentration": "Fixed chart30% viscosity curve, no tracked dissolved/free species or release; never reinterpret as measured mass fraction.",
        "pressure": "Each measured bar-absolute mean times100000 at both ends and groove.",
        "cavity": "Explicit zero-Pa-absolute oil-vapor-limit surrogate with refrigerant release suppressed. This is not measured mixture bubble pressure or R410A partial pressure.",
        "thermal_analogs": {"cp_J_kg_K": 2010., "conductivity_W_m_K": .13,
                           "beta_1_K": .7284/950., "steel_conductivity_W_m_K": 15.},
        "thermal_analog_basis": "Reuse BVC32-PVE cp/k and304 stainless analogs; pure-oil density-slope buoyancy analog separately declared under constant bulk rho950. None is a measured mixture thermal property.",
        "temperature": "TC0 supplies feed and liquid bath; film, shaft and housing temperatures solved. TC2 only referenceS and sensor comparison.",
        "rotation": "All selected motorRPM signs negative, same as pure-oil baseline; inherit its positive-theta-following-rotation frame and sensor coordinates. Motor sign alone does not establish an independent laboratory CW/CCW mapping.",
        "geometry": "Nominal only; no thermal expansion.",
        "approval": "User instructed full remote execution without stopping for approval after pilots; technical clarification supplied during preflight.",
        "limitations": "No solubility/release, composition evolution, thermal expansion, pressure-viscosity law or calibrated pressurized parasitic drag; condition stability does not prove steady friction or composition."}
    starts = {}
    for row in rows:
        key, rpm, load = row["selected_id"], float(row["rpm"]), float(row["force_n_mean"])
        temperature = float(row["tc0_c_mean"])
        assert float(row["speed_rpm_mean"]) < 0 and float(row["signed_requested_rpm"]) < 0
        nu = driver.viscosity(temperature, points)
        case = dict(case_id=key, rpm_nominal=rpm, load_mean_n=load,
                    tc0_boundary_temperature_c=temperature, tc0_boundary_temperature_k=temperature+273.15,
                    angular_speed_magnitude_rad_s=rpm*math.pi/30, angular_speed_rad_s=rpm*math.pi/30,
                    journal_diameter_m=.043, bearing_length_m=.04, radial_clearance_m=42.5e-6,
                    rho_kg_m3=950., nu_mm2_s=nu*1e6, mu_pa_s=950*nu,
                    sommerfeld_number=float(row["S"]), apparent_bearing_torque_nm=float(row["bearing_residual_nm"]),
                    friction_coefficient=float(row["f"]), measured=row)
        experiment["cases"].append(case)
        # Seeds affect convergence only; no measured torque is used to choose them.
        eps = float(np.interp(math.log(float(row["S"])), np.log([.004, .03, .1, 1.]), [.94, .85, .60, .15]))
        starts[key] = [eps, 300. if eps > .8 else 320.]
    for case in experiment["cases"]:
        key = case["case_id"]
        setup = {key: copy.deepcopy(reference[key]) for key in
                 ("force_balance", "material_settings", "acceptance", "comparison_sensor_locations")}
        setup["specified_boundaries"] = {
            "groove": copy.deepcopy(reference["specified_boundaries"]["groove"]),
            "housing": {"outer_radius_m": .0425}, "provenance": decisions["user_confirmation"]}
        setup["specified_boundaries"]["groove"]["pressure"] = "Measured chamber absolute pressure for this case"
        setup.update(status="conditional_fixed_chart30_surrogate", case_ids=list(starts), additional_cases=[])
        # The initial M02 pilot failed mass closure at the inherited pressure tolerance.
        setup["solver_tolerances"] = {"fbns_tol": 1e-10}
        setup["force_balance"]["initial_positions"] = starts
        setup["material_settings"].update(oil_density_kg_m3=950., oil_cp_J_kg_K=2010.,
            oil_conductivity_W_m_K=.13, oil_thermal_expansion_1_K=.7284/950.,
            bush_conductivity_W_m_K=15., cavitation_pressure_Pa=0.,
            reservoir_pressure_Pa=float(case["measured"]["pressure_bar_abs_mean"])*100000)
        setup["material_settings"]["pressure_reference"] = decisions["pressure"] + " " + decisions["cavity"]
        setup["material_settings"]["oil_thermal_analogue"] = decisions["thermal_analog_basis"]
        setup["scenario_decisions"] = decisions
        driver.save_json(destination / f"{key}_setup.json", setup)
    driver.save_json(destination / "experiment.json", experiment)
    driver.save_json(destination / "decisions.json", decisions)
    table = preflight / "solution_viscosity_surrogate.csv"
    (destination / table.name).write_bytes(table.read_bytes())
    driver.save_json(destination / "manifest.json", {p.name: driver.sha256(p) for p in sorted(destination.iterdir())})
    print(json.dumps({"cases": len(rows), "output": str(destination), "status": "prepared"}))


def options(inputs, destination, solver, rankfile):
    return SimpleNamespace(solver=str(solver) if solver else None, work=str(destination / "runs"),
        experiment=str(inputs / "experiment.json"), table=str(inputs / "solution_viscosity_surrogate.csv"),
        materials=str(destination / "setup.json"), ranks=4, rankfile=str(rankfile) if rankfile else None, trial=None,
        timeout=7200, max_thermal_iterations=30, coupling_max_iters=4000,
        coupling_relaxation=1., initial_state_root=None)


def run_case(args):
    destination = args.work / args.case
    destination.mkdir(parents=True, exist_ok=False)
    setup = json.loads((args.inputs / f"{args.case}_setup.json").read_text())
    setup["material_settings"]["cross_film_layers"] = args.layers
    if args.tight:
        setup["solver_tolerances"] = dict(steady_force_tolerance=2e-5, coupling_tolerance=1e-8, fbns_tol=1e-12)
        setup["cooling_relative_tolerance"] = 1e-5
        setup["cooling_temperature_tolerance_K"] = 1e-4
    driver.save_json(destination / "setup.json", setup)
    campaign = driver.Campaign(options(args.inputs, destination, args.solver, args.rankfile))
    campaign.mpi_version = subprocess.check_output([campaign.command()[0], "--version"], text=True).splitlines()[:2]
    case = campaign.cases[args.case]
    seed = setup["force_balance"]["initial_positions"][args.case]
    start = time.monotonic()
    try:
        result = campaign.balance(case, *seed)
        run = Path(result["run_dir"])
        checked = campaign.recheck(run, result, case)
        closed = campaign.find_completed(case, result["seed_eps"], result["seed_psi_deg"], closed=True)
        if checked is None or closed is None or Path(closed["run_dir"]) != run:
            raise ValueError("Result failed independent native and cooling recheck")
        tolerance = setup.get("solver_tolerances", {})
        thermal = tolerance.get("coupling_tolerance", 1e-6)
        if (result["force_error_rel"] >= tolerance.get("steady_force_tolerance", .002) or
                float(result["diagnostics"]["coupling_R_T"]) >= thermal or
                result["cross_film_profile_residual"] >= thermal or
                result["energy_balance_residual_rel"] >= thermal):
            raise ValueError("Result failed requested tighter numerical tolerances")
        config = parse_config((run / "config.txt").read_text())
        raw = driver.read_fields(run)
        assert float(config["omega"])*result["torque_Nm"] > 0
        observables = solid_observables(config, raw, .026)
        np.savez_compressed(destination / "fields.npz", **raw)
        accepted = {"case": case["measured"], "result": result, "config": config,
                    "observables": observables, "elapsed_seconds": time.monotonic()-start,
                    "worker_sha256": driver.sha256(__file__), "fields_sha256": driver.sha256(destination / "fields.npz"),
                    "tight": args.tight, "cross_film_layers": args.layers,
                    "status": "accepted_conditional_surrogate"}
        driver.save_json(destination / "accepted.json", accepted)
        print(json.dumps({"case": args.case, "torque_Nm": result["torque_Nm"], "TC_C": observables["TC_C"]}), flush=True)
    except Exception as error:
        driver.save_json(destination / "failure.json", {"case": args.case, "error": str(error),
            "elapsed_seconds": time.monotonic()-start, "status": "failed_preserved"})
        raise


def check_inputs(inputs, output):
    """Check all rendered physical controls and independent S/property identities."""
    verified = {}
    for key in (f"M{i:02d}" for i in range(1, 18)):
        args = options(inputs, output.parent / key, None, None)
        args.materials = str(inputs / f"{key}_setup.json")
        campaign = driver.Campaign(args)
        case = campaign.cases[key]
        original = case["measured"]
        seed = campaign.setup["force_balance"]["initial_positions"][key]
        config = campaign.config_for(case, *seed, 60., 100., output.parent / "unused")
        values = parse_config(config)
        pressure = float(original["pressure_bar_abs_mean"])*1e5
        assert float(values["bc_z_south_val"]) == float(values["bc_z_north_val"]) == pressure
        assert float(values["inlet_groove"].split()[2]) == pressure
        assert values["density_model"] == "CONSTANT" and float(values["rho"]) == 950.
        assert values["solution_table_clamp"] == "false"
        assert float(values["dissolved_gas_initial"]) == float(values["dissolved_gas_henry_coeff"]) == 0.
        assert values["fluid_property_model"] == "SINGLE_PHASE" and float(values["p_cav"]) == 0.
        assert float(values["R"]) == .0215 and float(values["c"]) == 42.5e-6 and float(values["L"]) == .04
        assert float(values["external_load_y"]) == float(original["force_n_mean"])
        assert float(values["density_temperature_slope"]) == -.7284
        assert float(values["liquid_specific_heat"]) == 2010. and float(values["thermal_conductivity"]) == .13
        nu = driver.viscosity(float(original["tc2_c_mean"]), campaign.points)
        s_value = 950*nu*float(original["rpm"])/60*.04*.043/float(original["force_n_mean"])*(.0215/42.5e-6)**2
        assert math.isclose(s_value, float(original["S"]), rel_tol=1e-13)
        assert math.isclose(float(values["mu"]), case["mu_pa_s"], rel_tol=1e-14)
        verified[key] = {"pressure_Pa_abs": pressure, "S_saved": float(original["S"]),
                         "feed_mu_Pa_s": float(values["mu"]), "reference_mu_Pa_s": 950*nu,
                         "omega_rad_s": float(values["omega"]),
                         "seed": seed, "inputs": campaign.input_hashes}
    with output.open("x", encoding="utf-8") as stream:
        json.dump({"status": "17_native_configurations_verified_unrun", "cases": verified}, stream, indent=2)
    print(json.dumps({"verified": len(verified), "output": str(output)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", type=Path, help="verified preflight directory")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--case", choices=[f"M{i:02d}" for i in range(1, 18)])
    parser.add_argument("--work", type=Path)
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--rankfile", type=Path)
    parser.add_argument("--tight", action="store_true")
    parser.add_argument("--layers", type=int, choices=(16, 32), default=16)
    parser.add_argument("--check-inputs", type=Path, help="write an exclusive validation report")
    args = parser.parse_args()
    if args.prepare:
        prepare(args.prepare, args.inputs)
    elif args.check_inputs:
        check_inputs(args.inputs, args.check_inputs)
    elif args.case and args.work and args.solver and args.rankfile:
        run_case(args)
    else:
        parser.error("run requires case, work, solver and rankfile")


if __name__ == "__main__":
    main()
