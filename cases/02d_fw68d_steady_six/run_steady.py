"""Run the six load-controlled steady FW68D cases, preserving every solver call.

Run ``python run_steady.py --self-check`` before a remote campaign. A fixed
geometry timing call is ``--trial A 0.85 120 --thermal-iterations 1``; it is
explicitly not a thermally closed bearing result.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

import numpy as np


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE / "experiment.json"
TABLE = HERE / "fw68d_experiment.csv"
SETUP = HERE / "steady_setup.json"
NUMBER = r"([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_sha256():
    root = HERE.parents[1]
    digest = hashlib.sha256()
    paths = [root / "CMakeLists.txt", *sorted((root / "src").rglob("*.cpp")),
             *sorted((root / "src").rglob("*.hpp"))]
    for path in paths:
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(bytes.fromhex(sha256(path)))
    return digest.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def density(temperature_c):
    return 954.9 - 0.7284 * temperature_c


def viscosity(temperature_c, points):
    """Return nu in m2/s using the experiment's ln(nu) interpolation."""
    for left, right in zip(points, points[1:]):
        a, b = left["temperature_c"], right["temperature_c"]
        if a <= temperature_c <= b:
            fraction = (temperature_c - a) / (b - a)
            return math.exp((1 - fraction) * math.log(left["nu_mm2_s"])
                            + fraction * math.log(right["nu_mm2_s"])) * 1e-6
    raise ValueError(f"FW68D viscosity table does not cover {temperature_c:g} C")


def oil_properties(temperature_c, points, materials):
    nu = viscosity(temperature_c, points)
    rho = density(temperature_c)
    cp = float(materials["oil_cp_J_kg_K"])
    k = float(materials["oil_conductivity_W_m_K"])
    if min(nu, rho, cp, k) <= 0:
        raise ValueError("nonpositive oil property")
    pr = nu * rho * cp / k
    return nu, rho, pr, k


def housing_coefficient(outer_temperature_k, ambient_k, points, materials):
    film_c = (outer_temperature_k + ambient_k) / 2 - 273.15
    nu, rho, pr, k = oil_properties(film_c, points, materials)
    diameter = 0.085
    alpha = k / (rho * float(materials["oil_cp_J_kg_K"]))
    rayleigh = (9.80665 * 0.7284 / rho * abs(outer_temperature_k - ambient_k)
                * diameter**3 / (nu * alpha))
    nusselt = (0.60 + 0.387 * rayleigh**(1 / 6)
                / (1 + (0.559 / pr)**(9 / 16))**(8 / 27))**2
    return nusselt * k / diameter, {"Ra_D": rayleigh, "Pr": pr, "Nu_D": nusselt}


def shaft_coefficient(shaft_temperature_k, ambient_k, omega, points, materials):
    if materials.get("shaft_end_cooling_model") != "rotating_disk":
        raise ValueError("shaft_end_cooling_model must be rotating_disk")
    coefficient = materials.get("shaft_end_nusselt_coefficient")
    if coefficient is None:
        raise ValueError("shaft-end Nusselt coefficient is unresolved")
    film_c = (shaft_temperature_k + ambient_k) / 2 - 273.15
    nu, _, pr, k = oil_properties(film_c, points, materials)
    radius = float(materials["shaft_exposed_end_radius_m"])
    reynolds = abs(omega) * radius**2 / nu
    nusselt = (float(coefficient)
                * reynolds**float(materials["shaft_end_re_exponent"])
                * pr**float(materials["shaft_end_pr_exponent"]))
    return nusselt * k / radius, {"Re_R": reynolds, "Pr": pr, "Nu_R": nusselt,
                                   "laminar_Re_under_1e5": reynolds < 1e5}


def outer_temperature(bore_mean_k, ambient_k, h_bush, radius, outer_radius, k_bush):
    resistance = 1 + h_bush * outer_radius / k_bush * math.log(outer_radius / radius)
    return ambient_k + (bore_mean_k - ambient_k) / resistance


def render_config(case, experiment, boundaries, materials, eps, psi, h_bush, h_shaft, output_dir, table):
    geometry = experiment["geometry"]
    radius = float(geometry["journal_radius_m"])
    clearance = float(geometry["radial_clearance_m"])
    temperature = float(case["tc0_boundary_temperature_k"])
    load = float(case["load_mean_n"])
    omega = float(case["angular_speed_magnitude_rad_s"])
    rho = density(temperature - 273.15)
    nu = viscosity(temperature - 273.15, experiment["viscosity"]["points"])
    values = {
        "R": radius, "c": clearance, "e": eps * clearance,
        "L": geometry["bearing_length_m"], "attitude_angle_deg": psi,
        "load_angle_deg": -90, "omega": omega, "n_theta_global": 720,
        "n_z_global": 192, "solution_mode": "STEADY_STATE", "dt": 2e-5,
        "end_t": 2e-5, "write_interval": 2e-5, "diagnostics_interval": 1,
        "motion_model": "STATIC", "motion_time_method": "EULER_IMPLICIT",
        "theta_convection_scheme": "UPWIND", "thermal_convection_scheme": "UPWIND",
        "gas_convection_scheme": "UPWIND", "external_load_x": 0,
        "external_load_y": -load, "external_load_z": 0,
        "cavitation_model": "FBNS", "cavitation_threshold": "SCALAR_PCAV",
        "cavitated_film_model": "STRIATED", "gas_pressure_coupling": "NONE",
        "p_cav": materials["cavitation_pressure_Pa"],
        "bc_z_south_type": "DIRICHLET", "bc_z_south_val": materials["reservoir_pressure_Pa"],
        "bc_z_north_type": "DIRICHLET", "bc_z_north_val": materials["reservoir_pressure_Pa"],
        "bc_z_south_thermal": "CONSTANT", "bc_z_north_thermal": "CONSTANT",
        "bc_z_south_temperature": temperature, "bc_z_north_temperature": temperature,
        "temperature_model": "ENERGY_EQUATION", "temperature_initial": temperature,
        "temperature_reference": temperature, "property_reference_temperature": temperature,
        "rho_cp": rho * float(materials["oil_cp_J_kg_K"]),
        "thermal_conductivity": materials["oil_conductivity_W_m_K"],
        "cross_film_thermal_model": "POLYNOMIAL", "cross_film_profile_order": 4,
        "cross_film_quadrature_points": 6, "bush_thermal_model": "CONDUCTION_2D",
        "bush_thermal_conductivity": materials["bush_conductivity_W_m_K"],
        "bush_outer_radius": boundaries["housing"]["outer_radius_m"],
        "bush_radial_cells": 16, "bush_outer_boundary": "CONVECTION",
        "bush_outer_heat_transfer": h_bush, "shaft_ambient_heat_transfer": h_shaft,
        "shaft_ambient_area_ratio": materials["shaft_end_area_ratio"],
        "thermal_ambient_temperature": temperature, "shaft_ambient_temperature": temperature,
        "fluid_property_model": "SINGLE_PHASE", "dissolved_gas_species": "AIR",
        "dissolved_gas_initial": 0, "solubility_model": "HENRY",
        "dissolved_gas_henry_coeff": 0, "density_model": "LINEAR_TEMPERATURE",
        "density_temperature_intercept": 954.9, "density_temperature_slope": -0.7284,
        "rho": rho, "liquid_viscosity_model": "KINEMATIC_TABLE",
        "solution_table_file": str(table.resolve()), "solution_table_clamp": "false",
        "kinematic_viscosity_log_interpolation": "true", "viscosity_pressure_coeff": 0,
        "mu": rho * nu, "groove_mixing_model": "PINNED",
        "inlet_groove": f"90 {boundaries['groove']['width_deg']} {materials['reservoir_pressure_Pa']} {temperature}",
        "coupling_max_iters": 200, "coupling_tolerance": 1e-6,
        "coupling_relaxation_mode": "ADAPTIVE", "coupling_relaxation": 1,
        "coupling_relaxation_min": 0.1, "coupling_relaxation_backoff": 0.5,
        "coupling_relaxation_growth": 1.05, "coupling_keep_best_iterate": "false",
        "log_outer_iters": "true", "log_coupling_history": "true",
        "output_dir": str(output_dir.resolve()),
        "output_write_3d": "false", "output_write_flat": "true",
        "output_verbosity": "DETAILED",
        "output_fields": "pressure,film_content,h,temperature,mu,temperature_journal,temperature_bearing,heat_flux_journal,heat_flux_bearing,fluid_force_x,fluid_force_y,fluid_force_z,pressure_force_x,pressure_force_y,friction_torque",
    }
    return "# Generated by run_steady.py; steady fixed-position load-balance evaluation.\n" + "\n".join(
        f"{key} = {value}" for key, value in values.items()) + "\n"


def read_fields(run_dir):
    pvd = run_dir / "results/flat/results.pvd"
    frames = list(ET.parse(pvd).getroot().iter("DataSet"))
    if not frames:
        raise ValueError("no flat VTK frame")
    frame = max(frames, key=lambda item: float(item.get("timestep")))
    pvts = pvd.parent / frame.get("file")
    arrays = {}
    for piece in ET.parse(pvts).getroot().iter("Piece"):
        source = pvts.parent / piece.get("Source")
        cell = next(ET.parse(source).getroot().iter("CellData"))
        for array in cell.findall("DataArray"):
            arrays.setdefault(array.get("Name"), []).append(
                np.fromstring(array.text or "", sep=" "))
    return {key: np.concatenate(values) for key, values in arrays.items()}


def last_diagnostic(run_dir):
    with (run_dir / "results/diagnostics.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("diagnostics.csv has no rows")
    return rows[-1]


def check_solution(run_dir, log, returncode, load, cavity_pressure):
    if returncode:
        raise ValueError(f"solver exit code {returncode}")
    diagnostic = last_diagnostic(run_dir)
    for key in ("liquid_residual", "cavitated_fraction", "energy_residual_W"):
        if not math.isfinite(float(diagnostic[key])):
            raise ValueError(f"nonfinite or missing {key} diagnostic")
    for key in ("coupling_R_p", "coupling_R_T"):
        if not math.isfinite(float(diagnostic[key])) or float(diagnostic[key]) >= 1e-6:
            raise ValueError(f"{key} failed: {diagnostic[key]}")
    if int(diagnostic["coupling_converged"]) != 1 or int(diagnostic["coupling_restored_best"]) != 0:
        raise ValueError("pressure/thermal coupling did not converge cleanly")
    if int(diagnostic["converged"]) != 1:
        raise ValueError("FBNS did not converge")
    fbns = re.findall(r"fbns iter=\d+ \|G\|=" + NUMBER + r" \|FB\|=" + NUMBER, log)
    if not fbns or max(map(float, fbns[-1])) > 1e-7:
        raise ValueError("FBNS complementarity residual missing or above 1e-7")
    force = re.findall(r"\| F=\(" + NUMBER + r",\s*" + NUMBER + r",\s*" + NUMBER + r"\) N\s*\| M=" + NUMBER, log)
    if not force:
        raise ValueError("final total-fluid-force/torque record missing")
    fx, fy, fz, torque = map(float, force[-1])
    fields = read_fields(run_dir)
    for field in ("p", "film_content", "h", "T", "T_bearing", "T_journal", "mu",
                  "q_bearing", "q_journal", "friction_torque"):
        if field not in fields or not len(fields[field]) or not np.isfinite(fields[field]).all():
            raise ValueError(f"missing/nonfinite {field} field")
    if float(np.min(fields["h"])) <= 0:
        raise ValueError("nonpositive film gap")
    if float(np.min(fields["p"])) < cavity_pressure - 1:
        raise ValueError("pressure below cavitation floor")
    if float(np.min(fields["film_content"])) < -1e-6 or float(np.max(fields["film_content"])) > 1 + 1e-6:
        raise ValueError("film content outside [0, 1]")
    if "F" not in fields or len(fields["F"]) % 3:
        raise ValueError("missing total fluid force vector F")
    vtk_force = fields["F"].reshape(-1, 3)
    if not np.isfinite(vtk_force).all():
        raise ValueError("nonfinite total fluid force vector F")
    if not np.allclose(vtk_force, vtk_force[0], rtol=1e-9, atol=1e-6):
        raise ValueError("total force differs across VTK cells")
    fx, fy, fz = map(float, vtk_force[0])
    torque = float(fields["friction_torque"][0])
    if not all(math.isfinite(value) for value in (fx, fy, fz, torque)):
        raise ValueError("nonfinite integrated force or torque")
    return {
        "Fx_N": fx, "Fy_N": fy, "Fz_N": fz, "torque_Nm": torque,
        "force_error_rel": math.hypot(fx, fy - load) / load,
        "h_min_m": float(np.min(fields["h"])),
        "p_min_Pa": float(np.min(fields["p"])), "p_max_Pa": float(np.max(fields["p"])),
        "cavitated_area_fraction": float(np.mean(fields["film_content"] < .999)),
        "film_content_min": float(np.min(fields["film_content"])),
        "T_min_K": float(np.min(fields["T"])), "T_max_K": float(np.max(fields["T"])),
        "T_bore_mean_K": float(np.mean(fields["T_bearing"])),
        "T_shaft_mean_K": float(np.mean(fields["T_journal"])),
        "diagnostics": diagnostic,
    }


class Campaign:
    def __init__(self, args):
        self.args = args
        self.experiment = json.loads(EXPERIMENT.read_text(encoding="utf-8"))
        self.materials_path = Path(args.materials).resolve()
        self.setup = json.loads(self.materials_path.read_text(encoding="utf-8"))
        self.materials = self.setup["material_settings"]
        self.boundaries = self.setup["specified_boundaries"]
        self.cases = {item["case_id"]: item for item in self.experiment["cases"]}
        for extra in self.setup.get("additional_cases", []):
            case = dict(self.cases[extra["base_case"]])
            case["case_id"] = extra["case_id"]
            case["tc0_boundary_temperature_c"] = extra["tc0_boundary_temperature_C"]
            case["tc0_boundary_temperature_k"] = extra["tc0_boundary_temperature_K"]
            case["nu_mm2_s"] = viscosity(case["tc0_boundary_temperature_c"],
                                            self.experiment["viscosity"]["points"]) * 1e6
            case["rho_kg_m3"] = density(case["tc0_boundary_temperature_c"])
            case["mu_pa_s"] = case["rho_kg_m3"] * case["nu_mm2_s"] * 1e-6
            self.cases[case["case_id"]] = case
        self.points = self.experiment["viscosity"]["points"]
        self.solver = Path(args.solver).resolve() if args.solver else None
        self.work = Path(args.work).resolve() if args.work else None
        self.sequence = 0
        self.h_warm = {}
        self.source_hash = source_sha256()
        self.mpi_version = None

    def command(self):
        mpi = shutil.which("mpiexec") or shutil.which("mpirun")
        if not mpi:
            raise FileNotFoundError("mpiexec/mpirun not found")
        command = [mpi]
        if self.args.cpu_list:
            if os.name == "nt":
                raise ValueError("--cpu-list requires OpenMPI on Linux")
            command += ["--bind-to", "core", "--cpu-list", self.args.cpu_list]
        command += ["-n", str(self.args.ranks), str(self.solver), "--config", "config.txt"]
        return command

    def single(self, case, eps, psi, h_bush, h_shaft, label):
        self.sequence += 1
        run_dir = self.work / case["case_id"] / f"{self.sequence:04d}_{label}_{uuid.uuid4().hex[:8]}"
        run_dir.mkdir(parents=True, exist_ok=False)
        config = render_config(case, self.experiment, self.boundaries, self.materials, eps, psi, h_bush, h_shaft,
                               run_dir / "results", TABLE)
        config_path = run_dir / "config.txt"
        config_path.write_text(config, encoding="utf-8")
        command = self.command()
        record = {
            "case_id": case["case_id"], "eps": eps, "psi_deg": psi,
            "h_bush_W_m2_K": h_bush, "h_shaft_W_m2_K": h_shaft,
            "run_dir": str(run_dir), "command": command,
            "mpi_ranks": self.args.ranks, "python_version": sys.version.split()[0],
            "mpi_version": self.mpi_version, "source_sha256": self.source_hash,
            "config_sha256": sha256(config_path), "executable_sha256": sha256(self.solver),
            "experiment_sha256": sha256(EXPERIMENT), "material_setup_sha256": sha256(self.materials_path),
            "property_table_sha256": sha256(TABLE), "status": "started",
        }
        save_json(run_dir / "metrics.json", record)
        started = time.monotonic()
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        try:
            with (run_dir / "run.log").open("w", encoding="utf-8") as stream:
                process = subprocess.run(command, cwd=run_dir, stdout=stream, stderr=subprocess.STDOUT,
                                         env=environment, timeout=self.args.timeout)
            record["returncode"] = process.returncode
            record["wall_seconds"] = time.monotonic() - started
            log = (run_dir / "run.log").read_text(errors="replace")
            record.update(check_solution(run_dir, log, process.returncode,
                                         float(case["load_mean_n"]), float(self.materials["cavitation_pressure_Pa"])))
            record["status"] = "valid_inner_solve"
        except Exception as error:
            record["wall_seconds"] = time.monotonic() - started
            record["status"] = "failed"
            record["error"] = str(error)
            save_json(run_dir / "metrics.json", record)
            raise RuntimeError(f"solver evaluation failed; retained {run_dir}: {error}") from error
        save_json(run_dir / "metrics.json", record)
        return record

    def evaluate(self, case, eps, psi, label, thermal_iterations=None):
        if not .0 < eps < .995:
            raise ValueError("eccentricity ratio must lie in (0, 0.995)")
        ambient = float(case["tc0_boundary_temperature_k"])
        omega = float(case["angular_speed_magnitude_rad_s"])
        if case["case_id"] in self.h_warm:
            h_bush, h_shaft = self.h_warm[case["case_id"]]
        else:
            h_bush, _ = housing_coefficient(ambient + 5, ambient, self.points, self.materials)
            h_shaft, _ = shaft_coefficient(ambient + 5, ambient, omega, self.points, self.materials)
        previous_outer = None
        previous_shaft = None
        count = thermal_iterations if thermal_iterations is not None else self.args.max_thermal_iterations
        for iteration in range(count):
            result = self.single(case, eps, psi, h_bush, h_shaft, f"{label}_thermal{iteration + 1}")
            outer = outer_temperature(result["T_bore_mean_K"], ambient, h_bush,
                                      self.experiment["geometry"]["journal_radius_m"],
                                      self.boundaries["housing"]["outer_radius_m"],
                                      float(self.materials["bush_conductivity_W_m_K"]))
            shaft = result["T_shaft_mean_K"]
            target_bush, housing_meta = housing_coefficient(outer, ambient, self.points, self.materials)
            target_shaft, shaft_meta = shaft_coefficient(shaft, ambient, omega, self.points, self.materials)
            rel = max(abs(target_bush / h_bush - 1), abs(target_shaft / h_shaft - 1))
            delta_temperature = (math.inf if previous_outer is None else
                                 max(abs(outer - previous_outer), abs(shaft - previous_shaft)))
            result.update(T_outer_mean_K=outer, T_shaft_K=shaft,
                          h_bush_target_W_m2_K=target_bush, h_shaft_target_W_m2_K=target_shaft,
                          h_relative_mismatch=rel, T_change_K=delta_temperature if math.isfinite(delta_temperature) else None,
                          housing_correlation=housing_meta, shaft_correlation=shaft_meta)
            if rel < 1e-3 and delta_temperature < .01:
                result["status"] = "thermally_closed"
                self.h_warm[case["case_id"]] = (h_bush, h_shaft)
                save_json(Path(result["run_dir"]) / "metrics.json", result)
                return result
            if thermal_iterations is not None and iteration == count - 1:
                result["status"] = "trial_only_not_thermally_closed"
                save_json(Path(result["run_dir"]) / "metrics.json", result)
                return result
            save_json(Path(result["run_dir"]) / "metrics.json", result)
            previous_outer, previous_shaft = outer, shaft
            h_bush = .5 * (h_bush + target_bush)
            h_shaft = .5 * (h_shaft + target_shaft)
        result["status"] = "thermal_closure_failed"
        save_json(Path(result["run_dir"]) / "metrics.json", result)
        raise RuntimeError(f"thermal closure cap reached for {case['case_id']}; last run {result['run_dir']}")

    def balance(self, case, eps, psi):
        load = float(case["load_mean_n"])
        history = []
        for iteration in range(self.args.max_force_iterations):
            result = self.evaluate(case, eps, psi, f"force{iteration}_base")
            residual = np.array([result["Fx_N"], result["Fy_N"] - load])
            history.append({"eps": eps, "psi_deg": psi, "force_error_rel": result["force_error_rel"],
                            "run_dir": result["run_dir"]})
            save_json(self.work / case["case_id"] / "force_history.json", history)
            if result["force_error_rel"] < .002:
                result["status"] = "load_and_thermal_closed"
                result["force_history"] = history
                save_json(self.work / case["case_id"] / f"result_{uuid.uuid4().hex[:8]}.json", result)
                return result
            de = min(.002, (.995 - eps) * .2)
            dp = .25
            e_result = self.evaluate(case, eps + de, psi, f"force{iteration}_de")
            p_result = self.evaluate(case, eps, psi + dp, f"force{iteration}_dp")
            force = np.array([result["Fx_N"], result["Fy_N"]])
            jacobian = np.column_stack(((np.array([e_result["Fx_N"], e_result["Fy_N"]]) - force) / de,
                                        (np.array([p_result["Fx_N"], p_result["Fy_N"]]) - force) / dp))
            step = np.linalg.solve(jacobian, -residual)
            scale = min(1., .08 / max(abs(step[0]), 1e-12), 12 / max(abs(step[1]), 1e-12))
            eps = float(np.clip(eps + scale * step[0], .02, .994))
            psi = float(psi + scale * step[1])
            save_json(self.work / case["case_id"] / "force_history.json", history)
        raise RuntimeError(f"force-balance cap reached for {case['case_id']}; history retained")


def self_check(campaign):
    expected = set("ABCDEF")
    if not set(campaign.cases).issuperset(expected):
        raise ValueError("expected exact A-F experimental cases")
    for case in campaign.cases.values():
        temperature = case["tc0_boundary_temperature_c"]
        nu = viscosity(temperature, campaign.points)
        mu = density(temperature) * nu
        if not math.isclose(mu, case["mu_pa_s"], rel_tol=2e-12):
            raise ValueError(f"case {case['case_id']} viscosity mismatch")
        if not math.isclose(case["angular_speed_magnitude_rad_s"],
                            case["rpm_nominal"] * math.pi / 30, rel_tol=2e-12):
            raise ValueError(f"case {case['case_id']} angular speed mismatch")
        if case["case_id"] in expected:
            pressure = case["load_mean_n"] / (case["journal_diameter_m"] * case["bearing_length_m"])
            sommerfeld = mu * (case["rpm_nominal"] / 60) / pressure * (
                case["journal_diameter_m"] / (2 * case["radial_clearance_m"]))**2
            if not math.isclose(sommerfeld, case["sommerfeld_number"], rel_tol=2e-12):
                raise ValueError(f"case {case['case_id']} Sommerfeld mismatch")
        amb = case["tc0_boundary_temperature_k"]
        config = render_config(case, campaign.experiment, campaign.boundaries, campaign.materials, .3, 140,
                               60, 100, Path("/tmp/steady-check"), TABLE)
        if "external_load_y = -" not in config or "solution_mode = STEADY_STATE" not in config:
            raise ValueError("steady or load config mismatch")
        outer = outer_temperature(amb + 10, amb, 60, .0215, .0425, 15)
        hb, _ = housing_coefficient(outer, amb, campaign.points, campaign.materials)
        hs, meta = shaft_coefficient(amb + 10, amb, case["angular_speed_magnitude_rad_s"],
                                     campaign.points, campaign.materials)
        if min(hb, hs) <= 0 or not meta["laminar_Re_under_1e5"]:
            raise ValueError("invalid external coefficient or laminar disk assumption")
    print("Self-check passed: six case mu/S/omega inputs, config force sign, thermal correlations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", help="compiled pancake executable")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--work", help="new or existing campaign output root")
    parser.add_argument("--cases", default="A,B,C,D,E,F", help="comma-separated A-F or E_low,E_high")
    parser.add_argument("--materials", default=str(SETUP), help="JSON with material_settings object")
    parser.add_argument("--cpu-list", help="OpenMPI physical-core list; Linux only")
    parser.add_argument("--timeout", type=int, default=7200, help="seconds per native solve")
    parser.add_argument("--max-thermal-iterations", type=int, default=30)
    parser.add_argument("--max-force-iterations", type=int, default=18)
    parser.add_argument("--trial", nargs=3, metavar=("CASE", "EPS", "PSI"),
                        help="one fixed-geometry benchmark evaluation")
    parser.add_argument("--thermal-iterations", type=int, default=1,
                        help="trial-only number of evaluations; result is unclosed")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    campaign = Campaign(args)
    if args.self_check:
        self_check(campaign)
        return
    if not args.solver or not args.work or args.ranks < 1:
        parser.error("physical runs require --solver, --work, and positive --ranks")
    if not campaign.solver.is_file():
        parser.error("solver executable does not exist")
    try:
        campaign.mpi_version = subprocess.run([campaign.command()[0], "--version"],
                                              capture_output=True, text=True, timeout=10).stdout.splitlines()[:2]
    except (OSError, subprocess.TimeoutExpired):
        campaign.mpi_version = ["unavailable"]
    campaign.work.mkdir(parents=True, exist_ok=True)
    starts = {"A": (.85, 120.), "B": (.6, 140.), "C": (.6, 140.),
              "D": (.35, 155.), "E": (.15, 165.), "F": (.15, 165.),
              "E_low": (.15, 165.), "E_high": (.15, 165.)}
    if args.trial:
        key, eps, psi = args.trial
        case = campaign.cases[key]
        result = campaign.evaluate(case, float(eps), float(psi), "trial", args.thermal_iterations)
        print(json.dumps({"status": result["status"], "run_dir": result["run_dir"],
                          "wall_seconds": result["wall_seconds"]}, indent=2))
        return
    for key in args.cases.split(","):
        key = key.strip()
        if key not in campaign.cases:
            raise ValueError(f"unknown case {key}; use A-F")
        result = campaign.balance(campaign.cases[key], *starts[key])
        print(f"{key}: force error {result['force_error_rel']:.5g}; result {result['run_dir']}", flush=True)


if __name__ == "__main__":
    main()
