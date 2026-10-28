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
import signal
import subprocess
import sys
import tempfile
import time
import uuid
import xml.etree.ElementTree as ET

import numpy as np


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE / "experiment.json"
TABLE = HERE / "fw68d_experiment.csv"
SETUP = HERE / "steady_setup.json"
NUMBER = r"([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)"
FBNS_TOL = 1e-6
TEMPERATURE_FLOOR_TOL_K = 1e-3
MASS_BALANCE_TOL = 1e-8
ENERGY_BALANCE_TOL = 1e-6
FORCE_BALANCE_TOL = .002


def applied_load(load):
    """The actuator below the bearing applies the measured load along +y."""
    return np.array([0.0, float(load)])


def force_residual(fx, fy, load):
    return np.array([fx, fy]) + applied_load(load)


def config_has_current_load(config, load):
    values = dict(line.split(" = ", 1) for line in config.splitlines() if " = " in line)
    force = applied_load(load)
    expected = {"external_load_x": force[0], "external_load_y": force[1],
                "external_load_z": 0.0, "load_angle_deg": 90.0}
    return all(float(values.get(key, "nan")) == value for key, value in expected.items())


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


def config_identity(config, ignore=()):
    """Hash physical settings while excluding named output/initial-guess paths."""
    lines = []
    for line in config.splitlines():
        key = line.split(" = ", 1)[0]
        if key == "steady_initial_state_file" and key in ignore:
            continue
        lines.append(f"{key} = <ignored>" if key in ignore else line)
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def warm_start_identity(config):
    """Compare physical inputs for an explicitly requested initial guess only."""
    values = dict(line.split(" = ", 1) for line in config.splitlines() if " = " in line)
    numerical = {"e", "attitude_angle_deg", "fbns_tol", "solution_table_file",
                 "bush_outer_heat_transfer", "shaft_ambient_heat_transfer",
                 "cross_film_thermal_model", "cross_film_layers", "cross_film_profile_order",
                 "cross_film_quadrature_points", "cross_film_relaxation"}
    return {key: value for key, value in values.items()
            if key not in numerical and not key.startswith(("output_", "coupling_", "steady_"))}


def run_solver(command, run_dir, environment, timeout, stop_grace=120):
    """Ask for a saved stop at the deadline; kill the process tree only after grace."""
    with (run_dir / "run.log").open("w", encoding="utf-8") as stream:
        process = subprocess.Popen(command, cwd=run_dir, stdout=stream, stderr=subprocess.STDOUT,
                                   env=environment, start_new_session=os.name != "nt")
        try:
            return process.wait(timeout=timeout), False
        except subprocess.TimeoutExpired:
            output = run_dir / "results"
            output.mkdir(exist_ok=True)
            (output / "STOP").write_text("Driver runtime deadline reached; save before stopping.\n")
            try:
                return process.wait(timeout=stop_grace), True
            except subprocess.TimeoutExpired:
                if os.name == "nt":
                    subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                                   stdout=stream, stderr=subprocess.STDOUT, check=False)
                else:
                    os.killpg(process.pid, signal.SIGKILL)
                return process.wait(), True


def mpi_command(mpi, solver, ranks, rankfile=None, platform=sys.platform):
    command = [str(mpi)]
    if rankfile:
        if not platform.startswith("linux"):
            raise ValueError("--rankfile requires OpenMPI on Linux")
        rankfile = Path(rankfile).resolve()
        if not rankfile.is_file():
            raise FileNotFoundError(f"OpenMPI rankfile not found: {rankfile}")
        command += ["--rankfile", str(rankfile), "--bind-to", "core", "--report-bindings"]
    return command + ["-n", str(ranks), str(solver), "--config", "config.txt"]


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


THERMAL_BOUNDARY_KEYS = (
    "bc_z_south_temperature", "bc_z_north_temperature", "groove_supply_temperature",
    "bush_outer_temperature", "shaft_ambient_temperature",
)


def boundary_temperatures(case):
    """Resolve physical reservoir temperatures in kelvin, leaving TC0 references alone."""
    values = dict.fromkeys(THERMAL_BOUNDARY_KEYS, float(case["tc0_boundary_temperature_k"]))
    overrides = case.get("thermal_boundaries_K", {})
    unknown = set(overrides) - set(values)
    if unknown:
        raise ValueError(f"unknown thermal boundary keys: {sorted(unknown)}")
    for key, value in overrides.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{key} must be a numeric temperature in kelvin")
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{key} must be finite and positive in kelvin")
        values[key] = float(value)
    return values


def render_config(case, experiment, boundaries, materials, eps, psi, h_bush, h_shaft,
                  output_dir, table, coupling_max_iters=2000, coupling_relaxation=1.0,
                  initial_state_file=None, native_equilibrium=True):
    geometry = experiment["geometry"]
    radius = float(geometry["journal_radius_m"])
    clearance = float(geometry["radial_clearance_m"])
    temperature = float(case["tc0_boundary_temperature_k"])
    temperatures = boundary_temperatures(case)
    load = float(case["load_mean_n"])
    load_vector = applied_load(load)
    omega = float(case["angular_speed_magnitude_rad_s"])
    rho = density(temperature - 273.15)
    nu = viscosity(temperature - 273.15, experiment["viscosity"]["points"])
    values = {
        "R": radius, "c": clearance, "e": eps * clearance,
        "L": geometry["bearing_length_m"], "attitude_angle_deg": psi,
        "load_angle_deg": 90, "omega": omega, "n_theta_global": 720,
        "n_z_global": 192, "solution_mode": "STEADY_STATE", "dt": 2e-5,
        "end_t": 2e-5, "write_interval": 2e-5, "diagnostics_interval": 1,
        "motion_model": "STEADY_EQUILIBRIUM" if native_equilibrium else "STATIC",
        "steady_force_tolerance": FORCE_BALANCE_TOL,
        "steady_pressure_inner_only": "true" if native_equilibrium else "false",
        "steady_output_interval": 50,
        "steady_position_step_limit": .05, "steady_position_fd_step": .001,
        "theta_convection_scheme": "UPWIND", "thermal_convection_scheme": "UPWIND",
        "gas_convection_scheme": "UPWIND", "external_load_x": load_vector[0],
        "external_load_y": load_vector[1], "external_load_z": 0,
        "cavitation_model": "FBNS", "fbns_tol": FBNS_TOL,
        "cavitation_threshold": "SCALAR_PCAV",
        "cavitated_film_model": "STRIATED", "gas_pressure_coupling": "NONE",
        "p_cav": materials["cavitation_pressure_Pa"],
        "bc_z_south_type": "DIRICHLET", "bc_z_south_val": materials["reservoir_pressure_Pa"],
        "bc_z_north_type": "DIRICHLET", "bc_z_north_val": materials["reservoir_pressure_Pa"],
        "bc_z_south_thermal": "CONSTANT", "bc_z_north_thermal": "CONSTANT",
        "bc_z_south_temperature": temperatures["bc_z_south_temperature"],
        "bc_z_north_temperature": temperatures["bc_z_north_temperature"],
        "temperature_model": "ENERGY_EQUATION", "temperature_initial": temperature + 1,
        "temperature_reference": temperature, "property_reference_temperature": temperature,
        "rho_cp": rho * float(materials["oil_cp_J_kg_K"]),
        "liquid_specific_heat": materials["oil_cp_J_kg_K"],
        "thermal_conductivity": materials["oil_conductivity_W_m_K"],
        "cross_film_thermal_model": "FINITE_VOLUME", "cross_film_layers": 16,
        "cross_film_profile_order": 4,
        "cross_film_quadrature_points": 6, "cross_film_relaxation": 0.3,
        "bush_thermal_model": "CONDUCTION_3D", "bush_axial_cells": 16,
        "bush_end_natural_convection": "true",
        "bush_thermal_conductivity": materials["bush_conductivity_W_m_K"],
        "bush_outer_radius": boundaries["housing"]["outer_radius_m"],
        "bush_radial_cells": 16, "bush_outer_boundary": "CONVECTION",
        "bush_outer_heat_transfer": h_bush, "shaft_ambient_heat_transfer": h_shaft,
        "shaft_ambient_area_ratio": materials["shaft_end_area_ratio"],
        "thermal_ambient_temperature": temperature,
        "shaft_ambient_temperature": temperatures["shaft_ambient_temperature"],
        "bush_outer_temperature": temperatures["bush_outer_temperature"],
        "fluid_property_model": "SINGLE_PHASE", "dissolved_gas_species": "AIR",
        "dissolved_gas_initial": 0, "solubility_model": "HENRY",
        "dissolved_gas_henry_coeff": 0, "density_model": "LINEAR_TEMPERATURE",
        "density_temperature_intercept": 954.9, "density_temperature_slope": -0.7284,
        "rho": rho, "liquid_viscosity_model": "KINEMATIC_TABLE",
        "solution_table_file": str(table.resolve()), "solution_table_clamp": "false",
        "kinematic_viscosity_log_interpolation": "true", "viscosity_pressure_coeff": 0,
        "mu": rho * nu, "groove_mixing_model": "FLUX_BALANCED",
        "inlet_groove": f"90 {boundaries['groove']['width_deg']} {materials['reservoir_pressure_Pa']} {temperatures['groove_supply_temperature']}",
        "coupling_max_iters": coupling_max_iters, "coupling_tolerance": 1e-6,
        "coupling_relaxation_mode": "FIXED", "coupling_relaxation": coupling_relaxation,
        "coupling_relaxation_min": 0.1, "coupling_relaxation_backoff": 0.5,
        "coupling_relaxation_growth": 1.05, "coupling_keep_best_iterate": "true",
        "log_outer_iters": "true", "log_coupling_history": "true",
        "output_dir": str(output_dir.resolve()),
        "output_write_3d": "false", "output_write_flat": "true",
        "output_verbosity": "DETAILED",
        "output_fields": "pressure,film_content,h,temperature,mu,temperature_journal,temperature_bearing,heat_flux_journal,heat_flux_bearing,fluid_force_x,fluid_force_y,fluid_force_z,pressure_force_x,pressure_force_y,friction_torque,bearing_x,bearing_y,bearing_z,temperature_housing_outer," + ",".join(f"cross_film_solid_{k}" for k in range(16)),
    }
    if case.get("sensitivity_diagnostics", False):
        values["output_fields"] += ",rho_liquid_solution,cross_film_f0,cross_film_f1,cross_film_f2,cross_film_couette_h"
    if initial_state_file is not None:
        values["steady_initial_state_file"] = str(Path(initial_state_file).resolve())
    return "# Generated by run_steady.py; steady CFD and optional native position equilibrium.\n" + "\n".join(
        f"{key} = {value}" for key, value in values.items()) + "\n"


def read_fields(run_dir, results_dir=None):
    pvd = (results_dir or run_dir / "results") / "flat/results.pvd"
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


def final_position(fields, clearance):
    """Read the constant native bearing vector; attitude uses e=-bearing position."""
    if "xB" not in fields or not len(fields["xB"]) or len(fields["xB"]) % 3:
        raise ValueError("missing bearing position vector xB")
    position = fields["xB"].reshape(-1, 3)
    if not np.isfinite(position).all() or not np.allclose(position, position[0], rtol=1e-9, atol=1e-12):
        raise ValueError("nonfinite or inconsistent bearing position vector")
    if clearance is None or not math.isfinite(clearance) or clearance <= 0:
        raise ValueError("invalid radial clearance for bearing position")
    x, y, z = map(float, position[0])
    eps = math.hypot(x, y) / clearance
    if not 0 <= eps < 1 or abs(z) > 1e-12:
        raise ValueError("inadmissible steady bearing position")
    return {"bearing_x_m": x, "bearing_y_m": y, "eps": eps,
            "psi_deg": math.degrees(math.atan2(-y, -x)) % 360}


def check_native_equilibrium(log, force_error):
    """Require the native final success marker and independent VTK load balance."""
    markers = re.findall(r"Equilibrium: (CONVERGED|NOT CONVERGED) R_F=" + NUMBER, log)
    if not markers or markers[-1][0] != "CONVERGED":
        raise ValueError("native steady equilibrium did not converge")
    native_error = float(markers[-1][1])
    if not math.isfinite(native_error) or not 0 <= native_error < FORCE_BALANCE_TOL:
        raise ValueError("native steady force residual above acceptance limit")
    if not math.isfinite(force_error) or not 0 <= force_error < FORCE_BALANCE_TOL:
        raise ValueError("independent VTK total-force residual above acceptance limit")
    return native_error


def check_solution(run_dir, log, returncode, load, cavity_pressure, minimum_boundary_k, viscosity_points,
                   native_equilibrium=False, clearance=None):
    if returncode:
        raise ValueError(f"solver exit code {returncode}")
    diagnostic = last_diagnostic(run_dir)
    for key in ("liquid_residual", "cavitated_fraction"):
        if not math.isfinite(float(diagnostic[key])):
            raise ValueError(f"nonfinite or missing {key} diagnostic")
    for key in ("liquid_residual_kg_s", "liquid_balance_scale_kg_s", "liquid_residual_rel"):
        if not diagnostic.get(key, "").strip() or not math.isfinite(float(diagnostic[key])):
            raise ValueError(f"missing or nonfinite dt-independent mass diagnostic {key}")
    mass_scale = float(diagnostic["liquid_balance_scale_kg_s"])
    mass_relative = float(diagnostic["liquid_residual_rel"])
    if mass_scale <= 0 or not math.isclose(
            mass_relative, float(diagnostic["liquid_residual_kg_s"]) / mass_scale,
            rel_tol=1e-8, abs_tol=1e-12):
        raise ValueError("inconsistent dt-independent liquid balance diagnostic")
    if abs(mass_relative) > MASS_BALANCE_TOL:
        raise ValueError(f"liquid balance relative residual above {MASS_BALANCE_TOL:g}")
    energy_keys = ("energy_balance_residual_W", "energy_balance_scale_W",
                   "energy_viscous_heat_W", "energy_wall_heat_W",
                   "energy_axial_enthalpy_out_W", "energy_groove_enthalpy_in_W",
                   "energy_mass_imbalance_correction_W", "energy_feed_mass_kg_s")
    for key in energy_keys:
        if not diagnostic.get(key, "").strip() or not math.isfinite(float(diagnostic[key])):
            raise ValueError(f"missing or nonfinite global energy diagnostic {key}")
    energy_scale = float(diagnostic["energy_balance_scale_W"])
    if energy_scale < 0:
        raise ValueError("negative global energy balance scale")
    energy_relative = abs(float(diagnostic["energy_balance_residual_W"])) / max(energy_scale, 1.)
    if energy_relative > ENERGY_BALANCE_TOL:
        raise ValueError(f"global energy balance relative residual above {ENERGY_BALANCE_TOL:g}")
    for key in ("cross_film_profile_t_min_K", "cross_film_profile_t_max_K",
                "cross_film_profile_residual"):
        if not diagnostic.get(key, "").strip() or not math.isfinite(float(diagnostic[key])):
            raise ValueError(f"missing or nonfinite cross-film profile diagnostic {key}")
    profile_min = float(diagnostic["cross_film_profile_t_min_K"])
    profile_max = float(diagnostic["cross_film_profile_t_max_K"])
    profile_residual = float(diagnostic["cross_film_profile_residual"])
    if profile_min < minimum_boundary_k - TEMPERATURE_FLOOR_TOL_K or profile_max < profile_min:
        raise ValueError("cross-film temperature profile below coldest physical reservoir or inverted")
    if not 0 <= profile_residual < 1e-6:
        raise ValueError("cross-film profile residual above 1e-6")
    table_min = float(viscosity_points[0]["temperature_c"]) + 273.15
    table_max = float(viscosity_points[-1]["temperature_c"]) + 273.15
    if profile_min < table_min or profile_max > table_max:
        raise ValueError("cross-film temperature profile outside FW68D viscosity table")
    if density(profile_max - 273.15) <= 0:
        raise ValueError("nonpositive density at hottest cross-film profile point")
    for key in ("coupling_R_p", "coupling_R_T"):
        if (not math.isfinite(float(diagnostic[key])) or float(diagnostic[key]) < 0 or
                ((key == "coupling_R_T" or not native_equilibrium) and float(diagnostic[key]) >= 1e-6)):
            raise ValueError(f"{key} failed: {diagnostic[key]}")
    if int(diagnostic["coupling_converged"]) != 1 or int(diagnostic["coupling_restored_best"]) != 0:
        raise ValueError("pressure/thermal coupling did not converge cleanly")
    if int(diagnostic["converged"]) != 1:
        raise ValueError("FBNS did not converge")
    fbns = re.findall(r"fbns iter=\d+ \|G\|=" + NUMBER + r" \|FB\|=" + NUMBER, log)
    if not fbns or max(map(float, fbns[-1])) >= FBNS_TOL:
        raise ValueError(f"FBNS residual missing or not below configured {FBNS_TOL:g}")
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
    if float(np.min(fields["mu"])) <= 0:
        raise ValueError("nonpositive dynamic viscosity")
    rho = 954.9 - 0.7284 * (fields["T"] - 273.15)
    if not np.isfinite(rho).all() or float(np.min(rho)) <= 0:
        raise ValueError("nonpositive or nonfinite density")
    for key in ("T", "T_bearing", "T_journal"):
        if float(np.min(fields[key])) < minimum_boundary_k - TEMPERATURE_FLOOR_TOL_K:
            raise ValueError(f"{key} below coldest physical reservoir by more than {TEMPERATURE_FLOOR_TOL_K:g} K")
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
    force_error = float(np.linalg.norm(force_residual(fx, fy, load))) / load
    position = final_position(fields, clearance) if native_equilibrium else {}
    if native_equilibrium:
        position["native_force_error_rel"] = check_native_equilibrium(log, force_error)
    return {
        **position,
        "Fx_N": fx, "Fy_N": fy, "Fz_N": fz, "torque_Nm": torque,
        "force_error_rel": force_error,
        "h_min_m": float(np.min(fields["h"])),
        "mu_min_Pa_s": float(np.min(fields["mu"])),
        "rho_min_kg_m3": float(np.min(rho)),
        "minimum_boundary_temperature_K": minimum_boundary_k,
        "temperature_floor_tolerance_K": TEMPERATURE_FLOOR_TOL_K,
        "mass_balance_tolerance": MASS_BALANCE_TOL,
        "energy_balance_tolerance": ENERGY_BALANCE_TOL,
        "energy_balance_residual_rel": energy_relative,
        "pressure_acceptance": ("converged_FBNS_below_1e-6; coupling_R_p_reported_only"
                                if native_equilibrium else "FBNS_and_coupling_R_p_below_1e-6"),
        "cross_film_profile_t_min_K": profile_min,
        "cross_film_profile_t_max_K": profile_max,
        "cross_film_profile_residual": profile_residual,
        "p_min_Pa": float(np.min(fields["p"])), "p_max_Pa": float(np.max(fields["p"])),
        "cavitated_area_fraction": float(np.mean(fields["film_content"] < .999)),
        "film_content_min": float(np.min(fields["film_content"])),
        "T_min_K": float(np.min(fields["T"])), "T_max_K": float(np.max(fields["T"])),
        "T_bore_min_K": float(np.min(fields["T_bearing"])),
        "T_bore_mean_K": float(np.mean(fields["T_bearing"])),
        "T_bore_max_K": float(np.max(fields["T_bearing"])),
        "T_outer_mean_K": float(np.mean(fields["temperature_housing_outer"])),
        "T_shaft_mean_K": float(np.mean(fields["T_journal"])),
        "energy_diagnostic_status": "global_balance_verified",
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
        balance_setup = self.setup["force_balance"]
        if (balance_setup["applied_load_unit_vector"] != [0, 1] or
                balance_setup["target_fluid_force_unit_vector"] != [0, -1] or
                balance_setup["load_direction_deg"] != 90 or
                self.boundaries["groove"]["theta_deg"] != 90):
            raise ValueError("setup must use upward actuator load, downward fluid target, and top groove")
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
            pressure = case["load_mean_n"] / (case["journal_diameter_m"] * case["bearing_length_m"])
            case["sommerfeld_number"] = (case["mu_pa_s"] * case["rpm_nominal"] / 60 / pressure
                                         * (case["journal_diameter_m"] / (2 * case["radial_clearance_m"]))**2)
            self.cases[case["case_id"]] = case
        overrides = self.setup.get("thermal_boundaries_by_case", {})
        if not isinstance(overrides, dict) or set(overrides) - set(self.cases):
            raise ValueError("thermal_boundaries_by_case must map known case IDs to kelvin values")
        for key, values in overrides.items():
            if not isinstance(values, dict):
                raise ValueError(f"case {key} thermal boundaries must be an object")
            self.cases[key]["thermal_boundaries_K"] = values
        for case in self.cases.values():
            case["sensitivity_diagnostics"] = bool(self.setup.get("sensitivity_scenario"))
            for temperature in boundary_temperatures(case).values():
                viscosity(temperature - 273.15, self.experiment["viscosity"]["points"])
        self.points = self.experiment["viscosity"]["points"]
        self.solver = Path(args.solver).resolve() if args.solver else None
        self.work = Path(args.work).resolve() if args.work else None
        self.native_equilibrium = not bool(args.trial)
        self.driver_hash = sha256(__file__)
        self.rankfile_hash = sha256(Path(args.rankfile).resolve()) if args.rankfile else None
        self.sequence = 0
        self.h_warm = {}
        self.latest_checkpoint = {}
        self.restart_positions = {}
        self.restart_provenance = {}
        self.source_hash = source_sha256()
        self.executable_hash = sha256(self.solver) if self.solver and self.solver.is_file() else None
        self.input_hashes = {"experiment_sha256": sha256(EXPERIMENT),
                             "material_setup_sha256": sha256(self.materials_path),
                             "property_table_sha256": sha256(TABLE)}
        self.completed = {key: [] for key in self.cases}
        if self.work and self.work.is_dir():
            for path in self.work.glob("*/*/metrics.json"):
                try:
                    record = json.loads(path.read_text())
                    if record.get("case_id") in self.completed:
                        self.completed[record["case_id"]].append((path.parent, record))
                except (OSError, ValueError):
                    continue
        for records in self.completed.values():
            records.sort(key=lambda item: item[0].stat().st_mtime)
        self.rechecked = {}
        self.mpi_version = None

    def matching_provenance(self, record):
        expected = {"source_sha256": self.source_hash,
                    "executable_sha256": self.executable_hash,
                    "mpi_ranks": self.args.ranks, "driver_sha256": self.driver_hash,
                    "rankfile_sha256": self.rankfile_hash,
                    "native_equilibrium": self.native_equilibrium, **self.input_hashes}
        return (record.get("returncode") == 0 and
                all(key in record and record[key] == value for key, value in expected.items()))

    def recheck(self, run_dir, record, case):
        if not self.matching_provenance(record):
            return None
        if run_dir in self.rechecked:
            return self.rechecked[run_dir]
        try:
            config_path = run_dir / "config.txt"
            if (record.get("config_sha256") != sha256(config_path) or
                    not config_has_current_load(config_path.read_text(), case["load_mean_n"])):
                return None
            checked = check_solution(run_dir, (run_dir / "run.log").read_text(errors="replace"),
                                     0, float(case["load_mean_n"]),
                                     float(self.materials["cavitation_pressure_Pa"]),
                                     min(boundary_temperatures(case).values()), self.points,
                                     self.native_equilibrium,
                                     float(self.experiment["geometry"]["radial_clearance_m"]))
            result = {**record, **checked, "run_dir": str(run_dir)}
            self.rechecked[run_dir] = result
            return result
        except (OSError, ValueError, KeyError, ET.ParseError):
            self.rechecked[run_dir] = None
            return None

    def config_for(self, case, eps, psi, h_bush, h_shaft, output_dir, initial_state_file=None):
        return render_config(case, self.experiment, self.boundaries, self.materials,
                             eps, psi, h_bush, h_shaft, output_dir, TABLE,
                             self.args.coupling_max_iters, self.args.coupling_relaxation,
                             initial_state_file, self.native_equilibrium)

    def find_completed(self, case, eps, psi, h_bush=None, h_shaft=None, closed=False):
        probe = self.config_for(case, eps, psi, h_bush or 60., h_shaft or 100.,
                                self.work / case["case_id"] / "_probe" / "results")
        ignore = {"output_dir", "steady_initial_state_file"}
        if closed:
            ignore |= {"bush_outer_heat_transfer", "shaft_ambient_heat_transfer"}
        identity = config_identity(probe, ignore)
        for run_dir, record in reversed(self.completed[case["case_id"]]):
            if not self.matching_provenance(record):
                continue
            if closed and record.get("status") not in ("thermally_closed", "load_and_thermal_closed"):
                continue
            if not closed and record.get("status") not in (
                    "valid_inner_solve", "thermally_closed", "load_and_thermal_closed",
                    "trial_only_not_thermally_closed"):
                continue
            config_path = run_dir / "config.txt"
            if (not config_path.is_file() or record.get("config_sha256") != sha256(config_path) or
                    config_identity(config_path.read_text(), ignore) != identity):
                continue
            if closed:
                values = dict(line.split(" = ", 1) for line in config_path.read_text().splitlines()
                              if " = " in line)
                if (float(values["bush_outer_heat_transfer"]) != record.get("h_bush_W_m2_K") or
                        float(values["shaft_ambient_heat_transfer"]) != record.get("h_shaft_W_m2_K")):
                    continue
            result = self.recheck(run_dir, record, case)
            if result is None:
                continue
            if closed:
                if result.get("h_relative_mismatch") is None or result.get("T_change_K") is None:
                    continue
                if not (float(result["h_relative_mismatch"]) < 1e-3 and
                        float(result["T_change_K"]) < .01):
                    continue
                temperatures = boundary_temperatures(case)
                h_b = float(result["h_bush_W_m2_K"])
                h_s = float(result["h_shaft_W_m2_K"])
                outer = result["T_outer_mean_K"]
                target_b, _ = housing_coefficient(outer, temperatures["bush_outer_temperature"], self.points, self.materials)
                target_s, _ = shaft_coefficient(result["T_shaft_mean_K"], temperatures["shaft_ambient_temperature"],
                                                case["angular_speed_magnitude_rad_s"],
                                                self.points, self.materials)
                if max(abs(target_b / h_b - 1), abs(target_s / h_s - 1)) >= 1e-3:
                    continue
                self.h_warm[case["case_id"]] = (h_b, h_s)
            checkpoint = run_dir / "results/steady_state.dat"
            if checkpoint.is_file():
                self.latest_checkpoint[case["case_id"]] = checkpoint
            return result
        return None

    def command(self):
        mpi = shutil.which("mpiexec") or shutil.which("mpirun")
        if not mpi:
            raise FileNotFoundError("mpiexec/mpirun not found")
        return mpi_command(mpi, self.solver, self.args.ranks, self.args.rankfile)

    def recover_checkpoint(self, case):
        key = case["case_id"]
        if key in self.latest_checkpoint:
            return self.latest_checkpoint[key]
        for run_dir, record in reversed(self.completed[key]):
            checkpoint = run_dir / "results/steady_state.dat"
            if (checkpoint.is_file() and self.recheck(run_dir, record, case) is not None and
                    (record.get("output_state_sha256") is None or
                     record["output_state_sha256"] == sha256(checkpoint))):
                self.latest_checkpoint[key] = checkpoint
                return checkpoint
        initial_root = getattr(self.args, "initial_state_root", None)
        if initial_root:
            root = Path(initial_root).resolve(strict=True)
            expected = self.config_for(case, .3, 270., 60., 100., self.work / key / "_probe")
            candidates = sorted((root / key).glob("*/results/**/steady_state.dat"),
                                key=lambda path: path.stat().st_mtime, reverse=True)
            for checkpoint in candidates:
                run_dir = next((parent for parent in checkpoint.parents if (parent / "config.txt").is_file()), None)
                if run_dir is None:
                    continue
                try:
                    record = json.loads((run_dir / "metrics.json").read_text())
                    config = (run_dir / "config.txt").read_text()
                    if (record.get("case_id") != key or
                            record.get("config_sha256") != sha256(run_dir / "config.txt") or
                            record.get("experiment_sha256") != self.input_hashes["experiment_sha256"] or
                            record.get("property_table_sha256") != self.input_hashes["property_table_sha256"] or
                            not config_has_current_load(config, case["load_mean_n"]) or
                            warm_start_identity(config) != warm_start_identity(expected)):
                        continue
                    if (checkpoint.parent == run_dir / "results" and record.get("output_state_sha256") and
                            record["output_state_sha256"] != sha256(checkpoint)):
                        continue
                    position = final_position(read_fields(run_dir, checkpoint.parent),
                                              self.experiment["geometry"]["radial_clearance_m"])
                    values = dict(line.split(" = ", 1) for line in config.splitlines() if " = " in line)
                    self.latest_checkpoint[key] = checkpoint
                    self.restart_positions[key] = position["eps"], position["psi_deg"]
                    self.h_warm[key] = float(values["bush_outer_heat_transfer"]), float(values["shaft_ambient_heat_transfer"])
                    self.restart_provenance[key] = {"run_dir": str(run_dir), "prior_status": record.get("status"),
                                                    "initial_guess_only": True, "source_sha256": record.get("source_sha256")}
                    return checkpoint
                except (OSError, ValueError, KeyError, ET.ParseError):
                    continue
        return None

    def single(self, case, eps, psi, h_bush, h_shaft, label):
        cached = self.find_completed(case, eps, psi, h_bush, h_shaft)
        if cached is not None:
            print(f"reusing completed inner solve {cached['run_dir']}", flush=True)
            return cached
        self.sequence += 1
        run_dir = self.work / case["case_id"] / f"{self.sequence:04d}_{label}_{uuid.uuid4().hex[:8]}"
        run_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = self.recover_checkpoint(case)
        eps, psi = self.restart_positions.pop(case["case_id"], (eps, psi))
        config = self.config_for(case, eps, psi, h_bush, h_shaft,
                                 run_dir / "results", checkpoint)
        config_path = run_dir / "config.txt"
        config_path.write_text(config, encoding="utf-8")
        command = self.command()
        record = {
            "case_id": case["case_id"], "eps": eps, "psi_deg": psi,
            "seed_eps": eps, "seed_psi_deg": psi, "native_equilibrium": self.native_equilibrium,
            "thermal_boundaries_K": boundary_temperatures(case),
            "h_bush_W_m2_K": h_bush, "h_shaft_W_m2_K": h_shaft,
            "run_dir": str(run_dir), "command": command,
            "mpi_ranks": self.args.ranks, "python_version": sys.version.split()[0],
            "mpi_version": self.mpi_version, "source_sha256": self.source_hash,
            "driver_sha256": self.driver_hash,
            "rankfile_sha256": self.rankfile_hash,
            "config_sha256": sha256(config_path), "executable_sha256": self.executable_hash,
            "initial_state_file": str(checkpoint) if checkpoint else None,
            "initial_state_sha256": sha256(checkpoint) if checkpoint else None,
            "restart_provenance": self.restart_provenance.get(case["case_id"]),
            **self.input_hashes, "status": "started",
        }
        save_json(run_dir / "metrics.json", record)
        started = time.monotonic()
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                           MKL_NUM_THREADS="1")
        try:
            returncode, deadline_reached = run_solver(command, run_dir, environment, self.args.timeout)
            record["returncode"] = returncode
            record["runtime_deadline_reached"] = deadline_reached
            record["wall_seconds"] = time.monotonic() - started
            state_path = run_dir / "results/steady_state.dat"
            if state_path.is_file():
                record["output_state_sha256"] = sha256(state_path)
            if deadline_reached:
                raise RuntimeError("runtime deadline reached; stop requested and saved outputs retained")
            log = (run_dir / "run.log").read_text(errors="replace")
            record.update(check_solution(run_dir, log, returncode,
                                         float(case["load_mean_n"]), float(self.materials["cavitation_pressure_Pa"]),
                                         min(boundary_temperatures(case).values()), self.points,
                                         self.native_equilibrium,
                                         float(self.experiment["geometry"]["radial_clearance_m"])))
            record["status"] = "valid_inner_solve"
        except Exception as error:
            record["wall_seconds"] = time.monotonic() - started
            record["status"] = "failed"
            record["error"] = str(error)
            save_json(run_dir / "metrics.json", record)
            raise RuntimeError(f"solver evaluation failed; retained {run_dir}: {error}") from error
        save_json(run_dir / "metrics.json", record)
        self.completed[case["case_id"]].append((run_dir, record))
        state_path = run_dir / "results/steady_state.dat"
        if state_path.is_file():
            record["output_state_sha256"] = sha256(state_path)
            self.latest_checkpoint[case["case_id"]] = state_path
            save_json(run_dir / "metrics.json", record)
        return record

    def evaluate(self, case, eps, psi, label, thermal_iterations=None):
        if not .0 < eps < .995:
            raise ValueError("eccentricity ratio must lie in (0, 0.995)")
        if not math.isfinite(psi):
            raise ValueError("attitude angle must be finite")
        if thermal_iterations is None:
            cached = self.find_completed(case, eps, psi, closed=True)
            if cached is not None:
                print(f"reusing thermally closed position {cached['run_dir']}", flush=True)
                return cached
        temperatures = boundary_temperatures(case)
        bush_bath = temperatures["bush_outer_temperature"]
        shaft_bath = temperatures["shaft_ambient_temperature"]
        omega = float(case["angular_speed_magnitude_rad_s"])
        if getattr(self.args, "initial_state_root", None):
            self.recover_checkpoint(case)
            eps, psi = self.restart_positions.pop(case["case_id"], (eps, psi))
        if case["case_id"] in self.h_warm:
            h_bush, h_shaft = self.h_warm[case["case_id"]]
        else:
            h_bush, _ = housing_coefficient(bush_bath + 5, bush_bath, self.points, self.materials)
            h_shaft, _ = shaft_coefficient(shaft_bath + 5, shaft_bath, omega, self.points, self.materials)
        previous_outer = None
        previous_shaft = None
        previous_mismatch = None
        count = thermal_iterations if thermal_iterations is not None else self.args.max_thermal_iterations
        for iteration in range(count):
            result = self.single(case, eps, psi, h_bush, h_shaft, f"{label}_thermal{iteration + 1}")
            outer = result["T_outer_mean_K"]
            shaft = result["T_shaft_mean_K"]
            target_bush, housing_meta = housing_coefficient(outer, temperatures["bush_outer_temperature"], self.points, self.materials)
            target_shaft, shaft_meta = shaft_coefficient(shaft, shaft_bath, omega, self.points, self.materials)
            rel = max(abs(target_bush / h_bush - 1), abs(target_shaft / h_shaft - 1))
            delta_temperature = (math.inf if previous_outer is None else
                                 max(abs(outer - previous_outer), abs(shaft - previous_shaft)))
            result.update(T_outer_mean_K=outer, T_shaft_K=shaft,
                          h_bush_target_W_m2_K=target_bush, h_shaft_target_W_m2_K=target_shaft,
                          h_relative_mismatch=rel, T_change_K=delta_temperature if math.isfinite(delta_temperature) else None,
                          housing_correlation=housing_meta, shaft_correlation=shaft_meta)
            if rel < 1e-3 and delta_temperature < .01:
                result["status"] = "load_and_thermal_closed" if self.native_equilibrium else "thermally_closed"
                self.h_warm[case["case_id"]] = (h_bush, h_shaft)
                save_json(Path(result["run_dir"]) / "metrics.json", result)
                return result
            if thermal_iterations is not None and iteration == count - 1:
                result["status"] = "trial_only_not_thermally_closed"
                save_json(Path(result["run_dir"]) / "metrics.json", result)
                return result
            save_json(Path(result["run_dir"]) / "metrics.json", result)
            if self.native_equilibrium:
                eps, psi = result["eps"], result["psi_deg"]
            previous_outer, previous_shaft = outer, shaft
            relaxation = .5 if previous_mismatch is not None and rel > 1.25 * previous_mismatch else 1.
            h_bush = float(np.clip(h_bush + relaxation * (target_bush - h_bush), h_bush / 2, 2 * h_bush))
            h_shaft = float(np.clip(h_shaft + relaxation * (target_shaft - h_shaft), h_shaft / 2, 2 * h_shaft))
            previous_mismatch = rel
        result["status"] = "thermal_closure_failed"
        save_json(Path(result["run_dir"]) / "metrics.json", result)
        raise RuntimeError(f"thermal closure cap reached for {case['case_id']}; last run {result['run_dir']}")

    def balance(self, case, eps, psi):
        """Native position updates share the pressure/energy loop; only h closes here."""
        result = self.evaluate(case, eps, psi, "equilibrium")
        save_json(self.work / case["case_id"] / f"result_{uuid.uuid4().hex[:8]}.json", result)
        return result


def self_check(campaign):
    expected = set("ABCDEF")
    if not set(campaign.cases).issuperset(expected):
        raise ValueError("expected exact A-F experimental cases")
    for case in campaign.cases.values():
        eps, psi = campaign.setup["force_balance"]["initial_positions"][case["case_id"]]
        if not 0 < eps < .995 or not math.isfinite(psi) or math.sin(math.radians(psi)) >= 0:
            raise ValueError("initial seed must put bearing centre above fixed journal centre")
        temperature = case["tc0_boundary_temperature_c"]
        nu = viscosity(temperature, campaign.points)
        mu = density(temperature) * nu
        if not math.isclose(mu, case["mu_pa_s"], rel_tol=2e-12):
            raise ValueError(f"case {case['case_id']} viscosity mismatch")
        if not math.isclose(case["angular_speed_magnitude_rad_s"],
                            case["rpm_nominal"] * math.pi / 30, rel_tol=2e-12):
            raise ValueError(f"case {case['case_id']} angular speed mismatch")
        if case["case_id"] in expected or case["case_id"] in ("E_low", "E_high"):
            pressure = case["load_mean_n"] / (case["journal_diameter_m"] * case["bearing_length_m"])
            sommerfeld = mu * (case["rpm_nominal"] / 60) / pressure * (
                case["journal_diameter_m"] / (2 * case["radial_clearance_m"]))**2
            if not math.isclose(sommerfeld, case["sommerfeld_number"], rel_tol=2e-12):
                raise ValueError(f"case {case['case_id']} Sommerfeld mismatch")
        amb = case["tc0_boundary_temperature_k"]
        temperatures = boundary_temperatures(case)
        config = render_config(case, campaign.experiment, campaign.boundaries, campaign.materials, .3, 140,
                               60, 100, Path("/tmp/steady-check"), TABLE)
        warmed = render_config(case, campaign.experiment, campaign.boundaries, campaign.materials, .3, 140,
                               61, 101, Path("/tmp/other-output"), TABLE,
                               initial_state_file=Path("/tmp/previous-state.dat"))
        path_keys = {"output_dir", "steady_initial_state_file"}
        if config_identity(config, path_keys) == config_identity(warmed, path_keys):
            raise ValueError("exact config identity discarded physical h settings")
        if config_identity(config, path_keys | {"bush_outer_heat_transfer", "shaft_ambient_heat_transfer"}) != (
                config_identity(warmed, path_keys | {"bush_outer_heat_transfer", "shaft_ambient_heat_transfer"})):
            raise ValueError("closed-position identity changed with h or warm-start path")
        values = dict(line.split(" = ", 1) for line in config.splitlines() if " = " in line)
        if (not config_has_current_load(config, case["load_mean_n"]) or
                "solution_mode = STEADY_STATE" not in config or
                values["motion_model"] != "STEADY_EQUILIBRIUM" or
                values["steady_pressure_inner_only"] != "true" or
                int(values["steady_output_interval"]) != 50 or
                float(values["steady_force_tolerance"]) != FORCE_BALANCE_TOL):
            raise ValueError("steady or load config mismatch")
        old_config = config.replace("external_load_y = " + str(float(case["load_mean_n"])),
                                    "external_load_y = " + str(-float(case["load_mean_n"])))
        if (config_has_current_load(old_config, case["load_mean_n"]) or
                config_identity(config, path_keys) == config_identity(old_config, path_keys)):
            raise ValueError("old downward-load config was not rejected")
        if not np.array_equal(force_residual(0.0, -case["load_mean_n"], case["load_mean_n"]), [0, 0]):
            raise ValueError("fluid target does not oppose upward applied load")
        if (float(values["temperature_initial"]) != amb + 1 or
                values["coupling_max_iters"] != "2000" or
                values["coupling_relaxation_mode"] != "FIXED" or
                values["coupling_relaxation"] != "1.0" or
                values["coupling_keep_best_iterate"] != "true" or
                float(values["fbns_tol"]) != FBNS_TOL or
                values["groove_mixing_model"] != "FLUX_BALANCED" or
                values["cross_film_thermal_model"] != "FINITE_VOLUME" or
                int(values["cross_film_layers"]) != 16 or
                float(values["cross_film_relaxation"]) != .3 or
                float(values["liquid_specific_heat"]) != float(campaign.materials["oil_cp_J_kg_K"])):
            raise ValueError("initial guess, boundary, or fixed-point defaults mismatch")
        for key in ("bc_z_south_temperature", "bc_z_north_temperature",
                    "shaft_ambient_temperature", "bush_outer_temperature"):
            if float(values[key]) != temperatures[key]:
                raise ValueError(f"{key} differs from the requested physical reservoir")
        for key in ("thermal_ambient_temperature", "temperature_reference", "property_reference_temperature"):
            if float(values[key]) != amb:
                raise ValueError(f"{key} numerical reference changed from TC0")
        groove = list(map(float, values["inlet_groove"].split()))
        if (len(groove) != 4 or groove[3] != temperatures["groove_supply_temperature"] or
                float(values["bc_z_south_val"]) != 101325 or
                float(values["bc_z_north_val"]) != 101325 or
                groove[2] != 101325 or float(values["p_cav"]) != 0):
            raise ValueError("groove/axial reservoir and cavity references mismatch")
        outer = outer_temperature(amb + 10, amb, 60, .0215, .0425, 15)
        hb, _ = housing_coefficient(outer, amb, campaign.points, campaign.materials)
        hs, meta = shaft_coefficient(amb + 10, amb, case["angular_speed_magnitude_rad_s"],
                                     campaign.points, campaign.materials)
        if min(hb, hs) <= 0 or not meta["laminar_Re_under_1e5"]:
            raise ValueError("invalid external coefficient or laminar disk assumption")
    position = final_position({"xB": np.tile([0., .6 * 4.25e-5, 0.], 2)}, 4.25e-5)
    if not math.isclose(position["eps"], .6) or position["psi_deg"] != 270:
        raise ValueError("native final position sign convention mismatch")
    check_native_equilibrium("Equilibrium: CONVERGED R_F=0.001 x=0 y=0", .001)
    for marker, error in (("Equilibrium: NOT CONVERGED R_F=0.001", .001),
                          ("Equilibrium: CONVERGED R_F=0.003", .001),
                          ("Equilibrium: CONVERGED R_F=0.001", .003), ("", .001)):
        try:
            check_native_equilibrium(marker, error)
        except ValueError:
            pass
        else:
            raise ValueError("invalid native or VTK force convergence accepted")
    trial = render_config(case, campaign.experiment, campaign.boundaries, campaign.materials,
                          .3, 140, 60, 100, Path("/tmp/trial"), TABLE, native_equilibrium=False)
    if "motion_model = STATIC" not in trial or config_identity(trial) == config_identity(config):
        raise ValueError("fixed-position trial was not isolated from native equilibrium")
    record = {"returncode": 0, "source_sha256": campaign.source_hash,
              "executable_sha256": campaign.executable_hash, "mpi_ranks": campaign.args.ranks,
              "driver_sha256": campaign.driver_hash, "rankfile_sha256": campaign.rankfile_hash,
              "native_equilibrium": campaign.native_equilibrium, **campaign.input_hashes}
    if not campaign.matching_provenance(record):
        raise ValueError("current provenance rejected")
    for key in ("driver_sha256", "rankfile_sha256", "native_equilibrium"):
        old = dict(record)
        old[key] = "old"
        if campaign.matching_provenance(old):
            raise ValueError(f"stale {key} accepted")
    for invalid in ({}, {"xB": np.array([0., 0.])},
                    {"xB": np.array([math.nan, 0., 0.])},
                    {"xB": np.array([4.25e-5, 0., 0.])},
                    {"xB": np.array([0., 0., 0., 1e-5, 0., 0.])}):
        try:
            final_position(invalid, 4.25e-5)
        except ValueError:
            pass
        else:
            raise ValueError("invalid native bearing position accepted")
    with tempfile.TemporaryDirectory() as temporary:
        probe = Campaign.__new__(Campaign)
        probe.args = argparse.Namespace(max_thermal_iterations=2)
        probe.native_equilibrium = True
        probe.h_warm = {}
        probe.points, probe.materials = campaign.points, campaign.materials
        probe.experiment, probe.boundaries = campaign.experiment, campaign.boundaries
        probe.find_completed = lambda *args, **kwargs: None
        seeds = []
        ambient = float(case["tc0_boundary_temperature_k"])
        radius = float(campaign.experiment["geometry"]["journal_radius_m"])
        outer_radius = float(campaign.boundaries["housing"]["outer_radius_m"])
        def fake_single(case, eps, psi, h_bush, h_shaft, label):
            seeds.append((eps, psi))
            directory = Path(temporary) / label
            directory.mkdir()
            resistance = 1 + h_bush * outer_radius / float(probe.materials["bush_conductivity_W_m_K"]) * math.log(outer_radius / radius)
            return {"run_dir": str(directory), "eps": .3, "psi_deg": 270.,
                    "T_bore_mean_K": ambient + 5 * resistance,
                    "T_outer_mean_K": ambient + 5,
                    "T_shaft_mean_K": ambient + 5}
        probe.single = fake_single
        result = probe.evaluate(case, .15, 345, "seed_check")
        if seeds != [(.15, 345), (.3, 270.)] or result["status"] != "load_and_thermal_closed":
            raise ValueError("external h closure did not seed from native solved position")
        rankfile = Path(temporary) / "rankfile"
        rankfile.write_text("rank 0=localhost slot=0\n", encoding="utf-8")
        command = mpi_command("mpirun", "/tmp/pancake", 1, rankfile, platform="linux")
        if command[1:6] != ["--rankfile", str(rankfile.resolve()), "--bind-to", "core", "--report-bindings"]:
            raise ValueError("rankfile MPI command mismatch")
        if command[-4:] != ["1", "/tmp/pancake", "--config", "config.txt"]:
            raise ValueError("rankfile MPI rank/executable command mismatch")
    print("Self-check passed: cases, boundaries, thermal correlations, native position/force acceptance, provenance and trial isolation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", help="compiled pancake executable")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--work", help="new or existing campaign output root")
    parser.add_argument("--cases", default="A,B,C,D,E,F", help="comma-separated A-F or E_low,E_high")
    parser.add_argument("--materials", default=str(SETUP), help="complete setup JSON; optional thermal_boundaries_by_case contains physical kelvin temperatures")
    parser.add_argument("--rankfile", help="existing OpenMPI rankfile; Linux only")
    parser.add_argument("--initial-state-root", help="prior campaign work root; compatible fields seed fresh solves only")
    parser.add_argument("--timeout", type=int, default=7200, help="seconds per native solve")
    parser.add_argument("--max-thermal-iterations", type=int, default=30)
    parser.add_argument("--coupling-max-iters", type=int, default=2000,
                        help="steady pressure/thermal fixed-point cap (default 2000)")
    parser.add_argument("--coupling-relaxation", type=float, default=1.0,
                        help="fixed-point relaxation in (0,1] (default 1)")
    parser.add_argument("--trial", nargs=3, metavar=("CASE", "EPS", "PSI"),
                        help="one fixed-geometry benchmark evaluation")
    parser.add_argument("--thermal-iterations", type=int, default=1,
                        help="trial-only number of evaluations; result is unclosed")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if (args.coupling_max_iters < 1 or args.max_thermal_iterations < 1 or args.timeout <= 0 or
            args.thermal_iterations < 1 or not 0 < args.coupling_relaxation <= 1):
        parser.error("coupling cap must be positive and fixed relaxation in (0,1]")
    campaign = Campaign(args)
    if args.self_check:
        self_check(campaign)
        return
    if not args.solver or not args.work or args.ranks < 1:
        parser.error("physical runs require --solver, --work, and positive --ranks")
    if not campaign.solver.is_file():
        parser.error("solver executable does not exist")
    requested = [args.trial[0]] if args.trial else [key.strip() for key in args.cases.split(",")]
    if campaign.setup.get("sensitivity_scenario") and not set(requested) <= set(campaign.setup["case_ids"]):
        parser.error("requested cases are outside this sensitivity scenario's case_ids")
    campaign.command()
    try:
        campaign.mpi_version = subprocess.run([campaign.command()[0], "--version"],
                                              capture_output=True, text=True, timeout=10).stdout.splitlines()[:2]
    except (OSError, subprocess.TimeoutExpired):
        campaign.mpi_version = ["unavailable"]
    campaign.work.mkdir(parents=True, exist_ok=True)
    starts = campaign.setup["force_balance"]["initial_positions"]
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
