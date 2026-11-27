"""Shared helpers of the 02h compressor OS-bearing campaign: config rendering, solver runs,
and readers for motion.csv and the flat VTK fields."""
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02d_fw68d_steady_six"))
import run_steady as steady_driver  # noqa: E402

OUTPUT_FIELDS = ("pressure,film_content,h,mu,inlet_indicator,fluid_force_x,fluid_force_y,fluid_force_z,"
                 "friction_torque,bearing_x,bearing_y,bearing_z,external_load_x,external_load_y")


def base_values(**overrides):
    """FBNS, isothermal, constant properties, no inlet, DIRICHLET ends. Every caller
    sets the geometry, speed, ends, cavitation datum and motion keys it needs."""
    values = dict(
        omega_ramp_time=0, min_film_thickness=1e-9, stop_on_nonpositive_film="true",
        theta_convection_scheme="UPWIND", thermal_convection_scheme="UPWIND", gas_convection_scheme="UPWIND",
        cavitation_model="FBNS", fbns_tol=1e-10, fbns_max_iters=200, theta_min=1e-6,
        max_outer_iters=300, outer_tol=1e-7, gas_pressure_coupling="NONE", bulk_modulus=1e9,
        bc_z_south_type="DIRICHLET", bc_z_north_type="DIRICHLET",
        bc_z_south_thermal="CONSTANT", bc_z_north_thermal="CONSTANT",
        temperature_model="ISOTHERMAL", property_reference_temperature=293.15,
        property_reference_pressure=1e5, rho_cp=1.8e6, thermal_conductivity=.13,
        cross_film_thermal_model="NONE", journal_heat_transfer=0, bearing_heat_transfer=0,
        fluid_property_model="SINGLE_PHASE", steady_gas_model="FROZEN",
        density_model="CONSTANT", liquid_viscosity_model="CONSTANT",
        viscosity_temperature_coeff=0, viscosity_pressure_coeff=0,
        solution_density_gas_coeff=0, solution_viscosity_gas_coeff=0,
        coupling_max_iters=200, coupling_tolerance=1e-6, coupling_relaxation=1.0,
        log_outer_iters="false", output_dir="results", filename_prefix="solution",
        output_verbosity="DETAILED", output_write_3d="false", output_write_flat="true",
        output_fields=OUTPUT_FIELDS, diagnostics_interval=1)
    values.update(overrides)
    return values


def render(values):
    lines = []
    for key, value in values.items():
        lines.append(f"{key} = {value!r}" if isinstance(value, float) else f"{key} = {value}")
    return "\n".join(lines) + "\n"


def run(run_dir, values, solver, ranks=4, mpi="mpiexec", timeout=24 * 3600, reuse=True, rankfile=None):
    """Write config.txt and run the solver in run_dir. A finished run with the same
    config (sha256) is reused. Returns a record dict."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    text = render(values)
    digest = hashlib.sha256(text.encode()).hexdigest()
    record_path = run_dir / "record.json"
    if reuse and record_path.exists():
        record = json.loads(record_path.read_text())
        if record.get("config_sha256") == digest and record.get("returncode") is not None:
            return record
    (run_dir / "config.txt").write_text(text)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    command = steady_driver.mpi_command(mpi, str(Path(solver).resolve()), ranks, rankfile)
    start = time.monotonic()
    code, timed_out = steady_driver.run_solver(command, run_dir, env, timeout)
    record = dict(config_sha256=digest, returncode=code, timed_out=timed_out,
                  wall_s=time.monotonic() - start, ranks=ranks, solver=str(solver),
                  solver_sha256=steady_driver.sha256(solver))
    record_path.write_text(json.dumps(record, indent=2))
    return record


def read_motion(run_dir):
    """motion.csv as a dict of float arrays."""
    with (Path(run_dir) / "results/motion.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{run_dir}: motion.csv has no rows")
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def read_fields(run_dir):
    return steady_driver.read_fields(Path(run_dir))


def last_log_line(run_dir, marker):
    lines = [line for line in (Path(run_dir) / "run.log").read_text(errors="replace").splitlines() if marker in line]
    return lines[-1] if lines else ""


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def wrap180(angle_deg):
    return (angle_deg + 180.0) % 360.0 - 180.0


def finite(*values):
    return all(math.isfinite(v) for v in values)
