"""Prepare and run matched native COUPLED refrigerant-release scenarios.

All source chart bases and caloric proxies are explicit conditional assumptions.
Preparation never launches a solver; failed runs remain available for inspection.
"""

import argparse
import bisect
import copy
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "cases/02e_lge_mixture"))
import run_mixture as baseline
driver = baseline.driver
parse_config = baseline.parse_config
BRANCHES = ("gas_off", "gas_on")
CALORIC_KEYS = ("dissolved_heat_capacity", "desorption_internal_energy",
                "gas_thermal_conductivity", "property_reference_temperature",
                "property_reference_pressure", "mu_gas")


def load_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def curve_value(curves, axis, temperature_c, logarithmic=False):
    """Mirror native exact-node and strict piecewise-linear curve interpolation."""
    axes = sorted(curves)
    if not axes[0] <= axis <= axes[-1]:
        raise ValueError("pressure/concentration outside digitized range")
    right = bisect.bisect_left(axes, axis)
    contributing = [axes[right]] if axes[right] == axis else axes[right-1:right+1]
    values = []
    for key in contributing:
        points = sorted(curves[key])
        temperatures = [point[0] for point in points]
        if not temperatures[0] <= temperature_c <= temperatures[-1]:
            raise ValueError(f"temperature {temperature_c:g} C outside curve {key:g} range")
        ordinate = [math.log(point[1]) if logarithmic else point[1] for point in points]
        values.append(float(np.interp(temperature_c, temperatures, ordinate)))
    value = values[0] if len(values) == 1 else float(np.interp(axis, contributing, values))
    return math.exp(value) if logarithmic else value


def property_curves(table):
    viscosity, solubility = {}, {}
    for row in load_rows(table):
        concentration, temperature = float(row["c_mass_fraction"]), float(row["T_C"])
        viscosity.setdefault(concentration, []).append((temperature, float(row["viscosity_mm2s"])))
        if row["P_MPa"]:
            solubility.setdefault(float(row["P_MPa"])*1e6, []).append((temperature, concentration))
    return viscosity, solubility


def bubble_pressure(curves, concentration, temperature_c):
    """Invert the native first-crossing isotherm, then verify strict chart support."""
    if concentration <= 0:
        return 0.
    pressures = sorted(curves)
    values = [float(np.interp(temperature_c, *zip(*sorted(curves[p])))) for p in pressures]
    pressure = pressures[-1]
    if concentration <= values[0]:
        pressure = pressures[0]
    elif concentration <= values[-1]:
        for index in range(1, len(pressures)):
            if values[index] < concentration:
                continue
            span = values[index]-values[index-1]
            pressure = (pressures[index] if span <= 0 else pressures[index-1] +
                        (concentration-values[index-1])/span*(pressures[index]-pressures[index-1]))
            break
    if abs(curve_value(curves, pressure, temperature_c)-concentration) > 1e-10:
        raise ValueError("bubble pressure outside strictly supported digitized range")
    return pressure


def dew_pressure(temperature_k):
    """Same Lemmon R410A dew ancillary and limits as the checked-in gas EOS."""
    if not 200 <= temperature_k <= 500:
        raise ValueError("R410A EOS requires 200 <= T <= 500 K")
    if temperature_k >= 344.494:
        return math.inf
    x = 1-temperature_k/344.494
    return 4901200*math.exp(344.494/temperature_k *
        (-7.4411*x+1.9883*x**1.6-2.4925*x**2.4-3.2633*x**5))


def pressure_roundoff(config):
    """Match native CoupledCavitation::pressure_admissible: slack [Pa] below p_v.

    The coupled residual is scaled by p_scale, so a state accepted at
    coupled_nonlinear_tolerance may undershoot p_v by that tolerance times p_scale.
    """
    tolerance = float(config.get("coupled_nonlinear_tolerance", "1e-10"))
    p_scale = max(1e5, float(config.get("bc_z_south_val", "0")), float(config.get("bc_z_north_val", "0")),
                  float(config["mu"])*abs(float(config["omega"]))*float(config["R"])**2/float(config["c"])**2)
    return max(1e-9, min(tolerance, 1e-8))*p_scale


def check_layer_eos(config, fields):
    """Require every resolved film layer and wall to support any local free gas."""
    pressure = fields["p"]
    free = fields["alpha_gas"] > 1e-10
    names = [f"cross_film_layer_{layer}" for layer in range(int(config["cross_film_layers"]))]
    names += ["T_journal", "T_bearing"]
    for name in names:
        temperature = fields.get(name)
        if temperature is None or temperature.shape != pressure.shape or not np.isfinite(temperature).all():
            raise ValueError(f"missing or nonfinite resolved film temperature: {name}")
        t, p = temperature[free], pressure[free]
        if np.any((t < 200) | (t > 500)) or np.any(p > 50e6):
            raise ValueError(f"free phase outside R410A EOS range in {name}")
        subcritical = t < 344.494
        t, p = t[subcritical], p[subcritical]
        x = 1-t/344.494
        dew = 4901200*np.exp(344.494/t*(-7.4411*x+1.9883*x**1.6-2.4925*x**2.4-3.2633*x**5))
        if np.any(p > dew*(1+1e-12)):
            raise ValueError(f"free phase outside stable R410A gas domain in {name}")


def validate_caloric(caloric):
    """Require recorded values and provenance, including a nonimplicit latent offset."""
    caloric = dict(caloric)
    if "mu_gas" not in caloric and "gas_viscosity" in caloric:
        caloric["mu_gas"] = caloric["gas_viscosity"]
    for key in CALORIC_KEYS:
        if isinstance(caloric.get(key), bool) or not isinstance(caloric.get(key), (int, float)):
            raise ValueError(f"explicit numeric caloric parameter required: {key}")
        if not math.isfinite(caloric[key]) or (key != "desorption_internal_energy" and caloric[key] <= 0):
            raise ValueError(f"invalid caloric parameter: {key}")
    if not caloric.get("basis") or not caloric.get("sources"):
        raise ValueError("caloric assumptions require basis and sources")
    if not 200 <= caloric["property_reference_temperature"] <= 500:
        raise ValueError("caloric reference outside R410A EOS temperature range")
    return caloric


def check_layer_viscosity(config, fields, curves):
    """Check used inverse viscosities against the unrelaxed current-state law."""
    if config["density_model"] != "LINEAR_TEMPERATURE":
        raise ValueError("layer viscosity check requires LINEAR_TEMPERATURE density")
    shape = fields["p"].shape
    for name in ("film_content", "dissolved_gas"):
        value = fields.get(name)
        if value is None or value.shape != shape or not value.size or not np.isfinite(value).all():
            raise ValueError(f"missing or nonfinite layer constitutive input: {name}")
    fill, concentration = fields["film_content"], fields["dissolved_gas"]
    mu_gas = float(config["mu_gas"])
    if not math.isfinite(mu_gas) or mu_gas <= 0:
        raise ValueError("nonpositive or nonfinite gas viscosity")
    axes = np.array(sorted(curves))
    if np.any((concentration < axes[0]) | (concentration > axes[-1])):
        raise ValueError("concentration outside digitized range")
    upper = np.searchsorted(axes, concentration)
    lower = np.where(axes[upper] == concentration, upper, upper-1)
    span = axes[upper]-axes[lower]
    fraction = np.divide(concentration-axes[lower], span, out=np.zeros_like(concentration), where=span > 0)
    prepared = []
    for axis in axes:
        points = np.array(sorted(curves[axis]))
        prepared.append((points[:, 0], np.log(points[:, 1])))
    maximum = 0.
    for layer in range(int(config["cross_film_layers"])):
        names = (f"cross_film_layer_{layer}", f"cross_film_inv_mu_{layer}")
        for name in names:
            value = fields.get(name)
            if (value is None or value.shape != shape or not value.size or
                    not np.isfinite(value).all() or np.any(value <= 0)):
                raise ValueError(f"missing, nonfinite or nonpositive layer constitutive field: {name}")
        temperature_c = fields[names[0]]-273.15
        density = (float(config["density_temperature_intercept"]) +
                   float(config["density_temperature_slope"])*temperature_c)
        if not np.isfinite(density).all() or np.any(density <= 0):
            raise ValueError("nonpositive or nonfinite layer liquid density")
        log_nu = np.zeros_like(temperature_c)
        for index, (temperatures, ordinates) in enumerate(prepared):
            selected = (lower == index) | (upper == index)
            t = temperature_c[selected]
            if np.any((t < temperatures[0]) | (t > temperatures[-1])):
                raise ValueError(f"temperature outside curve {axes[index]:g} range")
            weight = np.where(lower[selected] == index, 1-fraction[selected], fraction[selected])
            log_nu[selected] += weight*np.interp(t, temperatures, ordinates)
        nu = np.exp(log_nu)*1e-6
        physical = 1/(fill*density*nu+(1-fill)*mu_gas)
        if not np.isfinite(physical).all() or np.any(physical <= 0):
            raise ValueError("nonpositive or nonfinite physical layer inverse viscosity")
        maximum = max(maximum, float(np.max(np.abs(fields[names[1]]-physical)/physical)))
    if maximum >= float(config["coupling_tolerance"]):
        raise ValueError("raw layer viscosity constitutive mismatch exceeds requested tolerance")
    return maximum


def prepare(preflight, destination, caloric):
    """Freeze the selected matched inputs and independent measured-coordinate evidence."""
    caloric = validate_caloric(caloric)
    baseline.prepare(preflight, destination, "LINEAR_TEMPERATURE")
    table = ROOT / "data/R410A_FW68D/fw68d.csv"
    (destination / table.name).write_bytes(table.read_bytes())
    project = ROOT / "data/R410A_FW68D/source_project.json"
    metadata = json.loads(project.read_text())
    decisions = json.loads((destination / "decisions.json").read_text())
    decisions.update(status="conditional_native_coupled_comparison_prepared_unrun",
        concentration="Conditional 0.30 refrigerant mass / liquid-solution mass; source basis unspecified, in-run composition unmeasured.",
        source_chart_bases={key: metadata[key] for key in ("pressureBasis", "concentrationBasis")},
        source_chart_convention="User-authorized conditional absolute pressure and refrigerant mass-fraction interpretation.",
        cavity="Matched zero-Pa oil-vapor-limit surrogate; no invented refrigerant release/absorption rates.",
        matched_comparison={"common": "COUPLED, HOMOGENEOUS, DUKLER, full nu(c,T), rho(T), caloric law, boundary data, mesh, cooling and load balance",
            "gas_off": "SINGLE_PHASE; coupled_phase_equilibrium=false; composition transported without release",
            "gas_on": "TWO_PHASE; coupled_phase_equilibrium=true; instantaneous local equilibrium",
            "incoming_free_gas_mass_fraction": 0.,
            "attribution": "Gas-on minus matched gas-off quantifies release under this conditional model. Existing FBNS STRIATED fixed30 baseline differs structurally and cannot isolate gas effects."},
        caloric_assumptions=caloric,
        numerical_initialization="One gas-free FBNS pressure/film-content solve at unchanged geometry, composition and temperature initializes conservative COUPLED inventories. Only subsequent COUPLED convergence is accepted.",
        limitations="Conditional mass/absolute chart bases; unmeasured mixture density and caloric proxies; no calibrated pressurized parasitic drag; condition stability does not prove friction/composition steady state. Finite gas EOS and chart domains are hard limits.")
    decisions["thermal_analog_basis"] = "Pure-oil rho(T) proxy, explicit caloric law from caloric_assumptions, k_liquid=.13 and steel k=15 proxies; no measured FW68D mixture caloric data."
    experiment = json.loads((destination / "experiment.json").read_text())
    curves, saturation = property_curves(table)
    differences, domains = [], []
    old_points = experiment["viscosity"]["points"]
    full_points = [{"temperature_c": t, "nu_mm2_s": nu} for t, nu in sorted(curves[.3])]
    experiment["viscosity"]["points"] = full_points
    for case in experiment["cases"]:
        key, measured = case["case_id"], case["measured"]
        for sensor in ("tc0", "tc1", "tc2", "tc3"):
            temperature = float(measured[f"{sensor}_c_mean"])
            old = driver.viscosity(temperature, old_points)*1e6
            full = curve_value(curves, .3, temperature, True)
            differences.append(dict(case_id=key, temperature_role=sensor.upper(), temperature_c=temperature,
                old_fixed30_nu_mm2_s=old, full_c03_nu_mm2_s=full, difference_mm2_s=full-old,
                difference_percent=100*(full/old-1)))
        temperature = float(measured["tc0_c_mean"])
        pressure = float(measured["pressure_bar_abs_mean"])*1e5
        try:
            csat = curve_value(saturation, pressure, temperature)
            table_status, table_error = "within_table", ""
        except ValueError as error:
            csat, table_status, table_error = None, "outside_table", str(error)
        dew = dew_pressure(temperature+273.15)
        domains.append(dict(case_id=key, feed_temperature_c=temperature, feed_pressure_Pa_abs=pressure,
            feed_c_sat=csat, feed_table_status=table_status, feed_table_error=table_error,
            feed_dew_pressure_Pa_abs=dew, feed_p_over_dew=pressure/dew,
            hypothetical_feed_gas_EOS_valid=pressure <= min(dew*(1+1e-12), 50e6),
            interpretation="No imposed free gas at feed. EOS validity applies only where a free phase exists; local states must be checked during execution."))
        case["nu_mm2_s"] = curve_value(curves, .3, temperature, True)
        case["mu_pa_s"] = driver.density(temperature)*case["nu_mm2_s"]*1e-6
        original_setup = json.loads((destination / f"{key}_setup.json").read_text())
        for branch in BRANCHES:
            setup = copy.deepcopy(original_setup)
            setup.update(status="conditional_native_coupled_prepared_unrun", gas_branch=branch,
                         scenario_decisions=decisions)
            setup["caloric_assumptions"] = caloric
            setup["material_settings"]["pressure_reference"] = (
                "Measured bar-absolute mean times 100000 at both ends and groove; common zero-Pa "
                "oil-vapor-limit surrogate. " + ("Refrigerant release uses the conditional chart bubble pressure."
                if branch == "gas_on" else "Refrigerant release is disabled."))
            setup["force_balance"]["restart_policy"] = (
                "Fresh conservative species inventories for each cooling call; retain the last native "
                "equilibrium position as its seed. A fresh FBNS pressure/film-content initial guess "
                "is converted to conservative COUPLED inventories; legacy checkpoints are incompatible.")
            setup["acceptance"].update(
                pressure_and_temperature_coupling="Native COUPLED pressure/species residual below configured coupled_nonlinear_tolerance; volume, R_p, R_T, R_g, profile and unrelaxed viscosity residuals below configured coupling_tolerance; no exhausted cap or restored-best iterate.",
                fields="Finite resolved fields, positive gap and viscosity, nonnegative absolute pressure, admissible phase fractions, supported viscosity/bubble-pressure chart and stable free-gas EOS domain at all layers and walls.",
                external_convection="Both exterior coefficients close to relative 1e-3 and wall-temperature change 0.01 K normally, tightened to 1e-5 and 1e-4 K for focused checks. Each call starts fresh species inventories at the preceding native equilibrium position.",
                native_pressure="User-selected common COUPLED pressure/species tolerance 1e-8, applied to both the nonlinear solve and independent acceptance. Focused checks tighten load, thermal and exterior-cooling closure.",
                oil_and_species_balance_relative_residual_max="min(coupling_tolerance, 1e-8), using time-step-independent boundary/feed mass-flux scales",
                global_energy_relative_residual_max="configured coupling_tolerance",
                cross_film_raw_collocation_residual_max="configured coupling_tolerance",
                cross_film_raw_viscosity_relative_mismatch_max="configured coupling_tolerance, independently recomputed from current film content, composition and layer temperature",
                temperature_floor_K="No supply-temperature floor: desorption can cool the film. Enforce actual property domains.",
                property_domain="Strict local nu(c,T), chart-supported bubble pressure for gas-on, and free-gas EOS support at mean, layer and wall temperatures.")
            setup["acceptance"].pop("fbns", None)
            driver.save_json(destination / f"{key}_{branch}_setup.json", setup)
    driver.save_json(destination / "experiment.json", experiment)
    driver.save_json(destination / "decisions.json", decisions)
    write_rows(destination / "fixed30_vs_full_c03.csv", differences)
    write_rows(destination / "feed_domain_review.csv", domains)
    for filename, count in (("pure_coordinates.csv", 103), ("existing_cfd_coordinates.csv", 6)):
        source = preflight / filename
        if len(load_rows(source)) != count:
            raise ValueError(f"wrong preserved coordinate count: {filename}")
        (destination / filename).write_bytes(source.read_bytes())
    driver.save_json(destination / "property_review.json", {
        "status": "preflight_only_local_states_unrun", "source_table_sha256": driver.sha256(table),
        "source_project_sha256": driver.sha256(project), "source_chart_bases": decisions["source_chart_bases"],
        "maximum_fixed30_difference_percent": max(abs(r["difference_percent"]) for r in differences),
        "feed_table_failures": [r["case_id"] for r in domains if r["feed_table_status"] != "within_table"],
        "gas_EOS": "200--500 K; <=50 MPa; pressure <= dew below344.494 K; applies to free phase only",
        "table": "Unclamped nu(c,T) and strict chart-supported bubble-pressure inversion of c_sat(p,T), including all125 digitized points. A gas-free local pressure may exceed chart pressure support while its bubble point remains strictly supported; free phase must lie at that bubble pressure. Contributing curve temperature limits apply.",
        "viscosity_C03_temperature_range_C": [min(t for t, _ in curves[.3]), max(t for t, _ in curves[.3])],
        "matched_native_configurations": 2*len(experiment["cases"]), "preserved_coordinates": {"pure_oil": 103, "CFD": 6}})
    driver.save_json(destination / "manifest.json", {p.name: driver.sha256(p)
        for p in sorted(destination.iterdir()) if p.name != "manifest.json"})


def integrate_pressure(config, fields):
    """Integrate raw pressure over the actual bearing surface across all ranks."""
    nt, nz = int(config["n_theta_global"]), int(config["n_z_global"])
    theta, z = fields["theta_rad"], fields["z_m"]
    if len(theta) != nt*nz or len(z) != nt*nz:
        raise ValueError("incomplete raw pressure mesh")
    dt, dz = 2*math.pi/nt, float(config["L"])/nz
    order = np.lexsort((theta, z))
    tt, zz = theta[order].reshape(nz, nt), z[order].reshape(nz, nt)
    if not (np.allclose(tt, (np.arange(nt)+.5)*dt, rtol=0, atol=1e-12) and
            np.allclose(zz, (np.arange(nz)[:, None]+.5)*dz, rtol=0, atol=1e-12)):
        raise ValueError("duplicate or irregular raw pressure mesh")
    h = fields["h"][order].reshape(nz, nt)
    pressure = fields["p"][order].reshape(nz, nt)-float(config["p_cav"])
    ht = (np.roll(h, -1, axis=1)-np.roll(h, 1, axis=1))/(2*dt)
    hz = np.gradient(h, dz, axis=0) if nz > 1 else np.zeros_like(h)
    radius = float(config["R"])+h
    return np.array([np.sum(pressure*(radius*np.cos(tt)+ht*np.sin(tt))),
                     np.sum(pressure*(radius*np.sin(tt)-ht*np.cos(tt))),
                     np.sum(-pressure*radius*hz)])*dt*dz


def run_artifacts(run_dir):
    """Hash retained outputs, checkpoints and logs, including failed attempts."""
    return {p.relative_to(run_dir).as_posix(): driver.sha256(p)
            for p in sorted(run_dir.rglob("*")) if p.is_file() and p.name != "metrics.json"}


class Campaign(driver.Campaign):
    """Retain established cooling closure and restart checks, with native COUPLED gates."""

    def __init__(self, args):
        super().__init__(args)
        self.branch = self.setup["gas_branch"]
        if self.branch not in BRANCHES:
            raise ValueError("unknown matched gas branch")
        self.caloric = validate_caloric(self.setup["caloric_assumptions"])
        self.driver_hash = driver.sha256(__file__)
        manifest = self.experiment_path.parent / "manifest.json"
        for name, digest in json.loads(manifest.read_text()).items():
            if Path(name).name != name or driver.sha256(manifest.parent / name) != digest:
                raise ValueError(f"frozen input manifest mismatch: {name}")
        self.input_hashes.update(inputs_manifest_sha256=driver.sha256(manifest),
            base_driver_sha256=driver.sha256(driver.__file__),
            baseline_worker_sha256=driver.sha256(baseline.__file__))

    def recover_checkpoint(self, case):
        # Existing checkpoints omit conservative species inventories.
        return None

    def command(self):
        command = super().command()
        return command + ["-coupled_sub_pc_type", "lu", "-coupled_snes_monitor",
                          "-coupled_ksp_converged_reason"]

    def config_for(self, case, eps, psi, h_bush, h_shaft, output_dir, initial_state_file=None):
        if initial_state_file is not None:
            raise ValueError("COUPLED needs a complete inventory restart; legacy FBNS checkpoints are incompatible")
        values = parse_config(super().config_for(case, eps, psi, h_bush, h_shaft, output_dir, initial_state_file))
        values.update(cavitation_model="COUPLED", cavitated_film_model="HOMOGENEOUS",
            steady_pressure_inner_only="false",
            coupled_fbns_initialization="true",
            coupled_nonlinear_tolerance=1e-8,
            coupling_relaxation=.3,
            steady_position_fd_step=.0001,
            groove_mixing_model="PINNED",
            fluid_property_model="TWO_PHASE" if self.branch == "gas_on" else "SINGLE_PHASE",
            coupled_phase_equilibrium="true" if self.branch == "gas_on" else "false",
            dissolved_gas_species="R410A", dissolved_gas_initial=.3, dissolved_gas_max=1.,
            solubility_model="CURVE_TABLE", gas_mass_transfer_rate=0., gas_absorption_rate=0.,
            dissolved_gas_diffusivity=0., mixture_viscosity_model="DUKLER", gas_alpha_max=1.,
            bc_z_south_dissolved_gas=.3, bc_z_north_dissolved_gas=.3,
            bc_z_south_free_gas_mass_fraction=0., bc_z_north_free_gas_mass_fraction=0.,
            inlet_species="0 0.3 0", solution_table_file=str(self.table_path),
            gas_pressure_floor=1., **{key: self.caloric[key] for key in CALORIC_KEYS})
        values["output_fields"] += ",pressure_force_z,dissolved_gas,alpha_gas,gas_mass_transfer," + ",".join(
            name for layer in range(int(values["cross_film_layers"]))
            for name in (f"cross_film_layer_{layer}", f"cross_film_inv_mu_{layer}"))
        return "# Matched conditional native COUPLED branch: " + self.branch + "\n" + "\n".join(
            f"{key} = {value}" for key, value in values.items()) + "\n"

    def check_solution(self, run_dir, case, returncode=0):
        """Require native pressure, species, thermal and independently integrated load gates."""
        if returncode:
            raise ValueError(f"solver exit code {returncode}")
        config = parse_config((run_dir / "config.txt").read_text())
        summary_rows = load_rows(run_dir / "results/coupled_summary.csv")
        if len(summary_rows) != 1:
            raise ValueError("expected one final native COUPLED summary")
        summary = {key: float(value) for key, value in summary_rows[0].items()}
        if not all(math.isfinite(value) for value in summary.values()):
            raise ValueError("nonfinite native COUPLED summary")
        if any(summary[key] != value for key, value in
               (("converged", 1), ("pressure_converged", 1), ("cap_exhausted", 0), ("restored_best", 0))):
            raise ValueError("native COUPLED failed convergence, cap or restored-iterate gate")
        tolerance = float(config["coupling_tolerance"])
        pressure_tolerance = float(config.get("coupled_nonlinear_tolerance", "1e-10"))
        if not 0 < pressure_tolerance < 1 or not 0 <= summary["pressure_residual"] < pressure_tolerance:
            raise ValueError("native COUPLED pressure residual exceeds configured tolerance")
        for key in ("volume_error", "R_p", "R_T", "R_g", "profile_residual", "viscosity_residual"):
            if key not in summary or not 0 <= summary[key] < tolerance:
                raise ValueError(f"native COUPLED {key} fails requested tolerance")
        for component in ("oil", "species"):
            scale, residual = summary[f"{component}_balance_scale_kg_s"], summary[f"{component}_balance_kg_s"]
            relative = abs(residual)/max(scale, 1e-30)
            if scale < 0 or not math.isclose(relative, summary[f"{component}_balance_rel"], rel_tol=1e-8, abs_tol=1e-14):
                raise ValueError(f"inconsistent {component} balance scale")
            if relative >= min(tolerance, 1e-8):
                raise ValueError(f"{component} dt-independent balance exceeds acceptance")
        energy = abs(summary["energy_residual_W"])/max(summary["energy_scale_W"], 1.)
        if summary["energy_scale_W"] < 0 or energy >= tolerance:
            raise ValueError("native COUPLED energy balance exceeds requested tolerance")
        fields = driver.read_fields(run_dir)
        required = ("p", "film_content", "h", "T", "T_bearing", "T_journal", "mu", "q_bearing",
                    "q_journal", "friction_torque", "dissolved_gas", "alpha_gas", "temperature_housing_outer")
        for key in required:
            if key not in fields or not len(fields[key]) or not np.isfinite(fields[key]).all():
                raise ValueError(f"missing or nonfinite {key} field")
        if np.min(fields["h"]) <= 0 or np.min(fields["mu"]) <= 0:
            raise ValueError("nonpositive gap or viscosity")
        if np.min(fields["p"]) < -pressure_roundoff(config):
            raise ValueError("negative absolute pressure beyond converged round-off")
        for key in ("film_content", "alpha_gas", "dissolved_gas"):
            if np.min(fields[key]) < -1e-8 or np.max(fields[key]) > 1+1e-8:
                raise ValueError(f"{key} outside [0,1]")
        if self.branch == "gas_off" and np.max(fields["alpha_gas"]) > 1e-10:
            raise ValueError("matched gas-off contains released free gas")
        check_layer_eos(config, fields)
        # Evaporation may cool below supply temperature; require actual property domains.
        curves, saturation = property_curves(self.table_path)
        viscosity_mismatch = check_layer_viscosity(config, fields, curves)
        for concentration, temperature, pressure, alpha in zip(
                fields["dissolved_gas"], fields["T"]-273.15, fields["p"], fields["alpha_gas"]):
            curve_value(curves, float(concentration), float(temperature), True)
            if self.branch == "gas_on":
                bubble = bubble_pressure(saturation, float(concentration), float(temperature))
                if pressure < bubble-1. or (alpha > 1e-6 and abs(pressure-bubble) > 1.):
                    raise ValueError("raw phase fractions violate bubble-pressure equilibrium")
            if alpha > 1e-10 and (not 200 <= temperature+273.15 <= 500 or
                    pressure > min(50e6, dew_pressure(float(temperature)+273.15)*(1+1e-12))):
                raise ValueError("resolved free phase outside R410A stable gas EOS domain")
        if "F" not in fields or not len(fields["F"]) or len(fields["F"]) % 3:
            raise ValueError("missing integrated total fluid force")
        force = fields["F"].reshape(-1, 3)
        if not np.isfinite(force).all() or not np.allclose(force, force[0], rtol=1e-9, atol=1e-6):
            raise ValueError("inconsistent integrated force")
        pressure_force = integrate_pressure(config, fields)
        if "Fp" not in fields or len(fields["Fp"]) != len(fields["p"])*3:
            raise ValueError("missing native pressure force vector")
        if not np.allclose(fields["Fp"].reshape(-1, 3), pressure_force, rtol=1e-8, atol=1e-6):
            raise ValueError("raw pressure integration differs from native pressure force")
        fx, fy, fz = map(float, force[0])
        load = float(case["load_mean_n"])
        force_error = float(np.linalg.norm(driver.force_residual(fx, fy, load)))/load
        if force_error >= float(config["steady_force_tolerance"]):
            raise ValueError("independent VTK load balance fails")
        torque = float(fields["friction_torque"][0])
        if not np.allclose(fields["friction_torque"], torque, rtol=1e-9, atol=1e-12):
            raise ValueError("inconsistent integrated torque")
        diagnostic = driver.last_diagnostic(run_dir)
        profile_min = float(diagnostic["cross_film_profile_t_min_K"])
        profile_max = float(diagnostic["cross_film_profile_t_max_K"])
        if not math.isfinite(profile_min) or not math.isfinite(profile_max) or profile_min > profile_max:
            raise ValueError("invalid cross-film profile temperature envelope")
        return {**driver.final_position(fields, float(config["c"])),
            "Fx_N": fx, "Fy_N": fy, "Fz_N": fz, "torque_Nm": torque,
            "raw_pressure_force_N": pressure_force.tolist(),
            "force_error_rel": force_error, "native_force_error_rel": summary["force_residual"],
            "h_min_m": float(np.min(fields["h"])), "mu_min_Pa_s": float(np.min(fields["mu"])),
            "p_min_Pa": float(np.min(fields["p"])), "p_max_Pa": float(np.max(fields["p"])),
            "T_min_K": float(np.min(fields["T"])), "T_max_K": float(np.max(fields["T"])),
            "T_bore_min_K": float(np.min(fields["T_bearing"])), "T_bore_mean_K": float(np.mean(fields["T_bearing"])),
            "T_bore_max_K": float(np.max(fields["T_bearing"])),
            "T_outer_mean_K": float(np.mean(fields["temperature_housing_outer"])),
            "T_shaft_mean_K": float(np.mean(fields["T_journal"])),
            "energy_balance_residual_rel": energy, "cross_film_profile_residual": summary["profile_residual"],
            "cross_film_viscosity_residual": viscosity_mismatch,
            "cross_film_profile_t_min_K": profile_min, "cross_film_profile_t_max_K": profile_max,
            "oil_balance_rel": summary["oil_balance_rel"], "species_balance_rel": summary["species_balance_rel"],
            "alpha_gas_max": float(np.max(fields["alpha_gas"])), "alpha_gas_mean": float(np.mean(fields["alpha_gas"])),
            "dissolved_gas_min": float(np.min(fields["dissolved_gas"])),
            "dissolved_gas_max": float(np.max(fields["dissolved_gas"])),
            "pressure_acceptance": "native_COUPLED_pressure_species_volume_energy_force_profile_verified",
            "coupled_summary": summary, "diagnostics": diagnostic}

    def recheck(self, run_dir, record, case):
        if not self.matching_provenance(record):
            return None
        try:
            path = run_dir / "config.txt"
            if record.get("config_sha256") != driver.sha256(path) or not driver.config_has_current_load(path.read_text(), case["load_mean_n"]):
                return None
            if not record.get("artifact_sha256") or record["artifact_sha256"] != run_artifacts(run_dir):
                return None
            return {**record, **self.check_solution(run_dir, case), "run_dir": str(run_dir)}
        except (OSError, ValueError, KeyError, ET.ParseError):
            return None

    def single(self, case, eps, psi, h_bush, h_shaft, label):
        cached = self.find_completed(case, eps, psi, h_bush, h_shaft)
        if cached is not None:
            return cached
        self.sequence += 1
        run_dir = self.work / case["case_id"] / f"{self.sequence:04d}_{label}_{uuid.uuid4().hex[:8]}"
        run_dir.mkdir(parents=True, exist_ok=False)
        path = run_dir / "config.txt"
        path.write_text(self.config_for(case, eps, psi, h_bush, h_shaft, run_dir / "results"), encoding="utf-8")
        command = self.command()
        record = dict(case_id=case["case_id"], gas_branch=self.branch, eps=eps, psi_deg=psi,
            seed_eps=eps, seed_psi_deg=psi, native_equilibrium=True,
            h_bush_W_m2_K=h_bush, h_shaft_W_m2_K=h_shaft, run_dir=str(run_dir), command=command,
            mpi_ranks=self.args.ranks, python_version=sys.version.split()[0], mpi_version=self.mpi_version,
            source_sha256=self.source_hash, driver_sha256=self.driver_hash, rankfile_sha256=self.rankfile_hash,
            config_sha256=driver.sha256(path), executable_sha256=self.executable_hash,
            initial_state_file=None, initial_state_sha256=None,
            restart_note="Fresh species inventories: legacy checkpoint omits conservative species state.",
            status="started", runtime_timeout_seconds=self.args.timeout, **self.input_hashes)
        driver.save_json(run_dir / "metrics.json", record)
        started = time.monotonic()
        try:
            code, deadline = driver.run_solver(command, run_dir,
                dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"), self.args.timeout)
            record.update(returncode=code, runtime_deadline_reached=deadline, wall_seconds=time.monotonic()-started)
            if deadline:
                raise RuntimeError("runtime deadline reached; failed outputs preserved")
            record.update(self.check_solution(run_dir, case, code), status="valid_inner_solve")
        except Exception as error:
            record["artifact_sha256"] = run_artifacts(run_dir)
            record.update(status="failed", error=str(error), wall_seconds=time.monotonic()-started)
            driver.save_json(run_dir / "metrics.json", record)
            raise RuntimeError(f"failed native COUPLED solve preserved at {run_dir}: {error}") from error
        driver.save_json(run_dir / "metrics.json", record)
        record["artifact_sha256"] = run_artifacts(run_dir)
        driver.save_json(run_dir / "metrics.json", record)
        self.completed[case["case_id"]].append((run_dir, record))
        return record


def check_inputs(inputs, output):
    """Verify matching controls and preserved coordinates without solving the bearing."""
    verified = []
    experiment = json.loads((inputs / "experiment.json").read_text())
    for key in (case["case_id"] for case in experiment["cases"]):
        values = {}
        for branch in BRANCHES:
            options = baseline.options(inputs, output.parent / "unused", None, None)
            options.table = str(inputs / "fw68d.csv")
            options.materials = str(inputs / f"{key}_{branch}_setup.json")
            campaign = Campaign(options)
            case = campaign.cases[key]
            eps, psi = campaign.setup["force_balance"]["initial_positions"][key]
            values[branch] = parse_config(campaign.config_for(case, eps, psi, 60., 100., output.parent / "unused"))
            config = values[branch]
            pressure = float(case["measured"]["pressure_bar_abs_mean"])*1e5
            assert config["cavitation_model"] == "COUPLED"
            assert config["cavitated_film_model"] == "HOMOGENEOUS"
            assert config["steady_pressure_inner_only"] == "false"
            assert config["coupled_fbns_initialization"] == "true"
            assert float(config["coupled_nonlinear_tolerance"]) == 1e-8
            assert float(config["coupling_relaxation"]) == .3
            assert float(config["steady_position_fd_step"]) == .0001
            assert config["density_model"] == "LINEAR_TEMPERATURE"
            assert float(config["density_temperature_intercept"]) == 954.9
            assert float(config["density_temperature_slope"]) == -.7284
            output_fields = config["output_fields"].split(",")
            assert all(name in output_fields for layer in range(int(config["cross_film_layers"]))
                       for name in (f"cross_film_layer_{layer}", f"cross_film_inv_mu_{layer}"))
            assert config["solution_table_clamp"] == "false"
            assert config["solubility_model"] == "CURVE_TABLE"
            assert float(config["dissolved_gas_initial"]) == .3
            assert float(config["gas_mass_transfer_rate"]) == float(config["gas_absorption_rate"]) == 0.
            assert float(config["bc_z_south_val"]) == float(config["bc_z_north_val"]) == pressure
            assert float(config["inlet_groove"].split()[2]) == pressure
            assert config["inlet_species"] == "0 0.3 0"
            assert all(float(config[f"bc_z_{side}_free_gas_mass_fraction"]) == 0 for side in ("north", "south"))
            assert math.isclose(case["S_density_reference_TC2"], float(case["measured"]["S"])*
                                driver.density(float(case["measured"]["tc2_c_mean"]))/950., rel_tol=1e-13)
        different = [name for name in values["gas_off"] if values["gas_off"][name] != values["gas_on"][name]]
        assert set(different) == {"fluid_property_model", "coupled_phase_equilibrium"}, different
        verified.append({"case_id": key, "differing_controls": different})
    driver.save_json(output, {"status": f"{2*len(verified)}_matched_native_configs_verified_unrun", "cases": verified})


def run_case(args):
    """Run one branch/case with independent closure verification and exclusive outputs."""
    destination = args.work / args.branch / args.case
    destination.mkdir(parents=True, exist_ok=False)
    setup = json.loads((args.inputs / f"{args.case}_{args.branch}_setup.json").read_text())
    setup["material_settings"]["cross_film_layers"] = args.layers
    if args.tight:
        setup["solver_tolerances"] = dict(steady_force_tolerance=2e-5, coupling_tolerance=1e-8, fbns_tol=1e-12)
        setup["cooling_relative_tolerance"] = 1e-5
        setup["cooling_temperature_tolerance_K"] = 1e-4
    driver.save_json(destination / "setup.json", setup)
    options = baseline.options(args.inputs, destination, args.solver, args.rankfile)
    options.table = str(args.inputs / "fw68d.csv")
    options.ranks = args.ranks
    options.timeout = args.timeout or None
    campaign = Campaign(options)
    campaign.mpi_version = subprocess.check_output([campaign.command()[0], "--version"], text=True).splitlines()[:2]
    case = campaign.cases[args.case]
    seed = setup["force_balance"]["initial_positions"][args.case]
    started = time.monotonic()
    try:
        result = campaign.balance(case, *seed)
        run = Path(result["run_dir"])
        checked = campaign.recheck(run, result, case)
        closed = campaign.find_completed(case, result["seed_eps"], result["seed_psi_deg"], closed=True)
        if checked is None or closed is None or Path(closed["run_dir"]) != run:
            raise ValueError("failed independent native and cooling closure recheck")
        config = parse_config((run / "config.txt").read_text())
        fields = driver.read_fields(run)
        if float(config["omega"])*result["torque_Nm"] <= 0:
            raise ValueError("integrated journal torque has inconsistent rotation sign")
        observables = baseline.solid_observables(config, fields, .026)
        np.savez_compressed(destination / "fields.npz", **fields)
        driver.save_json(destination / "accepted.json", dict(case=case["measured"], result=result,
            config=config, observables=observables, elapsed_seconds=time.monotonic()-started,
            worker_sha256=driver.sha256(__file__), fields_sha256=driver.sha256(destination / "fields.npz"),
            tight=args.tight, gas_branch=args.branch, cross_film_layers=args.layers,
            status="accepted_conditional_native_coupled"))
    except Exception as error:
        driver.save_json(destination / "failure.json", dict(case=args.case, gas_branch=args.branch,
            error=str(error), elapsed_seconds=time.monotonic()-started, status="failed_preserved"))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", type=Path)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--caloric", type=Path)
    parser.add_argument("--check-inputs", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--branch", choices=BRANCHES)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--rankfile", type=Path)
    parser.add_argument("--ranks", type=int, choices=(1, 2, 4), default=4)
    parser.add_argument("--layers", type=int, choices=(16, 32), default=16)
    parser.add_argument("--tight", action="store_true")
    parser.add_argument("--timeout", type=int, default=0,
                        help="seconds per native solve; 0 (default) disables the wall-clock deadline")
    args = parser.parse_args()
    if args.timeout < 0:
        parser.error("--timeout must be nonnegative")
    if args.prepare and args.caloric:
        prepare(args.prepare, args.inputs, json.loads(args.caloric.read_text()))
    elif args.check_inputs:
        if args.check_inputs.exists():
            raise FileExistsError(args.check_inputs)
        check_inputs(args.inputs, args.check_inputs)
    elif args.case and args.branch and args.work and args.solver:
        cases = json.loads((args.inputs / "experiment.json").read_text())["cases"]
        if args.case not in [case["case_id"] for case in cases]:
            parser.error("--case must belong to the frozen experiment selection")
        run_case(args)
    else:
        parser.error("use --prepare PREFLIGHT --caloric JSON, or --check-inputs REPORT")


if __name__ == "__main__":
    main()
