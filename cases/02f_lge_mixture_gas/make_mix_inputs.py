"""Build 02f mixture inputs (gas_off and gas_on) for X1-X6 from selected_points.csv.

Usage:
  python3 make_mix_inputs.py <source_clone> <stribeck_mix_inputs> <selected_points.csv> <out_dir>
                             [--composition C] [--seeds seeds.json]

Templates (physics, caloric, tolerances) come from the Stribeck mix_thermal inputs (T3/P09 derived).
Per case: rpm, load, TC0 (supply/ambient) and chamber pressure from the nearest 30 s bin.
--composition sets the dissolved refrigerant mass fraction c (default 0.30, the
lge-thermal-v2-20261029 value). It sets the reference viscosity nu(c, TC0), the
experiment viscosity points, and the setup key read by run_gas.Campaign for the
initial field, both axial ends and the groove feed.
--seeds is a JSON map {branch: {case: [eps, psi_deg]}} of numerical seeds; a case
missing there gets the S-interpolated seed of the 0.30 campaign. Seeds affect
convergence only.
Outputs are written with LF line endings so the manifest is platform independent.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("source", type=Path)
parser.add_argument("template", type=Path)
parser.add_argument("points", type=Path)
parser.add_argument("out", type=Path)
parser.add_argument("--composition", type=float, default=None)
parser.add_argument("--seeds", type=Path)
args = parser.parse_args()
source, template, points, out = (p.resolve() for p in (args.source, args.template, args.points, args.out))
sys.path.insert(0, str(source / "cases/02f_lge_mixture_gas"))
import run_gas

composition = run_gas.validate_composition(
    run_gas.DEFAULT_REFRIGERANT_MASS_FRACTION if args.composition is None else args.composition)
default = composition == run_gas.DEFAULT_REFRIGERANT_MASS_FRACTION
seeds = json.loads(args.seeds.read_text()) if args.seeds else {}


def write_json(path, value):
    path.write_bytes((json.dumps(value, indent=2)).encode())


out.mkdir(parents=True, exist_ok=False)
for name in ("fw68d.csv", "property_review.json", "decisions.json"):
    shutil.copy2(template / name, out / name)
viscosity, _ = run_gas.property_curves(out / "fw68d.csv")
experiment = json.loads((template / "experiment.json").read_text())
experiment["viscosity"]["points"] = run_gas.composition_viscosity_points(viscosity, composition)
base = experiment["cases"][0]
setups = {b: json.loads((template / f"X01_{b}_setup.json").read_text()) for b in run_gas.BRANCHES}
rho = lambda t: 954.9 - 0.7284 * t

if not default:
    label = run_gas.composition_label(composition)
    decisions = json.loads((out / "decisions.json").read_text())
    decisions["concentration"] = (f"Measured {label} refrigerant mass / liquid-solution mass (weighing method, "
                                  "24 h at 14.4 bar); source chart basis unspecified, in-run composition unmeasured.")
    decisions["refrigerant_mass_fraction"] = composition
    write_json(out / "decisions.json", decisions)
    review = json.loads((out / "property_review.json").read_text())
    review["refrigerant_mass_fraction"] = composition
    review["viscosity_composition_temperature_range_C"] = [experiment["viscosity"]["points"][0]["temperature_c"],
                                                           experiment["viscosity"]["points"][-1]["temperature_c"]]
    review["viscosity_composition_note"] = (f"nu(c={label}, T): ln(nu) interpolated linearly in c between the "
                                            "bracketing digitized curves, unclamped, on their common temperature range.")
    write_json(out / "property_review.json", review)

cases = []
for row in csv.DictReader(points.open()):
    if row["fluid"] != "mix":
        continue
    key, rpm, load = row["case"], float(row["rpm"]), float(row["W"])
    t_c, p_pa, s = float(row["TC0"]), float(row["P_bar"]) * 1e5, float(row["S"])
    nu = run_gas.curve_value(viscosity, composition, t_c, True)
    case = copy.deepcopy(base)
    case.update(case_id=key, rpm_nominal=rpm, load_mean_n=load,
                angular_speed_magnitude_rad_s=rpm * math.pi / 30, angular_speed_rad_s=rpm * math.pi / 30,
                tc0_boundary_temperature_c=t_c, tc0_boundary_temperature_k=t_c + 273.15,
                rho_kg_m3=rho(t_c), nu_mm2_s=nu, mu_pa_s=rho(t_c) * nu * 1e-6, sommerfeld_number=s,
                apparent_bearing_torque_nm=None, friction_coefficient=float(row["f"]),
                measured=dict(row, pressure_bar_abs_mean=row["P_bar"], tc0_c_mean=row["TC0"],
                              tc2_c_mean=row["TC2"], bearing_torque_nm=row["bt_meas"]),
                S_density_reference_TC2=s * rho(float(row["TC2"])) / 950.,
                note="Thermal v2 mixture case; Oct 27 mixture rerun point, conditions from nearest 30 s bin")
    cases.append(case)
    eps = float(np.interp(math.log(s), [math.log(v) for v in (.004, .03, .1, 1.)], [.94, .85, .60, .15]))
    default_seed = [eps, 300. if eps > .8 else 320.]
    for branch, setup in setups.items():
        seed = [float(v) for v in seeds.get(branch, {}).get(key, default_seed)]
        setup = copy.deepcopy(setup)
        setup["material_settings"]["reservoir_pressure_Pa"] = p_pa
        setup["force_balance"]["initial_positions"] = {key: seed}
        setup["case_ids"] = [key]
        setup["additional_cases"] = []
        if not default:
            setup["refrigerant_mass_fraction"] = composition
            setup["scenario_decisions"] = json.loads((out / "decisions.json").read_text())
        assert setup["gas_branch"] == branch
        write_json(out / f"{key}_{branch}_setup.json", setup)
        print(f"  {key} {branch} seed eps={seed[0]:.6f} psi={seed[1]:.4f}")
    print(f"{key} rpm={rpm:.0f} W={load:.1f} TC0={t_c:.2f} P={p_pa / 1e5:.3f} bar c={composition:g} "
          f"nu(TC0)={nu:.3f} mm2/s S(TC2)={s:.4f}")
experiment["cases"] = cases
write_json(out / "experiment.json", experiment)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
write_json(out / "manifest.json", {p.name: sha(p) for p in sorted(out.iterdir())})
