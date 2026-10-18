"""Check the fixed-chart mixture driver without executing a CFD solve."""

import json
import math
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
driver = runpy.run_path(str(ROOT / "cases/02d_fw68d_steady_six/run_steady.py"))


def rejects(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError("invalid material or speed accepted")


def main():
    with TemporaryDirectory() as directory:
        work = Path(directory)
        experiment = json.loads(driver["EXPERIMENT"].read_text())
        setup = json.loads(driver["SETUP"].read_text())
        case = dict(experiment["cases"][0], case_id="M01", tc0_boundary_temperature_c=30.,
                    tc0_boundary_temperature_k=303.15, angular_speed_rad_s=-20.)
        case.pop("angular_speed_magnitude_rad_s")
        experiment["cases"] = [case]
        points = [{"temperature_c": 30., "nu_mm2_s": 8.2435},
                  {"temperature_c": 40., "nu_mm2_s": 6.6992}]
        experiment["viscosity"]["points"] = points
        setup["additional_cases"] = []
        setup["force_balance"]["initial_positions"] = {"M01": [.3, 270.]}
        materials = setup["material_settings"]
        materials.update(oil_density_kg_m3=950., oil_thermal_expansion_1_K=.7284 / 950,
                         reservoir_pressure_Pa=1500000.)
        experiment_path, setup_path, table_path = (work / name for name in
                                                  ("experiment.json", "setup.json", "fixed.csv"))
        experiment_path.write_text(json.dumps(experiment))
        setup_path.write_text(json.dumps(setup))
        table_path.write_text("c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n"
                              "0,30,,8.2435\n0,40,,6.6992\n")
        args = SimpleNamespace(experiment=experiment_path, table=table_path, materials=setup_path,
                               solver=None, work=str(work), ranks=1, coupling_max_iters=2000,
                               coupling_relaxation=1., trial=None, rankfile=None,
                               max_thermal_iterations=30)
        campaign = driver["Campaign"](args)
        config = campaign.config_for(case, .3, 270., 60., 100., work / "results")
        values = dict(line.split(" = ", 1) for line in config.splitlines() if " = " in line)
        assert values["solution_table_file"] == str(table_path.resolve())
        assert campaign.input_hashes["experiment_sha256"] == driver["sha256"](experiment_path)
        assert campaign.input_hashes["property_table_sha256"] == driver["sha256"](table_path)
        assert campaign.input_hashes["material_setup_sha256"] == driver["sha256"](setup_path)
        assert values["density_model"] == "CONSTANT" and float(values["rho"]) == 950
        assert float(values["rho_cp"]) == 950 * 2010
        assert math.isclose(float(values["density_temperature_slope"]), -.7284)
        assert math.isclose(float(values["mu"]), 950 * 8.2435e-6)
        assert values["solution_table_clamp"] == "false"
        assert values["kinematic_viscosity_log_interpolation"] == "true"
        assert values["cross_film_layers"] == "16"
        refined = driver["render_config"](case, experiment, campaign.boundaries,
                dict(materials, cross_film_layers=32), .3, 270., 60., 100., work / "results", table_path)
        assert refined == config.replace("cross_film_layers = 16\n", "cross_film_layers = 32\n")
        for value in (True, False, 16., 32., "32", 0, 8, 64, None):
            invalid = dict(materials, cross_film_layers=value)
            rejects(lambda: driver["render_config"](case, experiment, campaign.boundaries,
                    invalid, .3, 270., 60., 100., work, table_path))
            setup["material_settings"] = invalid
            setup_path.write_text(json.dumps(setup))
            rejects(lambda: driver["Campaign"](args))
        setup["material_settings"] = materials
        setup_path.write_text(json.dumps(setup))
        assert float(values["omega"]) == -20
        assert float(values["bc_z_south_val"]) == float(values["bc_z_north_val"]) == 1500000
        assert float(values["inlet_groove"].split()[2]) == 1500000
        assert float(values["p_cav"]) == 0
        assert campaign.cases["M01"]["sommerfeld_number"] == case["sommerfeld_number"]
        assert driver["density"](40.) == 925.764
        assert driver["density"](40., materials) == driver["density"](80., materials) == 950
        assert driver["angular_speed"]({"angular_speed_magnitude_rad_s": 20.}) == 20
        h_plus, _ = driver["shaft_coefficient"](310., 303.15, 20., points, materials)
        h_minus, _ = driver["shaft_coefficient"](310., 303.15, -20., points, materials)
        assert h_plus == h_minus
        _, metadata = driver["housing_coefficient"](310., 303.15, points, materials)
        nu, rho, _, k = driver["oil_properties"]((310. + 303.15) / 2 - 273.15, points, materials)
        expected_ra = 9.80665 * materials["oil_thermal_expansion_1_K"] * (310. - 303.15) * .085**3
        expected_ra /= nu * k / (rho * materials["oil_cp_J_kg_K"])
        assert math.isclose(metadata["Ra_D"], expected_ra)
        fixture = runpy.run_path(str(ROOT / "tests/test_fw68d_force_balance.py"))
        run = work / "fixture"
        log = fixture["write_solution"](run, case, 4.25e-5)
        result = driver["check_solution"](run, log, 0, case["load_mean_n"], 0., 303.15,
                                           points, True, 4.25e-5, materials=materials)
        assert result["rho_min_kg_m3"] == 950
        for key in ("oil_density_kg_m3", "oil_thermal_expansion_1_K",
                    "oil_cp_J_kg_K", "oil_conductivity_W_m_K"):
            for value in (0., -1., math.nan, math.inf):
                invalid = dict(materials, **{key: value})
                rejects(lambda: driver["render_config"](case, experiment, campaign.boundaries,
                        invalid, .3, 270., 60., 100., work, table_path))
        for value in (0., math.nan, math.inf):
            rejects(lambda: driver["angular_speed"]({"angular_speed_rad_s": value}))
    print("Mixture driver passed: actual input provenance, constant density, signed rotation, cooling and acceptance")


if __name__ == "__main__":
    main()
