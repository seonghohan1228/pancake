"""Verify conditional matched gas inputs without launching numerical campaigns."""

import csv
import json
import math
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
import unittest

ROOT = Path(__file__).resolve().parents[1]
worker = runpy.run_path(str(ROOT / "cases/02f_lge_mixture_gas/run_gas.py"))


def rows(path, values):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


class MatchedGasInputsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (ROOT / "tmp").mkdir(exist_ok=True)

    @staticmethod
    def caloric():
        return dict(dissolved_heat_capacity=1600., desorption_internal_energy=190000.,
            gas_thermal_conductivity=.013, gas_viscosity=1.3e-5,
            property_reference_temperature=293.15, property_reference_pressure=1e5,
            basis="explicit synthetic test values", sources=["test fixture"])

    def test_strict_interpolation_and_exact_node(self):
        interpolate = worker["curve_value"]
        curves = {0.: [(10., 4.), (20., 1.)], 1.: [(30., 8.), (40., 2.)]}
        self.assertAlmostEqual(interpolate(curves, 0., 15., True), 2.)
        with self.assertRaises(ValueError):
            interpolate(curves, .5, 15., True)
        with self.assertRaises(ValueError):
            interpolate(curves, -1., 15., True)
        with self.assertRaises(ValueError):
            worker["dew_pressure"](199.)

    def test_raw_pressure_integration_and_mesh_rejection(self):
        import numpy as np
        nt, nz, radius, length = 64, 3, .0215, .04
        theta = np.tile((np.arange(nt)+.5)*2*math.pi/nt, nz)
        fields = dict(theta_rad=theta, z_m=np.repeat((np.arange(nz)+.5)*length/nz, nt),
                      h=np.full(nt*nz, 40e-6), p=2e6+1e6*np.cos(theta))
        config = dict(n_theta_global=nt, n_z_global=nz, R=radius, L=length, p_cav=0)
        # Deliberately scramble rank ordering; coordinates determine the integration.
        order = np.arange(nt*nz)[::-1]
        actual = worker["integrate_pressure"](config, {key: value[order] for key, value in fields.items()})
        np.testing.assert_allclose(actual, [1e6*(radius+40e-6)*math.pi*length, 0, 0], atol=1e-9)
        fields["theta_rad"][0] = fields["theta_rad"][1]
        with self.assertRaises(ValueError):
            worker["integrate_pressure"](config, fields)

    def test_artifact_tampering_changes_manifest(self):
        with TemporaryDirectory(dir=ROOT / "tmp") as directory:
            root = Path(directory)
            (root / "run.log").write_text("success")
            first = worker["run_artifacts"](root)
            (root / "metrics.json").write_text("ignored mutable record")
            self.assertEqual(first, worker["run_artifacts"](root))
            (root / "run.log").write_text("tampered")
            self.assertNotEqual(first, worker["run_artifacts"](root))

    def test_resolved_layer_eos_catches_hidden_cold_gas(self):
        import numpy as np
        fields = {name: np.array([320.]) for name in
                  ("cross_film_layer_0", "cross_film_layer_1", "T_journal", "T_bearing")}
        fields.update(p=np.array([2e6]), alpha_gas=np.array([.1]))
        worker["check_layer_eos"]({"cross_film_layers": "2"}, fields)
        fields["cross_film_layer_1"][:] = 270.
        with self.assertRaises(ValueError):
            worker["check_layer_eos"]({"cross_film_layers": "2"}, fields)
        fields["alpha_gas"][:] = 0.
        worker["check_layer_eos"]({"cross_film_layers": "2"}, fields)
        del fields["cross_film_layer_1"]
        with self.assertRaises(ValueError):
            worker["check_layer_eos"]({"cross_film_layers": "2"}, fields)

    def test_caloric_provenance_required(self):
        with self.assertRaises(ValueError):
            worker["validate_caloric"]({})
        caloric = self.caloric()
        self.assertEqual(worker["validate_caloric"](caloric)["mu_gas"], 1.3e-5)

    def test_bubble_inverse_requires_strict_supported_state(self):
        inverse = worker["bubble_pressure"]
        curves = {1e5: [(10., .1), (30., .2)], 2e5: [(10., .3), (30., .4)]}
        self.assertAlmostEqual(inverse(curves, .25, 20.), 1.5e5)
        self.assertEqual(inverse(curves, 0., 100.), 0.)
        with self.assertRaises(ValueError):
            inverse(curves, .5, 20.)
        with self.assertRaises(ValueError):
            inverse(curves, .2, 5.)

    def test_matched_preparation(self):
        with TemporaryDirectory(dir=ROOT / "tmp") as directory:
            root = Path(directory)
            preflight = root / "preflight"
            preflight.mkdir()
            table = ROOT / "data/R410A_FW68D/fw68d.csv"
            curves, _ = worker["property_curves"](table)
            points = [{"temperature_c": t, "nu_mm2_s": nu} for t, nu in sorted(curves[.3])[1:]]
            rows(preflight / "viscosity_nodes.csv", points)
            (preflight / "solution_viscosity_surrogate.csv").write_text(
                "c_mass_fraction,T_C,P_MPa,viscosity_mm2s\n0,20,,10\n0,80,,2\n")
            measured = []
            for index in range(1, 18):
                temperature, rpm, load = 26.+index/10, 200., 1000.+index
                nu = worker["driver"].viscosity(temperature, points)
                measured.append(dict(selected_id=f"M{index:02d}", rpm=rpm, force_n_mean=load,
                    tc0_c_mean=temperature, tc1_c_mean=temperature+.1,
                    tc2_c_mean=temperature+.2, tc3_c_mean=temperature+.3,
                    speed_rpm_mean=-rpm, signed_requested_rpm=-rpm,
                    S=950*worker["driver"].viscosity(temperature+.2, points)*rpm/60*.04*.043/load*(.0215/42.5e-6)**2,
                    bearing_residual_nm=.2, f=.2/(load*.0215), pressure_bar_abs_mean=13.))
            rows(preflight / "case_matrix.csv", measured)
            rows(preflight / "pure_coordinates.csv", [dict(S=i, f=i*.01) for i in range(103)])
            rows(preflight / "existing_cfd_coordinates.csv", [dict(S=i, f=i*.01) for i in range(6)])
            inputs = root / "inputs"
            worker["prepare"](preflight, inputs, self.caloric())
            report = root / "check.json"
            worker["check_inputs"](inputs, report)
            checked = json.loads(report.read_text())
            self.assertEqual(checked["status"], "34_matched_native_configs_verified_unrun")
            self.assertEqual(len(checked["cases"]), 17)
            experiment = json.loads((inputs / "experiment.json").read_text())
            self.assertEqual([case["measured"] for case in experiment["cases"]],
                             [{key: str(value) for key, value in row.items()} for row in measured])
            review = json.loads((inputs / "property_review.json").read_text())
            self.assertLess(review["maximum_fixed30_difference_percent"], 1e-10)
            self.assertEqual(review["source_chart_bases"],
                             dict(pressureBasis="unspecified", concentrationBasis="unspecified"))
            self.assertEqual((inputs / "pure_coordinates.csv").read_bytes(),
                             (preflight / "pure_coordinates.csv").read_bytes())
            with self.assertRaises(FileExistsError):
                worker["prepare"](preflight, inputs, self.caloric())
            reporter = runpy.run_path(str(ROOT / "cases/02f_lge_mixture_gas/summarize_gas.py"))
            reporter["summarize"](root, root / "report")
            result = json.loads((root / "report/report.json").read_text())
            self.assertEqual(result["accepted_pairs"], 0)
            self.assertEqual(result["focused_checks_verified"], 0)
            pairs = worker["load_rows"](root / "report/matched_comparison.csv")
            self.assertEqual(len(pairs), 17)
            self.assertTrue(all(row["paired_status"] == "unavailable" and
                                row["gas_on_minus_off_torque_nm"] == "" for row in pairs))
            self.assertEqual(len(worker["load_rows"](root / "report/focused_checks.csv")), 6)
            namespace = reporter["summarize"].__globals__
            original_accept = namespace["accepted"]
            def broken_observables(*_):
                return (dict(config={"omega": "10"}, observables={"TC_C": [30., 31., 32., 33.]}),
                        dict(torque_Nm=.2, alpha_gas_max=.1, dissolved_gas_min=.25,
                             cross_film_profile_t_max_K=310.))
            namespace["accepted"] = broken_observables
            try:
                reporter["summarize"](root, root / "late_failure_report")
            finally:
                namespace["accepted"] = original_accept
            failed = worker["load_rows"](root / "late_failure_report/matched_comparison.csv")
            self.assertTrue(all(row["paired_status"] == "unavailable" and row["gas_on_torque_nm"] == ""
                                for row in failed))
            # COUPLED has no legacy FBNS checkpoint; shared cooling/cache paths must accept that.
            options = worker["baseline"].options(inputs, root / "cache", None, None)
            options.table = str(inputs / "fw68d.csv")
            options.materials = str(inputs / "M16_gas_off_setup.json")
            campaign = worker["Campaign"](options)
            case = campaign.cases["M16"]
            eps, psi = campaign.setup["force_balance"]["initial_positions"]["M16"]
            ambient = case["tc0_boundary_temperature_k"]
            material = worker["driver"].case_materials(case, campaign.materials)
            hb, _ = worker["driver"].housing_coefficient(ambient+5, ambient, campaign.points, material)
            hs, _ = worker["driver"].shaft_coefficient(ambient+5, ambient,
                worker["driver"].angular_speed(case), campaign.points, material)
            run = campaign.work / "M16" / "fixture"
            run.mkdir(parents=True)
            (run / "config.txt").write_text(campaign.config_for(case, eps, psi, hb, hs, run / "results"))
            (run / "run.log").write_text("synthetic shared-cache fixture")
            record = dict(returncode=0, source_sha256=campaign.source_hash, executable_sha256=None,
                mpi_ranks=4, driver_sha256=campaign.driver_hash, rankfile_sha256=None, native_equilibrium=True,
                config_sha256=worker["driver"].sha256(run / "config.txt"), status="load_and_thermal_closed",
                run_dir=str(run), seed_eps=eps, seed_psi_deg=psi, h_bush_W_m2_K=hb, h_shaft_W_m2_K=hs,
                h_relative_mismatch=0., T_change_K=0., T_outer_mean_K=ambient+5, T_shaft_mean_K=ambient+5,
                artifact_sha256=worker["run_artifacts"](run), **campaign.input_hashes)
            campaign.completed["M16"].append((run, record))
            campaign.check_solution = lambda *_: {}
            self.assertIsNotNone(campaign.recheck(run, record, case))
            self.assertIsNotNone(campaign.find_completed(case, eps, psi, closed=True))
            self.assertEqual(campaign.balance(case, eps, psi)["run_dir"], str(run))
            self.assertIsNone(campaign.recover_checkpoint(case))
            self.assertFalse((run / "results/steady_state.dat").exists())


if __name__ == "__main__":
    unittest.main()
