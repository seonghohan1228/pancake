"""Check final cooling-seed lookup and run identity without executing CFD."""

import json
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
worker = runpy.run_path(str(ROOT / "cases/02e_lge_mixture/run_mixture.py"))


class MixtureWorkerTest(unittest.TestCase):
    def exercise(self, closure):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            setup = {"material_settings": {}, "force_balance": {
                "initial_positions": {"M16": [.15, 345.]}}}
            (inputs / "M16_setup.json").write_text(json.dumps(setup))
            run = root / "native_final_cooling"
            run.mkdir()
            (run / "config.txt").write_text("omega = 20\n")
            case = {"measured": {"selected_id": "M16"}}
            result = {"run_dir": str(run), "seed_eps": .274, "seed_psi_deg": 318.2,
                      "force_error_rel": 1e-5, "diagnostics": {"coupling_R_T": 1e-9},
                      "cross_film_profile_residual": 1e-9,
                      "energy_balance_residual_rel": 1e-9, "torque_Nm": .1}
            campaign = Mock()
            campaign.cases = {"M16": case}
            campaign.command.return_value = ["fake_mpiexec"]
            campaign.balance.return_value = result
            campaign.recheck.return_value = result
            campaign.find_completed.return_value = (None if closure is None else
                {"run_dir": str(run if closure == "same" else root / "another_run")})
            args = SimpleNamespace(inputs=inputs, work=root / "work", case="M16",
                solver=root / "fake_solver", rankfile=root / "fake_rankfile", layers=16, tight=False)
            observables = Mock(return_value={"TC_C": [30., 31., 32.]})
            with patch.object(worker["driver"], "Campaign", return_value=campaign), \
                 patch.object(worker["driver"], "read_fields", return_value={"p": [1., 2.]}), \
                 patch.object(worker["subprocess"], "check_output", return_value="Fake MPI\n"), \
                 patch.dict(worker["run_case"].__globals__, {"solid_observables": observables}), \
                 patch("builtins.print"):
                if closure == "same":
                    worker["run_case"](args)
                else:
                    with self.assertRaisesRegex(ValueError, "independent native and cooling recheck"):
                        worker["run_case"](args)
            campaign.balance.assert_called_once_with(case, .15, 345.)
            campaign.recheck.assert_called_once_with(run, result, case)
            campaign.find_completed.assert_called_once_with(case, .274, 318.2, closed=True)
            destination = args.work / args.case
            if closure == "same":
                accepted = json.loads((destination / "accepted.json").read_text())
                self.assertEqual(accepted["status"], "accepted_conditional_surrogate")
                self.assertEqual(accepted["result"]["run_dir"], str(run))
                self.assertTrue((destination / "fields.npz").is_file())
                self.assertFalse((destination / "failure.json").exists())
                observables.assert_called_once()
            else:
                self.assertFalse((destination / "accepted.json").exists())
                self.assertFalse((destination / "fields.npz").exists())
                self.assertEqual(json.loads((destination / "failure.json").read_text())["status"],
                                 "failed_preserved")
                observables.assert_not_called()

    def test_accepts_final_cooling_seed_instead_of_initial_seed(self):
        self.exercise("same")

    def test_rejects_missing_closure_or_another_run(self):
        for closure in (None, "different"):
            with self.subTest(closure=closure):
                self.exercise(closure)


if __name__ == "__main__":
    unittest.main()
