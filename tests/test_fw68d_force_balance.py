"""Exercise the FW68D position search without running CFD."""

import json
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
driver = runpy.run_path(str(ROOT / "cases/02d_fw68d_steady_six/run_steady.py"))
Campaign = driver["Campaign"]
LOAD = 200.0
TARGET_EPS = 0.45
TARGET_PSI = 320.0


class MockCampaign:
    def __init__(self, work):
        self.work = work
        self.args = SimpleNamespace(max_force_iterations=18)
        self.failed_full_step = False
        (work / "A").mkdir()

    def best_previous_position(self, case, default):
        return default

    def evaluate(self, case, eps, psi, label):
        if label.endswith("_step1") and not self.failed_full_step:
            self.failed_full_step = True
            raise RuntimeError("mocked failed full step")
        fx = LOAD * (psi - TARGET_PSI) / 20.0
        fy = -LOAD * (1.0 + 4.0 * (eps - TARGET_EPS))
        return {"eps": eps, "psi_deg": psi, "Fx_N": fx, "Fy_N": fy,
                "force_error_rel": ((fx / LOAD) ** 2 + ((fy + LOAD) / LOAD) ** 2) ** 0.5,
                "run_dir": str(self.work / "A" / label)}


def main():
    outcomes = []
    assert driver["np"].array_equal(driver["applied_load"](LOAD), [0, LOAD])
    for start in ((0.2, 280.0), (0.7, 350.0)):
        with TemporaryDirectory() as directory:
            mock = MockCampaign(Path(directory))
            result = Campaign.balance(mock, {"case_id": "A", "load_mean_n": LOAD}, *start)
            history = json.loads(Path(result["force_history_path"]).read_text())
            assert mock.failed_full_step
            assert any("mocked failed full step" in item.get("rejected", "") for item in history)
            assert result["status"] == "load_and_thermal_closed"
            assert result["force_error_rel"] < 0.002
            assert abs(result["Fy_N"] + LOAD) / LOAD < 0.002
            assert abs(result["eps"] - TARGET_EPS) < 0.0005
            assert abs(result["psi_deg"] - TARGET_PSI) < 0.04
            outcomes.append((result["eps"], result["psi_deg"]))
    assert abs(outcomes[0][0] - outcomes[1][0]) < 0.0005
    assert abs(outcomes[0][1] - outcomes[1][1]) < 0.04
    with TemporaryDirectory() as directory:
        campaign = Campaign(SimpleNamespace(materials=driver["SETUP"], solver=None,
                            work=directory, ranks=1, coupling_max_iters=2000, coupling_relaxation=1.0))
        case = campaign.cases["E"]
        run = Path(directory) / "E" / "old_downward_load"
        (run / "results").mkdir(parents=True)
        config = campaign.config_for(case, .45, 140, 60, 100, run / "results")
        config = config.replace("external_load_y = " + str(float(case["load_mean_n"])),
                                "external_load_y = " + str(-float(case["load_mean_n"])))
        (run / "config.txt").write_text(config)
        (run / "run.log").write_text("")
        (run / "results/steady_state.dat").write_text("old orientation state")
        record = {"case_id": "E", "returncode": 0, "status": "thermally_closed",
                  "source_sha256": campaign.source_hash, "executable_sha256": None,
                  "mpi_ranks": 1, **campaign.input_hashes,
                  "config_sha256": driver["sha256"](run / "config.txt"),
                  "eps": .45, "psi_deg": 140, "force_error_rel": 0.0}
        campaign.completed["E"] = [(run, record)]
        assert campaign.matching_provenance(record)
        def forbid_old_validation(*args):
            raise AssertionError("old load reached output validation instead of orientation rejection")
        globals_ = Campaign.recheck.__globals__
        original_check = globals_["check_solution"]
        globals_["check_solution"] = forbid_old_validation
        try:
            assert campaign.recheck(run, record, case) is None
            assert campaign.find_completed(case, .45, 140, 60, 100) is None
            assert campaign.recover_checkpoint(case) is None
            assert campaign.best_previous_position(case, (.45, 320)) == (.45, 320)
        finally:
            globals_["check_solution"] = original_check
    print("FW68D force-balance passed: upward load, two starts/backtracking, downward equilibrium, old caches/restarts rejected")


if __name__ == "__main__":
    main()
