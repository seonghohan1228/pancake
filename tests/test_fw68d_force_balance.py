"""Exercise the FW68D position search without running CFD."""

import json
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
Campaign = runpy.run_path(str(ROOT / "cases/02d_fw68d_steady_six/run_steady.py"))["Campaign"]
LOAD = 200.0
TARGET_EPS = 0.45
TARGET_PSI = 140.0


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
        fy = LOAD * (1.0 + 4.0 * (eps - TARGET_EPS))
        return {"eps": eps, "psi_deg": psi, "Fx_N": fx, "Fy_N": fy,
                "force_error_rel": ((fx / LOAD) ** 2 + ((fy - LOAD) / LOAD) ** 2) ** 0.5,
                "run_dir": str(self.work / "A" / label)}


def main():
    outcomes = []
    for start in ((0.2, 100.0), (0.7, 170.0)):
        with TemporaryDirectory() as directory:
            mock = MockCampaign(Path(directory))
            result = Campaign.balance(mock, {"case_id": "A", "load_mean_n": LOAD}, *start)
            history = json.loads(Path(result["force_history_path"]).read_text())
            assert mock.failed_full_step
            assert any("mocked failed full step" in item.get("rejected", "") for item in history)
            assert result["status"] == "load_and_thermal_closed"
            assert result["force_error_rel"] < 0.002
            assert abs(result["eps"] - TARGET_EPS) < 0.0005
            assert abs(result["psi_deg"] - TARGET_PSI) < 0.04
            outcomes.append((result["eps"], result["psi_deg"]))
    assert abs(outcomes[0][0] - outcomes[1][0]) < 0.0005
    assert abs(outcomes[0][1] - outcomes[1][1]) < 0.04
    print("FW68D force-balance mock passed: two starts, failed full step, backtracking, same equilibrium")


if __name__ == "__main__":
    main()
