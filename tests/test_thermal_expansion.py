"""Check geometry limits, per-case isolation and signed comparison arithmetic."""
import copy
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "cases/02d_fw68d_steady_six"))
import run_steady as driver
import thermal_expansion as expansion


def main():
    expansion.self_check()
    with TemporaryDirectory() as directory:
        root = Path(directory)
        setup = json.loads(driver.SETUP.read_text())
        hot = expansion.expanded_geometry(65.7, 65.1)
        setup["geometry_by_case"] = {"B": hot}
        setup_path = root / "setup.json"
        driver.save_json(setup_path, setup)
        args = SimpleNamespace(materials=str(setup_path), solver=None, work=None, trial=None,
                               ranks=1, rankfile=None, coupling_max_iters=4000, coupling_relaxation=1.)
        campaign = driver.Campaign(args)
        b = campaign.cases["B"]
        nominal = copy.deepcopy(b)
        nominal.pop("geometry")
        cfg = campaign.config_for(b, .5, 320., 60., 100., root)
        old = campaign.config_for(nominal, .5, 320., 60., 100., root)
        assert driver.warm_start_identity(cfg) != driver.warm_start_identity(old)
        assert "geometry" not in campaign.cases["A"]
        assert driver.case_geometry(b, campaign.experiment)["radial_clearance_m"] == hot["radial_clearance_m"]
        hb, _ = driver.housing_coefficient(340., 335., campaign.points, driver.case_materials(b, campaign.materials))
        hb0, _ = driver.housing_coefficient(340., 335., campaign.points, campaign.materials)
        assert hb != hb0
        for bad in ({"coupling_tolerance": 1e-5}, {"fbns_tol": float("nan")}, {"steady_force_tolerance": 0.}, {"unknown": 1e-8}):
            setup["solver_tolerances"] = bad
            setup_path.write_text(json.dumps(setup))
            try:
                driver.Campaign(args)
            except ValueError:
                pass
            else:
                raise AssertionError(f"weakened or invalid tolerance accepted: {bad}")
    before, after, measured = .3, .301, .4
    signed_error = 100 * (after / measured - 1)
    removed = 100 * (1 - abs(after-measured)/abs(before-measured))
    assert abs(signed_error + 24.75) < 1e-12 and abs(removed-1.) < 1e-12
    print("thermal expansion geometry/isolation/error checks passed")


if __name__ == "__main__":
    main()
