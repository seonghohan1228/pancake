"""Re-close the Ferron conjugate equilibrium on G5/G6/G7 film grids."""

import argparse
import copy
import json
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
SOURCE_CASE = CASE_DIR.parent / "01_ferron_wallbc_fields"
sys.path.insert(0, str(SOURCE_CASE))

import force_balance as FB


PRODUCTION_GRIDS = {
    "G5": (480, 128),
    "G6": (720, 192),
    "G7": (1080, 288),
}
SMOKE_GRIDS = {"smoke": (180, 48)}


def configure_shared_driver():
    """Point the production driver at this case and require converged coupling."""
    FB.CASE_DIR = CASE_DIR
    FB.PROGRESS = CASE_DIR / "progress.log"
    FB.TEMPLATE = (SOURCE_CASE / "config_template.txt").read_text(encoding="utf-8")
    FB.COUPLING_MAX_ITERS = 200
    FB.COUPLING_RELAXATION = 1.0
    FB.COUPLING_EXTRA = """
coupling_relaxation_mode = ADAPTIVE
coupling_relaxation_min = 0.1
coupling_relaxation_backoff = 0.5
coupling_relaxation_growth = 1.05
coupling_keep_best_iterate = true
"""


def coupling_converged(result):
    for pass_name in ("pinned", "flow"):
        coupling = result[pass_name]["metrics"].get("coupling") or {}
        if (coupling.get("R_p", float("inf")) >= 1.0e-6
                or coupling.get("R_T", float("inf")) >= 1.0e-6):
            return False
    return True


def run(grids):
    configure_shared_driver()
    output_path = CASE_DIR / ("smoke_equilibrium.json" if "smoke" in grids
                              else "grid_equilibrium.json")
    results = json.loads(output_path.read_text(encoding="utf-8")) if output_path.exists() else {}
    production_spec = copy.deepcopy(FB.GROUPS["rerun_cht_robin10"])
    warm = production_spec["warm"]

    for tag, (n_theta, n_z) in grids.items():
        if tag in results:
            warm = (results[tag]["flow"]["eps"], results[tag]["flow"]["psi"])
            continue
        FB.N_THETA, FB.N_Z = n_theta, n_z
        spec = copy.deepcopy(production_spec)
        spec["warm"] = warm
        FB.GROUPS = {tag: spec}
        FB.S.progress(FB.PROGRESS,
                      f"{tag}: start {n_theta}x{n_z}, warm eps={warm[0]:.6f} psi={warm[1]:.4f}")
        result = FB.solve_point(tag, FB.W_TARGET, *warm)
        if not coupling_converged(result):
            raise RuntimeError(f"{tag}: force balance closed but pressure-temperature coupling did not")
        results[tag] = result
        output_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
        warm = (result["flow"]["eps"], result["flow"]["psi"])

    print(json.dumps({tag: {
        "grid": grids[tag],
        "eps": result["flow"]["eps"],
        "psi": result["flow"]["psi"],
        "run_dir": result["flow"]["dir"],
    } for tag, result in results.items()}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="run one unscored 180x48 harness check")
    args = parser.parse_args()
    run(SMOKE_GRIDS if args.smoke else PRODUCTION_GRIDS)


if __name__ == "__main__":
    main()
