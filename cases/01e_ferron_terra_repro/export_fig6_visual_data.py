"""Export compact profiles and a summary for the Ferron Figure 6 campaign."""

import json
import sys
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S

CASES = (
    "f6_depth_averaged",
    "f6_xf_adiabatic",
    "f6_xf_iso40",
    "f6_cht_fixed",
)


def main():
    equilibrium = json.loads((CASE_DIR / "fig6_equilibrium.json").read_text())
    profiles = {}
    summary = {}
    for name in CASES:
        result = equilibrium[name]
        run_dir = CASE_DIR / result["flow"]["dir"]
        pressure = S.midplane_profile(run_dir, "p")
        temperature = S.midplane_profile(run_dir, "T")
        bearing = S.midplane_profile(run_dir, "T_bearing")
        if pressure:
            profiles[name + "_theta_p_deg"] = (np.asarray(pressure[0]) - 90.0) % 360.0
            profiles[name + "_p_gauge_MPa"] = (np.asarray(pressure[1]) - 101325.0) / 1e6
        if temperature:
            profiles[name + "_theta_Tbar_deg"] = (np.asarray(temperature[0]) - 90.0) % 360.0
            profiles[name + "_T_bar_C"] = np.asarray(temperature[1]) - 273.15
        if bearing:
            profiles[name + "_theta_T_deg"] = (np.asarray(bearing[0]) - 90.0) % 360.0
            profiles[name + "_T_bearing_C"] = np.asarray(bearing[1]) - 273.15
        metrics = result["flow"]["metrics"]
        summary[name] = {
            "eps": result["flow"]["eps"],
            "psi_deg": result["flow"]["psi"],
            "load_N": metrics["W"],
            "p_max_gauge_MPa": (metrics["p_max"] - 101325.0) / 1e6,
            "T_film_max_C": metrics["T_max"] - 273.15,
            "configured_Qi_Lmin": result["q_supply_m3s"] * 60000.0,
            "diagnostic_source_Lmin": metrics["q_supply"] * 60000.0,
        }
    np.savez_compressed(CASE_DIR / "fig6_visual_data.npz", **profiles)
    (CASE_DIR / "fig6_analysis.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
