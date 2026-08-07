"""Export compact profile and G6 contour data from the remote VTK results."""

import json
import sys
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
SOURCE_CASE = CASE_DIR.parent / "01_ferron_wallbc_fields"
sys.path.insert(0, str(SOURCE_CASE))

import solver_io as S


def field_plane(run_dir, name):
    """Return one cell field on an ordered (z, groove-relative theta) grid."""
    arrays = S._read_pvts(S.final_pvts(run_dir))
    theta_raw = np.degrees(np.asarray(arrays["theta_rad"][0])) % 360.0
    theta_rel_raw = (theta_raw - 90.0) % 360.0
    z_raw = np.asarray(arrays["z_m"][0])
    values_raw = np.asarray(arrays[name][0])

    theta = np.unique(theta_rel_raw)
    z = np.unique(z_raw)
    values = np.full((len(z), len(theta)), np.nan)
    i_theta = np.searchsorted(theta, theta_rel_raw)
    i_z = np.searchsorted(z, z_raw)
    values[i_z, i_theta] = values_raw
    if np.isnan(values).any():
        raise RuntimeError(f"incomplete {name} plane")
    return theta, z, values


def main():
    equilibria = json.loads((CASE_DIR / "grid_equilibrium.json").read_text())
    output = {}
    for tag in ("G5", "G6", "G7"):
        run_dir = CASE_DIR / equilibria[tag]["flow"]["dir"]
        theta_t, temperature = S.midplane_profile(run_dir, "T_bearing")
        theta_p, pressure = S.midplane_profile(run_dir, "p")
        output[f"{tag}_theta_temperature_deg"] = (np.asarray(theta_t) - 90.0) % 360.0
        output[f"{tag}_temperature_bearing_C"] = np.asarray(temperature) - 273.15
        output[f"{tag}_theta_pressure_deg"] = (np.asarray(theta_p) - 90.0) % 360.0
        output[f"{tag}_pressure_gauge_MPa"] = (np.asarray(pressure) - 101325.0) / 1.0e6

    run_dir = CASE_DIR / equilibria["G6"]["flow"]["dir"]
    theta, z, pressure = field_plane(run_dir, "p")
    theta_t, z_t, temperature = field_plane(run_dir, "T")
    if not (np.array_equal(theta, theta_t) and np.array_equal(z, z_t)):
        raise RuntimeError("G6 pressure and temperature grids differ")
    output["G6_theta_deg"] = theta
    output["G6_z_over_L"] = z / 0.08 - 0.5
    output["G6_plane_pressure_gauge_MPa"] = (pressure - 101325.0) / 1.0e6
    output["G6_plane_temperature_C"] = temperature - 273.15

    np.savez_compressed(CASE_DIR / "visual_data.npz", **output)


if __name__ == "__main__":
    main()
