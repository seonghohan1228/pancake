"""Plot the repaired legacy GCT reproduction without a validation claim."""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import apply_style, save_figure, smart_legend


def convergence_plots(data):
    quantities = (
        ("W_N", r"$W$ (N)", "load"),
        ("p_max_MPa", r"$p_{\max}$ (MPa)", "pressure"),
        ("torque_Nm", r"$M_f$ (N m)", "torque"),
        ("T_max_C", r"$T_{\max}$ ($^\circ$C)", "temperature"),
    )
    for group, group_data in data["groups"].items():
        runs = group_data["runs"]
        for key, label, name in quantities:
            selected = [row for row in runs if row.get(key) is not None]
            if not selected:
                continue
            fig, axis = plt.subplots()
            axis.plot([row["grid"][0] for row in selected],
                      [row[key] for row in selected], marker="o",
                      label=group.replace("_", " "))
            axis.set_xlabel(r"$N_\theta$")
            axis.set_ylabel(label)
            smart_legend(axis, fontsize=11)
            save_figure(fig, CASE_DIR / f"terra_legacy_gct_{group}_{name}")
            plt.close(fig)


def contour_plots():
    visual_path = CASE_DIR / "legacy_gct_visual_data.npz"
    if not visual_path.exists():
        return
    import numpy as np
    with np.load(visual_path) as fields:
        for key in ("legacy_gct_g5_p", "legacy_gct_thd_g4_p", "legacy_gct_thd_g4_T"):
            if key + "_plane" not in fields:
                continue
            fig, axis = plt.subplots()
            level = axis.contourf(fields[key + "_theta_deg"],
                                  fields[key + "_z_over_L"],
                                  fields[key + "_plane"], levels=20,
                                  cmap="cividis" if key.endswith("_p") else "inferno")
            axis.set_xlabel(r"$\theta-\theta_g$ (deg)")
            axis.set_ylabel(r"$z/L$")
            label = r"$p-p_{\mathrm{amb}}$ (Pa)" if key.endswith("_p") else r"$T$ ($^\circ$C)"
            fig.colorbar(level, ax=axis, label=label)
            save_figure(fig, CASE_DIR / f"terra_{key}_contour")
            plt.close(fig)


def main():
    apply_style()
    data = json.loads((CASE_DIR / "legacy_gct_analysis.json").read_text(encoding="utf-8"))
    convergence_plots(data)
    contour_plots()


if __name__ == "__main__":
    main()
