"""Plot Ferron Figure 6 temperatures against four thermal closures."""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
CASES_DIR = CASE_DIR.parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import COLORS, MARKERS, apply_style, save_figure, smart_legend


def sorted_xy(x, y):
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order]


def reference(name):
    with (CASES_DIR / "ref_data" / name).open(encoding="utf-8") as stream:
        return np.asarray([[float(a), float(b)] for a, b in csv.reader(
            line for line in stream if not line.startswith("#"))])


def main():
    apply_style()
    experiment = reference("ferron1983_fig6_4000rpm_6kN_temperature.csv")
    cases = (
        ("f6_depth_averaged", "theta_Tbar_deg", "T_bar_C",
         r"depth-averaged THD ($\bar{T}$, not a surface)", "0.45", ":"),
        ("f6_xf_adiabatic", "theta_T_deg", "T_bearing_C",
         r"cross-film THD, both walls adiabatic", COLORS[0], "-"),
        ("f6_xf_iso40", "theta_T_deg", "T_bearing_C",
         r"cross-film THD, both walls $40\,^\circ$C", COLORS[1], "-."),
        ("f6_cht_fixed", "theta_T_deg", "T_bearing_C",
         r"cross-film + CHT, fixed outer $T_b=40\,^\circ$C", COLORS[4], "-"),
    )
    with np.load(CASE_DIR / "fig6_visual_data.npz") as data:
        fig, axis = plt.subplots(figsize=(7.2, 5.2))
        for tag, theta_key, value_key, label, color, linestyle in cases:
            theta, temperature = sorted_xy(data[tag + "_" + theta_key],
                                           data[tag + "_" + value_key])
            axis.plot(theta, temperature, color=color, linestyle=linestyle,
                      linewidth=1.8, label=label)
    axis.plot(experiment[:, 0], experiment[:, 1], linestyle="none",
              marker=MARKERS[2], markerfacecolor="white",
              label="Ferron et al. (1983)")
    axis.set_xlim(0.0, 360.0)
    axis.set_xlabel(r"$\theta-\theta_g$ (deg)")
    axis.set_ylabel(r"$T$ ($^\circ$C)")
    smart_legend(axis, fontsize=9)
    save_figure(fig, CASE_DIR / "terra_fig6_temperature_comparison")
    plt.close(fig)


if __name__ == "__main__":
    main()
