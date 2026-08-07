"""Compare Claude's archived G6 result with the corrected G6 result."""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
CASES_DIR = CASE_DIR.parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import COLORS, MARKERS, apply_style, save_figure, smart_legend


def read_reference(name):
    """Read one two-column Ferron reference trace."""
    rows = []
    with (CASES_DIR / "ref_data" / name).open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if row and not row[0].lstrip().startswith("#"):
                rows.append((float(row[0]), float(row[1])))
    return np.asarray(rows)


def sorted_profile(theta, values):
    """Return a profile sorted in groove-relative angle."""
    order = np.argsort(theta)
    return np.asarray(theta)[order], np.asarray(values)[order]


def periodic_interp(x, xp, values):
    """Interpolate a periodic profile without a seam at zero degrees."""
    xp, values = sorted_profile(xp, values)
    return np.interp(
        x,
        np.concatenate((xp - 360.0, xp, xp + 360.0)),
        np.tile(values, 3),
    )


def panel_label(ax, label, x=0.975, y=0.94, va="top", color="black"):
    ax.text(x, y, label, transform=ax.transAxes,
            ha="right", va=va, color=color)


def main():
    apply_style()
    archived = json.loads(
        (CASES_DIR / "01_ferron_wallbc_fields" / "extract.json").read_text()
    )["rerun_cht_robin10"]["profile"]
    temperature_ref = read_reference("ferron1983_fig4_2000rpm_4kN_temperature.csv")
    pressure_ref = read_reference("ferron1983_fig4_2000rpm_4kN_pressure.csv")

    with np.load(CASE_DIR / "visual_data.npz") as corrected:
        theta_new_t, temperature_new = sorted_profile(
            corrected["G6_theta_temperature_deg"],
            corrected["G6_temperature_bearing_C"],
        )
        theta_new_p, pressure_new = sorted_profile(
            corrected["G6_theta_pressure_deg"],
            corrected["G6_pressure_gauge_MPa"],
        )

    theta_old = np.asarray(archived["x_deg"])
    temperature_old = np.asarray(archived["T_bearing_C"])
    pressure_old = np.asarray(archived["p_MPa"])
    old_on_new_t = periodic_interp(theta_new_t, theta_old, temperature_old)
    old_on_new_p = periodic_interp(theta_new_p, theta_old, pressure_old)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.0), sharex="col")
    old_color = COLORS[1]
    new_color = COLORS[0]

    ax = axes[0, 0]
    ax.plot(temperature_ref[:, 0], temperature_ref[:, 1], linestyle="none",
            marker=MARKERS[1], color="black", label="Ferron et al. (1983)", zorder=5)
    ax.plot(theta_old, temperature_old, color=old_color,
            label=r"Claude, $Q_s=22.04$ L/min")
    ax.plot(theta_new_t, temperature_new, color=new_color,
            label=r"Corrected, $Q_s=5.06$ L/min")
    ax.set_ylabel(r"$T_b$ ($^\circ$C)")
    ax.set_ylim(39.0, 51.5)
    panel_label(ax, r"$\mathrm{(a)}$")
    smart_legend(ax, fontsize=11, loc="upper left")

    ax = axes[0, 1]
    ax.plot(pressure_ref[:, 0], pressure_ref[:, 1] / 1.0e6, linestyle="none",
            marker=MARKERS[1], color="black", label="Ferron et al. (1983)", zorder=5)
    ax.plot(theta_old, pressure_old, color=old_color, label="Claude")
    ax.plot(theta_new_p, pressure_new, color=new_color, label="Corrected")
    ax.set_ylabel(r"$p-p_{\mathrm{amb}}$ (MPa)")
    ax.set_ylim(-0.05, 1.50)
    panel_label(ax, r"$\mathrm{(b)}$")
    smart_legend(ax, fontsize=11, loc="upper left")

    ax = axes[1, 0]
    ax.plot(theta_new_t, temperature_new - old_on_new_t, color=new_color)
    ax.axhline(0.0, color="0.55", linewidth=0.9, linestyle="--")
    ax.set_ylabel(r"$T_{b,\mathrm{corr}}-T_{b,\mathrm{Claude}}$ (K)")
    panel_label(ax, r"$\mathrm{(c)}$", y=0.08, va="bottom")

    ax = axes[1, 1]
    ax.plot(theta_new_p, pressure_new - old_on_new_p, color=new_color)
    ax.axhline(0.0, color="0.55", linewidth=0.9, linestyle="--")
    ax.set_ylabel(r"$p_{\mathrm{corr}}-p_{\mathrm{Claude}}$ (MPa)")
    panel_label(ax, r"$\mathrm{(d)}$")

    for ax in axes.flat:
        ax.set_xlim(0.0, 360.0)
        ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    for ax in axes[1, :]:
        ax.set_xlabel(r"$\theta-\theta_g$ (deg)")
    fig.subplots_adjust(wspace=0.32, hspace=0.16)
    save_figure(fig, str(CASE_DIR / "ferron_claude_corrected_comparison"))
    plt.close(fig)


if __name__ == "__main__":
    main()
