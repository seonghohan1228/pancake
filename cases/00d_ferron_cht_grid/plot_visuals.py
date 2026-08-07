"""Plot the corrected Ferron grid profiles and the G6 film-plane fields."""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import MARKERS, apply_style, blues, save_figure, smart_legend


def read_reference(name):
    path = CASE_DIR.parent / "ref_data" / name
    rows = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if row and not row[0].lstrip().startswith("#"):
                rows.append((float(row[0]), float(row[1])))
    return np.asarray(rows)


def sorted_profile(theta, values):
    order = np.argsort(theta)
    return theta[order], values[order]


def plot_profiles(data):
    temperature_ref = read_reference("ferron1983_fig4_2000rpm_4kN_temperature.csv")
    pressure_ref = read_reference("ferron1983_fig4_2000rpm_4kN_pressure.csv")
    colors = blues(3)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), sharex=True)

    ax = axes[0]
    ax.plot(temperature_ref[:, 0], temperature_ref[:, 1], linestyle="none",
            marker=MARKERS[1], color="black", label="Ferron et al. (1983)", zorder=5)
    for color, tag in zip(colors, ("G5", "G6", "G7")):
        theta, values = sorted_profile(
            data[f"{tag}_theta_temperature_deg"],
            data[f"{tag}_temperature_bearing_C"],
        )
        ax.plot(theta, values, color=color, label=rf"$\mathrm{{{tag}}}$")
    ax.set_ylabel(r"$T_b$ ($^\circ$C)")
    ax.set_ylim(39.0, 51.5)
    ax.text(0.975, 0.94, r"$\mathrm{(a)}$", transform=ax.transAxes,
            ha="right", va="top")
    smart_legend(ax, fontsize=12, loc="upper left")

    ax = axes[1]
    ax.plot(pressure_ref[:, 0], pressure_ref[:, 1] / 1.0e6, linestyle="none",
            marker=MARKERS[1], color="black", label="Ferron et al. (1983)", zorder=5)
    for color, tag in zip(colors, ("G5", "G6", "G7")):
        theta, values = sorted_profile(
            data[f"{tag}_theta_pressure_deg"],
            data[f"{tag}_pressure_gauge_MPa"],
        )
        ax.plot(theta, values, color=color, label=rf"$\mathrm{{{tag}}}$")
    ax.set_ylabel(r"$p-p_{\mathrm{amb}}$ (MPa)")
    ax.set_ylim(-0.05, 1.50)
    ax.text(0.975, 0.94, r"$\mathrm{(b)}$", transform=ax.transAxes,
            ha="right", va="top")
    smart_legend(ax, fontsize=12, loc="upper left")

    for ax in axes:
        ax.set_xlabel(r"$\theta-\theta_g$ (deg)")
        ax.set_xlim(0.0, 360.0)
        ax.set_xticks(np.arange(0.0, 361.0, 60.0))
    fig.subplots_adjust(wspace=0.30)
    save_figure(fig, str(CASE_DIR / "ferron_cht_grid_profiles"))
    plt.close(fig)


def plot_contours(data):
    theta = data["G6_theta_deg"]
    z_over_l = data["G6_z_over_L"]
    pressure = data["G6_plane_pressure_gauge_MPa"]
    temperature = data["G6_plane_temperature_C"]
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.6), sharex=True)

    levels_p = np.linspace(-0.005, 1.36, 28)
    contour = axes[0].contourf(theta, z_over_l, pressure, levels=levels_p,
                               cmap="cividis", extend="both")
    colorbar = fig.colorbar(contour, ax=axes[0], pad=0.025)
    colorbar.set_label(r"$p-p_{\mathrm{amb}}$ (MPa)")
    axes[0].text(0.018, 0.90, r"$\mathrm{(a)}$", transform=axes[0].transAxes,
                 color="white", ha="left", va="top")

    levels_t = np.linspace(40.0, 50.7, 28)
    contour = axes[1].contourf(theta, z_over_l, temperature, levels=levels_t,
                               cmap="inferno", extend="both")
    colorbar = fig.colorbar(contour, ax=axes[1], pad=0.025)
    colorbar.set_label(r"$T$ ($^\circ$C)")
    axes[1].text(0.018, 0.90, r"$\mathrm{(b)}$", transform=axes[1].transAxes,
                 color="white", ha="left", va="top")

    for ax in axes:
        ax.set_ylabel(r"$z/L$")
        ax.set_ylim(-0.5, 0.5)
        ax.set_yticks((-0.5, -0.25, 0.0, 0.25, 0.5))
        ax.set_xlim(0.0, 360.0)
    axes[1].set_xlabel(r"$\theta-\theta_g$ (deg)")
    axes[1].set_xticks(np.arange(0.0, 361.0, 60.0))
    fig.subplots_adjust(hspace=0.20)
    save_figure(fig, str(CASE_DIR / "ferron_cht_g6_contours"))
    plt.close(fig)


def main():
    apply_style()
    with np.load(CASE_DIR / "visual_data.npz") as data:
        plot_profiles(data)
        plot_contours(data)


if __name__ == "__main__":
    main()
