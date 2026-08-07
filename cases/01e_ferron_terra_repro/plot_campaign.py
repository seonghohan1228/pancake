"""Create publication-style Terra/Ferron comparison figures."""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import MARKERS, apply_style, blues, save_figure, smart_legend


def reference(name):
    rows = []
    with (CASE_DIR.parent / "ref_data" / name).open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if row and not row[0].lstrip().startswith("#"):
                rows.append((float(row[0]), float(row[1])))
    return np.asarray(rows)


def sorted_xy(x, y):
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order]


def plot_profiles(data, tags):
    temperature_ref = reference("ferron1983_fig4_2000rpm_4kN_temperature.csv")
    pressure_ref = reference("ferron1983_fig4_2000rpm_4kN_pressure.csv")
    colors = blues(len(tags))
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), sharex=True)

    axes[0].plot(temperature_ref[:, 0], temperature_ref[:, 1], linestyle="none",
                 marker=MARKERS[1], color="black", label="Ferron et al. (1983)",
                 zorder=5)
    axes[1].plot(pressure_ref[:, 0], pressure_ref[:, 1] / 1.0e6,
                 linestyle="none", marker=MARKERS[1], color="black",
                 label="Ferron et al. (1983)", zorder=5)

    for color, tag in zip(colors, tags):
        label = tag.upper().replace("_", " ")
        x_t, t = sorted_xy(data[tag + "_theta_T_deg"], data[tag + "_T_bearing_C"])
        x_p, p = sorted_xy(data[tag + "_theta_p_deg"], data[tag + "_p_gauge_MPa"])
        axes[0].plot(x_t, t, color=color, label=label)
        axes[1].plot(x_p, p, color=color, label=label)

    axes[0].set_ylabel(r"$T_b$ ($^\circ$C)")
    axes[1].set_ylabel(r"$p-p_{\mathrm{amb}}$ (MPa)")
    axes[0].text(0.975, 0.94, r"$\mathrm{(a)}$", transform=axes[0].transAxes,
                 ha="right", va="top")
    axes[1].text(0.975, 0.94, r"$\mathrm{(b)}$", transform=axes[1].transAxes,
                 ha="right", va="top")
    for axis in axes:
        axis.set_xlabel(r"$\theta-\theta_g$ (deg)")
        axis.set_xlim(0.0, 360.0)
        axis.set_xticks(np.arange(0.0, 361.0, 60.0))
        smart_legend(axis, fontsize=10, loc="best")
    fig.subplots_adjust(wspace=0.30)
    save_figure(fig, str(CASE_DIR / "terra_ferron_profiles"))
    plt.close(fig)


def plot_contours(data, tag):
    theta = data[tag + "_theta_deg"]
    z_over_l = data[tag + "_z_over_L"]
    fields = (
        (tag + "_plane_p_gauge_MPa", r"$p-p_{\mathrm{amb}}$ (MPa)", "cividis"),
        (tag + "_plane_T_C", r"$T$ ($^\circ$C)", "inferno"),
    )
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.6), sharex=True)
    for panel, (key, label, cmap) in zip(axes, fields):
        contour = panel.contourf(theta, z_over_l, data[key], levels=24, cmap=cmap)
        colorbar = fig.colorbar(contour, ax=panel, pad=0.025)
        colorbar.set_label(label)
        panel.set_ylabel(r"$z/L$")
        panel.set_ylim(-0.5, 0.5)
    axes[0].text(0.018, 0.90, r"$\mathrm{(a)}$", transform=axes[0].transAxes,
                 color="white", ha="left", va="top")
    axes[1].text(0.018, 0.90, r"$\mathrm{(b)}$", transform=axes[1].transAxes,
                 color="white", ha="left", va="top")
    axes[1].set_xlabel(r"$\theta-\theta_g$ (deg)")
    axes[1].set_xlim(0.0, 360.0)
    axes[1].set_xticks(np.arange(0.0, 361.0, 60.0))
    fig.subplots_adjust(hspace=0.20)
    save_figure(fig, str(CASE_DIR / f"terra_{tag}_contours"))
    plt.close(fig)


def plot_peak_comparison(analysis):
    labels = []
    values = []
    colors = []
    for key, row in analysis["archive_comparison"].items():
        if key == "ferron_fig4":
            continue
        labels.append(key.replace("_", "\n"))
        values.append(row.get("peak_temperature_C", np.nan))
        colors.append("0.65")
    for key, row in analysis["runs"].items():
        if "peak_temperature_C" not in row:
            continue
        labels.append(key.replace("_", "\n"))
        values.append(row["peak_temperature_C"])
        colors.append("#0072B2")

    fig, axis = plt.subplots(figsize=(max(6.8, 0.62 * len(labels)), 4.5))
    x = np.arange(len(labels))
    axis.scatter(x, values, c=colors, marker=MARKERS[1])
    ferron_peak = analysis["archive_comparison"]["ferron_fig4"]["peak_temperature_C"]
    axis.axhline(ferron_peak, color="black", linestyle="--", linewidth=0.9,
                 label="Ferron et al. (1983)")
    axis.set_ylabel(r"$T_{b,\max}$ ($^\circ$C)")
    axis.set_xticks(x, labels, rotation=60, ha="right", fontsize=9)
    smart_legend(axis, fontsize=11)
    save_figure(fig, str(CASE_DIR / "terra_archive_ferron_peak_temperature"))
    plt.close(fig)


def main():
    apply_style()
    analysis = json.loads((CASE_DIR / "campaign_analysis.json").read_text(encoding="utf-8"))
    plot_peak_comparison(analysis)
    visual_path = CASE_DIR / "visual_data.npz"
    if not visual_path.exists():
        return
    with np.load(visual_path) as data:
        tags = [tag for tag in ("r0_g6", "r0_g7", "r1_g6", "r1_g7",
                                "r2_g6", "r2_g7")
                if tag + "_T_bearing_C" in data]
        if tags:
            plot_profiles(data, tags)
        for tag in ("r0_g6", "r1_g6", "r2_g6"):
            if tag + "_plane_T_C" in data:
                plot_contours(data, tag)


if __name__ == "__main__":
    main()
