"""Plot D-grid and D-dt convergence from local case-01e analyses."""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import apply_style, blues, save_figure, smart_legend


def sorted_xy(x, y):
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order]


def timestep_plot(data):
    rows = sorted(data["runs"].values(), key=lambda row: row["dt_s"], reverse=True)
    for key, label in (("eps_final", r"$\epsilon_f$"),
                       ("T_max_C", r"$T_{\max}$ ($^\circ$C)")):
        fig, axis = plt.subplots()
        axis.plot([row["dt_s"] for row in rows],
                  [row.get(key, float("nan")) for row in rows],
                  marker="o", label="D-dt")
        axis.set_xscale("log")
        axis.set_xlabel(r"$\Delta t$ (s)")
        axis.set_ylabel(label)
        smart_legend(axis, fontsize=11)
        save_figure(fig, CASE_DIR / f"terra_d_timestep_{key}")
        plt.close(fig)


def grid_plot(data):
    rows = sorted(
        (row for name, row in data.get("runs", {}).items()
         if name.startswith("d_grid_") and row.get("status") == "PASS"),
        key=lambda row: row["grid"][0],
    )
    if not rows:
        return
    n_theta = [row["grid"][0] for row in rows]
    quantities = (
        ("eps", r"$\epsilon$"),
        ("pressure_peak_gauge_MPa", r"$p_{\max}-p_{\mathrm{amb}}$ (MPa)"),
        ("peak_temperature_C", r"$T_{b,\max}$ ($^\circ$C)"),
        ("mean_error_K", r"$\overline{T_b-T_{b,\mathrm{exp}}}$ (K)"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)
    for axis, (key, label) in zip(axes.flat, quantities):
        axis.plot(n_theta, [row[key] for row in rows], marker="s")
        axis.set_xlabel(r"$N_\theta$")
        axis.set_ylabel(label)
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    save_figure(fig, CASE_DIR / "terra_d_grid_convergence")
    plt.close(fig)


def field_plots():
    visual_path = CASE_DIR / "d_grid_visual_data.npz"
    if not visual_path.exists():
        return
    with np.load(visual_path) as data:
        colors = blues(3)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
        for color, grid in zip(colors, ("g5", "g6", "g7")):
            tag = "d_grid_" + grid
            theta_t, temperature = sorted_xy(
                data[tag + "_theta_T_deg"], data[tag + "_T_bearing_C"])
            theta_p, pressure = sorted_xy(
                data[tag + "_theta_p_deg"], data[tag + "_p_gauge_MPa"])
            axes[0].plot(theta_t, temperature, color=color, label=grid.upper())
            axes[1].plot(theta_p, pressure, color=color, label=grid.upper())
        axes[0].set_xlabel(r"$\theta-\theta_g$ (deg)")
        axes[1].set_xlabel(r"$\theta-\theta_g$ (deg)")
        axes[0].set_ylabel(r"$T_b$ ($^\circ$C)")
        axes[1].set_ylabel(r"$p-p_{\mathrm{amb}}$ (MPa)")
        smart_legend(axes[0], fontsize=11)
        smart_legend(axes[1], fontsize=11)
        save_figure(fig, CASE_DIR / "terra_d_pressure_temperature_profiles")
        plt.close(fig)

        tag = "d_grid_g6"
        theta = data[tag + "_theta_deg"]
        z_over_l = data[tag + "_z_over_L"]
        fields = (
            (tag + "_plane_p_gauge_MPa", r"$p-p_{\mathrm{amb}}$ (MPa)", "cividis"),
            (tag + "_plane_T_C", r"$T$ ($^\circ$C)", "inferno"),
        )
        fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.6), sharex=True,
                                 constrained_layout=True)
        for axis, (key, label, cmap) in zip(axes, fields):
            contour = axis.contourf(theta, z_over_l, data[key], levels=24, cmap=cmap)
            fig.colorbar(contour, ax=axis, pad=0.025, label=label)
            axis.set_ylabel(r"$z/L$")
            axis.set_ylim(-0.5, 0.5)
        axes[1].set_xlabel(r"$\theta-\theta_g$ (deg)")
        save_figure(fig, CASE_DIR / "terra_d_g6_contours")
        plt.close(fig)


def main():
    apply_style()
    campaign_path = CASE_DIR / "d_grid_campaign_analysis.json"
    if not campaign_path.exists():
        campaign_path = CASE_DIR / "campaign_analysis.json"
    timestep_path = CASE_DIR / "d_discretization_analysis.json"
    if campaign_path.exists():
        grid_plot(json.loads(campaign_path.read_text(encoding="utf-8")))
    field_plots()
    if timestep_path.exists():
        timestep_plot(json.loads(timestep_path.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
