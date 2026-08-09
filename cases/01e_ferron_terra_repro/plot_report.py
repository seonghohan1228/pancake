"""Create the requested grid, timestep, pressure, and CHT report figures."""

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

from plot_style import MARKERS, apply_style, blues, save_figure, smart_legend


def convergence_limits(axis, values):
    values = np.asarray(values, dtype=float)
    span = float(np.max(values) - np.min(values))
    if span == 0.0:
        span = max(abs(float(values[0])) * 0.01, 1e-6)
    axis.set_ylim(float(np.min(values) - 0.40 * span),
                  float(np.max(values) + 0.10 * span))


def sorted_xy(x, y):
    order = np.argsort(x)
    return np.asarray(x)[order], np.asarray(y)[order]


def reference(name):
    path = CASES_DIR / "ref_data" / name
    with path.open(encoding="utf-8") as stream:
        return np.asarray([[float(a), float(b)] for a, b in csv.reader(
            line for line in stream if not line.startswith("#"))])


def grid_convergence():
    analysis = json.loads((CASE_DIR / "d_grid_campaign_analysis.json").read_text())
    rows = [analysis["runs"][f"d_grid_g{i}"] for i in range(1, 9)
            if f"d_grid_g{i}" in analysis["runs"]]
    cells = np.asarray([np.prod(row["grid"]) for row in rows])
    series = (
        ([row["eps"] for row in rows], r"$\epsilon$"),
        ([row["pressure_peak_gauge_MPa"] for row in rows],
         r"$p_{\max}-p_{\mathrm{amb}}$ (MPa)"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    for axis, (values, label) in zip(axes, series):
        axis.plot(cells, values, linestyle="none", marker=MARKERS[0])
        axis.set_xlabel(r"$N$")
        axis.set_ylabel(label)
        convergence_limits(axis, values)
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    save_figure(fig, CASE_DIR / "terra_report_grid_convergence")
    plt.close(fig)


def legacy_grid_convergence():
    analysis = json.loads((CASE_DIR / "legacy_gct_analysis.json").read_text())
    rows = analysis["groups"]["isothermal"]["runs"]
    cells = np.asarray([np.prod(row["grid"]) for row in rows])
    series = (
        ([row["W_N"] / 1e3 for row in rows], r"$W$ (kN)"),
        ([row["p_max_MPa"] for row in rows], r"$p_{\max}$ (MPa)"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    for axis, (values, label) in zip(axes, series):
        axis.plot(cells, values, linestyle="none", marker=MARKERS[0])
        axis.set_xlabel(r"$N$")
        axis.set_ylabel(label)
        convergence_limits(axis, values)
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    save_figure(fig, CASE_DIR / "terra_report_legacy_grid_convergence")
    plt.close(fig)


def timestep_convergence():
    old_data = json.loads((CASE_DIR / "d_timestep_history.json").read_text())
    extended_path = CASE_DIR / "d_timestep_100ms_history.json"
    data = json.loads(extended_path.read_text()) if extended_path.exists() else old_data
    extended_names = ("d_dt100_t1", "d_dt100_t2", "d_dt100_t3", "d_dt100_t4")
    fallback_names = ("d_dt_t1", "d_dt_t2", "d_dt_t3", "d_dt_t4")
    selected_names = extended_names if extended_path.exists() else fallback_names
    histories = [data[name] for name in selected_names if name in data]
    colors = blues(len(histories))
    means = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), constrained_layout=True)
    for color, run in zip(colors, histories):
        time_ms = np.asarray([row["time_s"] for row in run["history"]]) * 1e3
        eps = np.asarray([row["eps"] for row in run["history"]])
        axes[0].plot(time_ms, eps, color=color,
                     label=rf"$\Delta t={run['dt_us']:.0f}\,\mu\mathrm{{s}}$")
        means.append(float(np.mean(eps[time_ms >= time_ms[-1] - 50.0])))
    axes[0].set_xlabel(r"$t$ (ms)")
    axes[0].set_ylabel(r"$\epsilon$")
    smart_legend(axes[0], fontsize=11)

    dt_us = [run["dt_us"] for run in histories]
    axes[1].plot(dt_us, means, linestyle="none", marker=MARKERS[0])
    axes[1].set_xlabel(r"$\Delta t$ ($\mu$s)")
    axes[1].set_ylabel(r"$\overline{\epsilon}$")
    mean_max = max(means)
    axes[1].set_ylim(mean_max - 0.60 * 0.01, mean_max + 0.40 * 0.01)
    save_figure(fig, CASE_DIR / "terra_report_timestep_convergence")
    plt.close(fig)

    fig, axis = plt.subplots()
    stability_names = ("d_dt_tn2", "d_dt_tn1", "d_dt_t0")
    stability_runs = [old_data[name] for name in stability_names if name in old_data]
    stability_runs.extend(histories)
    dt_all = [run["dt_us"] for run in stability_runs]
    reached_ms = [run["history"][-1]["time_s"] * 1e3 for run in stability_runs]
    axis.plot(dt_all, reached_ms, linestyle="none", marker=MARKERS[0])
    axis.set_xlabel(r"$\Delta t$ ($\mu$s)")
    axis.set_ylabel(r"$t_{\mathrm{reached}}$ (ms)")
    convergence_limits(axis, reached_ms)
    save_figure(fig, CASE_DIR / "terra_report_timestep_stability")
    plt.close(fig)


def pressure_validation():
    experiment = reference("ferron1983_fig4_2000rpm_4kN_pressure.csv")
    with np.load(CASE_DIR / "d_grid_visual_data.npz") as current, \
            np.load(CASE_DIR / "visual_data.npz") as archived:
        theta_now, p_now = sorted_xy(current["d_grid_g7_theta_p_deg"],
                                     current["d_grid_g7_p_gauge_MPa"])
        theta_old, p_old = sorted_xy(archived["r0_g7_theta_p_deg"],
                                     archived["r0_g7_p_gauge_MPa"])
    fig, axis = plt.subplots()
    axis.plot(theta_now, p_now, label=r"R1: fixed outer $T_b$")
    axis.plot(theta_old, p_old, label=r"R0: outer $h=10$ W m$^{-2}$ K$^{-1}$")
    axis.plot(experiment[:, 0], experiment[:, 1] / 1e6,
              linestyle="none", marker=MARKERS[1], markerfacecolor="white",
              label="Ferron et al. (1983)")
    axis.set_xlim(0.0, 360.0)
    axis.set_xlabel(r"$\theta-\theta_g$ (deg)")
    axis.set_ylabel(r"$p-p_{\mathrm{amb}}$ (MPa)")
    smart_legend(axis, loc="upper left", fontsize=11)
    save_figure(fig, CASE_DIR / "terra_report_pressure_validation")
    plt.close(fig)


def cht_comparison():
    experiment = reference("ferron1983_fig4_2000rpm_4kN_temperature.csv")
    cases = (
        ("c_xf_adiabatic", "cross-film, both walls adiabatic"),
        ("c_xf_bush40_jad", r"cross-film, bush $40\,^\circ$C / journal adiabatic"),
        ("c_xf_bush40_j40", r"cross-film, both walls $40\,^\circ$C"),
        ("c_conj_fixed", r"conjugate bush, fixed outer $T_b=40\,^\circ$C"),
        ("c_conj_robin80", r"conjugate bush, outer $h=80$ W m$^{-2}$ K$^{-1}$"),
    )
    colors = blues(len(cases))
    with np.load(CASE_DIR / "visual_data.npz") as data:
        fig, axis = plt.subplots(figsize=(7.2, 5.2))
        for color, (tag, label) in zip(colors, cases):
            theta, temperature = sorted_xy(data[tag + "_theta_T_deg"],
                                           data[tag + "_T_bearing_C"])
            axis.plot(theta, temperature, color=color, label=label)
    axis.plot(experiment[:, 0], experiment[:, 1], linestyle="none",
              marker=MARKERS[2], markerfacecolor="white",
              label="Ferron et al. (1983)")
    axis.set_xlim(0.0, 360.0)
    axis.set_xlabel(r"$\theta-\theta_g$ (deg)")
    axis.set_ylabel(r"$T_b$ ($^\circ$C)")
    smart_legend(axis, fontsize=10)
    save_figure(fig, CASE_DIR / "terra_report_cht_comparison")
    plt.close(fig)


def main():
    apply_style()
    grid_convergence()
    legacy_grid_convergence()
    timestep_convergence()
    pressure_validation()
    cht_comparison()


if __name__ == "__main__":
    main()
