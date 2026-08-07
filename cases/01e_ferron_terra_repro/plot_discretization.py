"""Plot D-grid and D-dt convergence from local case-01e analyses."""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
STYLE_DIR = Path.home() / ".codex" / "skills" / "plot-data" / "scripts"
sys.path.insert(0, str(STYLE_DIR))

from plot_style import apply_style, save_figure, smart_legend


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
    fig, axis = plt.subplots()
    axis.plot([row["grid"][0] for row in rows],
              [row["peak_temperature_C"] for row in rows],
              marker="s", label="D-grid")
    axis.set_xlabel(r"$N_\theta$")
    axis.set_ylabel(r"$T_{b,\max}$ ($^\circ$C)")
    smart_legend(axis, fontsize=11)
    save_figure(fig, CASE_DIR / "terra_d_grid_temperature")
    plt.close(fig)


def main():
    apply_style()
    campaign_path = CASE_DIR / "campaign_analysis.json"
    timestep_path = CASE_DIR / "d_discretization_analysis.json"
    if campaign_path.exists():
        grid_plot(json.loads(campaign_path.read_text(encoding="utf-8")))
    if timestep_path.exists():
        timestep_plot(json.loads(timestep_path.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
