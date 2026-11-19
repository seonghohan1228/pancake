"""Figure of the last B2 revolution against crank angle (small multiples, one shared x axis).

Usage: python plot_part_b.py --work CAMPAIGN_DIR
Reads summary.json, per_crank_angle.csv and the last-revolution motion.csv; writes
CAMPAIGN_DIR/os_bearing_130hz_cycle.png.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import common  # noqa: E402

SERIES_1, SERIES_2 = "#2a78d6", "#eb6834"         # categorical slots 1 and 2 (reference palette)
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PERIOD = 1.0 / 130.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    args = parser.parse_args()
    summary = json.loads((args.work / "summary.json").read_text())
    motion = common.read_motion(Path(summary["b2"]["last_revolution_dir"].replace(
        "/home/user/pancake-runs/lge-compressor-transient-20261119-02a4b6c", str(args.work))))
    rows = list(csv.DictReader((args.work / "per_crank_angle.csv").open()))
    crank = (motion["time_s"] - motion["time_s"][0] + (motion["time_s"][1] - motion["time_s"][0])) / PERIOD * 360.0
    load = np.hypot(motion["load_x_N"], motion["load_y_N"])
    angle = np.array([float(r["crank_angle_deg"]) for r in rows])

    def column(key):
        return np.array([float(r[key]) if r.get(key) not in (None, "", "None") else np.nan for r in rows])

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "axes.facecolor": SURFACE,
                         "figure.facecolor": SURFACE})
    fig, axes = plt.subplots(5, 1, figsize=(7.0, 9.0), sharex=True, dpi=200)
    panels = [
        (axes[0], crank, load / 1e3, "Load on OS [kN]"),
        (axes[1], crank, motion["h_min_m"] * 1e6, "Minimum film [µm]"),
        (axes[2], crank, motion["friction_torque_Nm"], "Friction torque [N m]"),
    ]
    for ax, x, y, label in panels:
        ax.plot(x, y, color=SERIES_1, lw=1.6)
        ax.set_ylabel(label)
    # Quasi-steady reference: steady h_min at the instantaneous load (3 B0 points, linear).
    steady = summary["steady"]
    points = sorted((steady[k]["load_N"], steady[k]["h_min_um"]) for k in steady)
    quasi = np.interp(load, [p[0] for p in points], [p[1] for p in points])
    axes[1].lines[0].set_label("transient")
    axes[1].plot(crank, quasi, color=SERIES_2, lw=1.6, label="steady at the same load")
    axes[1].legend(frameon=False, loc="lower right", fontsize=8, ncol=2)
    ax = axes[3]
    start, end = column("bubble_start_from_hole_deg"), column("bubble_end_from_hole_deg")
    ax.plot(angle, start, color=SERIES_1, lw=1.6, marker="o", ms=3, label="bubble start")
    ax.plot(angle, end, color=SERIES_2, lw=1.6, marker="o", ms=3, label="bubble end")
    ax.axhline(0, color=INK_2, lw=0.8, ls="--")
    ax.set_ylabel("Bubble edge from\nhole [deg]")
    ax.legend(frameon=False, loc="upper right", fontsize=8, ncol=2)
    ax = axes[4]
    ax.plot(angle, column("bubble_gap_to_hole_edge_mm"), color=SERIES_1, lw=1.6, marker="o", ms=3)
    ax.axhline(0, color=INK_2, lw=0.8, ls="--")
    ax.set_ylabel("Bubble gap to\nhole edge [mm]")
    ax.set_xlabel("Crank angle [deg]")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 45))
    for ax in axes:
        ax.grid(color=GRID, lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("OS bearing, 130 Hz, GC+VC: last revolution", color=INK, fontsize=10)
    fig.tight_layout()
    out = args.work / "os_bearing_130hz_cycle.png"
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()
