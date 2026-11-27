"""Set T figures: minimum film thickness and load over the last (periodic) revolution.

  flooded_transient_hmin   h_min against crank angle for each fluid, steady h_min at the minimum
                           and maximum cycle load (dashed), contact criterion 0.5 um
  flooded_transient_load   load magnitude against crank angle (same x axis, for stacking)

Usage: python plot_flooded_transient.py --work CAMPAIGN
(run directories in summary.json are remote paths; they are mapped onto --work)
"""
import argparse
import json
from pathlib import Path

import numpy as np

import common
from plotstyle import MARKERS, MUTED, PRL, figure, save, setup

LABEL = {"F1": r"$c=0.50$", "F2": r"$c=0.60$"}
COLOR = {"F1": PRL[0], "F2": PRL[1]}


def last_motion(work, fluid):
    summary = json.loads((work / f"T_{fluid}" / "summary.json").read_text())
    run_dir = work / f"T_{fluid}" / Path(summary["last_revolution_dir"]).name
    m = common.read_motion(run_dir)
    t = np.asarray(m["time_s"])
    step = t[1] - t[0]
    phase = (t - t[0] + step) / (t[-1] - t[0] + step) * 360.0   # 0 < phase <= 360, no wrap
    return summary, phase, m


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    args = parser.parse_args()
    out = args.work / "figures"
    setup()
    data = {f: last_motion(args.work, f) for f in ("F1", "F2")}

    fig, ax = figure(w=13.0)
    for f, (s, phase, m) in data.items():
        ax.plot(phase, m["h_min_m"] * 1e6, color=COLOR[f], lw=1.0, label=LABEL[f])
        for key in ("min", "max"):
            ax.axhline(s["steady"][key]["h_min_um"], color=COLOR[f], lw=0.5, ls="--")
        ax.text(357, s["steady"]["min"]["h_min_um"] + 0.004, "steady, min load", fontsize=6.5, color=COLOR[f],
                ha="right", va="bottom")
        ax.text(357, s["steady"]["max"]["h_min_um"] - 0.004, "steady, max load", fontsize=6.5, color=COLOR[f],
                ha="right", va="top")
    ax.axhline(0.5, color=MUTED, lw=0.8)
    ax.text(30, 0.505, r"contact criterion, $h_\mathrm{min}=0.5\ \mathrm{µm}$", fontsize=6.5, color=MUTED)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 60))
    ax.set_ylim(0.4, 1.1)
    ax.set_xlabel(r"Crank angle, $\phi$ ($°$)")
    ax.set_ylabel(r"Minimum film thickness, $h_\mathrm{min}$ ($\mathrm{µm}$)")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62))
    save(fig, str(out / "flooded_transient_hmin"))

    fig, ax = figure(w=13.0)
    s, phase, m = data["F1"]
    ax.plot(phase, np.hypot(m["load_x_N"], m["load_y_N"]), color=PRL[0], lw=1.0)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 60))
    ax.set_xlabel(r"Crank angle, $\phi$ ($°$)")
    ax.set_ylabel(r"Bearing load, $W$ ($\mathrm{N}$)")
    save(fig, str(out / "flooded_transient_load"))
    print(out)


if __name__ == "__main__":
    main()
