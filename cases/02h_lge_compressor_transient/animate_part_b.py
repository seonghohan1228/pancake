"""Contours and GIF of the last B2 revolution (fields every 5 deg crank angle).

Panels: pressure, void fraction (1 - film content), film thickness, and pressure minus its
cycle mean (the time-dependent part) on the unwrapped film
(theta' in the solver frame, OS bore moving towards +theta' relative to the pin; z axial).
The feed hole is outlined. Colour limits are fixed over the cycle.

Usage: python animate_part_b.py --work CAMPAIGN_DIR [--snapshots 5 105 180]
Writes CAMPAIGN_DIR/contours/frame_XXX.png, os_bearing_130hz_cycle.gif, os_bearing_130hz_contours.png.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402
from PIL import Image  # noqa: E402

import common  # noqa: E402
import part_b  # noqa: E402

INK, INK_2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
# Sequential: one hue, light to dark (reference blue ramp steps 100 -> 700).
BLUE = LinearSegmentedColormap.from_list("blue", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
# Diverging: blue <-> orange with a neutral gray midpoint.
DIVERGING = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#e4e3df", "#f19a72", "#b8461b"])
ORANGE = LinearSegmentedColormap.from_list("orange", ["#fdf3ee", "#f7c9b3", "#f19a72", "#eb6834", "#b8461b", "#7a2c0f"])


def grids(fields):
    nt, nz, theta, z, shape = part_b.grid_of(fields)
    return theta, z * 1e3, shape("p") / 1e5, 1 - shape("film_content"), shape("h") * 1e6


def draw(fig, axes, frame, limits, title, p_mean):
    theta, z, p, void, h = frame
    extent = [0, 360, z[0] - (z[1] - z[0]) / 2, z[-1] + (z[1] - z[0]) / 2]
    specs = [(p, "Pressure [bar]", BLUE, limits["p"]), (void, "Void fraction [-]", ORANGE, (0, 1)),
             (h, "Film thickness [µm]", BLUE.reversed(), limits["h"]),
             (p - p_mean, "Pressure - cycle mean [bar]", DIVERGING, (-limits["dp"], limits["dp"]))]
    hole_w = np.degrees(part_b.HOLE_RADIUS / part_b.R) * 2
    for ax, (data, label, cmap, (lo, hi)) in zip(axes, specs):
        ax.clear()
        image = ax.imshow(data, origin="lower", aspect="auto", extent=extent, cmap=cmap, vmin=lo, vmax=hi,
                          interpolation="nearest")
        ax.add_patch(Ellipse((part_b.HOLE_THETA, part_b.HOLE_Z * 1e3), hole_w, part_b.HOLE_RADIUS * 2e3,
                             fill=False, ec=INK, lw=1.0))
        ax.set_ylabel("z [mm]")
        ax.set_xlim(0, 360)
        ax.set_xticks(range(0, 361, 45))
        ax.tick_params(colors=INK_2)
        if not hasattr(ax, "_bar"):
            ax._bar = fig.colorbar(image, ax=ax, pad=0.015)
            ax._bar.set_label(label, color=INK_2)
        else:
            ax._bar.update_normal(image)
        ax.text(part_b.HOLE_THETA, part_b.HOLE_Z * 1e3 + part_b.HOLE_RADIUS * 1e3 + 1.0, "hole", color=INK,
                ha="center", fontsize=7)
    axes[-1].set_xlabel("θ' [deg] (OS bore moves towards +θ' relative to the pin)")
    fig.suptitle(title, color=INK, fontsize=10)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--snapshots", type=float, nargs="+", default=[5.0, 105.0, 180.0])
    args = parser.parse_args()
    summary = json.loads((args.work / "summary.json").read_text())
    frames_dir = args.work / "b2" / "last_fields"
    motion = common.read_motion(frames_dir)
    frames = []
    for t, pvts in part_b.frame_list(frames_dir):
        if t <= 0:
            continue
        index = int(np.argmin(np.abs(motion["time_s"] - t)))
        frames.append(dict(crank=t / part_b.PERIOD * 360.0, load=float(np.hypot(motion["load_x_N"][index], motion["load_y_N"][index])),
                           h_min=float(motion["h_min_m"][index]) * 1e6, data=grids(part_b.read_frame(pvts))))
    p_mean = np.mean([f["data"][2] for f in frames], axis=0)
    limits = dict(p=(part_b.P_D / 1e5, max(f["data"][2].max() for f in frames)),
                  h=(0.0, max(f["data"][4].max() for f in frames)),
                  dp=max(np.abs(f["data"][2] - p_mean).max() for f in frames))
    plt.rcParams.update({"font.size": 8, "figure.facecolor": SURFACE, "axes.edgecolor": INK_2,
                         "axes.labelcolor": INK_2})
    out = args.work / "contours"
    out.mkdir(exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 8.0), sharex=True, dpi=110)
    images = []
    for k, f in enumerate(frames):
        draw(fig, axes, f["data"], limits,
             f"Crank {f['crank']:5.1f} deg   load {f['load'] / 1e3:.3f} kN   h_min {f['h_min']:.3f} µm", p_mean)
        fig.tight_layout()
        path = out / f"frame_{k:03d}.png"
        fig.savefig(path)
        images.append(Image.open(path).convert("P", palette=Image.ADAPTIVE, colors=128))
    plt.close(fig)
    gif = args.work / "os_bearing_130hz_cycle.gif"
    images[0].save(gif, save_all=True, append_images=images[1:], duration=120, loop=0, optimize=True)
    # Static snapshot sheet: rows = crank angles, columns = fields.
    picks = [min(frames, key=lambda f: abs(f["crank"] - a)) for a in args.snapshots]
    fig, axes = plt.subplots(4, len(picks), figsize=(4.0 * len(picks), 8.0), sharex=True, sharey=True, dpi=200,
                             squeeze=False)
    for col, f in enumerate(picks):
        draw(fig, axes[:, col], f["data"], limits, "", p_mean)
        axes[0, col].set_title(f"Crank {f['crank']:.0f} deg, load {f['load'] / 1e3:.2f} kN, h_min {f['h_min']:.3f} µm",
                               color=INK, fontsize=8)
        if col:
            for ax in axes[:, col]:
                ax.set_ylabel("")
    fig.suptitle("OS bearing, 130 Hz, GC+VC: pressure, void fraction, film thickness", color=INK, fontsize=10)
    fig.tight_layout()
    sheet = args.work / "os_bearing_130hz_contours.png"
    fig.savefig(sheet)
    print(gif, sheet, len(frames), "frames")


if __name__ == "__main__":
    main()
