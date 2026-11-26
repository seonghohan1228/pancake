"""GC+VC vs VC-only at the mean load: unwrapped film maps and centre-line pressure.

Files (in CAMPAIGN/figures):
  os_bearing_130hz_vc_vs_gcvc_mean.(pdf|svg|png)   2 x 2 maps: rows GC+VC, VC only; columns pressure,
                                                    film content; same colour scale per column
  os_bearing_130hz_centreline_pressure.(...)        centre-line pressure, both models
  os_bearing_130hz_centreline_pressure_zoom.(...)   same, 0-30 bar
Angle measured from the hole centre in the direction the OS bore moves relative to the pin.

Usage: python plot_vc_vs_gcvc.py --work VC_CAMPAIGN --ref B0_CAMPAIGN
(run directories in summary.json are remote paths; they are mapped onto --work/--ref)
"""
import argparse
import json
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.patches import Ellipse

import common
import part_b
from plotstyle import CM, CMAP, INK, MUTED, PRL, figure, note, save, setup
import matplotlib.pyplot as plt

X0 = 240.0   # angle window [-240, 120) deg keeps the load zone off the edge
REMOTE_REF = "/home/user/pancake-runs/lge-compressor-transient-20261119-02a4b6c"
REMOTE_WORK = "/home/user/pancake-runs/lge-compressor-vc-steady-20261123-00e5d1a"


def unwrapped(run_dir):
    fields = common.read_fields(Path(run_dir))
    nt, nz, theta, z, shape = part_b.grid_of(fields)
    from_hole = (theta - part_b.HOLE_THETA + X0) % 360.0 - X0   # window [-X0, 360 - X0)
    order = np.argsort(from_hole)
    p = shape("p")[:, order] / 1e5
    fc = shape("film_content")[:, order]
    h = shape("h")[:, order]
    hmin_angle = float(from_hole[order][np.argmin(h.min(axis=0))])
    return from_hole[order], z * 1e3, p, fc, hmin_angle


def local(path, remote, root):
    return Path(str(path).replace(remote, str(root)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--ref", required=True, type=Path)
    args = parser.parse_args()
    summary = json.loads((args.work / "summary.json").read_text())
    b0 = json.loads((args.ref / "b0" / "summary.json").read_text())
    gc_dir = local(b0["mean"]["accepted"]["run_dir"], REMOTE_REF, args.ref)
    vc_dir = local(summary["cases"]["V-mean"]["result"]["run_dir"], REMOTE_WORK, args.work)
    models = [("GC+VC", unwrapped(gc_dir)), ("VC only", unwrapped(vc_dir))]
    out = args.work / "figures"
    setup()

    # 2 x 2 maps, same colour scale per column.
    p_max = max(m[2].max() for _, m in models)
    left, bottom, gap_x, gap_y, w, h = 1.25, 1.0, 2.1, 0.75, 6.0, 3.0
    W, H = (left + 2 * w + gap_x + 1.5) * CM, (bottom + 2 * h + gap_y + 0.25) * CM
    fig = plt.figure(figsize=(W, H))
    hole_w = np.degrees(part_b.HOLE_RADIUS / part_b.R) * 2
    for row, (name, (angle, z, p, fc, hmin)) in enumerate(models):
        y0 = bottom + (1 - row) * (h + gap_y)
        for col, (data, lo, hi, label) in enumerate(((p, 0.0, p_max, r"Pressure, $p$ ($\mathrm{bar}$)"),
                                                     (fc, 0.0, 1.0, r"Film content, $\theta$"))):
            x0 = left + col * (w + gap_x)
            ax = fig.add_axes([x0 * CM / W, y0 * CM / H, w * CM / W, h * CM / H])
            image = ax.pcolormesh(angle, z, data, cmap=CMAP if col == 0 else CMAP.reversed(), vmin=lo, vmax=hi,
                                  shading="nearest", rasterized=True)
            ax.add_patch(Ellipse((0.0, part_b.HOLE_Z * 1e3), hole_w, part_b.HOLE_RADIUS * 2e3, fill=False,
                                 ec="white" if col == 0 else INK, lw=0.6))
            ax.axvline(hmin, color="white" if col == 0 else INK, lw=0.6, ls="--")
            ax.set_xlim(-X0, 360 - X0)
            ax.set_xticks(range(-240, 121, 60))
            ax.set_ylim(0, 20)
            ax.text(0.98 if col == 0 else 0.02, 0.90, name, transform=ax.transAxes, fontsize=6.5,
                    ha="right" if col == 0 else "left", color=INK)
            if row == 1:
                ax.set_xlabel(r"Angle from hole, $\theta-\theta_\mathrm{hole}$ ($°$)")
            else:
                ax.tick_params(labelbottom=False)
                ax.text(hmin + 4, 2.0, r"$h_\mathrm{min}$", fontsize=6.5, color="white" if col == 0 else INK)
            ax.set_ylabel(r"Axial position, $z$ ($\mathrm{mm}$)" if col == 0 else "")
            if row == 0:
                cax = fig.add_axes([(x0 + w + 0.15) * CM / W, (bottom) * CM / H, 0.2 * CM / W,
                                    (2 * h + gap_y) * CM / H])
                bar = fig.colorbar(image, cax=cax)
                bar.set_label(label)
                bar.outline.set_linewidth(0.4)
    save(fig, str(out / "os_bearing_130hz_vc_vs_gcvc_mean"))

    # Centre-line pressure, full range and zoom.
    j = None
    for zoom in (False, True):
        fig, ax = figure(w=13.0)
        for k, (name, (angle, z, p, fc, hmin)) in enumerate(models):
            j = int(np.argmin(np.abs(z - part_b.HOLE_Z * 1e3)))
            ax.plot(angle, p[j], color=PRL[k], lw=1.0, label=name)
        ax.axhline(part_b.P_D / 1e5, color=MUTED, lw=0.6, ls="--")
        ax.axvline(0.0, color=MUTED, lw=0.6, ls=":")
        ax.set_xlim(-X0, 360 - X0)
        ax.set_xticks(range(-240, 121, 30))
        ax.set_xlabel(r"Angle from hole, $\theta-\theta_\mathrm{hole}$ ($°$)")
        ax.set_ylabel(r"Centre-line pressure, $p$ ($\mathrm{bar}$)")
        hmin = models[0][1][4]
        if zoom:
            ax.set_ylim(0, 30)
            ax.legend(loc="lower right")
            ax.text(-236, part_b.P_D / 1e5 + 0.6, r"shell $=p_\mathrm{sat}=24.5\ \mathrm{bar}$", fontsize=6.5, color=MUTED)
            ax.text(2, 1.2, "hole", fontsize=6.5, color=MUTED)
            vc = models[1][1][2][j]
            low = models[1][1][0][np.argmin(vc)]
            note(ax, r"VC only: vapour cavity, $p=p_\mathrm{cav}=500\ \mathrm{Pa}$", xy=(low, 0.4), xytext=(low + 50, 6))
        else:
            ax.legend(loc="upper left")
            ax.text(2, 60, "hole", fontsize=6.5, color=MUTED)
            note(ax, r"$h_\mathrm{min}$", xy=(hmin, 5), xytext=(hmin + 25, 80))
        save(fig, str(out / ("os_bearing_130hz_centreline_pressure_zoom" if zoom else "os_bearing_130hz_centreline_pressure")))
    print(out)


if __name__ == "__main__":
    main()
