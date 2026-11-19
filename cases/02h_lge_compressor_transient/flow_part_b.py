"""Film mass flux (flow field) of one B2 frame: full film and a zoom around the feed hole.

q = rho [ U h theta / 2 - h^3 / (12 mu) grad p ]  per unit width (U = R omega, one moving surface);
in cavitated cells p = p_cav so only the Couette part remains. Streamlines over the void fraction,
with the void = 0.5 contour and pressure isobars near the hole. Also prints the 1D reformation
estimate of the liquid ring upstream of the hole, d = dp h^2 / (6 mu U (1 - theta_c)).

Usage: python flow_part_b.py --work CAMPAIGN_DIR [--crank 105]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

import animate_part_b as anim  # noqa: E402
import part_b  # noqa: E402

INK, INK_2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--crank", type=float, default=105.0)
    args = parser.parse_args()
    params = json.loads((args.work / "parameters.json").read_text())
    frames = [(t, f) for t, f in part_b.frame_list(args.work / "b2" / "last_fields") if t > 0]
    t, pvts = min(frames, key=lambda item: abs(item[0] / part_b.PERIOD * 360 - args.crank))
    fields = part_b.read_frame(pvts)
    theta, z_mm, p_bar, void, h_um = anim.grids(fields)
    mu, rho, U = params["mu_Pa_s"], params["rho_kg_m3"], part_b.R * part_b.OMEGA
    x = np.radians(theta) * part_b.R
    z = z_mm * 1e-3
    p, h, th = p_bar * 1e5, h_um * 1e-6, 1 - void
    dpdz, dpdx = np.gradient(p, z, x)
    qx = rho * (U * h * th / 2 - h**3 / (12 * mu) * dpdx)
    qz = rho * (-h**3 / (12 * mu) * dpdz)

    # 1D estimate of the ring upstream of the hole (mid-plane).
    j = int(np.argmin(np.abs(z - part_b.HOLE_Z)))
    i_hole = int(np.argmin(np.abs(theta - part_b.HOLE_THETA)))
    edge = part_b.HOLE_THETA - np.degrees(part_b.HOLE_RADIUS / part_b.R)
    upstream = np.where((theta < edge) & (void[j] >= 0.5))[0]
    i_cav = upstream.max()
    theta_c = th[j, i_cav]
    dp = params["p_feed_Pa"] - p[j, i_cav]
    d_est = dp * h[j, i_hole]**2 / (6 * mu * U * (1 - theta_c))
    d_sim = (np.radians(edge) - np.radians(theta[i_cav])) * part_b.R
    print(f"crank {t / part_b.PERIOD * 360:.1f} deg: dp = {dp / 1e5:.3f} bar, h_hole = {h[j, i_hole] * 1e6:.1f} um, "
          f"theta_c = {theta_c:.3f}, ring estimate {d_est * 1e3:.3f} mm, simulated {d_sim * 1e3:.3f} mm (cell edge)")

    plt.rcParams.update({"font.size": 8, "figure.facecolor": SURFACE, "axes.edgecolor": INK_2, "axes.labelcolor": INK_2})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), dpi=200, gridspec_kw=dict(width_ratios=[1.6, 1]))
    hole_w = np.degrees(part_b.HOLE_RADIUS / part_b.R) * 2
    speed = np.hypot(qx, qz) / rho * 1e6          # volume flux per width [mm^2/s]
    for ax, (t0, t1, z0, z1), density in ((axes[0], (0, 360, 0, 20), 1.6), (axes[1], (75, 105, 6.5, 13.5), 2.2)):
        sel_t = (theta >= t0) & (theta <= t1)
        sel_z = (z_mm >= z0) & (z_mm <= z1)
        tt, zz = theta[sel_t], z_mm[sel_z]
        sub = np.ix_(sel_z, sel_t)
        image = ax.pcolormesh(tt, zz, void[sub], cmap=anim.ORANGE, vmin=0, vmax=1, shading="nearest")
        # Streamlines in plotting units: theta [deg] per x, z [mm].
        u = qx[sub] / rho * (180 / np.pi / part_b.R)
        w = qz[sub] / rho * 1e3
        ax.streamplot(tt, zz, u, w, color="#1c5cab", density=density, linewidth=0.6, arrowsize=0.6)
        ax.contour(tt, zz, void[sub], levels=[0.5], colors=INK, linewidths=0.8)
        if ax is axes[1]:
            levels = np.linspace(params["p_shell_Pa"] / 1e5 + 0.02, params["p_feed_Pa"] / 1e5, 6)
            cs = ax.contour(tt, zz, p_bar[sub], levels=levels, colors="#52514e", linewidths=0.5, linestyles="dotted")
            ax.clabel(cs, fmt="%.2f", fontsize=6)
        ax.add_patch(Ellipse((part_b.HOLE_THETA, part_b.HOLE_Z * 1e3), hole_w, part_b.HOLE_RADIUS * 2e3,
                             fill=False, ec=INK, lw=1.2))
        ax.set_xlim(t0, t1)
        ax.set_ylim(z0, z1)
        ax.set_xlabel("θ' [deg] (OS bore moves towards +θ')")
        ax.set_ylabel("z [mm]")
    bar = fig.colorbar(image, ax=axes, pad=0.01, fraction=0.03)
    bar.set_label("Void fraction [-]", color=INK_2)
    axes[0].set_title("Film mass flux streamlines over void fraction (black: void = 0.5)", color=INK, fontsize=9)
    axes[1].set_title(f"Around the hole; dotted isobars [bar], feed {params['p_feed_Pa'] / 1e5:.2f} bar", color=INK, fontsize=9)
    fig.suptitle(f"OS bearing, 130 Hz, GC+VC, crank {t / part_b.PERIOD * 360:.0f} deg", color=INK, fontsize=10)
    out = args.work / f"os_bearing_130hz_flow_crank{args.crank:.0f}.png"
    fig.savefig(out, bbox_inches="tight")
    print(out, f"max |q|/rho {speed.max():.1f} mm^2/s")


if __name__ == "__main__":
    main()
