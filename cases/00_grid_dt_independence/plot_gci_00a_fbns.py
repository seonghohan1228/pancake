"""GCT plot for the FBNS 00a ladder (grid_G*_fbns)."""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, ".")
from plot_gci_00a import case_metrics
from plot_style import apply_style, save_figure, MARKERS, COLORS

GRIDS = {"grid_G1_fbns": 60, "grid_G2_fbns": 120, "grid_G3_fbns": 240,
         "grid_G4_fbns": 360, "grid_G5_fbns": 480, "grid_G6_fbns": 720,
         "grid_G7_fbns": 1080, "grid_G8_fbns": 1620}

results = {c: case_metrics(c) for c in GRIDS}
n_cells = np.array([n * (n * 16 // 60) for n in GRIDS.values()])
Ws = np.array([results[c][0] for c in GRIDS]) / 1e3
Ps = np.array([results[c][1] for c in GRIDS]) / 1e6


def convergence_ylim(ax, y):
    span = (y.max() - y.min()) / 0.5
    ax.set_ylim(y.min() - 0.4 * span, y.max() + 0.1 * span)


apply_style()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.0))
ax1.plot(n_cells, Ws, MARKERS[0], color=COLORS[0], ms=6)
ax1.set_xlabel(r"$N$")
ax1.set_ylabel(r"$W$ (kN)")
convergence_ylim(ax1, Ws)
ax2.plot(n_cells, Ps, MARKERS[0], color=COLORS[0], ms=6)
ax2.set_xlabel(r"$N$")
ax2.set_ylabel(r"$p_{\mathrm{max}}$ (MPa)")
convergence_ylim(ax2, Ps)
fig.tight_layout()
save_figure(fig, "gci_00a_fbns")
