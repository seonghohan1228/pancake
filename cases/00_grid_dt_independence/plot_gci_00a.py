"""Extract W and p_max from 00a flat VTK output, compute GCI (Celik 2008), plot."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from plot_style import apply_style, save_figure, MARKERS, COLORS, clean_log_ticks

CASE_DIR = Path(__file__).parent
GRIDS = {  # name -> n_theta
    "grid_G1": 60, "grid_G2": 120, "grid_G3": 240, "grid_G4": 360, "grid_G5": 480,
    "grid_G6": 720,
}
for name, n in (("grid_G7", 1080), ("grid_G8", 1620)):
    if list((Path(__file__).parent / name / "results" / "flat").glob("*.pvts")):
        GRIDS[name] = n


def read_cell_array(vts_path, name):
    root = ET.parse(vts_path).getroot()
    for da in root.iter("DataArray"):
        if da.get("Name") == name:
            vals = np.fromstring(da.text, sep=" ")
            ncomp = int(da.get("NumberOfComponents", "1"))
            return vals.reshape(-1, ncomp) if ncomp > 1 else vals
    raise KeyError(f"{name} not in {vts_path}")


def case_metrics(case):
    flat = CASE_DIR / case / "results" / "flat"
    pvts = sorted(flat.glob("*.pvts"))[-1]
    pieces = [flat / p.get("Source")
              for p in ET.parse(pvts).getroot().iter("Piece") if p.get("Source")]
    p_max = max(read_cell_array(pc, "p").max() for pc in pieces)
    fp = read_cell_array(pieces[0], "Fp")[0]  # constant field
    W = math.hypot(fp[0], fp[1])
    att = math.degrees(math.atan2(fp[0], -fp[1]))  # deg off the -y load line
    return W, p_max, att


def gci(f_coarse, f_mid, f_fine, r21, r32):
    """Celik et al. 2008; index 1 = fine. Returns observed order, extrapolated, GCI_fine."""
    e21, e32 = f_mid - f_fine, f_coarse - f_mid
    s = math.copysign(1.0, e32 / e21)
    p = 2.0
    for _ in range(50):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            break
        p = p_new
    f_ext = f_fine + (f_fine - f_mid) / (r21 ** p - 1.0)
    gci_fine = 1.25 * abs((f_mid - f_fine) / f_fine) / (r21 ** p - 1.0)
    return p, f_ext, gci_fine


results = {c: case_metrics(c) for c in GRIDS}
print(f"{'case':12s} {'n_theta':>7s} {'W [N]':>12s} {'p_max [MPa]':>12s} {'att [deg]':>9s} {'dW% vs coarser':>14s}")
prev_W = None
for c, n in GRIDS.items():
    W, pm, att = results[c]
    dW = "" if prev_W is None else f"{100 * (W - prev_W) / prev_W:+.3f}"
    print(f"{c:12s} {n:7d} {W:12.2f} {pm / 1e6:12.4f} {att:9.3f} {dW:>14s}")
    prev_W = W

triplet = list(GRIDS)[-3:]  # three finest available grids
n3, n4, n5 = (GRIDS[c] for c in triplet)
W3, W4, W5 = (results[c][0] for c in triplet)
P3, P4, P5 = (results[c][1] for c in triplet)
r21, r32 = n5 / n4, n4 / n3
print(f"GCI triplet: {triplet}")
pW, Wext, gciW = gci(W3, W4, W5, r21, r32)
pP, Pext, gciP = gci(P3, P4, P5, r21, r32)
print(f"\nW    : order p = {pW:.2f}, extrapolated = {Wext:.2f} N, GCI_fine = {100 * gciW:.3f}%")
print(f"p_max: order p = {pP:.2f}, extrapolated = {Pext / 1e6:.4f} MPa, GCI_fine = {100 * gciP:.3f}%")

apply_style()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 4.0))
n_cells = np.array([n * (n * 16 // 60) for n in GRIDS.values()])  # n_theta * n_z
Ws = np.array([results[c][0] for c in GRIDS]) / 1e3
Ps = np.array([results[c][1] for c in GRIDS]) / 1e6


def convergence_ylim(ax, y):
    """Data span occupies 40%-90% of the axis height (10% above max, 40% below min)."""
    span = (y.max() - y.min()) / 0.5
    ax.set_ylim(y.min() - 0.4 * span, y.max() + 0.1 * span)


ax1.plot(n_cells, Ws, MARKERS[0], color=COLORS[0], ms=6)
ax1.set_xlabel(r"$N$")
ax1.set_ylabel(r"$W$ (kN)")
convergence_ylim(ax1, Ws)

ax2.plot(n_cells, Ps, MARKERS[0], color=COLORS[0], ms=6)
ax2.set_xlabel(r"$N$")
ax2.set_ylabel(r"$p_{\mathrm{max}}$ (MPa)")
convergence_ylim(ax2, Ps)

fig.tight_layout()
save_figure(fig, str(CASE_DIR / "gci_00a"))
