"""Eccentricity-ratio history e(t)/c for the four dt cases (00b, FBNS, G6)."""
import sys
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, ".")
from plot_style import apply_style, save_figure, blues

C = 1.52e-4
CASES = [  # (label_dt_us, dir, results subdir)
    (40, "dt_T1_fbns", "results"),
    (20, "dt_T2_fbns", "results_remote"),
    (10, "dt_T3_fbns", "results_remote"),
    (5,  "dt_T4_fbns", "results_remote"),
]


def read_xb(vts_path):
    root = ET.parse(vts_path).getroot()
    for da in root.iter("DataArray"):
        if da.get("Name") == "xB":
            v = np.fromstring(da.text, sep=" ")
            return v[:3]
    raise KeyError(f"xB not in {vts_path}")


def history(case_dir, sub):
    flat = Path(case_dir) / sub / "flat"
    pvd = flat / "results.pvd"
    t, e = [], []
    for ds in ET.parse(pvd).getroot().iter("DataSet"):
        ts = float(ds.get("timestep"))
        if ts > 0.06 + 1e-12:
            continue
        pvts = flat / ds.get("file")
        piece = next(p.get("Source") for p in ET.parse(pvts).getroot().iter("Piece")
                     if p.get("Source"))
        xb = read_xb(flat / piece)
        t.append(ts)
        e.append(math.hypot(xb[0], xb[1]) / C)
    o = np.argsort(t)
    return np.array(t)[o], np.array(e)[o]


apply_style()
fig, ax = plt.subplots(figsize=(5.4, 4.2))
shades = blues(4)
for (dt_us, d, sub), col in zip(CASES, shades):
    t, e = history(d, sub)
    ax.plot(t * 1e3, e, "-", color=col, lw=1.4,
            label=rf"$\Delta t = {dt_us}\ \mu\mathrm{{s}}$")
ax.set_xlabel(r"$t$ (ms)")
ax.set_ylabel(r"$\varepsilon$")
ax.legend()
fig.tight_layout()
save_figure(fig, "e_vs_t_00b")
