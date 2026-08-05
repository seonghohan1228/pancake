"""Case 02 deliverable figure: eccentricity vs Sommerfeld number (log x) and
attitude angle vs eccentricity, for L/D = 1 and L/D = 0.5, from the flooded
(1 deg ambient groove at h_max) locus runs.
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
from plot_style import apply_style, save_figure, MARKERS, clean_log_ticks

SETS = [("locus_LD1.csv", r"$L/D$ = 1", MARKERS[0]),
        ("locus_LD05.csv", r"$L/D$ = 0.5", MARKERS[1])]


def read(name):
    with open(CASE_DIR / name, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main():
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    ax_s, ax_phi = axes

    for k, (fname, lab, mk) in enumerate(SETS):
        rows = read(fname)
        eps = [float(r["eps"]) for r in rows]
        somm = [float(r["S"]) for r in rows]
        phi = [float(r["attitude_deg"]) for r in rows]
        ax_s.plot(somm, eps, mk + "-", mfc="none", label=lab)
        ax_phi.plot(eps, phi, mk + "-", mfc="none", label=lab)

    ax_s.set_xscale("log")
    ax_s.set_xlabel(r"$S$")
    ax_s.set_ylabel(r"$\varepsilon$")
    ax_s.set_ylim(0, 1)
    clean_log_ticks(ax_s)
    ax_s.legend(loc="upper right")

    ax_phi.set_xlabel(r"$\varepsilon$")
    ax_phi.set_ylabel(r"$\phi$ (deg)")
    ax_phi.set_xlim(0, 1)
    ax_phi.set_ylim(0, 90)
    ax_phi.legend(loc="upper right")

    fig.subplots_adjust(wspace=0.32)
    save_figure(fig, str(CASE_DIR / "static_locus"))
    print("wrote static_locus.png")


if __name__ == "__main__":
    main()
