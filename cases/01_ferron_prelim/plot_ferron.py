"""Case 01 deliverable figure: mid-plane p(theta) at the force-balanced
equilibrium for both viscosities, overlaid with Ferron et al. (1983) Fig. 6.

Mapping used for the experiment (stated, not tuned): the CSV abscissa is
degrees from the supply-groove CENTRE measured in the direction of rotation;
the solver groove centre sits at theta = 90 deg and rotation is +theta, so
solver_theta = (90 + x) mod 360, i.e. both curves are plotted against
(theta - 90) mod 360. No extra shift is applied.
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S
from plot_style import apply_style, save_figure, MARKERS

REF_CSV = CASE_DIR.parent / "ref_data" / "ferron1983_fig6_4000rpm_6kN_pressure.csv"
GROOVE_DEG = 90.0
P_AMB = 101325.0


def read_ref():
    xs, ys = [], []
    with open(REF_CSV, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]) / 1.0e6)
    return xs, ys


def main():
    apply_style()
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    labels = [("mu40", r"$\mu$ = 0.0277 Pa s (40 $^{\circ}$C)"),
              ("mu55", r"$\mu$ = 0.0166 Pa s (55 $^{\circ}$C)")]

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for k, (key, lab) in enumerate(labels):
        run_dir = CASE_DIR / eq[key]["dir"]
        theta, p = S.midplane_pressure(run_dir)
        x = [(t - GROOVE_DEG) % 360.0 for t in theta]
        y = [(pv - P_AMB) / 1.0e6 for pv in p]
        order = sorted(range(len(x)), key=lambda i: x[i])
        ax.plot([x[i] for i in order], [y[i] for i in order], "-",
                lw=1.8, label=lab)

    xr, yr = read_ref()
    ax.plot(xr, yr, MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 6")

    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_ylim(-0.15, 2.85)
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4)
    save_figure(fig, str(CASE_DIR / "ferron_midplane_pressure"))
    print("wrote ferron_midplane_pressure.png")


if __name__ == "__main__":
    main()
