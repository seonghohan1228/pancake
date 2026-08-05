"""Case 01-VS-REYNSIM deliverable figures.

1. reynsim_midplane_pressure.png - mid-plane p(theta) at 4000 rpm / 6 kN:
   Ferron 1983 Fig. 6 experiment (points), previous-generation MATLAB/ReynSim
   (line), current pancake FBNS isothermal band (two thin lines).
2. reynsim_eccentricity.png - eps vs W: ReynSim's two numerical vintages, their
   own experimental digitisation (source unverified), our Fig. 5 digitisation
   rescaled 145/152, and our FBNS+THD sweep at 2000 rpm.

Frame conventions, all stated and none fitted:
  * abscissa = degrees from the supply-groove CENTRE in the direction of
    rotation.  Our mesh puts the groove at theta = 90 deg, so solver curves are
    drawn against (theta - 90) % 360, as in 01_ferron_prelim / 01_ferron_thd.
  * ReynSim's own frame is ours + 90 deg; see reynsim_io.py for the derivation.
  * ordinate = gauge pressure, p - 101325 Pa, for every curve.

Usage: python plot_vs_reynsim.py [pressure|eccentricity|all]
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S
import reynsim_io as R
from plot_style import apply_style, save_figure, MARKERS, COLORS

CASES = CASE_DIR.parent
REF = CASES / "ref_data"
PRELIM = CASES / "01_ferron_prelim"
THD = CASES / "01_ferron_thd"
GROOVE_DEG = 90.0
P_AMB = 101325.0
ECC_RESCALE = 145.0 / 152.0


def read_csv(path, yscale=1.0):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]) * yscale)
    return xs, ys


def solver_midplane(run_dir):
    theta, p = S.midplane_profile(run_dir, "p")
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [(p[i] - P_AMB) / 1.0e6 for i in order]


def plot_pressure():
    eq = json.loads((PRELIM / "equilibrium.json").read_text())
    fig, ax = plt.subplots(figsize=(6.6, 4.6))

    # Current solver: 01-prelim isothermal FBNS band at this exact condition.
    # One legend entry for the current solver: the isothermal band (40-55 C).
    band = [("mu40", "-"), ("mu55", "--")]
    curves = []
    for i, (key, ls) in enumerate(band):
        x, y = solver_midplane(PRELIM / eq[key]["dir"])
        curves.append((x, y))
        ax.plot(x, y, ls, lw=1.1, color=COLORS[2],
                label=("This year (JFO)" if i == 0 else None))
    if len(curves[0][0]) == len(curves[1][0]):
        ax.fill_between(curves[0][0], curves[0][1], curves[1][1],
                        color=COLORS[2], alpha=0.15, lw=0, zorder=0)

    # Previous generation.
    xr, pr = R.midplane_pressure()
    ax.plot(xr, [v / 1.0e6 for v in pr], "-", lw=1.8, color=COLORS[1],
            label="Last year")

    # Experiment.
    xe, ye = read_csv(REF / "ferron1983_fig6_4000rpm_6kN_pressure.csv", 1.0e-6)
    ax.plot(xe, ye, MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron (1983)")

    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.set_ylim(-0.2, 2.4)
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=11)
    save_figure(fig, str(CASE_DIR / "reynsim_midplane_pressure"))
    print("wrote reynsim_midplane_pressure.png")


def our_thd_sweep():
    eq = json.loads((THD / "equilibrium.json").read_text())
    pts = []
    for group, rec in eq.items():
        if group == "fig4_Jadiabatic":
            w = 4000.0
        elif group.startswith("fig5_W"):
            w = float(group[len("fig5_W"):])
        else:
            continue
        pts.append((w / 1000.0, rec["flow"]["eps"]))
    pts.sort()
    return [p[0] for p in pts], [p[1] for p in pts]


def plot_eccentricity():
    fig, ax = plt.subplots(figsize=(6.4, 4.6))

    # Experiments first (open black symbols).
    xe, ye = read_csv(REF / "ferron1983_fig5_2000rpm_eccentricity.csv")
    ax.plot([x / 1000.0 for x in xe], [y * ECC_RESCALE for y in ye],
            MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron Fig. 5, 2000 rpm ($C$ = 145 $\\mu$m $\\times$ 145/152)")
    ax.plot(R.VALIDATION_M["load_exp_kN"], R.VALIDATION_M["eps_exp"],
            MARKERS[2], mfc="none", mec="0.45", ls="none",
            label="experiment (ReynSim digitisation, source unverified)")

    # Previous code, two vintages.
    ax.plot([w / 1000.0 for w in R.VALIDATION_M["load_N"]], R.VALIDATION_M["eps_num"],
            "-", marker=MARKERS[3], mfc="none", color=COLORS[1], lw=1.6,
            label="ReynSim numerical (validation.m, speed unverified)")
    ax.plot([w / 1000.0 for w in R.PLOT_DATA_M["load_N"]], R.PLOT_DATA_M["eps_num"],
            "--", marker=MARKERS[4], mfc="none", color=COLORS[4], lw=1.4,
            label="ReynSim numerical (plot_data.m, earlier vintage)")

    # Current solver.
    xs, ys = our_thd_sweep()
    ax.plot(xs, ys, "-", marker=MARKERS[0], color=COLORS[0], lw=1.8,
            label="pancake FBNS + THD, 2000 rpm")

    ax.set_xlabel(r"$W$ (kN)")
    ax.set_ylabel(r"$\varepsilon = e/c$")
    ax.set_xlim(1, 11)
    ax.set_ylim(0.2, 0.95)
    ax.legend(loc="lower right", handlelength=1.8, borderaxespad=0.4, fontsize=10)
    save_figure(fig, str(CASE_DIR / "reynsim_eccentricity"))
    print("wrote reynsim_eccentricity.png")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    apply_style()
    if which in ("pressure", "all"):
        plot_pressure()
    if which in ("eccentricity", "all"):
        plot_eccentricity()


if __name__ == "__main__":
    main()
