"""Case 01-THD deliverable figures.

1. ferron_thd_pressure.png     - mid-plane p(theta), both journal-wall variants,
                                 vs Ferron 1983 Fig. 4 (2000 rpm, 4 kN).
2. ferron_thd_temperature.png  - mid-plane film T(theta), both variants,
                                 vs the same figure's bush temperature.
3. ferron_thd_eccentricity.png - eps vs load at 2000 rpm (J-adiabatic) vs the
                                 Fig. 5 experiment rescaled by 145/152.

Mapping used for the Fig. 4 experiment (stated, not tuned): the CSV abscissa is
degrees from the supply-groove CENTRE measured in the direction of rotation; the
solver groove centre sits at theta = 90 deg and rotation is +theta, so
solver_theta = (90 + x) mod 360 and both curves are drawn against
(theta - 90) mod 360.  No extra shift is applied.

Usage: python plot_thd.py [pressure|temperature|eccentricity|all]
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

REF = CASE_DIR.parent / "ref_data"
GROOVE_DEG = 90.0
P_AMB = 101325.0
ECC_RESCALE = 145.0 / 152.0   # Fig. 5 experiment is on a C = 145 um basis

VARIANTS = [("fig4_Jadiabatic", "J-adiabatic"),
            ("fig4_Jiso40", r"J-iso40 ($h_j$ = 10$^4$ W m$^{-2}$K$^{-1}$)")]


def read_csv(path, yscale=1.0):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]) * yscale)
    return xs, ys


def equilibria():
    return json.loads((CASE_DIR / "equilibrium.json").read_text())


def _profile(run_dir, field):
    theta, v = S.midplane_profile(run_dir, field)
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def plot_pressure(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for group, lab in VARIANTS:
        x, p = _profile(CASE_DIR / eq[group]["flow"]["dir"], "p")
        ax.plot(x, [(pv - P_AMB) / 1.0e6 for pv in p], "-", lw=1.8, label=lab)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv", 1.0e-6)
    ax.plot(xr, yr, MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4)
    save_figure(fig, str(CASE_DIR / "ferron_thd_pressure"))
    print("wrote ferron_thd_pressure.png")


def plot_temperature(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for group, lab in VARIANTS:
        x, t = _profile(CASE_DIR / eq[group]["flow"]["dir"], "T")
        ax.plot(x, [tv - 273.15 for tv in t], "-", lw=1.8, label=lab)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    ax.plot(xr, yr, MARKERS[2], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$T$ ($^{\circ}$C)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4)
    save_figure(fig, str(CASE_DIR / "ferron_thd_temperature"))
    print("wrote ferron_thd_temperature.png")


def plot_eccentricity(eq):
    loads, eps = [], []
    for group, rec in eq.items():
        if group == "fig4_Jadiabatic":
            loads.append(4000.0)
        elif group.startswith("fig5_W"):
            loads.append(float(group[len("fig5_W"):]))
        else:
            continue
        eps.append(rec["flow"]["eps"])
    order = sorted(range(len(loads)), key=lambda i: loads[i])
    loads = [loads[i] / 1000.0 for i in order]
    eps = [eps[i] for i in order]

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    ax.plot(loads, eps, "-", marker=MARKERS[0], lw=1.8,
            label="pancake, FBNS + THD (J-adiabatic)")
    xr, yr = read_csv(REF / "ferron1983_fig5_2000rpm_eccentricity.csv")
    ax.plot([x / 1000.0 for x in xr], [y * ECC_RESCALE for y in yr],
            MARKERS[1], mfc="none", mec="k", ls="none",
            label=r"Ferron et al. (1983), Fig. 5" "\n"
                  r"($C$ = 145 $\mu$m data $\times$ 145/152)")
    ax.set_xlabel(r"$W$ (kN)")
    ax.set_ylabel(r"$\varepsilon = e/c$")
    ax.set_xlim(0, 11)
    ax.set_ylim(0.3, 1.0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=12)
    save_figure(fig, str(CASE_DIR / "ferron_thd_eccentricity"))
    print("wrote ferron_thd_eccentricity.png")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    apply_style()
    eq = equilibria()
    if which in ("pressure", "all"):
        plot_pressure(eq)
    if which in ("temperature", "all"):
        plot_temperature(eq)
    if which in ("eccentricity", "all"):
        plot_eccentricity(eq)


if __name__ == "__main__":
    main()
