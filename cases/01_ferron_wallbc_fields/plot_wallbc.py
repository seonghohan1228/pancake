"""Case 01-WALLBC-FIELDS deliverable figures.

1. ferron_wallbc_temperature.png - mid-plane wall-side temperature for the four
   wall treatments vs Ferron 1983 Fig. 4 (2000 rpm, 4 kN). For the three
   cross-film groups the plotted quantity is T_bearing, the BUSH-SIDE FILM
   SURFACE temperature, which is the closest model quantity to Ferron's
   thermocouples 0.5 mm below the bush inner surface. favg_walls40 has no
   cross-film profile, so only the depth average T_bar exists; it is drawn
   dashed and LABELLED as a different quantity rather than silently compared.

2. ferron_wallbc_heatflux.png - mid-plane wall heat flux, positive INTO the
   fluid. For the cross-film groups these are the solver's own q_bearing /
   q_journal fields. favg_walls40 cannot write them (they are undefined without
   a cross-film profile and requesting them is a config error), so they are
   post-derived from its own Newton-cooling law, q = h (T_wall - T_bar), with
   h = 7200 (bearing) and 10000 (journal) W/(m^2 K).

Mapping for the experiment (stated, not tuned, identical to cases 01-THD,
01-THD-crossfilm and 01-CONJUGATE): solver_theta = (90 + x) mod 360; both curves
are drawn against (theta - 90) mod 360. No fitted shift.

Usage: python plot_wallbc.py [temperature|heatflux|all]
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
from plot_style import apply_style, save_figure, MARKERS

REF = CASE_DIR.parent / "ref_data"

GROUPS = [
    ("rerun_xf_adiabatic", "cross-film, both walls adiabatic", dict(lw=1.8, ls="-")),
    ("rerun_xf_bush40_J40", r"cross-film, bush 40 $^{\circ}$C / journal 40 $^{\circ}$C",
     dict(lw=1.8, ls="-.")),
    ("rerun_cht_robin10", r"conjugate 2-D bush, outer $h=10$", dict(lw=2.0, ls="--")),
    ("favg_walls40", r"depth-averaged, walls 40 $^{\circ}$C  ($\bar{T}$, not a surface)",
     dict(lw=1.4, ls=":", color="0.45")),
]


def read_csv(path, yscale=1.0):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]) * yscale)
    return xs, ys


def data():
    return json.loads((CASE_DIR / "extract.json").read_text())


def plot_temperature(d):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for g, lab, style in GROUPS:
        if g not in d:
            continue
        p = d[g]["profile"]
        key = "T_bearing_C" if d[g]["has_wall_fields"] else "T_bar_C"
        ax.plot(p["x_deg"], p[key], label=lab, **style)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    ax.plot(xr, yr, MARKERS[2], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$T$ ($^{\circ}$C)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.legend(loc="upper left", handlelength=1.9, borderaxespad=0.4, fontsize=9)
    save_figure(fig, str(CASE_DIR / "ferron_wallbc_temperature"))
    print("wrote ferron_wallbc_temperature.png")


def plot_heatflux(d):
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
    for ax, side, title in ((axes[0], "q_bearing", "bush side"),
                            (axes[1], "q_journal", "shaft side")):
        for g, lab, style in GROUPS:
            if g not in d:
                continue
            p = d[g]["profile"]
            ax.plot(p["x_deg"], [v / 1e3 for v in p[side]], label=lab, **style)
        ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
        ax.set_ylabel(r"$q$ (kW m$^{-2}$)")
        ax.text(0.985, 0.06, title, transform=ax.transAxes, ha="right", fontsize=11)
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    axes[1].set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    axes[0].legend(loc="lower left", handlelength=1.9, borderaxespad=0.4, fontsize=8.5)
    save_figure(fig, str(CASE_DIR / "ferron_wallbc_heatflux"))
    print("wrote ferron_wallbc_heatflux.png")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    apply_style()
    d = data()
    if which in ("temperature", "all"):
        plot_temperature(d)
    if which in ("heatflux", "all"):
        plot_heatflux(d)


if __name__ == "__main__":
    main()
