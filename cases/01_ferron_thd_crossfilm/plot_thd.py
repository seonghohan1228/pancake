"""Case 01-THD-crossfilm deliverable figures.

1. ferron_crossfilm_temperature.png - mid-plane BUSH-SIDE FILM SURFACE
   temperature for the three wall variants, vs Ferron 1983 Fig. 4 (2000 rpm,
   4 kN). The cross-film model resolves T(y), so the bush-side surface value
   exists for the first time; it is the closest model quantity to Ferron's
   thermocouples, which sit 0.5 mm BELOW the bush inner surface. The
   depth-averaged film temperature of the reference case is drawn as a thin grey
   line for scale, not as a competing prediction.
2. ferron_crossfilm_pressure.png    - mid-plane p(theta), same three variants,
   vs the same figure's pressure data.

Mapping used for the Fig. 4 experiment (stated, not tuned, identical to case
01-THD): the CSV abscissa is degrees from the supply-groove CENTRE measured in
the direction of rotation; the solver groove centre sits at theta = 90 deg and
rotation is +theta, so solver_theta = (90 + x) mod 360 and both curves are drawn
against (theta - 90) mod 360. No extra shift is applied and nothing is fitted.

Usage: python plot_thd.py [pressure|temperature|all]
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
DEPTH_AVG_CASE = CASE_DIR.parent / "01_ferron_thd"
GROOVE_DEG = 90.0
P_AMB = 101325.0

# The three wall variants are drawn with equal weight; none is preferred. The
# line styles only separate B and C, which coincide in the bush-side figure
# because a fixed bush wall makes that quantity the boundary condition itself.
VARIANTS = [
    ("walls_adiabatic",   "both walls adiabatic",
     dict(lw=1.8, ls="-")),
    ("bush40_Jadiabatic", r"bush 40 $^{\circ}$C, journal adiabatic",
     dict(lw=1.8, ls="--")),
    ("bush40_J40",        r"bush 40 $^{\circ}$C, journal 40 $^{\circ}$C",
     dict(lw=1.8, ls=":")),
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


def equilibria():
    return json.loads((CASE_DIR / "equilibrium.json").read_text())


def _profile(run_dir, field):
    theta, v = S.midplane_profile(run_dir, field)
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def _depth_averaged_reference():
    """Mid-plane depth-averaged T of the cross_film = NONE J-adiabatic run."""
    eq_path = DEPTH_AVG_CASE / "equilibrium.json"
    if not eq_path.exists():
        return None
    eq = json.loads(eq_path.read_text())
    rec = eq.get("fig4_Jadiabatic")
    if rec is None:
        return None
    return _profile(DEPTH_AVG_CASE / rec["flow"]["dir"], "T")


def plot_temperature(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ref = _depth_averaged_reference()
    if ref is not None:
        ax.plot(ref[0], [t - 273.15 for t in ref[1]], "-", lw=1.0, color="0.65",
                zorder=0, label="depth-averaged, no cross-film (case 01-THD)")
    for group, lab, style in VARIANTS:
        if group not in eq:
            continue
        x, t = _profile(CASE_DIR / eq[group]["flow"]["dir"], "cross_film_t_bush")
        ax.plot(x, [tv - 273.15 for tv in t], label=lab, **style)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    ax.plot(xr, yr, MARKERS[2], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$T_{\mathrm{bush\ surface}}$ ($^{\circ}$C)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=11)
    save_figure(fig, str(CASE_DIR / "ferron_crossfilm_temperature"))
    print("wrote ferron_crossfilm_temperature.png")


def plot_pressure(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    for group, lab, style in VARIANTS:
        if group not in eq:
            continue
        x, p = _profile(CASE_DIR / eq[group]["flow"]["dir"], "p")
        ax.plot(x, [(pv - P_AMB) / 1.0e6 for pv in p], label=lab, **style)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv", 1.0e-6)
    ax.plot(xr, yr, MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=11)
    save_figure(fig, str(CASE_DIR / "ferron_crossfilm_pressure"))
    print("wrote ferron_crossfilm_pressure.png")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    apply_style()
    eq = equilibria()
    if which in ("temperature", "all"):
        plot_temperature(eq)
    if which in ("pressure", "all"):
        plot_pressure(eq)


if __name__ == "__main__":
    main()
