"""Case 01-CONJUGATE deliverable figures.

1. ferron_conjugate_temperature.png - mid-plane BUSH-SIDE FILM SURFACE
   temperature (`T_bearing`) for the two outer-boundary variants of the 2-D
   conjugate bush, against the LOCAL 1-D bush closure (variant A of case
   01-THD-crossfilm, same model, same operating point) and Ferron 1983 Fig. 4
   (2000 rpm, 4 kN).
2. ferron_conjugate_pressure.png - mid-plane p(theta), same variants, vs the same
   figure's pressure data.

Mapping used for the Fig. 4 experiment (stated, not tuned, identical to cases
01-THD and 01-THD-crossfilm): the CSV abscissa is degrees from the supply-groove
CENTRE measured in the direction of rotation; the solver groove centre sits at
theta = 90 deg and rotation is +theta, so solver_theta = (90 + x) mod 360 and both
curves are drawn against (theta - 90) mod 360. No extra shift, nothing fitted.

Usage: python plot_conjugate.py [pressure|temperature|all]
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
LOCAL_CASE = CASE_DIR.parent / "01_ferron_thd_crossfilm"
GROOVE_DEG = 90.0
P_AMB = 101325.0

VARIANTS = [
    ("outer_fixed313", r"2-D bush, outer 40 $^{\circ}$C", dict(lw=1.8, ls="-")),
    ("outer_robin80", r"2-D bush, outer $h=80$ W m$^{-2}$K$^{-1}$",
     dict(lw=1.8, ls="--")),
]

# outer_robin80 did not close its force balance: the temperature KSP diverged at
# the third Newton step (pre-existing DIVERGED_DTOL defect, reproduced on BOTH
# hosts at different iterates, i.e. decomposition-dependent). Its best converged
# iterate is used instead and is labelled as such -- W = 4001.36 N (0.034 % from
# target) with the force angle 0.062 deg out against a 0.05 deg tolerance, so it
# is a near-equilibrium point, not a closed one.
FALLBACK_DIR = {"outer_robin80": "outer_robin80/fb_iter_03"}
FALLBACK_NOTE = {"outer_robin80": " (near-eq.)"}


def run_dir_for(eq, group):
    """Converged flow-pass directory, or the documented fallback iterate."""
    if group in eq:
        return CASE_DIR / eq[group]["flow"]["dir"], ""
    if group in FALLBACK_DIR:
        return CASE_DIR / FALLBACK_DIR[group], FALLBACK_NOTE[group]
    return None, ""


def read_csv(path, yscale=1.0):
    xs, ys = [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].lstrip().startswith("#"):
                continue
            xs.append(float(row[0]))
            ys.append(float(row[1]) * yscale)
    return xs, ys


def _profile(run_dir, field):
    got = S.midplane_profile(run_dir, field)
    if got is None:
        return None
    theta, v = got
    x = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in order], [v[i] for i in order]


def _local_bush(field):
    """Variant A of the LOCAL-bush case: same model, per-cell 1-D radial bush."""
    eq_path = LOCAL_CASE / "equilibrium.json"
    if not eq_path.exists():
        return None
    eq = json.loads(eq_path.read_text())
    rec = eq.get("walls_adiabatic")
    if rec is None:
        return None
    return _profile(LOCAL_CASE / rec["flow"]["dir"], field)


def plot_temperature(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ref = _local_bush("cross_film_t_bush")
    if ref is not None:
        ax.plot(ref[0], [t - 273.15 for t in ref[1]], "-", lw=1.1, color="0.6",
                zorder=0, label="LOCAL 1-D bush, adiabatic (case 01-THD-crossfilm)")
    for group, lab, style in VARIANTS:
        d, note = run_dir_for(eq, group)
        if d is None:
            continue
        got = _profile(d, "T_bearing")
        if got is None:
            continue
        ax.plot(got[0], [tv - 273.15 for tv in got[1]], label=lab + note, **style)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_temperature.csv")
    ax.plot(xr, yr, MARKERS[2], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$T_{\mathrm{bush\ surface}}$ ($^{\circ}$C)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=10)
    save_figure(fig, str(CASE_DIR / "ferron_conjugate_temperature"))
    print("wrote ferron_conjugate_temperature.png")


def plot_pressure(eq):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ref = _local_bush("p")
    if ref is not None:
        ax.plot(ref[0], [(pv - P_AMB) / 1.0e6 for pv in ref[1]], "-", lw=1.1,
                color="0.6", zorder=0,
                label="LOCAL 1-D bush, adiabatic (case 01-THD-crossfilm)")
    for group, lab, style in VARIANTS:
        d, note = run_dir_for(eq, group)
        if d is None:
            continue
        got = _profile(d, "p")
        if got is None:
            continue
        ax.plot(got[0], [(pv - P_AMB) / 1.0e6 for pv in got[1]], label=lab + note, **style)
    xr, yr = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv", 1.0e-6)
    ax.plot(xr, yr, MARKERS[1], mfc="none", mec="k", ls="none",
            label="Ferron et al. (1983), Fig. 4")
    ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
    ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4, fontsize=10)
    save_figure(fig, str(CASE_DIR / "ferron_conjugate_pressure"))
    print("wrote ferron_conjugate_pressure.png")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    apply_style()
    eq = json.loads((CASE_DIR / "equilibrium.json").read_text())
    if which in ("temperature", "all"):
        plot_temperature(eq)
    if which in ("pressure", "all"):
        plot_pressure(eq)


if __name__ == "__main__":
    main()
