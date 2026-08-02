"""Report figures for the 2026-08-02 verification campaign.

Produces the two figures that carry the findings not already covered by the
existing `validation_*.png` set:

  validation_cross_solver.png   -- pancake against an independent solver at
                                   Ferron's geometry, the evidence that the load
                                   over-prediction is not a discretisation error.
  validation_load_in_context.png -- the same over-prediction shown against the
                                   published models, which share it.

Usage:  $PY make_report_figures.py
"""
import json, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ferron as F
from plot_style import apply_style, save_figure, MARKERS, COLORS

HERE = os.path.dirname(os.path.abspath(__file__))

# Secondhand single points, digitised by a third party from figures in the open
# reproductions. Plotted as isolated markers, never joined into a curve, because
# that is all the provenance supports.
#   (load N, eccentricity ratio, label, source)
LITERATURE_POINTS = [
    (7200.0, 0.67, "Ferron's own theory",
     "Lubricants 11(7) 308 (2023) Fig. 8"),
    (4000.0, 0.56, "Alakhramsing FE model",
     "Lubricants 3(2) 256 (2015) Fig. 10a"),
]


def load_cross():
    path = os.path.join(HERE, "cross_solver_results.json")
    if not os.path.exists(path):
        raise SystemExit("run_cross_solver.py must be run first")
    data = json.load(open(path))
    return data["cases"] if isinstance(data, dict) else data


def figure_cross_solver(cases):
    """Load against eccentricity: two independent solvers and the measurement."""
    apply_style()
    fig, ax = plt.subplots(figsize=(4.8, 4.3))

    eps = [c["eps"] for c in cases]
    w_pancake = [c["pancake"]["W_N"] / 1e3 for c in cases]
    w_indep = [c["independent"]["W_N"] / 1e3 for c in cases]

    meas_eps = [e for _, e in F.FERRON_ECC_VS_LOAD]
    meas_w = [w / 1e3 for w, _ in F.FERRON_ECC_VS_LOAD]

    ax.errorbar(meas_eps, meas_w, xerr=F.ECC_DIGITISE_ERR, fmt=MARKERS[2],
                color=COLORS[6], markerfacecolor="white", markersize=6,
                linestyle="none", label="Ferron (1983), measured", zorder=3)
    ax.plot(eps, w_pancake, MARKERS[0] + "-", color=COLORS[0], markersize=6,
            label="pancake (Elrod-Adams)", zorder=5)
    ax.plot(eps, w_indep, MARKERS[1] + "--", color=COLORS[1], markersize=6,
            label="independent (Christopherson)", zorder=4)

    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$W$ (kN)")
    ax.set_xlim(0.40, 0.88)
    ax.set_ylim(0, 17)
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    save_figure(fig, "validation_cross_solver")
    plt.close(fig)


def figure_load_in_context(cases):
    """The over-prediction is shared with the published models."""
    apply_style()
    fig, ax = plt.subplots(figsize=(4.8, 4.3))

    meas_w = [w / 1e3 for w, _ in F.FERRON_ECC_VS_LOAD]
    meas_eps = [e for _, e in F.FERRON_ECC_VS_LOAD]
    ax.errorbar(meas_w, meas_eps, yerr=F.ECC_DIGITISE_ERR, fmt=MARKERS[2],
                color=COLORS[6], markerfacecolor="white", markersize=6,
                linestyle="none", label="Ferron (1983), measured", zorder=3)

    ax.plot([c["pancake"]["W_N"] / 1e3 for c in cases],
            [c["eps"] for c in cases], MARKERS[0] + "-", color=COLORS[0],
            markersize=6, label="pancake", zorder=5)

    for k, (w, e, label, _src) in enumerate(LITERATURE_POINTS):
        ax.plot(w / 1e3, e, MARKERS[3 + k], color=COLORS[2 + k], markersize=8,
                linestyle="none", label=label, zorder=4)

    ax.set_xlabel(r"$W$ (kN)")
    ax.set_ylabel(r"$\varepsilon$")
    ax.set_xlim(0, 17)
    ax.set_ylim(0.40, 0.90)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    save_figure(fig, "validation_load_in_context")
    plt.close(fig)


if __name__ == "__main__":
    cases = load_cross()
    figure_cross_solver(cases)
    figure_load_in_context(cases)
    print("wrote validation_cross_solver.png, validation_load_in_context.png")
