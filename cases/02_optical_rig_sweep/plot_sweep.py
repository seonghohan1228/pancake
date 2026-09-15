"""Case-02 figures: the P1 pressure overlay and the operating-window maps.

phi is the case-02 angle convention (see README): phi = 0 is the line of
centres, phi increases with journal rotation.
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent
# Case 01 also ships a plot_style; keep this case's copy ahead of it on the path.
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
sys.path.insert(0, str(CASE_DIR))
import solver_io as S
from plot_style import apply_style, save_figure, MARKERS, COLORS, smart_legend

BAR = 1.0e5


def midplane(run_dir):
    phi, p = S.midplane_profile(CASE_DIR / run_dir, "p")
    return np.asarray(phi), np.asarray(p) / BAR


def case_name(eps, rpm, key):
    return f"r1_e{round(eps * 100):03d}_n{rpm:04d}_{key}"


def plot_p1_overlay(ledger, eps, rpm):
    """Mid-span p(phi): degassed against p_sat = 0.8 bar at the P1 point."""
    apply_style()
    fig, ax = plt.subplots()
    for i, (key, label) in enumerate([("degassed", "degassed"),
                                      ("p080", r"$p_{sat} = 0.8$ bar")]):
        entry = ledger.get(case_name(eps, rpm, key))
        if not entry or entry.get("rc") != 0:
            continue
        phi, p = midplane(entry["dir"])
        ax.plot(phi, p, marker=MARKERS[i], color=COLORS[i], label=label,
                markevery=max(1, len(phi) // 24))
        ax.axhline(entry["p_cav"] / BAR, color=COLORS[i], ls=":", lw=0.9)
    ax.set_xlabel(r"$\varphi$ (deg)")
    ax.set_ylabel(r"$p$ (bar abs)")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 60))
    smart_legend(ax)
    save_figure(fig, CASE_DIR / f"rig_p1_pressure_e{round(eps * 100):03d}_n{rpm:04d}")
    plt.close(fig)


def plot_window(rows):
    """Cavity arc against speed, one line per eccentricity, at p_sat = 0.8 bar."""
    apply_style()
    fig, ax = plt.subplots()
    for i, eps in enumerate(sorted({r["eps"] for r in rows})):
        series = sorted([r for r in rows
                         if r["eps"] == eps and r["p_sat_setting"] == "p080"],
                        key=lambda r: r["rpm"])
        if not series:
            continue
        ax.plot([r["rpm"] for r in series], [r["cav_arc_deg"] for r in series],
                marker=MARKERS[i], color=COLORS[i], label=rf"$\varepsilon = {eps}$")
    ax.axhline(90.0, color="black", ls="--", lw=0.9)
    ax.set_xlabel(r"$N$ (rpm)")
    ax.set_ylabel(r"cavity arc (deg)")
    smart_legend(ax)
    save_figure(fig, CASE_DIR / "rig_cavity_arc_vs_speed")
    plt.close(fig)


def plot_exclusion(rows):
    """Degassed p_min against the vaporous floor -- the exclusion argument."""
    apply_style()
    fig, ax = plt.subplots()
    for i, eps in enumerate(sorted({r["eps"] for r in rows})):
        series = sorted([r for r in rows
                         if r["eps"] == eps and r["p_sat_setting"] == "degassed"],
                        key=lambda r: r["rpm"])
        if not series:
            continue
        ax.plot([r["rpm"] for r in series], [r["p_min_bar_abs"] for r in series],
                marker=MARKERS[i], color=COLORS[i], label=rf"$\varepsilon = {eps}$")
    ax.axhline(0.005, color="black", ls="--", lw=0.9)
    ax.axhline(0.1, color="black", ls=":", lw=0.9)
    ax.set_yscale("log")
    ax.set_xlabel(r"$N$ (rpm)")
    ax.set_ylabel(r"$p_{min}$ (bar abs)")
    smart_legend(ax)
    save_figure(fig, CASE_DIR / "rig_degassed_exclusion_margin")
    plt.close(fig)


def main():
    import csv
    ledger = json.loads((CASE_DIR / "sweep_runs.json").read_text())
    with (CASE_DIR / "sweep_results.csv").open(encoding="utf-8") as handle:
        rows = [{k: (float(v) if k in ("eps", "cav_arc_deg", "p_min_bar_abs", "W_N")
                     else int(v) if k == "rpm" else v)
                 for k, v in row.items()}
                for row in csv.DictReader(handle)]
    run1 = [r for r in rows if r["run"] == "RUN1"]
    analysis = json.loads((CASE_DIR / "sweep_analysis.json").read_text())
    plot_window(run1)
    plot_exclusion(run1)
    candidates = analysis.get("Q2_P1_candidates") or []
    if candidates:
        plot_p1_overlay(ledger, candidates[0]["eps"], candidates[0]["rpm"])
    else:
        print("no P1 candidate satisfies both criteria; skipping the overlay")


if __name__ == "__main__":
    main()
