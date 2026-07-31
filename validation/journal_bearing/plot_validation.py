"""Figures for the cavitation / energy-equation validation.

Consumes `raimondi_boyd_results.json` and `ferron_results.json` written by
`run_raimondi_boyd.py` and `run_ferron.py`.

Usage:  $PY plot_validation.py
"""
import json, math, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import apply_style, save_figure, MARKERS, COLORS, clean_log_ticks
from raimondi_boyd import RB_TABLE
import ferron as F

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        print(f"skipping: {path} not found")
        return None
    return json.load(open(path))


def fig_rb_sommerfeld(rb):
    """Sommerfeld number against eccentricity, computed vs the chart."""
    fig, ax = plt.subplots(figsize=(4.7, 4.2))
    for n, (ld_key, colour) in enumerate(zip(("0.5", "1"), (COLORS[0], COLORS[1]))):
        ref = RB_TABLE[ld_key]
        ax.plot([r[0] for r in ref], [r[1] for r in ref], "-", color=colour, lw=1.4,
                label=f"Raimondi \\& Boyd, $L/D = {ld_key}$" if False
                      else f"Raimondi & Boyd, $L/D$ = {ld_key}")
        pts = [r for r in rb["sweep"]
               if r["ld_key"] == ld_key and r["inlet"] == "flooded"
               and abs(r["mu"] - 0.01) < 1e-9]
        pts.sort(key=lambda r: r["eps"])
        ax.plot([r["eps"] for r in pts], [r["S"] for r in pts], MARKERS[n],
                color=colour, ms=6, ls="none", mfc=colour, mec="k", mew=0.6,
                label=f"pancake JFO, fed film, $L/D$ = {ld_key}")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$S$")
    ax.set_xlim(0.1, 1.0)
    clean_log_ticks(ax)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    save_figure(fig, "validation_raimondi_boyd_S")
    plt.close(fig)


def fig_rb_attitude(rb):
    fig, ax = plt.subplots(figsize=(4.7, 4.2))
    for n, (ld_key, colour) in enumerate(zip(("0.5", "1"), (COLORS[0], COLORS[1]))):
        ref = RB_TABLE[ld_key]
        ax.plot([r[0] for r in ref], [r[2] for r in ref], "-", color=colour, lw=1.4,
                label=f"Raimondi & Boyd, $L/D$ = {ld_key}")
        pts = [r for r in rb["sweep"]
               if r["ld_key"] == ld_key and r["inlet"] == "flooded"
               and abs(r["mu"] - 0.01) < 1e-9]
        pts.sort(key=lambda r: r["eps"])
        ax.plot([r["eps"] for r in pts], [r["phi"] for r in pts], MARKERS[n],
                color=colour, ms=6, ls="none", mfc=colour, mec="k", mew=0.6,
                label=f"pancake JFO, fed film, $L/D$ = {ld_key}")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$\phi$ (deg)")
    ax.set_xlim(0.1, 1.0)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    save_figure(fig, "validation_raimondi_boyd_attitude")
    plt.close(fig)


def fig_free_boundary(rb):
    """Rupture-angle convergence against the exact Swift-Stieber value.  This is
    the numerical prerequisite for every other cavitation result here."""
    cvg = rb.get("convergence") or []
    if not cvg:
        from raimondi_boyd import swift_stieber_infinite
        from run_raimondi_boyd import CONVERGENCE_MEASURED
        ex = swift_stieber_infinite(0.6)
        cvg = [dict(r, S_exact=ex["S"], phi_exact=ex["phi_deg"],
                    rupture_exact=ex["t2_deg"]) for r in CONVERGENCE_MEASURED]
    fig, ax = plt.subplots(figsize=(4.7, 4.2))
    exact = cvg[0]["rupture_exact"]
    ax.axhline(exact, color="k", ls="--", lw=1.2, label="exact Swift–Stieber")
    ticks = sorted({c["nth"] for c in cvg})
    for n, beta in enumerate(sorted({c["beta"] for c in cvg}, reverse=True)):
        rows = sorted((c for c in cvg if c["beta"] == beta), key=lambda c: c["nth"])
        ax.plot([c["nth"] for c in rows], [c["rupture"] for c in rows],
                MARKERS[n] + "-", color=COLORS[n], ms=6, mfc=COLORS[n],
                mec="k", mew=0.6,
                label=rf"$\beta = 10^{{{round(math.log10(beta))}}}$ Pa")
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel(r"$N_\theta$")
    ax.set_ylabel(r"$\theta_{\mathrm{rupture}}$ (deg)")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    save_figure(fig, "validation_free_boundary_convergence")
    plt.close(fig)


def fig_ferron_pressure(fr):
    fig, ax = plt.subplots(figsize=(4.9, 4.2))

    # Thermal bracket as a band: isothermal and THD at a COMMON eccentricity.
    # The load-matched cases cannot form a band, because matching the load forces
    # their peak pressures together by construction.
    br = fr.get("bracket")
    if br and "THD" in br and "ISOTHERMAL" in br:
        # The two bracket runs settle on slightly different groove angles, so
        # their X grids do not coincide; interpolate both onto a common grid.
        import numpy as np
        grid = np.linspace(-1.0, 1.0, 801)
        def resample(prof):
            xs_p = np.array([r[0] for r in prof])
            ys_p = np.array([r[1] for r in prof])
            order = np.argsort(xs_p)
            return np.interp(grid, xs_p[order], ys_p[order])
        lo = resample(br["THD"]["profile"])
        hi = resample(br["ISOTHERMAL"]["profile"])
        ax.fill_between(grid, np.minimum(lo, hi), np.maximum(lo, hi),
                        color=COLORS[2], alpha=0.30, lw=0,
                        label=r"thermal effect at fixed $\varepsilon$")

    xs = [p[0] for p in F.FERRON_PRESSURE_MPA]
    ps = [p[1] for p in F.FERRON_PRESSURE_MPA]
    ax.errorbar(xs, ps, yerr=F.PRESSURE_DIGITISE_ERR_MPA,
                xerr=F.PRESSURE_DIGITISE_ERR_XPI, fmt=MARKERS[0], color="k",
                ms=5, mfc="k", mec="k", lw=0, elinewidth=0.8, capsize=2,
                label="Ferron et al. 1983 (digitised)")
    # Only the THD case is drawn as a line.  The load-matched isothermal curve is
    # omitted deliberately: matching the 4 kN load forces the two peak pressures
    # together, so it lies on top of the THD curve and says nothing.  The shaded
    # band above carries the thermal effect instead.
    if "thd_adiabatic" in fr:
        prof = fr["thd_adiabatic"]["profile"]
        ax.plot([r[0] for r in prof], [r[1] for r in prof], "-",
                color=COLORS[0], lw=1.5, label="pancake THD, 4 kN")
    ax.set_xlabel(r"$X/\pi$")
    ax.set_ylabel(r"$p$ (MPa)")
    ax.set_xlim(-1.0, 1.0)
    # Add headroom above the peak so the legend has somewhere to sit that is not
    # on top of the data.
    ax.set_ylim(-0.05, 1.72)
    ax.legend(frameon=False, fontsize=9, loc="upper left", handlelength=1.6)
    save_figure(fig, "validation_ferron_pressure")
    plt.close(fig)


def fig_ferron_temperature(fr):
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    xs = [p[0] for p in F.FERRON_TEMPERATURE_C]
    Ts = [p[1] for p in F.FERRON_TEMPERATURE_C]
    ax.errorbar(xs, Ts, yerr=F.TEMPERATURE_DIGITISE_ERR_K, fmt=MARKERS[1],
                color="k", ms=5, mfc="none", mec="k", lw=0, elinewidth=0.8,
                capsize=2, label="Ferron et al. 1983 (digitised)")
    for n, (tag, label) in enumerate((("thd_adiabatic", "pancake THD, adiabatic"),
                                      ("thd_metalbush", "pancake THD, metal bush"))):
        if tag not in fr:
            continue
        prof = fr[tag]["profile"]
        ax.plot([r[0] for r in prof], [r[2] for r in prof], "-",
                color=COLORS[n], lw=1.4, label=label)
    ax.set_xlabel(r"$X/\pi$")
    ax.set_ylabel(r"$T$ ($^\circ$C)")
    ax.set_xlim(-1.0, 1.0)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    save_figure(fig, "validation_ferron_temperature")
    plt.close(fig)


def fig_ferron_load(fr):
    lc = fr.get("load_curve") or []
    if not lc:
        return
    fig, ax = plt.subplots(figsize=(4.7, 4.2))
    ax.errorbar([r["W_meas"] * 1e-3 for r in F_pairs()],
                [r["eps"] for r in F_pairs()], yerr=F.ECC_DIGITISE_ERR,
                fmt=MARKERS[0], color="k", ms=5, mfc="k", mec="k", lw=0,
                elinewidth=0.8, capsize=2, label="Ferron et al. 1983 (digitised)")
    lc = sorted(lc, key=lambda r: r["W_sim"])
    ax.plot([r["W_sim"] * 1e-3 for r in lc], [r["eps"] for r in lc],
            MARKERS[1] + "-", color=COLORS[0], ms=6, mfc=COLORS[0], mec="k",
            mew=0.6, lw=1.4, label="pancake THD")
    ax.set_xlabel(r"$W$ (kN)")
    ax.set_ylabel(r"$\varepsilon$")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    save_figure(fig, "validation_ferron_load")
    plt.close(fig)


def F_pairs():
    return [dict(W_meas=w, eps=e) for w, e in F.FERRON_ECC_VS_LOAD]


def main():
    apply_style()
    rb = load("raimondi_boyd_results.json")
    if rb:
        fig_rb_sommerfeld(rb)
        fig_rb_attitude(rb)
        fig_free_boundary(rb)
    fr = load("ferron_results.json")
    if fr:
        fig_ferron_pressure(fr)
        fig_ferron_temperature(fr)
        fig_ferron_load(fr)
    print("done")


if __name__ == "__main__":
    main()
