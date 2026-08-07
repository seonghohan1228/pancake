"""Case-05 comparison figures.

  etsion_extent_scatter.png  measured vs computed theta_st / theta_end, 1:1 line
  etsion_extent_vs_ps.png    theta_end and cavity extent vs ambient pressure
  etsion_cavity_map.png      computed theta-z cavity footprint, 1840 rpm sweep
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, save_figure, COLORS  # noqa: E402
import conditions as C  # noqa: E402
import extent as X  # noqa: E402

CASE_DIR = Path(__file__).resolve().parent
REF_BAND = 5.0  # deg, the paper's own max-clearance location uncertainty


def load_rows():
    with open(CASE_DIR / "comparison.csv", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))]


def unwrap_end(a):
    """Map theta_end onto a continuous axis through 0 deg (355 -> -5)."""
    return a - 360.0 if a > 180.0 else a


def fig_scatter(rows):
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
    style = {"fbns": ("o", COLORS[0], "FBNS"), "jfo": ("s", COLORS[1], "JFO / Elrod")}
    for ax, (key_e, key_s, label, lo, hi) in zip(
            axes,
            [("th_st_exp", "th_st_sim", r"$\vartheta_{\rm st}$", 170, 230),
             ("th_end_exp", "th_end_sim", r"$\vartheta_{\rm end}$", -20, 170)]):
        ax.plot([lo, hi], [lo, hi], "-", color="0.35", lw=1.2, zorder=1)
        ax.fill_between([lo, hi], [lo - REF_BAND, hi - REF_BAND],
                        [lo + REF_BAND, hi + REF_BAND],
                        color="0.75", alpha=0.45, lw=0, zorder=0,
                        label=r"reference $\pm5^\circ$")
        for model, (mk, col, lbl) in style.items():
            for rpm, fill in ((1840, col), (3000, "none")):
                xs, ys = [], []
                for r in rows:
                    if r["model"] != model or int(r["rpm"]) != rpm:
                        continue
                    e, s = float(r[key_e]), float(r[key_s])
                    if "end" in key_e:
                        e, s = unwrap_end(e), unwrap_end(s)
                    xs.append(e)
                    ys.append(s)
                ax.plot(xs, ys, mk, mfc=fill, mec=col, ms=9, mew=1.6, ls="none",
                        label=f"{lbl}, {rpm} rpm", zorder=3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"measured {label}  (deg)")
        ax.set_ylabel(f"pancake {label}  (deg)")
    axes[0].legend(loc="upper left", frameon=False, fontsize=11)
    save_figure(fig, str(CASE_DIR / "etsion_extent_scatter"))


def fig_vs_ps(rows):
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0), sharex=True)
    ref = C.read_reference()
    for col, rpm in enumerate(C.SPEEDS_RPM):
        ps = list(C.PS_KPA_GAUGE)
        # top row: theta_end ; bottom row: circumferential extent
        exp_end = [unwrap_end(ref[C.tag(rpm, p)]["theta_end"]) for p in ps]
        exp_ext = [ref[C.tag(rpm, p)]["extent"] for p in ps]
        axes[0][col].errorbar(ps, exp_end, yerr=REF_BAND, fmt="k*", ms=15,
                              capsize=4, lw=1.2, label="Etsion & Ludwig 1982")
        axes[1][col].errorbar(ps, exp_ext, yerr=REF_BAND * math.sqrt(2), fmt="k*",
                              ms=15, capsize=4, lw=1.2, label="Etsion & Ludwig 1982")
        for model, mk, colr, lbl in (("fbns", "o", COLORS[0], "pancake FBNS"),
                                     ("jfo", "s", COLORS[1], "pancake JFO / Elrod")):
            sel = {float(r["ps_kPa"]): r for r in rows
                   if r["model"] == model and int(r["rpm"]) == rpm}
            axes[0][col].plot(ps, [unwrap_end(float(sel[p]["th_end_sim"])) for p in ps],
                              mk + "-", color=colr, ms=8, label=lbl)
            axes[1][col].plot(ps, [float(sel[p]["ext_sim"]) for p in ps],
                              mk + "-", color=colr, ms=8, label=lbl)
        axes[0][col].axhline(0.0, color="0.6", ls=":", lw=1.0)
        axes[0][col].text(0.04, 0.06, f"{rpm} rpm", transform=axes[0][col].transAxes,
                          ha="left", fontsize=15)
        axes[1][col].text(0.04, 0.06, f"{rpm} rpm", transform=axes[1][col].transAxes,
                          ha="left", fontsize=15)
        axes[1][col].set_xlabel(r"ambient pressure $p_s$  (kPa gauge)")
    axes[0][0].set_ylabel(r"cavitation end $\vartheta_{\rm end}$  (deg)")
    axes[1][0].set_ylabel(r"circumferential extent  (deg)")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               fontsize=14, bbox_to_anchor=(0.5, -0.005))
    for a in axes.ravel():
        a.set_xlim(8, 60)
    fig.subplots_adjust(bottom=0.16)
    save_figure(fig, str(CASE_DIR / "etsion_extent_vs_ps"))


def fig_cavity_map():
    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.6), sharey=True)
    ref = C.read_reference()
    for ax, ps in zip(axes, C.PS_KPA_GAUGE):
        run = CASE_DIR / "fbns" / C.tag(1840, ps)
        if not run.exists():
            run = CASE_DIR / "results_remote" / "fbns" / C.tag(1840, ps)
        pts = X.cavity_boundary(run)
        if not pts:
            continue
        th = np.array([p[0] for p in pts])
        z = np.array([p[1] for p in pts])
        th = np.where(th > 180.0, th - 360.0, th)
        ax.plot(th, 1000.0 * (z - 0.5 * 0.0381), ".", ms=1.2, color=COLORS[0],
                rasterized=True)
        r = ref[C.tag(1840, ps)]
        for v, lbl, c in ((unwrap_end(r["theta_st"]), r"measured $\vartheta_{\rm st}$", COLORS[1]),
                          (unwrap_end(r["theta_end"]), r"measured $\vartheta_{\rm end}$", COLORS[2])):
            ax.axvline(v, color=c, ls="--", lw=1.6, label=lbl)
        ax.axvline(0.0, color="0.55", ls=":", lw=1.2, label="max. clearance")
        ax.set_xlim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_xlabel(r"$\vartheta$  (deg)")
        ax.set_title(rf"$p_s$ = {ps} kPa", fontsize=15)
    axes[0].set_ylabel(r"axial position $z$  (mm)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=13,
               bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(bottom=0.30)
    save_figure(fig, str(CASE_DIR / "etsion_cavity_map"))


if __name__ == "__main__":
    rows = load_rows()
    rows = [r for r in rows if r.get("th_st_sim")]
    fig_scatter(rows)
    fig_vs_ps(rows)
    fig_cavity_map()
    print("figures written")
