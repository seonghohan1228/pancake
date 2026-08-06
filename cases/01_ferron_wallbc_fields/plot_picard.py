"""Outer (Picard) coupling residual history of the conjugate Ferron CHT point.

Reads picard_study/history.json (one entry per relaxation setting, produced by
running the archived rerun_cht_robin10/fb_iter_07 config with
log_coupling_history = true on afml-r6) and draws max(R_p, R_T) against the
outer sweep index.

Output: ferron_picard_convergence.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import apply_style, save_figure

CASE_DIR = Path(__file__).resolve().parent
TOL = 1.0e-6

# label, key in history.json, style
SERIES = [
    (r"$\omega = 1$ (frozen protocol)", "base", dict(color="#D55E00", lw=2.2)),
    (r"$\omega = 0.95$", "r095", dict(color="#0072B2", lw=1.4)),
    (r"$\omega = 0.8$", "r08", dict(color="#009E73", lw=1.4)),
    (r"$\omega = 0.5$", "r05", dict(color="#56B4E9", lw=1.4)),
    (r"$\omega = 0.3$", "r03", dict(color="#CC79A7", lw=1.4)),
    ("adaptive", "adapt", dict(color="black", lw=2.2, ls="--")),
]


def main():
    apply_style()
    hist = json.loads((CASE_DIR / "picard_study" / "history.json").read_text())

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for label, key, style in SERIES:
        if key not in hist:
            continue
        rows = hist[key]
        n = [r[0] for r in rows]
        res = [max(r[1], r[2]) for r in rows]
        ax.semilogy(n, res, label=label, **style)

    ax.axhline(TOL, color="0.35", lw=1.0, ls=":")
    ax.annotate("coupling tolerance", xy=(126, TOL * 1.8), fontsize=13, color="0.35")

    # Mark the iterate the frozen-protocol run should have delivered.
    base = hist.get("base")
    if base:
        res = [max(r[1], r[2]) for r in base]
        i_best = min(range(len(res)), key=lambda i: res[i])
        ax.plot([base[i_best][0]], [res[i_best]], marker="o", ms=9, mfc="none",
                mec="#D55E00", mew=2.0, ls="none")
        ax.annotate(f"best iterate\n(outer {base[i_best][0]})",
                    xy=(base[i_best][0], res[i_best]),
                    xytext=(base[i_best][0] + 18, res[i_best] * 0.12),
                    fontsize=13, color="#D55E00",
                    arrowprops=dict(arrowstyle="-", color="#D55E00", lw=1.0))

    ax.set_xlabel("outer (Picard) sweep")
    ax.set_ylabel(r"$\max(R_p,\,R_T)$")
    ax.set_xlim(0, 205)
    ax.set_ylim(1e-9, 1e2)
    ax.legend(frameon=False, loc="upper right", ncol=2)
    save_figure(fig, str(CASE_DIR / "ferron_picard_convergence"))


if __name__ == "__main__":
    main()
