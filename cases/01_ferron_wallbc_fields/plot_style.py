"""Publication plot style for Seongho Han (KAIST Applied Fluid Mechanics Lab).

Usage:
    from plot_style import apply_style, save_figure, MARKERS, blues, clean_log_ticks, fit_slope_guide, slope_triangle
    apply_style()            # STIX/Times mathtext (works everywhere)
    apply_style(usetex=True) # full LaTeX with Times (newtx); falls back automatically
    ... build figure ...
    save_figure(fig, "jet_velocity")   # -> jet_velocity.png (300 dpi)
"""

import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker
from cycler import cycler

# Colorblind-friendly cycle (Wong palette), dark blue first.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
# Clear filled marker shapes, in preference order.
MARKERS = ["s", "o", "^", "D", "v", "<", ">"]

RC = {
    # Fonts: Times New Roman where installed (e.g. Windows); STIX fallback ships
    # with matplotlib and is visually near-identical.
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "Nimbus Roman", "Liberation Serif"],
    "mathtext.fontset": "stix",
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    # Box frame, but inward ticks only on the labeled sides (bottom/left).
    "axes.linewidth": 1.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,
    "ytick.right": False,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.grid": False,
    # Markers and lines.
    "lines.markersize": 8,
    "lines.markeredgewidth": 1.0,
    "lines.linewidth": 1.5,
    "errorbar.capsize": 3,
    "axes.prop_cycle": cycler(color=COLORS),
    # Legend: no frame, inside the axes.
    "legend.frameon": False,
    # Figure geometry and export.
    "figure.figsize": (4.5, 4.2),
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
}

_USETEX_PREAMBLE = r"\usepackage{newtxtext}\usepackage{newtxmath}"


def apply_style(usetex: bool = False) -> None:
    """Apply the lab style. usetex=True attempts full LaTeX (Times via newtx)
    and silently falls back to mathtext if the TeX toolchain is unavailable."""
    mpl.rcParams.update(RC)
    if usetex:
        mpl.rcParams.update({"text.usetex": True, "text.latex.preamble": _USETEX_PREAMBLE})
        try:  # probe: render a math string once
            fig = plt.figure(figsize=(1, 1))
            fig.text(0.5, 0.5, r"$\lambda_f$")
            fig.canvas.draw()
            plt.close(fig)
        except Exception:
            plt.close("all")
            mpl.rcParams.update({"text.usetex": False})
            warnings.warn("usetex unavailable (missing LaTeX packages); using STIX mathtext.")


def blues(n: int, dark_first: bool = True):
    """Sequential blues for a swept parameter (dark -> light), as in the lab's
    reference figures where marker shade encodes a parameter value."""
    lo, hi = 0.35, 0.95
    vals = np.linspace(hi, lo, n) if dark_first else np.linspace(lo, hi, n)
    return plt.cm.Blues(vals)


def save_figure(fig, basename: str, formats=("png",), dpi: int = 300):
    """Save the figure. Default is a single 300 dpi PNG; pass e.g.
    formats=("png", "eps") only when another format is explicitly requested.
    (EPS flattens transparency - keep marker faces opaque if exporting EPS.)"""
    paths = []
    for ext in formats:
        p = f"{basename}.{ext}"
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    return paths



class _PlainLogFormatter(mticker.LogFormatterSciNotation):
    """Log-axis formatter that writes 'obvious' numbers plainly (6, not
    6x10^0) and falls back to 10^n notation outside 0.01-999. It respects
    the parent class's decisions about WHICH ticks get labels."""

    def __call__(self, x, pos=None):
        s = super().__call__(x, pos)
        if s and 0.01 <= x < 1000:
            return f"{x:g}"
        return s


def clean_log_ticks(ax):
    """Apply plain-number tick labels to any log-scaled axis of `ax`.
    Call AFTER setting scales and limits on every log/log-log plot."""
    for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
        if scale == "log":
            axis.set_major_formatter(_PlainLogFormatter())
            axis.set_minor_formatter(_PlainLogFormatter())


def fit_slope_guide(ax, x, y, exponent=None, x_range=None,
                    color: str = "k", ls: str = "--", lw: float = 1.2):
    """Dashed power-law line FITTED to the data (log-log): the guide is a
    representation of the data, not an offset annotation. If `exponent` is
    given, only the prefactor is least-squares fitted; if None, both are.
    `x_range=(lo, hi)` restricts the fit and the drawn span to the region
    where the scaling holds. Returns (Line2D, prefactor, exponent)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = (x > 0) & (y > 0)
    if x_range is not None:
        m &= (x >= x_range[0]) & (x <= x_range[1])
    lx, ly = np.log(x[m]), np.log(y[m])
    if exponent is None:
        exponent, b = np.polyfit(lx, ly, 1)
    else:
        b = float(np.mean(ly - exponent * lx))
    pref = np.exp(b)
    xs = np.array([x[m].min(), x[m].max()])
    (line,) = ax.plot(xs, pref * xs ** exponent, ls=ls, color=color, lw=lw, zorder=1)
    return line, pref, exponent


def add_slope_guide(ax, exponent: float, x0: float, x1: float, anchor_y: float,
                    color: str = "k", ls: str = "--", lw: float = 1.2):
    """Manually anchored power-law guide y = anchor_y * (x/x0)^exponent from x0
    to x1. Prefer fit_slope_guide, which anchors to the data; use this only
    when Seongho asks for a specific placement. Returns the Line2D."""
    x = np.array([x0, x1], dtype=float)
    y = anchor_y * (x / x0) ** exponent
    (line,) = ax.plot(x, y, ls=ls, color=color, lw=lw, zorder=1)
    return line


def slope_triangle(ax, x0: float, y0: float, exponent: float, width_decades: float = 0.4,
                   color: str = "k", fontsize: int = 14, above: bool = True):
    """Right-triangle slope marker for log-log axes, corner at (x0, y0).
    Labels the run '1' and the rise with the exponent."""
    x1 = x0 * 10 ** width_decades
    y1 = y0 * (x1 / x0) ** exponent
    if above:
        ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color=color, lw=1.0, zorder=1)
        ax.text(np.sqrt(x0 * x1), y0, "1", ha="center", va="top", fontsize=fontsize)
        ax.text(x1, np.sqrt(y0 * y1), _fmt_exp(exponent), ha="left", va="center", fontsize=fontsize)
    else:
        ax.plot([x0, x1, x0, x0], [y0, y1, y1, y0], color=color, lw=1.0, zorder=1)
        ax.text(np.sqrt(x0 * x1), y1, "1", ha="center", va="bottom", fontsize=fontsize)
        ax.text(x0, np.sqrt(y0 * y1), _fmt_exp(exponent), ha="right", va="center", fontsize=fontsize)


def _fmt_exp(e: float) -> str:
    from fractions import Fraction
    fr = Fraction(e).limit_denominator(12)
    if fr.denominator == 1:
        return f"{fr.numerator}"
    return f"{fr.numerator}/{fr.denominator}"
