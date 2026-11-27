"""plot-data house style: PRL look (white, closed box, inward ticks all round), LaTeX variables and units,
PRL colours with Pit Wall accents for talks.

One style for papers and meetings; setup(journal=True) gives all-serif journal figures. Each plot is its own file with a fixed axes box,
so separate files line up when placed together. Scale all files by the same factor on slides.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex

CM = 1 / 2.54
INK, MUTED, CONTEXT = "#15151E", "#6B6E78", "#B9BCC4"   # text, secondary text, grey context
HIGHLIGHT = "#7A2BB8"                                    # Pit Wall purple: the series that carries the message
STATUS = {"good": "#00A651", "near": "#E0A100", "bad": "#E10600"}   # vs a target only
NEUTRAL = "#000000"                                     # single series, models, slope triangles
PRL = ("#000000", "#E00000", "#0000E0", "#008000")       # up to 4 series: black, red, blue, green
MARKERS = ("o", "s", "^", "D")                          # open markers, one per series, same order as PRL
_ANCH = ("#F25C54", "#7A2BB8", "#0B2A66")               # red -> purple -> blue, lightness falls 0.67 -> 0.30


# ---------- colour ----------
def _lin(c): c = np.asarray(c); return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
def _gam(c): c = np.clip(c, 0, 1); return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)
def _ok(h):
    r, g, b = _lin(to_rgb(h))
    l, m, s = np.cbrt([0.4122214708*r + 0.5363325363*g + 0.0514459929*b,
                       0.2119034982*r + 0.6806995451*g + 0.1073969566*b,
                       0.0883024619*r + 0.2817188376*g + 0.6299787005*b])
    return np.array([0.2104542553*l + 0.7936177850*m - 0.0040720468*s,
                     1.9779984951*l - 2.4285922050*m + 0.4505937099*s,
                     0.0259040371*l + 0.7827717662*m - 0.8086757660*s])
def _hex(L):
    l, m, s = (L[0] + 0.3963377774*L[1] + 0.2158037573*L[2]) ** 3, (L[0] - 0.1055613458*L[1] - 0.0638541728*L[2]) ** 3, \
              (L[0] - 0.0894841775*L[1] - 1.2914855480*L[2]) ** 3
    return to_hex(_gam(np.array([4.0767416621*l - 3.3077115913*m + 0.2309699292*s,
                                 -1.2684380046*l + 2.6097574011*m - 0.3413193965*s,
                                 -0.0041960863*l - 0.7034186147*m + 1.7076147010*s])))
def ramp(n):
    """n Pit Wall shades for an ordered parameter: red (low) -> purple -> blue (high)."""
    a, p, b = (_ok(h) for h in _ANCH)
    return [_hex(a + (p - a) * 2 * x if x <= .5 else p + (b - p) * (2 * x - 1)) for x in np.linspace(0, 1, n)]
CMAP = LinearSegmentedColormap.from_list("pitwall", ramp(256))


# ---------- style ----------
def _have(name):
    from matplotlib import font_manager as fm
    return any(f.name == name for f in fm.fontManager.ttflist)


def setup(lang="en", journal=False):
    """Call once before plotting.
    journal=False: sans words and tick numbers, serif math (symbols, units, numbers in math).
    journal=True:  all text serif, in the same face as the math (for journal submission).
    lang="ko" puts a CJK font first (else Hangul vanishes from math labels)."""
    mpl.rcdefaults()
    sans = next((f for f in ("TeX Gyre Heros", "Helvetica", "Arial", "Liberation Sans") if _have(f)), "DejaVu Sans")
    serif = "Latin Modern Roman" if _have("Latin Modern Roman") else "STIXGeneral"   # matches the math fontset below
    kr = next((f for f in ("Noto Sans CJK KR", "Pretendard", "Malgun Gothic", "NanumGothic") if _have(f)), sans)
    if _have("Latin Modern Roman"):   # LaTeX look with upright µ and ° when typed as characters
        math = {"mathtext.fontset": "custom", "mathtext.rm": "Latin Modern Roman",
                "mathtext.it": "Latin Modern Roman:italic", "mathtext.bf": "Latin Modern Roman:bold",
                "mathtext.fallback": "cm"}
    else:                              # STIX also keeps µ upright; 'cm' alone would slant it
        math = {"mathtext.fontset": "stix"}
    mpl.rcParams.update(math)
    mpl.rcParams.update({
        "font.family": [kr] if lang == "ko" else [serif if journal else sans],
        "font.size": 7, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
        "legend.frameon": False, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "lines.markersize": 3.8,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8, "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
        "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "svg.fonttype": "path", "pdf.fonttype": 42,   # SVG as outlines: looks the same on any machine
    })


def figure(w=6.0, h=4.5, left=1.25, bottom=1.0, right=0.25, top=0.25):
    """Fixed axes box (cm). Same defaults for every plot -> files align when placed together."""
    W, H = (left + w + right) * CM, (bottom + h + top) * CM
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([left * CM / W, bottom * CM / H, w * CM / W, h * CM / H])
    return fig, ax


def save(fig, path):
    """Write PDF + SVG (vector) and a 600 dpi PNG. No tight bbox: it would change the file size."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{path}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)


def slope_triangle(ax, x0, y0, decades, slope, color=NEUTRAL, label=None, fs=6.5):
    """Log-log slope triangle; colour black or the colour of the data it describes."""
    x1 = x0 * 10 ** decades; y1 = y0 * (x1 / x0) ** slope
    ax.plot([x0, x1, x1, x0], [y0, y0, y1, y0], color=color, lw=0.6)
    ax.text(np.sqrt(x0 * x1), y0 * (0.82 if slope > 0 else 1.1), "1", color=color, ha="center",
            va="top" if slope > 0 else "bottom", fontsize=fs)
    ax.text(x1 * 1.1, np.sqrt(y0 * y1), label or f"${slope:g}$", color=color, va="center", fontsize=fs)


def colorbar_inset(fig, ax, norm, label, loc=(0.60, 0.12, 0.32, 0.04), ticks=None, ticklabels=None):
    cax = ax.inset_axes(loc)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.4); cax.tick_params(labelsize=6, length=1.5, direction="out", top=False, right=False)
    if ticks is not None: cb.set_ticks(ticks)
    if ticklabels is not None: cb.set_ticklabels(ticklabels)
    cax.set_title(label, fontsize=6.5, pad=2)
    return cb


def note(ax, text, xy, xytext, color=INK, fs=6.5, **kw):
    """Arrow annotation in the house style. Always use this: matplotlib's default arrowhead colour is blue."""
    return ax.annotate(text, xy=xy, xytext=xytext, fontsize=fs, color=color,
                       arrowprops=dict(arrowstyle="-|>", lw=0.5, color=color, shrinkA=1, shrinkB=2), **kw)
