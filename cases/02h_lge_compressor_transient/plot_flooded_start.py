"""Flooded-start sweep figures (house style), one file per quantity and speed.

  flooded_hmin_<f>hz        minimum film thickness against refrigerant fraction, mean load,
                            bars = minimum to maximum load; Lambda = 3 levels of both roughnesses
  flooded_contact_<f>hz     share of the load carried by asperity contact
  flooded_power_<f>hz       friction power (film + contact)

Usage: python plot_flooded_start.py --work CAMPAIGN
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from plotstyle import MARKERS, MUTED, PRL, figure, note, save, setup

STYLE = {"N": (PRL[0], MARKERS[0], "no contact model"),
         "C2": (PRL[1], MARKERS[1], r"contact, $\sigma=0.2\ \mathrm{µm}$"),
         "C4": (PRL[2], MARKERS[2], r"contact, $\sigma=0.4\ \mathrm{µm}$")}
OFFSET = {"N": -0.008, "C2": 0.0, "C4": 0.008}


def series(rows, set_key, variant, key):
    out = []
    for fluid in ("F0", "F1", "F2", "F3", "F4"):
        sel = {r["load_N"]: r for r in rows if r["set"] == set_key and r["variant"] == variant and r["fluid"] == fluid}
        if not sel:
            continue
        loads = sorted(sel)
        mean = sel[loads[1]]
        vals = [sel[l].get(key) if sel[l]["converged"] else None for l in loads]
        good = [v for v in vals if v is not None]
        out.append(dict(c=mean.get("c") or rows_c(rows, fluid), mean=vals[1],
                        lo=min(good) if good else None, hi=max(good) if good else None, converged=mean["converged"],
                        last_h=mean.get("last_valid_h_min_um")))
    return out


def rows_c(rows, fluid):
    return next(r["c"] for r in rows if r["fluid"] == fluid and r.get("c") is not None)


def draw(ax, rows, set_key, key, scale=1.0):
    for variant, (color, marker, label) in STYLE.items():
        s = series(rows, set_key, variant, key)
        c = np.array([p["c"] for p in s]) + OFFSET[variant]
        ok = [i for i, p in enumerate(s) if p["mean"] is not None]
        ax.plot(c[ok], [s[i]["mean"] * scale for i in ok], marker, mfc="none", mec=color, color=color, ls="-", lw=0.6,
                label=label)
        for i in ok:
            if s[i]["lo"] is not None and s[i]["hi"] != s[i]["lo"]:
                ax.plot([c[i]] * 2, [s[i]["lo"] * scale, s[i]["hi"] * scale], color=color, lw=0.6)
        missing = [i for i, p in enumerate(s) if p["mean"] is None]
        yield variant, c, s, missing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True, type=Path)
    args = parser.parse_args()
    rows = json.loads((args.work / "summary.json").read_text())["rows"]
    out = args.work / "figures"
    setup()
    for set_key, freq in (("A", 130), ("B", 30)):
        fig, ax = figure()
        ax.set_yscale("log")
        for variant, c, s, missing in draw(ax, rows, set_key, "h_min_um"):
            color = STYLE[variant][0]
            for i in missing:
                ax.plot(c[i], 1.3e-2, "x", color=color)
        for sigma, ls in ((0.2, "--"), (0.4, ":")):
            ax.axhline(3 * sigma, color=MUTED, lw=0.6, ls=ls)
            ax.text(0.985, 3 * sigma * 1.07,  rf"$\Lambda=3,\ \sigma={sigma}\ \mathrm{{µm}}$", fontsize=6.5, color=MUTED, ha="right")
        ax.set_xlim(0.15, 1.0)
        ax.set_ylim(1e-2, 4)
        ax.set_xlabel(r"Refrigerant mass fraction, $c$")
        ax.set_ylabel(r"Minimum film thickness, $h_\mathrm{min}$ ($\mathrm{µm}$)")
        ax.legend(loc="lower left")
        save(fig, str(out / f"flooded_hmin_{freq}hz"))

        fig, ax = figure()
        list(draw(ax, rows, set_key, "contact_share_pct"))
        ax.axhline(10, color=MUTED, lw=0.6, ls="--")
        ax.set_xlim(0.15, 1.0)
        top = max(r.get("contact_share_pct") or 0 for r in rows if r["set"] == set_key)
        ax.set_ylim(-1, 5 * math.ceil(top / 5 + 1))
        ax.text(0.17, 10.4, r"$10\ \%$", fontsize=6.5, color=MUTED)
        ax.set_xlabel(r"Refrigerant mass fraction, $c$")
        ax.set_ylabel(r"Load carried by contact ($\%$)")
        ax.legend(loc="upper left")
        save(fig, str(out / f"flooded_contact_{freq}hz"))

        fig, ax = figure()
        ax.set_yscale("log")
        list(draw(ax, rows, set_key, "friction_power_W"))
        ax.set_xlim(0.15, 1.0)
        ax.set_xlabel(r"Refrigerant mass fraction, $c$")
        ax.set_ylabel(r"Friction power, $P_\mathrm{f}$ ($\mathrm{W}$)")
        ax.legend(loc="upper left")
        save(fig, str(out / f"flooded_power_{freq}hz"))
    print(out)


if __name__ == "__main__":
    main()
