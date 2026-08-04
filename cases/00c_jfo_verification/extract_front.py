"""Extract rupture/reformation angles, p_max and W' from a flat VTK run.

Reads every processor .vts of a `results/flat` tree, stitches the mid-z row in
circumferential order, and reports the cavitation front positions so a run can be
compared with `exact_1d_jfo.py`.

A cell is called cavitated when film_content < 0.999. The rupture angle is the
face between the last full-film cell and the first cavitated cell downstream of
the groove; the reformation angle is the face where the film becomes full again.
Reporting the face (not the cell centre) removes the half-cell bias, so the
remaining error is the scheme's.

Usage: python extract_front.py <case_dir> [<case_dir> ...]
"""
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

P_CAV = 1.0e5
R = 0.05
GROOVE_CENTER = 90.0
GROOVE_HALF_WIDTH = 9.0
FULL_FILM = 1.0 - 1.0e-6


def read_piece(path):
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    extent = [int(v) for v in piece.get("Extent").split()]
    nx = extent[1] - extent[0]
    ny = extent[3] - extent[2]
    arrays = {}
    for da in piece.iter("DataArray"):
        name = da.get("Name")
        if name == "Points":
            continue
        ncomp = int(da.get("NumberOfComponents", "1"))
        values = np.fromstring(da.text.strip(), sep=" ")
        if ncomp == 1:
            arrays[name] = values.reshape(ny, nx)
    return extent[0], nx, ny, arrays


def load_mid_row(case_dir):
    files = sorted(Path(case_dir).glob("results/flat/processor*/*.vts"))
    if not files:
        raise SystemExit(f"no vts files under {case_dir}")
    pieces = [read_piece(f) for f in files]
    pieces.sort(key=lambda p: p[0])
    n_theta = sum(p[1] for p in pieces)
    n_z = pieces[0][2]
    j = n_z // 2
    out = {}
    for name in pieces[0][3]:
        out[name] = np.concatenate([p[3][name][j] for p in pieces])
    assert len(out["p"]) == n_theta
    return n_theta, out


def analyse(case_dir):
    n_theta, row = load_mid_row(case_dir)
    d_theta = 360.0 / n_theta
    centers = (np.arange(n_theta) + 0.5) * d_theta
    film = row["film_content"]
    p = row["p"]

    cavitated = film < FULL_FILM

    # Walk downstream (increasing theta) from the groove trailing edge, wrapping
    # once, and take the first full->cavitated and cavitated->full transitions.
    # The front is reported at the face between the two cells, so the estimator
    # itself carries a +-0.5 cell quantization and no half-cell bias.
    start = int(np.argmax(centers > GROOVE_CENTER + GROOVE_HALF_WIDTH))
    order = [(start + k) % n_theta for k in range(n_theta)]
    rupture = reformation = None
    for k in range(1, len(order)):
        prev, cur = order[k - 1], order[k]
        face = (centers[prev] + 0.5 * d_theta) % 360.0
        if rupture is None and not cavitated[prev] and cavitated[cur]:
            rupture = face
        elif rupture is not None and reformation is None and cavitated[prev] and not cavitated[cur]:
            reformation = face

    # Midpoint rule over the cell-centred pressures: the finite-volume field's own
    # quadrature, so the load carries no extra interpolation error.
    theta_rad = np.radians(centers)
    gauge = p - P_CAV
    d_rad = math.radians(d_theta)
    fx = -float(np.sum(gauge * np.cos(theta_rad))) * d_rad * R
    fy = -float(np.sum(gauge * np.sin(theta_rad))) * d_rad * R

    residual = row.get("mass_residual")
    res_max = float(np.max(np.abs(residual))) if residual is not None else float("nan")
    return {
        "n": n_theta,
        "rupture": rupture,
        "reformation": reformation,
        "p_max": float(p.max()),
        "p_max_theta": float(centers[int(p.argmax())]),
        "load": float(math.hypot(fx, fy)),
        "cell_deg": d_theta,
        "mass_residual_max": res_max,
        "theta_min": float(film.min()),
    }


if __name__ == "__main__":
    exact_rupture = 302.766608
    exact_reformation = 79.006725
    exact_pmax = 7.597141e6
    exact_load = 495463.6
    print(f"{'case':<26}{'n':>6}{'dtheta':>8}{'rupture':>10}{'err':>9}{'err/cell':>9}"
          f"{'reform':>10}{'err':>9}{'err/cell':>9}{'p_max[MPa]':>12}{'W[N/m]':>12}")
    for case in sys.argv[1:]:
        r = analyse(case)
        er = r["rupture"] - exact_rupture if r["rupture"] is not None else float("nan")
        ef = r["reformation"] - exact_reformation if r["reformation"] is not None else float("nan")
        if ef > 180.0:
            ef -= 360.0
        print(f"{Path(case).name:<26}{r['n']:>6}{r['cell_deg']:>8.3f}"
              f"{r['rupture']:>10.3f}{er:>9.3f}{er / r['cell_deg']:>9.2f}"
              f"{r['reformation']:>10.3f}{ef:>9.3f}{ef / r['cell_deg']:>9.2f}"
              f"{r['p_max'] / 1e6:>12.4f}{r['load']:>12.1f}")
        print(f"{'':<26}  p_max err {100 * (r['p_max'] - exact_pmax) / exact_pmax:+.3f} %"
              f" | W' err {100 * (r['load'] - exact_load) / exact_load:+.3f} %"
              f" | min film_content {r['theta_min']:.4g}"
              f" | max |mass_residual| {r['mass_residual_max']:.3e} kg/s")
