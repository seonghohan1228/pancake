"""Cavity-extent extraction from a finished case-05 run.

Photographed bubbles must be compared with TOTAL VOID, 1 - liquid content.
The old test theta < 1 - alpha_gas counted only the additional vapor deficit;
it excluded gas-filled cells even though the experiment concerns gaseous
cavitation. Archived case-05 runs had alpha_gas = 0 and are unaffected.
The threshold is a numerical void threshold, not a measured optical detection
limit. Report threshold sensitivity when comparing diffuse gas fronts.

Two cavity definitions are produced:
  * `any_z`   -- a theta station counts as cavitated if ANY cell at that theta is.
                 This is the closer match to the reference: Etsion & Ludwig marked
                 START and END on AXIAL traces, i.e. a station was inside the
                 cavity if its axial scan showed cavitation anywhere.
  * `midplane`-- the same test restricted to the z = L/2 row.

Angles come straight out of the solver's own theta_rad field, and case 05 sets
attitude_angle_deg = 180 so solver theta already equals the paper's theta (from
maximum clearance, positive with rotation). No offset is applied here.
"""

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

CEILING_TOL = 1e-6


def _piece_arrays(vts_path):
    root = ET.parse(vts_path).getroot()
    out = {}
    for cd in root.iter("CellData"):
        for da in cd.findall("DataArray"):
            ncomp = int(da.get("NumberOfComponents", "1"))
            out[da.get("Name")] = ([float(v) for v in (da.text or "").split()], ncomp)
    return out


def read_pvts(pvts_path):
    root = ET.parse(pvts_path).getroot()
    merged = {}
    for piece in root.iter("Piece"):
        src = pvts_path.parent / piece.get("Source")
        for name, (vals, ncomp) in _piece_arrays(src).items():
            merged.setdefault(name, ([], ncomp))[0].extend(vals)
    return merged


def final_pvts(run_dir):
    pvd = Path(run_dir) / "results" / "flat" / "results.pvd"
    if not pvd.exists():
        return None
    root = ET.parse(pvd).getroot()
    items = [(float(ds.get("timestep")), pvd.parent / ds.get("file"))
             for ds in root.iter("DataSet")]
    if not items:
        return None
    items.sort(key=lambda x: x[0])
    return items[-1][1]


def _wrap(a):
    return a % 360.0


def _longest_cavitated_arc(theta_deg, flags):
    """Longest contiguous run of True on a periodic theta axis.

    Returns (theta_start, theta_end, extent_deg, n_cells) with the fronts placed
    at the FACE between the last full-film and the first cavitated cell, which is
    the convention PHYSICS.md 7.3.5 uses when it scores fronts against the exact
    1-D solution. Returns None when nothing is cavitated.
    """
    n = len(theta_deg)
    if n == 0 or not any(flags):
        return None
    if all(flags):
        return (0.0, 360.0, 360.0, n)
    # rotate so index 0 is the first cavitated cell whose predecessor is full film
    start = next(i for i in range(n) if flags[i] and not flags[(i - 1) % n])
    best = None
    i = start
    seen = 0
    while seen < n:
        if flags[i] and not flags[(i - 1) % n]:
            j, count = i, 0
            while flags[j] and count < n:
                j = (j + 1) % n
                count += 1
            if best is None or count > best[1]:
                best = (i, count)
        i = (i + 1) % n
        seen += 1
    i0, count = best
    d = 360.0 / n
    th_start = _wrap(theta_deg[i0] - 0.5 * d)
    th_end = _wrap(theta_deg[(i0 + count - 1) % n] + 0.5 * d)
    return (th_start, th_end, count * d, count)


def cavity_metrics(run_dir):
    """Extent metrics of the final frame. None if the run produced no output."""
    path = final_pvts(run_dir)
    if path is None:
        return None
    a = read_pvts(path)
    theta_rad = a["theta_rad"][0]
    z = a["z_m"][0]
    fc = a["film_content"][0]
    p = a["p"][0]
    h = a["h"][0]
    alpha = a["alpha_gas"][0] if "alpha_gas" in a else [0.0] * len(fc)

    cav = [fc[k] < 1.0 - CEILING_TOL for k in range(len(fc))]

    theta_deg = [_wrap(math.degrees(t)) for t in theta_rad]
    z_sorted = sorted(set(z))
    z_mid = z_sorted[len(z_sorted) // 2]

    by_theta_any = {}
    by_theta_mid = {}
    p_mid = {}
    width = {}          # cavitated axial width [m] per theta
    for k in range(len(fc)):
        t = round(theta_deg[k], 6)
        by_theta_any[t] = by_theta_any.get(t, False) or cav[k]
        if z[k] == z_mid:
            by_theta_mid[t] = cav[k]
            p_mid[t] = p[k]
    dz = (max(z_sorted) - min(z_sorted)) / (len(z_sorted) - 1) if len(z_sorted) > 1 else 0.0
    for k in range(len(fc)):
        t = round(theta_deg[k], 6)
        if cav[k]:
            width[t] = width.get(t, 0.0) + dz

    ths = sorted(by_theta_any)
    arc_any = _longest_cavitated_arc(ths, [by_theta_any[t] for t in ths])
    arc_mid = _longest_cavitated_arc(ths, [by_theta_mid.get(t, False) for t in ths])

    # theta_0 -- the reference's "location of minimum FILM pressure", i.e. of the
    # sub-cavity pressure loop, which lives in the FULL-FILM region. Taking a raw
    # argmin over the midplane would return an arbitrary cell of the p = p_cav
    # plateau instead. Restricted to full-film cells it is the rupture point
    # itself on any complementarity path -- see README.md 2.5.
    full_mid = {t: p_mid[t] for t in p_mid if not by_theta_mid.get(t, False)}
    theta_p_min = min(full_mid, key=lambda t: full_mid[t]) if full_mid else None
    p_min_full_mid = min(full_mid.values()) if full_mid else None

    m = {
        "n_cells": len(fc),
        "cav_fraction": sum(1 for c in cav if c) / len(cav),
        "p_max": max(p),
        "p_min": min(p),
        "p_min_midplane": min(p_mid.values()) if p_mid else None,
        "p_min_fullfilm_midplane": p_min_full_mid,
        "theta_0_deg": theta_p_min,
        "alpha_gas_max": max(alpha),
        "vapor_deficit_max": max(max(0.0, 1.0-fc[k]-alpha[k]) for k in range(len(fc))),
        "h_min": min(h),
        "z_mid": z_mid,
    }
    for label, arc in (("any_z", arc_any), ("midplane", arc_mid)):
        if arc is None:
            m[f"theta_st_{label}"] = None
            m[f"theta_end_{label}"] = None
            m[f"extent_{label}"] = 0.0
        else:
            m[f"theta_st_{label}"] = arc[0]
            m[f"theta_end_{label}"] = arc[1]
            m[f"extent_{label}"] = arc[2]
    if arc_any is not None:
        ws = [width.get(t, 0.0) for t in ths if by_theta_any[t]]
        m["axial_width_max_m"] = max(ws) if ws else 0.0
        # width at the angular midpoint of the cavity
        mid_ang = _wrap(arc_any[0] + 0.5 * arc_any[2])
        t_near = min(ths, key=lambda t: min(abs(t - mid_ang), 360.0 - abs(t - mid_ang)))
        m["axial_width_at_cavity_mid_m"] = width.get(t_near, 0.0)
    return m


def midplane_profile(run_dir, name="p"):
    path = final_pvts(run_dir)
    if path is None:
        return None
    a = read_pvts(path)
    if name not in a:
        return None
    theta, z, v = a["theta_rad"][0], a["z_m"][0], a[name][0]
    z_sorted = sorted(set(z))
    z_mid = z_sorted[len(z_sorted) // 2]
    rows = [(_wrap(math.degrees(t)), vv) for t, zv, vv in zip(theta, z, v) if zv == z_mid]
    rows.sort()
    return [r[0] for r in rows], [r[1] for r in rows]


def cavity_boundary(run_dir):
    """(theta_deg, z_m) traces of the cavitated cells, for the cavity map figure."""
    path = final_pvts(run_dir)
    if path is None:
        return None
    a = read_pvts(path)
    fc = a["film_content"][0]
    alpha = a["alpha_gas"][0] if "alpha_gas" in a else [0.0] * len(fc)
    theta = [_wrap(math.degrees(t)) for t in a["theta_rad"][0]]
    z = a["z_m"][0]
    pts = [(theta[k], z[k]) for k in range(len(fc))
           if fc[k] < 1.0 - CEILING_TOL]
    return pts


def diagnostics_tail(run_dir):
    csv_path = Path(run_dir) / "results" / "diagnostics.csv"
    if not csv_path.exists():
        return {}
    lines = [l for l in csv_path.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        return {}
    return dict(zip(lines[0].split(","), lines[-1].split(",")))


def coupling_health(run_dir, log_name="run.log"):
    path = Path(run_dir) / log_name
    if not path.exists():
        return None
    last = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Coupling:" in line:
            last = line
    if last is None:
        return None
    # Strip ANSI colour codes properly. The inherited helper used
    # rstrip("\x1b[0m"), which strips the characters {ESC,[,0,m} individually and
    # therefore eats the value itself whenever the residual is exactly "0".
    body = re.sub(r"\x1b\[[0-9;]*m", "", last).split("Coupling:", 1)[1]
    try:
        iters = int(body.split("outer iters")[0].strip())
        r_p = float(body.split("R_p=")[1].split()[0])
        r_t = float(body.split("R_T=")[1].split()[0])
    except (IndexError, ValueError):
        return None
    return {"iters": iters, "R_p": r_p, "R_T": r_t}
