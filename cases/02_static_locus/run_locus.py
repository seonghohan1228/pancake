"""Case 02: static eccentricity locus of a fully flooded plain journal bearing.

For each L/D and each eps = 0.1 ... 0.9 one STATIC + STEADY_STATE + FBNS run is
made at fixed attitude_angle_deg = -90 (minimum film at theta = 270 deg). The
load magnitude W, its direction, and hence the Sommerfeld number are outputs.

Feed variants:
  "groove"  (default, the reported locus) - a 1 deg ambient-pressure groove at
            theta = 90 deg, which is exactly the maximum-film-thickness line for
            attitude_angle_deg = -90. This is the fallback the task allowed and
            it is the configuration that actually reproduces a fully flooded
            film: the groove pins film_content = 1 at h_max, so the complete
            film starts where Raimondi & Boyd assume it starts.
  "noinlet" - no inlet keyword at all. The solver ACCEPTS this config, but the
            mass-conserving FBNS film then has no feed except axial diffusion
            through the ends, and it settles into a starved state with ~85 % of
            the domain cavitated at every eccentricity, carrying roughly 1/20 of
            the flooded load. Kept for the record, not for the R&B comparison.

Attitude angle convention (classical, journal-relative):
    the bearing moves and the journal is fixed, so the journal-relative line of
    centres points at psi = attitude_angle_deg and the equivalent journal load
    direction is the pressure-force direction atan2(Fpy, Fpx). The reported
    attitude angle is  phi = wrap180(psi - angle(Fp)), positive in the +theta
    (rotation) direction, which is the Raimondi & Boyd / Someya definition.

Sommerfeld number:  S = mu * N_revs * L * D * (R/c)^2 / W,  N_revs = omega/2pi.

Side flow is the total volumetric liquid outflow through the two axial ends,
recomputed from the final frame with the same FBNS Dirichlet boundary link the
solver assembles:  q_face = (h^3 / 12 mu) * R dtheta / (dz/2) * (p_cell - p_bc).

Usage: python run_locus.py [LD1|LD05|all] [groove|noinlet]
"""

import csv
import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_prelim"))
import solver_io as S

PROGRESS = CASE_DIR / "progress.log"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text()

R = 0.05
CLEARANCE = 1.52e-4
MU = 0.0277
RHO = 860.0
OMEGA = 418.879020479
PSI_DEG = -90.0
N_REVS = OMEGA / (2.0 * math.pi)
EPS_LIST = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

GEOMETRIES = {"LD1": {"L": 0.1, "n_z": 256}, "LD05": {"L": 0.05, "n_z": 128}}

NO_INLET_MARK = "# No inlet groove: fully flooded, fed only through the axial ends."
GROOVE_LINE = "inlet_groove = 90 1 101325"


def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0


def side_flow(run_dir, L, n_z, p_bc=101325.0):
    """(total outflow, net flux) through both axial ends, m^3/s."""
    a = S._read_pvts(S.final_pvts(run_dir))
    theta, z, p, h = a["theta_rad"][0], a["z_m"][0], a["p"][0], a["h"][0]
    z_lo, z_hi = min(z), max(z)
    d_theta = 2.0 * math.pi / 720
    d_z = L / n_z
    out = net = 0.0
    for t, zv, pv, hv in zip(theta, z, p, h):
        if zv not in (z_lo, z_hi):
            continue
        q = (hv ** 3 / (12.0 * MU)) * (R * d_theta) / (0.5 * d_z) * (pv - p_bc)
        net += q
        out += max(q, 0.0)
    return out, net


def run_geometry(key, variant="groove"):
    geom = GEOMETRIES[key]
    suffix = "" if variant == "groove" else f"_{variant}"
    rows = []
    for eps in EPS_LIST:
        run_dir = CASE_DIR / f"{key}{suffix}" / f"eps_{int(round(eps*100)):02d}"
        run_dir.mkdir(parents=True, exist_ok=True)
        cfg = TEMPLATE.format(e=repr(eps * CLEARANCE), L=repr(geom["L"]), n_z=geom["n_z"])
        if variant == "groove":
            cfg = cfg.replace(NO_INLET_MARK, GROOVE_LINE)
        (run_dir / "config.txt").write_text(cfg)
        rc, wall = S.run_case(run_dir)
        if rc != 0:
            S.progress(PROGRESS, f"{key}{suffix}/eps={eps:.1f} FAILED rc={rc}")
            continue
        m = S.frame_metrics(run_dir)
        W = m["W"]
        phi = wrap180(PSI_DEG - m["force_angle_deg"])
        somm = MU * N_REVS * geom["L"] * (2.0 * R) * (R / CLEARANCE) ** 2 / W
        q_out, q_net = side_flow(run_dir, geom["L"], geom["n_z"])
        rows.append({
            "eps": eps, "S": somm, "attitude_deg": phi, "W": W,
            "p_max": m["p_max"], "friction_torque": m["friction_torque"],
            "side_flow": q_out, "side_flow_net": q_net,
            "Fx": m["Fx"], "Fy": m["Fy"], "h_min": m["h_min"],
            "cav_frac": m["cav_frac"], "wall_s": round(wall, 1),
        })
        S.progress(PROGRESS, f"{key}{suffix}/eps={eps:.1f} W={W:.1f} N S={somm:.4f} "
                             f"phi={phi:.2f} deg p_max={m['p_max']:.4g} Pa "
                             f"Q={q_out:.4g} m3/s ({wall:.0f} s)")
    return rows


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    variant = sys.argv[2] if len(sys.argv) > 2 else "groove"
    keys = list(GEOMETRIES) if which == "all" else [which]
    suffix = "" if variant == "groove" else f"_{variant}"
    all_rows = {}
    for key in keys:
        rows = run_geometry(key, variant)
        all_rows[key] = rows
        cols = ["eps", "S", "attitude_deg", "W", "p_max", "friction_torque",
                "side_flow", "side_flow_net", "Fx", "Fy", "h_min", "cav_frac"]
        with open(CASE_DIR / f"locus_{key}{suffix}.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    (CASE_DIR / f"locus{suffix}.json").write_text(json.dumps(all_rows, indent=2))


if __name__ == "__main__":
    main()
