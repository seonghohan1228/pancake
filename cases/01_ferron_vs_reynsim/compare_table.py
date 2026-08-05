"""Build the 4000 rpm / 6 kN comparison table from data already on disk.

No solver is run.  Sources:
  experiment  cases/ref_data/ferron1983_fig6_4000rpm_6kN_pressure.csv
  ReynSim     reynsim_data/pressure_distribution_4000rpm_6000N.csv (+ locus, flagged)
  pancake     cases/01_ferron_prelim/equilibrium.json (isothermal FBNS, mu40/mu55)
"""

import csv
import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR))
import solver_io as S
import reynsim_io as R

CASES = CASE_DIR.parent
REF = CASES / "ref_data"
PRELIM = CASES / "01_ferron_prelim"
P_AMB = 101325.0
GROOVE_DEG = 90.0
CLEARANCE = 1.52e-4


def peak(xs, ys):
    i = max(range(len(ys)), key=lambda k: ys[k])
    return ys[i], xs[i]


def read_ref():
    xs, ys = [], []
    with open(REF / "ferron1983_fig6_4000rpm_6kN_pressure.csv", encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if row and not row[0].lstrip().startswith("#"):
                xs.append(float(row[0]))
                ys.append(float(row[1]))
    return xs, ys


def parabola_peak(xs, ys):
    """Three-point parabolic refinement around the discrete maximum."""
    i = max(range(len(ys)), key=lambda k: ys[k])
    if i == 0 or i == len(ys) - 1:
        return ys[i], xs[i]
    x0, x1, x2 = xs[i - 1], xs[i], xs[i + 1]
    y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d
    c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1
         + x0 * x1 * (x0 - x1) * y2) / d
    if a >= 0:
        return y1, x1
    xv = -b / (2 * a)
    return a * xv * xv + b * xv + c, xv


out = {}

# --- experiment -------------------------------------------------------------
xe, ye = read_ref()
pv, px = peak(xe, ye)
out["experiment"] = {
    "source": "Ferron 1983 Fig. 6, digitised (cases/ref_data)",
    "p_peak_MPa_gauge": pv / 1e6, "p_peak_theta_deg": px,
    "note": "peak is the largest of 15 digitised points, not a curve maximum",
    "eps": None, "attitude_deg": None,
}

# --- ReynSim ----------------------------------------------------------------
xr, pr = R.midplane_pressure()
pv, px = parabola_peak(xr, pr)
raw_t, raw_p = R.raw_midplane_pressure()
out["reynsim"] = {
    "source": "reynsim_data/pressure_distribution_4000rpm_6000N.csv",
    "p_peak_MPa_gauge": pv / 1e6,
    "p_peak_MPa_absolute": (pv + P_AMB) / 1e6,
    "p_peak_theta_deg": px,
    "p_raw_max_Pa": max(raw_p), "p_raw_min_Pa": min(raw_p),
    "groove_cell_Pa": max(v for t, v in zip(raw_t, raw_p) if 85 < t < 95),
    "n_theta": len(raw_t),
}

# ReynSim journal-centre locus found in the same KSME folder.  The filename and
# every script that reads it are silent about the operating point, so this is
# reported as UNVERIFIED.
loc = CASE_DIR / "reynsim_data" / "journal_center_locus.csv"
if loc.exists():
    rows = [r for r in csv.DictReader(open(loc, encoding="utf-8-sig"))]
    x, y = float(rows[-1]["avg(x)"]), float(rows[-1]["avg(y)"])
    out["reynsim"]["locus_unverified"] = {
        "x_m": x, "y_m": y, "e_m": math.hypot(x, y),
        "eps_at_C_152um": math.hypot(x, y) / CLEARANCE,
        "angle_from_plus_x_deg": math.degrees(math.atan2(y, x)),
        "note": "KSME/journal_center_locus.csv; operating point NOT stated anywhere",
    }

# --- pancake, 01-prelim isothermal FBNS ------------------------------------
eq = json.loads((PRELIM / "equilibrium.json").read_text())
for key in ("mu40", "mu55"):
    rec = eq[key]
    theta, p = S.midplane_profile(PRELIM / rec["dir"], "p")
    xs = [(t - GROOVE_DEG) % 360.0 for t in theta]
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in order]
    ys = [p[i] - P_AMB for i in order]
    pv, px = parabola_peak(xs, ys)
    out[f"pancake_{key}"] = {
        "source": f"cases/01_ferron_prelim/{rec['dir']} (FBNS, ISOTHERMAL)",
        "mu_Pa_s": rec["mu"],
        "p_peak_midplane_MPa_gauge": pv / 1e6, "p_peak_theta_deg": px,
        "p_max_field_MPa_gauge": (rec["metrics"]["p_max"] - P_AMB) / 1e6,
        "eps": rec["eps"], "attitude_deg": rec["psi"],
        "h_min_um": rec["metrics"]["h_min"] * 1e6,
        "W_N": rec["metrics"]["W"],
        "cav_area_frac": rec["metrics"]["cav_frac"],
        "friction_torque_Nm": rec["metrics"]["friction_torque"],
    }

(CASE_DIR / "comparison.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
