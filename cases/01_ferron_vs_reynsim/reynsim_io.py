"""Read the previous-generation MATLAB/ReynSim exports and put them in our frame.

Every convention decision below is derived from the ReynSim scripts or from an
independently known physical constant.  Nothing here is fitted to the Ferron
experiment.

Units
-----
`pressure_distribution_4000rpm_6000N.csv` column `p` is in Pa, ABSOLUTE.
  * `p_plot.m` divides the column by 1e6 and labels the axis MPa, so the column
    is Pa and not kPa (max = 1.6087e6 -> 1.61 MPa, the right order for this
    bearing).
  * The single groove cell reads p = 171320 Pa.  Ferron's supply is 0.07 MPa
    gauge, and 101325 + 70000 = 171325 Pa.  The match to 5 Pa fixes the datum
    as ABSOLUTE, independently of the experiment.
  * The cavitated plateau is 100 Pa absolute, i.e. ReynSim's rupture pressure
    is ~0 absolute (our solver uses p_cav = 1.0e5 Pa absolute).  Converting to
    gauge therefore puts the ReynSim cavitated zone at -0.101 MPa.  That is a
    real model difference, not a units artefact.

Theta origin
------------
ReynSim theta = (degrees from groove centre) + 90.
  * `p_plot.m` overlays the Ferron experiment after `mod(theta_exp + 90, 360)`,
    and its `theta_exp` array is the same groove-referenced digitisation as
    cases/ref_data/ferron1983_fig6_4000rpm_6kN_pressure.csv (15/15 points agree
    to <1 deg and <0.03 MPa).  So their frame is ours + 90 deg.
  * Confirmed structurally: the groove cell sits at ReynSim theta = 89.0 deg,
    i.e. their groove centre is at theta = 90 deg, exactly like our mesh.
  * `pressure_contour.m` applies the same +90 shift to a second Ferron file,
    so the convention is consistent across their scripts.
This is the identical mapping our own plots already use, x = (theta - 90) % 360.
"""

import csv
import math
from pathlib import Path

DATA = Path(__file__).resolve().parent / "reynsim_data"
P_AMB = 101325.0
GROOVE_DEG = 90.0


def _rows(path):
    with open(path, encoding="utf-8-sig") as fh:
        for row in csv.reader(fh):
            if row and not row[0].lstrip().startswith("#"):
                yield row


def midplane_pressure(name="pressure_distribution_4000rpm_6000N.csv"):
    """ReynSim mid-plane p(theta) -> (deg from groove, gauge Pa), sorted."""
    rows = list(_rows(DATA / name))
    header = rows[0]
    i_p, i_t = header.index("p"), header.index("theta")
    pts = []
    for r in rows[1:]:
        theta_deg = float(r[i_t]) * 180.0 / math.pi
        pts.append(((theta_deg - GROOVE_DEG) % 360.0, float(r[i_p]) - P_AMB))
    pts.sort()
    return [p[0] for p in pts], [p[1] for p in pts]


def raw_midplane_pressure(name="pressure_distribution_4000rpm_6000N.csv"):
    """Same file, untouched: (ReynSim theta deg, absolute Pa)."""
    rows = list(_rows(DATA / name))
    header = rows[0]
    i_p, i_t = header.index("p"), header.index("theta")
    pts = sorted((float(r[i_t]) * 180.0 / math.pi, float(r[i_p])) for r in rows[1:])
    return [p[0] for p in pts], [p[1] for p in pts]


# --- eccentricity vs load ---------------------------------------------------
# Inline arrays transcribed verbatim from the two ReynSim plotting scripts kept
# in reynsim_data/.  They disagree with each other; both are reported.

# KSFM/validation.m (2025-06-16).  Speed NOT stated anywhere in the script.
VALIDATION_M = {
    "load_N": [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0],
    "eps_num": [0.29, 0.42, 0.51, 0.57, 0.61, 0.65, 0.68, 0.70, 0.72],
    "load_exp_kN": [2.001499678, 3.001444845, 4.000304322, 4.99854693, 6.002785508,
                    7.000460596, 7.99801231, 9.001584669, 10.0049843],
    "eps_exp": [0.448205568, 0.535694636, 0.595467244, 0.639491864, 0.676586931,
                0.706123401, 0.732510274, 0.752597513, 0.768275316],
}

# LGEP/plot_data.m (2025-03-05).  An earlier ReynSim vintage.  Its experimental
# points reproduce our own Fig. 5 (C = 145 um) digitisation to <0.01 in eps; its
# numerical points are a different curve from validation.m's.
PLOT_DATA_M = {
    "load_N": [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0],
    "eps_num": [0.421, 0.548, 0.624, 0.676, 0.713, 0.742, 0.765, 0.784, 0.799],
    "load_exp_kN": [1.973871065, 2.973904431, 3.989421271, 5.009676943, 6.000911384,
                    7.008550739, 7.982026902, 9.0235415, 10.00816543],
    "eps_exp": [0.493289032, 0.595543222, 0.664413147, 0.713713817, 0.754386962,
                0.78009293, 0.808683151, 0.825749527, 0.848582213],
}
