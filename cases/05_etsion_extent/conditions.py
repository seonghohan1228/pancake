"""The eight Etsion & Ludwig (1982) conditions, and the config generator.

Single source of truth for the case-05 run matrix. Imported by run_case05.py
(which runs the solver) and by extent.py / the plotting scripts (which only need
the labels and the reference values), so the matrix is defined exactly once.

Ambient pressures are printed in the paper as GAUGE (differential transducer
referenced to atmosphere, "2 psig" ... "8 psig"). The solver boundary condition
needs ABSOLUTE, hence P_ATM is added here and nowhere else.
"""

import math
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
REF_CSV = CASE_DIR.parent / "ref_data" / "etsion1982_table2_extent.csv"

P_ATM = 101325.0  # Pa

# Cavitation datum: atmospheric air-saturation pressure, ABSOLUTE, identical for
# every ambient level. README.md section 2.5 gives the three reasons.
P_CAV = 101325.0

# The oil's own vapour pressure at 38 C: 0.1 mmHg, Table 1 of the paper. VERIFIED.
# Used only by the `vapour_datum` leg, which is the contrast run that shows what
# treating this bearing as VAPOROUS cavitation would predict -- and, incidentally,
# the only configuration in which pancake's TWO_PHASE outgassing does anything at
# all in steady state (see SUMMARY: with p_cav = p_sat the release band is empty).
P_VAPOUR = 13.33  # Pa absolute

# UNVERIFIED (not in the paper): paraffinic mineral oil near 40 C, 1 bar.
# Read by the JFO/Elrod leg only; FBNS ignores it (PHYSICS.md 7.3.8).
BULK_MODULUS = 1.7e9

N_THETA = 720   # frozen by case 00 (G6)
N_Z = 192

SPEEDS_RPM = (1840, 3000)
PS_KPA_GAUGE = (13.6, 27.2, 40.8, 54.4)

MODELS = {
    "fbns": {"cavitation_model": "FBNS", "fluid_property_model": "TWO_PHASE"},
    "jfo": {"cavitation_model": "JFO", "fluid_property_model": "TWO_PHASE"},
}


def omega_of(rpm):
    return 2.0 * math.pi * rpm / 60.0


def tag(rpm, ps_kpa):
    """Directory name for one condition, e.g. n1840_ps136."""
    return f"n{rpm}_ps{int(round(ps_kpa * 10)):03d}"


def conditions():
    """The eight (rpm, p_s gauge kPa) pairs, in Table 2 order."""
    return [(rpm, ps) for rpm in SPEEDS_RPM for ps in PS_KPA_GAUGE]


def read_reference():
    """Table 2 as a dict keyed by tag(). Values are floats, angles in degrees."""
    rows = {}
    with open(REF_CSV, "r", encoding="utf-8") as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split(",")
                continue
            vals = line.split(",")
            rec = dict(zip(header, vals))
            rpm = int(float(rec["N_rpm"]))
            ps = float(rec["p_s_kPa_gauge"])
            rows[tag(rpm, ps)] = {
                "rpm": rpm,
                "ps_kPa": ps,
                "theta_0": float(rec["theta_0_deg"]),
                "p_min_kPa_gauge": float(rec["p_min_kPa_gauge"]),
                "theta_st": float(rec["theta_st_deg"]),
                "theta_end": float(rec["theta_end_deg"]),
                "extent": float(rec["extent_deg"]),
            }
    return rows


def render_config(rpm, ps_kpa, model, *, n_theta=N_THETA, n_z=N_Z,
                  p_cav=P_CAV, bulk_modulus=BULK_MODULUS, inlet_line="",
                  cavitation_threshold="LOCAL_PSAT"):
    """Render config_template.txt for one run. Returns the config text.

    `cavitation_threshold` matters more than it looks. Under LOCAL_PSAT the solver
    uses max(p_cav, p_sat(c_d_init, T_init)) as the complementarity datum
    (PHYSICS.md 7.3.8), so a p_cav BELOW the air-saturation pressure is silently
    raised back to it and the outgassing band closes. SCALAR_PCAV is required to
    keep a datum below p_sat, which is what the `vapour_datum` leg needs.
    """
    template = (CASE_DIR / "config_template.txt").read_text(encoding="utf-8")
    spec = MODELS[model]
    template = template.replace("cavitation_threshold = LOCAL_PSAT",
                                f"cavitation_threshold = {cavitation_threshold}")
    return template.format(
        n_theta=n_theta,
        n_z=n_z,
        rpm=rpm,
        omega=repr(omega_of(rpm)),
        p_cav=repr(p_cav),
        bulk_modulus=repr(bulk_modulus),
        p_amb_abs=repr(P_ATM + ps_kpa * 1000.0),
        fluid_property_model=spec["fluid_property_model"],
        cavitation_model=spec["cavitation_model"],
        inlet_line=inlet_line,
    )
