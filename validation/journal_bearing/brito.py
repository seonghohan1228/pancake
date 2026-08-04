"""Brito, Bouyer, Fillon & Miranda (2006) two-axial-groove THD benchmark.

    Brito, F. P., Bouyer, J., Fillon, M. and Miranda, A. S., "Experimental
    Investigation of the Influence of Supply Temperature and Supply Pressure on
    the Performance of a Two-Axial Groove Hydrodynamic Journal Bearing",
    ASME J. Tribology 129(1), 2006, DOI 10.1115/1.2401206.

Open access at the Universidade do Minho repository:
https://repositorium.sdum.uminho.pt/bitstreams/3ec650eb-8054-416f-8efe-bedb3f79c816/download
(retrieved 2026-07-31; the /bitstream/1822/6136/... path 301-redirects to a
localhost URL and is broken).

Why this benchmark, alongside Ferron (see `ferron.py`)
------------------------------------------------------
Ferron is the only dataset in this validation set that tests temperature at all,
so the energy-equation verdict rests on a single experiment.  Brito is a second,
independent THD dataset on the *same rig lineage* -- the paper states it is the
LMS-Poitiers rig used by Ferron, Frene and Boncompain -- and it adds three things
Ferron does not have:

  * **oil-shaft interface temperature**, measured through a rotary joint.  That
    is the quantity pancake's fixed-wall energy model structurally cannot
    produce, and the suspected cause of the Ferron temperature-peak failure.
  * a **supply-temperature and supply-pressure sweep**, which is the axis where
    the Ferron comparison is weakest.
  * an **unambiguous clearance**: quoted diametrally, with a tolerance and a
    reference temperature.

ANGULAR CONVENTION -- differs from Ferron, check before comparing
-----------------------------------------------------------------
Brito's abscissa is "angle from centre of upstream groove", 0 to 360 deg, and
the paper states: "The angle origin is located at the center of the upstream
groove: the two feeding grooves here are called the upstream and downstream
grooves, depending on whether the groove is located upstream or downstream of
the minimum film thickness."  With the two grooves at +-90 deg to the load line
(abstract), the layout is therefore

    upstream groove at 0 deg,  load line at 90 deg,  downstream groove at 180 deg

with the angle increasing in the direction of rotation.  Ferron's X is measured
from the *anti-groove* position with the load at X = 0, so the two conventions
are offset and must not be interchanged.
"""
import math

# --- Table 1, p.4 -----------------------------------------------------------
D = 100e-3                 # bush inner diameter (nominal) [m]
R = 0.5 * D
L = 80e-3                  # bush length b [m]
C_DIAMETRAL = 171e-6       # diametral clearance at 20 C [m]
C = 0.5 * C_DIAMETRAL      # radial clearance, 85.5 um
GROOVE_AXIAL_LENGTH = 70e-3        # groove length a [m], 88% of the bush width
GROOVE_CIRC_EXTENSION = 16e-3      # circumferential extension w [m]
MU_40 = 0.0293             # Pa.s
MU_70 = 0.0111             # Pa.s
K_OIL = 0.13               # W/(m.K)
K_BUSH = 50.0              # W/(m.K) -- given, but pancake has no conjugate model
RHO = 870.0                # kg/m3
CP = 2000.0                # J/(kg.K)
RHO_CP = RHO * CP          # 1.74e6
T_AMBIENT_RANGE_C = (30.0, 55.0)

L_OVER_D = L / D                      # 0.800
C_OVER_R = C / R                      # 1.71e-3
GROOVE_ANGLE_DEG = GROOVE_CIRC_EXTENSION / (math.pi * D) * 360.0   # 18.33 deg

# Operating grid actually tested
SPEEDS_RPM = (1000.0, 3000.0, 4000.0)
LOADS_N = (2000.0, 10000.0)
SUPPLY_PRESSURES_PA = (70e3, 140e3, 210e3)     # gauge
SUPPLY_TEMPERATURES_C = (35.0, 40.0, 50.0)

# Geometry of the groove layout, from the convention note above.
LOAD_ANGLE_FROM_UPSTREAM_GROOVE_DEG = 90.0
GROOVE_ANGLES_FROM_UPSTREAM_DEG = (0.0, 180.0)


def andrade_E_mu():
    """Andrade slope from the two tabulated viscosity points.

    Unlike Ferron -- whose viscosity is given as an exponential-in-temperature
    law that has to be refitted onto pancake's Andrade form with a residual --
    Brito quotes mu at two temperatures, so the Andrade constant is determined
    exactly and reproduces both points to five digits.
    """
    T40, T70 = 313.15, 343.15
    return math.log(MU_40 / MU_70) / (1.0 / T40 - 1.0 / T70)      # 3476.8 K


E_MU = andrade_E_mu()

# --- Measured data ----------------------------------------------------------
# Fig. 5(c): maximum hydrodynamic pressure against supply pressure, W = 10 kN,
# Tf = 40 C.  Digitised from the high-DPI render; the paper states repeatability
# on maximum pressure of better than +-0.01 MPa in most cases, and the reading
# uncertainty here is about +-0.02 MPa.
#   {speed rpm: [(supply pressure Pa gauge, p_max MPa), ...]}
PMAX_VS_SUPPLY = {
    1000.0: [(70e3, 2.94), (140e3, 2.87), (210e3, 2.81)],
    3000.0: [(70e3, 2.75), (140e3, 2.66), (210e3, 2.55)],
    4000.0: [(70e3, 2.58), (140e3, 2.56), (210e3, 2.55)],
}
PMAX_DIGITISE_ERR_MPA = 0.02

# Fig. 6: locus of the shaft centre, W = 10 kN, Tf = 40 C.
#
# THIS IS A POLAR CHART, not a Cartesian one.  The origin (zero eccentricity) is
# the top-left corner; eccentricity is the RADIUS from it and attitude angle is
# the angle from the downward vertical, with the constant-eccentricity gridlines
# drawn as arcs and the constant-attitude lines as rays.  Reading eccentricity
# off the vertical axis and attitude off the horizontal axis -- as if it were
# Cartesian -- gives eccentricities about 0.15 too low and attitude angles about
# 15 deg too low, and produces the physically impossible result that this
# two-groove bearing carries its load at a *lower* eccentricity than a plain
# bearing would.
#
# Digitised numerically: the plot box was located from the image rules, giving
# the origin and 2193 px per unit eccentricity, and each marker position was
# taken as the centroid of dark pixels in a window around it.
#   {speed rpm: [(supply pressure Pa gauge, eccentricity ratio, attitude deg)]}
ECC_ATTITUDE = {
    1000.0: [(70e3, 0.761, 35.3), (210e3, 0.754, 36.6)],
    3000.0: [(70e3, 0.576, 39.4), (210e3, 0.563, 43.5)],
    4000.0: [(70e3, 0.558, 43.5), (210e3, 0.526, 45.8)],
}
# The 4000/70 and 3000/210 markers nearly touch, so their centroids can be
# contaminated by one another; those two carry the larger uncertainty.
ECC_DIGITISE_ERR = 0.02
ATTITUDE_DIGITISE_ERR_DEG = 2.0

# Independent re-digitisation of Fig. 6 by Seongho Han (digiplotter, 2026-08-04),
# done blind to the values above. Same layout as ECC_ATTITUDE.
#   Agreement with the polar reading: |d eps| <= 0.013 (mean 0.005), |d phi| <=
#   1.0 deg -- inside the stated uncertainties at every one of the six points.
#   The earlier CARTESIAN reading of this polar chart sat ~0.15 lower in eps and
#   is hereby REFUTED; the polar values stand and the load comparison is
#   scoreable again. Keep both digitisations: their spread is the measured
#   reading uncertainty, which no single reading provides.
ECC_ATTITUDE_INDEPENDENT = {
    1000.0: [(70e3, 0.7704, 34.26), (210e3, 0.7533, 36.28)],
    3000.0: [(70e3, 0.5769, 39.05), (210e3, 0.5658, 43.12)],
    4000.0: [(70e3, 0.5447, 42.58), (210e3, 0.5264, 45.54)],
}

# Consistency check on the reading, from the paper's own text: "the eccentricity
# decreases and the attitude angle increases with rises in value of those
# parameters [feeding pressure and rotational speed]".  The values above satisfy
# both monotonicities in both variables; the Cartesian misreading satisfied the
# eccentricity trend but inverted nothing about attitude, so it did not fail this
# check -- only the physical argument caught it.

# Fig. 8: oil flow rate against supply pressure, Tf = 40 C.  Flow meter accuracy
# is +-0.05 l/min; reading uncertainty about +-0.15 l/min.
#   {(speed rpm, load N): [(supply pressure Pa gauge, Q l/min), ...]}
FLOW_RATE = {
    (1000.0, 10000.0): [(70e3, 2.10), (140e3, 3.20), (210e3, 4.25)],
    (4000.0, 10000.0): [(70e3, 4.20), (140e3, 5.35), (210e3, 6.50)],
    (1000.0, 2000.0):  [(70e3, 1.70), (140e3, 2.55), (210e3, 3.35)],
    (4000.0, 2000.0):  [(70e3, 2.65), (140e3, 3.60), (210e3, 4.60)],
}
FLOW_DIGITISE_ERR_LPM = 0.15

# Fig. 5(a): mid-plane pressure profile, W = 10 kN, N = 4000 rpm, Tf = 40 C,
# against angle from the centre of the upstream groove.  Reading uncertainty
# about +-4 deg in angle and +-0.05 MPa in pressure.
#   {supply pressure Pa gauge: [(angle deg, p MPa), ...]}
PRESSURE_PROFILE_10kN_4000rpm = {
    210e3: [(20, 0.61), (50, 1.13), (70, 1.60), (85, 2.15), (100, 2.56),
            (112, 2.49), (130, 1.82), (148, 0.78), (168, 0.00), (200, -0.05),
            (250, 0.00), (270, 0.02), (310, 0.00), (350, 0.05)],
    70e3:  [(20, 0.53), (50, 1.12), (70, 1.63), (85, 2.19), (100, 2.61),
            (112, 2.48), (130, 1.77), (148, 0.72), (168, 0.00), (200, -0.06),
            (250, -0.02), (270, 0.00), (310, 0.00), (350, 0.03)],
}
PROFILE_DIGITISE_ERR_MPA = 0.05
PROFILE_DIGITISE_ERR_DEG = 4.0


def wrap_from_upstream_groove(theta_deg, upstream_groove_deg):
    """Map a pancake circumferential angle onto Brito's abscissa: degrees from
    the centre of the upstream groove, 0 to 360, increasing with rotation."""
    return (theta_deg - upstream_groove_deg) % 360.0


if __name__ == "__main__":
    print(f"Brito 2006 rig: D = {D*1e3:g} mm, L = {L*1e3:g} mm (L/D = {L_OVER_D:g}), "
          f"c = {C*1e6:.1f} um radial (Cd = {C_DIAMETRAL*1e6:g} um at 20 C), "
          f"c/R = {C_OVER_R:.3e}")
    print(f"two grooves at 0 and 180 deg from the upstream groove, "
          f"{GROOVE_ANGLE_DEG:.2f} deg wide, "
          f"{GROOVE_AXIAL_LENGTH/L*100:.0f}% of the bush width; load at "
          f"{LOAD_ANGLE_FROM_UPSTREAM_GROOVE_DEG:g} deg")
    print(f"oil: mu(40 C) = {MU_40:g}, mu(70 C) = {MU_70:g} Pa.s -> "
          f"E_mu = {E_MU:.1f} K (exact two-point fit)")
    print(f"     rho = {RHO:g}, rho_cp = {RHO_CP:.4g}, k = {K_OIL:g} W/(m.K)")
    print(f"conditions: N = {SPEEDS_RPM} rpm, W = {LOADS_N} N, "
          f"Pf = {tuple(p/1e3 for p in SUPPLY_PRESSURES_PA)} kPa gauge, "
          f"Tf = {SUPPLY_TEMPERATURES_C} C")
