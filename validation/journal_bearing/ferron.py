"""Ferron, Frene & Boncompain (1983) reference data and property conversions.

Primary reference (paywalled, ASME):

    Ferron, J., Frene, J. and Boncompain, R. (1983), "A Study of the
    Thermohydrodynamic Performance of a Plain Journal Bearing: Comparison
    Between Theory and Experiments", ASME J. Lubrication Technology 105(3),
    422-428, DOI 10.1115/1.3254632.

All numbers below are taken from the open-access reproduction

    Alakhramsing, S. S., van Ostayen, R. A. J. and Eling, R. (2015),
    "Thermo-Hydrodynamic Analysis of a Plain Journal Bearing on the Basis of a
    New Mass Conserving Cavitation Algorithm", Lubricants 3(2), 256-280,
    DOI 10.3390/lubricants3020256.

retrieved 2026-07-31 from the Delft University of Technology institutional
repository (resolver.tudelft.nl/uuid:b84a2854-b276-49e7-aa2e-a66bf9724717); the
MDPI copy is not reachable from this environment.  The per-number provenance is
recorded against each entry below.

NOTE ON THE GEOMETRY.  Secondary sources circulate L/D = 0.55 and C/R = 1.5e-3
for this rig.  Alakhramsing Table 1, which cites Ferron directly, gives
L = 80 mm (L/D = 0.80) and C = 145 um (C/R = 2.9e-3).  The values here follow
Table 1.
"""
import math

# --- Alakhramsing et al. (2015), Table 1, p.273 -- "Geometric details of the
# --- plain journal bearing and properties of the lubricant for the reference
# --- calculations [14]", where [14] is Ferron et al. (1983).
D = 100e-3              # shaft diameter [m]
R = 0.5 * D
L = 80e-3               # axial length of the bearing [m]
C = 145e-6              # radial clearance [m]
RHO = 860.0             # density [kg/m^3]
MU_INLET = 0.0277       # viscosity at the supply temperature [Pa.s]
CP = 2000.0             # heat capacity [J/(kg.K)]
RHO_CP = RHO * CP       # 1.72e6 J/(m^3.K)
K_THERMAL = 0.13        # thermal conductivity [W/(m.K)]
BETA_VISC = 3.3e-2      # temperature-viscosity index [1/K], mu = mu0 exp(-beta dT)
P_SUPPLY_GAUGE = 0.7e5  # supply pressure [Pa gauge]
T_SUPPLY_C = 40.0       # supply temperature [deg C]
T_SUPPLY = T_SUPPLY_C + 273.15
GROOVE_ANGLE_DEG = 18.0

L_OVER_D = L / D                     # 0.80
C_OVER_R = C / R                     # 2.9e-3

# --- Reference operating point of Alakhramsing Figure 9 (caption: w^a_x = 4 kN,
# --- omega_s = 2 krpm), whose panels (b) and (c) plot the Ferron measurements.
REF_LOAD_N = 4000.0
REF_RPM = 2000.0

# --- Measured mid-plane circumferential pressure, digitised from Alakhramsing
# --- Figure 9b (filled squares, "experiment").  X is the circumferential
# --- coordinate with X/pi = +-1 at the supply groove and X increasing in the
# --- direction of shaft rotation.  Digitising uncertainty is about +-0.02 in
# --- X/pi and +-0.03 MPa in p.
FERRON_PRESSURE_MPA = [
    (-0.90, 0.03), (-0.72, 0.01), (-0.50, 0.135), (-0.35, 0.30),
    (-0.17, 0.60), (-0.05, 1.15), (0.08, 1.29), (0.16, 1.21),
    (0.30, 0.42), (0.42, 0.02), (0.50, 0.03), (0.62, 0.02),
    (0.72, 0.03), (0.85, 0.03),
]
PRESSURE_DIGITISE_ERR_MPA = 0.03
PRESSURE_DIGITISE_ERR_XPI = 0.02

# --- Measured bearing-wall temperature, digitised from Alakhramsing Figure 9c
# --- (open circles, "experiment").  Ferron's thermocouples sat 0.5 mm below the
# --- bush surface.  Digitising uncertainty about +-0.3 K.
FERRON_TEMPERATURE_C = [
    (-0.90, 40.5), (-0.80, 40.5), (-0.60, 41.1), (-0.45, 41.3),
    (-0.32, 41.8), (-0.20, 42.8), (-0.13, 43.0), (0.05, 45.1),
    (0.22, 46.8), (0.32, 47.3), (0.35, 48.1), (0.55, 49.1),
    (0.75, 47.1), (0.92, 43.9),
]
TEMPERATURE_DIGITISE_ERR_K = 0.3

# --- Measured eccentricity ratio against applied load, digitised from
# --- Alakhramsing Figure 10a (open squares, "experiment Ferron et al, 1983").
# --- The panel does not label a speed; it sits in the same section as the
# --- 2000 rpm reference case.  Digitising uncertainty about +-0.02 in eps.
FERRON_ECC_VS_LOAD = [
    (2000.0, 0.47), (3000.0, 0.57), (4000.0, 0.67), (5000.0, 0.70),
    (6000.0, 0.74), (7000.0, 0.77), (8000.0, 0.78), (9000.0, 0.82),
]
ECC_DIGITISE_ERR = 0.02


def measured_peak_pressure_mpa():
    """(X/pi, p_max) of the largest digitised pressure point."""
    return max(FERRON_PRESSURE_MPA, key=lambda t: t[1])


def measured_peak_temperature_c():
    """(X/pi, T_max) of the largest digitised temperature point."""
    return max(FERRON_TEMPERATURE_C, key=lambda t: t[1])


def measured_rupture_xpi(threshold_mpa=0.10):
    """X/pi at which the measured pressure has fallen back to the cavity level,
    bracketed by the last point above `threshold_mpa` and the first point below.
    The measurement cannot resolve rupture more sharply than that bracket."""
    pts = sorted(FERRON_PRESSURE_MPA)
    peak_x = measured_peak_pressure_mpa()[0]
    last_hi = max(x for x, p in pts if x >= peak_x and p >= threshold_mpa)
    first_lo = min(x for x, p in pts if x > last_hi)
    return last_hi, first_lo


def andrade_from_exponential(beta=BETA_VISC, T_ref=T_SUPPLY, T_hi=None):
    """Least-squares Andrade slope E_mu [K] equivalent to Ferron's exponential
    viscosity law over the operating temperature range.

    Ferron/Alakhramsing use   mu = mu0 exp(-beta (T - T_ref)),
    pancake uses              mu = mu_ref exp(E_mu (1/T - 1/T_ref)).

    Matching the slope at T_ref alone gives E_mu = beta T_ref^2; fitting over the
    range the film actually spans is slightly better.  Returns (E_mu, max
    relative viscosity error over the range).
    """
    if T_hi is None:
        T_hi = T_ref + 15.0            # the measured film rise is about 9 K
    n = 61
    xs, ys = [], []
    for k in range(n):
        T = T_ref + (T_hi - T_ref) * k / (n - 1)
        xs.append(1.0 / T - 1.0 / T_ref)
        ys.append(-beta * (T - T_ref))
    num = sum(x * y for x, y in zip(xs, ys))
    den = sum(x * x for x in xs)
    E_mu = num / den
    err = max(abs(math.exp(E_mu * x) / math.exp(y) - 1.0) for x, y in zip(xs, ys))
    return E_mu, err


def wrap_xpi(theta_deg, groove_deg):
    """Map a pancake circumferential angle [deg, increasing with rotation] onto
    Ferron's X/pi, which is zero opposite the groove and +-1 at the groove."""
    x = (theta_deg - groove_deg - 180.0) % 360.0
    if x > 180.0:
        x -= 360.0
    return x / 180.0


if __name__ == "__main__":
    E_mu, err = andrade_from_exponential()
    print(f"geometry: D={D*1e3:g} mm, L={L*1e3:g} mm, L/D={L_OVER_D:g}, "
          f"c={C*1e6:g} um, c/R={C_OVER_R:.2e}")
    print(f"oil: rho={RHO:g}, mu({T_SUPPLY_C:g}C)={MU_INLET:g} Pa.s, "
          f"rho_cp={RHO_CP:.3g}, k={K_THERMAL:g}, beta={BETA_VISC:g} 1/K")
    print(f"Andrade refit: E_mu = {E_mu:.1f} K "
          f"(slope-match at T_ref would give {BETA_VISC*T_SUPPLY**2:.1f} K); "
          f"max viscosity error over 40-55 C = {err:.2%}")
    print(f"measured peak pressure    {measured_peak_pressure_mpa()[1]:.2f} MPa "
          f"at X/pi = {measured_peak_pressure_mpa()[0]:+.2f}")
    print(f"measured peak temperature {measured_peak_temperature_c()[1]:.1f} C "
          f"at X/pi = {measured_peak_temperature_c()[0]:+.2f} "
          f"(rise {measured_peak_temperature_c()[1]-T_SUPPLY_C:.1f} K)")
    lo, hi = measured_rupture_xpi()
    print(f"measured rupture bracketed by X/pi = {lo:+.2f} .. {hi:+.2f} "
          f"({(hi-lo)*180:.0f} deg wide)")
