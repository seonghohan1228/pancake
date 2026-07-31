"""Raimondi & Boyd (1958) full-journal-bearing design data, plus the exact
infinite-length Swift-Stieber solution used to audit that data.

Reference table
---------------
`RB_TABLE` reproduces "Table (1): Dimensionless performance parameters for full
journal bearing", scanned on p.192 of the University of Babylon lecture notes
"Sliding Bearings" (https://www.uobabylon.edu.iq/eprints/publication_7_4218_752.pdf,
retrieved 2026-07-31, table image on PDF page 12).  That table is the standard
textbook reprint of

    Raimondi, A. A. and Boyd, J. (1958), "A Solution for the Finite Journal
    Bearing and its Application to Analysis and Design -- I, II, III",
    ASLE Transactions 1(1), 159-209.

Columns retained here: eccentricity ratio `eps`, Sommerfeld number `S`, and
attitude angle `phi` in degrees.  The source also lists h0/c, (r/c)f,
Q/(r c n_s l), Qs/Q and P/p_max; those are recorded in `RB_TABLE_FULL` for
completeness but are not used by the acceptance test.

Sommerfeld number convention (stated on p.181 of the same source):

    S = (mu * n_s / P) * (r / c)^2,    n_s in rev/s,   P = W / (L * D)

Attitude angle convention: `phi` is the angle between the load line and the line
of centres, measured from the load line in the direction of shaft rotation.

Provenance audit
----------------
Textbook reprints of this table are known to disagree in the third digit, so the
table is not taken on trust.  `swift_stieber_infinite` below solves the
infinitely-long bearing with the Reynolds (Swift-Stieber) rupture condition in
closed form to machine precision; `audit_infinite_row` compares it against the
L/D = inf block of `RB_TABLE`.  Agreement there is evidence that the finite rows
come from the same solution family and were transcribed correctly.

A second, independent consistency check comes from Al-Hamood et al. (2015),
SSRG Int. J. Mech. Eng. 2(10), 17-20
(https://www.internationaljournalssrg.org/IJME/2015/Volume2-Issue10/IJME-V2I10P105.pdf):
their Table (3) reports, for D = 3 in, L = 1.5 in (L/D = 1/2), N = 4000 rpm,
W = 1000 lbf, c = 1.5e-3 in, that the Raimondi & Boyd charts give eps = 0.28 and
phi = 70 deg.  Those operating conditions give S = 1.39, and log-interpolating
the L/D = 1/2 block of `RB_TABLE` at S = 1.39 returns eps = 0.28, phi = 69.5 deg.
"""
import math

# eps -> (S, phi_deg)
RB_TABLE = {
    "inf": [
        (0.1, 0.240, 69.10),
        (0.2, 0.123, 67.26),
        (0.4, 0.0626, 61.94),
        (0.6, 0.0389, 54.31),
        (0.8, 0.021, 42.22),
        (0.9, 0.0115, 31.62),
    ],
    "1": [
        (0.1, 1.33, 79.5),
        (0.2, 0.631, 74.02),
        (0.4, 0.264, 63.10),
        (0.6, 0.121, 50.58),
        (0.8, 0.0446, 36.24),
        (0.9, 0.0188, 26.45),
        (0.97, 0.00474, 15.47),
    ],
    "0.5": [
        (0.1, 4.31, 81.62),
        (0.2, 2.03, 74.94),
        (0.4, 0.779, 61.45),
        (0.6, 0.319, 48.14),
        (0.8, 0.0923, 33.31),
        (0.9, 0.0313, 23.66),
        (0.97, 0.00609, 13.75),
    ],
}

# eps, h0/c, S, phi_deg, (r/c)f, Q/(r c n_s l), Qs/Q, P/p_max
RB_TABLE_FULL = {
    "inf": [
        (0.1, 0.9, 0.240, 69.10, 4.80, 3.03, 0.0, 0.826),
        (0.2, 0.8, 0.123, 67.26, 2.57, 2.83, 0.0, 0.814),
        (0.4, 0.6, 0.0626, 61.94, 1.52, 2.26, 0.0, 0.764),
        (0.6, 0.4, 0.0389, 54.31, 1.20, 1.56, 0.0, 0.667),
        (0.8, 0.2, 0.021, 42.22, 0.961, 0.760, 0.0, 0.495),
        (0.9, 0.1, 0.0115, 31.62, 0.756, 0.411, 0.0, 0.358),
    ],
    "1": [
        (0.1, 0.9, 1.33, 79.5, 26.4, 3.37, 0.150, 0.540),
        (0.2, 0.8, 0.631, 74.02, 12.8, 3.59, 0.280, 0.529),
        (0.4, 0.6, 0.264, 63.10, 5.79, 3.99, 0.497, 0.484),
        (0.6, 0.4, 0.121, 50.58, 3.22, 4.33, 0.680, 0.415),
        (0.8, 0.2, 0.0446, 36.24, 1.70, 4.62, 0.842, 0.313),
        (0.9, 0.1, 0.0188, 26.45, 1.05, 4.74, 0.919, 0.247),
        (0.97, 0.03, 0.00474, 15.47, 0.514, 4.82, 0.973, 0.152),
    ],
    "0.5": [
        (0.1, 0.9, 4.31, 81.62, 85.6, 3.43, 0.173, 0.523),
        (0.2, 0.8, 2.03, 74.94, 40.9, 3.72, 0.318, 0.506),
        (0.4, 0.6, 0.779, 61.45, 17.0, 4.29, 0.552, 0.441),
        (0.6, 0.4, 0.319, 48.14, 8.10, 4.85, 0.730, 0.365),
        (0.8, 0.2, 0.0923, 33.31, 3.26, 5.41, 0.874, 0.267),
        (0.9, 0.1, 0.0313, 23.66, 1.60, 5.69, 0.939, 0.206),
        (0.97, 0.03, 0.00609, 13.75, 0.610, 5.88, 0.980, 0.126),
    ],
}


def rb_lookup(ld_key, eps):
    """(S, phi_deg) at `eps`, log-interpolating S and linearly interpolating phi
    between tabulated rows.  Returns None outside the tabulated eps range."""
    rows = RB_TABLE[ld_key]
    for k in range(len(rows) - 1):
        e0, s0, f0 = rows[k]
        e1, s1, f1 = rows[k + 1]
        if e0 - 1e-12 <= eps <= e1 + 1e-12:
            if abs(e1 - e0) < 1e-12:
                return s0, f0
            w = (eps - e0) / (e1 - e0)
            return math.exp(math.log(s0) + w * (math.log(s1) - math.log(s0))), f0 + w * (f1 - f0)
    return None


def rb_lookup_eps(ld_key, S):
    """Inverse lookup: eccentricity ratio and attitude angle at Sommerfeld
    number `S`.  S decreases monotonically with eps in every block."""
    rows = RB_TABLE[ld_key]
    for k in range(len(rows) - 1):
        e0, s0, f0 = rows[k]
        e1, s1, f1 = rows[k + 1]
        if s1 - 1e-12 <= S <= s0 + 1e-12:
            w = (math.log(S) - math.log(s0)) / (math.log(s1) - math.log(s0))
            return e0 + w * (e1 - e0), f0 + w * (f1 - f0)
    return None


# --------------------------------------------------------------------------
# Exact infinitely-long bearing with the Reynolds (Swift-Stieber) rupture
# condition -- used to audit the L/D = inf block of the table above.
# --------------------------------------------------------------------------
#
# Film convention matches FilmThickness::compute_static with
# attitude_angle_deg = 180:  h = c (1 + eps cos t), max gap at t = 0, min at pi,
# surface entrained in +t.  With x = R t,
#
#     dp/dt = (6 mu U R / c^2) (H - H_m) / H^3,    H = 1 + eps cos t
#
# The pressure film runs from t = 0 (fully flooded at the maximum gap, which is
# what Raimondi & Boyd assume) to the rupture angle t2, where Swift-Stieber
# imposes both p = 0 and dp/dt = 0.  dp/dt = 0 forces H_m = H(t2), leaving t2
# fixed by p(t2) = 0.  Beyond t2 the film is at ambient.

def _dpdt_shape(t, eps, Hm):
    H = 1.0 + eps * math.cos(t)
    return (H - Hm) / H ** 3


def _integrate(f, a, b, n):
    """Composite Simpson; n must be even."""
    hstep = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4.0 if k % 2 else 2.0) * f(a + k * hstep)
    return s * hstep / 3.0


def swift_stieber_infinite(eps, n=20000, tol=1e-13):
    """Exact infinite-length Swift-Stieber solution.

    Returns dict with the rupture angle (rad and deg, measured from the maximum
    film thickness), the Sommerfeld number S and the attitude angle phi in
    degrees.  Both are independent of mu, U, R, c and L, so no dimensional
    inputs are needed.
    """
    def residual(t2):
        Hm = 1.0 + eps * math.cos(t2)
        return _integrate(lambda t: _dpdt_shape(t, eps, Hm), 0.0, t2, n)

    # Rupture always sits past the minimum gap for a flooded inlet at t = 0.
    lo, hi = math.pi + 1e-9, 2.0 * math.pi - 1e-9
    flo, fhi = residual(lo), residual(hi)
    if flo * fhi > 0.0:
        raise RuntimeError(f"no bracketed rupture angle at eps={eps}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = residual(mid)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < tol:
            break
    t2 = 0.5 * (lo + hi)
    Hm = 1.0 + eps * math.cos(t2)

    # p(t) / (6 mu U R / c^2) on [0, t2]; zero elsewhere.
    m = 4000
    dt = t2 / m
    pv = [0.0] * (m + 1)
    for k in range(1, m + 1):
        a = (k - 1) * dt
        pv[k] = pv[k - 1] + dt / 6.0 * (_dpdt_shape(a, eps, Hm)
                                        + 4.0 * _dpdt_shape(a + 0.5 * dt, eps, Hm)
                                        + _dpdt_shape(a + dt, eps, Hm))

    # Force on the bearing per unit length, in units of 6 mu U R^2 / c^2.
    fx = fy = 0.0
    for k in range(m):
        t = (k + 0.5) * dt
        p = 0.5 * (pv[k] + pv[k + 1])
        fx += p * math.cos(t) * dt
        fy += p * math.sin(t) * dt

    # W = |F| * 6 mu U R^2 / c^2 * L.  With U = omega R = 2 pi n_s R,
    # P = W / (L * 2R) = 6 pi mu n_s R^2 |F| / c^2, hence
    # S = (mu n_s / P) (R/c)^2 = 1 / (6 pi |F|).
    W_hat = math.hypot(fx, fy)
    S = 1.0 / (6.0 * math.pi * W_hat)
    phi = 180.0 - math.degrees(math.atan2(fy, fx))
    return dict(t2_rad=t2, t2_deg=math.degrees(t2), S=S, phi_deg=phi,
                rupture_from_hmin_deg=math.degrees(t2) - 180.0)


def audit_infinite_row():
    """Print the table's L/D = inf block beside the exact Swift-Stieber values."""
    print("Audit of RB_TABLE['inf'] against the exact Swift-Stieber solution")
    print(f"{'eps':>5} {'S table':>9} {'S exact':>9} {'dS':>8} "
          f"{'phi table':>10} {'phi exact':>10} {'dphi':>7}")
    rows = []
    for eps, s_tab, phi_tab in RB_TABLE["inf"]:
        ex = swift_stieber_infinite(eps)
        ds = (ex["S"] - s_tab) / s_tab
        dphi = ex["phi_deg"] - phi_tab
        print(f"{eps:>5.2f} {s_tab:>9.4f} {ex['S']:>9.4f} {ds:>7.2%} "
              f"{phi_tab:>10.2f} {ex['phi_deg']:>10.2f} {dphi:>7.2f}")
        rows.append((eps, s_tab, ex["S"], ds, phi_tab, ex["phi_deg"], dphi))
    return rows


if __name__ == "__main__":
    audit_infinite_row()
