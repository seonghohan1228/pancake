"""Analytical journal-bearing benchmarks for the pancake Reynolds solver.

Two exact limits of the isothermal, rigid, prescribed-eccentricity Reynolds
equation bracket any finite bearing, and between them they exercise both terms of
the operator:

  * Infinitely-long (Sommerfeld) limit, L/D -> inf: axial flow vanishes, so the
    circumferential term alone sets p(theta). Full-Sommerfeld (no cavitation)
    has a closed form, which isolates the theta discretization.

  * Short-bearing (Ocvirk) limit, L/D -> 0: the circumferential pressure flow is
    negligible against the axial one, so the axial term alone sets p(theta, z).

Film convention matches FilmThickness::compute_static with attitude_angle = 180:
    h(theta) = c (1 + eps cos theta),  max gap at theta = 0, min gap at theta = pi.
Surface speed U = omega R_journal, entrained in +theta.
"""
import math


def film(theta, c, eps):
    return c * (1.0 + eps * math.cos(theta))


def sommerfeld_long_pressure(theta, *, mu, U, R, c, eps, p_ref):
    """Full-Sommerfeld infinite-length pressure [Pa abs].

        p = p_ref + (6 mu U R / c^2) * eps sin(t) (2 + eps cos t)
                    / [ (2 + eps^2) (1 + eps cos t)^2 ]

    Antisymmetric about theta = 0/pi, so it integrates to zero net circumferential
    force and dips below p_ref on the divergent side -- the classic unphysical
    tension that cavitation models remove.
    """
    s, co = math.sin(theta), math.cos(theta)
    num = eps * s * (2.0 + eps * co)
    den = (2.0 + eps * eps) * (1.0 + eps * co) ** 2
    return p_ref + (6.0 * mu * U * R / (c * c)) * num / den


def ocvirk_short_pressure(theta, z, *, mu, U, R, c, eps, L, p_ref, gumbel=True):
    """Ocvirk short-bearing pressure [Pa abs].

    Integrating d/dz (h^3 dp/dz) = 6 mu U dh/dx with x = R theta and p = p_ref at
    z = 0, L gives, with zeta = z - L/2,

        p = p_ref + (3 mu U eps sin t) / (R c^2 (1 + eps cos t)^3) * (L^2/4 - zeta^2).

    Positive over 0 < t < pi (convergent side). Gumbel clips the negative lobe.
    """
    s, co = math.sin(theta), math.cos(theta)
    zeta = z - 0.5 * L
    amp = (3.0 * mu * U * eps * s) / (R * c * c * (1.0 + eps * co) ** 3)
    p = p_ref + amp * (0.25 * L * L - zeta * zeta)
    if gumbel and p < p_ref:
        return p_ref
    return p


def ocvirk_load(*, mu, U, R, c, eps, L, gumbel=True):
    """Short-bearing load capacity [N] by numerical quadrature of the pressure
    field. Reported alongside the Ocvirk closed form

        W = (mu U L^3 eps) / (4 c^2 R) * sqrt(pi^2 (1 - eps^2) + 16 eps^2)
            / (1 - eps^2)^2

    which is the standard Ocvirk result for the Gumbel (half-Sommerfeld) film.
    """
    closed = (mu * U * L ** 3 * eps) / (4.0 * c * c * R) * \
        math.sqrt(math.pi ** 2 * (1.0 - eps ** 2) + 16.0 * eps ** 2) / (1.0 - eps ** 2) ** 2

    n_t, n_z = 4000, 200
    fx = fy = 0.0
    dt = 2.0 * math.pi / n_t
    dz = L / n_z
    for i in range(n_t):
        t = (i + 0.5) * dt
        for j in range(n_z):
            z = (j + 0.5) * dz
            pg = ocvirk_short_pressure(t, z, mu=mu, U=U, R=R, c=c, eps=eps,
                                       L=L, p_ref=0.0, gumbel=gumbel)
            fx -= pg * math.cos(t) * R * dt * dz
            fy -= pg * math.sin(t) * R * dt * dz
    return dict(closed_form=closed, quadrature=math.hypot(fx, fy),
                fx=fx, fy=fy, attitude_deg=math.degrees(math.atan2(fy, fx)))


def sommerfeld_number(*, mu, N_rev_s, W, L, D, R, c):
    """S = (mu N / P) (R/c)^2 with P = W / (L D). The dimensionless group used by
    the Raimondi & Boyd charts, so it is comparable across clearance ratios."""
    P = W / (L * D)
    return (mu * N_rev_s / P) * (R / c) ** 2
