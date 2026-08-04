"""Exact incompressible-JFO solution, infinitely long journal bearing with one
fixed-pressure groove. Semi-analytic: steady 1D mass conservation makes the land
flux constant, q = U h_r / 2, so

    dp/dtheta = 6 mu U R (h - h_r) / h^3,   h(theta) = c + e sin(theta)

with p = p_s at both groove edges, p = p_cav and dp/dtheta = 0 at rupture theta_r,
and flux continuity (the JFO jump) automatic at reformation theta_f.
theta_r solves p(theta_r) = p_cav integrating from the groove trailing edge;
theta_f solves p(groove leading edge) = p_s integrating from theta_f.
Quadrature: composite Simpson on 2^20+1 points; bisection to 1e-12 rad.
"""
import math
import numpy as np

# Parameters matched to cases/00_grid_dt_independence grid template (Ferron-type)
R = 0.05
c = 1.52e-4
e = 9.12e-5          # eps = 0.6
mu = 0.0277
omega = 418.879020479
U = omega * R
p_s = 1.713e5
p_cav = 1.0e5
G1 = math.radians(81.0)    # groove leading edge (upstream side, receives film)
G2 = math.radians(99.0)    # groove trailing edge (film leaves groove here)
K = 6.0 * mu * U * R


def h_of(th):
    return c + e * np.sin(th)


def pressure_rise(th_a, th_b, h_r, n=1 << 20):
    """integral of K (h - h_r)/h^3 dtheta from th_a to th_b (Simpson)."""
    th = np.linspace(th_a, th_b, n + 1)
    f = K * (h_of(th) - h_r) / h_of(th) ** 3
    w = np.ones(n + 1); w[1:-1:2] = 4; w[2:-1:2] = 2
    return (th_b - th_a) / n / 3.0 * np.dot(w, f)


def solve_theta_r():
    """p(theta_r) = p_cav integrating from G2 with h_r = h(theta_r)."""
    def resid(th_r):
        return p_s + pressure_rise(G2, th_r, h_of(th_r)) - p_cav
    lo, hi = math.radians(270.001), math.radians(355.0)
    assert resid(lo) > 0 > resid(hi), (resid(lo), resid(hi))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def solve_theta_f(th_r):
    """p(G1 + 2 pi) = p_s integrating from theta_f with the same h_r."""
    h_r = h_of(th_r)
    end = G1 + 2.0 * math.pi
    def resid(th_f):
        return p_cav + pressure_rise(th_f, end, h_r) - p_s
    lo, hi = th_r + 1e-9, end - 1e-9
    if resid(lo) < 0:      # cannot reach p_s even from rupture: reforms at groove
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if resid(mid) > 0: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def pressure_profile(th_r, th_f, n=4096):
    """Return theta (deg) and absolute p over the full land for comparison."""
    h_r = h_of(th_r)
    end = G1 + 2.0 * math.pi
    thetas, ps = [], []
    # region A: G2 -> th_r
    ta = np.linspace(G2, th_r, n)
    pa = [p_s]
    for i in range(1, len(ta)):
        pa.append(pa[-1] + pressure_rise(ta[i-1], ta[i], h_r, n=256))
    thetas += list(ta); ps += pa
    # region B: cavitated
    tb = np.linspace(th_r, th_f, n // 4)
    thetas += list(tb); ps += [p_cav] * len(tb)
    # region C: th_f -> G1+2pi
    tc = np.linspace(th_f, end, n)
    pc = [p_cav]
    for i in range(1, len(tc)):
        pc.append(pc[-1] + pressure_rise(tc[i-1], tc[i], h_r, n=256))
    thetas += list(tc); ps += pc
    return np.array(thetas), np.array(ps)  # unwrapped radians


if __name__ == "__main__":
    th_r = solve_theta_r()
    th_f = solve_theta_f(th_r)
    print(f"exact rupture     theta_r = {math.degrees(th_r):.6f} deg")
    print(f"exact reformation theta_f = {math.degrees(th_f) if th_f else 'AT GROOVE'} deg"
          if th_f is None or True else "")
    if th_f is not None:
        print(f"                        = {math.degrees(th_f) % 360:.6f} deg (wrapped)")
    th_u, p = pressure_profile(th_r, th_f if th_f else G1 + 2*math.pi)
    th_deg = np.degrees(np.mod(th_u, 2 * math.pi))
    np.savetxt("exact_1d_profile.csv",
               np.column_stack([th_deg, p]), delimiter=",",
               header="theta_deg,p_abs_Pa", comments="")
    # per-unit-length load from gauge pressure, on unwrapped angle
    fx = -np.trapezoid((p - p_cav) * np.cos(th_u), th_u) * R
    fy = -np.trapezoid((p - p_cav) * np.sin(th_u), th_u) * R
    print(f"exact W' = {math.hypot(fx, fy):.2f} N/m, "
          f"p_max = {p.max()/1e6:.6f} MPa at {th_deg[p.argmax()]:.3f} deg")
