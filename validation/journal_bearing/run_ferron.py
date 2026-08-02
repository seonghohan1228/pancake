"""Target B -- Ferron, Frene & Boncompain (1983) thermohydrodynamic benchmark.

Runs the Ferron rig at its reference operating point with
`temperature_model = ENERGY_EQUATION` and Ferron's own oil, and compares the
predicted mid-plane circumferential pressure and bearing-wall temperature
against the measurements (see `ferron.py` for the data and its provenance).
An isothermal run at the inlet viscosity is included as the upper bracket: with
no thermal viscosity drop it must overpredict the peak pressure, and the
measurement should sit between the two.

Geometry convention.  `attitude_angle_deg = 180` puts the minimum film at
theta = 180 deg and the shaft entrains in +theta.  Ferron plots against X, which
is zero opposite the supply groove and +-pi at it.  The rig has a single axial
groove feeding the unloaded region, so the groove is taken to sit diametrically
opposite the load line; since the load direction is an output of the solve, the
groove angle is found by fixed-point iteration (groove = load angle + 180 deg).

That assumption is the one convention-dependent input in the comparison, so two
quantities are reported that do not depend on it at all -- the peak pressure
magnitude, and the arc from the pressure peak to the rupture point -- alongside
the peak position relative to the groove, which does.

Usage:  $PY run_ferron.py <path-to-pancake.exe> [workdir]
"""
import json, math, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpi_launch import solver_command, describe as describe_ranks
import ferron as F
from run_raimondi_boyd import read_field, read_vector3, OUTER_RE, COUPLING_RE

P_AMB = 1.0e5
# Elrod resolves pressure through theta - 1 = p/beta; the Raimondi & Boyd
# free-boundary audit (see run_raimondi_boyd.BULK) showed the rupture location
# needs p_max/beta of order 1e-2 on a 2880-cell mesh.  Ferron's peak is ~1.3 MPa,
# so beta = 1e8 puts p_max/beta at 1.3e-2.  The price is a ~1.3% density
# variation in the loaded region, which is reported per case.
BULK = 1.0e8
# 1440 x 60 rather than the 2880 x 80 used for Raimondi & Boyd.  The acceptance
# tolerances here are 15% on peak pressure and 10 deg on peak and rupture
# location, and the measured rupture is itself only bracketed to 22 deg by the
# spacing of Ferron's pressure taps, so a finer mesh would resolve nothing the
# measurement can test.  `resolution_check()` reports 720 vs 1440 so the claim is
# backed rather than asserted.
N_THETA = 1440
N_Z = 60
MAX_OUTER = 4000

E_MU, E_MU_ERR = F.andrade_from_exponential()

CFG = """R = {R}
c = {c}
e = {e}
L = {L}
attitude_angle_deg = 180.0
load_angle_deg = 0.0
tilt_slope_x = 0
tilt_slope_y = 0
n_theta_global = {nth}
n_z_global = {nz}
end_t = 0.0
dt = 1e-5
write_interval = 1e-5
solution_mode = STEADY_STATE
motion_time_method = EULER_IMPLICIT
theta_convection_scheme = UPWIND
thermal_convection_scheme = UPWIND
gas_convection_scheme = UPWIND
diagnostics_interval = 0
omega = {omega}
omega_ramp_time = 0.0
mu = {mu}
rho = {rho}
p_cav = {p_amb}
bulk_modulus = {bulk}
temperature_model = {tmodel}
temperature_initial = {t_sup}
temperature_reference = {t_sup}
journal_wall_temperature = {t_wall}
bearing_wall_temperature = {t_wall}
rho_cp = {rho_cp}
thermal_conductivity = {k_th}
journal_heat_transfer = {h_j}
bearing_heat_transfer = {h_b}
fluid_property_model = {fpm}
dissolved_gas_species = AIR
solubility_model = HENRY
density_model = CONSTANT
liquid_viscosity_model = {vmodel}
dissolved_gas_initial = 0
dissolved_gas_max = 1
# Ferron's oil carries no dissolved gas in this comparison.  The validator still
# insists on a positive Henry coefficient whenever a property model is active, so
# it is set to a placeholder; with dissolved_gas_initial = 0 and
# gas_mass_transfer_rate = 0 nothing is ever released and it has no effect.
dissolved_gas_henry_coeff = 1e-12
dissolved_gas_henry_temp_coeff = 0
dissolved_gas_bunsen_coeff = 0
dissolved_gas_liquid_density = 0
gas_mass_transfer_rate = 0
dissolved_gas_diffusivity = 0
solution_density_gas_coeff = 0
solution_viscosity_gas_coeff = {a_c}
viscosity_temperature_coeff = {e_mu}
viscosity_pressure_coeff = 0
gas_alpha_max = 1
mixture_viscosity_model = MCADAMS
mu_gas = 0
gas_pressure_floor = 1.0
property_reference_pressure = {p_amb}
property_reference_temperature = {t_sup}
solubility_table =
density_table =
viscosity_table =
motion_model = STATIC
min_film_thickness = 1e-09
stop_on_nonpositive_film = true
bc_z_south_type = DIRICHLET
bc_z_south_val = {p_amb}
bc_z_south_thermal = ZERO_GRADIENT
bc_z_north_type = DIRICHLET
bc_z_north_val = {p_amb}
bc_z_north_thermal = ZERO_GRADIENT
cavitation_model = JFO
cavitation_threshold = SCALAR_PCAV
gas_pressure_coupling = NONE
max_outer_iters = {max_outer}
outer_tol = 1e-10
theta_min = 1e-6
log_outer_iters = true
coupling_max_iters = 80
coupling_tolerance = 1e-8
coupling_relaxation = {relax}
steady_gas_model = FROZEN
allow_frozen_gas_steady = true
inlet_groove = {groove} {groove_size} {p_sup} {t_sup}
output_dir = out
filename_prefix = solution
output_verbosity = DETAILED
output_write_3d = false
output_write_flat = true
output_fields = pressure, film_content, h, temperature, mu, velocity_z, pressure_force_x, pressure_force_y
"""


def run_case(exe, work, tag, *, eps, rpm, groove_deg, thermal="ENERGY_EQUATION",
             h_wall=0.0, relax=0.7, nth=N_THETA, nz=N_Z):
    d = os.path.join(work, tag)
    os.makedirs(d, exist_ok=True)
    isothermal = (thermal == "ISOTHERMAL")
    cfg = CFG.format(
        R=F.R, c=F.C, e=eps * F.C, L=F.L, nth=nth, nz=nz,
        omega=2.0 * math.pi * rpm / 60.0, mu=F.MU_INLET, rho=F.RHO,
        p_amb=P_AMB, bulk=BULK, tmodel=thermal,
        t_sup=F.T_SUPPLY, t_wall=F.T_SUPPLY, rho_cp=F.RHO_CP, k_th=F.K_THERMAL,
        h_j=h_wall, h_b=h_wall,
        fpm="CONSTANT" if isothermal else "SINGLE_PHASE",
        vmodel="CONSTANT" if isothermal else "EMPIRICAL",
        e_mu=0.0 if isothermal else E_MU,
        # The EMPIRICAL viscosity law multiplies exp(a_c * c_dissolved).  The
        # config validator rejects a_c = 0, but Ferron's oil carries no dissolved
        # gas in this model (dissolved_gas_initial = 0), so c_dissolved is
        # identically zero and any non-zero a_c leaves the viscosity untouched.
        a_c=1.0e-12,
        max_outer=MAX_OUTER, relax=1.0 if isothermal else relax,
        groove=groove_deg, groove_size=F.GROOVE_ANGLE_DEG,
        p_sup=P_AMB + F.P_SUPPLY_GAUGE)
    open(os.path.join(d, "config.txt"), "w").write(cfg)
    r = subprocess.run(solver_command(exe), cwd=d,
                       capture_output=True, text=True, timeout=7200)

    outers = OUTER_RE.findall(r.stderr)
    coup = COUPLING_RE.findall(r.stdout) or COUPLING_RE.findall(r.stderr)
    if not outers:
        raise RuntimeError(f"{tag}: no outer log\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    n_it, mc, ns, _, _ = outers[-1]
    n_it, mc, ns = int(n_it), float(mc), int(ns)
    R_p = float(coup[-1][1]) if coup else float("nan")
    R_T = float(coup[-1][2]) if coup else float("nan")
    coup_iters = int(coup[-1][0]) if coup else -1
    # R_p is the largest pressure change between coupling iterations, normalised
    # by property_reference_pressure = 1e5 Pa, so the threshold below is "the
    # coupling loop has stopped moving the pressure by more than 100 Pa", i.e.
    # 0.01% of the ~1.2 MPa peak.  The configured coupling_tolerance of 1e-8
    # (1e-3 Pa) is far tighter than the comparison can use and the segregated
    # Picard loop does not reach it on the finer mesh within any sane budget;
    # measured behaviour is R_p = 2.6e-6 (0.26 Pa) at 720x40 and 1.6e-3 (160 Pa)
    # at 1440x60, and the two meshes agree on p_max to 1.3%.
    converged = (n_it < MAX_OUTER - 1) and mc < 1e-10 and ns == 0 and R_p < 1e-3

    f = os.path.join(d, "out", "flat", "processor0", "solution_0_0.vts")
    if not os.path.exists(f):
        raise RuntimeError(f"{tag}: no output\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    fp = read_vector3(f, "Fp")
    out = dict(p=read_field(f, nth, nz, "p"),
               theta=read_field(f, nth, nz, "film_content"),
               h=read_field(f, nth, nz, "h"),
               T=read_field(f, nth, nz, "T"),
               mu=read_field(f, nth, nz, "mu"),
               fx=fp[0], fy=fp[1], nth=nth, nz=nz,
               elrod_iters=n_it + 1, mc=mc, ns=ns, R_p=R_p, R_T=R_T,
               coupling_iters=coup_iters, converged=converged,
               groove_deg=groove_deg, eps=eps, rpm=rpm, thermal=thermal,
               h_wall=h_wall, dir=d)
    return out


def load_angle_deg(res):
    return math.degrees(math.atan2(res["fy"], res["fx"])) % 360.0


def solve_groove(exe, work, tag, *, eps, rpm, thermal, h_wall, tol=1.0, max_it=5,
                 nth=N_THETA, nz=N_Z):
    """Fixed point on 'the groove sits opposite the load line'.  Returns the
    converged case plus the iteration history.

    `tol` is 1 deg: the groove angle only sets the origin of Ferron's X
    coordinate, and the acceptance tolerance on peak and rupture location is
    10 deg, so resolving it further buys nothing and costs a full solve each."""
    groove, hist = 0.0, []
    res = None
    for k in range(max_it):
        res = run_case(exe, work, f"{tag}_g{k}", eps=eps, rpm=rpm,
                       groove_deg=groove, thermal=thermal, h_wall=h_wall,
                       nth=nth, nz=nz)
        target = (load_angle_deg(res) + 180.0) % 360.0
        delta = (target - groove + 180.0) % 360.0 - 180.0
        hist.append((groove, target, delta))
        groove = (groove + delta) % 360.0
        if abs(delta) < tol:
            break
    res["groove_history"] = hist
    return res


def profile(res):
    """Mid-plane pressure [MPa gauge] and temperature [deg C] against Ferron's
    X/pi, sorted."""
    nth, mid = res["nth"], res["nz"] // 2
    rows = []
    for i in range(nth):
        th = (i + 0.5) * 360.0 / nth
        rows.append((F.wrap_xpi(th, res["groove_deg"]),
                     (res["p"][mid][i] - P_AMB) / 1e6,
                     res["T"][mid][i] - 273.15))
    rows.sort()
    return rows


def metrics(res):
    rows = profile(res)
    xs = [r[0] for r in rows]
    ps = [r[1] for r in rows]
    Ts = [r[2] for r in rows]
    ip = max(range(len(ps)), key=lambda k: ps[k])
    iT = max(range(len(Ts)), key=lambda k: Ts[k])

    # Rupture: first angle past the pressure peak where the mid-plane film
    # content drops below full.
    nth, mid = res["nth"], res["nz"] // 2
    theta_peak_deg = (xs[ip] * 180.0 + res["groove_deg"] + 180.0) % 360.0
    i_peak = int(theta_peak_deg / 360.0 * nth)
    rupt_xpi = float("nan")
    for k in range(1, nth):
        i = (i_peak + k) % nth
        if res["theta"][mid][i] < 1.0 - 1e-5:
            rupt_xpi = F.wrap_xpi((i + 0.5) * 360.0 / nth, res["groove_deg"])
            break
    peak_to_rupture_deg = ((rupt_xpi - xs[ip]) * 180.0) % 360.0

    # Side leakage: axial Poiseuille flux through both end faces, evaluated from
    # the solved pressure with a one-sided gradient to the Dirichlet boundary.
    # Only full-film cells carry Poiseuille flow (Elrod switch function).
    Q = 0.0
    nz = res["nz"]
    dth = 2.0 * math.pi / nth
    dz = F.L / nz
    # Both ends leak outwards, so the two contributions add in magnitude; summing
    # them with opposite signs (as a signed axial flux would) cancels them to zero.
    for j in (0, nz - 1):
        for i in range(nth):
            if res["theta"][j][i] < 1.0 - 1e-5:
                continue
            dpdz = (res["p"][j][i] - P_AMB) / (0.5 * dz)
            hj = res["h"][j][i]
            Q += hj ** 3 / (12.0 * res["mu"][j][i]) * dpdz * F.R * dth

    W = math.hypot(res["fx"], res["fy"])
    return dict(eps=res["eps"], rpm=res["rpm"], thermal=res["thermal"],
                h_wall=res["h_wall"], groove_deg=res["groove_deg"],
                W_N=W, attitude_deg=(180.0 - load_angle_deg(res)) % 360.0,
                p_max_MPa=ps[ip], p_max_xpi=xs[ip],
                T_max_C=Ts[iT], T_max_xpi=xs[iT],
                T_rise_K=Ts[iT] - F.T_SUPPLY_C,
                rupture_xpi=rupt_xpi, peak_to_rupture_deg=peak_to_rupture_deg,
                Q_lpm=Q * 6.0e4, p_over_beta=ps[ip] * 1e6 / BULK,
                theta_min=min(min(r) for r in res["theta"]),
                converged=res["converged"], elrod_iters=res["elrod_iters"],
                coupling_iters=res["coupling_iters"],
                mc=res["mc"], ns=res["ns"], R_p=res["R_p"], R_T=res["R_T"])


def report(m, label):
    px, pm = F.measured_peak_pressure_mpa()
    tx, tm = F.measured_peak_temperature_c()
    lo, hi = F.measured_rupture_xpi()
    print(f"\n--- {label} ---")
    print(f"  converged={m['converged']}  elrod_it={m['elrod_iters']} "
          f"coupling_it={m['coupling_iters']} mc={m['mc']:.1e} ns={m['ns']} "
          f"R_p={m['R_p']:.1e} R_T={m['R_T']:.1e}")
    print(f"  groove at theta={m['groove_deg']:.1f} deg, attitude angle "
          f"{m['attitude_deg']:.1f} deg, load {m['W_N']:.0f} N "
          f"(measured 4000 N at eps=0.67)")
    print(f"  p_max      {m['p_max_MPa']:.3f} MPa vs measured {pm:.2f} MPa "
          f"({(m['p_max_MPa']-pm)/pm:+.1%})")
    print(f"  p_max loc  X/pi {m['p_max_xpi']:+.3f} vs measured {px:+.2f} "
          f"({(m['p_max_xpi']-px)*180:+.1f} deg)")
    print(f"  rupture    X/pi {m['rupture_xpi']:+.3f} vs measured bracket "
          f"{lo:+.2f}..{hi:+.2f}")
    print(f"  peak->rupture arc {m['peak_to_rupture_deg']:.1f} deg "
          f"(measurement bracket {(lo-px)*180:.0f}..{(hi-px)*180:.0f} deg)")
    print(f"  T_max      {m['T_max_C']:.2f} C vs measured {tm:.1f} C "
          f"({m['T_max_C']-tm:+.2f} K), rise {m['T_rise_K']:.2f} K "
          f"vs measured {tm-F.T_SUPPLY_C:.1f} K")
    print(f"  T_max loc  X/pi {m['T_max_xpi']:+.3f} vs measured {tx:+.2f}")
    print(f"  Q_side {m['Q_lpm']:.2f} l/min | theta_min {m['theta_min']:.4f} "
          f"| p_max/beta {m['p_over_beta']:.4f}")


def solve_eps_for_load(exe, work, tag, *, W_target, rpm, thermal, h_wall,
                       tol=0.02, max_it=6, nth=720, nz=40,
                       confirm_nth=N_THETA, confirm_nz=N_Z):
    """Find the eccentricity ratio that carries `W_target`, with the groove
    fixed point re-solved at every eccentricity.

    Alakhramsing Figure 9 is specified by load (4 kN), not by eccentricity, so
    matching the load is the like-for-like comparison for the pressure profile.
    The eccentricity that comes out is then itself a prediction, testable against
    the measured eps-load curve of Figure 10a.

    The search runs on a coarse mesh -- it only needs the load, which converges
    far faster than the free-boundary location -- and the accepted eccentricity is
    then re-run once at full resolution, which is the case that gets reported.
    A secant iteration on log(W) against log(eps/(1-eps)) converges in three or
    four solves; bisection would need twelve.
    """
    def evaluate(k, eps):
        res = solve_groove(exe, work, f"{tag}_e{k}", eps=eps, rpm=rpm,
                           thermal=thermal, h_wall=h_wall, nth=nth, nz=nz)
        W = math.hypot(res["fx"], res["fy"])
        print(f"    eps={eps:.4f} -> W={W:.0f} N  ({nth}x{nz})")
        sys.stdout.flush()
        return W

    def t_of(eps):
        return math.log(eps / (1.0 - eps))

    def eps_of(t):
        return 1.0 / (1.0 + math.exp(-t))

    t0, t1 = t_of(0.45), t_of(0.70)
    y0 = math.log(evaluate(0, eps_of(t0)) / W_target)
    y1 = math.log(evaluate(1, eps_of(t1)) / W_target)
    for k in range(2, max_it):
        if abs(y1) < math.log(1.0 + tol):
            break
        if abs(y1 - y0) < 1e-12:
            break
        t2 = t1 - y1 * (t1 - t0) / (y1 - y0)
        t2 = max(t_of(0.05), min(t_of(0.95), t2))
        t0, y0 = t1, y1
        t1 = t2
        y1 = math.log(evaluate(k, eps_of(t1)) / W_target)
    eps = eps_of(t1)
    print(f"    accepted eps={eps:.4f}; re-running at {confirm_nth}x{confirm_nz}")
    sys.stdout.flush()
    res = solve_groove(exe, work, f"{tag}_final", eps=eps, rpm=rpm,
                       thermal=thermal, h_wall=h_wall,
                       nth=confirm_nth, nz=confirm_nz)
    return eps, res, math.hypot(res["fx"], res["fy"])


def resolution_check(exe, work, eps, rpm=F.REF_RPM):
    """Same case at two meshes, to show the reported quantities are resolved at
    the tolerances this benchmark can actually test."""
    print("\n=== mesh resolution check (THD, adiabatic) ===")
    print(f"{'mesh':>10} {'p_max/MPa':>10} {'p_max X/pi':>11} {'rupture X/pi':>13} "
          f"{'T_max/C':>9} {'W/N':>8}")
    rows = []
    for nth, nz in ((720, 40), (1440, 60)):
        res = solve_groove(exe, work, f"ferron_res_{nth}", eps=eps, rpm=rpm,
                           thermal="ENERGY_EQUATION", h_wall=0.0, nth=nth, nz=nz)
        m = metrics(res)
        print(f"{nth:>5d}x{nz:<4d} {m['p_max_MPa']:>10.4f} {m['p_max_xpi']:>11.4f} "
              f"{m['rupture_xpi']:>13.4f} {m['T_max_C']:>9.3f} {m['W_N']:>8.0f}")
        rows.append(dict(nth=nth, nz=nz, p_max_MPa=m["p_max_MPa"],
                         p_max_xpi=m["p_max_xpi"], rupture_xpi=m["rupture_xpi"],
                         T_max_C=m["T_max_C"], W_N=m["W_N"]))
        sys.stdout.flush()
    return rows


def main(exe, work):
    os.makedirs(work, exist_ok=True)
    print(f"pancake: {exe}")
    print(f"Ferron rig: D={F.D*1e3:g} mm, L={F.L*1e3:g} mm (L/D={F.L_OVER_D:g}), "
          f"c={F.C*1e6:g} um (c/R={F.C_OVER_R:.2e}), groove {F.GROOVE_ANGLE_DEG:g} deg")
    print(f"oil: mu({F.T_SUPPLY_C:g}C)={F.MU_INLET:g} Pa.s, rho_cp={F.RHO_CP:.3g}, "
          f"k={F.K_THERMAL:g}, Ferron beta={F.BETA_VISC:g} 1/K")
    print(f"Andrade refit for pancake: E_mu = {E_MU:.1f} K "
          f"(max viscosity error {E_MU_ERR:.2%} over 40-55 C); "
          f"the optical-rig value 4803 K is NOT used")
    print(f"mesh {N_THETA}x{N_Z}, bulk_modulus {BULK:.0e} Pa")

    eps_meas = F.FERRON_ECC_VS_LOAD[2][1]     # 0.67 measured at 4 kN
    out = {}

    # Only the headline THD case is confirmed on the finer mesh; the wall-heat
    # sensitivity and the isothermal bracket are differences against it, and a
    # difference is resolved at a coarser mesh than an absolute value.
    cases = [
        ("thd_adiabatic", dict(thermal="ENERGY_EQUATION", h_wall=0.0),
         "THD, adiabatic walls (the assumption Alakhramsing justify for this rig)",
         (N_THETA, N_Z)),
        ("thd_metalbush", dict(thermal="ENERGY_EQUATION", h_wall=5000.0),
         "THD, metal bush + shaft, h = 5000 W/(m2.K)", (720, 40)),
        ("isothermal", dict(thermal="ISOTHERMAL", h_wall=0.0),
         "ISOTHERMAL at the 40 C inlet viscosity (upper bracket)", (720, 40)),
    ]
    for tag, kw, label, (cn, cz) in cases:
        print(f"\n### {label}")
        print("  root-finding the eccentricity that carries Figure 9's 4 kN load:")
        eps, res, W = solve_eps_for_load(exe, work, f"ferron_{tag}",
                                         W_target=F.REF_LOAD_N, rpm=F.REF_RPM,
                                         confirm_nth=cn, confirm_nz=cz, **kw)
        m = metrics(res)
        m["label"] = label
        m["profile"] = profile(res)
        m["eps_for_4kN"] = eps
        m["eps_measured"] = eps_meas
        out[tag] = m
        report(m, label)
        print(f"  eccentricity needed for 4 kN: {eps:.3f} vs measured {eps_meas:.2f} "
              f"({eps - eps_meas:+.3f})")
        sys.stdout.flush()

    # Load-carrying check across Ferron's measured eccentricity/load pairs.
    print("\n=== load capacity vs Ferron Figure 10a ===")
    print(f"{'W meas':>8} {'eps meas':>9} {'W sim':>9} {'dW':>8} {'conv':>5}")
    lc = []
    for W_meas, eps in F.FERRON_ECC_VS_LOAD[::3]:
        res = solve_groove(exe, work, f"ferron_load_eps{eps:g}", eps=eps,
                           rpm=F.REF_RPM, thermal="ENERGY_EQUATION", h_wall=0.0,
                           nth=720, nz=40)
        m = metrics(res)
        lc.append(dict(W_meas=W_meas, eps=eps, W_sim=m["W_N"],
                       converged=m["converged"], p_max_MPa=m["p_max_MPa"],
                       T_rise_K=m["T_rise_K"]))
        print(f"{W_meas:>8.0f} {eps:>9.2f} {m['W_N']:>9.0f} "
              f"{(m['W_N']-W_meas)/W_meas:>7.1%} {str(m['converged']):>5}")
        sys.stdout.flush()
    out["load_curve"] = lc
    out["resolution"] = resolution_check(exe, work, out["thd_adiabatic"]["eps_for_4kN"])

    here = os.path.dirname(os.path.abspath(__file__))
    dst = os.path.join(here, "ferron_results.json")
    json.dump(out, open(dst, "w"), indent=1)
    print(f"\nwrote {dst}")
    return out


if __name__ == "__main__":
    exe = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "_work")
    main(exe, work)
