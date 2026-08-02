"""Target D -- Brito et al. (2006) two-axial-groove THD benchmark.

A second, independent THD dataset for the energy equation, whose verdict
otherwise rests entirely on Ferron.  See `brito.py` for the rig, the oil and the
per-number provenance.

Geometry fixed point.  The two grooves sit at +-90 deg to the load line and
Brito's abscissa is measured from the *upstream* groove, so the layout is

    upstream groove at 0 deg, load at 90 deg, downstream groove at 180 deg

in Brito's frame.  In pancake the load direction is an output of the solve, so
the groove positions are found by iterating

    upstream_groove = angle(pressure_force) - 90 deg

which is the two-groove analogue of the `groove = load + 180` fixed point used
for Ferron's single groove.  Note the different offset: reusing Ferron's
convention here would rotate the whole comparison by 90 deg.

Known approximations, both stated rather than hidden:
  * Brito's grooves span 70 mm of the 80 mm bush width; `inlet_groove` is
    full-width, so the modelled grooves are 14% too long axially.
  * The measured pressure goes sub-ambient in the unloaded lobe (-0.05 MPa).
    With `cavitation_threshold = SCALAR_PCAV` at ambient the solver cannot go
    below ambient, so that part of the profile cannot be reproduced.

Usage:  $PY run_brito.py <path-to-pancake.exe> [workdir]
"""
import json, math, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import brito as B
from run_raimondi_boyd import read_field, read_vector3, OUTER_RE, COUPLING_RE

P_AMB = 1.0e5
# Same reasoning as the Ferron case: Elrod carries pressure inside theta, so the
# free boundary needs p_max/beta of order 1e-2.  Brito's peak is ~2.9 MPa, so
# beta = 3e8 gives p_max/beta ~ 1e-2 with under 1% density variation.
BULK = 3.0e8
N_THETA = 1440
N_Z = 60
MAX_OUTER = 4000
MAX_COUPLING = 80

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
temperature_reference = {t_ref}
journal_wall_temperature = {t_wall}
bearing_wall_temperature = {t_wall}
rho_cp = {rho_cp}
thermal_conductivity = {k_th}
journal_heat_transfer = {h_wall}
bearing_heat_transfer = {h_wall}
fluid_property_model = {fpm}
dissolved_gas_species = AIR
solubility_model = HENRY
density_model = CONSTANT
liquid_viscosity_model = {vmodel}
dissolved_gas_initial = 0
dissolved_gas_max = 1
dissolved_gas_henry_coeff = 1e-12
dissolved_gas_henry_temp_coeff = 0
dissolved_gas_bunsen_coeff = 0
dissolved_gas_liquid_density = 0
gas_mass_transfer_rate = 0
dissolved_gas_diffusivity = 0
solution_density_gas_coeff = 0
solution_viscosity_gas_coeff = 1e-12
viscosity_temperature_coeff = {e_mu}
viscosity_pressure_coeff = 0
gas_alpha_max = 1
mixture_viscosity_model = MCADAMS
mu_gas = 0
gas_pressure_floor = 1.0
property_reference_pressure = {p_amb}
property_reference_temperature = {t_ref}
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
coupling_max_iters = {max_coupling}
coupling_tolerance = 1e-8
coupling_relaxation = {relax}
steady_gas_model = FROZEN
allow_frozen_gas_steady = true
inlet_groove = {g_up} {g_width} {p_sup} {t_sup}
inlet_groove = {g_down} {g_width} {p_sup} {t_sup}
output_dir = out
filename_prefix = solution
output_verbosity = DETAILED
output_write_3d = false
output_write_flat = true
output_fields = pressure, film_content, h, temperature, mu, pressure_force_x, pressure_force_y
"""


def run_case(exe, work, tag, *, eps, rpm, upstream_groove_deg, p_supply_gauge,
             t_supply_c=40.0, thermal="ENERGY_EQUATION", h_wall=0.0, relax=0.7,
             nth=N_THETA, nz=N_Z):
    d = os.path.join(work, tag)
    os.makedirs(d, exist_ok=True)
    isothermal = (thermal == "ISOTHERMAL")
    t_sup = t_supply_c + 273.15
    cfg = CFG.format(
        R=B.R, c=B.C, e=eps * B.C, L=B.L, nth=nth, nz=nz,
        omega=2.0 * math.pi * rpm / 60.0, mu=B.MU_40, rho=B.RHO,
        p_amb=P_AMB, bulk=BULK, tmodel=thermal,
        t_sup=t_sup,
        # mu_40 is quoted at 40 C, so the viscosity reference stays at 40 C even
        # when the supply temperature is moved.  Tying it to t_sup would silently
        # re-define the oil at each new temperature.
        t_ref=313.15,
        t_wall=t_sup, rho_cp=B.RHO_CP, k_th=B.K_OIL, h_wall=h_wall,
        fpm="CONSTANT" if isothermal else "SINGLE_PHASE",
        vmodel="CONSTANT" if isothermal else "EMPIRICAL",
        e_mu=0.0 if isothermal else B.E_MU,
        max_outer=MAX_OUTER, max_coupling=MAX_COUPLING,
        relax=1.0 if isothermal else relax,
        g_up=upstream_groove_deg % 360.0,
        g_down=(upstream_groove_deg + 180.0) % 360.0,
        g_width=B.GROOVE_ANGLE_DEG,
        p_sup=P_AMB + p_supply_gauge)
    open(os.path.join(d, "config.txt"), "w").write(cfg)
    r = subprocess.run([exe, "--config", "config.txt"], cwd=d,
                       capture_output=True, text=True, timeout=7200)

    outers = OUTER_RE.findall(r.stderr)
    coup = COUPLING_RE.findall(r.stdout) or COUPLING_RE.findall(r.stderr)
    if not outers:
        raise RuntimeError(f"{tag}: no outer log\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    n_it, mc, ns, _, _ = outers[-1]
    n_it, mc, ns = int(n_it), float(mc), int(ns)
    R_p = float(coup[-1][1]) if coup else float("nan")
    R_T = float(coup[-1][2]) if coup else float("nan")
    # Record the iteration counts, not just the residuals: without them a run that
    # exhausted the coupling budget is indistinguishable in the JSON from one that
    # converged, and R_p alone does not say which.
    elrod_iters = n_it + 1
    coupling_iters = int(coup[-1][0]) if coup else -1
    coupling_capped = coupling_iters >= MAX_COUPLING
    converged = (n_it < MAX_OUTER - 1) and mc < 1e-10 and ns == 0 and R_p < 1e-3

    f = os.path.join(d, "out", "flat", "processor0", "solution_0_0.vts")
    if not os.path.exists(f):
        raise RuntimeError(f"{tag}: no output\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    fp = read_vector3(f, "Fp")
    return dict(p=read_field(f, nth, nz, "p"),
                theta=read_field(f, nth, nz, "film_content"),
                h=read_field(f, nth, nz, "h"),
                T=read_field(f, nth, nz, "T"),
                mu=read_field(f, nth, nz, "mu"),
                fx=fp[0], fy=fp[1], nth=nth, nz=nz,
                mc=mc, ns=ns, R_p=R_p, R_T=R_T, converged=converged,
                elrod_iters=elrod_iters, coupling_iters=coupling_iters,
                coupling_capped=coupling_capped,
                upstream_groove_deg=upstream_groove_deg % 360.0, eps=eps, rpm=rpm,
                p_supply_gauge=p_supply_gauge, t_supply_c=t_supply_c,
                thermal=thermal, h_wall=h_wall, dir=d)


def load_angle_deg(res):
    return math.degrees(math.atan2(res["fy"], res["fx"])) % 360.0


def solve_groove(exe, work, tag, *, eps, rpm, p_supply_gauge, tol=1.0, max_it=6,
                 groove0=None, **kw):
    """Fixed point on 'the load sits 90 deg downstream of the upstream groove'.

    The iteration contracts by only about 0.4 per step, so seeding it with the
    previous case's answer (`groove0`) cuts it from five solves to two -- the
    groove angle tracks the attitude angle, which varies slowly across the
    operating grid.
    """
    groove, hist = (0.0 if groove0 is None else groove0), []
    res = None
    for k in range(max_it):
        res = run_case(exe, work, f"{tag}_g{k}", eps=eps, rpm=rpm,
                       upstream_groove_deg=groove, p_supply_gauge=p_supply_gauge, **kw)
        target = (load_angle_deg(res) - B.LOAD_ANGLE_FROM_UPSTREAM_GROOVE_DEG) % 360.0
        delta = (target - groove + 180.0) % 360.0 - 180.0
        hist.append((groove, target, delta))
        groove = (groove + delta) % 360.0
        if abs(delta) < tol:
            break
    res["groove_history"] = hist
    return res


def profile(res):
    """Mid-plane pressure [MPa gauge] and temperature [C] against Brito's angle
    from the centre of the upstream groove."""
    nth, mid = res["nth"], res["nz"] // 2
    rows = [(B.wrap_from_upstream_groove((i + 0.5) * 360.0 / nth,
                                         res["upstream_groove_deg"]),
             (res["p"][mid][i] - P_AMB) / 1e6,
             res["T"][mid][i] - 273.15) for i in range(nth)]
    rows.sort()
    return rows


def metrics(res):
    rows = profile(res)
    ps = [r[1] for r in rows]
    Ts = [r[2] for r in rows]
    ip = max(range(len(ps)), key=lambda k: ps[k])
    iT = max(range(len(Ts)), key=lambda k: Ts[k])
    W = math.hypot(res["fx"], res["fy"])
    return dict(eps=res["eps"], rpm=res["rpm"],
                p_supply_gauge=res["p_supply_gauge"], t_supply_c=res["t_supply_c"],
                thermal=res["thermal"], h_wall=res["h_wall"],
                upstream_groove_deg=res["upstream_groove_deg"],
                W_N=W,
                # Attitude angle = angle from the load line to the line of
                # centres, measured with the rotation.  attitude_angle_deg = 180
                # puts the minimum film at theta = 180, so phi = 180 - load angle.
                # (Subtracting the 90 deg groove offset instead returns the
                # groove angle, which the fixed point drives to load - 90 and is
                # not the attitude angle at all.)
                attitude_deg=(180.0 - load_angle_deg(res)) % 360.0,
                p_max_MPa=ps[ip], p_max_angle_deg=rows[ip][0],
                p_min_MPa=min(ps),
                T_max_C=Ts[iT], T_max_angle_deg=rows[iT][0],
                T_rise_K=Ts[iT] - res["t_supply_c"],
                theta_min=min(min(r) for r in res["theta"]),
                p_over_beta=ps[ip] * 1e6 / BULK,
                mc=res["mc"], ns=res["ns"], R_p=res["R_p"],
                elrod_iters=res["elrod_iters"],
                coupling_iters=res["coupling_iters"],
                coupling_capped=res["coupling_capped"],
                converged=res["converged"])


def solve_eps_for_load(exe, work, tag, *, W_target, rpm, p_supply_gauge,
                       tol=0.02, max_it=6, nth=720, nz=40,
                       confirm_nth=N_THETA, confirm_nz=N_Z, **kw):
    """Eccentricity that carries `W_target`, by secant on log W against the
    logit of eps, then one confirmation run at full resolution."""
    def evaluate(k, eps):
        res = solve_groove(exe, work, f"{tag}_e{k}", eps=eps, rpm=rpm,
                           p_supply_gauge=p_supply_gauge, nth=nth, nz=nz, **kw)
        W = math.hypot(res["fx"], res["fy"])
        print(f"      eps={eps:.4f} -> W={W:.0f} N")
        sys.stdout.flush()
        return W

    t_of = lambda e: math.log(e / (1.0 - e))
    eps_of = lambda t: 1.0 / (1.0 + math.exp(-t))
    t0, t1 = t_of(0.35), t_of(0.60)
    y0 = math.log(evaluate(0, eps_of(t0)) / W_target)
    y1 = math.log(evaluate(1, eps_of(t1)) / W_target)
    for k in range(2, max_it):
        if abs(y1) < math.log(1.0 + tol) or abs(y1 - y0) < 1e-12:
            break
        t2 = max(t_of(0.05), min(t_of(0.95), t1 - y1 * (t1 - t0) / (y1 - y0)))
        t0, y0, t1 = t1, y1, t2
        y1 = math.log(evaluate(k, eps_of(t1)) / W_target)
    eps = eps_of(t1)
    res = solve_groove(exe, work, f"{tag}_final", eps=eps, rpm=rpm,
                       p_supply_gauge=p_supply_gauge,
                       nth=confirm_nth, nz=confirm_nz, **kw)
    return eps, res


def main(exe, work):
    os.makedirs(work, exist_ok=True)
    print(f"pancake: {exe}")
    print(f"Brito rig: D={B.D*1e3:g} mm, L={B.L*1e3:g} mm (L/D={B.L_OVER_D:g}), "
          f"c={B.C*1e6:.1f} um radial, c/R={B.C_OVER_R:.3e}")
    print(f"two grooves {B.GROOVE_ANGLE_DEG:.2f} deg wide at 0 and 180 deg from the "
          f"upstream groove; load at {B.LOAD_ANGLE_FROM_UPSTREAM_GROOVE_DEG:g} deg")
    print(f"oil: mu(40C)={B.MU_40:g}, mu(70C)={B.MU_70:g} -> E_mu={B.E_MU:.1f} K "
          f"(exact two-point fit); rho_cp={B.RHO_CP:.4g}, k={B.K_OIL:g}")
    print(f"mesh {N_THETA}x{N_Z}, bulk_modulus {BULK:.0e} Pa\n")

    W_APPLIED = 10000.0
    out = []
    # Each case is run at the MEASURED eccentricity (Fig. 6) rather than
    # root-finding to the applied load.  Brito measures eccentricity, attitude
    # angle and peak pressure at the same operating point, so prescribing the
    # measured eccentricity turns all three into independent predictions from a
    # single solve -- and costs one solve per case instead of six.
    print("=== W = 10 kN, Tf = 40 C, run at the measured eccentricity ===")
    print("    eps and attitude from Fig 6; p_max from Fig 5(c); "
          "load is then a prediction against the applied 10 kN")
    print(f"{'N rpm':>6} {'Pf kPa':>7} {'eps':>6} {'W sim':>8} {'dW':>8} "
          f"{'p_max sim':>10} {'p_max meas':>11} {'dp':>8} "
          f"{'phi sim':>8} {'phi meas':>9} {'ang':>7} {'T_max':>7} {'conv':>5}")
    groove_seed = 32.0        # from the smoke test; updated as the sweep runs
    for rpm in B.SPEEDS_RPM:
        meas_p = dict(B.PMAX_VS_SUPPLY[rpm])
        meas_e = {p: (e, a) for p, e, a in B.ECC_ATTITUDE[rpm]}
        for pf in (70e3, 210e3):
            eps, phi_meas = meas_e[pf]
            tag = f"brito_N{int(rpm)}_Pf{int(pf/1e3)}"
            res = solve_groove(exe, work, tag, eps=eps, rpm=rpm,
                               p_supply_gauge=pf, groove0=groove_seed)
            m = metrics(res)
            m["profile"] = profile(res)
            m["p_max_meas_MPa"] = meas_p[pf]
            m["eps_meas"], m["attitude_meas_deg"] = eps, phi_meas
            m["W_applied_N"] = W_APPLIED
            out.append(m)
            groove_seed = m["upstream_groove_deg"]
            print(f"{rpm:>6.0f} {pf/1e3:>7.0f} {eps:>6.3f} {m['W_N']:>8.0f} "
                  f"{(m['W_N']-W_APPLIED)/W_APPLIED:>7.1%} "
                  f"{m['p_max_MPa']:>10.3f} {meas_p[pf]:>11.2f} "
                  f"{(m['p_max_MPa']-meas_p[pf])/meas_p[pf]:>7.1%} "
                  f"{m['attitude_deg']:>8.1f} {phi_meas:>9.1f} "
                  f"{m['p_max_angle_deg']:>7.1f} {m['T_max_C']:>7.2f} "
                  f"{str(m['converged']):>5}")
            sys.stdout.flush()

    dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), "brito_results.json")
    json.dump(out, open(dst, "w"), indent=1)
    print(f"\nwrote {dst}")
    return out


if __name__ == "__main__":
    exe = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "_work")
    main(exe, work)
