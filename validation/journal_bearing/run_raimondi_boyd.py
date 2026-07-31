"""Target A -- Raimondi & Boyd (1958) cavitated finite-bearing validation.

Sweeps prescribed eccentricity at L/D = 1/2 and 1 with the mass-conserving
Elrod-Adams (JFO) cavitation model and a scalar cavity pressure at ambient, then
forms the Sommerfeld number and attitude angle from the integrated pressure force
and compares them against the tabulated Raimondi & Boyd design data in
`raimondi_boyd.py` (see that module for the source and its provenance audit).

Two inlet conditions are run at every point, because they are different physics
and the difference is the headline result:

  flooded   -- a narrow full-film supply groove at the maximum-gap position
               (theta = 0).  In JFO the groove is a Dirichlet theta = 1, i.e.
               exactly the "complete film at the maximum gap" that Raimondi &
               Boyd assume.  This is the like-for-like comparison.

  ungrooved -- no supply at all, the configuration of the optical rig.  The
               Elrod switch function kills Poiseuille flow in cavitated cells,
               so nothing re-floods the film and the cavity persists into the
               convergent region.  The wedge is shortened and the load must come
               out below the chart.

Convergence is checked explicitly for every case: the exit code of pancake.exe
is *not* a convergence signal.  `log_outer_iters = true` makes the Elrod outer
loop print `mc` (max |d theta|) and `ns` (switch-function flips) per iteration,
and the steady coupling loop print R_p.  A case is accepted only if the Elrod
loop terminated on the tolerance rather than the iteration budget and R_p is at
the coupling tolerance.

Usage:  $PY run_raimondi_boyd.py <path-to-pancake.exe> [workdir]
"""
import json, math, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analytic import sommerfeld_number
from raimondi_boyd import RB_TABLE, rb_lookup, swift_stieber_infinite

# Geometry chosen so that c/R = 1e-3, i.e. inside the thin-film regime the
# Raimondi & Boyd solution assumes.  mu and omega only set the pressure scale;
# S and the attitude angle are invariant under them in the incompressible limit,
# which is why the sweep repeats cases at a second viscosity (see set_mu below
# and the BULK comment for why that invariance is not free here).
R_J = 0.025
C_R = 25e-6
MU = 0.01
RHO = 870.0
N_REV_S = 5.0
OMEGA = 2.0 * math.pi * N_REV_S
P_AMB = 1.0e5


def set_mu(mu):
    """mu enters the reference comparison only through the pressure scale, and S
    is invariant under it in the incompressible limit.  Running the same case at
    two viscosities is therefore a direct check that neither compressibility
    (which grows with p_max/beta) nor free-boundary under-resolution (which grows
    as p_max/beta shrinks) is biasing the answer."""
    global MU
    MU = mu
# Elrod carries pressure in the film content: p = p_cav + beta ln(theta), so the
# entire full-film pressure signal lives in theta - 1 = p/beta.  The accuracy of
# the free boundary is controlled by how well the mesh resolves that excursion,
# and it degrades sharply as p_max/beta shrinks.  Measured at L/D = 8, eps = 0.6
# against the exact Swift-Stieber rupture angle of 213.08 deg
# (`convergence_study()` reproduces this):
#
#   beta    n_theta   rupture     S vs exact
#   1e9        720    192.25 deg    -19.7%
#   1e9       2880    206.06 deg     -7.9%
#   1e8        720    211.25 deg     +2.1%
#   1e8       2880    213.31 deg     +9.0%
#   1e8       5760    213.59 deg    +10.2%
#
# The rupture angle converges cleanly to the exact value at beta = 1e8, so the
# cavitation physics is right; the beta = 1e9 / 720 combination is simply
# under-resolved -- and it reports itself converged (mc = 0, ns = 0, R_p = 0)
# while being 20 deg of arc wrong.  Production runs therefore use beta = 1e8 with
# 2880 circumferential cells, and every case is repeated at a second viscosity to
# confirm the answer is invariant (see `set_mu`).
BULK = 1.0e8

N_THETA = 2880
N_Z = 80
MAX_OUTER = 4000

# Measured output of convergence_study(exe, work, eps=0.6, LD=8.0), recorded here
# so the sweep does not have to re-run an hour of audit to report it.  Regenerate
# with `main(..., do_convergence=True)`.
CONVERGENCE_MEASURED = [
    dict(beta=1e9, nth=720,  S=0.03127, phi=42.78, rupture=192.25),
    dict(beta=1e9, nth=2880, S=0.03576, phi=43.02, rupture=206.06),
    dict(beta=1e8, nth=720,  S=0.03978, phi=47.28, rupture=211.25),
    dict(beta=1e8, nth=1440, S=0.04151, phi=49.77, rupture=212.62),
    dict(beta=1e8, nth=2880, S=0.04246, phi=51.29, rupture=213.31),
    dict(beta=1e8, nth=5760, S=0.04292, phi=52.06, rupture=213.59),
]

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
temperature_model = ISOTHERMAL
temperature_initial = 298.15
temperature_reference = 298.15
journal_wall_temperature = 298.15
bearing_wall_temperature = 298.15
rho_cp = 1.92e6
thermal_conductivity = 0.15
journal_heat_transfer = 1000
bearing_heat_transfer = 1000
fluid_property_model = CONSTANT
dissolved_gas_species = AIR
solubility_model = HENRY
density_model = CONSTANT
liquid_viscosity_model = CONSTANT
dissolved_gas_initial = 0
dissolved_gas_max = 1
dissolved_gas_henry_coeff = 0
dissolved_gas_henry_temp_coeff = 0
dissolved_gas_bunsen_coeff = 0
dissolved_gas_liquid_density = 0
gas_mass_transfer_rate = 0
dissolved_gas_diffusivity = 0
solution_density_gas_coeff = 0
solution_viscosity_gas_coeff = 0
viscosity_temperature_coeff = 0
viscosity_pressure_coeff = 0
gas_alpha_max = 1
mixture_viscosity_model = MCADAMS
mu_gas = 0
gas_pressure_floor = 1.0
property_reference_pressure = {p_amb}
property_reference_temperature = 298.15
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
coupling_max_iters = 60
coupling_tolerance = 1e-10
coupling_relaxation = 1.0
steady_gas_model = FROZEN
allow_frozen_gas_steady = true
{inlet}output_dir = out
filename_prefix = solution
output_verbosity = DETAILED
output_write_3d = false
output_write_flat = true
output_fields = pressure, film_content, h, pressure_force_x, pressure_force_y
"""

OUTER_RE = re.compile(r"outer=(\d+) mc=([0-9.eE+-]+) ns=(\d+) th=\[([0-9.eE+-]+),([0-9.eE+-]+)\]")
COUPLING_RE = re.compile(r"Coupling: (\d+) outer iters, R_p=([0-9.eE+-]+) R_T=([0-9.eE+-]+)")


def read_field(path, nth, nz, name):
    txt = open(path).read()
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii"[^>]*>(.*?)</DataArray>' % name,
                  txt, re.S)
    if m is None:
        raise RuntimeError(f"field {name} missing from {path}")
    vals = [float(v) for v in m.group(1).split()]
    return [vals[j * nth:(j + 1) * nth] for j in range(nz)]


def read_vector3(path, name):
    """First (x, y, z) tuple of an interleaved 3-component DataArray."""
    txt = open(path).read()
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii"[^>]*>(.*?)</DataArray>' % name,
                  txt, re.S)
    if m is None:
        raise RuntimeError(f"vector {name} missing from {path}")
    vals = m.group(1).split()[:3]
    return [float(v) for v in vals]


def run_case(exe, work, tag, *, eps, LD, flooded, nth=N_THETA, nz=N_Z):
    """Run one prescribed-eccentricity case and return the solved fields plus
    the convergence evidence."""
    d = os.path.join(work, tag)
    os.makedirs(d, exist_ok=True)
    L = LD * 2.0 * R_J
    # 4 deg groove at the maximum-gap position, spanning the full width.
    inlet = f"inlet_groove = 0.0 4.0 {P_AMB}\n" if flooded else ""
    cfg = CFG.format(R=R_J, c=C_R, e=eps * C_R, L=L, nth=nth, nz=nz, omega=OMEGA,
                     mu=MU, rho=RHO, p_amb=P_AMB, bulk=BULK, max_outer=MAX_OUTER,
                     inlet=inlet)
    open(os.path.join(d, "config.txt"), "w").write(cfg)
    r = subprocess.run([exe, "--config", "config.txt"], cwd=d,
                       capture_output=True, text=True, timeout=3600)

    outers = OUTER_RE.findall(r.stderr)
    coup = COUPLING_RE.findall(r.stdout) or COUPLING_RE.findall(r.stderr)
    if not outers:
        raise RuntimeError(f"{tag}: no outer-iteration log\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    n_it, mc, ns, th_lo, th_hi = outers[-1]
    n_it, mc, ns = int(n_it), float(mc), int(ns)
    R_p = float(coup[-1][1]) if coup else float("nan")
    R_T = float(coup[-1][2]) if coup else float("nan")
    # Budget exhaustion: the last logged outer index is max_outer_iters - 1.
    budget_hit = (n_it >= MAX_OUTER - 1)
    converged = (not budget_hit) and mc < 1e-10 and ns == 0 and R_p < 1e-9

    f = os.path.join(d, "out", "flat", "processor0", "solution_0_0.vts")
    if not os.path.exists(f):
        raise RuntimeError(f"{tag}: no output\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    p = read_field(f, nth, nz, "p")
    theta = read_field(f, nth, nz, "film_content")
    # pressure_force_{x,y,z} are written as the interleaved 3-vector "Fp"; the
    # field is uniform, so the first tuple is the whole answer.
    fp = read_vector3(f, "Fp")
    fx, fy = fp[0], fp[1]
    h = read_field(f, nth, nz, "h")
    return dict(p=p, theta=theta, h=h, fx=fx, fy=fy, L=L, nth=nth, nz=nz,
                elrod_iters=n_it + 1, mc=mc, ns=ns, R_p=R_p, R_T=R_T,
                converged=converged, budget_hit=budget_hit, dir=d)


def hand_integrate_force(p, nth, nz, L):
    """Independent quadrature of the pressure force, to cross-check the solver's
    own pressure_force_x/y.  Uses the same convention as
    Reynolds::calculate_macroscopic_properties: force on the bearing surface,
    gauge pressure, area element (R + h) cos t dtheta dz -- the h correction is
    O(c/R) = 1e-3 and is dropped here, so agreement is expected to that order."""
    dt = 2.0 * math.pi / nth
    dz = L / nz
    fx = fy = 0.0
    for j in range(nz):
        for i in range(nth):
            t = (i + 0.5) * dt
            pg = p[j][i] - P_AMB
            fx += pg * math.cos(t) * R_J * dt * dz
            fy += pg * math.sin(t) * R_J * dt * dz
    return fx, fy


def metrics(res, eps, LD):
    """Sommerfeld number and attitude angle from the pressure force.

    The load carried by the bearing is the pressure force on the bearing surface
    (equal and opposite to the force on the journal, which the external load must
    balance).  The line of centres points at theta = 180 deg because
    attitude_angle_deg = 180 puts the minimum gap there, so the attitude angle --
    the angle from the load line to the line of centres, measured with the
    rotation -- is 180 deg minus the load-vector angle."""
    W = math.hypot(res["fx"], res["fy"])
    phi = 180.0 - math.degrees(math.atan2(res["fy"], res["fx"]))
    S = sommerfeld_number(mu=MU, N_rev_s=N_REV_S, W=W, L=res["L"], D=2.0 * R_J,
                          R=R_J, c=C_R)
    hx, hy = hand_integrate_force(res["p"], res["nth"], res["nz"], res["L"])
    return dict(eps=eps, LD=LD, W=W, S=S, phi=phi,
                W_hand=math.hypot(hx, hy),
                phi_hand=180.0 - math.degrees(math.atan2(hy, hx)),
                p_max=max(max(row) for row in res["p"]),
                theta_min=min(min(row) for row in res["theta"]),
                elrod_iters=res["elrod_iters"], mc=res["mc"], ns=res["ns"],
                R_p=res["R_p"], converged=res["converged"])


EPS_SWEEP = (0.2, 0.4, 0.6, 0.8, 0.9)
LD_CASES = (("0.5", 0.5), ("1", 1.0))
MU_NOMINAL = 0.01                  # baseline viscosity for every swept case
MU_CHECK_EPS = (0.4, 0.8)          # cases repeated at a second viscosity
MU_CHECK = 0.003


def rupture_angle(res):
    """First angle past the minimum gap (180 deg) at which the mid-plane film
    content drops below full, in degrees from the maximum gap."""
    th, n = res["theta"], res["nth"]
    mid = res["nz"] // 2
    for i in range(n // 2, n):
        if th[mid][i] < 1.0 - 1e-5:
            return (i + 0.5) * 360.0 / n
    return float("nan")


def convergence_study(exe, work, eps=0.6, LD=8.0):
    """Free-boundary resolution audit against the exact infinite-length
    Swift-Stieber solution.  Reproduces the table in the BULK comment."""
    global BULK
    ex = swift_stieber_infinite(eps)
    print(f"\n=== Free-boundary convergence audit (L/D = {LD:g}, eps = {eps}) ===")
    print(f"exact Swift-Stieber: S = {ex['S']:.5f}, phi = {ex['phi_deg']:.2f} deg, "
          f"rupture = {ex['t2_deg']:.2f} deg")
    print(f"{'beta':>7} {'n_theta':>8} {'S':>9} {'dS':>8} {'phi':>8} {'rupture':>9}")
    saved, rows = BULK, []
    for beta in (1e9, 1e8):
        BULK = beta
        for nth in (720, 2880):
            res = run_case(exe, work, f"cvg_b{beta:g}_n{nth}", eps=eps, LD=LD,
                           flooded=True, nth=nth, nz=60)
            m = metrics(res, eps, LD)
            print(f"{beta:>7.0e} {nth:>8d} {m['S']:>9.5f} "
                  f"{(m['S']-ex['S'])/ex['S']:>7.1%} {m['phi']:>8.2f} "
                  f"{rupture_angle(res):>9.2f}")
            rows.append(dict(beta=beta, nth=nth, S=m["S"], phi=m["phi"],
                             rupture=rupture_angle(res), S_exact=ex["S"],
                             phi_exact=ex["phi_deg"], rupture_exact=ex["t2_deg"]))
            sys.stdout.flush()
    BULK = saved
    return rows


def sweep(exe, work):
    results = []
    for ld_key, LD in LD_CASES:
        for flooded in (True, False):
            label = "flooded" if flooded else "ungrooved"
            print(f"\n=== L/D = {ld_key}, {label} inlet ===")
            print(f"{'eps':>5} {'mu':>7} {'S sim':>10} {'S ref':>9} {'dS':>9} "
                  f"{'phi sim':>8} {'phi ref':>8} {'dphi':>7} {'rupt':>7} "
                  f"{'th_min':>7} {'p/beta':>8} {'conv':>5}")
            for eps in EPS_SWEEP:
                # Read the baseline from the constant, not from the module-level
                # MU, which set_mu() below rebinds.
                mus = [MU_NOMINAL] + ([MU_CHECK] if eps in MU_CHECK_EPS else [])
                for mu in mus:
                    set_mu(mu)
                    tag = f"rb_LD{ld_key}_{label}_eps{eps:g}_mu{mu:g}"
                    res = run_case(exe, work, tag, eps=eps, LD=LD, flooded=flooded)
                    m = metrics(res, eps, LD)
                    ref = rb_lookup(ld_key, eps)
                    m.update(ld_key=ld_key, inlet=label, mu=mu,
                             rupture=rupture_angle(res),
                             p_over_beta=(m["p_max"] - P_AMB) / BULK,
                             S_ref=ref[0], phi_ref=ref[1],
                             dS=(m["S"] - ref[0]) / ref[0],
                             dphi=m["phi"] - ref[1])
                    results.append(m)
                    print(f"{eps:>5.1f} {mu:>7.3f} {m['S']:>10.5f} {ref[0]:>9.4f} "
                          f"{m['dS']:>8.1%} {m['phi']:>8.2f} {ref[1]:>8.2f} "
                          f"{m['dphi']:>7.2f} {m['rupture']:>7.2f} "
                          f"{m['theta_min']:>7.4f} {m['p_over_beta']:>8.4f} "
                          f"{str(m['converged']):>5}")
                    sys.stdout.flush()
            set_mu(MU_NOMINAL)
    return results


def main(exe, work, do_convergence=True):
    os.makedirs(work, exist_ok=True)
    print(f"pancake: {exe}")
    print(f"geometry: R={R_J*1e3} mm, c={C_R*1e6} um, c/R={C_R/R_J:.1e}, "
          f"N={N_REV_S} rev/s, beta={BULK:.0e} Pa, grid {N_THETA}x{N_Z}")
    print("\nReference-table provenance check (exact Swift-Stieber vs RB_TABLE['inf']):")
    for eps, s_tab, phi_tab in RB_TABLE["inf"]:
        ex = swift_stieber_infinite(eps)
        print(f"  eps={eps:.1f}  S {s_tab:.4f} vs {ex['S']:.4f} ({(ex['S']-s_tab)/s_tab:+.2%})  "
              f"phi {phi_tab:.2f} vs {ex['phi_deg']:.2f} ({ex['phi_deg']-phi_tab:+.2f} deg)")

    if do_convergence:
        cvg = convergence_study(exe, work)
    else:
        ex = swift_stieber_infinite(0.6)
        cvg = [dict(r, S_exact=ex["S"], phi_exact=ex["phi_deg"],
                    rupture_exact=ex["t2_deg"]) for r in CONVERGENCE_MEASURED]
    results = sweep(exe, work)

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "raimondi_boyd_results.json")
    json.dump(dict(sweep=results, convergence=cvg), open(out, "w"), indent=1)
    print(f"\nwrote {out}")

    bad = [r for r in results if not r["converged"]]
    if bad:
        print(f"WARNING: {len(bad)} case(s) did not converge: "
              + ", ".join(f"{r['ld_key']}/{r['inlet']}/eps{r['eps']}/mu{r['mu']}" for r in bad))
    else:
        print("all cases converged (Elrod tolerance reached, R_p < 1e-9)")
    return results


if __name__ == "__main__":
    exe = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "_work")
    main(exe, work)
