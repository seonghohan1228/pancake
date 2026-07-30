"""Run the analytical journal-bearing benchmarks against pancake.

Usage:  python run_benchmarks.py <path-to-pancake.exe> [workdir]

Case A -- infinitely-long (Sommerfeld) limit.  L/D swept upward with
FULL_SOMMERFELD (no cavitation), comparing the mid-plane p(theta) against the
closed form.  Isolates the circumferential operator; the residual is the axial
end effect, which must fall as L/D grows.

Case B -- short-bearing (Ocvirk) limit.  L/D swept downward with GUMBEL,
comparing p(theta, z) against the closed form.  Isolates the axial operator; the
residual is the neglected circumferential pressure flow, which must fall as
L/D shrinks (Ocvirk error is O((L/D)^2)).

Both are verification against exact solutions of the same PDE the solver
discretizes, so a passing result bounds discretization error but says nothing
about whether the Reynolds equation itself models the experiment.
"""
import math, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analytic import sommerfeld_long_pressure, ocvirk_short_pressure, film

R_J = 0.02475
C_R = 250e-6
MU = 0.0865
RHO = 868.8
EPS = 0.6
RPM = 1000.0
OMEGA = 2.0 * math.pi * RPM / 60.0
U_SURF = OMEGA * R_J

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
p_cav = {p_cav}
bulk_modulus = 1e9
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
property_reference_pressure = 1.0e5
property_reference_temperature = 298.15
solubility_table =
density_table =
viscosity_table =
motion_model = STATIC
min_film_thickness = 1e-09
stop_on_nonpositive_film = true
bc_z_south_type = DIRICHLET
bc_z_south_val = {p_ref}
bc_z_south_thermal = ZERO_GRADIENT
bc_z_north_type = DIRICHLET
bc_z_north_val = {p_ref}
bc_z_north_thermal = ZERO_GRADIENT
cavitation_model = {cav}
cavitation_threshold = SCALAR_PCAV
gas_pressure_coupling = NONE
max_outer_iters = 300
outer_tol = 1e-10
theta_min = 1e-6
log_outer_iters = false
coupling_max_iters = 60
coupling_tolerance = 1e-10
steady_gas_model = FROZEN
allow_frozen_gas_steady = true
output_dir = out
filename_prefix = solution
output_verbosity = SIMPLIFIED
output_write_3d = false
output_write_flat = true
output_fields = pressure, film_content, h
"""


def read_field(path, nth, nz, name):
    txt = open(path).read()
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="ascii"[^>]*>(.*?)</DataArray>' % name,
                  txt, re.S)
    vals = [float(v) for v in m.group(1).split()]
    return [vals[j * nth:(j + 1) * nth] for j in range(nz)]


def run(exe, work, tag, **kw):
    d = os.path.join(work, tag)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "config.txt"), "w").write(CFG.format(**kw))
    r = subprocess.run([exe, "--config", "config.txt"], cwd=d,
                       capture_output=True, text=True, timeout=600)
    f = os.path.join(d, "out", "flat", "processor0", "solution_0_0.vts")
    if not os.path.exists(f):
        raise RuntimeError(f"{tag} produced no output:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return read_field(f, kw["nth"], kw["nz"], "p")


def case_a(exe, work):
    print("\n=== Case A: infinitely-long (Sommerfeld) limit, FULL_SOMMERFELD ===")
    print("  mid-plane p(theta) vs closed form; residual = axial end effect")
    p_ref = 1.0e6
    nth, nz = 720, 40
    print(f"{'L/D':>6} {'max|err|/amp':>13} {'rms/amp':>10} {'p_max sim':>12} {'p_max exact':>12}")
    out = []
    for LD in (1.0, 2.0, 4.0, 8.0):
        L = LD * 2.0 * R_J
        p = run(exe, work, f"sommerfeld_LD{LD:g}", R=R_J, c=C_R, e=EPS * C_R, L=L,
                nth=nth, nz=nz, omega=OMEGA, mu=MU, rho=RHO, p_cav=1.0,
                p_ref=p_ref, cav="FULL_SOMMERFELD")
        mid = p[nz // 2]
        exact = [sommerfeld_long_pressure((i + 0.5) * 2 * math.pi / nth, mu=MU, U=U_SURF,
                                          R=R_J, c=C_R, eps=EPS, p_ref=p_ref)
                 for i in range(nth)]
        amp = max(exact) - min(exact)
        err = [abs(a - b) for a, b in zip(mid, exact)]
        rms = math.sqrt(sum(e * e for e in err) / len(err))
        print(f"{LD:>6.1f} {max(err)/amp:>13.4%} {rms/amp:>10.4%} "
              f"{max(mid)/1e5:>12.4f} {max(exact)/1e5:>12.4f}")
        out.append((LD, max(err) / amp, rms / amp))
    return out


def case_b(exe, work):
    print("\n=== Case B: short-bearing (Ocvirk) limit, GUMBEL ===")
    print("  full p(theta,z) vs closed form; residual = neglected circumferential flow, O((L/D)^2)")
    p_ref = 1.0e5
    nth, nz = 720, 60
    print(f"{'L/D':>6} {'max|err|/amp':>13} {'rms/amp':>10} {'p_max sim':>12} {'p_max exact':>12}")
    out = []
    for LD in (0.5, 0.25, 0.125, 0.0625):
        L = LD * 2.0 * R_J
        p = run(exe, work, f"ocvirk_LD{LD:g}", R=R_J, c=C_R, e=EPS * C_R, L=L,
                nth=nth, nz=nz, omega=OMEGA, mu=MU, rho=RHO, p_cav=p_ref,
                p_ref=p_ref, cav="GUMBEL")
        errs, exact_all, sim_all = [], [], []
        for j in range(nz):
            z = (j + 0.5) * L / nz
            for i in range(nth):
                t = (i + 0.5) * 2 * math.pi / nth
                ex = ocvirk_short_pressure(t, z, mu=MU, U=U_SURF, R=R_J, c=C_R,
                                           eps=EPS, L=L, p_ref=p_ref, gumbel=True)
                errs.append(abs(p[j][i] - ex)); exact_all.append(ex); sim_all.append(p[j][i])
        amp = max(exact_all) - min(exact_all)
        rms = math.sqrt(sum(e * e for e in errs) / len(errs))
        print(f"{LD:>6.4g} {max(errs)/amp:>13.4%} {rms/amp:>10.4%} "
              f"{max(sim_all)/1e5:>12.4f} {max(exact_all)/1e5:>12.4f}")
        out.append((LD, max(errs) / amp, rms / amp))
    return out


if __name__ == "__main__":
    exe = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "_work")
    os.makedirs(work, exist_ok=True)
    print(f"pancake: {exe}")
    print(f"geometry: R={R_J*1e3} mm, c={C_R*1e6} um, eps={EPS}, N={RPM:g} rpm, "
          f"mu={MU} Pa.s, U={U_SURF:.4f} m/s")
    case_a(exe, work)
    case_b(exe, work)
