"""FBNS gas-effect pairs at the 500 RPM LGE tester conditions (c = 0.274).

Each pair (P4v/P4g at 768.1 N, P5v/P5g at 2610.3 N) differs only in the gas model:
  v: SINGLE_PHASE, SCALAR_PCAV, p_cav = 500 Pa (vaporous cavitation only)
  g: TWO_PHASE, steady_gas_model = EQUILIBRIUM, LOCAL_PSAT (uniform bubble pressure
     of the feed composition at the feed temperature, CURVE_TABLE), floor p_cav = 500 Pa
FBNS as on the optical rig: isothermal at TC0, no cross-film or housing conduction,
constant liquid viscosity. The solver has no native load balance for FBNS with dissolved
gas, so the driver balances the load with a Newton loop of STATIC solves:
  v run: unknowns (mu, psi) at the eccentricity of the accepted COUPLED gas-off h_min
  g run: unknowns (eps, psi) at the viscosity of its v run
Converged when |F|/W is within FORCE_TOL of 1 and the force direction within
ANGLE_TOL_DEG of the load direction. Every solve is kept under <work>/<case>/runs.

Usage: python3 run_fbns.py --case P4v --work DIR --solver EXE --table fw68d.csv
       [--ranks N --rankfile FILE --mpi mpirun]
"""
import argparse
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02d_fw68d_steady_six"))
sys.path.insert(0, str(HERE.parent / "02f_lge_mixture_gas"))
import run_steady as driver  # noqa: E402
import run_gas as gas_driver  # noqa: E402

R, CLEARANCE, LENGTH = .0215, 42.5e-6, .040
N_THETA, N_Z = 720, 192
OMEGA = 52.35987755982988                      # 500 RPM
TC0_K = 297.5390033315531                      # X4/X5 TC0 (24.389 C), feed temperature
P_CHAMBER = 1657756.124499086                  # 16.58 bar abs, both ends and the groove
GROOVE_DEG, GROOVE_WIDTH_DEG = 90.0, 13.298312524795712
COMPOSITION = .274
RHO = 954.9 - .7284*(TC0_K - 273.15)           # LINEAR_TEMPERATURE density of the campaign at TC0
MU_GAS = 1.2652527835974629e-05                # R410A vapour viscosity of the accepted campaign
P_CAV = 500.0
FORCE_TOL, ANGLE_TOL_DEG = 1e-3, .1
FBNS_TOL, MASS_TOL = 1e-6, 1e-8
MAX_NEWTON = 30
# The solver accepts R410A only under COUPLED. The free-gas density is the only use of the
# species here, so PROPANE stands in; check_static records whether any free gas formed
# (it does not when the cavity sits at the uniform bubble pressure), i.e. whether the
# stand-in, and the table clamp needed for out-of-chart flash queries, touch the solution.
GAS_SPECIES = "PROPANE"

# Accepted COUPLED gas-off references (lge-thermal-v3-c274-20261114, mix/gas_off/<X>/accepted.json).
REFERENCE = {
    "X4": dict(load=768.0791197054391, eps=0.7887024018305416, psi=356.148545274,
               h_min=8.980200472028191e-06, torque=0.06863233503655095),
    "X5": dict(load=2610.276612252356, eps=0.9304830787327196, psi=316.55515572741024,
               h_min=2.954697817450341e-06, torque=0.12406863169122363),
}
CASES = {
    "P4v": dict(ref="X4", gas=False), "P4g": dict(ref="X4", gas=True, pair="P4v"),
    "P5v": dict(ref="X5", gas=False), "P5g": dict(ref="X5", gas=True, pair="P5v"),
}
OUTPUT_FIELDS = ("pressure,film_content,h,temperature,mu,alpha_gas,dissolved_gas,mass_residual,"
                 "fluid_force_x,fluid_force_y,fluid_force_z,pressure_force_x,pressure_force_y,"
                 "friction_torque,bearing_x,bearing_y,bearing_z")


def render_config(gas, eps, psi, mu, table):
    """Full solver config; the v and g runs differ only in the four gas-model keys."""
    values = dict(
        R=R, c=CLEARANCE, e=eps*CLEARANCE, L=LENGTH, attitude_angle_deg=psi % 360,
        n_theta_global=N_THETA, n_z_global=N_Z, omega=OMEGA, omega_ramp_time=0,
        solution_mode="STEADY_STATE", dt=2e-5, end_t=2e-5, write_interval=2e-5, diagnostics_interval=1,
        motion_model="STATIC", bearing_initial_from_attitude="true", min_film_thickness=1e-9,
        stop_on_nonpositive_film="true",
        theta_convection_scheme="UPWIND", thermal_convection_scheme="UPWIND", gas_convection_scheme="UPWIND",
        cavitation_model="FBNS", fbns_tol=1e-10, fbns_max_iters=200, theta_min=1e-6,
        max_outer_iters=300, outer_tol=1e-7,
        cavitation_threshold="LOCAL_PSAT" if gas else "SCALAR_PCAV", p_cav=P_CAV,
        gas_pressure_coupling="NONE", bulk_modulus=1e9,
        bc_z_south_type="DIRICHLET", bc_z_south_val=P_CHAMBER, bc_z_north_type="DIRICHLET",
        bc_z_north_val=P_CHAMBER, bc_z_south_thermal="CONSTANT", bc_z_north_thermal="CONSTANT",
        bc_z_south_temperature=TC0_K, bc_z_north_temperature=TC0_K,
        temperature_model="ISOTHERMAL", temperature_initial=TC0_K, temperature_reference=TC0_K,
        property_reference_temperature=293.15, property_reference_pressure=1e5,
        rho_cp=1883641.4504463263, thermal_conductivity=.13, cross_film_thermal_model="NONE",
        journal_heat_transfer=0, bearing_heat_transfer=0,
        fluid_property_model="TWO_PHASE" if gas else "SINGLE_PHASE",
        steady_gas_model="EQUILIBRIUM" if gas else "FROZEN",
        dissolved_gas_species=GAS_SPECIES, dissolved_gas_initial=COMPOSITION, dissolved_gas_max=1.0,
        solubility_model="CURVE_TABLE", solution_table_file=table, solution_table_clamp="true",
        dissolved_gas_henry_coeff=0, dissolved_gas_diffusivity=0, gas_mass_transfer_rate=1e4,
        solution_density_gas_coeff=0, solution_viscosity_gas_coeff=0,
        density_model="CONSTANT", rho=RHO, liquid_viscosity_model="CONSTANT", mu=mu,
        viscosity_temperature_coeff=0, viscosity_pressure_coeff=0,
        mixture_viscosity_model="MCADAMS", mu_gas=MU_GAS, gas_alpha_max=.9, gas_pressure_floor=1e3,
        bc_z_south_dissolved_gas=COMPOSITION, bc_z_north_dissolved_gas=COMPOSITION,
        bc_z_south_free_gas_mass_fraction=0, bc_z_north_free_gas_mass_fraction=0,
        groove_mixing_model="PINNED",
        inlet_groove=f"{GROOVE_DEG:g} {GROOVE_WIDTH_DEG!r} {P_CHAMBER!r} {TC0_K!r}",
        inlet_species=f"0 {COMPOSITION} 0",
        coupling_max_iters=200, coupling_tolerance=1e-6, coupling_relaxation=1.0,
        coupling_relaxation_mode="ADAPTIVE", coupling_relaxation_min=.1, coupling_relaxation_backoff=.5,
        coupling_relaxation_growth=1.05, coupling_keep_best_iterate="true", log_outer_iters="true",
        output_dir="results", filename_prefix="solution", output_verbosity="DETAILED",
        output_write_3d="false", output_write_flat="true", output_fields=OUTPUT_FIELDS)
    return "\n".join(f"{k} = {v!r}" if isinstance(v, float) else f"{k} = {v}" for k, v in values.items()) + "\n"


def check_static(run_dir, returncode, load, gas):
    """Gates of one STATIC solve. Load balance is the driver's job and is not checked here."""
    if returncode:
        raise ValueError(f"solver exit code {returncode}")
    log = (run_dir / "run.log").read_text(errors="replace")
    diag = driver.last_diagnostic(run_dir)
    if int(diag["converged"]) != 1:
        raise ValueError("FBNS did not converge")
    fbns = re.findall(r"fbns iter=\d+ \|G\|=" + driver.NUMBER + r" \|FB\|=" + driver.NUMBER, log)
    if not fbns or max(map(float, fbns[-1])) >= FBNS_TOL:
        raise ValueError(f"FBNS residual missing or not below {FBNS_TOL:g}")
    if gas and (int(diag["coupling_converged"]) != 1 or int(diag["coupling_restored_best"]) != 0):
        raise ValueError("FBNS/gas-equilibrium outer loop did not converge cleanly")
    mass_rel = float(diag["liquid_residual_rel"])
    if not math.isfinite(mass_rel) or abs(mass_rel) > MASS_TOL:
        raise ValueError(f"liquid balance relative residual {mass_rel} above {MASS_TOL:g}")
    f = driver.read_fields(run_dir)
    for key in ("p", "film_content", "h", "mu", "theta_rad", "z_m", "F", "friction_torque"):
        if key not in f or not len(f[key]) or not np.isfinite(f[key]).all():
            raise ValueError(f"missing/nonfinite {key} field")
    if float(np.min(f["h"])) <= 0 or float(np.min(f["mu"])) <= 0:
        raise ValueError("nonpositive film gap or viscosity")
    if float(np.min(f["p"])) < P_CAV - 1:
        raise ValueError("pressure below the cavitation floor")
    if float(np.min(f["film_content"])) < -1e-6 or float(np.max(f["film_content"])) > 1 + 1e-6:
        raise ValueError("film content outside [0, 1]")
    alpha = f.get("alpha_gas", np.zeros_like(f["p"]))
    if not gas and float(np.max(np.abs(alpha))) > 0:
        raise ValueError("free gas in a VC-only run")
    if gas and float(np.max(alpha + f["film_content"])) > 1 + 1e-6:
        raise ValueError("liquid plus free gas exceeds the gap")
    force = f["F"].reshape(-1, 3)
    if not np.allclose(force, force[0], rtol=1e-9, atol=1e-6):
        raise ValueError("total force differs across VTK cells")
    fx, fy, _ = map(float, force[0])
    torque = float(f["friction_torque"][0])
    if OMEGA*torque <= 0:
        raise ValueError("journal torque has inconsistent rotation sign")
    dissolved = f.get("dissolved_gas", np.full_like(f["p"], COMPOSITION))
    proxy = dict(alpha_gas_max=float(np.max(alpha)), dissolved_min=float(np.min(dissolved)),
                 dissolved_max=float(np.max(dissolved)))
    proxy["species_and_clamp_inert"] = bool(proxy["alpha_gas_max"] == 0 and
                                            np.allclose(dissolved, COMPOSITION, rtol=0, atol=1e-12))
    return dict(**proxy, Fx_N=fx, Fy_N=fy, torque_Nm=torque, h_min_m=float(np.min(f["h"])),
                p_min_Pa=float(np.min(f["p"])), p_max_Pa=float(np.max(f["p"])),
                mass_residual_rel=mass_rel, diagnostics=diag), f


def force_errors(result, load):
    """Relative magnitude error and direction error (deg) against the load (fluid force = -W e_y)."""
    fx, fy = result["Fx_N"], result["Fy_N"]
    magnitude = math.hypot(fx, fy)/load - 1
    angle = math.degrees(math.atan2(fy, fx)) + 90
    return magnitude, (angle + 180) % 360 - 180


class Case:
    def __init__(self, args):
        self.args, self.key = args, args.case
        self.spec = CASES[args.case]
        self.ref = REFERENCE[self.spec["ref"]]
        self.load = self.ref["load"]
        self.work = Path(args.work).resolve() / args.case
        (self.work / "runs").mkdir(parents=True, exist_ok=True)
        self.table = str(Path(args.table).resolve())
        self.count = len(list((self.work / "runs").iterdir()))
        self.history = []

    def solve(self, eps, psi, mu, label):
        config = render_config(self.spec["gas"], eps, psi, mu, self.table)
        digest = hashlib.sha256(config.encode()).hexdigest()[:8]
        for old in (self.work / "runs").glob(f"*_{digest}/record.json"):
            record = json.loads(old.read_text())
            if record.get("config_sha256_8") == digest:
                return record
        self.count += 1
        run = self.work / "runs" / f"{self.count:04d}_{label}_{digest}"
        run.mkdir()
        (run / "config.txt").write_text(config)
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        command = driver.mpi_command(self.args.mpi, self.args.solver, self.args.ranks, self.args.rankfile)
        start = time.monotonic()
        code, timed_out = driver.run_solver(command, run, env, self.args.timeout)
        record = dict(run_dir=str(run), config_sha256_8=digest, eps=eps, psi_deg=psi % 360, mu_Pa_s=mu,
                      wall_s=time.monotonic() - start, returncode=code, timed_out=timed_out)
        try:
            result, _ = check_static(run, code, self.load, self.spec["gas"])
            record.update(result, valid=True)
            record["force_mag_err"], record["force_angle_err_deg"] = force_errors(record, self.load)
        except Exception as error:  # a failed solve is kept and reported, never used as a result
            record.update(valid=False, error=str(error))
        driver.save_json(run / "record.json", record)
        print(f"[{self.key}] {run.name} eps={eps:.6f} psi={psi % 360:.4f} mu={mu:.6e} "
              + (f"dF={record['force_mag_err']:+.2e} dang={record['force_angle_err_deg']:+.3f} "
                 f"hmin={record['h_min_m']*1e6:.4f}um M={record['torque_Nm']:.5f} {record['wall_s']:.0f}s"
                 if record["valid"] else f"FAILED {record['error']}"), flush=True)
        return record

    def residual(self, record):
        return np.array([record["force_mag_err"], math.radians(record["force_angle_err_deg"])])

    def balance(self, a0, psi0):
        """Newton on (a, psi): a = mu (VC run, eps fixed) or eps (GC run, mu fixed)."""
        gas = self.spec["gas"]
        eps_fixed = None if gas else 1 - self.ref["h_min"]/CLEARANCE
        mu_fixed = self.pair_mu if gas else None

        def run(a, psi, label):
            eps, mu = (a, mu_fixed) if gas else (eps_fixed, a)
            if not 0 < eps < 1 or mu <= 0:
                raise ValueError(f"inadmissible iterate eps={eps} mu={mu}")
            record = self.solve(eps, psi, mu, label)
            self.history.append({k: record.get(k) for k in ("run_dir", "eps", "psi_deg", "mu_Pa_s", "valid",
                                                            "force_mag_err", "force_angle_err_deg", "error")})
            return record

        step_a = 2e-4 if gas else 1e-3*a0           # FD steps: eps or relative mu
        limit_a = .02 if gas else .2*a0
        x = np.array([a0, psi0])
        current = run(*x, "newton0")
        if not current["valid"]:
            raise ValueError(f"initial solve failed: {current['error']}")
        jacobian = None
        for iteration in range(1, MAX_NEWTON + 1):
            r = self.residual(current)
            if abs(current["force_mag_err"]) < FORCE_TOL and abs(current["force_angle_err_deg"]) < ANGLE_TOL_DEG:
                return current, iteration - 1
            if jacobian is None:
                columns = []
                for index, step in ((0, step_a), (1, .05)):
                    probe = x.copy()
                    probe[index] += step
                    record = run(*probe, f"fd{iteration}_{index}")
                    if not record["valid"]:
                        raise ValueError(f"finite-difference solve failed: {record['error']}")
                    columns.append((self.residual(record) - r)/step)
                jacobian = np.column_stack(columns)
            delta = -np.linalg.solve(jacobian, r)
            scale = min(1., limit_a/abs(delta[0]) if delta[0] else 1., 5./abs(delta[1]) if delta[1] else 1.)
            delta *= scale
            for _ in range(5):
                candidate = run(*(x + delta), f"newton{iteration}")
                if candidate["valid"] and np.linalg.norm(self.residual(candidate)) < np.linalg.norm(r):
                    break
                delta *= .5
            else:
                jacobian = None                     # fresh Jacobian at the same point
                continue
            dr = self.residual(candidate) - r
            jacobian = jacobian + np.outer(dr - jacobian @ delta, delta)/float(delta @ delta)  # Broyden
            x, current = x + delta, candidate
        raise ValueError(f"load balance did not converge in {MAX_NEWTON} Newton iterations")


def cavity_report(fields):
    """Void fraction 1 - film_content >= 0.5: angular extent, area fraction and groove distance."""
    nt, nz = len(np.unique(fields["theta_rad"])), len(np.unique(fields["z_m"]))
    order = np.lexsort((fields["theta_rad"], fields["z_m"]))
    void = (1 - fields["film_content"][order]).reshape(nz, nt)
    h = fields["h"][order].reshape(nz, nt)
    theta_axis = np.degrees(fields["theta_rad"][order].reshape(nz, nt)[0])
    theta_hmin = float(theta_axis[np.argmin(h[nz//2])])
    columns = (void >= .5).any(axis=0)
    out = dict(theta_hmin_deg=theta_hmin, bubble_area_fraction=float(np.mean(void >= .5)),
               void_max=float(void.max()))
    if not columns.any():
        return out
    if columns.all():
        out.update(bubble_start_deg=None, bubble_end_deg=None, note="bubble region spans the full circumference")
        return out
    # Longest circular run of columns, in the rotation direction (theta increasing, omega > 0).
    shift = int(np.argmin(columns))
    rolled = np.roll(columns, -shift)
    runs, start = [], None
    for index, flag in enumerate(np.append(rolled, False)):
        if flag and start is None:
            start = index
        elif not flag and start is not None:
            runs.append((start, index - 1))
            start = None
    first, last = max(runs, key=lambda item: item[1] - item[0])
    dtheta = 360./nt
    start_deg = float((theta_axis[(first + shift) % nt] - dtheta/2) % 360)
    end_deg = float((theta_axis[(last + shift) % nt] + dtheta/2) % 360)
    groove_edge = GROOVE_DEG - GROOVE_WIDTH_DEG/2
    out.update(bubble_start_deg=start_deg, bubble_end_deg=end_deg,
               bubble_start_from_hmin_deg=(start_deg - theta_hmin) % 360,
               bubble_end_from_hmin_deg=(end_deg - theta_hmin) % 360,
               bubble_end_to_groove_edge_deg=(groove_edge - end_deg) % 360,
               bubble_regions=len(runs), void_max=float(void.max()))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", required=True, choices=sorted(CASES))
    parser.add_argument("--work", required=True)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--mpi", default="mpirun")
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--rankfile")
    parser.add_argument("--timeout", type=float, default=6*3600)
    parser.add_argument("--mu0", type=float, default=.0112, help="VC seed viscosity [Pa s]")
    args = parser.parse_args()
    case = Case(args)
    started = time.monotonic()
    destination = case.work
    try:
        if case.spec["gas"]:
            pair = json.loads((case.work.parent / case.spec["pair"] / "accepted.json").read_text())["result"]
            case.pair_mu = pair["mu_Pa_s"]
            seed = (pair["eps"], pair["psi_deg"])
        else:
            seed = (args.mu0, case.ref["psi"])
        result, iterations = case.balance(*seed)
        run = Path(result["run_dir"])
        checked, fields = check_static(run, 0, case.load, case.spec["gas"])
        save = {k: v for k, v in fields.items() if k in ("theta_rad", "z_m", "p", "h", "film_content", "T", "mu",
                                                          "alpha_gas", "dissolved_gas", "F", "friction_torque")}
        np.savez_compressed(destination / "fields.npz", **save)
        sat = gas_driver.property_curves(Path(args.table))[1]
        p_sat = gas_driver.bubble_pressure(sat, COMPOSITION, TC0_K - 273.15)
        driver.save_json(destination / "accepted.json", dict(
            case=args.case, spec=case.spec, reference=case.ref, result=result, newton_iterations=iterations,
            cavity=cavity_report(fields), uniform_p_sat_Pa=p_sat if case.spec["gas"] else None,
            elapsed_seconds=time.monotonic() - started, history=case.history,
            driver_sha256=driver.sha256(__file__), fields_sha256=driver.sha256(destination / "fields.npz"),
            solver_sha256=driver.sha256(args.solver), ranks=args.ranks))
        print(f"[{args.case}] accepted {run.name}")
    except Exception as error:
        driver.save_json(destination / "failure.json", dict(case=args.case, error=str(error), history=case.history,
                                                            elapsed_seconds=time.monotonic() - started))
        print(f"[{args.case}] FAILED: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
