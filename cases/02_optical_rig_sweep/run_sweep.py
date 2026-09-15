"""Optical-rig operating-window sweep (case 02).

Pinned-eccentricity steady solves: eps is prescribed and the film load W is an
output, so each sweep point is one solve rather than a force-balance ladder.
Resumable: completed points are keyed in sweep_runs.json and never re-run.

Usage (from repo root):
    python3 cases/02_optical_rig_sweep/run_sweep.py RUN1
    python3 cases/02_optical_rig_sweep/run_sweep.py RUN2 RUN4 THD
    python3 cases/02_optical_rig_sweep/run_sweep.py BAND      # eps +/- 0.04
"""

import json
import math
import sys
import time
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
RUN_ROOT = CASE_DIR / "runs" / "production"
TEMPLATE = (CASE_DIR / "config_template.txt").read_text(encoding="utf-8")
LEDGER = CASE_DIR / "sweep_runs.json"
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S

# Rev B geometry, 2026-08-10. D = 49.3 mm journal in a 50.0 mm PMMA bore.
R_JOURNAL = 0.02465
C_NOMINAL = 350.0e-6
C_ALTERNATE = 250.0e-6
L_NOMINAL = 0.050             # L/D = 1.014
P_AMBIENT = 101325.0
GRID = (720, 192)             # G6, frozen campaign resolution

# Lab oil: white mineral oil, CAS 8020-83-5, ISO VG 46. nu(40) = 46 cSt gives
# nu(20) ~ 80 cSt at VI ~ 100, so mu(20) = 850 * 80e-6 = 0.068 Pa s. This is a
# grade-based conversion, not a measurement -- the MU suite brackets it.
MU_20 = 0.068
ANDRADE_E = 3556.0            # [K], fitted to 0.061 Pa s at 25 degC / 0.011 at 75 degC.
                              # No measured pair for this oil: THD results are
                              # qualitative only.
T_AMBIENT = 293.15

# Cavitation threshold settings [Pa abs]. "degassed" is the vaporous floor used
# in the exclusion argument; "iso" is the isopentane-doped charge (3.6 wt% in
# paraffin of MW ~350 gives x = 0.153, and Raoult against p_v = 77 kPa at 20 degC
# gives 12 kPa); the rest are prescribed gaseous release pressures.
P_CAV = {"degassed": 500.0, "iso": 1.2e4, "p020": 2.0e4, "p050": 5.0e4,
         "p080": 8.0e4, "p150": 1.5e5}

ISOTHERMAL = """temperature_model = ISOTHERMAL
fluid_property_model = CONSTANT
liquid_viscosity_model = CONSTANT
viscosity_temperature_coeff = 0
cross_film_thermal_model = NONE
journal_wall_temperature = 293.15
bearing_wall_temperature = 293.15
journal_heat_transfer = 0
bearing_heat_transfer = 0
"""

# THD corner check: Andrade mu(T), depth-averaged film, no conjugate bush.
# Adiabatic walls bound the heating from above; PMMA (k = 0.19) is close to it,
# so the steel journal is the only realistic heat path (h_journal > 0).
def thd_block(h_journal):
    return f"""temperature_model = ENERGY_EQUATION
fluid_property_model = SINGLE_PHASE
liquid_viscosity_model = EMPIRICAL
viscosity_temperature_coeff = {ANDRADE_E}
cross_film_thermal_model = NONE
journal_wall_temperature = 293.15
bearing_wall_temperature = 293.15
journal_heat_transfer = {h_journal}
bearing_heat_transfer = 0
"""


def omega_of(rpm):
    return 2.0 * math.pi * rpm / 60.0


def point(name, eps, rpm, p_cav_key, c=C_NOMINAL, L=L_NOMINAL,
          p_end=P_AMBIENT, thermal=ISOTHERMAL, mu=MU_20, run=""):
    return dict(name=name, run=run, eps=eps, rpm=rpm, p_cav_key=p_cav_key,
                p_cav=P_CAV[p_cav_key], c=c, L=L, p_end=p_end,
                thermal=thermal, mu=mu, omega=omega_of(rpm))


def tag(eps, rpm, key):
    return f"e{round(eps * 100):03d}_n{rpm:04d}_{key}"


def run1():
    out = []
    for eps in (0.5, 0.6, 0.7):
        for rpm in (500, 1000, 1500, 2000, 3000):
            for key in ("degassed", "iso", "p020", "p050", "p080"):
                out.append(point(f"r1_{tag(eps, rpm, key)}", eps, rpm, key, run="RUN1"))
    return out


def run2():
    """Clearance contrast: identical points at c = 250 um."""
    out = []
    for eps in (0.5, 0.7):
        for rpm in (1000, 2000):
            for key in ("degassed", "p080"):
                out.append(point(f"r2_{tag(eps, rpm, key)}", eps, rpm, key,
                                 c=C_ALTERNATE, run="RUN2"))
    return out


def run4():
    """Back-pressure lever: raise both ends and see whether the cavity survives."""
    out = []
    for p_end in (1.5e5, 2.0e5):
        for key in ("p080", "p150"):
            out.append(point(f"r4_p{int(p_end / 1e3):03d}_{key}", 0.7, 2000, key,
                             p_end=p_end, run="RUN4"))
    return out


def band():
    """eps +/- 0.02 uncertainty band: shim +/-3 um, bore +/-5 um, touch-zero
    +/-3 um, so sqrt(3^2 + 5^2 + 3^2) / 350 = +/-0.019."""
    out = []
    for eps in (0.48, 0.52, 0.58, 0.62, 0.68, 0.72):
        for rpm in (500, 1000, 1500, 2000, 3000):
            for key in ("degassed", "p080"):
                out.append(point(f"bd_{tag(eps, rpm, key)}", eps, rpm, key, run="BAND"))
    return out


def mu_sensitivity():
    """Cavity geometry is nonlinear in mu even though p and W are linear: the
    threshold is an absolute pressure, so changing mu slides the whole field
    past it. VG 46 is a grade label, not a measurement, so bracket it."""
    out = []
    for mu in (0.05, 0.068, 0.10, 0.15):
        for key in ("degassed", "p080"):
            out.append(point(f"mu_{round(mu * 1000):03d}_{key}", 0.6, 1500, key,
                             mu=mu, run="MU"))
    return out


def thd():
    """Worst-corner thermal check against the isothermal assumption."""
    out = []
    for h_journal in (0, 10000):
        block = thd_block(h_journal)
        for key in ("degassed", "p080"):
            out.append(point(f"th_h{h_journal}_{key}", 0.7, 3000, key,
                             thermal=block, run="THD"))
    return out


SUITES = {"RUN1": run1, "RUN2": run2, "RUN4": run4, "BAND": band, "THD": thd,
          "MU": mu_sensitivity}


def progress(message):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {message}"
    with (CASE_DIR / "progress.log").open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    print(line, flush=True)


def load_ledger():
    return json.loads(LEDGER.read_text()) if LEDGER.exists() else {}


def save_ledger(ledger):
    LEDGER.write_text(json.dumps(ledger, indent=2))


def evaluate(spec, ledger):
    if spec["name"] in ledger and ledger[spec["name"]].get("rc") == 0:
        return ledger[spec["name"]]
    run_dir = RUN_ROOT / spec["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.txt").write_text(TEMPLATE.format(
        c=spec["c"], e=spec["eps"] * spec["c"], L=spec["L"],
        n_theta=GRID[0], n_z=GRID[1], omega=spec["omega"], mu=spec["mu"],
        p_cav=spec["p_cav"], p_end=spec["p_end"], thermal_block=spec["thermal"]))
    rc, wall = S.run_case(run_dir)
    metrics = S.frame_metrics(run_dir) if rc == 0 else None
    row = {k: spec[k] for k in ("run", "eps", "rpm", "p_cav_key", "p_cav", "c", "L", "p_end", "mu", "omega")}
    row.update(dir=run_dir.relative_to(CASE_DIR).as_posix(), rc=rc, wall_s=round(wall, 1),
               metrics=metrics)
    ledger[spec["name"]] = row
    save_ledger(ledger)
    progress(f"{spec['name']}: rc={rc} wall={row['wall_s']}s "
             f"W={(metrics or {}).get('W', float('nan')):.1f} N")
    return row


def main():
    suites = sys.argv[1:] or ["RUN1"]
    unknown = [s for s in suites if s not in SUITES]
    if unknown:
        raise SystemExit(f"unknown suite(s): {unknown}; choose from {list(SUITES)}")
    specs = [s for name in suites for s in SUITES[name]()]
    progress(f"start {suites}: {len(specs)} points, grid {GRID}, mu_20={MU_20} Pa s")
    ledger = load_ledger()
    for spec in specs:
        evaluate(spec, ledger)
    failed = [k for k, v in ledger.items() if v.get("rc") != 0]
    progress(f"done {suites}: {len(specs)} points, {len(failed)} failed")
    if failed:
        progress(f"failed: {failed}")


if __name__ == "__main__":
    main()
