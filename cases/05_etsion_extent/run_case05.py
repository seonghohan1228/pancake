"""Run the case-05 matrix. Works on Windows (local, 4 ranks) and afml-r6 (24).

Usage:
    python run_case05.py smoke              # coarse 180x48, one condition, fbns
    python run_case05.py matrix             # the 16 production runs
    python run_case05.py gate               # the one cross-platform gate condition
    python run_case05.py sens               # pcav-tracking + bulk-modulus sensitivity
    python run_case05.py matrix --ranks 24  # explicit rank count

Every run gets its own directory holding the exact config.txt that produced it,
its run.log, and its results/ tree -- nothing is overwritten (campaign rule).
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import conditions as C  # noqa: E402
import extent as X      # noqa: E402

CASE_DIR = Path(__file__).resolve().parent
PROGRESS = CASE_DIR / "progress.log"

IS_WINDOWS = platform.system() == "Windows"
if IS_WINDOWS:
    EXE = Path(r"C:\Dev\pancake\build-windows-mingw\pancake.exe")
    MINGW_BIN = r"C:\msys64\mingw64\bin"
    DEFAULT_RANKS = 4
else:
    EXE = Path.home() / "pancake" / "build-linux" / "pancake"
    MINGW_BIN = None
    DEFAULT_RANKS = 24


def progress(message):
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp}  {message}"
    with open(PROGRESS, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def run_solver(work_dir, ranks, log_name="run.log"):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    if MINGW_BIN:
        env["PATH"] = MINGW_BIN + os.pathsep + env["PATH"]
        cmd = ["mpiexec", "-n", str(ranks), str(EXE), "--config", "config.txt"]
    else:
        cmd = ["mpirun", "-n", str(ranks), "--bind-to", "core",
               str(EXE), "--config", "config.txt"]
    t0 = time.perf_counter()
    with open(work_dir / log_name, "w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.run(cmd, cwd=str(work_dir), stdout=log,
                              stderr=subprocess.STDOUT, env=env, text=True,
                              encoding="utf-8", errors="replace")
    return proc.returncode, time.perf_counter() - t0


def do_run(run_dir, cfg_text, ranks, label):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "config.txt").write_text(cfg_text, encoding="utf-8")
    rc, secs = run_solver(run_dir, ranks)
    if rc != 0:
        progress(f"{label} FAILED rc={rc} ({secs:.0f}s)")
        return None
    m = X.cavity_metrics(run_dir)
    if m is None:
        progress(f"{label} produced no output ({secs:.0f}s)")
        return None
    m["wall_s"] = secs
    m["ranks"] = ranks
    m["coupling"] = X.coupling_health(run_dir)
    d = X.diagnostics_tail(run_dir)
    for k in ("liquid_residual", "gas_residual", "mass_residual"):
        if k in d:
            m[k] = float(d[k])
    progress(f"{label} ok {secs:.0f}s  cav={m['cav_fraction']:.3f} "
             f"th_st={m['theta_st_any_z']} th_end={m['theta_end_any_z']} "
             f"ext={m['extent_any_z']:.1f} p_max={m['p_max']:.4g} "
             f"alpha_max={m['alpha_gas_max']:.4g}")
    return m


def save(results, name):
    path = CASE_DIR / name
    if path.exists():
        stamp = time.strftime("%Y%m%d_%H%M%S")
        path.rename(path.with_name(f"{path.stem}_superseded_{stamp}{path.suffix}"))
    path.write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    progress(f"wrote {path.name}")


def cmd_smoke(args):
    rpm, ps = 1840, 13.6
    cfg = C.render_config(rpm, ps, "fbns", n_theta=180, n_z=48)
    m = do_run(CASE_DIR / "smoke" / f"{C.tag(rpm, ps)}_fbns", cfg, args.ranks, "smoke/fbns")
    cfg = C.render_config(rpm, ps, "jfo", n_theta=180, n_z=48)
    m2 = do_run(CASE_DIR / "smoke" / f"{C.tag(rpm, ps)}_jfo", cfg, args.ranks, "smoke/jfo")
    save({"fbns": m, "jfo": m2}, "smoke.json")


def cmd_matrix(args):
    out = {}
    for model in ("fbns", "jfo"):
        for rpm, ps in C.conditions():
            key = f"{model}/{C.tag(rpm, ps)}"
            run_dir = CASE_DIR / model / C.tag(rpm, ps)
            cfg = C.render_config(rpm, ps, model)
            out[key] = do_run(run_dir, cfg, args.ranks, key)
            save(out, "matrix.json")


def cmd_vapour(args):
    """Contrast leg: cavitation datum at the OIL VAPOUR PRESSURE (Table 1, 13.3 Pa).

    Two things at once. (i) It is the "what if this were vaporous cavitation"
    control, with no fitted number -- 13.3 Pa is printed in the paper. (ii) It is
    the only steady configuration in which pancake's TWO_PHASE outgassing is not
    inert: with the datum at the air-saturation pressure, complementarity keeps
    p >= p_sat everywhere and the equilibrium flash releases nothing, so alpha_gas
    is identically zero. Dropping the datum to p_vap opens a 101 kPa band in which
    dissolved air genuinely comes out of solution while the film is still full,
    which is where FBNS (moved ceiling, gas volumetrically ACTIVE) and JFO
    (fixed ceiling, gas volumetrically INERT) can finally differ by gas.
    """
    out = {}
    for model in ("fbns", "jfo"):
        for rpm, ps in C.conditions():
            if model == "jfo" and ps not in (13.6, 54.4):
                continue  # JFO spot-checks only; FBNS carries the full sweep
            key = f"vapour_{model}/{C.tag(rpm, ps)}"
            cfg = C.render_config(rpm, ps, model, p_cav=C.P_VAPOUR,
                                  cavitation_threshold="SCALAR_PCAV")
            # The moved-ceiling feedback loop alpha_gas -> theta_full -> theta -> p
            # -> alpha_gas is stiff on an open release band (PHYSICS.md 7.3.7,
            # "Where this is stiff"). At the default damping this leg exhausts
            # coupling_max_iters with R_p ~ 0.27 and the best-iterate rollback
            # fires. These three settings are the documented lever and they turn
            # it into a clean fixed point (R_p ~ 1e-7 in ~220 outer iterations).
            cfg = (cfg.replace("void_coupling_relaxation = 0.3",
                               "void_coupling_relaxation = 0.05")
                      .replace("coupling_relaxation = 1.0",
                               "coupling_relaxation = 0.3")
                      .replace("coupling_max_iters = 300",
                               "coupling_max_iters = 2000"))
            out[key] = do_run(CASE_DIR / f"vapour_{model}" / C.tag(rpm, ps),
                              cfg, args.ranks, key)
            save(out, "vapour.json")


def cmd_gate(args):
    rpm, ps = 1840, 27.2
    cfg = C.render_config(rpm, ps, "fbns")
    host = "r6" if not IS_WINDOWS else "win"
    key = f"gate_{host}"
    m = do_run(CASE_DIR / "gate" / key, cfg, args.ranks, key)
    save({key: m}, f"gate_{host}.json")


def cmd_sens(args):
    out = {}
    # (a) p_cav tracking the ambient: the datum degeneracy demonstration (README 2.5)
    for rpm, ps in ((1840, 13.6), (1840, 54.4)):
        p_amb = C.P_ATM + ps * 1000.0
        cfg = C.render_config(rpm, ps, "fbns", p_cav=p_amb)
        key = f"sens_pcav_tracks_ambient/{C.tag(rpm, ps)}"
        out[key] = do_run(CASE_DIR / "sens_pcav_tracks_ambient" / C.tag(rpm, ps),
                          cfg, args.ranks, key)
        save(out, "sens.json")
    # (b) bulk modulus, JFO leg only (the one UNVERIFIED number it depends on)
    for beta in (1.0e9, 1.7e9):
        rpm, ps = 1840, 27.2
        cfg = C.render_config(rpm, ps, "jfo", bulk_modulus=beta)
        key = f"sens_beta/beta_{beta:.1e}"
        out[key] = do_run(CASE_DIR / "sens_beta" / f"beta_{beta:.1e}",
                          cfg, args.ranks, key)
        save(out, "sens.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["smoke", "matrix", "gate", "sens", "vapour"])
    ap.add_argument("--ranks", type=int, default=DEFAULT_RANKS)
    args = ap.parse_args()
    progress(f"=== {args.what} on {platform.node()} ({platform.system()}), "
             f"{args.ranks} ranks, exe={EXE} ===")
    {"smoke": cmd_smoke, "matrix": cmd_matrix, "gate": cmd_gate,
     "sens": cmd_sens, "vapour": cmd_vapour}[args.what](args)


if __name__ == "__main__":
    main()
