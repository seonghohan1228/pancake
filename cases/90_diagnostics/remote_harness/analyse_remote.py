"""Case 00b dt-sweep tail (T2/T3/T4) analysis for the afml-r6 remote runs.

Follows the protocol fixed by the local T1 gate:
  * the bearing settles onto a LIMIT CYCLE, not a fixed point, so equilibrium
    metrics are post-settling MEANS over t in [0.03, 0.06] s, never final-frame
    samples;
  * the key ladder question is whether the cycle amplitude is dt-independent
    (front/coupling artifact) or shrinks with dt (temporal resolution), so the
    cycle is characterised per-step from diagnostics.csv cavitated_fraction,
    which is the only per-step signal available (VTK frames are 25/50/100 steps
    apart and cannot resolve a 2*dt oscillation).

Reads <run_dir>/<results_sub>/flat/results.pvd and <run_dir>/<results_sub>/diagnostics.csv.
"""

import json
import math
import re
import statistics
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

CLEARANCE = 1.52e-4
LOAD_N = 6000.0
WINDOW = (0.03, 0.06)   # post-settling comparison window (coordinator)
FULL_END = 0.06


def cfg_value(cfg_path, key, cast=float):
    m = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", cfg_path.read_text(), re.M)
    return cast(m.group(1)) if m else None


def uniform_vector(vts_path, name, ncomp=3):
    """Leading tuple of a cell field the solver fills uniformly (xB, Fp, ...).

    Only the first piece and its first values are parsed: these fields are
    constant over the mesh and a full parse is ~20 MB of ASCII per frame.
    """
    blob = vts_path.read_bytes()
    m = re.search(rb'Name="%s"[^>]*>' % name.encode(), blob)
    if not m:
        return None
    chunk = blob[m.end():m.end() + 400].split(b"<")[0].split()
    return [float(v) for v in chunk[:ncomp]]


def frame(pvts_path):
    root = ET.parse(pvts_path).getroot()
    piece = next(root.iter("Piece"), None)
    if piece is None:
        return {}
    src = pvts_path.parent / piece.get("Source")
    out = {}
    xb = uniform_vector(src, "xB")
    if xb:
        out["e"] = math.hypot(xb[0], xb[1]) / CLEARANCE
        out["angle"] = math.degrees(math.atan2(xb[1], xb[0]))
    fp = uniform_vector(src, "Fp")
    if fp:
        out["W"] = math.hypot(fp[0], fp[1])
    att = uniform_vector(src, "bearing_attitude_angle", ncomp=1)
    if att:
        out["attitude"] = att[0]
    return out


def series(res_dir):
    pvd = res_dir / "flat" / "results.pvd"
    if not pvd.exists():
        return []
    root = ET.parse(pvd).getroot()
    items = [(float(ds.get("timestep")), pvd.parent / ds.get("file"))
             for ds in root.iter("DataSet")]
    items.sort(key=lambda x: x[0])
    return items


def read_diagnostics(res_dir):
    path = res_dir / "diagnostics.csv"
    if not path.exists():
        return None
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    head = lines[0].split(",")
    rows = [dict(zip(head, l.split(","))) for l in lines[1:]]
    return rows


def cycle_spectrum(rows, dt, t_lo):
    """Dominant frequency and amplitude of the per-step cavitated_fraction.

    Returns the peak of the one-sided amplitude spectrum over t >= t_lo, with
    the mean removed. f_nyquist = 1/(2*dt) is the step-rate Nyquist; a peak
    sitting at f_nyquist is the 2*dt numerical signature.
    """
    sub = [r for r in rows if float(r["time"]) >= t_lo]
    if len(sub) < 32:
        return None
    x = np.array([float(r["cavitated_fraction"]) for r in sub])
    n = len(x)
    x = x - x.mean()
    win = np.hanning(n)
    amp = np.abs(np.fft.rfft(x * win)) * 2.0 / win.sum()
    freq = np.fft.rfftfreq(n, d=dt)
    k = int(np.argmax(amp[1:]) + 1)
    order = np.argsort(amp[1:])[::-1][:5] + 1
    return {
        "n_samples": n,
        "t_window": [float(sub[0]["time"]), float(sub[-1]["time"])],
        "f_nyquist_hz": 1.0 / (2.0 * dt),
        "dominant_f_hz": float(freq[k]),
        "dominant_amp": float(amp[k]),
        "f_over_nyquist": float(freq[k]) / (1.0 / (2.0 * dt)),
        "cav_frac_mean": float(np.mean([float(r["cavitated_fraction"]) for r in sub])),
        "cav_frac_p2p": float(np.ptp([float(r["cavitated_fraction"]) for r in sub])),
        "cav_frac_std": float(np.std([float(r["cavitated_fraction"]) for r in sub])),
        "top5_hz": [float(freq[i]) for i in order],
        "top5_amp": [float(amp[i]) for i in order],
    }


def analyse(run_dir, results_sub="results"):
    run_dir = Path(run_dir)
    res = run_dir / results_sub
    cfg = run_dir / "config.txt"
    dt = cfg_value(cfg, "dt")
    entry = {
        "run": run_dir.name,
        "dt": dt,
        "end_t": cfg_value(cfg, "end_t"),
        "n_theta": cfg_value(cfg, "n_theta_global", int),
        "n_z": cfg_value(cfg, "n_z_global", int),
        "write_interval_steps": cfg_value(cfg, "write_interval") / dt,
    }

    hist = []
    for t, path in series(res):
        if t > FULL_END + 1e-12:
            continue
        f = frame(path)
        if f.get("e") is not None:
            hist.append({"t": t, **f})
    entry["n_frames"] = len(hist)
    if not hist:
        return entry

    t = [h["t"] for h in hist]
    e = [h["e"] for h in hist]
    e_final = e[-1]

    # settling: first output time after which |e - e_final| stays < 1% of c
    settling = None
    for k in range(len(e)):
        if all(abs(v - e_final) < 0.01 for v in e[k:]):
            settling = t[k]
            break
    de = [e[i + 1] - e[i] for i in range(len(e) - 1)]
    nz = [d for d in de if abs(d) > 1e-9]
    flips = sum(1 for i in range(len(nz) - 1) if nz[i] * nz[i + 1] < 0)

    entry["transient"] = {
        "t_last": t[-1],
        "e_final": e_final,
        "e_max": max(e),
        "overshoot": max(e) - e_final,
        "settling_time_s": settling,
        "slope_sign_changes": flips,
        "character": "monotone" if flips == 0 else "oscillatory",
    }

    win = [h for h in hist if WINDOW[0] - 1e-12 <= h["t"] <= WINDOW[1] + 1e-12]
    if win:
        ew = [h["e"] for h in win]
        ww = [h["W"] for h in win if h.get("W") is not None]
        aw = [h["angle"] for h in win if h.get("angle") is not None]
        tw = [h["attitude"] for h in win if h.get("attitude") is not None]
        s = {
            "window_s": list(WINDOW),
            "n_frames": len(win),
            "e_mean": statistics.mean(ew),
            "e_std": statistics.pstdev(ew),
            "e_peak_to_peak": max(ew) - min(ew),
        }
        if ww:
            s["W_mean"] = statistics.mean(ww)
            s["W_std"] = statistics.pstdev(ww)
            s["W_min"], s["W_max"] = min(ww), max(ww)
            s["W_mean_err_pct"] = 100.0 * (statistics.mean(ww) - LOAD_N) / LOAD_N
        if aw:
            s["angle_mean"] = statistics.mean(aw)
        if tw:
            s["attitude_mean"] = statistics.mean(tw)
            s["attitude_std"] = statistics.pstdev(tw)
        entry["w3060"] = s

    rows = read_diagnostics(res)
    if rows:
        outer = [int(r["outer_iters"]) for r in rows]
        conv = [int(r["converged"]) for r in rows]
        entry["diagnostics"] = {
            "n_steps": len(rows),
            "final_time": float(rows[-1]["time"]),
            "outer_iters_mean": statistics.mean(outer),
            "outer_iters_max": max(outer),
            "n_not_converged": sum(1 for c in conv if c == 0),
            "liquid_residual_max": max(abs(float(r["liquid_residual"])) for r in rows),
        }
        entry["cycle"] = cycle_spectrum(rows, dt, WINDOW[0])
    return entry


if __name__ == "__main__":
    out = {}
    for arg in sys.argv[1:]:
        run_dir, _, sub = arg.partition(":")
        entry = analyse(run_dir, sub or "results")
        out[entry["run"]] = entry
        c = entry.get("cycle") or {}
        w = entry.get("w3060") or {}
        tr = entry.get("transient") or {}
        d = entry.get("diagnostics") or {}
        print(f"\n=== {entry['run']}  dt={entry['dt']:.1e}  frames={entry['n_frames']} "
              f"steps={d.get('n_steps')} ===")
        print(f"  settling {tr.get('settling_time_s')} s   {tr.get('character')}   "
              f"overshoot {tr.get('overshoot'):.4f}   e_final {tr.get('e_final'):.5f}")
        print(f"  [0.03,0.06]: e_mean {w.get('e_mean'):.5f}  e_p2p {w.get('e_peak_to_peak'):.5f}  "
              f"attitude {w.get('attitude_mean')}")
        print(f"  W_mean {w.get('W_mean'):.1f} N  ({w.get('W_mean_err_pct'):+.2f}% vs 6000)  "
              f"W_std {w.get('W_std'):.1f}")
        if c:
            print(f"  cycle: f_dom {c['dominant_f_hz']:.1f} Hz  = {c['f_over_nyquist']:.3f} x Nyquist"
                  f" ({c['f_nyquist_hz']:.0f} Hz)   cav_frac p2p {c['cav_frac_p2p']:.4f}")
        print(f"  newton mean {d.get('outer_iters_mean'):.1f} max {d.get('outer_iters_max')}  "
              f"not_converged {d.get('n_not_converged')}")
    print()
    print(json.dumps(out, indent=1, default=float)[:0])  # keep stdout clean
    Path(sys.argv[0]).with_name("dt_tail_analysis.json").write_text(
        json.dumps(out, indent=1, default=float))
