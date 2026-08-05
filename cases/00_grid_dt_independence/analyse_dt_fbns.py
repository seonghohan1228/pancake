"""Case 00b (FBNS rerun) post-processing: equilibrium, settling time, force balance.

Reuses extract.py's VTK/diagnostics readers but targets the dt_T*_fbns runs on
the G6 grid (720x192). Writes dt_fbns_analysis.json in this directory.

Stats are reported twice: "full" over the whole record and "w060" over the
common 0 to 0.06 s window used to compare timesteps.

Per run:
  e_final        |xB| / c at the last written frame                       [-]
  angle_final    attitude angle of the bearing centre offset              [deg]
  settling_time  first output time after which |e - e_final| < 0.01       [s]
                 (e is normalised by c, so 0.01 == 1% of the clearance)
  overshoot      max(e) - e_final                                         [-]
  character      monotone if de/dt never changes sign, else oscillatory
  W_final        |pressure-force resultant| at the last frame             [N]
  W_err_pct      100 * (W_final - 6000) / 6000, the force-balance check    [%]

The T1 record shows the bearing does not reach a fixed point: it settles onto
a limit cycle (see SUMMARY.md), so the single-frame e_final / W_final are one
arbitrary sample of that cycle and are not reproducible between runs. The
post-settling averages are the comparable metrics:
  e_mean, e_std, e_peak_to_peak      eccentricity ratio over the cycle     [-]
  W_mean, W_std, W_mean_err_pct      force balance against the 6 kN load
  angle_mean, attitude_mean          offset and attitude angle             [deg]
The dt-to-dt relative change is computed from w060.e_mean, not from e_final.
"""

import json
import math
import re
import statistics
from pathlib import Path

import extract

CASE_DIR = Path(__file__).resolve().parent
RUNS = ["dt_T0_fbns", "dt_T1_fbns", "dt_T2_fbns",
        "dt_T3_fbns", "dt_T4_fbns", "dt_T5_fbns"]
LOAD_N = 6000.0
WINDOW_S = 0.06  # common analysis window agreed after the T1 settling check


def _uniform_vector(vts_path, name, ncomp=3, _cache={}):
    """First tuple of a cell field the solver fills uniformly (xB, Fp, F).

    extract._read_pvts parses every float of every piece, which is ~20 MB of
    ASCII per frame and far too slow over a 300-frame series. These fields are
    constant over the mesh, so one piece and its leading values are enough.
    """
    key = (str(vts_path), name)
    if key in _cache:
        return _cache[key]
    blob = vts_path.read_bytes()
    m = re.search(rb'Name="%s"[^>]*>' % name.encode(), blob)
    if not m:
        _cache[key] = None
        return None
    chunk = blob[m.end():m.end() + 400].split(b"<")[0].split()
    vals = [float(v) for v in chunk[:ncomp]]
    _cache[key] = vals
    return vals


def _fast_frame(pvts_path):
    """(xB, Fp) for one frame, read from the first piece only."""
    root = extract.ET.parse(pvts_path).getroot()
    piece = next(root.iter("Piece"), None)
    if piece is None:
        return {}
    src = pvts_path.parent / piece.get("Source")
    out = {}
    xb = _uniform_vector(src, "xB")
    if xb:
        out["bearing_x"], out["bearing_y"] = xb[0], xb[1]
        out["ecc_ratio"] = math.hypot(xb[0], xb[1]) / extract.CLEARANCE
        out["bearing_angle_deg"] = math.degrees(math.atan2(xb[1], xb[0]))
    fp = _uniform_vector(src, "Fp")
    if fp:
        out["W"] = math.hypot(fp[0], fp[1])
        out["force_angle_deg"] = math.degrees(math.atan2(fp[1], fp[0]))
    ang = _uniform_vector(src, "bearing_attitude_angle", ncomp=1)
    if ang:
        out["bearing_attitude_angle"] = ang[0]
    return out


def collect(name):
    run_dir = CASE_DIR / name
    if not run_dir.exists():
        return None
    cfg = run_dir / "config.txt"
    entry = {
        "dt": extract._cfg_value(cfg, "dt"),
        "end_t": extract._cfg_value(cfg, "end_t"),
        "n_theta": extract._cfg_value(cfg, "n_theta_global", int),
        "n_z": extract._cfg_value(cfg, "n_z_global", int),
        "cavitation_model": extract._cfg_value(cfg, "cavitation_model", str),
    }
    series = extract._series(run_dir)
    hist = []
    for t, path in series:
        fm = _fast_frame(path)
        hist.append({
            "t": t,
            "e": fm.get("ecc_ratio"),
            "angle": fm.get("bearing_angle_deg"),
            "attitude": fm.get("bearing_attitude_angle"),
            "W": fm.get("W"),
            "force_angle": fm.get("force_angle_deg"),
        })
    # the last frame is worth the full parse (p_max, h_min, cavitated area)
    if series:
        entry["final_full"] = extract._frame_metrics(series[-1][1])
    entry["n_frames"] = len(hist)
    entry["history"] = hist

    for tag, t_end in (("full", None), ("w060", WINDOW_S)):
        sub = [h for h in hist if h["e"] is not None
               and (t_end is None or h["t"] <= t_end + 1e-12)]
        if not sub:
            continue
        e = [h["e"] for h in sub]
        t = [h["t"] for h in sub]
        e_final = e[-1]
        settling = None
        for k in range(len(e)):
            if all(abs(v - e_final) < 0.01 for v in e[k:]):
                settling = t[k]
                break
        de = [e[i + 1] - e[i] for i in range(len(e) - 1)]
        nz = [d for d in de if abs(d) > 1e-9]
        flips = sum(1 for i in range(len(nz) - 1) if nz[i] * nz[i + 1] < 0)
        stats = {
            "t_last": t[-1],
            "n_frames": len(sub),
            "e_final": e_final,
            "e_max": max(e),
            "overshoot": max(e) - e_final,
            "settling_time_s": settling,
            "slope_sign_changes": flips,
            "character": "monotone" if flips == 0 else "oscillatory",
            "angle_final": sub[-1].get("angle"),
            "attitude_final": sub[-1].get("attitude"),
            "W_final": sub[-1].get("W"),
        }
        if sub[-1].get("W") is not None:
            stats["W_err_pct"] = 100.0 * (sub[-1]["W"] - LOAD_N) / LOAD_N

        # The run settles onto a limit cycle rather than a fixed point, so the
        # single-frame e_final / W_final are one arbitrary sample of that cycle.
        # Averages over the post-settling part of the window are the metrics
        # that can actually be compared between timesteps.
        if settling is not None:
            post = [h for h in sub if h["t"] >= settling]
            ep = [h["e"] for h in post]
            wp = [h["W"] for h in post if h.get("W") is not None]
            ap = [h["angle"] for h in post if h.get("angle") is not None]
            if ep:
                stats["e_mean"] = statistics.mean(ep)
                stats["e_std"] = statistics.pstdev(ep)
                stats["e_peak_to_peak"] = max(ep) - min(ep)
            if wp:
                stats["W_mean"] = statistics.mean(wp)
                stats["W_std"] = statistics.pstdev(wp)
                stats["W_mean_err_pct"] = 100.0 * (statistics.mean(wp) - LOAD_N) / LOAD_N
            if ap:
                stats["angle_mean"] = statistics.mean(ap)
                stats["attitude_mean"] = statistics.mean(ap) - 180.0
        entry[tag] = stats
    if "full" in entry:
        entry.update({k: entry["full"][k] for k in
                      ("e_final", "settling_time_s", "character", "W_final",
                       "overshoot", "angle_final", "t_last")})
    entry["diagnostics"] = extract._diagnostics(run_dir)
    return entry


def main():
    out = {}
    for name in RUNS:
        entry = collect(name)
        if entry is None:
            continue
        out[name] = entry
        d = entry.get("diagnostics", {})
        print(f"{name}: dt={entry['dt']:.2e} frames={entry['n_frames']} "
              f"t_last={entry.get('t_last')} e_final={entry.get('e_final')} "
              f"settle={entry.get('settling_time_s')} {entry.get('character')} "
              f"W={entry.get('W_final')} steps={d.get('final_step')} "
              f"iters_max={d.get('outer_iters_max')} "
              f"n_not_conv={d.get('n_not_converged')}", flush=True)

    # relative change of e_final versus the next-larger dt
    order = [r for r in RUNS if r in out and out[r].get("w060", {}).get("e_mean") is not None]
    for i, r in enumerate(order):
        if i == 0:
            out[r]["e_rel_change_vs_larger_dt_pct"] = None
        else:
            prev = out[order[i - 1]]["w060"]["e_mean"]
            cur = out[r]["w060"]["e_mean"]
            out[r]["e_rel_change_vs_larger_dt_pct"] = (
                100.0 * (cur - prev) / prev if prev else None)
    (CASE_DIR / "dt_fbns_analysis.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
