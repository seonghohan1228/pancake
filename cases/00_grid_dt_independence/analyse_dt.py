"""00b timestep-study post-processing: equilibrium state, settling time, monotonicity.

Reads metrics.json (written by extract.py) and reports, per dt run:
  - final eccentricity ratio and bearing-offset angle
  - settling time: first VTK output time after which |e(t) - e_final| < 0.01
    (1% of the clearance, since e is already normalised by c) for the rest of
    the record. Resolution is the config write_interval.
  - trajectory character: number of sign changes of de/dt over the record
    (0 -> monotone, >0 -> oscillatory) and the overshoot above e_final.
"""

import json
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
RUNS = ["dt_T1", "dt_T2", "dt_T3", "dt_T4", "dt_T5"]
T_WINDOW = 0.15  # common analysis window [s]; runs made before end_t was frozen
                 # at 0.15 s carry extra frames that are discarded here.
TOL = 1e-9


def analyse(entry):
    hist = [h for h in entry.get("history", [])
            if h.get("ecc_ratio") is not None and h["t"] <= T_WINDOW + TOL]
    if not hist:
        return None
    t = [h["t"] for h in hist]
    e = [h["ecc_ratio"] for h in hist]
    e_final = e[-1]

    settling = None
    for k in range(len(e)):
        if all(abs(v - e_final) < 0.01 for v in e[k:]):
            settling = t[k]
            break

    de = [e[i + 1] - e[i] for i in range(len(e) - 1)]
    nz = [d for d in de if abs(d) > 1e-9]
    flips = sum(1 for i in range(len(nz) - 1) if nz[i] * nz[i + 1] < 0)
    overshoot = max(e) - e_final

    return {
        "n_frames": len(hist),
        "t_last": t[-1],
        "e_final": e_final,
        "angle_final": hist[-1].get("bearing_angle_deg"),
        "e_max": max(e),
        "overshoot": overshoot,
        "settling_time_s": settling,
        "slope_sign_changes": flips,
        "character": "monotone" if flips == 0 else "oscillatory",
        "history_t": t,
        "history_e": e,
    }


def trajectory_deviation(out):
    """max |e_run(t) - e_ref(t)| over the frame times shared with the finest dt.

    Sharper than the settling time, whose resolution is limited by
    write_interval (5 ms). Reference is the smallest available dt.
    """
    avail = [r for r in RUNS if r in out and out[r].get("history_t")]
    if not avail:
        return
    ref = avail[-1]
    ref_map = dict(zip(out[ref]["history_t"], out[ref]["history_e"]))
    for r in avail:
        common = [(t, e) for t, e in zip(out[r]["history_t"], out[r]["history_e"])
                  if any(abs(t - tr) < 1e-9 for tr in ref_map)]
        dev = [abs(e - ref_map[min(ref_map, key=lambda x: abs(x - t))]) for t, e in common]
        out[r]["ref_run"] = ref
        out[r]["max_traj_dev_vs_ref"] = max(dev) if dev else None
        out[r]["e_final_dev_vs_ref"] = abs(out[r]["e_final"] - out[ref]["e_final"])


def main():
    m = json.loads((CASE_DIR / "metrics.json").read_text())
    out = {}
    for name in RUNS:
        if name not in m:
            continue
        entry = m[name]
        a = analyse(entry) or {}
        a["dt"] = entry["dt"]
        a["end_t"] = entry["end_t"]
        a["diagnostics"] = entry.get("diagnostics", {})
        a["final_frame"] = entry.get("final", {})
        out[name] = a
        d = a.get("diagnostics", {})
        print(f"{name}: dt={entry['dt']:.2e} frames={a.get('n_frames')} "
              f"t_last={a.get('t_last')} e_final={a.get('e_final')} "
              f"settle={a.get('settling_time_s')} {a.get('character')} "
              f"steps={d.get('final_step')} t_diag={d.get('final_time')}")
    trajectory_deviation(out)
    (CASE_DIR / "dt_analysis.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
