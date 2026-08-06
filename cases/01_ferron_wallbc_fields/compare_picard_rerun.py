"""Compare a best-iterate-fix rerun against its archived iterate-200 predecessor.

Both sides are force-balance logs written by force_balance.py, so the comparison
is iterate by iterate at identical (eps, psi) and identical config. Anything that
moves is the best-iterate rollback, because nothing else in the protocol changed.

Usage:
    python compare_picard_rerun.py <archived_fb_log.json> <rerun_fb_log.json>

Prints one row per force-balance call and a worst-case summary. Exit status is
always 0; this is a reporting tool, not a gate.
"""

import json
import sys
from pathlib import Path

# metric -> (label, absolute tolerance below which the difference is "printed digits")
METRICS = {
    "W": ("W [N]", 5.0e-2),
    "force_angle_deg": ("angle [deg]", 1.0e-3),
    "p_max": ("p_max [Pa]", 1.0e2),
    "T_max": ("T_max [K]", 1.0e-3),
    "h_min": ("h_min [m]", 1.0e-10),
    "cav_frac": ("cav frac", 1.0e-6),
    "friction_torque": ("torque [N m]", 1.0e-4),
}


def load(path):
    return {r["dir"].split("/")[-1]: r for r in json.loads(Path(path).read_text())["iters"]}


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    old, new = load(sys.argv[1]), load(sys.argv[2])

    worst = {k: 0.0 for k in METRICS}
    n_paired = 0
    drifted = []
    print(f"{'iterate':12s} {'picard':>14s} " + " ".join(f"{lbl:>16s}" for lbl, _ in METRICS.values()))
    for name in sorted(set(old) & set(new)):
        a, b = old[name], new[name]
        if not (a.get("ok") and b.get("ok")):
            print(f"{name:12s} {'(failed run)':>14s}")
            continue
        # The Newton driver chooses the next (eps, psi) from the previous
        # iterate's forces, so once anything moves at all the search points
        # diverge. Only iterates evaluated at the SAME geometry are comparable.
        if abs(a["eps"] - b["eps"]) > 1e-9 or abs(a["psi"] - b["psi"]) > 1e-9:
            drifted.append((name, b["eps"] - a["eps"], b["psi"] - a["psi"]))
            continue
        n_paired += 1
        ca = (a.get("coupling") or {}).get("iters")
        cb = (b.get("coupling") or {}).get("iters")
        cells = []
        for k in METRICS:
            if k in a and k in b:
                d = b[k] - a[k]
                worst[k] = max(worst[k], abs(d))
                cells.append(f"{d:16.3e}")
            else:
                cells.append(f"{'-':>16s}")
        print(f"{name:12s} {str(ca) + '->' + str(cb):>14s} " + " ".join(cells))

    print(f"\nworst |rerun - archive| over the {n_paired} iterates at identical (eps, psi):")
    for k, (lbl, tol) in METRICS.items():
        flag = "within printed digits" if worst[k] <= tol else "*** MOVED ***"
        print(f"  {lbl:16s} {worst[k]:12.4e}   (report tol {tol:g})  {flag}")

    if drifted:
        print("\niterates the Newton search placed elsewhere (d = rerun - archive):")
        for name, de, dp in drifted:
            print(f"  {name:12s} d_eps = {de:+.3e}   d_psi = {dp:+.3e} deg")


if __name__ == "__main__":
    main()
