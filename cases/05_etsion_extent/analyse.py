"""Score case 05 against Etsion & Ludwig Table 2 and emit the comparison table.

Reads matrix.json (the 16 production runs) plus the supporting legs, writes
comparison.csv and prints the table that goes into the report.

Angle errors are wrapped to (-180, 180] before being reported, because theta_end
crosses 0 deg in half the conditions and a naive difference would read 355 vs 5
as a 350 deg error instead of a 10 deg one.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import conditions as C  # noqa: E402

CASE_DIR = Path(__file__).resolve().parent
REF_UNCERTAINTY_DEG = 5.0  # the paper's own stated max-clearance location error


def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0


def load(name):
    p = CASE_DIR / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def rows():
    ref = C.read_reference()
    matrix = load("matrix_r6.json") or load("matrix.json")
    out = []
    for model in ("fbns", "jfo"):
        for rpm, ps in C.conditions():
            t = C.tag(rpm, ps)
            m = matrix.get(f"{model}/{t}")
            r = ref[t]
            if m is None:
                out.append({"model": model, "tag": t, "rpm": rpm, "ps_kPa": ps,
                            "status": "FAILED"})
                continue
            st, en = m["theta_st_any_z"], m["theta_end_any_z"]
            out.append({
                "model": model, "tag": t, "rpm": rpm, "ps_kPa": ps, "status": "ok",
                "th_st_exp": r["theta_st"], "th_st_sim": st,
                "th_st_err": wrap180(st - r["theta_st"]),
                "th_end_exp": r["theta_end"], "th_end_sim": en,
                "th_end_err": wrap180(en - r["theta_end"]),
                "ext_exp": r["extent"], "ext_sim": m["extent_any_z"],
                "ext_err": m["extent_any_z"] - r["extent"],
                "th_st_mid": m["theta_st_midplane"], "th_end_mid": m["theta_end_midplane"],
                "theta_0_sim": m["theta_0_deg"], "theta_0_exp": r["theta_0"],
                "p_min_exp_abs_kPa": 101.325 + r["p_min_kPa_gauge"],
                "p_min_sim_abs_kPa": m["p_min"] / 1000.0,
                "cav_frac": m["cav_fraction"], "p_max_kPa": m["p_max"] / 1000.0,
                "alpha_gas_max": m["alpha_gas_max"],
                "axial_width_mid_mm": 1000.0 * m.get("axial_width_at_cavity_mid_m", 0.0),
                "axial_width_max_mm": 1000.0 * m.get("axial_width_max_m", 0.0),
                "wall_s": m.get("wall_s"),
            })
    return out


def fmt(v, w=8, p=1):
    if v is None:
        return " " * (w - 1) + "-"
    if isinstance(v, str):
        return f"{v:>{w}}"
    return f"{v:{w}.{p}f}"


def main():
    rs = rows()
    lines = []
    hdr = (f"{'model':5} {'rpm':>5} {'p_s':>6} | {'th_st':>7}{'exp':>6}{'err':>7} |"
           f" {'th_end':>7}{'exp':>6}{'err':>7} | {'ext':>7}{'exp':>6}{'err':>7} |"
           f" {'cav':>6} {'pmax':>7} {'a_gas':>7}")
    lines.append(hdr)
    lines.append("-" * len(hdr))
    for r in rs:
        if r["status"] != "ok":
            lines.append(f"{r['model']:5} {r['rpm']:5d} {r['ps_kPa']:6.1f} | FAILED")
            continue
        lines.append(
            f"{r['model']:5} {r['rpm']:5d} {r['ps_kPa']:6.1f} |"
            f"{fmt(r['th_st_sim'],8)}{fmt(r['th_st_exp'],6,0)}{fmt(r['th_st_err'],7)} |"
            f"{fmt(r['th_end_sim'],8)}{fmt(r['th_end_exp'],6,0)}{fmt(r['th_end_err'],7)} |"
            f"{fmt(r['ext_sim'],8)}{fmt(r['ext_exp'],6,0)}{fmt(r['ext_err'],7)} |"
            f"{fmt(r['cav_frac'],7,3)}{fmt(r['p_max_kPa'],8,0)}{fmt(r['alpha_gas_max'],8,4)}")
    # per-model error summary
    lines.append("")
    for model in ("fbns", "jfo"):
        sel = [r for r in rs if r["model"] == model and r["status"] == "ok"]
        if not sel:
            continue
        def stats(key):
            v = [abs(r[key]) for r in sel]
            return sum(v) / len(v), max(v)
        a1, m1 = stats("th_st_err")
        a2, m2 = stats("th_end_err")
        a3, m3 = stats("ext_err")
        within = sum(1 for r in sel if abs(r["th_st_err"]) <= REF_UNCERTAINTY_DEG)
        within2 = sum(1 for r in sel if abs(r["th_end_err"]) <= REF_UNCERTAINTY_DEG)
        lines.append(f"{model.upper():5} mean|err| deg  th_st {a1:5.1f} (max {m1:5.1f})   "
                     f"th_end {a2:5.1f} (max {m2:5.1f})   extent {a3:5.1f} (max {m3:5.1f})   "
                     f"within +/-5 deg: th_st {within}/8, th_end {within2}/8")
    text = "\n".join(lines)
    print(text)
    (CASE_DIR / "comparison_table.txt").write_text(text + "\n", encoding="utf-8")

    keys = [k for k in rs[0] if k != "status"]
    with open(CASE_DIR / "comparison.csv", "w", encoding="utf-8") as fh:
        fh.write("# case 05 vs Etsion & Ludwig 1982 Table 2. Angles in deg, "
                 "errors wrapped to (-180,180]. Reference carries a common +/-5 deg\n"
                 "# systematic uncertainty on every angle (paper's own statement).\n")
        fh.write(",".join(keys) + "\n")
        for r in rs:
            fh.write(",".join("" if r.get(k) is None else str(r.get(k)) for k in keys) + "\n")


if __name__ == "__main__":
    main()
