"""Score case-01e D-dt against pre-registered limits."""
import json
import math
from pathlib import Path
CASE_DIR = Path(__file__).resolve().parent
LIMITS = {"eps_final_rel_pct": .5, "W_rel_pct": 1., "attitude_final_deg": .25, "T_max_K": .10}
def rel(a,b): return abs(a-b) / max(abs(b), 1e-300) * 100
def main():
    path=CASE_DIR / "d_timestep_runs.json"
    if not path.exists(): raise SystemExit("no D-dt results")
    rows={}
    for name, run in json.loads(path.read_text()).items():
        m=run.get("metrics") or {}
        rows[name]={"dt_s":run["dt_s"],"rc":run["rc"],"eps_final":m.get("eps_final"),"W_N":m.get("W"),"attitude_final_deg":m.get("attitude_final_deg"),"T_max_C":m.get("T_max",float("nan"))-273.15,"final_time_s":m.get("final_time_s"),"flow_converged_all_steps":m.get("flow_converged_all_steps",False),"complete":run["rc"]==0 and m.get("final_time_s",0.0)>=0.06 and m.get("flow_converged_all_steps",False)}
    checks={}
    for coarse,fine,label in (("d_dt_t2","d_dt_t3","T2_to_T3"),("d_dt_t3","d_dt_t4","T3_to_T4")):
        a,b=rows.get(coarse),rows.get(fine)
        if not a or not b or not a["complete"] or not b["complete"]: checks[label]={"pass":False,"reason":"missing or unconverged run"}; continue
        v={"eps_final_rel_pct":rel(a["eps_final"],b["eps_final"]),"W_rel_pct":rel(a["W_N"],b["W_N"]),"attitude_final_deg":abs(a["attitude_final_deg"]-b["attitude_final_deg"]),"T_max_K":abs(a["T_max_C"]-b["T_max_C"])}
        checks[label]={"values":v,"pass":all(v[k]<=LIMITS[k] for k in v)}
    out={"limits":LIMITS,"runs":rows,"checks":checks,"pass":all(v["pass"] for v in checks.values())}
    (CASE_DIR / "d_discretization_analysis.json").write_text(json.dumps(out,indent=2))
    print("VERDICT:", "PASS" if out["pass"] else "FAIL")
if __name__ == "__main__": main()
