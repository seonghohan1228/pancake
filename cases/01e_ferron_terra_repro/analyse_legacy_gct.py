"""Report-only analysis of the legacy GCT reproduction."""
import json
from pathlib import Path
CASE_DIR=Path(__file__).resolve().parent
LIMITS={"W_rel_pct":1.,"torque_rel_pct":1.,"T_max_K":.1}
def change(a,b,key): return abs(a[key]-b[key])/max(abs(b[key]),1e-300)*100
def main():
    raw=json.loads((CASE_DIR/"legacy_gct_runs.json").read_text()); groups={"isothermal":[f"legacy_gct_g{i}" for i in range(1,6)],"nominal_thd":["legacy_gct_thd_g3","legacy_gct_thd_g4"]}; out={"classification":"REPORT_ONLY_REPAIRED_LEGACY","limits":LIMITS,"groups":{},"repair_ledger":json.loads((CASE_DIR/"legacy_repair_ledger.json").read_text())}
    for group,names in groups.items():
        rows=[]
        for name in names:
            if name not in raw: continue
            m=raw[name].get("metrics") or {}; ok=raw[name]["rc"]==0 and m.get("flow_converged",False) and m.get("coupling_converged",False); rows.append({"name":name,"grid":raw[name]["grid"],"rc":raw[name]["rc"],"status":"PASS" if ok else "NOT_CONVERGED","W_N":m.get("W"),"p_max_MPa":m.get("p_max",0)/1e6,"torque_Nm":m.get("friction_torque"),"T_max_C":m.get("T_max",0)-273.15 if "T_max" in m else None,"coupling_R_p":m.get("coupling_R_p"),"coupling_R_T":m.get("coupling_R_T")})
        comparisons=[]
        for a,b in zip(rows,rows[1:]):
            item={"coarse":a["name"],"fine":b["name"],"W_rel_pct":change(a,b,"W_N"),"torque_rel_pct":change(a,b,"torque_Nm")}
            if a["T_max_C"] is not None and b["T_max_C"] is not None: item["T_max_K"]=abs(a["T_max_C"]-b["T_max_C"])
            comparisons.append(item)
        out["groups"][group]={"runs":rows,"consecutive_changes":comparisons}
    (CASE_DIR/"legacy_gct_analysis.json").write_text(json.dumps(out,indent=2)); print("LEGACY GCT: REPORT ONLY")
if __name__=="__main__": main()
