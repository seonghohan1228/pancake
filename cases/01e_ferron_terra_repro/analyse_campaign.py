"""Create the case-01e verdict ledger without reading prior case results."""

import csv
import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parent / "01_ferron_wallbc_fields"))
import solver_io as S

LIMITS = {"eps_rel_pct": .5, "psi_deg": .25, "p_max_rel_pct": 1., "mean_error_K": .1,
          "rms_error_K": .1, "peak_temperature_K": .25, "peak_location_deg": 2.}

def reference():
    rows=[]
    with (CASE_DIR.parent / "ref_data" / "ferron1983_fig4_2000rpm_4kN_temperature.csv").open() as f:
        for r in csv.reader(f):
            if r and not r[0].lstrip().startswith("#"): rows.append((float(r[0]),float(r[1])+273.15))
    return rows

def pressure_reference():
    rows=[]
    path=CASE_DIR.parent / "ref_data" / "ferron1983_fig4_2000rpm_4kN_pressure.csv"
    with path.open() as f:
        for r in csv.reader(f):
            if r and not r[0].lstrip().startswith("#"): rows.append((float(r[0]),float(r[1])))
    return rows

def interp(theta, values, x):
    pairs=sorted((a%360,b) for a,b in zip(theta,values)); xs=[pairs[-1][0]-360,*[p[0] for p in pairs],pairs[0][0]+360]; ys=[pairs[-1][1],*[p[1] for p in pairs],pairs[0][1]]; x%=360
    for a,b,u,v in zip(xs,xs[1:],ys,ys[1:]):
        if a<=x<=b: return u+(v-u)*(x-a)/(b-a)
    raise RuntimeError("periodic interpolation failed")

def circular(a,b): return abs((a-b+180)%360-180)
def pct(a,b): return abs(a-b)/max(abs(b),1e-300)*100

def score(name, result, ref):
    if "failed" in result: return {"label":result["label"],"status":"FAILED", "reason":result["failed"]}
    flow=result["flow"]; m=flow["metrics"]; row={"label":result["label"],"kind":result["kind"],"grid":result["grid"],"run_dir":flow["dir"],"q_supply_Lmin":result["q_supply_m3s"]*6e4,"eps":flow["eps"],"psi_deg":flow["psi"],"W_N":m["W"],"force_angle_deg":m["force_angle_deg"],"p_max_MPa_abs":m["p_max"]/1e6,"coupling":m.get("coupling"),"status":"PASS"}
    c=m.get("coupling") or {}
    if c.get("R_p",math.inf)>=1e-6 or c.get("R_T",math.inf)>=1e-6: row["status"]="FAILED_COUPLING"; return row
    pressure=S.midplane_profile(CASE_DIR/flow["dir"],"p")
    if pressure:
        theta,p=pressure; i=max(range(len(p)),key=p.__getitem__)
        row["pressure_peak_gauge_MPa"]=(p[i]-101325.0)/1e6
        row["pressure_peak_location_deg"]=(theta[i]-90)%360
    profile=S.midplane_profile(CASE_DIR/flow["dir"],"T_bearing")
    if profile:
        theta,t=profile; errors=[interp(theta,t,90+x)-y for x,y in ref]; i=max(range(len(t)),key=t.__getitem__)
        row.update(mean_error_K=sum(errors)/len(errors), rms_error_K=math.sqrt(sum(e*e for e in errors)/len(errors)), max_abs_error_K=max(map(abs,errors)), peak_temperature_C=t[i]-273.15, peak_location_deg=(theta[i]-90)%360)
    shaft=S.midplane_profile(CASE_DIR/flow["dir"],"T_journal")
    if shaft:
        row["shaft_temperature_C"]=sum(shaft[1])/len(shaft[1])-273.15
    return row

def grid_checks(a,b):
    checks={"eps_rel_pct":pct(a["eps"],b["eps"]),"psi_deg":circular(a["psi_deg"],b["psi_deg"]),"p_max_rel_pct":pct(a["p_max_MPa_abs"],b["p_max_MPa_abs"])}
    for key in ("mean_error_K","rms_error_K"):
        if key in a and key in b: checks[key]=abs(a[key]-b[key])
    if "peak_temperature_C" in a and "peak_temperature_C" in b:
        checks["peak_temperature_K"]=abs(a["peak_temperature_C"]-b["peak_temperature_C"]); checks["peak_location_deg"]=circular(a["peak_location_deg"],b["peak_location_deg"])
    return {k:{"value":v,"limit":LIMITS[k],"pass":v<=LIMITS[k]} for k,v in checks.items()}

def main():
    p=CASE_DIR/"campaign_equilibrium.json"
    if not p.exists(): raise SystemExit("no production campaign_equilibrium.json")
    ref=reference(); p_ref=pressure_reference(); data=json.loads(p.read_text()); runs={name:score(name,r,ref) for name,r in data.items()}
    exp_p_x,exp_p=max(p_ref,key=lambda item:item[1])
    for row in runs.values():
        if "pressure_peak_gauge_MPa" in row:
            row["pressure_peak_error_pct"]=100.0*(row["pressure_peak_gauge_MPa"]*1e6-exp_p)/exp_p
            row["pressure_peak_location_error_deg"]=circular(row["pressure_peak_location_deg"],exp_p_x)
    grids={}
    for label in ("r0","r1","r2"):
        a,b=runs.get(label+"_g6"),runs.get(label+"_g7")
        if a and b and a["status"]==b["status"]=="PASS":
            checks=grid_checks(a,b); grids[label]={"checks":checks,"pass":all(v["pass"] for v in checks.values())}
    fidelity={
        "outer_radius_fixed_temperature": "R1/R2 only",
        "bush_axial_hb_80": "NOT_MODELLED",
        "shaft_r_z_with_43C_ends": "NOT_MODELLED; R2 is a uniform scalar surrogate",
        "fresh_supply_Qi": "UNRESOLVED; FLOW_WEIGHTED uses a model-derived groove source",
        "temperature_validation_eligible": False,
    }
    out={"limits":LIMITS,"runs":runs,"grid_verdicts":grids,"ferron_fidelity":fidelity,"archive_comparison":json.loads((CASE_DIR/"baseline.json").read_text())}
    (CASE_DIR/"campaign_analysis.json").write_text(json.dumps(out,indent=2))
    print("case status q_Lmin eps psi_deg mean_K rms_K peak_C peak_deg")
    for name,row in runs.items(): print(name,row["status"],f"{row.get('q_supply_Lmin',float('nan')):.4f}",f"{row.get('eps',float('nan')):.6f}",f"{row.get('psi_deg',float('nan')):.3f}",f"{row.get('mean_error_K',float('nan')):+.3f}",f"{row.get('rms_error_K',float('nan')):.3f}",f"{row.get('peak_temperature_C',float('nan')):.3f}",f"{row.get('peak_location_deg',float('nan')):.3f}")
    for name,row in grids.items(): print(f"{name} GRID {'PASS' if row['pass'] else 'FAIL'}")

if __name__ == "__main__": main()
