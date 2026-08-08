"""Sequential REPAIRED_LEGACY reproduction of Claude's Case-00a GCT."""
import argparse
import json
import sys
from pathlib import Path

CASE_DIR=Path(__file__).resolve().parent
SOURCE=CASE_DIR.parent/"00_grid_dt_independence"
sys.path.insert(0,str(CASE_DIR.parent/"01_ferron_wallbc_fields"))
import solver_io as S

MATRIX={"legacy_gct_g1":("grid_G1",(60,16)),"legacy_gct_g2":("grid_G2",(120,32)),"legacy_gct_g3":("grid_G3",(240,64)),"legacy_gct_g4":("grid_G4",(360,96)),"legacy_gct_g5":("grid_G5",(480,128)),"legacy_gct_thd_g3":("grid_G3_thd",(240,64)),"legacy_gct_thd_g4":("grid_G4_thd",(360,96))}

def repaired_text(source_name):
    text=(SOURCE/source_name/"config.txt").read_text(encoding="utf-8")
    if source_name.endswith("_thd"):
        text=text.replace("dissolved_gas_henry_coeff = 0", "dissolved_gas_henry_coeff = 1e-12")
        text=text.replace("solution_viscosity_gas_coeff = 0", "solution_viscosity_gas_coeff = 1e-12")
    return text

def main():
    p=argparse.ArgumentParser(); p.add_argument("mode",nargs="?",choices=("ISOTHERMAL","THD","ALL"),default="ALL"); p.add_argument("--smoke",action="store_true"); args=p.parse_args()
    selected=[name for name in MATRIX if args.mode=="ALL" or (args.mode=="THD") == ("_thd_" in name)]
    if args.smoke:
        selected=[name for name in selected if name in ("legacy_gct_g1","legacy_gct_thd_g3")]
    output=CASE_DIR/("smoke_legacy_gct_runs.json" if args.smoke else "legacy_gct_runs.json"); runs=json.loads(output.read_text()) if output.exists() else {}
    run_root=CASE_DIR/"runs"/("smoke" if args.smoke else "production")
    for name in selected:
        if (runs.get(name, {}).get("rc") == 0 and runs[name].get("metrics")
                and "coupling_converged" in runs[name]["metrics"]):
            print(f"{name}: already archived, skipped"); continue
        runs.pop(name, None)
        source,grid=MATRIX[name]; work=run_root/name; work.mkdir(parents=True,exist_ok=True)
        (work/"config.txt").write_text(repaired_text(source),encoding="utf-8")
        rc,wall=S.run_case(work); metrics=S.frame_metrics(work) if rc==0 else None
        if metrics is not None:
            diagnostics=S.diagnostics_tail(work)
            metrics["flow_converged"]=diagnostics.get("converged")=="1"
            metrics["coupling_converged"]=diagnostics.get("coupling_converged")=="1"
            metrics["coupling_R_p"]=float(diagnostics.get("coupling_R_p","nan"))
            metrics["coupling_R_T"]=float(diagnostics.get("coupling_R_T","nan"))
        runs[name]={"source_config":str(SOURCE/source/"config.txt"),"run_dir":work.relative_to(CASE_DIR).as_posix(),"grid":grid,"classification":"REPAIRED_LEGACY" if source.endswith("_thd") else "LEGACY_EXACT","rc":rc,"wall_s":round(wall,1),"metrics":metrics}
        output.write_text(json.dumps(runs,indent=2),encoding="utf-8")
        if rc!=0 or metrics is None: raise RuntimeError(f"{name}: failed")

if __name__=="__main__": main()
