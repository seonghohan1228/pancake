"""Export local legacy GCT pressure and nominal-THD temperature planes."""
import json
import sys
from pathlib import Path
import numpy as np
CASE_DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(CASE_DIR.parent/"01_ferron_wallbc_fields"))
import solver_io as S
def plane(run_dir,name):
    a=S._read_pvts(S.final_pvts(run_dir)); theta=(np.degrees(np.asarray(a["theta_rad"][0]))-90)%360; z=np.asarray(a["z_m"][0]); v=np.asarray(a[name][0]); xs=np.unique(theta); zs=np.unique(z); out=np.full((len(zs),len(xs)),np.nan); out[np.searchsorted(zs,z),np.searchsorted(xs,theta)]=v
    if np.isnan(out).any(): raise RuntimeError("incomplete field plane")
    return xs,zs/.08-.5,out
def main():
    runs=json.loads((CASE_DIR/"legacy_gct_runs.json").read_text(encoding="utf-8")); out={}
    for tag,field in (("legacy_gct_g5","p"),("legacy_gct_thd_g4","p"),("legacy_gct_thd_g4","T")):
        if tag not in runs: continue
        run=CASE_DIR/runs[tag]["run_dir"]
        if not S.final_pvts(run): continue
        theta,z,v=plane(run,field); key=f"{tag}_{field}"; out[key+"_theta_deg"]=theta; out[key+"_z_over_L"]=z; out[key+"_plane"]=v-101325 if field=="p" else v-273.15
    np.savez_compressed(CASE_DIR/"legacy_gct_visual_data.npz",**out)
if __name__=="__main__": main()
