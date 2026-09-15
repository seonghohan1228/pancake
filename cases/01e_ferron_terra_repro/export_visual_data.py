"""Export compact local profiles and fields from completed case-01e runs."""

import json
import sys
from pathlib import Path
import numpy as np

CASE_DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(CASE_DIR.parent/"01_ferron_wallbc_fields"))
import solver_io as S

def plane(run_dir, field):
    a=S._read_pvts(S.final_pvts(run_dir)); theta=(np.degrees(np.asarray(a["theta_rad"][0]))-90)%360; z=np.asarray(a["z_m"][0]); value=np.asarray(a[field][0]); xs=np.unique(theta); zs=np.unique(z); out=np.full((len(zs),len(xs)),np.nan); out[np.searchsorted(zs,z),np.searchsorted(xs,theta)]=value
    if np.isnan(out).any(): raise RuntimeError(f"incomplete {field} plane")
    return xs,zs,out

def main():
    eq=json.loads((CASE_DIR/"campaign_equilibrium.json").read_text()); out={}
    for name,r in eq.items():
        if "flow" not in r: continue
        d=CASE_DIR/r["flow"]["dir"]; t=S.midplane_profile(d,"T_bearing"); t_bar=S.midplane_profile(d,"T"); p=S.midplane_profile(d,"p")
        if t: out[name+"_theta_T_deg"]=(np.asarray(t[0])-90)%360; out[name+"_T_bearing_C"]=np.asarray(t[1])-273.15
        if t_bar: out[name+"_theta_Tbar_deg"]=(np.asarray(t_bar[0])-90)%360; out[name+"_T_bar_C"]=np.asarray(t_bar[1])-273.15
        if p: out[name+"_theta_p_deg"]=(np.asarray(p[0])-90)%360; out[name+"_p_gauge_MPa"]=(np.asarray(p[1])-101325)/1e6
    for name in ("r0_g6","r1_g6","r2_g6","d_grid_g6"):
        if name not in eq or "flow" not in eq[name]: continue
        d=CASE_DIR/eq[name]["flow"]["dir"]; theta,z,p=plane(d,"p"); _,_,T=plane(d,"T"); out[name+"_theta_deg"]=theta; out[name+"_z_over_L"]=z/.08-.5; out[name+"_plane_p_gauge_MPa"]=(p-101325)/1e6; out[name+"_plane_T_C"]=T-273.15
    np.savez_compressed(CASE_DIR/"visual_data.npz",**out)

if __name__=="__main__": main()
