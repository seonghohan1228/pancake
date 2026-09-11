"""Audit archived solves, force closure, mass balance and selected grid repeats."""
import json
import re
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def main():
    """Fail on invalid completed scan points; keep numerical sensitivity explicit."""
    rows=[json.loads(p.read_text()) for p in ROOT.glob("scan_n*_w*.json")]
    checks=[]
    for row in rows:
        run=ROOT/row["dir"]
        log=(run/"run.log").read_text()
        residual=re.findall(r"fbns iter=\d+ \|G\|=([\deE+.\-]+) \|FB\|=([\deE+.\-]+)",log)
        g,fb=map(float,residual[-1])
        diag=row["diagnostics"]
        source=float(diag["inlet_mass_source"])
        boundary=float(diag["liquid_boundary_flux"])
        mass_error=abs(source+boundary)/max(abs(source),abs(boundary),1e-30)
        a=np.load(run/"fields.npz")
        pressure=a["p"]
        fraction=a["film_content"]
        physical=bool(np.isfinite(pressure).all() and np.isfinite(fraction).all()
                      and pressure.min()>=row["p_cav"]-1 and fraction.min()>=0
                      and fraction.max()<=1.000001 and np.min(a["h"])>0)
        checks.append(dict(load=row["load_target"],rpm=row["rpm"],G=g,FB=fb,
                           force_error_rel=row["force_error_rel"],mass_error_rel=mass_error,
                           physical=physical,passed=bool(g<1e-7 and fb<1e-7 and
                           row["force_error_rel"]<.002 and mass_error<1e-5 and physical)))
    fine=[]
    for path in ROOT.glob("fine_n*_w*.json"):
        b=json.loads(path.read_text())
        a=json.loads((ROOT/path.name.replace("fine_","scan_")).read_text())
        fine.append(dict(load=b["load_target"],rpm=b["rpm"],
            arc_coarse=a["mid_arc_deg"],arc_fine=b["mid_arc_deg"],
            area_coarse=a["cav_area_pct"],area_fine=b["cav_area_pct"],
            h_min_coarse=a["h_min_exact_um"],h_min_fine=b["h_min_exact_um"],
            p_max_coarse=a["p_max_bar"],p_max_fine=b["p_max_bar"],
            force_error_fine=b["force_error_rel"]))
    result=dict(points=len(checks),all_passed=all(r["passed"] for r in checks),
                checks=checks,grid_comparison=fine)
    (ROOT/"validation.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    if not result["all_passed"]:
        raise SystemExit("At least one result fails validation")


if __name__=="__main__":
    main()
