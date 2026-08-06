"""Score a conjugate rerun T_bearing profile against the 16 digitised Fig. 4 points."""
import sys, csv, math
from pathlib import Path
sys.path.insert(0, str(Path.home()/"pancake_picard/cases/01_ferron_wallbc_fields"))
import solver_io as S

REF = Path.home()/"pancake_picard/cases/ref_data"
cands = sorted(REF.glob("*temperature*.csv"))
print("reference candidates:", [c.name for c in cands])
if not cands: sys.exit(0)
rows=[]
with open(cands[0]) as fh:
    for line in fh:
        line=line.strip()
        if not line or line.startswith("#"): continue
        p=line.replace(",", " ").split()
        try: rows.append((float(p[0]), float(p[1])))
        except (ValueError, IndexError): pass
print("n ref points:", len(rows))

def score(d, label):
    th, v = S.midplane_profile(d, "T_bearing")
    # solver_theta = (90 + x) mod 360
    sol = sorted(((t, x) for t, x in zip(th, v)))
    def interp(xdeg):
        st = (90.0 + xdeg) % 360.0
        for i in range(len(sol)-1):
            if sol[i][0] <= st <= sol[i+1][0]:
                a,b = sol[i], sol[i+1]
                w = (st-a[0])/(b[0]-a[0]) if b[0]!=a[0] else 0.0
                return a[1] + w*(b[1]-a[1])
        return sol[-1][1]
    errs=[interp(x)-273.15-y for x,y in rows]
    n=len(errs)
    mean=sum(errs)/n
    rms=math.sqrt(sum(e*e for e in errs)/n)
    pts=[(t,x) for t,x in zip(th,v) if abs(((t-90+180)%360)-180) > 12]
    tpk,vpk=max(pts,key=lambda a:a[1])
    print("%-42s n=%2d mean=%+.3f K rms=%.3f K maxabs=%.3f K   peak=%.3f C @ %.2f deg  (exp 49.04 @ 258.2 -> %+.2f K)" %
          (label, n, mean, rms, max(abs(e) for e in errs), vpk-273.15, (tpk-90.0)%360.0, vpk-273.15-49.04))

score(Path.home()/"picard_rerun/adaptive/01_ferron_conjugate/outer_robin80/fb_iter_10", "outer_robin80 h=80 ADAPTIVE converged")
score(Path.home()/"picard_rerun/frozen/01_ferron_conjugate/outer_fixed313/fb_iter_10", "outer_fixed313 (control, matches archive)")
score(Path.home()/"picard_rerun/frozen/01_ferron_wallbc_fields/rerun_cht_robin10/fb_iter_07", "rerun_cht_robin10 h=10 (best iterate)")
score(Path.home()/"picard_rerun/adaptive/01_ferron_wallbc_fields/rerun_cht_robin10/fb_iter_07", "rerun_cht_robin10 h=10 (ADAPTIVE, converged)")
score(Path.home()/"pancake/cases/01_ferron_wallbc_fields/rerun_cht_robin10/fb_iter_07", "rerun_cht_robin10 h=10 ARCHIVE (iterate 200)")
