import sys, pathlib
sys.path.insert(0, '/home/user/cases_run/01_ferron_conjugate')
import solver_io as S
d = pathlib.Path('/home/user/diag_h10')
m = S.frame_metrics(d)
th, tb = S.midplane_profile(d, 'T_bearing')
x = [(t-90.0) % 360.0 for t in th]
o = sorted(range(len(x)), key=lambda i: x[i])
xs = [x[i] for i in o]; vs = [tb[i]-273.15 for i in o]
keep = [(a,b) for a,b in zip(xs,vs) if 12 <= a <= 348]
pk = max(keep, key=lambda r: r[1])
print("W        = %.2f N"%m['W'])
print("T_max    = %.2f C"%(m['T_max']-273.15))
print("Tb_peak  = %.2f C @ %.1f deg   (min %.2f C)"%(pk[1], pk[0], min(b for _,b in keep)))
