"""Sommerfeld-collapse discriminator: Ferron 1983, Brito 2006 and pancake on one
S-eps axis.

Question. The load-capacity comparison is an open FAIL (docs/UNFINISHED_BUSINESS.md
sec 2.1). Both experiments report eccentricity as a proximity-probe displacement
divided by an ASSUMED radial clearance, so the clearance basis enters the
"measurement" itself. Does a single clearance rescale per rig collapse experiment
onto solver (=> instrumentation), or does a load- or speed-dependent residual
survive every clearance choice (=> physics gap)?

Method (the clearance-free invariant).
    S    = mu N L D (R/c)^2 / W          (N in rev/s, D = 2R)
    eps  = e / c                          (e = probe displacement, c = assumed clearance)
    =>  Lambda := eps^2 / S = e^2 W / (mu N L D R^2)
Lambda contains no clearance at all. It is therefore fixed by what the rig
actually measured (e, W) and by nominal fluid/geometry constants, and is
invariant under any rescale of the assumed clearance. Under c -> k c a data point
slides along the parabola S ∝ eps^2 in the (S, eps) plane; Lambda labels the
parabola.

Consequently each measured point has an EXACT per-point effective clearance, no
fitting required: invert the solver's own Lambda_solver(eps) = eps^2/S_solver(eps)
(monotone increasing) at Lambda_measured to get eps*, then

    c_eff = e / eps* = c_assumed * (eps_measured / eps*).

A single clearance rescale explains the dataset if and only if c_eff is CONSTANT
across the dataset's independent variable (load for Ferron, speed for Brito).
The scatter of c_eff is the discriminator; its level is not, because a uniform
multiplicative solver load error is exactly degenerate with a uniform clearance
error (see COLLAPSE_ANALYSIS.md sec 6).

Note mu, N, L, D, R cancel between Lambda_measured and Lambda_solver whenever the
solver run uses the same nominal values as the experiment, which it does here.
c_eff/c_assumed is therefore insensitive to the viscosity convention.

Run:  python collapse_analysis.py
Out:  s_eps_collapse.png, collapse_analysis.json (stdout carries every table)
"""
import json
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from plot_style import apply_style, save_figure, clean_log_ticks, MARKERS, COLORS

# --------------------------------------------------------------------------
# Rig constants
# --------------------------------------------------------------------------
R = 0.05
D = 2.0 * R
L = 0.08                      # both rigs: L/D = 0.8

FERRON = dict(
    mu=0.0277,                # Pa s at 40 C, the supply temperature
    rpm=2000.0,
    c_digitised=145e-6,       # basis of cases/ref_data/ferron1983_fig5_*.csv
    c_solver=152e-6,          # basis of the case 01 THD sweep
    L=0.08,
)
BRITO = dict(
    mu=0.0293,                # Pa s at 40 C = Tf
    c_quoted=85.5e-6,         # 171 um diametral at 20 C
    c_hot=93.024e-6,          # +8.8% thermal-growth bracket
    L=0.08,
)


def sommerfeld(mu, rpm, c, W, Lb=L):
    """S = mu N L D (R/c)^2 / W, N in rev/s."""
    return mu * (rpm / 60.0) * Lb * D * (R / c) ** 2 / W


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
# Ferron Fig. 5 experiment, C = 145 um basis (cases/ref_data/...csv)
FERRON_EXP = np.array([
    (2018.3076478185085, 0.49206769529848343),
    (3006.0606611162334, 0.5967814755431203),
    (3966.4854373047838, 0.6610061636151398),
    (5025.2540163424355, 0.716115559835213),
    (6058.019134308065, 0.7513966924505936),
    (7016.536233003122, 0.7858697873848016),
    (8015.795872452074, 0.8128947643336984),
    (9047.395187505012, 0.8299943678093049),
    (10062.718682672989, 0.8504037917361227),
])
FERRON_EXP_SIGMA = 0.015      # cluster scatter (+-0.01) plus digitising error

# pancake case 01 THD Fig. 5 sweep (L/D = 0.8, 2000 rpm, c = 152 um, J-adiabatic,
# FLOW_WEIGHTED groove). cases/01_ferron_thd/SUMMARY.md Part B.
FERRON_SOLVER = np.array([
    (2000.0, 0.4061),
    (3000.0, 0.5051),
    (4000.0, 0.5738),
    (6000.0, 0.6636),
    (8000.0, 0.7207),
    (10000.0, 0.7607),
])

# Brito Fig. 6, W = 10 kN, Tf = 40 C, c basis 85.5 um.
# Two independent digitisations (validation/journal_bearing/brito.py).
BRITO_EXP = {   # rpm: (eps_polar, eps_independent) at Pf = 70 kPa gauge
    1000.0: (0.761, 0.7704),
    3000.0: (0.576, 0.5769),
    4000.0: (0.558, 0.5447),
}
BRITO_EXP_210 = {   # Pf = 210 kPa gauge; no matching solver run
    1000.0: (0.754, 0.7533),
    3000.0: (0.563, 0.5658),
    4000.0: (0.526, 0.5264),
}
# pancake eccentricity carrying the applied 10 kN (validation/journal_bearing/
# brito_eps_for_load.json, Pf = 70 kPa, ENERGY_EQUATION).
BRITO_SOLVER = {    # rpm: (eps at c = 85.5 um, eps at c = 93.02 um, converged?)
    1000.0: (0.6483792075937971, 0.6782876440997776, False),
    3000.0: (0.42600428203299634, 0.46070792466203914, True),
    4000.0: (0.36654964476767854, 0.4009827081812275, True),
}
BRITO_W = 10000.0

# pancake case 02 isothermal plain-bearing loci (this directory).
LOCUS_LD1 = np.array([
    (0.1, 1.3159536988125662), (0.2, 0.6277873928379692), (0.3, 0.3869194498963979),
    (0.4, 0.25935978192631015), (0.5, 0.17802370857761493), (0.6, 0.12060378210282854),
    (0.7, 0.07762911087983446), (0.8, 0.04448089506183887), (0.9, 0.018797544075144074),
])
LOCUS_LD05 = np.array([
    (0.1, 4.201489872376309), (0.2, 2.000959102985472), (0.3, 1.2075818529119708),
    (0.4, 0.7798334854213217), (0.5, 0.5066604999874068), (0.6, 0.3180904548125792),
    (0.7, 0.18456340503684138), (0.8, 0.09154870567902457), (0.9, 0.031205215450644787),
])
# Raimondi & Boyd 1958 (validation/journal_bearing/raimondi_boyd.py), eps -> S.
RB_LD1 = np.array([(0.1, 1.33), (0.2, 0.631), (0.4, 0.264), (0.6, 0.121),
                   (0.8, 0.0446), (0.9, 0.0188)])
RB_LD05 = np.array([(0.1, 4.31), (0.2, 2.03), (0.4, 0.779), (0.6, 0.319),
                    (0.8, 0.0923), (0.9, 0.0313)])

# L/D interpolation weight, linear in ln(L/D) between the 1 and 0.5 blocks.
W_LD08 = math.log(0.8) / math.log(0.5)


def interp_ld08(a1, a05):
    """Interpolate an (eps, S) locus to L/D = 0.8, linear in ln S and ln(L/D),
    on the common eps grid of the two inputs."""
    eps = a1[:, 0]
    assert np.allclose(eps, a05[:, 0])
    lnS = np.log(a1[:, 1]) + W_LD08 * (np.log(a05[:, 1]) - np.log(a1[:, 1]))
    return np.column_stack([eps, np.exp(lnS)])


# --------------------------------------------------------------------------
# The clearance-free inversion
# --------------------------------------------------------------------------
class SolverCurve:
    """Lambda_solver(eps) = eps^2 / S_solver(eps), and its inverse.

    Stored as ln(Lambda) against eps, which is very close to linear for every
    curve here, so linear interpolation and a short linear extrapolation off
    either end are both well behaved."""

    def __init__(self, eps, S, label):
        order = np.argsort(eps)
        self.eps = np.asarray(eps, float)[order]
        self.S = np.asarray(S, float)[order]
        self.lnlam = np.log(self.eps ** 2 / self.S)
        self.label = label
        if np.any(np.diff(self.lnlam) <= 0):
            raise ValueError(f"{label}: Lambda(eps) not monotone")

    def lam(self, eps):
        return np.exp(np.interp(eps, self.eps, self.lnlam))

    def eps_of_lambda(self, lam):
        """Invert, with linear extrapolation in (eps, ln Lambda) off either end.
        Returns (eps*, extrapolated_flag)."""
        y = math.log(lam)
        if y < self.lnlam[0]:
            s = (self.lnlam[1] - self.lnlam[0]) / (self.eps[1] - self.eps[0])
            return self.eps[0] + (y - self.lnlam[0]) / s, True
        if y > self.lnlam[-1]:
            s = (self.lnlam[-1] - self.lnlam[-2]) / (self.eps[-1] - self.eps[-2])
            return self.eps[-1] + (y - self.lnlam[-1]) / s, True
        return float(np.interp(y, self.lnlam, self.eps)), False

    def S_of_eps(self, eps):
        """Solver Sommerfeld number at eps, extrapolating through ln Lambda."""
        eps = np.atleast_1d(np.asarray(eps, float))
        lo = (self.lnlam[1] - self.lnlam[0]) / (self.eps[1] - self.eps[0])
        hi = (self.lnlam[-1] - self.lnlam[-2]) / (self.eps[-1] - self.eps[-2])
        y = np.interp(eps, self.eps, self.lnlam)
        y = np.where(eps < self.eps[0], self.lnlam[0] + lo * (eps - self.eps[0]), y)
        y = np.where(eps > self.eps[-1], self.lnlam[-1] + hi * (eps - self.eps[-1]), y)
        out = eps ** 2 / np.exp(y)
        return out if out.size > 1 else float(out[0])

    def eps_of_S(self, S):
        """Invert S_of_eps by bisection on the monotone branch (S decreasing)."""
        lo, hi = 0.02, 0.995
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if self.S_of_eps(mid) > S:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def c_eff_point(curve, eps_meas, c_assumed, S_assumed, sigma_eps):
    """Per-point effective clearance [m] and its 1-sigma [m] from the eps reading
    error. S_assumed is fixed by (c_assumed, W, speed) and carries no eps error."""
    lam = eps_meas ** 2 / S_assumed
    eps_star, extrap = curve.eps_of_lambda(lam)
    c = c_assumed * eps_meas / eps_star
    out = []
    for d in (-sigma_eps, sigma_eps):
        em = eps_meas + d
        es, _ = curve.eps_of_lambda(em ** 2 / S_assumed)
        out.append(c_assumed * em / es)
    return c, eps_star, 0.5 * abs(out[1] - out[0]), extrap


# --------------------------------------------------------------------------
# Build the curves
# --------------------------------------------------------------------------
apply_style()

ferron_curve = SolverCurve(
    FERRON_SOLVER[:, 1],
    sommerfeld(FERRON["mu"], FERRON["rpm"], FERRON["c_solver"], FERRON_SOLVER[:, 0]),
    "pancake THD, L/D=0.8 (Ferron rig)")

brito_rpm = np.array(sorted(BRITO_SOLVER))
brito_eps_solv = np.array([BRITO_SOLVER[r][0] for r in brito_rpm])
brito_S_solv = np.array([sommerfeld(BRITO["mu"], r, BRITO["c_quoted"], BRITO_W)
                         for r in brito_rpm])
brito_curve = SolverCurve(brito_eps_solv, brito_S_solv, "pancake THD (Brito rig)")

pan_ld08 = interp_ld08(LOCUS_LD1, LOCUS_LD05)
pan08_curve = SolverCurve(pan_ld08[:, 0], pan_ld08[:, 1], "pancake isothermal L/D=0.8")
rb_ld08 = interp_ld08(RB_LD1, RB_LD05)
rb08_curve = SolverCurve(rb_ld08[:, 0], rb_ld08[:, 1], "Raimondi & Boyd L/D=0.8")

results = {"method": "Lambda = eps^2/S is clearance-free; c_eff = c_assumed * eps_meas/eps*"}

# --- bonus: score case 02 directly against the transcribed R&B table ---------
# docs/UNFINISHED_BUSINESS.md sec 2.1 lists this transcription as the pending
# first step; it is already in validation/journal_bearing/raimondi_boyd.py.
PAN_PHI = {"1": {0.1: 79.43101915804715, 0.2: 73.8575686908869, 0.4: 62.555866486586524,
                 0.6: 50.50445734073219, 0.8: 36.189850762585195, 0.9: 26.432211381269497},
           "0.5": {0.1: 81.92645793220572, 0.2: 75.2802725639686, 0.4: 61.87328205982004,
                   0.6: 48.081056295304165, 0.8: 32.913023961453916, 0.9: 23.52624356375344}}
RB_PHI = {"1": {0.1: 79.5, 0.2: 74.02, 0.4: 63.10, 0.6: 50.58, 0.8: 36.24, 0.9: 26.45},
          "0.5": {0.1: 81.62, 0.2: 74.94, 0.4: 61.45, 0.6: 48.14, 0.8: 33.31, 0.9: 23.66}}
print("\n=== case 02 locus vs the transcribed Raimondi & Boyd table ===")
rb_score = {}
for key, pan, rb in (("1", LOCUS_LD1, RB_LD1), ("0.5", LOCUS_LD05, RB_LD05)):
    print(f"\n  L/D = {key}:  {'eps':>5} {'S (pancake)':>12} {'S (R&B)':>9} {'dS %':>7} "
          f"{'phi (pan)':>10} {'phi (R&B)':>10} {'d phi':>7}")
    rows = []
    for e, S_rb in rb:
        S_pan = float(pan[np.isclose(pan[:, 0], e), 1][0])
        dS = 100 * (S_pan - S_rb) / S_rb
        dphi = PAN_PHI[key][e] - RB_PHI[key][e]
        rows.append(dict(eps=e, S_pancake=S_pan, S_rb=S_rb, dS_pct=dS, dphi_deg=dphi))
        print(f"        {e:5.1f} {S_pan:12.5f} {S_rb:9.4f} {dS:+7.2f} "
              f"{PAN_PHI[key][e]:10.2f} {RB_PHI[key][e]:10.2f} {dphi:+7.2f}")
    rb_score[key] = rows
    print(f"        max |dS| = {max(abs(r['dS_pct']) for r in rows):.2f} %, "
          f"max |d phi| = {max(abs(r['dphi_deg']) for r in rows):.2f} deg")
results["case02_vs_raimondi_boyd"] = rb_score

# --- similarity audit: does the solver's own eps(S) survive a clearance change?
sim = []
for r in brito_rpm:
    e20, ehot, _ = BRITO_SOLVER[r]
    S_hot = sommerfeld(BRITO["mu"], r, BRITO["c_hot"], BRITO_W)
    eps_pred, _ = brito_curve.eps_of_lambda(ehot ** 2 / S_hot)
    sim.append((r, ehot, eps_pred, ehot - eps_pred))
results["similarity_audit"] = [
    dict(rpm=r, eps_hot_run=a, eps_from_20C_curve=b, d_eps=d) for r, a, b, d in sim]

print("\n=== similarity audit: 20 C curve predicting the +8.8% hot-clearance runs ===")
print(f"{'rpm':>6} {'eps hot run':>12} {'from 20C curve':>15} {'d eps':>8}")
for r, a, b, d in sim:
    print(f"{r:6.0f} {a:12.4f} {b:15.4f} {d:+8.4f}")

# --------------------------------------------------------------------------
# Ferron
# --------------------------------------------------------------------------
def ferron_table(curve, label):
    rows = []
    for Wexp, eps145 in FERRON_EXP:
        S145 = sommerfeld(FERRON["mu"], FERRON["rpm"], FERRON["c_digitised"], Wexp)
        c, eps_star, sig, ex = c_eff_point(curve, eps145, FERRON["c_digitised"],
                                           S145, FERRON_EXP_SIGMA)
        rows.append(dict(W=Wexp, eps145=eps145, Lambda=eps145 ** 2 / S145,
                         eps_star=eps_star, c_eff_um=c * 1e6, sigma_um=sig * 1e6,
                         extrapolated=ex))
    return rows


print("\n=== FERRON 1983 Fig. 5, 2000 rpm: per-point effective clearance ===")
for label, curve in (("pancake THD L/D=0.8 (case 01)", ferron_curve),
                     ("pancake isothermal L/D=0.8 (case 02 interp)", pan08_curve),
                     ("Raimondi & Boyd L/D=0.8", rb08_curve)):
    rows = ferron_table(curve, label)
    ce = np.array([r["c_eff_um"] for r in rows])
    results.setdefault("ferron", {})[label] = dict(
        rows=rows, mean_c_eff_um=float(ce.mean()), sd_c_eff_um=float(ce.std(ddof=1)),
        spread_pct=float(100 * (ce.max() - ce.min()) / ce.mean()))
    print(f"\n  reference: {label}")
    print(f"  {'W [N]':>8} {'eps(145)':>9} {'Lambda':>8} {'eps*':>7} "
          f"{'c_eff [um]':>11} {'+-':>6} {'extrap':>7}")
    for r in rows:
        print(f"  {r['W']:8.0f} {r['eps145']:9.4f} {r['Lambda']:8.4f} {r['eps_star']:7.4f} "
              f"{r['c_eff_um']:11.1f} {r['sigma_um']:6.1f} {str(r['extrapolated']):>7}")
    p = np.polyfit(np.log([r["W"] for r in rows]), ce, 1)
    print(f"  mean c_eff = {ce.mean():.1f} um, sd = {ce.std(ddof=1):.1f} um, "
          f"peak-to-peak {100*(ce.max()-ce.min())/ce.mean():.1f} %, "
          f"trend d c_eff/d ln W = {p[0]:+.2f} um/decade-e")
    results["ferron"][label]["trend_um_per_lnW"] = float(p[0])

# --------------------------------------------------------------------------
# Brito
# --------------------------------------------------------------------------
def brito_table(curve, exp_dict, tag):
    rows = []
    for rpm in sorted(exp_dict):
        e_pair = exp_dict[rpm]
        eps_m = float(np.mean(e_pair))
        sig = max(0.5 * abs(e_pair[0] - e_pair[1]), 0.02)   # spread or stated 0.02
        Sq = sommerfeld(BRITO["mu"], rpm, BRITO["c_quoted"], BRITO_W)
        c, eps_star, sc, ex = c_eff_point(curve, eps_m, BRITO["c_quoted"], Sq, sig)
        rows.append(dict(rpm=rpm, eps_meas=eps_m, sigma_eps=sig, Lambda=eps_m ** 2 / Sq,
                         eps_star=eps_star, c_eff_um=c * 1e6, sigma_um=sc * 1e6,
                         extrapolated=ex, tag=tag))
    return rows


print("\n=== BRITO 2006 Fig. 6, W = 10 kN: per-point effective clearance ===")
for label, curve in (("pancake THD (Brito rig, Pf=70 kPa)", brito_curve),
                     ("pancake isothermal L/D=0.8 (case 02 interp)", pan08_curve),
                     ("Raimondi & Boyd L/D=0.8", rb08_curve)):
    rows = brito_table(curve, BRITO_EXP, "70kPa")
    ce = np.array([r["c_eff_um"] for r in rows])
    results.setdefault("brito", {})[label] = dict(
        rows=rows, mean_c_eff_um=float(ce.mean()), sd_c_eff_um=float(ce.std(ddof=1)),
        spread_pct=float(100 * (ce.max() - ce.min()) / ce.mean()))
    print(f"\n  reference: {label}")
    print(f"  {'rpm':>6} {'eps meas':>9} {'Lambda':>8} {'eps*':>7} "
          f"{'c_eff [um]':>11} {'+-':>6} {'extrap':>7}")
    for r in rows:
        print(f"  {r['rpm']:6.0f} {r['eps_meas']:9.4f} {r['Lambda']:8.4f} {r['eps_star']:7.4f} "
              f"{r['c_eff_um']:11.1f} {r['sigma_um']:6.1f} {str(r['extrapolated']):>7}")
    print(f"  mean c_eff = {ce.mean():.1f} um, sd = {ce.std(ddof=1):.1f} um, "
          f"peak-to-peak {100*(ce.max()-ce.min())/ce.mean():.1f} %, "
          f"1000->4000 rpm change {100*(ce[-1]-ce[0])/ce[0]:+.1f} %")
    results["brito"][label]["change_1000_to_4000_pct"] = float(100 * (ce[-1] - ce[0]) / ce[0])

rows210 = brito_table(brito_curve, BRITO_EXP_210, "210kPa")
results["brito_210kPa"] = rows210
print("\n  Pf = 210 kPa (no matching solver run; same 70 kPa reference curve)")
for r in rows210:
    print(f"  {r['rpm']:6.0f} {r['eps_meas']:9.4f} c_eff = {r['c_eff_um']:.1f} um")

# --------------------------------------------------------------------------
# Residual after imposing ONE clearance per rig
# --------------------------------------------------------------------------
print("\n=== residual in eps after collapsing each rig on its own single best c_eff ===")
resid = {}
for tag, curve, c_best, pts in (
        ("ferron", ferron_curve,
         np.mean([r["c_eff_um"] for r in results["ferron"]["pancake THD L/D=0.8 (case 01)"]["rows"]]) * 1e-6,
         [(W, e) for W, e in FERRON_EXP]),
        ("brito", brito_curve,
         np.mean([r["c_eff_um"] for r in results["brito"]["pancake THD (Brito rig, Pf=70 kPa)"]["rows"]]) * 1e-6,
         [(r, float(np.mean(BRITO_EXP[r]))) for r in sorted(BRITO_EXP)])):
    print(f"\n  {tag}: single c_eff = {c_best*1e6:.1f} um")
    rr = []
    for x, eps_assumed in pts:
        if tag == "ferron":
            c_as, rpm, Wv = FERRON["c_digitised"], FERRON["rpm"], x
        else:
            c_as, rpm, Wv = BRITO["c_quoted"], x, BRITO_W
        k = c_best / c_as
        eps_new = eps_assumed / k
        S_new = sommerfeld(FERRON["mu"] if tag == "ferron" else BRITO["mu"],
                           rpm, c_best, Wv)
        eps_solver_at_S = curve.eps_of_S(S_new)
        rr.append(dict(x=x, eps_rescaled=eps_new, S=S_new,
                       eps_solver=eps_solver_at_S, d_eps=eps_new - eps_solver_at_S))
        print(f"    x = {x:7.0f}  S = {S_new:7.4f}  eps(exp) = {eps_new:.4f}  "
              f"eps(solver) = {eps_solver_at_S:.4f}  d = {eps_new-eps_solver_at_S:+.4f}")
    resid[tag] = rr
results["residual_single_c"] = resid

# --------------------------------------------------------------------------
# What the clearance error is worth in load: solver load at the measured
# displacement, divided by the applied load, as a function of the clearance basis.
# --------------------------------------------------------------------------
print("\n=== load over-prediction W_solver(at measured e) / W_applied vs clearance basis ===")
lo = {}
print("\n  Ferron 2000 rpm")
hdr = [145, 152, 160, 171]
print("      W [N]  " + "".join(f"c={h:>4} um " for h in hdr))
for Wexp, eps145 in FERRON_EXP:
    e_disp = eps145 * FERRON["c_digitised"]
    row = []
    for cc in hdr:
        eps_here = e_disp / (cc * 1e-6)
        S_here = ferron_curve.S_of_eps(eps_here)
        W_solv = sommerfeld(FERRON["mu"], FERRON["rpm"], cc * 1e-6, 1.0) / S_here
        row.append(W_solv / Wexp)
    lo[f"ferron_{Wexp:.0f}"] = row
    print(f"  {Wexp:9.0f}  " + "".join(f"{v:10.2f} " for v in row))

print("\n  Brito 10 kN")
hdrb = [85.5, 93.0, 99.7, 104.3]
print("      rpm  " + "".join(f"c={h:>5} um " for h in hdrb))
for rpm in sorted(BRITO_EXP):
    eps_m = float(np.mean(BRITO_EXP[rpm]))
    e_disp = eps_m * BRITO["c_quoted"]
    row = []
    for cc in hdrb:
        eps_here = e_disp / (cc * 1e-6)
        S_here = brito_curve.S_of_eps(eps_here)
        W_solv = sommerfeld(BRITO["mu"], rpm, cc * 1e-6, 1.0) / S_here
        row.append(W_solv / BRITO_W)
    lo[f"brito_{rpm:.0f}"] = row
    print(f"  {rpm:8.0f}  " + "".join(f"{v:11.2f} " for v in row))
results["load_overprediction"] = dict(ferron_bases_um=hdr, brito_bases_um=hdrb, ratios=lo)

with open("collapse_analysis.json", "w") as fh:
    json.dump(results, fh, indent=1)

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4))
ax0, ax1, ax2 = axes

c_f = results["ferron"]["pancake THD L/D=0.8 (case 01)"]["mean_c_eff_um"] * 1e-6
c_b = results["brito"]["pancake THD (Brito rig, Pf=70 kPa)"]["mean_c_eff_um"] * 1e-6

eps_fine = np.linspace(0.30, 0.92, 200)


def draw_solver(ax, legend=False):
    ax.plot(pan08_curve.S_of_eps(eps_fine), eps_fine, "-", color="0.55", lw=1.4,
            label="pancake isothermal, $L/D=0.8$" if legend else None, zorder=2)
    ax.plot(rb08_curve.S_of_eps(eps_fine), eps_fine, ":", color="0.35", lw=1.4,
            label="Raimondi & Boyd, $L/D=0.8$" if legend else None, zorder=2)
    ef = np.linspace(ferron_curve.eps[0], ferron_curve.eps[-1], 120)
    ax.plot(ferron_curve.S_of_eps(ef), ef, "-", color=COLORS[0], lw=2.0,
            label="pancake THD, Ferron rig" if legend else None, zorder=3)
    eb = np.linspace(brito_curve.eps[0], brito_curve.eps[-1], 120)
    ax.plot(brito_curve.S_of_eps(eb), eb, "--", color=COLORS[1], lw=2.0,
            label="pancake THD, Brito rig" if legend else None, zorder=3)


# ---- panel (a): as published -------------------------------------------
draw_solver(ax0, legend=True)
Sf145 = sommerfeld(FERRON["mu"], FERRON["rpm"], FERRON["c_digitised"], FERRON_EXP[:, 0])
ax0.plot(Sf145, FERRON_EXP[:, 1], MARKERS[0], mfc="none", mec=COLORS[0], ms=7,
         ls="none", label=r"Ferron, $C=145\,\mu$m", zorder=4)
k171 = FERRON["c_digitised"] / 171e-6
ax0.plot(Sf145 * k171 ** 2, FERRON_EXP[:, 1] * k171, MARKERS[3], mfc=COLORS[0],
         mec=COLORS[0], ms=6, ls="none", alpha=0.55,
         label=r"Ferron, $C_{\rm eff}=171\,\mu$m", zorder=4)
Sb = np.array([sommerfeld(BRITO["mu"], r, BRITO["c_quoted"], BRITO_W) for r in brito_rpm])
eb_meas = np.array([np.mean(BRITO_EXP[r]) for r in brito_rpm])
eb_sig = np.array([max(0.5 * abs(BRITO_EXP[r][0] - BRITO_EXP[r][1]), 0.02) for r in brito_rpm])
ax0.errorbar(Sb, eb_meas, yerr=eb_sig, fmt=MARKERS[1], mfc="none", mec=COLORS[1],
             ecolor=COLORS[1], ms=8, ls="none", label=r"Brito, $c=85.5\,\mu$m", zorder=4)

# rescale trajectory: a clearance change slides a point along S ∝ eps^2
lam_guide = eb_meas[0] ** 2 / Sb[0]
eg = np.linspace(0.44, 0.94, 50)
ax0.plot(eg ** 2 / lam_guide, eg, "-", color="k", lw=1.0, alpha=0.75, zorder=1)
ax0.annotate("$c$ rescale:\n$S\\propto\\varepsilon^{2}$", xy=(0.155, 0.905),
             xycoords="data", fontsize=12, ha="center", va="center")

ax0.set_xscale("log")
ax0.set_xlim(0.035, 1.6)
ax0.set_ylim(0.28, 0.98)
ax0.set_xlabel(r"$S=\mu N L D (R/c)^{2}/W$")
ax0.set_ylabel(r"$\varepsilon$")
clean_log_ticks(ax0)
ax0.legend(loc="lower left", fontsize=9.5, handletextpad=0.4, labelspacing=0.22,
           borderpad=0.2, frameon=True, framealpha=0.9, edgecolor="none")
ax0.text(0.97, 0.94, "(a) as published", transform=ax0.transAxes, fontsize=14,
         ha="right", va="top")

# ---- panel (b): after the per-rig collapse -----------------------------
draw_solver(ax1)
kf = FERRON["c_digitised"] / c_f
ax1.plot(Sf145 * kf ** 2, FERRON_EXP[:, 1] * kf, MARKERS[0], mfc="none", mec=COLORS[0],
         ms=7, ls="none", label=rf"Ferron $\times\,c_{{\rm eff}}={c_f*1e6:.0f}\,\mu$m",
         zorder=4)
kb = BRITO["c_quoted"] / c_b
ax1.errorbar(Sb * kb ** 2, eb_meas * kb, yerr=eb_sig * kb, fmt=MARKERS[1], mfc="none",
             mec=COLORS[1], ecolor=COLORS[1], ms=8, ls="none",
             label=rf"Brito $\times\,c_{{\rm eff}}={c_b*1e6:.0f}\,\mu$m", zorder=4)
for r, x, y in zip(brito_rpm, Sb * kb ** 2, eb_meas * kb):
    ax1.annotate(f"{r/1000:.0f}k", xy=(x, y), xytext=(4, -13), textcoords="offset points",
                 fontsize=11, color=COLORS[1])
ax1.set_xscale("log")
ax1.set_xlim(0.035, 1.6)
ax1.set_ylim(0.28, 0.98)
ax1.set_xlabel(r"$S=\mu N L D (R/c_{\rm eff})^{2}/W$")
ax1.set_ylabel(r"$\varepsilon$")
clean_log_ticks(ax1)
ax1.legend(loc="lower left", fontsize=10.5, handletextpad=0.4, labelspacing=0.22,
           borderpad=0.2, frameon=True, framealpha=0.9, edgecolor="none")
ax1.text(0.97, 0.94, "(b) one $c_{\\rm eff}$ per rig", transform=ax1.transAxes,
         fontsize=14, ha="right", va="top")

# ---- panel (c): per-point c_eff ----------------------------------------
fr = results["ferron"]["pancake THD L/D=0.8 (case 01)"]["rows"]
br = results["brito"]["pancake THD (Brito rig, Pf=70 kPa)"]["rows"]
ax2.errorbar([r["W"] / 1e3 for r in fr], [r["c_eff_um"] for r in fr],
             yerr=[r["sigma_um"] for r in fr], fmt=MARKERS[0], mfc="none",
             mec=COLORS[0], ecolor=COLORS[0], ms=7, ls="none",
             label="Ferron (2 kN$\\to$10 kN, 2000 rpm)")
ax2.axhline(c_f * 1e6, color=COLORS[0], lw=1.2, ls="-", alpha=0.6)
ax2b = ax2.twiny()
ax2b.grid(False)
ax2b.errorbar([r["rpm"] / 1e3 for r in br], [r["c_eff_um"] for r in br],
              yerr=[r["sigma_um"] for r in br], fmt=MARKERS[1], mfc="none",
              mec=COLORS[1], ecolor=COLORS[1], ms=8, ls="none",
              label="Brito (10 kN, 1000$\\to$4000 rpm)")
ax2b.plot([r["rpm"] / 1e3 for r in br], [r["c_eff_um"] for r in br], "--",
          color=COLORS[1], lw=1.2, alpha=0.7)
ax2b.set_xlim(0.3, 4.7)
ax2b.set_xticks([1, 2, 3, 4])
ax2b.xaxis.set_minor_locator(mticker.NullLocator())
ax2b.set_xlabel(r"Brito: $N$ [krpm]", color=COLORS[1], fontsize=18, labelpad=6)
ax2b.tick_params(axis="x", colors=COLORS[1], direction="in", length=6, width=1.0)
for yv, tx, ls in ((145, "145", ":"), (152, "152", ":"), (171, "171", ":"),
                   (85.5, "85.5", "-.")):
    ax2.axhline(yv, color="0.65", lw=1.0, ls=ls, zorder=0)
    ax2.text(1.45, yv + 1.5, tx, fontsize=11, color="0.45", va="bottom")
ax2.set_xlim(1.2, 11.0)
ax2.set_ylim(72, 192)
ax2.set_xticks([2, 4, 6, 8, 10])
ax2.set_xlabel(r"Ferron: $W$ [kN]", color=COLORS[0], fontsize=18)
ax2.tick_params(axis="x", colors=COLORS[0])
ax2.set_ylabel(r"$c_{\rm eff}$ [$\mu$m]")
h0, l0 = ax2.get_legend_handles_labels()
h1, l1 = ax2b.get_legend_handles_labels()
ax2.legend(h0 + h1, l0 + l1, loc="upper right", fontsize=10.5, handletextpad=0.4,
           labelspacing=0.22, borderpad=0.2)
ax2.text(0.04, 0.05, "(c) per-point $c_{\\rm eff}$", transform=ax2.transAxes,
         fontsize=14)

fig.tight_layout(w_pad=1.6)
save_figure(fig, "s_eps_collapse")
print("\nwrote s_eps_collapse.png and collapse_analysis.json")
