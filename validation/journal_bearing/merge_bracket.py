"""Fold the fixed-eccentricity thermal bracket into `ferron_results.json`.

The three headline Ferron cases are each run at the eccentricity that carries the
measured 4 kN load, which is the right comparison for the pressure *profile* but
makes the isothermal-vs-THD bracket vacuous: matching the load forces the two
peak pressures together.  The bracket therefore comes from a separate pair of
runs at a common eccentricity, and this script attaches their profiles so the
figure can shade between them.

Usage:  $PY merge_bracket.py
"""
import json, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ferron as F
import run_ferron as RF

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 0.5345          # the eccentricity the THD case needs for 4 kN
NTH, NZ = 720, 40


def load_profile(tag):
    """Read the last groove-iteration case directory for `tag` and reduce it."""
    dirs = sorted(d for d in os.listdir(os.path.join(HERE, "_work"))
                  if d.startswith(f"ferron_bracket_{tag}_g"))
    if not dirs:
        return None
    d = os.path.join(HERE, "_work", dirs[-1])
    f = os.path.join(d, "out", "flat", "processor0", "solution_0_0.vts")
    if not os.path.exists(f):
        return None
    groove = None
    for line in open(os.path.join(d, "config.txt")):
        if line.startswith("inlet_groove"):
            groove = float(line.split("=", 1)[1].split()[0])
    fp = RF.read_vector3(f, "Fp")
    res = dict(p=RF.read_field(f, NTH, NZ, "p"),
               theta=RF.read_field(f, NTH, NZ, "film_content"),
               h=RF.read_field(f, NTH, NZ, "h"),
               T=RF.read_field(f, NTH, NZ, "T"),
               mu=RF.read_field(f, NTH, NZ, "mu"),
               fx=fp[0], fy=fp[1], nth=NTH, nz=NZ,
               elrod_iters=-1, mc=0.0, ns=0, R_p=0.0, R_T=0.0,
               coupling_iters=-1, converged=True, groove_deg=groove,
               eps=EPS, rpm=F.REF_RPM,
               thermal="ISOTHERMAL" if tag == "ISOTHERMAL" else "ENERGY_EQUATION",
               h_wall=0.0)
    m = RF.metrics(res)
    m["profile"] = RF.profile(res)
    return m


def main():
    dst = os.path.join(HERE, "ferron_results.json")
    fr = json.load(open(dst)) if os.path.exists(dst) else {}
    bracket = {}
    for tag in ("THD", "ISOTHERMAL"):
        m = load_profile(tag)
        if m:
            bracket[tag] = m
            print(f"  {tag}: p_max = {m['p_max_MPa']:.4f} MPa, W = {m['W_N']:.0f} N, "
                  f"T_max = {m['T_max_C']:.2f} C")
    if not bracket:
        print("no bracket runs found")
        return
    bracket["eps"] = EPS
    fr["bracket"] = bracket
    json.dump(fr, open(dst, "w"), indent=1)
    thd, iso = bracket["THD"]["p_max_MPa"], bracket["ISOTHERMAL"]["p_max_MPa"]
    meas = F.measured_peak_pressure_mpa()[1]
    print(f"\nat eps = {EPS}: isothermal {iso:.4f} MPa is {(iso-thd)/thd:+.1%} above THD "
          f"{thd:.4f} MPa")
    print(f"measured peak {meas:.2f} +/- {F.PRESSURE_DIGITISE_ERR_MPA} MPa "
          f"{'sits inside' if thd <= meas <= iso else 'is NOT enclosed by'} the bracket")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
