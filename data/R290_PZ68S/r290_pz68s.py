"""
R290 (propane) + PZ68S oil — liquid-solution properties from (pressure, temperature).

    from r290_pz68s import properties
    viscosity_mm2s, solubility_pct = properties(P_MPa=1.0, T_C=40)   # -> (6.843, 21.746)

properties(P_MPa, T_C) -> (viscosity_mm2s, solubility_pct)
    viscosity_mm2s : kinematic viscosity of the oil + dissolved-R290 solution [mm^2/s]
    solubility_pct : dissolved R290 mass fraction [%]

Input range: P 0.1-4.0 MPa, T -60..160 C (inputs outside are clamped to the box).
Beyond the saturation boundary the values are held at the saturated-solution value
(excess R290 leaves the liquid as free gas, handled separately).

Data: r290_pz68s.csv (same folder). Interpolation: linear in T, linear in ln(P),
viscosity in log10.
"""
import csv
import os
import math
import bisect

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "r290_pz68s.csv")


def _load():
    grid = {}
    with open(_DATA, newline="") as f:
        for r in csv.DictReader(f):
            grid[(round(float(r["P_MPa"]), 3), int(r["T_C"]))] = (
                float(r["solubility_pct"]),
                float(r["viscosity_mm2s"]),
            )
    pressures = sorted({k[0] for k in grid})
    temperatures = sorted({k[1] for k in grid})
    return grid, pressures, temperatures


_GRID, _PS, _TS = _load()


def _bracket(values, x):
    if x <= values[0]:
        return values[0], values[0], 0.0
    if x >= values[-1]:
        return values[-1], values[-1], 0.0
    i = bisect.bisect_right(values, x)
    return values[i - 1], values[i], (x - values[i - 1]) / (values[i] - values[i - 1])


def properties(P_MPa, T_C):
    """Return (viscosity_mm2s, solubility_pct) at pressure [MPa] and temperature [C]."""
    P = min(max(P_MPa, _PS[0]), _PS[-1])
    T = min(max(T_C, _TS[0]), _TS[-1])
    plo, phi, _ = _bracket(_PS, P)
    tlo, thi, ft = _bracket(_TS, T)
    fp = 0.0 if phi == plo else (math.log(P) - math.log(plo)) / (math.log(phi) - math.log(plo))

    def bilinear(idx, log=False):
        def val(p, t):
            v = _GRID[(round(p, 3), t)][idx]
            return math.log10(v) if log else v
        a = val(plo, tlo) * (1 - ft) + val(plo, thi) * ft
        b = val(phi, tlo) * (1 - ft) + val(phi, thi) * ft
        v = a * (1 - fp) + b * fp
        return 10 ** v if log else v

    viscosity = bilinear(1, log=True)
    solubility = bilinear(0, log=False)
    return round(viscosity, 4), round(max(0.0, min(100.0, solubility)), 4)


if __name__ == "__main__":
    for P, T in [(1.0, 40), (0.5, 60), (2.0, 100), (0.3, 25)]:
        v, s = properties(P, T)
        print(f"P={P} MPa, T={T} C -> viscosity={v} mm2/s, solubility={s} %")
