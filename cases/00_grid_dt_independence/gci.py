"""Richardson extrapolation and GCI (Celik et al. 2008, ASME JFE 130:078001).

Grids are ordered 1 = finest, 2 = medium, 3 = coarsest. Representative cell
size h_i = sqrt(A / N_i) with A the (identical) unwrapped bearing area, so the
refinement ratios reduce to r = sqrt(N_fine / N_coarse).
"""

import json
import math
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent


def apparent_order(phi1, phi2, phi3, r21, r32, tol=1e-12, itmax=500):
    """Celik eq. (2): solve p = |ln|e32/e21| + q(p)| / ln(r21) by fixed point."""
    e21 = phi2 - phi1
    e32 = phi3 - phi2
    if e21 == 0.0 or e32 == 0.0:
        return None, None
    ratio = e32 / e21
    s = 1.0 if ratio > 0 else -1.0
    p = 1.0
    for _ in range(itmax):
        try:
            q = math.log((r21 ** p - s) / (r32 ** p - s))
        except (ValueError, ZeroDivisionError):
            return None, s
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if not math.isfinite(p_new):
            return None, s
        if abs(p_new - p) < tol:
            p = p_new
            break
        p = 0.5 * p + 0.5 * p_new  # damped, the raw map can oscillate
    return p, s


def gci_triplet(phi1, phi2, phi3, r21, r32, p_fixed=None):
    """Full Celik table for one quantity. phi1 = finest grid."""
    p_obs, s = apparent_order(phi1, phi2, phi3, r21, r32)
    p = p_fixed if p_fixed is not None else p_obs
    out = {"phi1": phi1, "phi2": phi2, "phi3": phi3,
           "r21": r21, "r32": r32, "p_observed": p_obs, "sign": s, "p_used": p}
    if p is None:
        return out
    phi_ext21 = (r21 ** p * phi1 - phi2) / (r21 ** p - 1.0)
    phi_ext32 = (r32 ** p * phi2 - phi3) / (r32 ** p - 1.0)
    ea21 = abs((phi1 - phi2) / phi1)
    ea32 = abs((phi2 - phi3) / phi2)
    eext21 = abs((phi_ext21 - phi1) / phi_ext21)
    gci21 = 1.25 * ea21 / (r21 ** p - 1.0)
    gci32 = 1.25 * ea32 / (r32 ** p - 1.0)
    out.update({
        "phi_ext21": phi_ext21, "phi_ext32": phi_ext32,
        "ea21": ea21, "ea32": ea32, "eext21": eext21,
        "gci21_fine": gci21, "gci32_medium": gci32,
    })
    return out


def main():
    m = json.loads((CASE_DIR / "metrics.json").read_text())
    # 1 = finest (G5), 2 = medium (G4), 3 = coarse (G3)
    n = {k: m[k]["n_theta"] * m[k]["n_z"] for k in ("grid_G3", "grid_G4", "grid_G5")}
    r21 = math.sqrt(n["grid_G5"] / n["grid_G4"])
    r32 = math.sqrt(n["grid_G4"] / n["grid_G3"])

    results = {"r21": r21, "r32": r32, "N": n, "quantities": {}}
    for label, key in (("W", "W"), ("p_max", "p_max"),
                       ("friction_torque", "friction_torque"),
                       ("force_angle_deg", "force_angle_deg"),
                       ("cav_frac_0999", "cav_frac_0999")):
        phi1 = m["grid_G5"]["final"][key]
        phi2 = m["grid_G4"]["final"][key]
        phi3 = m["grid_G3"]["final"][key]
        entry = gci_triplet(phi1, phi2, phi3, r21, r32)
        entry["fixed_p2"] = gci_triplet(phi1, phi2, phi3, r21, r32, p_fixed=2.0)
        results["quantities"][label] = entry

    (CASE_DIR / "gci.json").write_text(json.dumps(results, indent=2))
    print(f"r21 = {r21:.6f}  r32 = {r32:.6f}")
    for label, e in results["quantities"].items():
        print(f"\n{label}")
        print(f"  phi (G3,G4,G5) = {e['phi3']:.8g}, {e['phi2']:.8g}, {e['phi1']:.8g}")
        print(f"  p_observed = {e['p_observed']}")
        if "gci21_fine" in e:
            print(f"  phi_ext21  = {e['phi_ext21']:.8g}")
            print(f"  ea21 = {e['ea21']*100:.4f}%  eext21 = {e['eext21']*100:.4f}%"
                  f"  GCI21_fine = {e['gci21_fine']*100:.4f}%")
            print(f"  ea32 = {e['ea32']*100:.4f}%  GCI32_medium = {e['gci32_medium']*100:.4f}%")
        f = e["fixed_p2"]
        if "gci21_fine" in f:
            print(f"  [p=2] phi_ext21 = {f['phi_ext21']:.8g}"
                  f"  GCI21_fine = {f['gci21_fine']*100:.4f}%"
                  f"  GCI32_medium = {f['gci32_medium']*100:.4f}%")


if __name__ == "__main__":
    main()
