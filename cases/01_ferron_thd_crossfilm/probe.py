"""Config-validation + wall-BC-mapping probe at a coarse grid.

Confirms (a) the three wall variants parse and run, and (b) that the large-HTC
limit really does pin the film-surface temperature onto
thermal_ambient_temperature, which is how a "fixed 313.15 K wall" is expressed
in the conjugate cross-film model (it has no prescribed-wall-temperature option).
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import solver_io as S

CASE_DIR = Path(__file__).resolve().parent
TEMPLATE = (CASE_DIR / "config_template.txt").read_text()
VARIANTS = {
    "walls_adiabatic":   dict(shaft_htc=0.0,  bush_htc=0.0,  bush_k=50.0),
    "bush40_Jadiabatic": dict(shaft_htc=0.0,  bush_htc=1e7,  bush_k=1e7),
    "bush40_J40":        dict(shaft_htc=1e7,  bush_htc=1e7,  bush_k=1e7),
}
out = {}
for name, w in VARIANTS.items():
    d = CASE_DIR / "probe_cfg" / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "config.txt").write_text(TEMPLATE.format(
        e=repr(0.574 * 1.52e-4), attitude=repr(-36.3),
        n_theta=180, n_z=48, coupling_max_iters=200, coupling_relaxation=repr(1.0),
        q_supply=repr(1.0e-4), groove_mixing_model="PINNED", **{k: repr(v) for k, v in w.items()}))
    rc, wall = S.run_case(d)
    m = S.frame_metrics(d) if rc == 0 else None
    tb = S.midplane_profile(d, "cross_film_t_bush")
    ts = S.midplane_profile(d, "cross_film_t_shaft")
    out[name] = {"rc": rc, "wall_s": round(wall, 1),
                 "T_max_C": (m or {}).get("T_max", 0) - 273.15 if m else None,
                 "T_bush_range_C": [min(tb[1]) - 273.15, max(tb[1]) - 273.15] if tb else None,
                 "T_shaft_C": ts[1][0] - 273.15 if ts else None,
                 "picard": (m or {}).get("coupling"),
                 "liquid_residual": (m or {}).get("liquid_residual")}
    print(name, json.dumps(out[name]), flush=True)
(CASE_DIR / "probe_cfg" / "probe.json").write_text(json.dumps(out, indent=2))
