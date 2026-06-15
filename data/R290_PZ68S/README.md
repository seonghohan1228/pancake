# R290 / PZ68S property lookup

Get oil-solution **viscosity** and **solubility** from **pressure** and **temperature**.

```python
from r290_pz68s import properties
viscosity_mm2s, solubility_pct = properties(P_MPa=1.0, T_C=40)
# -> (6.843, 21.746)
```

**Function**

`properties(P_MPa, T_C) -> (viscosity_mm2s, solubility_pct)`
- `viscosity_mm2s` — kinematic viscosity of the oil + dissolved-R290 solution [mm²/s]
- `solubility_pct` — dissolved R290 mass fraction [%]

**Files**
- `r290_pz68s.csv` — data grid: `P_MPa, T_C, solubility_pct, viscosity_mm2s`
- `r290_pz68s.py` — the lookup function above

**Range:** P 0.1–4.0 MPa, T −60…160 °C (inputs are clamped to this box). Beyond the
saturation boundary, values are held at the saturated-solution value.
