# Fluid directory

Oil-refrigerant pairs for the pancake GUI. Each pair is one folder:

```
fluids/
  <pair id>/
    pair.txt       header: names, source, density law, ...
    solution.csv   saturation pressure and viscosity table
```

The GUI lists every folder that contains `pair.txt`. Pairs can be copied,
edited and created in the GUI (Fluid > Open fluid library) or in a text
editor. The folder location is shown in the GUI and can be changed in
Tools > Settings.

## pair.txt

`key = value` lines; `#` starts a comment.

| Key | Unit | Meaning |
|---|---|---|
| `name` | | Display name |
| `oil`, `refrigerant` | | Names |
| `source` | | Where the data comes from (citation, supplier chart, "synthetic") |
| `note` | | Shown with the pair, e.g. "example only, not for design" |
| `solver_species` | | Gas model the solver uses for released gas: `R410A` (real-fluid equation of state), `PROPANE` or `AIR` (ideal gas) |
| `density_intercept_kg_m3` | kg/m^3 | Liquid density at 0 degC |
| `density_slope_kg_m3_K` | kg/(m^3 K) | Density change per kelvin (negative) |
| `oil_vapour_pressure_Pa` | Pa | Pressure at which the oil itself cavitates (vapour cavitation) |
| `specific_heat_J_kgK` | J/(kg K) | Liquid specific heat |
| `thermal_conductivity_W_mK` | W/(m K) | Liquid thermal conductivity |
| `table` | | Table file name, normally `solution.csv` |
| `gas_viscosity_Pa_s` | Pa s | Optional. Free-gas viscosity; needed for gas transport |
| `dissolved_heat_capacity_J_kgK` | J/(kg K) | Optional. Needed for thermal cases with gas release |
| `desorption_energy_J_kg` | J/kg | Optional. Energy taken up when dissolved refrigerant is released |
| `gas_thermal_conductivity_W_mK` | W/(m K) | Optional. Free-gas conductivity |

## solution.csv

The solver's `solution_table_file` format. Header and columns:

```
c_mass_fraction,T_C,P_MPa,viscosity_mm2s
```

| Column | Unit | Meaning |
|---|---|---|
| `c_mass_fraction` | kg/kg | Refrigerant content of the oil |
| `T_C` | degC | Temperature |
| `P_MPa` | MPa, absolute | Saturation pressure of this content at this temperature; blank for pure oil (content 0) |
| `viscosity_mm2s` | mm^2/s | Kinematic viscosity of the solution at this content and temperature |

Rows with the same content form one curve. One table gives both:

- the viscosity at any content and temperature (interpolated in log viscosity);
- the saturation pressure at any content and temperature, and the content
  that is saturated at a given pressure and temperature.

Rules checked when a pair is loaded or saved:

- at least two content curves with a saturation pressure, each with at least
  two temperatures, and pure-oil rows without a pressure;
- at fixed pressure the saturated content falls as temperature rises;
- every value is a finite positive number.

**Range:** a case whose ambient state (content, temperature, pressure) lies
outside the table turns the Fluid step yellow and marks the results. The solver
holds the end values of the table there (no extrapolation), so such results are
not reliable.

## Shipped example

`example_vg68_rx` is **synthetic, example only, not for design**. It is generated
by `tools/make_example_fluid.py` from the ISO VG 68 grade point, an assumed
viscosity index, Raoult's law and an invented vapour pressure for a fictitious
refrigerant "RX". The synthetic load table `examples/load_table_example.csv`
comes from the same script.

LGE material data (FW68D + R410A) is not in this folder. It is test data only,
under `tests/data/lge/`.
