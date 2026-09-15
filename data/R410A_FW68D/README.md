# FW68D / R410A property table

The [table](fw68d.csv) contains all 125 digitized points from the user's
[project JSON](source_project.json). The [original chart](source_chart.png) is
archived separately from the JSON's image field. No numerical points were
discarded or replaced by a fitted correlation.

Pressure is **absolute**, and concentration is **refrigerant mass divided by
total liquid-solution mass**, as requested by the user. The archived project's
`unspecified` basis labels are provenance, not the model's operating convention.

## Columns and range

| CSV column | Meaning | Solver units |
|---|---|---|
| `c_mass_fraction` | R410A mass fraction, 0 = pure oil, 1 = pure refrigerant | kg/kg |
| `T_C` | Calibrated chart temperature | K after adding 273.15 |
| `P_MPa` | Absolute pressure at a concentration/isobar intersection | Pa after multiplying by 10⁶ |
| `viscosity_mm2s` | Kinematic liquid viscosity | m²/s after multiplying by 10⁻⁶ |

Pure-oil rows have blank pressure because the oil curve has no pressure labels.
The other curves each contain the 11 isobars at 0.1, 0.2, 0.4, 0.6, 0.8, 1.0,
1.5, 2.0, 2.5, 3.0 and 4.0 MPa absolute.

| Refrigerant mass fraction | Points | Digitized temperature range, °C |
|---|---|---|
| 0 | 15 | -19.96 to 120.15 |
| 0.1 | 11 | -22.71 to 130.01 |
| 0.2 | 11 | -36.17 to 98.03 |
| 0.3 | 11 | -42.28 to 84.02 |
| 0.4 | 11 | -45.76 to 75.92 |
| 0.5 | 11 | -48.01 to 71.03 |
| 0.6 | 11 | -49.41 to 67.37 |
| 0.7 | 11 | -50.46 to 65.41 |
| 0.8 | 11 | -51.04 to 63.61 |
| 0.9 | 11 | -51.62 to 62.68 |
| 1 | 11 | -51.74 to 62.07 |

The chart calibration uses piecewise linear interpolation of `ln(T_C + 273.15)`
against x-pixel and `ln(nu)` against y-pixel, between adjacent saved anchors.
Pixels outside the calibration anchors use the nearest anchor segment. This is
needed for several cold, refrigerant-rich points left of the -40 °C anchor.
`target_T_C` is a digitization target; the recorded x-pixel determines the exported
temperature. For example, the pure-oil point near 40 °C decodes to 40.017739 °C
and 74.427078 mm²/s. The saved calibration is retained even where coarse anchors
limit accuracy. Exported decimal precision does not imply measurement accuracy.

## Configuration

Start with [config_r410a_fw68d.txt](../../config_r410a_fw68d.txt), or set:

```ini
fluid_property_model = SINGLE_PHASE
solubility_model = CURVE_TABLE
liquid_viscosity_model = KINEMATIC_TABLE
solution_table_file = data/R410A_FW68D/fw68d.csv
solution_table_clamp = true
kinematic_viscosity_log_interpolation = true
density_model = CONSTANT
rho = 940
density_temperature_intercept = 954.9
density_temperature_slope = -0.7284
```

| Option | Behavior |
|---|---|
| `solubility_model = CURVE_TABLE` | Interpolate saturated mass fraction from absolute pressure and temperature |
| `liquid_viscosity_model = KINEMATIC_TABLE` | Interpolate kinematic viscosity at actual concentration and temperature, then return μ = ρν |
| `kinematic_viscosity_log_interpolation = true` | Interpolate ln(ν) in temperature and concentration |
| `kinematic_viscosity_log_interpolation = false` | Interpolate ν linearly in temperature and concentration |
| `solution_table_clamp = true` | Hold endpoints outside each contributing curve's temperature range or the outer axis range |
| `solution_table_clamp = false` | Reject a query outside either contributing curve's range |
| `density_model = CONSTANT` | Use `rho`, with 940 kg/m³ in this example |
| `density_model = LINEAR_TEMPERATURE` | Use intercept + slope × T in °C, default 954.9 − 0.7284 T |

The GUI exposes these choices under Fluid property models. The existing HENRY,
TABLE, CONSTANT, EMPIRICAL and COMPOSITION_TABLE models remain available.
`fluid_property_model = CONSTANT` bypasses all property laws, including the
linear density law. Use SINGLE_PHASE or TWO_PHASE to enable them.

For R410A simulations, the existing solver requires `dissolved_gas_species = R410A`
and `cavitation_model = COUPLED`. The example supplies a small isothermal smoke
case. Its cavity pressure and free-gas viscosity are labelled placeholders; the
chart supplies neither of those properties nor heat of solution or exchange rates.

## Interpolation and limits

Viscosity is interpolated within each constant-concentration curve in temperature,
then between the two enclosing concentrations. Pressure does not independently
modify this ν(c,T) model. Treating equilibrium measurements as an actual-composition
law assumes that pressure effects at fixed c,T are negligible.

Solubility is interpolated within each isobar in temperature, then linearly between
absolute pressures. Its inverse returns the first pressure crossing of the same
interpolated isotherm. No zero-pressure or dilute-concentration data are invented.
With clamping enabled, the solubility table has a 0.1 mass-fraction lower endpoint
and a 1.0 upper endpoint; endpoint holding outside the measured region is a model
extension, not measured solubility. For zero dissolved charge, the existing bubble
pressure helper returns zero. Disable clamping when range violations should stop
a run. Exact curve nodes query only their own temperature range.

The two requested density laws apply to the whole liquid solution independently
of concentration and pressure. At 40 °C they give 940 and 925.764 kg/m³ respectively.
Neither law is a measured R410A liquid density model. The c=1 chart endpoint is
pure **liquid** refrigerant viscosity, distinct from free-gas viscosity and the
existing R410A gas EOS. A pure endpoint lookup does not establish the validity of
the bearing solver or these density approximations for a pure-refrigerant run.

## Regenerate and verify

From the repository root:

```sh
python data/R410A_FW68D/import_project.py
python -m unittest discover -s data/R410A_FW68D
ctest --test-dir build-linux --output-on-failure -R 'fw68d|composition_viscosity|property_tables'
```

An optional JSON argument imports a revised digitization. Recompile is unnecessary
after editing the CSV. The loader validates numeric fields, positive viscosity,
absolute pressure, mass fractions, unique curve temperatures and decreasing
solubility with temperature on each isobar. Relative file paths resolve from the
configuration file's directory, then the working directory.
