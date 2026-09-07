# FVC68B–R410A viscosity reference data

`fvc68b_r410a_1999_viscosity.csv` transcribes Table 35, printed page 139
(PDF page 156), in Gunsel and Pozebanchuk, *Elastohydrodynamic Lubrication
with Polyolester Lubricants and HFC Refrigerants*, final report volume I,
DOE/CE/23810-102, ARTI 670-54400, April 1999.
[Original DOE/OSTI report](https://www.osti.gov/servlets/purl/8212).

The report identifies its ISO-68 PVE as **Idemitsu FVC68B**. It does not
identify that oil as FW68D. Table 35 attributes the property data to its
reference 77, a November 1997 communication from L. Homolish at
Idemitsu-Kosan. The report's bearing-scale measurements concern EHL film
thickness; these viscosity data do not validate a journal-bearing cavity.

The source's refrigerant weight percentages 0, 15, 20 and 30 are converted
to mass fractions 0, 0.15, 0.20 and 0.30. Temperatures remain degrees Celsius.
Dynamic viscosity is converted from cP to Pa s by multiplying by 0.001.
All 12 values were checked visually against the rendered source table.

The corresponding configuration is:

```text
liquid_viscosity_model = COMPOSITION_TABLE
viscosity_composition_file = cases/ref_data/fvc68b_r410a_1999_viscosity.csv
```

The CSV has a rectangular grid with columns `c_mass_fraction,T_C,mu_Pa_s`.
The solver interpolates log viscosity in actual dissolved mass fraction
and temperature. It rejects extrapolation beyond 23–65°C or 0–30 wt%.
This table therefore does **not** cover the user's full 20–80°C range, and
must not be selected as an FW68D material model. It provides a traceable
data-import/interpolation verification fixture. The concentration and
temperature dependence of FW68D must be supplied separately.
