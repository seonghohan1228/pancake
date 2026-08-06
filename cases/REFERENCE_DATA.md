# Reference data to digitize (user-provided)

> Campaign map: `cases/INDEX.md`. Ladder and decision gates:
> `cases/VALIDATION_PLAN.md`.

Digitized experimental/numerical curves live in `cases/ref_data/`. A plotting agent
reads these and overlays solver results, so follow the format exactly.

`ref_data/` is shared by every case and is **not** a case directory - it has no
stage number, no runs and no verdict. Every case reads it by relative path
(`CASE_DIR.parent / "ref_data" / ...`), so it must stay where it is. Delivered
files as of 2026-08-05: Ferron 1983 Fig. 4 (pressure + temperature + source PNG),
Fig. 5 (eccentricity + PNG), Fig. 6 (pressure + temperature + PNG).

## File format

- One CSV per curve: `cases/ref_data/ferron1983_fig6_4000rpm_6kN_pressure.csv`
  (paper_figure_condition_quantity; lowercase, underscores).
- Header lines start with `#`: paper, figure, curve label, axis quantities and UNITS
  exactly as printed, operating condition, and the abscissa convention (theta origin
  and direction). Flag any uncertainty in the header - the delivered Fig. 6 files
  are the model to copy.
- Two columns `x,y`, `.` decimal, no thousands separators.
- Keep the paper's axis units; the plotting agent converts.
- Contours/boundaries: same format, columns are the boundary trace (theta, z).
- Tables (not figures): one CSV per table, named columns in the header comment.

## Status and needs, by priority

### DELIVERED
- ferron1983_fig6_4000rpm_6kN_pressure.csv + _temperature.csv (+ source png). Used
  by Case 01-prelim; peak location agreed to ~3 deg with the stated theta origin.
- **etsion1982_table2_extent.csv** (2026-08-05, agent). Etsion & Ludwig 1982
  Table 2: theta_0, p_min, theta_st, theta_end for all 8 conditions.
  **TRANSCRIBED from the PDF text layer, not digitised** - the extraction emits the
  table column-major with the headers detached, so every value was re-attached to
  its column by cross-checking against seven explicit statements in the Results
  text; the checks are listed in the CSV header. Also records the rig, the oil
  properties, the +/-5 deg angular uncertainty the paper states, and the fact that
  **Table 2 carries no axial cavity width** (so the axial comparison of
  `cases/05_etsion_extent/` is N/A, not pending). Used by Case 05.

### 1. Case 01-THD - blocks the THD validation run (see 01_ferron_thd/README.md)

| # | What | Where in Ferron 1983 | Form |
|---|---|---|---|
| 1 | Oil viscosity at TWO or more temperatures, or the stated mu(T) law | oil properties / test conditions section (their oil table) | short CSV or just type the numbers in chat |
| 2 | Thermal boundary treatment used in their theory + any measured shaft/journal temperature | theory section + results text | numbers in chat is fine |
| 3 | Measured eccentricity ratio and attitude angle at 4000 rpm / 6 kN, if stated | text or eccentricity figure | numbers in chat |
| 4 | Eccentricity vs load points (each speed) | eccentricity figure | one CSV per speed |

### 2. Case 02 - solver side done, waiting on reference values

| # | What | Source | Form |
|---|---|---|---|
| 1 | S vs eps, attitude angle (min); friction/flow variables (nice to have) at L/D = 1 and 0.5 | Raimondi & Boyd 1958 tables/charts, and/or Someya Databook cylindrical-bearing tables | one CSV per L/D per source; note the source's cavitation BC (Reynolds/half-Sommerfeld) in the header |

Solver curves already produced with a 1-deg ambient groove at h_max (flooded
framing) - if your table's bearing has a different groove arrangement, say so and
the sweep is rerun to match.

### 3. Case 05 - gaseous cavitation (solver support landed, codex/fbns-gas)

| # | What | Where | Form | status |
|---|---|---|---|---|
| 1 | Pressure profiles inside the cavitated zone | Etsion & Ludwig 1982, Figs. 6/7 (4 panels each, p/p_s vs Z/L, 8-9 theta stations per panel) | one CSV per profile, note theta or z station in header | still needed - blocks sub-case 05b, which is separately blocked on `p_void` convergence |
| 2 | Cavity extent boundary (theta-z trace) per condition | same paper, Figs. 5(a)-(d) photographs | one CSV per condition | still needed, and **optional**: Table 2 already gives theta_st / theta_end numerically, which is what Case 05 scored |
| 3 | Operating/fluid parameters of their rig (geometry, speed, oil, air content if stated) | same paper | numbers in chat | **CLOSED** - all printed in the paper and transcribed into `etsion1982_table2_extent.csv`. Air content is NOT stated; Case 05 uses Henry saturation at atmospheric with an Ostwald coefficient of 0.09 and records it as UNVERIFIED |
| 4 | (optional) cavity pressure/temperature maps | Braun & Hendricks 1984, ASLE Trans. 27:1 | profiles as CSVs | not in the library (`UNFINISHED_BUSINESS.md` §3.6) |

### 4. Case 00c #2 - 2D contour verification (optional but recommended)

| # | What | Source | Form |
|---|---|---|---|
| 1 | 2D pressure contours + cavitation boundary + midplane profile, grooved journal bearing case | Ausas et al. 2007, ASME J. Tribol. 129:868 (or say "use Woloszynski 2015" - open access, I digitize it myself) | contour traces + profile CSVs |

Cases 00/00b (internal) and 03 (endpoint vs 01/02) need nothing from you.
