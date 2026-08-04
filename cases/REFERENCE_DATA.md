# Reference data to digitize (user-provided)

Digitized experimental/numerical curves live in `cases/ref_data/`. A plotting agent
reads these and overlays solver results, so follow the format exactly.

## File format

- One CSV per curve: `cases/ref_data/case01/ferron1983_fig6_3000rpm_6kN.csv`
  (paper_figure_condition; lowercase, underscores).
- Header lines start with `#`: paper, figure, curve label, axis quantities and UNITS
  exactly as printed in the paper, operating condition (speed, load), and the
  abscissa convention (e.g., "theta measured from groove center, direction of rotation").
- Then two columns `x,y`, no thousands separators, `.` decimal.
- Keep the paper's axis units; the plotting agent converts (units must be in the header).
- Contours/boundaries (e.g., cavity extent): same format, columns are the boundary
  trace (theta, z), one file per boundary.

Example:

```
# paper: Ferron, Frene, Boncompain 1983, J. Lubr. Technol. 105:422
# figure: 6
# curve: experiment, N = 2000 rpm, W = 4000 N
# x: theta [deg], from supply groove, in rotation direction
# y: pressure [MPa], gauge
0.0,0.02
10.0,0.11
...
```

## What to provide, by case

### Case 01 — Ferron 1983 (priority 1; verify figure numbers against the PDF, they are from memory)

| # | Quantity | Where in paper | Conditions wanted |
|---|---|---|---|
| 1 | Midplane pressure vs theta, experiment points | pressure-distribution figure (~Fig. 6) | 2+ speed/load combos, note each |
| 2 | Bush inner-surface midplane temperature vs theta | temperature figure (~Fig. 9) | same combos as #1 |
| 3 | Eccentricity ratio vs load | eccentricity figure | each speed shown |
| 4 | Oil flow rate vs speed or load (optional) | flow figure | as available |

### Case 02 — static locus (priority 2; tables, not figures)

| # | Quantity | Source | Format |
|---|---|---|---|
| 1 | S (Sommerfeld) vs eps, attitude angle, friction var., flow | Someya Databook, plain cylindrical tables (pick L/D nearest 0.5 and 1.0) | one CSV per L/D, columns: S,eps,phi_deg,... (name columns in header) |
| 2 | Same from Raimondi & Boyd 1958, L/D = 1 and 1/2 | charts/tables | same |

Solver will be run at the table L/D values, not at 0.8.

### Case 05 — gaseous cavitation (priority 3)

| # | Quantity | Where | Notes |
|---|---|---|---|
| 1 | Pressure profiles inside cavitated zone (axial and/or circumferential) | Etsion & Ludwig 1982, J. Lubr. Technol. 104:157 | one CSV per profile, note theta or z station |
| 2 | Cavity extent boundary (theta-z trace) | Etsion & Ludwig 1982 photos/sketches | one CSV per condition |
| 3 | Cavity pressure/temperature maps (optional) | Braun & Hendricks 1984, ASLE Trans. 27:1 | profiles as CSVs |

Cases 00 (internal) and 03 (endpoint checked against 01/02 data) need no digitization.
Case 04 (secondary THD) will get its own list if we pursue it.
