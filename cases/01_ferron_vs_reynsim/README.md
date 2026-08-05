# Case 01-VS-REYNSIM

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Three-way overlay at the Ferron et al. (1983) conditions: the **experiment**, the
**previous-generation MATLAB/ReynSim solver** the group used before pancake, and
**pancake**. The point of the case is not another physics comparison - it is to
show, on one axis and in one frame, where the new solver sits relative to the old
one against the same measurements.

Results are in `comparison.json` (every reported number) and the three figures;
this case has no separate `SUMMARY.md`.

## Frame conventions, all stated and none fitted

- Abscissa: degrees from the supply-groove **centre**, in the direction of
  rotation. Our mesh puts the groove at theta = 90 deg, so solver curves are
  drawn against `(theta - 90) % 360`, exactly as in `01_ferron_prelim` and
  `01_ferron_thd`.
- ReynSim's own frame is ours + 90 deg; the derivation is in `reynsim_io.py`.
- Ordinate: gauge pressure, `p - 101325` Pa, for every curve.
- The ReynSim eccentricity locus in `journal_center_locus.csv` carries **no
  stated operating point** anywhere in the delivered data. It is flagged
  `locus_unverified` in `comparison.json` and must stay flagged.

## Deliverables

| figure | what it shows |
|---|---|
| `reynsim_midplane_pressure.png` | mid-plane p(theta) at 4000 rpm / 6 kN: experiment (points), ReynSim (line), pancake FBNS isothermal band (two thin lines) |
| `reynsim_eccentricity.png` | eps vs W: ReynSim's two numerical vintages, their own experimental digitisation (source unverified), our Fig. 5 digitisation rescaled 145/152, and our FBNS+THD sweep at 2000 rpm |
| `fig4_three_way.png` | the same three-way overlay at the Fig. 4 point (2000 rpm, 4 kN) |

Headline: ReynSim's peak is 1.507 MPa gauge against the experiment's 1.828 MPa,
i.e. **-17.6 %**, while pancake's two isothermal viscosities bracket the
measurement at -5.0 % and +6.9 %. Peak locations agree within about 3 deg for all
three.

## Cross-references (fix these before renaming anything)

- `compare_table.py` and `plot_vs_reynsim.py` read
  `cases/01_ferron_prelim/equilibrium.json`.
- `plot_vs_reynsim.py` and `plot_fig4_comparison.py` read `cases/01_ferron_thd/`.
- All three read `cases/ref_data/ferron1983_*.csv`.

## What the subdirectories are

- `reynsim_data/` - the delivered MATLAB/ReynSim material: pressure
  distributions, centreline profiles, the journal-centre locus, and the `.m`
  files that produced them. Read-only input to this case.
- `bin/` - the pinned Windows `pancake.exe` plus its MinGW runtime DLLs. This is
  the binary that produced the case's numbers, kept with the case so the run is
  reproducible; untracked (`*.exe` / `*.dll` are gitignored at the repo root).
- `smoke/` - the coarse-grid config smoke run.
- `diag/` - the four-run linear-solver diagnostic set (`serial`, `repeat`,
  `relax_0.5`, `relax_0.7`). `CHANGES.md` quotes this set for the block-Jacobi
  vs ASM Krylov-iteration table, so it has to stay reachable.
- `fig6_Jadiabatic/` - the Fig. 6 force-balance attempt; `fb_iter_03` failed
  (`rc=1`) and the run stops there. See `progress.log`.

## Files

`config_template.txt`, `force_balance.py`, `solver_io.py`, `reynsim_io.py`,
`plot_style.py`, `compare_table.py` -> `comparison.json`, `plot_vs_reynsim.py`,
`plot_fig4_comparison.py`, `progress.log`, `fb_run.out`.
