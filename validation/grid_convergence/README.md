# Grid-convergence study (WP-3 / AUDIT §C1)

Runs the steady eccentric journal case (`config_base.txt`, ε = 0.8,
Elrod-Adams) on three grid levels — 60×20, 120×40, 240×80 — and tabulates per
grid:

- `cavitated_fraction` — cavitation extent (fraction of cells with θ < 1),
  from the diagnostics CSV;
- `liquid_residual` — global mass-balance closure (should sit at
  `linear_rtol`);
- fluid force and friction torque (PowerShell script, parsed from the
  DETAILED step log).

Usage from the build directory (executables and `mpiexec` on PATH):

```powershell
powershell -File ..\validation\grid_convergence\run.ps1 -Ranks 2
```

```bash
../validation/grid_convergence/run.sh 2
```

Outputs `grid_convergence_summary.csv` in the working directory. Cavitation
extent and load should converge monotonically under refinement; archive the
summary alongside result reviews. To study scheme sensitivity, change
`theta_convection_scheme` in `config_base.txt` (UPWIND / LINEAR / VANLEER /
TYPE_DIFFERENCING) and rerun.
