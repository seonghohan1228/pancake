# Case 01-PRELIM

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Ferron et al. (1983) Fig. 6 point (4000 rpm, 6 kN), **isothermal**. The first
stage of the Ferron family: place the solver at the operating point with no
thermal physics at all and see how far an isothermal film gets on its own.

Two viscosities are run to bracket the missing THD effect - 0.0277 Pa s (the
40 C supply value) and 0.0166 Pa s (a rough 55 C effective value). The experiment
should sit between them if the only thing missing is the viscosity drop. It does;
`SUMMARY.md` has the numbers.

Superseded as a *result* by `01_ferron_thd`, which removes the hand-picked
viscosity. Kept because two things in it are still load-bearing for the whole
family:

1. **The sign convention**, derived empirically here rather than assumed - the
   pressure force sits in the -y half plane while the bearing is pushed in +y, so
   the STATIC target is `Fp = (0, -6000)` N. Every later case inherits it.
2. **The theta mapping**, `solver_theta = (90 + x) mod 360`, supported by the
   fact that the peak location is reproduced to about 3 deg and is
   viscosity-insensitive. Every later case plots against `(theta - 90) mod 360`
   with no fitted shift on this basis.

`02_static_locus/run_locus.py` imports `solver_io.py` from this directory, and
`01_ferron_vs_reynsim` reads this case's `equilibrium.json` for its pancake
column. Do not move or rename this directory without fixing those two.

## Configuration

Base config `cases/00_grid_dt_independence/grid_G6_fbns/config.txt`: R = 0.05 m,
c = 1.52e-4 m, L = 0.08 m, G6 720 x 192, omega = 418.879 rad/s,
`STEADY_STATE` + `STATIC` + `ISOTHERMAL`, `cavitation_model = FBNS`,
p_cav = 1.0e5 Pa, `inlet_groove = 90 18 1.713e5`, Dirichlet ambient axial ends,
`mpiexec -n 2`.

## Files

- `config_template.txt` - parametrised config (e, attitude, mu).
- `solver_io.py` - run driver and VTK/pvts readers (shared with case 02).
- `force_balance.py` - the 2-D damped-Newton driver on (eps, attitude).
- `fb_iter_00 .. fb_iter_13` - one archived run per solver call, mu40 in
  `00..06`, mu55 in `07..13`; `probe/` is the initial orientation run.
- `fb_log.json`, `equilibrium.json` - every iterate, and the two converged points.
- `plot_ferron.py` -> `ferron_midplane_pressure.png`.
- `progress.log` - running record.
