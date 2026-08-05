# Case 01-WALLBC-FIELDS

> Stage 01 (Ferron validation family). Campaign index: `../INDEX.md`.

Ferron et al. (1983) Fig. 4 point (2000 rpm, 4 kN), the **final four-way wall
treatment comparison**, run with the film-surface temperatures and the wall heat
fluxes written out as fields so the mechanism can be shown and not just inferred.

This is the production batch of the wall-BC thread. It closes the sequence
`01_ferron_thd` (no cross-film) -> `01_ferron_thd_crossfilm` (cross-film, three
local-bush variants) -> `01_ferron_conjugate` (2-D bush, verdict campaign).
Results in `SUMMARY.md`.

First production use of the rank-robust linear solve (`f0e772d`): the
`rerun_cht_robin10` group closed its force balance here where it had previously
diverged on both hosts.

## Groups

Everything except the thermal wall block and the `output_fields` line is frozen
from `01_ferron_thd_crossfilm`; both varying pieces are supplied by
`force_balance.py` so the rest of the config cannot drift.

| group | film model | walls |
|---|---|---|
| `rerun_xf_adiabatic` | cross-film T(y), LOCAL 1-D bush | both adiabatic |
| `rerun_xf_bush40_J40` | cross-film T(y), LOCAL 1-D bush | both driven onto 40 C by the large-conductance limit (h = k = 1e7) |
| `favg_walls40` | **no** cross-film; depth-averaged T | Newton cooling at 40 C, journal h = 10000, bearing h = 7200 (the 01-THD values) |
| `rerun_cht_robin10` | cross-film T(y), conjugate 2-D bush (k = 250 W/(m K), R2 = 0.10 m, 16 radial cells) | Robin outer surface, h = 10 W/(m^2 K) to 313.15 K; shaft globally adiabatic |

Two things that travel with the numbers:

- `favg_walls40` has no film-surface fields to write - they do not exist under
  `cross_film_thermal_model = NONE`, and requesting them is a config error. Its
  wall fluxes are therefore **post-derived** as `q = h (T_wall - T_bar)`,
  positive into the fluid. Its temperature column is `T-bar`, not `T_bearing`,
  and is labelled as such everywhere.
- `rerun_cht_robin10` uses the **user-selected h = 10 W/(m^2 K)**, whereas the
  verdict campaign in `01_ferron_conjugate` used the paper's own h = 80. Both
  stay quotable; the difference between them was measured rather than assumed,
  and the frozen-geometry version of that measurement was later corrected. See
  `01_ferron_conjugate/README.md`.

## Warm starts

Each group starts from the converged equilibrium of its closest archived twin
rather than the cold (0.60, -40) guess. That removes one or two Newton steps per
pass and cannot change the answer - it is the same fixed point, approached from
closer in. Reruns reproduce their archived predecessors (eps exact to 6 digits
where conditions are identical), which is the check that this is true.

## Protocol

Frozen from `01_ferron_thd_crossfilm`: G6 720 x 192, FBNS, p_cav = 1.0e5 Pa,
Ferron polynomial viscosity, 18 deg fed groove at theta = 90 deg, damped-Newton
force balance to 0.1 % in |W| and 0.05 deg in force direction, two-pass
PINNED -> FLOW_WEIGHTED groove closure with Q_supply from the first pass's own
steady mass balance. Run on afml-r6, 24 pinned ranks, remote HEAD hash-verified,
with a cross-platform gate that matched Windows to all printed digits.

## Files

- `config_template.txt`, `force_balance.py` - the parametrised config and the
  four-group driver (`python force_balance.py [group ...]`).
- `solver_io.py`, `plot_style.py` - harness copies, unchanged.
- `equilibrium.json`, `extract.json` - converged points and every extracted
  profile; `driver.log`, `progress.log` - the running record.
- `plot_wallbc.py` -> `ferron_wallbc_temperature.png`,
  `ferron_wallbc_heatflux.png`.
- `<group>/fb_iter_NN/` - one archived solver call each, with its `config.txt`.

## Where the bulk data is

983 MB of VTK stays on afml-r6 at
`~/pancake/cases/01_ferron_wallbc_fields/`; only provenance and the extracted
profiles are mirrored locally. The coarse 180 x 24 config smoke for all four
groups is remote-only, at `.../01_ferron_wallbc_fields/smoke_coarse_180x24/`.
The extraction script that produced `extract.json` is archived at
`cases/90_diagnostics/remote_harness/extract.py`.
