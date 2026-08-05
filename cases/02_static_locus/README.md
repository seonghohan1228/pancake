# Case 02 - static eccentricity locus

> Stage 02 (static locus). Campaign index: `../INDEX.md`.

Solver side of an eccentricity-vs-Sommerfeld and attitude-angle locus for a plain
journal bearing at L/D = 1 and L/D = 0.5, produced so it can be laid over the
Raimondi & Boyd (1958) or Someya tables.

**POSTPONED by the user (2026-08-05).** The reference tables are not in the
library: the Raimondi & Boyd part held there is Part III / the 120-deg partial
bearing, and Someya is a book. Load-capacity validation was folded into Case 01
instead, via Ferron Fig. 5 - experimental, and on the same bearing. The curves
here are complete and are kept as supporting material; `SUMMARY.md` has them in
full, and `REFERENCE_DATA.md` records what would be needed to close the case.

## Configuration

R = 0.05 m, c = 1.52e-4 m; L = 0.1 m (L/D = 1, n_z = 256) and L = 0.05 m
(L/D = 0.5, n_z = 128); n_theta = 720 throughout; mu = 0.0277 Pa s, rho = 860;
omega = 418.879 rad/s; `STEADY_STATE` + `STATIC` + `ISOTHERMAL`,
`cavitation_model = FBNS`, p_cav = 1.0e5 Pa; Dirichlet 101325 Pa axial ends;
`mpiexec -n 2`. The sweep fixes `attitude_angle_deg = -90` and steps
eps = 0.1 ... 0.9 - load magnitude and direction are **outputs**, not targets, so
no force balance is needed.

## The feed anomaly, and why the reported locus has a groove

The task asked for the `inlet_groove` line to be deleted. The solver accepts a
config with no inlets - it does not refuse, warn or fail - but the resulting
bearing is **starved, not flooded**: with mass-conserving FBNS and no feed, the
film settles at a cavitated area fraction of 0.76-0.86 at every eccentricity and
carries roughly 1/20 of the flooded load. Net axial boundary flux is about
4e-14 kg/s, so it is a genuine converged solution of the mass-conserving problem
- just not the film the classical tables assume.

The reported locus therefore uses the documented fallback: a **1-degree
ambient-pressure groove at theta = 90 deg** (`inlet_groove = 90 1 101325`). With
`attitude_angle_deg = -90` the maximum film thickness sits exactly at 90 deg, so
this pins a complete film at h_max - where the classical solution assumes the
full film starts - while adding no supply overpressure. The starved runs are kept
in `LD1_noinlet/`, `LD05_noinlet/` and their CSVs.

If the eventual reference table uses a different groove arrangement, say so and
the sweep is re-run to match.

## Definitions used in the CSVs

- `S = mu N L D (R/c)^2 / W`, N in rev/s, D = 2R. For this geometry the numerator
  group is 1998.5 N at L/D = 1 and 999.3 N at L/D = 0.5.
- Attitude angle `phi = wrap180(psi - angle(Fp))`, positive in the +theta
  (rotation) direction - the Raimondi & Boyd / Someya definition. The bearing
  moves and the journal is fixed, hence the subtraction.
- Side flow is recomputed from the final frame with exactly the FBNS Dirichlet
  boundary link the solver assembles, summed over outflow faces and divided by
  rho; `side_flow_net` is also stored.

## Files

- `config_template.txt`, `run_locus.py` - sweep driver;
  `python run_locus.py [LD1|LD05|all] [groove|noinlet]`. It imports
  `solver_io.py` from `../01_ferron_prelim`.
- `locus_LD1.csv`, `locus_LD05.csv`, `locus.json` - the reported flooded locus.
- `locus_LD1_noinlet.csv`, `locus_LD05_noinlet.csv`, `locus_noinlet.json` - the
  starved runs.
- `LD1/`, `LD05/`, `LD1_noinlet/`, `LD05_noinlet/` - one directory per run.
- `plot_locus.py` -> `static_locus.png`; `progress.log`, `sweep*.out`.
