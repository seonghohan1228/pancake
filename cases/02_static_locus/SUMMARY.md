# Case 02 - Static eccentricity locus, plain journal bearing

Purpose: produce the solver side of an eccentricity-vs-Sommerfeld and
attitude-angle locus that can be laid over the Raimondi & Boyd (1958) / Someya
tables once the user supplies them.

## Configuration

| item | value |
|---|---|
| geometry | R = 0.05 m, c = 1.52e-4 m; L = 0.1 m (L/D = 1) and L = 0.05 m (L/D = 0.5) |
| grid | n_theta = 720; n_z = 256 (L/D = 1), n_z = 128 (L/D = 0.5) |
| fluid | mu = 0.0277 Pa s, rho = 860 kg/m3 |
| speed | omega = 418.879 rad/s (4000 rpm), N = 66.667 rev/s |
| mode | `STEADY_STATE` + `STATIC` + `ISOTHERMAL`, `cavitation_model = FBNS`, p_cav = 1.0e5 Pa |
| axial ends | Dirichlet 101325 Pa (fully flooded / submerged) |
| sweep | `attitude_angle_deg = -90` fixed, eps = 0.1 ... 0.9; load magnitude and direction are OUTPUTS |
| ranks | `mpiexec -n 2` |

Sommerfeld number: `S = mu * N * L * D * (R/c)^2 / W`, with N = omega / 2pi in
rev/s and D = 2R. For this geometry `mu * N * L * D * (R/c)^2` = 1998.5 N at
L/D = 1 and 999.3 N at L/D = 0.5.

Attitude angle: the bearing moves while the journal is fixed, so the
journal-relative line of centres points at `psi = attitude_angle_deg` and the
equivalent journal load direction is the pressure-force direction
`atan2(Fpy, Fpx)`. The reported angle is `phi = wrap180(psi - angle(Fp))`,
positive in the +theta (rotation) direction - the Raimondi & Boyd / Someya
definition. Verified against case 01: at eps = 0.482 the balanced Ferron point
gives phi = 58.7 deg, and the eps = 0.6 / psi = -90 orientation probe gives
phi = 50.2 deg, both in the expected plain-bearing range.

Side flow is the total volumetric liquid outflow through the two axial ends,
recomputed from the final frame with exactly the FBNS Dirichlet boundary link
the solver assembles:
`q_face = (h^3 / 12 mu) * (R dtheta) / (dz/2) * (p_cell - p_bc)`, summed over
faces with `q_face > 0` and divided by rho. `side_flow_net` (net of inflow) is
also stored in the CSVs.

## Solver refusal check and the feed anomaly

The task asked to delete the `inlet_groove` line. **The solver accepts a config
with no inlets** - it does not refuse, warn, or fail; `cfg.inlets.empty()` is an
early-return in every consumer (`reynolds.cpp:192`, `film_thickness.cpp:58`,
`energy.cpp:433`, `diagnostics.cpp:44`).

**But the resulting bearing is starved, not flooded.** With mass-conserving FBNS
and no feed except axial diffusion through the ends, the film settles into a
state with a cavitated area fraction of 0.76-0.86 at EVERY eccentricity, and it
carries roughly 1/20 of the flooded load:

| eps | W (no inlet) | W (flooded) | cav frac (no inlet) |
|---|---|---|---|
| 0.1 | 52 N | 1518 N | 0.764 |
| 0.6 | 813 N | 16568 N | 0.852 |
| 0.9 | 36773 N | 106302 N | 0.854 |

Net axial boundary flux is ~4e-14 kg/s, so that starved state is a genuine
converged steady solution of the mass-conserving problem, not a failure. It is
simply not the fully flooded film the Raimondi & Boyd tables assume.

The documented fallback was therefore used for the reported locus: a
**1-degree ambient-pressure groove at theta = 90 deg**, `inlet_groove = 90 1
101325`. With `attitude_angle_deg = -90` the maximum film thickness is exactly
at theta = 90 deg, so this groove pins a complete film at h_max - precisely
where the classical solution assumes the full film starts - while adding no
supply overpressure (the groove is at ambient). The no-inlet runs are kept for
the record in `LD1_noinlet/`, `LD05_noinlet/`, `locus_LD1_noinlet.csv`,
`locus_LD05_noinlet.csv`.

## Locus, L/D = 1 (`locus_LD1.csv`, runs in `LD1/eps_XX/`)

| eps | S | phi (deg) | W (N) | p_max (Pa, abs) | torque (N m) | side flow (m3/s) |
|---|---|---|---|---|---|---|
| 0.1 | 1.3160 | 79.43 | 1518.5 | 3.795e5 | 6.037 | 2.559e-5 |
| 0.2 | 0.6278 | 73.86 | 3182.9 | 7.000e5 | 6.165 | 5.086e-5 |
| 0.3 | 0.3869 | 68.25 | 5164.4 | 1.112e6 | 6.394 | 7.582e-5 |
| 0.4 | 0.2594 | 62.56 | 7704.4 | 1.692e6 | 6.748 | 1.005e-4 |
| 0.5 | 0.1780 | 56.68 | 11224.4 | 2.581e6 | 7.278 | 1.249e-4 |
| 0.6 | 0.1206 | 50.50 | 16568.4 | 4.087e6 | 8.076 | 1.491e-4 |
| 0.7 | 0.0776 | 43.82 | 25740.5 | 7.009e6 | 9.342 | 1.732e-4 |
| 0.8 | 0.0445 | 36.19 | 44922.9 | 1.413e7 | 11.604 | 1.971e-4 |
| 0.9 | 0.0188 | 26.43 | 106301.7 | 4.291e7 | 16.993 | 2.210e-4 |

## Locus, L/D = 0.5 (`locus_LD05.csv`, runs in `LD05/eps_XX/`)

| eps | S | phi (deg) | W (N) | p_max (Pa, abs) | torque (N m) | side flow (m3/s) |
|---|---|---|---|---|---|---|
| 0.1 | 4.2015 | 81.93 | 237.8 | 1.897e5 | 3.014 | 1.494e-5 |
| 0.2 | 2.0010 | 75.28 | 499.3 | 2.966e5 | 3.067 | 2.988e-5 |
| 0.3 | 1.2076 | 68.60 | 827.4 | 4.458e5 | 3.160 | 4.481e-5 |
| 0.4 | 0.7798 | 61.87 | 1281.2 | 6.765e5 | 3.305 | 5.974e-5 |
| 0.5 | 0.5067 | 55.06 | 1971.9 | 1.068e6 | 3.523 | 7.468e-5 |
| 0.6 | 0.3181 | 48.08 | 3140.9 | 1.809e6 | 3.853 | 8.965e-5 |
| 0.7 | 0.1846 | 40.80 | 5413.3 | 3.440e6 | 4.385 | 1.047e-4 |
| 0.8 | 0.0915 | 32.91 | 10913.4 | 8.061e6 | 5.356 | 1.198e-4 |
| 0.9 | 0.0312 | 23.53 | 32017.3 | 3.071e7 | 7.752 | 1.351e-4 |

## Figure

`static_locus.png` (single 300 dpi PNG, lab style, `plot_locus.py`): left panel
eps vs S with a logarithmic S axis, right panel phi vs eps, both L/D on the same
figure with distinct open markers (squares for L/D = 1, circles for L/D = 0.5).

## Notes for the table comparison

- Both loci are monotone and smooth over the whole eps range; no run failed and
  no outer-iteration limit was hit.
- The two attitude-angle curves nearly collapse, crossing near eps = 0.3, with
  L/D = 0.5 slightly above at low eps and slightly below at high eps. That is
  the qualitative shape of the Raimondi & Boyd locus family.
- Side flow is close to linear in eps and, at equal eps, scales roughly with the
  number of end faces per unit length; the axial cell size dz = 3.906e-4 m is
  identical for the two geometries by construction (n_z was halved with L).
- Cavitated area fraction runs 0.39-0.46 across the sweep, i.e. slightly less
  than the classical half-Sommerfeld 0.5, as expected for a mass-conserving
  model with a finite reformation zone.
- A quick sanity check against remembered textbook Raimondi & Boyd values (NOT a
  verified digitisation; replace with the user's tables) puts every point within
  about 2 % in S and 0.5 deg in phi - e.g. L/D = 1, eps = 0.6: solver
  S = 0.1206 / phi = 50.50 deg vs the commonly quoted 0.121 / 50.6 deg.

## Files

- `config_template.txt` - parametrised config (e, L, n_z).
- `run_locus.py` - sweep driver; `python run_locus.py [LD1|LD05|all] [groove|noinlet]`.
- `locus_LD1.csv`, `locus_LD05.csv`, `locus.json` - flooded locus (reported).
- `locus_LD1_noinlet.csv`, `locus_LD05_noinlet.csv`, `locus_noinlet.json` - starved runs.
- `LD1/`, `LD05/`, `LD1_noinlet/`, `LD05_noinlet/` - one directory per run with config and results.
- `plot_locus.py`, `static_locus.png` - the deliverable figure.
- `progress.log` - running record.
