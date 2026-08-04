# Case 00: Grid and timestep independence

First step of the validation campaign. Establishes the grid and dt used by all later
validation cases (Ferron THD benchmark geometry: R = 50 mm, L = 80 mm, c = 152 um hot clearance
per the paper (145 um cold + thermal expansion),
single 18-deg axial groove at 90 deg, 4000 rpm / 6 kN = Ferron Fig. 6 condition,
ISO-VG-32-like oil mu = 0.0277 Pa s).

## Runner instructions

1. Build: `cmake --preset windows-native-mingw` then build; binary `pancake.exe`
   (see BUILDING.md). Tests/ctest need mingw64 on PATH.
2. For each run directory below, from inside that directory:
   `mpiexec -n 4 <path-to>/pancake.exe --config config.txt`
   Use at most 4 ranks. n_theta of every case is divisible by 4.
3. Do NOT delete anything. Leave `results/` and the diagnostics CSV inside each run dir.
   Tee solver stdout/stderr to `run.log` in each run directory, and append one line per
   event (case started / finished / failed, wall-clock) to `progress.log` in this
   directory so the user can follow along by opening those files.
4. If STEADY_STATE fails to converge for a grid run, note it, then rerun with
   `solution_mode = TRANSIENT`, `end_t = 0.2`, `dt = 1e-5` (save as config_transient.txt
   in the same dir) and use the final state.

## Run matrix

### 00a grid (steady, fixed eccentricity 0.6, isothermal, JFO)

| run | n_theta x n_z | model |
|---|---|---|
| grid_G1 | 60 x 16 | isothermal |
| grid_G2 | 120 x 32 | isothermal |
| grid_G3 | 240 x 64 | isothermal |
| grid_G4 | 360 x 96 | isothermal |
| grid_G5 | 480 x 128 | isothermal |
| grid_G3_thd | 240 x 64 | THD (energy eq, empirical viscosity) |
| grid_G4_thd | 360 x 96 | THD confirmation of chosen grid |

### 00b timestep (transient, moving bearing, 6 kN load step, grid 240 x 64)

Bearing mass 10 kg (user-chosen test mass: slower, more oscillatory transient makes
the harsher dt test; not a Ferron rig value). Grid runs are STATIC and unaffected.

| run | dt [s] | deg/step at 4000 rpm |
|---|---|---|
| dt_T1 | 4e-5 | 0.96 |
| dt_T2 | 2e-5 | 0.48 |
| dt_T3 | 1e-5 | 0.24 |
| dt_T4 | 5e-6 | 0.12 |
| dt_T5 | 2.5e-6 | 0.06 |

If 00a selects a grid other than G3, update `n_theta_global`/`n_z_global` in all dt_T*
configs before running (note the change in the summary).

## Metrics to extract per run

00a: load W = sqrt(Fx^2 + Fy^2) from pressure_force_x/y (final diagnostics CSV row),
attitude angle of the force vector, peak pressure p_max (flat VTK), cavitation area
fraction (cells with film_content < 0.999), friction_torque.
THD runs additionally: max and mean temperature.

00b: equilibrium eccentricity ratio and attitude angle from bearing_x/bearing_y
(final state), settling time (time when |e(t) - e_final| < 1% of c), and whether the
trajectory is monotone or oscillatory.

## Reporting

Write `SUMMARY.md` in this directory:
- One table per sub-study: run vs metrics above, plus relative change from the
  next-coarser grid / next-larger dt.
- Richardson extrapolation and GCI on W and p_max using G3/G4/G5 (procedure:
  Celik et al. 2008, ASME J. Fluids Eng. 130:078001).
- Wall-clock time per run.
- Flag any run that hit max_outer_iters, stopped on non-positive film, or produced
  mass-balance drift in the diagnostics CSV.

Pass criteria are decided by the user, target guidance: <1% change in W and p_max
between consecutive grids, <1% change in equilibrium eccentricity between consecutive dts.
