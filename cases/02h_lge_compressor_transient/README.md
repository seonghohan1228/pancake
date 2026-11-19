# 02h LGE compressor OS bearing at 130 Hz (transient, GC+VC)

Campaign `lge-compressor-transient-20261119-<commit7>`. Part A adds a time-dependent load
input and implicit motion coupling to the solver and checks them. Part B runs the orbiting-scroll
(OS) bearing of the LGE compressor at 130 Hz with the crank-angle load table, gaseous plus
vaporous cavitation (GC+VC), and reports where the bubble region sits relative to the feed
hole over one revolution, with the minimum film thickness and the friction torque.

| File | Purpose |
|---|---|
| `common.py` | config rendering, solver runs (config sha256 reuse), `motion.csv` and VTK readers |
| `part_a_checks.py` | Part A solver checks on a coarse plain bearing |
| `part_b.py` | Part B stages: `prepare`, `b0`, `b1`, `b2`, `report` |

## Part A: load input and implicit coupling (solver)

Keys (PHYSICS.md 12.1-12.3):

- `load_model = CONSTANT | ROTATING | TABLE`. `CONSTANT` is `external_load_x|y|z`, unchanged.
  `ROTATING`: `load_magnitude`, `load_rotation_start_deg`, `load_rotation_speed` (one constant
  rate). `TABLE`: `load_table_file`, `load_table_abscissa = PHASE_DEG` (periodic, phase advances
  at `load_frequency`, offset `load_phase_offset_deg`) or `TIME_S` (not periodic),
  `load_table_format = POLAR | CARTESIAN`, `load_table_scale`, `load_table_angle_offset_deg`.
  Rows are interpolated linearly in their Cartesian parts; jumps are kept.
- `motion_coupling = EXPLICIT | IMPLICIT`, `motion_coupling_tolerance`,
  `motion_coupling_position_tolerance`, `motion_coupling_max_iters`, `motion_coupling_fd_step`,
  `motion_coupling_step_limit`.
- Output: `results/motion.csv` every step of a `MOVING_BEARING` transient; the load of the step
  in the `external_load_*` fields; `results/final_state.dat` at the end of a single-phase FBNS
  transient (restart through `steady_initial_state_file`).

Checks (`part_a_checks.py`): see the results table below.

## Part B: frame and mapping

The crank pin turns inside the OS boss; the OS does not rotate. The solver runs in the frame
that turns with the crank:

- The pin, with the feed hole, is at rest. It is the solver's non-rotating mesh surface.
- The OS bore turns at -omega in that frame. The angle axis is mirrored (theta' = -theta) so the
  bore surface moves at +omega, the solver's rotation sense.
- In the solver, the moving body (mass, external load) is the outer surface and the rotating
  surface is the inner one. The OS bore is the outer surface and the body that moves relative
  to the pin, so the solver's moving body is the OS directly: no sign flip of displacement or
  load. The rotating surface being the inner one in the solver and the outer one in the machine
  does not change the film (the Couette flux depends on the sum of the surface speeds). The
  journal friction torque the solver reports is the torque on its rotating surface; the
  friction power is torque x omega in both descriptions.
- Load: the crank drives the OS forward against the tangential gas force, so the film force
  on the OS points along the orbit velocity and the external (gas) load on the OS points
  against it. In the mirrored crank frame this is a fixed direction, theta' = 90 deg; only its
  size changes with crank angle. The pin side of the film is loaded on the leading side, and
  the feed hole sits 180 deg from it, at theta' = 90 deg (unloaded side, assumption).
- Fictitious forces of the rotating frame on the bearing's relative motion are dropped:
  m omega^2 = 0.46 N per um, Coriolis about 1 N, against a film stiffness of order 1000 N per um.
  The OS centrifugal force on the 3 mm orbit (1,373 N) is carried by the wrap flanks (LGE sheet).

## Part B: inputs

| Item | Value | Source |
|---|---|---|
| Bore, length, radial clearance | 21 mm, 20 mm, 30 um (loosest fit; drawing 17.5-30, mean 23.75) | LGE drawing, Seongho |
| Speed | 130 Hz | LGE |
| Shell (end) pressure | 25 kgf/cm2 A = 24.517 bar, both ends | LGE, Seongho |
| Feed pressure | shell + 0.3 bar | forced-vortex estimate (assumption) |
| Hole | radial hole on the pin, 2.7 mm, mid-length, theta' = 90 deg | LGE (diameter), assumption (position) |
| Load | F_OS_N of `Sweep, 130 Hz`, rows 0-359 deg, wrap 359 -> 0 | LGE sheet export |
| OS mass | 0.686 kg | LGE sheet |
| Oil | FW68D + R410A saturated at the shell pressure, 70 C | Seongho |
| Model | FBNS, LOCAL_PSAT, SINGLE_PHASE, isothermal, constant mu and rho | prompt |

## Results (campaign `lge-compressor-transient-20261119-02a4b6c`, D:/Seongho/Sim/pancake)

### Part A checks (plain bearing, 180 x 48, W = 1000 N, part_a/part_a_summary.json)

| Check | Result | Pass |
|---|---|---|
| Rotating load at +omega vs static (eps) | 0.84008 vs 0.84224 (-0.26 %) | < 1 % |
| Same, journal torque vs static bearing torque T_j - W e sin(phi) | 0.02972 vs 0.03000 N m (-0.93 %); vs static journal torque -31.7 % | < 1 % |
| Same, attitude from load | +147.31 vs -147.63 deg (mirrored) | within 1 deg |
| IMPLICIT vs EXPLICIT, m = 10 kg, dt = 2 us, 1000 steps | position 0.21 %, torque 0.27 % | < 0.5 % |
| m = 0.686 kg, dt = 20 us (10 x explicit limit), 40 ms | IMPLICIT settles at the static point (eps 3e-5); EXPLICIT diverges at step 364 | yes |

### Part B (GC+VC, 720 x 192, 0.25 deg per step)

Steady (B0):

| Load [N] | eps | h_min [um] | Torque [N m] | p_max [bar] | Cavitated area |
|---|---|---|---|---|---|
| 3,024 | 0.9370 | 1.891 | 0.0726 | 373 | 0.57 |
| 3,450 | 0.9437 | 1.688 | 0.0773 | 439 | 0.57 |
| 3,604 | 0.9458 | 1.625 | 0.0790 | 464 | 0.57 |

Effective p_cav = p_sat(c_feed, 70 C) = 24.517 bar = end pressure; no margin was needed.

Transient: B1 drift <= 6.1e-6 (eps), 0.0014 deg (attitude) per revolution at 1, 0.5 and 0.25 deg
per step. B2 repeats after 3 revolutions. Half dt: cycle-minimum h_min -0.06 %, mean torque +0.03 %.
Film force minus load: 49.6 N (1.4 %) at most, at the 180 deg load jump (OS inertia).

| Last revolution | min | mean | max |
|---|---|---|---|
| h_min [um] | 1.641 (at 105 deg crank) | 1.688 | 1.762 |
| Torque [N m] | 0.0749 | 0.0773 | 0.0784 |
| Friction power [W] | 61.2 | 63.2 | 64.0 |
| p_max [bar] | 382 | 440 | 461 |
| Cavitated area | 0.569 | 0.571 | 0.572 |
| Hole mass flow [g/s] | 1.207 | 1.217 | 1.220 |

Bubble region (void >= 0.5), angles in the direction the OS bore moves relative to the pin: it
starts 20 deg after the thinnest film and ends 263 deg after it, at every crank angle. Measured
from the hole it runs from -140 to +102-103 deg, so it wraps around the hole at other axial
positions, but it never reaches the hole: a liquid ring of 0.35 mm (about 4 cells) remains around
the hole edge at all crank angles, and the void in the cells around the hole is 0. The squeeze
film narrows the h_min swing (1.64-1.76 um) against the steady values at the same loads
(1.62-1.89 um) and delays the thinnest film to about 100 deg after the peak load.

Limits: no case before this checked a time-varying load; the Part A checks and B1a, B2a-c are the
only evidence. GC is instant release and re-dissolution at one uniform pressure. Isothermal, no
asperity contact, no shaft tilt, fictitious frame forces dropped. The 0.35 mm ring is grid-scale
(no grid check). Hole position and angle, feed pressure, no bush, the load direction and the
saturated oil state (a flood-back condition) are assumptions until LGE confirms them.
