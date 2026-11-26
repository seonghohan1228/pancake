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

### Part C: turning inlet, Part B repeated in the OS frame (part_c/summary.json)

Frame: non-rotating, with the pin centre fixed. OS bore at rest, pin turning at +omega, hole and
load turning with the crank (`inlet_rotation_speed` = omega, TABLE load with direction
-90 deg + crank angle). Started from the Part B state at the start of its accepted revolution;
2 revolutions.

| Cycle quantity | Part B (crank frame) | Part C (OS frame) | Difference | Pass (< 2 %) |
|---|---|---|---|---|
| h_min minimum [um] | 1.6406 | 1.6713 | +1.87 % | yes |
| Torque, dissipation / omega [N m] | 0.07738 | 0.07703 | -0.46 % | yes |
| Raw solver torque (rotating surface) [N m] | 0.07734 (OS bore) | 0.04399 (pin) | -43 % | not comparable |
| Mean eps | 0.94375 | 0.94264 | -0.12 % | |

Liquid balance residual <= 8.7e-8 per step with the hole moving half a cell per step.
The h_min difference is a 0.12 % difference in eps (33 nm). In Part C the whole pattern moves
through the mesh at 0.25 deg per step, so its time-discretization error is larger than in Part B,
where the pattern is nearly stationary; a half-dt Part C run would show which way it converges
(not run).

Limits: no case before this checked a time-varying load; the Part A checks and B1a, B2a-c are the
only evidence. GC is instant release and re-dissolution at one uniform pressure. Isothermal, no
asperity contact, no shaft tilt, fictitious frame forces dropped. The 0.35 mm ring is grid-scale
(no grid check). Hole position and angle, feed pressure, no bush, the load direction and the
saturated oil state (a flood-back condition) are assumptions until LGE confirms them.

### VC-only comparison (`vc_steady.py`, campaign `lge-compressor-vc-steady-20261123-00e5d1a`)

The three B0 loads with `SCALAR_PCAV`, p_cav = 500 Pa ("the refrigerant stays dissolved").
(GC+VC - VC)/VC: h_min +6.2 to +6.7 %, torque -10.4 to -11.1 %, p_max -7.6 to -8.0 %, hole flow
+15 to +23 %. VC-only films reach 500 Pa (vapour cavity, 13 % of the area); GC+VC cavity at
shell pressure covers 57 %. Plots: `plot_vc_vs_gcvc.py`.

### Flooded start (`flooded_start.py`, campaign `lge-compressor-flooded-start-20261124-c6715cb`)

Refrigerant-rich feed at a cold start, steady at fixed feed concentration (the film renews in
~15 ms, the sump in seconds). Fluids F0 (c 0.225, 70 C, 5.085 mPa s), F1-F4 (c 0.5/0.6/0.8/0.95,
20 C, 2.62/1.28/0.417/0.159 mPa s), 130 Hz (B0 loads) and 30 Hz ("P2, 30 Hz": 1,432/2,969/4,106 N),
variants N (no contact), C2/C4 (Greenwood-Tripp + Patir-Cheng, composite roughness 0.2/0.4 um,
modulus 107 GPa as the 10/23 study). p_cav = shell pressure (SCALAR_PCAV, no dissolved gas). Feed
overpressure 0.3 bar x (f / 130 Hz)^2. Load balance on the solver's `fluid_force` (film +
contact). Newton on (ln h_geo, attitude). Plots: `plot_flooded_start.py`.

Mean load (h_min in um; Lambda = h_min / sigma; share = load carried by contact):

| f | Fluid | N h_min | C2 h_min, Lambda, share | C4 h_min, Lambda, share |
|---|---|---|---|---|
| 130 Hz | F0 | 1.688 | 1.689, 8.4, 0 % | 1.708, 4.3, 0 % |
| | F1 | 0.940 | 0.948, 4.7, 0 % | 0.986, 2.5, 0.02 % |
| | F2 | 0.489 | 0.511, 2.6, 0.01 % | 0.562, 1.4, 0.7 % |
| | F3 | 0.170 | 0.211, 1.1, 1.4 % | 0.258, 0.65, 4.9 % |
| | F4 | 0.068 | 0.109, 0.55, 4.4 % | 0.139, 0.35, 9.3 % |
| 30 Hz | F0 | 0.516 | 0.537, 2.7, 0 % | 0.588, 1.5, 0.7 % |
| | F1 | 0.279 | 0.312, 1.6, 0.4 % | 0.367, 0.92, 3.1 % |
| | F2 | 0.142 | 0.184, 0.92, 2.2 % | 0.230, 0.57, 6.7 % |
| | F3 | 0.048 | 0.087, 0.44, 6.4 % | 0.111, 0.28, 12.4 % |
| | F4 | 0.019 | 0.046, 0.23, 9.7 % | 0.055, 0.14, 16.2 % |

Checks: B0 reproduction exact; SCALAR_PCAV = LOCAL_PSAT for F0 (difference 4e-14); shell pressure
1.44 MPa changes nothing (4e-14); at F0 the C2/C4 changes (+0.08/+1.2 % h_min) come from the flow
factors, contact is < 1e-4; grid 1440 x 384 at the thinnest case (30 Hz, N, F4, max load, 14 nm):
h_min -3.2 %, torque +1.0 %. Contact enters the moving-OS equation of motion (it is part of
`fluid_force`); a 30 deg table transient (F2, C4) runs at 0.25 and 0.125 deg per step, h_min
differs 0.05 %.

Limits: one homogeneous liquid, isothermal at the feed temperature (no boiling by friction heat),
p_cav = shell pressure, oil density at every c, fixed load direction (the 30 Hz load swings 32 deg),
assumed roughness, modulus and asperity parameters (load share and friction only, no wear).

### Set T: flooded-start transients (`flooded_transient.py`, campaign `lge-compressor-flooded-transient-20261125-55808eb`)

130 Hz load table, moving OS, contact model off, steady mean-load start, 0.25 deg per step,
chained revolutions until h_min repeats (3 each), last revolution again at half dt. Contact
criterion: h_min < 0.5 um (fixed value, not a roughness).

| Fluid | Steady h_min, max / mean / min load (um) | Cycle h_min (um) | Cycle share < 0.5 um | Half-dt change (h_min / torque) |
|---|---|---|---|---|
| F1, c = 0.50 | 0.904 / 0.940 / 1.058 | 0.909-0.995 | 0 | -0.10 % / +0.04 % |
| F2, c = 0.60 | 0.470 / 0.489 / 0.552 | 0.470-0.525 | 72 % (24.5-153, 207-337 deg) | -0.28 % / +0.06 % |

The cycle minimum matches the steady max-load value (F1 +0.5 %, F2 +0.1 %), and the cycle
maximum stays below the steady min-load value: the film follows the slow load decline and
damps the jumps. Steady runs at the max load bracket the transient minimum, so the steady
criterion c(0.5 um) = 0.59-0.62 at 130 Hz holds for the load cycle.

