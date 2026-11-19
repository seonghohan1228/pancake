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

## Results

Filled from the campaign `summary.json` after the runs.
