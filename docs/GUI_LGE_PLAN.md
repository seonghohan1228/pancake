# GUI for LGE bearing engineers: plan (Step 0)

Status: rev. 2 approved 2026-11-26 (D1-D8 as recommended); **rev. 3 for approval** (section 0).
Branch: `lge-gui` (from `main` a11f131), worktree `C:\Dev\pancake\worktrees\lge-gui`.
Prompt: `C:\Workspace\02. 2025 LG전자\07. 최종발표\2026-11-26 pancake prompt (GUI for LGE engineers).md` (rev. 2), plus the rev. 3 update in chat (2026-11-26).

Sections:
0. Rev. 3 changes (supersede the older sections where they differ)
1. Decisions needed before the build
2. Solver base
3. Inventory of the current GUI (Basic / Expanded / Auto / Removed)
4. Model selection (what the GUI sets automatically)
5. Screen map
6. Live feedback and output cadence
7. GUI inputs and `SimulationConfig`
8. Fluid directory
9. Code structure, files to delete or move
10. Acceptance tests: sources and recorded values
11. Draft glossary
12. Open items and risks

---

## 0. Rev. 3 changes (2026-11-26)

State when rev. 3 arrived: D1-D3 done and committed; the case model, fluid library and all screens are written and build; steady acceptance runs pass (section 10). Nothing below is built yet. Measurements made for this revision are marked *measured*.

### 0.1 Decisions for rev. 3

| # | Question | Proposal |
|---|---|---|
| R1 | Default model with dissolved gas on | **FBNS with the cavity pressure = saturation pressure** (the pressure plateau): fast (5-15 core-min steady), the form of the deck's compressor results. "Gas transport (detailed, slower)" in the Advanced view switches to COUPLED. **You decide.** |
| R2 | Thermal model for a thermal *transient* case | The deck's thermal form (finite-volume film energy + 3-D housing conduction) is refused by the solver in transient (*measured*: "FINITE_VOLUME requires steady"). Offer thermal transient with **polynomial cross-film energy + 2-D housing conduction**, the most complete form the solver runs in transient, labelled "simplified thermal model (transient)". Alternative: film-averaged energy with wall heat-transfer coefficients (no housing). |
| R3 | Cavity shear for all FBNS cases | **STRIATED**, fixed, read-only in Advanced (0.6). It is the shear that matches the FBNS cavity (liquid streaks at Θ, massless void). FULL_FILM shears the whole gap and overstates friction loss (solver documentation). *Measured* on the compressor GC+VC case: torque −18.2 %, h_min +0.02 %, cavitated area equal, p_max −0.03 %. The deck's compressor torque, power loss and the tester P5g torque were FULL_FILM, so GUI torque will read about 18 % lower than the deck at 57 % cavitated area. Acceptance test 2 compares h_min and cavitated area (unchanged); torque is checked against a FULL_FILM run of the same config. |
| R4 | Gas phase of a pair | Optional gas block (0.4). The solver has gas data only for R410A (real-fluid EOS), propane and air (ideal gas, fixed molar mass), so the block names one of these. A free ideal gas with its own molar mass needs a new solver key (`gas_molar_mass`); proposed, not built. |
| R5 | Shipped example pair | **No gas block** (FBNS only; COUPLED red with the reason). It stays non-LGE: synthetic VG 68 oil + fictitious refrigerant RX. Alternative: call the example refrigerant propane with synthetic solubility, which enables COUPLED demos but attaches invented solubility to a real refrigerant. |

### 0.2 Run-time estimate (before Run)

Reference core-time C_ref at 720 × 192 cells on afml-r6 (core-time = wall time × processes), from the rev. 3 measurements:

| Form | Low | Typical | High |
|---|---|---|---|
| Isothermal steady, FBNS | 5 core-min | 10 core-min | 15 core-min |
| Isothermal steady, COUPLED | 0.8 core-h | 9.6 core-h | 30 core-h |
| Isothermal transient, FBNS, per revolution | C(n) × 0.7 | C(n) = 0.5 core-h × (n/360)^1.44 | C(n) × 1.4 |
| Thermal steady, FBNS (pure oil) | 0.8 core-h | 2.5 core-h | 8 core-h |
| Thermal steady, COUPLED | 10 core-h | 22 core-h | 50 core-h |
| Thermal transient (simplified, R2), per revolution | 1.0 × isothermal | 1.5 × isothermal | 2.0 × isothermal |

- **COUPLED isothermal steady** was given as 2-74 min wall (median 24 min). I assumed 24 processes (afml-r6), which gives the core-times above. **Please confirm the process count.**
- **Transient:** n is the number of steps per revolution. The two measured points (3.7 core-h at 1440 steps, 0.5 core-h at 360) fix the exponent 1.44. Total = per revolution × the revolutions set (2-3 to settle; the GUI says so).
- **Thermal transient factor:** *measured* on 180 × 48 cells, 4 processes, time per step from 40 minus 20 steps:

  | Form | Time per step | Relative to isothermal |
  |---|---|---|
  | Isothermal | 0.18 s | 1.0 |
  | Film-averaged energy | 0.26 s | 1.5 |
  | Polynomial cross-film + 2-D housing | 0.18 s | 1.0 |

  Scaled to 720 × 192 at 1440 steps per revolution, the isothermal case gives 4.6 core-h per revolution, against 3.7 measured, so the cell scaling holds within 25 %. The user's 25-70 core-h per revolution referred to the deck thermal form, which cannot run transient.
- **Scaling:**
  - **Cells:** linear in n_θ × n_z.
  - **Processes:** wall = core-time / (p · η(p)) with η(p) = 1 / (1 + 0.02 (p − 1)), about 68 % at 24. This is an assumption, checked by acceptance test 10.
  - **Machine:** factor 1. The compressor B0 run on 8 local processes took 13.5 core-min, inside the afml-r6 FBNS range.
- **Display:** "about 25 min on 4 cores (range 12-40 min)". The headline is the typical value at the chosen core count.

### 0.3 Thermal transient (an option, not blocked)

- *Measured:* the solver runs energy + transient together in two forms: film-averaged energy, and polynomial cross-film + 2-D housing.
  - Time per step is in 0.2.
  - The deck thermal form is steady-only: it is red with the solver's reason if forced.
- **Before Run, a confirmation dialog shows:**
  - the time per revolution and the total for the chosen revolutions and cores;
  - the line: "The housing warms over seconds to minutes; a few revolutions (milliseconds each) do not reach the thermal periodic state.";
  - the faster two-step route as a second button: (1) thermal steady at the mean load of the table, with the deck thermal form; (2) isothermal transient with the film temperature set to the mean film temperature of step 1. The GUI writes both cases and queues them.
- Run starts only after the user confirms.

### 0.4 Fluid file: optional gas block

```
# Gas phase (optional). Needed for gas transport (COUPLED).
gas_model = R410A | PROPANE | AIR     # R410A: real-fluid EOS; PROPANE, AIR: ideal gas
gas_viscosity_Pa_s = ...
gas_heat_capacity_cv_J_kgK = ...      # optional; ideal gases only (overrides the default)
gas_thermal_conductivity_W_mK = ...
dissolved_heat_capacity_J_kgK = ...   # thermal + gas release
desorption_energy_J_kg = ...          # thermal + gas release
```

- `solver_species` is replaced by `gas_model`.
- A pair without the block runs FBNS only. Gas transport (and any thermal case with refrigerant in the oil, which needs COUPLED) turns the Solution node red: "This pair has no gas-phase data; gas transport is not available."
- The FBNS plateau needs no gas data: it uses only the saturation pressure from the table.

### 0.5 Showing FBNS results honestly

- **FBNS:** the cavity is liquid at fraction Θ plus a massless bubble that carries no shear, at the cavity pressure. No gas mass is tracked, and GC and VC differ only in the cavity pressure.
  - Field maps offer "Void fraction 1 − Θ" (outline at 0.5). If density is shown, it is "Effective density Θ · ρ_liquid".
  - `alpha_gas`, `dissolved_gas` and `rho` are hidden for FBNS. If a user opens a file view of them, they are labelled "not computed in this model".
  - The results card says "Cavitated area (void fraction > 0)".
- **COUPLED:** free-gas volume fraction (`alpha_gas`), dissolved content (`dissolved_gas`), gas mass transfer (`gas_mass_transfer`), mixture density (`rho`) and mixture viscosity (`mu`). Tooltip on the viscosity: "Void-weighted (Dukler) mixture viscosity; not validated for oil-refrigerant mixtures."
- **Test 11:** an FBNS result never shows a nonzero gas field. This is checked on the compressor case G_eps0.9300 (`D:\Seongho\Sim\pancake\lge-gc-sign-eccentricity-20261125-16b9796`) and on a GUI-run FBNS case.

### 0.6 Fixed choices shown read-only in Advanced

- **Cavity shear:** "No shear from the bubble; only the liquid shears" (solver key STRIATED) for every FBNS case (R3). COUPLED keeps its own full-film treatment (the solver accepts STRIATED only with JFO/FBNS).
- **Solver form, gas model, thermal form:** read-only lines (section 4).

### 0.7 Layout rule: no overlapping text

- Every label is measured (`CalcTextSize`) before it is placed. Long labels wrap or are cut with "…" and a tooltip.
- Plot annotations that would collide are stacked, as for the h_max/h_min/feed marks in the live views.
- **Check:** the screenshot series at 1600 × 1000 and at 1280 × 800. Each screenshot is inspected for overlapping or clipped text before acceptance (test 9).

### 0.8 New acceptance tests

| # | Test | Pass |
|---|---|---|
| 10 | Run-time estimate vs measured, this machine: one isothermal steady case (compressor GC+VC) and one thermal steady case | estimate within a factor of 2 of the measured wall time |
| 11 | FBNS result never shows a nonzero gas field | `alpha_gas`, `dissolved_gas` hidden or labelled; void fraction from Θ only |
| 12 | Thermal transient small test | *measured* (0.2): the two simplified forms run; time per step reported; the deck thermal form is red with the solver's reason |

### 0.9 Other findings while preparing rev. 3

- **X4 config:** it says `cavitated_film_model = HOMOGENEOUS`. The parser silently ignores unknown values of this key, so the run used FULL_FILM. The same silent fallback affects other enum keys. Proposal: make unknown enum values a parse error (solver change, small). Not built.
- **Steady result snapshots:** they all share step 0, so the Results view did not reload between snapshots or between two steady records. To be fixed in `poll_loader` by comparing the step file path, not the step number.

---

## 1. Decisions needed before the build

| # | Question | Recommendation |
|---|---|---|
| D1 | Solver base | Merge `lge-compressor-transient` into `lge-gui`. It contains the `src/` work of `lge-onset-chart` and `lge-mix-c274`. The `src/` merge is clean; only `CHANGES.md` and `PLAN.md` conflict. See section 2. |
| D2 | Load balance with dissolved gas on. `STEADY_EQUILIBRIUM` rejects FBNS with dissolved gas (`config.hpp:1927-1938`, branch version). The deck's GC+VC steady runs (compressor B0 mean, tester P5g) are `STATIC` runs, load-balanced by Python drivers (`02h part_b.py balance`, `02g run_fbns.py balance`). The GUI cannot reproduce them with `pancake.exe -c` alone. | **(a)** A narrow solver change: allow `STEADY_EQUILIBRIUM` for FBNS with `LOCAL_PSAT` (SINGLE_PHASE, and TWO_PHASE with `steady_gas_model = EQUILIBRIUM`). First find out why the guard exists. Accept the change only if it reproduces the driver results within 0.5 % (B0 mean: h_min 1.6877 µm, 57.09 %; P5g: h_min 3.242 µm). **(b)** Fallback: the GUI runs the Newton load balance itself over `STATIC` runs and saves the balanced config. The terminal then reproduces only the final point. |
| D3 | `to_config_text()` drops keys. It never writes `coupling_max_iters`, `coupling_tolerance`, `coupling_relaxation*`, `coupling_keep_best_iterate`, `steady_gas_model`, `allow_frozen_gas_steady`, `cavitated_film_model`, `consistent_boundary_flux`, `void_coupling_relaxation`, `log_coupling_history` or `bc_z_*_theta`. The current GUI saves before every run (`gui_app.cpp:399-431`), so a loaded config silently loses them. The bug is on `main` too. | Fix the writer (serialization only, no numerics) and add a test: parse, write, parse, then compare every known key. Separately, the GUI saves loaded configs by patching the original text (section 7.3). |
| D4 | Speed as a function of time. The solver has a load table but no speed table (only `omega` + `omega_ramp_time`). | Constant speed in transient runs for now; the load table is supported. A speed table would be a solver change, proposed separately if LGE needs it. |
| D5 | Form of "dissolved gas on" in the FBNS path. The compressor deck run used SINGLE_PHASE + `LOCAL_PSAT` (pressure plateau at p_sat, no gas-fraction field). P5g (deck p45) used TWO_PHASE + `LOCAL_PSAT` + `steady_gas_model = EQUILIBRIUM`. | SINGLE_PHASE + `LOCAL_PSAT`, the compressor form (`MODEL_SELECTION.md` calls the FBNS free-gas extension lagged). The bubble region is drawn from 1 − film content. Test 2 then also runs the tester case in this form and reports the difference from P5g. |
| D6 | Onset chart Φ(S). The 7 gas-free points come from tester runs G1–G8 (`lge-onset-chart-20261116`, FW68D + R410A, tester geometry). No code in the repository computes Φ yet. | Ship the curve as a nondimensional data file with its source and limits stated: "tester bearing, L/D = 0.93, groove feed; approximate for other bearings". No LGE property data goes with it. After a run, the GUI also gives the direct check (film p_min against p_sat), which holds for any bearing. |
| D7 | Example fluid. | A synthetic pair labelled "synthetic, example only, not for design". Viscosity follows ISO VG 68 (68 mm²/s at 40 °C) through a Walther fit, with log-mixing thinning by a fictitious refrigerant. Solubility comes from an ideal-solution law with a fictitious vapour-pressure curve. Published data can replace it later. |
| D8 | Default mesh level. | Medium (360 × 96) for new cases. The acceptance tests use Fine (720 × 192), the deck grid. |

---

## 2. Solver base

All three branches fork from `acc20fd` (2026-10-19). `main` has gained 5 commits since then; one touches `src/` (`1e034f6 fix(io)`).

| Branch | `src/` added relative to `main` | Needed by the GUI |
|---|---|---|
| `lge-mix-c274` (6 `src/` commits) | COUPLED inner tolerance (`coupled_nonlinear_tolerance`), relaxed outer coupling, θ ≤ 1 and p ≥ p_v admissibility fixes, mixed lubrication (Greenwood-Tripp, Patir-Cheng; `mixed_lubrication.hpp`), resolved cross-film updates. 17 files, +942/−108. | Yes: the X-runs and the onset chart ran on it. Mixed lubrication comes along but stays hidden in the GUI. |
| `lge-onset-chart` | Same `src/` as c274 (`3e3d1cc` changes only drivers). | Yes, through c274. |
| `lge-compressor-transient` (8 `src/` commits, contains both above) | Adds `load_model` (CONSTANT/ROTATING/TABLE), load table by crank angle or time, implicit motion coupling, `inlet_rotation_speed`. 19 files, +1719/−147. | Yes: transient load table, compressor runs. |

Proposal: `git merge lge-compressor-transient` into `lge-gui`, then resolve `CHANGES.md` and `PLAN.md` by keeping both sides. `git merge-tree` shows no `src/` conflict. Build and run `ctest` before any GUI work.

---

## 3. Inventory of the current GUI

Legend: **B** Basic (always shown) · **E** Expanded ("Advanced" toggle in the panel) · **A** Auto (set by a model choice, read-only in Expanded) · **R** Removed from GUI (the key stays valid in files; a loaded value is kept and listed as "not shown").

### 3.1 Parameters (`schema.cpp`: 62 double, 4 int, 4 bool specs, plus combos in `panels_setup.cpp`)

| Group | Key | Now | New | New label / note |
|---|---|---|---|---|
| Geometry | `R` | common | B | Bearing diameter D (= 2R) |
| | `c` | common | B | Radial clearance |
| | `L` | common | B | Bearing length |
| | `e`, `attitude_angle_deg` | common | E / A | Fixed shaft position mode only; otherwise the initial guess (A) |
| | `load_angle_deg` | common | R | Legacy; load direction comes from the load vector |
| | `tilt_slope_x`, `tilt_slope_y` | adv. | E | Shown as tilt angle γ and tilt direction θ_γ (deck p18) |
| Operating | `omega` | common | B | Speed (rpm / Hz / rad/s) |
| | `omega_ramp_time` | adv. | R | |
| Lubricant | `mu`, `rho` | common | A | From the fluid pair at the film temperature |
| | `p_cav` | common | A | Oil vapour pressure, from the fluid pair |
| | `bulk_modulus` | adv. | R | Elrod only |
| Time | `end_t` | common | E | Number of revolutions (transient) |
| | `dt` | common | E | Steps per revolution |
| | `write_interval` | common | E | Field outputs per revolution |
| Cavitation | `outer_tol`, `theta_min`, `max_outer_iters` | adv. | R | Elrod only |
| Boundaries | `bc_z_south_val`, `bc_z_north_val` | common | A | = ambient pressure; editable in the bearing-ends zone (E) |
| Thermal | `temperature_initial`, `temperature_reference` | common/adv. | A | = ambient temperature |
| | `journal_wall_temperature`, `bearing_wall_temperature` | adv. | E | Shaft / bearing surface zones |
| | `journal_heat_transfer`, `bearing_heat_transfer` | adv. | E | Shaft / bearing surface zones |
| | `rho_cp`, `thermal_conductivity` | adv. | A | From the fluid pair |
| | `gas_heat_capacity_cv`, `dissolved_heat_capacity`, `desorption_internal_energy`, `gas_thermal_conductivity`, `oil_vapor_temperature_coeff` | adv. | R | COUPLED caloric internals |
| Fluid models | `dissolved_gas_initial` | adv. | B | Refrigerant content (wt %), or "saturated at ambient pressure" |
| | `dissolved_gas_max`, Henry (2), `dissolved_gas_liquid_density`, `gas_mass_transfer_rate`, `gas_absorption_rate`, `dissolved_gas_diffusivity`, `solution_density_gas_coeff`, `solution_viscosity_gas_coeff`, `viscosity_temperature_coeff`, `viscosity_pressure_coeff`, `viscosity_poly_k0..k2` | adv. | R | Replaced by the fluid-pair tables |
| | `density_temperature_intercept/slope` | common | A | From the fluid pair |
| | `gas_alpha_max`, `gas_pressure_floor`, `mu_gas` | adv. | R | Gas equation-of-state internals |
| | `property_reference_pressure/temperature` | adv. | A | |
| Motion | `external_load_x/y` | common | B | Load magnitude and direction |
| | `external_load_z` | adv. | R | |
| | `bearing_mass`, `bearing_stiffness_x`, `bearing_damping_x` | adv. | R | Raw support stiffness and damping |
| | `min_film_thickness`, `stop_on_nonpositive_film` | adv. | R | Guards; kept at defaults |
| Numerics | `linear_rtol`, `coupled_nonlinear_tolerance` | adv. | R / E | `coupled_nonlinear_tolerance` in Expanded convergence (COUPLED only) |
| | `diagnostics_interval` | adv. | A | = 1 |
| | `log_outer_iters` | adv. | R | |
| Output | `output_write_3d`, `output_write_flat` | adv. | A | Both on (flat for the GUI, 3D for ParaView) |
| Combos | `solution_mode` | common | B | Analysis type: steady / transient |
| | `cavitation_model` (5 options) | common | A | FBNS or COUPLED (section 4). Full Sommerfeld, Gumbel, JFO: R |
| | `cavitation_threshold` | common | A | From "dissolved gas on/off" |
| | `coupled_phase_equilibrium` | common | A | |
| | `gas_pressure_coupling` | common | R | Always NONE |
| | `bc_z_*_type` | adv. | R | Always fixed pressure (open ends) |
| | `bc_z_*_thermal` | adv. | E | Bearing-ends zone |
| | `temperature_model` | adv. | B | Thermal analysis off / on |
| | `fluid_property_model`, `dissolved_gas_species`, `solubility_model`, `density_model`, `liquid_viscosity_model`, `mixture_viscosity_model` | adv. | A | From the fluid pair + dissolved gas setting |
| | `*_table_file`, `solution_table_file`, `viscosity_composition_file`, `solution_table_clamp`, `kinematic_viscosity_log_interpolation` | adv. | A | From the fluid pair |
| | `motion_model` | adv. | A | From "shaft position" and analysis type |
| | `output_verbosity` | adv. | A | DETAILED |
| | `output_dir`, `filename_prefix` | adv. | A | `<case>_results/`, `solution` |
| | `motion_time_method` | adv. | A | Implicit Euler |
| | `theta_convection_scheme` | adv. | R | |
| Inlets | inlet editor (type, θ, z, size, p_supply, t_supply) | common | B / E | One feed in Basic (type, angle, size); supply pressure and temperature in the feed zone (E). Further feeds in a config are kept and listed |

### 3.2 Keys the GUI does not show today but the new GUI needs

| Key(s) | New | Why |
|---|---|---|
| `motion_model = STEADY_EQUILIBRIUM`, `steady_force_tolerance` | A / E | Load-balanced steady case (the current combo lacks this option) |
| `coupling_max_iters`, `coupling_tolerance` | E | Outer iteration limit and tolerance |
| `fbns_tol`, `fbns_max_iters` | E | Inner solution tolerance (FBNS) |
| `load_model`, `load_table_file`, `load_table_abscissa`, `load_table_format`, `load_frequency`, `load_phase_offset_deg`, `load_table_scale` | B (transient) | Load table |
| `motion_coupling`, `motion_coupling_tolerance`, `motion_coupling_max_iters` | A / E | Transient load-to-position coupling |
| `cross_film_thermal_model`, `cross_film_layers` | A / E | Thermal on = FINITE_VOLUME; layers in Expanded mesh |
| `bush_thermal_model`, `bush_outer_radius`, `bush_thermal_conductivity`, `bush_outer_heat_transfer`, `bush_outer_temperature`, `bush_end_*`, `bush_radial_cells`, `bush_axial_cells` | B / E | Housing (thermal on): outer radius in Basic, the rest in the housing zone |
| `shaft_ambient_*`, `thermal_ambient_temperature` | E / A | Shaft zone; ambient = operating ambient |
| `steady_gas_model` | A | |
| `steady_output_interval` | A | Live field views in steady runs (section 6) |
| `inlet_rotation_speed` | E | "Feed hole turns with the shaft" (transient FBNS only) |
| `roughness_rq_shaft`, `roughness_rq_bush` | B | Surface roughness for Λ. Checked: read only when contact or flow factors are on (`config.hpp:2355`), so they do not change results |

Stay removed (as the prompt lists): contact / mixed lubrication (`asperity_*`, `flow_factor_model`, `contact_modulus`, `boundary_friction_coefficient`), PETSc linear options, convection schemes, gas EOS internals, `output_fields`, debug logging. **Note for the user guide:** the flooded-start cases (deck p48) reach Λ < 3. The GUI shows Λ in red below 3 and says that contact is not modelled. It does not add a contact model.

### 3.3 Panels and selectors today

| Item | Now | New |
|---|---|---|
| Case Setup window: run bar, MPI count, validation line | one window | Run bar stays (top), cores always visible (default 1, per case). Validation moves to the tree lights |
| Search filter, Advanced checkbox, Presets combo | toolbar | R (an "Advanced" toggle per panel replaces them; example cases replace presets) |
| 12 collapsible parameter groups | schema-generated | R: replaced by seven workflow panels |
| Species presets (`apply_dissolved_gas_species_preset`) | code | R: the fluid directory replaces them |
| Inlet editor | Boundaries group | Feed in Geometry; feed values in the feed zone |
| Results: Field map (heatmap, cavitation and inlet overlays, colormap, range, CSV/PNG) | tab | Kept and rebuilt: angle reference, 1:1, α = 0.5 contour, h_min line, four fields |
| Results: Profiles (p, h at a z slice) | tab | Kept; adds the polar plot |
| Results: Summary | tab | Key-results card (section 5.6) |
| Monitors (h_min, T, |F|, torque, cavitated %, residual) | window | Kept; adds force residual from `equilibrium_history.csv` and the transient histories from `motion.csv` |
| Solver Log | window | Kept, under View |
| Menus File / Run / Results / Tools / Help | menu bar | Kept; File gains Duplicate case; Tools > Settings gains the fluid directory |

---

## 4. Model selection (Auto)

The engineer picks three things: dissolved gas on/off, thermal on/off, steady/transient. Expanded adds one choice: the gas model (plateau or transport). The GUI derives the rest, following `docs/MODEL_SELECTION.md` and the deck runs.

| Engineer choice | `cavitation_model` | Key settings | Deck runs of this form |
|---|---|---|---|
| Dissolved gas **off** (VC only) | FBNS | `fluid_property_model = SINGLE_PHASE`, `cavitation_threshold = SCALAR_PCAV`, `p_cav` = oil vapour pressure | compressor VC steady (R-mean), O1–O4 |
| Dissolved gas **on** (GC+VC), default plateau model | FBNS | SINGLE_PHASE, `LOCAL_PSAT`, `dissolved_gas_initial` = content, `steady_gas_model = FROZEN` (see D5) | compressor B0 mean, rev03 transient |
| Dissolved gas on, Expanded "gas transport" | COUPLED | TWO_PHASE, `coupled_phase_equilibrium = true`, `SCALAR_PCAV`, DUKLER mixture viscosity | X1–X5, onset chart, optical rig |
| Thermal on | (either) | `ENERGY_EQUATION`, `cross_film_thermal_model = FINITE_VOLUME`, `bush_thermal_model = CONDUCTION_3D` | O/X runs |
| Steady + load | | `STEADY_STATE`, `motion_model = STEADY_EQUILIBRIUM` (needs D2 for gas on) | O/X runs |
| Steady, fixed shaft position (Expanded) | | `STATIC`, `e`, `attitude_angle_deg` | optical rig |
| Transient | | `TRANSIENT`, `MOVING_BEARING`, `motion_coupling = IMPLICIT`, `EULER_IMPLICIT`; load table → `load_model = TABLE` | rev03 |

**Found during the build (2026-11-26):** the solver allows the film energy model in steady runs only, and FBNS film energy only for pure oil (`FLUX_BALANCED` requires content 0, as in the O-runs). A thermal case with refrigerant in the oil therefore runs COUPLED, gas off as SINGLE_PHASE and gas on as equilibrium TWO_PHASE, as in the X-runs. Gas release in COUPLED needs the free-gas viscosity, and with heat also the caloric data (dissolved heat capacity, desorption energy). These are optional fields of a fluid pair; a pair without them turns the case red with that reason. Thermal + transient is red. The cavitated area is the fraction of cells with film content below 1: it reproduces the deck values (57.09 %, 11.83 %); the void >= 0.5 measure gives 52.8 % and 7.4 %.

Unsupported combinations turn the Solution node red with one sentence. Examples: COUPLED + transient + turning feed; COUPLED with a fluid pair whose refrigerant has no equation of state in the solver (only R410A real-fluid, air and propane ideal gas exist).

Expanded shows the chosen form read-only, e.g. "Solver form: FBNS (complementarity), steady load balance (Newton), isothermal".

---

## 5. Screen map

### 5.1 Main window

```
+--------------------------------------------------------------------------------------+
| File  Run  Results  View  Tools  Help                                                |
+--------------------------------------------------------------------------------------+
| [> Run F5] [x Cancel]  CPU cores [ 1 v]   Case: compressor_example.txt *  ...status  |
+-------------------------+------------------------------------------------------------+
| WORKFLOW                |  PANEL (selected node)                [ ] Advanced         |
|  (G) Geometry           |  Geometry                                                  |
|  (G) Fluid              |   Bearing diameter D   [ 21.0 ] [mm v]                     |
|  (Y) Operating cond.  ! |   Bearing length L     [ 20.0 ] [mm v]                     |
|  (G) Boundary cond.     |   Radial clearance c   [ 30.0 ] [um v]                     |
|  (G) Mesh               |   Oil feed             (o) Hole  ( ) Groove                |
|  (G) Solution           |   Feed angle           [ 90.0 ] [deg v]  (from fixed frame)|
|  ( ) Results            |   Hole diameter        [ 2.7  ] [mm v]                     |
|                         |   -- Advanced --  tilt angle, tilt direction               |
| ESTIMATES (live)        +------------------------------------------------------------+
|  S        = 0.016       |  LIVE VIEW   [Geometry | Mesh]   Angle zero [h_max v]      |
|  p_sat    = 24.5 bar    |        feed o                                              |
|  margin   = 0.0 bar     |      .-------.         clearance x 50 (exaggerated)        |
|  margin/P = 0.00        |     /  .---.  \        rotation (arrow), load W (arrow)    |
|  GC: expected (onset)   |    |  ( o-->)  |       theta = 0 at h_max, theta_in = 0    |
|  run time ~ 3 min       |     \  '---'  /        at feed                             |
|  [onset chart mini]     |      '-------'                                             |
+-------------------------+------------------------------------------------------------+
| (Y) Operating conditions: enter the load.                                            |
+--------------------------------------------------------------------------------------+
```

Lights: green complete; yellow input needed; red conflict, the case will not run. The reason appears as one sentence in the panel, in the bottom bar and as a tooltip on the tree node. Run is enabled only when every node is green. Red uses the solver's `SimulationConfig::validate()` on the mapped config; yellow is a GUI check (missing input).

### 5.2 Panels (one per tree node)

```
Fluid                                         Operating conditions
 Oil-refrigerant pair [Example VG68 + RX v]    Analysis type   (o) Steady  ( ) Transient
   [Open fluid library...]                     Speed           [ 7800 ] [rpm v]
 Refrigerant content  (o) Saturated at         Load            [ 3450 ] [N v]  dir [270] deg
                          ambient pressure     Ambient pressure[ 24.5 ] [bar v]
                      ( ) Given  [22.5] [wt%]  Ambient temp.   [ 70.0 ] [degC v]
   -> 22.5 wt% at 24.5 bar, 70 degC            (transient) Load table [os_load.csv] [Import CSV...]
 Surface roughness Rq shaft [0.2] bush [0.2] um   crank-angle / time, magnitude+direction / x,y,
 Within table range: yes (green)                  frequency [130] Hz, revolutions [3]
                                                  [plot of load vs crank angle]

Boundary conditions (zones; double-click opens the pop-up, the zone lights up in the views)
 Bearing ends ............ pressure = ambient (24.5 bar); inflow temperature
 Oil feed ................ supply pressure = ambient + [0.3] bar; supply temperature
 Shaft surface ........... (thermal) heat transfer to ambient
 Bearing surface ......... (thermal) housing conduction
 Housing outer surface ... (thermal) heat transfer coeff., outer temperature

Mesh                                          Solution
 Level (o) Coarse 180 x 48                     Dissolved gas   [x] on (GC+VC)  [ ] off (VC only)
       ( ) Medium 360 x 96                     Thermal         [ ] on (film, housing, shaft)
       ( ) Fine   720 x 192                    Comparison mode [ ] run gas on and off together
 cells: 2.0 deg x 0.10 mm; hole spans 13 cells -- Advanced --
 -- Advanced -- n_theta, n_z, cross-film layers  Convergence (default, controls):
                                                  load balance  2e-3   |F+W|/|W|
                                                  outer         1e-6   pressure/temperature change
                                                  inner         1e-10  film equations (FBNS)
                                                  max outer it. 200
                                                 Time stepping: steps/rev, revolutions, outputs/rev
                                                 Solver form (read-only): FBNS, Newton load balance
```

### 5.3 Fluid library dialog

```
+---------------------------------------------------------------------------------+
| Fluid directory: C:\...\pancake\fluids   [Change...]                            |
+--------------------------+------------------------------------------------------+
| Pairs                    | Header: oil, refrigerant, source, valid range,       |
|  example_vg68_rx  (syn.) |         density law, oil vapour pressure             |
|  my_pair                 | Table (c, T, p_sat, nu)  [paste from Excel]          |
| [New] [Copy] [Delete]    |  c [wt%] | T [degC] | p_sat [bar] | nu [mm2/s]        |
| [Open folder]            |  ...                                                 |
|                          | Plots: nu vs T per content | p_sat vs content at T_amb|
|                          | Warnings: "content 95 wt% above table (max 90)"      |
|                          |                       [Save] [Revert]                |
+--------------------------+------------------------------------------------------+
```

### 5.4 Run monitor (while running)

```
Stage: Load balance, iteration 12 of max 200 (force residual 3e-2)   [=====      ] 40%
 h_min [um]   |----\___         force residual (log) |\__
 torque [N m] |  __/----        cavitated area [%]   |  _/--
 Field view refreshed every 10 outer iterations (steady) / every output step (transient)
 Command: mpiexec -n 4 pancake.exe -c "C:\cases\comp.txt"           [Copy command]
```

### 5.5 Batch / sweep

```
Sweep  input [Load v]  from [2000] to [4000] N  steps [5]   or list [....]  or [Add case files...]
       cores per case [1]   [Write batch script]   [Run batch]
+----+--------+--------+---------+-----------+----------+--------+
| #  | load N | status | h_min um| torque Nm | cav. %   | onset  |
| 1  | 2000   | done   | 2.41    | 0.071     | 52       | GC     |
| 2  | 2500   | running| ...     |           |          |        |
+----+--------+--------+---------+-----------+----------+--------+
[h_min vs load, grows as cases finish]   [onset chart: S vs margin/P, case points]
```

### 5.6 Results (steady)

```
Key results  (step final, converged: yes)            Angle zero [h_max v]  Units [...]
 h_min 1.69 um | eps 0.944 | phi 289.6 deg | torque 0.0773 N m | power 63.2 W
 p_max 439 bar | cavitated area 57 % | Lambda 6.0 | onset: GC expected (p_min < p_sat)
 shear share of film force: 3 %          (yellow/red mark if a table range was left
                                          or the run did not converge)
[Field map: p | h | gas fraction | T]  x = angle (zero: h_max), y = z/L, 1:1
   alpha = 0.5 contour, feed outline, dashed h_min line, labels above the map
[Polar: p at mid-plane, ambient and p_sat circles, over/under-ambient fills]   [Profiles p, h vs angle]
[Export: PNG | CSV | Case report]
```

Transient adds time histories with a cursor (h_min, eps, torque, cavitated area). The field maps follow the cursor, and a play/pause control animates them. The moving zeros (h_max, h_min) follow the instant shown.

Comparison mode: one table (gas on, gas off, Δ %) and overlaid profiles. VC is red, GC+VC green and highlights purple, as in the deck.

---

## 6. Live feedback and output cadence

No solver change is needed for live views:

| Run | Monitors | Field views |
|---|---|---|
| Steady, load balance | `equilibrium_history.csv` (force residual, position, R_p, R_T), `diagnostics.csv` | The GUI sets `steady_output_interval = 10`; the solver writes `iterations/NNNNNN/` snapshots. Views refresh every 10 outer iterations. |
| Steady, fixed position | `diagnostics.csv`, `Coupling outer=` log lines | Same, via `steady_output_interval` |
| Transient | `motion.csv` (every step: h_min, ε, φ, torque, p_max, cavitated fraction, load) | VTK at every `write_interval` (default 36 per revolution) |

The stage text comes from the log ("Load balance, iteration k", "Revolution 2 of 3, crank 140°").

Instant estimates (no run):
- S = μ N L D (R/c)² / W, N in rev/s, μ from the fluid table at the ambient temperature and the content.
- p_sat(T_amb, content) by the solver's own `saturation_pressure` (header-only, linked into the GUI).
- Pressure margin p_amb − p_sat, its ratio to P = W/(LD), the chart point and the verdict: "gaseous cavitation expected / borderline (within ±20 % of the curve) / not expected".
- **Colours mean design risk, not model state** (Seongho, 2026-11-26: this is a design tool, and GC expected is not a good outcome). Gas release expected = red, borderline = yellow, not expected = green. The same rule applies to the results card: Λ < 3 red, a table range left yellow, not converged red. Model colours (VC red, GC+VC green, as in the deck) are used only for curves that compare the two models, never for a verdict.
- The run-time estimate is cells × typical iterations × time per cell-iteration. It is measured on this machine and scaled by cores; it is labelled "rough".

After the run, the onset check is direct: GC is expected if the gas-off film p_min < p_sat. In a gas-on run the plateau at p_sat shows it as a cavitated region at p_sat.

---

## 7. GUI inputs and `SimulationConfig`

### 7.1 Case model

`CaseModel` (new, `case_model.hpp`) holds only engineer-level inputs, in SI units. Panels edit only this.

| CaseModel field | Config key(s) | Mapping |
|---|---|---|
| `diameter`, `length`, `radial_clearance` | `R`, `L`, `c` | R = D/2 |
| `feed.type`, `feed.angle_deg`, `feed.size`, `feed.z` | `inlet_circular` / `inlet_groove` (first entry) | Hole: size = diameter, radius = size/2, z = L/2 by default |
| `feed.supply_overpressure`, `feed.supply_temperature` | 3rd/4th inlet tokens | p_supply = p_amb + Δp; T_supply default = T_amb |
| `feed.turns_with_shaft` | `inlet_rotation_speed` = ω | Transient FBNS only |
| `tilt_angle`, `tilt_direction_deg` | `tilt_slope_x/y` | slope = γ (cos θ_γ, sin θ_γ) in the solver frame; checked against `film_thickness.cpp:21` |
| `speed` | `omega` | |
| `load`, `load_direction_deg` | `external_load_x/y` | Polar to Cartesian in the solver frame (load "down" = deck −Y) |
| `load_table` (file, abscissa, format, frequency, revolutions) | `load_model = TABLE`, `load_table_*`, `load_frequency`, `end_t` | end_t = revolutions / f |
| `ambient_pressure` | `bc_z_south_val`, `bc_z_north_val` | Both ends = p_amb unless edited in the zone |
| `ambient_temperature` | `temperature_initial`, `temperature_reference`, `thermal_ambient_temperature`, inflow T | |
| `fluid_pair` | `solution_table_file`, `solubility_model = CURVE_TABLE`, `liquid_viscosity_model`, `density_model` + coefficients, `p_cav`, `rho_cp`, `thermal_conductivity`, `dissolved_gas_species` | From `fluids/<pair>/pair.txt` (section 8). Isothermal FBNS: `liquid_viscosity_model = CONSTANT`, `mu` from the table (as in the deck runs; the exact rule is checked against `part_b.py`) |
| `content` / `saturated` | `dissolved_gas_initial` | Saturated: c_sat(p_amb, T_amb) from the table |
| `roughness_shaft`, `roughness_bush` | `roughness_rq_shaft/bush` | Inert unless contact is on |
| `dissolved_gas_on`, `gas_model`, `thermal_on`, `analysis_type`, `shaft_position_mode` | section 4 table | |
| `housing_outer_diameter` and housing zone values | `bush_*` | Thermal on only |
| `mesh_level` / `n_theta`, `n_z`, `cross_film_layers` | `n_theta_global`, `n_z_global`, `cross_film_layers` | Coarse 180×48, Medium 360×96, Fine 720×192 |
| `steps_per_revolution`, `outputs_per_revolution` | `dt`, `write_interval` | dt = 2π/(ω n) |
| convergence values | `steady_force_tolerance`, `coupling_tolerance`, `coupling_max_iters`, `fbns_tol`, `fbns_max_iters`, `coupled_nonlinear_tolerance`, `motion_coupling_tolerance`, `motion_coupling_max_iters` | Direct |
| `cores` | not in config: command line `mpiexec -n N` | Stored as the comment `#@gui cores = N`; default 1 |

GUI-only view settings (angle zero, units, mesh level label, pair name, "saturated") go into `#@gui key = value` comment lines. The solver ignores comments, so the GUI needs no case format of its own.

### 7.2 Two tested mappings

- `to_config(case, base) -> SimulationConfig`: starts from `base` (the loaded config, or defaults for a new case) and writes only the keys in 7.1.
- `from_config(config) -> {case, hidden}`: `hidden` lists the keys whose value differs from what `to_config` would write. The GUI says "This case uses N settings the GUI does not show; they are kept on save" and lists them.

Unit test: for each acceptance config, `to_config(from_config(cfg), cfg)` equals `cfg` in every key in `known_config_keys`.

### 7.3 Saving

- Loaded file: patch the original text. Lines of unchanged keys stay byte-identical (CRLF, comments, blank lines, afml-r6 paths with Korean characters). Changed keys are rewritten in place, and new keys are appended under `# GUI`.
- New case: `to_config_text()` (after the D3 fix).

### 7.4 Terminal reproducibility

- The results window shows the exact command: `pancake.exe -c "<file>"`, or `mpiexec -n N pancake.exe -c "<file>"`, with a Copy button.
- Sweeps write `<sweep>/case_001.txt ...` plus `run_all.ps1` and `run_all.bat`.

---

## 8. Fluid directory

`fluids/` next to `pancake_gui.exe`; the location can be changed in Settings. One folder per pair:

```
fluids/
  README.md                      format, units, how to add a pair
  example_vg68_rx/               SYNTHETIC, example only, not for design
    pair.txt                     key = value header
    solution.csv                 c_mass_fraction,T_C,P_MPa,viscosity_mm2s
```

`pair.txt` holds: `name`, `oil`, `refrigerant`, `source`, `note = example only, not for design`, `species` (solver gas model), `density_model` + `density_intercept_kg_m3`, `density_slope_kg_m3_K`, `oil_vapour_pressure_Pa`, `heat_capacity_J_kgK`, `thermal_conductivity_W_mK`, `valid_T_C = lo hi`, `valid_c = lo hi`.

`solution.csv` is the format the solver already reads (`solution_table_file`, FW68D layout). One row = one content curve point, so one table gives both views the dialog draws: viscosity vs T per content, and p_sat vs content at T. No solver change is needed. The dialog keeps the table as one grid (content × temperature), with saturation pressure and viscosity per cell, and pastes from Excel.

Range: the solver clamps with `solution_table_clamp = true` (the default in deck configs). The GUI checks the case (T_amb, content, and thermal T range after the run) against the table range before Run. Outside it, the Fluid node turns yellow and says so; the results card carries a yellow mark. There is no silent extrapolation.

Shipped examples (synthetic, labelled): `example_vg68_rx` and `examples/load_table_example.csv` (crank angle, magnitude, direction).

LGE data for tests only: `tests/data/lge/fw68d.csv` (copy of `data/R410A_FW68D/fw68d.csv`) and `tests/data/lge/os_load_130hz.csv`. Both stay out of `fluids/` and out of the release folder.

---

## 9. Code structure; files to delete or move

New files (each ≤ ~800 lines):

| File | Role |
|---|---|
| `case_model.hpp/.cpp` | CaseModel, `to_config`, `from_config`, hidden-key report |
| `config_text_patch.hpp/.cpp` | Text-preserving save |
| `fluid_library.hpp/.cpp` | Reads/writes `fluids/`, range checks |
| `estimates.hpp/.cpp` | S, p_sat, margin, onset verdict, run-time estimate |
| `onset_chart.hpp/.cpp` + `onset_curve.csv` | Φ(S) data and plot |
| `angle_reference.hpp` | Angle-zero transforms (fixed, h_max, h_min, feed, opposite feed) |
| `workflow_tree.cpp` | Tree, status lights, reasons |
| `panel_geometry.cpp`, `panel_fluid.cpp`, `panel_operating.cpp`, `panel_boundaries.cpp`, `panel_mesh.cpp`, `panel_solution.cpp` | One per node |
| `view_geometry.cpp`, `view_mesh.cpp` | Live views |
| `dialog_fluid_library.cpp` | Fluid dialog |
| `results_card.cpp`, `results_maps.cpp`, `results_polar.cpp`, `results_transient.cpp`, `results_compare.cpp` | Results |
| `run_monitor.cpp` | Monitors, stage, command line |
| `batch_runner.hpp/.cpp`, `panel_batch.cpp` | Sweep queue, scripts, live table |
| `case_report.cpp` | One-page report (self-contained HTML with the figures; prints to PDF) |

Delete: `panels_setup.cpp` (989 lines), the generic form in `schema.cpp` (`double_specs` and friends; the unit families move to `units.hpp/.cpp`), `apply_dissolved_gas_species_preset`, the presets combo, the support-stiffness mirroring (`normalize_motion_config`, `mirror_isotropic_support`, `warn_unrepresentable_motion_config`; replaced by the hidden-key report).

Split: `panels_results.cpp` (952 lines) into the results files above. Keep `solver_runner`, `results` (VTK reader), `exports`, `settings`, `theme`, `main_gui`.

Tests: new `tests/test_gui_case_model.cpp` (round trip, mappings, angle references, estimates), registered in CTest. The GUI build stays Windows-only; the case model and estimates are plain C++ and also build on Linux for the tests.

---

## 10. Acceptance tests: sources and recorded values

| Test | Source | Recorded value |
|---|---|---|
| 1 Round trip, O-run | `D:\Seongho\Sim\pancake\lge-thermal-v2-20261029\oil\O4\O4\0003_equilibrium_thermal3_806c7579\config.txt` | every solver key equal |
| 1 X-run | `lge-thermal-v3-c274-20261114\mix\gas_on\X4\runs\X4\0003_equilibrium_thermal3_474e01a5\config.txt` | |
| 1 Optical rig | `silicone-coupled-20261114-b271002-local\runs\A-G_002_pinned_720x192\config.txt` (CRLF) | |
| 1 Compressor steady | `lge-compressor-vc-steady-20261123-00e5d1a\runs\R-mean\001_newton0\config.txt` | |
| 1 Compressor transient | `lge-compressor-transient-20261119-02a4b6c\b2\rev03\config.txt` | |
| 2 Compressor GC+VC, 3450 N | B0 mean (`...02a4b6c\b0\mean\009_newton6`) | cavitated area 57.09 %, h_min 1.6877 µm, ε 0.94375, φ 289.57°, torque 0.07734 N m, power 63.17 W, p_max 439.5 bar |
| 2 Tester 500 rpm, 2610 N | P5g (`lge-fbns-tester-20261117\runs\P5g\runs\0007_newton4_0755bbea`) | h_min 3.242 µm, ε 0.9237, φ 307.19°, torque 0.1196 N m, p_max 90.58 bar, cavitated area 11.83 % |
| 2 Compressor transient | rev03 | h_min 1.6406–1.7619 µm over the cycle |
| 5 Estimates | deck p45, p43 | S = 0.016 (computed 0.0155–0.0159), p_sat 9.65 bar at 27.4 wt % |

The deck's 1.62–1.89 µm range is the steady h_min from maximum to minimum load (3,604 to 3,024 N), not one run. Test 2 at 3,450 N therefore compares with 1.6877 µm.

The cavitated-area metric must match the deck. `part_b.py` has two: `cavitated_fraction` (film content < 1 − 1e-8) and `bubble_area_fraction` (void ≥ 0.5). The build checks which one gives 57.09 % and uses that definition, with the definition stated in the glossary.

Path note: the stored configs carry afml-r6 paths (`/media/user/새 볼륨/...`). Round trip keeps them verbatim. Reproduction runs use local copies, which have matching sha256 per the survey.

---

## 11. Draft glossary (becomes `docs/GUI_GLOSSARY.md`)

| UI term | Definition | Symbol, unit | Deck term |
|---|---|---|---|
| Bearing diameter | Nominal bore (shaft) diameter | D = 2R, mm | 내경 |
| Bearing length | Axial length of the film | L, mm | 길이 |
| Radial clearance | Bore radius minus shaft radius | c, µm | 반경 틈새 |
| Oil feed (hole / groove) | Where oil enters the film | – | 급유 구멍 / 급유 그루브 |
| Feed angle | Feed centre, measured in the direction of rotation | θ_in = 0 at the feed, deg | 급유부 기준 각도 θ_in |
| Angle from maximum film | Angle from h_max in the direction of rotation | θ, deg | 최대 유막두께 기준 θ |
| Attitude angle | Angle from the load direction to the line of centres (h_min) | φ, deg | attitude angle φ |
| Tilt angle, tilt direction | Shaft axis relative to the bearing axis | γ (mrad), θ_γ (deg) | 기울어짐 γ |
| Speed | Shaft rotational speed | N, rpm | 회전 속도 |
| Load | Force on the bearing, magnitude and direction | W, N | 하중 |
| Specific load | Load per projected area | P = W/(LD), bar | 면압 |
| Ambient pressure | Pressure around the bearing (chamber, shell) | p_amb, bar | 주변 압력 (챔버 압력) |
| Ambient temperature | Temperature of the surroundings and the supplied oil | T_amb, °C | 주변 온도 |
| Oil-refrigerant pair | Oil and dissolved refrigerant with their property tables | – | 오일-냉매 |
| Refrigerant content | Mass fraction of refrigerant dissolved in the oil | w, wt % | 냉매 함량 |
| Saturated at ambient pressure | Content at which the oil is in equilibrium with the refrigerant at p_amb, T_amb | – | 포화 오일 |
| Saturation pressure | Pressure below which dissolved refrigerant comes out of solution at the given content and temperature | p_sat, bar | 포화압력 |
| Pressure margin | Ambient minus saturation pressure | p_amb − p_sat, bar | 압력 여유 |
| Pressure drop ratio | How far the film pressure falls below ambient, per specific load | Φ = (p_amb − p_min)/P | 압력 강하 Φ |
| Sommerfeld number | Bearing duty parameter | S = μ N L D (R/c)² / W | Sommerfeld 수 S |
| Vapour cavitation (VC) | Film ruptures when the oil itself boils at its vapour pressure | VC | VC (증기 캐비테이션) |
| Gaseous cavitation (GC) | Dissolved refrigerant comes out of solution when the film pressure falls below p_sat | GC | GC (기체 캐비테이션) |
| Dissolved-gas release | Model setting: GC on (GC+VC) or off (VC only) | – | GC+VC 모델 / VC 모델 |
| Cavitated region | Part of the film filled partly by gas or vapour (definition fixed in the build, section 10) | area % | 버블 영역 |
| Gas fraction | Volume fraction of gas or vapour in the gap | α | 기체 부피분율 α |
| Minimum film thickness | Smallest gap | h_min, µm | 최소유막두께 |
| Eccentricity ratio | Shaft centre offset per clearance | ε = e/c | 편심률 |
| Friction torque | Viscous torque on the shaft | M, N m | 유막 토크 / 마찰 토크 |
| Friction power loss | Torque times angular speed | P_f = M ω, W | 마찰 손실 |
| Peak pressure | Highest film pressure | p_max, bar | 최고 압력 |
| Film thickness ratio | Minimum film per combined roughness | Λ = h_min/σ, σ = (Rq₁² + Rq₂²)^½ | 유막두께비 Λ |
| Shear share of film force | Viscous part of the total film force | F_v/F, % | – |
| Load balance | Iteration that moves the shaft until the film force carries the load | – | 하중 평형 |
| Iterative solution | Repeated solving until the residuals fall below the tolerances (Expanded: FBNS, Newton, Picard, Krylov) | – | 반복 해석 |
| Mesh level | Number of cells around and along the film | n_θ × n_z | 격자 |
| Thermal analysis | Film, housing and shaft temperatures solved | – | 열 해석 |
| Steady / transient | Fixed operating point / load varying in time | – | 정상 / 과도 해석 |

---

## 12. Open items and risks

1. D2 determines whether test 2 runs from `pancake.exe -c` alone. Until it is settled, the load-balanced GC+VC steady case is the main risk.
2. The tester case in the plateau form (D5) may differ from P5g (TWO_PHASE). The test reports the difference; it does not hide it.
3. Viscosity rule for isothermal FBNS: the deck configs use `liquid_viscosity_model = CONSTANT` with μ set by the driver. The GUI must use the same rule (table at T_amb and content); it is checked against `part_b.py` and `run_fbns.py`.
4. Onset curve: the Φ values are computed from the G1–G8 run outputs (p_min). If those outputs are incomplete, I stop on this item and say so.
5. Field refresh in steady runs costs I/O (one VTK set per 10 outer iterations; 720 × 192 ASCII). If it slows runs noticeably, the interval becomes adaptive.
6. Screenshots (test 9): the agent session cannot capture the GUI window (noted in memory). The plan is an in-app "Save screenshot of every screen" command that uses the existing PNG capture path, so the screenshots are produced by the app itself.
7. Release folder: `pancake.exe`, `pancake_gui.exe`, the MS-MPI runtime (installer, `msmpi.dll`), the MinGW runtime DLLs if not static, PETSc DLLs if dynamic, `fluids/`, `examples/`. The exact list comes from a clean-machine check (test 8).
8. Build effort is large (about 25 new files). Proposed order: D1 merge → D3 writer fix → case model + tests → fluid directory → panels and tree → live views and estimates → results → batch and comparison → report → screenshots and docs. A commit after each step.
