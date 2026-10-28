# LGE thermal sensitivity campaign

Prepared from the confirmed 40 mm bearing and the corrected experimental torque comparison. **Production is not launched.** The first group contains eight additional load-balanced equilibria. The six September 28 results supply the baseline.

The implementation is [run_steady.py](run_steady.py). [prepare_sensitivity.py](prepare_sensitivity.py) generates complete setup files and verifies their actual native configuration differences. [summarize_sensitivity.py](summarize_sensitivity.py) rechecks accepted results and creates the comparison figures. Preparation evidence and generated inputs are in [the campaign directory](../../results/lge_sensitivity_20261011).

## Baseline and source match

The reference is [the October 8 discrepancy review](<F:/Seongho/Research/chats/2026-10-08 LGE experiment and simulation error review.md>) and [the corrected experimental table](<F:/Seongho/Research/attachments/LGE-2026-10-02-stribeck-rpm-drag-cases.csv>). Use `new_exp_torque_nm`, not `exp_torque_nm` or the older report error columns. Error is `100*(CFD/new_exp_torque_nm - 1)`; temperature errors are CFD minus experiment in kelvin.

| Case | Corrected experimental torque (N m) | Baseline CFD (N m) | Torque error (%) | TC1 / TC2 / TC3 error (K) |
|---|---:|---:|---:|---|
| A | 0.373168 | 0.268886 | -27.945 | +0.207 / +0.164 / +0.357 |
| B | 0.461914 | 0.308961 | -33.113 | +0.237 / +0.301 / +1.913 |
| C | 0.358950 | 0.197388 | -45.010 | +0.046 / +0.016 / +0.157 |
| D | 0.161938 | 0.148962 | -8.013 | +0.481 / +0.483 / +0.550 |
| E | 0.209897 | 0.289616 | +37.980 | +4.119 / +4.092 / +4.821 |
| F | 0.109124 | 0.143498 | +31.500 | +0.674 / +0.678 / +0.719 |

The initial archive comparison covered 84 source/case/build files. The only differences from `resolved_snapshot/source.tar.gz` were the case README and a comment in `src/resolved_cross_film.hpp`. Native equations, driver, setup, experimental inputs and viscosity table matched. The archived aggregate source hash differs because the header comment participates in that hash. This is not a changed physical model. The revised driver's no-scenario rendering is byte-identical to the archived driver for all six cases, including reference temperatures, supplied h coefficients and initial positions. The baseline accepted-state metadata and input hashes were checked. No matched-baseline reruns are required for this preparation. Any later prediction-changing source/configuration edit requires a fresh match before attributing a difference to a physical input.

Geometry stays R=21.5 mm, L=40 mm, c=42.5 micrometres. Load is upward, groove centre 90 degrees, groove width 5 mm. Both axial pressures and groove pressure stay 101325 Pa absolute; p_cav stays 0 Pa absolute. Every sensitivity uses native force equilibrium. Eccentricity and attitude remain solved outputs. The 720 by 192 pressure plane, 16 film layers and 16 by 16 housing radial/axial resolution stay fixed.

## Implemented first group

| Setup | Cases | Axial incoming oil and groove supply | Housing and shaft baths | New equilibria |
|---|---|---|---|---:|
| `supply_plus_1K.json` | E | TC0 + 1 K | original TC0 | 1 |
| `cooling_plus_1K.json` | E | original TC0 | TC0 + 1 K | 1 |
| `supply_TC2.json` | A-F | each measured TC2 | original TC0 | 6 |

The TC2 inputs below retain the experimental table's precision. Kelvin values in generated JSON are calculated as Celsius + 273.15.

| Case | TC2 (degrees C) | TC2 - TC0 (K) |
|---|---:|---:|
| A | 31.48782181362495 | 1.814274693107193 |
| B | 65.8447735548189 | 4.35786116068954 |
| C | 30.24139772640574 | 1.601269041987745 |
| D | 29.39117452990355 | 1.086833693448031 |
| E | 54.967458032366274 | 3.402198369050514 |
| F | 29.31858253758852 | 0.946947600226996 |

For E, original TC0 is 324.71525966331575 K. Supply +1 K changes only `bc_z_south_temperature`, `bc_z_north_temperature`, and the fourth `inlet_groove` token to 325.71525966331575 K. Cooling +1 K changes only `bush_outer_temperature` and `shaft_ambient_temperature` to that value. Full-precision machine-readable changes are in `prepared_verified/*_config_diff.json`. The additional `output_fields` entries are diagnostic output only: native density and F0/F1/F2/Couette-thickness fields permit exact torque and flux decomposition without changing native equations.

The new top-level setup object `thermal_boundaries_by_case` maps case IDs to any of these five kelvin temperatures. `groove_supply_temperature` maps explicitly to the fourth native groove token. Missing values retain that case's TC0. Unknown cases/keys, nonnumeric/nonfinite/nonpositive values and temperatures outside the viscosity table are rejected. `common_solver_settings` and `per_case_mapping` remain descriptive; editing those sections alone is not a scenario implementation.

Incoming axial thermal data remain inflow-only in the unchanged native solver. Groove mixing remains `FLUX_BALANCED`; the groove film temperature and bearing wall temperature remain solved. `temperature_initial`, `temperature_reference`, `property_reference_temperature`, `thermal_ambient_temperature`, reference rho/cp/mu and the viscosity law retain their baseline definitions. Local viscosity still responds to solved temperature. Housing and shaft correlations use their own actual bath temperatures for both driving temperature and oil properties, including when rechecking cached coefficient closure. The housing bath also controls the native annular-end calculation.

The physical temperature floor is the coldest actual supply/cooling reservoir minus 0.001 K. It replaces the common-TC0 assumption without weakening the original case. The full setup hash and rendered configuration isolate scenario caches. Changed temperatures remain incompatible with prior-physics checkpoint matching.

TC2 is an embedded housing measurement, not an inlet measurement. TC2-supply cases are hypotheses, and TC2 residuals are marked dependent diagnostics rather than independent validation. E_low/E_high changed supply, cooling and reference properties together. Their derivative is not used as a supply-only sensitivity.

## Run commands and cost

Local preparation, already performed from `C:/Dev/pancake`:

```powershell
python .\cases\02d_fw68d_steady_six\run_steady.py --help
python .\cases\02d_fw68d_steady_six\run_steady.py --self-check
python .\tests\test_fw68d_force_balance.py
cmake --fresh --preset windows-native-mingw
cmake --build --preset windows-native-mingw-release
ctest --test-dir build-windows-mingw --output-on-failure
```

The preparation generator requires a new output directory. This preserves earlier input packages:

```powershell
python .\cases\02d_fw68d_steady_six\prepare_sensitivity.py `
  --experiment-csv 'F:\Seongho\Research\attachments\LGE-2026-10-02-stribeck-rpm-drag-cases.csv' `
  --baseline .\results\fw68d_physical_20260928\resolved_results `
  --archive-source .\results\fw68d_physical_20260928\resolved_snapshot\source.tar.gz `
  --out .\results\lge_sensitivity_20261011\prepared_new
```

The selected snapshot is on `afml-r6`, reached through `user@143.248.99.6`, at `/home/user/pancake-runs/lge-sensitivity-20261011-thermal`. It does not modify an existing remote checkout. Four MPI ranks use the previously verified physical cores 0-3. There are 28 physical cores, so this conservative sequential schedule leaves room for other work. Recheck occupancy before launch. Build/test commands in that directory are:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cmake -S . -B build-linux -DCMAKE_BUILD_TYPE=Release
cmake --build build-linux -j8
ctest --test-dir build-linux --output-on-failure
python3 cases/02d_fw68d_steady_six/run_steady.py --self-check
python3 tests/test_fw68d_force_balance.py
```

After approval, run these three commands sequentially from that remote root. No unresolved placeholders remain:

```bash
python3 cases/02d_fw68d_steady_six/run_steady.py --solver "$PWD/build-linux/pancake" --ranks 4 --rankfile "$PWD/inputs/four_cores.rankfile" --cases E --materials "$PWD/inputs/supply_plus_1K.json" --work "$PWD/results/supply_plus_1K" --timeout 7200
python3 cases/02d_fw68d_steady_six/run_steady.py --solver "$PWD/build-linux/pancake" --ranks 4 --rankfile "$PWD/inputs/four_cores.rankfile" --cases E --materials "$PWD/inputs/cooling_plus_1K.json" --work "$PWD/results/cooling_plus_1K" --timeout 7200
python3 cases/02d_fw68d_steady_six/run_steady.py --solver "$PWD/build-linux/pancake" --ranks 4 --rankfile "$PWD/inputs/four_cores.rankfile" --cases A,B,C,D,E,F --materials "$PWD/inputs/supply_TC2.json" --work "$PWD/results/supply_TC2" --timeout 7200
```

The prepared `run_first_group.sh` retains these commands and per-group driver logs. It checks the selected snapshot and input hashes, aborts at the first failed command, and refuses existing production result directories.

| Quantity | Estimate from full archived solve histories |
|---|---:|
| Additional equilibria | 8 |
| Matched baseline repeats | 0 |
| Expected native solves across h-closure passes | about 27 |
| Sequential four-rank wall time | 2.384 h |
| CPU allocation | 9.534 core-hours |
| Planning allowance, 1-2 times historical duration | 2.4-4.8 h, 9.5-19.1 core-hours |

This estimate sums every A-F archived native call, then adds two more E-equivalent costs. It is not based only on each final call. Changed thermal conditions may alter iteration counts; this range is a planning allowance, not a measured confidence interval. The default allows up to 30 cooling passes per equilibrium, hence 240 native calls across eight successful cases at the cap. `--timeout 7200` applies to each native solve, with 120 seconds cooperative-stop grace. It is not a two-hour campaign deadline. Stop and diagnose a failed/capped solve; do not loosen acceptance gates. No billing rate was supplied, so monetary cost is unknown.

Every changed-physics case starts fresh. Within one unchanged scenario the driver may use its own verified previous h-closure checkpoint. Do not pass `--initial-state-root` for this group. `resolved_results` is not a compatible checkpoint work root, and an accepted baseline is never reused as a changed-physics prediction. Do not treat a `--trial` or deliberately capped preflight as an equilibrium.

## Acceptance and outputs

The driver retains configuration, source, executable, driver, setup, experimental-table, property-table, rankfile and checkpoint hashes, exact native commands, elapsed solve time, MPI/Python versions, logs and return codes for each solve. The transfer manifest and workstation build/environment logs identify the complete selected snapshot.

Acceptance is unchanged except for the correctly generalized temperature floor. Require native final CONVERGED and independent force residual <0.002; inner FBNS pressure/complementarity <1e-6; normalized temperature and profile residual <1e-6; clean coupling convergence without restored best or cap; relative mass residual <=1e-8; global energy residual <=1e-6; valid viscosity-table domain, finite positive gap/viscosity/density and admissible filling. External h mismatch must be <1e-3 and mean surface change <0.01 K. Pressure-field change is reported, not substituted for these gates.

Postprocessing, after all eight equilibria pass:

```bash
python3 cases/02d_fw68d_steady_six/summarize_sensitivity.py --prepared "$PWD/inputs" --work-root "$PWD/results" --solver "$PWD/build-linux/pancake" --out "$PWD/analysis"
```

The reporter independently rechecks provenance, the exact scenario configuration, physical acceptance and coefficient closure before plotting. It retains total torque, Couette/pressure-shear split and corrected torque errors; TC1-3 errors; film temperature and viscosity; minimum gap, eccentricity, attitude, peak pressure and cavity area/arc; incoming/outgoing mass at each axial end and groove; viscous heat, wall/advective energy, separate cylinder/end/disk cooling and closure diagnostics. Positive split mass rates denote the named incoming or outgoing direction. Heat removal is positive away from the film. Native signed energy fields are retained unchanged.

The figure set compares baseline/perturbed torque and temperature residuals, per-kelvin torque/temperature responses and changes in absolute errors. The TC2-supply slope is a secant over TC2-TC0, not a small-step derivative. Improvement means smaller absolute error for a particular output; reductions in torque error do not automatically validate the temperature or cavity prediction. TC2 dependence must appear in the figure caption and report.

Retrieve the raw production roots together with analysis, inputs and verification, preserving failed runs and source metadata:

```powershell
scp -r user@143.248.99.6:/home/user/pancake-runs/lge-sensitivity-20261011-thermal/results .\results\lge_sensitivity_20261011\remote_results
scp -r user@143.248.99.6:/home/user/pancake-runs/lge-sensitivity-20261011-thermal/analysis .\results\lge_sensitivity_20261011\remote_analysis
```

## Separate optional groups

These are selectable follow-up groups, not included in the eight-run script. Start with Case E, then a Case C control where noted; extend to A-F only after seeing direction and magnitude. Each row changes one physical input at a time and recloses load and thermal balance. The thermal-boundary interface does not silently implement these other options. Narrow driver/native wiring and a fresh configuration-diff check are required wherever stated.

| Group | Proposed first comparison / extra equilibria | Implementation and evidence needed |
|---|---|---|
| Housing outer-cylinder cooling | E at 0.5x and 2x h correlation, 2 | Multiply only `housing_coefficient` results in both iteration and cache recheck. Does not multiply annular ends. Illustrative amplitudes, not measured uncertainty. |
| Annular-end cooling | E end cooling disabled, 1; later separate 0.5x/2x, 2 | Native `bush_end_natural_convection=false`, `bush_end_heat_transfer=0` is an endpoint hypothesis. A multiplier on the local natural-convection correlation needs a small solver extension and tests. A constant coefficient is a different closure. |
| Shaft cooling | E at 0.5x/2x, 2 | Existing `material_settings.shaft_end_nusselt_coefficient` 0.31/1.24 changes the actual shaft correlation, holding baths and housing fixed. Not measured bounds. |
| Clearance | E and C at 41.5/43.5 micrometres, 4 | Explicit geometry override needed. Update `c`, seed `e=eps*c`, physical bore R+c, wall thickness and derived metadata; retain fixed shaft R, L and OD. Native CHT inner radius remains R under the existing thin-film approximation. Use measured clearance bounds when available. |
| Pure-oil viscosity law | E and C for two independently supported tables, 4 when available | Copy the actual nu(T) table and wire its path and points into native and Python property/cooling calculations, hashes and domain checks. A reference mu change alone is ineffective as a physical-law sensitivity. No arbitrary numeric law bounds proposed. |
| Thermal properties | E, one at a time: cp and oil k at +/-20%, bush k at +/-20%, 6 | `material_settings.oil_cp_J_kg_K`, `oil_conductivity_W_m_K`, `bush_conductivity_W_m_K` already enter native configs. The first two also feed Python correlations and native annular cooling. The amplitudes are illustrative. Driver self-check assumptions must match any selected physical scenario. |
| Load and pressure | One low/high pair per independently supported scalar input, 2 each | Wire applied force, each axial pressure, and groove pressure separately. The existing material reservoir-pressure setting ties all three reservoirs, so it cannot isolate them. Maintain absolute datum and pressure >= p_cav. Logged chamber pressure is not automatically groove pressure. |
| Cavitation and cavity shear | One supported alternative per selected case | `material_settings.cavitation_pressure_Pa` reaches native p_cav, but common self-check defaults must be scenario-aware. Threshold/shear alternatives need explicit wiring and their own validity checks. Gas release/starvation may require model and measurement work. |
| Numerical error | B: theta-only 720->1440 and z-only 192->384, 2; then coupled solid radial-only and axial-only 16->32, 2 | Explicit mesh/output wiring and fresh starts; never reuse incompatible layer/mesh checkpoints. Tighten force/thermal criteria separately afterward. Existing through-film and fixed-flux solid checks do not certify these directions. |
| Omitted energy | No repeat of existing B pressure-work pair | Reuse [the saved pressure-work comparison](../../results/fw68d_pressure_work_20260929/analysis/comparison_report.md). Any additional term needs an independent physical justification and one-term paired comparison. |

The optional counts are menu entries, not an approved or recommended total. Load, pressure, oil law and cavity-model counts cannot be finalized without supported alternative inputs. Misalignment, shaped clearance, thermoelastic deformation, mixed lubrication and resolved groove recirculation are model extensions. Sensor-position uncertainty should first be explored by resampling existing solid fields. Motor calibration, parasitic subtraction and unknown historical inputs require experimental evidence; bearing CFD alone cannot repair them.

Reuse [B film-layer evidence](../../results/fw68d_physical_20260928/resolved_results/B_film32_verification.json): 16->32 layers changes torque +0.0047626%, peak pressure +0.0225613%, and TC1-3 by at most 0.0020015 K. Existing A/B fixed-interface-flux housing checks change sensor temperatures by no more than about 0.0024 K, but do not reclose the coupled fluid response. B pressure work changes torque by about -0.023%, and TC1/TC2/TC3 by +0.051/+0.041/-0.013 K. It cannot explain B's roughly -33% corrected torque discrepancy. The pressure-work report's old stationary-journal prose is superseded by the archived source convention: rotating journal, translated bearing equilibrium.

## Interpretation before new predictions

Existing fields and the corrected torque comparison require torque increases in A-D and decreases in E/F. Warmer incoming oil has the expected viscosity-driven sign for E/F and the opposite screening sign for A-D. The coupled runs must determine magnitude and any compensating load/cavity/heat-flow response. No new thermal sensitivity direction or magnitude is claimed yet. Common-boundary E_low/high results are not substituted for the separate supply/cooling tests.

After the eight runs, compare every scenario against the same-case baseline, report absolute and per-kelvin changes, and flag worsened independent sensor errors. Classify mechanisms as sufficient, directionally useful but too small, adverse to other outputs, or unidentifiable. Do not infer unique inlet temperature, clearance, viscosity or parasitic torque by fitting away the torque residual. The [validation ledger](../VALIDATION_STATUS.md) continues to limit physical-validation claims.

## Verified preparation outcome, October 12

- Native Windows configure/build and all 59 CTest tests passed. The extended Python driver/reporting regression and all three scenario self-checks passed.
- The isolated workstation build passed 58 tests initially; the composition-viscosity test lacked its transferred CSV fixture. After transferring the original fixture, that test passed. The complete set of 59 tests is covered. The original failed log is retained.
- The first two preflight launch attempts encountered CRLF in the Windows-written OpenMPI rankfile. The generator now writes LF explicitly. Failed launches, the old rankfile and repaired input hashes are retained.
- The two repaired four-rank checks each performed one full-mesh native coupled iteration: supply +1 K, 41.63 s; cooling baths +1 K, 41.53 s. Both returned native code 4 at the intended cap and were rejected by the driver. They are diagnostic checks, not physical predictions. Output fields were present and the torque decomposition matched the native total within 2.1e-15 N m.
- Remote executable SHA-256 `62ddd6a0906943e2a7942330742757042626036602e8ab1b2073e0e563969551` is identical to the archived resolved-campaign executable. This provides stronger baseline matching than the source-comment audit alone.
- Historical production storage is 3.7 GiB for the original batch. Reserve approximately 5-10 GiB for the eight-run group and retained intermediate outputs; this is a planning allowance. The workstation had about 68 GiB free before preparation.
- Final inputs are in [prepared_verified](../../results/lge_sensitivity_20261011/prepared_verified). Remote test logs, manifest, environment and capped outputs are retrieved under the same local campaign directory. Synthetic plot-layout checks were inspected separately and are explicitly labelled test-only. Actual comparison plots and physical conclusions await the eight approved equilibria.

At the preparation checkpoint, no production scenario had been launched. The eight runs were subsequently authorized and started as recorded below.


## Parallel production launch, October 12

The user authorized all eight prepared equilibria and concurrent four-rank runs.
Production started at **09:57:50 KST (00:57:50 UTC)** on afml-r6. Seven processes
use disjoint four-core allocations over all 28 physical cores; MPI affinity was
independently probed and checked before launching any solver. Logical sibling
CPUs share a physical-core allocation and are not counted as extra cores.

Initial cases: E supply +1 K, E baths +1 K, and TC2-supply E/F/B/D/C.
TC2-supply A takes the first free allocation. Every case has a separate result
subdirectory and driver log. The validated solver, physical inputs and driver
are unchanged. Failed runs remain preserved and do not block unrelated jobs.

The exact launcher, commands, rankfile hashes, core allocations and PIDs are in
[the launch record](../../results/lge_sensitivity_20261011/PARALLEL_LAUNCH.json)
and [parallel launcher](../../results/lge_sensitivity_20261011/run_parallel.py).
Remote `parallel/status.json` updates each case's state, exit code and timing.
The historical serial estimate of 2.4-4.8 hours is not the elapsed estimate for
this parallel launch; shared memory and I/O mean the speedup must be measured.

A detached postprocessor waits for all eight successful driver exits, rechecks
every equilibrium using its exact per-case rankfile, generates the comparisons,
and archives all raw successful and failed attempts for retrieval. Its state is
`postprocess_status.json`. The added reporter option is `--execution-manifest`;
missing, failed, duplicate or hash-mismatched execution records are rejected by
regression checks. The original snapshot reporter remains preserved remotely;
the parallel-aware copy is `summarize_sensitivity_parallel.py`.

No new physical prediction is accepted merely because a job is running or its
inner pressure residual has converged. Full force/thermal/conservation checks
remain required before interpretation.


## Completed thermal sensitivity results October 12

All eight equilibria passed the unchanged native and independent acceptance
checks. The seven-slot parallel campaign took **25 min 56.6 s**, with 26 native
solves and 9.877 allocated core-hours. The complete 864,594,919-byte archive was
retrieved, SHA-256 verified and extracted locally; all eight final configuration
and checkpoint hashes were independently checked after extraction.

**Kinematic scope:** these results use a shaft with a fixed centre and rotating
surface, and a nonrotating bearing whose centre is adjusted for steady load
balance. Position adjustment introduces no physical translational velocity.
On October 12, the user confirmed that the shaft centre is fixed while the
shaft rotates and the bearing centre moves. This matches the completed steady
equilibria and resolves the motion clarification. No kinematic change or rerun
is required.

### Does changing temperature boundary conditions explain torque

Within the implemented model, **no tested temperature boundary condition fixes
the torque discrepancy**. Every tested warming lowers torque. A-D already
underpredict torque, so their discrepancy worsens. E and F have the helpful
sign, but the TC2-supply hypothesis removes only 7.90% and 1.44% of their original
absolute torque gaps. Their errors remain +34.98% and +31.05%, respectively.
Every scenario also raises TC1-TC3, worsening the positive temperature errors,
including independent TC1 and TC3 in the TC2-input cases.

For E, supply +1 K changes torque by -1.845 mN m and TC1-TC3 by +0.130 to
+0.134 K. Baths +1 K changes torque by -6.192 mN m and sensors by +0.593 to
+0.597 K. The bath-temperature torque response is 3.36 times the supply response
for these one-kelvin tests. Only 2.31% and 7.77% of E's original torque gap are
removed, respectively, while the already excessive predicted temperatures rise.

TC2 is a housing measurement used here as a hypothetical inlet value. Its
residual is descriptive, not independent validation, for A-F TC2-supply cases.
The results do not identify an actual inlet temperature, clearance, viscosity
law or parasitic torque. They do not establish the effects of cooling strength,
which was not varied, or justify extrapolating these slopes to fit the residual.

### Torque comparison

Torque error is 100(M_CFD/M_exp - 1), using only corrected `new_exp_torque_nm`.
Positive error means excess CFD torque. Contributions sum to total torque;
all torque values in this table are N m.

| Scenario | Baseline torque | Perturbed torque | Couette | Pressure shear | Experiment | Error baseline (%) | Error perturbed (%) |
|---|---|---|---|---|---|---|---|
| E supply +1 K | 0.289616 | 0.287770 | 0.286220 | 0.001550 | 0.209897 | +37.98 | +37.10 |
| E baths +1 K | 0.289616 | 0.283424 | 0.281849 | 0.001575 | 0.209897 | +37.98 | +35.03 |
| A TC2 supply | 0.268886 | 0.263833 | 0.212743 | 0.051090 | 0.373168 | -27.95 | -29.30 |
| B TC2 supply | 0.308961 | 0.291478 | 0.256448 | 0.035030 | 0.461914 | -33.11 | -36.90 |
| C TC2 supply | 0.197388 | 0.194006 | 0.172144 | 0.021862 | 0.358950 | -45.01 | -45.95 |
| D TC2 supply | 0.148962 | 0.147638 | 0.145608 | 0.002029 | 0.161938 | -8.01 | -8.83 |
| E TC2 supply | 0.289616 | 0.283318 | 0.281741 | 0.001577 | 0.209897 | +37.98 | +34.98 |
| F TC2 supply | 0.143498 | 0.143003 | 0.142481 | 0.000522 | 0.109124 | +31.50 | +31.05 |

### Thermocouple errors and temperature response

Errors are CFD minus experiment in K. Baseline and perturbed values appear as `baseline -> perturbed`. The last column is the TC1/TC2/TC3 temperature change per kelvin of imposed boundary change; TC2-supply values are finite-interval secants.

| Scenario | TC1 error (K) | TC2 error (K) | TC3 error (K) | TC response (K/K) |
|---|---|---|---|---|
| E supply +1 K | +4.119 -> +4.253 | +4.092 -> +4.223 | +4.821 -> +4.951 | 0.134/0.131/0.130 |
| E baths +1 K | +4.119 -> +4.713 | +4.092 -> +4.688 | +4.821 -> +5.418 | 0.593/0.596/0.597 |
| A TC2 supply | +0.207 -> +0.806 | +0.164 -> +0.769 | +0.357 -> +0.965 | 0.330/0.333/0.335 |
| B TC2 supply | +0.237 -> +2.077 | +0.301 -> +2.111 | +1.913 -> +3.719 | 0.422/0.415/0.415 |
| C TC2 supply | +0.046 -> +0.501 | +0.016 -> +0.474 | +0.157 -> +0.619 | 0.284/0.286/0.288 |
| D TC2 supply | +0.481 -> +0.610 | +0.483 -> +0.612 | +0.550 -> +0.679 | 0.119/0.119/0.119 |
| E TC2 supply | +4.119 -> +4.581 | +4.092 -> +4.542 | +4.821 -> +5.269 | 0.136/0.133/0.132 |
| F TC2 supply | +0.674 -> +0.736 | +0.678 -> +0.740 | +0.719 -> +0.781 | 0.065/0.065/0.065 |

### Film temperatures and viscosity

Film values are cell-centre statistics; the resolved cross-film maximum additionally includes film layers. Temperature is in degrees C and dynamic viscosity in Pa s.

| Scenario | Film min/mean/max (C) | Layer max (C) | Viscosity min/mean/max (Pa s) |
|---|---|---|---|
| E supply +1 K | 55.865/59.214/59.482 | 59.790 | 0.030393/0.030724/0.035137 |
| E baths +1 K | 55.537/59.604/59.909 | 60.220 | 0.029877/0.030248/0.035601 |
| A TC2 supply | 31.688/32.250/32.353 | 32.390 | 0.102169/0.102713/0.105729 |
| B TC2 supply | 66.194/67.775/68.440 | 68.843 | 0.021222/0.021800/0.023221 |
| C TC2 supply | 30.405/30.735/30.792 | 30.824 | 0.110725/0.111046/0.112955 |
| D TC2 supply | 29.896/30.036/30.053 | 30.071 | 0.115019/0.115122/0.116215 |
| E TC2 supply | 57.246/59.618/59.838 | 60.134 | 0.029962/0.030229/0.033244 |
| F TC2 supply | 29.990/30.092/30.102 | 30.117 | 0.114727/0.114787/0.115419 |

### Solved geometry pressure and cavity extent

Eccentricity and attitude are solved outputs. Cavity area and midplane extent use film content below 0.999; the extent is the total angular measure, not necessarily one contiguous arc.

| Scenario | Minimum gap (um) | Eccentricity ratio | Attitude (deg) | Peak pressure (MPa abs) | Cavity area (%) | Midplane extent (deg) |
|---|---|---|---|---|---|---|
| E supply +1 K | 35.667 | 0.160771 | 354.155 | 0.507803 | 6.413 | 45.5 |
| E baths +1 K | 35.559 | 0.163307 | 353.970 | 0.507911 | 6.463 | 46.0 |
| A TC2 supply | 5.762 | 0.864413 | 301.823 | 11.324682 | 29.565 | 122.5 |
| B TC2 supply | 11.674 | 0.725318 | 314.917 | 5.351527 | 26.117 | 106.5 |
| C TC2 supply | 12.734 | 0.700381 | 318.187 | 3.177774 | 23.738 | 100.5 |
| D TC2 supply | 31.740 | 0.253182 | 352.990 | 0.452751 | 6.525 | 47.5 |
| E TC2 supply | 35.562 | 0.163238 | 354.233 | 0.509084 | 6.454 | 45.5 |
| F TC2 supply | 37.016 | 0.129032 | 358.794 | 0.260197 | 0.000 | 0.0 |

### Mass transport

Positive incoming/outgoing magnitudes in g/s. South and north figures exclude groove reservoir faces, which are accounted for in the groove columns. Split fluxes reproduce the native axial and groove balances.

| Scenario | South in | South out | North in | North out | Groove in | Groove out |
|---|---|---|---|---|---|---|
| E supply +1 K | 0.132472 | 0.346075 | 0.132472 | 0.346075 | 0.430017 | 0.002810 |
| E baths +1 K | 0.134617 | 0.350875 | 0.134617 | 0.350875 | 0.436650 | 0.004134 |
| A TC2 supply | 0.082202 | 0.272166 | 0.082202 | 0.272166 | 0.412575 | 0.032647 |
| B TC2 supply | 0.346728 | 1.657990 | 0.346728 | 1.657990 | 2.622524 | 0.000000 |
| C TC2 supply | 0.052578 | 0.220586 | 0.052578 | 0.220586 | 0.336017 | 0.000000 |
| D TC2 supply | 0.028288 | 0.075124 | 0.028288 | 0.075124 | 0.093671 | 0.000000 |
| E TC2 supply | 0.133902 | 0.352176 | 0.133902 | 0.352176 | 0.437564 | 0.001017 |
| F TC2 supply | 0.016275 | 0.037669 | 0.016275 | 0.037669 | 0.042789 | 0.000000 |

### Heat generation and removal

All values are W. Wall removal equals heat entering the housing plus heat entering the shaft. Net axial enthalpy out and groove enthalpy in are the native channels; their reference is the unchanged numerical temperature reference. The full energy ledger and negligible mass-imbalance correction are retained in the data.

| Scenario | Viscous generation | Into housing | Into shaft | Outer cylinder out | Annular ends out | Shaft disk out | Axial enthalpy out | Groove enthalpy in |
|---|---|---|---|---|---|---|---|---|
| E supply +1 K | 45.20402 | 8.14180 | 27.66089 | 4.32267 | 3.81910 | 27.66089 | 10.22559 | 0.82426 |
| E baths +1 K | 44.52119 | 7.39260 | 25.57484 | 3.92833 | 3.46425 | 25.57484 | 11.49255 | -0.06120 |
| A TC2 supply | 5.52408 | 1.52934 | 3.07908 | 0.80827 | 0.72106 | 3.07908 | 2.26091 | 1.34524 |
| B TC2 supply | 45.77691 | 6.75440 | 25.40703 | 3.57043 | 3.18398 | 25.40703 | 36.58695 | 22.97148 |
| C TC2 supply | 4.06175 | 1.17538 | 2.43312 | 0.62108 | 0.55430 | 2.43312 | 1.53474 | 1.08149 |
| D TC2 supply | 3.09222 | 0.92597 | 1.97087 | 0.48914 | 0.43683 | 1.97087 | 0.40000 | 0.20463 |
| E TC2 supply | 44.50429 | 8.69516 | 29.14182 | 4.61231 | 4.08283 | 29.14182 | 9.64393 | 2.97662 |
| F TC2 supply | 2.99504 | 0.92035 | 1.95753 | 0.48618 | 0.43417 | 1.95753 | 0.19860 | 0.08144 |

### Convergence and conservation

All values below are absolute relative residuals, except the last column in K. Gates are force <0.002, thermal and layer residuals <1e-6, mass <1e-8, energy <1e-6, cooling coefficient mismatch <0.001 and cooling temperature change <0.01 K. Native inner FBNS residuals also pass 1e-6. No capped or restored-best iterate is accepted.

| Scenario | Force | Thermal | Film layers | Mass | Energy | Cooling h | Cooling dT (K) |
|---|---|---|---|---|---|---|---|
| E supply +1 K | 1.298e-03 | 5.077e-07 | 5.289e-07 | 1.026e-10 | 3.404e-12 | 3.741e-04 | 0.007221 |
| E baths +1 K | 1.252e-03 | 4.405e-07 | 4.599e-07 | 7.001e-11 | 1.112e-11 | 2.542e-04 | 0.004572 |
| A TC2 supply | 8.030e-04 | 4.009e-07 | 5.502e-07 | 1.592e-10 | 5.249e-09 | 7.251e-04 | 0.005438 |
| B TC2 supply | 1.340e-03 | 3.469e-07 | 5.621e-07 | 9.682e-11 | 2.227e-09 | 5.842e-05 | 0.000917 |
| C TC2 supply | 8.608e-04 | 1.219e-07 | 1.438e-07 | 4.544e-10 | 1.695e-09 | 9.599e-04 | 0.005938 |
| D TC2 supply | 1.430e-03 | 1.786e-07 | 1.839e-07 | 1.235e-10 | 4.128e-11 | 7.887e-05 | 0.000414 |
| E TC2 supply | 1.114e-03 | 3.954e-07 | 4.117e-07 | 4.128e-11 | 1.435e-11 | 4.604e-04 | 0.009267 |
| F TC2 supply | 1.686e-03 | 2.460e-07 | 2.518e-07 | 1.007e-14 | 2.310e-11 | 9.255e-05 | 0.000484 |


The maximum torque reconstruction error is 8.9e-16 N m. Maximum separately
reconstructed housing and shaft heat imbalances are 2.76e-5 W and 7.90e-8 W.
All changed scenarios started fresh; subsequent cooling iterations reused only
states belonging to that identical physical scenario.

### Figures and reproducible evidence

- [Baseline and perturbed torque and sensor errors](../../results/lge_sensitivity_20261011/analysis/baseline_vs_perturbed.png)
- [Sensitivity slopes and changes in absolute error](../../results/lge_sensitivity_20261011/analysis/sensitivities_and_error_changes.png)
- [Full accepted numerical records](../../results/lge_sensitivity_20261011/analysis/results.json)
- [Completion and transfer audit](../../results/lge_sensitivity_20261011/COMPLETION_AUDIT.json)
- [Final commands core bindings and runtimes](../../results/lge_sensitivity_20261011/PARALLEL_COMPLETION.json)
- [Independent remote verification and archive checksum](../../results/lge_sensitivity_20261011/POSTPROCESS_COMPLETION.json)
- [Retrieved raw runs inputs and original figures](../../results/lge_sensitivity_20261011/production)

Both final 300 dpi PNG figures were visually inspected. Filled markers and small
categorical x offsets separate overlapping series; offsets do not encode a
physical x quantity. Sensitivity torque changes are displayed in mN m, while
raw results remain in SI. The figure captions identify hypothetical inputs,
the lack of measured uncertainty bars, and TC2's input dependence. The original
remote plots and the exact final plotting code/style are retained separately.


## Before and after temperature and viscosity

Matched same-case baseline and perturbed fields are aligned on the original
720 by 192 theta-z grid, with no spatial interpolation. Means below are
equal-area means on that grid, not volume- or mass-weighted averages. Film
temperature is the depth-averaged temperature field. Viscosity is the exported
local dynamic viscosity field. Changes include the reclosed equilibrium
position; the viscosity law itself is unchanged. These are the completed
rotating-shaft cases with the user-confirmed motion convention above.

| Case and perturbation | Film T before (C) | After (C) | Change (K) | Mean viscosity before (mPa s) | After (mPa s) | Change (%) | Spatial comparison |
|---|---|---|---|---|---|---|---|
| E Supply +1 K | 59.049 | 59.214 | +0.165 | 30.928 | 30.724 | -0.662 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/E_supply_plus_1K_fields.png) |
| E Baths +1 K | 59.049 | 59.604 | +0.555 | 30.928 | 30.248 | -2.199 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/E_cooling_plus_1K_fields.png) |
| A TC2 supply | 31.588 | 32.250 | +0.662 | 106.279 | 102.713 | -3.355 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/A_supply_TC2_fields.png) |
| B TC2 supply | 65.529 | 67.775 | +2.245 | 23.866 | 21.800 | -8.656 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/B_supply_TC2_fields.png) |
| C TC2 supply | 30.235 | 30.735 | +0.500 | 113.961 | 111.046 | -2.558 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/C_supply_TC2_fields.png) |
| D TC2 supply | 29.895 | 30.036 | +0.141 | 116.227 | 115.122 | -0.950 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/D_supply_TC2_fields.png) |
| E TC2 supply | 59.049 | 59.618 | +0.569 | 30.928 | 30.229 | -2.260 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/E_supply_TC2_fields.png) |
| F TC2 supply | 30.025 | 30.092 | +0.067 | 115.189 | 114.787 | -0.349 | [Before, after and difference](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/F_supply_TC2_fields.png) |

[All-case temperature and viscosity overview](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/temperature_viscosity_before_after.png) compares both states and highlights absolute temperature and relative viscosity changes. [Sensor temperatures and their changes](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/sensor_temperature_changes.png) separately compare TC1-TC3.

Each spatial comparison places before, after and difference from left to right,
with temperature above viscosity. Before/after colour limits match within each
row; difference colour limits are symmetric about zero. Different scenarios
have separate difference scales, so compare the colourbar numbers. Viscosity
difference maps show pointwise relative change, while the overview compares
changes in the area means. Full-precision numbers and input hashes are in
[the comparison data](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/summary.json)
and [field provenance](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/provenance.json).


## Revised contours and experimental errors, October 12

Widened all eight contour comparisons and reduced the panel gaps. The panels
share theta and z coordinates, show tick labels only on the outside, and share
one before/after colourbar per row. Temperature retains inferno, viscosity
retains viridis, and differences retain the red-blue scale centred at zero.
The earlier layout is preserved in the analysis directory. Numerical values
and before/after colour limits are unchanged.

[Experimental errors for cases A-F](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/experimental_errors.png)
show signed torque and TC1-TC3 errors before and after the selected TC2-supply
perturbation for A-F. The user requested the last E row only; its supply +1 K
and baths +1 K trials remain in the full table below and the archived
`error_layout_eight_comparisons` figure. Repeated scenario text is omitted from
row labels. Each sensor has its own framed Before/After block, with TC1, TC2
or TC3 written once below the pair.
Torque uses the corrected `new_exp_torque_nm`: `100 * (CFD / experiment - 1)`.
Sensor errors are `CFD - experiment` in K. The figure follows the user-provided
reference: paired horizontal torque-error bars on the left, with an annotated
temperature-error heatmap on the right and aligned case rows. Light blue is
before and dark blue is after. Each sensor has paired before/after heatmap
columns on a common Oranges scale of 0-5.5 K. Daggers mark TC2-input residuals
that are not independent validation. Every plotted value was recomputed from
the experimental measurements and checked against the accepted results.
Zero means agreement. The earlier scatter layout is preserved alongside the
analysis in `error_layout_before_reference`.

| Case and perturbation | Torque error before -> after (%) | TC1 error before -> after (K) | TC2 error before -> after (K) | TC3 error before -> after (K) |
|---|---:|---:|---:|---:|
| E Supply +1 K | +37.98 -> +37.10 | +4.119 -> +4.253 | +4.092 -> +4.223 | +4.821 -> +4.951 |
| E Baths +1 K | +37.98 -> +35.03 | +4.119 -> +4.713 | +4.092 -> +4.688 | +4.821 -> +5.418 |
| A TC2 supply | -27.95 -> -29.30 | +0.207 -> +0.806 | +0.164 -> +0.769 (input) | +0.357 -> +0.965 |
| B TC2 supply | -33.11 -> -36.90 | +0.237 -> +2.077 | +0.301 -> +2.111 (input) | +1.913 -> +3.719 |
| C TC2 supply | -45.01 -> -45.95 | +0.046 -> +0.501 | +0.016 -> +0.474 (input) | +0.157 -> +0.619 |
| D TC2 supply | -8.01 -> -8.83 | +0.481 -> +0.610 | +0.483 -> +0.612 (input) | +0.550 -> +0.679 |
| E TC2 supply | +37.98 -> +34.98 | +4.119 -> +4.581 | +4.092 -> +4.542 (input) | +4.821 -> +5.269 |
| F TC2 supply | +31.50 -> +31.05 | +0.674 -> +0.736 | +0.678 -> +0.740 (input) | +0.719 -> +0.781 |

TC2 errors marked input are not independent validation because the measured
TC2 temperature was used as the assumed inlet. No measured spatial temperature
field or film viscosity is available, so no experimental contour or viscosity
error is inferred. Measurement uncertainty is not supplied in the comparison
table. [Exact errors and source hashes](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/experimental_errors.json)
and [the plotting script](../../results/lge_sensitivity_20261011/analysis/thermal_viscosity/compare_experiment_errors.py)
are retained with the figures.
