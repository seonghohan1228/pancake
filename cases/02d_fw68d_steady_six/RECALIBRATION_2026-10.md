# Measurement recalibration for the A–F cases (October 2026)

Two experimental inputs of the six FW68D cases were re-measured in October
2026. On 2026-10-23 Seongho asked to apply both and rerun the CFD. This note
records what changed, why, and what it does to the inputs and the comparison.
Original values stay in `experiment.json`; the recalibrated set is
`experiment_recalibrated_20261023.json`, built by
`recalibration_2026_10/make_recalibrated_experiment.py`.

| Quantity | Enters CFD? | Before | After |
|---|---|---|---|
| Radial load $W$ | Yes (`external_load_y`) | $F = 2752.246\,(V - 0.0130312)$ N | $F = 2068.1\,V$ N (bounds 1844.7–2291.5 N/V) |
| Parasitic torque $M_p$ | No (comparison target) | Sep 11 sweep, test bearing installed dry | Oct 21 sweep, test bearing removed; up branch primary, down branch kept |

## 1. Radial load (load-cell recalibration, 2026-10-20)

**Measurement.** Record `live-daq_20261020_213732_047904700.csv`. The load
cylinder was stepped from 1 to 9.4 bar abs and back with the rod side at
1 bar abs; the final 0 V plateau is the no-contact zero. Force on the load cell
is the cylinder force $A\,\Delta p$ ($A$ for a 100 mm bore) minus seal friction,
which produces a 300–900 N hysteresis loop. The loading branch bounds the
calibration from above, the unloading branch from below. Bracketing gives

$$F = K V_{CTT},\qquad K = 2068.1\ \mathrm{N/V}\ (1844.7\text{–}2291.5\ \mathrm{N/V}),\qquad V_0 = 0\ \mathrm{V}.$$

The provisional law ($2752.2$ N/V, from one 5 bar point on 2 September) lies
above the loading-branch upper bound.

**Application.** Each saved force is converted back to its CTT voltage with the
provisional law and then recalibrated:
$F_{new} = K\,(F_{old}/2752.246 + 0.0130312)$.

| Case | rpm | $W$ provisional (N) | $W$ recalibrated (N) | Bounds (N) | Ratio |
|---|---:|---:|---:|---|---:|
| A | 200 | 5276.3 | 3991.7 | 3560.6–4422.9 | 0.757 |
| B | 1500 | 3237.8 | 2460.0 | 2194.2–2725.7 | 0.760 |
| C | 200 | 1973.0 | 1509.5 | 1346.5–1672.6 | 0.765 |
| D | 200 | 382.5 | 314.4 | 280.4–348.4 | 0.822 |
| E | 1500 | 461.3 | 373.6 | 333.2–413.9 | 0.810 |
| F | 200 | 190.6 | 170.2 | 151.8–188.6 | 0.893 |

**Limits.** The cylinder part number is not confirmed (100 mm bore assumed).
Seal friction is bracketed, not removed; the bounds column is the resulting
range. A known-mass check of the load cell has not been done.

Evidence: `recalibration_2026_10/load_cell_calibration_20261020.json` and
`..._analysis.py` (copied from the research vault,
`tmp/lge_recal_20261020/`).

## 2. Parasitic torque (re-measured 2026-10-21)

**Measurement.** Record `live-daq_20261021_204036_717183600.csv`
(SHA-256 `8fbce5ce…1f075`), test bearing removed, dry, unloaded, atmospheric,
no DAQ sensors. Sweep 100→2000→100 rpm in 100 rpm steps. Hold means use samples
≥30 s after each step within 3 % of the target speed.

| rpm | Sep 11, bearing in, dry (N·m) | Oct 21 up (N·m) | Oct 21 down (N·m) |
|---:|---:|---:|---:|
| 200 | 0.643 | 0.438 | 0.326 |
| 1500 | 0.650 | 0.482 | 0.444 |

Full table: `recalibration_2026_10/parasitic_torque_20261021_by_rpm.csv`.

**Findings.** The new reference is 0.15–0.40 N·m lower than Sep 11; the removed
dry test-bearing rub accounts for most of that. Below 1300 rpm the down branch
is 0.10–0.17 N·m lower than the up branch (one-way state change near
1400 rpm, likely seal/grease warming). The up branch is the primary reference
for the long steady runs.

**Limits.** The reference is dry, unloaded and atmospheric. With it, the 200 rpm
high-$S$ residual torque is 2–4× Petroff, so the test runs carry extra losses
the reference does not contain (support bearings under load, wetted seals in
cold oil, oil drag). The torque comparison remains conditional on this.

## 3. What changes in the CFD and the comparison

- CFD input: only `load_mean_n` (hence `external_load_y`). Geometry,
  viscosity, temperatures, boundaries and solver settings are unchanged.
- Experimental Sommerfeld number rises by $1/\text{ratio}$ (×1.12–1.32).
- Experimental bearing torque: $M_b = M_{motor} - M_p$ with the Oct 21 up
  branch (`apparent_bearing_torque_20261021_up_nm`); down branch and Sep 11
  values are kept for comparison.
- Friction coefficient uses the recalibrated load.

## Rerun record

Campaign 2026-10-23, `run_steady.py --experiment experiment_recalibrated_20261023.json`,
solver built from branch `codex/lge-recalibrated-inputs` (source identical to
`main` `acc20fd`), Linux, GCC 13, OpenMPI, PETSc 3.19.6, 2 MPI ranks.
Results are summarized in `RECALIBRATION_RESULTS_2026-10-23.md` when complete.
