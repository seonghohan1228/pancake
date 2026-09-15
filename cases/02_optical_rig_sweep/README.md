# Case 02 — Optical rig operating-window sweep (Rev B geometry)

Maps the cavitation operating window of the PMMA optical journal-bearing rig so
that the servo speed, camera stations, transducer ranges and the P1 test point
can be fixed before machining. Supersedes the discarded 96-case sweep, which was
built on the previous 50.0 mm journal / 250 µm clearance geometry.

## Geometry and physics (frozen 2026-08-11)

| Quantity | Value |
|---|---|
| Journal diameter `D` | 49.300 mm (`R` = 24.65 mm) |
| Sleeve bore | 50.000 mm, PMMA, honed |
| Radial clearance `c` | 350 µm (RUN 2 contrast at 250 µm) |
| Relative clearance `ψ = c/R` | 1.42 % |
| Bearing length `L` | **50.0 mm (`L/D` = 1.014)** |
| Supply | none — ends-fed only, both ends 1 bar abs |
| Oil | paraffin, ρ = 850 kg/m³, µ(20 °C) = 0.068 Pa·s (grade estimate) |
| Temperature | isothermal 20 °C (THD suite checks this) |
| Cavitation | FBNS complementarity, `SCALAR_PCAV` |
| Grid | G6, 720 × 192 |

## Angle convention

`attitude_angle_deg = 0` throughout, so with the solver's film law
`h = c − e cos θ` ([film_thickness.cpp:27](../../src/film_thickness.cpp:27)):

- **`φ = 0` is the line of centres**, i.e. the minimum-film line.
- `φ` increases in the direction of journal rotation (`ω > 0`).
- The film diverges over `0 < φ < 180°`, so the **cavity lives there**; it
  converges over `180° < φ < 360°`, where the **pressure peak** sits.
- `attitude_deg` in the CSV is `(force_angle_deg + 180) mod 360`: the angle from
  the line of centres to the load line, positive in the direction of rotation.

Case 01e reports angles offset by −90° because its supply groove sits at θ = 90°.
This case has no groove, so no offset is applied. **Convert to rig azimuth before
drilling camera ports** — the mapping from `φ = 0` to a physical clock position
depends on how the loaded journal is oriented in the housing, which is a rig
drawing decision, not a solver output.

## Suites

| Suite | Points | Content |
|---|---|---|
| `RUN1` | 60 | ε ∈ {0.5, 0.6, 0.7} × N ∈ {500, 1000, 1500, 2000, 3000} rpm × p_sat ∈ {degassed, 0.2, 0.5, 0.8} bar abs |
| `RUN2` | 8 | same points at `c` = 250 µm, ε ∈ {0.5, 0.7}, N ∈ {1000, 2000}, p_sat ∈ {degassed, 0.8} |
| `RUN4` | 4 | back-pressure lever: both ends at 1.5 / 2.0 bar abs, p_sat ∈ {0.8, 1.5}, ε = 0.7, 2000 rpm |
| `BAND` | 60 | ε ± 0.04 uncertainty band on the RUN1 grid at p_sat ∈ {degassed, 0.8} |
| `THD` | 4 | worst corner (ε = 0.7, 3000 rpm) with the energy equation and Andrade µ(T) |

```bash
python3 cases/02_optical_rig_sweep/run_sweep.py RUN1 RUN2 RUN4 THD
```

Then `analyse_sweep.py` (writes `sweep_results.csv`, `sweep_analysis.json` with
the Q1–Q5 answers) and `plot_sweep.py`.

Runs are resumable: `sweep_runs.json` keys every completed point and re-running a
suite skips them.

## Decisions taken, and why

**RUN 3 (L/D variation) is removed.** It existed to test whether raising L/D to
0.75 could recover cavitation signal lost to the larger clearance. With the rig
now at L/D = 1, 0.75 is a step *down*; the lever is already spent. If a length
study is still wanted it should sweep *downward* from 1.0, which is a different
question — flag it and it can be added as a new suite.

**`LOCAL_PSAT` replaced by `SCALAR_PCAV`.** The original spec asked for
`cavitation_threshold = LOCAL_PSAT` with p_sat prescribed per case. In this
solver `p_sat` is not an input: `effective_p_cav()`
([config.hpp:1957](../../src/config.hpp:1957)) derives it from a Henry-law
solubility model, and under `fluid_property_model = CONSTANT` it returns 0 and
falls back to `p_cav` — all four "p_sat" cases would have been identical.
`LOCAL_PSAT` is only `max(p_cav, p_sat)`, a uniform datum, so setting `p_cav`
directly is the same physics with no silent-failure mode. "degassed" is
`p_cav` = 500 Pa (0.005 bar abs), matching the vaporous floor used in the
exclusion argument; the config rejects `p_cav ≤ 0`.

**Pinned ε, not pinned load.** Every case-01e production run solves a force
balance to a target load. Here ε is prescribed and `W` is read out, so each point
is a single steady solve. The case-01e `force_balance.py` machinery is therefore
not reused, only `solver_io`.

**Viscosity.** The lab oil is white mineral oil, CAS 8020-83-5, ISO VG 46. µ is
now a **grade-based conversion, not a fit**: ν(40 °C) = 46 cSt by grade
definition, so ν(20 °C) ≈ 80 cSt at VI ≈ 100, and µ = ρν = 850 × 80e-6 =
**0.068 Pa·s**. Density 0.85 g/mL at 20 °C is confirmed from the datasheet.

`ANDRADE_E` = 3556 K is still the old two-point fit to a different oil and is
used only by the `THD` suite; with no measured pair for this oil, **THD results
are qualitative only**.

A grade label is not a measurement, and the `MU` suite shows why that matters:
across µ ∈ {0.05, 0.068, 0.10, 0.15} the load and peak pressure scale linearly,
but the degassed cavity arc runs 25° → 103°, a factor of four. Measuring the
oil remains the highest-value input refinement available.

**Isothermal main matrix.** The `THD` suite runs the worst corner with the energy
equation and Andrade µ(T), adiabatic and cooled-journal bracketing, to quantify
what the isothermal assumption costs. If the µ shift there is material the whole
matrix should be re-run thermally.

## Open items

- **Lengthening the bearing costs transducer headroom.** A 180 × 48 probe at the
  worst corner (ε = 0.7, 3000 rpm, µ = 0.075) gives W = 679 N, p_max = 7.28 bar abs,
  p_min pinned at the 500 Pa floor with a 112° cavity, and 160 W dissipated. Shaft
  deflection is then ~6 % of `c` and PMMA bore expansion ~4 % of `c`, both
  tolerable, but p_max sits just above the planned 0–7 bar channel — size the
  high-range transducer from `sweep_analysis.json → transducer` once RUN1 lands.
- **The degassed exclusion argument does not hold everywhere.** At that corner the
  degassed case cavitates on its own, so it is not a clean control. Q3 in the
  analysis and `rig_degassed_exclusion_margin.png` map where the margin survives;
  the P1 point must be chosen inside that region.
- **µ(20 °C) = 0.068 Pa·s (grade estimate) is fitted, not measured** — no datasheet available for
  the rig oil. Film pressures and `W` scale linearly with µ, but the cavity
  geometry does **not**, because the threshold is an absolute pressure: a different
  µ moves the whole pressure field relative to a fixed `p_cav`. If the real oil
  turns out more than ~20 % away from 0.075 Pa·s the matrix has to be re-run,
  which costs a couple of hours. Measuring the oil, or back-calculating µ from the
  rig's own friction torque on first spin-up, is worth doing early.
- **Resolved 2026-08-12: ε band.** Now ±0.02, from the measured stack-up
  sqrt(3² + 5² + 3²) / 350 = ±0.019 (shim ±3 µm, bore ±5 µm, touch-zero ±3 µm).
  The earlier ±0.04 was twice the real resolution.
- **Resolved 2026-08-12: grooveless ends-fed anchor (V1).** `verify_setup.py`
  compares ε = 0.5 / 1000 rpm against the Ocvirk short-bearing solution and
  passes; see `verify_setup.json`. The solver sits at ~0.47 of Ocvirk because
  Ocvirk's load integral keeps the sub-ambient lobe that JFO truncates.
- **Resolved 2026-08-12: grid independence (V2).** The worst corner at
  720 × 192 against 1080 × 288 agrees to 0.004 % in load, 0.003 % in peak
  pressure and 0.157 % in cavity arc, against a 2 % gate. G6 transfers to this
  geometry.
- **Rigid-bore correction is post-processing only.** `analyse_sweep.py` applies
  the quoted 23 µm/MPa PMMA compliance to `p_max` and reports it as a percentage
  of `c`; the solver itself assumes a rigid bore. Once that percentage stops being
  small, the assumption — not the correction — is what fails.
