# Case 05 — Etsion & Ludwig 1982: cavity EXTENT in a submerged journal bearing

**Scope: EXTENT ONLY.** This case scores the circumferential start and end of the
cavitated zone against Table 2 of Etsion & Ludwig (1982). It does **not** attempt
the cavity-*pressure* comparison of their Figs. 6/7 — that is case 05b in
`docs/UNFINISHED_BUSINESS.md` §3.9 and is blocked on the `p_void` convergence work
(§4.3 there). See "What this case deliberately does not claim" below.

Written before the runs, per `cases/VALIDATION_PLAN.md` rules.

---

## 1. Reference

Etsion, I. & Ludwig, L. P., *Observation of Pressure Variation in the Cavitation
Region of Submerged Journal Bearings*, ASME J. Lubrication Technology **104**(2),
April 1982, pp. 157–163 (ASLE-ASME Joint Lubrication Conf., New Orleans, Oct 1981,
Paper 81-Lub-11).

PDF: `F:\Seongho\Research\literature\pdfs\Etsion and Ludwig - 1982 - Observation of
Pressure Variation in the Cavitation Region of Submerged Journal Bearings.pdf`

Transcribed reference table: `cases/ref_data/etsion1982_table2_extent.csv`
(**transcribed from the text layer, not digitised**; the column-to-value assignment
was cross-checked against seven explicit statements in the paper's Results section
— the checks are listed in the CSV header).

**Why this paper.** The cavity contents were pulled into an evacuated tube and run
through a mass spectrometer: **air, no oil vapour**. The oil's own vapour pressure
is 13.3 Pa, three orders below every measured cavity pressure. That is in-paper,
measurement-based proof that the mechanism is **gaseous (outgassing) cavitation**,
which is the regime `docs/CAVITATION_OIL_REFRIGERANT.md` argues for by analogy.
And the ambient pressure `p_s` is the **control variable** — a direct dissolved-air
lever, swept over four levels at each of two speeds.

---

## 2. Parameter table and verification status

`VERIFIED` = printed in the paper. `DERIVED` = arithmetic on verified values.
`UNVERIFIED` = not in the paper; a stated default with its justification.

### 2.1 Geometry — all VERIFIED (load-bearing)

| Parameter | Value | Config | Status | Source |
|---|---|---|---|---|
| Journal radius R | 25.4 mm | `R = 0.0254` | **VERIFIED** | "A rotor 5.08 cm dia" (p.157) |
| Bearing length L | 38.1 mm | `L = 0.0381` | **VERIFIED** | "3.81 cm long (L/D = 0.75)" |
| L/D | 0.75 | — | **VERIFIED** | stated |
| Mean radial clearance c | 114 µm | `c = 1.14e-4` | **VERIFIED** | "the average radial clearance was 0.0114 cm" |
| Eccentricity ratio ε | 0.4, **fixed, static** | `e = 4.56e-5` | **VERIFIED** | "an eccentricity of 0.4 was used for all the data reported herein"; set by shims |
| h_min | 68.4 µm | — | **VERIFIED** | paper states 0.00684 cm — independently confirms ε = 0.4 |
| Misalignment / tilt | not modelled | `tilt_slope_* = 0` | **VERIFIED as a known deviation** | paper attributes the slight axial asymmetry of its traces to imperfect alignment; magnitude not quantified, so it is not modelled |
| Journal runout | < 12.7 µm | not modelled | **VERIFIED** | "runouts of less than 0.00127 cm" |
| Grooves / inlets | **none** | no `inlet_groove` | **VERIFIED** | oil in/out ports on **both end faces**; no circumferential groove anywhere in the description |

### 2.2 Operating conditions — all VERIFIED (load-bearing)

| Parameter | Value | Status | Source |
|---|---|---|---|
| Speeds | **1840 and 3000 rpm** | **VERIFIED** | "Tests were conducted at two different shaft speeds of 1840 and 3000 rpm" |
| ω | 192.68281 / 314.15927 rad/s | DERIVED | 2πN/60 |
| Ambient levels p_s | **13.6 / 27.2 / 40.8 / 54.4 kPa GAUGE** | **VERIFIED — the full set, as instructed** | "four levels of the supply pressure, p_s, were examined ranging from 13.6 KPa (2 psig) to 54.4 KPa (8 psig) at intervals of 13.6 KPa (2 psi)"; the concluding remarks repeat the same range; individual panels of Figs. 5–7 are labelled 13.6, 27.2, 40.8, 54.4 |
| p_s reference | **gauge** (differential transducer referenced to atmosphere) | **VERIFIED** | "The differential pressure was measured with respect to the atmospheric pressure"; "2 psig"…"8 psig" |
| p_s absolute (used for the axial Dirichlet ends) | 114.925 / 128.525 / 142.125 / 155.725 kPa | DERIVED | p_s + 101.325 kPa |
| Both ends flooded at the same p_s | yes | **VERIFIED** | "Both ends of the bearing were flooded with oil… oil in and oil out ports on both ends"; nomenclature: "p_s = ambient supply pressure on both sides of bearing" |
| Oil-in temperature | 30–38 °C, controlled | **VERIFIED** | heat exchanger; run to thermal equilibrium before every scan. Fig. 7(b) is annotated T = 36 °C |
| Temperature used here | **38 °C = 311.15 K, isothermal** | DERIVED choice | the only temperature at which the paper's viscosity is tabulated; the rig is temperature-controlled and every scan is taken at thermal equilibrium, so an isothermal run at the property temperature is the consistent choice. Using 38 °C rather than 30 °C is the conservative end (lowest µ ⇒ weakest hydrodynamics) |

### 2.3 Fluid — VERIFIED where the paper prints it (load-bearing)

| Parameter | Value | Config | Status | Source |
|---|---|---|---|---|
| Oil | paraffinic, **unnamed** | — | **VERIFIED (that it is unnamed)** | "A paraffinic oil was used as the test fluid" |
| Dynamic viscosity µ(38 °C) | **0.026 Pa·s** | `mu = 0.026` | **VERIFIED** | Table 1 |
| µ at other T | 2.1 (−29 °C), 0.84 (−18 °C), 0.0047 (99 °C), 0.0017 (204 °C) Pa·s | not used (isothermal) | **VERIFIED** | Table 1 — enough to fit µ(T) if a THD leg is ever added |
| Specific gravity | 0.8285 → ρ = **828.5 kg/m³** | `rho = 828.5` | **VERIFIED** | Table 1 |
| Oil vapour pressure at 38 °C | 0.1 mmHg = **13.3 Pa** | not used | **VERIFIED** | Table 1. Three orders below every measured cavity pressure ⇒ vaporous cavitation is ruled out by the reference itself |
| Bulk modulus β | **1.7 GPa assumed** | `bulk_modulus = 1.7e9` | **UNVERIFIED — not in the paper** | literature value for a paraffinic mineral oil near 40 °C / 1 bar. **Used only by the JFO/Elrod leg** — FBNS does not read `bulk_modulus` at all (`PHYSICS.md` §7.3.8), so the primary (FBNS) result carries no dependence on this number. A β sensitivity leg (1.0 vs 1.7 GPa) is run on one JFO condition to bound it. Flagged for user sign-off per the case-00 rule. |
| Specific heat / conductivity | not needed | — | n/a | isothermal case |

### 2.4 Dissolved-air state — UNVERIFIED, stated default (this is the modelling choice)

The paper does **not** report the dissolved-air content of its oil. Per the task
brief, the default is **Henry-law air saturation**, with the Bunsen/Ostwald
coefficient the campaign already records for air in mineral oil.

| Parameter | Value | Config | Status |
|---|---|---|---|
| Species | air | `dissolved_gas_species = AIR` | **VERIFIED by the paper's own mass spectrometry** |
| Solubility law | Henry, c_sat = H(T)·p | `solubility_model = HENRY` | UNVERIFIED (standard closure; `CAVITATION_OIL_REFRIGERANT.md` §5.1) |
| Ostwald/Bunsen coefficient L | **0.09** (≈ 9 vol% air at 1 atm) | — | UNVERIFIED — the classical 8–12 vol% figure for mineral oil (Pinkus 1990, via [R1] in `docs/CAVITATION_OIL_REFRIGERANT.md`); the same 0.09 the campaign records in `UNFINISHED_BUSINESS.md` §3.9 (05c) |
| Henry slope H | **1.2164e-9 1/Pa** | `dissolved_gas_henry_coeff` | DERIVED: H = L·ρ_gas(1 atm, 38 °C)/(ρ_liq·p) so that `ostwald_coefficient()` returns exactly 0.09 |
| Feed dissolved mass fraction c_d | **1.23252e-4** | `dissolved_gas_initial` | DERIVED: H × 101325 Pa = air saturation at **atmospheric** exposure |
| Free-gas viscosity µ_g | 1.9e-5 Pa·s (air, 38 °C) | `mu_gas = 1.9e-5` | UNVERIFIED (handbook) |
| Two-phase µ closure | McAdams | `mixture_viscosity_model = MCADAMS` | UNVERIFIED — open question §7.3 of `CAVITATION_OIL_REFRIGERANT.md`; kept identical to the campaign's other configs so this case introduces no new free parameter |

### 2.5 THE cavitation datum, and why it is what it is

This is the single most consequential choice in the case, so it is stated in full.

**Chosen: `p_cav` = 101 325 Pa ABSOLUTE**, identical for all eight conditions,
realised as `p_cav = 101325` with `cavitation_threshold = LOCAL_PSAT` and
`dissolved_gas_initial` set to Henry saturation at 1 atm, so that
`effective_p_cav = max(p_cav, p_sat(c_d,init, T_init)) = 101 325 Pa` exactly
(`SimulationConfig::effective_p_cav`, `PHYSICS.md` §7.3.8).

**Three reasons.**

1. **The paper postulates exactly this.** Its closing mechanism paragraph:
   the oil "is saturated with air under atmospheric pressure" and "becomes
   oversaturated in the subcavity pressure loop". The air-release datum is
   therefore atmospheric saturation, not saturation at the elevated end pressure.

2. **It is the Henry-law default read correctly.** Air enters solution where the
   oil meets air — the reservoir, at atmosphere. The 13.6–54.4 kPa gauge applied at
   the bearing ends is a pump/seal pressure downstream of that exposure, and gas
   solubility does not equilibrate to it on the residence time of a once-through
   feed line.

3. **It is the only choice under which p_s is a physical variable rather than a
   datum shift.** With constant properties and Dirichlet ends, the incompressible
   Reynolds/complementarity system is invariant under
   `p → p + Δ`, `p_cav → p_cav + Δ`, `p_bc → p_bc + Δ`. Setting the datum to
   "Henry saturation at each ambient" would make `p_cav = p_ambient,abs` and
   would therefore produce the **identical cavity at all four ambient levels** —
   destroying the very sweep this case exists to score. That degeneracy is
   demonstrated numerically, not asserted: run leg `sens_pcav_tracks_ambient`
   repeats two conditions with the tracking datum (see §4.3).

**Tier-0 limitation, stated prominently.** With a single scalar `p_cav`, the model
**cannot** reproduce the measured minimum film pressure. Etsion & Ludwig measure
p_min = −13.6 kPa gauge at 1840 rpm and −27.2 to −34.0 kPa gauge at 3000 rpm, i.e.
**67.3 to 87.7 kPa absolute — below the 101.3 kPa release datum.** Complementarity
forbids p < p_cav anywhere, so every pancake run will pin its minimum at exactly
101 325 Pa and will report θ₀ (the location of the minimum) as the rupture point
itself. The measured 5–13° gap between θ₀ and θ_st, and the 14–34 kPa of
sub-datum tension, are a **finite-rate supersaturation effect** (Tier 2 in
`CAVITATION_OIL_REFRIGERANT.md` §6) that the paper's own authors name and that no
equilibrium model of any kind can produce. This case therefore scores **θ_st and
θ_end**, reports the p_min failure as a known and expected Tier-0 outcome, and does
not tune anything to hide it.

### 2.6 Angle convention

Etsion's θ is measured **from the line of maximum clearance, positive in the
direction of shaft rotation** (their Nomenclature and Fig. 4). pancake's static
film is `h = c − e·cos(θ − ψ)` with `ψ = attitude_angle_deg`, so h is minimum at
θ = ψ and maximum at θ = ψ + 180°, and ω > 0 drives Couette flow in +θ.

**Setting `attitude_angle_deg = 180` makes pancake's θ identical to Etsion's θ**:
max clearance at 0°, min clearance at 180°, converging arc 0–180°, diverging arc
180–360°. No post-hoc angle offset is applied anywhere in this case. Under `STATIC`
motion the attitude angle is an input that fixes the frame; the load magnitude and
direction are outputs.

**The reference's own angular uncertainty is ±5° on every angle**, systematic and
common to all four columns (the paper states it twice). Agreement is scored
against that band.

### 2.7 Discretisation — frozen by case 00

`n_theta_global = 720` (0.5°/cell), `n_z_global = 192` (0.198 mm/cell),
`solution_mode = STEADY_STATE`. Case 00 froze G6 = 720×192 and dt = 2e-5 s for
every subsequent case; this is a steady case, so only the grid applies. 0.5° cells
are ten times finer than the reference's 5° scan increment.

---

## 3. Protocol

- `motion_model = STATIC`, ε = 0.4 prescribed, load is an output.
- `solution_mode = STEADY_STATE`; `steady_gas_model = EQUILIBRIUM` (the equilibrium
  flash; `FROZEN` is refused by the validator in steady state and would be wrong).
- `temperature_model = ISOTHERMAL` at 311.15 K — the rig is temperature-controlled
  and every scan is taken after thermal equilibrium is reached.
- **Submerged/flooded**: `bc_z_south_type = bc_z_north_type = DIRICHLET` at the
  **absolute** ambient pressure of that level. No inlet groove.
- **The case-02 no-inlet warning is checked, not assumed.** Case 02 found that with
  no inlet, mass-conserving FBNS settles into a *starved* branch (cavitated
  fraction 0.76–0.86 at every ε, ~1/20 of the flooded load). Etsion's rig is
  genuinely fed through both end faces, so the axial Dirichlet ends *are* the feed
  and the flooded branch should be reached — but this is verified explicitly before
  any extent is scored (§4.2). **Documented fallback if it starves:** add a 1°
  ambient-pressure inlet line at maximum clearance (θ = 0), which is the flooded
  framing case 02 already used for its non-starved locus, and label every affected
  result with the deviation.

### 3.1 Run matrix (16 production runs)

2 speeds × 4 ambient levels × 2 cavitation models:

| leg | `cavitation_model` | `fluid_property_model` | runs |
|---|---|---|---|
| `fbns/` | `FBNS` | `TWO_PHASE` (air) | 8 |
| `jfo/` | `JFO` | `TWO_PHASE` (air) | 8 |

Both legs use `gas_pressure_coupling = NONE` (FBNS rejects `VOID_COUPLED` at
validation; JFO is kept on the same setting so the two legs differ in exactly one
thing — the cavitation solver).

**Why both models is a result, not redundancy.** FBNS and JFO treat the gas volume
differently by construction:

- **FBNS + TWO_PHASE** moves the full-film ceiling to θ_full = 1 − α_gas
  (`PHYSICS.md` §7.3.7), so released air is **volumetrically active**: it displaces
  liquid and the complementarity is written against the moved ceiling.
- **JFO/Elrod** keeps the ceiling at θ_full = 1 and reaches the pressure only
  through the mixture viscosity/density; with `gas_pressure_coupling = NONE` the
  released air is **volumetrically inert**.

At L = 0.09 the released air can occupy at most ≈9 % of the gap, so the two legs
should differ by a bounded, predictable amount. Whether they do — and whether that
difference points the right way relative to the measurement — is the cross-model
finding this case is asked to report.

### 3.2 Supporting legs (not part of the 16)

| leg | purpose |
|---|---|
| `smoke/` | one coarse 180×48 condition, local, to validate the config before dispatch |
| `gate/` | one production condition run **locally on 4 ranks** and on **afml-r6 on 24 ranks**; must match to 0.01 % (mandatory cross-platform/rank gate, `VALIDATION_PLAN.md`) |
| `sens_pcav_tracks_ambient/` | 2 conditions with `p_cav` = ambient absolute — demonstrates the datum degeneracy argued in §2.5 |
| `sens_beta/` | 1 JFO condition at β = 1.0e9 vs 1.7e9 — bounds the one UNVERIFIED number the JFO leg depends on |

**Added during the batch (2026-08-05), after the smoke run showed `alpha_gas ≡ 0`:**

| leg | purpose |
|---|---|
| `sens_shift_invariance/` | 1 condition with Δ = 40.8 kPa added to **both** the datum and both axial ends. The clean way to demonstrate §2.5's degeneracy claim, because `sens_pcav_tracks_ambient/` turned out not to be solvable (it puts the whole flooded boundary on the Fischer–Burmeister kink `a = b = 0`) |
| `vapour_fbns/`, `vapour_jfo/` | the same conditions with the datum at the **oil vapour pressure, 13.3 Pa** (Table 1 of the paper — no fitted number). Two purposes: the "what if this were vaporous cavitation" control, and the only steady configuration in which the `TWO_PHASE` outgassing is not inert. Under `LOCAL_PSAT` the datum is raised to `p_sat` by construction, complementarity then pins `p ≥ p_sat` everywhere, and the equilibrium flash releases nothing — so `cavitation_threshold = SCALAR_PCAV` is required to open a release band at all |

---

## 4. What is scored, and how

### 4.1 Metrics

| Metric | Reference | Solver definition |
|---|---|---|
| θ_st (rupture) | Table 2 | first θ, scanning downstream from minimum clearance, at which the film is cavitated |
| θ_end (reformation) | Table 2 | last θ of the same contiguous cavitated arc (wraps through 0°) |
| circumferential extent | derived, (θ_end − θ_st) mod 360 | same |
| θ₀ (min-pressure location) | Table 2 | argmin of the midplane pressure — **reported, expected to coincide with θ_st on any Tier-0 path** (§2.5) |
| p_min | Table 2 | **reported as a known failure**: pinned at p_cav by complementarity |
| axial cavity width at midplane | **NOT AVAILABLE** | Table 2 has four data columns and none of them is an axial width; Figs. 6/7 carry START/END marks on axial traces but no cavity width in z. **Scored: N/A.** The solver's axial width is still reported, unscored, as supporting material |

A cell counts as cavitated when `film_content < FullFilm::bound − tol`, i.e.
`theta < (1 − alpha_gas) − 1e-6` on the FBNS+TWO_PHASE leg and `theta < 1 − 1e-6`
on the JFO leg — the ceiling is model-dependent and must not be hard-coded to 1
(the trap `PHYSICS.md` §7.3.7 documents). Both a **midplane** (z = L/2) and an
**any-z** definition are extracted; the reference's START/END marks come from axial
traces, so the any-z definition is the closer match and is the one scored, with
midplane reported alongside.

### 4.2 Gate conditions before any extent is scored

1. **Flooded, not starved.** Global cavitated fraction must be well below the
   0.76–0.86 that case 02 saw on its starved no-inlet branch, and the load must be
   O(the flooded value), not ~1/20 of it.
2. **Converged.** FBNS Newton converged on every outer iteration; steady Picard
   residual below `coupling_tolerance`; liquid mass residual at machine-noise level.
3. **Complementarity invariants hold** (`PHYSICS.md` §7.3.7): min(p − p_cav) ≥ 0 and
   max(θ − (1 − α_gas)) ≤ 0 in the returned state.
4. **Rank/platform gate**: the `gate/` condition matches local 4-rank to remote
   24-rank within 0.01 %.

### 4.3 Figures

- `etsion_extent_scatter.png` — measured vs computed θ_st and θ_end, 8 conditions ×
  2 models, with the 1:1 line and the reference's ±5° band.
- `etsion_extent_vs_ps.png` — cavity extent and θ_end vs ambient pressure, one
  panel per speed, experiment and both models overlaid. This is where the paper's
  stated trend lives ("as the supply pressure is increased, the downstream end of
  the cavitation retreats toward the maximum clearance line").
- `etsion_cavity_map.png` — computed θ–z cavity boundary for the four 1840 rpm
  levels, i.e. the solver's analogue of their Figs. 5(a)–(d) photographs.

---

## 5. What this case deliberately does not claim

- **Not the cavity pressure.** The paper's headline result is that cavity pressure
  varies by ~50 kPa and exceeds p_s where the cavity penetrates the converging
  clearance. Every converged pancake path predicts a flat plateau at p_cav. That is
  a **refutation of the plateau assumption by the reference**, and it is case 05b's
  job to quantify it; here it is stated, not measured.
- **Not p_min.** See §2.5 — structurally out of reach for an equilibrium
  complementarity model, and the reference's own postulated mechanism
  (supersaturation) is Tier 2.
- **Not the reverse flow.** The paper reports a low-frequency backward/forward
  oscillation of oil streamers at the downstream cavity end. A 2-D Reynolds solver
  with a homogeneous cavitated film has no representation of it.
- **Not the axial cavity width.** Not in Table 2 (§4.1).
- **Not a claim about R290/PZ68S.** This is air-in-oil. It validates the *machinery*
  the refrigerant application will use, on the only regime for which a measured
  bearing dataset exists at all (`UNFINISHED_BUSINESS.md` §3.8).

---

## 6. Status

| | |
|---|---|
| Parameter sign-off | all load-bearing parameters (ambient levels, viscosity, geometry, speeds, ε, density) **VERIFIED from the paper** — no invented numbers, so the case is cleared to run under the task brief |
| Open for user sign-off | bulk modulus β = 1.7 GPa (JFO leg only, bounded by `sens_beta/`); Ostwald coefficient 0.09; `mixture_viscosity_model = MCADAMS` |
| Execution host | afml-r6, 24 pinned ranks (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`) |
| Results | `SUMMARY.md` (or `ETSION_REPORT.md` if the harness blocks that filename) |
