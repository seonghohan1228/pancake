# The Sommerfeld-collapse discriminator

Date: 2026-08-05. Branch `codex/validation-cavitation`. Analysis only: no solver
run, no source change. Every number below comes from data already archived in this
repository, recomputed by `collapse_analysis.py` in this directory.

Answers option (3) of `docs/UNFINISHED_BUSINESS.md` §2.1: *"Plot S vs ε for Ferron,
Brito and pancake on one axis. If the three collapse under a common clearance
rescale, the discrepancy is instrumentation; if they do not, it is physics."*

---

## 0. Verdict

**A single clearance rescale collapses Ferron completely and Brito only partly.**

| | collapse | best-fit `c_eff` | residual after collapse |
|---|---|---|---|
| **Ferron 1983** Fig. 5, 2–10 kN at 2000 rpm | **yes** | **160 ± 4 µm** (nominal 145; the figure's own second legend entry 152) | ±0.010 in ε, no load structure beyond the digitising systematic |
| **Brito 2006** Fig. 6, 10 kN at 1000–4000 rpm | **no** | **95.4 / 99.5 / 104.3 µm** at 1000 / 3000 / 4000 rpm (quoted 85.5) | −0.045 / −0.002 / +0.042 in ε, **monotone in speed** |

Read together:

1. **There is no load-dependent residual anywhere.** Over a 5× load range at fixed
   speed, one clearance number removes the whole of the Ferron discrepancy —
   +54 % to +85 % in load at the nominal 145 µm — and leaves a residual smaller
   than the digitising noise. This rules out a **load-dependent** error over the
   tested range. It does not exonerate the absolute load calculation: a uniform
   multiplicative load bias is exactly degenerate with a uniform clearance bias.
2. **There is a speed-dependent residual that no constant clearance removes.**
   Brito's required `c_eff` grows +9.3 % from 1000 to 4000 rpm at fixed load. The
   trend survives three independent reference curves, is replicated at a second
   supply pressure, and is 3.5σ against the reading uncertainty over the full
   range (1.8σ on the two converged points alone).
3. **That surviving residual is thermal, not hydrodynamic.** The pancake *THD*
   reference needs the smallest speed trend (+9.3 %); the isothermal and Raimondi
   & Boyd references need +14.7 % and +13.4 %. The energy equation therefore
   already removes about a third of it. What is left is a genuine gap, and its
   two candidate names are both already on the open list: no thermo-elastic
   clearance growth (`validation/VALIDATION_REPORT.md` §5.3) and an under-predicted
   film temperature (`docs/UNFINISHED_BUSINESS.md` §2.4).

**The honest limit of this result.** The *level* of `c_eff` is not certified by
this analysis and cannot be — see §8. A uniform multiplicative solver load error
is exactly degenerate with a uniform clearance error. What is certified is the
*absence of structure*: constant in load, not constant in speed.

---

## 1. Method — the clearance-free invariant

Both experiments report ε as a proximity-probe displacement `e` divided by an
**assumed** radial clearance. The assumed clearance therefore enters the
"measurement" itself, and it enters `S` as well:

```
S   = mu N L D (R/c)^2 / W          N in rev/s, D = 2R
eps = e / c
```

A rescale `c -> k c` maps `(S, eps) -> (S/k^2, eps/k)`, so a data point slides
along the parabola `S ∝ eps^2`. The label of that parabola,

```
Lambda := eps^2 / S = e^2 W / (mu N L D R^2),
```

contains **no clearance at all**. It is fixed by what the rig actually measured
(`e`, `W`) and by nominal constants. This gives an exact per-point effective
clearance with no fitting: invert the solver's own `Lambda_solver(eps) =
eps^2/S_solver(eps)` (monotone increasing) at the measured `Lambda`, obtaining
`eps*`, then

```
c_eff = e / eps* = c_assumed * (eps_measured / eps*).
```

**A single clearance rescale explains a dataset if and only if `c_eff` is constant
across that dataset's independent variable** — load for Ferron, speed for Brito.
The *scatter* of `c_eff` is the discriminator. Its *level* is not (§8).

Two properties worth stating because they remove obvious objections:

- **The viscosity convention cancels.** `mu`, `N`, `L`, `D`, `R` appear identically
  in `Lambda_measured` and in `Lambda_solver` whenever the solver run uses the same
  nominal values as the experiment, which it does here. `c_eff/c_assumed` is a
  pure function of the dimensionless data. Choosing `mu(40 °C)` versus any other
  reference temperature does not move it.
- **The rescale is multiplicative, not additive.** The Lubricants 11(7) 308 (2023)
  re-reduction is recorded in the Ferron CSV header as "reduce the published ε by
  0.10". An additive shift is *not* a clearance rescale. Applied correctly, as
  ×145/171 = 0.848, that re-reduction over-corrects (§6).

---

## 2. Inputs and provenance

| dataset | source in repo | notes |
|---|---|---|
| Ferron Fig. 5 experiment, 9 points, `C` = 145 µm basis | `cases/ref_data/ferron1983_fig5_2000rpm_eccentricity.csv` | header states ±0.01 cluster scatter; σ_ε = 0.015 used (scatter + digitising) |
| pancake THD, Ferron rig, `L/D` = 0.8, 2000 rpm, `c` = 152 µm | `cases/01_ferron_thd/SUMMARY.md` Part B (6 loads, `flow` pass, J-adiabatic) | every point closed \|W\| to <0.1 %, Picard converged |
| Brito Fig. 6 experiment, 10 kN, `Pf` = 70 and 210 kPa, `c` = 85.5 µm basis | `validation/journal_bearing/brito.py`, `ECC_ATTITUDE` and `ECC_ATTITUDE_INDEPENDENT` | mean of the two independent digitisations (commit `81dd953`); σ_ε = 0.02 |
| pancake THD, Brito rig, ε carrying the applied 10 kN | `validation/journal_bearing/brito_eps_for_load.json` | **1000 rpm row is not converged** (coupling capped at 80, `R_p` = 5.23) |
| pancake isothermal plain-bearing loci, `L/D` = 1 and 0.5 | `cases/02_static_locus/locus_LD1.csv`, `locus_LD05.csv` | FBNS, flooded, 1° ambient groove at `h_max` |
| Raimondi & Boyd 1958 design data, `L/D` = 1 and 0.5 | `validation/journal_bearing/raimondi_boyd.py` `RB_TABLE` | provenance in §4 |

Both rigs are `R` = 50 mm, `L` = 80 mm, `L/D` = 0.8 — the LMS-Poitiers lineage
(`brito.py` records that Brito states it is the rig used by Ferron, Frene and
Boncompain). The `L/D` = 1 and 0.5 loci are interpolated to `L/D` = 0.8 linearly
in `ln S` and `ln(L/D)`.

---

## 3. Method audit

**(a) Reynolds similarity, tested on the solver's own runs.** The method assumes
the solver's `eps(S)` relation is clearance-independent. That is exact for the
isothermal Reynolds equation and only approximate under THD, because the
temperature rise itself depends on `c`. The Brito set contains a direct test: the
same operating points were solved at `c` = 85.5 µm and at the +8.8 % hot-clearance
bracket `c` = 93.02 µm. Predicting the hot runs from the 20 °C curve alone:

| N (rpm) | ε, hot run | ε, from the 20 °C curve | Δε |
|---|---|---|---|
| 1000 | 0.6783 | 0.6781 | +0.0002 |
| 3000 | 0.4607 | 0.4633 | −0.0026 |
| 4000 | 0.4010 | 0.4017 | −0.0008 |

Similarity holds to \|Δε\| ≤ 0.0026 across an 8.8 % clearance change under the full
energy equation — an order of magnitude below every effect discussed here. The
method's central assumption is therefore not a weak point.

**(b) Reference-curve sensitivity.** Every `c_eff` was recomputed against three
independent reference curves: the matching pancake THD run, the pancake isothermal
plain-bearing locus interpolated to `L/D` = 0.8, and the Raimondi & Boyd table
interpolated the same way. Results in §5 and §6. The *conclusions* — flat for
Ferron, sloped for Brito — are identical under all three.

---

## 4. Case 02 versus the transcribed Raimondi & Boyd table

`docs/UNFINISHED_BUSINESS.md` §2.1 lists "transcribe R&B Table 1 … and overlay on
the existing `locus_LD1.csv`" as the pending first step, and case 02's own SUMMARY
scores itself against "remembered textbook values (NOT a verified digitisation)".
**A transcription is already in the repository** and the comparison can be done now.

**Provenance.** `validation/journal_bearing/raimondi_boyd.py`, `RB_TABLE` /
`RB_TABLE_FULL`, transcribed by Seongho Han in commit `fb204b1` (2026-07-31,
"Validate the cavitation model and energy equation against published data"). The
source is *not* Part III directly: it is "Table (1): Dimensionless performance
parameters for full journal bearing" on p.192 of the University of Babylon lecture
notes *Sliding Bearings*
(`https://www.uobabylon.edu.iq/eprints/publication_7_4218_752.pdf`, retrieved
2026-07-31), the standard textbook reprint of Raimondi, A. A. and Boyd, J. (1958),
ASLE Transactions **1**(1), 159–209. The module records that textbook reprints of
this table are known to disagree in the third digit and audits it two ways: against
a closed-form Swift–Stieber infinite-length solution computed in the same module
(the `L/D` = ∞ block), and against an independent worked example in Al-Hamood et
al. (2015). The `S` convention there — `S = (mu n_s / P)(r/c)^2` with `P = W/(LD)`
— is identical to case 02's.

**Score.** `attitude_angle_deg = -90` sweep, flooded, FBNS, `n_theta` = 720.

| `L/D` = 1 | ε | `S` pancake | `S` R&B | Δ`S` | φ pancake | φ R&B | Δφ |
|---|---|---|---|---|---|---|---|
| | 0.1 | 1.31595 | 1.3300 | −1.06 % | 79.43° | 79.50° | −0.07° |
| | 0.2 | 0.62779 | 0.6310 | −0.51 % | 73.86° | 74.02° | −0.16° |
| | 0.4 | 0.25936 | 0.2640 | −1.76 % | 62.56° | 63.10° | −0.54° |
| | 0.6 | 0.12060 | 0.1210 | −0.33 % | 50.50° | 50.58° | −0.08° |
| | 0.8 | 0.04448 | 0.0446 | −0.27 % | 36.19° | 36.24° | −0.05° |
| | 0.9 | 0.01880 | 0.0188 | −0.01 % | 26.43° | 26.45° | −0.02° |

| `L/D` = 0.5 | ε | `S` pancake | `S` R&B | Δ`S` | φ pancake | φ R&B | Δφ |
|---|---|---|---|---|---|---|---|
| | 0.1 | 4.20149 | 4.3100 | −2.52 % | 81.93° | 81.62° | +0.31° |
| | 0.2 | 2.00096 | 2.0300 | −1.43 % | 75.28° | 74.94° | +0.34° |
| | 0.4 | 0.77983 | 0.7790 | +0.11 % | 61.87° | 61.45° | +0.42° |
| | 0.6 | 0.31809 | 0.3190 | −0.29 % | 48.08° | 48.14° | −0.06° |
| | 0.8 | 0.09155 | 0.0923 | −0.81 % | 32.91° | 33.31° | −0.40° |
| | 0.9 | 0.03121 | 0.0313 | −0.30 % | 23.53° | 23.66° | −0.13° |

Max \|Δ`S`\| = 1.76 % (`L/D` = 1) and 2.52 % (`L/D` = 0.5); max \|Δφ\| = 0.54° and
0.42°. Case 02's remembered-value sanity check ("within about 2 % in `S` and 0.5°
in φ") is **confirmed against an actual transcription**, and the worst point is the
ε = 0.1 / `L/D` = 0.5 corner where the reprint's own third digit is least reliable.

Two things this is worth beyond the score:

- Raimondi & Boyd impose the **Swift–Stieber** rupture condition with no mass
  conservation in the cavity; case 02 runs **mass-conserving FBNS/JFO**. Agreeing
  to ≤2.5 % in `S` over the whole ε range is also a statement about how little the
  cavitation model changes the *load locus* of a flooded plain bearing — which is
  exactly why the load discrepancy against the experiments could never have been a
  cavitation-model problem.
- **It does not fully close §2.1 option 1.** This is the textbook reprint of the
  design data, not Part III, and it carries no `L/D` = 1/4 block. A Part III
  transcription is still worth having for the third digit; it is no longer worth
  *blocking* on.

---

## 5. Ferron 1983 — per-point effective clearance

Reference: pancake THD, `L/D` = 0.8, 2000 rpm (`cases/01_ferron_thd` Part B).
`Lambda` is clearance-free and is the same in every variant.

| `W` [N] | ε (145 µm basis) | `Lambda` | ε* | `c_eff` [µm] | ±σ [µm] |
|---|---|---|---|---|---|
| 2018 | 0.4921 | 0.5564 | 0.4412 | 161.7 | 2.3 |
| 3006 | 0.5968 | 1.2189 | 0.5356 | 161.5 | 2.1 |
| 3966 | 0.6610 | 1.9732 | 0.5970 | 160.5 | 2.1 |
| 5025 | 0.7161 | 2.9341 | 0.6482 | 160.2 | 2.0 |
| 6058 | 0.7514 | 3.8942 | 0.6843 | 159.2 | 2.0 |
| 7017 | 0.7859 | 4.9337 | 0.7141 | 159.6 | 2.0 |
| 8016 | 0.8129 | 6.0306 | 0.7386 | 159.6 | 2.0 |
| 9047 | 0.8300 | 7.0961 | 0.7583 | 158.7 | 2.0 |
| 10063 | 0.8504 | 8.2854 | 0.7770 | 158.7 | 1.9 (extrapolated) |

**mean 160.0 µm, sd 1.1 µm, peak-to-peak 1.9 % over a 5× load range.**

Reference-curve sensitivity:

| reference | mean `c_eff` | sd | peak-to-peak |
|---|---|---|---|
| pancake THD, `L/D` = 0.8 | 160.0 µm | 1.1 µm | 1.9 % |
| pancake isothermal, `L/D` = 0.8 | 162.9 µm | 1.6 µm | 2.8 % |
| Raimondi & Boyd, `L/D` = 0.8 | 163.0 µm | 1.0 µm | 2.0 % |

Quoted result: **`c_eff` (Ferron) = 160 ± 4 µm**, where ±4 combines the ±1–2 µm
point scatter, the ±2 µm per-point digitising systematic and the ±3 µm spread
across reference curves.

Residual in ε after imposing the single value 160.0 µm:

| `W` [N] | `S` | ε (exp, rescaled) | ε (solver) | Δε |
|---|---|---|---|---|
| 2018 | 0.3575 | 0.4460 | 0.4358 | +0.0102 |
| 3006 | 0.2400 | 0.5409 | 0.5309 | +0.0100 |
| 3966 | 0.1819 | 0.5991 | 0.5955 | +0.0037 |
| 5025 | 0.1436 | 0.6491 | 0.6476 | +0.0014 |
| 6058 | 0.1191 | 0.6811 | 0.6861 | −0.0051 |
| 7017 | 0.1028 | 0.7123 | 0.7151 | −0.0028 |
| 8016 | 0.0900 | 0.7368 | 0.7395 | −0.0027 |
| 9047 | 0.0798 | 0.7523 | 0.7611 | −0.0088 |
| 10063 | 0.0717 | 0.7708 | 0.7798 | −0.0090 |

rms 0.0066 against a digitising σ_ε of 0.015. **The collapse is complete: the
residual is below the reading noise.**

*The one structure that is there, stated honestly.* `c_eff` drifts −2.0 µm per
e-fold in load, i.e. −3.2 µm (−2.0 %) across the whole sweep, and the sign is
resolvable against the ±0.33 µm point-to-point scatter about the fit (t ≈ 9). It
is nonetheless smaller than the ±2 µm digitising systematic, which is correlated
across points and so cannot be beaten down by having nine of them. Its direction
— experiment above the solver at low load, below at high load — matches the tilt
case 01 already reported ("the measured points are softer than the solver at low
load"), and is consistent with the solver running too cold and therefore too
viscous where the film is thickest. It is **1/5 the size of the Brito speed
trend** and should not be over-read.

---

## 6. Ferron: which published clearance is right

The value that works, 160 µm, is *inside* the interval already spanned by the
source's own numbers. Ferron's Fig. 5 legend carries both 145 µm and 152 µm; the
Lubricants 11(7) 308 (2023) re-reduction argues for ≈171 µm. No new hypothesis is
needed. The three candidates score as follows, in terms of the solver load carried
at the measured probe displacement divided by the applied load:

| `W` applied [N] | `c` = 145 µm | 152 µm | **160 µm** | 171 µm |
|---|---|---|---|---|
| 2018 | 1.54 | 1.27 | **1.04** | 0.82 |
| 3006 | 1.62 | 1.30 | **1.04** | 0.79 |
| 3966 | 1.64 | 1.30 | **1.02** | 0.75 |
| 5025 | 1.71 | 1.32 | **1.01** | 0.73 |
| 6058 | 1.72 | 1.29 | **0.97** | 0.69 |
| 7017 | 1.81 | 1.34 | **0.98** | 0.68 |
| 8016 | 1.85 | 1.36 | **0.98** | 0.67 |
| 9047 | 1.81 | 1.32 | **0.95** | 0.64 |
| 10063 | 1.84 | 1.33 | **0.95** | 0.63 |

- 145 µm: +54 % to +85 %. The 2–4 kN end reproduces the campaign's reported
  +52–68 %.
- 171 µm applied as a proper multiplicative rescale: **−18 % to −37 %, an
  over-correction**. The "+0.011 ± 0.003" that case 01 obtained for this
  re-reduction used the additive "ε − 0.10", which is not a clearance rescale.
  This analysis therefore *sharpens* case 01's "both readings cannot be right and
  this case cannot resolve which" into a point estimate that rejects both endpoints.
- 160 µm: within ±5 % at every load.

**Sensitivity, and why the load FAIL looked so large.** At fixed probe
displacement the tables above give `W ∝ c^n` with `n` between −3.8 and −4.8
(−4.8 for Ferron 145 → 160 µm, −4.3 to −4.7 for Brito, −3.8 for the long
extrapolation to 171 µm). A 1 % error in the assumed clearance moves the predicted
load by roughly 4–5 %. The entire +52–85 % Ferron load discrepancy is what a
**10 % clearance error** looks like. Any load score taken at a probe-derived
eccentricity is, to first order, a measurement of the clearance basis.

---

## 7. Brito 2006 — per-point effective clearance

Reference: pancake THD on the Brito rig, `Pf` = 70 kPa gauge, ε carrying the
applied 10 kN. Measured ε is the mean of the two independent digitisations.

| `N` [rpm] | ε measured | `Lambda` | ε* | `c_eff` [µm] | ±σ [µm] |
|---|---|---|---|---|---|
| 1000 | 0.7657 | 4.3884 | 0.6865 | 95.4 | 1.7 (extrapolated) |
| 3000 | 0.5764 | 0.8291 | 0.4954 | 99.5 | 1.9 |
| 4000 | 0.5514 | 0.5688 | 0.4522 | 104.3 | 1.9 |

**mean 99.7 µm, sd 4.5 µm, +9.3 % from 1000 to 4000 rpm.** Quoted clearance is
85.5 µm radial (171 µm diametral at 20 °C).

Reference-curve sensitivity — the trend survives, and is *smallest* against the
THD reference:

| reference | mean `c_eff` | sd | 1000 → 4000 rpm |
|---|---|---|---|
| pancake THD, Brito rig | 99.7 µm | 4.5 µm | **+9.3 %** |
| pancake isothermal, `L/D` = 0.8 | 102.4 µm | 7.0 µm | +14.7 % |
| Raimondi & Boyd, `L/D` = 0.8 | 102.1 µm | 6.4 µm | +13.4 % |

Independent replication at the second supply pressure (`Pf` = 210 kPa; scored
against the same 70 kPa solver curve, so read the trend, not the level):
94.4 / 98.4 / 101.9 µm, **+7.9 %**.

Residual in ε after imposing the single value 99.7 µm:

| `N` [rpm] | `S` | ε (exp, rescaled) | ε (solver) | Δε |
|---|---|---|---|---|
| 1000 | 0.0983 | 0.6566 | 0.7018 | −0.0452 |
| 3000 | 0.2948 | 0.4943 | 0.4963 | −0.0019 |
| 4000 | 0.3930 | 0.4728 | 0.4309 | +0.0420 |

rms 0.037 against a digitising σ_ε of 0.02: **about twice the reading noise, and
monotone in speed.** Compare Ferron's 0.0066 against 0.015.

Load over-prediction at the measured displacement, for reference:

| `N` [rpm] | `c` = 85.5 µm | 93.0 µm | 99.7 µm | 104.3 µm |
|---|---|---|---|---|
| 1000 | 1.99 | 1.16 | 0.77 | 0.60 |
| 3000 | 2.03 | 1.35 | 0.99 | 0.82 |
| 4000 | 2.37 | 1.61 | 1.20 | 1.00 |

**Two caveats on this table specifically.** (i) The Brito solver reference is only
three points, and the 1000 rpm column needs extrapolation past its top; the entries
are indicative. (ii) These ratios (+99 % to +137 % at the quoted clearance) are
larger than the "+65–106 %" in `validation/VALIDATION_REPORT.md` §4.1, because that
band mixes the 20 °C and hot-clearance rows. The hot-clearance bracket there is
also **basis-inconsistent**: it changes the solver's `c` while leaving the measured
ε on the 85.5 µm basis, so it changes `S` without changing ε. The columns above
change both together, which is the only self-consistent way to vary a clearance
basis. None of §0's conclusions depend on this table.

**Weaknesses of the Brito evidence, stated plainly.**

- The 1000 rpm solver row **did not converge** (coupling capped at 80 iterations,
  `R_p` = 5.23) and also requires extrapolating a three-point curve. It is the
  weakest of the three.
- Dropping it entirely, the 3000 → 4000 rpm step alone is +4.8 µm against a
  combined ±2.7 µm — **1.8σ**, positive and monotone but not on its own decisive.
  Over the full range it is 3.5σ.
- Brito's bearing has two axial grooves; its solver reference is that geometry, but
  the two clearance-free cross-checks (isothermal plain, R&B) are plain-bearing
  curves. They are cross-checks on the trend, not substitutes for the reference.
- The tolerance Brito quotes on the 171 µm diametral clearance was not available
  to this analysis, so "is 199 µm diametral within tolerance?" cannot be answered
  here.

---

## 8. What this cannot decide, and why

**The level of `c_eff` is degenerate with a uniform solver load error.** Scaling
the solver's load by a constant factor `f` at fixed ε multiplies `Lambda_solver`
by `f`; so does mis-stating the clearance. The two are indistinguishable in the
`S`–ε plane, by construction, for any dataset whose ε is probe-derived. Therefore:

- The fact that **both** rigs require a clearance larger than quoted — Ferron
  +10 % over 145 µm, Brito +17 % over 85.5 µm at the slowest speed — is equally
  consistent with (a) a shared clearance-calibration systematic on the LMS-Poitiers
  rig lineage that both papers inherit, and (b) pancake over-predicting load by a
  roughly constant factor of 1.6–2.0. Nothing here separates them.
- What *is* separated is structure. Constant-in-load (Ferron) rules out a
  load-dependent physics error. Rising-in-speed (Brito) rules out "one clearance
  number explains everything".

**The surviving speed trend has two live explanations and this data cannot choose.**

1. **Real speed-dependent clearance.** Brito quotes `c` at 20 °C; the film runs
   40–55 °C and hotter at 4000 rpm. Thermal (and centrifugal) growth of the bush
   is physically required to be positive and speed-increasing. The required +9 %
   over 1000 → 4000 rpm is the same order as the +8.8 % uniform hot bracket already
   in the record. pancake has **no thermo-elastic clearance model at all**
   (`validation/VALIDATION_REPORT.md` §5.3), so this is a modelling gap in the
   *geometry*. The data do not prove that interpretation over a uniform load bias.
2. **Under-predicted film temperature.** Too little viscosity drop leaves too much
   load capacity, hence too low an ε — the observed sign, and the same direction as
   the Ferron −2.3 K peak residual and the 5–7 K cold bias in
   `cases/01_ferron_thd/SUMMARY.md`.

The evidence tilts toward the mechanism being thermal in general — the THD
reference needs a 9.3 % trend where the isothermal reference needs 14.7 %, so the
energy equation is already absorbing roughly a third of it — but it does not
distinguish "the metal moved" from "the oil was hotter than we computed".

---

## 9. The single measurement that would settle it

**An absolute, in-situ film-thickness measurement at one operating point** —
ultrasonic through the bush, in the manner of Kasolang & Dwyer-Joyce,
*Tribology Transactions* **51**(2):231–242 (2008), already noted in
`docs/UNFINISHED_BUSINESS.md` §3.6 with a free accepted version at White Rose.

Why this one and not another:

- It returns `h(θ)` in micrometres **without dividing by an assumed clearance**.
  From `h_max` and `h_min` at a single condition you get `c = (h_max + h_min)/2`
  and `e = (h_max − h_min)/2` **separately**. That is precisely the degeneracy of
  §8: it fixes the level of `c_eff` independently, and so decides whether 160 µm
  and ~100 µm are real clearances (instrumentation) or whether the solver's load
  level is high by the corresponding factor.
- Repeated at two speeds on the same rig it *also* settles §8's second question,
  by measuring `c(N)` directly and saying whether the required +9 % is metal
  motion or a thermal-model error.
- Nothing cheaper does either job. Re-digitising the figures cannot (commit
  `81dd953` already proved the digitisation is not the problem). More solver runs
  cannot (they are one side of the comparison). A Raimondi & Boyd or Someya
  transcription cannot either — §4 shows pancake already matches R&B to ≤2.5 %,
  so a table comparison validates the solver but says nothing about either rig's
  clearance.

A cheaper partial substitute, if only an existing rig is available: perform the
centring/"clearance circle" calibration at the **operating** temperature and speed
rather than cold. That yields `c(N, T)` and settles the trend, but not the level.

---

## 10. Consequences for the record

Recommended re-scoring of `validation/VALIDATION_REPORT.md` §3 and
`docs/UNFINISHED_BUSINESS.md` §2.1 (proposed here, not applied):

| row | current | proposed |
|---|---|---|
| Ferron, load at measured ε | **OPEN**, +52–68 % | **NOT SCOREABLE as a load test** — the full discrepancy is absorbed by a single `c_eff` = 160 ± 4 µm, a value bracketed by the source's own two legend clearances and the 2023 re-reduction; residual after collapse is below the digitising noise at every load from 2 to 10 kN |
| Brito, load via ε-for-load | **FAIL / OPEN**, Δε = −0.113/−0.150/−0.191 | **OPEN, re-framed** — not a load-capacity failure; a *speed-dependent* residual equivalent to `c_eff` rising +9.3 % over 1000–4000 rpm, of which the energy equation already removes ≈1/3. Re-run the 1000 rpm point (non-converged) before quoting |
| Raimondi & Boyd, `S` and φ | **PASS** (from `validation/journal_bearing`) | extend to case 02: ≤1.8 % in `S` and ≤0.54° in φ at `L/D` = 1, ≤2.5 % and ≤0.42° at `L/D` = 0.5, against the transcription in `raimondi_boyd.py` (§4) |

Three follow-ups this analysis makes cheap:

1. **Re-run the Brito 1000 rpm point to convergence** with the ASM/fallback ladder
   from `f0e772d`. It is the weakest leg of the only surviving residual.
2. **Report load discrepancies as `c_eff`, not as a percentage.** With
   `W ∝ c^(−4 … −5)`, a percentage in load is a misleadingly dramatic restatement
   of a small clearance uncertainty.
3. **Any future load score must use a source whose independent variable is not a
   probe-derived eccentricity** — either R&B/Someya (dimensionless), or a rig with
   an absolute film-thickness measurement. `docs/UNFINISHED_BUSINESS.md` §2.1's
   ordering of preferences was right; §4 above partially discharges its option 1.

---

## 11. Figure

`s_eps_collapse.png` (single 300 dpi PNG, lab style per `plot_style.py`), three
panels, ε against a logarithmic `S`, datasets as markers and solver families as
lines:

- **(a) as published.** Ferron on its 145 µm basis and on the 171 µm re-reduction,
  Brito on its 85.5 µm basis, against pancake THD (Ferron rig, Brito rig), pancake
  isothermal `L/D` = 0.8 and Raimondi & Boyd `L/D` = 0.8. The black guide is the
  clearance-rescale trajectory `S ∝ ε²` — the direction, and the *only* direction,
  in which a clearance change can move a point.
- **(b) one `c_eff` per rig.** Ferron ×160 µm collapses onto the solver curve over
  the whole load range. Brito ×100 µm cannot: the 1000 rpm point sits above the
  curve and the 4000 rpm point below it.
- **(c) per-point `c_eff`.** Ferron against load (lower blue axis) is flat at
  160 µm; Brito against speed (upper orange axis) climbs from 95 to 104 µm. Grey
  lines mark the published clearances 145, 152, 171 µm (Ferron) and 85.5 µm
  (Brito). **This panel is the discriminator.**

## 12. Files

- `collapse_analysis.py` — the whole analysis; prints every table above.
  `python collapse_analysis.py` from this directory.
- `collapse_analysis.json` — every number, machine-readable.
- `s_eps_collapse.png` — the figure.
