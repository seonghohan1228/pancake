# Unfinished business — validation campaign inventory

> **Snapshot, not current status.** Written before case 05a, the Picard-cap fix,
> and the Sommerfeld-collapse analysis. Keep it as the detailed backlog and audit
> trail; use `cases/VALIDATION_STATUS.md` for current verdicts and priorities.

Date: 2026-08-05. Written against `codex/thermal-ksp-robustness` @ `f0e772d`.

What this is: a prioritised, actionable list of what the validation campaign has
**not** closed, why each item matters, what it needs, and the first concrete step.
Companion to `cases/VALIDATION_PLAN.md` (the case ladder), `HANDOFF.md` (verified
defects D1–D16), `docs/CODE_AUDIT_2026-08.md` (code-level legacy/maintainability
findings — §1.3 and §2.1 there are cited below as live physics risks) and
`validation/VALIDATION_REPORT.md` (the previous campaign's verdict table, now
partly stale — see §4.9).

Where the campaign stands: grid/dt independence **PASSED** and frozen (FBNS, G6
720×192, dt = 2e-5 s); FBNS verified against the exact 1-D JFO solution
(`cases/00c_jfo_verification/`); Ferron pressure **PASS** (+0.8 %, 1.25°),
eccentricity consistent within the published clearance ambiguity, THD temperature
with 2-D conjugate bush at mean +0.02 K / rms 1.17 K with a **−2.3 K peak
residual** and the user verdict pending.

---

## 1. Top five, with reasoning

| # | Item | Why it is first | §|
|---|---|---|---|
| 1 | **Load capacity is an open FAIL, and the postponement of R&B/Someya left it un-adjudicated** | The campaign's own headline quantity. Ferron: **OPEN, +52–68 %** at measured ε. Brito: **FAIL/OPEN**, solver reaches the applied 10 kN at ε 0.12–0.18 below measurement, growing with speed. The digitisation excuse is gone (Brito Fig. 6 independently re-read, `81dd953`). Nothing in the current ladder scores load against a source whose clearance is not in dispute. | §2.1 |
| 2 | **Gaseous cavitation — the campaign's stated goal — has no experimental comparison, and the model cannot yet make the comparison the best dataset asks for** | Two gaps, not one. (a) No case 05 exists. (b) Etsion & Ludwig measured ~50 kPa of *variation inside the cavity*, and **no converged pancake path predicts a spatially varying cavity pressure**: FBNS pins `p = p_cav` (scalar), `VOID_COUPLED` has the per-cell field but does not converge. Split the case: extent is scoreable today, pressure is blocked. | §3 |
| 3 | **Etsion & Ludwig 1982 is on disk, complete, and closes items 1–3 of `cases/REFERENCE_DATA.md` §3** | Every parameter pancake needs is printed in the paper (R, L, c, ε, speeds, µ, ρ, four ambient levels), the cavity contents were mass-spec'd as **air with no oil vapour**, and the oil vapour pressure is 13 Pa — an independent proof the mechanism is gaseous. This is the single highest-value unopened box in the whole campaign. | §3.1 |
| 4 | **Case 03 (transient → equilibrium) is ~1 day of machine time and closes a ladder row** | The machinery already ran in case 00b (`MOVING_BEARING`, 6 kN load step, G6, FBNS). The only new thing is comparing the transient endpoint to case 01's static force-balanced ε. It is the cheapest remaining PASS/FAIL on the board and it independently checks the force-balance driver. | §2.2 |
| 5 | **The −2.3 K Ferron peak residual has four named suspects and none has been run** | The temperature verdict is *pending*, not *closed*, and a pending verdict blocks the THD claim. Four levers exist (groove λ, collocation saturation, global-adiabatic shaft, STRIATED). Each is a config change on a frozen protocol — hours, not days. | §2.4 |

Honest counterweight: item 1 is the one that could still be a *measurement*
problem rather than a solver problem — three independent lines of evidence in
`validation/VALIDATION_REPORT.md` §4 say so. That is exactly why it needs a
source whose clearance is not derived from a proximity probe.

**One finding from §5 that is not a speed issue and should be read before the
rest.** The conjugate Ferron runs (`outer_robin80`, `rerun_cht_robin10`) hit the
200-iteration Picard cap on every run, and the coupling residual **rises after
outer ≈ 51 and then oscillates**. The loop returns iterate 200, which is worse
than iterate 51. The archived conjugate results are therefore **not the best
iterate the solver found** — including the −2.3 K peak in item 5. Check whether
that moves any reported number before investing in §2.4's four suspects. Both
relaxation levers already exist in the code and are switched off (set to 1.0).

Sections §5 (acceleration), §6 (why finite volume) and §7 (Elrod vs FBNS vs LCP)
are reference material for the write-up rather than campaign tasks.

---

## 2. Open validation targets

### 2.1 Load-capacity locus — R&B Part III / Someya (POSTPONED), and the open load FAIL

**What.** Case 02 produced both solver loci (L/D = 1 and 0.5, ε = 0.1…0.9,
`cases/02_static_locus/SUMMARY.md`) and they are waiting on reference tables that
the user postponed on 08-05 (R&B Part I in the library is the 120° partial
bearing; Someya is a book, not in the PDF set). Load validation was folded into
case 01 via Ferron Fig. 5.

**Why it matters.** Folding it into Ferron does not resolve it, because Ferron's
ε is a proximity-probe displacement over an *assumed* clearance — the very
quantity in dispute. The record currently reads:

| Source | Verdict | Size |
|---|---|---|
| Ferron 1983, load at measured ε | **OPEN** | +52–68 % |
| Brito 2006, ε carrying the applied 10 kN | **FAIL/OPEN** | Δε = −0.113 / −0.150 / −0.191 at 1000/3000/4000 rpm |
| Case 01 THD Fig. 5 sweep | consistent | flat −0.055 ± 0.005 under the 145/152 rescale; +0.011 ± 0.003 under the Lubricants 2023 re-reduction |

Two independent readings of Ferron Fig. 5 differ by 0.10 in ε — an order of
magnitude more than anything the solver contributes. Brito's applied load is
known exactly and the deficit still grows with speed. Case 02's own sanity check
(unverified from memory) put every point within ~2 % in S and 0.5° in φ of
remembered R&B values.

**What it needs.** A load reference whose independent variable is not a
clearance-derived eccentricity. In order of preference:

1. **Raimondi & Boyd 1958 Part III tables** (full journal bearing) or the **Someya
   databook** — dimensionless, no clearance ambiguity, and case 02's solver curves
   are already computed and archived. Pure digitisation/transcription work.
2. **Ferron's own theory curve** rather than his measurements — the framing
   `validation/VALIDATION_REPORT.md` §4 already recommends. His theory
   over-predicts his own data by ~75–80 %, worse than pancake; reproducing his
   theory is a code-to-code check that isolates the solver from the clearance
   question.
3. **A dimensionless collapse.** Plot S vs ε for Ferron, Brito and pancake on one
   axis. If the three collapse under a common clearance rescale, the discrepancy
   is instrumentation; if they do not, it is physics.

**Effort.** (1) is a few hours of transcription plus a plotting pass — the solver
side is done. (2) is 1–2 days (implement Ferron's own THD closure as a variant,
or accept the published curve digitised). (3) is half a day of plotting.

**First step.** Transcribe R&B Part III Table 1 (L/D = 1, and L/D = 1/2) for
S, ε, φ into `cases/ref_data/raimondi_boyd_1958_full_LD1.csv` per
`cases/REFERENCE_DATA.md`, noting the source's cavitation BC in the header, and
overlay on the existing `cases/02_static_locus/locus_LD1.csv`. If the tables
cannot be sourced, run option (3) first — it costs half a day and it tells you
whether to keep looking.

### 2.2 Case 03 — transient to equilibrium (trivial now)

**What.** Release the bearing under a constant external load from a cold start
and check that the trajectory settles on the same equilibrium the case 01 static
force-balance found.

**Why it matters.** It is a cross-check of two independent code paths onto one
number: the `MOVING_BEARING` ODE + transient FBNS versus the Python damped-Newton
force balance in `cases/01_ferron_thd/force_balance.py`. If they disagree, one of
them is wrong and every case-01 equilibrium is suspect. It is also the only
remaining ladder row that needs no external data (`cases/REFERENCE_DATA.md`:
"Cases 00/00b (internal) and 03 (endpoint vs 01/02) need nothing from you").

**What it needs.** Nothing new. `cases/00_grid_dt_independence/dt_T3_fbns/config.txt`
is already a G6 720×192 FBNS transient with `motion_model = MOVING_BEARING`,
`bearing_mass = 10`, `external_load_y = 6000`, `end_t = 0.06`. Swap dt to the
frozen 2e-5 s, put in the Ferron parameter set (R = 50 mm, L = 80 mm,
c = 152 µm, groove at +90°, µ(40 °C) = 0.0277, thermoviscosity), and run at
2 kN / 4 kN / 6 kN.

**Effort.** Half a day of setup, a few hours of machine time on afml-r6.

**First step.** `cp cases/00_grid_dt_independence/dt_T3_fbns/config.txt
cases/03_transient_equilibrium/config_template.txt`, retarget the geometry, and
score `|ε_transient(t→∞) − ε_static|` against the force-balance tolerance
(0.1 % in |W|). Caution already recorded in case 00b: **detrend before any DFT**
of the trajectory — a naive transform shows a fake 33 Hz "half-whirl" that is a
window bin-1 artefact.

### 2.3 2-D contour verification vs Ausas 2007 / Woloszynski 2015 (POSTPONED)

**What.** Case 00c test #2: compare the 2-D pressure field and the cavitation
boundary against a published grooved-journal-bearing benchmark.

**Why it matters.** Everything verified so far about the free boundary is
**1-D** (`cases/00c_jfo_verification` test #1, the exact infinitely-long
solution). The FBNS front is ≤ 0.92 cells there — but nothing has tested the
*axial* shape of the rupture boundary, which is where a 2-D grooved bearing
differs from a 1-D one and where the axial Dirichlet half-cell link (§7.3.4) is
exercised. Case 05 will need exactly this capability (Etsion & Ludwig's data are
axial scans, §3.1).

**What it needs.** Digitised contours. Woloszynski, Podsiadlo & Stachowiak 2015
(Tribol. Lett. 58:18) is the preferred source — the user has said it is open
access and they would digitise it themselves. Ausas et al. 2007 (ASME J. Tribol.
129:868) is the alternative. Note Ausas, Jai & Buscaglia 2009 (ASME JT
131(3):031702) **ships reference source code** (`PLAN.md`:1831), so it is a
code-to-code option needing no digitisation at all.

**Effort.** 1 day once contours exist; 2–3 days for the Ausas 2009 code-to-code
route including building their code.

**First step.** Ask the user to confirm Woloszynski 2015 is reachable and to
digitise one cavitation-boundary trace and one centreline profile. In parallel,
add an axial-profile extractor to `cases/*/solver_io.py` — it is needed for case
05 regardless (§3.1 caveat).

### 2.4 The Ferron −2.3 K peak residual — four untried suspects

**Current best.** `cases/01_ferron_wallbc_fields/SUMMARY.md`, group
`rerun_cht_robin10` (2-D conjugate bush, Robin h = 10 to 40 °C, shaft globally
adiabatic, force-balanced at ε = 0.6033): peak 46.71 °C @ 265.2° against a
measured 49.04 °C @ 258.2°; mean **+0.02 K**, rms **1.17 K** over 16 points; peak
value **−2.3 K**, early arc ~+1 K hot.

So the aggregate is excellent and the *shape* has a systematic tilt: too cold at
the peak, too hot early. That is the signature of insufficient circumferential
transport, not of a level error.

Suspects, in the order they should be tried:

| Suspect | Mechanism | Cost | Expected sign |
|---|---|---|---|
| **Groove mixing λ** | `FLOW_WEIGHTED` is running with `Q_supply` derived from the solver's own mass balance (case 01 Deviation 1), giving λ ≈ 0.8 — i.e. 80 % fresh oil. If the real rig recirculates more, the whole arc starts hotter and the peak rises. `PINNED` (λ = 1) is the other bracket and moves ε only +0.5 %. | hours | raises the whole curve, incl. the early arc that is already hot — likely **wrong sign** |
| **Collocation saturation** | `PHYSICS.md` §13.1.2: the convective term is saturated against the dissipative one because unbounded it produced several-hundred-K cores. Reduced Péclet on this rig is 10–30, so **the saturation binds**, and where it binds "the quartic degenerates towards a fixed flattened shape" (`PLAN.md`:142–146). A flattened profile under-reports the bush-side surface temperature. | 1–2 days | **raises the peak**; the leading candidate |
| **Shaft global-adiabatic** | `shaft_ambient_heat_transfer = 0` makes the shaft carry heat from the loaded arc back upstream, which is *why* the peak moved correctly. But a real shaft loses heat to ambient; a non-zero value would cool the early arc (which is +1 K hot) more than the peak. | hours | **fixes the tilt**, small effect on peak |
| **STRIATED cavitated thermal** | `cavitated_film_model = STRIATED` vs the `FULL_FILM` default — **never run on any Ferron case**. D14 measured *no steady effect* on this rig (T_max moves 0.003 K; θ cancels exactly between source and mass flux) so this is close to a formality, but it is cheap and it removes a named suspect. | hours | ~zero, per the D14 measurement |

**Also unchecked:** the Ferron polynomial viscosity law is temperature-only — no
Barus piezo-viscous term at the rig's ~70 bar peak (`PLAN.md`:190–196 flags that
no validation case has been re-run on `POLYNOMIAL` with that gap stated).

**Effort.** A/B/C/D as a single batch: ~2 days including the force-balance
re-closure at each variant.

**First step.** Run the STRIATED and shaft-Robin variants first (both are
one-token config changes on the frozen `01_ferron_wallbc_fields` protocol) to
clear two suspects cheaply, then invest in the collocation study. **Re-close the
force balance for every variant** — the 01-conjugate CORRECTION is explicit that
a frozen-geometry probe is not comparable to a re-closed equilibrium (it is what
made an earlier h = 10 conclusion wrong by 1.3 K).

### 2.5 Case 04 — THD on a second dataset (optional, on request)

Boncompain et al. 1986 / Dowson 1966–67. Not started. Worth doing only after
§2.4 closes; a second dataset is how you tell a Ferron-specific artefact from a
model error, and the Brito rig already gives a partial second opinion (its 5 of 6
committed runs are unconverged and need re-running — `validation/VALIDATION_REPORT.md`
line 141, "**INCOMPLETE**").

---

## 3. Gaseous cavitation — the stated goal, currently unvalidated

The campaign's goal statement includes gaseous cavitation. There is **no case 05
directory**, no reference data in `cases/ref_data/`, and no experimental
comparison of any kind. This section ranks the concrete options.

### 3.0 First, what pancake can and cannot model today

This determines which comparison figures are even askable.

| Capability | Status | Where |
|---|---|---|
| Cavity **extent** under gaseous cavitation | **Runnable and claimed trustworthy** | `PLAN.md`:337–340 — "Predicted cavity *extent* is trustworthy; predicted cavity *pressure* is nominal" |
| Uniform onset threshold at the bubble point | `cavitation_threshold = LOCAL_PSAT` | uniform **by construction** — `effective_p_cav()` evaluates the bubble point at the *feed* state (config scalars), so it is a datum, not a field (`PHYSICS.md` §7.3.8) |
| Per-cell bubble-point field `p_sat(T, c_d)` | **NOT IMPLEMENTED** (HANDOFF D4) | this is the Tier-0 → Tier-1 step that `docs/CAVITATION_OIL_REFRIGERANT.md` §9.2 calls "the #1 recommended change" |
| Per-cell cavity plateau `p_void` | exists on the Elrod path only, **does not converge** | `R_p` 0.57–7.2 against a 1e-9 tolerance, relaxation- and grid-dependent (`PLAN.md`:261–274) |
| `VOID_COUPLED` under FBNS | **rejected at config validation** | `src/config.hpp:1474`; `docs/FBNS_GAS_COUPLING.md` §8.3 |
| FBNS + `TWO_PHASE` (moved ceiling θ ≤ 1 − α_gas) | **landed and verified** — complementarity invariants exactly 0, per-cell mass residual 5.1e-15 of peak Couette flux | `PHYSICS.md` §7.3.7 |
| Two-phase µ closures (McAdams / Dukler / Einstein / K–D / Grando) | present, **uncalibrated for R290** | `fluid_properties.cpp`; open question §7 Q3 of the cavitation doc |
| PTSV solubility/viscosity table ingestion | **landed**, with real measured data | `data/R290_PZ68S/r290_pz68s.csv`, `config_r290_pz68_table.txt` |

**The consequence that matters.** On every *converged* path, the cavitated cell
sits at a **single scalar pressure**. That is precisely the assumption the best
available experiment refutes (§3.1). So Case 05 must be split:

- **05a — cavity extent.** Runnable today with FBNS or JFO. Scoreable. Do this.
- **05b — cavity pressure.** Blocked on the `p_void` convergence work (§4.3).
  Frame it honestly as a falsification target for a model the code already
  documents as nominal, not as a claim.

### 3.1 RANK 1 — Etsion & Ludwig 1982 (data is on disk, complete)

`F:\Seongho\Research\literature\pdfs\Etsion and Ludwig - 1982 - Observation of
Pressure Variation in the Cavitation Region of Submerged Journal Bearings.pdf`
— J. Lub. Tech. **104**(2):157–163. (The title in the old plan notes,
"dynamic behavior of flooded journal bearings", is wrong; volume/page match.)

**Rig, complete:** R = 25.4 mm, L = 38.1 mm, **L/D = 0.75**; mean radial clearance
**114 µm**; **ε = 0.4 for every reported condition**; acrylic housing, a single
0.2 mm tap on an indexing fixture scanning the whole cavity; differential
transducer referenced to atmosphere; **no groove — both ends flooded** at the
same supply pressure; runout < 12.7 µm; tests at **1840 and 3000 rpm**.

**Oil (the datum that was missing):** paraffinic, unnamed. Table 1 gives
µ = **0.026 Pa·s at 38 °C** (and 2.1 at −29 °C, 0.84 at −18 °C, 0.0047 at 99 °C,
0.0017 at 204 °C — enough to fit a µ(T) law), specific gravity **0.8285**
(ρ ≈ 828.5 kg/m³), and **vapour pressure 0.1 mmHg ≈ 13 Pa at 38 °C**. Oil-in
temperature 30–38 °C, annotated per panel.

**Why it is the right target.** They pulled cavity contents into an evacuated tube
and ran a **mass spectrometer: air, no oil vapour**. The oil vapour pressure is
13 Pa, three orders below the measured cavity pressures. That is an independent,
in-paper proof that the mechanism is gaseous — the exact regime
`docs/CAVITATION_OIL_REFRIGERANT.md` argues for, established by measurement
rather than by analogy. And **ambient pressure is the control variable**:
p_s = 13.6 / 27.2 / 40.8 / 54.4 kPa gauge, which is a direct dissolved-air lever.

**The figures.**
- **Table 2** — tabulated θ₀, p_min, θ_st, θ_end for all 8 conditions. *Transcribe,
  do not digitise.* **This is the case-05a comparison.**
- **Fig. 6 (1840 rpm) and Fig. 7 (3000 rpm)** — cavity-zone pressure maps, four
  panels each (one per p_s). x = dimensionless axial position Z/L ∈ [−0.5, 0.5];
  y = p/p_s (differential, gauge); 8–9 curves per panel labelled by circumferential
  θ. Digitisable at 170 dpi. **This is case 05b.**

**The headline result, and why it is a falsification target.** Cavity pressure is
**not constant**: ~50 kPa of variation inside the bubbles, and **p/p_s > 1**
where the cavity penetrates the converging clearance. pancake predicts a flat
plateau on every converged path (§3.0). Reporting that honestly — "extent
matches, the plateau assumption is refuted by the reference, here is the
magnitude of the error it costs" — is a stronger result than a fitted agreement.

**Comparison figures.**
- 05a: measured θ_st, θ_end, θ₀ vs solver, 8 conditions, one scatter with the
  1:1 line. Plus cavity extent vs p_s (their stated trend: raising p_s retreats
  the downstream cavity end toward max clearance).
- 05b: p/p_s vs Z/L at 3–4 of their labelled θ stations, one panel per p_s, solver
  overlaid.

**What pancake needs configured.** Submerged/flooded axial `DIRICHLET` ends at
p_s + ambient; no inlet groove — **but see the case 02 finding**: with no inlet,
mass-conserving FBNS settles into a *starved* state (cavitated fraction 0.76–0.86
at every ε, ~1/20 of the flooded load). Etsion & Ludwig's rig is genuinely
flooded through both end faces, so the axial Dirichlet ends *are* the feed and
this should be fine — but verify it, because the case 02 no-inlet runs are the
warning.

**Effort.** 05a: 1 day (transcribe Table 2, build 8 configs, run, score).
05b: 2–3 days of digitisation plus whatever §4.3 costs.

**Caveat flagged by the source hunt.** Their data are **axial** scans (p vs Z/L at
fixed θ). The existing harness (`cases/*/solver_io.py`) is built around midplane
θ-profiles. Add an axial extractor first — it is a small script change and §2.3
needs it too.

**First step.** Transcribe Table 2 into
`cases/ref_data/etsion_ludwig_1982_table2_cavity_extent.csv` with a header per
`cases/REFERENCE_DATA.md`, and stand up `cases/05_gaseous_cavitation/` with the
parameter table for user sign-off (all values are printed in the paper, so this
should sign off cleanly — no invented numbers).

### 3.2 RANK 2 — Ochiai et al. 2019 (open access, and the only geometry pancake can actually represent)

Ochiai, Sakai & Hashimoto, *Lubricants* **7**(9):74 (2019), DOI
10.3390/lubricants7090074. **CC-BY, fully open access.**

Photographic cavity extent in a **transparent plain journal bearing**:
D = 25 mm, radial clearance 125 µm, L = 14.5 mm (**L/D = 0.58**), 3500 rpm,
single feed hole, run **both flooded (2.6 cm³/s) and starved (0.5 cm³/s)**.

**Why it ranks this high despite giving less than Etsion & Ludwig.** It is a
plain cylindrical bearing with one feed hole — the exact geometry pancake meshes.
No pads, no pockets, no offset halves. The photographs give a θ–z cavity boundary
that maps straight onto `cases/REFERENCE_DATA.md`'s "contours/boundaries" format,
and the flooded/starved pair is a second, independent axis to score against
(pancake's case 02 already showed the solver has a genuine starved solution
branch, so this tests it deliberately rather than by accident).

**What it does not give:** no void fraction, no cavity pressure. Extent only.

**Effort.** 1–2 days including the boundary digitisation. Pairs naturally with
case 05a (§3.1) — same comparison figure, second data point.

**First step.** Download the PDF (open access) and confirm the photographs carry
a usable scale and a stated θ origin before committing to digitisation.

### 3.3 RANK 3 — Engels et al. 2024 (the only measured in-film void fraction; a physics claim, not a benchmark)

Engels, Wettmarshausen, Stottrop, Hagemann, Weißbacher, Schwarze & Bender,
*"Experimental Identification of the Void Fraction in a Large Hydrodynamic Offset
Halves Bearing"*, *Lubricants* **13**(1):7, published 29 Dec 2024, DOI
10.3390/lubricants13010007. **CC-BY.** This is the "Engels 2024" the old plans
refer to and that `validation/VALIDATION_REPORT.md` lists as **NOT RUN — the
direct Elrod observable; nothing currently tests it**. It is **experimental**, not
CFD.

**What they measured.** Fractional film content θ ∈ [0,1] — i.e. **void fraction
1 − θ, the direct Elrod/JFO state variable** — resolved over the full
circumference and the full 350 mm bearing width, simultaneously with film
pressure maps, film thickness, 123 thermocouples and calorimetric power loss.
Method: two co-located sensors in the rotating shaft, eddy-current for true gap
height and capacitive for effective permittivity, with θ recovered through a
Bruggeman effective-medium model. They claim it as the first quantitative in-film
void-fraction measurement.

**The result that matters for pancake.** **θ < 1 in regions where p exceeds
10 bar** (θ ≈ 0.8–0.9 there), with void fraction highest at pad leading/trailing
edges and at the axial edges. That falsifies a vapour-pressure cavitation
criterion outright and is direct measured support for a gas-transport model —
which is precisely the assumption `docs/CAVITATION_OIL_REFRIGERANT.md` argues for
from analogy.

**Why it is nonetheless not a benchmark, and this must be said explicitly.** The
bearing is an offset-halves design with two pads, leading-edge lube pockets, deep
grooves and a hydrostatic jacking pocket — unreproducible on a cylindrical θ×z
mesh. Worse, the companion CFD paper concludes the observed cavitation requires
**fluid inertia**, which the Reynolds equation drops by construction, and
attributes the mechanism to **pseudo-cavitation** (already-undissolved gas
expanding) rather than fresh desorption, having found Osterland's gaseous model
gave "insufficient desorption rate". **Do not promise a quantitative match.**

**Companion, for the record:** Wettmarshausen, Engels, Hagemann, Stottrop,
Weißbacher, Schwarze & Bender, *"An Experimentally Validated Cavitation Model for
Hydrodynamic Bearings Using Non-Condensable Gas"*, *Lubricants* **13**(4):140
(2025), DOI 10.3390/lubricants13040140, CC-BY. CFD only (VOF + non-condensable
gas, TEHD). Its Figs. 4b/5 and 8 are **line plots along θ at the bearing centre**,
much easier to digitise than the 2024 colour maps. **Note:
`docs/CAVITATION_OIL_REFRIGERANT.md` [R7] cites this paper under a wrong title
("Modeling Cavitation in Hydrodynamic Bearings with VOF and a Non-Condensable Gas
Phase") — fix it.**

**Recommended use.** Cite as measured evidence that gaseous cavitation with
θ < 1 at p ≫ p_vap is real, which is the load-bearing assumption of pancake's gas
model. State plainly that it is not a geometry pancake can reproduce.

**Effort.** Hours (a citation and a paragraph), not a run.

### 3.4 RANK 4 — gas-transport model constants (the calibration data nobody has)

`docs/CAVITATION_OIL_REFRIGERANT.md` §7 Q1/Q3 and `PLAN.md`:329–336 both flag
that the desorption/resorption rate constants and the two-phase µ closure are
**uncalibrated**. Two open-access sources supply numbers:

- **Reinke et al., *Fluids* 2020 **5**:208 / 2022 **7**:159 / 2023 **8**:44** —
  a scaled Couette journal-bearing model (R₁ = 146.44 mm, gap 3.62 mm, paraffin
  µ = 190 mPa·s), run at **100 % vs 25 % air saturation**, reporting time-resolved
  **gas release rate** as Schnerr–Sauer constants C_R (release) and C_A
  (absorption) versus rotation angle. Line plots, easy to digitise. This is the
  best available source for calibrating `gas_mass_transfer_rate`. Their finding
  **C_A = 0 — no reabsorption once bubbles form** — bears directly on
  `PLAN.md`:329–336, which records pancake's symmetric release/resorption rate as
  a known limitation. If C_A = 0 holds, pancake's symmetry is not a small
  approximation.
- **Wettmarshausen 2025** uses **Bunsen coefficient β_B = 0.09** (literature range
  0.08–0.1 for mineral oil), tuned by matching θ. pancake's Bunsen path
  (`fluid_properties.cpp`) can take this directly.

**Effort.** 1 day to digitise Reinke's C_R(θ) and set up a single-cell
release-rate comparison against `FluidProperties::gas_mass_transfer_rate`.

**First step.** This is the cheapest way to stop `gas_mass_transfer_rate` being an
unjustified number, and it is independent of every other item here. Do it
alongside §3.1.

### 3.5 RANK 5 — the LG R290 / PZ68S application (most relevant, no bearing data)

**What exists.** `data/R290_PZ68S/r290_pz68s.csv` — a real measured PTSV surface
(P 0.1–4.0 MPa, T −60…160 °C, solubility % and kinematic viscosity), ingested and
runnable via `config_r290_pz68_table.txt` with `solubility_model = TABLE`,
`liquid_viscosity_model = TABLE`, `cavitation_threshold = LOCAL_PSAT`. The shipped
constants put p_sat = 1.0 MPa at exactly the 10 bar feed pressure, so **the feed
enters saturated** and flashes on any pressure drop.

**What it can validate today: nothing, and that should be said plainly.** There is
no measured bearing observable. The PTSV table is the *input*, so reproducing it
is not validation. What is available is a physical *finding* rather than a
comparison: the saturated feed releases enough propane to fill ~23 % of the gap in
one 5 µs step, the ceiling drops by that much at once, and the load swings between
5.8e2 N and 2.4e5 N on alternate steps — **timestep-independent** at
dt = 5e-6, 5e-8 and 5e-9 s, so it is the model's displacement rate, not a temporal
artefact (`PHYSICS.md` §7.3.7). Every Newton solve converges; the load history is
transient noise for tens of steps. The Elrod path with `VOID_COUPLED` does not
resolve it either (it drifts monotonically instead).

**Why it matters anyway.** This is the application the whole gaseous-cavitation
line exists for, and the flash oscillation is a genuine, reproducible, honestly
characterised model behaviour. It is publishable as a *finding about
volumetrically-active gas models on a saturated feed* even without a comparison.

**What would unlock validation.** Any in-house LG measurement with a bearing
observable: friction power loss vs suction pressure, minimum film thickness, a
bearing-temperature trace, or a bubble-point measurement on the actual PZ68S
charge. **`[USER INPUT NEEDED]`** — this is the one item on this list that no
amount of literature work can supply.

**Effort.** Zero until data arrives; then 2–3 days.

**First step.** Ask what LG measured. Even a single (p_suction, power loss) pair at
two speeds turns this from a demonstration into a validation.

### 3.6 Cavity-pressure magnitude, and Braun & Hendricks (needs sourcing)

ASLE Trans. **27**:1, cavity pressure/temperature maps. **Not present** in
`F:\Seongho\Research\literature\pdfs` (236 PDFs, all checked). `cases/REFERENCE_DATA.md`
lists it as optional. It is currently cited in
`docs/CAVITATION_OIL_REFRIGERANT.md`:94 **second-hand** through San Andrés' TAMU
notes ("steady-loaded flooded cavities formed by gas coming out of solution") —
that characterisation is unverified against the primary source, which is a
citation-hygiene problem in its own right.

**Value.** A second, independent cavity-pressure dataset would let 05b
discriminate between "pancake's plateau is wrong" and "Etsion & Ludwig's tap
disturbs the cavity". Worth sourcing, not worth blocking on.

**First step.** Interlibrary request through KAIST (confirmed to exist: *ASLE
Trans.* **27**(1):1–14, DOI 10.1080/05698198408981539); meanwhile, downgrade the
`CAVITATION_OIL_REFRIGERANT.md`:94 citation to explicitly say "as reported by
San Andrés" until the primary is in hand.

**Two other cavity-pressure leads, both cheap to check:**

- **Shen & Khonsari, *Tribol. Lett.* **51**:153–160 (2013)**, DOI
  10.1007/s11249-013-0158-2 — reports **p_cav ≈ 30 kPa absolute**, far below
  ambient. Directly usable as a p_cav datum. **Caveat: paywalled and unread — it
  is not confirmed whether the experiments are theirs or a reanalysis of the
  literature. Verify before citing.**
- **Kasolang & Dwyer-Joyce, *Tribology Transactions* **51**(2):231–242 (2008)**,
  DOI 10.1080/10402000801947717 — ultrasonic film thickness around the full
  circumference plus cavitation onset location; D = 75 mm, c = 25 µm, Shell Turbo
  T68, 200–800 rpm, 3/6/10 kN. Paywalled at T&F but the **accepted version is free
  at eprints.whiterose.ac.uk/id/eprint/94787**. Line plots, easy to digitise.
  Gives an h(θ) and onset-angle target with no z-resolution and no void fraction.
  Their explicit finding that the cavitation boundary is **gradual, not sharp**
  ("air bubbles forming gently") is a useful counterweight to the ≤1-cell front
  precision pancake reports against the exact 1-D solution.

### 3.7 REJECT — the Jakobsson MATLAB CSV: do not use

`F:\Seongho\MATLAB\LGEP\jakobsson-experiment-data.csv` (18 rows, θ in deg 1.86→343.55,
p in Pa 19178→357808; duplicate under `case_data/data/`).

**Provenance: unestablishable.** It is provably a hand digitisation from a raster
plot — all 36 numbers are exact integer multiples of a fixed pixel scale (θ in
units of 1/29 deg, p in units of 10000/73 Pa), which fabricated or transcribed
data would not carry. But **which figure of which source is nowhere recorded**: a
full scan of `F:\Seongho` found only two MATLAB `readtable` calls and a legend
string `"Jakobsson"`. No Jakobsson or Floberg PDF exists anywhere in the library,
and `validation/journal_bearing/README.md`:334 already records that J&F 1957 is a
Chalmers technical report that could not be obtained — so the digitised source was
almost certainly a secondary reproduction that would have to be named before it
could be cited.

**Two independent reasons to drop it:** (i) an uncitable source cannot appear in a
validation report; (ii) **Jakobsson & Floberg is a vaporous-cavitation reference**,
so it is the wrong physics for this goal even if provenance were established.

Incidental finding worth recording: `plot_experimental_cfd.m` draws
`yline(20, ... 'p_{sat}')`, implying the companion simulation assumed
p_sat = 20 kPa. That modelling choice is documented nowhere else.

**First step.** Add a one-line note to `cases/REFERENCE_DATA.md` marking this file
as rejected-with-reason, so nobody re-discovers and re-uses it.

### 3.8 Dead ends, recorded so they are not searched again

- **No modern re-measurement of the classics exists.** Jakobsson & Floberg 1957,
  Etsion & Ludwig 1982, Braun & Hendricks 1984, Sun & Brewe 1991 and
  Heshmat/Pinkus are *still* cited as primary data in 2025 papers. Braun & Hannon,
  *Proc. IMechE Part J* **224**:839–863 (2010) is the standard review, not new
  data. Engels 2024 is the closest thing to a replacement and deliberately uses a
  different bearing type. **Conclusion: the 1982 data is not superseded, and using
  it is not a compromise.**
- **No experimental void fraction, outgassing pressure or cavity pressure for
  R290-in-oil in a bearing film could be found.** This is a genuine gap in the
  literature, not a search failure — which makes the LG in-house question (§3.5)
  more valuable, not less.
- **The Int. J. Refrigeration refrigerant-bearing paper is simulation-only.**
  Correct citation: Wang, Zhong, Lei & Wu, *Int. J. Refrigeration* **158**:253–263
  (2024), DOI 10.1016/j.ijrefrig.2023.12.001; R134a/RL10H at 0.2 MPa, homogeneous
  two-phase, equilibrium vs non-equilibrium desorption. It is [R4] in
  `docs/CAVITATION_OIL_REFRIGERANT.md`, cited there as "(2023)" with a bare
  ScienceDirect pii and described as "most directly on-target". **It contains no
  experimental data**, so the ~3 %→21 % load-capacity figure the cavitation doc
  quotes from it is a simulation result, not a measurement. Fix the citation and
  the framing.
- **X-ray / neutron radiography void-fraction work exists** (Argonne APS on
  nozzles; neutron imaging of orifices, 2024) but **nothing applied to a bearing
  film**. Not a near-term option.
- **"Kobayashi & Murata (1986)"**, cited secondhand in the refrigerant-bearing
  literature as having measured oil-film pressure vs groove type and
  refrigerant/oil mixing ratio, **could not be confirmed to exist** — no hit in
  Purdue e-Pubs or Crossref. **Do not cite without a copy.** A different, real
  Purdue paper (Kobayashi, Hayashi & Yoshimura, ICEC 1998, Paper 1258) is a
  friction/wear study with no pressure or void-fraction fields.
- Layton et al., *Lubricants* **14**(2):76 (2026) is a clean open-access TEHD
  half-journal dataset (in-film pressure, temperatures, ultrasonic film thickness,
  5–20 kN, 1000–4000 rpm) — but cavitation appears in it only as elevated
  measurement uncertainty. Useful for §2.5 (case 04, a second THD dataset), not
  for case 05.

### 3.9 Recommended shape of case 05

Do not look for one paper that does everything. Split it:

| Sub-case | Target | Source | Blocked? |
|---|---|---|---|
| 05a | cavity **extent** (θ_st, θ_end, θ₀) vs ambient pressure | Etsion & Ludwig Table 2, 8 conditions | **no — run it** |
| 05b | cavity **extent**, second geometry, flooded + starved | Ochiai 2019 photographs | no |
| 05c | **gas-release rate** constants | Reinke 2022 C_R/C_A; β_B = 0.09 | no |
| 05d | cavity **pressure** p/p_s vs Z/L | Etsion & Ludwig Figs. 6/7 | **yes — §4.3** |
| 05e | void fraction, physics claim only | Engels 2024 | no (not a benchmark) |
| 05f | R290/PZ68S application | LG in-house | **yes — no data exists** |

---

## 4. Known defects and limitations still open

Verified against `HANDOFF.md` (D1–D16), `PLAN.md`, `CHANGES.md` and the config
validator. Of the sixteen numbered defects, **only D6 and D14 have fix entries**;
D2 is partially fixed. Note there are **no TODO/FIXME comments anywhere in
`src/`** — all deferred work lives in the markdown, which is disciplined but means
the markdown must be kept honest.

### 4.1 Thermal runaway with thermally-transparent grooves — guarded, not fixed

`PHYSICS.md` §7.3.6 ("Known instability, both cavitation models"); `PLAN.md`:60–65;
`CHANGES.md`:699–709.

A groove that pins pressure but is left thermally transparent
(`feeds_fresh_oil = false`) diverges under a strong Andrade slope: Poiseuille
dissipation ∝ 1/µ with ∇p held by the pressure pin gives hotter → thinner →
hotter without bound. Measured **T ~ 1e6 K** on `grid_G3_fbns` at E_µ = 3700 K,
under FBNS *and* under JFO, at every under-relaxation from 1.0 down to 0.15.

**Status.** Workaround only: give the groove a supply temperature. **The
enforcement is documentation, not a config error** — the solver will happily run
the divergent configuration. Every case README repeats the warning by hand.

**Risk.** Someone runs a THD case with a pressure-only groove and gets a
plausible-looking answer before it blows up, or a mildly wrong one that does not.

**Fix.** Either (a) a validation error when `ENERGY_EQUATION` is on and an inlet
pins pressure without `feeds_fresh_oil`, or (b) the real fix in `energy.cpp`
(bound the Poiseuille dissipation at a pressure-pinned cell). (a) is 30 minutes
and removes the whole risk class; (b) is its own branch.

**First step.** Do (a). It is a strictly-better guard than a README sentence and
it costs nothing.

### 4.2 FBNS Schur-complement conditioning near the complementarity switch

`PLAN.md`:23–26; `CHANGES.md`:55–68.

The theory bound `cond(A_F) ≤ 4 + 2√2` holds by construction, but the **Schur
complement** goes near-singular whenever a complementarity pivot sits exactly on
its threshold — a last-bit difference flips the pivot and leaves a near-singular
row. The Elrod film-content operator has the same failure mode.

**Status.** Mitigated by the ASM preconditioner + `linear_solver_fallbacks` ladder
(`f0e772d`), which is in first production use and closed `rerun_cht_robin10`'s
force balance where it previously diverged on both hosts. The **conditioning fix
itself is open**: row scaling, or a regularised pivot near the switch.

**Live evidence it is not fully covered:** `test_fbns` fails at `-n 4` (Elrod
film-content KSP diverges) and **ctest only runs `-n 2` and `-n 6`**, so the
failure is unpoliced (`cases/01_ferron_thd_crossfilm/SUMMARY.md`:58).

**First step.** Add `-n 4` to the ctest matrix so the failure is visible, then
decide whether the ladder is sufficient. Making a known failure visible costs one
line and is a prerequisite for judging the rest.

### 4.3 `VOID_COUPLED` — non-convergent on Elrod, rejected under FBNS

Two separate problems that are easy to conflate.

**(a) On the Elrod path, `VOID_COUPLED` does not converge.** `R_p` 0.57–7.2
against a 1e-9 tolerance; answers are relaxation- and grid-dependent, which
proves they are not fixed points; the solver reports success on iteration-budget
exhaustion (`PLAN.md`:261–274, `PHYSICS.md`:1969 "**must not be used for…**").

**(b) Under FBNS it is rejected at config validation** (`src/config.hpp:1474`),
because `p_void` is advanced and under-relaxed **only** inside
`Reynolds::solve_elrod`'s outer flag-update loop, and FBNS has no such loop — so
it would run the plateau frozen at `p_cav` while reporting `VOID_COUPLED`.
Delivered as a documented limitation with a specific error message rather than a
half-implementation (`docs/FBNS_GAS_COUPLING.md` §8.3). That was the right call.

**Why this is now on the critical path.** §3.1 makes it load-bearing: a spatially
varying cavity pressure is exactly what the best gaseous-cavitation dataset
measures, and (a)+(b) together mean pancake cannot currently produce one.

**What it needs.** Not a bug fix — a numerical-methods change. Two candidates are
already recorded: fold the flash into the complementarity condition (i.e. make
`p_cav` a genuine per-cell field in the FB argument, which the memo shows is the
*easy* half), or Newton on the coupled `(p, θ, p_void)` residual. The hard half is
the plateau update site and its relaxation schedule.

**Effort.** 1–2 weeks, and it needs its own validation — not a side task.

**First step.** Decide whether 05b is worth it before starting. If yes, prototype
the per-cell `p_cav` array in `FBNS::cell_outflow`, the residual and
`axial_boundary_pressure` (mechanically easy per §4.2 of the memo) with `p_void`
held at a *prescribed analytic* field, to separate the solver work from the
closure work.

### 4.4 `LOCAL_PSAT` is a uniform scalar, not the per-cell field the physics calls for

HANDOFF **D4**. `docs/CAVITATION_OIL_REFRIGERANT.md` §9.2 item 1 calls promoting
`p_cav` to a field `p_sat(T, c_d)` "**the #1 recommended change**" and the
Tier-0 → Tier-1 step. `PHYSICS.md` §7.3.8 argues the scalar is *correct on the
FBNS path* — because `effective_p_cav()` evaluates the bubble point at the feed
state, so it is uniform by construction and injects no net force.

**Both are right, and the tension is the finding:** the scalar is the correct
*datum* for the model as implemented, and the model as implemented is Tier 0, not
Tier 1. **pancake does not currently implement the recommended gaseous-cavitation
model.** That should be stated in any write-up rather than left for a reviewer.

**First step.** One sentence in `PHYSICS.md` §16.5 (model provenance and honest
validity) saying the shipped threshold is a feed-state scalar and the per-cell
bubble-point surface is future work. Costs nothing, removes an overclaim.

### 4.5 The Elrod ceiling merge (`VoidParams::theta_full` vs `FullFilm::bound`)

`docs/FBNS_GAS_COUPLING.md` §8.4; HANDOFF **D2** (partially fixed).

`FullFilm::bound` is now the single definition of θ_full and moves only under
FBNS + `TWO_PHASE`. The Elrod operator still carries its own copy in
`VoidParams::theta_full`. Deliberately not merged: under `gas_pressure_coupling =
NONE` the Elrod free gas is volumetrically inert **by design** (the shipped model
of case 05), and under `VOID_COUPLED` merging would be a second independent
modelling change to a path that is already non-convergent — and it would break the
byte-identity that made the phase-3 change checkable.

**Status.** Correct call, correctly documented. Whoever merges them "should expect
to re-validate case 05, not just re-run it". Leave it until §4.3 is decided.

### 4.6 Bulk modulus β — the Elrod trap that does not announce itself

`validation/VALIDATION_REPORT.md` §2.4; `PLAN.md`:311–315, 235–237.

The accuracy of the Elrod rupture boundary is governed by **`p_max/β`, not by any
residual the solver reports**. At β = 1e9 on a 720-cell mesh the rupture angle is
**20.8° wrong** while `mc = 0`, `ns = 0`, `R_p = 0` and the process exits
successfully. Requirement: keep `p_max/β ≥ 4e-3`, and refine β *with* the mesh in
any grid study or the rupture error masquerades as grid error. Raising β toward
incompressibility makes agreement **worse** and stops the outer loop converging
(200+ iterations at β = 1e10 vs FBNS's 13).

**Status.** FBNS retires this entirely (`bulk_modulus` is unused, `PHYSICS.md`:716),
which is the main practical argument for making FBNS the default everywhere. The
trap remains live for any Elrod/JFO run — and case 05's shipped gaseous path *is*
JFO.

**First step.** Add the `p_max/β ≥ 4e-3` check to `SimulationConfig::validate()`
as a post-run warning, or at minimum to the case-05 protocol.

### 4.7 Gas transport is solved in non-conservative (material) form

`PHYSICS.md` §15.1, "two deliberate visibility limits": the dissolved-gas
transport is in material form, so `r_gas` in advecting cases includes that
discretisation choice rather than measuring a defect. Also: variable-property runs
refresh ρ_l after the segregated solve, so `r_l` then measures the property lag of
the sweep.

**Why it matters here.** This is a real conservation hole inside a solver whose
central verification claim is per-cell conservation (§6). For gaseous cavitation
the transported dissolved fraction is the state variable that sets the threshold,
so a non-conservative form is not cosmetic.

**First step.** Quantify before fixing: run a closed-domain gas advection case and
report the drift in `M_gas,tot`. If it is at round-off in practice, document it and
move on; if not, it is a prerequisite for §3.

### 4.8 Other open items worth naming

| Item | Where | Note |
|---|---|---|
| Groove hot-oil carry-over **off by default** | `PHYSICS.md` §13.1.1 | `groove_mixing_model = PINNED` (λ = 1) is the default and bit-identical to the old hard pin. The original defect (sized at "+8 K inlet = −22 % load") is still the shipped behaviour unless configured. Every case must state which it used. |
| D3 — initial θ uses raw `p_cav`, not `effective_p_cav()` | `src/main.cpp:185` | verified still live; matters exactly on the `LOCAL_PSAT` gaseous path |
| D1 — van Leer wrong for reversed flow; no z-face TVD | `HANDOFF.md`:186; `src/fvm.cpp:174-193` | **Reclassify: no longer latent.** `CHANGES.md`:1412 calls it "latent for omega > 0", which holds only for the Elrod θ equation (pure Couette face flux). On **energy** (`energy.cpp:542`) and **gas transport** (`gas_transport.cpp:283, :309`) the Poiseuille term can flip the face flux negative **with omega > 0**, and `VANLEER` passes straight through. `omega < 0` is also supported, which puts every Elrod face on the buggy branch. The campaign is safe only because it runs `UPWIND` everywhere. See `docs/CODE_AUDIT_2026-08.md` §1.3. |
| **Coefficient unification is incomplete** | `src/energy.cpp:245, :270`; `src/gas_transport.cpp:247-248` | `reynolds.cpp`, `diagnostics.cpp` and `fbns.cpp` route through `reynolds_coefficients.hpp`; **`energy.cpp`'s axial boundary flux helpers and all of `gas_transport.cpp` do not.** Under `cross_film_thermal_model` the energy equation transports heat through interior faces with Dowson's `F2` and through **axial boundary faces with the classical `h³/12µ`** — and the Ferron cross-film runs have Dirichlet axial ends, so the mismatched faces are the oil outflow. Gas is advected on a velocity field inconsistent with the pressure solve that produced it. See `docs/CODE_AUDIT_2026-08.md` §2.1. |
| D5 — stale thermal-property ghosts between gas and energy solves | `HANDOFF.md`:266 | MPI correctness; matters at 24 ranks |
| D7 — zero initial guess every `KSPSolve` | `src/linear_system.cpp:121` | see §5.1 — this is an acceleration item as much as a defect |
| D8 — Picard relaxation acts on p/T while the Elrod state is θ | `HANDOFF.md`:305 | `coupling_relaxation < 1` under-damps the variable that matters |
| No axial re-flooding of a cavitated film | `PLAN.md`:341–348 | "a **model limitation, not a physical prediction**". Ungrooved geometries need an explicit inlet as a numerical stand-in — **directly relevant to §3.1**, whose rig has no groove |
| Cavity pressure assumed uniform | `PLAN.md`:337–340 | see §3.0 |
| No vaporous-cavitation physics | `PLAN.md`:321–328 | "Do not quote the vaporous case as a validated prediction of vapour behaviour" |
| Reabsorption kinetics symmetric | `PLAN.md`:329–336 | one rate for release and resorption |
| Brito: 5 of 6 committed runs unconverged | `validation/VALIDATION_REPORT.md`:141 | must be re-run before anything Brito is quoted |
| No thermal/elastic clearance deformation | `validation/VALIDATION_REPORT.md` §5.3 | relevant to the clearance ambiguity in §2.1 |
| WP-8 dynamic coefficients (K_ij, C_ij) | `PLAN.md`:1299 | TODO; not required by the current ladder |

### 4.9 Record-keeping defects (cheap, and they cost credibility)

These are not physics, but a reviewer will find them.

1. **`validation/VALIDATION_REPORT.md` is stale.** It reports Ferron temperature
   peak *location* as "**FAIL** — 71° out, and **unchanged** by the conjugate wall
   model" and lists D14 and hot-oil carry-over as open gaps. All three have since
   changed: the 2-D conjugate bush puts the peak at 265.2° vs 258.2° measured,
   D14 is fixed, and `FLOW_WEIGHTED` carry-over landed. **Re-score it.**
2. **`PLAN.md`:2094–2108 (the Phase F table) is stale** — F1 (D6), F3 (D2) and
   F8 (D14) are marked TODO against landed fixes.
3. **`VALIDATION_PLAN.md`:5 says "D1–D16, all still open"** — contradicted by the
   D6 and D14 fix entries.
4. **B-number collision:** `HANDOFF.md`:417–442 and `VALIDATION_PLAN.md`:64–74 use
   two different B1/B2/B3/B4 sets, making "D14/B3" and "B4" citations ambiguous.
5. **Braun & Hendricks is cited second-hand** as if primary
   (`docs/CAVITATION_OIL_REFRIGERANT.md`:94) — see §3.3.
6. **`main` is 31 commits behind the tip of the stack.** Nine `codex/*` branches
   are unmerged: `fbns-cavitation` (24) → `full-film-predicate` (25) →
   `fbns-gas` (26) → `viscosity-polynomial` (27) → `fbns-crossfilm` (29) →
   `conjugate-bush` (30) → `thermal-ksp-robustness` (31), plus
   `validation-cavitation` (22) and `void-coupled-plateau` (4). Every result in
   this campaign — FBNS, the exact-1-D verification, the conjugate bush, the KSP
   robustness fix — lives on unmerged work, and the case archive is on a
   *different* branch from the solver that produced the later results. The repo
   rule is "one task, one branch, one verified merge, then delete the branch".
   **This is the largest process risk on the list.**

**First step for §4.9.** Items 1–4 are an afternoon and they are the difference
between a document a reviewer trusts and one they audit. Item 6 needs a user
decision on merge order.

---

## 5. Speeding up the campaign

Measured baseline: the archived force-balance work is **196 solver runs,
6227 s total, mean 31.8 s** (summed from the `fb_log.json` files). Ordered by
payoff ÷ effort.

**Framing fact that changes how you read the rest.** The force-balance loop is
**not in the solver**. It is a Python driver
(`cases/01_ferron_thd/force_balance.py:126-158`) that spawns a fresh `mpiexec`
per iterate (`solver_io.py:31-45`), and every iterate initialises fields to
uniform constants (`src/main.cpp:186-208`). There is **no restart, resume or
field-load path anywhere** in `io.cpp`, `config.hpp` or `main.cpp`. The existing
"warm start" is *parameter* warm start (reuse of identical `(cfgkey, eps, psi)`
records, pass 2 starting from pass 1's converged point, per-load seeding) —
**field warm start does not exist**.

### 5.1 P1 — the Picard limit cycle in the conjugate cases (~35 % of the batch, levers already in the code)

`coupling_max_iters` defaults to 30 (`config.hpp:918`) but the driver overrides it
to **200** (`force_balance.py:94`). Across 180 records: min 6, median 11, mean
26.4, max 200. Twelve runs sit exactly at the cap. The split is perfectly clean by
wall boundary condition:

| group | bush outer BC | Picard mean | max |
|---|---|---|---|
| `outer_fixed313` | Dirichlet | 28.0 | 31 |
| `outer_robin80` | Robin h = 80 | **196.2** | 200 |
| `rerun_cht_robin10` | Robin h = 10 | **200.0** | 200 (all 8 runs) |
| all non-conjugate THD | — | 6–16 | 16 |

**It is a genuine cycle, not slow convergence.** Sampling `|G|` at each outer's
first Newton evaluation in `cases/01_ferron_conjugate/outer_robin80/fb_iter_00/run.log`:
`1.12e-4` → `1.2e-7` at outer ≈ 51 → then **rises** to 3–4.5e-7 and oscillates for
the remaining 150 iterations. **The loop returns iterate 200, which is worse than
iterate 51.** `cases/01_ferron_conjugate/progress.log` independently records
`R_p = 1.06e-3 / 9.83e-4 / 9.73e-4` at 100 / 200 / 400 outers — flat — and that it
vanishes at 360×96, so it is grid-dependent.

**Both relaxation levers already exist and are set to 1.0, i.e. off:**
`coupling_relaxation` (`config.hpp:920`, applied at `main.cpp:99-104, 388-389`)
and `cross_film_relaxation` (`config.hpp:755`). Both cycling configs run at 1.0.

- **Payoff:** ~2170 s of the 6227 s batch (35 %) is spent past the useful iterate.
  Matching the `outer_fixed313` count would be ~6× on those runs.
- **Effort:** hours. Sweep `coupling_relaxation ≈ 0.5–0.7`; add stagnation
  detection and keep-best-iterate at `main.cpp:384-392`.
- **Risk:** low. **Relaxation before Anderson** — the signature is loop gain
  slightly above 1 with a reversing trend, which is exactly what a relaxation
  factor fixes. Anderson needs new code and joint p/T mixing; do not reach for it
  until the free lever has been tried.
- **Correctness note, not just speed:** returning iterate 200 when iterate 51 was
  better means **the archived conjugate results are not the best iterate the
  solver found**. Worth checking whether it moves any reported number.

### 5.2 P2 — `Fields` string lookups inside per-cell loops

`Fields` is `std::map<std::string, unique_ptr<Field>>` (`field.hpp:33`); every
call site passes a `const char*` literal, so each lookup constructs a temporary
`std::string` and walks a red-black tree. libstdc++ SSO is 15 chars, so these
names **heap-allocate on every lookup**: `rho_cp_liquid_solution` (22),
`cavitation_clamp_mass` (21), `cross_film_couette_h` (20),
`rho_liquid_solution` (19), `mu_/cp_liquid_solution` (18), `k_liquid_solution` (17).

Worst body: **`FBNS::cell_outflow` (`fbns.cpp:137-191`)**. Its `gamma_p` lambda
resolves through `ReynoldsCoeff::fbns_mobility` to 6 map lookups plus 2 mallocs,
and is called 8× per cell, plus ~10 direct lookups — roughly **66 map lookups and
30 malloc/free pairs per cell**, inside the Newton loop at `fbns.cpp:436`. The
file's own comment at `fbns.cpp:240-247` states the coefficients are frozen
through the whole Newton solve, so hoisting the mobility into a precomputed
`Field` is **provably safe by its own argument**. Second worst:
`FluidProperties::update_solution_fields` (`fluid_properties.cpp:274-330`), called
**three times per outer iteration** (`main.cpp:299, 339, 364`).

- **Payoff:** structurally the largest raw win in the codebase — but **no profile
  exists and I will not put a number on it.** A malloc/free pair plus a 5-deep
  string compare is order 50–100 ns against maybe 100–200 ns of actual arithmetic
  per cell, which *suggests* lookups dominate `cell_outflow`, possibly by an order
  of magnitude. Treat as a hypothesis worth 30 minutes of `perf` before committing
  days.
- **Effort:** 1 day for the two hot routines (hoist `Field&` references above the
  loops); 2–3 days for the structural fix (a handle type so the coefficient
  helpers take `const Field*` rather than `const char*`), which removes the class.
- **Risk:** low. Pure mechanical hoisting, bit-identity preserved. **Survives
  `-O3`** — the lambdas inline, the map lookups and mallocs do not.

### 5.3 P3 — field warm starting / restart files (does not exist)

The driver's Newton makes 3 solver calls per step: base, a `J_eps` probe at
`eps + 5e-3`, and a `J_psi` probe at `psi + 0.5°` (`force_balance.py:142-143`).
The probes are **~0.8 % perturbations of an already-converged state, and each one
restarts from a uniform field.**

- **Payoff:** the two Jacobian probes are **2/3 of all 196 solver calls**. Starting
  a probe from the base iterate's fields should drop its Picard count from ~8–12 to
  ~2–3 — roughly **2× on the force-balance batch**, more on the conjugate cases.
- **Effort:** 2–3 days. Write physical cells of
  `pressure`/`theta`/`temperature`/`dissolved_gas`/`free_gas_mass` to a binary
  restart; add an `initial_fields_from` config key; load after field registration
  and before `store_old_time` (`main.cpp:253`).
- **Risk:** low–moderate. Main hazard is a stale restart silently seeding a
  different configuration — guard with a config-hash header in the file.
- **Bonus:** the same file format is the prerequisite for §5.5, and it would also
  make the 983 MB of afml-r6 VTK archives re-startable rather than re-runnable.

### 5.4 P4 — over-solved linear systems, and the zero initial guess (defect D7)

`linear_rtol = 1e-8` (`config.hpp:701`), never overridden by any case template.
`LinearSystem::solve` re-does everything on every call: `MatZeroEntries` plus
per-entry `MatSetValue` re-assembly (`linear_system.cpp:154-195`),
`KSPSetOperators` (`:213`, which forces PC re-setup so ILU(0) is refactorised
every solve), and `VecZeroEntries(x)` (`:214`). **`KSPSetInitialGuessNonzero`,
`PCSetReusePreconditioner` and batched `MatSetValues` appear nowhere in `src/`.**

Direct evidence of waste, repeated ~200× per run in
`outer_robin80/fb_iter_00/run.log`:

```
fbns iter=0 |G|=4.489276e-07
fbns iter=1 |G|=3.150233e-14
```

One Newton step drives the mass residual **8 orders below what the outer loop can
consume** (`coupling_tolerance = 1e-6`), and a second Newton iteration runs purely
to confirm it.

- **Effort:** hours for the nonzero initial guess (this is **defect D7**, already
  flagged "MEDIUM (cheap, big win)"); 1 day for an inexact-Newton `linear_rtol`
  schedule plus validation that converged answers are unchanged.
- **Risk:** moderate, and specifically: loosening `rtol` interacts with
  `reynolds.cpp:454`, where `snap_tol = 5e-7` was chosen *because* `rtol = 1e-8`
  gives errors ~1e-5 near the free boundary. **Loosen one and the other must move,
  or flag flips increase.** Freezing the preconditioner is the riskier sibling:
  `linear_system.hpp:24-30` documents that these operators go near-singular at
  complementarity pivots (§4.2), which is why the fallback ladder exists — a frozen
  ILU would raise the fallback rate.

**Prerequisite, and it is 5 lines.** KSP iteration counts on *successful* solves
are **never printed** — `linear_system.cpp` reports `n_iters` only in the failure
and fallback branches (`:242, :251, :259`). There is no distribution of normal
Krylov effort anywhere in `cases/` or `validation/`. **You cannot tune this
without adding that instrumentation first.** Do it before anything else in §5.

### 5.5 P5 — coarse-to-fine grid cascading

No multigrid, `PCMG`, `DMDA`, prolongation or inter-mesh interpolation exists.
But the grids are uniform and **exactly nested**: G3 240×64 → G6 720×192 is exactly
3× in both directions, so prolongation is trivial to write. Cost ratio from the
00a wall-clock table: `grid_G3_thd` 3.6 s vs `grid_G4_thd` 9.5 s, against 14–52 s
for the G6 THD force-balance iterates — a coarse balance costs roughly 1/10.

**Important caveat from the 00a sweep.** Load is converged from G3 onward
(W = 11762.8 / 11759.1 / 11752.3 across G3/G4/G5, 0.09 % spread) but **force angle
is not**: −134.564 / −134.130 / −133.746, drifting ~0.4° per refinement against a
force-balance tolerance of 0.05° (`force_balance.py:48`). So a coarse balance
lands ε well but ψ about 1° off — a good initial guess given the Newton damping
caps (0.15 in ε, 20° in ψ), but **never the final answer**.

- **Payoff:** ~2× on the batch, but **largely overlapping with §5.3 rather than
  stacking with it.**
- **Effort:** 3–4 days, and it **requires §5.3 first** (needs field I/O).
- **Risk:** moderate — interpolating a cavitation front between grids can seed an
  inconsistent θ/p complementarity state. Safest: interpolate **pressure only** and
  let FBNS rebuild θ.

### 5.6 P6 — the loop nest order is backwards nearly everywhere

`Field::operator()` indexes `data[(j + n_ghost) * stride_theta + (i + n_ghost)]`
(`field.hpp:23`), so **`i` is the contiguous axis**. Nearly every solver loop is
`for i { for j {`, so the inner loop strides by `stride_theta` doubles: one cache
line per element, no vectorisation. ~50 sites across `reynolds.cpp`, `energy.cpp`,
`fbns.cpp`, `fvm.cpp`, `gas_transport.cpp`, `cross_film.cpp`, `diagnostics.cpp`,
`fluid_properties.cpp`. Only `film_thickness.cpp` and `io.cpp` have the
cache-friendly order.

- **Payoff:** unquantified. Memory-bound kernels typically gain 2–4×, but here the
  map lookups (§5.2) currently dominate and mask it. **Do §5.2 first, then
  re-measure** — fixing order before lookups will show almost nothing.
- **Effort:** 1–2 days, mechanical.
- **Risk:** low, with one real catch: floating-point **reduction order changes** in
  any accumulating loop, which breaks the bit-identity guarantees asserted at
  `reynolds_coefficients.hpp:31-34`. Check the regression suite.

### 5.7 P7 — thermal timestep relaxation (weakest item; deprioritise)

Answering the question as posed: in `STEADY_STATE` the energy equation has **no
time term at all** (`energy.cpp:592-593` gates `ddt_weighted` on non-steady), so
every force-balance run has no thermal dt — thermal cost is Picard iterations,
i.e. §5.1. In `TRANSIENT`, energy uses the same `cfg.dt` as Reynolds with no
sub-cycling anywhere.

**The counterpoint against pursuing it.** The 00b cost table shows per-step cost
*falls* as dt falls, because Newton iterations per step drop (mean 20.18 → 2.03 →
1.81 → 1.62 across T1–T4), so T4 did 4× T2's steps in 3.6× the wall clock. At the
production dt = 2e-5 the FBNS Newton mean is **2.03** — the flow solve is already
cheap and the single energy solve is a comparable share. Sub-cycling would at best
remove ~1 of ~3 linear solves per step, so the ceiling is ~30 %.

**And the premise is untested:** no run in the archive varied thermal dt
independently, so *"thermal transients may tolerate 5–10× the mechanical dt"* is a
hypothesis, not a finding. It needs an experiment (~half a day) before the 2–3 days
of implementation — and the implementation's hard part is time-averaging
`heat_generation` over skipped steps, which is easy to get subtly wrong.

### 5.8 P8 — the 24-rank ceiling is hardware, not the decomposition

Already measured on afml-r6 (2× Xeon E5-2660 v4, **28 physical cores**), marginal
s/step: n=4 1.783 (100 %), n=8 0.973 (92 %), n=16 0.614 (73 %), **n=24 0.410
(73 %)**, n=32 0.548 (41 %), n=48 0.416 (36 %), n=56 0.449 (28 %). **The knee is
the 28-physical-core boundary** — n=32 is slower in absolute wall clock than n=24.
At n=24 each rank owns 30 θ cells × 192 z, nowhere near the algorithmic limit.

Hard algorithmic limit of the current code: `communicator.cpp:31-33` exchanges only
with immediate neighbours, so a rank needs `n_theta_local >= n_ghost = 2`, i.e.
**360 ranks at 720 θ**. Note `mesh.cpp:14-23` has **no guard** against
`n_theta_local < n_ghost` or `== 0` — a latent silent-corruption bug worth a
one-line assert regardless of any scaling work.

**Cost of 2-D decomposition:** `mesh.cpp:10-11` (`n_z_local`/`offset_z` hardcoded
to global), `communicator.cpp:47-54` (2-neighbour `Sendrecv` → 4-neighbour plus
corners), `io.cpp` (per-rank VTK z extents), and critically
`linear_system.cpp:138-140`, where `global_row = global_theta * n_z_global + j`
**assumes each rank owns complete z columns** — which is what makes its rows
contiguous and satisfies the `MatSetSizes(n_local, ...)` ownership contract at
`:52`. **Yes, PETSc's matrix distribution is tied to the θ-only layout**, and
breaking it means a real row remapping or a `DMDA`.

- **Payoff: zero on current hardware.** Only pays above ~28 cores.
- **Effort:** 1–2 weeks. **Risk:** high.
- **Verdict: do not do this until §5.1–§5.4 are exhausted.** More ranks is the
  wrong answer to a 35 % limit cycle.

### 5.9 A build-configuration footgun (not a speedup, but fix it)

`CMakeLists.txt` never sets or defaults `CMAKE_BUILD_TYPE` — no `-O`, no `-march`,
no LTO anywhere in the build system. The documented Linux command in `CLAUDE.md:24`
and `AGENTS.md:24` is `cmake -S . -B build-linux` with **no build type, which
compiles at `-O0`**; `BUILDING.md:116`, the only Linux command that sets a type,
sets `Debug`.

**Be honest about what this is and is not.** The production Windows binary *was*
built Release with `-O3 -DNDEBUG` (`build-windows-mingw/CMakeCache.txt`,
confirmed in `compile_commands.json`), so **every measured number in this
repository came from an optimised build.** This is a trap waiting for the next
Linux user, not 5–15× sitting unclaimed. Fix: one line plus a doc correction. Also
note presets `release` and `release-mingw` set `PANCAKE_BUILD_SOLVER=OFF` —
anyone benchmarking "the release preset" gets the GUI and no solver.

### 5.10 Suggested order

Add KSP iteration logging (minutes, and a prerequisite for §5.4) → **§5.1**
(hours, 35 % of the batch, levers already in the code) → §5.2 hot routines with a
profile first (1 day) → §5.4 initial guess (hours) → §5.3 restart (2–3 days) →
§5.5 cascading → §5.6 loop order → stop. §5.7 and §5.8 do not earn their effort
yet.

---

## 6. Is finite volume required for mass conservation?

Short answer: **no — a conservative flux form is required, and finite volume is
the cheapest way to be sure you have one.** pancake should keep it.

### 6.1 What pancake actually does

Finite volume, cell-centred, on the structured N_θ × N_z grid
(`PHYSICS.md` §5, `ARCHITECTURE.md`:12 "finite-volume operators"). Each cell is a
control volume; the divergence theorem is applied on the cylindrical surface;
face mobilities are **harmonically averaged**, `γ_e = 2γ_iγ_{i+1}/(γ_i+γ_{i+1})`;
the assembled stencil is `a_P p_P − Σ a_f p_nbr = source` with
`a_P = Σ a_f`. Time is backward Euler. Both cavitation solvers and the
diagnostics read their face coefficients from **one** header,
`src/reynolds_coefficients.hpp`.

### 6.2 Why conservative fluxes matter *specifically* for JFO/FBNS

Three reasons, in increasing order of how much they should change your mind.

1. **The free boundary is a mass-flux statement.** JFO says the Couette flux
   entering a rupture front equals the flux leaving it. A non-conservative
   discretisation puts a spurious source at exactly the cell where the physics
   lives, so the front lands wherever the artificial source balances — and it
   will still converge, and still report a small residual.

2. **Telescoping is what makes §15.1 a *check* rather than a claim.** Because
   interior faces cancel pairwise between neighbours, summing the per-cell
   residual over the domain returns the global closure *exactly*:
   `Σ_P r_P = S_storage − Φ_boundary − S_inlet − M_clamp/Δt`. That identity is
   the whole reason `mass_residual` can distinguish a cell-by-cell conservative
   scheme from large local errors that cancel on summation.
   `tests/test_mass_residual.cpp` demonstrates the point by perturbing one
   converged cell: the **global `r_l` does not move at all** while the field rises
   to 2e3× the Couette flux in the five affected cells. Without a face-flux form
   you cannot construct that test.

3. **It is how the FBNS-vs-Elrod comparison was made trustworthy.** The
   `mass_residual` work ruled out a non-conservative flux or boundary treatment
   *first* (front cells at 8e-11 of peak Couette flux under Elrod), which is what
   left the free-boundary treatment as the only remaining candidate for the 7–13
   cell front error. Conservation was the instrument, not the conclusion.

### 6.3 Is FVM *required*?

No. What is required is that **each interior face flux is computed once and
applied to its two adjacent cells with opposite signs.** Finite volume delivers
that by construction. A finite-difference scheme written in *conservative*
(flux-difference) form on a uniform grid is algebraically identical to FVM and
conserves equally well. A finite-difference scheme in *non-conservative* form
(expanding `∂(ρhθ)/∂θ` by the chain rule) does not, and will misplace the JFO
front regardless of which cavitation algorithm sits on top.

**FBNS in particular is discretisation-agnostic.** The Fischer–Burmeister
reformulation is pointwise and algebraic: it acts on the discrete residual `G`
and the Jacobian blocks `A` (Poiseuille) and `B` (upwind Couette + backward-Euler
storage), and it never inspects how `A` and `B` were assembled. So FBNS does not
require FVM — it requires `G` to be a genuine discrete mass balance, which is the
same conservative-flux requirement stated a different way.

*Unverified:* which spatial discretisation Woloszynski, Podsiadlo & Stachowiak
(Tribol. Lett. 58:18, 2015) themselves used. The paper is paywalled and the
abstract describes only the reformulation and the solve. **First step:** check
whether the PDF is in `F:\Seongho\Research\literature\pdfs` and read §2; it is a
one-line citation fix either way, and worth getting right before publishing a
claim about it. Note there is a sibling Tribology Letters paper the same year,
*A General Finite Volume Method for the Solution of the Reynolds Lubrication
Equation with a Mass-Conserving Cavitation Model* (10.1007/s11249-015-0588-0),
which suggests the FV framing was being developed alongside.

### 6.4 The practical answer for pancake — keep it, for three specific reasons

1. **It is what makes the verification chain real** (§6.2, point 2). Drop the flux
   form and `mass_residual` becomes uninterpretable.
2. **Non-uniform grids are where FD and FVM stop being equivalent.** Any future
   grid stretching toward the groove or the axial ends — likely, given §2.3 and
   §3.1 both need axial resolution — keeps conservation under FVM and loses it
   under naive FD.
3. **One flux definition is enforceable; four are not.** The 2026-08 finding is
   the empirical argument: four drifted copies of the face coefficients had
   accumulated in the diagnostics, and on exactly the cross-film and
   `VOID_COUPLED` runs the field exists to police, they reported a **fabricated
   ~10 % imbalance**. Unifying onto `src/reynolds_coefficients.hpp` dropped the
   `VOID_COUPLED` mass residual **9.4×** with the solved fields byte-identical —
   i.e. 90 % of the previously reported imbalance was the diagnostics mismatch,
   not the solver. A flux-based formulation is what makes a single shared
   definition possible at all.

   **But the unification is not finished, and that qualifies this point.**
   `reynolds.cpp`, `diagnostics.cpp` and `fbns.cpp` honour the header;
   `energy.cpp`'s axial boundary flux helpers and the whole of
   `gas_transport.cpp` still hand-roll `h³/(12µ)` (§4.8). So the guarantee
   currently holds for the *mass* balance and not for the enthalpy or gas
   fluxes — which is exactly the drift the header exists to prevent, reappearing
   on the two paths the campaign most recently extended. Closing it is a
   prerequisite for trusting a cross-film energy balance, not a tidy-up.

### 6.5 Two honest caveats

- **Conservation is necessary, not sufficient.** The Elrod path is fully
  conservative to round-off and *still* put the rupture front 7–13 cell widths
  wrong. Anyone reading "mass-conserving" as "accurate" is making the mistake this
  campaign spent case 00c disproving.
- **pancake is not conservative everywhere.** The dissolved-gas transport is
  solved in **material (non-conservative) form**, so `r_gas` in advecting cases
  includes that discretisation choice (`PHYSICS.md` §15.1). That is a real hole in
  the middle of the gaseous-cavitation path — see §4.7.

---

## 7. Elrod–Adams vs FBNS vs LCP — recap for the docs

All numbers below are measured in this repository unless marked *(lit.)*.

| Axis | Elrod–Adams (JFO), `cavitation_model = JFO` | FBNS, `= FBNS` | LCP / Giacopini *(lit.)* |
|---|---|---|---|
| Free-boundary treatment | switch function `g ∈ {0,1}` read through `p − p_cav`, which is `O(β⁻¹)` small at the front because `p = p_cav + β ln θ` | Fischer–Burmeister `FB(a,b) = a+b−√(a²+b²)` imposed pointwise; no flag array, no active-set loop | complementarity as a linear complementarity problem, solved by pivoting (Lemke) |
| Mass conservation | per-cell, at round-off (`mass_residual` 1.1e-7 of peak Couette flux worst case, = the Krylov tolerance) | per-cell, at round-off (**3.3e-14** of peak Couette flux on the front cells themselves) | mass-conserving by construction |
| **Rupture-front error vs the exact 1-D solution** | **first-order in degrees, ~7–13 *cell widths*** — does not converge relative to the mesh | **≤ 0.92 cells at every resolution**; reformation ≤ 0.66 cells | not measured here |
| p_max error @ 480 / 960 cells | +24.0 % / +12.2 % | +0.14 % / +0.07 % | — |
| W′ error @ 480 / 960 cells | +31.1 % / +17.8 % | +0.07 % / −0.04 % | — |
| Free parameters | **bulk modulus β**, and the front accuracy is governed by `p_max/β`, not by any reported residual (20.8° wrong at β = 1e9, 720 cells, all indicators green). Raising β makes it **worse** and stops convergence | **none** — `bulk_modulus` is unused | none |
| Nonlinear iterations | outer active-set flag loop: 13 at β = 1e9, **200+ at β = 1e10** | **4 Newton iterations, grid-independent** from N_θ = 120 to 960; last step drops 7 orders (finite termination once the active set stops changing) | pivot count grows with problem size |
| Linear-solve cost per iteration | one 5-point sparse solve | one 5-point sparse solve on the Schur complement, which **equals `B_G` at the solution** — i.e. the cost of an uncavitated pressure solve | pivoting, not a single sparse solve |
| Conditioning | β-conditioned | `cond(A_F) ≤ 4 + 2√2` by construction; **but** the Schur complement is near-singular when a pivot sits on the switch (§4.2, open) | — |
| Boundary re-flooding | needs the `consistent_boundary_flux` patch — the barotropic link overstates re-flooding by ~`β/(p_bc − p_cav)` | physical half-cell Poiseuille link, correct for free | — |
| Enthalpy transport in the cavitated arc | Poiseuille enthalpy flux must be **gated** by the full-film switch, because the barotropic `p` there is not the physical pressure | never gated — a ruptured cell is pinned at `p_cav`, so every gradient it produces is physical and the energy equation convects on the same velocity field that carried mass | — |
| Higher-order θ scheme | `TYPE_DIFFERENCING` (Vijayaraghavan & Keith), `VANLEER`, `LINEAR` available | ignored — upwind only (a higher-order stencil would add θ-second-neighbour columns to `B`) | — |
| Per-cell cavity pressure | `VOID_COUPLED` field exists, **does not converge** | rejected at validation (§4.3) | — |
| Two-phase ceiling θ ≤ 1 − α_gas | not moved (Elrod free gas is volumetrically inert under `NONE`, by design) | **landed**; complementarity invariants exactly 0 | — |

**What the campaign measured on the Ferron case.** Two distinct numbers that are
easy to conflate:

- On the **exact 1-D benchmark** (`cases/00c_jfo_verification/`), JFO is still
  **+12.2 % in p_max and +17.8 % in W′ at 960 cells**, where FBNS is +0.07 % and
  −0.04 %. This is a verification result against a known answer.
- On the **Ferron-like bearing**, Elrod's rupture front sits **12–18° from the
  exact 1-D solution**, and that gap is a property of the barotropic smearing, not
  of the grid: it stays near 13° from N_θ = 240 to 720 (`PHYSICS.md` §7.3.6).
  The residual FBNS-vs-Elrod temperature difference on the shared full-film arc is
  **0.44 K** at ε = 0.3, scaling with the thermal rise (0.24 K at ε = 0.2, 1.42 K
  at ε = 0.6) — that is Elrod's misplaced front leaking into its pressure field and
  hence into the Poiseuille share of enthalpy transport (~25 % of the total).

**Why FBNS was adopted, in one sentence.** The Elrod error is not a tuning problem
— `mass_residual` proved the flux and boundary treatment are conservative to
round-off, which left the free-boundary treatment itself as the only candidate,
and FBNS removes the mechanism rather than tuning it.

**Fronts are deliberately not cross-compared between the two solvers.** Elrod's
front is 12–18° from exact on this bearing, so agreeing with it would be evidence
of a defect. The exact 1-D benchmark is the only front reference that means
anything.
