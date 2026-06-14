# AUDIT.md — Physics, Boundary-Condition & CFD-Methods Audit of `pancake`

> **What this is.** A literature-grounded audit of the physical models, boundary
> conditions, and numerical methods in `pancake` (2D FVM Reynolds journal-bearing
> solver with Elrod–Adams JFO cavitation, oil+refrigerant mixture model, film
> energy equation, moving bearing). It answers three questions per issue: **what
> the codebase does** (verified `file:line`), **what the research says** (primary
> literature, quoted), and **what needs to be done**.
>
> **Inputs.** (1) The two-agent review thread `REVIEW.md` (Reviewer A: W1–W11,
> M1–M4; Reviewer B: N1–N9) and the consolidated `REVIEW_VERDICT.md`; every
> load-bearing code claim from that thread was re-verified directly against the
> source for this audit. (2) A multi-agent deep-research pass over the tribology /
> refrigeration / CFD literature: 21 sources fetched (ASME *J. Tribology*, STLE
> *Tribology Transactions*, *Wear*, *Tribology International*, *Int. J.
> Refrigeration*, *Tribology Letters*, Proc. IMechE Part J, San Andrés TAMU
> notes), 97 claims extracted with direct quotes.
>
> **Citation integrity note.** Every literature statement below traces to a
> **direct quote** pulled from the named source during the research pass (21
> sources fetched, 97 claims extracted with verbatim quotes and source URLs). An
> adversarial 3-vote verification was attempted on the top 25 claims across two
> runs; both were cut short by session rate limits, so the adversarial pass is
> **incomplete**. Three claims were formally confirmed before the limits hit —
> Braun & Hannon's gaseous/vaporous/pseudo regime taxonomy (2-0), Etsion &
> Ludwig's sub-cavity pressure measurements (3-0), and Braun & Hannon's
> "thermodynamic aspects" framing (3-0). The remaining claims rest on their
> extracted source quotes, not on model memory; the verifier "refutations" in the
> run logs are rate-limit abstentions (0-0), **not** substantive rejections, and
> none touches a finding here. Where a number sits behind a paywall (only the
> abstract was readable), that is flagged inline. Treat the unverified quotes as
> high-confidence-but-unaudited: spot-check the load-bearing numbers (the −3%/+21%
> load shift, the 3–40% β-softening error, the 7.13% vs 21.75% dissolved-fraction
> comparison) against the primary PDFs before citing them externally.

---

## 0. Executive summary

The Reynolds/Elrod–Adams core is canonical and is not re-litigated here. The
audit finds **two model-level defects that contradict (not merely simplify) the
literature**, both in the gaseous-cavitation path:

1. **The free-gas field is volumetrically inert** — it neither displaces liquid
   nor relieves pressure, while every published refrigerant-oil and
   non-condensable-gas bearing model couples the gas into film continuity
   (§ A2). The published magnitude of the effect this suppresses is ~**−3% load
   at ε = 0.3 and +21% at ε = 0.9** (Int. J. Refrigeration 2023).
2. **Cavitation onset is pinned to a fixed vaporous `p_cav` (0.3 bar) while the
   film is gas-saturated at 10 bar** — fifty years of experiments (Etsion &
   Ludwig 1982; Braun & Hendricks 1984) and every refrigerant-bearing model
   (Grando–Priest–Prata 2006 →  Int. J. Refrig. 2023) put gas-laden void onset
   at the **mixture saturation condition**, with measured cavity pressures
   near-ambient and operating-condition-dependent, not at a fixed sub-bar
   threshold (§ A1).

Below those, in decreasing severity: a **single-sweep "steady state"** that no
published THD/two-phase code would report as an operating point (§ C4); an
**Einstein dilute viscosity law driven to α_g → 1**, where it is directionally
wrong (§ A3); a **THD level below the 1983–84 validated baseline** while the
default config presents an energy solve (§ A5); a **β-conditioning problem the
literature has documented and solved** three different ways (§ C2–C3); **dead
numerics** (TVD, time-method selectors) advertised by docs/GUI (§ C1, C5);
**doc↔config drift on the flagship R290 case**, whose 10-bar sump sits exactly
at the model's own saturation pressure (§ B1); and a **validation suite with
zero quantitative benchmarks** despite the literature shipping ready-made ones,
including one with downloadable reference source code (§ E).

One review recommendation is **corrected** by the literature: the finite-rate
outgassing kinetics `k_dg(c_d − c_eq)` is *not* a defect — non-equilibrium
release/absorption has explicit experimental motivation (Etsion & Ludwig's
p/p_s > 1 observation) and peer-reviewed precedent (Hao & Gu 2014). Keep the
kinetics; fix what the released gas is allowed to do (§ A2.3).

---

# Part I — Physical models

## A1. Cavitation onset: fixed vaporous `p_cav` for a gas-saturated film — **UNPHYSICAL for the R290 case**

**What pancake does.** JFO rupture occurs at a single constant `p_cav`
(`config_r290_pz68.txt:31`: `p_cav = 3.0e4` Pa absolute = 0.3 bar), independent
of temperature and dissolved-gas content. The cavitated-zone pressure is
identically `p_cav` (`PHYSICS.md` EOS: `p = p_cav + β ln θ`, g = 0 ⇒ p = p_cav).
Meanwhile the model's own Henry law (`fluid_properties.cpp:76–85`) puts the
saturation pressure of the shipped mixture at
`p_sat = c_d/H = 0.2175 / 2.175e-7 = 1.0e6 Pa = 10 bar` — so the film is
supersaturated and releasing gas across the entire 10 → 0.3 bar band while JFO
still reports "full film."

**What the research says.**

- Cavitation in bearings is **three distinct regimes** — gaseous, vaporous,
  pseudo — not one phenomenon with one pressure. *(Braun & Hannon 2010, Proc.
  IMechE Part J, DOI 10.1243/13506501JET772 — review keywords list exactly these
  regimes; adversarially confirmed 2-0.)* The review explicitly frames cavitation
  modelling as requiring "attention to the thermodynamic aspects" grounded in
  experiment, not just algorithmics.
- **Etsion & Ludwig 1982** (*J. Lubrication Technology* 104(2):157–163;
  confirmed 3-0): direct pressure mapping in a submerged journal bearing shows
  measured sub-cavity pressures of only **−13.6 kPa gauge at 1840 rpm and at
  most −34 kPa at 3000 rpm** — tens of kPa of tension, "far above the oil's
  vapor pressure (~0.1 mmHg at 38 °C)". Mass-spectrometer sampling of the cavity
  found **air and no oil vapors** — the mechanism was "release of dissolved air
  as the oil pressure fell below the saturation pressure." The cavity pressure
  is also **not constant**: it "gradually rises until it smoothly joins the full
  film pressure" over the last ~45–50° of the cavity.
- **Braun & Hendricks 1984** (ASLE Transactions 27(1):1–14): in a flooded
  bearing (L/D = 0.75, ε = 0.4) the minimum film pressure was
  **operating-condition-dependent**: no sub-atmospheric pressure at all at
  2000 rpm (0.13–0.137 MPa), falling to ~0.08 MPa at 4000 rpm and ~0.034 MPa at
  6000 rpm. Cavity contents were "more like air than hydrocarbon vapors"
  (spectrographic analysis); with **CO₂-saturated oil the cavitation zone was
  larger than with air-saturated oil**, "due to the larger solubility in oil of
  CO₂" — i.e. **cavitation extent tracks dissolved-gas solubility**, exactly the
  dependence a fixed `p_cav` cannot represent. Their theoretical frame:
  `P_c = P_v + P_g`; "the more gas in the oil, the closer the cavity pressure
  will be to the atmospheric pressure and the sooner the rupture of oil film
  will occur."
- San Andrés (TAMU Notes 06) summarizes the same experiments: phase change (oil
  vaporization) "requires a source of energy not readily available in actual
  operation" — steady flooded-bearing cavitation is gaseous.
- **Every published refrigerant-oil bearing model puts onset at the mixture
  saturation pressure.** Grando, Priest & Prata (Tribology & Interface Eng.
  Series, 2005, S0167892205800501): "the release of refrigerant gas from the
  lubricant mixture **when saturation pressure is reached** in the divergent
  region … is considered directly, with a two-phase flow thereafter." Int. J.
  Refrigeration 2023 (S0140700723004577): pressure variation "results in the
  variation in solubility of refrigerant in oil, and causes gaseous cavitation."
- **Timescale guidance for when each regime applies.** Sun & Brewe 1992 (ASME
  *J. Tribology* 114, via San Andrés notes): vaporization/condensation is fast
  versus machine periods (> 1 ms) while **gas diffusion is orders of magnitude
  slower** — dynamically formed bubbles are vapor; steadily loaded gas-rich films
  outgas. A refrigerant compressor bearing under quasi-steady load is squarely in
  the gaseous regime.

**Verdict.** For pure-oil cases a constant low `p_cav` is the accepted JFO
idealization (with the caveat — Etsion & Ludwig — that even then enclosed-cavity
pressure is not strictly constant). For the **R290/PZ68 case it is unphysical**:
the dominant void mechanism in the shipped configuration is outgassing at
~10 bar, 30× above the configured rupture threshold. The code computes the
correct supersaturation driver (`c_d − c_eq(p,T)`) and then lets it do nothing
to the void or the pressure (§ A2).

**What to do.**
1. Make gaseous void onset track `p_sat(T, c_d)` (the Henry/solubility
   inversion already in the code) rather than `p_cav`; reserve `p_cav` for true
   vapor rupture. The clean way to get this is not a second threshold but the
   coupled closure of § A2 — when released gas occupies volume in the film
   continuity, pressure caps near `p_sat` automatically (this is the mechanism
   in Grando's and the 2023 models).
2. Log `p_sat(T, c_d)` and the boundary saturation ratio at startup so a config
   like the current one (sump exactly at saturation) is visible (§ B1).
3. Document, per Braun & Hannon's taxonomy, which regime each config is meant to
   represent.

## A2. Two-phase closure: released gas is volumetrically inert — **UNPHYSICAL, no published precedent**

**What pancake does (all re-verified).**
- Mixture density fed to Reynolds is liquid-only: `ρ = ρ_l·θ`
  (`fluid_properties.cpp:216`).
- `reynolds.cpp` contains **zero references** to `alpha_gas` or `free_gas_mass`
  (grep-verified); the Elrod solve is untouched by gas.
- The outgassing exchange writes `dissolved_gas`, `free_gas_mass`, `alpha_gas`,
  and the rate `gas_mass_transfer` as **output fields only**
  (`fluid_properties.cpp:274–282`); the source enters no solved continuity
  equation. Gas transport (`gas_transport.cpp`) advects the scalars but its
  result feeds back only through the Einstein viscosity bump (§ A3).
- `α_g = free_mass/(ρ_g h)` is defined against the **full gap**
  (`fluid_properties.cpp:274–276`) while θ is liquid volume per gap — nothing
  enforces `θ + α_g ≤ 1`. With `gas_alpha_max = 1.0` (`config_r290_pz68.txt:69`,
  allowed by `config.hpp:767–769`) the release cap is an entire gap of gas on
  top of the liquid bookkeeping.

**What the research says.**

- **The Elrod variable already carries an implicit Σα = 1 closure.** San Andrés
  Notes 06: in the universal algorithm "θ is termed the fractional film content
  and **(1−θ) represents the void (gas volume) fraction**." A second,
  independently transported void variable breaks the closure the formulation is
  built on.
- **Homogeneous-mixture coupling is the published standard for exactly this
  problem.** Grando, Priest & Prata 2006 (*Tribology Letters*, DOI
  10.1007/s11249-006-9027-6): "the liquid–gas mixture is replaced … by a
  monophasic pseudofluid, whose density and dynamic viscosity are given by
  `ρ̄ = φρ_g + (1−φ)ρ_l` and `μ̄ = χμ_g + (1−χ)μ_l`" — and these mixture
  properties **enter the Reynolds equation directly**. The companion study
  (S0167892205800501) emphasizes the payoff: "the methodology presented
  **eliminates the use of intermediate boundary conditions** in determining
  cavitation."
- Int. J. Refrigeration 2023 (S0140700723004577): "a numerical model is
  developed to **simultaneously solve** the gaseous cavitation and two-phase
  lubrication problems." Quantified consequence of the coupling: load capacity
  "**reduced by about 3% at a lower eccentricity of 0.3, but … greatly enhanced
  by about 21% at a higher eccentricity of 0.9** due to the gaseous cavitation."
  That is the error band pancake's inert-gas treatment leaves on the table.
- ASME *J. Tribology* 147(2):024101 (2025, DOI 10.1115/1.4066414): the gas-oil
  refrigerant bearing problem is solved through a "**solubility-based
  compressible Reynolds equation**" — dissolved/free gas in the pressure
  equation itself, to the point that even the *dynamic coefficients* are derived
  by perturbing it.
- **No fetched source — across 21 — transports a free-gas fraction that does not
  feed back into film continuity.** The closest concept, JFO's own (1−θ), is by
  construction inside the mass balance.

**A2.3 — Where the review thread needs correcting: the kinetics are fine.**
`REVIEW_VERDICT.md` treats the finite-rate exchange as part of the problem
surface; the literature says otherwise:
- Hao & Gu 2014 (*Tribology International* 78:14–26, DOI
  10.1016/j.triboint.2014.04.028): a gaseous cavitation model "considering
  probable **non-equilibrium dissolution effects**… the rate of gas absorption or
  release process is modeled to be related to local variables" — direct
  peer-reviewed precedent for `k_dg(c_d − c_eq)`. (Caution for citation use:
  the claim that the non-equilibrium model "performs in better accordance with
  experiment" belongs to the group's later tilting-pad paper — Ding et al.,
  *Science Progress* 2021, DOI 10.1177/00368504211029431 — not the 2014 paper.)
- Grando et al. 2006 bracket kinetics with equilibrium vs. never-redissolve
  limits and find the choice **materially changes the answer**: "Under
  non-equilibrium conditions, the behaviour of the bearing is significantly
  altered. Gas is present along all the bearing … the pressure profile spreads
  across a wider extent."
- Etsion & Ludwig 1982 observed the physical loop the kinetics represent:
  pressures **above supply** (p/p_s > 1) deep in the convergent zone, explained
  by oil "over-saturated in the sub-cavity pressure loop … air which is
  subsequently pressurized … forces the air back into solution at the other end
  of the cavity" — finite-rate release *and* resorption, which pancake's
  exchange term already models (`fluid_properties.cpp:250–268`).

**Verdict.** The transported-scalar + finite-rate-exchange machinery is
publishable; its **decoupling from the pressure equation is not**. As shipped,
the "mass-conserving" label is true for liquid only, `θ + α_g` is unconstrained,
and the single physical mechanism of gaseous cavitation — gas displacing liquid
and capping pressure near saturation — cannot occur.

**What to do (one of, explicitly chosen and documented).**
1. **Homogeneous closure (Grando route, recommended):** one void fraction
   combining vapor and released gas; mixture density `ρ̄(p, T, α)` (and § A3
   viscosity) into the Elrod/Reynolds solve. JFO θ and α_g merge; Σα = 1 by
   construction.
2. **Two-field closure:** keep θ and α_g but let α_g consume part of (1−θ)
   (enforce `α_l + α_g ≤ 1`) and add the gas-continuity source
   `∂(ρ_g α_g h)/∂t + ∇·(ρ_g α_g h u) = Δṁ_g` into the coupled pressure/void
   system so outgassing relieves pressure toward `p_sat`.
3. Until either lands: relabel `alpha_gas`/`free_gas_mass`/`gas_mass_transfer`
   in docs and output as **non-conservative diagnostics**, and gate any R290
   load/extent claims accordingly (the 2023 paper's −3%/+21% numbers bound the
   suppressed physics).

## A3. Mixture viscosity: Einstein dilute law to α_g → 1, and the sign of the dominant effect — **UNPHYSICAL at high void; dissolved-gas thinning is the first-order effect**

**What pancake does.** Free-gas correction
`μ = μ_l(1 + 2.5 α_g)` with α_g clamped only at 0.99
(`fluid_properties.cpp:205–206`) ⇒ up to ×3.475 *thickening*;
`gas_alpha_max = 1.0` permitted (`config.hpp:767–769`). Dissolved-gas thinning
*is* implemented (`exp(a_c·c_d)`, `fluid_properties.cpp:154–165` — Andrade T ×
gas dilution × Barus p), which is the right structure.

**What the research says.**
- Measured R-290/synthetic-oil viscosity (Int. J. Refrigeration 2023,
  S0140700723002451; double-capillary apparatus, PAG 68 / **POE 75** / PVE 68,
  303–348 K — bracketing pancake's 40 °C calibration point): "the absorbed R-290
  has a great impact on the oil viscosity in the low temperature region" —
  dissolved propane **reduces** viscosity, more strongly as solubility rises.
  The validated correlation is the **Eyring-MTSM activity model**, "i.e., a
  non-ideal activity-coefficient mixture formulation — not a linear mixing rule
  or dilute-suspension (Einstein-type) correction."
- Grando et al. (S0167892205800501): "a reduction in the load carrying capacity
  occurs for the lubricant mixture **due to its lower viscosity**." Their free-gas
  viscosity is quality-weighted toward μ_g (homogeneous model), which **limits to
  the gas viscosity as gas dominates** — the physically correct asymptote.
- The Einstein 2.5-coefficient correction is a dilute-suspension result; no
  fetched bearing model applies it beyond dilute void, and pancake's own
  PHYSICS.md flags it as first-order only.

**Verdict.** Directionally wrong for α_g ≳ 0.1–0.3 and grossly wrong as
α_g → 1 (a gas-filled gap modeled as 3.5× *thicker* than oil). The dissolved
side of the model is structurally right but calibrated at one point (§ A4).

**What to do.** Replace the Einstein bump with a mixture law that limits to μ_g
(quality-weighted McAdams/Dukler-type as in Grando, or at minimum cap
`gas_alpha_max ≈ 0.5–0.6` and document the validity bound as a runtime guard).
For the liquid solution, prefer the existing `TABLE`/correlation path seeded
from measured R290/POE data (the 2023 dataset covers exactly this pair class);
keep the one-point Henry+exponential as a documented local tangent (§ A4).

## A4. Solubility: one-point linear Henry at 21.75 wt% — **NON-STANDARD at this concentration**

**What pancake does.** `c_eq = H(T)·p`, H calibrated at a single (10 bar, 40 °C)
point, van't Hoff T-factor (`fluid_properties.cpp:76–85`);
`dissolved_gas_initial = 0.2175` (`config_r290_pz68.txt:57`). `BUNSEN` and
`TABLE` alternatives exist (`fluid_properties.cpp:87–97`) but are unused by the
shipped configs.

**What the research says.** Grando et al. characterize the mixture "through
**correlations** for solubility, density and viscosity" (full isotherms, not a
one-point tangent); their reference case runs at **2 bar with w_sat = 7.13%**
HFC-134a/POE — pancake's flagship case (10 bar, 21.75%) is far deeper into the
non-dilute regime where the 2023 R-290 paper required an activity model
(Eyring-MTSM) rather than any linear law. Hao & Gu's baseline is Bunsen
solubility (pressure-dependent), the same family as pancake's unused `BUNSEN`
path.

**What to do.** Treat linear Henry as a documented local tangent around the
calibration point; for the R290/PZ68 case populate `solubility_table` (and the
viscosity/density tables) from measured isotherms; state the validity window in
PHYSICS.md. This is config/data work, not new code — the `TABLE` plumbing
already exists.

## A5. THD fidelity: depth-averaged single-T film, lumped HTCs, θ-blind energy terms — **SIMPLIFICATION below the 1983–84 validated baseline; one internal inconsistency**

**What pancake does (verified).** One depth-averaged film temperature; lumped
wall loss `h_wall(T − T_wall_blend)` with fixed wall temperatures
(`energy.cpp:124–140`); no solid conduction; energy capacity `ρc_p·h` with full
gap and **no θ** (`energy.cpp:301`); enthalpy fluxes with full `h_face`,
`ρc_p_face` (`energy.cpp:316–321, 327–339`); Couette dissipation and friction
shear `μωR/h` un-gated by θ (`reynolds.cpp:519,543`; only the Poiseuille part is
g(θ)-gated, `reynolds.cpp:516`). Meanwhile the gas transport *does* θ-weight its
capacity (`gas_transport.cpp:31,193`) — the inconsistency flagged as N5.

**What the research says.**
- **Ferron, Frêne & Boncompain 1983** (ASME *J. Lub. Technology* 105(3):422,
  DOI 10.1115/1.3253647) — the canonical THD benchmark — "takes into account
  heat transfer between the film and **both the shaft and the bush**," plus
  cavitation and lubricant recirculation (hot-oil carryover), and ships measured
  wall pressure/temperature distributions, eccentricity, and flow rates across
  speeds and loads. Even with all of that, agreement with experiment is only
  "satisfactory." The experiments also show **differential thermal dilatation
  measurably shifts the operating eccentricity** — a clearance-feedback channel
  pancake does not have.
- **Lund & Tonnesen 1984** (ASME *J. Tribology* 106(2):237–244, DOI
  10.1115/1.3260891): even the classic "**approximate**" THD level retains a
  **fourth-order polynomial cross-film temperature profile** coupled to **sleeve
  conduction** (Part I, 106(2):228–236), and still reports "significant
  discrepancies … overall agreement is satisfactory" against their own test
  data.
- Both papers are exactly the validation targets a THD bearing code is expected
  to hit (§ E).

**Verdict.** A depth-averaged film with lumped HTCs is a legitimate *first*
model but sits **below the fidelity that was already standard-and-validated in
1983–84**; quantitative peak-temperature or thermal-load claims from it are not
defensible. The θ-blind energy capacity/advection/shear in cavitated cells
overstates thermal mass, transported enthalpy, and friction power by up to ~1/θ
in striated zones and is internally inconsistent with the code's own gas solve.

**What to do.** (1) θ-weight energy capacity, enthalpy fluxes, and
cavitated-zone Couette shear (or document the full-film assumption and bound
it). (2) Keep lumped HTCs only with a calibration note; plan the assumed
cross-film profile + 1D conjugate wall model before any temperature validation
(§ E: Ferron, Lund–Tonnesen). (3) Ship the default config `ISOTHERMAL` or label
the energy output diagnostic until § C4's coupling iteration exists (the default
`config.txt` currently pairs `ENERGY_EQUATION` with `CONSTANT` properties — a
passive scalar).

---

# Part II — Boundary conditions & lubricant feed

## B1. The flagship R290 case: doc↔config drift, and a sump sitting exactly at saturation — **DOC-CODE MISMATCH + undeclared modelling choice**

**What pancake does (verified).** `PHYSICS.md:700–701` documents the R290
example at `p_cav = 1e5`, `bc_z = 1e5` ("1 bar sump"); the shipped
`config_r290_pz68.txt` uses `p_cav = 3.0e4` and `bc_z = 1.0e6`
(`config_r290_pz68.txt:31,107,110`). `PHYSICS.md:148` separately documents axial
ends as "Dirichlet p = p_cav … (open bearing ends)" while both shipped configs
are submerged (1.5 MPa default, 1.0 MPa R290). And the 10-bar sump **equals the
model's own saturation pressure** (`0.2175/2.175e-7 = 1.0e6 Pa`), so boundary
oil enters exactly saturated and nearly the whole interior is supersaturated —
all of it inert to pressure per § A2.

**What the research says.** Submerged/flooded ends are a legitimate, classical
configuration — Etsion & Ludwig's canonical experiment ran "both ends … flooded
with oil at equal supply pressure" (13.6–54.4 kPa gauge) and showed supply
pressure directly modulates cavity extent. High-pressure refrigerant ambients
are studied territory (ASME JT 147(2) 2025: attitude angle and K/C shift
materially under pressurized refrigerant vs. atmospheric). But the published
refrigerant-POE reference case (Grando 2006) runs at **2 bar, 7.13% dissolved**
— pancake's 10 bar / 21.75% is a much more aggressive operating point and needs
to be declared, not implied. The 2023 Int. J. Refrig. paper adds a flooded-end
caveat pancake cannot currently represent: sub-ambient film pressure (to
−1 MPa at ε = 0.8) drives **reverse end flow (−0.057 dimensionless) that ingests
ambient gas** — "ambient gas might be sucked into the bearing" — a mechanism
that requires the § A2 coupling to mean anything.

**What to do.** Reconcile `PHYSICS.md` §6/§13.2.1 with the shipped configs (pick
the documented case and regenerate the other); state whether sump = p_sat is
deliberate; log `p_sat`, boundary saturation ratio, and open-vs-submerged intent
at startup; note gas ingestion at ends as out-of-scope until A2 lands.

## B2. Feed model: pressure-pinned flooded inlets only — **SIMPLIFICATION; starved operation has been inside the Elrod framework since 1989**

**What pancake does.** Inlets are penalty-pinned to a supply pressure
(`reynolds.cpp:79–88`, penalty 1e20, `θ_supply = exp((p_supply−p_cav)/β)`), with
flooded initialization (`reynolds.cpp:244–250`); no mass-flux or
inlet-film-fraction option.

**What the research says.** Vijayaraghavan & Keith 1989 (*Wear* 134(2):377–397,
DOI 10.1016/0043-1648(89)90137-3) ran the Elrod algorithm on a grooved bearing
"for **both flooded and starved inlet conditions**," with parametric "degrees of
starvation" and supply pressure, validated "with available experimental and
theoretical results." ASME *J. Tribology* 144(12):121801 (2022) shows the
modern design expectation: groove **length, width, and angular location are
model inputs**, and optimal placement is load-dependent (100–140° to the load
line, L_g/L > 0.6 for light load; 80–115° for heavy load) — i.e., an Elrod-type
solver intended for design must parameterize feed geometry and admit starvation,
not hard-wire a flooded pressure patch.

**What to do.** Add a prescribed-`θ_inlet` (or feed mass-flux) inlet mode beside
the pressure inlet; document the current model as flooded-only until then.
(Groove geometry is already parameterized in `theta_in`/`z`/`size` terms —
the missing piece is the starved supply condition, not the geometry.)

## B3. Constant cavity pressure as a boundary idealization — **ACCEPTABLE, but say so**

Etsion & Ludwig's measured cavity pressure varies over the cavity's last
~45–50°, "contrary to common assumptions of constant cavitation pressure …
incorrect in the case of enclosed cavitations." All Elrod-type codes share this
idealization; pancake should state it as a known model limit in PHYSICS.md
rather than silently inheriting it — particularly because the submerged-end
configs are exactly the "enclosed cavity" case the experiments criticize.

## B4. Implementation-level axial-boundary flux findings (line-level review, 2026-06-11) — **one live inconsistency, three latent**

A dedicated cross-solver review of how Reynolds/Elrod, energy, and gas
transport each assemble and gate their axial-boundary fluxes (current working
tree, after the WP-3/WP-4/WP-5 changes). What **checked out**: the FVM contract
is coherent (laplacian skips z-boundary faces, `FVM::divergence` zeroes
z-boundary advective fluxes, every boundary flux is added explicitly once — no
double counting); in **full-film** boundary cells the Elrod Dirichlet link
`Γ_base·Δθ/(Δz/2)` is first-order equivalent to the physical Poiseuille flux
`ρh³Δp/(12μ·Δz/2)` via the EOS chain rule `β·dθ ≈ θ·dp`, matching the
energy/gas boundary mass flux; energy and gas use the same g-gated
film-averaged velocity for interior and boundary faces, with their
non-conservative corrections including the boundary fluxes consistently
([energy.cpp:264](src/energy.cpp:264), [gas_transport.cpp:108](src/gas_transport.cpp:108));
the dissolved/free-gas reservoir composition contracts and supply-cell
overwrites are coherent and clamp-mass-accounted. The findings:

**BC-3 (live; medium): the Dirichlet θ-link is unphysical in *cavitated*
boundary cells, and its mass flux carries no enthalpy or composition.** The
boundary link uses ungated `Γ_base` ([reynolds.cpp:380–406](src/reynolds.cpp:380),
comment "flooded-bearing condition"). In a cavitated boundary cell the EOS
`p = p_cav + β ln θ` does not hold (p ≡ p_cav), so the physical reformation
inflow is `ρ₀h³(p_bc − p_cav)/(12μ·Δz/2)` — but the assembled link delivers
`ρ₀βh³(θ_bc − θ_P)/(12μ·Δz/2)`, an overestimate by ≈ `β(1−θ_P)/(p_bc−p_cav)`
(~10²–10³× for θ_P ≈ 0.9 at a 10-bar sump): the boundary row re-floods
unphysically fast. Simultaneously the energy and gas boundary fluxes are
g-gated to **zero** in the same cells ([energy.cpp:149](src/energy.cpp:149),
[gas_transport.cpp:44](src/gas_transport.cpp:44)), so the re-flooding liquid
arrives with no enthalpy and silently adopts the cell's dissolved-gas content
instead of the reservoir's. Live only when the cavitation zone reaches the
axial ends (shipped 10–15 bar submerged ends usually keep boundary rows
full-film; short-bearing / high-ε / low-sump cases hit it).

**BC-1 (latent; medium): Elrod `INLET_OUTLET` outflow is a liquid-mass wall.**
In the outflow state no boundary term is assembled (zero-gradient θ ⇒ zero
flux, [reynolds.cpp:388–414](src/reynolds.cpp:388)) — liquid cannot leave —
while energy and gas compute outflow through the same face from the pressure
gradient. Enthalpy and gas exit a boundary that mass cannot cross.

**BC-2 (latent; medium): stale, cross-solver-inconsistent direction switching.**
The Elrod inflow check reads `fields["pressure"]`, which is recovered only
*after* the outer loop ([reynolds.cpp:482](src/reynolds.cpp:482)) — i.e. the
previous timestep's pressure, frozen across outer iterations — while energy/gas
switch on current-step pressure. The three solvers can disagree about flow
direction at the same face in the same step; chatter risk near `p ≈ bc_val`.

**BC-4 (latent; low): `INLET_OUTLET` means different physics per cavitation
model.** The Gumbel solver treats it as unconditional Dirichlet
([reynolds.cpp:159–168](src/reynolds.cpp:159)); Elrod gates it inflow-only.

**BC-5 (consistent; document): `NEUMANN` ends are sealed/symmetry planes** in
all three solvers (`uses_pressure_boundary` false everywhere) — correct as a
z-symmetry model, wrong if read as "open end"; document the intent.

All shipped configs use `DIRICHLET`, so BC-1/2/4 are opt-in paths. **Fix:**
one shared boundary mass-flux function consumed by all three solvers plus the
diagnostics balance, with the cavitated-cell reformation flux driven by
`(p_bc − p_cav)` — see `AUDIT_PLAN.md` WP-11.

---

# Part III — CFD methods & numerics

## C1. Convection discretization: first-order upwind everywhere; TVD machinery dead — **PARTLY STANDARD, PARTLY NON-STANDARD, plus a documentation bug**

**What pancake does (verified).** Every production `FVM::divergence` call is
`UPWIND` — Elrod Couette (`reynolds.cpp:259`), energy (`energy.cpp:344`), both
gas scalars (`gas_transport.cpp:241,260`). The van Leer/minmod limiters and
deferred correction in `fvm.cpp:13–27,118–162` are dead code; `PHYSICS.md §5.5`
advertises them. The Gumbel path central-differences the same Couette term the
Elrod path upwinds — the documented 1.2% A.1↔Gumbel residual is a
scheme-mismatch artifact, not physics.

**What the research says — a nuance the review thread missed.**
- Upwinding the **cavitated-zone** θ-transport is *standard*: even the
  state-of-the-art FBNS solver "in the cavitation region … a **first-order
  upwind formula was employed**" (Woloszynski, Podsiadlo & Stachowiak 2015,
  *Tribology Letters*, DOI 10.1007/s11249-015-0487-4).
- What is *not* standard is upwinding **everywhere**: the canonical
  Vijayaraghavan & Keith 1989 scheme (STLE *Tribology Transactions*
  32(2):225–233, DOI 10.1080/10402008908981882) "makes use of **type
  differencing** … to automatically switch from central to upwind differences,
  and vice versa, **at cavitation boundaries**" — central (elliptic) in the full
  film, upwind (hyperbolic) only in/across the cavitated zone (San Andrés Notes
  06 describes the same zone-dependence).
- First-order upwind's false diffusion (~|u|Δx/2) is the textbook objection
  (Patankar; Versteeg & Malalasekera) and lands precisely on the rupture /
  reformation / thermal / gas fronts — the outputs the project wants to
  validate, on a 60×20 R290 grid.

**What to do.** Either implement V&K-style type differencing (central in full
film) or wire the existing TVD into θ/T/c solves; in all cases delete-or-honor
the PHYSICS.md §5.5 advertisement, make the Gumbel and Elrod paths use the same
Couette scheme so the regression limit is clean, and publish a 1×/2×/4×
grid-convergence study for cavitation extent, h_min, and load (§ E).

## C2. The β = 1e9 conditioning problem — **KNOWN, DOCUMENTED, AND SOLVED IN THE LITERATURE**

**What pancake does.** `p = p_cav + β ln θ` with β = 1e9 Pa (both configs);
supply band maps to θ ∈ [1, ~1.001]; PETSc noise O(1e−5) in θ ⇒ O(10 kPa)
pressure noise; `snap_tol = 5e-7` flag-chatter hack with its own admission
comment (`reynolds.cpp:328–337`).

**What the research says.**
- The dependence is documented: "The solution of cavitation problem is shown to
  **strongly depend on the specific values chosen for the lubricant bulk
  modulus**," and traditional Elrod schemes yield "a stiff discrete numerical
  system" whose **stability is dominated by β** ("the impact of static load and
  mesh size is negligible on numerical stability compared to dominant
  significance of varying bulk modulus") — ASME *J. Tribology* 139(3):031703
  (2017), which proposes a finite-volume modification stable "for any chosen
  value" of β.
- The artificial-softening tradeoff is quantified there: reduction is "tolerable
  with maximum two order of magnitudes reduction of β to avoid large errors of
  more than **3–40%** in calculated results."
- Practitioner reality (San Andrés Notes 06): "analysts use an artificial low
  value of the bulk modulus … **1/100 to 1/10 of the actual physical
  magnitude**" (~2.41 GPa for mineral oil), and "the relevant literature is yet
  to report the details on a robust numerical procedure for the universal
  lubrication model" (as of those notes). Pancake's β = 1 GPa ≈ 0.4× physical
  is *praiseworthy physically* and *maximally stiff numerically* — the exact
  corner the 2017 paper targets.
- The modern way out removes β entirely: complementarity. Giacopini, Fowell,
  Dini & Strozzi 2010 (ASME *J. Tribology* 132(4):041702, DOI 10.1115/1.4002215)
  formulate mass-conserving cavitation as an LCP ("p − p_cav ≥ 0, θ ≥ 0,
  (p − p_cav)·θ = 0") that handles reformation; Woloszynski et al. 2015 (FBNS)
  reformulate Elrod–Adams/JFO as an unconstrained smooth system solved by
  Newton, with "roughly **two orders of magnitude reduction in computational
  time**" versus established algorithms (incl. Ausas 2009 and Bertocchi 2013)
  and **mesh-independent iteration counts** (900 → 360,000 DOF). Ausas, Jai &
  Buscaglia 2009 (ASME *J. Tribology* 131(3):031702, DOI 10.1115/1.3142903)
  enforce p–θ complementarity with "a simple but effective relaxation scheme" —
  no β at all.

**What to do.** Near-term: report a pressure-space residual, tie KSP `rtol` to
`Δp_target/β`, keep the snap as a stopgap with its limits documented. Real fix:
migrate the cavitation kernel to a complementarity formulation (FBNS is the
efficiency reference; Giacopini 2010 the formulation reference; Ausas 2009 the
transient reference with downloadable source). This one change retires § C2 and
most of § C3.

## C3. Switch-function iteration fragility — **KNOWN failure mode of the classic algorithm**

**What pancake does.** Outer flag loop with flooded re-initialization on
iteration 0 (`reynolds.cpp:237–253`), θ-clamp + snap, convergence on flag count
and Δθ (`reynolds.cpp:339–344`).

**What the research says.** Pre-FBNS, "there is no approach that allows for an
efficient and stable search of the pressure distribution and the cavity
fraction" — the benchmark Elrod–Adams implementation needed **empirically tuned
relaxation (ω_p = 1.8, ω_θ = 0.5)** (Woloszynski 2015). San Andrés: line solvers
"produce numerical instabilities at the nodes where the cavitation zone starts
and ends," forcing point-wise under-relaxation and slow convergence (his 1D demo
took 600 iterations). The 2017 ASME paper exists precisely because "convergence
speed … considerably faster than the widely used techniques" was worth a paper.

**Verdict & what to do.** Pancake's heuristics (flooded init, snap, flag-count
convergence) are recognizable members of this fragile family, neither better nor
worse than the literature's baseline — but the field has moved on. Same remedy
as § C2: complementarity/FBNS. Until then, log iteration counts and flag-flips
per step so fragility is visible in runs, and add the relaxation knobs the
literature found necessary rather than relying on snap alone.

## C4. Coupling architecture: one segregated sweep is not an operating point — **NON-STANDARD as "steady state"**

**What pancake does (verified).** `STEADY_STATE` executes the time loop once
(`main.cpp:205–208`); per step the sequence is properties → Reynolds →
properties → gas transport/exchange → properties → forces → energy
(`main.cpp:220–260`), so Reynolds always consumes last-sweep μ, ρ and the energy
solve never feeds back within the step. No outer iteration, no coupled residual.
Additionally `STEADY_STATE` zeroes `gas_dt` (`main.cpp:243`), and
`update_gas_state` requires `dt > 0` (`fluid_properties.cpp:222–223`) — **a
steady GAS_CAVITATION_MIXTURE run silently freezes the gas model** at
`dissolved_gas_initial` with α_g ≡ 0 and no warning from `validate()`.

**What the research says.** Published coupled solvers iterate within the step:
Ausas 2009 discretizes the journal's equations of motion with Newmark "and the
dynamical variables are **updated within the same relaxation process**" as the
Reynolds/cavitation solve. The THD benchmarks (§ A5) are converged
multiphysics solutions by construction; the refrigerant models (§ A2) solve gas
and lubrication "simultaneously." No fetched source reports a one-sweep
segregated state as a steady operating point.

**What to do.** Add an under-relaxed outer (Picard) iteration over
{Reynolds, properties, gas, energy} with a coupled residual and convergence
log; make `STEADY_STATE` iterate that loop to tolerance. Error (or hard-warn)
on `GAS_CAVITATION_MIXTURE ∧ STEADY_STATE` until the frozen-gas behavior is
replaced by an equilibrium gas state. Document that one transient step ≠
operating point.

## C5. Dead and dormant numerics switches — **CORRECTNESS-OF-EXPECTATIONS bugs**

Verified facts:
- `pressure_time_method` / `temperature_time_method` are parsed, stored,
  GUI-exposed, and written back (`config.hpp:310–312,467–469,901–903`;
  `gui_win32.cpp:1081–1193,1376–1378`) but **read by no solver**; both shipped
  configs select `EULER_EXPLICIT`, silently ignored — the field solves hardwire
  backward Euler (`FVM::ddt_weighted`).
- The Elrod capacity/squeeze terms assemble **only when**
  `MOTION_MODEL == MOVING_BEARING` and the run is transient
  (`reynolds.cpp:261–280`): a **fixed-bearing transient is a quasi-steady θ
  sequence** with no `∂(ρθh)/∂t` at all — e.g. an `omega_at_time` speed ramp on
  a fixed bearing responds instantaneously. (New finding of this audit; the
  review thread's S4 praise implicitly assumed the term was always on.)
- Energy `ddt` is correctly gated on transient (`energy.cpp:345–347`), making
  the Reynolds-side gap more surprising, not less.

**What to do.** Delete the dead selectors (or implement them); assemble the θ
capacity term for all transient runs, not just moving-bearing ones; document
backward-Euler-only status in PHYSICS.md and the GUI.

## C6. FSI time coupling and motion integrator — **SIMPLIFICATION with a known sharper alternative**

**What pancake does.** Fluid force frozen per step and fed to the motion
integrator (`journal_motion.cpp:131–145`); with the shipped `k = c = 0` the
selected RK4 degenerates to constant-acceleration; no sub-iteration.

**What the research says.** Ausas 2009 (above): motion variables updated within
the fluid relaxation loop — an implicit/iterative FSI coupling in the reference
transient cavitation code. Squeeze-film FSI with explicit force lag is
`dt`-fragile in stiff (small-clearance) regimes; pancake's PLAN already
anticipates a Picard loop.

**What to do.** Fold bearing-motion update into the § C4 outer iteration (or at
minimum document a `dt` stability guideline); see § D for the coefficients the
dynamics work actually needs.

---

# Part IV — Rotor-bearing dynamics deliverables

## D1. No linearized K/C extraction — **GAP vs. standard practice, including for this exact application**

**What the research says.**
- Lund 1987 (ASME *J. Tribology* 109(1):37–41, DOI 10.1115/1.3261324) is the
  canonical formalization: stiffness/damping coefficients via small
  perturbation, used for "unbalance response, stability" — whirl margins are
  conventionally judged from K_ij/C_ij, with Lund himself delimiting where
  linearization holds (small motion about equilibrium) and where nonlinear
  transient orbits are needed. The two are complementary; pancake currently
  offers only the latter.
- ASME *J. Tribology* 147(2):024101 (2025) does it for **gas-oil refrigerant
  bearings specifically**, perturbing the solubility-based compressible Reynolds
  equation, and finds two-phase operation "significantly" raises direct
  stiffness and coupled damping at light load, with attitude angle elevated
  under pressurized refrigerant — i.e., the K/C numbers pancake cannot produce
  are exactly where the refrigerant physics shows up for rotordynamics.

**What to do.** Add a small-perturbation post-process at a converged operating
point (requires § C4): perturb x, y, ẋ, ẏ; difference the fluid force; report
K_ij, C_ij and a whirl-frequency-ratio/stability margin. Validate the
perturbation against the nonlinear transient for one case (Lund's
applicability check).

---

# Part V — Validation: what the literature treats as the price of admission

**What pancake has.** Unit/regression tests (sign conventions, EOS round-trips,
ghost exchange, decoupled-limit equivalences). No quantitative comparison to any
published curve or experiment; no conservation monitor (`PLAN.md`'s O(1e−12)
mass goal is asserted, never logged).

**Benchmarks the fetched literature supplies, in suggested order:**

1. **Ausas et al. 2009** (ASME JT 131(3):031702): oscillatory squeeze flow
   **with exact analytical solution** + dynamically loaded journal bearing, and
   the article "is accompanied by the **ready-to-compile source code**" —
   code-to-code and code-to-exact in one shot, directly exercising the transient
   Elrod path.
2. **Vijayaraghavan & Keith 1989** (STLE Trib. Trans. 32(2):225–233): slider +
   journal cases vs. Elrod's algorithm (code-to-code; note their own caveat that
   the 1989 paper's evaluation is *not* experimental — the grooved experimental
   comparison is the companion work). **V&K 1989, Wear 134(2):377–397**: grooved
   (misaligned) bearing **flooded and starved**, "compared with available
   experimental and theoretical results."
3. **Giacopini et al. 2010** (ASME JT 132(4):041702): analytical/semi-analytical
   suite (textured bearings, squeeze-film dampers) for the mass-conserving
   kernel; also the natural target if the § C2 complementarity migration
   happens; FBNS (Trib. Lett. 2015) timing/iteration data as the efficiency
   yardstick.
4. **THD:** Ferron 1983 (measured wall T/p distributions, eccentricity, flow)
   and Lund–Tonnesen 1984 — gate any quantitative temperature claim on these
   (§ A5 upgrades first).
5. **Refrigerant two-phase:** Grando 2006's verification pattern — recover the
   classical Reynolds-BC solution in the equilibrium limit, then bracket
   equilibrium vs. non-equilibrium kinetics; trend-compare against the 2023
   Int. J. Refrig. case (−3%/+21% load shift, sub-ambient end pressures).
6. **Cavity-physics sanity:** Etsion & Ludwig / Braun & Hendricks regime checks
   — near-ambient gaseous cavity pressures at the measured operating points, and
   cavitation-extent growth with dissolved-gas solubility.
7. **Conservation monitor (prerequisite to all of it):** per-step global liquid
   and total-propane balances (domain integral + boundary fluxes + exchange),
   logged, with the θ-clamp/snap and gas clamps accounted — the review thread's
   N9, which the operator-split + clamps currently make invisible.
8. **Reported metrics:** include total power loss and h_min alongside load and
   extent — the standard design outputs per ASME JT 144(12):121801 (2022).

---

# Part VI — Priority order (merged with REVIEW_VERDICT, updated by the literature)

| # | Action | Closes | Literature anchor |
|---|--------|--------|-------------------|
| 1 | Couple released gas into film continuity (homogeneous ρ̄/μ̄ **or** α_g-consumes-(1−θ) + gas continuity); onset at `p_sat(T,c_d)`; keep finite-rate kinetics | A1, A2, A3(sign), B1(physics) | Grando 2006; IJR 2023 (−3%/+21%); ASME JT 147(2); Hao & Gu 2014 |
| 2 | Outer Picard iteration + coupled residual; steady = converged loop; error on steady+gas | C4, A5(coupling), C6 | Ausas 2009 (in-loop updates); all THD benchmarks |
| 3 | Numerics honesty: type-differencing or TVD for θ/T/c; same Couette scheme both solvers; delete/implement time-method selectors; assemble θ-capacity for all transients; grid-convergence study | C1, C5 | V&K 1989; FBNS 2015 (upwind-in-cavity OK); Patankar |
| 4 | Conservation + regime logging: global mass residuals, `p_sat` & saturation ratio at startup, Re guard, iteration/flag-flip counts | E.7, B1, C3 | San Andrés notes; N9/M3 |
| 5 | Property data: solubility/viscosity tables from measured R290/POE isotherms; mixture-μ law with gas-limit asymptote; cap α_g | A3, A4 | IJR 2023 (Eyring-MTSM); Grando μ̄ |
| 6 | Docs reconciliation: PHYSICS.md §6/§13.2.1 vs shipped configs; declare submerged-end + sump-at-saturation intent; constant-p_cav caveat | B1, B3 | Etsion & Ludwig 1982 |
| 7 | Validation campaign in § E order (Ausas first — exact solution + reference code) | E | Ausas 2009; V&K; Giacopini; Ferron; Lund–Tonnesen |
| 8 | K/C small-perturbation post-process + whirl margin (after #2) | D1 | Lund 1987; ASME JT 147(2) |
| 9 | Starved-inlet option (`θ_inlet` / mass flux) | B2 | V&K Wear 1989; ASME JT 144(12) 2022 |
| 10 | Strategic: complementarity/FBNS cavitation kernel (retires β-stiffness, switch-iteration fragility; ~100× speedup reported) | C2, C3 | Giacopini 2010; Woloszynski 2015; ASME JT 139(3) 2017 |

---

## Appendix — Source list (fetched & quoted in this audit)

1. Braun & Hannon (2010), *Cavitation formation and modelling for fluid film
   bearings: a review*, Proc. IMechE Part J, DOI 10.1243/13506501JET772.
2. Etsion & Ludwig (1982), *Observation of pressure variation in the cavitation
   region of submerged journal bearings*, ASME J. Lubr. Technol. 104(2):157–163.
3. Braun & Hendricks (1984), *An experimental investigation of the
   vaporous/gaseous cavity characteristics of an eccentric journal bearing*,
   ASLE Trans. 27(1):1–14, DOI 10.1080/05698198408981539.
4. Vijayaraghavan & Keith (1989), *Development and evaluation of a cavitation
   algorithm*, STLE Tribology Trans. 32(2):225–233, DOI 10.1080/10402008908981882.
5. Vijayaraghavan & Keith (1989), *Effect of cavitation on the performance of a
   grooved misaligned journal bearing*, Wear 134(2):377–397,
   DOI 10.1016/0043-1648(89)90137-3.
6. Ausas, Jai & Buscaglia (2009), *A mass-conserving algorithm for dynamical
   lubrication problems with cavitation*, ASME J. Tribology 131(3):031702,
   DOI 10.1115/1.3142903.
7. Giacopini, Fowell, Dini & Strozzi (2010), *A mass-conserving complementarity
   formulation to study lubricant films in the presence of cavitation*, ASME
   J. Tribology 132(4):041702, DOI 10.1115/1.4002215.
8. Woloszynski, Podsiadlo & Stachowiak (2015), *Efficient solution to the
   cavitation problem in hydrodynamic lubrication* (FBNS), Tribology Letters,
   DOI 10.1007/s11249-015-0487-4.
9. (2017), *A robust modification to the universal cavitation algorithm*, ASME
   J. Tribology 139(3):031703 (β-sensitivity; 3–40% softening errors;
   eigenvalue stability analysis).
10. (2022), *Implementation of the Elrod algorithm in practical journal bearing
    groove design*, ASME J. Tribology 144(12):121801.
11. San Andrés, *Modern Lubrication Theory, Notes 06: Liquid cavitation in
    bearings*, TAMU rotorlab (β practice; type differencing; Sun & Brewe 1992
    timescales; Braun & Hendricks summary).
12. Grando, Priest & Prata (2006), *A two-phase flow approach to cavitation
    modelling in journal bearings*, Tribology Letters,
    DOI 10.1007/s11249-006-9027-6; and companion chapter, Tribology & Interface
    Engineering Series (S0167892205800501).
13. Hao & Gu (2014), *Non-equilibrium gaseous cavitation model*, Tribology
    International 78:14–26, DOI 10.1016/j.triboint.2014.04.028 (see also Ding et
    al. 2021, Science Progress, DOI 10.1177/00368504211029431 for the
    experimental-accordance claim).
14. (2023), *Two-phase lubrication characteristics of journal bearing in
    refrigerant–oil system under high-pressure environment considering gaseous
    cavitation*, Int. J. Refrigeration, S0140700723004577.
15. (2023), *Viscosity of saturated R-290/synthetic-oil mixtures (PAG 68, POE
    75, PVE 68)*, Int. J. Refrigeration, S0140700723002451 (Eyring-MTSM).
16. (2025), *Static characteristics and dynamic coefficients of gas–oil-
    lubricated journal bearings in high-pressure refrigerant environment*, ASME
    J. Tribology 147(2):024101, DOI 10.1115/1.4066414.
17. Ferron, Frêne & Boncompain (1983), *A study of the thermohydrodynamic
    performance of a plain journal bearing — comparison between theory and
    experiments*, ASME J. Lubr. Technol. 105(3):422.
18. Lund & Tonnesen (1984), *An approximate analysis of the temperature
    conditions in a journal bearing, Part II*, ASME J. Tribology 106(2):237–244,
    DOI 10.1115/1.3260891 (Part I: 106(2):228–236).
19. Lund (1987), *Review of the concept of dynamic coefficients for fluid film
    journal bearings*, ASME J. Tribology 109(1):37–41, DOI 10.1115/1.3261324.

Code facts in this audit were verified against the working tree on branch
`feature/journal_bearing` (see REVIEW.md Round 2 and the re-verification pass
for `file:line` anchors).
