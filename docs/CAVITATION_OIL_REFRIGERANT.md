# Cavitation Modeling of Oil–Refrigerant Mixtures in Journal Bearings

**Subject:** PAG-class compressor oil (PZ68S, ISO VG 68) + R290 (propane), modeled in the `pancake` 2D Reynolds solver.
**Scope requested:** physics understanding first; survey all modeling frameworks then recommend; treat properties and cavitation regime across the *whole* pressure–temperature–solubility–viscosity (PTSV) map, not a single operating point.
**Author:** deep-research synthesis (multi-agent web research + adversarial verification, 24/25 claims confirmed by 3-vote check) cross-referenced against the current `pancake` source tree.
**Date:** 2026-06-14

---

## 0. Document status & handoff guide

This document is a **research synthesis meant to be taken over and validated by another agent or by you**. Every substantive section carries a status tag:

| Tag | Meaning |
|-----|---------|
| `[VERIFIED]` | Backed by a primary source that survived 3-vote adversarial verification (see ledger in §11). Vote shown as `3-0`. |
| `[SYNTHESIS]` | My reasoning that *combines* verified facts. The inputs are cited; the combination is mine and should be sanity-checked. |
| `[NEEDS VERIFICATION]` | Asserted from general thermophysics or single-source/abstract evidence; pull the full text or run numbers before relying on it. |
| `[CODE]` | Statement about the current `pancake` implementation, with `file:line`. Verified by reading the tree on 2026-06-14; re-check if the code has moved. |
| `[USER INPUT NEEDED]` | Data or a decision only you can supply. |

**How to take this over:** start at §11 (verification ledger) to see exactly which claims are load-bearing and how strongly. The single most important caveat: **almost all quantitative bearing-cavitation evidence in the literature is for _air-in-oil_, not R290-in-PAG.** Only one source ([R4]) directly treats refrigerant-in-oil journal-bearing cavitation. The air–oil results transfer *by physical analogy* (dissolved gas leaving solution at a solubility-set pressure), but **R290's much higher solubility and very different gas-phase properties mean you must not carry over numerical coefficients uncritically.**

---

## 1. Executive summary — the combined physical picture

`[SYNTHESIS]` For a PAG oil saturated with dissolved R290 in a hydrodynamic journal bearing, the literature converges on a single coherent picture:

1. **It is gaseous (outgassing) cavitation, not vaporous boiling of the oil.** `[VERIFIED 3-0]` A lubricant film has *two* distinct low-pressure thresholds. As pressure falls in the divergent wedge, the **dissolved refrigerant comes out of solution first**, at the local **solution saturation (bubble-point) pressure** — which for an R290/PAG mixture is far above the oil's own vapor pressure. True vaporous cavitation of the base oil requires reaching the oil vapor pressure, which is essentially never attained for a nonvolatile PAG. [R1][R3][R4]

2. **The controlling threshold is a _field_, not a constant.** `[SYNTHESIS]` The outgassing pressure is the **bubble-point of the local solution**, which depends on local temperature *and* local dissolved fraction. Across the PTSV map this is a surface `p_sat(T, c_d)`, not the single `p_cav` used in classical JFO. [R4][R11]

3. **The correct base description is "gas dissolved in liquid oil" governed by solubility (Henry / Bunsen / VLE), _not_ a liquid–liquid mixture** — at least up to saturation. `[VERIFIED 3-0]` Below saturation the film is a single-phase liquid solution whose viscosity and density depend on how much R290 is dissolved. At/after saturation the desorbed R290 forms a **homogeneous two-phase (gas+liquid) film** whose void fraction, mixture viscosity, and mixture density feed back into the pressure solution. [R4]

4. **Both viscosity _and_ density of the two-phase film reshape the pressure profile** — this is not a density-only effect. `[VERIFIED 3-0]` Reported load-capacity changes from including the two-phase mixture: ~3% at low eccentricity rising to ~21% at high eccentricity. [R4]

5. **Equilibrium vs. finite-rate desorption is the main subtlety.** `[VERIFIED 3-0]` Instantaneous-equilibrium solubility (a Henry/Bunsen closure evaluated at local `p,T`) is a defensible *steady-loaded* baseline, but at high speed / dynamic loading it is quantitatively wrong because gas diffusion out of solution is orders of magnitude slower than the rotor period. Non-equilibrium (finite-rate, history-carrying) models match measured minimum film thickness and force balance better. [R1][R2]

6. **The whole thing bolts onto your existing JFO/Elrod–Adams machinery.** `[VERIFIED 3-0]` A compressible Reynolds equation with a barotropic density–pressure closure is a *formal approximation of JFO/Elrod–Adams* [R5]. So `pancake` keeps its θ / switch-function structure; the conceptual change is to **reinterpret `p_cav` as the local bubble-point `p_sat(T,c_d)`** and to make density/viscosity composition- and void-fraction-dependent. **Much of this is already implemented in your working tree** (see §9).

**One-line recommendation:** adopt a **solubility-based compressible (Elrod–Adams-type) Reynolds equation** — equilibrium Bunsen/Henry closure as the baseline, with a finite-rate desorption option for dynamic regimes — ingest the supplier PTSV chart as the solubility/viscosity/density surfaces, and let the desorbed gas drive a homogeneous two-phase density+viscosity coupling. This is exactly the architecture `pancake`'s new `fluid_properties` + `gas_transport` modules are reaching for; §9 validates what is right and flags the gaps.

---

## 2. Physics — cavitation taxonomy and what actually happens in an oil + R290 film

### 2.1 Classical taxonomy for a pure liquid `[VERIFIED 3-0]`

Three named mechanisms are conventionally distinguished:

- **Vaporous cavitation** — the liquid itself flashes to vapor when local pressure drops to its **vapor pressure**. Fast (vaporization/condensation timescale ≪ rotor period). [R1]
- **Gaseous cavitation** — **dissolved gas comes out of solution** when pressure drops to the gas **release/saturation pressure**, which is *higher* than the liquid's vapor pressure. Governed by solubility, and rate-limited by slow diffusion. [R1]
- **Pseudo-cavitation** — entrained gas bubbles already present simply **expand** as pressure drops (no phase change, no mass transfer). [R7]

> [R1] (San Andrés, TAMU Notes 06): "for pressures below ambient (~1 bar), the dissolved gases in a lubricant (air for example) are released… Most mineral oils contain between 8 and 12% in volume of dissolved air (Pinkus, 1990)."

### 2.2 The two-threshold rule — which species cavitates `[VERIFIED 3-0]`

The key fact that resolves the user's "which fluid is cavitating?" question:

> A lubricant has **two release thresholds**. The species with the **higher release pressure** leaves first. Dissolved gas/refrigerant outgasses near or below ambient; the base oil would vaporize only at its (much lower) vapor pressure, which is rarely reached. [R1][R3]

> [R3] (Li, Song, Hao, Gu, ASME J. Tribol. 134:031701): "only when the film pressure decreases, the gaseous cavitation will occur and the dissolved air will be released from the oil," and this is **compatible with the JFO condition**.

`[SYNTHESIS]` For **R290 in PAG**: replace "dissolved air" with "dissolved propane." Propane is far more volatile and far more soluble than air, so the effect is *stronger and the saturation pressure is higher and more composition-sensitive* than the classic air-in-mineral-oil case. The base PAG is essentially nonvolatile, so vaporous cavitation of the oil is not the operative mechanism anywhere on the realistic PTSV map. → **gaseous/outgassing cavitation dominates.**

### 2.3 Liquid–liquid mixture vs. gas-dissolved-in-liquid — when each picture applies `[VERIFIED 3-0 / SYNTHESIS]`

This is the user's central modeling-frame question. The answer is **both, separated by the saturation surface**:

| Region of the PTSV map | Physical state | Correct frame |
|---|---|---|
| `p > p_sat(T,c_d)` (compressed, full film) | R290 fully dissolved in PAG → **single-phase liquid solution** | **Miscible liquid mixture.** Use mixture (solution) viscosity & density as functions of dissolved fraction, T, p. |
| `p ≈ p_sat` (the threshold) | onset of outgassing | **Bubble-point / VLE.** `p_sat` is the local saturation pressure set by solubility. |
| `p < p_sat` (divergent wedge / cavitation) | desorbed R290 gas + liquid solution → **two-phase film** | **Homogeneous gas–liquid two-phase.** Void fraction α; mixture viscosity & density from two-phase closures. |

So it is **not** a pure liquid–liquid problem and **not** a pure gas-in-liquid problem — it is **gas-dissolved-in-liquid that becomes a two-phase film once the solubility limit is crossed.** [R4] implements exactly this transition for a refrigerant–oil journal bearing.

> [R4] (Int. J. Refrigeration 2023, S0140700723004577): "Pressure changes in refrigerant compressors result in variation in solubility of refrigerant in oil, causing gaseous cavitation of the refrigerant-oil mixture in the bearing"; "The release of refrigerant gas from the lubricant mixture **when saturation pressure is reached in the divergent region** of the bearing is considered directly, with a two-phase flow thereafter," via a "**solubility-based compressible Reynolds equation (CRE).**"

### 2.4 The bubble-point as the controlling threshold `[VERIFIED 3-0]`

In JFO/Elrod the cavitation threshold is a constant `p_cav`. For the mixture it must become the **solution bubble-point** — the pressure at which the *first* bubble of R290 appears for the local composition and temperature. A first-principles route exists:

> [R11] (Restrepo, Bent & Michels, AIChE J. 2009): a model "based on applicable theory for the excess Gibbs energy of nonideal solutions," using "Wohl three-suffix equations … to describe nonideal liquid behavior," with code that "predicts **dew and bubble point** conditions while calculating enthalpy and entropy across varying compositions and temperatures" (reported < 4% error for most mixtures except very low T).

`[SYNTHESIS]` Practically, the **supplier PTSV chart already encodes the bubble-point surface**: the solubility curve "mass fraction of R290 dissolved vs (T, p)" *is* the saturation locus. Reading it backwards — "at temperature T, what pressure saturates the current dissolved fraction `c_d`?" — gives `p_sat(T, c_d)`. See §6.

### 2.5 Equilibrium vs. non-equilibrium (the timescale subtlety) `[VERIFIED 3-0]`

The single most important physical caveat for accuracy:

> [R1] (citing Sun & Brewe 1992): "the characteristic time for liquid vaporization … is very small when compared to the typical period of rotating machinery (>1 ms), while … the characteristic time for gas diffusion is **orders of magnitude larger**. Hence a dynamic cavitation bubble must contain fluid vapor." Braun & Hendricks (1984) measured *steady-loaded* flooded cavities formed by gas coming out of solution.

> [R2] (Hao & Gu, Tribology International 78:14-26, 2014): the equilibrium baseline "is a generalization … of an equilibrium dissolution model upon **Bunsen solubility**" (Henry-type, local-pressure-dependent); the non-equilibrium extension makes "the rate of gas absorption or release … related to local variables" and "considers the historical effects on the local dissolution status." Result: "the minimum film thickness of the **non-equilibrium** gaseous cavitation model is closer to the experimental data and its bearing force balance is also better." Ding et al. 2021 (tilting-pad): min-film-thickness error **8.2% (non-eq) vs 17.5% (eq)**; force-unbalance **3.8% vs 9.5%**.

**Important scope limit** `[VERIFIED 3-0]`: [R6] (Groper & Etsion 2001) found dissolved-gas diffusion *inadequate* to explain cavity-zone pressure build-up in a submerged bearing — but that finding is specific to **pressure build-up within an established cavity**, and does **not** mean outgassing is irrelevant to cavitation *onset*, void fraction, or foaming. Do not over-extend it.

**Takeaway for `pancake`** `[SYNTHESIS]`: an **instantaneous-equilibrium** closure (the default) is reasonable for **steady-loaded** operation; for **high-speed or dynamically-loaded** operation you need a **finite-rate transport** equation for the dissolved fraction. Your tree already has both paths (`gas_mass_transfer_rate`, `gas_transport.cpp`) — see §9.

---

## 3. The combined physical model (the synthesis you asked for)

`[SYNTHESIS]` Pulling §2 together into one consistent model. State the film with three local fields beyond the usual `h`:

- `θ` — fractional film content / liquid content (the JFO universal variable). `θ ≥ 1`: full film, liquid compressed; `θ < 1`: cavitated, with **void fraction `α = 1 − θ`**.
- `c_d` — dissolved-R290 mass fraction in the liquid solution (kg R290 / kg liquid).
- `m_g` (or `α`) — desorbed free-gas content.

**Governing logic across the PTSV map:**

```
At each cell, given local T and the local solution composition c_d:
  p_sat = bubble-point pressure of the (oil + c_d·R290) solution at T      [from PTSV chart / VLE]

  if  p ≥ p_sat:   FULL FILM
      single-phase liquid solution
      μ = μ_solution(T, p, c_d)         (mixture viscosity, dissolved gas thins the oil)
      ρ = ρ_solution(T, p, c_d)         (mixture density)
      Reynolds equation is ELLIPTIC; p solved from mass conservation
      compressibility set by liquid bulk modulus β

  if  p < p_sat:   CAVITATED / TWO-PHASE
      pressure pinned at p = p_sat  (the JFO complementarity, generalized)
      R290 desorbs → void fraction α = 1 − θ grows
      μ = μ_2φ(α, μ_solution, μ_gas)    (two-phase mixture viscosity)
      ρ = ρ_2φ(α, ρ_solution, ρ_gas)    (two-phase mixture density)
      Reynolds equation is HYPERBOLIC/parabolic in θ (Couette transport of liquid content)
```

**The generalized JFO complementarity** `[SYNTHESIS, from R4+R5]`:

$$
p \ge p_\mathrm{sat}(T, c_d), \qquad \theta \le 1, \qquad (p - p_\mathrm{sat})\,(1-\theta) = 0
$$

i.e. either full film (`p ≥ p_sat`, `θ = 1`) or cavitated (`p = p_sat`, `θ < 1`). This is **exactly the classical JFO complementarity with `p_cav` promoted to the local field `p_sat(T,c_d)`.** That single substitution is the conceptual heart of the whole extension.

**Regime map across the supplier PTSV map** `[SYNTHESIS]`:

| Local condition | Dominant mechanism | What the model must do |
|---|---|---|
| High p, any T (loaded arc) | none (fully dissolved) | liquid-solution μ, ρ; classic Reynolds |
| p falling toward `p_sat`, steady load | equilibrium outgassing | equilibrium Bunsen/Henry `c_sat(p,T)`; θ-switch at `p_sat` |
| p < `p_sat`, **high speed / dynamic** | **non-equilibrium** outgassing | finite-rate desorption (history term); cavity may be vapor-rich, not gas-saturated |
| pre-entrained bubbles, mild Δp | pseudo-cavitation | bubble expansion (EOS of free gas), no mass transfer |
| p → oil vapor pressure (rare) | vaporous | only if `p_sat` mechanism exhausted; usually negligible for PAG |

---

## 4. Survey of modeling frameworks (then recommendation in §6)

### 4.1 JFO / Elrod–Adams switch-function, mass-conserving `[VERIFIED 3-0]`

The mass-conserving classic. One universal variable `θ`; a switch `g ∈ {0,1}` toggles the Reynolds equation between elliptic (full film) and parabolic/hyperbolic (cavitation):

> [R1]: "In the cavitation zone … θ is termed the **fractional film content** and (1−θ) represents the **void (gas volume) fraction**"; "P = P_cav + g·Λ·ln(θ)"; "A switch function (g=1 or 0) … switches the character of the … (Reynolds) equation from elliptic to parabolic in the full film and cavitation regions."

- **Assumptions:** single fixed cavitation pressure `p_cav`; cavity is a constant-pressure liquid–gas striation; instantaneous switching (equilibrium).
- **Strengths:** mass-conserving (correct film reformation), robust, cheap, well understood; **this is what `pancake` already implements.**
- **Limits for our problem:** `p_cav` is a *constant* — it cannot represent a composition/temperature-dependent bubble-point, nor void-fraction-dependent mixture properties, without extension.

### 4.2 Homogeneous-equilibrium / compressible-mixture (barotropic) `[VERIFIED 3-0]`

Treat the film as a single compressible medium with a density–pressure law (variable bulk modulus). The key theoretical result:

> [R5] (Bayada & Chupin, Tribology International 71 (2014) 38-49): "This new cavitation model is based on a **compressible Reynolds equation in which the density-pressure relation is obtained from a barotropic-isentropic assumption. It can be viewed as an approximation of the Jakobson-Floberg-Olsson/Elrod–Adams cavitation model.**"

- **Assumptions:** local thermodynamic equilibrium; a barotropic `ρ(p)` closure; the JFO saturation pressure and bulk modulus map onto the compressible-model parameters.
- **Strengths:** smooth, single PDE, naturally accommodates **variable density and mixture compressibility**; the cleanest home for a solubility-driven `ρ(p, c_d, T)` closure; recovers JFO in the stiff-`β` limit.
- **Limits (stated by the authors):** equivalence to JFO is **approximate, not exact** — discrepancies under light loading, starvation, strong viscosity variation; the compressible model *permits* sub-saturation (negative) pressure that strict JFO forbids; recovering JFO numerically needs a **high bulk modulus (β up to ~1000 GPa)**.

### 4.3 Bubbly / Rayleigh–Plesset / non-equilibrium outgassing / VOF two-phase `[VERIFIED 3-0]`

Resolve the gas phase explicitly — either as a transported void fraction with finite-rate mass transfer, or (most detailed) full VOF CFD with a non-condensable-gas phase.

> [R7] (Wettmarshausen, Engels et al., *Lubricants* 2025, 13(4):140): "a new approach for modeling cavitation in hydrodynamic bearings by using **CFD with the volume of fluid method and a phase of non-condensable gas** in the lubrication oil," experimentally validated, concluding the dominant mechanism is "pseudo-cavitation (a form of gaseous cavitation)." Also: "cavitation in hydrodynamic journal bearings is still a not completely understood phenomenon … it is unclear which proportions of different cavitation types are present."

> [R2]: finite-rate non-equilibrium outgassing/redissolution carrying local dissolution history (see §2.5).

- **Assumptions:** explicit void-fraction transport; rate constants for desorption/redissolution (calibrated, **not** first-principles); optionally bubble dynamics (needs interfacial tension, nucleation model).
- **Strengths:** most physically faithful; required for dynamic loading and for genuine foam/bubble effects.
- **Limits:** rate coefficients are fitted; full VOF CFD is far heavier than a Reynolds solver; **R290-specific rate data is essentially unavailable** — the validated cases are air-in-oil [R7] or air-in-oil thrust bearings [R2].

### 4.4 Side-by-side comparison `[SYNTHESIS]`

| Axis | JFO / Elrod–Adams (§4.1) | Compressible-mixture / barotropic (§4.2) | Two-phase / non-equilibrium (§4.3) |
|---|---|---|---|
| Mass-conserving | ✔ | ✔ | ✔ |
| Variable bubble-point `p_sat(T,c_d)` | ✗ (constant `p_cav`) | ✔ natural | ✔ |
| Void-fraction-dependent μ, ρ | ✗ (needs bolt-on) | ✔ via `ρ(p)`, μ-closure | ✔ explicit |
| Non-equilibrium desorption | ✗ | ✗ (equilibrium) | ✔ |
| Cost | lowest | low–moderate | high (esp. VOF) |
| Recovers JFO | — | ✔ (stiff β) [R5] | ✔ (fast-rate limit) |
| Best for | steady, classic | **steady mixture across PTSV** | dynamic / high-speed / foam |
| Maps to `pancake` today | **yes (current)** | **yes (1-step extension)** | partial (`gas_transport` exists) |

---

## 5. Property closures across the PTSV map

### 5.1 Solubility → the saturation threshold `[VERIFIED 3-0]`

The solubility surface does double duty: it sets the dissolved fraction in the liquid region **and** defines `p_sat` (outgassing threshold). Closure options, in increasing fidelity:

- **Henry's law:** `c_sat = H(T)·p`, with `H(T) = H_ref·exp[E_H(1/T − 1/T_ref)]`. Simple, good for dilute gas.
- **Bunsen coefficient:** `c_sat = β_Bunsen·(p/p_ref)·(T_ref/T)·ρ_g,ref/ρ_liq`. The form used by the air-in-oil equilibrium baselines [R1][R2].
- **Tabulated PTSV / VLE model:** the supplier chart, or a thermodynamic VLE model [R9][R11]. **Required for R290/PAG**, because:

> [R9] (Yokozeki, Purdue ICEC 1994): "Often, the correlations are unable to represent all data with temperature, pressure, and compositions at the same time," and "General and practical models, based on **thermodynamic theory**, have been developed for refrigerant-oil solubility." → Use the measured PTSV surface or a VLE model, not a one-parameter fit.

`[VERIFIED 3-0]` **R290/PAG-specific magnitude:** [R10] (IIR FRIDOC #22122): "POE and PAG indicated smaller amount of R290 dissolution than mineral oil. Especially, **PAG showed almost 50% of solubility compared with mineral oil**"; corroborated "R290 was least soluble in PAG VG68 (~5%)… Partially miscible PAG is used for R290 applications." So your PZ68S dissolves *less* R290 than a mineral oil would — but enough (single-digit-% scale, condition-dependent) to drive gaseous cavitation. The "~5%" / "half" figures are single-condition; **this is exactly why a full PTSV surface is needed, not a constant.**

### 5.2 Liquid-solution (mixture) viscosity — full film `[VERIFIED 3-0 + REFUTED]`

Dissolved R290 **thins** the oil. Computing this:

- **Preferred: tabulated PTSV viscosity surface** `μ_solution(T, p, c_d)`, or an empirical correlation fitted to it: e.g. `μ = μ_oil·exp[E_μ(1/T − 1/T_ref)]·exp[a_c·c_d]·exp[α_p(p − p_ref)]` (Andrade temperature × dissolved-gas thinning × Barus piezo-viscosity).
- **REFUTED `0-3`:** "a simple composition-based mixing rule suffices for refrigerant-oil mixture viscosity." **Do not** rely on a single Grunberg–Nissan / log-mixing rule alone for R290/PAG — use the measured surface or a VLE-anchored correlation. [R9]

`[NEEDS VERIFICATION]` Log-mixing (`ln μ = Σ x_i ln μ_i`) / Walther–ASTM remain acceptable *interpolants* of the supplier data, but not as standalone predictors.

### 5.3 Two-phase effective viscosity — once gas is present `[VERIFIED 3-0]`

When `α > 0`, the mixture viscosity must interpolate between liquid and gas with the correct limits (`μ = μ_liq` at quality 0, `μ = μ_gas` at quality 1):

> [R8] (Awad & Muzychka, Exp. Thermal & Fluid Sci. 33:106-113, 2008): "Using an analogy between thermal conductivity of porous media and viscosity in two-phase flow, new definitions for two-phase viscosity are proposed," validated against R-12, R-22, R134a, R410A, R717, **argon, and propane** data.

Standard homogeneous closures (all share the correct limits) `[NEEDS VERIFICATION on best choice for R290]`:
- **McAdams (harmonic / quality):** `1/μ = x/μ_g + (1−x)/μ_l`
- **Cicchitti (mass-weighted):** `μ = x·μ_g + (1−x)·μ_l`
- **Dukler / Beattie–Whalley (void-fraction-weighted):** `μ = α·μ_g + (1−α)·μ_l`
- **Einstein/Taylor dilute-suspension:** `μ = μ_l(1 + 2.5α)` (low void only)

Which one best matches **oil + R290** across the void range is an open calibration question (see §10).

> `[VERIFIED 3-0]` Why this matters: [R4] reports "variations of viscosity and density of the gas-liquid mixture influence the pressure profile **significantly**"; air-in-oil work [R2] finds "the air volume fraction … affects the mixture viscosity significantly, eventually influencing the shear stress … and bearing mechanical loss," with **viscosity dominating shear in the cavitated zone.** Viscosity coupling is not optional.

### 5.4 Density and void fraction `[VERIFIED 3-0 / SYNTHESIS]`

- **Liquid-solution density** `ρ_solution(T,p,c_d)`: mass–volume mixing of oil and dissolved-R290 partial volumes, or from the PTSV/VLE data.
- **Two-phase mixture density:** `ρ_2φ = (1−α)ρ_solution + α·ρ_gas`, with **void fraction `α = 1 − θ`** linking directly to the Elrod content variable.
- **Coupling to compressibility/cavitation index** `[VERIFIED 3-0]`: in the Elrod formulation, `θ` carries a *dual* meaning — `θ > 1` is liquid compression via bulk modulus β; `θ < 1` is fractional film content with `(1−θ)` the gas void fraction [R1]. The mixture density therefore enters the conserved quantity `ρh` and the gas EOS sets how `α` responds to pressure in the cavitated zone.

### 5.5 Ingesting the supplier PTSV chart — concrete recipe `[SYNTHESIS]`

The Daniel/Idemitsu-type PTSV chart (viscosity + solubility vs T, p; "100% = pure refrigerant") is **exactly the liquid-region closure surface.** Recommended ingestion:

1. **Solubility surface → `c_sat(p, T)`** and its inverse **`p_sat(T, c_d)`** (the outgassing threshold). This is the most important extraction — it replaces the constant `p_cav`.
2. **Viscosity surface → `μ_solution(T, p, c_d)`** for the full-film region (table interpolation).
3. **Density:** if the chart gives mixture density, table it; otherwise derive `ρ_solution` from mass–volume mixing using R290 liquid-phase partial density `[USER INPUT NEEDED]`.
4. **Beyond saturation (gas region):** the PTSV chart typically stops at the saturation line. For the **desorbed gas phase** you need a **separate R290 gas EOS** (ideal-gas floor, or a real-gas/REFPROP EOS) for `ρ_gas`, `μ_gas` `[USER INPUT NEEDED]`.
5. Feed `ρ_2φ`, `μ_2φ` back into the Reynolds coefficients (§3).

---

## 6. Recommendation

`[SYNTHESIS — primary recommendation]` **Adopt a solubility-based compressible Elrod–Adams Reynolds equation**, structured as follows:

1. **Keep the Elrod–Adams θ / switch-function backbone** (`pancake` already has it). It is mass-conserving and it is a *formal approximation of the very compressible-mixture model you want* [R5], so you are not throwing anything away.
2. **Promote `p_cav` → `p_sat(T, c_d)`**, the local solution bubble-point read from the PTSV solubility surface (§5.1). This is the single most important physics change. The complementarity in §3 is unchanged in form.
3. **Make μ and ρ composition- and void-fraction-dependent** (§5.2–5.4): liquid-solution properties below `p_sat`; two-phase mixture properties above it. Couple both back into the Reynolds coefficients.
4. **Default to an equilibrium (Bunsen/Henry) closure**; provide a **finite-rate desorption** option (transport equation for `c_d` with a release/resorption rate) for high-speed/dynamic regimes (§2.5). [R2]
5. **Ingest the supplier PTSV chart** as the solubility/viscosity/(density) surfaces; supply an **R290 gas-phase EOS** separately for the desorbed phase.

**Why this and not the alternatives:**
- Pure JFO with constant `p_cav` (§4.1) cannot represent the composition/temperature-dependent threshold that *is* the physics here.
- Full VOF two-phase CFD (§4.3) is the most faithful but is far heavier than a Reynolds solver and has **no R290-specific calibration data** — overkill and under-supported for now. Keep it as a validation reference, not the production model.
- The barotropic compressible-mixture view (§4.2) is the *theoretical justification* for steps 1–3, not a separate code path.

**Fidelity ladder** (pick per regime / as validation allows):
- **Tier 0 (baseline):** Elrod–Adams, constant `p_cav`. *(current `pancake` default)*
- **Tier 1 (recommended):** + `p_sat(T,c_d)` from PTSV + equilibrium solubility + two-phase μ,ρ. Steady-loaded.
- **Tier 2:** + finite-rate desorption transport for `c_d`. Dynamic / high-speed.
- **Tier 3 (research/validation):** explicit bubble dynamics or VOF; needs interfacial tension, nucleation, rate constants.

---

## 7. Open physical questions (regime-determining)

These are not yet answerable from the literature and drive which tier you need:

1. **`[USER INPUT NEEDED]` What are R290's desorption/diffusion timescales in PZ68S relative to the rotor period?** This decides equilibrium (Tier 1) vs. non-equilibrium (Tier 2). High speed / dynamic load → almost certainly Tier 2 [R2].
2. **`[USER INPUT NEEDED]` Does the PTSV chart extend past the saturation line, or must an R290 gas-phase EOS be supplied separately?** (See §5.5 step 4.)
3. **`[NEEDS VERIFICATION]` Which two-phase μ-closure (Awad–Muzychka / McAdams / Beattie–Whalley / Dukler / Einstein) best matches oil+R290 across the void range,** and how to reconcile the *fractional-film content* `θ` (volume-based) with the *mass-quality*-based two-phase property models?
4. **`[USER INPUT NEEDED]` Interfacial/surface tension and a nucleation criterion** — only needed if you go to Tier 3 bubble dynamics.

---

## 8. Data you must supply (PZ68S + R290)

| Quantity | Why | Likely source |
|---|---|---|
| Solubility surface `c_sat(p,T)` (the PTSV "S" curves) | sets `c_d` and `p_sat` threshold | supplier PTSV chart |
| Viscosity surface `μ(p,T,c_d)` (the PTSV "V" curves) | full-film mixture viscosity | supplier PTSV chart |
| Mixture/solution density vs (p,T,c_d) | `ρ_solution`; conserved `ρh` | supplier data or mass–volume mixing + R290 partial density |
| R290 gas-phase EOS (`ρ_gas`, `μ_gas` vs p,T) | desorbed-phase density/viscosity beyond saturation | REFPROP / CoolProp / ideal-gas floor |
| Desorption/redissolution rate constants | Tier 2 non-equilibrium | calibration (no R290 literature values) |
| Interfacial tension, nucleation | Tier 3 only | measurement / literature |
| Operating envelope (speed, load, suction/discharge p, T) | picks the tier and the active region of the PTSV map | your machine |

---

## 9. Where this fits in `pancake` — mapping to the existing code

`[CODE — read 2026-06-14]` **Your working tree has already built most of Tier 1–2.** This section validates the choices and flags gaps. Files: `src/equation_of_state.hpp`, `src/reynolds.cpp`, `src/fluid_properties.{hpp,cpp}`, `src/gas_transport.{hpp,cpp}`, `src/config.hpp`, `src/diagnostics.cpp`.

### 9.1 What is already correct and matches the literature

| Literature concept | `pancake` implementation | Status |
|---|---|---|
| Elrod–Adams θ / switch (§4.1) | `EOS::switch_function` (`equation_of_state.hpp:33`), barotropic `p = p_cav + β ln θ` (`:14`), `ρ = ρ_0 θ` (`:22`); applied as Γ multiplier in `reynolds.cpp:324` | ✔ faithful [R1] |
| Void fraction `α = 1−θ` | diagnostics count `θ<1` → `cavitated_fraction` (`diagnostics.cpp:282`) | ✔ [R1] |
| Equilibrium Bunsen/Henry solubility (§5.1) | `fluid_properties.cpp:70-100` (Henry `c_sat=H(T)p`, Bunsen, table) | ✔ [R2] |
| Solution density mass–volume mixing (§5.4) | `fluid_properties.cpp:102-135` (`MASS_VOLUME_MIXING`) | ✔ |
| Solution viscosity Andrade×thinning×Barus (§5.2) | `fluid_properties.cpp:154-165` (`EMPIRICAL_CORRELATION`) | ✔ — but see REFUTED note below |
| Two-phase μ closures (§5.3) | `fluid_properties.cpp:174-200`: Einstein, **Dukler-void**, **McAdams-quality**, Krieger–Dougherty, **Grando linear-quality** | ✔ [R8] menu present |
| Finite-rate desorption (§2.5, Tier 2) | `fluid_properties.cpp:258-322` (`update_gas_state`, `gas_mass_transfer_rate`); `gas_transport.cpp` transports `dissolved_gas` + `free_gas_mass` | ✔ [R2] |
| Refrigerant species + mixture model | `config.hpp` `DissolvedGasSpecies::PROPANE`, `FluidPropertyModel`, `GAS_CAVITATION_MIXTURE` | ✔ |
| PTSV table ingestion (§5.5) | `config.hpp` `solubility_table`, `viscosity_table`, `density_table`; `interpolate_property_table` | ✔ scaffold present |

### 9.2 Gaps / mismatches to address `[CODE + SYNTHESIS]`

1. **`p_cav` is still a single global constant** (`config.hpp:396`, used in `reynolds.cpp` and `EOS::pressure_from_theta`). **This is the #1 recommended change (§6 step 2):** replace the scalar with a **field `p_sat(i,j) = p_sat(T_ij, c_d,ij)`** derived from the solubility surface. The complementarity and switch logic are otherwise unchanged. This is the gap between Tier 0 and Tier 1.
2. **θ normalization uses `ρ_0` (reference density), not `ρ_cav` (liquid density at cavitation pressure).** `equation_of_state.hpp:22` defines `ρ = ρ_0·θ`; the canonical Elrod definition normalizes by `ρ_cav` [R1]. The difference is `O(p/β)` — *numerically tiny* (since `ρ_cav ≈ ρ_0`) but **definitionally distinct**; document it so a reviewer doesn't read it as a bug. `[VERIFIED 3-0 nuance]`
3. **REFUTED-claim guard (§5.2):** the `EMPIRICAL_CORRELATION` and `LOG_MIXING` viscosity paths are fine *as interpolants of measured PTSV data*, but a single mixing rule was **refuted 0-3** as a standalone predictor for refrigerant–oil viscosity. Prefer `TABLE` (PTSV) where data exists; treat the correlation coefficients as fitted-to-data, not first-principles.
4. **Two-phase μ-closure choice is unvalidated for R290** (§7 Q3). The five models in `fluid_properties.cpp:174-200` give materially different μ(α); pick/calibrate against data, and reconcile the **volume-based `θ`** with the **mass-quality-based** closures (`mass_quality = m_g/(m_g+ρ_l θ h)`, `gas_transport`/`fluid_properties`).
5. **Gas-phase EOS** (`gas_pressure_floor`, `mu_gas` in `config.hpp`) is currently a simple ideal-gas floor; for R290 near its critical region a real-gas EOS may be needed (§5.5 step 4). `[NEEDS VERIFICATION]`
6. **Equilibrium-vs-non-equilibrium selection should be tied to the operating regime** (§7 Q1): expose `gas_mass_transfer_rate → ∞` (equilibrium) vs finite (non-equilibrium) as the deliberate Tier-1/Tier-2 switch, and document which regime the bearing runs in.

### 9.3 Suggested validation path `[SYNTHESIS]`

- **Sanity:** with `c_d = 0` and constant properties, Tier-1 must reduce *exactly* to the current Tier-0 Elrod–Adams result (regression test).
- **Stiff-β check:** confirm the barotropic compressible model recovers JFO as `β` grows [R5] (your default `β = 1e9` Pa is modest; [R5] needed up to ~1000 GPa for tight JFO recovery — check sensitivity).
- **Cross-check against [R4]:** reproduce the qualitative trend (load-capacity reduction growing with eccentricity, ~3%→~21%) as a model-credibility test.
- **Mass conservation:** `diagnostics.cpp` already tracks liquid + gas mass balance and residuals — use `gas_residual`/`liquid_residual` to confirm the two-phase coupling conserves mass.

---

## 10. References

Sources that survived 3-vote adversarial verification are marked ✔. Several primary PDFs returned HTTP 403 during automated verification; their abstract quotes were confirmed via multiple independent retrievals but **you should pull the full texts for exact model equations before implementing.**

- **[R1]** San Andrés, L. *Modern Lubrication Theory, Notes 06: Liquid Cavitation Model.* Texas A&M. ✔ — https://rotorlab.tamu.edu/me626/Notes_pdf/Notes06%20Liquid%20cavitation%20model.pdf *(cites Sun & Brewe 1992 J. Tribol. 114:612; Braun & Hendricks 1984; Pinkus 1990)*
- **[R2]** Hao, Song & Gu. *Numerical modeling for gaseous cavitation of oil film and non-equilibrium dissolution effects in thrust bearings.* Tribology International 78:14-26 (2014). ✔ — https://www.sciencedirect.com/science/article/abs/pii/S0301679X1400173X
- **[R3]** Li, Song, Hao, Gu. *Cavitation Mechanism of Oil-Film Bearing and Development of a New Gaseous Cavitation Model Based on Air Solubility.* ASME J. Tribol. 134:031701. ✔ — https://www.researchgate.net/publication/268989053
- **[R4]** *Two-phase lubrication characteristics of journal bearing in refrigerant-oil system under high-pressure environment considering gaseous cavitation.* Int. J. Refrigeration (2023). ✔ **[most directly on-target — refrigerant-in-oil journal bearing]** — https://www.sciencedirect.com/science/article/abs/pii/S0140700723004577
- **[R5]** Bayada & Chupin. *A new model for cavitation … compressible Reynolds equation as an approximation of JFO/Elrod–Adams.* Tribology International 71:38-49 (2014) / ASME J. Tribol. 135:041702. ✔ — https://www.sciencedirect.com/science/article/abs/pii/S0301679X13003575 ; open manuscript https://lmbp.uca.fr/~chupin/FICHIERS-RECHERCHE/Bayada-Chupin13.pdf
- **[R6]** Groper & Etsion. *The Effect of Shear Flow and Dissolved Gas Diffusion on the Cavitation in a Submerged Journal Bearing.* ASME J. Tribol. 123:494 (2001). ✔ — https://asmedigitalcollection.asme.org/tribology/article-abstract/123/3/494
- **[R7]** Wettmarshausen, Engels et al. *Modeling Cavitation in Hydrodynamic Bearings with VOF and a Non-Condensable Gas Phase.* Lubricants 13(4):140 (2025). ✔ — https://doi.org/10.3390/lubricants13040140
- **[R8]** Awad & Muzychka. *Effective property models for homogeneous two-phase flows (two-phase viscosity).* Experimental Thermal and Fluid Science 33:106-113 (2008). ✔ *(validated incl. propane)* — https://www.sciencedirect.com/science/article/abs/pii/S0894177708001027
- **[R9]** Yokozeki. *Solubility and Viscosity of Refrigerant-Oil Mixtures.* Purdue Int. Compressor Eng. Conf. (1994). ✔ — https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2001&context=icec
- **[R10]** *Characteristics of refrigerant-oil mixture: R290-mineral oil, POE and PAG.* IIR FRIDOC #22122. ✔ **[R290/PAG solubility magnitude]** — https://iifiir.org/en/fridoc/characteristics-of-refrigerant-oil-mixture-r290-mineral-oil-poe-and-pag-22122
- **[R11]** Restrepo, Bent & Michels. *A general model for the thermophysical properties of refrigerant/lubricant mixtures (excess Gibbs energy, Wohl equations; dew/bubble points).* AIChE Journal (2009), doi 10.1002/aic.11944. ✔ — https://aiche.onlinelibrary.wiley.com/doi/10.1002/aic.11944
- **[R12]** Ding et al. *Non-equilibrium gaseous cavitation in a tilting-pad bearing* (2021), PMC10306131. *(non-eq vs eq error figures)* — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10306131
- *Supporting / secondary (lower weight, not all individually vote-verified):* IntechOpen ch. 51347 (mass-conserving cavitation review); ASME J. Tribol. 147(2):024101 (refrigerant bearing static/dynamic coefficients); RG 313256837 (*high-pressure viscosity of refrigerant/oil systems*); Int. J. Refrigeration S0140700724004523 (2024); Oxford IJLCT 6(4):243 (two-phase viscosity).

---

## 11. Verification ledger (for the next agent)

From the deep-research run: **5 search angles → 21 sources fetched → 83 claims extracted → top 25 adversarially verified (3 independent skeptics each, 2/3 refutes kills a claim) → 24 confirmed, 1 killed → 12 synthesized findings.**

| # | Load-bearing claim | Vote | Used in |
|---|---|---|---|
| 1 | Two pressure thresholds; dissolved gas outgasses first; gaseous cavitation dominant, oil vaporization rare | 3-0 ✔ | §1, §2.2 |
| 2 | Refrigerant-oil cavitation is solubility-controlled outgassing at the saturation pressure, not oil boiling | 3-0 ✔ | §1, §2.3 |
| 3 | Desorbed gas → homogeneous two-phase film; μ AND ρ reshape pressure (3%→21% load effect) | 3-0 ✔ | §1, §5.3 |
| 4 | Compressible/barotropic Reynolds ≈ formal approximation of JFO/Elrod–Adams | 3-0 ✔ | §4.2, §6 |
| 5 | Elrod θ dual meaning (compression / fractional film content), switch g∈{0,1}; matches `pancake` | 3-0 ✔ | §4.1, §9 |
| 6 | Equilibrium closure insufficient at high speed; non-equilibrium matches experiment better | 3-0 ✔ | §2.5 |
| 7 | Gas-diffusion timescale ≫ rotor period ≫ vaporization; regime-dependent; Groper-Etsion scope limit | 3-0 ✔ | §2.5 |
| 8 | Homogeneous two-phase (VOF + NCG) validated alternative; cavitation-type proportions unsettled | 3-0 ✔ | §4.3 |
| 9 | Awad–Muzychka two-phase μ with correct limits (incl. propane data); McAdams/Cicchitti/Dukler/Beattie–Whalley menu | 3-0 ✔ | §5.3 |
| 10 | Refrigerant-oil solubility/μ/ρ are simultaneous f(T,p,composition); no single empirical rule captures all | 3-0 ✔ | §5.1 |
| 11 | R290 least soluble in PAG (~5%, ~half of mineral oil); partially miscible | 3-0 ✔ | §5.1 |
| 12 | Excess-Gibbs / Wohl VLE model predicts dew/bubble points vs composition & T (<4% err) | 3-0 ✔ | §2.4 |
| — | **REFUTED:** "a simple mixing rule suffices for refrigerant-oil mixture viscosity" | **0-3 ✗** | §5.2 guard |

**Global caveats carried forward:** (a) most quantitative evidence is *air-in-oil*, only [R4] is refrigerant-in-oil — transfer by analogy, recalibrate for R290; (b) several PDFs were 403-blocked during verification — pull full texts before coding equations; (c) non-equilibrium rate constants are fitted, no R290 values exist; (d) `pancake`'s θ uses `ρ_0` not `ρ_cav` normalization (tiny but distinct); (e) the refrigerant-oil two-phase bearing literature is *active (2014-2025)* and still evolving.

---

*End of document. Next agent: see §0 for the status-tag legend, §9.2 for the concrete code gaps, and §11 for exactly how strong each claim is.*
