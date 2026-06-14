# REVIEW_VERDICT.md — Combined verdict & action plan

> Synthesis of the two-agent thread in `REVIEW.md` (Reviewer A — Round 1, W1–W11/M1–M4;
> Reviewer B — Round 2, verdicts + N1–N9). This document keeps **only the issues that
> need action** and states **how to address each**. Strengths (S1–S6) are omitted by
> design. IDs in brackets trace back to `REVIEW.md`. Items are ordered by priority, not
> by original numbering.

**Consensus status:** both reviewers green-light the architecture (correct Reynolds /
Elrod–Adams core, sound phased plan). All action items below are agreed by both
reviewers except where a divergence is explicitly flagged. Reviewer B re-prioritised A's
two headline concerns (W1, W2) *behind* the new finding N2.

---

## Two corrections that change the fix (read before acting)

Reviewer B's code audit refined two of Reviewer A's claims. These change *how* you fix
them, so they come first:

1. **W1 onset pressure — use the model's own `p_sat`, not propane's vapour pressure.**
   A implied `p_cav` (0.3 bar) should rise toward propane's ~13.7 bar saturation. B
   showed the model already drives gas release off `c_d − c_{d,eq}(p,T)`
   (`fluid_properties.cpp:246`), and its *own* outgassing pressure is the Henry inversion
   `p_sat = c_d/H = 0.2175 / 2.175e-7 = 1.0e6 Pa = 10 bar` — the mixture bubble-point, not
   the pure-fluid vapour pressure. **Fix the coupling (N2), not the numeric value of
   `p_cav`.**
2. **W5 thermoviscous feedback — it's lagged, not absent.** μ(T,p,c) *is* implemented for
   `EMPIRICAL_CORRELATION` (`fluid_properties.cpp:154–165`) and fed back to Reynolds one
   timestep late. It's only truly absent in the default `CONSTANT` config. **The THD
   issue is the missing outer iteration (N1), not unimplemented feedback.**

---

## Tier 1 — MAJOR: blockers before any quantitative two-phase or THD number is trustworthy

### A. The two-phase mixture is not closed — released gas is inert to pressure and void  [N2, W2, W1, M1]
**Issue (consensus, strengthened by B at code level).** The "mass-conserving" solve
conserves **liquid only**. `reynolds.cpp` never reads `alpha_gas`/`free_gas_mass`; the JFO
density is `ρ = ρ_l·θ` (`fluid_properties.cpp:216`), which excludes the gas volume. The
outgassing source `Δm_g` (`fluid_properties.cpp:250–258`) updates *diagnostic* fields only
— it enters no solved continuity equation. Consequently:
- Released gas can neither open void (θ is decoupled, W1) nor relieve/cap pressure at
  `p_sat`, which is *the* mechanism of gaseous cavitation.
- There is no `Σα = 1` closure, so `α_l(=θ) + α_g` is unconstrained and can exceed 1
  (W2). Mixture density and the Reynolds source `∂(ρh)/∂θ` are therefore inconsistent.
- Net effect: the R290 case will **over-predict load and under-predict the gassy/cavitated
  extent**.

**How to fix.** Pick **one** closed formulation and state it explicitly:
- **(i) Single homogeneous void fraction** combining vapour + released gas, feeding one
  mixture density into Reynolds (the approach in the cited refrigerant / non-condensable-gas
  papers); **or**
- **(ii) Let `α_g` consume part of `(1−θ)`** and add a gas-continuity source
  `∂(ρ_g α_g h)/∂t + ∇·(…)` to the pressure equation so outgassing actually caps `p` at
  `p_sat(T,c_d)`.
- Tie gaseous-void onset to `p_sat(T,c_d)` (the Henry inversion, per Correction 1) and let
  `p_cav` represent only true vapour rupture; document that the two fronts differ (M1).
- Until a closed form lands, **label `alpha_gas`/`free_gas_mass` as non-conservative
  post-processing**, not solved physics.

### B. No outer coupling iteration — `STEADY_STATE` is a single segregated sweep  [N1]
**Issue (consensus).** The only iterated nonlinearity is the Elrod flag loop
(`reynolds.cpp:237`). Inter-equation coupling is one Gauss–Seidel pass per step
(`main.cpp:220–260`): Reynolds is assembled with the **previous** step's μ,ρ. For
`STEADY_STATE` the time loop runs exactly once (`main.cpp:208`), so the steady "solution"
is one sweep — pressure, temperature, viscosity and gas state are never mutually
converged, and there is no residual to reveal the splitting error. The constant-property
regressions (S6) can't catch this because there is nothing to couple in that limit.

**How to fix.** Add an **outer Picard / under-relaxed iteration** over
{Reynolds, properties, gas, energy} with a coupled residual and a convergence log. At
minimum, loop the segregated sweep to convergence in `STEADY_STATE` and document that one
transient step ≠ a converged operating point.

---

## Tier 2 — medium: correctness/fidelity defects that bias the reported numbers

### C. First-order upwind is the only convection scheme; TVD limiters are dead code  [N3]
**Issue.** Every solver passes `ConvectionScheme::UPWIND` (Elrod Couette
`reynolds.cpp:259`, energy `energy.cpp:344`, gas `gas_transport.cpp:241,259`). The
van-Leer/MINMOD limiters in `fvm.cpp:13–27,138–162` are never exercised. First-order
upwind injects false diffusion `~|u|Δx/2` that smears exactly the JFO rupture/reformation,
gas, and thermal fronts you want validated — and makes them mesh-dependent. (Also explains
the A.1↔Gumbel 1.2% residual: the two "equivalent" solvers don't even use the same Couette
scheme — Gumbel central-differences it, `reynolds.cpp:146–147`.)
**How to fix.** Wire the existing (deferred-correction-ready) TVD into the θ/T/c solves,
**or** explicitly document upwind-only. Either way, publish a grid-convergence study of
cavitation extent and `h_min`/load at 1×/2×/4× resolution before any quantitative claim.

### D. θ→p recovery is ill-conditioned in the full-film zone  [N4, caveat on S3]
**Issue.** With `p = p_cav + β ln θ` and β = 1e9, the whole load-bearing range maps to
θ ∈ [1, ~1.001]. Solver noise of O(1e-5) in θ becomes O(1e4 Pa ≈ 0.1 bar) of pressure
noise, which the load integral `∫(p−p_cav)dA` accumulates. The `snap_tol=5e-7` hack
(`reynolds.cpp:333`) suppresses flag chatter, not this amplification. This is the numerical
price of the (physically correct) β — interacts with N1 (no residual to detect it) and N3
(smeared front sits where amplification is worst).
**How to fix.** Report a **pressure-space residual** (not just θ-change) in the outer log;
scale KSP `rtol` so θ is resolved to ≪ Δp_target/β; consider relaxing in a pressure-like
variable in the full-film region.

### E. Energy capacity, convected enthalpy, and Couette shear ignore θ in cavitated cells  [N5, W6]
**Issue (B flags as partly a modelling choice — worth a deliberate decision).** Energy
capacity uses full gap `ρ·cp·h` with no θ (`energy.cpp:301`); enthalpy flux uses full
`h_face`/`ρcp_face` (`energy.cpp:316–321`); friction torque and viscous force use full
`μωR/h` in every cell (`reynolds.cpp:519,543`). The gas solve *does* weight by `ρ=ρ_l·θ`
(`gas_transport.cpp:31,193`), so the treatment is internally inconsistent. In the striated
zone (θ down to `theta_min=1e-6`) this overstates thermal mass, advected enthalpy, and
shear drag by up to ~1/θ — inflating the reported temperature field and
friction-torque/power-loss (first-order THD deliverables).
**How to fix.** Weight energy capacity/advection and cavitated-zone Couette shear by θ (or
`g(θ)` plus a striation model), **or** explicitly document and bound the full-film
thermal/shear assumption in the cavitated region.

### F. Thermoviscous coupling is lagged and depth-averaged  [W5, W6]
**Issue.** Per Correction 2, μ(T,p,c) feedback exists for `EMPIRICAL_CORRELATION` but is
explicit/lagged and (in steady state) never iterated — folds into the N1 fix. Separately,
the THD model is a single depth-averaged film temperature with lumped wall HTCs
(`energy.cpp:124`) and fixed wall temps — it loses the wall-to-film drop, the
journal/bush heat split, and inlet hot-oil mixing.
**How to fix.** (a) Resolve the lag via the N1 outer iteration. (b) Default `config.txt`
should ship `ISOTHERMAL` until coupling is iterated, **or** clearly label the energy solve
as diagnostic-only. (c) Plan an assumed across-film (parabolic/Φ-weighted) profile and at
least a 1D conjugate wall model before quantitative temperature claims.

### G. Two-phase property laws used outside their valid range  [W3, W4]
**Issue.** (W3) Einstein `μ = μ_l(1 + 2.5α_g)` (`fluid_properties.cpp:206`) is a dilute
result, but `gas_alpha_max = 1.0` lets `μ → 3.5μ_l` — gas *thickens* the film, which is
directionally backwards. (W4) Linear Henry `c_{d,eq}=H·p` calibrated to one (10 bar,
40 °C) point is used at `dissolved_gas_initial = 0.2175` (~22% mass), far from dilute;
extrapolation is unreliable.
**How to fix.** Cap `gas_alpha_max` well below 1 (≈0.5–0.6, Krieger–Dougherty), or adopt a
mixture viscosity that limits to gas viscosity as `α_g→1`; make the dilute bound a hard
validation gate. For solubility, keep linear Henry only as a local tangent and prefer the
`TABLE` path seeded from measured R290/POE isotherms (Purdue datasets); state the validity
window.

### H. Lagged FSI coupling and no linearized K/C coefficients  [W9]
**Issue.** The bearing equation of motion uses a lagged (explicit) fluid force; with
`k=c=0` and a force frozen per step (`journal_motion.cpp:136–145`), RK4 degenerates to
constant-acceleration, so `motion_time_method` currently buys nothing. Explicit lag in a
stiff thin-film FSI is instability-prone and `dt`-sensitive. No linearized stiffness/damping
(Kᵢⱼ, Cᵢⱼ) are extracted — the standard rotordynamic deliverable for oil-whirl assessment.
**How to fix.** Document a `dt` stability guideline or default to a Picard/implicit-force
coupling for stiff supports; add a small-perturbation post-process to report Kᵢⱼ/Cᵢⱼ and a
whirl/stability margin.

---

## Tier 3 — correctness-of-expectations, dead code, and config/doc drift (low effort, high honesty value)

### I. Inert time-scheme selectors  [N6]
`pressure_time_method` and `temperature_time_method` are parsed, stored, GUI-exposed, and
saved, but never read — the field solves hardwire backward Euler. Both shipped configs
select `EULER_EXPLICIT`, silently ignored; `CRANK_NICOLSON`/`RK2`/`RK4` have no effect.
**Fix:** implement the selected schemes for the field equations, **or** delete the dead
selectors and document backward-Euler-only.

### J. `STEADY_STATE` silently disables the gaseous-cavitation model  [N7]
With `STEADY_STATE`, `gas_dt=0` so both gas transports and `update_gas_state` no-op; a
steady `GAS_CAVITATION_MIXTURE` run is frozen at `dissolved_gas_initial` with `α_g≡0` and
**no warning**. A user could publish "R290 steady" numbers with the gas model off.
**Fix:** emit a hard warning/error for `GAS_CAVITATION_MIXTURE ∧ STEADY_STATE`, or compute
a frozen local-equilibrium gas state so a steady run carries a representative gas field.

### K. Config ↔ PHYSICS.md drift on the flagship R290 case  [N8, W7]
`PHYSICS.md §6` documents open `p=p_cav` axial ends, but `config.txt` ships submerged
`bc_z_*_val=1.5e6`. Worse, `PHYSICS.md §13.2.1` lists the R290 case as `p_cav=1e5`,
`bc=1e5` (1 bar sump), while the shipped `config_r290_pz68.txt` has `p_cav=3.0e4` and
`bc=1.0e6` (10 bar — 10× higher), and that 10-bar sump sits *exactly* at the model's
`p_sat`, so boundary oil enters saturated and most interior cells are supersaturated
(which, by N2, never affects pressure).
**Fix:** reconcile the `PHYSICS.md` table with the shipped configs; state the intended
sump-vs-saturation relationship; add a startup log of derived quantities (`p_sat`, boundary
saturation ratio, open vs submerged end).

### L. Pressure-only inlets cannot represent feed starvation  [W8]
Inlets impose a pressure (`reynolds.cpp:85`) with flooded init, so the leading edge is
always full-film; starvation (θ<1 at inlet, common in grooved/jet-fed/refrigerant bearings)
can't be modelled.
**Fix:** offer a mass-flux / inlet-film-fraction inlet option (prescribe `θ_inlet` or feed
`ṁ`); document the current model as flooded-inlet only.

---

## Tier 4 — validation & runtime guards (the overarching open gap)

### M. No quantitative validation archived  [W11]
Repo evidence is unit/regression tests only; no match to a published pressure/
cavitation-extent curve or measured bearing.
**Fix — run and archive:** (1) 1D slider JFO profile vs analytic; (2) a published
journal-bearing cavitation-extent comparison (e.g. Vijayaraghavan & Keith 1989); (3) a
logged global mass-residual time history; (4) R290 case vs the refrigerant-oil two-phase
literature trends.

### N. No global mass-conservation diagnostic across the operator split  [N9]
Liquid (θ + clamps/snap), dissolved gas, and free gas advance in three segregated solves
with independent open-boundary fluxes and non-mass-neutral clamps; nothing reconciles total
propane, and the PLAN's "O(1e-12)" goal is never logged.
**Fix:** log per-step global liquid-mass and total-gas residuals (domain integral +
boundary flux + exchange) and gate cavitation/gas results on them.

### O. No `Re` / thin-film validity guard at runtime  [M3]
`validate()` (`config.hpp:606`) has no regime check; high-speed/low-viscosity cases could
silently leave the laminar thin-film regime.
**Fix:** compute and report `Re` at startup and warn when the thin-film/laminar assumptions
are exceeded.

---

## Acknowledged future work — note, don't block

- **Conjugate heat transfer to the solids; T-dependent gas density in energy capacity**
  [W6, M2] — deferred; required before quantitative temperature claims.
- **EHD / thermal deformation (Phase D.3/D.4)** [M4] — Barus piezo-viscosity is present;
  elastic/thermal deflection is not. Can bias `h_min`/peak pressure for highly loaded or
  high-pressure refrigerant bearings.
- **Gumbel/half-Sommerfeld fallback** [W10] — **no action needed**: it is non-mass-conserving
  but the guard `GAS_CAVITATION_MIXTURE ⇒ ELROD_ADAMS` already exists (`config.hpp:758`).
  Just keep it out of any quantitative cavitation study.

---

## Suggested order of work

1. **Decide and implement the mixture closure (A / N2)** — everything two-phase downstream
   depends on whether gas shares void with θ or is post-processing.
2. **Add the outer coupling iteration + coupled residual log (B / N1)** — unblocks
   trustworthy THD/gas steady states and exposes D/N4.
3. **Numerics honesty pass (C, I, J / N3, N6, N7)** — cheap: wire-or-document TVD, fix the
   time-scheme selectors, warn on steady+gas.
4. **Config/doc reconciliation + diagnostics (K, N, O / N8, W7, N9, M3)** — restores the
   validation story's credibility.
5. **Validation campaign (M / W11)** and the remaining fidelity items (E, F, G, H).
