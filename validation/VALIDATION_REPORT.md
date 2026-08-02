# pancake — verification and validation report

**Date:** 2026-08-02 · **Branch:** `codex/validation-cavitation` · **Build:** MinGW64 / PETSc, 23/23 tests pass

This is the consolidated record of what has been verified, what has been
validated, what is known to be wrong, and what remains open. Every number is
traceable to a committed script and a results JSON; the "Reproduce" column names
the command.

Terminology used strictly throughout: **verification** = is the code solving the
equations right (against analytic solutions and an independent solver);
**validation** = are those the right equations (against experiment).

---

## 1. Headline verdicts

| Question | Verdict |
|---|---|
| Is the Reynolds/cavitation discretisation correct? | **Yes.** Two independent checks: Raimondi & Boyd charts (S within 5.8%, attitude within 4.08°) and a from-scratch independent solver at the benchmark geometry (loads within 1.3–6.2%) |
| Is the 52–68% load over-prediction a solver defect? | **No.** An independent solver reproduces it, and Ferron's *own* published theory over-predicts his own data by ~75–80% — worse than pancake |
| Is the energy equation implemented correctly? | **Dimensionally yes**, formulation is a correct depth-averaged form. But it is **structurally incomplete** — see §5 |
| Can the 71° temperature-peak error be fixed by tuning? | **No.** With wall temperature a prescribed constant there is no mechanism to move heat upstream. It needs a conjugate model |
| Is the gaseous-cavitation (VOID_COUPLED) path usable? | **No.** The outer coupling does not converge; answers are relaxation- and grid-dependent |
| Are temperature results quotable as-is? | **Only with the wall condition attached**, and noting the measurement is a *bush* temperature, not a film temperature |

---

## 2. Verification — is the code solving the equations right?

### 2.1 Analytic limits

| Benchmark | Result | Reproduce |
|---|---|---|
| Sommerfeld long-bearing limit | Converges to closed form at expected rate | `run_benchmarks.py` |
| Ocvirk short-bearing limit | Converges to closed form at expected rate | `run_benchmarks.py` |
| Exact Swift-Stieber free boundary | Rupture angle within 0.23° at β = 1e8 / 2880 cells | `run_raimondi_boyd.convergence_study()` |

### 2.2 Raimondi & Boyd finite-width charts — **PASS**

Fed film (supply groove at maximum gap), ε = 0.2–0.9, L/D = 1/2 and 1.

| Metric | Result | Criterion |
|---|---|---|
| Sommerfeld number | within **5.8%** | 10% |
| Attitude angle | within **4.08°** | 5° |

Reference data itself audited: the transcribed chart reproduces an independently
derived exact Swift-Stieber solution to 0.6% in S and 0.4° in attitude.

**Caveat that matters operationally:** without a supply groove the same solver has
essentially *no* load capacity (0.007–1.8 N against chart loads of order kN),
because the Elrod switch suppresses Poiseuille flow in cavitated cells and
nothing re-floods the film. The charts assume a fed film; so must you.

*Reproduce:* `python run_raimondi_boyd.py <pancake.exe>` → `raimondi_boyd_results.json`

### 2.3 Independent-solver cross-check at the benchmark geometry — **PASS** (new, 2026-08-02)

R&B only covers L/D = 1/2 and 1 at c/R = 1e-3. The THD benchmarks run at
L/D = 0.80, c/R = 2.9e-3 (Ferron) — outside any previously verified regime.
`christopherson.py` is a from-scratch Python finite-volume Reynolds solver with
Christopherson clamping, sharing **no code** with pancake and having no
`p_max/β` failure mode.

Self-audit of the reference solver against the R&B charts: **1.3% in S, 0.80° in
attitude** — tighter than pancake's own, so it is a legitimate reference.

Isothermal, Ferron geometry, 2000 rpm:

| ε (measured load) | pancake | independent | Δ load | Δ peak pressure |
|---|---|---|---|---|
| 0.47 (2 kN) | 3248 N | 3206 N | +1.3% | +3.1% |
| 0.67 (4 kN) | 7933 N | 7591 N | +4.5% | +3.5% |
| 0.78 (8 kN) | 15010 N | 14129 N | +6.2% | +4.1% |

*Reproduce:* `python run_cross_solver.py <pancake.exe>` → `cross_solver_results.json`

**Is the residual few-percent difference a defect? No — it is mesh convergence.**
The headline table ran pancake at 720×40 against the reference at 1440×60. At a
fixed β = 1e8, refining pancake's mesh moves it monotonically onto the reference:

| grid | pancake | independent | Δ |
|---|---|---|---|
| 720×40 | 15010 N | 14139 N | +6.2% |
| 1440×60 | 14717 N | 14129 N | +4.2% |
| 2880×60 | 14540 N | 14129 N | **+2.9%** |

The reference is already grid-converged at 1440×60 (unchanged at 2880). pancake
converges more slowly in mesh because its Couette term is an **upwinded flux on
the film content — first order**, while the reference solves for pressure
directly with the wedge as a source and has no convective discretisation at all.
Monotone convergence onto an independent solution under refinement is the
signature of discretisation error, not a defect. **No defect was found.**

**A finding that inverts an intuition, recorded so it is not re-derived:**
raising β to approach the incompressible limit the reference solves makes
agreement *worse*, not better —

| β | Δ vs reference | `p_max/β` |
|---|---|---|
| 1e8 | +2.9% | 0.058 |
| 3e8 | +3.5% | 0.019 |
| 1e9 | +9.4% | 0.0059 |

Raising β shrinks `p_max/β` toward the 4e-3 free-boundary resolution floor
(§2.4), so it trades a small compressibility difference for a larger
under-resolution error. The two effects oppose each other and β = 1e8 sits near
the optimum for this case. **Do not raise β to flatter a comparison.**

*Practical consequence:* 1440×60 (the production Ferron mesh) carries ~4% of
discretisation error in load. That is inside the R&B validation tolerance and an
order of magnitude below the 52–68% discrepancy under investigation, but it
should be stated whenever a load number is quoted to better than ~5%.


### 2.4 Free-boundary resolution — a trap that does not announce itself

Elrod carries pressure inside the film content (`p = p_cav + β ln θ`), so the
accuracy of the rupture boundary is governed by **`p_max/β`, not by any residual
the solver reports**. At β = 1e9 on a 720-cell mesh the rupture angle is **20.8°
wrong** while `mc = 0`, `ns = 0`, `R_p = 0` and the process exits successfully.

**Requirement: keep `p_max/β ≥ 4e-3`.** This also means β must be refined
together with the mesh in any grid-convergence study, or the rupture error
masquerades as grid error.

---

## 3. Validation — are they the right equations?

| Target | Quantity | Verdict |
|---|---|---|
| Raimondi & Boyd (1958) | S, attitude angle | **PASS** (§2.2) |
| Ferron et al. (1983) | peak pressure | **PASS** — −9.9% |
| Ferron et al. (1983) | peak pressure location | **PASS** — +4.2° |
| Ferron et al. (1983) | rupture position | **PASS** — inside measured bracket |
| Ferron et al. (1983) | peak temperature | **PASS adiabatic** (−1.62 K); **FAIL** metal bush (−6.10 K) |
| Ferron et al. (1983) | temperature peak *location* | **FAIL** — 71° out (§5.1) |
| Ferron et al. (1983) | load at measured ε | **OPEN** — +52–68%, shared with the literature (§4) |
| Brito et al. (2006) | load, via eccentricity-for-load | **NOT SCOREABLE** — the digitisation ambiguity spans the discrepancy (§4.1) |
| Brito et al. (2006) | pressure, flow, interface T | **INCOMPLETE** — 5 of 6 committed runs unconverged; re-run required |
| Cole & Hughes (1956) | reformation | **NOT SCOREABLE** — not sourced, and not a prediction of this solver in either inlet configuration |
| Engels et al. (2024) | measured void fraction | **NOT RUN** — the direct Elrod observable; nothing currently tests it |
| Grando et al. | refrigerant-oil equilibrium | **NOT RUN** — blocked on §5.3 |

**Do not quote** any run with `R_p > 1e-2` (units of 1e5 Pa), `mc > 1e-10`,
`ns > 0`, or the coupling-iteration cap hit. Five of six committed Brito rows
fail this; one moves the pressure 6.7% per iteration.

---

## 4. The load over-prediction — assessed, not a defect

Three independent lines of evidence:

1. **An independent solver reproduces it** (§2.3) — agreement 1.3–6.2% between
   two unrelated discretisations, while both exceed the measured load by 60–90%.
2. **Ferron's own theory over-predicts his own measurements by ~75–80%** — worse
   than pancake. His theory reaches ε = 0.67 at ~7.2 kN where he measured it at
   4 kN. Alakhramsing's FE model is displaced the same way.
3. **The "measured" eccentricity is not measured directly.** It is a
   proximity-probe displacement divided by an *assumed* clearance, so a clearance
   error enters the measurement itself. Lubricants 11(7) 308 (2023) reduces the
   published ε by 0.10 on exactly this argument (c_eff ≈ 171 µm vs 145 µm
   nominal). Independent reproductions of the same 4 kN point disagree by 0.05
   in ε.

**Ruled out by measurement, not argument:** mesh (1.4% in peak pressure between
720×40 and 1440×60), groove placement (4.8% over a 90° sweep), groove width (3%
over 9→36°), thermal treatment (adiabatic / metal-bush / isothermal agree on the
ε carrying 4 kN to within 0.02), and clearance thermal growth (needs the bush
56 K hotter than the shaft to supply the required factor).

**Recommended framing:** validate on peak pressure, its position, and the
temperature profile — the quantities all three published models match. Quote
W(ε) as an open discrepancy shared with the literature, with the
clearance-derived eccentricity as the leading suspect. If a load criterion is
wanted, reproduce Ferron's *theory* curve, not his measurements.

---

### 4.1 Brito — the load comparison cannot currently discriminate

Same lesson as Ferron, applied before the fact: rather than prescribing a
digitised eccentricity and calling the load a prediction, solve for the
eccentricity that carries the load Brito **applied** and therefore knew exactly.
The uncertain quantity becomes a check instead of an input.

Eccentricity carrying the applied 10 kN, `Pf` = 70 kPa gauge, `Tf` = 40 °C:

| N (rpm) | ε for 10 kN, c at 20 °C | ε hot (+8.8%) | ε digitised (polar) | Δε |
|---|---|---|---|---|
| 1000 | 0.648 | 0.678 | 0.761 | −0.113 |
| 3000 | 0.426 | 0.461 | 0.576 | −0.150 |
| 4000 | 0.367 | 0.401 | 0.558 | −0.191 |

Three readings, in order of how much they should change your mind:

1. **The digitisation ambiguity is larger than the discrepancy.** `brito.py`
   records that the earlier *Cartesian* reading of Fig. 6 — a polar chart — gave
   eccentricities about 0.15 lower than the polar re-digitisation now in use.
   That puts the Cartesian values near 0.61 / 0.43 / 0.41, and the solver answers
   0.648 / 0.426 / 0.367. **The solver tracks the Cartesian reading to within
   ~0.04 at every speed and sits 0.11–0.19 below the polar one.** This does not
   prove the polar reading wrong — using the solver to pick the digitisation
   would be circular — but it does mean **the Brito load comparison cannot
   discriminate until Fig. 6 is re-read from the source.** Report it as
   unscoreable, not as a 65–106% error.
2. **The hot-clearance bracket closes 0.030–0.035**, roughly a quarter of the gap
   at 1000 rpm and a fifth at 4000 rpm. Real, insufficient, and it rests on
   *assumed* metal temperatures (bush ~55 °C, shaft ~60 °C) that Brito did not
   measure. It is a bracket, never a correction.
3. **The deficit grows monotonically with speed** (−0.113, −0.150, −0.191),
   which is the qualitative signature of a thermal mechanism. It is *not* a U²
   growth, so it does not by itself confirm the cross-film hypothesis of §5.2 —
   and eccentricity is a compressed measure of a load error, so the scaling
   should be re-tested on load at fixed ε after the Phase 2 energy fixes.

Convergence: the two 1000 rpm rows did not meet the convergence criterion and are
marked accordingly in `brito_eps_for_load.json`; the 3000 and 4000 rpm rows did.

*Reproduce:* `python -c "import run_brito; run_brito.eps_for_load_sweep(EXE, WORK)"`
→ `brito_eps_for_load.json`

## 5. Known limitations of the physics

### 5.1 No conjugate wall model — cause of the 71° temperature-peak error

`Energy::apply_wall_transfer` is Newton cooling to a **prescribed constant** wall
temperature. No bush or shaft conduction; no circumferential heat spreading by
the rotating shaft. With no path for heat to travel upstream against the flow,
the temperature can only climb monotonically from groove to groove — which is
exactly the failure shape observed. **No value of `h_wall` fixes this**; the two
runnable settings are unphysical limits (perfect insulator, or infinite sink
pinned at supply temperature).

### 5.2 Uniform temperature across the film is invalid on these rigs

The energy equation is depth-averaged with viscosity at the mean temperature.
That is valid only when Br ≪ 1, reduced Péclet ≪ 1, and the walls are at equal
temperature. **All three fail here:**

| Criterion | Value | Requirement |
|---|---|---|
| Brinkman `μU²/(k ΔT)` | ≈ 2.6 | ≪ 1 |
| Reduced Péclet `ρcUh²/(kR)` | ≈ 10–30 | ≪ 1 |
| Wall temperatures | unequal | equal |

Dissipation-driven core superheat `μU²/(8k)` is ≈ 3 K (Ferron, 2000 rpm) to ≈ 6 K
(Brito, 4000 rpm) and **scales as U²**, worth ~10–20% of load at high speed
through the ~3.4%/K viscosity sensitivity.

Wall heat flux and effective viscosity both follow from the same missing `T(z)`
profile, so they cannot be fixed independently. This yields a **discriminator**:
a residual error that grows as U² is cross-film physics; a speed-independent one
is geometry or property data.

### 5.3 Other recorded gaps

- **Cavitated-zone thermal inventory is not θ-weighted** (defect D14): capacity,
  conduction and convected flux use the full film in ruptured cells, overstating
  the liquid mass flux through the ~100° cavitated arc by ~2.5×.
- **No hot-oil carry-over**: groove cells are hard-pinned to supply temperature.
  Sized in-repo at ~6% of load; a real groove mixes 60–80% recirculated oil.
- **VOID_COUPLED does not converge**: `R_p` 0.57–7.2 against a 1e-9 tolerance,
  relaxation- and grid-dependent. Fixing it is a numerical-methods change (fold
  the flash into the Elrod switch, or Newton on the coupled residual), not a bug
  fix.
- **16 verified defects (D1–D16) remain open** — see `HANDOFF.md`. D6 was fixed
  2026-08-02; the rest stand.
- **No thermal or elastic deformation** of the clearance.

---

## 6. Reference-data provenance and its limits

Errors here are as consequential as solver errors, so provenance is recorded per
number in `ferron.py` and `brito.py`. Three findings worth carrying into any
report:

1. **The competing Ferron geometry is a different experiment.** `L/D = 0.55,
   C/R = 1.5e-3` in circulation belongs to Tonnesen & Hansen (1981) — also a
   100 mm bearing loaded to ~9 kN. Ferron is L/D = 0.80, c = 145 µm **radial**,
   confirmed against three independent open sources. Note 1.5e-3 is also what you
   get by misreading 145 µm as diametral, so the wrong number has two routes in.
2. **The eccentricity-vs-load data is 2000 rpm.** The source panel is unlabelled;
   confirmed via an independent reproduction captioned "At 2000 rpm", agreeing
   within 0.03 in ε at every load.
3. **The measurement is a bush temperature 0.5 mm below the surface**, not a film
   temperature — not the quantity pancake reports. The two open reproductions
   also disagree by ~30° on where the peak sits within a 250–285° plateau, which
   bounds how precisely a peak-location error can be scored against this data at
   all. (The observed 71° error still exceeds that bound, so it is real.)

Every digitised value carries its reading uncertainty next to it in the data
modules; `brito.py` is the template.

---

## 7. Figures

| File | Shows |
|---|---|
| `validation_cross_solver.png` | pancake vs independent solver vs measurement, load against eccentricity |
| `validation_load_in_context.png` | the over-prediction alongside the published models that share it |
| `validation_raimondi_boyd_S.png` | Sommerfeld number against the charts |
| `validation_raimondi_boyd_attitude.png` | attitude angle against the charts |
| `validation_ferron_pressure.png` | mid-plane pressure profile against measurement |
| `validation_ferron_temperature.png` | temperature profile against measurement |
| `validation_ferron_load.png` | load against measured eccentricity |
| `validation_free_boundary_convergence.png` | rupture angle against `p_max/β` |

*Regenerate:* `python make_report_figures.py` and `python plot_validation.py`

---

## 8. What to do next

Priority order, with full detail in `../VALIDATION_PLAN.md`:

1. **Re-run Brito** — five of six committed cases are unconverged, and the
   eccentricity rests on an unverified polar re-digitisation.
2. **Phase 2 energy fixes** — θ-weight the cavitated thermal terms, scale Couette
   dissipation by θ, add groove carry-over, then the cross-film `T(z)` profile
   with conjugate wall coupling and a generalized Reynolds equation.
3. **Grid and time-step independence** with GCI, holding `p_max/β ≥ 4e-3`.
4. **Engels et al. (2024) void fraction** — the direct Elrod observable, and the
   one class of measurement nothing here currently tests.
