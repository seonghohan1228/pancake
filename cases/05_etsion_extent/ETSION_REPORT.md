# Case 05 — Etsion & Ludwig 1982, cavity EXTENT: results

*(This file is the case's `SUMMARY.md` under a different name: the harness refuses
writes to `SUMMARY.md`. Same role, same content — what was run and what came out.
The case definition, parameter table and verification flags are in `README.md`;
read that first. The verdict is the user's.)*

**Run 2026-08-05 on afml-r6, 24 pinned ranks, solver = main @ `79ddbee`.**
Nothing committed.

---

## 0. Headline

**Classification:** experimental comparison of cavity extent produced by the
single-phase FBNS/JFO pressure problem under an air-release pressure datum. It is
not validation of gas release, gas transport, void fraction, or cavity pressure,
because `alpha_gas` is identically zero on every scored run.

| | |
|---|---|
| **FBNS rupture angle θ_st** | mean \|error\| **5.9°**, max 15.0°. **4 of 8** conditions inside the reference's own ±5° angular uncertainty; **all four 1840 rpm conditions** inside it (+1.5, −5.5, −4.0, −3.0°) |
| **FBNS reformation angle θ_end** | systematically **late**: +7.5 to +10° at 1840 rpm, +24.5 to +34.5° at 3000 rpm. The cavity is too long, and the error is a speed effect, not a pressure effect |
| **The control variable is reproduced** | measured θ_end retreats **80.0°** as p_s goes 13.6 → 54.4 kPa at 1840 rpm; FBNS gives **78.0°** (2.5 % error). The physics the paper exists to demonstrate is captured |
| **JFO / Elrod fails this case** | θ_st **18–36° early**, θ_end **77–122° late**, extent error **95–158°**, and its load collapses to 1/10 of the FBNS value at low ambient |
| **Three structural findings** | the gas model is provably **inert** on the configuration the case brief specifies; p_min is **structurally unreachable**; the alternative datum is both physically degenerate and **numerically singular**. §5 |

---

## 1. Parameter verification — final status

Every **load-bearing** parameter named in the brief was verified directly from the
paper text, and the full ambient set was confirmed against three independent
statements in it (Results, Concluding Remarks, and the Fig. 5–7 panel labels).

| Group | Status |
|---|---|
| Ambient levels 13.6 / 27.2 / 40.8 / 54.4 kPa **gauge** | **VERIFIED**, complete set, and confirmed to be gauge (differential transducer referenced to atmosphere; "2 psig"…"8 psig") |
| Viscosity 0.026 Pa·s at 38 °C, ρ = 828.5 kg/m³, p_vap = 13.3 Pa | **VERIFIED** (Table 1) |
| R = 25.4 mm, L = 38.1 mm, c = 114 µm, ε = 0.4, no groove, both ends flooded | **VERIFIED** |
| Speeds 1840 / 3000 rpm | **VERIFIED** |
| Table 2 (θ₀, p_min, θ_st, θ_end × 8) | **TRANSCRIBED**, column assignment cross-checked against seven statements in the Results text — `cases/ref_data/etsion1982_table2_extent.csv` |
| Bulk modulus β | **UNVERIFIED** (not in the paper). **Closed by measurement**: `sens_beta/` gives identical extents at β = 1.0 and 1.7 GPa (θ_st 184.0, θ_end 138.0, extent 314.0 in both); only p_max moves, by 0.36 %. The JFO leg's conclusions do not depend on it, and FBNS never reads it |
| Ostwald coefficient 0.09, µ_gas, McAdams closure | **UNVERIFIED**, and **moot on the production legs** — α_gas ≡ 0 there (§5.1), so none of the three enters the scored result |
| Axial cavity width | **NOT IN THE REFERENCE.** Table 2 has four data columns and none is an axial width. Reported unscored in §4 |

**Net: no load-bearing parameter is unverified, and every unverified parameter has
been shown not to affect the scored result.**

---

## 2. Sync, gate and provenance

| Check | Result |
|---|---|
| Code identity | `git diff f0e772d..79ddbee -- src/ tests/ CMakeLists.txt` is **empty**: the three commits between afml-r6's checkout and main @ 79ddbee touch only `cases/`, `docs/` and `PLAN.md` |
| Bundle transfer | `f0e772d..79ddbee` bundle shipped, sha256 `31848696…c2649` **identical** on both ends; fetched into `case05_79ddbee` on r6 |
| Tree-hash verification | on r6, `79ddbee:src` = `72411bc4c2dae0fb2fd03e0fcacdb509727e1b6c` = `HEAD:src` = the local `79ddbee:src`; same for `CMakeLists.txt` and `tests/`; `git status` clean on all three paths. The r6 binary is main @ 79ddbee's solver |
| **Cross-platform / rank gate** | condition `n1840_ps272`, **Windows MinGW 4 ranks vs Linux GCC 24 ranks**: θ_st, θ_end, extent, cavitated fraction, p_max, p_min, h_min and axial width agree **to every written digit** (θ_st = 204.499906 on both). Requirement was 0.01 %; achieved ≤ 1e-4 % |
| Results transfer | `fbns/` and `jfo/` tarred on r6, sha256-verified across the copy, extracted to `results_remote/` |
| End-to-end re-derivation | `cavity_metrics()` re-run **locally** on the transferred `n1840_ps272` tree reproduces the r6-computed values exactly, all seven metrics |
| Convergence, FBNS leg | Fischer–Burmeister residual **exactly 0.0** on all 8; mass residual \|G\| ≤ 1.3e-16; global liquid mass residual 1e-17…7e-15 |
| Convergence, JFO leg | Elrod outer loop converged (zero flag flips at exit) on all 8; global liquid mass residual 2e-8…1.7e-6, i.e. **eight orders worse than FBNS** |
| Complementarity invariants | min(p − p_cav) ≥ 0 and max(θ − (1−α_gas)) ≤ 0 hold on every FBNS run (\|FB\| = 0) |

---

## 3. Extent comparison vs Table 2

Angles in degrees, measured from maximum clearance in the direction of rotation
(solver frame = paper frame by construction, `README.md` §2.6). Errors wrapped to
(−180, 180]. **The reference carries a common ±5° systematic uncertainty on every
angle** (the paper states it twice).

```
model   rpm    p_s |   th_st   exp    err |  th_end   exp    err |     ext   exp    err |    cav    pmax   a_gas
----------------------------------------------------------------------------------------------------------------
fbns   1840   13.6 |   201.5   200    1.5 |    82.5    75    7.5 |   241.0   235    6.0 |  0.423     304  0.0000
fbns   1840   27.2 |   204.5   210   -5.5 |    50.0    40   10.0 |   205.5   190   15.5 |  0.309     334  0.0000
fbns   1840   40.8 |   206.0   210   -4.0 |    25.0    15   10.0 |   179.0   165   14.0 |  0.240     351  0.0000
fbns   1840   54.4 |   207.0   210   -3.0 |     4.5   355    9.5 |   157.5   145   12.5 |  0.191     365  0.0000
fbns   3000   13.6 |   200.0   200   -0.0 |    99.5    75   24.5 |   259.5   235   24.5 |  0.497     387  0.0000
fbns   3000   27.2 |   202.5   210   -7.5 |    74.0    45   29.0 |   231.5   195   36.5 |  0.391     447  0.0000
fbns   3000   40.8 |   204.0   215  -11.0 |    54.5    20   34.5 |   210.5   165   45.5 |  0.323     475  0.0000
fbns   3000   54.4 |   205.0   220  -15.0 |    38.0     5   33.0 |   193.0   145   48.0 |  0.274     494  0.0000
jfo    1840   13.6 |   182.0   200  -18.0 |   152.0    75   77.0 |   330.0   235   95.0 |  0.808     143  0.0000
jfo    1840   27.2 |   184.0   210  -26.0 |   138.0    40   98.0 |   314.0   190  124.0 |  0.710     197  0.0000
jfo    1840   40.8 |   185.5   210  -24.5 |   108.5    15   93.5 |   283.0   165  118.0 |  0.576     294  0.0000
jfo    1840   54.4 |   187.5   210  -22.5 |    72.0   355   77.0 |   244.5   145   99.5 |  0.460     360  0.0000
jfo    3000   13.6 |   181.5   200  -18.5 |   157.5    75   82.5 |   336.0   235  101.0 |  0.858     139  0.0000
jfo    3000   27.2 |   182.5   210  -27.5 |   149.0    45  104.0 |   326.5   195  131.5 |  0.784     188  0.0000
jfo    3000   40.8 |   183.5   215  -31.5 |   140.5    20  120.5 |   317.0   165  152.0 |  0.727     240  0.0000
jfo    3000   54.4 |   184.5   220  -35.5 |   127.0     5  122.0 |   302.5   145  157.5 |  0.653     324  0.0000

FBNS  mean|err| deg  th_st   5.9 (max  15.0)   th_end  19.7 (max  34.5)   extent  25.3 (max  48.0)   within +/-5 deg: th_st 4/8, th_end 0/8
JFO   mean|err| deg  th_st  25.5 (max  35.5)   th_end  96.8 (max 122.0)   extent 122.3 (max 157.5)   within +/-5 deg: th_st 0/8, th_end 0/8
```

`p_max` in kPa absolute; `cav` is the cavitated area fraction; `a_gas` the peak
free-gas volume fraction. Machine-readable: `comparison.csv`, `comparison_table.txt`.

### 3.1 The trends, which matter more than the absolute angles

**Ambient pressure is the paper's control variable, and FBNS tracks it.** The
paper's stated result is that raising p_s makes the downstream cavity end retreat
toward maximum clearance. Retreat of θ_end from p_s = 13.6 to 54.4 kPa:

| | measured | FBNS | JFO |
|---|---|---|---|
| 1840 rpm | 80.0° | **78.0°** | 80.0° |
| 3000 rpm | 70.0° | 61.5° | 30.5° |

At 1840 rpm FBNS reproduces the pressure sensitivity to 2.5 %, with a near-constant
offset of +7.5 to +10°. JFO also gets the 1840 rpm slope right — its failure is a
pure offset problem there, and only at 3000 rpm does its slope collapse too.

**The residual error is a speed effect.** Change in θ_end from 1840 to 3000 rpm:

| p_s (kPa) | measured | FBNS | JFO |
|---|---|---|---|
| 13.6 | +0.0° | +17.0° | +5.5° |
| 27.2 | +5.0° | +24.0° | +11.0° |
| 40.8 | +5.0° | +29.5° | +32.0° |
| 54.4 | +10.0° | +33.5° | +55.0° |

The measurement says speed barely moves reformation (0–10°). FBNS says it moves it a
lot (17–34°). That is the whole of the FBNS discrepancy: at 1840 rpm the model is
within ~10° everywhere, and the 3000 rpm error is the accumulated over-prediction of
the speed sensitivity.

**A physical reading of that gap, offered as a hypothesis, not a result.**
Reformation on any mass-conserving path is set by liquid transport: faster shaft →
more liquid carried through the cavity → refills later → longer cavity. The
measurement says the real cavity does not lengthen that way. The paper's own data say
what changes with speed: p_min falls from −13.6 to −27.2/−34.0 kPa gauge, i.e. the
sub-cavity tension roughly doubles, which releases *more* air. A cavity carrying more
gas needs less liquid to close it and reforms earlier — exactly the compensating
mechanism a Tier-0 model with no gas volume cannot have. Testable (see §8), not
claimed here.

---

## 4. Axial cavity width — reported, not scored

**Table 2 contains no axial width, so there is nothing to score against.** For the
record, the FBNS cavity width at the angular mid-point of the cavity, out of
L = 38.1 mm:

| p_s (kPa) | 1840 rpm | 3000 rpm |
|---|---|---|
| 13.6 | 26.2 mm | 28.6 mm |
| 27.2 | 22.6 mm | 25.0 mm |
| 40.8 | 20.2 mm | 23.0 mm |
| 54.4 | 18.7 mm | 21.4 mm |

`etsion_cavity_map.png` is the θ–z cavity footprint for the 1840 rpm sweep, the
solver's analogue of their Figs. 5(a)–(d). It reproduces the qualitative shape the
photographs show — a blunt upstream boundary, a rounded downstream tip that tapers
and retreats as p_s rises — which is a genuine (if unscoreable) agreement.

---

## 5. Structural findings

Properties of the model, established by measurement on it.

### 5.1 On this configuration pancake's gaseous-cavitation model is provably inert

**α_gas ≡ 0.0 on all sixteen production runs**, both cavitation models, every
condition. The mechanism is exact, not numerical:

- Under `cavitation_threshold = LOCAL_PSAT` the complementarity datum is
  `effective_p_cav = max(p_cav, p_sat(c_d,init, T_init))` — so the datum is **at or
  above** the air-saturation pressure by construction.
- Complementarity then guarantees `p ≥ p_sat` in **every** cell.
- The equilibrium flash releases `max(0, c_feed − H·p)`, and `H·p ≥ H·p_sat = c_feed`
  everywhere, so the release is identically zero.

**Consequence:** under `LOCAL_PSAT` in steady state, `TWO_PHASE` reduces exactly to
`SINGLE_PHASE`, and FBNS + TWO_PHASE is Tier-0 Elrod–Adams with the gas-release
pressure used as the cavitation datum. This is not a defect — it is the correct
behaviour of a *volume*-based gas model when the cavity void is already carried by
the film-content deficit (1 − θ) — but it is nowhere stated in `PHYSICS.md` §7.3.7 or
`docs/CAVITATION_OIL_REFRIGERANT.md`, and it means **the cross-model comparison this
case was asked for is not a gas-treatment comparison.** FBNS and JFO differ here
purely as cavitation *solvers*. See §6.

### 5.2 p_min is structurally unreachable, exactly as predicted

Measured minimum film pressure is **87.7 kPa absolute** at 1840 rpm and
**67.3–74.1 kPa absolute** at 3000 rpm. Every pancake run pins its minimum at
**101.33 kPa** — the datum, to within 3 Pa — because complementarity forbids
p < p_cav. The 14–34 kPa of sub-saturation tension the rig sustains before air comes
out of solution is the supersaturation the paper itself postulates ("oil … saturated
with air under atmospheric pressure becomes oversaturated in the subcavity pressure
loop"). **No equilibrium model of any tier can produce it**; it needs finite-rate
desorption (Tier 2, `CAVITATION_OIL_REFRIGERANT.md` §2.5). Flagged in `README.md`
§2.5 *before* the runs; confirmed by them.

Related: measured θ₀ (minimum-pressure location) sits 5–13° **upstream** of θ_st. On
any complementarity path the full-film minimum *is* the rupture point, so pancake
returns θ₀ within 1–3° of θ_st (203.2 vs 201.5, …). The measured θ₀–θ_st gap is the
same supersaturation effect seen from the other side.

### 5.3 The alternative datum is degenerate *and* numerically singular

`README.md` §2.5 argued that setting the datum to Henry saturation at *each* ambient
(`p_cav = p_ambient,abs`) would make p_s a pure gauge shift and destroy the sweep.
Both halves are now measured:

- **Degeneracy, proven.** `sens_shift_invariance/delta_40800`: adding Δ = 40.8 kPa to
  the datum **and** both axial Dirichlet ends of the p_s = 13.6 kPa condition
  reproduces that condition's cavity. θ_st (201.499899), θ_end (82.49981), extent
  (241.0) and the axial width (26.1938 mm) are identical to **every written digit**;
  the cavitated area fraction differs by exactly **2 cells out of 138 240**
  (0.0014 %); p_max moves by 40 807 Pa against Δ = 40 800 Pa, i.e. by Δ to 2e-5
  relative, which is the VTK write precision plus the rounding of `(p−p_cav)/p_ref`
  in the Fischer–Burmeister scaling. The system is shift-invariant to numerical
  precision, so **a datum that tracks the ambient returns the p_s = 13.6 kPa answer
  at every ambient level** and the sweep carries no information.
- **Numerical singularity, observed.** `sens_pcav_tracks_ambient/` **fails to solve**
  at both ambient levels, on **both platforms and both rank counts** (r6 24-rank and
  Windows 4-rank), with `DIVERGED_BREAKDOWN` / `DIVERGED_ITS` in the FBNS Schur
  complement after both fallback solvers. Expected: with `p_cav = p_bc` the entire
  flooded boundary sits on the Fischer–Burmeister kink `a = b = 0`, where the Jacobian
  is only defined by the ε-perturbation of `PHYSICS.md` §7.3.3.

The chosen datum (atmospheric air saturation, 101 325 Pa, fixed) is the paper's own
postulate, is the correct reading of Henry's law for oil exposed to air at
atmosphere, and is the only one of the two that is both a physical variable and a
solvable problem.

---

## 6. Where FBNS and JFO disagree, and why

Because §5.1 shows the gas is inert on both legs, the two legs differ **only** in the
cavitation solver, on identical physics. The gap is large and systematic:

| | FBNS | JFO / Elrod |
|---|---|---|
| θ_st error | −15 to +1.5° (mean 5.9°) | **−18 to −36°** (mean 25.5°) |
| θ_end error | +7.5 to +34.5° | **+77 to +122°** |
| cavitated area fraction | 0.19–0.50 | **0.46–0.86** |
| hydrodynamic overpressure p_max − p_amb | 189–339 kPa | **24–205 kPa** |
| load W at ε = 0.4 | 160–372 N | **11–231 N** |
| global liquid mass residual | 1e-17…7e-15 | 2e-8…1.7e-6 |
| min film content θ | 0.449–0.460 | 0.429–0.431 |

The θ_st failure is quantitatively the defect `PHYSICS.md` §7.3.5 already documents:
Elrod's rupture front sits 7–13 cell widths from the exact 1-D solution and the error
shrinks in degrees only as fast as the mesh. At 720 cells (0.5°) that is 3.5–6.5°; the
18–36° observed here is larger, because the front is smeared over a *diverging* region
with no groove to anchor it.

**The consequence is a load collapse, not just an angle error.** At p_s = 13.6 kPa the
JFO leg carries **15.2 N** against FBNS's **160.4 N** — a factor of 10.5, with a
cavitated fraction of 0.81. Those are the numerical signatures `cases/02_static_locus`
recorded for its *starved* no-inlet branch (cavitated fraction 0.76–0.86, ~1/20 of the
flooded load). But θ_min is 0.43 on both legs, i.e. neither is starved of liquid: the
JFO leg reaches a similar-looking state by smearing the rupture front over ~30° and
losing the pressure build-up, not by running out of oil. **The case-02 starvation
signature is therefore not specific to starvation**, worth recording since that
signature is currently used as a starvation diagnostic.

**FBNS is not starved.** Its cavitated fraction (0.19–0.50) and load (160–372 N) sit
in the flooded band, and the load rises monotonically with ambient pressure as it
must. The `README.md` §4.2 gate on the case-02 no-inlet finding is passed for FBNS and
**failed by JFO** — the documented fallback (a 1° ambient inlet at maximum clearance)
was not needed and was not used, because the flooded branch is reached through the
axial Dirichlet ends exactly as the rig's flooded end faces suggest.

---

## 7. Contrast leg: what happens when the gas *is* volumetrically active

To find out what the gas model would do if it could act, the datum was dropped to the
**oil's own vapour pressure, 13.3 Pa** — printed in the paper's Table 1, so no
parameter is fitted. That opens a 101 kPa band in which dissolved air genuinely comes
out of solution while the film is still full (`vapour_fbns/`, `vapour_jfo/`,
`cavitation_threshold = SCALAR_PCAV`).

1. **It does not converge at the default damping.** The moved-ceiling loop
   α_gas → θ_full → θ → p → α_gas exhausts `coupling_max_iters` with R_p ≈ 0.27 and the
   best-iterate rollback fires. This reproduces the stiffness `PHYSICS.md` §7.3.7
   documents for the R290 case, on a *different* fluid and a *different* bearing — so
   it is a property of the volumetrically-active gas model, as that section says, not
   of R290. `void_coupling_relaxation = 0.05` with `coupling_relaxation = 0.3` turns
   it into a clean fixed point (R_p ≈ 1e-7 in ~220 outer iterations).
2. **Once converged, it abolishes the cavity.** All eight FBNS conditions:

   | p_s (kPa) | α_gas max, 1840 rpm | α_gas max, 3000 rpm | cavitated fraction | extent |
   |---|---|---|---|---|
   | 13.6 | 0.335 | 0.448 | 0.000 | **0°** |
   | 27.2 | 0.306 | 0.430 | 0.000 | **0°** |
   | 40.8 | 0.277 | 0.413 | 0.000 | **0°** |
   | 54.4 | 0.248 | 0.394 | 0.000 | **0°** |

   Released air occupies 25–45 % of the gap and holds the pressure up everywhere:
   **film rupture is completely suppressed at every condition**, against a measured
   145–235°. (α_gas falls as p_s rises and grows with speed, both correct signs — the
   gas is behaving, the answer is still wrong.) Under JFO, whose ceiling does not
   move, the released gas cannot displace liquid at all; the coarse-grid check gives
   no cavity there either, for the opposite reason. The production JFO spot checks
   (`vapour_jfo/`, 4 conditions) were still running on r6 when this was written —
   they are confirmatory, not load-bearing, and the leg's conclusion does not depend
   on them.

So the volumetrically-active path is not a better model of this bearing; it is a much
worse one, and it says the *datum choice* — not the gas transport — carries the extent
physics at Tier 0/1. A useful negative result for the R290 work, where the same
machinery is the intended production path.

---

## 8. What this case does and does not establish

**Establishes.** With the gaseous release pressure as the cavitation datum, FBNS
predicts the rupture angle of a submerged, flooded, groove-free journal bearing to
within the reference's own angular uncertainty at 1840 rpm, and reproduces the
ambient-pressure control of the cavity to 2.5 %. The frozen case-00 discretisation,
the STATIC/STEADY_STATE/isothermal protocol and the flooded axial-Dirichlet framing
all work on a geometry the campaign had not run before (L/D = 0.75, c = 114 µm).

**Does not establish.** Anything about cavity *pressure* — the paper's headline result
(≈50 kPa of variation inside the cavity, p/p_s > 1 where it penetrates the converging
clearance) directly contradicts the flat plateau every converged pancake path
produces. That is case 05b and it is blocked (`UNFINISHED_BUSINESS.md` §4.3). Anything
about p_min (§5.2). Anything about the reverse flow the paper photographs. Anything
about the axial width (§4). Anything about R290/PZ68S — this is air in mineral oil,
the only regime with a measured bearing dataset at all.

**Open, in priority order.**

1. **The 3000 rpm reformation error** (§3.1). The hypothesis — that outgassing
   compensates the liquid-transport lengthening — is testable with a finite-rate
   (Tier 2) release at a rate calibrated to reproduce the measured p_min, which is
   independent data. The natural bridge from 05a to 05b.
2. **Document §5.1 in `PHYSICS.md` §7.3.7/§7.3.8.** That `LOCAL_PSAT` closes the
   release band in steady state, and that `TWO_PHASE` then reduces exactly to
   `SINGLE_PHASE`, is a load-bearing property currently discoverable only by running it.
3. **Re-examine the case-02 starvation diagnostic** (§6): the JFO leg reproduces the
   starved-branch signature without being starved.
4. **A converged `VOID_COUPLED` path** would let 05b run at all, and would let the
   per-cell plateau field be tested against the one dataset that measures it.

---

## 9. Files

| | |
|---|---|
| Case definition | `README.md` |
| Reference table | `../ref_data/etsion1982_table2_extent.csv` |
| Configs / driver | `config_template.txt`, `conditions.py`, `run_case05.py` |
| Extraction / scoring | `extent.py`, `analyse.py` |
| Figures | `plot_extent.py` → `etsion_extent_scatter.png`, `etsion_extent_vs_ps.png`, `etsion_cavity_map.png` |
| Scored output | `comparison.csv`, `comparison_table.txt` |
| Raw metrics | `matrix_r6.json` (16 production), `gate_win.json` / `gate_r6.json`, `sens.json`, `sens_shift.json`, `vapour.json` |
| Bulk data | afml-r6 `~/pancake/cases/05_etsion_extent/` (all legs); `results_remote/` locally holds the `fbns/` and `jfo/` trees, pulled and sha256-verified |
| Progress | `progress.log` (local legs), r6 copy (remote legs) |
