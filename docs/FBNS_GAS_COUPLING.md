# FBNS + gaseous cavitation: design memo

Status: **implemented** (2026-08-05, branch `codex/fbns-gas`).
`cavitation_model = FBNS` with `fluid_property_model = TWO_PHASE` now runs; see
PHYSICS.md §7.3.7 for the landed model and §8 below for what this memo got right,
what it got wrong, and what is still guarded. The analysis in §§1-4 is kept as
written, unamended, so it can be checked against the outcome.

- Step 1 of §4, the shared full-film predicate — **DONE** 2026-08-05, branch
  `codex/full-film-predicate` (28ac88a): `FullFilm` in `src/full_film.hpp`, all
  fifteen sites routed, byte-identical across five regression cases.
- Steps 2-4 — **DONE** 2026-08-05, branch `codex/fbns-gas`, except the
  `VOID_COUPLED` half of step 2, which stays guarded (§8).

Written while landing FBNS phase 2 (energy equation + `SINGLE_PHASE` properties),
branch `codex/fbns-cavitation`.

## 1. The question

FBNS imposes the JFO conditions pointwise as a complementarity problem:

$$p - p_{cav} \ge 0, \qquad 1 - \theta \ge 0, \qquad (p - p_{cav})(1 - \theta) = 0,$$

solved through the Fischer-Burmeister function $FB(a,b) = a + b - \sqrt{a^2+b^2}$
with $a = (p - p_{cav})/p_{ref}$ and $b = 1 - \theta$.

With a resolved free-gas phase occupying a volume fraction $\alpha_{gas}$ of the
gap, the liquid can no longer fill it. The full-film bound becomes

$$\theta \le 1 - \alpha_{gas}, \qquad\text{i.e.}\qquad b = 1 - \alpha_{gas} - \theta \ge 0 .$$

Does that compose with the FB function and the Schur-complement solve?

## 2. Finding 1 — the constraint composes cleanly

Yes, provided $\alpha_{gas}$ is **lagged** (frozen through one Newton solve, which
is what the existing operator-split gas path already does: `GasTransport::solve`
then `FluidProperties::update_gas_state` run *after* the flow solve in the outer
loop, so $\alpha_{gas}$ is a coefficient, not an unknown).

Under that lag the new $b$ is affine in $\theta$ with the *same* derivative as
before:

$$\frac{\partial b}{\partial \theta} = -1 .$$

Therefore:

- $J_{G,b} = B$ is **unchanged**. The mass-balance block already differentiates
  with respect to the cavity fraction $\alpha = 1-\theta$ and only ever uses
  $\partial(\text{outflow})/\partial\alpha = -\partial(\text{outflow})/\partial\theta$.
- $J_{F,p}$ and $J_{F,b}$ keep their FB forms, so the column swap, the
  $\mathrm{cond}(A_F) \le 4 + 2\sqrt2$ bound (Woloszynski et al. 2015, Thm. 1) and
  the Schur complement are untouched.
- Only three lines of `src/fbns.cpp` change in substance:
  - `alpha_max`: `1 - theta_min` becomes `max(1 - alpha_gas(i,j) - theta_min, 0)`;
  - the residual: `b = 1 - alpha_gas(i,j) - theta(i,j)`;
  - the update and the cold start: `theta = 1 - alpha_gas - b`, and the cold-start
    fill becomes `1 - alpha_gas` instead of `1`.

The Poiseuille mobility needs nothing new: the free gas already reaches it through
the mixture viscosity `FluidProperties::gas_mixture_viscosity` writing the `mu`
field, which the face mobilities read per cell (§7.3.6 of `PHYSICS.md`).

So the answer to the design question posed in the phase-2 brief is: **the bound
composes; Newton is not destabilised by it.**

## 3. Finding 2 — the blocker is `theta >= 1`, not the constraint

The whole codebase identifies "this cell holds a full film" by the literal test
`theta(i,j) >= 1.0`. As of this writing that test appears at 15 sites in four
files:

| File | Sites | What the switch gates |
|---|---|---|
| `src/reynolds.cpp` | 5 | Poiseuille flux in the Elrod operator and in `calculate_velocities` |
| `src/energy.cpp` | 4 | Poiseuille enthalpy flux, dissipation, STRIATED liquid weighting |
| `src/gas_transport.cpp` | 3 | Poiseuille transport of dissolved gas, axial boundary velocity |
| `src/diagnostics.cpp` | 3 | face flux and boundary flux in the mass balance |

Move the bound to $\theta \le 1 - \alpha_{gas}$ and, the instant any gas comes out
of solution, **every cell in the domain fails `theta >= 1`**. Every one of those
switches turns off at once. The solver would then report a fully cavitated bearing
that carries no Poiseuille flow, no Poiseuille enthalpy, and no dissolved-gas
transport anywhere — silently, with a converged Newton residual, because the FBNS
solve itself would be perfectly happy.

That is the reason for the guard. It is not a numerical-stability finding; it is a
representation finding. The film-content field is currently overloaded: it is both
"how much liquid is in this cell" and "is this cell ruptured", and free gas
separates those two meanings for the first time.

## 4. What lifting the guard requires

In order:

1. **Introduce an explicit rupture indicator.** — **DONE 2026-08-05.** Either a
   `full_film` field written by whichever cavitation solver ran, or a single shared
   predicate `is_full_film(fields, cfg, i, j)` that reads `theta >= 1 - alpha_gas - eps`
   under TWO_PHASE and `theta >= 1` otherwise. A predicate is preferable: it is one
   definition, it cannot go stale, and it makes the Elrod path byte-identical by
   construction. Replace all 15 sites.

   Landed as `FullFilm` (`src/full_film.hpp`). The ceiling lives in one function,
   `FullFilm::bound(fields, cfg, i, j)`, which today returns exactly `1.0`; the
   TWO_PHASE branch is the single line written out in its docstring. All fifteen
   sites (plus two adjacent ones in the same files) call
   `is_full_film` / `is_full_film_face` or their 1/0 `gate` / `face_gate` forms.
   Verified byte-identical on JFO steady, FBNS steady, JFO+energy, JFO+energy
   +STRIATED and a transient TWO_PHASE run: 41 output files, 41 sha256 matches.
2. **Decide the cavity pressure with gas present.** In the vapour cavity FBNS pins
   $p = p_{cav}$. With a released species the plateau is the local bubble point,
   which is what `cavitation_threshold = LOCAL_PSAT` already supplies as a scalar
   datum (`effective_p_cav()`), and what `gas_pressure_coupling = VOID_COUPLED`
   supplies as a per-cell lagged field on the Elrod path. FBNS would need
   `p_cav` to become a per-cell array in `cell_outflow`, the residual, and
   `axial_boundary_pressure`. Mechanically easy; the modelling question of which
   closure is right is settled in `docs/CAVITATION_OIL_REFRIGERANT.md`.
3. **Expect the outer loop to need under-relaxation.** The Elrod path already
   limit-cycles on the $\alpha_{gas} \leftrightarrow p \leftrightarrow \theta$
   feedback and carries `void_coupling_relaxation = 0.3` for exactly that reason.
   FBNS removes the *inner* iteration but not this outer one, so the same lever is
   likely to be needed; nothing about FBNS makes the Picard loop better
   conditioned.
4. **Steady state has no gas transport.** `GasTransport::solve` returns
   immediately for $\Delta t \le 0$, so a `STEADY_STATE` FBNS + TWO_PHASE run only
   makes sense with `steady_gas_model = EQUILIBRIUM` (local flash). That should be
   enforced, not warned about, for the FBNS path.

Steps 1 and 2 are the real work and they are cross-cutting: step 1 touches four
files that the Elrod, Gumbel and diagnostics paths all share, so it needs its own
byte-identity regression before anything gas-related is attempted on top of it.
Bundling it into the FBNS phase-2 change would have put a repo-wide refactor of the
full-film test behind a feature nobody can run yet.

## 5. Recommendation

Land step 1 (the shared `is_full_film` predicate) as a standalone, byte-identical
refactor on its own branch, with the existing Elrod/Gumbel regression as its
acceptance test. Steps 2–4 then become a small change confined to `fbns.cpp` plus
config validation, and the gas guard can be lifted with the 20-line diff sketched
in §2.

Step 1 is now landed (see the status line). What remains for whoever picks up
step 2: change `FullFilm::bound` to return `1 - alpha_gas` under TWO_PHASE, make
`p_cav` a per-cell array in `cell_outflow`, the FBNS residual and
`axial_boundary_pressure`, and enforce `steady_gas_model = EQUILIBRIUM` for
`STEADY_STATE` FBNS + TWO_PHASE. Note that the Elrod operator carries its own
copy of the same ceiling in `Reynolds`'s `VoidParams::theta_full`; the two should
merge onto `FullFilm::bound` in that change rather than diverge further.

Until then `cavitation_model = JFO` remains the path for gaseous cavitation
(validation case 05), and FBNS covers cases 00, 00c, 01 (isothermal and THD) and
02–03.

## 8. Outcome (2026-08-05, branch `codex/fbns-gas`)

### 8.1 What the memo got right

Finding 1 held exactly. With `alpha_gas` lagged the bound is affine in `theta`
with the same derivative, both Jacobian blocks are untouched, and the landed diff
in `src/fbns.cpp` is the three substantive changes §2 predicted plus the ceiling
lambda. The complementarity invariants come out at **0 exactly** on the test case
(`tests/test_fbns.cpp` test 11) and per-cell mass conservation stays at Krylov
precision, `5.1e-15` of the peak Couette flux.

Finding 2 held too: `FullFilm::bound` was the whole of the representation work.
Moving the ceiling was one line there and zero call-site edits, and the Elrod and
Gumbel paths came out byte-identical across six regression cases / 54 output
files against 28ac88a, including JFO + TWO_PHASE and JFO + VOID_COUPLED +
LOCAL_PSAT.

§4.3 held: the outer loop does need under-relaxation. The lever already existed —
`FluidProperties::gas_state_relaxation` — and now fires whenever the gas displaces
liquid rather than only under `VOID_COUPLED`.

### 8.2 What the memo got wrong

**The cold start.** §2 says "the cold-start fill becomes `1 - alpha_gas` instead
of `1`". The *fill* generalizes, but the cold-start **test** must not: it is a
sentinel comparison against the value `main.cpp` initialises `theta` to, not a
full-film test. Generalizing it to `theta >= theta_full` declares a cold start on
every gassy step (an uncavitated gassy cell sits at `theta = 1 - alpha_gas` by
definition), throws the warm start away, and refills `theta` down to the new
ceiling in one go — a spurious `d(theta)/dt` sink that drove the pressure to
`2.9e8` Pa on the R290/PZ68 case. The landed code keeps the literal `1` and says
why, at length, at the site.

**Two config.hpp helpers were missed.** `elrod_boundary_outflow` and
`is_reformation_boundary` take film-content *values*, not cell indices, and live in
`config.hpp`, which `full_film.hpp` includes — so they could not be routed through
the predicate during step 1. Both now take the ceiling as an explicit parameter
with no default (a silent `1.0` default would reintroduce exactly the stale-ceiling
bug the predicate exists to prevent). Seven call sites updated; a no-op today,
since every caller is on a path where the bound is 1.

**Step 4 was already enforced.** `STEADY_STATE` + `steady_gas_model = FROZEN` was
already a validation error for any TWO_PHASE run, FBNS included, with an explicit
`allow_frozen_gas_steady` opt-out. Nothing to add.

**§4.2's "mechanically easy" is optimistic** — see §8.3.

### 8.3 What stays guarded, and why

`gas_pressure_coupling = VOID_COUPLED` is still rejected under FBNS.

The FB-argument side is indeed easy: `a = (p - p_cav)/p_ref` would read a per-cell
datum, the Poiseuille Jacobian `A` is unchanged because `dp/dp* = p_ref` whatever
constant offset `p_cav` carries, and the `source = A p* - G` assembly already
absorbs the extra `grad(p_cav)` term because `G` is evaluated from the real
pressure field. That is the part §4.2 saw.

What §4.2 missed is that the plateau field has no update site outside the Elrod
solver. `p_void` is advanced — and under-relaxed — inside
`Reynolds::solve_elrod`'s outer flag-update loop (`reynolds.cpp`, `if
(fields.has("p_void") && outer > 0)`), because recomputing the closure on the
instantaneous `theta` at every assembly makes the `p_void <-> theta` feedback
limit-cycle. FBNS has no flag-update loop to hang that on. Wiring it in means a new
update site with its own relaxation schedule and its own validation, on top of a
feedback loop the Elrod path already needs damping for. Delivered as a documented
limitation with a specific error message rather than a half-implementation.

Note also that `p_cav` under `cavitation_threshold = LOCAL_PSAT` is **not** a
field: `effective_p_cav()` evaluates the bubble point at the *feed* state
(`dissolved_gas_initial`, `temperature_initial`), both config scalars, so it is
uniform by construction — deliberately, so the datum injects no net force. FBNS's
scalar already reads the correct local value there. `p_void` is the only genuinely
per-cell cavitation pressure in the codebase.

### 8.4 The Elrod ceiling merge is still open

§5 recommends merging `Reynolds`'s `VoidParams::theta_full` onto
`FullFilm::bound`. That did **not** land, on purpose. `FullFilm::bound` moves only
under FBNS + TWO_PHASE (`FullFilm::gas_displaces_liquid`), because:

- under `gas_pressure_coupling = NONE` the Elrod path's free gas is volumetrically
  inert *by design* — that is the shipped, validated model of validation case 05,
  not an oversight to correct in passing;
- under `VOID_COUPLED` the Elrod operator already carries the moved ceiling in
  `VoidParams::theta_full`, where it feeds the barotropic EOS and the switch
  function. Dropping the fifteen Poiseuille gates onto it as well is a second,
  independent modelling change to a validated path, and it would break the
  byte-identity that made this change checkable.

Whoever merges them should expect to re-validate case 05, not just re-run it.

### 8.5 Known stiffness

`config_r290_pz68_quick.txt` with `cavitation_model = FBNS` converges every Newton
solve (61/61 steps, 2–11 iterations, mean 6.3) and keeps the gas mass bounded, but
its load history is transient noise for tens of steps: the saturated feed flashes
enough propane to fill ~23% of the gap per step, and the resulting ceiling drop
swings the load between 5.8e2 N and 2.4e5 N on alternate steps, decaying ~30% per
step pair. The swing is timestep-independent (identical at dt = 5e-6, 5e-8 and
5e-9 s), so it is the model's displacement rate, not a temporal-resolution
artefact. The Elrod path with `VOID_COUPLED` does not resolve this case either —
it drifts monotonically instead. Lower `void_coupling_relaxation` to damp it.
Details in PHYSICS.md §7.3.7.
