# FBNS + gaseous cavitation: design memo

Status: **not implemented, deliberately.** `cavitation_model = FBNS` with
`fluid_property_model = TWO_PHASE` is rejected at config validation. This memo
records why, what was actually checked, and what it would take to lift the guard.

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

1. **Introduce an explicit rupture indicator.** Either a `full_film` field written
   by whichever cavitation solver ran, or a single shared predicate
   `is_full_film(fields, cfg, i, j)` that reads `theta >= 1 - alpha_gas - eps`
   under TWO_PHASE and `theta >= 1` otherwise. A predicate is preferable: it is one
   definition, it cannot go stale, and it makes the Elrod path byte-identical by
   construction. Replace all 15 sites.
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

Until then `cavitation_model = JFO` remains the path for gaseous cavitation
(validation case 05), and FBNS covers cases 00, 00c, 01 (isothermal and THD) and
02–03.
