# The Picard cap defect: fix, measurement, and what it did (and did not) change

> **Scope correction (2026-08-07):** the Picard conclusions remain valid for the
> runs tested, but their favorable Ferron temperature score is superseded. The
> FLOW_WEIGHTED closure used a grid-divergent inlet diagnostic; corrected case
> 00d gives 5.06 L/min and +3.82 K mean / 4.02 K RMS on G6.

Branch `codex/picard-best-iterate`, from `main` @ `79ddbee`.
Solver commit shipped to afml-r6 and hash-verified: see §6.
Date: 2026-08-05.

---

## 0. One-paragraph answer

The outer (Picard) coupling loop really did deliver the wrong iterate: on every
conjugate Ferron run it exhausted `coupling_max_iters = 200` and returned a state
**6.3x worse in residual** than one it had already passed through around outer 56.
That is now fixed — the loop tracks the lowest-residual iterate, returns it when
the cap is exhausted, and says `NOT CONVERGED` out loud. The cause of the cycle is
a loop gain marginally above one: `coupling_relaxation = 0.95` already converges
in 48 sweeps where 1.0 never does, and a new adaptive schedule gets there in 52
with no tuned value. **But the delivered numbers do not move.** Across the
iterate-200 state, the rolled-back best iterate and nine converged relaxation
settings, the bush-side wall temperature profile agrees to <= 1e-3 K, `W` to
0.11 N in 4002 N, `p_max` to 40 Pa in 1.405 MPa and `T_max` to 1 mK. **The -2.3 K
CHT peak residual, the +0.02 K mean and the 1.17 K rms stand unchanged, and the
four-case temperature figure did not need regenerating.**

---

## 1. The defect

`src/main.cpp`, STEADY_STATE branch, before this branch:

```cpp
coupling_R_p = normalized_max_change(fields["pressure"], p_outer_prev, mesh, p_ref);
coupling_R_T = normalized_max_change(fields["temperature"], T_outer_prev, mesh, T_ref);
relax_field(fields["pressure"],    p_outer_prev, cfg.coupling_relaxation, mesh);
relax_field(fields["temperature"], T_outer_prev, cfg.coupling_relaxation, mesh);
if (coupling_R_p < cfg.coupling_tolerance && coupling_R_T < cfg.coupling_tolerance)
    break;
```

Two problems, both in what happens when the `for (outer < max_coupling)` loop
falls off its end instead of hitting `break`:

1. **It returns the last iterate.** A Picard iteration is not a descent method;
   nothing makes `max(R_p, R_T)` monotone. If the loop is cycling, the last sweep
   is an arbitrary point on the cycle.
2. **It reports nothing.** The `Coupling: N outer iters, R_p=... R_T=...` line was
   printed only under `log_outer_iters`, in cyan, with no verdict attached. A run
   that never met its own tolerance looked exactly like one that did.

The archived force-balance logs show how systematic this was. Every one of the
eight `rerun_cht_robin10` solver calls, and four of the six `outer_robin80`
calls that ran to completion, ended at exactly 200 sweeps with `R_p` between
2.6e-5 and 3.4e-4 against a `coupling_tolerance` of 1e-6. By contrast every `outer_fixed313` call converged in
26-31 sweeps with `R_p = 0` exactly, and every non-conjugate THD run in 6-16.

---

## 2. The fix

New module `src/coupling.hpp` / `src/coupling.cpp`, so the loop *control* is
separable from the loop *body* and testable without a solver.

### 2.1 Best-iterate rollback

`Coupling::Controller` ranks iterates by `max(R_p, R_T)` — the single scalar that
actually has to fall below `coupling_tolerance`. `Coupling::Snapshot` deep-copies
every registered field, and it copies the raw storage, so ghost layers and
`old_data` come along and the restored state needs no ghost exchange to be usable.

Three decisions worth defending:

| decision | why |
|---|---|
| rollback fires **only** on cap exhaustion without meeting the tolerance | a run that meets the tolerance has an acceptable last iterate by definition, so every converging case stays bit-identical to the pre-fix solver. This is the property the regression test asserts. |
| snapshot taken **before** under-relaxation | the unrelaxed end-of-sweep state is what the solvers produced and what the recorded residual describes. A relaxed blend of two sweeps is neither, and under FBNS it is not a complementarity solution of anything. |
| snapshot covers **all** fields, not just p and T | the outputs are p, theta, T, wall temperatures, wall fluxes, mu, and the integrated forces. Restoring p and T alone would deliver a state whose pressure and film content came from different sweeps. |

`coupling_keep_best_iterate = false` restores the historical behaviour exactly.

### 2.2 Non-convergence is now loud

A steady run that ends without meeting the tolerance prints, in red, **regardless
of `log_outer_iters`**:

```
Coupling: 200 outer iters, R_p=0.000269524 R_T=1.61411e-05 | NOT CONVERGED (tol=1e-06) -- hit coupling_max_iters=200
Coupling rollback: delivering outer 56 (residual 4.29799e-05, R_p=4.29799e-05 R_T=3.06711e-06) instead of
  outer 200 (residual 0.000269524); the final iterate is 6.27093x worse. The loop is not converging -
  treat this result as NON-CONVERGED.
```

`results/diagnostics.csv` gains eight columns, **appended** so header-keyed
readers (`solver_io.diagnostics_tail`) keep working:
`coupling_iters, coupling_R_p, coupling_R_T, coupling_converged,
coupling_best_iter, coupling_best_residual, coupling_restored_best,
coupling_relaxation`. In TRANSIENT they are written empty rather than as a
misleading "converged", because there is no outer loop.

The existing `Coupling: ... R_p=... R_T=...` prefix is preserved verbatim so the
case harness parser (`solver_io.coupling_health`) is unaffected.

### 2.3 Adaptive under-relaxation, and `log_coupling_history`

`coupling_relaxation_mode = ADAPTIVE` (default `FIXED`, i.e. inert):

```
omega <- max(omega_min, backoff * omega)   if max(R_p,R_T) rose
omega <- min(omega_max, growth  * omega)   otherwise
```

with `coupling_relaxation_min` 0.1, `coupling_relaxation_backoff` 0.5,
`coupling_relaxation_growth` 1.1, and `omega_max = coupling_relaxation`. It
targets the limit-cycle signature specifically: a gain slightly above 1 with a
reversing trend.

`log_coupling_history = true` prints `R_p`, `R_T`, `omega` and a `*best` marker
for every sweep. Off by default (it is up to `coupling_max_iters` lines per run);
it is what produced every trace in this report.

`coupling_max_iters`, `coupling_tolerance` and `coupling_relaxation` are now
validated (they were not), and the ADAPTIVE bounds are validated when that mode is
selected.

---

## 3. Measurement: why the loop cycles, and what stops it

All runs below: the archived `rerun_cht_robin10/fb_iter_07` configuration (the
converged FLOW_WEIGHTED equilibrium, eps = 0.603314, psi = -37.7183), G6
720 x 192, FBNS + ENERGY_EQUATION + cross-film polynomial + `CONDUCTION_2D` bush
with `bush_outer_heat_transfer = 10`, 24 pinned MPI ranks on afml-r6.

### 3.1 The `omega = 1` trace

`max(R_p, R_T)` decays geometrically for ~55 sweeps and then locks into a
sustained oscillation. It never comes within 40x of the tolerance again:

| outer | 1 | 11 | 21 | 31 | 41 | 51 | **56** | 61 | 71 | 101 | 151 | 200 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| max(R_p,R_T) | 1.6e+1 | 8.7e-2 | 1.1e-2 | 1.3e-3 | 1.7e-4 | 7.7e-5 | **4.3e-5** | 5.3e-5 | 1.6e-4 | 2.7e-4 | 2.0e-4 | 2.7e-4 |

In physical units the cycle is about **27 Pa** of pressure ripple against a
1.405 MPa peak, and about **5 mK** of temperature ripple. That is the number that
decides §4: the cycle is small.

(The earlier finding in `docs/UNFINISHED_BUSINESS.md` §5.1 sampled the FBNS mass
residual `|G|` per sweep and read the minimum at outer ~51. Measured directly on
the coupling residual it is outer 56. Same phenomenon, different probe.)

### 3.2 Relaxation sweep

Sweeps to `coupling_tolerance = 1e-6`, `coupling_max_iters = 200`:

| `coupling_relaxation` | outer sweeps | converged | final max(R_p,R_T) |
|---|---|---|---|
| **1.0** (frozen protocol) | **200** | **no** | 2.7e-4 |
| 0.95 | 48 | yes | 4.3e-7 |
| 0.9 | 53 | yes | 2.9e-7 |
| 0.8 | 64 | yes | 1.5e-7 |
| 0.7 | 76 | yes | 1.1e-7 |
| 0.5 | 113 | yes | 7.1e-8 |
| 0.3 | 197 | yes | 5.9e-8 |
| ADAPTIVE (0.1 / 0.5 / 1.1) | **52** | yes | 1.8e-7 |
| ADAPTIVE (0.1 / 0.7 / 1.1) | 70 | yes | 2.6e-8 |
| ADAPTIVE (0.1 / 0.8 / 1.05) | 69 | yes | 8.7e-7 |

Three things follow.

1. **The gain is only marginally above 1.** `omega = 0.95` is enough. This is
   exactly the diagnosis `UNFINISHED_BUSINESS.md` §5.1 predicted and the reason
   it said "relaxation before Anderson". Anderson acceleration was not needed and
   was not written.
2. **There is an interior optimum near 0.95**, and over-damping is expensive:
   0.3 costs 197 sweeps, four times 0.95.
3. **ADAPTIVE lands within 4 sweeps of the best fixed value without being told
   it.** It ends at `omega = 0.974` having dipped to a minimum of 0.5 — i.e. it
   spends one brief excursion buying stability and then returns to nearly full
   steps. The shipped defaults (backoff 0.5, growth 1.1) were chosen over
   backoff 0.7 and 0.8 by this measurement, not by taste.

Figure: `ferron_picard_convergence.png`. Raw per-sweep traces:
`picard_study/history.json`; regenerate the figure with `python plot_picard.py`.

### 3.3 Cost of the fix

Same configuration, 200 sweeps, all three ~equivalent within host noise:

| run | best-iterate tracking | per-sweep history log | wall (s) |
|---|---|---|---|
| `ovh_nofix` | off | off | 99.2 |
| `ovh_fixonly` | on | off | 95.0 |
| `base` | on | on | 85.8 |

The three differ by more than any ordering the overhead could explain (the
*slowest* is the one with the fix switched off), so **the snapshot cost is below
run-to-run scatter**. That is expected: a capture is ~6 MB of `memcpy` per rank
and only fires on a new residual minimum, which on a cycling run happens ~56
times, not 200.

Host caveat, stated because the wall times look wrong against the archive: afml-r6
was **shared** during this session — a second 24-rank job from another session was
running for part of it, putting 72 ranks on 56 cores. Absolute wall times here are
therefore not comparable with the archived batch (which recorded ~52 s for the
same run). **Outer-sweep counts are contention-free and are the number to quote**;
every wall time in this report is indicative only.

---

## 4. What the fix changed in the delivered answers: nothing beyond printed digits

This is the question the defect actually raised, and the answer is measured, not
argued.

### 4.1 Integral metrics, 12 runs, one equilibrium

| run | W [N] | force angle [deg] | p_max [Pa] | T_max [K] | h_min [m] | cav. frac. | torque [N m] |
|---|---|---|---|---|---|---|---|
| `ovh_nofix` (iterate 200, no rollback) | 4002.37011 | -90.0132613 | 1405340 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61956 |
| `base` (rolled back to outer 56) | 4002.38011 | -90.0132477 | 1405340 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61957 |
| `adapt` (converged, 52) | 4002.41011 | -90.0131774 | 1405360 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61959 |
| `r095` (converged, 48) | 4002.48010 | -90.0130553 | 1405380 | 320.154 | 6.02963e-05 | 0.264510995 | 2.61962 |
| `r08` (converged, 64) | 4002.41011 | -90.0131894 | 1405350 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61958 |
| `r05` (converged, 113) | 4002.39011 | -90.0132309 | 1405350 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61957 |
| `r03` (converged, 197) | 4002.38011 | -90.0132403 | 1405350 | 320.155 | 6.02963e-05 | 0.264510995 | 2.61957 |
| **spread, all 12 runs** | **0.11 N (2.7e-5)** | **2.1e-4 deg** | **40 Pa (2.8e-5)** | **1 mK** | **0** | **0** | **6e-5 (2.3e-5)** |

Rollback vs no rollback, in isolation: **W moves by 0.01 N** (2.5e-6 relative) and
nothing else moves at all. The force-balance tolerance is 0.1 % of 4000 N = 4 N,
i.e. 400x larger.

### 4.2 The quantity the -2.3 K claim is made from

`SUMMARY.md` scores the **bush-side film-surface temperature** `T_bearing`, not
`T_max`. Midplane profile, groove band +-12 deg excluded, over the eight runs whose
fields were extracted:

| run | `T_bearing` peak [C] | location [deg from groove] | max abs deviation from `base` over the whole profile |
|---|---|---|---|
| `base` (rollback) | 46.712 | 265.25 | — |
| `ovh_fixonly` | 46.712 | 265.25 | **0** (bit-identical) |
| `ovh_nofix` (iterate 200) | 46.712 | 265.25 | 1e-3 K |
| `adapt` | 46.712 | 265.75 | 1e-3 K |
| `r095` | 46.711 | 264.75 | 1e-3 K |
| `r08` / `r05` / `r03` | 46.712 | 265.75 / 265.25 / 265.25 | 1e-3 K |

1e-3 K is the last written digit of the VTK output, so this is an upper bound, not
a measurement of a real difference.

**Therefore:** peak 46.71 C @ 265.2 deg against a measured 49.04 C @ 258.2 deg;
mean +0.02 K; rms 1.17 K; **peak residual -2.3 K**. Unchanged. The four-case
figure `ferron_wallbc_temperature.png` was **not** regenerated because there is
nothing in it that would move by a drawable amount.

### 4.3 The one thing that does move, and why it is not a result

`T_max_theta_deg` — the theta of the hottest *cell* of the film temperature field
— takes the values 291.25 / 292.25 / 292.75 / 293.25 / 294.25 across the twelve
runs, a spread of 3 deg (6 cells at the 0.5 deg mesh spacing). The maximum is
nearly flat there, so the argmax hops between adjacent cells under perturbations
of order 1 mK. This is quantisation of an argmax, not a physical shift, and it is
a reason to prefer the fitted/peak-band metrics the case analysis already uses
over a raw argmax. It affects `T_max_theta_deg`, which is a diagnostic column; it
does not affect the `T_bearing` peak location in §4.2, which moves by at most one
cell.

---

## 5. Production reruns

All from the hash-verified commit in §6, 24 pinned ranks on afml-r6. The frozen
lanes keep `coupling_max_iters = 200`, `coupling_relaxation = 1.0`, so the
comparison against the archive is like-for-like: **the only difference from the
runs that produced the archive is the best-iterate rollback itself.**

Comparison tool: `compare_picard_rerun.py <archived_fb_log> <rerun_fb_log>`. It
pairs iterates by force-balance slot and only compares those evaluated at
identical `(eps, psi)`, because once anything moves at all the Newton driver
chooses different search points downstream.

### 5.1 `01_ferron_wallbc_fields / rerun_cht_robin10`, frozen protocol

The rollback fired on all eight calls — `coupling_restored_best = 1`, best iterate
56 to 153, final iterate 6.3x to 8.1x worse than the best — and every call is now
logged `NOT CONVERGED`.

Paired iterates (identical `eps`, `psi`), worst difference over the three that
stayed paired:

| W [N] | angle [deg] | p_max [Pa] | T_max [K] | h_min [m] | cav. frac. | torque [N m] |
|---|---|---|---|---|---|---|
| 9.1e-7 | 1.4e-6 | **0** | **0** | **0** | **0** | **0** |

The Newton search then drifts by `d_eps <= 8.6e-10` and `d_psi <= 5.3e-6 deg`.
Converged equilibrium:

| | eps | psi [deg] | Q_supply [m^3/s] |
|---|---|---|---|
| archive (iterate 200) | 0.603313784 | -37.718344987 | 3.6734606e-04 |
| rerun (best iterate) | 0.603313785 | -37.718350290 | 3.6734584e-04 |

i.e. eps agrees to 9 significant figures, psi to 5.3e-6 deg against a 0.05 deg
tolerance, and Q_supply to 6.1e-7 relative.

Delivered fields at the equilibrium call (`fb_iter_07`):

| | W [N] | angle [deg] | p_max [Pa] | `T_bearing` peak | h_min [m] | cav. frac. |
|---|---|---|---|---|---|---|
| archive (iterate 200) | 4002.37011 | -90.013261 | 1405340 | 46.7120 C @ 265.25 deg | 6.029630e-05 | 0.264510995 |
| rerun (best iterate) | 4002.38011 | -90.013254 | 1405340 | 46.7120 C @ 265.25 deg | 6.029630e-05 | 0.264510995 |
| rerun (ADAPTIVE, fully converged) | 4002.38011 | -90.013215 | 1405340 | 46.7120 C @ 265.75 deg | 6.029660e-05 | 0.264510995 |

`max |dT_bearing|` over the whole midplane profile: **1e-3 K** for both reruns
against the archive — the last written digit of the VTK output.

### 5.2 The same case, fully converged (ADAPTIVE)

Same driver, `PANCAKE_COUPLING_EXTRA="coupling_relaxation_mode = ADAPTIVE;..."`.
Every call converged:

| | sweeps per call | total sweeps | total wall [s] |
|---|---|---|---|
| frozen protocol | 200 x 8 (none converged) | 1600 | 1384 (692 of it during the host-contention episode; ~790 at the clean rate) |
| ADAPTIVE | 40, 42, 42, 43, 52, 56, 56, 57 (all converged) | **388** | **245** |

**4.1x fewer sweeps, and the answer is a converged one instead of a
non-converged one.** Equilibrium: eps = 0.603311807, psi = -37.718242158, against
the archived 0.603313784 / -37.718344987 — 2.0e-6 in eps (3.3e-6 relative) and
1.0e-4 deg in psi. Both are far below the force-balance tolerance (0.1 % in |W|,
0.05 deg), and the wall temperature profile is the same to 1e-3 K.

### 5.3 `01_ferron_conjugate / outer_fixed313`, frozen protocol

This is the control: the loop converges here (26-31 sweeps, `R_p = 0` exactly), so
no rollback fires and nothing should change.

| paired iterates | W [N] | angle [deg] | p_max, T_max, h_min, cav. frac., torque |
|---|---|---|---|
| rerun vs 24-rank archive | 7.9e-7 | 1.4e-6 | **all exactly 0** |
| **control**: 4-rank archive vs 24-rank archive, both pre-fix | 6.0e-8 | 5.7e-7 | all exactly 0 |

The residual 1e-6-level scatter in `W` and the force angle is the same order as the
scatter between two *pre-fix* archives of the same case at different rank counts,
so it is reduction-order noise, not the fix. Equilibrium reproduces the archive
exactly at every printed digit: **eps = 0.586599 / psi = -36.9839** (PINNED) and
**eps = 0.593397 / psi = -37.5036** (FLOW_WEIGHTED).

### 5.4 `01_ferron_conjugate / outer_robin80`, frozen protocol — and a new result

This variant had **no archived equilibrium** — the original batch aborted at
`fb_iter_06` on the rank-dependent KSP divergence fixed in `f0e772d`, and both
`SUMMARY.md` and `README.md` record it as abandoned. It now closes.

Against the partial archive, over the three iterates that stayed paired
(4-rank pre-KSP-fix archive vs 24-rank post-fix rerun, so this comparison
carries the rank change and the KSP fix as well as the rollback):

| W [N] | angle [deg] | p_max [Pa] | T_max [K] | h_min [m] | cav. frac. | torque [N m] |
|---|---|---|---|---|---|---|
| 4.7e-5 | 1.4e-5 | **0** | 1e-3 | **0** | **0** | **0** |

with downstream drift `d_eps <= 1.1e-8`, `d_psi <= 4.5e-5 deg`. `fb_iter_06`, the
call that used to abort, now runs.

**First force-balanced h = 80 conjugate equilibrium:** eps = 0.592143 / psi =
-36.9500 (PINNED), **eps = 0.602770 / psi = -37.7046** (FLOW_WEIGHTED), W =
4002.17 N at -90.0133 deg, p peak 1.3034 MPa gauge, h_min 60.38 um, cavitated
fraction 0.2645, torque 2.6266 N m.

Re-closed a second time under ADAPTIVE, where every call converges: eps =
0.602769 / psi = -37.7045 — the same equilibrium to six digits, from 501 outer
sweeps instead of 2178 (**4.3x fewer**; 173 s against 1097 s on an unshared host).

### 5.5 Scoring against the 16 digitised Fig. 4 points

The quantity the campaign actually reports. `T_bearing` interpolated to the 16
thermocouple angles, `solver_theta = (90 + x) mod 360`, groove band +-12 deg
excluded from the peak:

| run | mean [K] | rms [K] | max abs [K] | peak | peak residual |
|---|---|---|---|---|---|
| `rerun_cht_robin10` h = 10, **archive (iterate 200)** | +0.017 | 1.173 | 2.358 | 46.712 C @ 265.25 deg | **-2.33 K** |
| `rerun_cht_robin10` h = 10, **rerun (best iterate)** | +0.017 | 1.173 | 2.358 | 46.712 C @ 265.25 deg | **-2.33 K** |
| `rerun_cht_robin10` h = 10, **rerun (ADAPTIVE, converged)** | +0.017 | 1.173 | 2.358 | 46.712 C @ 265.75 deg | **-2.33 K** |
| `outer_fixed313` (control, converges anyway) | -2.971 | 3.581 | 6.586 | 42.488 C @ 247.75 deg | -6.55 K |
| `outer_robin80` h = 80, **force-balanced for the first time** | -0.201 | 1.214 | 2.614 | 46.450 C @ 264.75 deg | -2.59 K |
| `outer_robin80` h = 80, ADAPTIVE (fully converged) | -0.201 | 1.214 | 2.615 | 46.449 C @ 263.75 deg | -2.59 K |

The first three rows are **identical to three decimal places**, which is the
answer to the question this defect raised: the reported +0.02 K mean, 1.17 K rms
and -2.3 K peak residual are the same whether you take the iterate the cap landed
on, the best iterate the loop visited, or a fully converged solve. `outer_fixed313`
reproduces its archived -2.97 / 3.58 exactly, which is the control.

**By-product worth flagging.** The h = 80 variant, once its force balance closes,
scores mean **-0.20 K** / rms **1.21 K** — essentially as good as h = 10 (+0.02 /
1.17) and far better than the **-1.31 K / 1.76 K** recorded for it in
`01_ferron_conjugate/README.md`, which was measured on a *near-equilibrium*, never
force-balanced state because the run used to abort. So the outer-surface
coefficient matters much less than the record implies: **closing the force
balance is what supplies the heat, not the value of h** (0.26 K separates h = 80
from h = 10; 4.0 K separates either from `FIXED_TEMPERATURE`). That refines, and
does not contradict, the CORRECTION already recorded in that README.

---

## 6. Provenance

| item | value |
|---|---|
| branch | `codex/picard-best-iterate` (from `main` @ `79ddbee`) |
| solver commit run on afml-r6 | `78b310c5c8a630152971dd3cbadc0839d65002ed` (solver code identical to `282d61bd2725faa6d73554f7b2a3661ffb597ef4`; the later commit only adds the driver env hooks) |
| remote HEAD hash-verified | yes, `git rev-parse HEAD` on afml-r6 matched the local hash before every batch |
| remote tree | `~/pancake_picard` (a fresh clone; the shared `~/pancake` was **not** touched) |
| build | `cmake -B build-linux -DCMAKE_BUILD_TYPE=Release`, GCC, PETSc via pkg-config |
| gate | `ctest` 34/34 on afml-r6 (GCC) and 34/34 on Windows (MinGW/MS-MPI) |
| ranks | 24, `-x OMP_NUM_THREADS=1 -x OPENBLAS_NUM_THREADS=1` |
| study runs (12 single-point) | `r6:~/pancake/cases/01_ferron_wallbc_fields/picard_study_r6/` (331 MB) |
| wall-BC reruns | `r6:.../01_ferron_wallbc_fields/picard_rerun_frozen/` (219 MB), `picard_rerun_adaptive/` (221 MB) |
| conjugate reruns | `r6:.../01_ferron_conjugate/picard_rerun_frozen/` (602 MB), `picard_rerun_adaptive/` (302 MB) |
| local artefacts | `picard_study/*.json` (every `fb_log`, `equilibrium`, residual history and metric set), `ferron_picard_convergence.png`, `plot_picard.py`, `compare_picard_rerun.py` |

Bulk VTK stays on afml-r6 by the standing decision for this case; only provenance
and the extracted numbers are mirrored into the repo.

**Host was shared for part of this work.** A second session's 24-rank case-05
batch (`~/pancake/cases/05_etsion_extent/`) ran alongside the first rerun attempt,
putting 72 ranks on 56 cores. With MPI busy-wait that is not a 1.3x slowdown but a
**45x** one: the first `rerun_cht_robin10` call took 692 s against ~99 s once the
lanes were serialized. Two consequences: **(a)** every wall time in this report is
indicative only, and outer-sweep counts are the number to quote; **(b)** the
standing policy is one 24-rank job at a time — the lanes here were serialized by
`chain.sh` after that was measured. Contention cannot affect any solved value.

---

## 7. What is still open

1. **The best iterate is not a converged iterate.** On this case it sits at
   `max(R_p,R_T) = 4.3e-5`, 43x above `coupling_tolerance`. The rollback is an
   honesty and robustness measure. **Under-relaxation is the fix**, and whether
   the frozen case protocol should adopt `coupling_relaxation = 0.95` or
   `ADAPTIVE` is a user decision, because it changes the protocol under which
   every archived conjugate number was produced — even though, per §4, the numbers
   themselves do not move.
2. **HANDOFF D8 is untouched.** The relaxation acts on p and T while the Elrod
   primary state is theta. Under FBNS this matters less (FBNS re-solves
   `(p, theta)` from scratch every sweep and the relaxed p enters only through the
   property update and the dissipation term), but the Elrod path is still
   under-damped in the variable that matters.
3. **Why the gain crosses 1 at all** is not answered here. The suspect named in
   `PHYSICS.md` §13.1.2 is the frozen-bush Newton linearisation
   (`dT_b/dT_bar = 0`) combined with a nearly adiabatic outer surface, which
   removes the damping the film energy equation would otherwise get from the wall
   term. The measurement that would test it — refine `bush_radial_cells` and the
   theta grid and watch the critical omega — was not run.
4. **`coupling_tolerance` is a max-change criterion, not a residual.** It measures
   how much the iterate moved, not how far the state is from satisfying the
   coupled system. A run can stall with a small change and a large residual. That
   is a pre-existing property of the loop, unchanged here, and it is why "the
   residual bottomed out at 4.3e-5" should be read as "the sweep stopped moving by
   more than 4.3e-5", not as an error bound.
5. **Merge note.** This branch is based on `main` @ `79ddbee`; `main` has since
   moved to `16b7c55`, which adds a hard error for unrecognised config keys. The
   six keys introduced here (`coupling_keep_best_iterate`,
   `coupling_relaxation_mode`, `coupling_relaxation_min`,
   `coupling_relaxation_backoff`, `coupling_relaxation_growth`,
   `log_coupling_history`) are `else if` branches in the same parse chain, so
   they land ahead of that terminal `else` and need no change — but whoever
   rebases should confirm the new key list in `tests/test_config_keys.cpp`
   covers them, and that `src/config.hpp` merges cleanly (both sides edit the
   parse chain and `validate()`). No case config on this branch uses a key that
   `16b7c55` renamed: the two sets of touched `cases/` files do not overlap at
   all, which is why only the 27 files new or changed here were mirrored to
   afml-r6 rather than running `cases/sync_to_r6.sh` — a full sync from this
   branch would have pushed 106 stale `config.txt` files over `main`'s corrected
   ones. **Run the full sync from `main` after the merge.**
6. **Edge case in the ranking, stated so it is not discovered.** Because the
   metric is a per-sweep change, a genuinely *diverging* loop (residual rising
   monotonically from the first sweep) has its minimum at sweep 1, and the
   rollback would deliver the first sweep — the state closest to the uniform
   initial guess. In practice sweep 1 has by far the largest change of any sweep
   (16.1 on the trace in §3.1) so this does not arise, and such a run is loudly
   flagged NON-CONVERGED either way; but the rollback is a "best of what was
   visited" guarantee, not a correctness guarantee.
