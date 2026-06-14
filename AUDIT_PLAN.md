# AUDIT_PLAN.md — Remediation Work Packages for Implementing Agents

> **Purpose.** Executable work plan derived from `AUDIT.md` (literature-grounded
> audit), `REVIEW_VERDICT.md` (consolidated punch list), and `REVIEW.md`
> (two-agent review, findings W1–W11/M1–M4/N1–N9). Each work package (WP) is
> written so an agent with **no prior context** can implement it: physics with
> equations and citations, exact code touch-points, and test cases with
> acceptance criteria.
>
> **Relationship to PLAN.md.** PLAN.md remains the phase tracker. Note that
> PLAN.md §B.4 already *specified* the coupled mixture closure and Picard loop
> that WP-1/WP-2 below deliver — the current implementation diverged from that
> plan (gas was left volumetrically inert). WP-1/WP-2 are therefore B.4
> completed as originally designed, informed by the audit.

---

## Ground rules for every agent (read first)

1. Read `AUDIT.md` for the *why* of your WP; this file is the *how*.
2. Honor `CLAUDE.md` conventions: snake_case, minimal comments (headers
   self-documenting), modular/swappable models, **new physics must be
   disablable at the highest level** (config enum/flag, default = old
   behavior), `std::filesystem::path`, GCC/Clang/MSVC-buildable.
3. **Decoupled-limit regressions are sacred.** Every WP keeps the full existing
   `ctest` suite green. When you add physics, the "off" setting must reproduce
   the previous results bit-for-bit (or to solver tolerance — state which).
4. Test harness: add `tests/test_<name>.cpp`, register with
   `add_pancake_test(...)` in `tests/CMakeLists.txt` (runs under `mpiexec -n 2`
   by default; add an `-n 6` variant if your change touches ghost exchange, cf.
   `test_elrod_a2_mpi6`).
5. On completion: append a dated entry to `CHANGES.md`, update the relevant
   `PHYSICS.md` section with the model equations actually implemented, mark the
   WP status in `PLAN.md` Phase E, and update `README.md` if user-facing
   options changed.
6. Notation matches `PHYSICS.md`: θ = fractional film content, β = bulk
   modulus, `p_cav` = cavitation/plateau pressure, c_d = dissolved gas mass
   fraction, α_g = free-gas volume fraction, m_g = free gas mass per unit area
   (`free_gas_mass`), ρ_g = gas density, h = gap, ω = shaft speed, R = radius.

**Dependency graph (do in this order unless parallelizing):**

```
WP-4 (diagnostics)  ──────────────┐  (no deps; gives observability for all)
WP-3 (numerics honesty) ──────────┤  (no deps)
WP-5 (property models) ──┐        │
WP-2 (outer coupling) ───┼──> WP-1 (two-phase closure) ──> WP-7 (validation, gas cases)
WP-6 (θ-weighted energy/shear)────┘        WP-2 ──> WP-8 (K/C extraction)
WP-9 (starved inlet)  (independent)
WP-10 (FBNS kernel)   (strategic; after WP-1 stabilizes)
WP-7 (validation, pure-oil cases V1–V4) can start immediately.
```

---

## WP-1 — Couple released gas into the film continuity / pressure problem

**Closes:** AUDIT §A1, §A2 (the two "contradicts-the-literature" defects).
**Severity:** MAJOR. **Depends on:** WP-2 (outer loop), WP-5 (mixture μ); WP-4
strongly recommended first (conservation logs make debugging this possible).

### Physics

**The defect.** Today `ρ = ρ_l·θ` (`fluid_properties.cpp:216`) and
`reynolds.cpp` never reads `alpha_gas`/`free_gas_mass`; the outgassing source
`Δṁ_g` updates diagnostics only. Released gas can neither open void nor relieve
pressure — the defining mechanism of gaseous cavitation
(Braun & Hannon 2010, DOI 10.1243/13506501JET772; Etsion & Ludwig 1982, ASME
J. Lubr. Technol. 104(2):157–163). Published magnitude of what this suppresses:
load −3% at ε = 0.3, **+21% at ε = 0.9** (Int. J. Refrigeration 2023,
S0140700723004577).

**Target closure (Σα = 1).** Volume fractions liquid/gas must satisfy

$$\alpha_l + \alpha_g = 1,\qquad \alpha_l \equiv \theta,\qquad
\alpha_g = \frac{m_g}{\rho_g(p,T)\,h},\qquad
\rho_g = \frac{p}{R_g T}.$$

Mixture density and the film mass per unit area become

$$\bar\rho = \rho_l\,\theta + \rho_g\,\alpha_g, \qquad
M = \bar\rho\,h = \rho_l\theta h + m_g .$$

(Homogeneous-mixture precedent: Grando, Priest & Prata 2006, Tribology Letters,
DOI 10.1007/s11249-006-9027-6 — "ρ̄ = φρ_g + (1−φ)ρ_l … enter the Reynolds
equation directly"; same closure in Int. J. Refrig. 2023 and in the
solubility-based compressible Reynolds equation of ASME J. Tribology
147(2):024101, DOI 10.1115/1.4066414.)

**Saturation pressure (gaseous onset).** Already implicit in the code's Henry
law — make it explicit and *use* it:

$$p_{sat}(T, c_d) = \frac{c_d}{H(T)},\qquad
H(T) = H_{ref}\exp\!\left[E_H\!\left(\tfrac1T - \tfrac1{T_{ref}}\right)\right]$$

(inversion of `equilibrium_dissolved_gas`, `fluid_properties.cpp:76–85`). For
the shipped R290/PZ68 config: p_sat = 0.2175/2.175e-7 = 1.0 MPa. Void onset for
gas-laden oil must track p_sat, not the vaporous `p_cav` (Grando 2006: gas
released "when saturation pressure is reached … two-phase flow thereafter";
experimental: Braun & Hendricks 1984, ASLE Trans. 27(1):1–14 — cavitation
extent grows with gas solubility; cavity contents are released gas, not vapor).

**Stage 1 — constrained-θ void coupling (incremental, keeps Elrod machinery).**
Implement behind `gas_pressure_coupling = NONE | VOID_COUPLED`
(default `NONE` = current behavior):

1. **Gas displaces liquid (void ceiling).** The full-film ceiling for θ is no
   longer 1 but

   $$\theta_{full}(i,j) = 1 - \alpha_g(i,j),$$

   used in the switch function `g(θ) = [θ ≥ θ_full]`, in the post-solve snap
   (`reynolds.cpp:333–337` snaps to θ_full, not 1.0), and in the supply/BC film
   content. A full-film cell that releases gas at p < p_sat acquires α_g > 0,
   its ceiling drops, it cavitates — **gaseous rupture at p_sat emerges
   mechanistically** instead of being imposed.
2. **Local void pressure.** The cavitated-zone plateau becomes cell-local:

   $$p_{void}(i,j) = \max\!\left(p_{cav},\; \frac{m_g R_g T}{(1-\theta)\,h}\right),$$

   i.e. the gas filling the void sets the plateau; at release equilibrium it
   approaches p_sat. Pressure recovery uses
   `p = p_void + β̄ ln(θ/θ_full)` in pressurized cells, `p = p_void` in
   cavitated cells. (Experimental support for non-fixed, near-saturation cavity
   pressures: Etsion & Ludwig 1982; Braun & Hendricks 1984 — cavity pressure
   0.137 → 0.034 MPa as speed rises.)
3. **Mixture compressibility.** In pressurized cells containing gas, the
   effective bulk modulus is the harmonic mixture (gas isothermal
   compressibility = 1/p):

   $$\frac{1}{\bar\beta} = \frac{\theta}{\beta} + \frac{\alpha_g}{p},$$

   so Γ in those cells uses β̄ — gas presence makes the film soft and caps
   pressure near p_void. This *is* the capping mechanism; do not add an ad-hoc
   clamp.
4. **Mass bookkeeping.** Liquid-solution continuity keeps θ as unknown with
   release sink, gas continuity is the existing `free_gas_mass` transport with
   the release source — both per unit area:

   $$\partial_t(\rho_l\theta h) + \nabla\!\cdot\!(\text{Couette/Poiseuille fluxes}) = -\Delta\dot m_g,
   \qquad
   \partial_t m_g + \nabla\!\cdot\!(m_g \mathbf{u}) = +\Delta\dot m_g,$$

   with the **existing** finite-rate kinetics (keep them — peer-reviewed
   precedent: Hao & Gu 2014, Tribology Int. 78:14–26, DOI
   10.1016/j.triboint.2014.04.028; experimental motivation: Etsion & Ludwig's
   p/p_s > 1 release/resorption loop):

   $$\Delta\dot m_g = k_{dg}\,(c_d - c_{eq}(p,T))\,\rho_l\,h\,\max(\theta,\theta_{min}).$$

5. **Forces/outputs** use ρ̄ and the WP-5 mixture viscosity.

**Stage 2 (strategic, pairs with WP-10).** Full homogeneous compressible
mixture EOS ρ̄(p, T, w) with equilibrium flash + kinetic relaxation, replacing
the hard switch entirely (Grando 2006; Bayada & Chupin 2013, ASME J. Tribology
135(4):041702). Do not start Stage 2 until Stage 1 is validated.

> **Design caution for the implementing agent.** Stage 1's assembly
> (θ_full ceiling + local p_void + β̄) is a synthesis consistent with the cited
> closures, not a transcription of a single published scheme. The 0-D and 1-D
> tests below are therefore *mandatory acceptance gates*, and the equilibrium
> limit must reproduce Grando's behavior (pressure plateau at p_sat; classical
> Reynolds-BC solution recovered when gas is disabled).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Add `enum class GasPressureCoupling { NONE, VOID_COUPLED }`, key `gas_pressure_coupling` (default NONE). Validation: `VOID_COUPLED` requires `GAS_CAVITATION_MIXTURE` and TRANSIENT (until WP-2's steady loop lands). Add `gas_constant_specific` if not present (R_g). |
| `fluid_properties.cpp/hpp` | Export `saturation_pressure(c_d, T, cfg)` (Henry inversion incl. TABLE path); set `rho = ρ_l θ + ρ_g α_g` when coupling on; keep `rho = ρ_l θ` when NONE. |
| `reynolds.cpp` (`solve_elrod`) | Per-cell `theta_full` field from α_g (lagged one outer iteration); switch, snap, Γ (with β̄), and Couette flux base density use ρ̄/θ_full; pressure recovery via local `p_void`. Release sink in the θ equation source. |
| `equation_of_state.hpp` | Generalize `switch_function(theta, theta_full)`, `pressure_from_theta(theta, theta_full, p_void, beta_bar)`; keep old signatures as θ_full=1, p_void=p_cav wrappers. |
| `gas_transport.cpp` | No structural change (transport already exists); ensure the release source is consistent with the θ-equation sink (one function, used by both). |
| `main.cpp` | Order inside WP-2's outer loop: properties → Reynolds(θ, p) → gas transport → exchange → properties → … (already close; the difference is iteration, WP-2). |
| `io.cpp` | Output `p_void`, `theta_full`, `p_sat` as optional fields. |

### Test cases (`tests/test_gas_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Decoupled regression**: `gas_pressure_coupling = NONE` on `config_r290_pz68_quick.txt` | Fields identical to current master to solver tolerance (≤ 1e-10 relative) |
| 2 | **Σα invariant**: random release/transport steps | `θ + α_g ≤ 1 + 1e-12` every cell, every step |
| 3 | **0-D closed-cell equilibrium**: single cell, no flux, p < p_sat initially, kinetics on | p → p_sat asymptotically; total propane mass (c_d ρ_l θ h + m_g) constant to 1e-12; equilibrium c_d = H p_sat |
| 4 | **0-D resorption**: start gassy at p > p_sat | m_g → 0, c_d → c_eq(p); no negative masses |
| 5 | **1-D pressure cap**: converging-diverging channel, supersaturated feed | In the gassy zone, p ∈ [p_void, p_sat·(1+5%)]; no cell sustains p ≪ p_sat with θ = θ_full while supersaturated |
| 6 | **Mass conservation (periodic θ, no-flux z)** | Global liquid and total-propane residual (WP-4 logger) ≤ 1e-10 per step |
| 7 | **Load-shift direction** (journal, ε = 0.9, R290 config) | Load with coupling ON > load with NONE (direction per IJR 2023 +21%; log the ratio, no hard % gate) |
| 8 | **Vaporous limit**: c_d = 0, no gas | Identical to plain Elrod (p_void ≡ p_cav, θ_full ≡ 1) |

---

## WP-2 — Outer coupling (Picard) iteration; real steady state; steady+gas guard

**Closes:** AUDIT §C4 (REVIEW N1, N7). **Severity:** MAJOR. **Depends on:** none
(but unblocks WP-1, WP-6 fidelity, WP-8).

### Physics

A segregated one-sweep "solution" is not a fixed point of the coupled system;
published transient cavitation codes iterate the coupling inside the step —
Ausas, Jai & Buscaglia 2009 (ASME J. Tribology 131(3):031702, DOI
10.1115/1.3142903) update even the journal motion "within the same relaxation
process." Fixed-point (Picard) iteration with under-relaxation:

$$\phi^{(k+1)} \leftarrow \omega_\phi\,\phi^{new} + (1-\omega_\phi)\,\phi^{(k)},
\qquad \phi \in \{\mu, \rho, T, c_d, m_g\},\ \omega_\phi \in (0,1].$$

Convergence on normalized residuals each outer iteration k:

$$R_p = \frac{\max_{ij}|p^{(k)}-p^{(k-1)}|}{p_{ref}},\quad
R_T = \frac{\max_{ij}|T^{(k)}-T^{(k-1)}|}{T_{ref}},\quad
R_c = \max_{ij}|c_d^{(k)}-c_d^{(k-1)}|,$$

stop when all < `coupling_tolerance` (default 1e-6) or `coupling_max_iters`
(default 30; warn on hitting the cap). **Report R_p in pressure space**, not
θ-space — β amplifies θ-noise by ~1e9 (AUDIT §C2/N4), so a θ-residual hides
pressure-level error.

`STEADY_STATE` = run this loop to convergence once (this is the operating
point). Steady + gas: currently `gas_dt = 0` silently freezes the gas model
(`main.cpp:243`, `fluid_properties.cpp:222`) — replace with either a hard
config error or `steady_gas_model = FROZEN | EQUILIBRIUM` where EQUILIBRIUM
sets c_d = c_eq(p,T) and α_g from the equilibrium flash each outer iteration.

### Code implementation

| File | Change |
|------|--------|
| `main.cpp` | Extract the per-step block (`main.cpp:220–260`) into `solve_coupled_step(fields, sys, mesh, cfg, bearing_state, dt)`; wrap in the outer loop with relaxation + residual log (rank-0, `Utils::log`). Steady mode: iterate to tolerance instead of `step == 0` single pass. |
| `config.hpp` | Keys: `coupling_max_iters`, `coupling_tolerance`, `coupling_relaxation` (one factor; per-field override later if needed), `steady_gas_model`. `validate()`: error on `GAS_CAVITATION_MIXTURE && STEADY_STATE && steady_gas_model == FROZEN` unless user sets an explicit `allow_frozen_gas_steady = true`. |
| `journal_motion.cpp` | (With WP-1/WP-8 in mind) optionally update bearing state inside the outer loop (Ausas-style) behind `motion_coupling = LAGGED | IN_LOOP`, default LAGGED. |
| `io.cpp`/logging | Per-step line: outer iters, final R_p/R_T/R_c, Elrod inner iters, flag flips (WP-4 shares this). |

### Test cases (`tests/test_coupling.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Constant-property limit**: `fluid_property_model = CONSTANT`, ISOTHERMAL | Converges in exactly 1 outer iteration; results identical to master |
| 2 | **Monotone convergence**: EMPIRICAL_CORRELATION μ(T) with strong dissipation | R_p, R_T decrease monotonically (allowing first-iter transient); converged < max iters |
| 3 | **Ordering independence**: permute sweep order (Reynolds↔energy) | Converged steady fields agree to `coupling_tolerance` (this is the definition of a fixed point — fails on master) |
| 4 | **Steady ≡ long transient**: steady solve vs transient run to t→∞, same config | max|p_steady − p_transient|/p_ref < 5·tolerance |
| 5 | **Guard**: GAS_CAVITATION_MIXTURE + STEADY_STATE + FROZEN | `validate()` throws/errors with actionable message |
| 6 | **Relaxation robustness**: ω = 0.5 and 1.0 on test 2 | Same converged answer |

---

## WP-3 — Numerics honesty: convection schemes, dead selectors, capacity-term gap

**Closes:** AUDIT §C1, §C5 (REVIEW N3, N6 + audit's new finding). **Severity:**
medium, cheap. **Depends on:** none.

### Physics

- First-order upwind *inside the cavitated zone* is standard — even FBNS uses
  it (Woloszynski et al. 2015, Tribology Letters, DOI 10.1007/s11249-015-0487-4:
  "in the cavitation region … a first-order upwind formula was employed").
  Upwinding *everywhere* is not: the canonical scheme is **type differencing** —
  central in the full film, upwind across/inside cavitation (Vijayaraghavan &
  Keith 1989, Tribology Trans. 32(2):225–233, DOI 10.1080/10402008908981882;
  San Andrés Notes 06). False diffusion of full upwind ≈ |u|Δx/2 (Patankar
  1980, ch. 5) lands exactly on the rupture/reformation/thermal/gas fronts.
- TVD face reconstruction (already coded, dead in `fvm.cpp:13–27,118–162`):
  $$\phi_f = \phi_U + \tfrac12\,\psi(r)\,(\phi_D-\phi_U),\quad
  r = \frac{\phi_U-\phi_{UU}}{\phi_D-\phi_U},\quad
  \psi_{vanLeer}(r)=\frac{r+|r|}{1+|r|}.$$
- Transient capacity: $\partial(\rho_0\theta h)/\partial t$ must be assembled
  for **all** transient runs; currently gated on `MOVING_BEARING`
  (`reynolds.cpp:261–262`), so a fixed-bearing transient (e.g. an
  `omega_ramp_time` startup) is quasi-steady in θ — film-history physics
  silently absent.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | Key `theta_convection_scheme = UPWIND \| TVD_VANLEER \| TVD_MINMOD \| TYPE_DIFFERENCING` (default UPWIND for back-compat), similarly `thermal_convection_scheme`, `gas_convection_scheme`. **Delete** `pressure_time_method`/`temperature_time_method` members + GUI combos; parser accepts the old keys with a one-line warning ("ignored: field solves are backward-Euler"). |
| `reynolds.cpp` | Pass configured scheme to `FVM::divergence` (`:259`); TYPE_DIFFERENCING = per-face: central if `g=1` on both sides, upwind otherwise. Change the transient gate at `:261` to `cfg.solution_mode != STEADY_STATE` for the `ddt_weighted` capacity; squeeze source remains conditional on `dh_dt` presence. Align the Gumbel Couette discretization (`:146–147`) with the Elrod path's scheme so the A.1↔Gumbel limit compares like-for-like. |
| `energy.cpp`, `gas_transport.cpp` | Pass configured schemes (`energy.cpp:344`, `gas_transport.cpp:241,260`). |
| `gui_win32.cpp` | Remove dead time-method combos; add scheme dropdowns. |
| `fvm.cpp` | None expected — limiters/deferred correction exist; verify ghost-width (n_ghost = 2) suffices for the UU stencil at partition boundaries (it does; add the MPI test below anyway). |

### Test cases (extend `tests/test_fvm.cpp`, new `tests/test_schemes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **TVD boundedness**: 1-D step advection of θ | No new extrema (θ ∈ [min, max] of initial data); sharper front than upwind (L1 error strictly smaller) |
| 2 | **Order of accuracy**: smooth 1-D profile | Upwind ~1st order, TVD ≥ ~1.5 observed order on refinement |
| 3 | **A.1↔Gumbel match**: matched schemes both paths | Residual strictly below the current 1.2% (record value in test) |
| 4 | **Capacity-term fix**: fixed bearing, ω step at t=0 | θ relaxes over finite time (exponential-like), not instantaneously; steady final state matches steady solve |
| 5 | **Regression**: all schemes default UPWIND | Existing ctest suite untouched |
| 6 | **MPI**: TVD run, 2 vs 6 ranks | Rank-count independence of fields to 1e-12 (cf. `test_elrod_a2_mpi6`) |
| 7 | **Config back-compat**: old config with `pressure_time_method = CRANK_NICOLSON` | Parses, warns, runs |

Also deliverable: `validation/grid_convergence/` script running 60×20 / 120×40 /
240×80 on the journal case, reporting cavitation extent, h_min, load (archived
table — AUDIT §C1 requirement; feeds WP-7).

---

## WP-4 — Conservation, regime, and convergence diagnostics

**Closes:** AUDIT Part V item 7, §B1 logging (REVIEW N9, M3, N8-logging).
**Severity:** low effort / high value. **Depends on:** none. **Do first.**

### Physics

Per-step global balances on the physical domain (Ω = bearing surface,
dA = R dθ dz):

$$M_l = \int_\Omega \rho_l\,\theta\,h\,dA,\qquad
M_{gas,tot} = \int_\Omega \big(c_d\,\rho_l\,\theta\,h + m_g\big)\,dA,$$

$$r_l^{n+1} = \frac{M_l^{n+1}-M_l^{n} - \Delta t\,(\Phi_{in}-\Phi_{out}+S_{exch})}{\max(M_l^{n},\epsilon)},$$

likewise r_gas (exchange terms cancel in M_gas,tot — only boundary fluxes
remain). Boundary fluxes from the assembled face fluxes (Couette + Poiseuille at
z-boundaries and inlet penalties), not re-derived. Clamp/snap operations
(`reynolds.cpp:324–337`, gas clamps) must be **accounted as explicit source
terms** in the balance so their mass effect is visible (REVIEW N9: clamps are
not mass-neutral). Target per PLAN A.2: O(1e-12) in closed domains.

Regime guards at startup (config validate + rank-0 log):

- Couette Reynolds number and Taylor criterion (laminar validity):
  $$Re_c = \frac{\rho\,\omega R\,c}{\mu},\qquad
  Ta = Re_c\sqrt{c/R} \lesssim 41,$$
  warn above (Taylor-vortex onset; standard journal-bearing criterion, cf.
  Hamrock, *Fundamentals of Fluid Film Lubrication*; San Andrés Notes).
- Saturation context (AUDIT §B1): log `p_sat(T_init, c_d_init)`, ratio
  `bc_z_*_val / p_sat`, and `p_cav / p_sat`; warn when boundary oil enters
  saturated (ratio ≥ 1) so the current R290 configuration is visible.
- Per-step: Elrod outer iterations, flag flips, (after WP-2) coupling residuals.

### Code implementation

| File | Change |
|------|--------|
| `src/diagnostics.hpp/cpp` (new) | `Diagnostics::mass_balance(fields, mesh, cfg, dt)` returning {M_l, M_gas, r_l, r_gas, clamp_mass}; MPI_Allreduce inside; rank-0 CSV append (`results/<run>/diagnostics.csv`) + log line. |
| `reynolds.cpp` | Accumulate clamp/snap mass deltas into a fields-resident scalar (or returned struct) for the balance. |
| `config.hpp::validate()` | Re/Ta computation + warnings; p_sat ratio warnings; key `diagnostics_interval` (default 1 = every step). |
| `main.cpp` | Call after each step; write CSV. |

### Test cases (`tests/test_diagnostics.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Closed-domain conservation**: periodic θ, zero-gradient z, no inlets, transient | \|r_l\| ≤ 1e-12 per step (the PLAN goal, now actually checked) |
| 2 | **Open-domain accounting**: Dirichlet z ends | Balance closes: residual after flux accounting ≤ 1e-10 |
| 3 | **Clamp visibility**: initialize θ < θ_min somewhere | clamp_mass ≠ 0 reported, balance still closes |
| 4 | **Gas total**: WP-1 test 3 config | r_gas ≤ 1e-12 (exchange cancels) |
| 5 | **Re guard**: synthetic config with μ tiny | validate() warns (Ta > 41) |
| 6 | **Saturation log**: R290 config | startup log contains p_sat = 1.0e6 and ratio = 1.0 warning |

---

## WP-5 — Property models: mixture viscosity with gas-limit asymptote; solubility tables

**Closes:** AUDIT §A3, §A4 (REVIEW W3, W4). **Severity:** medium. **Depends
on:** none (WP-1 consumes the result).

### Physics

**Free-gas (bubbly/foam) viscosity.** Einstein μ = μ_l(1 + 2.5α_g) is a dilute
result (α ≲ 0.1) and is directionally wrong as α_g → 1 (gas must *thin* the
film). Replace with selectable models, all satisfying μ̄(0) = μ_l:

- `DUKLER_VOID` (linear void-weighted; the Grando 2006 density-side analogue):
  $$\bar\mu = \alpha_g\,\mu_g + (1-\alpha_g)\,\mu_l \;\to\; \mu_g \text{ as } \alpha_g\to1.$$
- `MCADAMS_QUALITY` (homogeneous two-phase standard; quality
  x_g = m_g/(m_g + ρ_lθh)):
  $$\frac{1}{\bar\mu} = \frac{x_g}{\mu_g} + \frac{1-x_g}{\mu_l}.$$
- `EINSTEIN_DILUTE` (keep, for back-compat) — but `validate()` errors when
  selected with `gas_alpha_max > 0.6` (Krieger–Dougherty packing bound):
  $$\bar\mu = \mu_l\,(1-\alpha_g/\alpha_{max})^{-2.5\,\alpha_{max}}\ \text{(KD, optional fourth model)}.$$

(Grando 2006 quote: "μ̄ = χμ_g + (1−χ)μ_l" — quality-weighted linear; offer it
as `LINEAR_QUALITY` if exact Grando replication is wanted for WP-7 V5.)

**Dissolved-gas (solution) viscosity.** The validated form for R290/oil pairs
is the Eyring-MTSM activity model (Int. J. Refrigeration 2023,
S0140700723002451: "viscosities of R-290/oil mixtures could be well represented
by the Eyring-MTSM model"; measured 303–348 K on PAG 68 / POE 75 / PVE 68).
Do **not** implement the activity model — feed measured isotherms through the
existing `TABLE` path (`viscosity_table`, `solubility_table`,
`density_table`), and document linear Henry + `exp(a_c c_d)` as a local tangent
valid near the (10 bar, 40 °C) calibration point only. Grando's reference case
(2 bar, w_sat = 7.13 %) shows how far the shipped 21.75 % case extrapolates.

> **Data task (flag to user, do not fabricate):** populating the R290/PZ68
> tables requires the actual isotherm data (e.g. the IJR 2023 dataset or
> supplier data for PZ68). Ship the *plumbing* + a clearly-labeled placeholder
> table; mark the case "table data pending" in config comments and README.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `enum class GasMixtureViscosityModel { EINSTEIN_DILUTE, DUKLER_VOID, MCADAMS_QUALITY, KRIEGER_DOUGHERTY, LINEAR_QUALITY }`, key `gas_mixture_viscosity_model` (default EINSTEIN_DILUTE = current behavior); `mu_gas` config key; validation pairing with `gas_alpha_max` as above. |
| `fluid_properties.cpp` | Replace hardcoded Einstein at `:202–207` with the dispatch; α_g clamp from `gas_alpha_max`, not the hardcoded 0.99. |
| `config_r290_pz68*.txt` | Set `gas_alpha_max = 0.6`, select `MCADAMS_QUALITY`; comment the validity statement. |

### Test cases (extend `tests/test_fluid_properties.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Asymptotes**: each model | μ̄(α=0) = μ_l exactly; DUKLER/MCADAMS/KD: μ̄ → μ_g (or finite KD limit) as α→α_max; monotone decreasing in α (except EINSTEIN which is documented-increasing) |
| 2 | **Einstein guard**: EINSTEIN + α_max = 1.0 | validate() error |
| 3 | **Regression**: default EINSTEIN, α_max as shipped | Existing test_fluid_properties values unchanged |
| 4 | **Table round-trip**: synthetic 3-point viscosity/solubility tables | Interpolation exact at nodes, linear between, clamped outside |
| 5 | **Calibration point**: R290 case (existing test) | Still passes with new defaults documented |

---

## WP-6 — θ-weighted energy capacity, enthalpy flux, and cavitated-zone shear

**Closes:** AUDIT §A5 inconsistency (REVIEW N5, part of W6). **Severity:**
medium. **Depends on:** none (interacts with WP-2 for converged THD).

### Physics

In a striated cavitated zone only the liquid fraction θ carries thermal mass,
enthalpy, and (predominantly) shear. Replace, behind
`cavitated_film_model = FULL_FILM | STRIATED` (default FULL_FILM = current):

$$\text{capacity: } \rho c_p\,\theta\,h\,\frac{\partial T}{\partial t},\qquad
\text{flux: } \rho c_p\,\theta_f\,h_f\,\bar u\,T,\qquad
\text{Couette shear: } \tau_\theta = \theta\,\frac{\mu\,\omega R}{h}\ (\text{cavitated cells}),$$

with full-film cells (θ = 1) unchanged. This restores internal consistency: the
gas transport already θ-weights its capacity (`gas_transport.cpp:31,193`).
Friction torque/power-loss are first-order THD deliverables (groove-design
practice reports power loss + h_min: ASME J. Tribology 144(12):121801, 2022) —
the FULL_FILM model overstates them by up to ~1/θ where θ → θ_min = 1e-6.
Cross-film profile and conjugate conduction remain Phase D scope; record that
the validated 1983–84 baseline includes them (Ferron 1983, ASME 105(3):422;
Lund & Tonnesen 1984, ASME 106(2):237–244 — 4th-order cross-film polynomial +
sleeve conduction) so temperature claims stay gated on WP-7 V6.

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `cavitated_film_model` enum + key; validation: STRIATED requires ELROD_ADAMS. |
| `energy.cpp` | Capacity `:301` ×θ; flux loops `:308–340` ×θ_face (face value: upwind θ consistent with flux direction); wall-loss term: decide and document whether Q_w scales with θ (liquid contact) — recommended ×θ for consistency. |
| `reynolds.cpp` | `:519` and `:543`: Couette term ×(g_f + (1−g_f)·θ) i.e. full shear in full film, θ-scaled in cavitated cells. |

### Test cases (extend `tests/test_energy.cpp`, `tests/test_forces.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Full-film limit**: θ ≡ 1, STRIATED | Bit-identical to FULL_FILM |
| 2 | **Scaling**: uniform cavitated strip θ = 0.5 | Friction torque contribution of that strip = 0.5× FULL_FILM value; energy capacity likewise |
| 3 | **Energy balance**: WP-4 style closed thermal balance with θ-weighting | Heat generated − wall loss = dE/dt to 1e-10 |
| 4 | **Direction**: journal case with cavitation, STRIATED vs FULL_FILM | Power loss strictly lower with STRIATED; T_max not higher |
| 5 | **Regression**: default FULL_FILM | Suite unchanged |

---

## WP-7 — Validation campaign (archived, quantitative)

**Closes:** AUDIT Part V (REVIEW W11). **Depends on:** V1–V4 none; V5 needs
WP-1; V6 needs WP-2+WP-6.

### Structure

`validation/<case_id>/` containing `config*.txt`, `reference/` (digitized CSV +
provenance note: paper, figure/table number, digitization method), `run.ps1` +
`run.sh`, `compare.py` (tolerance-gated, exit code), `README.md`. Each case also
gets a coarse "smoke" variant wired into ctest (`add_test(NAME validation_<id>_smoke ...)`).

### Cases, references, acceptance

| ID | Case | Reference | Acceptance gate |
|----|------|-----------|-----------------|
| V1 | 1-D JFO slider (rupture + reformation), analytic | Elrod & Adams 1975; PLAN A.2 derivation | Rupture/reformation locations within 1 cell at 200 cells; p_peak ≤ 1 %; θ-profile L2 ≤ 1 % |
| V2 | Oscillatory squeeze film (exact solution) + dynamically loaded JB | Ausas, Jai & Buscaglia 2009, ASME JT 131(3):031702 (paper ships reference source code — use it for code-to-code) | p(t) L2 ≤ 2 % vs exact; JB orbit qualitative match + mass residual ≤ 1e-10 |
| V3 | Grooved journal bearing, flooded + (post-WP-9) starved | Vijayaraghavan & Keith 1989: Trib. Trans. 32(2):225–233 (algorithm cases, code-to-code) and Wear 134(2):377–397 (grooved, vs experiment) | Cavitation extent ±1 cell @ published grid; load ±5 % |
| V4 | LCP analytic suite (pocketed slider, squeeze damper) | Giacopini et al. 2010, ASME JT 132(4):041702 | p-profile L2 ≤ 2 % on the 1-D cases |
| V5 | Refrigerant equilibrium limit | Grando, Priest & Prata 2006 (Trib. Lett.) | With gas coupling ON and equilibrium kinetics (k_dg → large): recover classical Reynolds-BC solution for the moderately loaded case (their validation); plateau at p_sat |
| V6 | THD test bearing | Ferron, Frêne & Boncompain 1983, ASME 105(3):422; Lund & Tonnesen 1984, ASME 106(2):237 | Phase 1: qualitative T(θ) shape + documented gaps (no conjugate conduction). Hard gates deferred to Phase D; record deltas. |
| V7 | Grid convergence + conservation history | WP-3/WP-4 outputs | Monotone convergence of extent/load/h_min; archived plots |

**Citation hygiene for the implementing agent:** the V&K 1989 Trib. Trans.
paper is code-to-code only (its own abstract: comparisons "with the results
obtained using Elrod's algorithm"); the experimental comparison is the Wear
companion. Do not cite the former as experimental validation.

---

## WP-8 — Linearized dynamic coefficients (K_ij, C_ij) and whirl margin

**Closes:** AUDIT §D1 (REVIEW W9b). **Depends on:** WP-2 (needs a converged
operating point).

### Physics

Small-perturbation method (Lund 1987, ASME J. Tribology 109(1):37–41, DOI
10.1115/1.3261324 — the canonical reference; applied to gas-oil refrigerant
bearings via the perturbed compressible Reynolds equation in ASME JT
147(2):024101, DOI 10.1115/1.4066414):

$$K_{ij} = -\left.\frac{\partial F_i}{\partial x_j}\right|_{eq}
\approx -\frac{F_i(x_j^{eq}+\Delta x_j)-F_i(x_j^{eq}-\Delta x_j)}{2\,\Delta x_j},
\qquad
C_{ij} = -\left.\frac{\partial F_i}{\partial \dot x_j}\right|_{eq},$$

i, j ∈ {x, y}; force = fluid force on the journal/bearing per the existing
convention. Perturbation sizes: Δx = ε_K·c with ε_K ~ 1e-2 (linearity check:
halving Δx changes K by < 1 %); velocity perturbations Δẋ = ε_C·c·ω applied as
prescribed `dh/dt` with frozen position. Each force evaluation is a **converged
WP-2 operating point** (transient effects off; squeeze term from the prescribed
ẋ only).

Stability deliverable: for a rigid rotor of mass m per bearing, the
characteristic equation `det(m s² I + C s + K) = 0`; report eigenvalues, onset
speed (Re(s) = 0), and whirl frequency ratio Im(s)/ω (≈ 0.5 for plain bearings
— sanity anchor; San Andrés Notes 05). Lund's validity caveat (small motion
about equilibrium) goes in the output header.

### Code implementation

| File | Change |
|------|--------|
| `src/dynamic_coefficients.hpp/cpp` (new) | `compute_dynamic_coefficients(fields, sys, mesh, cfg, state)` → struct {K[2][2], C[2][2], whirl_ratio, onset_mass_or_speed}; loops ±perturbations calling the WP-2 steady driver. |
| `config.hpp` | `compute_dynamic_coefficients` (bool), `perturbation_position_fraction`, `perturbation_velocity_fraction`. |
| `main.cpp` | Post-run hook when enabled; CSV + log output. |
| `io.cpp` | Append K/C block to a summary file. |

### Test cases (`tests/test_dynamic_coefficients.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Linearity**: Δx and Δx/2 | K_ij relative change < 1 % |
| 2 | **Sign structure** (plain journal, moderate ε) | K_xx, K_yy > 0; K_xy·K_yx < 0 (cross-coupling signs per Lund 1987 / San Andrés Notes 05) |
| 3 | **Symmetry of C** (aligned, no gas) | C_xy ≈ C_yx within 5 % |
| 4 | **Consistency vs transient**: small free decay orbit vs linear model prediction | Frequency/decay within 10 % |
| 5 | **Whirl anchor**: light-load plain bearing | whirl ratio ∈ [0.45, 0.55] |

---

## WP-9 — Starved / mass-flux inlet option

**Closes:** AUDIT §B2 (REVIEW W8). **Depends on:** none.

### Physics

Starved feed = the supply cannot maintain a full film at the inlet; prescribe
inlet film fraction θ_in < 1 (or feed flow). JFO reformation downstream is
already handled by Elrod machinery. Precedent: starved cases inside the Elrod
framework since Vijayaraghavan & Keith 1989 (Wear 134(2):377–397 — "both
flooded and starved inlet conditions … various degrees of starvation"); groove
geometry as a design input per ASME JT 144(12):121801 (2022).

Two modes on `InletConfig`:

- `STARVED_THETA`: penalty-pin θ = θ_in (θ_in ∈ (θ_min, 1)); cavitated inlet
  cell (g = 0), Couette flux carries ρ_l θ_in h u into the domain.
- `MASS_FLUX`: prescribe feed ṁ [kg/s]; convert to equivalent θ_in over the
  inlet footprint from the local Couette flux:
  θ_in = ṁ / (ρ_l (ωR/2) h L_z,inlet), clamped to (θ_min, θ_supply(p)).

### Code implementation

| File | Change |
|------|--------|
| `config.hpp` | `InletConfig::Type` extends with `STARVED_THETA`, `MASS_FLUX` (parse `inlet_N_mode`, `inlet_N_theta`, `inlet_N_mdot`); validation: starved modes require ELROD_ADAMS. |
| `reynolds.cpp` (`apply_inlet_conditions`, `:60–91`) | Branch: pressure mode unchanged; starved modes pin θ_in via the same penalty; MASS_FLUX computes θ_in per cell from local h. |
| `energy.cpp`/`gas_transport.cpp` | Inlet thermal/composition contracts treat starved inlets like supply cells (existing reservoir-composition path). |
| `gui_win32.cpp` | Inlet mode selector + value field. |

### Test cases (`tests/test_inlets.cpp` extension)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Flooded equivalence**: θ_in = θ_supply(p_supply) | Matches pressure-inlet result to 1e-10 |
| 2 | **Monotonicity**: θ_in ∈ {1.0, 0.8, 0.5} | Load strictly decreasing; cavitation extent strictly increasing |
| 3 | **Flux fidelity** (MASS_FLUX) | Integrated inlet mass flux = prescribed ṁ within 1 % |
| 4 | **Conservation**: starved case under WP-4 logger | Balance closes ≤ 1e-10 |

---

## WP-10 — Strategic: complementarity (FBNS) cavitation kernel

**Closes:** AUDIT §C2, §C3 (β-stiffness, switch-iteration fragility; REVIEW
N4 + S3 caveat). **Depends on:** stabilize WP-1 first (formulation choice
interacts). **This retires `snap_tol` and the β-conditioning class of bugs.**

### Physics

JFO mass conservation as a complementarity problem (Giacopini et al. 2010, ASME
JT 132(4):041702, DOI 10.1115/1.4002215):

$$\bar p \ge 0,\qquad \vartheta \ge 0,\qquad \bar p\,\vartheta = 0,
\qquad \bar p = p - p_{cav},\ \ \vartheta = 1-\theta,$$

reformulated as a smooth unconstrained system via the Fischer–Burmeister
function (Woloszynski, Podsiadlo & Stachowiak 2015, Tribology Letters, DOI
10.1007/s11249-015-0487-4):

$$\phi_{FB}(\bar p_{ij}, \vartheta_{ij}) = \bar p_{ij} + \vartheta_{ij}
- \sqrt{\bar p_{ij}^2 + \vartheta_{ij}^2} = 0,$$

solved with Newton + Schur-complement reduction (FBNS) coupled to the
discretized Reynolds residual R(p, θ) = 0. Literature-reported performance: ~2
orders of magnitude faster than Elrod-Adams/LCP-pivoting baselines,
mesh-independent Newton counts (900 → 360 000 DOF), upwind θ-flux in the
cavitated region retained. Motivation for leaving the β-switch: solution
"strongly depend[s]" on β, β dominates the discrete system's stability
(eigenvalue analysis), softening tolerable only ~2 orders (3–40 % error beyond)
— ASME JT 139(3):031703 (2017); practitioner β-reduction practice per San
Andrés Notes 06. No β appears in the FBNS formulation; physical β returns only
through an optional compressible full-film EOS.

### Code implementation

| File | Change |
|------|--------|
| `src/reynolds_fbns.cpp/hpp` (new) | `Reynolds::solve_fbns(...)` parallel to `solve_elrod`; Newton loop with analytic Jacobian blocks (Reynolds w.r.t. p and θ; φ_FB diagonal blocks); reuse `LinearSystem`/PETSc (KSP inside Newton, or SNES if available). |
| `config.hpp` | `cavitation_model = GUMBEL \| ELROD_ADAMS \| FBNS`; `fbns_max_newton`, `fbns_tolerance`. |
| `main.cpp` | Dispatch. WP-1's mixture terms enter through ρ̄, p_void analogues — coordinate with whoever owns WP-1 Stage 2. |

### Test cases (`tests/test_fbns.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Equivalence**: V1 slider + journal case vs ELROD_ADAMS | p, θ agree within discretization tolerance (≤ 0.5 %) |
| 2 | **Complementarity**: converged fields | max(p̄_ij·ϑ_ij) ≤ 1e-10·p_ref |
| 3 | **Mesh-independent iterations** | Newton counts within ±2 across 60×20 → 240×80 |
| 4 | **No snap needed** | Flag-chatter test (the scenario behind `snap_tol`, reynolds.cpp:328–337 comment) converges without snapping |
| 5 | **Timing** | ≥ 5× faster than ELROD_ADAMS outer loop on the 240×80 case (log; literature says ~100×, don't gate on that) |

---

## Phase E status table (mirror in PLAN.md; update as you complete)

| WP | Title | Status | Owner | Tests green | Docs updated |
|----|-------|--------|-------|-------------|--------------|
| WP-4 | Diagnostics & guards | **DONE** 2026-06-11 | Claude | `test_diagnostics` (6 cases) | CHANGES, PHYSICS §6/§13.2.1/§15, README |
| WP-3 | Numerics honesty | **DONE** 2026-06-11 | Claude | `test_schemes` + `test_schemes_mpi6` | CHANGES, PHYSICS §5.5/§7, README |
| WP-5 | Property models | **DONE** 2026-06-11 (isotherm table data pending from user) | Claude | `test_fluid_properties` extended | CHANGES, PHYSICS §13.2.1, README, configs |
| WP-2 | Outer coupling / steady | TODO | — | — | — |
| WP-1 | Two-phase closure | TODO | — | — | — |
| WP-6 | θ-weighted energy/shear | TODO | — | — | — |
| WP-7 | Validation campaign | TODO (V7 grid-convergence harness seeded: `validation/grid_convergence/`) | — | — | — |
| WP-9 | Starved inlet | TODO | — | — | — |
| WP-8 | K/C coefficients | TODO | — | — | — |
| WP-10 | FBNS kernel | TODO | — | — | — |

Implementation deltas vs the WP specs (2026-06-11):
- WP-4: balances are formed against the assembled discrete operators (storage
  from `theta` vs `theta.old`), so no Tracker state is needed; a `linear_rtol`
  config key was added because the 1e-12 closed-domain gate is only reachable
  when the Krylov tolerance is tightened (default stays 1e-8).
- WP-3: the Gumbel alignment factorizes `rho*h` into central `rho0*h` times
  scheme-interpolated `theta = rho/rho0` — upwinding the whole product had
  *increased* the A.1 mismatch to 8.4%; the factorized form reaches ~0.01%.
  `test_elrod_a1`/`test_elrod_a2` now declare STEADY_STATE (their quasi-steady
  intent) since the capacity term is no longer silently absent.
- WP-5: `LINEAR_QUALITY` (Grando replication) shipped alongside the four main
  models; the mixture-viscosity dispatch lives in
  `FluidProperties::gas_mixture_viscosity` for direct unit testing.

Also in scope, fold into the first WP that touches the file (tracked here so
nothing is lost):

- **Docs reconciliation (AUDIT §B1/B3):** PHYSICS.md §6 ("open ends p = p_cav")
  vs shipped submerged configs; PHYSICS.md §13.2.1 R290 table (1e5/1e5) vs
  shipped (3e4/1e6); add the constant-cavity-pressure caveat (Etsion & Ludwig).
  → do with WP-4 (it adds the p_sat logging that makes the choice explicit).

## WP-11 — Axial-boundary flux consistency (single source of truth)

**Closes:** AUDIT §B4 (line-level BC review, 2026-06-11). **Severity:** medium
(BC-3 live when cavitation reaches the ends; BC-1/2/4 latent — no shipped
config uses `INLET_OUTLET`). **Depends on:** none; coordinate with WP-1 (the
shared flux function is where `p_void` will later replace `p_cav`).

### Physics

The three solvers currently assemble the *same* axial boundary three different
ways. Required invariant: **one liquid mass flux per boundary face, used by
Elrod assembly, energy convection, gas convection, and the diagnostics
balance**, with direction-dependent temperature/composition layered on top.

Per face (south shown; north mirrored), with `p_bc = elrod_boundary_pressure(bc_val)`:

$$\dot m''_{bc} =
\begin{cases}
\dfrac{\rho_0\,\theta_P\,h^3}{12\mu}\,\dfrac{p_{bc}-p_P}{\Delta z/2} & \theta_P \ge 1 \ \text{(full film — EOS chain rule } \beta\,d\theta=\theta\,dp\text{)}\\[2ex]
\max\!\left(0,\ \dfrac{\rho_0\,h^3}{12\mu}\,\dfrac{p_{bc}-p_{cav}}{\Delta z/2}\right) & \theta_P < 1 \ \text{(cavitated — physical reformation rate)}
\end{cases}$$

The cavitated branch replaces the current ungated `Γ_base` Dirichlet link,
which overestimates re-flooding by ≈ `β(1−θ_P)/(p_bc−p_cav)` (~10²–10³×)
because the EOS `p = p_cav + β ln θ` is invalid where p ≡ p_cav. After WP-1,
`p_cav` in the cavitated branch becomes `p_void(i,j)`.

`INLET_OUTLET` semantics unified: for the **mass** equation it behaves as
DIRICHLET (bidirectional, like the Gumbel path already does); direction
dependence applies only to *what the flux carries* — inflow brings
`temperature_reference` (RESERVOIR mode) and `dissolved_gas_initial` / zero
free gas; outflow convects cell values (energy/gas already do exactly this).
The lagged `fields["pressure"] < bc_val` gate (stale by one timestep, frozen
across outer iterations) is deleted; direction comes from the sign of the
shared flux, evaluated with the current iterate.

Cavitated-cell reformation inflow must also carry reservoir enthalpy and
composition (today the g-gate zeroes both, so re-flooding liquid silently
adopts the cell's T and c_d).

### Code implementation

| File | Change |
|------|--------|
| `src/boundary_flux.hpp/cpp` (new) | `axial_boundary_mass_flux(fields, mesh, cfg, side, i)` implementing the two-branch formula; plus `axial_boundary_volume_flux` (free gas) and a struct return {mass_flux, is_inflow}. |
| `reynolds.cpp` (`solve_elrod`) | Full-film boundary cells: keep the Dirichlet link (equivalent). Cavitated boundary cells: drop the `Γ_base` link; add the reformation flux as an explicit θ-equation source `+ṁ''_{bc}·RΔθ/(ρ_0 h ...)` consistent with the assembled units. Delete the `INLET_OUTLET` lagged gate (`:388–414`); INLET_OUTLET assembles like DIRICHLET. |
| `energy.cpp` | `south/north_boundary_flux` call the shared function (×`c_p` face value) instead of the local g-gated expression; RESERVOIR/OPEN logic unchanged. |
| `gas_transport.cpp` | `south/north_boundary_mass_flux`/`_volume_flux` call the shared function; composition logic unchanged. |
| `diagnostics.cpp` | Balance uses the same shared function (removes any residual reimplementation). |
| `config.hpp` | Document NEUMANN = sealed/symmetry plane in the key comment + PHYSICS.md §6. |

### Test cases (`tests/test_boundary_fluxes.cpp`)

| # | Test | Acceptance |
|---|------|-----------|
| 1 | **Triple consistency**: random full-film states | Elrod-assembled boundary flux, energy flux ÷ ρc_p-face, gas mass flux agree cell-by-cell to 1e-12 (by construction once shared) |
| 2 | **Cavitated reformation rate**: 1-D z-strip, cavitated interior, submerged end | Re-flood mass rate = `ρ₀h³(p_bc−p_cav)/(12μ·Δz/2)` analytic within 1 %; assert the old `Γ_base` link would give > 10× (regression demonstrating the fix) |
| 3 | **Outflow continuity**: interior p > p_bc with INLET_OUTLET | Liquid mass leaves; WP-4 balance closes ≤ 1e-10; energy/gas convect the same direction in the same step |
| 4 | **Reformation carries reservoir state**: cavitated end cell re-floods (RESERVOIR mode) | Cell c_d moves toward `dissolved_gas_initial`, T toward `temperature_reference` — not frozen |
| 5 | **Direction-flip robustness**: p_bc oscillating about interior pressure | Flux continuous through zero; no chatter; conservation holds |
| 6 | **Full-film regression**: shipped configs (ends full-film) | Results identical to current master to solver tolerance |
| 7 | **Gumbel parity**: INLET_OUTLET vs DIRICHLET in Gumbel | Identical (documents BC-4 resolution) |

## Reference list (canonical for all WPs — copy citations from here)

Same as `AUDIT.md` appendix. Key DOIs: Elrod-Adams numerics — 10.1080/10402008908981882,
10.1115/1.3142903, 10.1115/1.4002215, 10.1007/s11249-015-0487-4; refrigerant
two-phase — 10.1007/s11249-006-9027-6, S0140700723004577, 10.1115/1.4066414,
10.1016/j.triboint.2014.04.028; experiments — Etsion & Ludwig (ASME JLT
104(2):157), Braun & Hendricks (10.1080/05698198408981539); THD — Ferron 1983
(ASME 105(3):422), Lund & Tonnesen (10.1115/1.3260891); dynamics — Lund 1987
(10.1115/1.3261324); review — Braun & Hannon (10.1243/13506501JET772).
