# Case 01-THD-crossfilm - wall-BC comparison with cross-film T(y), FBNS

Ferron Fig. 4 point (2000 rpm, 4 kN), branch codex/fbns-crossfilm (e487c25):
cross-film polynomial T(y) active, mu(T(y)) inside the generalized (Dowson)
Reynolds mobilities (Poiseuille rho_l*F2, Couette thickness 2F1/F0). Frozen
parameter set and force-balance protocol of cases/01_ferron_thd. 33 archived
solver calls; every number in thd_analysis.json.

Wall-BC mapping (the conjugate model has no fixed-T option): "fixed 40 C" =
large-conductance limit onto thermal_ambient_temperature, realised to 1 mK
(probe_cfg/probe.json). Variants:
A walls_adiabatic (shaft_ambient_heat_transfer = 0, bush_outer_heat_transfer = 0),
B bush40_Jadiabatic, C bush40_J40.

## Results (flow pass, W closure <= 0.003%)

| quantity | A adiabatic | B bush40/J-ad | C bush40/J40 | depth-avg ref |
|---|---|---|---|---|
| eps | 0.6086 | 0.5860 | 0.5646 | 0.5738 |
| phi (classical) deg | 52.45 | 52.65 | 53.86 | 53.68 |
| T_max depth-avg C | 48.41 | 43.70 | 41.81 | 43.92 |
| bush-side T peak C | 50.25 @ 347.7 deg | 40.00 (=BC) | 40.00 (=BC) | n/a |
| T_shaft C | 45.13 | 44.05 | 40.00 | n/a |
| p peak gauge MPa | 1.3097 (+3.9%) | 1.2876 (+2.2%) | 1.2643 (+0.3%) | 1.2702 (+0.8%) |
| h_min um | 59.49 | 62.93 | 66.18 | 64.78 |
| torque N m | 2.596 | 2.786 | 2.890 | 2.808 |

Experiment: p peak 1.2602 MPa @ 201.0 deg; T peak 49.04 C @ 258.2 deg.

Bush-side film-surface T vs the 16 digitised Fig. 4 points:
- A: mean +0.18 K, rms 2.07 K; by segment: -0.96 K (0-140 deg), -0.47 K
  (140-270 deg), +3.80 K (>270 deg). Magnitude of the rise captured; the
  adiabatic film keeps heating to the groove, so the peak lands at 347.7 deg
  instead of the measured 258.2 deg - an adiabatic wall has no mechanism for
  the measured mid-arc turndown (that is bush conduction physics).
- B, C: bush-side surface = BC = flat 40.00 C, structurally uninformative in
  this comparison; depth-averaged values carry their information (too cold).

Deviations and caveats (full list in the run report): groove-edge one-cell
spike excluded (+-12 deg band) from peaks; adiabatic bush side can sit below
feed T (order-4 closure property, pre-existing); PHYSICS.md 13.1.2 collocation
caveats apply (reduced Peclet 10-30 binds the saturation); Q_supply still the
self-consistent two-pass value; nothing tuned, no fitted shifts.

Figures: ferron_crossfilm_temperature.png (bush-side T, three variants +
depth-averaged reference + experiment; model quantity = bush-side film surface
temperature, closest available to Ferron's 0.5 mm sub-surface thermocouples;
mapping solver_theta = (90 + x) mod 360, no fitted shift),
ferron_crossfilm_pressure.png.

Stage-1 code notes (commit e487c25): shared src/reynolds_coefficients.hpp is
now the single definition of all face coefficients (4 drifted diagnostics
copies found and unified); FBNS+cross-film vs Elrod+cross-film agree to 0.48 K
on the shared full-film arc; cross-film off is byte-identical to 6dacec0;
VOID_COUPLED mass_residual drops 9.4x with solved fields byte-identical (90%
of the previously reported imbalance was the diagnostics mismatch, not the
solver). Pre-existing unpoliced failure: test_fbns at -n 4 (Elrod film-content
KSP diverges; ctest only runs -n 2/-n 6) - same BCGS/block-Jacobi robustness
family as the temperature-KSP defect already ticketed.
