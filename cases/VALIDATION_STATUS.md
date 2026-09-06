# Validation status

This is the authoritative ledger of claims for the current validation campaign.
`VALIDATION_PLAN.md` defines planned work, `INDEX.md` locates evidence, and the
case reports contain detail. Root `VALIDATION_PLAN.md`,
`validation/VALIDATION_REPORT.md`, and `docs/UNFINISHED_BUSINESS.md` are dated
snapshots and do not override this file.

Snapshot: 2026-08-07, solver commit `df0802b` on
`codex/ferron-cht-grid-check`.

## 2026-09-06 audit addendum

The August campaign evidence below is retained as a dated snapshot. Current
code/physics repairs and later branch evidence are reconciled in
[the September audit](../docs/PHYSICS_AUDIT_2026-09.md).

- Full native Windows and Ubuntu WSL CTest suites: **40/40 pass** on the repaired
  source, including conservation on MPI-1/2/4 and a forced solver failure.
- The existing exact 1-D FBNS cases were rerun: both fronts within one cell on
  120/240/480/960 grids; finest peak error +0.0701%, load error -0.0361%.
- Synthetic oil-air transport/reaction agrees across 1/2/4 ranks and conserves
  gas to less than 2e-11 relative step residual. It also exhibits about 2.7e-4
  excess gap volume after the split update. Numerical conservation therefore
  does **not** validate its pressure-volume model.
- Later local D-series Ferron campaign artifacts report passing grid and dt
  comparisons for their revised setups. They remain in the separate dirty
  campaign checkout and explicitly deny temperature-validation eligibility.
  This does not replace the earlier failed grid check with a passing verdict.
- No experimental gas-release or vaporous/gaseous mixture validation is added.
  R410A is unsupported and now rejected explicitly. The old scored Etsion cases
  still have zero released gas and cannot validate finite-rate cavitation.

Reproduction: [audit case and evidence](90_diagnostics/physics_audit/README.md).

## 2026-09-06 paired-result addendum

[Paired runs and five figures](90_diagnostics/physics_comparison/SUMMARY.md)
compare `b9dcaf3` with `2bc4d3c` on identical inputs. The tested pure-oil steady
Ferron and optical-rig fields are unchanged at output precision. This does not
change their experimental-validation status. A closed reaction probe directly
shows the old 11.548% gas loss and repaired conservation to roundoff.

**The extended gas transient fails timestep acceptance for pressure/load.**
The 1 ms synthetic air case has alternating pressure and force at 10, 5 and
2.5 us; refinement does not suppress the oscillation. All inner solves
converge, but transient outer coupling is not iterated. Maximum gas fraction
is comparatively insensitive while excess gap volume remains nonzero.
Do not quote the terminal force difference as a physical load change or infer
experimental readiness from the earlier 0.1 ms conservation check.

## Reproducible baseline

| Host | Toolchain | Result |
|---|---|---|
| Windows development PC | MinGW GCC 16.1, MS-MPI, distributed PETSc, Release | fresh configure + full rebuild; 35/35 CTest tests pass |
| Ubuntu workstation `afml-r6` | GCC 11.4, OpenMPI, PETSc, Release | isolated checkout of the same commit; full build; 35/35 CTest tests pass |

The local and remote commit hashes and transferred Git-bundle SHA-256 were
checked before the Ubuntu build. Existing remote campaign trees were not used or
modified.

## Current claims

| Area | Defensible claim | Limit |
|---|---|---|
| Solver regression tests | Cross-platform baseline passes, including MPI-4/6, FBNS, mass residual, energy, cross-film, gas coupling, and steady coupling | Unit/regression tests verify implementations, not physical adequacy |
| FBNS free boundary | Verified against the exact 1-D JFO solution: rupture and reformation within one cell, integrated errors converge, and per-cell mass residual is at solver precision | Published 2-D contour verification remains open |
| Campaign resolution | G6 720x192 and `dt = 2e-5 s` remain operational settings for the original campaign | The focused Ferron CHT check fails its G6-to-G7 mean/RMS temperature limits (0.111/0.107 K versus 0.100 K); G6 is not certified for those quoted thermal metrics |
| Raimondi-Boyd locus | Numerical benchmark agrees within about 2.5% in Sommerfeld number and 0.54 deg in attitude over the transcribed points | This verifies the hydrodynamic locus against tabulated theory, not experiment |
| Ferron pressure | Corrected G6 gives 1.353 MPa gauge at 199.75 deg: +7.31% and 1.25 deg from the digitised peak; G6-to-G7 `p_max` changes 0.078% | One operating point does not validate every pressure field; the digitisation has about 0.1 MPa peak scatter |
| Ferron temperature | No validation claim. Corrected G6 gives mean +3.82 K, RMS 4.02 K and peak 50.38 C; G7 gives +3.71 K, 3.91 K and 50.27 C | The old +0.02 K / 1.17 K result used a diagnostic source that grew as `1/dz`; corrected G6-to-G7 mean/RMS changes also fail the registered grid limits |
| Absolute load capacity | No load-dependent residual remains in the Ferron collapse over the tested 5x load range | Uniform clearance bias and uniform solver load bias are exactly degenerate; absolute load accuracy is unresolved |
| Brito speed trend | No constant clearance rescale removes the remaining speed-dependent discrepancy | Thermal clearance growth and under-predicted temperature are plausible but not uniquely identified |
| Etsion cavity extent | FBNS partly reproduces rupture and the ambient-pressure trend; JFO performs poorly | Every scored run has `alpha_gas = 0`; this is not validation of gas release, void fraction, or cavity pressure |
| Gaseous cavitation | No current experimental validation claim | Spatial cavity pressure and finite-rate desorption are outside the converged model |

## Language controls

Allowed:

- "FBNS is numerically verified for the tested 1-D free-boundary problem."
- "Ferron pressure agrees within 7.3% in peak and 1.3 deg in location at the
  scored condition; the conjugate temperature comparison currently fails."
- "The Etsion extent comparison identifies where the equilibrium plateau model
  succeeds and fails."

Not allowed yet:

- "The CFD code is fully validated."
- "The load calculation is exonerated."
- "The gaseous-cavitation model is experimentally validated."
- "G6 makes every reported local pressure or coupled THD field grid-independent."

## Next decision gates

1. Replace or justify the inferred fresh-oil mixing flow using Ferron's measured
   inlet-flow evidence before tuning wall or cavitation models to temperature.
2. Extend the corrected CHT grid sequence before certifying a production grid.
3. Break the absolute-load degeneracy with independent clearance or film-thickness
   evidence; more rescaling of the same proximity-probe data cannot do it.
4. Specify and verify a finite-rate gas-release/cavity-pressure formulation before
   attempting case 05b or claiming gaseous-cavitation validation.
