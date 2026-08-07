# Validation status

This is the authoritative ledger of claims for the current validation campaign.
`VALIDATION_PLAN.md` defines planned work, `INDEX.md` locates evidence, and the
case reports contain detail. Root `VALIDATION_PLAN.md`,
`validation/VALIDATION_REPORT.md`, and `docs/UNFINISHED_BUSINESS.md` are dated
snapshots and do not override this file.

Snapshot: 2026-08-07, solver commit `ad53384` on
`codex/validation-reconciliation`.

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
| Campaign resolution | G6 720x192 and `dt = 2e-5 s` are accepted operational settings for the current campaign | Cell-centred `p_max` did not demonstrate grid independence; the later Ferron coupled THD/CHT model was not covered by the original grid pair |
| Raimondi-Boyd locus | Numerical benchmark agrees within about 2.5% in Sommerfeld number and 0.54 deg in attitude over the transcribed points | This verifies the hydrodynamic locus against tabulated theory, not experiment |
| Ferron pressure | At the scored 2000 rpm / 4 kN condition, peak and location agree closely with the digitised experiment | One operating point and unresolved `p_max` discretization uncertainty do not validate every pressure field |
| Ferron temperature | The final conjugate comparison gives mean +0.02 K and RMS 1.17 K over 16 bush-temperature points; Picard correction does not move the quoted values | Peak remains about 2.3 K low; user verdict and coupled-model grid check remain open |
| Absolute load capacity | No load-dependent residual remains in the Ferron collapse over the tested 5x load range | Uniform clearance bias and uniform solver load bias are exactly degenerate; absolute load accuracy is unresolved |
| Brito speed trend | No constant clearance rescale removes the remaining speed-dependent discrepancy | Thermal clearance growth and under-predicted temperature are plausible but not uniquely identified |
| Etsion cavity extent | FBNS partly reproduces rupture and the ambient-pressure trend; JFO performs poorly | Every scored run has `alpha_gas = 0`; this is not validation of gas release, void fraction, or cavity pressure |
| Gaseous cavitation | No current experimental validation claim | Spatial cavity pressure and finite-rate desorption are outside the converged model |

## Language controls

Allowed:

- “FBNS is numerically verified for the tested 1-D free-boundary problem.”
- “Ferron pressure and conjugate bush-temperature comparisons are promising at
  the scored condition.”
- “The Etsion extent comparison identifies where the equilibrium plateau model
  succeeds and fails.”

Not allowed yet:

- “The CFD code is fully validated.”
- “The load calculation is exonerated.”
- “The gaseous-cavitation model is experimentally validated.”
- “G6 makes every reported local pressure or coupled THD field grid-independent.”

## Next decision gates

1. Obtain the user's verdict on the final Ferron temperature comparison and the
   Etsion extent comparison without broadening either claim.
2. Run a focused coupled THD/CHT grid check for the quantities actually quoted.
3. Break the absolute-load degeneracy with independent clearance or film-thickness
   evidence; more rescaling of the same proximity-probe data cannot do it.
4. Specify and verify a finite-rate gas-release/cavity-pressure formulation before
   attempting case 05b or claiming gaseous-cavitation validation.
