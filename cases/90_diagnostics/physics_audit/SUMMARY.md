# Physics audit results, 2026-09-06

All 40 CTest cases pass on native Windows and Ubuntu WSL. The tests include
closed oil/gas inventory conservation, changing carrier capacity, reservoir
composition, gas-volume error detection, failed-solver exit status, and exact
transient heating/rigid-body timing. [Test logs and reproduction](README.md).

| Exact 1-D FBNS grid | Peak-pressure error | Load error | Rupture / reformation error (cells) |
|---|---|---|---|
| 120 | +0.540% | +0.695% | 0.922 / 0.664 |
| 240 | +0.280% | +0.273% | 0.844 / 0.329 |
| 480 | +0.142% | +0.066% | 0.689 / 0.658 |
| 960 | +0.070% | -0.036% | 0.378 / 0.315 |

The synthetic oil-air case advances ten intervals of 10 microseconds. Across
1/2/4 MPI ranks, the normalized gas conservation residual stays below 2e-11;
saved final fields agree within 6e-10 normalized difference. The final maximum
gas fraction is 0.00170326. These are synthetic constants, not oil-air validation.

**The gas-volume closure fails:** the maximum raw occupied-volume excess is
0.000269731 of the gap (0.0270%). The alpha/EOS clipping defect is zero in this
mild case, so the measured excess comes from split pressure/phase updates rather
than the alpha cap. The diagnosis remains visible even when mass conservation
and Newton convergence pass.

R410A and separate vapor/noncondensable inventories remain unsupported.
[Required physics and validation work](../../../docs/PHYSICS_AUDIT_2026-09.md).
