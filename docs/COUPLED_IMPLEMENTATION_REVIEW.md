# Cavitation and heating: what was corrected and what the evidence means

This extends the [history audit](PHYSICS_AUDIT_2026-09.md) and
[FW68D review](FW68D_CAVITATION_REVIEW.md). The target remains FW68D with
air or R410A, reservoirs up to 20 bar absolute, reservoir temperatures
20–80°C, and refrigerant charging near 10 bar. Film temperature is now
predicted rather than prescribed when thermal mode is selected.

The result is a conservative, homogeneous thin-film solver with a native
R410A free-fluid EOS. It is not yet a validated prediction of either FW68D
experiment. The distinction matters: correct conservation and thermodynamic
identities can be demonstrated without knowing the material coefficients,
whereas experimental pressure, load, cavity extent and heating cannot.

## 1. The original pressure and gas updates could describe different films

The pressure solver used one gas volume, after which transport and phase
exchange changed the dissolved and free inventories. Pressure was not
reconciled with the new gas volume before the step was accepted. Consequently,
a small pressure residual could coexist with liquid plus gas occupying more
than the geometrical gap. In the archived 1-ms synthetic oil-air case,
halving the timestep did not produce a stable load history. This was a
coupling failure, not evidence of physical bearing oscillation.

The earlier audit repaired an additional inventory error: a closed reaction
test lost 11.548% of its total dissolved-plus-free gas. It now conserves that
inventory to roundoff. Conservation alone did not repair the pressure/volume
inconsistency, which is why this implementation replaces the split path for
the target mixture simulations.

COUPLED solves pressure and three mass inventories simultaneously. Oil mass
and total dissolved/free species mass are conserved. The dissolved balance
supplies finite-rate phase exchange, and gas density comes from the accepted
pressure and temperature. The volume equation requires that liquid and the
gas-containing void fill precisely the available gap. No subsequent gas cap
deletes mass or hides excess volume. See the
[complete residual and flux equations](COUPLED_CAVITATION.md).

The existing FBNS solver remains useful as a gas-free JFO reference. Its old
segregated TWO_PHASE combination must not be used as the replacement for
this simultaneous mixture solve. A successful regression of legacy behavior
does not establish that behavior's physical validity.

## 2. Air release pressure is not the oil vapor pressure

Vaporous cavitation occurs when the liquid cannot sustain pressure above
its vapor-equilibrium limit. Gaseous cavitation instead depends on dissolved
gas inventory, solubility, time available for transfer and bubble compression.
Oil equilibrated with a gas at 10 bar does not impose a universal 10-bar
pressure floor throughout the bearing. Raising an atmospheric outlet to that
bubble point changes the experiment rather than modeling its cavities.

The repaired boundary validation rejects incompatible requested pressure
floors. In COUPLED the oil vapor pressure is $p_v(T)$, while the free-species
partial pressure is $p-p_v(T)$. Gas-containing cavities can therefore have
spatially varying pressure. At zero free-species inventory, the equations
recover the negligible-vapor-mass JFO limit.

Release and absorption have separate rates. Setting absorption to zero
does not freeze bubble volume: the bubble still compresses through its EOS.
Conversely, compressing a bubble must not be mistaken for immediate dissolution.
Published oil measurements demonstrate rapid release and delayed absorption,
but they do not supply universal rate constants for FW68D
(Reinke et al., *[Fluids](https://doi.org/10.3390/fluids7050159)*, 2022).

The optional local-equilibrium formulation is explicitly a limiting model.
It permits direct comparison with equilibrium solubility treatments, such as
the air-solubility journal-bearing model
(Li et al., *[J. Tribol.](https://doi.org/10.1115/1.4006702)*, 2012).
Agreement with that assumption would not demonstrate correct transient
desorption or absorption times.

## 3. Reservoirs must follow the actual flow direction

Both axial ends are pressure boundaries. At each face, the pressure solution
determines whether the end supplies fluid or receives it. Incoming fluid uses
that reservoir's temperature, dissolved fraction and entrained free-gas mass
fraction. Outgoing fluid carries the local film state. The corrected boundary
mobility also extrapolates outgoing properties: an arbitrary cold downstream
reservoir must not impose a fictitious cold, viscous resistance without
reverse flow.

A groove or circular hole is an internal pressure patch. Its required net
mass source follows from storage and neighboring fluxes. Positive source
uses its specified supply composition and enthalpy; negative source drains
the local mixture. A draining inlet is therefore not a fresh-oil reset or a
fixed-temperature cell. An end-fed bearing needs no internal inlet.

These conditions are verified with independent one-dimensional pressure,
species and thermal profiles in both directions, plus inflow/drainage tests
for both patch shapes. The current groove extends across the axial length;
finite-length grooves and feed-pipe jets require additional geometry/modeling.

## 4. Heating changes the phase state as well as viscosity

The earlier mixture workflow could append a liquid-only temperature update
after the phase solve. That does not conserve the energy of the accepted
oil/dissolved/free mixture. It also misses the work needed to compress a
bubble and the energy difference between dissolved and free refrigerant.

Temperature is now a fifth coupled unknown. Stored energy is the sum of
component mass times component internal energy. The balance includes
Couette and pressure-flow dissipation, compression work, in-plane conduction,
wall heat removal and reservoir enthalpy. A specified positive desorption
energy produces cooling during release and heating during absorption through
the changing inventories. Adding a second phase-transfer heat source would
double-count that contribution.

The tests cover exact uniform shear heating, wall cooling, independently
calculated adiabatic bubble compression, closed desorption energy, heated
pure-vapor complementarity and reservoir enthalpy. Repeating a nonlinear
solve at the same physical time does not heat or release gas a second time.
The [energy derivation](COUPLED_ENERGY.md) defines every input and unit.

This predicts a common film-average temperature. It does not yet resolve
temperature variation across the clearance or conduction in the shaft and
bush. A measured wall temperature cannot be compared directly to this fluid
temperature. Oil-vapor mass and evaporation enthalpy are neglected in the
JFO limit; strong thermal evaporation requires an expanded phase model.

## 5. R410A cannot be represented by inherited propane or ideal-air constants

R410A now has its own fixed-blend Helmholtz density and caloric EOS
(Lemmon, *[Int. J. Thermophys.](https://doi.org/10.1023/A:1025048800563)*, 2003).
At 40°C and 20 bar, ideal-gas density is about 26.4% below the real-fluid
value; the corresponding fixed-mass gas volume is about 35.8% too large.
That directly distorts available liquid volume, compressibility and load.
The implementation is checked at 15 density states and four caloric states
against independently generated CoolProp evaluations of the published EOS.

The 20-bar reservoir limit must not cap the EOS at 20 bar. Hydrodynamic film
pressure can be higher, so the free-fluid API permits higher pressure within
its stated thermodynamic domain. Conversely, at 20°C a free R410A phase
above its dew pressure cannot simply remain a stable gas. Gas-free dissolved
refrigerant is a different state and can exist in pressurized oil.

The solver rejects a resolved subcritical free gas outside the stable gas
domain. Equilibrium dissolution can remove it if the solution has capacity;
finite-rate transport may require a separate condensed refrigerant phase.
Replacing the absorption rate with an arbitrary large value to hide that
failure would change the physics. Blend fractionation and a refrigerant-rich
liquid phase are not included by the fixed-blend EOS.

The remote 20°C run also exposed a distinct implementation defect: a tiny
outgoing gas-mass perturbation used to form a numerical derivative was
tested as if it were an accepted bubble. The resulting domain failure
prevented construction of the Jacobian. The repaired check permits outgoing
trial continuation and retains the final physical-domain test. A new
regression reproduces the old failure and passes after the repair. With
the corrected code, all four synthetic 20-bar equilibrium cases at
20/40/60/80°C complete; finite-rate 20/40°C cases still expose the missing
condensed-phase treatment.

## 6. Saturated property tables do not describe arbitrary dissolved charge

A saturated pressure/temperature viscosity table implicitly links composition
to equilibrium pressure. During rapid pressure changes, the transported
liquid can be under- or supersaturated. Using that table as if it described
the actual liquid can apply the wrong viscosity before any gas has transferred.

The new composition table is indexed by actual dissolved mass fraction and
temperature and interpolates log viscosity. It rejects values outside the
supplied data range. The source-backed
[FVC68B–R410A dataset](../cases/ref_data/fvc68b_r410a_1999_viscosity.md) verifies
this import/interpolation path. It is not FW68D and spans only 23–65°C.
Its measured viscosity dependence illustrates why a nominal ISO viscosity
grade or a pure-refrigerant EOS cannot replace mixture-property measurements.

The recent R410A journal-bearing study uses FVS50, also a different oil,
and does not expose all coefficients needed to reproduce its material model
(Li et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2026.106916)*, 2026).
No unavailable FW68D coefficient has been inferred from its plots.

## 7. Verification and experimental validation are separate gates

The current native Windows, local WSL and remote workstation builds pass all
51 regression tests. Final bearing results are in the
[workstation report](../cases/90_diagnostics/coupled_validation/SUMMARY.md).
The [campaign record](../cases/90_diagnostics/coupled_validation/README.md)
preserves inputs, source/binary hashes, solver return codes and machine data.

On the workstation the coupled 1-ms oil-air load is 598.6536, 598.6482 and
598.6455 N at 10, 5 and 2.5 µs, while the paired segregated loads remain
erratic. With heating enabled, the 48 × 8 oil-air case reaches 45.80°C and
its load drops 7.12% relative to its isothermal counterpart. The last spatial
refinement, 192 × 32 to 384 × 64, changes the thermal load by 0.234%.
These values establish numerical behavior for the disclosed synthetic inputs.

The independent one-dimensional JFO solution tests pressure and moving
full-film/cavity boundaries in the zero-gas limit. A separate dynamic vapor
benchmark is available in the literature and remains a future full-case
reproduction (Brewe, *[J. Tribol.](https://doi.org/10.1115/1.3261288)*, 1986).

The Etsion/Ludwig experiment tests a different regime: gas cavities with two
flooded pressure ends and no internal inlet. Its selected photograph gives a
235° cavity extent. The final coupled runs predict substantially shorter
equilibrium cavities (165° on both tested grids at a 0.001 liquid-deficit
threshold); disabling absorption leaves weak gas around the entire
circumference at low detection thresholds. Neither result reproduces the
photograph. This mismatch must remain visible rather than being removed by
choosing a convenient plotting threshold
(Etsion and Ludwig, *[J. Lubr. Technol.](https://doi.org/10.1115/1.3253174)*, 1982).

For FW68D, the outstanding material inputs are solubility versus pressure and
temperature, viscosity versus actual dissolved concentration and temperature,
dissolved partial volume, consistent caloric/desorption data, preparation
temperature and thermal reservoir/wall conditions. Final experimental claims
also need mesh/timestep independence and comparison of like observables.
