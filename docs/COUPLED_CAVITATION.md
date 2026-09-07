# Coupled conservation and gas volume

This document describes the implemented isothermal and thermal model, its tests, and its
limits. It does not certify FW68D material properties or experimental accuracy.

## What changed

The segregated TWO_PHASE path first solved pressure using the previous gas
volume, then transported and released gas. The subsequent gas update changed
the liquid inventory and void volume without resolving pressure. A small
pressure residual therefore did not imply a conserved, volume-compatible
accepted state. Lowering the pressure datum to a representative oil vapor
pressure made this failure particularly severe.

`cavitation_model = COUPLED` instead solves pressure, oil mass, dissolved
species mass and free species mass simultaneously, with temperature as a fifth
unknown when `temperature_model = ENERGY_EQUATION`. Both axial pressure ends
are bidirectional reservoirs. A groove or circular hole constrains pressure
and exchanges the composition selected by the actual flow direction.

## State and equations

All inventories below are per bearing surface area, in kg/m²:

\[
M_o,\qquad M_d,\qquad M_g,\qquad
M_l=M_o+M_d,\quad c_d={M_d\over M_l},\quad
\theta_l={M_l\over\rho_l(c_d,p,T)h}.
\]

Oil and the dissolved/free species have separate conserved inventories:

\[
\partial_tM_o+\nabla_s\cdot F_o=0,
\qquad
\partial_t(M_d+M_g)+\nabla_s\cdot(F_d+F_g)=0.
\]

The dissolved balance closes finite-rate exchange:

\[
\partial_tM_d+\nabla_s\cdot F_d=-\Gamma,
\qquad
\Gamma=k(1-c_{\rm eq})(M_d-M_{d,\rm eq}),
\]

\[
M_{d,\rm eq}=\min\left(M_d+M_g,
 {c_{\rm eq}\over1-c_{\rm eq}}M_o\right).
\]

The release and absorption rates can differ. Setting absorption to zero
preserves free gas during compression; it does not prevent the bubble volume
from shrinking. The inventory-limited equilibrium mass prevents absorption
from demanding gas that is absent. Backward Euler integrates this model;
the frozen-pressure exponential used by the old split model is not reused.

Alternatively, `coupled_phase_equilibrium = true` replaces the kinetic
equation with complementarity between free mass and undersaturation:

\[
c_{\rm eq}-c_d\ge0,\quad M_g\ge0,\quad
(c_{\rm eq}-c_d)M_g=0.
\]

This is instantaneous phase equilibrium, not a prediction of bubble growth
time. In steady mode the transport equations retain their fluxes and omit
storage. They do not reset each cell to a prescribed fresh-oil charge.

## Pressure and volume

`p_cav` is the oil vapor pressure in this model, not the air saturation
pressure. Oil vapor is in equilibrium and its mass is neglected. Gas partial
pressure is \(p_g=p-p_v(T)\). The free gas occupies
the remaining gap volume:

\[
M_g=\rho_g(p_g,T)h(1-\theta_l),\qquad
p\ge p_v,\quad 0\le\theta_l\le1.
\]

For a noncondensable gas mixed with a small amount of oil vapor this uses
the partial-pressure approximation. At high vapor concentrations a
multicomponent fugacity model would be needed. The oil vapor pressure is
\(p_v(T)=p_{\rm cav}\exp[E_v(1/T_* - 1/T)]\), with
`oil_vapor_temperature_coeff = E_v` and the property reference temperature
\(T_*\). Zero slope deliberately holds the vapor pressure constant. Fit these
inputs to the actual oil; no FW68D vapor-pressure curve is supplied.

The nonlinear residual is

\[
a+b-\sqrt{a^2+b^2+2g}=0,\quad
a={\rho_g\over\rho_{g,*}},\quad b=1-\theta_l,\quad
g={M_g\over\rho_{g,*}h}.
\]

For nonnegative free mass this enforces \(ab=g\) and the pressure/liquid-volume
bounds. With zero gas it reduces to the JFO vapor-cavity complementarity
condition. With gas it enforces the equation of state and gap volume within
the same solve as conservation. No subsequent gas-volume cap repairs the
answer by deleting or hiding mass.

`alpha_gas` reports the void occupied by the gas-containing bubble mixture.
Oil vapor and noncondensable gas coexist in that bubble; they are not two
spatially disjoint volumes. In the pure vapor limit the liquid deficit is
`1 - film_content`, while `alpha_gas` is zero.

## Fluxes and viscosity

The model assumes a common depth-averaged phase velocity. The mobility uses
\(\rho_m=(M_o+M_d+M_g)/h\) and the volume-weighted viscosity
\(\mu_m=\theta_l\mu_l+(1-\theta_l)\mu_g\):

\[
F_{p,f}=-\left({\rho_m h^3\over12\mu_m}\right)_f
{p_R-p_L\over\Delta s}\,\ell_f,
\qquad
F_{i,f}=Y_{i,\rm up}F_{p,f}+{U\over2}h_f{M_{i,\rm up}\over h_{\rm up}}\ell_f.
\]

The mobility is harmonically interpolated. Pressure and Couette contributions
have their own upstream donors; they can oppose one another. The known face
gap is the mean of neighboring gaps. Upwinding the geometry itself introduced
avoidable first-order bias in the load and JFO free boundary. Axial faces
have no Couette term. Dissolved diffusion is balanced by the opposite oil
diffusion flux, preserving total mixture mass. With no gas and negligible
vapor viscosity, the pressure mobility reduces to \(\rho_lh^3/(12\mu_l)\),
the classical mass-conserving Reynolds limit.

This homogeneous constitutive approximation does not resolve individual
bubbles, bubble slip, surface tension, nucleation sites, collapse shocks or
erosion. Thin-film validity and experimental cavity topology must be assessed
separately. The volume-weighted viscosity is an explicit modeling choice,
not a universal experimentally validated oil/refrigerant mixing law.

## Boundary conditions

An axial `DIRICHLET` condition fixes absolute reservoir pressure and allows
both inflow and outflow. The pressure distance is half a cell. Incoming fluid
uses that end's dissolved fraction and free mass fraction. Outgoing fluid
uses the adjacent interior composition and temperature, including the
properties used in the boundary mobility. Otherwise changing the temperature
of a downstream reservoir changes the calculated discharge resistance even
when no fluid returns from that reservoir. `NEUMANN` represents a sealed end
in the new model. No internal inlet is necessary.

For a pressure patch, the net mixture source is obtained from its otherwise
unconstrained storage and face-flux balance. A positive source brings in
reservoir composition; a negative source removes local composition. The
species equations include that source. `inlet_species = index c q` assigns
the dissolved fraction in incoming liquid and free mass fraction in the
incoming mixture to an already declared, zero-based groove/hole index.
Unspecified `c = -1` inherits the initial dissolved charge; `q = 0` supplies
liquid without pre-existing free gas.

This pressure-patch approximation does not resolve the flow through a feed
pipe or the three-dimensional entrance jet. It represents a reservoir whose
specified pressure is attained over the selected film cells.

## Numerical acceptance

PETSc solves the four or five unknowns together with a sparse colored finite-difference
Jacobian. The matrix pattern and coloring persist between solves. The
unscaled mass balances are also evaluated using the actual discrete fluxes.
The normalized residual must be below the acceptance tolerance; a solver
termination reason alone is insufficient. State bounds and the free-gas
thermodynamic domain are checked before publishing a new state.

The Fischer–Burmeister residual is evaluated in rationalized form when its
sum is positive, avoiding subtraction of nearly equal numbers. Physical
inventories are not clipped afterward. Numerical zero tolerances permit free
mass down to \(-10^{-9}M_{g,*}\) and liquid fraction up to \(1+10^{-8}\);
these are acceptance tolerances, not allowed physical negative phases. A
resolved R410A free phase must also pass the stable EOS-domain check.

The previous physical inventories remain fixed throughout a transient solve.
Calling the nonlinear solver twice without accepting a new time step must
return the same state, not advance phase exchange twice. Steady solves may
use conservative pseudo-time continuation to obtain an initial guess; only
the final equations with zero storage establish steady convergence.

## Evidence obtained so far

The unit checks below pass on the current source. Final workstation results,
including before/after figures and the failed experimental comparison, are
in [the campaign report](../cases/90_diagnostics/coupled_validation/SUMMARY.md).
Intermediate local implementations remain archived under their own hashes.

* Closed-cell oil and total-species conservation; pressure rises with gas
  release while gas and liquid occupy exactly the prescribed gap.
* Gas compression at zero absorption satisfies Boyle's law.
* An independent quadratic equilibrium pressure solution agrees with the
  coupled constant-volume Henry/ideal-gas flash.
* Axial full-film pressure matches the exact linear profile, and transported
  concentration matches an independently evaluated implicit upwind solution
  in both flow directions.
* Circular-hole and groove inflow/outflow cases conserve species. Drainage
  does not reset the local composition to the nominal supply composition.
* A rotating pure-vapor case actually cavitates and satisfies the JFO
  pressure/film-content complementarity condition.
* R410A density matches independent CoolProp calls, and the real EOS is used
  within closed-cell release tests at 20, 40, 60 and 80 °C.
* The controlled 48 × 8 oil–air workstation comparison at 1 ms gives loads
  of 598.6536, 598.6482 and 598.6455 N at 10, 5 and 2.5 µs. The segregated
  comparator gives 1258.1, 3259.9 and 1401.5 N and violates gap volume. All material constants
  in this comparison are synthetic. See the preserved inputs and manifest in
  `cases/90_diagnostics/coupled_validation/final_evidence`.

These checks establish numerical behavior and limiting-case consistency.
They are not an experimental validation of FW68D–air or FW68D–R410A.

## Literature comparisons under evaluation

The equilibrium air model provides a constitutive/steady-bearing comparison
(Li et al., *[J. Tribol.](https://doi.org/10.1115/1.4006702)*, 2012).
The submerged bearing provides cavity pressure and rupture/reformation
measurements with two flooded ends and no internal supply groove
(Etsion and Ludwig, *[J. Lubr. Technol.](https://doi.org/10.1115/1.3253174)*, 1982).
Transient photographs distinguish release, compression and delayed absorption
(Reinke et al., *[Fluids](https://doi.org/10.3390/fluids7050159)*, 2022).
Dynamic vapor-cavity results supply a separate JFO transient benchmark
(Brewe, *[J. Tribol.](https://doi.org/10.1115/1.3261288)*, 1986).

The available recent R410A thin-film study uses FVS50, not FW68D, and its
constitutive coefficients are not all supplied in the article
(Li et al., *[Int. J. Refrig.](https://doi.org/10.1016/j.ijrefrig.2026.106916)*, 2026).
It is therefore not a source of an identified FW68D calibration. Quantitative
reproduction must state which inputs were measured, inferred or unavailable.
