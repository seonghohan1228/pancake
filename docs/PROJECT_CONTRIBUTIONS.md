# Project Contribution Summary

This summary lists changes that materially advance the solver, model coverage,
or maintainability of `pancake`. It intentionally excludes one-off error fixes,
temporary build repairs, and the work of creating or redesigning the Windows GUI
itself.

## Numerical Models

- Added a mass-conserving Elrod-Adams/JFO cavitation path using `theta` as the
  universal film-content variable, while keeping the older Gumbel pressure-clamp
  path selectable.
- Added axial boundary controls that distinguish pressure constraints from
  film-content constraints, so pressure and JFO liquid content are not
  accidentally conflated at z boundaries.
- Added structured inlet support for supply pressure regions, including groove
  and circular inlet definitions that can be represented in `config.txt`.
- Added moving-bearing dynamics with external load, support stiffness/damping,
  initial displacement/velocity, squeeze-film terms, and film-thickness safety
  controls.
- Added a film-averaged energy-equation model with temperature, wall heat
  transfer, thermal advection, and viscous heat-generation post-processing.
- Added a fluid-property layer for pure oil and dissolved/free-gas mixtures,
  with configurable density, viscosity, solubility, and gas release/resorption
  behavior.

## Post-Processing and Output

- Added Cartesian pressure, viscous, and total fluid-force integration using a
  consistent bearing-side force convention.
- Added velocity reconstruction for output and downstream force, torque, and
  heat-generation calculations.
- Standardized VTK output names toward compact OpenFOAM-style fields such as
  `p`, `T`, `U`, `Fp`, `Fv`, `F`, `film_content`, and `h`.
- Added output-selection controls so expensive or noisy visualization fields can
  be enabled or disabled without changing the simulation state.
- Added per-step logging with a configurable simplified/detailed verbosity mode
  so run progress and diagnostics are visible independently of VTK write cadence.

## Testing and Validation

- Expanded regression coverage for Elrod/JFO boundary behavior, inlet pressure
  handling, force integration, fluid-property models, output-name compatibility,
  and thermal heat-generation behavior near cavitation fronts.
- Added documentation explaining the current governing equations, numerical
  sequence, cavitation active set, energy equation, moving-bearing model, and
  dissolved/free-gas property coupling.

## Current Model Boundaries

- The gaseous-cavitation implementation is currently a segregated property and
  local release/resorption model. It is not yet a full dissolved/free-gas
  advection-diffusion transport solve.
- The energy model is film-averaged and uses wall heat-transfer closures. It is
  not yet a fully coupled ETHD deformation/thermal-viscosity iteration.
- The nominal Courant number shown in the live summary is based on Couette
  advection only; pressure-driven velocities require a solved pressure field and
  remain a post-solve diagnostic rather than an input-summary value.
