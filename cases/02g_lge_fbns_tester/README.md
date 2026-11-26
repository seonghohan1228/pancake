# 02g LGE tester, FBNS gas-effect pairs (500 RPM, c = 0.274)

Campaign `lge-fbns-tester-20261117`. Gives the change in film torque and h_min caused by
gaseous cavitation at the X4 (768.1 N) and X5 (2610.3 N) tester conditions. The two runs
of a pair differ only in the gas model. Absolute torque and temperature stay with the
accepted thermal COUPLED campaign (`lge-thermal-v3-c274-20261114`).

| Run | Load (N) | Gas model |
|---|---|---|
| P4v / P5v | 768.1 / 2610.3 | SINGLE_PHASE, SCALAR_PCAV, p_cav = 500 Pa (VC only) |
| P4g / P5g | 768.1 / 2610.3 | TWO_PHASE, EQUILIBRIUM, LOCAL_PSAT (uniform bubble pressure at feed c and T) |

Model, as on the optical rig: FBNS, isothermal at TC0 = 297.539 K, no cross-film or
housing conduction, constant liquid viscosity and density, groove and both ends at
16.58 bar, 720 x 192.

Driver `run_fbns.py`: the solver has no native load balance for FBNS with dissolved gas,
so the driver runs a Newton loop (FD Jacobian, Broyden updates, step limits, backtracking)
over STATIC solves. Converged when |F|/W is within 0.1 % of 1 and the force direction within
0.1 deg of the load direction.
- v run: eps fixed at the accepted COUPLED gas-off h_min, unknowns (mu, psi). The mu found
  is the pair's viscosity.
- g run: mu of its v run, unknowns (eps, psi).

Driver-level stand-ins (solver unchanged):
- `dissolved_gas_species = PROPANE`: the solver accepts R410A only under COUPLED. The
  species sets only the free-gas density. Each run records `species_and_clamp_inert`
  (no free gas and dissolved fraction 0.274 everywhere), i.e. whether it matters.
- `solution_table_clamp = true`: the equilibrium flash queries the FW68D chart outside its
  temperature support (it aborts otherwise) even though the cavity sits at the bubble
  pressure. Inert under the same check.

Usage: `python3 run_fbns.py --case P4v --work DIR --solver EXE --table fw68d.csv --ranks N
[--rankfile FILE]`. Run P4v/P5v before P4g/P5g (the g run reads its pair's accepted.json).
Outputs per case: `accepted.json` or `failure.json`, `fields.npz`, every solve under `runs/`.
