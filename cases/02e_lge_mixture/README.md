# LGE undersaturated FW68D–R410A comparison, 18 October 2026

This campaign retains all 17 selected M01–M17 experimental windows and the
existing 103 pure-FW68D experiment points and six CFD markers. The controlling
input is `LGE-2026-10-18-mixture-pancake-run-prompt.txt`. The campaign uses nominal
geometry and preserves existing results.

The user confirmed that the absolute-pressure sensor measures the chamber,
that both bearing ends and the groove share this pressure, and that exposed
housing and shaft surfaces contact liquid at approximately TC0. The user
describes the fluid as undersaturated, with measured refrigerant content of
approximately 30%. Conditioning at 15 bar and 20 °C does not establish the in-run
composition or identify the concentration basis.

## Physical scope and assumptions

| Input or model | Selected value | Evidence and limitation |
|---|---|---|
| Load and speed | Each saved mean force and nominal RPM | All 17 selected IDs retained; no torque-based selection |
| Rotation | Positive native theta follows rotation, as in the accepted pure-oil frame | All selected measured/requested motor speeds have the same negative sign as the baseline; laboratory CW/CCW mapping is not independently established |
| Reservoir pressure | Each measured mean in bar absolute × 100000 Pa/bar | User confirms common chamber/end/groove connection |
| Composition | Fixed chart-30% viscosity scenario | Percent basis remains unverified; not assigned mass fraction 0.3 |
| Viscosity | Log nu linear in Celsius between decoded nodes | Local solved film temperature; no extrapolation or clamping |
| Density | Constant 950 kg/m³ | Conditional plotting assumption, not a measurement |
| Species/release | None tracked; release omitted | SINGLE_PHASE effective-fluid surrogate, not native coupled refrigerant validation |
| Cavity floor | 0 Pa absolute | Explicit oil-vapor-limit idealization with refrigerant release suppressed; not measured mixture bubble pressure or refrigerant partial pressure |
| Fluid heat capacity | 2010 J/(kg K) | Existing BVC32-PVE analogue, not measured mixture cp |
| Fluid conductivity | 0.13 W/(m K) | Existing BVC32-PVE analogue, not measured mixture k |
| Exterior buoyancy beta | 0.7284/950 = 0.0007667368421 K⁻¹ | Declared pure-oil density-slope analogue; independent Boussinesq approximation under constant bulk density |
| Steel conductivity | 15 W/(m K) | Existing 304 stainless analogue |
| Feed and liquid bath | Each measured TC0 | Proxy boundary assumption; film, shaft and housing temperatures remain solved |
| Geometry | R = 21.5 mm, c = 42.5 µm, L = 40 mm, housing OD = 85 mm | Nominal; no thermal expansion |

The native R410A EOS does not supply mixture solubility, transport properties,
or dissolution energy. The chart contains isobars, but its original pressure
and concentration bases are explicitly unspecified. Moreover, the current
COUPLED path rejects the baseline finite-volume cross-film thermal treatment.
It is therefore not enabled. The surrogate cannot predict degassing, changing
composition, or a physically established mixture cavity threshold.

The ten decoded viscosity nodes cover −26.0429574434 to 84.0172809983 °C.
Checks give 8.2435157220 mm²/s at 30 °C and 6.6992248806 mm²/s at 40 °C.
The viscosity-only CSV uses `c_mass_fraction=0` as the solver's no-tracked-species
coordinate and blank pressure; this bookkeeping does not describe the physical
mixture as refrigerant-free and supplies no solubility law. Full legacy solution
tables retain their prior validation. All out-of-domain temperatures fail.

Cooling reuses the horizontal-cylinder, annular vertical-plate and rotating-disk
correlations from the [accepted campaign](../02d_fw68d_steady_six/README.md).
All use the same surrogate liquid viscosity, density, cp and k at their film
temperature. The cylinder and annular faces use the stated beta. Constant bulk
density affects Reynolds mass flux and energy capacity; the configured negative
density slope supplies only the separate exterior buoyancy approximation.

## Source and numerical plan

The accepted pure-oil source is the verified uncommitted September 28 snapshot
based on Git commit `d197ad594e651db440c4a246ef272efdcf1b3018`, not that commit
alone. Its source archive SHA-256 is
`e414b50fcde4cce55a1426173cd28cbfeecd548b7fa56675a547a198f2b1c8f3`.
The baseline executable hash is
`62ddd6a0906943e2a7942330742757042626036602e8ab1b2073e0e563969551`.
The task started from local main `3f75685b19d7a451021cd6749f51a2e8b4b30e2c` in
an isolated worktree. All six original driver renders match the archived
configurations byte-for-byte. The comparison preserves historical CFD markers;
they are not recalculated under mixture assumptions.

Use 720 × 192 film cells, 16 cross-film layers and 16 × 16 radial/axial housing cells.
Keep the top full-length groove, its 13.2983125248° width, upward applied load
and TC1/TC2/TC3 probes at r = 26 mm, z = 20 mm, theta = 240/270/300°.
Solve eccentricity and attitude numerically; the shaft center is fixed.

Pilots are M02 (lowest saved S), M11 (median S), and M16 (highest S).
After all pass, run the other fourteen automatically. Then run M02/M16 with
tighter tolerances and M02 with 32 film layers at those tighter tolerances.
The 20-equilibrium plan is one 17-case campaign plus three focused numerical
checks. No physical parameter is fitted to torque or temperature.

Accept only clean native convergence and independent force residual below
0.002, FBNS pressure/complementarity residuals below 1e-10 and temperature/profile
residuals below 1e-6. Require absolute relative mass imbalance at most 1e-8,
global energy residual at most 1e-6, cooling-coefficient mismatch below 1e-3
and wall-temperature change below 0.01 K. Gap, viscosity and density must be
finite and positive, with admissible film content and temperatures inside the
strict table domain. Caps, restored best iterates and nonzero
return codes are rejected. Tight checks use force tolerance 2e-5, FBNS tolerance
1e-12, thermal/profile and global energy tolerances 1e-8, cooling mismatch
1e-5 and temperature change 1e-4 K. The mass-balance limit remains 1e-8.
A tight check does not prove grid independence; the 32-layer comparison tests
only the film-normal direction.

The first M02 pilot failed the unchanged mass-balance gate. Its signed relative
mass residual was -4.553754639e-8, exceeding the 1e-8 magnitude limit despite
pressure-equation residual G = 6.023062e-9. The three initial pilot attempts
remain preserved under `/home/user/pancake-runs/lge-mixture-20261018`.
The retry uses the isolated
`/home/user/pancake-runs/lge-mixture-20261018-pressure-refined` directory and
the normal/tight FBNS tolerances of 1e-10/1e-12 above. This changes numerical
pressure accuracy only. No physical input changes, no acceptance limit is
relaxed, and no failed trial is counted as a converged result.

The initial M11/M16 solves passed native and cooling checks but exposed a
worker lookup bug: final verification used the original position seed instead
of the last cooling call's recorded seed. Both saved runs pass an independent
recheck with the recorded seed. The worker and reporter now use that seed and
require the same run directory; a regression test covers accepted, missing and
wrong-run closure. The brief pressure-refined attempt was interrupted and
preserved before completing its pilots. The corrected snapshot runs under
`/home/user/pancake-runs/lge-mixture-20261018-corrected`.

Each full remote case uses four MPI ranks. Six disjoint allocations use 24 of
afml-r6's 28 physical cores, leaving four free. Pilot timing supersedes the
initial 2–6 h wall-time allowance. Historical eight-case sensitivity cost was
9.877 allocated core-hours; scaling to 20 equilibria gives 24.69 core-hours before
allowance for changed viscosity, pressure, tight convergence and the 32-layer
case. Reserve 15–25 GB for retained raw files and retrieval. These are planning
estimates, not measured runtimes. All failed attempts remain preserved.

## Reproduction and outputs

All 17 normal cases passed independent raw-field and provenance verification.
Predicted bearing torque is 46.68-86.31% below the saved corrected experimental
torque. The separate motor diagnostic is 6.01-37.41% low. TC1-3 errors span
-0.13355 to +0.70534 K, with mean absolute error 0.35923 K over all 51 probes.
The largest normal force error is 0.09902%, relative mass imbalance magnitude
3.478e-10 and relative energy residual 1.780e-10.

Tighter M02/M16 checks change torque by +0.033291%/+0.009700%, respectively.
Their largest absolute probe-temperature change is 0.000124 K. Doubling M02's
film layers at tight tolerances changes torque by +3.0671e-6% and the three
probes by -1.1628e-5, -7.8264e-6 and +6.9052e-6 K. These observed differences
are much smaller than the experiment/model discrepancies; they are not a
discretization uncertainty bound or a full mesh-independence claim. All raw
evidence and both prior attempts have been retrieved and verified. Shared-driver
metrics retain its
baseline admission labels; case settings and the campaign worker and reporter
enforce the tighter criteria stated above.

```text
python cases/02e_lge_mixture/prepare_inputs.py --inputs EVIDENCE --out PREFLIGHT
python cases/02e_lge_mixture/run_mixture.py --prepare PREFLIGHT --inputs PREPARED
python cases/02e_lge_mixture/run_mixture.py --inputs PREPARED --check-inputs AUDIT.json
python cases/02e_lge_mixture/launch_campaign.py --root REMOTE_ROOT
python cases/02e_lge_mixture/summarize_mixture.py --root REMOTE_ROOT --execution-manifest REMOTE_ROOT/execution.json --plot-style-dir REMOTE_ROOT/style --out NEW_ANALYSIS
```

Full runs occur only after a local build/check, an isolated remote source/input
transfer, a remote rebuild and tests. Per-solve records retain source, input,
driver, worker, executable, configuration, checkpoint and rankfile hashes,
commands, runtime versions, MPI ranks, return codes and convergence histories.

The reporter verifies raw output, load/thermal/cooling closure and provenance.
It writes one `case_results.csv`, a methods/results note, a verification JSON
and the Stribeck PNG using the user's plot-data style. Failed or missing cases
remain explicit rows. Experimental S is retained exactly, based on measured
TC2; CFD uses the same x coordinate and $f=M/(WR)$. Bearing torque error is
$100(M_{\rm CFD}-M_{\rm residual})/M_{\rm residual}$. The separate motor diagnostic
adds the RPM-specific parasitic reference to predicted bearing torque and
compares the sum with measured directional motor torque, using measured motor
torque as its percentage denominator. TC errors are CFD minus experiment
in kelvin at the three housing probes; TC0 is a boundary proxy, not an
independent validation point.

Native torque integrates Couette and pressure shear over the cylindrical
journal bearing surface. It excludes exposed shaft-disk drag, seals, motor
losses and separately resolved groove sidewalls. The transferred RPM-specific
parasitic correction is unverified for pressurized mixture. Temporal min/max
bars remain distinct from unknown composition/calibration uncertainty and
the numerical checks. Agreement alone does not validate the mixture model.

## Verified completion (2026-10-19)

The 20 accepted jobs completed in 95 min 25.7 s. Their 80 native calls returned
zero; all 59 remote native tests passed. The archive also preserves 12 earlier
calls, including three intentionally interrupted return-code-1 calls. Total
allocated native-solver time across these 92 calls was 20.466822
core-hours. All corrected final cases pass independent raw-field, conservation,
pressure/complementarity, thermal/profile and cooling checks.

The 2,702,683,177-byte evidence archive has SHA-256
`c6752b66cb5a1dc6a418326514781798411a0e90a3bf9c1ec753239cf60678f2`.
Local retrieval additionally verified 4,463 evidence-file hashes. Full raw
files are in `results/lge_mixture_20261018/remote_evidence` in the parent checkout;
the four final analysis products are in its sibling `analysis` directory.
`retrieval_audit.json` records the external completed archive manifest; the
status snapshot inside the archive was necessarily written before archiving
finished. Research attachments and the existing discussion/project log were
updated without replacing earlier evidence.

The executed source bundle SHA-256 is
`f70ea67f76f7e31013b6867521ce4d4124612d03dcd77eea35d072847c0de30e`;
the remote binary SHA-256 is
`fd841b121a6e7db184801d4380f6aa2126a91a4ff8ebe71403ff2e4ca695ce92`.
The final reporter overlay adds signed/absolute dimensional torque errors; its
hash is `3860c5a04f1dc2e918eab62bf7bc359940260500ee0981817d9ce5f29721488e`.
The original reporter and overlay manifests are preserved. Solver, driver,
worker and physical inputs remained frozen during the campaign.
