# LGE thermal expansion, October 17, 2026

This campaign compares the original A-F reservoirs, loads and speeds with
homogeneous free expansion of the shaft and monolithic bearing/housing.
It uses the specified 304 analogue, alpha = 16e-6 K^-1, and dimensions at 20 C.
The exact grade and mechanical restraint are not inferred from this test.

## Inputs and implementation

The controlling prompt and full-precision geometry/experiment tables are copied
unchanged into `results/lge_thermal_expansion_20261017/inputs` in the parent
checkout. Experimental torque is exclusively `new_exp_torque_nm` from the
October 2 corrected comparison. The original measured TC1-3 remain independent
outputs. No TC2-supply, E_low or E_high perturbation is combined with expansion.

`thermal_expansion.py` computes each new geometry from the 20 C reference,
Ds0 = 0.043 m, Db0 = 0.043085 m, L0 = 0.04 m, Do0 = 0.085 m:

    Ds = Ds0 * (1 + alpha * (Ts - 20))
    Db = Db0 * (1 + alpha * (Tb - 20))
    c = (Db - Ds) / 2
    L = L0 * (1 + alpha * (Tb - 20))
    Do = Do0 * (1 + alpha * (Tb - 20))

Temperatures in these expressions are Celsius. The isothermal shaft supplies
Ts. Tb is the volume mean over the actual solid grid of each solve, using shell
weights `(r[k+1]^2-r[k]^2)/(Ro^2-R^2)`. All 16 radial fields are used. Each of the
16 native solid axial slabs repeats exactly into 12 of the 192 output rows.
The postprocessor verifies that equality before averaging. It preserves the
existing conduction inner radius R, rather than changing it to the physical
bore R+c. Equal shaft/bearing warming expands the original clearance by the
common scale; zero warming returns the reference dimensions.

`run_steady.py` accepts optional `geometry_by_case` entries containing
`journal_radius_m`, `radial_clearance_m`, `bearing_length_m`,
`bush_outer_radius_m` and `sensor_radius_m`. Unknown cases, incomplete entries,
nonfinite/nonpositive dimensions and sensors outside the physical wall fail
before native execution. Absence of this object preserves baseline rendering.

Geometry propagation was traced through these paths:

| Quantity | Implementation |
|---|---|
| Mesh, surface speed and torque | Native configuration R, c, L; existing mesh and postprocessing consume them. No native equations change. |
| Position seed | Valid epsilon times each case's c; native load balance re-solves epsilon and attitude. |
| Housing lateral cooling | `case_materials` passes actual Do into the existing housing correlation, including cache closure checks. |
| Annular ends and solid mesh | Native R, Ro and L drive the existing finite-volume conduction and end cooling. |
| Shaft disks | Actual R enters the disk correlation; two exposed disks give area ratio R/L. |
| Groove | Original 90-degree centre and angular width retained, with full axial extent. The reference arc scales with the bearing. |
| TC1-3 | Radius 0.026 times the bearing scale; z=L/2; angles 240/270/300 degrees. Log-radius interpolation matches the archived reporter. |
| Sommerfeld/friction | Reference coordinates remain available. Torque is compared directly; no radius-dependent friction conversion is needed. |

The shaft centre stays fixed with surface rotation. Bearing translation is a
steady position adjustment and introduces no translational surface velocity.

## Verification and execution

The task has its own branch and checkout under `tmp/lge-thermal-expansion`.
The parent's staged/unstaged October 12 sensitivity work is untouched. The
baseline archive SHA-256 is
`e414b50fcde4cce55a1426173cd28cbfeecd548b7fa56675a547a198f2b1c8f3`.
Native source differs only by a header comment. All six no-geometry rendered
inputs are byte-identical to the archived driver using the same arguments.
Final archived inputs match after output/checkpoint/property paths and the
initial position are accounted for. The archived 4000-iteration cap is retained.
The six supplied starting geometries and temperatures were independently
reconstructed from archived NPZ fields.

Local checks compile the native Windows solver and exercise solid energy,
native equilibrium and FW68D properties. Python checks cover zero rise, equal
warming, clearance consistency, cooling/sensor propagation, invalid inputs,
solid-volume weighting, axial mapping and the existing restart/acceptance tests.
The Windows runtime requires `C:/msys64/mingw64/bin` and Microsoft MPI on PATH.
The first test launch without those runtime paths failed; both the failed log
and the successful corrected launch are retained.

The remote root is `/home/user/pancake-runs/lge-thermal-expansion-20261017` on
afml-r6. Six independent four-rank workers occupy socket 0 cores 0-11 and socket
1 cores 0-11. Four of the 28 physical cores remain free. OMP, BLAS and MKL thread
counts are one. Source/config/executable hashes, MPI/Python/compiler versions,
rankfiles, commands, return codes, fields and failed attempts are retained.
Each complete six-case round was budgeted at 20-35 minutes; the planned
first/feedback/tight-repeat campaign at 1-3 hours and 25-70 core-hours.

Each worker first solves a fresh tighter nominal baseline. It then solves the
supplied geometry with the original tolerances, retaining its first-pass
prediction. Subsequent geometry updates use temperatures from the new fields.
Feedback stops at |delta c| <= 1e-9 m and all of |delta Ds|, |delta Db|,
|delta L|, |delta Do| <= 1e-8 m. At most five geometry updates are allowed.
A fresh tighter solve of the settled geometry checks numerical resolution;
if its temperatures require another permitted update, that update uses the
tighter criteria. Capped feedback is reported as unresolved.

Original acceptance is force-vector error <0.002; FBNS, temperature, layer and
global energy criteria 1e-6; relative mass closure 1e-8; positive gap and valid
temperature/property bounds; exterior h mismatch <1e-3 and surface change
<0.01 K. Tight repeats request force 2e-5, FBNS/thermal 1e-8, h mismatch 1e-5
and surface change 1e-4 K. These are numerical choices, not experimental
uncertainties. Geometry changes always start fresh, without old checkpoints.
Only a verified same-geometry cooling iteration may reuse its own state.

The comparison reports signed torque error `100*(M_CFD/M_exp-1)` and sensor
error `T_CFD-T_exp`. Improvement is judged separately for every observable.
Any expansion response comparable to stopping/repeat differences is labeled
unresolved. No coefficient, bath temperature, clearance or torque correction is
fitted. The model excludes nonuniform distortion, mounting restraint, elastic
load deformation, contact and different stainless grades.

The original setup's geometry descriptions remain reference metadata. The optional per-case geometry and rendered native configuration identify each actual solve. The final analysis explicitly reconstructs physical bore/wall thickness, scaled groove arc, disk area, sensor coordinates and shaft surface speed for every final geometry. A capped 24 by 16 expanded-B native Windows smoke test also verified field output and moving sensor sampling; it is retained as a nonproduction check.

## Completed comparison

All six cases closed their force, thermal, mass, energy and exterior-cooling
checks and their geometry feedback. B and E required one geometry update;
A, C, D and F required none. The first pass is retained separately from the
updated and tight predictions. The final column below isolates geometry using
the fresh tight nominal/expanded pair; comparing archived baseline directly
with final expansion also includes the small tightening correction.

| Case | Archived torque (N m) | First expansion (N m) | Final expansion (N m) | Archived error (%) | Final error (%) | Matched expansion change (%) |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.268886284 | 0.269050846 | 0.269116780 | -27.94502 | -27.88325 | +0.061278 |
| B | 0.308960678 | 0.309796906 | 0.309884739 | -33.11301 | -32.91296 | +0.271758 |
| C | 0.197387852 | 0.197475736 | 0.197520637 | -45.00972 | -44.97272 | +0.047152 |
| D | 0.148962350 | 0.149076586 | 0.149074470 | -8.01280 | -7.94357 | +0.076815 |
| E | 0.289615693 | 0.290397360 | 0.290433106 | +37.98010 | +38.36954 | +0.270900 |
| F | 0.143498190 | 0.143619788 | 0.143619569 | +31.50014 | +31.61137 | +0.084735 |

Expansion removes respectively 0.1580%, 0.5491%, 0.0576% and 0.8818% of the
original torque residual in A-D. It increases the residual in E and F by
0.9843% and 0.3537%. All 18 absolute TC errors worsen. Matched sensor warming
ranges from 0.000937 to 0.019726 K, far smaller than the existing disagreement.
The free-expansion sensitivity therefore does not explain the torque or
temperature discrepancies.

All torque and TC response signs are resolved against the normal/tight
response difference. The largest torque response difference is 5.59% of the
small expansion response (case C); the largest sensor response difference is
1.48e-5 K. These are stopping-sensitivity results on the existing mesh, not
statistical uncertainties or a new mesh-independence study.

The maximum final relative residuals across A-F are: force 1.175e-5,
mass 2.141e-11, global energy 2.910e-10, temperature 4.840e-9 and
layer/solid profile 6.723e-9. Final FBNS is at most 1.135e-16. Exterior h
mismatch is at most 5.263e-6 and surface-temperature change at most
2.655e-5 K. All satisfy the tighter requested gates. Actual solved dimensions,
proposed next dimensions, Ts/Tb, all sensor errors, minimum gap, attitude,
eccentricity, maximum pressure, cavity fraction and every geometry stage are
in the full machine-readable and Markdown results.

The successful production attempt ran from 05:16:14 to 07:25:02 UTC on
October 17, 2026 (2 h 8 min 48 s). It contains 20 accepted stages and
87 successful native MPI invocations. Including the six preserved startup
failures, native wall time times rank count totals 31.291 core-hours; this
is allocated MPI time, not measured CPU utilization. The Linux executable
is byte-identical to the archived baseline executable. Windows passed eight
focused CTests; all 59 Linux CTests passed after supplying the omitted test
fixture, with failed and successful logs both retained.

The local evidence root is
`C:/Dev/pancake/results/lge_thermal_expansion_20261017`.
`analysis/RESULTS.md`, `analysis/results.json` and `analysis/comparison.png`
contain the comparison. `remote_evidence` preserves inputs, raw VTK/NPZ fields,
checkpoints, source, executables, commands, rankfiles, histories and failed
attempts. `retrieval_audit.json` records the archive SHA-256 and the checks
of all accepted configuration/field/checkpoint files after extraction.
The archive contains the collector's status at packaging time; the final
retrieval manifest records completion. Original September 28 raw baseline
fields remain in their original directory with separate audited hashes.
`analysis/report_render_audit.json` records the display-only clarification
that nominal baseline rows do not apply geometry feedback; the original
archived report and all numerical results remain preserved.
