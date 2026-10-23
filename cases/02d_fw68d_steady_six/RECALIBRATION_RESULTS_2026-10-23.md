# A–F rerun with recalibrated loads (2026-10-23)

Inputs: [RECALIBRATION_2026-10.md](RECALIBRATION_2026-10.md). Generated data, run
records and scripts (not in git): `results/fw68d_recalibrated_20261023/` in the
local checkout, also archived in the research vault
(`attachments/LGE-2026-10-23-*`).

## 1. Smooth surfaces, full driver closure

`run_steady.py --experiment experiment_recalibrated_20261023.json`, 720×192,
16 film layers, load, thermal and energy closure as in the Sep 28 campaign.
Solver built from `main` `acc20fd` source (executable SHA-256 `4acdca3f…e25`),
Linux, GCC 13, OpenMPI, PETSc 3.19.6, 2 ranks. `fbns_tol` tightened from
1e-6 to 1e-8 (`steady_setup_fbns1e-8.json` in the results folder) for D–A after D's first attempt
failed the 1e-8 liquid-balance acceptance at −1.46e-8; F used 1e-6. All six
cases reached `load_and_thermal_closed`, force error 0.02–0.10 %.

| Case | rpm | $W$ old → new (N) | $h_{min}$ old → new (µm) | CFD torque old → new (N·m) | Change |
|---|---:|---|---|---|---:|
| A | 200 | 5276 → 3992 | 5.91 → 7.40 | 0.2689 → 0.2418 | −10.1 % |
| B | 1500 | 3238 → 2460 | 12.33 → 14.89 | 0.3090 → 0.2865 | −7.3 % |
| C | 200 | 1973 → 1510 | 12.95 → 15.52 | 0.1974 → 0.1834 | −7.1 % |
| D | 200 | 383 → 314 | 31.81 → 33.68 | 0.1490 → 0.1470 | −1.3 % |
| E | 1500 | 461 → 374 | 35.71 → 37.01 | 0.2896 → 0.2857 | −1.4 % |
| F | 200 | 191 → 170 | 37.02 → 37.59 | 0.1435 → 0.1431 | −0.3 % |

Torque error $100(M_{CFD}/M_{exp}-1)$, recalibrated load:

| Case | vs Oct 21 up | vs Oct 21 down | vs Sep 11 (rpm-specific) |
|---|---:|---:|---:|
| A | −58.2 % | −65.0 % | −35.2 % |
| B | −54.4 % | −57.1 % | −38.0 % |
| C | −67.5 % | −72.9 % | −48.9 % |
| D | −59.9 % | −69.3 % | −9.2 % |
| E | −24.2 % | −31.3 % | +36.1 % |
| F | −54.5 % | −66.4 % | +31.1 % |

The load recalibration lowers CFD torque by 7–10 % at high load and by about
1 % at high $S$. It widens every gap slightly; it does not change the picture.

## 2. Mixed lubrication on A and B (`codex/lge-af-recalibrated-contact`)

Greenwood–Tripp contact plus Patir–Cheng flow factors (commit `75299d8`),
starting from each closed smooth configuration with its exterior cooling held
fixed. $\sigma$ split equally between shaft and bush; $\sigma/\beta = 0.01$ and
$\eta\beta\sigma$ chosen to give the measured $K$; $E^* = 1.099\times10^{11}$ Pa,
$\mu_b = 0.07$. Tilt 0.416 mrad in the load plane.

| Case | Variant | $h_{min}$ (µm) | Torque (N·m) | Contact torque (N·m) |
|---|---|---:|---:|---:|
| A | smooth aligned (control) | 7.40 | 0.24180 | — |
| A | smooth, tilt | 1.92 | 0.24296 | — |
| A | mixed, tilt, $\sigma$ 0.4 µm, $K$ 0.001 | 1.92 | 0.24357 | 0 |
| A | mixed, tilt, $\sigma$ 0.5 µm, $K$ 0.002 | 1.93 | 0.24390 | 0 |
| B | smooth aligned (control) | 14.89 | 0.28653 | — |
| B | smooth, tilt | 9.51 | 0.28849 | — |
| B | mixed, tilt, nominal | 9.51 | 0.28859 | 0 |
| B | mixed, tilt, upper | 9.51 | 0.28865 | 0 |

The control reproduces the closed smooth torque to all printed digits on the
second binary. At the recalibrated loads the tilted minimum gap is 3.8–4.8
roughness heights in A, at or above the $H<4$ cut-off of the $F_{5/2}$ fit, so
the asperity load is zero; flow factors add 0.3 % (A) and 0.04 % (B). Contact
does not contribute to the A–F torque gap at these inputs.

Note: `mixed_lubrication.hpp` attributes the $F_{5/2}$ fit
$4.4086\times10^{-5}(4-H)^{6.804}$ to Arcoumanis et al. (1997). A secondary
check (Jedynak, Tribol. Int. 2018, Table B.1) attributes it to Hu et al.
(1994); Arcoumanis used a different polynomial. The fit is 11 % low at $H=0$.
