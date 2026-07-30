# Journal-bearing analytical benchmarks

Verification of the isothermal, rigid, prescribed-eccentricity Reynolds solve
against the two exact limits that bracket any finite bearing.

```bash
python run_benchmarks.py ../../build-windows-mingw/pancake.exe
```

`analytic.py` holds the closed forms; `run_benchmarks.py` generates configs, runs
`pancake.exe`, and reports the error normalised by the pressure amplitude.

## Results (2026-07-31, 720 x 40/60 mesh, eps = 0.6, 1000 rpm)

Case A -- infinitely-long (Sommerfeld) limit, `FULL_SOMMERFELD`, mid-plane
p(theta). Residual is the axial end effect and must vanish as L/D grows.

| L/D | max err | rms err |
|-----|---------|---------|
| 1   | 22.86%  | 14.94%  |
| 2   | 7.32%   | 4.94%   |
| 4   | 0.76%   | 0.51%   |
| 8   | 0.10%   | 0.04%   |

At L/D = 8 the peak pressure is 14.5944 bar against an exact 14.5923 bar (0.014%).

Case B -- short-bearing (Ocvirk) limit, `GUMBEL`, full p(theta, z). Residual is
the neglected circumferential pressure flow, which is O((L/D)^2).

| L/D    | max err | rms err | error ratio |
|--------|---------|---------|-------------|
| 0.5    | 37.12%  | 7.87%   | --          |
| 0.25   | 14.63%  | 2.95%   | 2.54        |
| 0.125  | 4.40%   | 0.86%   | 3.32        |
| 0.0625 | 1.17%   | 0.22%   | 3.77        |

The error ratio approaching 4 under successive halving of L/D confirms the
expected second-order convergence to the Ocvirk limit.

## Scope

These are verification against exact solutions of the same PDE the solver
discretizes. They bound discretization error in both operator terms; they say
nothing about whether the Reynolds equation models a given experiment, and they
do not exercise the cavitation or gas models (Case A has no cavitation, Case B
uses Gumbel clipping).
