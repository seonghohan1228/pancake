"""Shared MPI launch policy for the validation harnesses.

The harnesses ran the solver on a single rank. pancake decomposes theta across
ranks, so the benchmark meshes (720-2880 circumferential cells) parallelise well
and the sweeps are the long pole in this work.

Rank count comes from `PANCAKE_VALIDATION_RANKS`, defaulting to 1 so existing
behaviour is unchanged unless asked for. Set it to the cores you want to spend:

    PANCAKE_VALIDATION_RANKS=4 python run_ferron.py <pancake.exe>

Rank-independence, measured 2026-08-02 (Ferron geometry, eps = 0.6, isothermal,
720x40):

    linear_rtol = 1e-8    W = 5711.2731 (1 rank) vs 5711.5528 (4)   4.9e-5
    linear_rtol = 1e-12   W = 5711.252515 both                      BIT-IDENTICAL

So the partitioning is exact and the residual difference at the default tolerance
is the linear solve, not the ghost exchange. Two things follow: 4 ranks may be
used for any of these sweeps, and load carries ~5e-5 of relative solver noise at
`linear_rtol = 1e-8`, which is far below the percent-level differences these
benchmarks resolve but is not zero.

Correctness note: rank count must not change the answer, and the ghost exchange
is covered by MPI tests (`test_elrod_a2_mpi6`, `test_schemes_mpi6`). Defect D5 is
the known exception -- `rho_cp` and `k` ghosts go stale between the gas solve and
the energy solve, which makes the ENERGY_EQUATION path rank-count dependent while
a gas model is active. The THD benchmarks here run with `steady_gas_model =
FROZEN` and no dissolved gas, so that path is not exercised, but any run that
turns a gas model on should be verified at one rank before its numbers are
quoted.
"""
import os
import shutil


def rank_count():
    """Ranks to launch, from PANCAKE_VALIDATION_RANKS (default 1)."""
    try:
        n = int(os.environ.get("PANCAKE_VALIDATION_RANKS", "1"))
    except ValueError:
        return 1
    return max(1, n)


def launcher():
    """`mpiexec` if it is on PATH, else None."""
    for name in ("mpiexec", "mpirun"):
        path = shutil.which(name)
        if path:
            return path
    return None


def solver_command(exe, config="config.txt"):
    """Argv to run one case, parallel when both requested and available."""
    n = rank_count()
    if n <= 1:
        return [exe, "--config", config]
    mpi = launcher()
    if mpi is None:
        raise RuntimeError(
            f"PANCAKE_VALIDATION_RANKS={n} but neither mpiexec nor mpirun is on "
            "PATH. Add the MPI runtime to PATH, or unset the variable to run "
            "serially. Falling back silently would make the sweep look fast "
            "while it was not, so this is an error.")
    return [mpi, "-n", str(n), exe, "--config", config]


def describe():
    n = rank_count()
    return f"{n} rank" + ("" if n == 1 else "s") + (" (serial)" if n == 1 else "")
