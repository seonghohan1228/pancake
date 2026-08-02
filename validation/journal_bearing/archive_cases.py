"""Archive the case data behind each validation verdict.

The results JSONs hold scalar metrics and mid-plane profiles, which is what the
comparison plots need, but two things were only ever in the scratch tree:

  * the exact `config.txt` each case ran with.  Regenerating it from the harness
    is not the same guarantee -- the harness changes, the archived input does not.
  * the 2-D (theta, z) fields, without which no contour or colormap of pressure,
    temperature or film content can be made without re-running the solver.

`_work/` is gitignored and runs to over a gigabyte of VTK, so it is not the
answer.  This writes a compact, committable archive instead: per case, the
verbatim config plus the four fields as compressed float32, about 300 kB each
against 2.6 MB of .vts.

Usage:
    $PY archive_cases.py                 # archive every case found under _work/
    $PY archive_cases.py <workdir> ...   # archive only these trees
"""
import json, os, re, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_raimondi_boyd import read_field

HERE = os.path.dirname(os.path.abspath(__file__))
ARCHIVE = os.path.join(HERE, "case_archive")
FIELDS = ("p", "film_content", "h", "T")


def grid_from_config(cfg_text):
    """n_theta, n_z as the solver actually ran them."""
    values = {}
    for line in cfg_text.splitlines():
        if "=" in line and not line.strip().startswith("#"):
            key, _, value = line.partition("=")
            values[key.strip()] = value.strip()
    return int(values["n_theta_global"]), int(values["n_z_global"])


# Two tiers, because committing every 2-D field is 50 MB of which 20 MB is
# root-finding iterations that back no verdict, and git keeps it forever.
#
#   configs  -- EVERY case, verbatim, ~2 kB each. This is what makes a case
#               exactly reproducible, and it is cheap enough to keep all of.
#   fields   -- only cases a figure or a verdict actually draws on.
#
# Regenerate the full 2-D archive from a scratch tree at any time with
# `archive_cases.py --all <workdir>`; it is derived data, not a source of truth.
FIELD_TIER = (
    "thd_adiabatic", "thd_metalbush", "isothermal",   # Ferron headline verdicts
    "x_eps",                                          # cross-solver comparison
    "d_g", "d_b",                                     # grid / beta diagnosis
    "cvg_",                                           # free-boundary convergence
    "ss_LD",                                          # analytic Swift-Stieber
    "rb_LD1_flooded_eps0.6_mu0.01",                   # a passing R&B case
    "rb_LD1_ungrooved_eps0.6_mu0.01",                 # the ungrooved collapse
)
FIELD_TIER = FIELD_TIER + ("brito_N",)   # Brito operating grid

# Secant and fixed-point iterates are tagged `_e<k>` / `_g<k>` with a SINGLE digit
# (both loops cap at 6), so the pattern is anchored and single-digit on purpose:
# a substring test wrongly ate `d_g1440` and `d_g2880`, which are grid labels, and
# a multi-digit pattern would eat them again.
#
# Within an iterate family the LAST index is the converged answer, not a
# throwaway -- Brito's groove fixed point reports `..._g2` as its result. So keep
# the highest index per family and drop the rest, rather than assuming a `_final`
# tag that only the Ferron harness writes.
ITERATE_RE = re.compile(r"^(?P<base>.*_[eg])(?P<idx>\d)$")


def superseded_iterates(names):
    """Names that a later iterate of the same family replaces."""
    best = {}
    for name in names:
        m = ITERATE_RE.match(name)
        if m:
            base, idx = m.group("base"), int(m.group("idx"))
            if idx > best.get(base, (-1, None))[0]:
                best[base] = (idx, name)
    keep = {entry[1] for entry in best.values()}
    return {n for n in names if ITERATE_RE.match(n) and n not in keep}


def wants_fields(name, dropped):
    if name in dropped:
        return False
    return any(tag in name for tag in FIELD_TIER)


def find_cases(roots):
    for root in roots:
        for dirpath, _dirnames, filenames in os.walk(root):
            if "config.txt" not in filenames:
                continue
            vts = os.path.join(dirpath, "out", "flat", "processor0",
                               "solution_0_0.vts")
            if os.path.exists(vts):
                yield dirpath, vts


def archive(roots, verbose=True, keep_all=False):
    os.makedirs(ARCHIVE, exist_ok=True)
    index, total = [], 0
    cases = sorted(find_cases(roots))
    names = [os.path.relpath(d, HERE).replace(os.sep, "__")
             .replace("_work__", "").replace("..__", "") for d, _ in cases]
    dropped = superseded_iterates(names)
    for (case_dir, vts), name in zip(cases, names):
        cfg_text = open(os.path.join(case_dir, "config.txt")).read()
        try:
            nth, nz = grid_from_config(cfg_text)
        except (KeyError, ValueError):
            if verbose:
                print(f"  skip {name}: unreadable grid")
            continue

        arrays = {}
        if keep_all or wants_fields(name, dropped):
            for field in FIELDS:
                try:
                    arrays[field] = np.asarray(read_field(vts, nth, nz, field),
                                               dtype=np.float32)
                except Exception:
                    pass      # not every case writes every field

        out = os.path.join(ARCHIVE, name + ".npz")
        np.savez_compressed(out, config=np.array(cfg_text), n_theta=nth,
                            n_z=nz, **arrays)
        size = os.path.getsize(out)
        total += size
        index.append(dict(case=name, file=os.path.basename(out),
                          n_theta=nth, n_z=nz,
                          fields=sorted(arrays), bytes=size))
        if verbose:
            tier = "config+fields" if arrays else "config only"
            print(f"  {name:<44} {nth}x{nz} {size/1e3:>8.0f} kB  {tier}", flush=True)

    json.dump(index, open(os.path.join(ARCHIVE, "index.json"), "w"), indent=1)
    if verbose:
        print(f"\n{len(index)} cases, {total/1e6:.1f} MB total -> {ARCHIVE}")
        print("Load one with:")
        print("  d = np.load('case_archive/<case>.npz')")
        print("  p = d['p']            # (n_z, n_theta) absolute pressure [Pa]")
        print("  print(str(d['config']))")
    return index


if __name__ == "__main__":
    args = sys.argv[1:]
    keep_all = "--all" in args
    roots = [a for a in args if not a.startswith("--")] or [os.path.join(HERE, "_work")]
    archive(roots, keep_all=keep_all)
