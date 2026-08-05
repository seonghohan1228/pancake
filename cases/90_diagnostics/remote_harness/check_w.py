"""Verify that Fp / xB really do vary frame-to-frame in a finished run.

analyse_remote.py reported W_std = 0 over [0.03, 0.06] while the eccentricity
still moved, which is inconsistent; this dumps the raw per-frame values (and
cross-checks the first against the last MPI piece) to find out which is wrong.
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

NAME_RE = 'Name="{}"[^>]*>'


def uni(path, name, ncomp=3):
    b = Path(path).read_bytes()
    m = re.search(NAME_RE.format(name).encode(), b)
    if not m:
        return None
    return [float(v) for v in b[m.end():m.end() + 400].split(b"<")[0].split()[:ncomp]]


run = Path(sys.argv[1])
pvd = run / "flat" / "results.pvd"
root = ET.parse(pvd).getroot()
items = sorted((float(d.get("timestep")), pvd.parent / d.get("file"))
               for d in root.iter("DataSet"))
print(f"{len(items)} frames")
for t, f in items[::10]:
    r = ET.parse(f).getroot()
    pieces = [f.parent / p.get("Source") for p in r.iter("Piece")]
    fp0 = uni(pieces[0], "Fp")
    fpN = uni(pieces[-1], "Fp")
    xb = uni(pieces[0], "xB")
    print(f"t={t:.4f}  Fp[piece0]={fp0}  Fp[piece{len(pieces)-1}]={fpN}  xB={xb}")
