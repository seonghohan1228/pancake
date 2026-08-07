"""Report the resolved pressure jump immediately outside the Ferron groove."""

import json
import math
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
SOURCE_CASE = CASE_DIR.parent / "01_ferron_wallbc_fields"
sys.path.insert(0, str(SOURCE_CASE))

import solver_io as S


def periodic_distance_deg(a, b):
    """Smallest absolute angular distance in degrees."""
    return abs((a - b + 180.0) % 360.0 - 180.0)


def main():
    equilibria = json.loads((CASE_DIR / "grid_equilibrium.json").read_text())
    for tag, result in equilibria.items():
        run_dir = CASE_DIR / result["pinned"]["dir"]
        angles, pressure = S.midplane_profile(run_dir, "p")
        inlet = [i for i, angle in enumerate(angles)
                 if periodic_distance_deg(angle, 90.0) <= 9.0]
        west, east = inlet[0], inlet[-1]
        d_theta = 360.0 / len(angles)
        q_supply = result["q_supply"]
        print(
            f"{tag} n_theta={len(angles)} dtheta={d_theta:.6f}deg "
            f"q={q_supply * 6.0e4:.6f}L/min q*dtheta={q_supply*d_theta:.9e}\n"
            f"  west: groove {angles[west]:.6f}deg {pressure[west]:.3f}Pa; "
            f"film {angles[west-1]:.6f}deg {pressure[west-1]:.3f}Pa; "
            f"jump={pressure[west]-pressure[west-1]:.3f}Pa\n"
            f"  east: groove {angles[east]:.6f}deg {pressure[east]:.3f}Pa; "
            f"film {angles[(east+1)%len(angles)]:.6f}deg "
            f"{pressure[(east+1)%len(angles)]:.3f}Pa; "
            f"jump={pressure[east]-pressure[(east+1)%len(angles)]:.3f}Pa"
        )


if __name__ == "__main__":
    main()
