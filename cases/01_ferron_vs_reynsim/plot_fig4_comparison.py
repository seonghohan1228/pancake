"""Fig. 4 (2000 rpm, 4 kN) three-way comparison, per user instruction:
Ferron (1983) experiment / Last year (ReynSim case_00 centerline, condition
inferred - see SUMMARY section 8) / This year (JFO): FBNS+THD J-adiabatic.
Case_00 datum caveat: its cavitated plateau is exactly 0.0; plotted as gauge
as-is. If that vintage was absolute like the KSME export, the curve shifts
down by 0.101 MPa."""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

from plot_style import apply_style, save_figure, MARKERS, COLORS
from plot_vs_reynsim import solver_midplane, read_csv

CASE_DIR = Path(__file__).parent
THD = CASE_DIR / ".." / "01_ferron_thd"
REF = CASE_DIR / ".." / "ref_data"


def read_case00():
    th, p = [], []
    with open(CASE_DIR / "reynsim_data" / "case_00_pressure_centerline.csv") as f:
        rd = csv.reader(f)
        next(rd)
        for row in rd:
            p.append(float(row[0]))
            th.append(float(row[1]))
    import math
    # ReynSim frame = ours + 90 deg (proven in SUMMARY 2.2, same convention
    # across their scripts); plateau is 0.0 -> plotted as gauge as-is.
    x = [(math.degrees(t) - 90.0) % 360.0 for t in th]
    o = sorted(range(len(x)), key=lambda i: x[i])
    return [x[i] for i in o], [p[i] / 1e6 for i in o]


apply_style()
eq = json.loads((THD / "equilibrium.json").read_text())
run_dir = THD / eq["fig4_Jadiabatic"]["flow"]["dir"]

fig, ax = plt.subplots(figsize=(6.6, 4.6))

xs, ys = solver_midplane(run_dir)
ax.plot(xs, ys, "-", lw=1.6, color=COLORS[2], label="This year (JFO)")

xr, pr = read_case00()
ax.plot(xr, pr, "-", lw=1.8, color=COLORS[1], label="Last year")

xe, ye = read_csv(REF / "ferron1983_fig4_2000rpm_4kN_pressure.csv", 1.0e-6)
ax.plot(xe, ye, MARKERS[1], mfc="none", mec="k", ls="none",
        label="Ferron (1983)")

ax.set_xlabel(r"$\theta - \theta_{\mathrm{groove}}$ (deg)")
ax.set_ylabel(r"$p - p_{\mathrm{amb}}$ (MPa)")
ax.set_xlim(0, 360)
ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
ax.legend(loc="upper left")
save_figure(fig, str(CASE_DIR / "fig4_three_way"))
print("wrote fig4_three_way.png")
