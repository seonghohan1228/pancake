"""Build inspectable load-speed maps and unwrapped steady-field contours."""
import csv
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np
from plot_style import apply_style, save_figure

ROOT = Path(__file__).resolve().parent
CFG = json.loads((ROOT/"settings.json").read_text())
apply_style()


def load_rows():
    """Load converged points and compute explicitly defined dimensionless groups."""
    rows = [json.loads(p.read_text()) for p in sorted(ROOT.glob("scan_n*_w*.json"))]
    for row in rows:
        speed = row["rpm"]/60
        pressure_scale = CFG["mu"]*2*math.pi*speed*(CFG["R"]/CFG["c"])**2
        row["S"] = CFG["mu"]*speed/(row["load_target"]/(2*CFG["R"]*CFG["L"]))*(CFG["R"]/CFG["c"])**2
        row["Pi_c"] = (CFG["p_end"]-row["p_cav"])/pressure_scale
        row["power_W"] = row["torque_Nm"]*2*math.pi*speed
        row["void_max_pct"] = 100*row["void_max"]
    return rows


def maps(rows):
    """Interpolate only within the sampled rectangular load-speed matrix."""
    rpms = sorted(set(r["rpm"] for r in rows))
    loads = sorted(set(r["load_target"] for r in rows))
    lookup = {(r["rpm"],r["load_target"]):r for r in rows}
    fig, axes = plt.subplots(1,3,figsize=(14,4.3),layout="constrained")
    specs = [("cav_area_pct",r"$A_{\mathrm{cav}}/A$ (%)","Blues"),
             ("void_max_pct",r"$(1-\vartheta)_{\max}$ (%)","Blues"),
             ("h_min_exact_um",r"$h_{\min}$ ($\mu$m)","viridis")]
    for ax,(key,label,cmap) in zip(axes,specs):
        values = np.array([[lookup[(n,w)][key] for n in rpms] for w in loads])
        plot = ax.contourf(rpms,np.array(loads)/1000,values,levels=12,cmap=cmap)
        lines = ax.contour(rpms,np.array(loads)/1000,values,levels=5,colors="black",linewidths=.5)
        ax.clabel(lines,fontsize=10,fmt="%.0f")
        xx,yy = np.meshgrid(rpms,np.array(loads)/1000)
        ax.scatter(xx,yy,s=13,c="black",zorder=3)
        ax.set_xlabel(r"$N$ (rpm)")
        ax.set_ylabel(r"$W$ (kN)")
        ax.set_xticks([500,1000,1500,2000])
        ax.tick_params(axis="x",labelsize=12)
        fig.colorbar(plot,ax=ax,label=label,shrink=.92)
    save_figure(fig,str(ROOT/"load_rpm_map"))
    plt.close(fig)


def field_mosaic(rows,field,basename):
    """Use identical color bounds and housing coordinates across all panels."""
    lookup = {(r["rpm"],r["load_target"]):r for r in rows}
    rpms,loads = [500,1000,2000],[1000,3000,5000]
    fig,axes = plt.subplots(3,3,figsize=(13,8.7),sharex=True,sharey=True)
    fig.subplots_adjust(left=.09,right=.86,bottom=.10,top=.98,wspace=.22,hspace=.10)
    vmax = max(r["p_max_bar"] for r in rows if r["rpm"] in rpms and r["load_target"] in loads)
    void_max = math.ceil(10*max(r["void_max"] for r in rows if r["rpm"] in rpms and r["load_target"] in loads))/10
    for i,w in enumerate(loads):
        for j,n in enumerate(rpms):
            ax = axes[i,j]
            row = lookup[(n,w)]
            a = np.load(ROOT/row["dir"]/"fields.npz")
            theta = (a["theta_deg"]+180)%360-180
            order = np.argsort(theta)
            theta = theta[order]
            fc = a["film_content"][:,order]
            value = np.maximum(0,1-fc) if field=="void" else a["p"][:,order]/1e5
            norm = None if field=="void" else PowerNorm(.5,vmin=0,vmax=vmax)
            levels = np.linspace(0,void_max,21) if field=="void" else np.linspace(0,math.sqrt(math.ceil(vmax/5)*5),25)**2
            plot = ax.contourf(theta,a["z_over_L"],value,levels=levels,
                              cmap="Blues" if field=="void" else "magma",norm=norm)
            ax.contour(theta,a["z_over_L"],fc,levels=[.999],colors="#D55E00",linewidths=.7)
            # Mark the full-length feed groove in housing coordinates.
            ax.axvspan(CFG["groove_angle"]-CFG["groove_width"]/2,
                       CFG["groove_angle"]+CFG["groove_width"]/2,color=".65",alpha=.65)
            label = rf"$W={w/1000:.0f}$ kN"+"\n"+rf"$N={n}$ rpm"
            ax.text(.02,.96,label,transform=ax.transAxes,ha="left",va="top",fontsize=12,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
            ax.set_xlim(-180,180)
            ax.set_ylim(-.5,.5)
            ax.set_xticks([-180,-90,0,90,180])
            ax.set_yticks([-.5,0,.5])
            if i==2: ax.set_xlabel(r"$\theta$ (deg)")
            if j==0: ax.set_ylabel(r"$z/L$")
    label = r"$1-\vartheta$" if field=="void" else r"$p_{\mathrm{abs}}$ (bar)"
    cax=fig.add_axes([.89,.14,.025,.78])
    cb=fig.colorbar(plot,cax=cax,label=label)
    cb.minorticks_off()
    cb.set_ticks(np.linspace(0,void_max,6) if field=="void" else np.linspace(0,math.ceil(vmax/5)*5,6))
    save_figure(fig,str(ROOT/basename))
    plt.close(fig)


def main():
    rows=load_rows()
    if len(rows)!=len(CFG["loads"])*len(CFG["rpms"]):
        raise RuntimeError(f"Matrix incomplete: {len(rows)} points")
    columns=["load_target","rpm","eps","psi","W","force_error_rel","p_max_bar",
             "cav_area_pct","mid_arc_deg","void_max","void_volume_mm3","h_min_exact_um",
             "S","Pi_c","power_W","R_p","dir"]
    with (ROOT/"sweep_summary.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,extrasaction="ignore")
        writer.writeheader(); writer.writerows(sorted(rows,key=lambda r:(r["load_target"],r["rpm"])))
    (ROOT/"sweep_summary.json").write_text(json.dumps(rows,indent=2))
    maps(rows)
    field_mosaic(rows,"void","cavity_contours")
    field_mosaic(rows,"pressure","pressure_contours")


if __name__=="__main__":
    main()
