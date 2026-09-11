"""Plot provisional steady liquid-temperature rises with common geometry axes."""
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from plot_style import apply_style,save_figure,blues,smart_legend

ROOT=Path(__file__).resolve().parent
apply_style()


def mosaic(rows,mode):
    """Use one scale within each cooling scenario; label temperature rise."""
    selected=[r for r in rows if r["mode"]==mode]
    lookup={(r["load_target"],r["rpm"]):r for r in selected}
    vmax=max(r["rise_max_K"] for r in selected)
    step=2 if vmax<15 else 10
    vmax=math.ceil(vmax/step)*step
    fig,axes=plt.subplots(3,3,figsize=(13,8.7),sharex=True,sharey=True)
    fig.subplots_adjust(left=.09,right=.86,bottom=.10,top=.98,wspace=.22,hspace=.10)
    for i,w in enumerate((1000,3000,5000)):
        for j,n in enumerate((500,1000,2000)):
            ax=axes[i,j]
            row=lookup[(w,n)]
            fields=np.load(ROOT/row["dir"]/"fields.npz")
            theta=(fields["theta_deg"]+180)%360-180
            order=np.argsort(theta)
            rise=fields["T"][:,order]-row["thermal_inputs"]["inlet_temperature_K"]
            plot=ax.contourf(theta[order],fields["z_over_L"],rise,
                             levels=np.linspace(0,vmax,25),cmap="inferno",norm=Normalize(0,vmax))
            ax.contour(theta[order],fields["z_over_L"],fields["film_content"][:,order],
                       levels=[.999],colors="#56B4E9",linewidths=.7)
            ax.axvspan(85,95,color=".65",alpha=.65)
            label=rf"$W={w/1000:.0f}$ kN"+"\n"+rf"$N={n}$ rpm"
            ax.text(.02,.96,label,transform=ax.transAxes,ha="left",va="top",fontsize=12,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
            ax.set_xlim(-180,180); ax.set_ylim(-.5,.5)
            ax.set_xticks([-180,-90,0,90,180]); ax.set_yticks([-.5,0,.5])
            if i==2: ax.set_xlabel(r"$\theta$ (deg)")
            if j==0: ax.set_ylabel(r"$z/L$")
    cb=fig.colorbar(plot,cax=fig.add_axes([.89,.14,.025,.78]),label=r"$T-T_{\mathrm{in}}$ (K)")
    cb.minorticks_off(); cb.set_ticks(np.linspace(0,vmax,6))
    save_figure(fig,str(ROOT/f"temperature_{mode}_contours"))
    plt.close(fig)


def contrast(rows):
    """Show peak temperature rise versus speed for each assumed cooling condition."""
    fig,ax=plt.subplots(figsize=(7,4.6),layout="constrained")
    colors=blues(3)
    for mode,marker,ls in (("cooled","o","-"),("adiabatic","s","--")):
        for w,color in zip((1000,3000,5000),colors):
            subset=sorted((r for r in rows if r["mode"]==mode and r["load_target"]==w),key=lambda r:r["rpm"])
            ax.plot([r["rpm"] for r in subset],[r["rise_max_K"] for r in subset],
                    marker=marker,linestyle=ls,color=color,label=rf"$W={w/1000:.0f}$ kN, {mode}")
    ax.set_xlabel(r"$N$ (rpm)"); ax.set_ylabel(r"$T_{\max}-T_{\mathrm{in}}$ (K)")
    smart_legend(ax,fontsize=11)
    save_figure(fig,str(ROOT/"temperature_cooling_comparison"))
    plt.close(fig)


def cooling_mosaic(rows):
    """Compare cooling scenarios at 3 kN with a labelled scale for each row."""
    fig,axes=plt.subplots(2,3,figsize=(13,6.4),sharex=True,sharey=True)
    fig.subplots_adjust(left=.09,right=.86,bottom=.13,top=.98,wspace=.22,hspace=.14)
    for i,mode in enumerate(("cooled","adiabatic")):
        selected=[r for r in rows if r["mode"]==mode and r["load_target"]==3000]
        lookup={r["rpm"]:r for r in selected}
        step=2 if mode=="cooled" else 10
        vmax=math.ceil(max(r["rise_max_K"] for r in selected)/step)*step
        for j,n in enumerate((500,1000,2000)):
            ax=axes[i,j]; row=lookup[n]
            fields=np.load(ROOT/row["dir"]/"fields.npz")
            theta=(fields["theta_deg"]+180)%360-180; order=np.argsort(theta)
            rise=fields["T"][:,order]-row["thermal_inputs"]["inlet_temperature_K"]
            plot=ax.contourf(theta[order],fields["z_over_L"],rise,
                             levels=np.linspace(0,vmax,25),cmap="inferno",norm=Normalize(0,vmax))
            ax.contour(theta[order],fields["z_over_L"],fields["film_content"][:,order],
                       levels=[.999],colors="#56B4E9",linewidths=.7)
            ax.axvspan(85,95,color=".65",alpha=.65)
            ax.text(.02,.96,mode+"\n"+rf"$N={n}$ rpm",transform=ax.transAxes,
                    ha="left",va="top",fontsize=12,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
            ax.set_xlim(-180,180); ax.set_ylim(-.5,.5)
            ax.set_xticks([-180,-90,0,90,180]); ax.set_yticks([-.5,0,.5])
            if i==1: ax.set_xlabel(r"$\theta$ (deg)")
            if j==0: ax.set_ylabel(r"$z/L$")
        position=axes[i,0].get_position()
        cb=fig.colorbar(plot,cax=fig.add_axes([.89,position.y0,.025,position.height]),
                        label=r"$T-T_{\mathrm{in}}$ (K)")
        cb.minorticks_off(); cb.set_ticks(np.linspace(0,vmax,5))
    save_figure(fig,str(ROOT/"temperature_comparison_contours"))
    plt.close(fig)


def main():
    rows=json.loads((ROOT/"temperature_summary.json").read_text())
    if len(rows)!=18:
        raise RuntimeError("Expected 18 completed thermal screening points")
    for mode in ("cooled","adiabatic"):
        mosaic(rows,mode)
    contrast(rows)
    cooling_mosaic(rows)


if __name__=="__main__":
    main()
