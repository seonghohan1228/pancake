"""Plot archived controlled comparison data, using the lab's scientific style."""
import csv
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'physics_comparison'))
from plot_style import apply_style, blues, clean_log_ticks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence',type=Path,default=HERE/'evidence')
    parser.add_argument('--output',type=Path,default=HERE/'figures')
    args=parser.parse_args()
    apply_style()
    figures=args.output
    figures.mkdir(parents=True,exist_ok=True)
    colors=blues(3)
    manifest=json.loads((args.evidence/'manifest.json').read_text())
    data={r['name']:r for r in manifest['runs']}
    fig,axes=plt.subplots(1,2,figsize=(13,5.3))
    for ax,model in zip(axes,['fbns','coupled']):
        for color,dt in zip(colors,[2.5,5,10]):
            rows=list(csv.DictReader((args.evidence/f'{model}_{dt:g}us_history.csv').open()))
            t=np.array([float(r['time_s']) for r in rows])*1e3
            load=np.array([float(r['load_N']) for r in rows])
            keep=t>=.01-1e-10
            ax.plot(t[keep],load[keep]/(1000 if model=='fbns' else 1),color=color,label=rf'$\Delta t={dt:g}\,\mathrm{{\mu s}}$')
        ax.set_xlabel(r'$t$ (ms)')
        ax.set_ylabel(r'$W$ (kN)' if model=='fbns' else r'$W$ (N)')
        ax.set_xlim(.01,1)
        ax.legend(frameon=False,fontsize=14)
    axes[0].text(.02,.98,'(a) Segregated',transform=axes[0].transAxes,va='top',fontsize=17)
    axes[1].text(.02,.98,'(b) Coupled',transform=axes[1].transAxes,va='top',fontsize=17)
    axes[1].margins(y=.12)
    fig.tight_layout()
    fig.savefig(figures/'air_load_timestep.png',dpi=300)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7,5.3))
    for model,color,marker in [('fbns','#969696','s'),('coupled',colors[0],'o')]:
        dt=np.array([2.5,5,10])
        err=[max(data[f'{model}_{v:g}us']['max_abs_gas_volume_excess'],1e-16) for v in dt]
        ax.plot(dt,err,marker=marker,color=color,label='Segregated' if model=='fbns' else 'Coupled')
    ax.set_xlabel(r'$\Delta t$ ($\mu$s)')
    ax.set_ylabel('Maximum gap-volume defect')
    ax.set_yscale('log')
    clean_log_ticks(ax)
    ax.legend(frameon=False,fontsize=15)
    fig.tight_layout()
    fig.savefig(figures/'air_volume_closure.png',dpi=300)
    plt.close(fig)


if __name__=='__main__': main()
