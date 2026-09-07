"""Render traceable thermal, JFO and gas-cavity comparisons from saved evidence."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'physics_comparison'))
from plot_style import apply_style, blues, save_figure


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def column(data, name):
    return np.array([float(row[name]) for row in data])


def records(folder):
    data = json.loads((folder / 'manifest.json').read_text())
    failed = [r['name'] for r in data['runs'] if r['return_code']]
    if failed:
        print('Excluded failed cases from successful-result curves:', ', '.join(failed))
    return {r['name']: r for r in data['runs'] if r['return_code'] == 0}


def thermal_figures(folder, output, refinement=None):
    runs = records(folder)
    if refinement:
        runs.update(records(refinement))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for j, species in enumerate(['air', 'r410a']):
        for mode, color, style in [('isothermal', '0.45', '--'),
                                   ('energy_equation', '#0072B2', '-')]:
            name = f'{species}_{mode}'
            if name not in runs:
                raise ValueError(f'Missing successful comparison: {name}')
            data = rows(folder / f'{name}_history.csv')
            t = column(data, 'time_s') * 1000
            keep = t > 0
            label = 'Isothermal' if mode == 'isothermal' else 'Coupled heating'
            axes[0, j].plot(t[keep], column(data, 'load_N')[keep], style,
                            color=color, label=label)
            axes[1, j].plot(t[keep], column(data, 'temperature_max_K')[keep] - 273.15,
                            style, color=color)
        axes[0, j].text(.03, .97, '(a) Oil–air' if j == 0 else '(b) Oil–R410A',
                        transform=axes[0, j].transAxes, va='top')
        axes[0, j].legend(loc='upper center',bbox_to_anchor=(.5,1.2),ncol=2,fontsize=12)
        axes[0, j].set_ylabel(r'$W$ (N)')
        axes[1, j].set_ylabel(r'$T_{\max}$ (°C)')
        for ax in axes[:, j]:
            ax.set_xlabel(r'$t$ (ms)')
    fig.tight_layout()
    save_figure(fig, str(output / 'thermal_effect'))
    plt.close(fig)

    names = sorted([name for name in runs if name.startswith('air_thermal_grid')]
                   + ['air_energy_equation'], key=lambda name: runs[name]['n_theta'])
    time_names = ['air_energy_equation', 'air_thermal_dt5e-05', 'air_thermal_dt2.5e-05']
    if all(n in runs for n in names + time_names):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        for ax, selected, key, factor, label in [
            (axes[0], names, 'n_theta', 1, r'$N_\theta$'),
            (axes[1], time_names, 'dt_s', 1e6, r'$\Delta t$ ($\mu$s)')]:
            ax.plot([runs[n][key] * factor for n in selected],
                    [runs[n]['final']['load_N'] for n in selected], '-o')
            ax.set_xlabel(label)
            ax.set_ylabel(r'$W(t=20\,\mathrm{ms})$ (N)')
            ax.ticklabel_format(axis='y',style='plain',useOffset=False)
        fig.tight_layout()
        save_figure(fig, str(output / 'thermal_resolution'))
        plt.close(fig)


def benchmark_figures(folder, output, jfo_folder=None):
    runs = records(folder)
    jfo_folder = jfo_folder or folder
    jfo_runs = records(jfo_folder)
    exact = rows(jfo_folder / 'jfo_exact.csv')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(column(exact, 'theta_deg'), column(exact, 'pressure_Pa') / 1e5,
                 'k-', label='Independent JFO')
    grids = sorted(int(name.split('_')[1]) for name in jfo_runs if name.startswith('jfo_'))
    for n, color in zip(grids, blues(len(grids), dark_first=False)):
        name = f'jfo_{n}'
        if name not in jfo_runs:
            raise ValueError(f'Missing successful JFO case: {name}')
        data = rows(jfo_folder / f'{name}_profile.csv')
        angle = column(data, 'theta_deg')
        angle = np.where(angle < 99, angle + 360, angle)
        order = np.argsort(angle)
        axes[0].plot(angle[order], column(data, 'pressure_Pa')[order] / 1e5,
                     '--', color=color, label=rf'$N_\theta={n}$')
    axes[0].set_xlabel(r'$\theta$ (deg)')
    axes[0].set_ylabel(r'$p$ (bar, absolute)')
    axes[0].legend(fontsize=12)
    axes[1].loglog(grids, [jfo_runs[f'jfo_{n}']['pressure_relative_l2'] for n in grids], '-o')
    axes[1].set_xlabel(r'$N_\theta$')
    axes[1].set_ylabel(r'Relative pressure $L_2$ error')
    axes[1].set_xticks(grids,labels=[str(n) for n in grids])
    axes[1].xaxis.set_minor_formatter(plt.NullFormatter())
    fig.tight_layout()
    save_figure(fig, str(output / 'jfo_verification'))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.axhspan(225, 245, color='0.88', label='Photograph: 235° ± 10° extent bound')
    for name, label, color, marker in [
        ('etsion_120_equilibrium', 'Equilibrium, 120 × 32', '#56B4E9', 's'),
        ('etsion_240_equilibrium', 'Equilibrium, 240 × 64', '#0072B2', 'o'),
        ('etsion_120_retained', 'No absorption, 120 × 32', '#D55E00', '^')]:
        if name not in runs:
            continue
        angles = runs[name]['cavity_angles']
        thresholds = sorted(float(key) for key in angles)
        values = [next(value[2] for key, value in angles.items() if float(key) == t)
                  for t in thresholds]
        ax.semilogx(thresholds, values, marker=marker, color=color, label=label)
    ax.set_xlabel(r'Liquid-deficit threshold, $1-\theta_l$')
    ax.set_ylabel('Circumferential cavity extent (deg)')
    ax.set_ylim(0, 375)
    ax.legend(fontsize=12, loc='lower left')
    fig.tight_layout()
    save_figure(fig, str(output / 'etsion_gas_comparison'))
    plt.close(fig)


def material_figure(output):
    data = rows(ROOT / 'cases/ref_data/fvc68b_r410a_1999_viscosity.csv')
    fig, ax = plt.subplots(figsize=(7, 5.3))
    for t, color, marker in zip([23, 45, 65], blues(3), ['s', 'o', '^']):
        selected = [row for row in data if float(row['T_C']) == t]
        selected.sort(key=lambda row: float(row['c_mass_fraction']))
        ax.semilogy(column(selected, 'c_mass_fraction'), column(selected, 'mu_Pa_s') * 1000,
                    marker=marker, color=color, label=rf'$T={t}\,^\circ\mathrm{{C}}$')
    ax.set_xlabel(r'Dissolved R410A mass fraction, $c_d$')
    ax.set_ylabel(r'$\mu_l$ (mPa s)')
    ax.text(.97, .97, 'FVC68B PVE oil\nDOE report, Table 35',
            transform=ax.transAxes, ha='right', va='top')
    ax.legend(loc='lower left')
    fig.tight_layout()
    save_figure(fig, str(output / 'fvc68b_viscosity_reference'))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--thermal', type=Path)
    parser.add_argument('--benchmarks', type=Path)
    parser.add_argument('--refinement', type=Path)
    parser.add_argument('--jfo', type=Path)
    parser.add_argument('--output', type=Path, default=HERE / 'figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    material_figure(args.output)
    if args.thermal:
        thermal_figures(args.thermal, args.output, args.refinement)
    if args.benchmarks:
        benchmark_figures(args.benchmarks, args.output, args.jfo)


if __name__ == '__main__':
    main()
