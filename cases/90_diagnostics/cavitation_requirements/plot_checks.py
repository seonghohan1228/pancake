"""Reproduce the kinetic and pure-R410A property diagnostics; no FW68D fit."""
import csv
import json
import sys
from pathlib import Path

import CoolProp
from CoolProp.CoolProp import PropsSI, PhaseSI
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'physics_comparison'))
from plot_style import apply_style, save_figure, COLORS, blues


def main():
    evidence = HERE / 'evidence'
    figures = HERE / 'figures'
    evidence.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    apply_style()
    data = np.genfromtxt(evidence / 'kinetics.csv', delimiter=',', names=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), layout='constrained')
    invariants = {}
    for rate, color, label in [(20, COLORS[0], r'$k_a=k_r$'),
                               (0, COLORS[1], r'$k_a=0$')]:
        d = data[data['absorption_rate'] == rate]
        for ax, y in zip(axes, [d['free_mass'] * 1e6, d['alpha_gas'] * 100]):
            ax.plot(d['time'] * 1e3, y, color=color, label=label)
            ax.set_xlabel(r'$t$ (ms)')
        invariants[str(rate)] = {
            'max_relative_oil_drift': float(np.max(np.abs(d['oil_mass'] / d['oil_mass'][0] - 1))),
            'max_relative_gas_drift': float(np.max(np.abs(d['total_gas_mass'] / d['total_gas_mass'][0] - 1))),
            'final_free_mass_kg_m2': float(d['free_mass'][-1]),
            'final_alpha': float(d['alpha_gas'][-1])}
    axes[0].set_ylabel(r'$M_g$ (mg/m$^2$)')
    axes[1].set_ylabel(r'$\alpha_g$ (%)')
    axes[0].legend(loc='upper left', fontsize=14)
    save_figure(fig, str(figures / 'absorption_vs_compression'))
    plt.close(fig)

    gas_constant = PropsSI('GAS_CONSTANT', 'R410A') / PropsSI('MOLAR_MASS', 'R410A')
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), layout='constrained')
    shades = blues(4)
    for idx, temp_c in enumerate([20, 40, 60, 80]):
        temp = temp_c + 273.15
        curve = []
        for pressure in np.linspace(1e5, 2e6, 191):
            phase = PhaseSI('T', temp, 'P', pressure, 'R410A')
            rho = PropsSI('D', 'T', temp, 'P', pressure, 'R410A')
            gas = phase in ('gas', 'supercritical_gas')
            ideal = pressure / (gas_constant * temp)
            propane = pressure / (188.55 * temp)
            rows.append([temp_c, pressure, phase, rho, ideal, propane])
            curve.append([pressure / 1e5, rho if gas else np.nan,
                          100 * (ideal / rho - 1) if gas else np.nan])
        curve = np.array(curve)
        for ax, y in zip(axes, [curve[:, 1], curve[:, 2]]):
            ax.plot(curve[:, 0], y, label=f'{temp_c} °C', color=shades[idx])
            ax.set_xlabel(r'$p$ (bar, absolute)')
    axes[0].set_ylabel(r'$\rho_{\mathrm{R410A}}$ (kg/m$^3$)')
    axes[1].set_ylabel('Ideal-gas density error (%)')
    axes[0].legend(loc='upper left', fontsize=13)
    save_figure(fig, str(figures / 'r410a_eos'))
    plt.close(fig)
    with (evidence / 'r410a_eos.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['temperature_C', 'pressure_Pa', 'phase', 'density_kg_m3',
                         'ideal_R410A_kg_m3', 'wrong_propane_kg_m3'])
        writer.writerows(rows)
    metadata = {'CoolProp_version': CoolProp.__version__,
                'CoolProp_git_revision': CoolProp.__gitrevision__,
                'EOS': 'HEOS::R410A (fixed nominal blend, pseudo-pure)',
                'EOS_reference': 'https://doi.org/10.1023/A:1025048800563',
                'oil': 'none; these are pure-refrigerant properties, not FW68D/R410A VLE',
                'kinetics': 'prescribed-pressure parcel; rates illustrative, not calibrated',
                'gas_constant_J_kg_K': gas_constant,
                'kinetic_invariants': invariants}
    (evidence / 'checks.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    # Replot the archived study comparison, explicitly without a new solver run.
    with (HERE.parent.parent / '05_etsion_extent' / 'comparison.csv').open(encoding='utf-8') as stream:
        archive = list(csv.DictReader(line for line in stream if not line.startswith('#')))
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), layout='constrained')
    for ax, rpm in zip(axes, [1840, 3000]):
        selected = [r for r in archive if r['model'] == 'fbns' and int(r['rpm']) == rpm]
        x = [float(r['ps_kPa']) for r in selected]
        # Unwrap downstream angles across 360 degrees for a continuous curve.
        unwrap = lambda v: v-360 if v > 180 else v
        measured = [unwrap(float(r['th_end_exp'])) for r in selected]
        simulated = [unwrap(float(r['th_end_sim'])) for r in selected]
        ax.errorbar(x, measured, yerr=5, fmt='o', ls='none', label='Measured', color=COLORS[0])
        ax.plot(x, simulated, 's-', label='Archived FBNS', color=COLORS[1])
        ax.set_xlabel(r'$p_s$ (kPa, gauge)')
        ax.set_ylabel(r'$\vartheta_{\mathrm{end}}$ (°)')
        ax.text(0.95, 0.94, f'{rpm} rpm', transform=ax.transAxes, ha='right', va='top', fontsize=14)
    axes[0].legend(loc='lower left', fontsize=13)
    save_figure(fig, str(figures / 'etsion_reformation'))
    plt.close(fig)
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
