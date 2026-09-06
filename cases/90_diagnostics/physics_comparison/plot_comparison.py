"""Render saved comparison CSVs with the lab's publication style (300 dpi PNG)."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from plot_style import apply_style, save_figure, COLORS

HERE = Path(__file__).resolve().parent
REVISIONS = ('previous', 'current')
LABELS = ('Previous', 'Repaired')


def read(path):
    return np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))


def draw(ax, x, y, revision, markers=False, lines=True):
    i = REVISIONS.index(revision)
    ax.plot(x, y, color=COLORS[i], ls=('-' if i == 0 else '--') if lines else 'none',
            label=LABELS[i], lw=2.0, marker=('s', 'o')[i] if markers else None,
            markevery=max(1, len(x) // 12), markersize=5, mec='black', mew=.5)
    ax.ticklabel_format(axis='y', style='plain', useOffset=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=HERE / 'evidence')
    parser.add_argument('--output', type=Path, default=HERE / 'figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.0), layout='constrained')
    for revision in REVISIONS:
        data = read(args.evidence / f'exchange_{revision}.csv')
        t = data['time'] * 1000
        draw(axes[0], t, 100 * (data['total_gas_area'] / data['total_gas_area'][0] - 1), revision)
        draw(axes[1], t, data['alpha_gas'] * 100, revision)
    axes[0].set_ylabel('Gas mass change (%)')
    axes[1].set_ylabel(r'$\alpha_g$ (%)')
    for ax in axes:
        ax.set_xlabel(r'$t$ (ms)')
        ax.set_xlim(0, 20)
    axes[0].legend(loc='lower left', handlelength=2)
    save_figure(fig, str(args.output / 'closed_exchange'))
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.3, 7.8), layout='constrained')
    quantities = [('pressure_max_Pa', 1e-5, r'$p_{\max}$ (bar)'),
                  ('load_N', 1, r'$W$ (N)'),
                  ('alpha_gas_max', 100, r'$\alpha_{g,\max}$ (%)'),
                  ('temperature_max_K', 1, r'$T_{\max}-T_0$ (K)')]
    for revision in REVISIONS:
        data = read(args.evidence / f'air_{revision}_history.csv')
        # Show only shared advanced states, 10 us <= t <= 1 ms.
        data = data[(data['comparison_time_s'] >= 1e-5 - 1e-12) &
                    (data['comparison_time_s'] <= .001 + 1e-12)]
        for ax, (key, factor, label) in zip(axes.flat, quantities):
            y = data[key] - (313.15 if key == 'temperature_max_K' else 0)
            draw(ax, data['comparison_time_s'] * 1000, y * factor, revision)
            ax.set_ylabel(label)
            ax.set_xlabel(r'$t$ (ms)')
            ax.set_xlim(0, 1)
    axes[0, 1].legend(loc='best', handlelength=2)
    save_figure(fig, str(args.output / 'air_transient'))
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.1), layout='constrained')
    for revision in REVISIONS:
        ferron = read(args.evidence / f'ferron_{revision}_profile.csv')
        angle = (ferron['angle_deg'] - 90) % 360
        order = np.argsort(angle)
        draw(axes[0], angle[order], (ferron['p'][order] - 101325) / 1e6, revision, True)
        draw(axes[1], angle[order], ferron['T_bearing'][order] - 273.15, revision, True)
        optical = read(args.evidence / f'optical_{revision}_profile.csv')
        draw(axes[2], optical['angle_deg'], optical['p'] / 1e5, revision, True)
    axes[0].set_ylabel(r'$p-p_{\mathrm{end}}$ (MPa)')
    axes[1].set_ylabel(r'$T_{\mathrm{bearing}}$ (°C)')
    axes[2].set_ylabel(r'$p$ (bar)')
    for ax in axes[:2]:
        ax.set_xlabel(r'$\varphi_{\mathrm{groove}}$ (deg)')
    axes[2].set_xlabel(r'$\varphi$ (deg)')
    for ax in axes:
        ax.set_xlim(0, 360)
        ax.set_xticks([0, 90, 180, 270, 360])
    axes[0].legend(loc='best', handlelength=2)
    save_figure(fig, str(args.output / 'steady_profiles'))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), layout='constrained')
    t = np.linspace(0, .015, 100)
    axes[0].plot(t * 1000, .05 * (100 * .02 / .0001)**2 * t / 1.9e6,
                 color=COLORS[2], lw=1.8, label='Exact at labelled time')
    axes[1].plot(t * 1000, t**2 * 1e6, color=COLORS[2], lw=1.8)
    for revision in REVISIONS:
        data = read(args.evidence / f'timing_{revision}_history.csv')
        draw(axes[0], data['reported_time_s'] * 1000, data['temperature_max_K'] - 300, revision, True, False)
        draw(axes[1], data['reported_time_s'] * 1000, data['bearing_z_m'] * 1e6, revision, True, False)
    axes[0].set_ylabel(r'$T-T_0$ (K)')
    axes[1].set_ylabel(r'$z_B$ (µm)')
    for ax in axes:
        ax.set_xlabel('Labelled time (ms)')
        ax.set_xlim(-.4, 15.4)
    axes[0].legend(loc='upper left', fontsize=12, handlelength=2)
    save_figure(fig, str(args.output / 'time_labels'))
    plt.close(fig)

    if (args.evidence / 'air_dt5us_history.csv').exists():
        fig, axes = plt.subplots(1, 2, figsize=(10.3, 4), layout='constrained')
        for i, (name, label) in enumerate((('current', r'$\Delta t=10$ µs'),
                                         ('dt5us', r'$\Delta t=5$ µs'),
                                         ('dt2p5us', r'$\Delta t=2.5$ µs'))):
            data = read(args.evidence / f'air_{name}_history.csv')
            time_key = 'comparison_time_s' if name == 'current' else 'time_s'
            data = data[data[time_key] >= .0008 - 1e-12]
            for ax, key, factor in ((axes[0], 'load_N', 1), (axes[1], 'alpha_gas_max', 100)):
                ax.plot(data[time_key] * 1000, data[key] * factor,
                        color=COLORS[i], ls=('-', '--', ':')[i], lw=1.6, label=label)
        axes[0].set_ylabel(r'$W$ (N)')
        axes[1].set_ylabel(r'$\alpha_{g,\max}$ (%)')
        axes[1].legend(loc='best', handlelength=2)
        for ax in axes:
            ax.set_xlabel(r'$t$ (ms)')
            ax.set_xlim(.8, 1)
            ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        save_figure(fig, str(args.output / 'air_timestep_sensitivity'))
        plt.close(fig)


if __name__ == '__main__':
    main()
