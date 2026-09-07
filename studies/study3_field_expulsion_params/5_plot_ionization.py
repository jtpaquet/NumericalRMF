#!/usr/bin/env python3
"""Study 3 extension - plots for 4_ionization_fraction.py.

    figures/study3_ionization_tau.pdf     tau_penetration vs f_ion, one panel
                                           per T_e, the 10-... wait: 1-50 ms
                                           target band shaded, the 7 requested
                                           f_ion marked
    figures/study3_ionization_window.pdf  critical f_ion(p) per T_e and the
                                           band that gives 1-50 ms, against
                                           the requested f_ion grid
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import paths
import device as dev
from utils import plot_options

TARGET_MS = (1.0, 50.0)


def load(fname):
    with open(paths.results(fname)) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            r[k] = float(v)
    return rows


def main():
    plot_options['font.size'] = 15
    mpl.rcParams.update(plot_options)

    grid = load('study3_ionization_grid.csv')
    window = load('study3_ionization_window.csv')
    Tes = sorted(set(r['Te_eV'] for r in grid))
    ps = sorted(set(r['p_mtorr'] for r in grid))
    cmap = plt.get_cmap('viridis')

    # ---------------------------------------------------- tau vs f_ion, per Te
    fig, axes = plt.subplots(1, len(Tes), figsize=(4.2*len(Tes), 5), sharey=True)
    for ax, Te in zip(axes, Tes):
        ax.axhspan(*TARGET_MS, color='crimson', alpha=0.15,
                   label='observed 1-50 ms' if Te == Tes[0] else None)
        for j, p in enumerate(ps):
            rows = sorted((r for r in grid if r['Te_eV'] == Te and r['p_mtorr'] == p),
                          key=lambda r: r['f_ion'])
            f = [r['f_ion'] for r in rows]
            tau = [r['t_pen_ms'] for r in rows]
            ax.plot(f, tau, 'o-', color=cmap(j/max(1, len(ps) - 1)),
                    label=f'{p:.0f} mTorr' if Te == Tes[0] else None)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$f_{ion}$')
        ax.set_title(fr'$T_e$={Te:.1f} eV')
        ax.grid(True, alpha=0.3, which='both')
    axes[0].set_ylabel(r'$t_{pen}$ (ms), Milroy Eq. (17)')
    fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    fig.suptitle('Nonlinear penetration time vs. ionisation fraction: none of the '
                 '7 requested values land in the observed band')
    plt.tight_layout()
    out = paths.figure('study3_ionization_tau.pdf')
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')

    # --------------------------------------------------- critical f_ion window
    fig, axes = plt.subplots(1, len(Tes), figsize=(4.2*len(Tes), 5), sharey=True)
    for ax, Te in zip(axes, Tes):
        rows = sorted((r for r in window if r['Te_eV'] == Te), key=lambda r: r['p_mtorr'])
        p = [r['p_mtorr'] for r in rows]
        fc = [r['f_crit'] for r in rows]
        f1 = [r['f_1ms'] for r in rows]
        f50 = [r['f_50ms'] for r in rows]
        ax.fill_between(p, f50, f1, color='crimson', alpha=0.3,
                        label='band -> 1-50 ms' if Te == Tes[0] else None)
        ax.plot(p, fc, 'k--', lw=1.5, label=r'$f_{crit}$ ($\gamma=\gamma_c$)'
                if Te == Tes[0] else None)
        for fv in dev.F_ION:
            ax.axhline(fv, color='gray', lw=0.6, alpha=0.6)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('fill pressure (mTorr)')
        ax.set_title(fr'$T_e$={Te:.1f} eV')
        ax.set_ylim(1e-6, 2)
        ax.grid(True, alpha=0.3, which='major')
    axes[0].set_ylabel(r'$f_{ion}$')
    fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    fig.suptitle('Ionisation fraction at threshold and the narrow band above it '
                 'giving 1-50 ms\n(thin gray lines: the 7 requested $f_{ion}$ values)')
    plt.tight_layout()
    out = paths.figure('study3_ionization_window.pdf')
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')
    plt.close('all')


if __name__ == '__main__':
    main()
