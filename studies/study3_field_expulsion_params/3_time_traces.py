#!/usr/bin/env python3
"""Study 3 - the actual "outside minus centre" field traces vs real time,
for each T_e in the reference runs, against the observed 10/25 ms timescale.

    figures/study3_time_traces.pdf
"""
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import paths
import device as dev
from utils import plot_options

TARGET_TIMES_MS = (10.0, 25.0)
MEASURED_GAP_RATIO = (dev.MEASURED_B_WALL_G - dev.MEASURED_B_CENTER_G)/dev.B_W_G


def main():
    plot_options['font.size'] = 16
    mpl.rcParams.update(plot_options)

    with open(paths.results('study3_reference_by_Te.pkl'), 'rb') as f:
        reference = pickle.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    cmap = plt.get_cmap('plasma')
    Tes = sorted(reference)
    for i, Te in enumerate(Tes):
        res = reference[Te]
        c = cmap(i/max(1, len(Tes) - 1))
        gap_r = np.abs(res['BrR'] - res['Br0'])
        gap_th = np.abs(res['BthR'] - res['Bth0'])
        axes[0].plot(res['t_ms'], gap_r, color=c, lw=2,
                     label=fr'$T_e$={Te:.1f} eV ($\lambda$={res["lam"]:.1f})')
        axes[1].plot(res['t_ms'], gap_th, color=c, lw=2,
                     label=fr'$T_e$={Te:.1f} eV ($\lambda$={res["lam"]:.1f})')

    for ax, ylabel in zip(axes, [r'$|B_r(R) - B_r(0)|\,/\,B_w$',
                                  r'$|B_\theta(R) - B_\theta(0)|\,/\,B_w$']):
        for tt in TARGET_TIMES_MS:
            ax.axvline(tt, color='k', ls='--', alpha=0.5)
        ax.axhline(MEASURED_GAP_RATIO, color='r', ls=':', lw=2,
                   label='measured gap (with plasma)' if ax is axes[0] else None)
        ax.set_xlabel('time (ms)')
        ax.set_ylabel(ylabel)
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
    axes[0].legend(fontsize=11, loc='lower right')
    fig.suptitle('Field-expulsion transient vs. the observed 10 ms / 25 ms marks\n'
                 f'(dashed vertical: 10, 25 ms; dotted red: measured '
                 f'({dev.MEASURED_B_WALL_G:.1f}-{dev.MEASURED_B_CENTER_G:.1f}) G '
                 f'/ {dev.B_W_G:.0f} G = {MEASURED_GAP_RATIO:.2f})')
    plt.tight_layout()
    out = paths.figure('study3_time_traces.pdf')
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')

    print(f'\n{"Te(eV)":>7} {"lambda":>7} {"t_settle(ms)":>13} {"final gap":>10}')
    for Te in Tes:
        res = reference[Te]
        gap = abs(res['BrR'][-1] - res['Br0'][-1])
        print(f'{Te:7.2f} {res["lam"]:7.2f} {res["t_settle_ms"]:13.4f} {gap:10.4f}')
    plt.close('all')


if __name__ == '__main__':
    main()
