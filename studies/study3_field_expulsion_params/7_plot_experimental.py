#!/usr/bin/env python3
"""Study 3 - plots for 6_load_experimental_data.py.

    figures/study3_experimental_Te.pdf    inferred T_e vs pressure and vs
                                           added DC-coil resistance (weaker
                                           bias field at higher resistance)
    figures/study3_experimental_time.pdf  the 90%-of-plateau coupling time,
                                           same two slices
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import matplotlib as mpl
import matplotlib.pyplot as plt

import paths
from utils import plot_options


def load(fname='study3_experimental_summary.csv'):
    with open(paths.results(fname)) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in ('t_on', 'ratio_x', 'ratio_z', 'dx_final', 'dz_final',
                 't_couple90_ms', 'Te_x_eV', 'Te_z_eV'):
            r[k] = float(r[k])
        r['ohm'] = None if r['ohm'] == '' else int(r['ohm'])
        r['p_mtorr'] = int(r['p_mtorr'])
        r['reversed'] = r['reversed'] == 'True'
        r['noDC'] = r['noDC'] == 'True'
        r['Te_x_at_floor'] = r['Te_x_at_floor'] == 'True'
    return rows


def main():
    plot_options['font.size'] = 15
    mpl.rcParams.update(plot_options)
    rows = load()

    pscan = sorted((r for r in rows if r['ohm'] == 7 and not r['reversed'] and not r['noDC']),
                   key=lambda r: r['p_mtorr'])
    oscan = sorted((r for r in rows if r['p_mtorr'] == 22), key=lambda r: (r['ohm'] is None, r['ohm']))

    for metric, ylabel, fname in [
        ('Te_x_eV', r'inferred $T_e$ (eV), linear-regime model', 'study3_experimental_Te.pdf'),
        ('t_couple90_ms', 'time to 90% of plateau (ms)', 'study3_experimental_time.pdf'),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        floor = [r.get('Te_x_at_floor', False) for r in pscan] if metric == 'Te_x_eV' else [False]*len(pscan)
        axes[0].plot([r['p_mtorr'] for r in pscan], [r[metric] for r in pscan], 'o-', color='#2b6cb0')
        for r, fl in zip(pscan, floor):
            if fl:
                axes[0].plot(r['p_mtorr'], r[metric], 'o', color='crimson', zorder=5)
        axes[0].set_xlabel('fill pressure (mTorr)')
        axes[0].set_title('ohm = 7 (fixed bias field)')
        ohm_x = [(-1 if r['ohm'] is None else r['ohm']) for r in oscan]
        axes[1].plot(ohm_x, [r[metric] for r in oscan], 'o-', color='#2b6cb0')
        axes[1].set_xlabel('added DC-coil resistance ($\\Omega$, -1 = noDC)')
        axes[1].set_title('p = 22 mTorr (fixed pressure)')
        for ax in axes:
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)
        fig.suptitle('Experimental shots: ' + ylabel)
        plt.tight_layout()
        out = paths.figure(fname)
        plt.savefig(out, bbox_inches='tight')
        print(f'saved {out}')
    plt.close('all')


if __name__ == '__main__':
    main()
