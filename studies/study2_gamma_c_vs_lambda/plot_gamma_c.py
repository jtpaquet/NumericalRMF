#!/usr/bin/env python3
"""Plot results/gamma_c_vs_lambda.csv against Milroy 1999 Fig. 4.

Two panels: gamma_c vs lambda with Milroy's Eqs. (14) and (15) overlaid, and the
ratio to Eq. (15), which is where any systematic disagreement shows up.

    python plot_gamma_c.py
    python plot_gamma_c.py --csv gamma_c_vs_lambda.csv --log
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import paths
from rmf_solver import gamma_c as gamma_c_milroy
from utils import plot_options


def load(path):
    with open(path) as f:
        rows = [r for r in csv.DictReader(f)]
    rows.sort(key=lambda r: float(r['lam']))
    lam = np.array([float(r['lam']) for r in rows])
    gc = np.array([float(r['gamma_c']) for r in rows])
    lo = np.array([float(r['gamma_c_lo']) for r in rows])
    hi = np.array([float(r['gamma_c_hi']) for r in rows])
    ok = np.isfinite(gc)
    return lam[ok], gc[ok], lo[ok], hi[ok]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--csv', default='gamma_c_vs_lambda.csv')
    ap.add_argument('--log', action='store_true',
                    help='log-log axes for the left panel instead of the '
                         'linear gamma_c/lambda vs lambda axes Milroy uses '
                         'in his Fig. 4')
    ap.add_argument('--out', default='gamma_c_vs_lambda')
    args = ap.parse_args()

    path = args.csv if os.path.isabs(args.csv) else paths.results(args.csv)
    if not os.path.exists(path):
        sys.exit(f'no {path} - run gamma_c_scan.py first')
    lam, gc, lo, hi = load(path)
    if len(lam) == 0:
        sys.exit(f'{path} has no usable rows')

    plot_options['font.size'] = 22
    mpl.rcParams.update(plot_options)

    if args.log:
        lam_f = np.logspace(np.log10(max(0.3, lam.min()*0.8)),
                            np.log10(lam.max()*1.25), 400)
    else:
        lam_f = np.linspace(max(0.0, lam.min()*0.8), lam.max()*1.05, 400)
    eq = np.array([gamma_c_milroy(x) for x in lam_f])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(18, 7.5))

    if args.log:
        # gamma_c vs lambda, log-log
        ax.plot(lam_f, 1.12*lam_f, 'k--', lw=1.8,
                label=r'Milroy Eq. (14), $1.12\lambda$ (expulsion)')
        ax.plot(lam_f, eq, 'k-', lw=2.2, label='Milroy Eq. (15) (penetration)')
        ax.errorbar(lam, gc, yerr=[gc - lo, hi - gc], fmt='o', color='tab:blue',
                    ms=8, lw=1.5, capsize=4, label='this code')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$\lambda = R/\delta$')
        ax.set_ylabel(r'$\gamma_c$')
    else:
        # gamma_c/lambda vs lambda, linear - matches Milroy 1999 Fig. 4
        with np.errstate(divide='ignore', invalid='ignore'):
            eq_ratio = np.where(lam_f > 0, eq/np.where(lam_f > 0, lam_f, 1), np.nan)
            eq14_ratio = np.full_like(lam_f, 1.12)
        ax.plot(lam_f, eq14_ratio, 'k--', lw=1.8,
                label=r'Milroy Eq. (14), $1.12$ (expulsion)')
        ax.plot(lam_f, eq_ratio, 'k-', lw=2.2, label='Milroy Eq. (15) (penetration)')
        ax.errorbar(lam, gc/lam, yerr=[(gc - lo)/lam, (hi - gc)/lam],
                    fmt='D-', color='tab:blue', ms=8, lw=1.5, capsize=4,
                    label='this code')
        ax.set_xlim(0, lam.max()*1.05)
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$\gamma_c/\lambda$')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=17, loc='upper left')

    ratio = gc/np.array([gamma_c_milroy(x) for x in lam])
    ax2.axhline(1.0, color='k', lw=2)
    ax2.axvline(6.5, color='gray', ls=':', lw=1.5)
    ax2.plot(lam, ratio, 'o-', color='tab:blue', ms=8, lw=1.5)
    if args.log:
        ax2.set_xscale('log')
    else:
        ax2.set_xlim(0, lam.max()*1.05)
    ax2.set_xlabel(r'$\lambda = R/\delta$')
    ax2.set_ylabel(r'$\gamma_c$ / Milroy Eq. (15)')
    ax2.grid(True, which='both', alpha=0.3)
    ax2.text(6.6, ax2.get_ylim()[1], r'  $\lambda=6.5$', va='top', fontsize=16,
             color='gray')

    fig.suptitle(r'Critical $\gamma$ for full RMF penetration '
                 '(Milroy 1999, Fig. 4)', fontsize=24)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        out = paths.figure(f'{args.out}.{ext}')
        plt.savefig(out, bbox_inches='tight')
        print(f'saved {out}')

    print(f'\n{"lambda":>8} {"gamma_c":>10} {"Milroy":>9} {"ratio":>7}')
    for x, g in zip(lam, gc):
        e = gamma_c_milroy(x)
        print(f'{x:8.1f} {g:10.3f} {e:9.3f} {g/e:7.3f}')
    print(f'\nmean ratio {ratio.mean():.3f}, spread '
          f'{ratio.min():.3f} to {ratio.max():.3f}')


if __name__ == '__main__':
    main()
