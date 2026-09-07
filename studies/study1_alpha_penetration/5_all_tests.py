#!/usr/bin/env python3
"""Run all four scheme x ramp combinations and overlay them, one panel per gamma.

Same settings for all four, so the panels isolate the two effects:

    down a column   change of turn-on ramp at fixed discretisation
    across a row    change of discretisation at fixed ramp

Milroy's Sec. III.A points are marked on each panel.

    python 5_all_tests.py
    python 5_all_tests.py --nr 64 --periods 120
"""
import argparse
import os
import pickle
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                  # for _stability_run
sys.path.insert(0, os.path.dirname(os.path.dirname(_HERE)))  # repo root

import numpy as np

import paths
from _stability_run import (GAMMAS, LAM, MILROY, TWO_PI, add_common_args,
                            execute, gamma_c, penetration_time)

CASES = [
    ('legacy_ramp0', True,  0.0,        'original scheme, no ramp',   'tab:orange', '--'),
    ('legacy_ramp3', True,  3*TWO_PI,   'original scheme, 3T ramp',   'tab:red',    ':'),
    ('fixed_ramp0',  False, 0.0,        'corrected scheme, no ramp',  'tab:blue',   '-'),
    ('fixed_ramp3',  False, 3*TWO_PI,   'corrected scheme, 3T ramp',  'tab:green',  '-.'),
]


def figure(all_results, nr, dt, fname='alpha_all_tests'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from utils import plot_options

    plot_options['font.size'] = 20
    mpl.rcParams.update(plot_options)

    fig, axes = plt.subplots(1, len(GAMMAS), figsize=(19, 7), sharey=True)
    for ax, gam in zip(np.atleast_1d(axes), GAMMAS):
        for label, _, _, desc, color, ls in CASES:
            res = all_results[label][gam]
            ax.plot(res['times'], res['alphas'], color=color, ls=ls, lw=2.2,
                    label=desc if gam == GAMMAS[0] else None)
            if res['blew_up_at'] is not None:
                ax.axvline(res['blew_up_at'], color=color, ls=':', lw=1.5)
        t_pen, a_s = MILROY[gam]
        ax.scatter([t_pen if t_pen else 78], [a_s], color='k', s=260, marker='*',
                   zorder=5, label='Milroy 1999' if gam == GAMMAS[0] else None)
        ax.axhline(1.0, color='k', ls='--', alpha=0.35)
        ax.set_title(fr'$\gamma$={gam} ($\gamma/\gamma_c$={gam/gamma_c(LAM):.2f})',
                     fontsize=20)
        ax.set_xlabel('Time (RMF periods)')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 80)
        ax.set_ylim(0, 1.2)
    np.atleast_1d(axes)[0].set_ylabel(r'Penetration Factor $\alpha$')
    handles, labels = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=17,
               bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(fr'$\lambda$={LAM}, $N_r$={nr}, $dt$={dt}', fontsize=22)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        out = paths.figure(f'{fname}.{ext}')
        plt.savefig(out, bbox_inches='tight')
        print(f'  saved {out}')
    plt.close('all')


def main():
    ap = add_common_args(argparse.ArgumentParser(description=__doc__.splitlines()[0]))
    args = ap.parse_args()

    print('All four scheme x ramp combinations\n')
    print(f'  lambda  : {LAM}   gamma_c = {gamma_c(LAM):.3f} (Milroy Eq. 15)')
    print(f'  Nr      : {args.nr}')
    print(f'  dt      : {args.dt} ({int(round(TWO_PI/args.dt))} steps/period)')
    print(f'  periods : {args.periods}')

    all_results = {}
    for label, legacy, rise_time, desc, _, _ in CASES:
        print(f'\n--- {desc} ---')
        all_results[label] = execute(label, legacy, rise_time, args.nr, args.dt,
                                     args.periods, make_figure=False, save=True)

    out = paths.results('alpha_all_tests.pkl')
    with open(out, 'wb') as f:
        pickle.dump(all_results, f)
    print(f'\n  wrote {out}')

    # summary table: t_pen for the two penetrating gammas, against Milroy
    print('\nsummary (t_pen at 95% of alpha_s, alpha_s in brackets)\n')
    head = f'{"case":28}' + ''.join(f'{f"gamma={g}":>22}' for g in GAMMAS)
    print(head)
    for label, _, _, desc, _, _ in CASES:
        row = f'{desc:28}'
        for gam in GAMMAS:
            res = all_results[label][gam]
            a = res['alphas']
            if res['blew_up_at'] is not None:
                row += f'{"blew up":>22}'
            elif a[-1] > 0.9:
                tp = penetration_time(res['times'], a, 0.95)
                row += f'{f"{tp:.1f} T ({a[-1]:.3f})":>22}'
            else:
                row += f'{f"-- ({a[-1]:.3f})":>22}'
        print(row)
    row = f'{"Milroy 1999":28}'
    for gam in GAMMAS:
        t_pen, a_s = MILROY[gam]
        row += f'{f"{t_pen} T ({a_s})" if t_pen else f"-- ({a_s})":>22}'
    print(row)

    if not args.no_figure:
        figure(all_results, args.nr, args.dt)


if __name__ == '__main__':
    main()
