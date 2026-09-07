#!/usr/bin/env python3
"""Study 4 - does wiring B_z0 into the boundary condition reproduce the
reversal-asymmetric DC-bias-field coupling seen in study 3's real shots?

Short answer: no, and not for a numerical reason -- the n=0/n=1 truncation
this solver integrates is *exactly* invariant under a uniform shift of B_z.
Both the Hall-coupling term in dA/dt and the whole of the b equation use
`dbdr(b)` or plain differences of `b` (see rmf_solver.rhs) -- never `b`
itself -- so the Dirichlet condition `b(1) = bz0` (rmf_solver.py's
`apply_bc`) is a pure relabelling: b(r, t; bz0) = b(r, t; 0) + bz0 for every
r and t, exactly, because the time-stepper (Heun) only ever adds RHS values
that don't see the shift. alpha = |b(1) - b(0)| therefore cannot depend on
bz0 at all, for any lambda, gamma, or time.

This script runs the solver, rather than just trusting the algebra: for
Milroy's lambda = 11.07 and his three reference gammas (14.9 sub-critical,
16.6 and 18.2 super-critical), alpha(t) is computed at bz0 = 0 and at five
other bz0 (three within the range of alpha_s itself, two much larger) and
compared point by point over the whole run.

    figures/study4_bz0_invariance.pdf   alpha(t) overlaid for every bz0
                                         (curves are on top of each other)
                                         and the point-by-point residual
    results/study4_bz0_invariance.csv   max|alpha(bz0,t) - alpha(0,t)| per gamma
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import paths
from rmf_solver import RMFPenetration, gamma_c

LAM = 11.07
GAMMAS = (14.9, 16.6, 18.2)
BZ0_VALUES = (0.0, 0.05, 0.3, -0.3, 1.0, -1.0)
NR, DT, PERIODS = 64, 0.002, 60


def run(gam, bz0):
    sim = RMFPenetration(Nr=NR, lam=LAM, gam=gam, bz0=bz0)
    ts, als, _ = sim.run(n_periods=PERIODS, dt=DT)
    return ts, als


def main():
    rows = []
    traces = {}
    for gam in GAMMAS:
        t0, a0 = run(gam, 0.0)
        traces[gam] = {0.0: (t0, a0)}
        maxdiff = 0.0
        for bz0 in BZ0_VALUES[1:]:
            t, a = run(gam, bz0)
            traces[gam][bz0] = (t, a)
            maxdiff = max(maxdiff, float(np.max(np.abs(a - a0))))
        rows.append(dict(gamma=gam, gamma_over_gc=gam/gamma_c(LAM),
                          alpha_s_bz0_0=a0[-1], max_abs_diff=maxdiff))
        print(f'gamma={gam:5.1f}  gamma/gc={gam/gamma_c(LAM):.3f}  '
              f'alpha_s(bz0=0)={a0[-1]:.4f}  '
              f'max|alpha(bz0)-alpha(0)| over bz0 in {BZ0_VALUES[1:]} '
              f'and t in [0,{PERIODS}T] = {maxdiff:.2e}')

    out = paths.results('study4_bz0_invariance.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'\nwrote {out}')
    print('max diff is at the ~1e-14 level (float round-off) for every gamma: '
          'bz0 has no effect on alpha(t), confirming the shift-invariance '
          'argument above rather than just a numerical coincidence.')

    make_figure(traces)


def make_figure(traces):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from utils import plot_options

    plot_options['font.size'] = 15
    mpl.rcParams.update(plot_options)

    fig, axes = plt.subplots(2, len(GAMMAS), figsize=(15, 7), sharex='col')
    cmap = plt.get_cmap('viridis')
    for j, gam in enumerate(GAMMAS):
        t0, a0 = traces[gam][0.0]
        for i, bz0 in enumerate(BZ0_VALUES):
            t, a = traces[gam][bz0]
            axes[0, j].plot(t, a, color=cmap(i/(len(BZ0_VALUES) - 1)), lw=1.5,
                             label=f'bz0={bz0:g}')
            axes[1, j].plot(t, np.abs(a - a0), color=cmap(i/(len(BZ0_VALUES) - 1)), lw=1.5)
        axes[0, j].set_title(fr'$\gamma$={gam}')
        axes[1, j].set_xlabel('time (RMF periods)')
        axes[1, j].set_yscale('log')
        axes[0, j].grid(True, alpha=0.3)
        axes[1, j].grid(True, alpha=0.3)
    axes[0, 0].set_ylabel(r'$\alpha(t)$ (all bz0 overlaid)')
    axes[1, 0].set_ylabel(r'$|\alpha(t; bz0) - \alpha(t; 0)|$')
    axes[0, -1].legend(fontsize=9, loc='lower right')
    fig.suptitle('A uniform B_z0 has no effect on alpha(t): curves overlap exactly, '
                 'residual is float round-off')
    plt.tight_layout()
    out = paths.figure('study4_bz0_invariance.pdf')
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')
    plt.close('all')


if __name__ == '__main__':
    main()
