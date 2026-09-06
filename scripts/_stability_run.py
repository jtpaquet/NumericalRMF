"""Shared driver for the four alpha(t) / stability runs in this directory.

The four scripts differ only in two switches:

    scheme      'legacy' = the discretisation of nr32_fix.py
                'fixed'  = the second-order discretisation of rmf_solver.py
    rise_time   0 or 3 RMF periods

Everything else (lambda, the three gammas, Nr, dt, the number of periods) is
common and settable from the command line, so the four are directly comparable.

Each run writes results/alpha_<label>.pkl and figures/alpha_<label>.{pdf,png}.
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import paths
from rmf_solver import RMFPenetration, gamma_c, TWO_PI
from run_milroy import penetration_time

LAM = 11.07
GAMMAS = (14.9, 16.6, 18.2)
SNAPSHOTS = (5, 10, 20, 30, 40, 45, 50, 60)
COLORS = {14.9: 'blue', 16.6: 'green', 18.2: 'red'}
# Milroy 1999 Sec. III.A: (t of full penetration, alpha_s)
MILROY = {14.9: (None, 0.42), 16.6: (40, 0.98), 18.2: (27, 0.98)}


def run_case(gam, legacy, rise_time, Nr, dt, n_periods):
    """One gamma. Stops at the first non-finite alpha and reports the period."""
    sim = RMFPenetration(Nr=Nr, lam=LAM, gam=gam, legacy=legacy, rise_time=rise_time)
    nstep = int(round(TWO_PI/dt))
    snaps = {0: (sim.A.copy(), sim.b.copy(), sim.r.copy(), sim.tau)}
    ts, als, blew_up_at = [], [], None
    for p in range(n_periods):
        for _ in range(nstep):
            sim.step(dt)
        a = sim.alpha()
        if not np.isfinite(a):
            blew_up_at = p + 1
            break
        ts.append(p + 1.0)
        als.append(a)
        if (p + 1) in SNAPSHOTS:
            snaps[p + 1] = (sim.A.copy(), sim.b.copy(), sim.r.copy(), sim.tau)
    return np.array(ts), np.array(als), snaps, blew_up_at


def figure(results, label, title, fname):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from utils import plot_options

    plot_options['font.size'] = 22
    mpl.rcParams.update(plot_options)
    gc = gamma_c(LAM)

    fig, ax = plt.subplots(figsize=(13, 8))
    for gam in GAMMAS:
        res = results[gam]
        ax.plot(res['times'], res['alphas'], color=COLORS[gam], lw=2.2,
                label=fr'$\gamma$={gam} ($\gamma/\gamma_c$={gam/gc:.2f})')
        if res['blew_up_at'] is not None:
            ax.axvline(res['blew_up_at'], color=COLORS[gam], ls=':', lw=2)
        t_pen, a_s = MILROY[gam]
        ax.scatter([t_pen if t_pen else 80], [a_s], color=COLORS[gam], s=200,
                   marker='*', zorder=5, edgecolors='k', linewidths=0.6,
                   label='Milroy 1999' if gam == GAMMAS[0] else None)
    ax.axhline(1.0, color='k', ls='--', alpha=0.4)
    ax.set_xlabel('Time (RMF periods)')
    ax.set_ylabel(r'Penetration Factor $\alpha$')
    ax.set_title(title, fontsize=20)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.0)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 80)
    ax.set_ylim(0, 1.2)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        out = paths.figure(f'{fname}.{ext}')
        plt.savefig(out, bbox_inches='tight')
        print(f'  saved {out}')
    plt.close('all')


def main(label, legacy, rise_time, description):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--nr', type=int, default=64, help='radial points (default 64)')
    ap.add_argument('--dt', type=float, default=0.002, help='time step (default 0.002)')
    ap.add_argument('--periods', type=int, default=200,
                    help='RMF periods to run (default 200; gamma=14.9 needs ~200 '
                         'to settle)')
    ap.add_argument('--dt-scan', type=float, nargs='+', metavar='DT',
                    help='stability sweep: run each dt and report where it blows '
                         'up. Note the explicit scheme goes unstable *after* the '
                         'RMF has penetrated, not at start-up, so the sweep needs '
                         'enough periods to get past penetration - use at least '
                         '--periods 60. No pickle or figure is written.')
    ap.add_argument('--no-figure', action='store_true')
    args = ap.parse_args()

    scheme = 'legacy (nr32_fix.py)' if legacy else 'corrected (rmf_solver.py)'
    print(f'{description}\n')
    print(f'  scheme     : {scheme}')
    print(f'  rise_time  : {rise_time:.4g} '
          f'({rise_time/TWO_PI:.3g} RMF periods)')
    print(f'  lambda     : {LAM}   gamma_c = {gamma_c(LAM):.3f} (Milroy Eq. 15)')
    print(f'  Nr         : {args.nr}')

    if args.dt_scan:
        print(f'  periods    : {args.periods}\n')
        print('  (the explicit scheme fails after penetration, not at start-up,')
        print('   so a short sweep will report everything as stable)\n')
        print(f'{"dt":>9} {"gamma":>7} {"stable":>8} {"blew up at":>12} {"alpha_end":>10}')
        for dt in args.dt_scan:
            for gam in GAMMAS:
                _, a, _, bad = run_case(gam, legacy, rise_time, args.nr, dt,
                                        args.periods)
                print(f'{dt:9.5f} {gam:7.1f} {str(bad is None):>8} '
                      f'{"-" if bad is None else bad:>12} '
                      f'{a[-1] if len(a) else float("nan"):10.4f}', flush=True)
        return

    print(f'  dt         : {args.dt} ({int(round(TWO_PI/args.dt))} steps/period)')
    print(f'  periods    : {args.periods}\n')

    results, unstable = {}, False
    print(f'{"gamma":>7} {"gamma/gc":>9} {"alpha_s":>9} {"t_pen":>8} '
          f'{"Milroy alpha_s":>15} {"Milroy t_pen":>13}  status')
    for gam in GAMMAS:
        t, a, snaps, bad = run_case(gam, legacy, rise_time, args.nr, args.dt,
                                    args.periods)
        results[gam] = dict(times=t, alphas=a, snapshots=snaps, lam=LAM, gam=gam,
                            Nr=args.nr, dt=args.dt, rise_time=rise_time,
                            legacy=legacy, blew_up_at=bad)
        if bad is not None:
            unstable = True
            print(f'{gam:7.1f} {gam/gamma_c(LAM):9.3f} {"-":>9} {"-":>8} '
                  f'{"-":>15} {"-":>13}  BLEW UP at period {bad}')
            continue
        a_s = a[-1]
        t_pen = penetration_time(t, a, 0.95) if a_s > 0.9 else float('nan')
        m_t, m_a = MILROY[gam]
        print(f'{gam:7.1f} {gam/gamma_c(LAM):9.3f} {a_s:9.4f} {t_pen:8.1f} '
              f'{m_a:15.2f} {m_t if m_t else float("nan"):13.1f}  ok')

    out = paths.results(f'alpha_{label}.pkl')
    with open(out, 'wb') as f:
        pickle.dump(results, f)
    print(f'\n  wrote {out}')

    if not args.no_figure:
        title = (fr'{"original" if legacy else "corrected"} scheme, '
                 fr'$\tau_r$={rise_time/TWO_PI:.3g}$T$, '
                 fr'$N_r$={args.nr}, $dt$={args.dt}')
        figure(results, label, title, f'alpha_{label}')

    if unstable:
        print('\n  at least one gamma went unstable; try a smaller --dt, or '
              '--dt-scan to find the limit')
