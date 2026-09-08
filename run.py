#!/usr/bin/env python3
"""Run one RMF-penetration case at any gamma/lambda/Nr/dt.

    python run.py --gamma 16.6 --lam 11.07 --nr 32 --dt 0.001

See run_milroy.py for Milroy's three reference gammas together at his
production resolution, and studies/ for the comparison and validation work
built on top of this solver.
"""
import argparse
import pickle

import numpy as np

import paths
from rmf_solver import RMFPenetration, gamma_c, TWO_PI
from run_milroy import penetration_time


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gamma', type=float, default=16.6, help='RMF strength (default 16.6)')
    ap.add_argument('--lam', type=float, default=11.07, help='R/skin depth (default 11.07)')
    ap.add_argument('--nr', type=int, default=32, help='radial points (default 32)')
    ap.add_argument('--dt', type=float, default=0.002, help='time step (default 0.002)')
    ap.add_argument('--periods', type=int, default=200,
                     help='RMF periods to run (default 200)')
    ap.add_argument('--rise-time', type=float, default=0.0,
                     help='turn-on ramp, in RMF periods (default 0, Milroy step; '
                          'Hugrass & Grimm 1981 use 0.32)')
    ap.add_argument('--bz0', type=float, default=0.0, help='external DC bias field (default 0)')
    ap.add_argument('--save', metavar='NAME', help='also write results/NAME.pkl')
    args = ap.parse_args()

    gc = gamma_c(args.lam)
    print(f'lambda = {args.lam}   gamma = {args.gamma}   '
          f'gamma/gamma_c = {args.gamma/gc:.3f} (gamma_c = {gc:.3f}, Milroy Eq. 15)')
    print(f'Nr = {args.nr}   dt = {args.dt} ({int(round(TWO_PI/args.dt))} steps/period)   '
          f'periods = {args.periods}\n')

    sim = RMFPenetration(Nr=args.nr, lam=args.lam, gam=args.gamma,
                          rise_time=args.rise_time*TWO_PI, bz0=args.bz0)
    t, a, snaps = sim.run(n_periods=args.periods, dt=args.dt)

    if not np.all(np.isfinite(a)):
        raise RuntimeError(f'went unstable; reduce --dt below {args.dt}')

    tp = penetration_time(t, a) if a[-1] > 0.9 else np.nan
    print(f'alpha_s = {a[-1]:.4f}   t_pen (95%) = {tp:.1f} T')

    if args.save:
        out = paths.results(f'{args.save}.pkl')
        with open(out, 'wb') as f:
            pickle.dump(dict(times=t, alphas=a, snapshots=snaps, lam=args.lam,
                              gam=args.gamma, Nr=args.nr, dt=args.dt), f)
        print(f'wrote {out}')


if __name__ == '__main__':
    main()
