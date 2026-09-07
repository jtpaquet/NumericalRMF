#!/usr/bin/env python3
"""Study 2 - reproduce Milroy 1999 Fig. 4: the critical gamma for full RMF
penetration as a function of lambda.

Milroy's answer, for comparison:

    gamma_c = 1.12 lambda                                    lambda <= 6.5   Eq. (14)
    gamma_c = 1.12 lambda (1 + 0.12 (lambda - 6.5)^0.4)       lambda >  6.5   Eq. (15)

Method
------
For each lambda, bisect on gamma. A gamma "penetrates" if a cold start
(A = b = 0, RMF switched on as a step at t = 0 - the corrected scheme with
rise_time = 0) drives alpha above --alpha-pen within the time budget.

Milroy finds gamma_c with a staircase: hold gamma until alpha settles, step up,
repeat. Bisection from cold starts gives the same threshold far more cheaply,
because below threshold the cold start and the staircase land on the same low-alpha
branch - the hysteresis only matters coming *down* from the penetrated branch,
which is the expulsion threshold (Eq. 14), not this one.

The time budget is the whole difficulty. Milroy's Eq. (17), tau_P =
lambda^2 / (2 sqrt(gamma_N)) with gamma_N = (gamma - gamma_c)/gamma_c, says the
penetration time diverges at threshold, so a finite budget always reads gamma_c
slightly high. Setting the budget to f * lambda^2 RMF periods makes that bias
lambda-independent:

    gamma_N_min = (1 / (4 pi f))^2      ->   f = 1 gives ~0.6% high

Cost
----
Cost per lambda scales as roughly lambda^4, which is brutal:

    Nr    ~ 8 lambda        (to resolve the skin depth delta = 1/lambda)
    dt    ~ dr^2            (the Hall/whistler limit; see TODO 7b)
    T_max ~ lambda^2        (the penetration time)
    steps ~ T_max/dt ~ lambda^4

lambda <= 20 is comfortable, 30-40 is hours, and lambda >= 50 is out of reach
with this explicit integrator - that is exactly why Milroy has a semi-implicit
option and why Hugrass 1985 says his approximate equations let him "obtain
solutions for large values of R/delta in a reasonable amount of computing time".
Run with --dry-run first to see the plan and the step counts before committing.

Results are appended to results/gamma_c_vs_lambda.csv as they are produced, and
--resume skips lambdas already in the file, so the scan can be stopped and
restarted, or run small-lambda-first and abandoned when it gets too slow.

    python gamma_c_scan.py --dry-run
    python gamma_c_scan.py --lam-max 20
    python gamma_c_scan.py --resume
"""
import argparse
import csv
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import paths
from rmf_solver import RMFPenetration, gamma_c as gamma_c_milroy, TWO_PI

# lambda <= 6.5 in steps of 0.5, then the coarser list above it
LAMBDAS = [round(0.5*i, 1) for i in range(2, 14)] + [
    7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
    25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

CSV_FIELDS = ['lam', 'gamma_c', 'gamma_c_lo', 'gamma_c_hi', 'milroy_eq',
              'ratio', 'alpha_s', 't_pen', 'Nr', 'dt', 'max_periods',
              'n_evals', 'seconds']


# ----------------------------------------------------------------- resolution
def choose_Nr(lam, ppd, nr_min, nr_max):
    """Radial points: ppd cells across the skin depth delta = 1/lambda."""
    return int(min(nr_max, max(nr_min, math.ceil(ppd*lam) + 1)))


def choose_dt(lam, Nr, hall_c, diff_c):
    """Explicit stability limit, the tighter of two constraints.

    Diffusion:  dt <= lambda^2 dr^2.
    Hall:       the whistler branch has omega = k^2 b / 2 in these units, with
                b -> 1 once penetrated, so dt ~ dr^2 independent of lambda.
                Measured: lambda = 11.07, Nr = 64 (dr^2 = 2.5e-4) is stable at
                dt = 0.002 (~8 dr^2) and fails at 0.003, hence hall_c = 5.
    """
    dr = 1.0/(Nr - 1)
    return min(hall_c*dr**2, diff_c*lam**2*dr**2)


def choose_periods(lam, factor, floor):
    return int(max(floor, math.ceil(factor*lam**2)))


# --------------------------------------------------------------- one gamma
def penetrates(lam, gam, Nr, dt, max_periods, alpha_pen, plateau, quiet=True):
    """Cold start. Returns (verdict, alpha_s, t_pen_periods, periods_run).

    verdict is True (reached alpha_pen), False (did not), or None (went
    non-finite). Exits early on penetration, and on a flat sub-threshold
    plateau, so most evaluations cost far less than the full budget.
    """
    sim = RMFPenetration(Nr=Nr, lam=lam, gam=gam)      # corrected, rise_time = 0
    nstep = int(round(TWO_PI/dt))
    alphas = []
    for p in range(max_periods):
        for _ in range(nstep):
            sim.step(dt)
        a = sim.alpha()
        if not np.isfinite(a):
            return None, float('nan'), float('nan'), p + 1
        alphas.append(a)
        if a >= alpha_pen:
            return True, a, float(p + 1), p + 1
        # settled well below threshold: no point running the rest of the budget
        if len(alphas) > plateau and a < 0.7*alpha_pen:
            window = alphas[-plateau:]
            if max(window) - min(window) < 1e-5*max(1e-12, abs(a)) * plateau:
                return False, a, float('nan'), p + 1
    return False, alphas[-1], float('nan'), max_periods


def refine_t_pen(lam, gam, Nr, dt, max_periods, pen_frac):
    """Re-run a penetrating case to get t_pen at pen_frac of the final alpha.

    Done separately because the bisection exits the moment alpha crosses
    alpha_pen, which is before alpha_s is known.
    """
    sim = RMFPenetration(Nr=Nr, lam=lam, gam=gam)
    nstep = int(round(TWO_PI/dt))
    ts, als = [], []
    for p in range(max_periods):
        for _ in range(nstep):
            sim.step(dt)
        a = sim.alpha()
        if not np.isfinite(a):
            break
        ts.append(p + 1.0)
        als.append(a)
    if not als:
        return float('nan'), float('nan')
    t, a = np.array(ts), np.array(als)
    thr = pen_frac*a[-1]
    if a.max() < thr:
        return a[-1], float('nan')
    i = int(np.argmax(a >= thr))
    t_pen = t[0] if i == 0 else t[i-1] + (thr - a[i-1])*(t[i]-t[i-1])/(a[i]-a[i-1])
    return a[-1], t_pen


# ------------------------------------------------------------------ bisection
def find_gamma_c(lam, Nr, dt, max_periods, args):
    """Bracket then bisect. Returns (lo, hi, n_evals) or (nan, nan, n) on failure."""
    g0 = gamma_c_milroy(lam)
    n = 0

    hi = None
    for mult in (1.15, 1.4, 1.8, 2.5, 4.0):
        v, _, _, _ = penetrates(lam, g0*mult, Nr, dt, max_periods,
                                args.alpha_pen, args.plateau)
        n += 1
        if v is None:
            print(f'    unstable at gamma={g0*mult:.3f}; reduce --hall-c')
            return float('nan'), float('nan'), n
        if v:
            hi = g0*mult
            break
    if hi is None:
        print(f'    no penetration up to gamma={g0*4:.3f}; budget too short?')
        return float('nan'), float('nan'), n

    lo = None
    for mult in (0.85, 0.7, 0.5, 0.3, 0.15):
        if g0*mult >= hi:
            continue
        v, _, _, _ = penetrates(lam, g0*mult, Nr, dt, max_periods,
                                args.alpha_pen, args.plateau)
        n += 1
        if v is None:
            print(f'    unstable at gamma={g0*mult:.3f}; reduce --hall-c')
            return float('nan'), float('nan'), n
        if not v:
            lo = g0*mult
            break
    if lo is None:
        lo = 1e-3*g0

    while (hi - lo)/lo > args.tol and n < args.max_evals:
        mid = 0.5*(lo + hi)
        v, _, _, _ = penetrates(lam, mid, Nr, dt, max_periods, args.alpha_pen,
                                args.plateau)
        n += 1
        if v is None:
            print(f'    unstable at gamma={mid:.3f}; reduce --hall-c')
            return float('nan'), float('nan'), n
        if v:
            hi = mid
        else:
            lo = mid
    return lo, hi, n


# ----------------------------------------------------------------------- main
def plan(lam, args):
    Nr = choose_Nr(lam, args.ppd, args.nr_min, args.nr_max)
    dt = choose_dt(lam, Nr, args.hall_c, args.diff_c)
    periods = choose_periods(lam, args.tmax_factor, args.min_periods)
    steps = periods*int(round(TWO_PI/dt))
    return Nr, dt, periods, steps


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam', type=float, nargs='+',
                    help='only these lambdas (default: the full list)')
    ap.add_argument('--lam-max', type=float, default=None,
                    help='skip lambdas above this (the cost goes as lambda^4)')
    ap.add_argument('--out', default='gamma_c_vs_lambda.csv')
    ap.add_argument('--resume', action='store_true',
                    help='skip lambdas already present in the output file')
    ap.add_argument('--dry-run', action='store_true',
                    help='print the plan and the step counts, run nothing')

    ap.add_argument('--alpha-pen', type=float, default=0.9,
                    help='alpha counted as "fully penetrated" (default 0.9). The '
                         'transition is sharp for large lambda but gradual below '
                         '~6.5, where gamma_c depends on this choice - worth a '
                         'sensitivity check there')
    ap.add_argument('--pen-frac', type=float, default=0.99,
                    help='t_pen is reported at this fraction of the final alpha '
                         '(default 0.99)')
    ap.add_argument('--tol', type=float, default=0.005,
                    help='relative width of the final gamma bracket (default 0.5%%)')
    ap.add_argument('--max-evals', type=int, default=20)
    ap.add_argument('--plateau', type=int, default=40,
                    help='periods of a flat sub-threshold alpha before giving up '
                         'on a gamma (default 40)')

    ap.add_argument('--ppd', type=float, default=8.0,
                    help='radial points per skin depth (default 8)')
    ap.add_argument('--nr-min', type=int, default=24)
    ap.add_argument('--nr-max', type=int, default=1024)
    ap.add_argument('--hall-c', type=float, default=5.0,
                    help='dt <= hall_c*dr^2 (default 5; measured limit is ~8)')
    ap.add_argument('--diff-c', type=float, default=0.4,
                    help='dt <= diff_c*lambda^2*dr^2 (default 0.4)')
    ap.add_argument('--tmax-factor', type=float, default=1.0,
                    help='time budget in units of lambda^2 RMF periods '
                         '(default 1.0, a ~0.6%% high bias in gamma_c)')
    ap.add_argument('--min-periods', type=int, default=60)
    args = ap.parse_args()

    lams = args.lam if args.lam else LAMBDAS
    if args.lam_max is not None:
        lams = [x for x in lams if x <= args.lam_max]

    out = paths.results(args.out)
    done = set()
    if args.resume and os.path.exists(out):
        with open(out) as f:
            done = {float(row['lam']) for row in csv.DictReader(f)}
        lams = [x for x in lams if x not in done]
        print(f'resuming: {len(done)} lambdas already in {out}')

    bias = (1.0/(4*math.pi*args.tmax_factor))**2
    print(f'\ngamma_c vs lambda - corrected scheme, rise_time = 0, cold starts')
    print(f'budget {args.tmax_factor:g}*lambda^2 periods '
          f'(floor {args.min_periods}) -> gamma_c biased ~{100*bias:.1f}% high')
    print(f'penetration when alpha >= {args.alpha_pen}, '
          f'bracket to {100*args.tol:g}%\n')
    print(f'{"lambda":>8} {"Nr":>6} {"dt":>10} {"periods":>9} '
          f'{"steps/eval":>12} {"Milroy":>9}')
    total = 0
    for lam in lams:
        Nr, dt, periods, steps = plan(lam, args)
        total += steps
        flag = '  <- Nr capped' if Nr >= args.nr_max else ''
        print(f'{lam:8.1f} {Nr:6d} {dt:10.2e} {periods:9d} {steps:12,d} '
              f'{gamma_c_milroy(lam):9.3f}{flag}')
    print(f'\n~{total:,d} steps per bisection pass; a full scan is roughly '
          f'8-12x that (early exits make penetrating evaluations cheaper).')
    if args.dry_run:
        print('\n--dry-run: nothing executed.')
        return

    new_file = not os.path.exists(out)
    with open(out, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if new_file:
            w.writeheader()
            f.flush()
        print(f'\n{"lambda":>8} {"gamma_c":>10} {"Milroy":>9} {"ratio":>7} '
              f'{"alpha_s":>8} {"t_pen":>8} {"evals":>6} {"time":>9}')
        for lam in lams:
            Nr, dt, periods, _ = plan(lam, args)
            t0 = time.time()
            lo, hi, n = find_gamma_c(lam, Nr, dt, periods, args)
            if math.isnan(lo):
                gc = a_s = t_pen = float('nan')
            else:
                gc = 0.5*(lo + hi)
                a_s, t_pen = refine_t_pen(lam, hi, Nr, dt, periods, args.pen_frac)
            secs = time.time() - t0
            eq = gamma_c_milroy(lam)
            w.writerow(dict(lam=lam, gamma_c=f'{gc:.6g}', gamma_c_lo=f'{lo:.6g}',
                            gamma_c_hi=f'{hi:.6g}', milroy_eq=f'{eq:.6g}',
                            ratio=f'{gc/eq:.6g}', alpha_s=f'{a_s:.6g}',
                            t_pen=f'{t_pen:.6g}', Nr=Nr, dt=f'{dt:.6g}',
                            max_periods=periods, n_evals=n, seconds=f'{secs:.1f}'))
            f.flush()
            print(f'{lam:8.1f} {gc:10.3f} {eq:9.3f} {gc/eq:7.3f} {a_s:8.3f} '
                  f'{t_pen:8.1f} {n:6d} {secs:8.1f}s', flush=True)
    print(f'\nwrote {out}   (plot with: python plot_gamma_c.py)')


if __name__ == '__main__':
    main()
