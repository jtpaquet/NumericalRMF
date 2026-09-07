#!/usr/bin/env python3
"""Study 3 extension - does an ionisation fraction < 1 let the *nonlinear*
Milroy penetration threshold explain the observed 1-50 ms, pressure-dependent
coupling time?

Motivation
----------
Part 1 of this study (1_scan_grid.py) assumed full single ionisation
(n_e = n_gas) and found gamma << gamma_c everywhere: purely classical,
sub-ms field diffusion, no trace of the observed timescale. But the actual
ionisation fraction right after a cosmic-ray-seeded avalanche is whatever it
is -- it is not an input, it is the thing breakdown produces. Since
gamma = B_w/(e n eta) grows as n shrinks, a lower ionisation fraction pushes
gamma up, and if it crosses gamma_c(lambda) the *nonlinear* penetration time,
Milroy Eq. (17)

    tau_P(lambda, gamma) = lambda^2 / (2 sqrt(gamma_N)),   gamma_N = (gamma-gamma_c)/gamma_c

diverges at threshold -- in principle capable of producing a slow, strongly
gamma-dependent (hence density- and pressure-dependent) coupling time of the
kind observed, unlike the fast, threshold-free classical diffusion of part 1.

This script:
  1. validates Eq. (17) against the real time-dependent solver in this
     study's own lambda range (previously only checked at Milroy's
     lambda = 11.07, see validate_numerics.py and run_milroy.py);
  2. evaluates it at the seven requested ionisation fractions
     f_ion in {1e-6 ... 1} across the (p, T_e) grid;
  3. solves for the ionisation fraction that would put each (p, T_e) point
     exactly at gamma_c, and the narrow band around it that reproduces
     1-50 ms via the Eq. (17) divergence -- to see whether that band is
     anywhere near a value a stochastic avalanche could plausibly land on.

Outputs
-------
    results/study3_ionization_grid.csv     tau_penetration at the 7 f_ion values
    results/study3_ionization_window.csv   critical f_ion and the 1/50 ms band
"""
import csv
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import paths
import device as dev
from rmf_solver import RMFPenetration, gamma_c, tau_penetration
from run_milroy import penetration_time

VALIDATION_LAMBDAS = (9.57, 15.57)     # this study's Te=2.24, 4.73 cases
VALIDATION_MULTS = (1.3, 2.0)          # gamma/gamma_c to test the formula at


# --------------------------------------------------------- root-finding (1D)
def solve_n_for_gamma(gamma_target, Te, n_guess, iters=6):
    """n such that gam_of(n, Te) = gamma_target, iterating eta's weak (log)
    dependence on n through lnLambda to self-consistency."""
    n = n_guess
    for _ in range(iters):
        eta = dev.eta_spitzer(n, Te)
        n = dev.B_W/(dev.E_CHARGE*gamma_target*eta)
    return n


def n_at_gamma_over_gc(Te, ratio, n_guess, iters=6):
    """n such that gamma(n, Te)/gamma_c(lambda(Te, n)) = ratio."""
    n = n_guess
    for _ in range(iters):
        lam = dev.lam_of(Te, n)
        gc = gamma_c(lam)
        n = solve_n_for_gamma(ratio*gc, Te, n, iters=3)
    return n, dev.lam_of(Te, n), gamma_c(dev.lam_of(Te, n))


def n_for_tau(Te, tau_ms, n_guess, iters=6):
    """n such that tau_penetration(lambda(Te,n), gamma(n,Te)) = tau_ms."""
    n = n_guess
    tau_periods = tau_ms/dev.period_ms()
    for _ in range(iters):
        lam = dev.lam_of(Te, n)
        gc = gamma_c(lam)
        gN = (lam**2/(2.0*tau_periods))**2
        n = solve_n_for_gamma(gc*(1.0 + gN), Te, n, iters=3)
    return n, dev.lam_of(Te, n), gamma_c(dev.lam_of(Te, n))


# --------------------------------------------------------------- validation
def validate_eq17():
    print('validating Milroy Eq. (17) against the real solver at this '
          "study's lambdas (Milroy only checked lambda = 11.07):")
    print(f'{"lambda":>7} {"gamma":>7} {"gam/gc":>7} {"Nr":>5} {"periods":>8} '
          f'{"alpha_s":>8} {"t_pen(num)":>11} {"Eq17":>7} {"ratio":>6} {"sec":>6}')
    rows = []
    for lam in VALIDATION_LAMBDAS:
        gc = gamma_c(lam)
        for mult in VALIDATION_MULTS:
            gam = mult*gc
            Nr = dev.choose_Nr(lam, ppd=8.0)
            dr = 1.0/(Nr - 1)
            dt = min(0.002, 5.0*dr**2, 0.4*lam**2*dr**2)
            eq17 = tau_penetration(lam, gam)
            periods = int(max(80, np.ceil(3*eq17)))
            sim = RMFPenetration(Nr=Nr, lam=lam, gam=gam)
            t0 = time.time()
            t, a, _ = sim.run(n_periods=periods, dt=dt)
            secs = time.time() - t0
            tp_num = penetration_time(t, a, 0.95)
            ratio = tp_num/eq17 if np.isfinite(tp_num) else np.nan
            print(f'{lam:7.2f} {gam:7.2f} {mult:7.2f} {Nr:5d} {periods:8d} '
                  f'{a[-1]:8.3f} {tp_num:11.2f} {eq17:7.2f} {ratio:6.2f} {secs:6.1f}s')
            rows.append(dict(lam=lam, gam=gam, gam_over_gc=mult, Nr=Nr,
                              periods=periods, alpha_s=a[-1], t_pen_num=tp_num,
                              t_pen_eq17=eq17, ratio=ratio))
    print()
    return rows


def main():
    val_rows = validate_eq17()

    # ------------------------------------------------- the 7 requested f_ion
    print('tau_penetration (Eq. 17) at the requested ionisation fractions:')
    grid_rows = []
    for p in dev.P_MTORR_GRID:
        n_gas = dev.density_from_pressure(p)
        for Te in dev.TE_EV_GRID:
            for f in dev.F_ION:
                n = f*n_gas
                lam = dev.lam_of(Te, n)
                gam = dev.gam_of(n, Te)
                gc = gamma_c(lam)
                tau_p = tau_penetration(lam, gam)
                tau_ms = tau_p*dev.period_ms()
                grid_rows.append(dict(p_mtorr=p, Te_eV=Te, f_ion=f, n_m3=n,
                                      lam=lam, gam=gam, gamma_c=gc,
                                      gam_over_gc=gam/gc, t_pen_periods=tau_p,
                                      t_pen_ms=tau_ms))
    out = paths.results('study3_ionization_grid.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(grid_rows[0].keys()))
        w.writeheader()
        w.writerows(grid_rows)
    n_in_band = sum(1 for r in grid_rows if 1.0 <= r['t_pen_ms'] <= 50.0)
    print(f'wrote {out}  ({len(grid_rows)} rows, {n_in_band} land in 1-50 ms)\n')

    # ---------------------------------------- critical f_ion and the window
    print('critical ionisation fraction (gamma = gamma_c) and the band that '
          'gives 1-50 ms via the Eq. (17) divergence:')
    print(f'{"p(mTorr)":>9} {"Te(eV)":>7} {"f_crit":>9} {"f(50ms)":>9} '
          f'{"f(1ms)":>9} {"band/f_crit":>11}')
    window_rows = []
    for p in dev.P_MTORR_GRID:
        n_gas = dev.density_from_pressure(p)
        for Te in dev.TE_EV_GRID:
            n_crit, lam, gc = n_at_gamma_over_gc(Te, 1.0, n_gas)
            n_50, _, _ = n_for_tau(Te, 50.0, n_crit)
            n_1, _, _ = n_for_tau(Te, 1.0, n_crit)
            f_crit, f_50, f_1 = n_crit/n_gas, n_50/n_gas, n_1/n_gas
            band = abs(f_1 - f_50)/f_crit
            window_rows.append(dict(p_mtorr=p, Te_eV=Te, lam=lam, gamma_c=gc,
                                    f_crit=f_crit, f_50ms=f_50, f_1ms=f_1,
                                    band_over_fcrit=band))
            print(f'{p:9.1f} {Te:7.2f} {f_crit:9.2e} {f_50:9.2e} {f_1:9.2e} '
                  f'{band:11.2e}')
    out = paths.results('study3_ionization_window.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(window_rows[0].keys()))
        w.writeheader()
        w.writerows(window_rows)
    print(f'\nwrote {out}   (plot with: python 5_plot_ionization.py)')


if __name__ == '__main__':
    main()
