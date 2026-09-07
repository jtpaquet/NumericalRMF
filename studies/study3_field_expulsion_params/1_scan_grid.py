#!/usr/bin/env python3
"""Study 3 - can the fixed-ion RMF model reproduce the observed field
expulsion (t_pen ~ 10 ms, full expulsion ~ 25 ms after breakdown) for
plausible pre-ionization argon fill conditions?

Pre-ionization state
--------------------
Argon at room temperature (25 C), fill pressure 10-100 mTorr -> a plasma
density n_e via the ideal-gas density of the neutral fill (i.e. assuming the
breakdown fully singly-ionises the fill; see device.py). Plasma electron
temperature T_e = 0.5-10 eV is scanned independently, since it is set by the
discharge physics, not the fill pressure.

Method
------
n and T_e enter the solver only through the Spitzer resistivity eta(n, T_e)
(device.py), which sets both dimensionless numbers:

    lambda(T_e)   = R sqrt(mu0 w / (2 eta))       (n enters only via lnLambda)
    gamma(n, T_e) = B_w / (e n eta)

At this device's drive amplitude (B_w = 6 G) gamma comes out << gamma_c
everywhere in the (n, T_e) box (see the printed table) - i.e. this is deep in
the *linear* regime, where the Hall term is negligible and the transverse
field profile is the analytic small-gamma solution (Hugrass 1985 Eq. 21,
checked in validate_numerics.linear_limit): B_r(r)/B_w, B_theta(r)/B_w do not
depend on gamma to leading order, only on lambda. So instead of a brute-force
run per (n, T_e) pair, this script runs the real time-dependent solver *once
per T_e* at a small reference gamma to get the lambda-dependent screening
transient, then checks that assumption directly with one extra run at the
largest gamma actually reached in the grid.

Outputs
-------
    results/study3_reference_by_Te.pkl   time-dependent Br,Btheta at r=0,1
                                          for each T_e (the reference runs)
    results/study3_grid.csv              the full (p, T_e) grid: n, lambda,
                                          gamma, gamma/gamma_c, and the
                                          reference-run predictions for the
                                          settling time and steady screening
"""
import csv
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import paths
import device as dev
from rmf_solver import RMFPenetration, gamma_c, TWO_PI

P_MTORR = dev.P_MTORR_GRID
TE_EV = dev.TE_EV_GRID
GAM_REF = 0.01                        # nominal small gamma for the reference runs


def run_profiles(sim, n_periods, dt):
    """Like RMFPenetration.run(), but also records B_r, B_theta at r=0 and r=1."""
    nstep = int(round(TWO_PI/dt))
    t, Br0, BrR, Bth0, BthR, alphas = [np.empty(n_periods) for _ in range(6)]
    for p in range(n_periods):
        for _ in range(nstep):
            sim.step(dt)
        r, Br, Bth, Bz = sim.profiles()
        t[p] = p + 1.0
        Br0[p], BrR[p] = Br[0], Br[-1]
        Bth0[p], BthR[p] = Bth[0], Bth[-1]
        alphas[p] = sim.alpha()
    return t, Br0, BrR, Bth0, BthR, alphas


def settling_time_ms(t_ms, gap, frac=0.9):
    """First time |B_r(R)-B_r(0)| (or theta) reaches frac of its final value."""
    thr = frac*gap[-1]
    if gap.max() < thr:
        return np.nan
    i = int(np.argmax(gap >= thr))
    if i == 0:
        return t_ms[0]
    return t_ms[i-1] + (thr - gap[i-1])*(t_ms[i] - t_ms[i-1])/(gap[i] - gap[i-1])


def main():
    print(f'device: f_RMF = {dev.F_RMF/1e3:.0f} kHz (period = {dev.period_ms():.4f} ms), '
          f'R = {dev.R*100:.0f} cm, B_w = {dev.B_W_G:.1f} G '
          f'(calibration cross-check: {dev.B_W_FROM_CALIBRATION_G:.1f} G)')
    print(f'measured with plasma coupled: B(wall) ~ {dev.MEASURED_B_WALL_G:.1f} G, '
          f'B(centre) ~ {dev.MEASURED_B_CENTER_G:.1f} G\n')

    # ---------------------------------------------------- reference runs, per Te
    n_mid = dev.density_from_pressure(float(np.median(P_MTORR)))
    reference = {}
    print(f'{"Te(eV)":>7} {"lambda":>7} {"Nr":>5} {"dt":>9} {"periods":>8} '
          f'{"t_settle(ms)":>13} {"seconds":>8}')
    for Te in TE_EV:
        lam = dev.lam_of(Te, n_mid)
        Nr = dev.choose_Nr(lam)
        dt = dev.choose_dt(lam, Nr)
        tau = dev.skin_time_ms(Te, n_mid)
        periods = dev.choose_periods(tau)
        sim = RMFPenetration(Nr=Nr, lam=lam, gam=GAM_REF)
        t0 = time.time()
        t, Br0, BrR, Bth0, BthR, alphas = run_profiles(sim, periods, dt)
        secs = time.time() - t0
        t_ms = t*dev.period_ms()
        t_settle = settling_time_ms(t_ms, np.abs(BrR - Br0))
        reference[Te] = dict(Te=Te, lam=lam, gam=GAM_REF, Nr=Nr, dt=dt,
                              periods=periods, t_ms=t_ms, Br0=Br0, BrR=BrR,
                              Bth0=Bth0, BthR=BthR, alphas=alphas,
                              t_settle_ms=t_settle)
        print(f'{Te:7.2f} {lam:7.2f} {Nr:5d} {dt:9.2e} {periods:8d} '
              f'{t_settle:13.4f} {secs:8.1f}s')
    out = paths.results('study3_reference_by_Te.pkl')
    with open(out, 'wb') as f:
        pickle.dump(reference, f)
    print(f'wrote {out}\n')

    # ---------------------------- validate the gamma-independence assumption
    Te_hi, p_lo = TE_EV.max(), P_MTORR.min()
    n_lo = dev.density_from_pressure(p_lo)
    gam_max = dev.gam_of(n_lo, Te_hi)
    ref = reference[Te_hi]
    sim = RMFPenetration(Nr=ref['Nr'], lam=ref['lam'], gam=gam_max)
    _, Br0, BrR, Bth0, BthR, _ = run_profiles(sim, ref['periods'], ref['dt'])
    d_ref = abs(ref['BrR'][-1] - ref['Br0'][-1])
    d_hi = abs(BrR[-1] - Br0[-1])
    print(f'linear-regime check at the largest gamma in the grid '
          f'(p={p_lo:.0f} mTorr, Te={Te_hi:.1f} eV, gamma={gam_max:.3f} '
          f'vs reference gamma={GAM_REF}):')
    print(f'  |B_r(R)-B_r(0)|/B_w: reference {d_ref:.4f}, actual gamma {d_hi:.4f} '
          f'({100*abs(d_hi-d_ref)/d_ref:.2f}% difference)\n')

    # --------------------------------------------------------------- the grid
    rows = []
    print(f'{"p(mTorr)":>9} {"Te(eV)":>7} {"n(m^-3)":>10} {"lambda":>7} '
          f'{"gamma":>9} {"gamma_c":>8} {"gam/gamc":>9} {"t_settle(ms)":>13} '
          f'{"Br0/Bw":>7} {"BrR/Bw":>7}')
    for p in P_MTORR:
        n = dev.density_from_pressure(p)
        for Te in TE_EV:
            lam = dev.lam_of(Te, n)
            gam = dev.gam_of(n, Te)
            gc = gamma_c(lam)
            ref = reference[Te]
            row = dict(p_mtorr=p, Te_eV=Te, n_m3=n, lam=lam, gam=gam,
                       gamma_c=gc, gam_over_gc=gam/gc,
                       t_settle_ms=ref['t_settle_ms'],
                       Br0_over_Bw=ref['Br0'][-1], BrR_over_Bw=ref['BrR'][-1],
                       Bth0_over_Bw=ref['Bth0'][-1], BthR_over_Bw=ref['BthR'][-1])
            rows.append(row)
            print(f'{p:9.1f} {Te:7.2f} {n:10.3e} {lam:7.2f} {gam:9.3e} '
                  f'{gc:8.3f} {gam/gc:9.2e} {ref["t_settle_ms"]:13.4f} '
                  f'{ref["Br0"][-1]:7.3f} {ref["BrR"][-1]:7.3f}')
    out = paths.results('study3_grid.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'\nwrote {out}   (plot with: python 2_plot_grid.py)')


if __name__ == '__main__':
    main()
