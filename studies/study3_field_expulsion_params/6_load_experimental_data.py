#!/usr/bin/env python3
"""Study 3 - real shots: does the classical linear-regime picture (parts 1-2
of this study) match the actual measured diamagnetic response?

Data
----
`data/env_diff/*.csv`: bdot probe + Rogowski-inferred envelope pairs at the
chamber centre, 250 kHz RMF, 43 A DC field coil circuit, various added series
resistance (weakens the DC bias field) and argon fill pressure. Provenance
(from the user, verbatim):

    w = np.ones(100)/100                                  # 2.5 MHz raw -> envelope
    bdot_x_env  = convolve(abs(bdot_x_G),  w, 'same')*pi/2
    B_rog_x_env = I_rog_2_env * 0.163 G/A                  # coil geometry factor
    Env_diff_x  = bdot_x_env - B_rog_x_env
    (further boxcar-averaged over 5000 samples and downsampled to 0.8 ms)

`bdot_x/z_env` is the actually-measured field at the probe (chamber centre).
`B_rog_x/z_env` is what that field *would* be from the actual driven antenna
current alone, i.e. the vacuum-coupled prediction at the same point -- so it
plays the role of a shot-by-shot, dynamically measured B_w, not the fixed
6 G "no plasma" value used in parts 1-2 (which only holds before breakdown;
once the plasma loads the antenna, the driven current itself drops, and
B_rog tracks that). Env_diff = bdot - B_rog is therefore this shot's own
B_r(0) - B_w(t), and bdot/B_rog is exactly the model's `br0_over_bw(lambda)`
(device.py), evaluated at the plasma's actual (unknown) T_e.

RMF is on from t = 45 ms to (nominally) 145 ms; a hardware issue means it
doesn't fully cut off until t = 180 ms, which shows up as a glitch/recovery
near 145 ms in every shot (not real plasma physics) -- 4_ionization etc. are
unaffected since this is a separate part of the study; the summary window
used below (110-135 ms) stays clear of it.

Output
------
    results/study3_experimental_summary.csv
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.optimize import brentq

import paths
import device as dev

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'env_diff')
T_ON_SEARCH = (10.0, 200.0)          # search window for the RF turn-on step
STEADY_WINDOW = (110.0, 135.0)       # quasi-steady window, clear of the ~145 ms glitch
ONSET_THRESHOLD = 0.05               # |Env_diff| above this = "RF is on"


def parse_shot(fn):
    base = os.path.basename(fn)[:-4]
    parts = [p for p in base.split('_') if p]
    shot = parts[2]
    p_mtorr = ohm = None
    for tok in parts[5:]:
        if tok.endswith('mTorr'):
            p_mtorr = int(tok[:-5])
        elif tok.endswith('ohm'):
            ohm = int(tok[:-3])
    return dict(shot=shot, ohm=ohm, noDC=('noDC' in parts),
               reversed=('reversed' in parts), p_mtorr=p_mtorr, fn=fn)


def load_shot(fn):
    t, dx, dz, bx, bz, rx, rz = ([] for _ in range(7))
    with open(fn) as f:
        for row in csv.DictReader(f):
            t.append(float(row['t_ms']))
            dx.append(float(row['Env_diff_x'])); dz.append(float(row['Env_diff_z']))
            bx.append(float(row['bdot_x_env_LP'])); bz.append(float(row['bdot_z_env_LP']))
            rx.append(float(row['B_rog_x_env_LP'])); rz.append(float(row['B_rog_z_env_LP']))
    return {k: np.array(v) for k, v in
            dict(t=t, dx=dx, dz=dz, bx=bx, bz=bz, rx=rx, rz=rz).items()}


def Te_floor(n_guess):
    """Lowest T_e where lnLambda(n, T_e) > 0, i.e. where Spitzer/weak-coupling
    theory is still valid at this density -- the real physical floor, not an
    arbitrary bracket edge."""
    Te = 0.3
    while dev.coulomb_log_ei(n_guess, Te) > 0.5 and Te > 1e-3:
        Te *= 0.7
    return Te/0.7


def Te_from_ratio(ratio, n_guess, Te_hi=15.0):
    """Invert br0_over_bw(lambda(Te)) = ratio for T_e, at fixed n (weak lnLambda
    dependence). Returns (Te, floored): floored=True means ratio is weaker
    (closer to 1, less screening) than the model can produce anywhere the
    Spitzer resistivity is still valid -- i.e. no T_e >= Te_floor fits."""
    if not (0.0 < ratio < 1.0):
        return np.nan, False
    lo = Te_floor(n_guess)
    r_lo, r_hi = dev.br0_over_bw(dev.lam_of(lo, n_guess)), dev.br0_over_bw(dev.lam_of(Te_hi, n_guess))
    if ratio >= r_lo:
        return lo, True
    if ratio <= r_hi:
        return Te_hi, True
    f = lambda Te: dev.br0_over_bw(dev.lam_of(Te, n_guess)) - ratio
    return brentq(f, lo, Te_hi), False


def analyze(fn):
    info = parse_shot(fn)
    d = load_shot(fn)
    t = d['t']

    mag = np.sqrt(d['dx']**2 + d['dz']**2)
    on_idx = np.argmax((mag > ONSET_THRESHOLD) & (t > T_ON_SEARCH[0]) & (t < T_ON_SEARCH[1]))
    t_on = t[on_idx]

    steady = (t >= STEADY_WINDOW[0]) & (t <= STEADY_WINDOW[1])
    ratio_x = np.median(d['bx'][steady]/d['rx'][steady]) if steady.sum() else np.nan
    ratio_z = np.median(d['bz'][steady]/d['rz'][steady]) if steady.sum() else np.nan
    dx_final = np.median(d['dx'][steady]) if steady.sum() else np.nan
    dz_final = np.median(d['dz'][steady]) if steady.sum() else np.nan

    thr = 0.9*dx_final
    post = (t >= t_on) & (t <= 135.0)
    idx = np.where(d['dx'][post] <= thr)[0]           # dx_final is negative
    t_couple90 = (t[post][idx[0]] - t_on) if len(idx) else np.nan

    n_gas = dev.density_from_pressure(info['p_mtorr']) if info['p_mtorr'] else np.nan
    if info['p_mtorr']:
        Te_x, floor_x = Te_from_ratio(ratio_x, n_gas)
        Te_z, floor_z = Te_from_ratio(ratio_z, n_gas)
    else:
        Te_x = Te_z = np.nan
        floor_x = floor_z = False

    info.update(t_on=t_on, ratio_x=ratio_x, ratio_z=ratio_z, dx_final=dx_final,
               dz_final=dz_final, t_couple90_ms=t_couple90, Te_x_eV=Te_x, Te_z_eV=Te_z,
               Te_x_at_floor=floor_x, Te_z_at_floor=floor_z)
    return info


def main():
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.csv')))
    if not files:
        raise SystemExit(f'no shots found under {DATA_DIR}')
    rows = [analyze(fn) for fn in files]
    for r in rows:
        r.pop('fn')

    print(f'{"shot":>6} {"p(mTorr)":>9} {"ohm":>5} {"rev":>5} {"t_on":>6} '
          f'{"ratio_x":>8} {"ratio_z":>8} {"t_c90(ms)":>10} {"Te_x(eV)":>10} {"Te_z(eV)":>10}')
    for r in sorted(rows, key=lambda r: (r['reversed'], r['ohm'] or -1, r['p_mtorr'])):
        ohm_s = 'noDC' if r['noDC'] else r['ohm']
        tex = f"{r['Te_x_eV']:.3f}{'*' if r['Te_x_at_floor'] else ''}"
        tez = f"{r['Te_z_eV']:.3f}{'*' if r['Te_z_at_floor'] else ''}"
        print(f"{r['shot']:>6} {r['p_mtorr']:9d} {str(ohm_s):>5} {str(r['reversed']):>5} "
              f"{r['t_on']:6.1f} {r['ratio_x']:8.3f} {r['ratio_z']:8.3f} "
              f"{r['t_couple90_ms']:10.2f} {tex:>10} {tez:>10}")
    print("(* = ratio weaker than the model gives anywhere lnLambda > 0; Te is that floor, a lower bound)")

    out = paths.results('study3_experimental_summary.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'\nwrote {out}   (plot with: python 7_plot_experimental.py)')


if __name__ == '__main__':
    main()
