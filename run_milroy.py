"""Reproduce Sec. III.A of Milroy 1999 (Figs. 1 and 2) with the corrected solver.

Milroy's reference numbers for R = 10 cm, n = 0.333e20 m^-3, w = 2.2e5 s^-1,
T_e = 5 eV  ->  lam = 11.07:

    gamma = 14.9 (B_w =  90 G) : never fully penetrates, alpha_s = 0.42
    gamma = 16.6 (B_w = 100 G) : alpha_s = 0.98, full penetration at t = 40 T
    gamma = 18.2 (B_w = 110 G) : full penetration at t = 27 T
"""
import pickle
import numpy as np
import paths
from rmf_solver import RMFPenetration, gamma_c, tau_penetration, TWO_PI

LAM = 11.07
NR = 64
DT = 0.002
N_PERIODS = 200                       # gamma = 14.9 needs ~200 T to reach alpha_s
GAMMAS = [14.9, 16.6, 18.2]
SNAPSHOTS = (5, 10, 20, 30, 40, 45, 50, 60)

MILROY = {14.9: dict(alpha_s=0.42, t_pen=None),
          16.6: dict(alpha_s=0.98, t_pen=40),
          18.2: dict(alpha_s=0.98, t_pen=27)}


def penetration_time(t, a, frac=0.95):
    """First time alpha reaches frac of alpha_s, with alpha_s = the last value."""
    thr = frac*a[-1]
    if a.max() < thr:
        return np.nan
    i = int(np.argmax(a >= thr))
    return t[0] if i == 0 else t[i-1] + (thr - a[i-1])*(t[i]-t[i-1])/(a[i]-a[i-1])


def main(out='all_results_corrected.pkl'):
    out = paths.results(out)
    gc = gamma_c(LAM)
    print(f'lambda = {LAM}   gamma_c = {gc:.3f} (Milroy Eq. 15); '
          f'expulsion threshold 1.12*lambda = {1.12*LAM:.3f} (Eq. 14)')
    results = {}
    for gam in GAMMAS:
        sim = RMFPenetration(Nr=NR, lam=LAM, gam=gam)
        t, a, snaps = sim.run(n_periods=N_PERIODS, dt=DT, snapshot_periods=SNAPSHOTS)
        if not np.all(np.isfinite(a)):
            raise RuntimeError(f'gamma={gam} went unstable; reduce dt below {DT}')
        results[gam] = dict(times=t, alphas=a, snapshots=snaps,
                            lam=LAM, gam=gam, Nr=NR, dt=DT)
        ref = MILROY[gam]
        # only meaningful once the RMF actually penetrates
        tp = penetration_time(t, a) if a[-1] > 0.9 else np.nan  # 95% of alpha_s
        eq17 = tau_penetration(LAM, gam)
        print(f'gamma={gam:5.1f}  gamma/gamma_c={gam/gc:5.3f}   '
              f'alpha_s={a[-1]:.3f} (Milroy {ref["alpha_s"]:.2f})   '
              f't_pen={tp:5.1f} T (Milroy {ref["t_pen"] or float("nan"):5.1f}, '
              f'Eq.17 {eq17:5.1f})')
    with open(out, 'wb') as f:
        pickle.dump(results, f)
    print(f'wrote {out}')
    return results


if __name__ == '__main__':
    main()
