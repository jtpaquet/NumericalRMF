"""Milroy's own threshold procedure (1999, Sec. III.A, Figs. 4 and 6).

Hold gamma fixed until alpha settles, step it up, repeat.  The jump in alpha_s
locates gamma_c (Eq. 15).  Then step gamma back down from a penetrated state:
the RMF is expelled at the lower threshold 1.12*lambda (Eq. 14), so the two
branches bracket the hysteresis loop of his Fig. 9.

This is the sharpest quantitative check available against the paper that does
not need any extra physics: gamma_c is a single number, and Milroy gives a
closed-form fit for it.
"""
import numpy as np
from rmf_solver import RMFPenetration, gamma_c, TWO_PI

LAM = 11.07
NR = 64
DT = 0.0015           # the explicit scheme tightens as gamma grows; 0.002 fails ~gamma=20
HOLD = 80             # RMF periods held at each gamma


def alpha_s_fit(lam, gam):
    """Milroy Eq. (18): alpha_s = 1.6/sqrt(lam) * exp(-4*((gc-gam)/gc)^p).

    The exponent is printed as 8 in the paper, which cannot be right (it makes
    the exponential indistinguishable from 1 over the whole sub-critical range,
    so alpha_s would be flat at 1.6/sqrt(lam) = 0.48 all the way down).  p ~ 0.8
    reproduces the alpha_s = 0.42 he quotes at gamma = 14.9; both are plotted so
    the discrepancy is visible rather than hidden.
    """
    gc = gamma_c(lam)
    x = max(0.0, (gc - gam)/gc)
    return 1.6/np.sqrt(lam)*np.exp(-4.0*x**0.8)


def sweep(sim, gammas, label):
    rows = []
    for gam in gammas:
        sim.gam = float(gam)
        for _ in range(HOLD):
            for _ in range(int(round(TWO_PI/DT))):
                sim.step(DT)
        a = sim.alpha()
        if not np.isfinite(a):
            print(f'{label}: went unstable at gamma={gam}; reduce DT')
            break
        rows.append((float(gam), a))
        print(f'{label}  gamma={gam:5.2f}  alpha_s={a:.4f}  (Eq.18 fit {alpha_s_fit(LAM, gam):.3f})',
              flush=True)
    return np.array(rows)


def main():
    gc = gamma_c(LAM)
    print(f'lambda={LAM}: Eq.(15) gamma_c={gc:.3f}, Eq.(14) expulsion at {1.12*LAM:.3f}\n')

    up = sweep(RMFPenetration(Nr=NR, lam=LAM, gam=0.0),
               np.arange(10.0, 16.51, 0.25), 'up  ')
    np.savetxt('alpha_s_up.txt', up, header='gamma alpha_s')

    down_sim = RMFPenetration(Nr=NR, lam=LAM, gam=17.0)
    for _ in range(120):                       # get fully penetrated first
        for _ in range(int(round(TWO_PI/DT))):
            down_sim.step(DT)
    print(f'\npenetrated at gamma=17.0: alpha_s={down_sim.alpha():.4f}\n')
    down = sweep(down_sim, np.arange(16.75, 9.99, -0.25), 'down')
    np.savetxt('alpha_s_down.txt', down, header='gamma alpha_s')

    def threshold(rows, rising):
        """gamma at which alpha_s crosses 0.75, bracketed by the sweep step."""
        a = rows[:, 1]
        idx = np.flatnonzero((a[:-1] < 0.75) & (a[1:] >= 0.75) if rising
                             else (a[:-1] >= 0.75) & (a[1:] < 0.75))
        return (rows[idx[0], 0], rows[idx[0]+1, 0]) if len(idx) else (np.nan, np.nan)

    lo, hi = threshold(up, True)
    print(f'\nmeasured penetration threshold in ({lo}, {hi}]   Milroy Eq.(15): {gc:.3f}')
    lo, hi = threshold(down, False)
    print(f'measured expulsion   threshold in [{hi}, {lo})   Milroy Eq.(14): {1.12*LAM:.3f}')


if __name__ == '__main__':
    main()
