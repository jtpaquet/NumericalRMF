"""Ablation + convergence study behind the numerical choices in rmf_solver.py.

Each block isolates one difference between the original nr32_fix.py scheme and
the corrected one, so the effect on the comparison with Milroy 1999 Sec. III.A
can be attributed rather than guessed at.
"""
import numpy as np
from scipy.special import j1, iv
from rmf_solver import RMFPenetration, gamma_c, TWO_PI
from run_milroy import penetration_time

LEGACY = dict(A_op='naive', bc_order=1, b_scheme='legacy', rise_time=3*TWO_PI)
FIXED = dict(A_op='flux', bc_order=2, b_scheme='fv', rise_time=0.0)


def linear_limit(lam=11.07, gam=1e-3, n_periods=40, dt=0.002):
    """Check against the analytic small-gamma solution.

    Hugrass 1985 Eq. (20)-(21) (also Jones & Hugrass 1981): as gamma -> 0 the
    Hall term drops out and the steady n=1 flux obeys
        A'' + A'/r - A/r^2 + 2 i lam^2 A = 0,     A(1) + A'(1) = 2 gamma
    whose solution is  A = 2 gamma I1(k r) / (k I0(k)),  k = (1 - i) lam.
    (Hugrass writes -2i lam^2 because his phasor is e^{+i(wt-theta)}; ours is
    the complex conjugate convention, hence k = (1-i) lam rather than (1+i) lam.)

    This exercises the radial operator and the r = 1 boundary condition together,
    with no free parameters.
    """
    k = (1 - 1j)*lam
    print(f'\n== linear limit vs Hugrass 1985 Eq. (21), lambda = {lam} ==')
    print(f'{"Nr":>6} {"legacy L2":>12} {"|A(1)|/exact":>13} '
          f'{"corrected L2":>13} {"|A(1)|/exact":>13}')
    for Nr in (16, 32, 64, 128):
        row = []
        for legacy in (True, False):
            s = RMFPenetration(Nr=Nr, lam=lam, gam=gam, legacy=legacy, rise_time=0.0)
            s.run(n_periods=n_periods, dt=dt)
            A = s.A*np.exp(1j*s.tau)                # strip the e^{-i t} of the drive
            ex = 2*gam*iv(1, k*s.r)/(k*iv(0, k))
            row += [np.sqrt(np.mean(np.abs(A - ex)**2))/np.abs(ex).max(),
                    abs(A[-1])/abs(ex[-1])]
        print(f'{Nr:6d} {row[0]:12.3e} {row[1]:13.5f} {row[2]:13.3e} {row[3]:13.5f}')


def ramp(gam=16.6, Nr=64, n_periods=70, dt=0.002):
    """How much the RF turn-on ramp costs.

    Hugrass & Grimm 1981 Eq. (3) applies B_w(t) = (1 - e^{-t/tau_r}) B_w with
    tau_r = 0.4 us at omega = 5e6 /s, i.e. omega*tau_r = 2.0 = 0.32 RMF periods.
    Milroy 1999 applies the RMF as a step at t = 0 instead.
    """
    print(f'\n== gamma = {gam}: cost of the turn-on ramp ==')
    print(f'{"rise_time":24} {"t(0.9 alpha_s)":>15} {"alpha_s":>9}')
    for name, rt in (('0 (Milroy, step)', 0.0),
                     ('0.32 T (H&G 1981)', 2.0),
                     ('1 T', TWO_PI),
                     ('3 T (original code)', 3*TWO_PI)):
        s = RMFPenetration(Nr=Nr, gam=gam, rise_time=rt)
        t, a, _ = s.run(n_periods=n_periods, dt=dt)
        print(f'{name:24} {penetration_time(t, a, 0.9):15.1f} {a[-1]:9.4f}')


def operator_order():
    """L[A] on A = J1(kr): the 3-point form is only 1st order, the flux form 2nd."""
    k = 5.0
    print('\n== truncation error of L[A] for A = J1(5r) ==')
    print(f'{"Nr":>6} {"3-point":>12} {"flux":>12}')
    for Nr in (32, 64, 128, 256):
        r = np.linspace(0, 1, Nr)
        dr = r[1] - r[0]
        A, ex = j1(k*r), -k**2*j1(k*r)
        Ln = ((A[2:] - 2*A[1:-1] + A[:-2])/dr**2
              + (A[2:] - A[:-2])/(2*dr)/r[1:-1] - A[1:-1]/r[1:-1]**2)
        rh, rA = r[:-1] + dr/2, r*A
        w = (rA[1:] - rA[:-1])/(rh*dr)
        Lf = (w[1:] - w[:-1])/dr
        print(f'{Nr:6d} {np.abs(Ln-ex[1:-1]).max():12.3e} '
              f'{np.abs(Lf-ex[1:-1]).max():12.3e}')


def ablation(gam=16.6, Nr=64, n_periods=70, dt=0.002):
    """One change at a time, starting from the legacy scheme."""
    print(f'\n== gamma = {gam}, Nr = {Nr}: effect of each numerical choice ==')
    cases = [('legacy (nr32_fix.py)', {}),
             ('  + 2nd-order Robin BC', dict(bc_order=2)),
             ('  + finite-volume b eq.', dict(b_scheme='fv')),
             ('  + flux-form L[A]', dict(A_op='flux')),
             ('  + no gamma ramp', dict(rise_time=0.0)),
             ('all of the above', FIXED)]
    print(f'{"scheme":26} {"alpha_final":>11} {"t(0.95 alpha_s)":>16}')
    for name, over in cases:
        kw = dict(LEGACY)
        kw.update(over)
        s = RMFPenetration(Nr=Nr, gam=gam, **kw)
        t, a, _ = s.run(n_periods=n_periods, dt=dt)
        print(f'{name:26} {a[-1]:11.4f} {penetration_time(t, a, 0.95):16.1f}')


def grid_convergence(gam=14.9, n_periods=250, dt=0.002):
    """Subcritical steady state: the corrected scheme is converged by Nr = 32."""
    print(f'\n== gamma = {gam} (subcritical): alpha_s vs Nr, {n_periods} T ==')
    print(f'{"Nr":>6} {"legacy":>10} {"corrected":>10}')
    for Nr in (32, 64, 128):
        row = []
        for kw in (LEGACY, FIXED):
            s = RMFPenetration(Nr=Nr, gam=gam, **kw)
            _, a, _ = s.run(n_periods=n_periods, dt=dt)
            row.append(a[-1])
        print(f'{Nr:6d} {row[0]:10.4f} {row[1]:10.4f}')
    print(f'  Milroy 1999 quotes alpha_s = 0.42 for gamma = 14.9 '
          f'(gamma/gamma_c = {gam/gamma_c(11.07):.3f})')


if __name__ == '__main__':
    linear_limit()
    operator_order()
    ramp()
    ablation()
    grid_convergence()
