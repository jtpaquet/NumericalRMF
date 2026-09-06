"""Figures for the Milroy 1999 Sec. III.A comparison.

Corrections relative to the first version of this script:
  * gamma_c is Milroy's Eq. (15) penetration threshold (15.13 at lambda=11.07),
    not 1.12*lambda = 12.40, which is the *expulsion* threshold (Eq. 14).
    With the right value gamma = 14.9 sits just *below* threshold
    (gamma/gamma_c = 0.985), which is why it never penetrates.
  * the alpha_s target quoted by Milroy for gamma = 14.9 is 0.42, not 0.48
    (0.48 = 1.6/sqrt(lambda) is the gamma -> gamma_c limit of his Eq. 18 fit).
  * B_z/B_omega = lambda^2 * b / gamma.  The stored B array is b = B_z/lambda^2,
    so dividing it by gamma alone understates the axial field by lambda^2 = 123.
  * the field-line plots are drawn in the frame co-rotating with the RMF, which
    means Re[A e^{i(theta + tau)}]: A already carries the e^{-i tau} of the
    drive, so e^{i(theta - tau)} leaves the picture spinning backwards at 2*omega.
  * |B_theta(R)| comes from the boundary condition (A'(1) = 2 gamma e^{-i tau}
    - A(1)) instead of a one-sided difference.
"""
import pickle
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from rmf_solver import gamma_c
from utils import plot_options

PICKLE = sys.argv[1] if len(sys.argv) > 1 else 'all_results_corrected.pkl'
LAM = 11.07
GAMMA_C = gamma_c(LAM)
COLORS = {14.9: 'blue', 16.6: 'green', 18.2: 'red'}
PLOT_TIMES = [10, 30, 45, 60]

# Milroy 1999, Sec. III.A: (t of full penetration, alpha_s)
MILROY = {14.9: (None, 0.42), 16.6: (40, 0.98), 18.2: (27, 0.98)}
XMAX = 80          # x-limit of the alpha(t) figure


def load(fname):
    with open(fname, 'rb') as f:
        return pickle.load(f)


def gamma_of(res, gam):
    """Drive amplitude at the snapshot time (constant unless a ramp was used)."""
    return res.get('gam', gam)


def b_theta_wall(A, tau, gam):
    """|B_theta(R)|/B_omega straight from the Robin boundary condition."""
    return abs(2*gam*np.exp(-1j*tau) - A[-1])/gam


def dA_dr(A, r, tau, gam):
    dr = r[1] - r[0]
    d = np.empty_like(A)
    d[1:-1] = (A[2:] - A[:-2])/(2*dr)
    d[0] = (-3*A[0] + 4*A[1] - A[2])/(2*dr)
    d[-1] = 2*gam*np.exp(-1j*tau) - A[-1]          # exact
    return d


def fig_alpha(all_results, fname='alpha_vs_time.pdf'):
    fig, ax = plt.subplots(figsize=(14, 8))
    for gam in sorted(all_results):
        res = all_results[gam]
        ax.plot(res['times'], res['alphas'], color=COLORS.get(gam), lw=2,
                label=fr'$\gamma$={gam} ($\gamma/\gamma_c$={gam/GAMMA_C:.2f})')
        t_pen, a_s = MILROY[gam]
        # subcritical cases have no penetration time; park the marker at the
        # right-hand edge, where the curve has settled on alpha_s
        ax.scatter([t_pen if t_pen else XMAX], [a_s],
                   color=COLORS.get(gam), s=200, marker='*', zorder=5,
                   label='Milroy 1999' if gam == min(all_results) else None)
        print(f'gamma={gam}: alpha_s={res["alphas"][-1]:.4f}  (Milroy {a_s})')
    ax.axhline(1.0, color='k', ls='--', alpha=0.5)
    ax.set_xlabel('Time (RMF periods)')
    ax.set_ylabel(r'Penetration Factor $\alpha$')
    Nr = all_results[min(all_results)].get('Nr', '?')
    ax.set_title(rf'$\alpha$ vs Time - $\lambda$={LAM}, $N_r$={Nr}')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.2)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    print(f'Saved {fname}')


def fig_alpha_compare(new, old, fname='alpha_vs_time_comparison.pdf'):
    """alpha(t) before and after the numerical corrections, against Milroy."""
    fig, ax = plt.subplots(figsize=(14, 8))
    for gam in sorted(new):
        c = COLORS.get(gam)
        if gam in old:
            ax.plot(old[gam]['times'], old[gam]['alphas'], color=c, lw=1.6, ls='--',
                    alpha=0.65,
                    label='original scheme' if gam == min(new) else None)
        ax.plot(new[gam]['times'], new[gam]['alphas'], color=c, lw=2.4,
                label=fr'$\gamma$={gam} ($\gamma/\gamma_c$={gam/GAMMA_C:.2f})')
        t_pen, a_s = MILROY[gam]
        ax.scatter([t_pen if t_pen else XMAX], [a_s], color=c, s=220, marker='*',
                   zorder=5, edgecolors='k', linewidths=0.6,
                   label='Milroy 1999' if gam == min(new) else None)
    ax.axhline(1.0, color='k', ls=':', alpha=0.5)
    ax.set_xlabel('Time (RMF periods)')
    ax.set_ylabel(r'Penetration Factor $\alpha$')
    ax.set_title(r'$\lambda$=11.07, $N_r$=64; dashed = original scheme', fontsize=22)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.2)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    print(f'Saved {fname}')


def fig_field_lines(all_results, gam):
    res = all_results[gam]
    r_plot = np.linspace(0.0, 1.0, 200)
    theta = np.linspace(0, 2*np.pi, 200)
    Xm = np.outer(np.cos(theta), r_plot)
    Ym = np.outer(np.sin(theta), r_plot)
    circ = np.linspace(0, 2*np.pi, 400)

    fig, axes = plt.subplots(1, len(PLOT_TIMES), figsize=(18, 8))
    for ax, t in zip(np.atleast_1d(axes), PLOT_TIMES):
        A, b, r, tau = res['snapshots'][t]
        Ai = np.interp(r_plot, r, A.real) + 1j*np.interp(r_plot, r, A.imag)
        # co-rotating frame: A carries e^{-i tau}, so undo it with e^{+i tau}
        psi = np.real(Ai[None, :]*np.exp(1j*(theta[:, None] + tau)))
        ax.contour(Xm, Ym, psi, levels=30, colors='blue', linewidths=0.8)
        ax.plot(np.cos(circ), np.sin(circ), 'k-', lw=2)
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.set_aspect('equal')
        ax.set_title(fr'$t$={t}T, $\alpha$={abs(b[-1]-b[0]):.3f}', fontsize=24)
        ax.set_xlabel('$x/R$')
        ax.set_ylabel('$y/R$')
    fig.suptitle(fr'RMF Penetration - $\gamma$={gam}, $\lambda$={LAM}, '
                 fr'$N_r$={res.get("Nr", "?")}', fontsize=36)
    plt.tight_layout()
    fname = f'field_lines_gam{gam}.pdf'
    plt.savefig(fname, bbox_inches='tight')
    print(f'Saved {fname}')


def fig_profiles(all_results, gam):
    res = all_results[gam]
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    axes = axes.flatten()
    for ax, t in zip(axes, PLOT_TIMES):
        A, b, r, tau = res['snapshots'][t]
        g = gamma_of(res, gam)
        dA = dA_dr(A, r, tau, g)
        with np.errstate(divide='ignore', invalid='ignore'):
            Br = np.abs(A)/r
        Br[0] = abs(dA[0])                       # |B_r| = |B_theta| on axis for n=1
        ax.plot(r, Br/g, 'b-', lw=2, label=r'$|B_r|/B_\omega$')
        ax.plot(r, np.abs(dA)/g, 'g--', lw=2, label=r'$|B_\theta|/B_\omega$')
        # B_z is an order of magnitude larger than the RMF once penetrated
        # (lambda^2/gamma = 7.4 here), so it needs its own axis
        axz = ax.twinx()
        axz.plot(r, LAM**2*b/g, 'r:', lw=2, label=r'$B_z/B_\omega$')
        axz.set_ylim(-LAM**2/g*1.05, 0.05*LAM**2/g)
        axz.tick_params(axis='y', colors='r', labelsize=20)
        ax.set_xlabel(r'$r/R$')
        ax.set_title(fr'$t$={t}T, $\alpha$={abs(b[-1]-b[0]):.3f}')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 2.0)
        ax.tick_params(axis='both', labelsize=20)
    handles, labels = axes[0].get_legend_handles_labels()
    hz, lz = axz.get_legend_handles_labels()
    handles, labels = handles + hz, labels + lz
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5),
               fontsize=24)
    fig.supylabel(r'$|B_\perp|/B_\omega$   (right axis: $B_z/B_\omega$)', fontsize=22)
    fig.suptitle(fr'Radial Profiles - $\gamma$={gam} '
                 fr'($\gamma/\gamma_c$={gam/GAMMA_C:.2f}), $\lambda$={LAM}')
    plt.tight_layout()
    fname = f'radial_profiles_gam{gam}.pdf'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    print(f'Saved {fname}')


if __name__ == '__main__':
    np.set_printoptions(precision=2)
    plot_options['font.size'] = 24
    mpl.rcParams.update(plot_options)

    all_results = load(PICKLE)
    print(f'gamma_c({LAM}) = {GAMMA_C:.3f}')
    fig_alpha(all_results)
    try:
        fig_alpha_compare(all_results, load('all_results_legacy_Nr64.pkl'))
    except FileNotFoundError:
        print('no all_results_legacy_Nr64.pkl; skipping the before/after figure')
    for gam in sorted(all_results):
        fig_field_lines(all_results, gam)
        fig_profiles(all_results, gam)
    plt.close('all')
    print('Done.')
