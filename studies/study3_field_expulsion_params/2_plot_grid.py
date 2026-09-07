#!/usr/bin/env python3
"""Study 3 - heatmaps of the (pressure, T_e) grid from 1_scan_grid.py.

    figures/study3_gamma_ratio.pdf     how deep into the linear regime we are
    figures/study3_settling_time.pdf   predicted settling time vs the 10/25 ms
                                        target
    figures/study3_screening.pdf       steady B_r(0)/B_w, B_r(R)/B_w vs the
                                        measured ~0.17, ~0.27-0.33
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import paths
import device as dev
from utils import plot_options

MEASURED_WALL_RATIO = dev.MEASURED_B_WALL_G/dev.B_W_G
MEASURED_CENTER_RATIO = dev.MEASURED_B_CENTER_G/dev.B_W_G
TARGET_TIMES_MS = (10.0, 25.0)


def load_grid(fname='study3_grid.csv'):
    with open(paths.results(fname)) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            r[k] = float(v)
    ps = sorted(set(r['p_mtorr'] for r in rows))
    Tes = sorted(set(r['Te_eV'] for r in rows))
    grid = {(r['p_mtorr'], r['Te_eV']): r for r in rows}
    return ps, Tes, grid


def as_array(ps, Tes, grid, key):
    Z = np.empty((len(Tes), len(ps)))
    for i, Te in enumerate(Tes):
        for j, p in enumerate(ps):
            Z[i, j] = grid[(p, Te)][key]
    return Z


def heatmap(ax, ps, Tes, Z, cmap='viridis', log=False, contours=(), clabel_fmt='%g',
            vmin=None, vmax=None):
    Zp = np.log10(Z) if log else Z
    im = ax.pcolormesh(ps, Tes, Zp, shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('fill pressure (mTorr)')
    ax.set_ylabel(r'$T_e$ (eV)')
    for c in contours:
        cval = np.log10(c) if log else c
        cs = ax.contour(ps, Tes, Zp, levels=[cval], colors='white', linewidths=2.0)
        ax.clabel(cs, fmt=clabel_fmt)
    return im


def main():
    plot_options['font.size'] = 16
    mpl.rcParams.update(plot_options)
    ps, Tes, grid = load_grid()

    # ------------------------------------------------------------- gamma/gamma_c
    fig, ax = plt.subplots(figsize=(8, 6))
    ratio = as_array(ps, Tes, grid, 'gam_over_gc')
    im = heatmap(ax, ps, Tes, ratio, cmap='magma', log=True)
    cb = fig.colorbar(im, ax=ax, label=r'$\log_{10}(\gamma/\gamma_c)$')
    ax.set_title(r'Depth into the linear regime ($\gamma \ll \gamma_c$ everywhere)'
                 '\n' r'$B_w$ = %.0f G, f = %.0f kHz' % (dev.B_W_G, dev.F_RMF/1e3))
    out = paths.figure('study3_gamma_ratio.pdf')
    plt.tight_layout()
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')

    # -------------------------------------------------------------- settling time
    fig, ax = plt.subplots(figsize=(8, 6))
    t_settle = as_array(ps, Tes, grid, 't_settle_ms')
    im = heatmap(ax, ps, Tes, t_settle, cmap='viridis', log=True,
                 contours=[10.0, 25.0], clabel_fmt='%g ms')
    cb = fig.colorbar(im, ax=ax, label=r'$\log_{10}(t_{settle}$ / ms$)$')
    ax.set_title('Predicted field-settling time vs. the observed 10-25 ms')
    out = paths.figure('study3_settling_time.pdf')
    plt.tight_layout()
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')
    print(f'  model range: {t_settle.min():.4f} - {t_settle.max():.4f} ms '
          f'over the whole (p, Te) box  (observed: 10-25 ms)')

    # ---------------------------------------------------------------- screening
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    Br0 = as_array(ps, Tes, grid, 'Br0_over_Bw')
    BrR = as_array(ps, Tes, grid, 'BrR_over_Bw')
    im0 = heatmap(axes[0], ps, Tes, Br0, cmap='cividis',
                  contours=[MEASURED_CENTER_RATIO], clabel_fmt='measured %.2f')
    axes[0].set_title(r'steady $|B_r(0)|/B_w$')
    fig.colorbar(im0, ax=axes[0])
    imR = heatmap(axes[1], ps, Tes, BrR, cmap='cividis',
                  contours=[MEASURED_WALL_RATIO], clabel_fmt='measured %.2f')
    axes[1].set_title(r'steady $|B_r(R)|/B_w$')
    fig.colorbar(imR, ax=axes[1])
    fig.suptitle('Steady screening vs. the measured center/wall fields with plasma')
    out = paths.figure('study3_screening.pdf')
    plt.tight_layout()
    plt.savefig(out, bbox_inches='tight')
    print(f'saved {out}')
    plt.close('all')


if __name__ == '__main__':
    main()
