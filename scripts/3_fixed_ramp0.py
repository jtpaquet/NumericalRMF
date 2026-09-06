#!/usr/bin/env python3
"""alpha(t) for the three Milroy gammas with the CORRECTED scheme, no ramp.

Scheme is rmf_solver.py's default:
  * Robin BC via a ghost node at r = 1+dr with a centred derivative, and the
    boundary node evolved with the PDE (second order)
  * L[A] in the flux form d/dr[(1/r) d/dr (rA)] (second order, exact for r and r^3)
  * b equation as the conservative finite-volume flux G = r b' + P/(2 lam^2),
    which gets the axial half-cell and its Hall source right by construction

rise_time = 0 matches Milroy, who applies alpha_1 = (R B_w/2) e^{-i omega t} as a
step at t = 0. This is the configuration behind the numbers in COMPARISON.md.

    python 3_fixed_ramp0.py
"""
from _stability_run import main

if __name__ == '__main__':
    main(label='fixed_ramp0', legacy=False, rise_time=0.0,
         description='Corrected scheme, no turn-on ramp (matches Milroy)')
