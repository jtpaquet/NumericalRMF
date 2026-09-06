#!/usr/bin/env python3
"""alpha(t) for the three Milroy gammas with the ORIGINAL scheme, no ramp.

Scheme is exactly nr32_fix.py:
  * Robin BC by backward difference:  A[-1] = (2 A_ext + A[-2]/dr) / (1 + 1/dr)
  * L[A] as the plain 3-point  A'' + A'/r - A/r^2
  * b equation Hall term as the live line of nr32_fix.py,
        hall_B = 0.5 * d_product / (2 * lam4 * r)     i.e.  dP/dr / (4 lam4 r)
    (the line above it in that file, dP/dr / (lam4 r), is overwritten and never
    runs; 1/(4 lam4) is the coefficient consistent with A = 2 A_z1)
  * b on axis as  rhs_B[0] = 2 (B[1]-B[0]) / dr^2 / (2 lam2), no Hall source

rise_time = 0: the RMF is switched on as a step at t = 0, as Milroy does.
This is the case to check for the instability that motivated the ramp.

    python 1_legacy_ramp0.py
    python 1_legacy_ramp0.py --dt-scan 0.004 0.003 0.002 0.001 --periods 60
"""
from _stability_run import main

if __name__ == '__main__':
    main(label='legacy_ramp0', legacy=True, rise_time=0.0,
         description='Original scheme, no turn-on ramp')
