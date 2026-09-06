#!/usr/bin/env python3
"""alpha(t) for the three Milroy gammas with the CORRECTED scheme, 3-period ramp.

Same discretisation as 3_fixed_ramp0.py, with tau_r = 3 T. Comparing this against
3_fixed_ramp0.py isolates the cost of the ramp from the cost of the
discretisation, and against 2_legacy_ramp3.py isolates the discretisation at
fixed ramp.

    python 4_fixed_ramp3.py
"""
from _stability_run import main
from rmf_solver import TWO_PI

if __name__ == '__main__':
    main(label='fixed_ramp3', legacy=False, rise_time=3*TWO_PI,
         description='Corrected scheme, 3-period turn-on ramp')
