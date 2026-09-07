#!/usr/bin/env python3
"""alpha(t) for the three Milroy gammas with the ORIGINAL scheme, 3-period ramp.

Same discretisation as 1_legacy_ramp0.py (see its docstring), but with the
turn-on ramp gamma(t) = gamma (1 - exp(-t/tau_r)) at tau_r = 3 T, which is what
nr32_fix.py used. This reproduces the stored all_results_Nr32_6283steps.pkl and
all_results_Nr64_6283steps_fine.pkl.

For reference, Hugrass & Grimm 1981 Eq. (3) is where the ramp comes from, but
their tau_r = 0.4 us at omega = 5e6 /s is omega*tau_r = 2.0, i.e. 0.32 T.

    python 2_legacy_ramp3.py
"""
from _stability_run import main
from rmf_solver import TWO_PI

if __name__ == '__main__':
    main(label='legacy_ramp3', legacy=True, rise_time=3*TWO_PI,
         description='Original scheme, 3-period turn-on ramp')
