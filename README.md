# NumericalRMF

Rotating-magnetic-field current drive in a plasma column: the fixed-ion,
n = 0 / n = 1 model of Hugrass & Grimm, J. Plasma Phys. **26**, 455 (1981) and
Milroy, Phys. Plasmas **6**, 2771 (1999), Sec. III.A.

```
python run_milroy.py        # the three reference gammas -> all_results_corrected.pkl
python plot.py              # alpha(t), before/after, field lines, radial profiles
python validate_numerics.py # analytic linear limit, ramp cost, convergence, ablation
python study_thresholds.py  # gamma_c and the expulsion threshold
```

At lambda = 11.07, against Milroy's Sec. III.A numbers:

| gamma | alpha_s (orig / now / Milroy) | t_pen (orig / now / Milroy) |
|---|---|---|
| 14.9 | 0.387 / **0.400** / 0.42 | never |
| 16.6 | 0.990 / **0.992** / 0.98 | 46.8 T / **40.1 T** / 40 T |
| 18.2 | 0.992 / **0.994** / 0.98 | 29.8 T / **23.8 T** / 27 T |

The measured penetration threshold falls between gamma = 15.25 and 15.50,
against 15.13 from Milroy's Eq. (15). In the gamma -> 0 limit the solver matches
the analytic solution of Hugrass 1985 Eq. (21) to 0.04% at Nr = 64.

See [COMPARISON.md](COMPARISON.md) for the full comparison and for the
prioritised list of where this code departed from the published models.
