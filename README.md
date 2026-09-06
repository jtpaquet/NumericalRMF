# NumericalRMF

Rotating-magnetic-field current drive in a plasma column: the fixed-ion,
n = 0 / n = 1 model of Milroy, Phys. Plasmas **6**, 2771 (1999), Sec. III.A
(following Hugrass & Grimm, J. Plasma Phys. **26**, 455 (1981)).

```
python run_milroy.py        # the three reference gammas -> all_results_corrected.pkl
python plot.py              # alpha(t), field lines, radial profiles
python validate_numerics.py # convergence and ablation study
python study_thresholds.py  # gamma_c and the expulsion threshold
```

At lambda = 11.07 the code gives alpha_s = 0.40 / 0.99 / 0.99 and full
penetration at -- / 40 T / 24 T for gamma = 14.9 / 16.6 / 18.2, against
Milroy's 0.42 / 0.98 / 0.98 and -- / 40 T / 27 T. The measured penetration
threshold falls between gamma = 15.25 and 15.50, against 15.13 from his
Eq. (15). See [COMPARISON.md](COMPARISON.md).
