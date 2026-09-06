# NumericalRMF

Rotating magnetic field (RMF) current drive in a cylindrical plasma column, in
the fixed-ion model of Milroy 1999 Sec. III.A.

The goal of the project is to reproduce Milroy's Sec. III.A results (his Figs. 1,
2, 4 and 6) with an independent implementation, and then apply the same solver to
our own device. **Milroy 1999 is the reference for this project**; the other
papers in `articles/` are context — Hugrass & Grimm 1981 is the original
initial-value calculation, Jones & Hugrass 1981 and Hugrass 1985 give the
steady-state analysis and an analytic limit worth testing against.

## The model

An infinitely long plasma cylinder of radius `R`, surrounded by vacuum. Ions are
immobile with uniform density `n` (`u = 0`), electrons are massless, resistivity
is uniform and isotropic, `∇P_e = 0`. Antennae impose an external rotating field

```
B_rmf = B_w cos(wt) x̂ + B_w sin(wt) ŷ,        w_ci << w << w_ce
```

so the electrons are magnetised by the RMF while the ions are not, and an
azimuthal electron current is driven. With those assumptions the system is
Ohm's law plus Maxwell,

```
dA/dt = u x B - eta J - (1/en)(J x B - grad P_e),   B = curl A,   J = curl B / mu_0
```

### Normalisation

Milroy's scaling (his Sec. II.B): `t_0 = 1/w`, `r_0 = R`, `A_0 = R B_w / gamma`,
`B_0 = A_0/R`, `J_0 = B_0/(mu_0 R)`. Two dimensionless numbers control everything:

```
lambda = R/delta = R sqrt(mu_0 w / 2 eta)     delta = classical skin depth
gamma  = w_ce/nu_ei = (1/e) B_w/(n eta)       RMF strength
```

and the equations collapse to Milroy Eq. (10),

```
dA/dt = -(1/2 lambda^2) ( J + J x B )
```

### Reduction to 1-D

Keeping the `n = 0` component of `B_z` and the `n = 1` component of `A_z` — the
truncation Milroy uses throughout Sec. III.A — and writing

```
A = 2 A_z1        (Hugrass's convention: the BC becomes A + A' = 2 gamma e^{-it})
b = B_z/lambda^2  (so alpha = b(1) - b(0) with no extra factors)
```

gives the closed pair actually integrated in `rmf_solver.py`:

```
dA/dt = L[A]/(2 lambda^2)  -  i A b' / (2 r)

db/dt = (1/(2 lambda^2 r)) d/dr [ r b' + P/(2 lambda^2) ],   P = Im( L[A] conj(A) )

L[A] = A'' + A'/r - A/r^2 = d/dr[ (1/r) d/dr (r A) ]
```

Boundary and initial conditions:

```
A(0) = 0            n = 1 flux vanishes on axis
A(1) + A'(1) = 2 gamma e^{-it}    vacuum matching, Milroy Eq. (6)/(13)
b'(0) = 0           symmetry
b(1) = B_z0 = 0     applied axial bias field
A = b = 0 at t = 0
```

The Robin condition at `r = 1` is the one that eliminates the unknown internal
screening currents: outside the plasma `A_z,n = alpha_n (r/R)^n + beta_n
(r/R)^-n`, and `A_z,n(R) + (R/n) dA_z,n/dr|_R = 2 alpha_n` cancels `beta_n`,
leaving only the known antenna drive `alpha_1 = (R B_w/2) e^{-iwt}`.

### Diagnostic

```
alpha = (2/(mu_0 n e w R^2)) ( B_z(R) - B_z(0) )  =  b(1) - b(0)
```

the driven azimuthal current normalised to what it would be if every electron
rotated synchronously with the RMF. `alpha_s` is its steady value (the last value
of a converged run); `t_pen` is when `alpha` first reaches 95% of `alpha_s`.

## Numerics

- Uniform radial mesh, `r_i = i dr`, `dr = 1/(Nr-1)`; axis included.
- `L[A]` in the flux form `d/dr[(1/r) d/dr (rA)]`, second order and exact for
  both `A ~ r` and `A ~ r^3`. The plain 3-point form is only first order in the
  max norm because of the point next to the axis.
- Robin condition at `r = 1` via a ghost node at `r = 1+dr` with a centred
  derivative, and the boundary node evolved with the PDE. Second order.
- `b` equation in finite-volume flux form `G = r b' + P/(2 lambda^2)`, which is
  conservative and handles the axial half-cell without a special limit — the
  treatment Hugrass & Grimm describe as "assuming A_z, B_z are constant over the
  first radial zone and using Gauss's theorem".
- Time integration: explicit Heun (RK2 predictor-corrector), second order.
  Stability, not accuracy, sets the step: `dt = 0.002` at `Nr = 64` matches
  `dt = 0.001` to five digits, `dt = 0.003` blows up, and the limit tightens as
  `gamma` grows (`0.002` fails near `gamma = 20`).
- Production: `lambda = 11.07`, `Nr = 64`, `dt = 0.002`, 200 RMF periods,
  ~30 s per `gamma`.

## Layout

```
rmf_solver.py          the solver
paths.py               where outputs go (results/ and figures/)
run_milroy.py          the three reference gammas
plot.py                alpha(t), before/after, field lines, radial profiles
validate_numerics.py   analytic linear limit, ramp cost, operator order,
                       ablation, grid convergence
study_thresholds.py    ascending/descending gamma sweeps for gamma_c
utils.py               matplotlib style
nr32_fix.py            the original script, kept as the provenance of
                       results/all_results_Nr*.pkl

studies/               numbered studies (see below)
animation/             mp4 generation
articles/              the source papers
results/               pickles and tabulated data       (all outputs land here)
figures/               pdf and png
figures/animations/    mp4
```

Outputs are written through `paths.py`, which anchors them to the repo root, so
scripts can be run from any directory.

```
python run_milroy.py
python plot.py
python validate_numerics.py
python study_thresholds.py
```

## Studies

### 1 — `studies/study1_alpha_penetration/` : alpha vs time, penetration tests

alpha(t) for gamma = 14.9, 16.6, 18.2, in four configurations that differ only
in the discretisation and the turn-on ramp, so each effect can be isolated:

| script | discretisation | `rise_time` |
|---|---|---|
| `1_legacy_ramp0.py` | original (`nr32_fix.py`) | 0 |
| `2_legacy_ramp3.py` | original | 3 T |
| `3_fixed_ramp0.py` | corrected (`rmf_solver.py`) | 0 |
| `4_fixed_ramp3.py` | corrected | 3 T |
| `5_all_tests.py` | runs all four and overlays them, one panel per gamma | |

Each writes `results/alpha_<label>.pkl` and `figures/alpha_<label>.{pdf,png}`,
and prints alpha_s and t_pen against Milroy's numbers. Common options `--nr`,
`--dt`, `--periods`, `--no-figure`, plus on the first four

```
python 1_legacy_ramp0.py --dt-scan 0.004 0.003 0.002 0.001 --periods 60
```

which reports, per gamma, whether the run stayed finite and at which period it
blew up. The explicit scheme fails *after* the RMF has penetrated, not at
start-up, so a scan shorter than ~60 periods reports everything as stable.

### 2 — `studies/study2_gamma_c_vs_lambda/` : Milroy Fig. 4

`gamma_c_scan.py` finds the critical gamma for full penetration at each lambda
by bisection on cold starts (corrected scheme, `rise_time = 0`), and
`plot_gamma_c.py` plots it against Milroy's Eqs. (14) and (15).

Resolution and time step are chosen per lambda (`Nr ~ 8 lambda` to resolve the
skin depth, `dt` from the Hall/whistler limit `~dr^2`), and the time budget is
`--tmax-factor * lambda^2` RMF periods, which makes the finite-budget bias in
gamma_c independent of lambda — about 0.6% high at the default.

The cost goes as **lambda^4** (`Nr ~ lambda`, `dt ~ dr^2`, `T_max ~ lambda^2`),
so start with

```
python gamma_c_scan.py --dry-run
```

which prints the plan and the step count per lambda without running anything.
lambda <= 20 is comfortable, 25-40 runs into hours, and lambda >= 50 is out of
reach with the explicit integrator — that is what the semi-implicit scheme
(TODO 7b) would buy. Results are appended to `results/gamma_c_vs_lambda.csv` as
they are produced and `--resume` skips lambdas already there, so the scan can be
stopped and restarted:

```
python gamma_c_scan.py --lam-max 20
python plot_gamma_c.py
```

## Where it stands

`lambda = 11.07`, `Nr = 64`, both schemes run to 200 T so the comparison is like
for like:

| gamma | gamma/gamma_c | alpha_s (orig / now / Milroy) | t_pen (orig / now / Milroy) |
|---|---|---|---|
| 14.9 | 0.985 | 0.387 / **0.400** / 0.42 | never |
| 16.6 | 1.097 | 0.990 / **0.992** / 0.98 | 46.8 T / **40.9 T** / 40 T |
| 18.2 | 1.203 | 0.992 / **0.994** / 0.98 | 29.8 T / **24.4 T** / 27 T |

(`t_pen` at 95% of `alpha_s`. Figure: `figures/alpha_vs_time_comparison.pdf`.)

Note that `gamma_c = 15.13` here is Milroy's Eq. (15) penetration threshold, not
`1.12 lambda = 12.40`, which is his Eq. (14) value for *expulsion* of an already
penetrated RMF. On that scale `gamma = 14.9` sits 1.5% *below* threshold, which
is why it never penetrates.

## Validation

- **Analytic linear limit.** As `gamma -> 0` the Hall term drops and the steady
  n=1 flux is `A = 2 gamma I1(kr)/(k I0(k))` with `k = (1-i) lambda` (Hugrass
  1985 Eqs. 20-21). The solver matches it to 0.04% at `Nr = 64` and converges at
  second order.
- **Threshold.** Milroy's staircase procedure puts the penetration threshold
  between `gamma = 15.25` and `15.50`, against 15.13 from his Eq. (15).
- **alpha_s(gamma).** The ascending sweep (`results/alpha_s_up.txt`) tracks his Eq. (18)
  fit within the ~20% he claims for it.
- **Grid convergence.** `alpha_s` at `gamma = 14.9` is converged to four digits
  by `Nr = 32`.

## References

- R. D. Milroy, Phys. Plasmas **6**, 2771 (1999) — the reference for this project
- W. N. Hugrass and R. C. Grimm, J. Plasma Phys. **26**, 455 (1981)
- I. R. Jones and W. N. Hugrass, J. Plasma Phys. **26**, 441 (1981)
- W. N. Hugrass, Aust. J. Phys. **38**, 157 (1985)
- R. D. Milroy, Phys. Plasmas **7**, 4135 (2000) — 2-D MHD successor, moving ions
