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

## Files

| file | what it is |
|---|---|
| `rmf_solver.py` | the solver, plus `gamma_c()` and `tau_penetration()` from Milroy's Eqs. (15) and (17). `legacy=True` reproduces the original `nr32_fix.py` scheme for comparison |
| `run_milroy.py` | the three reference gammas -> `all_results_corrected.pkl` |
| `plot.py` | alpha(t), before/after, field lines, radial profiles |
| `validate_numerics.py` | analytic linear limit, ramp cost, operator order, ablation, grid convergence |
| `study_thresholds.py` | ascending/descending gamma sweeps for `gamma_c` and the expulsion threshold |
| `utils.py` | matplotlib style |
| `nr32_fix.py` | the original script, kept as the provenance of `all_results_Nr*.pkl` |
| `COMPARISON.md` | full comparison with Milroy, and where this code departed from the published models |
| `TODO.md` | what is left |
| `articles/` | the source papers |

```
python run_milroy.py
python plot.py
python validate_numerics.py
python study_thresholds.py
```

## Where it stands

`lambda = 11.07`, `Nr = 64`, both schemes run to 200 T so the comparison is like
for like:

| gamma | gamma/gamma_c | alpha_s (orig / now / Milroy) | t_pen (orig / now / Milroy) |
|---|---|---|---|
| 14.9 | 0.985 | 0.387 / **0.400** / 0.42 | never |
| 16.6 | 1.097 | 0.990 / **0.992** / 0.98 | 46.8 T / **40.9 T** / 40 T |
| 18.2 | 1.203 | 0.992 / **0.994** / 0.98 | 29.8 T / **24.4 T** / 27 T |

(`t_pen` at 95% of `alpha_s`. Figure: `alpha_vs_time_comparison.pdf`.)

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
- **alpha_s(gamma).** The ascending sweep (`alpha_s_up.txt`) tracks his Eq. (18)
  fit within the ~20% he claims for it.
- **Grid convergence.** `alpha_s` at `gamma = 14.9` is converged to four digits
  by `Nr = 32`.

## References

- R. D. Milroy, Phys. Plasmas **6**, 2771 (1999) — the reference for this project
- W. N. Hugrass and R. C. Grimm, J. Plasma Phys. **26**, 455 (1981)
- I. R. Jones and W. N. Hugrass, J. Plasma Phys. **26**, 441 (1981)
- W. N. Hugrass, Aust. J. Phys. **38**, 157 (1985)
- R. D. Milroy, Phys. Plasmas **7**, 4135 (2000) — 2-D MHD successor, moving ions
