# Comparison with Milroy 1999, Sec. III.A

R. D. Milroy, *A numerical study of rotating magnetic fields as a current drive
for field reversed configurations*, Phys. Plasmas **6**, 2771 (1999).

Reference case: R = 10 cm, n = 0.333e20 m^-3, B_w = 90/100/110 G,
omega = 2.2e5 s^-1, T_e = 5 eV  ->  lambda = 11.07, gamma = 14.9/16.6/18.2.

## What the code actually integrates

Milroy's Eqs. (10)-(13), truncated to the n = 0 component of B_z and the n = 1
component of A_z (the truncation he uses for all of Sec. III.A). Writing

    A = 2 A_z1        so the vacuum matching condition is A + A' = 2 gamma e^{-i t}
    b = B_z / lambda^2    so alpha = b(1) - b(0) with no extra factors

the closed pair is

    dA/dt = L[A]/(2 lambda^2)  -  i A b' / (2 r)
    db/dt = (1/(2 lambda^2 r)) d/dr [ r b' + P/(2 lambda^2) ],   P = Im(L[A] conj(A))
    L[A]  = A'' + A'/r - A/r^2 = d/dr[(1/r) d/dr (r A)]

with A(0) = 0, b'(0) = 0, b(1) = B_z0 = 0 and A = b = 0 at t = 0.

The factor 2 in `A` and the compensating 1/4 in the Hall term of the b equation
were already present (but undocumented) in `nr32_fix.py`, where the boundary
condition reads `2*gamma` and the Hall term carries `0.5*dP/(2*lam4*r)`. Those
two cancel, so the original scheme was self-consistent; the differences below
are all discretisation, not model.

## Results

`lambda = 11.07`, `Nr = 64`, `dt = 0.002`. Penetration time is the first time
alpha reaches 90% of its final value.

| gamma | gamma/gamma_c | Milroy alpha_s | this code alpha_s | Milroy t_pen | old code t_pen | this code t_pen |
|------:|--------------:|---------------:|------------------:|-------------:|---------------:|----------------:|
| 14.9  | 0.985         | 0.42           | 0.400             | never        | never          | never           |
| 16.6  | 1.097         | 0.98           | 0.992             | 40 T         | 47 T           | 40.1 T          |
| 18.2  | 1.203         | 0.98           | 0.994             | 27 T         | 30 T           | 23.8 T          |

Milroy's own empirical fit, Eq. (17) `tau_P = lambda^2 / (2 sqrt(gamma_N))`,
gives 31 T and 22 T for these two cases, so both his numbers and ours sit ~10-30%
above his own fit; that fit spans lambda = 6.5 to 100 and is only claimed to be
"a reasonable approximation".

The `gamma = 18.2` case is the one that now overshoots (23.8 T vs 27 T). Its
alpha(t) rises almost vertically through 0.9, so the number is sensitive to
exactly what "full penetration" means when read off Fig. 2 by eye; the same
threshold applied to `gamma = 16.6` lands on 40 T exactly.

### Threshold

Repeating Milroy's staircase procedure (hold gamma until alpha settles, step up,
repeat) puts the penetration threshold between gamma = 15.25 (alpha_s = 0.52)
and gamma = 15.50 (alpha_s = 0.99). Milroy's Eq. (15) gives

    gamma_c = 1.12 lambda (1 + 0.12 (lambda - 6.5)^0.4) = 15.13

so the code reproduces his threshold to within one sweep step (~1-2%). See
`study_thresholds.py`.

## What changed, and how much each piece was worth

Ablation at `gamma = 16.6`, `Nr = 64`, one change at a time from the original
scheme (`validate_numerics.py`):

| scheme | alpha_final | t(0.95 alpha_s) |
|---|---:|---:|
| original (`nr32_fix.py`) | 0.9898 | 47.7 T |
| + 2nd-order Robin BC at r = 1 | 0.9915 | 46.9 T |
| + finite-volume b equation (axis) | 0.9902 | 47.7 T |
| + flux-form `L[A]` | 0.9898 | 47.7 T |
| + no gamma ramp | 0.9900 | 41.7 T |
| all of the above | 0.9918 | 40.9 T |

**1. The `gamma` ramp is the whole timing discrepancy.** `nr32_fix.py` switches
the RMF on as `gamma (1 - exp(-t/3T))`. Milroy applies `alpha_1 = (R B_w/2)
e^{-i omega t}` as a step at t = 0 with B = 0 everywhere. The ramp costs about
6 T on both gamma = 16.6 and gamma = 18.2 — which is the entire 47 vs 40 gap.
It is also visible in the old figure: the alpha(t) curves are flat for the
first ~5 T before they start to climb, and Milroy's rise from t = 0.

**2. First-order boundary condition at r = 1.** The Robin condition was applied
as `A[-1] + (A[-1]-A[-2])/dr = 2 gamma e^{-i t}`, which is O(dr). The skin depth
is 1/lambda = 0.09, so at Nr = 64 that is a ~3% error in the drive amplitude at
the wall — checking the snapshots, `|B_theta(R)|/B_w` computed one-sided is
1.725 where the second-order value is 1.777. Because the penetration time
diverges as (gamma - gamma_c)^(-1/2), a percent-level error in the effective
drive is a several-percent error in t_pen. Fixed with a ghost node at r = 1+dr
and a centred derivative.

**3. The axis in the b equation.** `rhs_B[0] = 2 (b1-b0)/dr^2 / (2 lambda^2)`
is half the correct axial limit of the cylindrical Laplacian (`4 (b1-b0)/dr^2`),
and the Hall source was dropped at r = 0 although `(1/r) dP/dr -> P''(0)` is
finite there. Since alpha is measured *at* r = 0, both errors bias alpha low.
Recast as a finite-volume flux `G = r b' + P/(2 lambda^2)`, which is exactly
conservative and gets the axial half-cell right by construction.

**4. `L[A]` is only first order as written.** The plain 3-point form
`A'' + A'/r - A/r^2` is exact for A ~ r but not for A ~ r^3, and the O(dr) error
at the point next to the axis does not average out:

| Nr | 3-point | flux form |
|---:|---:|---:|
| 32 | 2.46e-01 | 3.15e-02 |
| 64 | 1.23e-01 | 7.63e-03 |
| 128 | 6.14e-02 | 1.88e-03 |
| 256 | 3.06e-02 | 4.66e-04 |

(max error of L[A] for A = J1(5r)). The flux form `d/dr[(1/r) d/dr (rA)]` is
second order and exact for both r and r^3.

Items 2-4 barely move the penetration *time* but they dominate the subcritical
steady state, which is where the old numbers were furthest from Milroy:

| Nr | original scheme | corrected |
|---:|---:|---:|
| 32 | 0.3682 | 0.4005 |
| 64 | 0.3867 | 0.4004 |
| 128 | 0.3941 | 0.4004 |

(alpha_s at gamma = 14.9, run to 250 T.) Both schemes head for the same answer,
but the corrected one is grid-converged to four digits at Nr = 32 while the
original is still moving at Nr = 128 and would need Nr of a few hundred to get
there. The stored Nr = 32 and Nr = 64 runs agreeing with each other on the
*penetrating* cases is therefore not evidence of convergence on this one.

**5. alpha_s was not run to steady state.** At gamma = 14.9 alpha is still
climbing at 101 T; it settles by ~200 T. The stored `all_results_Nr32_*.pkl`
value of 0.367 is a 101 T value, not alpha_s.

## Plotting corrections

- `gamma_c` was taken as `1.12*lambda = 12.40`. That is Milroy's Eq. (14), the
  threshold at which an already-penetrated RMF is *expelled*. The penetration
  threshold for lambda > 6.5 is Eq. (15), 15.13. With the right value the legend
  reads gamma/gamma_c = 0.98 / 1.10 / 1.20 instead of 1.20 / 1.34 / 1.47, and
  the reason gamma = 14.9 never penetrates becomes obvious: it is *below*
  threshold, by 1.5%.
- the Milroy target marker for gamma = 14.9 was 0.48. The paper quotes 0.42;
  0.48 = 1.6/sqrt(lambda) is the gamma -> gamma_c limit of his Eq. (18) fit.
- `B_z/B_w` was plotted as `b/gamma`, but `b` is stored as `B_z/lambda^2`, so
  the axial field came out lambda^2 = 123 times too small (a flat line at the
  bottom of the old profile figures). It should be `lambda^2 b / gamma`, which
  reaches -7.3 at full penetration — the physically interesting statement that
  the driven axial field is 7x the RMF that drives it.
- the field-line plots use `Re[A e^{i(theta - tau)}]`. `A` already carries the
  `e^{-i tau}` of the drive, so this leaves the picture rotating backwards at
  2 omega instead of standing still in the RMF frame. Milroy plots in the
  co-rotating frame, which is `Re[A e^{i(theta + tau)}]`.
- `|B_theta(R)|` was a one-sided difference. The boundary condition gives it
  exactly: `A'(1) = 2 gamma e^{-i tau} - A(1)`.

## Cost

- vectorising the two `for i in range(1, Nr-1)` loops: 5.7x at Nr = 64
  (536 -> 94 us/step).
- `dt = 0.002` gives answers identical to `dt = 0.001` to five digits at
  Nr = 64. `dt = 0.003` is unstable for the corrected scheme (the corrected
  axial coefficient is twice the old one, so it is slightly more restrictive),
  and the limit tightens as gamma grows — `dt = 0.002` fails around gamma = 20.
- together: ~11x, i.e. ~30 s per gamma rather than ~5 min.

## What is left

The remaining differences (alpha_s = 0.400 vs 0.42 at gamma = 14.9, 23.8 T vs
27 T at gamma = 18.2) are at the level of what can be read off a 1999 figure,
and gamma = 14.9 sits at gamma/gamma_c = 0.985 where alpha_s is on the steepest
part of the curve — our own measured gamma_c differs from his Eq. (15) fit by
about a step of 0.25, and that alone accounts for the 0.02.

The one modelling knob still closed is the n = 1 truncation. The quadratic Hall
term generates an n = 2 component of B_z, which feeds back on n = 1 at second
order; Milroy keeps it available (up to 42 modes for the anisotropic and
finite-coil cases) but states that "for most of the problems studied in this
paper, only the n = 1 mode needs to be included". Opening it up is a model
extension rather than a numerical fix, so it is deliberately not done here.
