# Comparison with Milroy 1999, Sec. III.A

R. D. Milroy, *A numerical study of rotating magnetic fields as a current drive
for field reversed configurations*, Phys. Plasmas **6**, 2771 (1999), which
follows Hugrass & Grimm, J. Plasma Phys. **26**, 455 (1981) and Jones & Hugrass,
J. Plasma Phys. **26**, 441 (1981); the same reduced n = 0 / n = 1 system is
written out in Hugrass, Aust. J. Phys. **38**, 157 (1985).

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

## Where this code departs from the published models

In priority order, by effect on alpha(t).

**1. The turn-on ramp is real but ~10x too slow.** `nr32_fix.py` used
`gamma (1 - exp(-t/3T))`. That functional form is Hugrass & Grimm 1981 Eq. (3),
`B_t(t) = (1 - e^{-t/tau_r}) B_w`, where tau_r is "characteristic of the RF
source used to generate the rotating field" -- so the ramp is not invented. But
their value is `tau_r = 0.4 us` at `omega = 5e6 /s`, i.e. `omega*tau_r = 2.0`,
which is **0.32 of an RMF period**, not 3 periods. Milroy drops the ramp
entirely and applies `alpha_1 = (R B_w / 2) e^{-i omega t}` as a step at t = 0.
Measured (corrected scheme, Nr = 64):

| tau_r | t_pen (gamma = 16.6) | t_pen (gamma = 18.2) | alpha_s |
|---|---:|---:|---:|
| 0 (Milroy) | 40.1 T | 23.8 T | 0.9918 |
| 0.32 T (Hugrass & Grimm's value) | 40.5 T | 24.2 T | 0.9918 |
| 1 T | 41.7 T | 25.4 T | 0.9919 |
| 3 T (original code) | 46.0 T | 29.6 T | 0.9918 |
| *Milroy 1999* | *40 T* | *27 T* | *0.98* |

At the published value the ramp costs 0.4 T and is irrelevant; at 3 T it costs
6 T. It has no effect at all on alpha_s, as a steady state should not. The run
is perfectly stable with `rise_time = 0`, so nothing is lost by dropping it;
`RISE_TIME_HG1981 = 2.0` is defined in `rmf_solver.py` if you want to keep it.

**2. First-order outer boundary condition.** Now measurable against an exact
solution rather than by argument. Hugrass 1985 Eqs. (20)-(21) give the
`gamma -> 0` limit in closed form: `A = 2 gamma I1(k r) / (k I0(k))` with
`k = (1 - i) lambda`, satisfying the same `A(1) + A'(1) = 2 gamma` of his
Eq. (32). Running the code at `gamma = 1e-3`:

| Nr | original L2 | original \|A(1)\|/exact | corrected L2 | corrected \|A(1)\|/exact |
|---:|---:|---:|---:|---:|
| 16 | 1.64e-01 | 1.396 | 3.23e-02 | 0.996 |
| 32 | 6.39e-02 | 1.170 | 6.43e-03 | 1.002 |
| 64 | 2.82e-02 | 1.080 | 1.43e-03 | 1.000 |
| 128 | 1.32e-02 | 1.039 | 4.34e-04 | 1.000 |

At Nr = 64 -- the production grid -- the original scheme gets the flux at the
wall **8% wrong** and converges at first order. The corrected one is 0.04% and
second order.

**3. The axis.** Hugrass & Grimm say exactly what to do: "The point r = 0
requires special treatment to remove the apparent singularity ... This is
performed by assuming that A_z, B_z are constant over the first radial zone and
using Gauss's theorem" -- i.e. a finite volume over the axial half cell, which
is what `b_scheme='fv'` now does. The original `2 (b1-b0)/dr^2` with no Hall
source is neither that nor the pointwise limit.

**4. alpha_s read before it settled** (101 T; it needs ~200 T at gamma = 14.9).

**5. gamma_c taken as 1.12*lambda**, which is the expulsion threshold (Milroy
Eq. 14), not the penetration threshold (Eq. 15).

**6. Harmonic truncation -- the one real modelling difference left.** The code
keeps only n = 0 in B_z and n = 1 in A_z. That is Hugrass 1985's *steady-state*
approximation ("the second and higher harmonics will be neglected"), and Milroy
says "for most of the problems studied in this paper, only the n = 1 mode needs
to be included". But Hugrass & Grimm's *initial-value* calculations were done on
a full 16x16 polar (r, theta) mesh, so they carried every harmonic. B_z picks up
even harmonics and A_z odd ones, so what is dropped is B_z2 feeding back on A_z1
at second order -- an error that grows with gamma. That is consistent with the
residual: at gamma = 16.6 we land on Milroy's 40 T, at gamma = 18.2 we are 3 T
fast. Opening this up is a model extension, so it is deliberately not done.

**7. Time integration.** Hugrass & Grimm used a time-centred (Crank-Nicolson)
scheme, semi-implicit in the linear terms with SOR iteration on the nonlinear
Hall part, and note they could run "about 10 times the Courant condition at the
origin". Milroy offers the same choice. This code is explicit Heun, so dt is
capped near 0.002 at Nr = 64. Not an accuracy problem, but it is why the runs
cost what they cost.

**8. Outer boundary on B_z.** Hugrass & Grimm close the system with
flux-conserving rings at r = b, so `B_z(a,t)` follows from axial flux
conservation (their Eq. 17), and they show smaller b gives *faster* penetration.
Milroy and Hugrass 1985 instead hold `B_z(R) = B_z0`. This code matches
Milroy -- correct as it stands, but it is the wrong condition if you ever
compare against Hugrass & Grimm's figures directly.

**9. Definition of the penetration time.** Hugrass & Grimm define it by
`alpha(tau) = 0.5 alpha_s` (their Eq. 31), and quote the far-above-threshold
scaling `tau = mu_0 a^2 sigma / 8`, which in these units is `lambda^2/4` = 4.9 T.
Milroy says "full penetration achieved at t = 40T" with no stated criterion.
On the corrected runs `alpha = 0.5 alpha_s` is reached at 24 T (gamma = 16.6)
and 16 T (gamma = 18.2). Comparisons need the same definition on both sides.

**10. Electron inertia** `(m/ne^2) dJ/dt`, which Hugrass & Grimm retain when
`omega <~ nu_ei`, is correctly dropped here: `nu_ei = omega_ce/gamma = 1.1e8 /s`
against `omega = 2.2e5 /s`, so `omega/nu_ei ~ 2e-3`.

## Results

`lambda = 11.07`, `Nr = 64`, `dt = 0.002`. Penetration time is the first time
alpha reaches 90% of its final value.

Both schemes at the same Nr = 64 and the same 200 T, so the comparison is like
for like (`all_results_legacy_Nr64.pkl` vs `all_results_corrected.pkl`; the
figure is `alpha_vs_time_comparison.pdf`).

| gamma | gamma/gamma_c | alpha_s: original / corrected / Milroy | t_pen: original / corrected / Milroy |
|------:|--------------:|---------------------------------------:|-------------------------------------:|
| 14.9  | 0.985         | 0.387 / **0.400** / 0.42               | never / never / never                |
| 16.6  | 1.097         | 0.990 / **0.992** / 0.98               | 46.8 T / **40.1 T** / 40 T           |
| 18.2  | 1.203         | 0.992 / **0.994** / 0.98               | 29.8 T / **23.8 T** / 27 T           |

Milroy's own empirical fit, Eq. (17) `tau_P = lambda^2 / (2 sqrt(gamma_N))`,
gives 31 T and 22 T for these two cases, so both his numbers and ours sit ~10-30%
above his own fit; that fit spans lambda = 6.5 to 100 and is only claimed to be
"a reasonable approximation".

The `gamma = 18.2` case is the one that now overshoots (23.8 T vs 27 T), and
this is not a matter of where the threshold is drawn: alpha goes from 0.833 at
23 T to 0.986 at 25 T, so every reasonable criterion gives 24-25 T. What it
means is that our gamma dependence is slightly too steep --
t_pen(16.6)/t_pen(18.2) = 1.71 for us against 1.48 for Milroy and 1.44 from his
own Eq. (17). Departure 6 above (the harmonic truncation, whose error grows with
gamma) is the most likely cause; his 27 T is also a number read off a 1999
figure. Either way both cases moved substantially towards him, and neither is
worse than before.

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

**1. The `gamma` ramp is the whole timing discrepancy** — see item 1 above.
It costs about 6 T on both gamma = 16.6 and gamma = 18.2, which is the entire
47 vs 40 gap. It is visible in the old figure too: the alpha(t) curves are flat
for the first ~5 T before they start to climb, where Milroy's rise from t = 0.

**2. First-order boundary condition at r = 1.** The Robin condition was applied
as `A[-1] + (A[-1]-A[-2])/dr = 2 gamma e^{-i t}`, which is O(dr). The skin depth
is 1/lambda = 0.09, so at Nr = 64 that is an 8% error in the flux at the wall
against the analytic solution (table in item 2 above). Because the penetration time
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
