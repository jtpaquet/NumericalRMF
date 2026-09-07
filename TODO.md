# TODO

Working list from the audit of the solver against Milroy 1999 and the source
papers. `COMPARISON.md` has the evidence for each item.

## Done

- [x] **1. RF turn-on ramp.** `gamma(t) = gamma (1 - exp(-t/tau_r))` is real —
  Hugrass & Grimm 1981 Eq. (3) — but `tau_r = 3 T` was ~10x their value
  (`tau_r = 0.4 us` at `omega = 5e6 /s`, i.e. `omega tau_r = 2.0`, or 0.32 T).
  Milroy applies the RMF as a step at `t = 0`. Cost at `gamma = 16.6`: 6 T at
  `3 T`, 0.4 T at their value. Default is now 0; `RISE_TIME_HG1981` holds their
  value if it is ever wanted. It was *not* needed for stability — see the note
  at the bottom.

- [x] **2. Outer boundary condition.** The Robin condition was discretised with
  a backward difference, which applies the vacuum matching half a cell inside the
  plasma: `O(dr)`, and 8% wrong at `Nr = 64` because the skin depth is only
  `5.7 dr`. Now a ghost node at `r = 1+dr` with a centred derivative, verified
  against the analytic `gamma -> 0` solution (Hugrass 1985 Eq. 21).

- [x] **3. Axis of the b equation.** `2(b1-b0)/dr^2` was half the correct
  `4(b1-b0)/dr^2` (the `b'/r -> b''(0)` contribution was dropped), and the Hall
  source was missing entirely although `(1/r) dP/dr -> P''(0)` is finite. Recast
  as the finite-volume flux `G = r b' + P/(2 lambda^2)`, which is what Hugrass &
  Grimm describe ("A_z, B_z constant over the first radial zone, Gauss's
  theorem"). A factor-2 error in the axial Hall face (`P ~ r^2` there, so the
  arithmetic mean is 2x the true `P(dr/2) = P1/4`) was found and fixed while
  writing this up; it changes the answers by nothing measurable.

- [x] **4/9. Definitions.** `alpha_s` = the last value of a converged run,
  `t_pen` = when `alpha` first reaches 95% of `alpha_s`. Note the literature is
  not consistent here: Hugrass & Grimm use `alpha(tau) = 0.5 alpha_s`, Milroy
  says "full penetration" with no stated criterion.

- [x] **5a. gamma_c.** Milroy's Eq. (15) `gamma_c = 1.12 lambda (1 + 0.12
  (lambda-6.5)^0.4) = 15.13`, not `1.12 lambda = 12.40` (that is his Eq. 14
  expulsion threshold). Measured by the staircase procedure: between 15.25 and
  15.50.

- [x] **7a. Integrator check.** It is RK2 — Heun's method, the explicit
  trapezoidal predictor-corrector — and it is correct and second order. Nothing
  to fix.

- [x] **10. Electron inertia.** Correctly dropped: `omega/nu_ei ~ 2e-3` for
  Milroy's parameters. Hugrass & Grimm only retain it when `omega <~ nu_ei`.

## Next

- [ ] **5b. gamma_c vs lambda (Milroy Fig. 4).** Study 2 is written
  (`studies/study2_gamma_c_vs_lambda/`) but not yet run. Cost goes as lambda^4,
  so lambda <= 20 is comfortable, 25-40 is hours, and lambda >= 50 needs item 7b
  first. Run `--dry-run` for the plan before committing to a scan.

- [ ] **5c. Hysteresis loop (Milroy Fig. 9).** `study_thresholds.py` has the
  ascending branch saved in `results/alpha_s_up.txt` (gamma = 10 to 16); the
  descending branch died partway on a loaded machine. One clean run of both gives
  the expulsion threshold, expected near `1.12 lambda = 12.40` at
  lambda = 11.07.

- [ ] **7b. Semi-implicit time integration.** Crank-Nicolson on the linear
  diffusion (tridiagonal solve) with the Hall term explicit. Both Hugrass & Grimm
  and Milroy do this; H&G report running "about 10 times the Courant condition at
  the origin". Would lift the `dt = 0.002` cap, which is a stability limit, not
  an accuracy one — `dt = 0.002` already matches `dt = 0.001` to five digits.
  **Do this before item 6**, since more harmonics will make the explicit limit
  worse.

- [ ] **6. More azimuthal modes.** Stay at n = 0 / n = 1 for now — that is what
  Milroy uses in Sec. III.A. Later: `B_z` carries even harmonics and `A_z` odd
  ones, so the first thing dropped is `B_z2` feeding back on `A_z1` at second
  order, an error that grows with `gamma`. Hugrass & Grimm carried all harmonics
  on a 16x16 polar mesh; Milroy allows up to 42 modes but says n = 1 suffices
  here. This is the most likely explanation for the residual below.

## Open questions

- [ ] **`gamma = 18.2` penetrates ~2-3 T too fast** (24.4 T against Milroy's
  27 T) while `gamma = 16.6` lands on his 40 T. Not a threshold-definition
  artefact: `alpha` goes from 0.833 at 23 T to 0.986 at 25 T, so any criterion
  gives 24-25 T. Our gamma dependence is slightly too steep —
  `t_pen(16.6)/t_pen(18.2)` is 1.71 for us against 1.48 for Milroy and 1.44 from
  his own Eq. (17). Item 6 is the leading suspect; his 27 T is also a number read
  off a 1999 figure.

- [ ] **`alpha_s = 0.400` against Milroy's 0.42** at `gamma = 14.9`. Grid- and
  time-converged, so it is not numerics. `gamma/gamma_c = 0.985` puts it on the
  steepest part of the `alpha_s(gamma)` curve, where a 1-2% difference in
  `gamma_c` moves `alpha_s` by this much — our measured threshold is one 0.25
  step above his Eq. (15) fit.

- [ ] **Side-by-side with Milroy's Fig. 1** (field-line evolution). The field
  lines are plotted in the co-rotating frame now, but nobody has put the two
  figures next to each other at his four times.

## Not doing

- **8. Flux-conserving rings.** Hugrass & Grimm close the system with rings at
  `r = b`, so `B_z(a,t)` follows from axial flux conservation (their Eq. 17), and
  they show smaller `b` gives faster penetration. Milroy instead holds
  `B_z(R) = B_z0`, which is what this code does — and our device has no flux
  conserver, so Milroy's condition is also the physically right one for the
  eventual application.

## Housekeeping

- [x] Outputs moved out of the repo root into `results/` and `figures/`, routed
  through `paths.py`. The duplicate 21 MB pickle under `animation/` is gone.

- [x] Scripts grouped under `studies/study1_alpha_penetration/`, with
  `5_all_tests.py` added to run all four and overlay them.

- [ ] Fold the study 1 runs into `COMPARISON.md` now that they have been run at
  full length (Nr = 32, dt = 0.002, nothing unstable; fixed_ramp0 closest to
  Milroy).

---

Note on item 1: the ramp was originally introduced because the run seemed to be
unstable without it. Retested directly - the original scheme with `rise_time = 0`
runs fine at `Nr = 32` and `Nr = 64` for every `dt` that is stable *with* the
ramp, so the step start is not what was failing.

What *does* fail is the explicit scheme once the RMF has penetrated. At
`Nr = 64`, `dt = 0.004` survives the whole approach and then blows up at period
57 for `gamma = 16.6` and period 32 for `gamma = 18.2` - in both cases just after
`alpha` reaches 1 - while `gamma = 14.9`, which never penetrates, never blows up.
So the limit is set by the nonlinear Hall coupling at large `b`, not by the
start-up transient, and the ramp only ever postponed it by postponing
penetration. `scripts/*.py --dt-scan` reproduces this, and it is the argument for
item 7b.
