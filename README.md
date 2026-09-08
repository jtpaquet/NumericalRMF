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

An infinitely long plasma cylinder of radius $R$, surrounded by vacuum. Ions are
immobile with uniform density $n$ ($u = 0$), electrons are massless, resistivity
is uniform and isotropic, $\nabla P_e = 0$. Antennae impose an external rotating
field

$$
\mathbf{B}_{rmf} = B_w \cos(\omega t)\,\hat{x} + B_w \sin(\omega t)\,\hat{y},
\qquad \omega_{ci} \ll \omega \ll \omega_{ce}
$$

so the electrons are magnetised by the RMF while the ions are not, and an
azimuthal electron current is driven. With those assumptions the system is
Ohm's law plus Maxwell,

$$
\frac{\partial \mathbf{A}}{\partial t} = \mathbf{u}\times\mathbf{B} - \eta \mathbf{J}
- \frac{1}{en}\left(\mathbf{J}\times\mathbf{B} - \nabla P_e\right),
\qquad \mathbf{B} = \nabla\times\mathbf{A}, \qquad \mathbf{J} = \frac{\nabla\times\mathbf{B}}{\mu_0}
$$

### Normalisation

Milroy's scaling (his Sec. II.B): $t_0 = 1/\omega$, $r_0 = R$,
$A_0 = R B_w/\gamma$, $B_0 = A_0/R$, $J_0 = B_0/(\mu_0 R)$. Two dimensionless
numbers control everything:

$$
\lambda = \frac{R}{\delta} = R\sqrt{\frac{\mu_0 \omega}{2\eta}}
\qquad (\delta = \text{classical skin depth})
$$

$$
\gamma = \frac{\omega_{ce}}{\nu_{ei}} = \frac{1}{e}\frac{B_w}{n\eta}
\qquad (\text{RMF strength})
$$

and the equations collapse to Milroy Eq. (10),

$$
\frac{\partial A}{\partial t} = -\frac{1}{2\lambda^2}\left(J + \mathbf{J}\times\mathbf{B}\right)
$$

### Reduction to 1-D

Keeping the $n = 0$ component of $B_z$ and the $n = 1$ component of $A_z$ — the
truncation Milroy uses throughout Sec. III.A — and writing

$$
A = 2 A_{z1} \qquad \text{(Hugrass's convention: the BC becomes } A + A' = 2\gamma e^{-it}\text{)}
$$

$$
b = \frac{B_z}{\lambda^2} \qquad \text{(so } \alpha = b(1) - b(0) \text{ with no extra factors)}
$$

gives the closed pair actually integrated in `rmf_solver.py`:

$$
\frac{\partial A}{\partial t} = \frac{L[A]}{2\lambda^2} - \frac{i A b'}{2r}
$$

$$
\frac{\partial b}{\partial t} = \frac{1}{2\lambda^2 r}\frac{d}{dr}\left[r b' + \frac{P}{2\lambda^2}\right],
\qquad P = \Im\left(L[A]\,\bar{A}\right)
$$

$$
L[A] = A'' + \frac{A'}{r} - \frac{A}{r^2} = \frac{d}{dr}\left[\frac{1}{r}\frac{d}{dr}(rA)\right]
$$

Boundary and initial conditions:

$$
\begin{aligned}
A(0) &= 0 &&\text{n = 1 flux vanishes on axis} \\
A(1) + A'(1) &= 2\gamma e^{-it} &&\text{vacuum matching, Milroy Eq. (6)/(13)} \\
b'(0) &= 0 &&\text{symmetry} \\
b(1) &= B_{z0} = 0 &&\text{applied axial bias field} \\
A = b &= 0 \text{ at } t=0
\end{aligned}
$$

The Robin condition at $r = 1$ is the one that eliminates the unknown internal
screening currents: outside the plasma
$A_{z,n} = \alpha_n (r/R)^n + \beta_n (r/R)^{-n}$, and
$A_{z,n}(R) + \frac{R}{n}\frac{dA_{z,n}}{dr}\Big|_R = 2\alpha_n$ cancels
$\beta_n$, leaving only the known antenna drive
$\alpha_1 = \frac{R B_w}{2} e^{-i\omega t}$.

### Diagnostic

$$
\alpha = \frac{2}{\mu_0 n e \omega R^2}\left(B_z(R) - B_z(0)\right) = b(1) - b(0)
$$

the driven azimuthal current normalised to what it would be if every electron
rotated synchronously with the RMF. $\alpha_s$ is its steady value (the last
value of a converged run); $t_{pen}$ is when $\alpha$ first reaches 95% of
$\alpha_s$.

## Numerics

- Uniform radial mesh, $r_i = i\,dr$, $dr = 1/(N_r-1)$; axis included.
- $L[A]$ in the flux form $\frac{d}{dr}\left[\frac{1}{r}\frac{d}{dr}(rA)\right]$,
  second order and exact for both $A \sim r$ and $A \sim r^3$. The plain
  3-point form is only first order in the max norm because of the point next
  to the axis.
- Robin condition at $r = 1$ via a ghost node at $r = 1+dr$ with a centred
  derivative, and the boundary node evolved with the PDE. Second order.
- $b$ equation in finite-volume flux form $G = r b' + P/(2\lambda^2)$, which is
  conservative and handles the axial half-cell without a special limit — the
  treatment Hugrass & Grimm describe as "assuming $A_z$, $B_z$ are constant
  over the first radial zone and using Gauss's theorem".
- Time integration: explicit Heun (RK2 predictor-corrector), second order.
  Stability, not accuracy, sets the step: $dt = 0.002$ at $N_r = 64$ matches
  $dt = 0.001$ to five digits, $dt = 0.003$ blows up, and the limit tightens
  as $\gamma$ grows ($0.002$ fails near $\gamma = 20$).
- Production: $\lambda = 11.07$, $N_r = 64$, $dt = 0.002$, 200 RMF periods,
  ~30 s per $\gamma$.

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

### 1 — `studies/study1_alpha_penetration/` : $\alpha$ vs time, penetration tests

$\alpha(t)$ for $\gamma = 14.9, 16.6, 18.2$, in four configurations that differ
only in the discretisation and the turn-on ramp, so each effect can be isolated:

| script | discretisation | `rise_time` |
|---|---|---|
| `1_legacy_ramp0.py` | original (`nr32_fix.py`) | 0 |
| `2_legacy_ramp3.py` | original | 3 T |
| `3_fixed_ramp0.py` | corrected (`rmf_solver.py`) | 0 |
| `4_fixed_ramp3.py` | corrected | 3 T |
| `5_all_tests.py` | runs all four and overlays them, one panel per $\gamma$ | |

Each writes `results/alpha_<label>.pkl` and `figures/alpha_<label>.{pdf,png}`,
and prints $\alpha_s$ and $t_{pen}$ against Milroy's numbers. Common options
`--nr`, `--dt`, `--periods`, `--no-figure`, plus on the first four

```
python 1_legacy_ramp0.py --dt-scan 0.004 0.003 0.002 0.001 --periods 60
```

which reports, per $\gamma$, whether the run stayed finite and at which period
it blew up. The explicit scheme fails *after* the RMF has penetrated, not at
start-up, so a scan shorter than ~60 periods reports everything as stable.

### 2 — `studies/study2_gamma_c_vs_lambda/` : Milroy Fig. 4

`gamma_c_scan.py` finds the critical $\gamma$ for full penetration at each
$\lambda$ by bisection on cold starts (corrected scheme, `rise_time = 0`), and
`plot_gamma_c.py` plots it against Milroy's Eqs. (14) and (15).

Resolution and time step are chosen per $\lambda$ ($N_r \sim 8\lambda$ to
resolve the skin depth, $dt$ from the Hall/whistler limit $\sim dr^2$), and the
time budget is `--tmax-factor` $\times\ \lambda^2$ RMF periods, which makes the
finite-budget bias in $\gamma_c$ independent of $\lambda$ — about 0.6% high at
the default.

The cost goes as $\lambda^4$ ($N_r \sim \lambda$, $dt \sim dr^2$,
$T_{max} \sim \lambda^2$), so start with

```
python gamma_c_scan.py --dry-run
```

which prints the plan and the step count per $\lambda$ without running
anything. $\lambda \le 20$ is comfortable, 25-40 runs into hours, and
$\lambda \ge 50$ is out of reach with the explicit integrator — that is what
the semi-implicit scheme (TODO 7b) would buy. Results are appended to
`results/gamma_c_vs_lambda.csv` as they are produced and `--resume` skips
$\lambda$s already there, so the scan can be stopped and restarted:

```
python gamma_c_scan.py --lam-max 20
python plot_gamma_c.py
```

### 3 — `studies/study3_field_expulsion_params/` : field expulsion vs. fill pressure and $T_e$

Maps the pre-ionization state -- argon at 25 C, 10-100 mTorr fill, giving
$n_e$ via the ideal-gas density of the fill (full single ionisation assumed)
-- and an independently scanned $T_e = 0.5\text{-}10$ eV onto $(\lambda, \gamma)$
through Spitzer resistivity (`device.py`), using this device's own numbers:
$f_{RMF} = 250$ kHz, $R = 10$ cm, $B_w = 6$ G vacuum antenna field (cross-checked
against $0.163\ \text{G/A} \times 35\ \text{A} = 5.7$ G from the coil calibration).

$\gamma/\gamma_c$ stays below 2% everywhere in the box (`figures/study3_gamma_ratio.pdf`)
-- at $B_w = 6$ G this device never leaves the linear (classical skin-effect)
regime, so `1_scan_grid.py` runs the real solver once per $T_e$ at a small
reference $\gamma$ and checks $\gamma$-independence directly against the
largest $\gamma$ actually reached in the grid (0.02% difference).

|                                     | model                                    | observed          |
|---|---|---|
| field-settling time                | 0.008-0.03 ms, set by $T_e$ alone        | 10-25 ms          |
| steady wall - centre $\lvert B_r \rvert$ split | matches measurement at $T_e \sim 3\text{-}4$ eV    | (1.6-1.0)/6 = 0.11 |
| steady wall $\lvert B_r \rvert$ alone          | matches measurement at $T_e \sim 0.7\text{-}0.8$ eV | 1.6/6 = 0.27      |

```
python 1_scan_grid.py    # the solver runs + results/study3_grid.csv
python 2_plot_grid.py    # gamma/gamma_c, settling time, screening heatmaps
python 3_time_traces.py  # |B_r(R)-B_r(0)| vs real time, against 10/25 ms
```

Classical Spitzer-resistivity diffusion reproduces the *steady* wall and
centre fields, but only at two different, mutually inconsistent values of
$T_e$ -- and it misses the *timescale* by ~1000x everywhere in the box,
because the pressure (density) barely matters here: with $\gamma$ this small,
only $T_e$ (through $\eta$ and $\lambda$) sets the field's own diffusion time,
and that time is always sub-millisecond for $T_e = 0.5\text{-}10$ eV. The
observed 25 ms field-expulsion time is therefore more likely tracking the
post-breakdown density (ionisation) build-up itself than classical field
diffusion into an already-formed plasma; reproducing it would need $n_e(t)$,
or a cross-field (magnetised) resistivity that this isotropic-$\eta$ model
does not have.

**Ionisation-fraction extension.** Part 1 above assumed full single ionisation
($n_e = n_{gas}$); `4_ionization_fraction.py` instead scans the actual unknown
-- the fraction of the fill that a stochastic, cosmic-ray-seeded avalanche has
ionised, $f_{ion} = n_e/n_{gas} \in \{10^{-6}, \ldots, 10^{-1}, 1\}$ -- since a
lower $f_{ion}$ raises $\gamma$ and can cross $\gamma_c$, where Milroy's
Eq. (17) penetration time *diverges*: in principle a slow, pressure-dependent
mechanism, unlike part 1's threshold-free classical diffusion. Eq. (17) is
validated against the real solver at this study's own $\lambda$s first
(previously only checked at Milroy's $\lambda = 11.07$): matches to 8-10%, the
same residual already noted above for Milroy's own case.

```
python 4_ionization_fraction.py   # Eq. (17) validation + the f_ion grid/window
python 5_plot_ionization.py       # t_pen vs f_ion; critical f_ion(p) per T_e
```

None of the 7 requested $f_{ion}$ land in the observed 1-50 ms window across
the whole 5x5 (pressure, $T_e$) grid (`figures/study3_ionization_tau.pdf`) --
$f_{ion}$ is either deep subcritical (never penetrates; classical diffusion
only, part 1's sub-ms result) or deep supercritical (penetrates in << 1 ms).
What *would* work is a narrow band of $f_{ion}$ straddling $\gamma_c$ where
Eq. (17)'s divergence supplies 1-50 ms, but that band is razor-thin relative
to the critical value itself (`figures/study3_ionization_window.pdf`): 0.03%
of $f_{crit}$ at $T_e = 0.5$ eV, widening to ~50% only by $T_e = 10$ eV.
Landing in it at every pressure, as the data requires, needs either $T_e$ on
the high end of the stated range (eV, where the band is an order-1 fraction
of $f_{crit}$) or a mechanism that pins $f_{ion}$ near threshold rather than
letting a stochastic avalanche land wherever it lands. Combined with part 1's
~1000x timescale
miss in the sub-critical regime, both mechanisms this model can produce point
the same way: **the 1-50 ms coupling time this device sees is most likely the
avalanche's own statistical growth time**, not RMF field diffusion (linear or
nonlinear) into an already-formed plasma.

**Real shots.** `data/env_diff/*.csv` (30 shots, provenance in
`6_load_experimental_data.py`) are bdot-probe vs. Rogowski-inferred envelope
pairs at the chamber centre: `bdot_x/z_env` is the actually-measured field,
`B_rog_x/z_env` is what that field would be from the driven antenna current
alone (a shot-by-shot, dynamically measured $B_w$ -- it tracks the antenna
loading down once the plasma couples, unlike the fixed 6 G no-plasma value
used above). Their ratio is exactly the model's `br0_over_bw(lambda)`
(the $\gamma \to 0$ limit, `device.py`) at the plasma's real, unknown $T_e$.

```
python 6_load_experimental_data.py   # per-shot ratio, coupling time, inferred Te
python 7_plot_experimental.py        # inferred Te and coupling time vs p and ohm
```

The 90%-of-plateau coupling time rises from ~1.6 ms at 7-22 mTorr to ~10-18 ms
by 60-99 mTorr (`figures/study3_experimental_time.pdf`) -- the same pressure
trend the 1-50 ms estimate was based on, now quantified shot by shot. Inverting
the steady-state ratio through `br0_over_bw` gives $T_e = 0.03\text{-}0.2$ eV
(`figures/study3_experimental_Te.pdf`) -- **an order of magnitude colder**
than the 0.5-10 eV range parts 1-2 scanned, rising with both pressure and
added DC-coil resistance (a weaker bias field). Take this $T_e$ as an
effective, single-point-model number rather than a literal measurement: the
added-resistance dependence alone (weak bias field couples fast and deep,
strong bias field couples slow and partial, and the effect disappears when
the bias field is reversed) is a first-order effect this radial, bias-field-
free model cannot represent at all, so some of what gets read as "$T_e$" here
is really that missing physics.

### 4 — `studies/study4_dc_bias_field/` : does an external $B_{z0}$ explain the reversal asymmetry?

Study 3's real shots showed a first-order effect its linear-regime model couldn't
address at all: weaker DC bias field (higher added-resistance) coupled fast and
deep, stronger bias field (0-2 Ω) coupled slow and only partially (sometimes
relaxing back toward vacuum mid-shot), and the whole effect vanished when the
bias polarity was reversed -- classic FRC field-reversal phenomenology. The
natural next step was to stop hardcoding $b(1) = B_{z0} = 0$ (TODO item 8) and
see whether a nonzero boundary value reproduces it.

It doesn't, and not for a numerical reason. `rmf_solver.py` now takes a `bz0`
argument and applies it at the boundary, but every term that touches $b$ in
`rhs()` -- the Hall-coupling $b'$ in $\partial A/\partial t$, and the
differences that make up the whole of $\partial b/\partial t$ -- uses a
*gradient* of $b$, never $b$ itself. So $b(1) = b_{z0}$ is a pure relabelling:
the solution shifts uniformly by $b_{z0}$ and $\alpha = |b(1) - b(0)|$ cannot
depend on it, for any $\lambda$, $\gamma$, or time. `1_bz0_invariance.py`
checks this against the actual solver rather than just the algebra: at
Milroy's $\lambda = 11.07$ and his three reference $\gamma$s (one
sub-critical, two super-critical), $\alpha(t)$ at $b_{z0} \in
\{0.05, 0.3, -0.3, 1, -1\}$ matches $b_{z0} = 0$ to float round-off
($\sim 10^{-14}$) at every point of a 60-period run
(`figures/study4_bz0_invariance.pdf`).

```
python 1_bz0_invariance.py
```

There is a related but genuine asymmetry buried in the equations: conjugating
$A$ and negating $b$ maps a solution at $(b_{z0}, \gamma e^{-i\tau})$ to one at
$(-b_{z0}, \gamma e^{+i\tau})$, so flipping the bias field *without* also
reversing the RMF's rotation sense is not a symmetry of this system in
general -- but it's moot here, since $b_{z0}$ already drops out of the
dynamics entirely regardless of sign.

**This makes the study 3 caveat more precise, not just confirmed.** It isn't
that this radial model *probably* lacks the physics; wiring $B_{z0}$ into the
existing boundary condition literally cannot change this model's output, by
construction. Real FRC field reversal is a statement about magnetic-null
topology along the axis -- something an infinite, translationally-invariant
cylinder with only the $n=0$ (axial, $r$-only) and $n=1$ (transverse) harmonics
structurally cannot represent, independent of what boundary value is chosen.
Testing the reversal-asymmetric-coupling hypothesis for real would need actual
axial structure (finite length, end effects) or a force-balance equation where
the *absolute* value of $B_z$ enters, not just its gradient -- both well beyond
the fixed-ion, Ohm's-law-only reduction used throughout this project. Left as
future work; TODO.md is updated accordingly.

## Where it stands

$\lambda = 11.07$, $N_r = 64$, both schemes run to 200 T so the comparison is
like for like:

| $\gamma$ | $\gamma/\gamma_c$ | $\alpha_s$ (orig / now / Milroy) | $t_{pen}$ (orig / now / Milroy) |
|---|---|---|---|
| 14.9 | 0.985 | 0.387 / **0.400** / 0.42 | never |
| 16.6 | 1.097 | 0.990 / **0.992** / 0.98 | 46.8 T / **40.9 T** / 40 T |
| 18.2 | 1.203 | 0.992 / **0.994** / 0.98 | 29.8 T / **24.4 T** / 27 T |

($t_{pen}$ at 95% of $\alpha_s$. Figure: `figures/alpha_vs_time_comparison.pdf`.)

Note that $\gamma_c = 15.13$ here is Milroy's Eq. (15) penetration threshold,
not $1.12\lambda = 12.40$, which is his Eq. (14) value for *expulsion* of an
already penetrated RMF. On that scale $\gamma = 14.9$ sits 1.5% *below*
threshold, which is why it never penetrates.

## Validation

- **Analytic linear limit.** As $\gamma \to 0$ the Hall term drops and the
  steady $n=1$ flux is $A = 2\gamma I_1(kr)/(k I_0(k))$ with
  $k = (1-i)\lambda$ (Hugrass 1985 Eqs. 20-21). The solver matches it to 0.04%
  at $N_r = 64$ and converges at second order.
- **Threshold.** Milroy's staircase procedure puts the penetration threshold
  between $\gamma = 15.25$ and $15.50$, against 15.13 from his Eq. (15).
- **$\alpha_s(\gamma)$.** The ascending sweep (`results/alpha_s_up.txt`) tracks
  his Eq. (18) fit within the ~20% he claims for it.
- **Grid convergence.** $\alpha_s$ at $\gamma = 14.9$ is converged to four
  digits by $N_r = 32$.

## References

- R. D. Milroy, Phys. Plasmas **6**, 2771 (1999) — the reference for this project
- W. N. Hugrass and R. C. Grimm, J. Plasma Phys. **26**, 455 (1981)
- I. R. Jones and W. N. Hugrass, J. Plasma Phys. **26**, 441 (1981)
- W. N. Hugrass, Aust. J. Phys. **38**, 157 (1985)
- R. D. Milroy, Phys. Plasmas **7**, 4135 (2000) — 2-D MHD successor, moving ions
