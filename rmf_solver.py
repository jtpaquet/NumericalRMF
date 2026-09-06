"""Rotating-magnetic-field penetration into a plasma column (fixed-ion model).

Reproduces the calculations of Sec. III.A of

    R. D. Milroy, "A numerical study of rotating magnetic fields as a current
    drive for field reversed configurations", Phys. Plasmas 6, 2771 (1999)

which in turn follows Hugrass & Grimm, J. Plasma Phys. 26, 455 (1981).

Model
-----
Infinite cylinder, fixed ions (u = 0), massless electrons, uniform n, no
electron-pressure gradient.  In Milroy's dimensionless variables (t0 = 1/w,
r0 = R, A0 = R*B_w/gamma, B0 = A0/R) Eqs. (10)-(13) read

    dA/dt = -(1/2 lam^2) ( J + J x B ),   B = curl A,   J = curl B

with the vacuum-matching boundary condition on the n = 1 component of A_z

    A_z1 + dA_z1/dr = gamma e^{-i t}      at r = 1.

Keeping only the n = 0 component of B_z and the n = 1 component of A_z (the
truncation Milroy uses for all of Sec. III.A) and writing

    A = 2 A_z1        (so the boundary condition becomes A + A' = 2 gamma e^{-i t})
    b = B_z / lam^2   (so alpha = b(1) - b(0) directly)

gives the closed pair actually integrated here

    dA/dt = L[A] / (2 lam^2)  -  i A b' / (2 r)
    db/dt = (1 / (2 lam^2 r)) d/dr [ r b' + P / (2 lam^2) ],   P = Im( L[A] conj(A) )

    L[A] = A'' + A'/r - A/r^2 = d/dr [ (1/r) d/dr (r A) ]

with A(0) = 0, b'(0) = 0 and b(1) = B_z0 = 0.

Numerics
--------
The defaults are second-order accurate in r; ``legacy=True`` restores the
first-order choices of the original ``nr32_fix.py`` so the two can be compared
(see ``validate_numerics.py``).  The three differences are

  * ``A_op``    : ``'flux'`` writes L[A] as d/dr[(1/r) d/dr (rA)] and differences
                  it conservatively.  The plain 3-point form ``'naive'`` is only
                  *first* order in the max norm because of the point next to the
                  axis (it is exact for A ~ r but not for A ~ r^3).
  * ``bc_order``: the Robin condition at r = 1 is applied with a ghost node and a
                  centred derivative (2nd order) instead of a one-sided
                  difference (1st order).  This matters more than its formal
                  order suggests: the RMF skin depth is 1/lam ~ 0.09, so a first
                  order boundary error is a percent-level error in the *drive*,
                  and the penetration time diverges as (gamma - gamma_c)^(-1/2).
  * ``b_scheme``: finite-volume form of the b equation.  The legacy version uses
                  2 (b1-b0)/dr^2 on axis where the correct axial limit of the
                  cylindrical Laplacian is 4 (b1-b0)/dr^2, and drops the Hall
                  source there entirely although (1/r) dP/dr -> P''(0) is finite.
                  alpha is measured at r = 0, so both errors bias it low.

``rise_time`` switches on the RF turn-on ramp gamma(t) = gamma (1 - exp(-t/tau_r))
of Hugrass & Grimm 1981 Eq. (3).  It is a real feature of their model -- tau_r is
"characteristic of the RF source" -- but their value is small: tau_r = 0.4 us at
omega = 5e6 /s, i.e. omega*tau_r = 2.0, or 0.32 of an RMF period (RISE_TIME_HG1981
below).  Milroy 1999 drops it and applies the RMF as a step at t = 0, so the
default here is 0; at the Hugrass & Grimm value the penetration time moves by
0.4 T, while at 3 T it moves by 6 T.  See ``validate_numerics.ramp()``.
"""
import numpy as np

TWO_PI = 2.0*np.pi

# Hugrass & Grimm 1981: tau_r = 0.4 us at omega = 5e6 /s, in units of 1/omega
RISE_TIME_HG1981 = 2.0


def gamma_c(lam):
    """Critical gamma for full penetration, Milroy 1999 Eqs. (14)-(15).

    Note this is *not* 1.12*lam for lam > 6.5; 1.12*lam is the (lower) value at
    which an already-penetrated RMF is expelled again.
    """
    if lam <= 6.5:
        return 1.12*lam
    return 1.12*lam*(1.0 + 0.12*(lam - 6.5)**0.4)


def tau_penetration(lam, gam):
    """Empirical penetration time, Milroy 1999 Eq. (17), in RMF periods."""
    gN = (gam - gamma_c(lam))/gamma_c(lam)
    if gN <= 0:
        return np.inf
    return lam**2/(2.0*np.sqrt(gN))/TWO_PI


class RMFPenetration:
    def __init__(self, Nr=64, lam=11.07, gam=16.6, rise_time=0.0,
                 legacy=False, A_op=None, bc_order=None, b_scheme=None):
        self.Nr, self.lam, self.gam = Nr, lam, gam
        self.lam2, self.lam4 = lam**2, lam**4
        self.rise_time = rise_time
        self.A_op = A_op if A_op else ('naive' if legacy else 'flux')
        self.bc_order = bc_order if bc_order else (1 if legacy else 2)
        self.b_scheme = b_scheme if b_scheme else ('legacy' if legacy else 'fv')

        self.dr = 1.0/(Nr - 1)
        self.r = np.linspace(0.0, 1.0, Nr)
        self.rh = self.r[:-1] + 0.5*self.dr          # r_{i+1/2}
        self.V = self.r*self.dr                      # finite-volume weights
        self.V[0] = self.dr**2/8.0                   # axial half cell
        self.A = np.zeros(Nr, dtype=complex)
        self.b = np.zeros(Nr)
        self.tau = 0.0

    # ------------------------------------------------------------------ drive
    def gamma_t(self, tau):
        if self.rise_time <= 0.0:
            return self.gam
        return self.gam*(1.0 - np.exp(-tau/self.rise_time))

    def drive(self, tau):
        return self.gamma_t(tau)*np.exp(-1j*tau)

    def _ghost(self, A, tau):
        """A at r = 1 + dr, from A + A' = 2 gamma e^{-i t} with A' centred at r = 1."""
        return A[-2] + 2.0*self.dr*(2.0*self.drive(tau) - A[-1])

    # -------------------------------------------------------------- operators
    def lapA(self, A, tau):
        r, dr = self.r, self.dr
        L = np.zeros(self.Nr, dtype=complex)
        if self.A_op == 'naive':
            L[1:-1] = ((A[2:] - 2*A[1:-1] + A[:-2])/dr**2
                       + (A[2:] - A[:-2])/(2*dr)/r[1:-1] - A[1:-1]/r[1:-1]**2)
            if self.bc_order == 2:
                Ag = self._ghost(A, tau)
                L[-1] = ((Ag - 2*A[-1] + A[-2])/dr**2
                         + (Ag - A[-2])/(2*dr)/r[-1] - A[-1]/r[-1]**2)
        else:
            rA = r*A
            w = (rA[1:] - rA[:-1])/(self.rh*dr)      # (1/r) d(rA)/dr at r_{i+1/2}
            L[1:-1] = (w[1:] - w[:-1])/dr
            if self.bc_order == 2:
                Ag = self._ghost(A, tau)
                wN = ((1.0 + dr)*Ag - rA[-1])/((1.0 + 0.5*dr)*dr)
                L[-1] = (wN - w[-1])/dr
        return L

    def dbdr(self, b):
        dr = self.dr
        g = np.zeros_like(b)
        g[1:-1] = (b[2:] - b[:-2])/(2*dr)
        g[-1] = (3*b[-1] - 4*b[-2] + b[-3])/(2*dr)   # 2nd-order one-sided
        return g

    def rhs(self, A, b, tau):
        r, dr, N = self.r, self.dr, self.Nr
        L = self.lapA(A, tau)
        P = np.imag(L*np.conj(A))

        rA = np.zeros(N, dtype=complex)
        s = slice(1, N if self.bc_order == 2 else N - 1)
        rA[s] = L[s]/(2*self.lam2) - 1j*A[s]*self.dbdr(b)[s]/(2*r[s])

        rb = np.zeros(N)
        if self.b_scheme == 'fv':
            # flux G = r b' + P/(2 lam^2) ;  db/dt = (1/(2 lam^2)) div G
            G = self.rh*(b[1:] - b[:-1])/dr + 0.5*(P[1:] + P[:-1])/(2*self.lam2)
            rb[0] = G[0]/(2*self.lam2*self.V[0])     # G(-1/2) = 0 by symmetry
            rb[1:-1] = (G[1:] - G[:-1])/(2*self.lam2*self.V[1:-1])
        else:
            rb[1:-1] = ((b[2:] - 2*b[1:-1] + b[:-2])/dr**2
                        + (b[2:] - b[:-2])/(2*dr)/r[1:-1])/(2*self.lam2) \
                       + (P[2:] - P[:-2])/(2*dr)/(4*self.lam4*r[1:-1])
            rb[0] = 2*(b[1] - b[0])/dr**2/(2*self.lam2)
        return rA, rb

    def apply_bc(self, A, b, tau):
        A[0] = 0.0
        if self.bc_order == 1:
            A[-1] = (2*self.drive(tau) + A[-2]/self.dr)/(1.0 + 1.0/self.dr)
        b[-1] = 0.0                                   # B_z(R) = B_z0 = 0
        return A, b

    # ---------------------------------------------------------- time stepping
    def step(self, dt):
        """Heun (explicit predictor-corrector), as in Milroy's first option."""
        rA1, rb1 = self.rhs(self.A, self.b, self.tau)
        As, bs = self.A + dt*rA1, self.b + dt*rb1
        As, bs = self.apply_bc(As, bs, self.tau + dt)
        rA2, rb2 = self.rhs(As, bs, self.tau + dt)
        self.A += 0.5*dt*(rA1 + rA2)
        self.b += 0.5*dt*(rb1 + rb2)
        self.tau += dt
        self.A, self.b = self.apply_bc(self.A, self.b, self.tau)

    # ------------------------------------------------------------ diagnostics
    def alpha(self):
        """Driven current normalised to synchronous rotation (Milroy Sec. III.A)."""
        return abs(self.b[-1] - self.b[0])

    def B_theta_wall(self):
        """|B_theta(R)|/B_omega, taken from the boundary condition itself.

        A'(1) = 2 gamma e^{-i t} - A(1) exactly, so no one-sided difference is
        needed for the quantity Milroy plots in his Fig. 10.
        """
        return abs(2*self.drive(self.tau) - self.A[-1])/self.gam

    def profiles(self):
        """(r, |B_r|/B_w, |B_theta|/B_w, B_z/B_w) for the n=1 transverse field.

        B_z is the n=0 axial field; note B_z/B_w = lam^2 * b / gamma - the lam^2
        is easy to drop, and without it the axial field looks ~100x too small.
        """
        r, dr = self.r, self.dr
        dA = np.empty_like(self.A)
        dA[1:-1] = (self.A[2:] - self.A[:-2])/(2*dr)
        dA[0] = (-3*self.A[0] + 4*self.A[1] - self.A[2])/(2*dr)
        dA[-1] = 2*self.drive(self.tau) - self.A[-1]          # exact, from the BC
        with np.errstate(divide='ignore', invalid='ignore'):
            Br = np.abs(self.A)/r
        Br[0] = abs(dA[0])                                    # |B_r| = |B_theta| on axis
        return r, Br/self.gam, np.abs(dA)/self.gam, self.lam2*self.b/self.gam

    def psi(self, r_plot, theta_plot):
        """Flux function in the frame co-rotating with the RMF (Milroy's Figs. 1, 3).

        The physical flux is Re[A e^{i theta}] and A carries the e^{-i t} of the
        drive, so the RMF is held still by evaluating Re[A e^{i(theta + t)}].
        """
        Ai = (np.interp(r_plot, self.r, self.A.real)
              + 1j*np.interp(r_plot, self.r, self.A.imag))
        return np.real(Ai[None, :]*np.exp(1j*(theta_plot[:, None] + self.tau)))

    # ------------------------------------------------------------------- runs
    def run(self, n_periods=101, dt=0.002, snapshot_periods=()):
        nstep = int(round(TWO_PI/dt))
        snaps = {0: (self.A.copy(), self.b.copy(), self.r.copy(), self.tau)}
        ts, als = [], []
        want = set(snapshot_periods)
        for p in range(n_periods):
            for _ in range(nstep):
                self.step(dt)
            ts.append(p + 1.0)
            als.append(self.alpha())
            if (p + 1) in want:
                snaps[p + 1] = (self.A.copy(), self.b.copy(), self.r.copy(), self.tau)
        return np.array(ts), np.array(als), snaps
